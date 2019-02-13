/** ==========================================================================
#    This file is part of the finite element software ParMooN.
# 
#    ParMooN (cmg.cds.iisc.ac.in/parmoon) is a free finite element software  
#    developed by the research groups of Prof. Sashikumaar Ganesan (IISc, Bangalore),
#    Prof. Volker John (WIAS Berlin) and Prof. Gunar Matthies (TU-Dresden):
#
#    ParMooN is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    ParMooN is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
#
#    If your company is selling a software using ParMooN, please consider 
#    the option to obtain a commercial license for a fee. Please send 
#    corresponding requests to sashi@iisc.ac.in

# =========================================================================*/ 
   
// =======================================================================
// @(#)SquareStructure2D.C        1.6 09/17/99
//
// Class:       TSquareStructure2D
//
// Purpose:     build and store a structure for a square matrix in 2d
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <DefineParams.h>
#include <FEDatabase2D.h>
#include <SquareStructure2D.h>
#include <Database.h>
#include <string.h>

#include <Joint.h>

#ifdef __MORTAR__
#ifdef __ADD_LINK__
#include <Database.h>
#include <It_Mortar.h>
#include <MortarJoint.h>
#include <stdlib.h>
#endif                                            // __ADD_LINK__
#endif                                            // __MORTAR__

#include <Database.h>
/** dummy constructor, needed only derives classes */
TSquareStructure2D::TSquareStructure2D()
: TSquareStructure()
{
}


/** generate the matrix structure, both space are 2D */
TSquareStructure2D::TSquareStructure2D(TFESpace2D *Space)
: TSquareStructure()
{
  int i,j,k,l,n,N_, n1,n2,m,p,q,r;
  int *GlobalNumbers;
  int *BeginIndex;
  int *Numbers, *NumbersNeighbour;

  int N_Dirichlet;
  int N_Inner;
  int N_Hanging;
  int *BoundNodeBounds, N_NonDiri, N_BoundaryNodeTypes;
  int HangingBound;

  int *AuxPtr, *HangingAuxPtr;
  int *KColAux, *HangingKColAux;
  int index, oldindex;

  int Offset, *DOF, end, begin;
  THangingNode **HangingNodes;
  THangingNode *hn;

  TCollection *Coll;

  TBaseCell *cell, *neigh;
  FE2D CurrentElement, CurrentNeighbour;

  FESpace = Space;

  ActiveBound=Space->GetN_ActiveDegrees();
  HangingBound=Space->GetHangingBound();

  N_BoundaryNodeTypes=Space->GetN_DiffBoundaryNodeTypes();
  BoundNodeBounds=Space->GetN_BoundaryNodes();
  N_NonDiri = 0;
  for(i=0;i<N_BoundaryNodeTypes;i++)
    N_NonDiri += BoundNodeBounds[i];

  N_Dirichlet=Space->GetN_Dirichlet();
  N_Inner=Space->GetN_Inner();
  N_Hanging=Space->GetN_Hanging();

  N_Rows = ActiveBound+N_Hanging+N_Dirichlet;
  N_Columns = N_Rows;
  // assembles matrices without shorter rows for non-active dof
   if (TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
  {
    N_NonDiri += N_Dirichlet;
    N_Dirichlet = 0;
    ActiveBound = N_Inner + N_NonDiri;
  }
  ColOrder = 0;
  // AuxPtr[i] will contain an upper bound for the number of
  // matrix entries in row i
  l=N_Rows+1;
  AuxPtr=new int[l];
  memset(AuxPtr, 0, l*sizeof(int));

  l=N_Hanging+1;
  HangingAuxPtr=new int[l];
  memset(HangingAuxPtr, 0, l*sizeof(int));

  GlobalNumbers=Space->GetGlobalNumbers();

  BeginIndex=Space->GetBeginIndex();

  Offset=ActiveBound;

  // loop over all elements
  N_=Space->GetN_Cells();
  Coll = FESpace->GetCollection();

  // associate each cell with her number in the collection
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);                        // set clipboard to number of the cell in collection
  }

  // loop over the mesh cells
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = Space->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    // loop over the local degrees of freedom
    for(j=0;j<n1;j++)
    {
      k=Numbers[j];
      if(k<ActiveBound)
      {
        AuxPtr[k]+=n1;

        // additional entries necessary if integrals over the edges
        // appear in the discretization
        if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
        {
          l=cell->GetN_Edges();                   // # edges
          for(m=0;m<l;m++)                        // for all edges
          {
                                                  // neighbour cell
            neigh = cell->GetJoint(m)->GetNeighbour(cell);
            if(neigh)
            {
              n = neigh->GetClipBoard();
              CurrentNeighbour = Space->GetFE2D(n, neigh);
              n2 = TFEDatabase2D::GetFE2D(CurrentNeighbour)->GetN_DOF();
              AuxPtr[k]+=n2;
            }                                     //endif
          }                                       //endfor m
        }
      }
      else
      {
                                                  // not adjusted for edge stabilization because hanging nodes are no more used
        if(k<HangingBound) HangingAuxPtr[k-Offset]+=n1;
      }                                           // endif
    }                                             // endfor j
  }                                               // endfor i

  // add rows for Dirichlet nodes in  space
  N_Entries=N_Dirichlet;
  Offset=ActiveBound+N_Hanging;
  for(i=0,j=Offset;i<N_Dirichlet;i++,j++)
  {
    AuxPtr[j]=1;
  }

  // add couplings for hanging nodes of  space
  HangingNodes=Space->GetHangingNodes();
  Offset=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    // cout << hn;
    // there is the additional entry in diagonal
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes() + 1;
    AuxPtr[j]=n;
    N_Entries+=n;
    // cout << "AuxPtr[j]: " << AuxPtr[j] << endl;
  }

  // additional space for storing the columns caused by the
  // hanging nodes of  space => some new columns
  Offset=ActiveBound;
  HangingNodes=Space->GetHangingNodes();
  // cout << "N_Hanging: " << N_Hanging << endl;
  for(i=0;i<N_Hanging;i++)
  {
    hn=HangingNodes[i];
    m=HangingAuxPtr[i];
    // cout << "HangingAuxPtr[i]: " << m << endl;
    DOF=hn->GetDOF();
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    for(j=0;j<n;j++)
    {
      k=DOF[j];
      if(k<ActiveBound)
        AuxPtr[k] += m;
    }
  }

#ifdef __MORTAR__
#ifdef __ADD_LINK__
  // reserve space for additional link term
  TIt_Mortar *It1 = (TIt_Mortar *) TDatabase::IteratorDB[It_Mortar1];
  TIt_Mortar *It2 = (TIt_Mortar *) TDatabase::IteratorDB[It_Mortar2];
  TBaseCell *CellM, *CellNM;
  FE2D CurrElementID;
  TFE2D *CurrElement;
  TFEDesc2D *CurrDesc;
  int lev, infoNM, infoM;
  double lamM0, lamM1 = 0, lamNM0, lamNM1, startX, startY, endX, endY;
  double delX, delY;
  bool NewMortarSideEle;
  int **J_DOF_NM, N_J_DOF_NM, **J_DOF_M, N_J_DOF_M;
  int indexM = 0, indexNM;

  // set ClipBoard to number in collection
  N_ = Space2D->GetN_Cells();
  Coll = Space2D->GetCollection();
  for(i=0;i<N_;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  // loop over all mortar edges
  N_ = It1->GetN_MortarFace();
  for (i=0;i<N_;i++)
  {
    NewMortarSideEle = true;
    lev = MAX_ItLevel + (i << 8);
    It1->Init(lev);
    It2->Init(-lev);

    It1->GetPoint(startX, startY);
    It2->GetPoint(endX, endY);
    delX = endX - startX;
    delY = endY - startY;

    while (CellNM = It1->Next(infoNM))
    {
      // get a new element on non-mortar side
      lamNM0 = GetLambda(startX, startY, CellNM->GetVertex(infoNM),
        delX, delY);
      lamNM1 = GetLambda(startX, startY, CellNM->GetVertex((infoNM+1)
        % CellNM->GetN_Vertices()), delX, delY);
      if (CellNM->GetClipBoard() == -1)
      {
        cerr << "Error in SquareStructure2D: cell out of collection!!!"
          << endl;
        exit(-3);
      }

      indexNM = CellNM->GetClipBoard();
      CurrElementID = Space2D->GetFE2D(indexNM, CellNM);
      CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
      CurrDesc = CurrElement->GetFEDesc2D();

      J_DOF_NM = CurrDesc->GetJointDOF();
      N_J_DOF_NM = CurrDesc->GetN_JointDOF();

      if (!NewMortarSideEle)
      {
        // add space for DOFs if the next while clause wont be
        // entered
        Numbers = GlobalNumbers + BeginIndex[indexNM];
        for (m=0;m<N_J_DOF_NM;m++)
          AuxPtr[Numbers[J_DOF_NM[infoNM][m]]] += N_J_DOF_M;

        Numbers = GlobalNumbers + BeginIndex[indexM];
        for (m=0;m<N_J_DOF_M;m++)
          AuxPtr[Numbers[J_DOF_M[infoM][m]]] += N_J_DOF_NM;
      }

      // check which side is next
      if (lamM1 < lamNM1) NewMortarSideEle = true;

      while (NewMortarSideEle)
      {
        // get a new element on mortar side
        if (!(CellM = It2->Next(infoM))) break;
        lamM0 = GetLambda(startX, startY, CellM->GetVertex((infoM+1)
          % CellM->GetN_Vertices()), delX, delY);
        lamM1 = GetLambda(startX, startY, CellM->GetVertex(infoM),
          delX, delY);
        if (CellM->GetClipBoard() == -1)
        {
          cerr << "Error in MatrixStructure: cell out of collection!!!"
            << endl;
          exit(-3);
        }

        indexM = CellM->GetClipBoard();
        CurrElementID = Space2D->GetFE2D(indexM, CellM);
        CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc2D();

        J_DOF_M = CurrDesc->GetJointDOF();
        N_J_DOF_M = CurrDesc->GetN_JointDOF();

        // add space for DOFs
        Numbers = GlobalNumbers + BeginIndex[indexNM];
        for (m=0;m<N_J_DOF_NM;m++)
          AuxPtr[Numbers[J_DOF_NM[infoNM][m]]] += N_J_DOF_M;

        Numbers = GlobalNumbers + BeginIndex[indexM];
        for (m=0;m<N_J_DOF_M;m++)
          AuxPtr[Numbers[J_DOF_M[infoM][m]]] += N_J_DOF_NM;

        // check which side is next
        if (lamM0 >= lamNM1 || lamM1 >= lamNM1) NewMortarSideEle = false;
      }
    }
  }
#endif                                          // __ADD_LINK__
#endif                                          // __MORTAR__

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=N_Rows;
  l=AuxPtr[0];
  AuxPtr[0]=0;

  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
  }

  /*
  cout << endl;
  cout << "AuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << AuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound: " << AuxPtr[N_Rows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows];                               // upper bound for number of matrix entries
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  // sum up the array HangingAuxPtr, HangingAuxPtr[i] will now contain
  // the index for Hanging KColAux array, the column numbers for row i
  // are in the intervall // [ HangingAuxPtr[i], HangingAuxPtr[i+1] )
  N_=N_Hanging;
  l=HangingAuxPtr[0];
  HangingAuxPtr[0]=0;
  // cout << HangingAuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=HangingAuxPtr[i+1];
    HangingAuxPtr[i+1]=HangingAuxPtr[i]+l;
    l=k;
    // cout << HangingAuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "HangingAuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << HangingAuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound for hanging nodes: ";
  // cout << HangingAuxPtr[N_Hanging] << endl;

  // get memory for HangingKColAux array, initialize it with -1
  l=HangingAuxPtr[N_Hanging];                     //upper bound for number of matrix entries
  // cout << "l= " << l << endl;
  HangingKColAux=new int[l];
  memset(HangingKColAux, -1, sizeof(int)*l);
  HangingN_Entries=0;

  N_=Space->GetN_Cells();
  Coll = FESpace->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = Space->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      for(k=0;k<n1;k++)
      {
        m=Numbers[k];
        n=Numbers[j];
        if(n<ActiveBound)
        {
          // this  node is a real node (inner or Neumann)
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            N_Entries++;
          }
        }
        else
        {
          if(n<HangingBound)
          {
            // this node is a hanging node in  space
            index=HangingAuxPtr[n-ActiveBound];
            l=HangingKColAux[index];
            // check whether this column is already in this row
            while(l!=-1 && l!=m)
            {
              index++; l=HangingKColAux[index];
            }
            if(l==-1)
            {
              // this is a new column for this row
              HangingKColAux[index]=m;
              HangingN_Entries++;
            }
          }                                       // endif
        }                                         // endif
      }                                           // endfor k

      // additional entries necessary if integrals over the edges
      // appear in the discretization
      if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
      {
        r=cell->GetN_Edges();
        for(p=0;p<r;p++)                          // for all edges
        {
                                                  // neighbour cell
          neigh=cell->GetJoint(p)->GetNeighbour(cell);
          if(neigh)
          {
            q = neigh->GetClipBoard();
            NumbersNeighbour=GlobalNumbers+BeginIndex[q];

            CurrentNeighbour = Space->GetFE2D(q, neigh);
            n2 = TFEDatabase2D::GetFE2D(CurrentNeighbour)->GetN_DOF();

            for(k=0;k<n2;k++)
            {
              m=NumbersNeighbour[k];
              n=Numbers[j];
              if(n<ActiveBound)
              {
                // this  node is a real node (inner or Neumann)
                index=AuxPtr[n];
                l=KColAux[index];
                // check whether this column is already in this row
                while(l!=-1 && l!=m)
                {
                  index++; l=KColAux[index];
                }
                if(l==-1)
                {
                  // this is a new column for this row
                  KColAux[index]=m;
                  N_Entries++;
                }
              }
              else
              {
                if(n<HangingBound)
                {
                  // this node is a hanging node in  space
                  index=HangingAuxPtr[n-ActiveBound];
                  l=HangingKColAux[index];
                  // check whether this column is already in this row
                  while(l!=-1 && l!=m)
                  {
                    index++; l=HangingKColAux[index];
                  }
                  if(l==-1)
                  {
                    // this is a new column for this row
                    HangingKColAux[index]=m;
                    HangingN_Entries++;
                  }
                }                                 // endif
              }                                   // endif
            }                                     // endfor k
          }                                       // endif
        }                                         // endfor p
      }
    }                                             // endfor j
  }                                               // endfor i

  // for(i=0;i<AuxPtr[N_Rows]_;i++)
  // {
  // cout << "KColAux: " << KColAux[i] << endl;
  // }

#ifdef __MORTAR__
#ifdef __ADD_LINK__
  // put additional link DOFs into structure
  // loop over all mortar edges
  int an, tn, *Numbers1, *Numbers2;
  int NeumannBound = N_Inner + BoundNodeBounds[0];

  N_ = It1->GetN_MortarFace();
  for (i=0;i<N_;i++)
  {
    NewMortarSideEle = true;
    lev = MAX_ItLevel + (i << 8);
    It1->Init(lev);
    It2->Init(-lev);

    It1->GetPoint(startX, startY);
    It2->GetPoint(endX, endY);
    delX = endX - startX;
    delY = endY - startY;

    while (CellNM = It1->Next(infoNM))
    {
      // get a new element on non-mortar side
      lamNM0 = GetLambda(startX, startY, CellNM->GetVertex(infoNM),
        delX, delY);
      lamNM1 = GetLambda(startX, startY, CellNM->GetVertex((infoNM+1)
        % CellNM->GetN_Vertices()), delX, delY);
      if (CellNM->GetClipBoard() == -1)
      {
        cerr << "Error in MatrixStructure: cell out of collection!!!"
          << endl;
        exit(-3);
      }

      indexNM = CellNM->GetClipBoard();
      CurrElementID = Space2D->GetFE2D(indexNM, CellNM);
      CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
      CurrDesc = CurrElement->GetFEDesc2D();

      J_DOF_NM = CurrDesc->GetJointDOF();
      N_J_DOF_NM = CurrDesc->GetN_JointDOF();

      if (!NewMortarSideEle)
      {
        // add DOFs if the next while clause wont be entered
        Numbers1 = GlobalNumbers + BeginIndex[indexNM];
        Numbers2 = GlobalNumbers + BeginIndex[indexM];
        for (n=0;n<N_J_DOF_NM;n++)
        {
          for (m=0;m<N_J_DOF_M;m++)
          {
            tn = Numbers1[J_DOF_NM[infoNM][n]];
            an = Numbers2[J_DOF_M[infoM][m]];

            if (tn < NeumannBound)
            {
              index = AuxPtr[tn];
              l = KColAux[index];

              // check whether this column is already in this row
              while (l != -1 && l != an)
                l = KColAux[++index];

              if (l == -1)
              {
                // this is a new column for this row
                KColAux[index] = an;
                N_Entries++;
              }
            }

            if (an < NeumannBound)
            {
              index = AuxPtr[an];
              l = KColAux[index];

              // check whether this column is already in this row
              while (l != -1 && l != tn)
                l = KColAux[++index];

              if (l == -1)
              {
                // this is a new column for this row
                KColAux[index] = tn;
                N_Entries++;
              }
            }
          }
        }
      }

      // check which side is next
      if (lamM1 < lamNM1) NewMortarSideEle = true;

      while (NewMortarSideEle)
      {
        // get a new element on mortar side
        if (!(CellM = It2->Next(infoM))) break;
        lamM0 = GetLambda(startX, startY, CellM->GetVertex((infoM+1)
          % CellM->GetN_Vertices()), delX, delY);
        lamM1 = GetLambda(startX, startY, CellM->GetVertex(infoM),
          delX, delY);
        if (CellM->GetClipBoard() == -1)
        {
          cerr << "Error in MatrixStructure: cell out of collection!!!"
            << endl;
          exit(-3);
        }

        indexM = CellM->GetClipBoard();
        CurrElementID = Space2D->GetFE2D(indexM, CellM);
        CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc2D();

        J_DOF_M = CurrDesc->GetJointDOF();
        N_J_DOF_M = CurrDesc->GetN_JointDOF();

        // add DOFs
        Numbers1 = GlobalNumbers + BeginIndex[indexNM];
        Numbers2 = GlobalNumbers + BeginIndex[indexM];
        for (n=0;n<N_J_DOF_NM;n++)
        {
          for (m=0;m<N_J_DOF_M;m++)
          {
            tn = Numbers1[J_DOF_NM[infoNM][n]];
            an = Numbers2[J_DOF_M[infoM][m]];

            if (tn < NeumannBound)
            {
              index = AuxPtr[tn];
              l = KColAux[index];

              // check whether this column is already in this row
              while (l != -1 && l != an)
                l = KColAux[++index];

              if (l == -1)
              {
                // this is a new column for this row
                KColAux[index] = an;
                N_Entries++;
              }
            }

            if (an < NeumannBound)
            {
              index = AuxPtr[an];
              l = KColAux[index];

              // check whether this column is already in this row
              while (l != -1 && l != tn)
                l = KColAux[++index];

              if (l == -1)
              {
                // this is a new column for this row
                KColAux[index] = tn;
                N_Entries++;
              }
            }
          }
        }

        // check which side is next
        if (lamM0 >= lamNM1 || lamM1 >= lamNM1) NewMortarSideEle = false;
      }
    }
  }
#endif                                          // __ADD_LINK__
#endif                                          // __MORTAR__

  /*
  // check hanging node data
  cout << endl;
  cout << "check hanging node data" << endl;
  N_=N_Hanging;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<index && HangingKColAux[j]!=-1;j++)
      cout << setw(4) << HangingKColAux[j];
  cout << endl;
  }
  */

  /*
  cout << "Number of matrix entries (hanging nodes): ";
  cout << HangingN_Entries << endl;
  cout << endl;
  */

  // compress HangingKColAux array to HangingKCol by deleting all -1's
  // build the HangingRowPtr array
  N_=N_Hanging;
  HangingKCol=new int[HangingN_Entries];
  HangingRowPtr=HangingAuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<m && HangingKColAux[j]!=-1;j++)
    {
      HangingKCol[index]=HangingKColAux[j];
      index++;
    }                                             // endfor j
    HangingRowPtr[i]=oldindex;
    // cout << HangingRowPtr[i] << endl;
  }                                               // endfor i
  HangingRowPtr[N_]=index;

  // free HangingKColAux
  delete [] HangingKColAux;

  // add the additional columns from hanging nodes to other nodes
  Offset=ActiveBound;
  HangingNodes=Space->GetHangingNodes();
  for(i=0;i<N_Hanging;i++)
  {
    // cout << "hanging node: " << i << endl;
    hn=HangingNodes[i];
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();

    for(j=0;j<n;j++)                              // loop over all nodes in coupling
    {
      k=DOF[j];
      // cout << "k= " << k << endl;
      if(k<ActiveBound)
      {
        // node is either inner or Neumann node

        end=HangingAuxPtr[i+1];
        for(oldindex=HangingAuxPtr[i];oldindex<end;oldindex++)
        {
          m=HangingKCol[oldindex];
          // cout << "m= " << m << endl;
          index=AuxPtr[k];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            N_Entries++;
          }
        }                                         // endfor
      }                                           // endif
    }                                             // endfor j
  }                                               // endfor i

  /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=ActiveBound;
  N_=N_Rows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++)
  cout << setw(4) << KColAux[j];
  cout << endl;
  }
  */

  /*
  cout << "Number of matrix entries: ";
  cout << N_Entries << endl;
  cout << endl;
  */

  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=ActiveBound;
  KCol=new int[N_Entries];
  RowPtr=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      KCol[index]=KColAux[j];
      index++;
    }                                             // endfor j
    RowPtr[i]=oldindex;
    // cout << setw(4) << i << " RowPtr[i]: " << RowPtr[i] << endl;
  }                                               // endfor i

  // cout << "index: " << index << endl;
  // cout << "RowPtr[N_]: " << RowPtr[N_] << endl;
  Offset=index-RowPtr[N_];
  for(i=0,j=ActiveBound;i<=N_Hanging;i++,j++)
  {
    // cout << "HangingRowPtr[i]: " << HangingRowPtr[i] << endl;
    RowPtr[j]+=Offset;
    // cout << setw(4) << j << " RowPtr[j]: " << RowPtr[j] << endl;
  }

  j=ActiveBound+N_Hanging;
  Offset=RowPtr[j];
  for(i=0;i<N_Dirichlet;i++,j++)
  {
    RowPtr[j+1]=RowPtr[j]+1;
    // cout << setw(4) << j+1 << " RowPtr[j+1]: " << RowPtr[j+1] << endl;
  }

  // add information for hanging and Dirichlet nodes into matrix
  HangingNodes=Space->GetHangingNodes();
  Offset=ActiveBound;
  m=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();
    index=AuxPtr[j];
    KCol[index]=m;
    index++;
    m++;
    for(k=0;k<n;k++)
    {
      // cout << "index: " << index << " DOF[k]" << DOF[k] << endl;
      KCol[index]=DOF[k];
      index++;
    }
  }

  // add Dirichlet rows
  j=HangingBound;
  index=RowPtr[ActiveBound+N_Hanging];
  for(i=0;i<N_Dirichlet;i++)
  {
    // cout << "index: " << index << endl;
    KCol[index]=j;
    j++;
    index++;
  }

  /*    // print out the whole matrix structure
      cout << endl;
      N_=N_Rows;
      for(i=0;i<N_;i++)
      {
        cout << RowPtr[i] << "---" << RowPtr[i+1]-1 << endl;
        cout << "Rows: " << setw(4) << i << ": ";
        end=RowPtr[i+1];
        for(j=RowPtr[i];j<end;j++)
          cout << setw(4) << KCol[j];
        cout << endl;
  }
  */

  // free KColAux
  delete [] KColAux;

#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif
  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
    cout << endl;
    cout << "Information on the stored matrix structure" << endl;
    cout << "Number of rows: " << N_Rows << endl;
    cout << "Number of columns: " << N_Columns << endl;
    cout << "Number of matrix entries: " << N_Entries << endl;
  }

}



/** generate the matrix structure, all arrays are already defined */
TSquareStructure2D::TSquareStructure2D(int n, int N_entries, int *col_ptr,
int *row_ptr)
: TSquareStructure()
{
  N_Rows = n;
  N_Columns = n;
  ActiveBound = n;
  N_Entries = N_entries;
  HangingN_Entries = 0;
  KCol = col_ptr;
  HangingKCol = NULL;
  RowPtr = row_ptr;
  HangingRowPtr = NULL;
  FESpace = NULL;
  ColOrder = 0;
}

TSquareStructure2D::TSquareStructure2D(int n)
 : TSquareStructure(n), FESpace(NULL)
{
}

TSquareStructure2D::~TSquareStructure2D()
{
}


/** check if SquareStructure lhs and rhs are equal*/
bool operator==(const TSquareStructure2D &lhs, const TSquareStructure2D &rhs)
{
  if(&lhs == &rhs)
    return true;
  if(lhs.FESpace == rhs.FESpace)
  {
    const TSquareStructure* lhsTmp = &lhs;
    const TSquareStructure* rhsTmp = &rhs;
    return *lhsTmp==*rhsTmp;
  }
  return false;
}

/** check if SquareStructure lhs and rhs are not equal*/
bool operator!=(const TSquareStructure2D &lhs, const TSquareStructure2D &rhs)
{
  return !(rhs == lhs);
}
