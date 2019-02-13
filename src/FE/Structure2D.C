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
// @(#)Structure2D.C        1.10 11/24/99
// 
// Class:       TStructure2D
//
// Purpose:     build and store a matrix structure
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation (Gunar Matthies)
//
//              04.08.1998 start reimplementation (Gunar Matthies)
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <DefineParams.h>

#include <FEDatabase2D.h>
#include <MooNMD_Io.h>
#include <Structure2D.h>
#include <string.h>
#include <stdlib.h>

#ifdef __MORTAR__
#ifdef __ADD_LINK__
  #include <Database.h>
  #include <It_Mortar.h>
  #include <MortarBaseJoint.h>
  #include <stdlib.h>
#endif // __ADD_LINK__
#endif // __MORTAR__
#include <Database.h>

/** generate the matrix structure, both spaces are 2D */
TStructure2D::TStructure2D(TFESpace2D *testspace, TFESpace2D *ansatzspace)
{
  TCollection *coll;
  TBaseCell *cell;
  int i,j,k,l,m,n, N_, n1, n2, index, oldindex;
  int TestN_BoundNodeTypes, AnsatzN_BoundNodeTypes;
  int TestN_Dirichlet, AnsatzN_Dirichlet;
  int TestN_Inner, AnsatzN_Inner;
  int *TestN_BoundNodes, *AnsatzN_BoundNodes;
  int TestSumBoundNodes, AnsatzSumBoundNodes;
  int *TestGlobalNumbers, *AnsatzGlobalNumbers;
  int *TestNumbers, *AnsatzNumbers;
  int *TestBeginIndex, *AnsatzBeginIndex;
  int *AuxPtr, *KColAux, end;
  FE2D CurrentElement; 

  int TestN_Hanging, AnsatzN_Hanging;
  int *HangingAuxPtr;
  int TestActiveBound, TestHangingBound;
  int AnsatzActiveBound, AnsatzHangingBound;
  int Offset, *DOF, begin;
  THangingNode **TestHangingNodes;
  THangingNode **AnsatzHangingNodes;
  THangingNode *hn;
  int *HangingKColAux;

  int N, N1, N2;

  // test if both spaces are defined on the same triangulation
  if(testspace->GetCollection() != ansatzspace->GetCollection())
  {
    OutPut("Structure2D.C : grid for test and ansatz space is not the same !" << endl); 
    exit(4711); 
    return;
  }

  // get collection of mesh cells
  coll = testspace->GetCollection();

  // test space and ansatz space differ
  TestSpace2D = testspace;
  AnsatzSpace2D = ansatzspace;
  TestSpace1D = NULL;
  AnsatzSpace1D = NULL;

  // get information from the spaces
  AnsatzN_BoundNodeTypes = ansatzspace->GetN_DiffBoundaryNodeTypes();
  TestN_BoundNodeTypes = testspace->GetN_DiffBoundaryNodeTypes();

  AnsatzN_Dirichlet = ansatzspace->GetN_Dirichlet();
  TestN_Dirichlet = testspace->GetN_Dirichlet();

  AnsatzN_Inner = ansatzspace->GetN_Inner();
  TestN_Inner = testspace->GetN_Inner();

  AnsatzN_BoundNodes = ansatzspace->GetN_BoundaryNodes();
  TestN_BoundNodes = testspace->GetN_BoundaryNodes();

  AnsatzSumBoundNodes = 0;
  for(i=0;i<AnsatzN_BoundNodeTypes;i++)
    AnsatzSumBoundNodes += AnsatzN_BoundNodes[i];

  TestSumBoundNodes = 0;
  for(i=0;i<TestN_BoundNodeTypes;i++)
    TestSumBoundNodes += TestN_BoundNodes[i];

  AnsatzGlobalNumbers = ansatzspace->GetGlobalNumbers();
  TestGlobalNumbers = testspace->GetGlobalNumbers();

  AnsatzBeginIndex = ansatzspace->GetBeginIndex();
  TestBeginIndex = testspace->GetBeginIndex();

  N_Rows = TestSpace2D->GetN_DegreesOfFreedom();
  N_Columns = AnsatzSpace2D->GetN_DegreesOfFreedom();

  TestN_Hanging = TestSpace2D->GetN_Hanging();
  TestActiveBound = TestSpace2D->GetN_ActiveDegrees();
  TestHangingBound= TestSpace2D->GetHangingBound();
  TestHangingNodes = TestSpace2D->GetHangingNodes();

  AnsatzN_Hanging = AnsatzSpace2D->GetN_Hanging();
  AnsatzActiveBound = AnsatzSpace2D->GetN_ActiveDegrees();
  AnsatzHangingBound= AnsatzSpace2D->GetHangingBound();
  AnsatzHangingNodes = AnsatzSpace2D->GetHangingNodes();

  // AuxPtr[i] will contain an upper bound for the number of 
  // matrix entries in row i
  l = N_Rows + 1;
  AuxPtr = new int[l];
  memset(AuxPtr, 0, l*SizeOfInt);

  l=TestN_Hanging+1;
  HangingAuxPtr=new int[l];
  memset(HangingAuxPtr, 0, l*sizeof(int));

  Offset = TestActiveBound;

  // number of mesh cells
  N_ = coll->GetN_Cells();
  // loop over all mesh cells
  for(i=0;i<N_;i++)
  {
    // get cell i 
    cell = coll->GetCell(i);
    // get global number of d.o.f. of test and ansatz function 
    // which live on this mesh cell
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];
    
    // get fe spaces on the mesh cell
    CurrentElement = TestSpace2D->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    CurrentElement = AnsatzSpace2D->GetFE2D(i, cell);
    n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
    // update vector which stores information on the length of the rows
    for(j=0;j<n1;j++)
    {
      k = TestNumbers[j];
      // real dof or Dirichlet node
      if(k<TestActiveBound || k>=TestHangingBound)
      {
        AuxPtr[k] += n2;
      }
      else
      {
        HangingAuxPtr[k-Offset]+=n2;
      } // endif

      if(AnsatzN_Hanging)
      {
        for(l=0;l<n2;l++)
        {
          m = AnsatzNumbers[l];
          if(m>=AnsatzActiveBound && m<AnsatzHangingBound)
          {
            // hanging node in ansatz space
            hn = AnsatzHangingNodes[m-AnsatzActiveBound];
            n = TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
            AuxPtr[k] += n;
          }
        }
      } // endif AnsatzN_Hanging
    } // endfor j
  } // endfor i

  // additional space for storing the columns caused by the 
  // hanging nodes of  space => some new columns
  Offset=TestActiveBound;
  // cout << "TestN_Hanging: " << TestN_Hanging << endl;
  for(i=0;i<TestN_Hanging;i++)
  {
    hn=TestHangingNodes[i];
    m=HangingAuxPtr[i];
    // cout << "HangingAuxPtr[i]: " << m << endl;
    DOF=hn->GetDOF();
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    for(j=0;j<n;j++)
    {
      k=DOF[j];
      if(k<TestActiveBound)
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
  N_ = TestSpace2D->GetN_Cells();
  coll = TestSpace2D->GetCollection();
  for(i=0;i<N_;i++)
    coll->GetCell(i)->SetClipBoard(i);

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
        cerr << "Error in MatrixStructure2D: cell out of collection!!!"
             << endl;
        exit(-3);
      }

      indexNM = CellNM->GetClipBoard();
      CurrElementID = AnsatzSpace2D->GetFE2D(indexNM, CellNM);
      CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
      CurrDesc = CurrElement->GetFEDesc2D();

      J_DOF_NM = CurrDesc->GetJointDOF();
      N_J_DOF_NM = CurrDesc->GetN_JointDOF();

      if (!NewMortarSideEle)
      {
        // add space for DOFs if the next while clause wont be
        // entered
        TestNumbers = TestGlobalNumbers + TestBeginIndex[indexNM];
        for (m=0;m<N_J_DOF_NM;m++)
          AuxPtr[TestNumbers[J_DOF_NM[infoNM][m]]] += N_J_DOF_M;

        TestNumbers = TestGlobalNumbers + TestBeginIndex[indexM];
        for (m=0;m<N_J_DOF_M;m++)
          AuxPtr[TestNumbers[J_DOF_M[infoM][m]]] += N_J_DOF_NM;
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
          cerr << "Error in MatrixStructure2D: cell out of collection!!!"
               << endl;
          exit(-3);
        }

        indexM = CellM->GetClipBoard();
        CurrElementID = AnsatzSpace2D->GetFE2D(indexM, CellM);
        CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc2D();

        J_DOF_M = CurrDesc->GetJointDOF();
        N_J_DOF_M = CurrDesc->GetN_JointDOF();

        // add space for DOFs
        TestNumbers = TestGlobalNumbers + TestBeginIndex[indexNM];
        for (m=0;m<N_J_DOF_NM;m++)
          AuxPtr[TestNumbers[J_DOF_NM[infoNM][m]]] += N_J_DOF_M;

        TestNumbers = TestGlobalNumbers + TestBeginIndex[indexM];
        for (m=0;m<N_J_DOF_M;m++)
          AuxPtr[TestNumbers[J_DOF_M[infoM][m]]] += N_J_DOF_NM;

        // check which side is next
        if (lamM0 >= lamNM1 || lamM1 >= lamNM1) NewMortarSideEle = false;
      }
    }
  }
#endif // __ADD_LINK__
#endif // __MORTAR__

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=N_Rows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  // cout << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    // cout << "i: " << i << "  l: " << l << endl;
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    // cout << AuxPtr[i+1] << endl;
  }

  // cout << "Upper bound: " << AuxPtr[N_Rows] << endl;
  
  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows]; // upper bound for number of matrix entries 
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  // sum up the array HangingAuxPtr, HangingAuxPtr[i] will now contain 
  // the index for Hanging KColAux array, the column numbers for row i 
  // are in the intervall // [ HangingAuxPtr[i], HangingAuxPtr[i+1] )
  N_=TestN_Hanging;
  l=HangingAuxPtr[0];
  HangingAuxPtr[0]=0;
//   cout << HangingAuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=HangingAuxPtr[i+1];
    HangingAuxPtr[i+1]=HangingAuxPtr[i]+l;
    l=k;
    // cout << HangingAuxPtr[i+1] << endl;
  }

  // get memory for HangingKColAux array, initialize it with -1
  l=HangingAuxPtr[TestN_Hanging]; //upper bound for number of matrix entries 
  // cout << "l= " << l << endl;
  HangingKColAux=new int[l];
  memset(HangingKColAux, -1, sizeof(int)*l);
  HangingN_Entries=0;

  N_Entries=0;
  N_ = coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];

    cell = coll->GetCell(i);

    CurrentElement = TestSpace2D->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    CurrentElement = AnsatzSpace2D->GetFE2D(i, cell);
    n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;
    
    for(j=0;j<n1;j++)
    {
      for(k=0;k<n2;k++)
      {
        m=AnsatzNumbers[k];
        n=TestNumbers[j];

        if(n<TestActiveBound || n>=TestHangingBound)
        {
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m && index<AuxPtr[n+1])
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
          // this node is a hanging node in  space
          index=HangingAuxPtr[n-TestActiveBound];
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
        } // endif

        // hanging nodes in ansatz space
        if(m>=AnsatzActiveBound && m<AnsatzHangingBound)
        {
          hn = AnsatzHangingNodes[m-AnsatzActiveBound];
          N2 = TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
          DOF = hn->GetDOF();
          for(N1=0;N1<N2;N1++)
          {
            index=AuxPtr[n];
            l = KColAux[index];
            N = DOF[N1];
            while(l!=-1 && l!=N && index<AuxPtr[n+1])
            {
              index++; l=KColAux[index];
            }
            if(l==-1)
            {
              // this is a new column for this row
              KColAux[index]=N;
              N_Entries++;
            }
          } // endfor N1
        } // endif
      } // endfor k
    } // endfor j
  } // endfor i
  
#ifdef __MORTAR__
#ifdef __ADD_LINK__
  // put additional link DOFs into structure
  // loop over all mortar edges
  int an, tn;
  int TestNeumannBound = TestN_Inner + TestN_BoundNodes[0];
  int AnsatzNeumannBound = AnsatzN_Inner + AnsatzN_BoundNodes[0];

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
        cerr << "Error in MatrixStructure2D: cell out of collection!!!"
             << endl;
        exit(-3);
      }

      indexNM = CellNM->GetClipBoard();
      CurrElementID = AnsatzSpace2D->GetFE2D(indexNM, CellNM);
      CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
      CurrDesc = CurrElement->GetFEDesc2D();

      J_DOF_NM = CurrDesc->GetJointDOF();
      N_J_DOF_NM = CurrDesc->GetN_JointDOF();

      if (!NewMortarSideEle)
      {
        // add DOFs if the next while clause wont be entered
        TestNumbers = TestGlobalNumbers + TestBeginIndex[indexNM];
        AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[indexM];
        for (n=0;n<N_J_DOF_NM;n++)
        {
          for (m=0;m<N_J_DOF_M;m++)
          {
            an = AnsatzNumbers[J_DOF_M[infoM][m]];
            tn = TestNumbers[J_DOF_NM[infoNM][n]];

            if (tn < TestNeumannBound)
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

            if (an < AnsatzNeumannBound)
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
          cerr << "Error in MatrixStructure2D: cell out of collection!!!"
               << endl;
          exit(-3);
        }

        indexM = CellM->GetClipBoard();
        CurrElementID = AnsatzSpace2D->GetFE2D(indexM, CellM);
        CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc2D();

        J_DOF_M = CurrDesc->GetJointDOF();
        N_J_DOF_M = CurrDesc->GetN_JointDOF();

        // add DOFs
        TestNumbers = TestGlobalNumbers + TestBeginIndex[indexNM];
        AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[indexM];
        for (n=0;n<N_J_DOF_NM;n++)
        {
          for (m=0;m<N_J_DOF_M;m++)
          {
            an = AnsatzNumbers[J_DOF_M[infoM][m]];
            tn = TestNumbers[J_DOF_NM[infoNM][n]];

            if (tn < TestNeumannBound)
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

            if (an < AnsatzNeumannBound)
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
#endif // __ADD_LINK__
#endif // __MORTAR__

/*
  // check
  cout << endl;
  cout << "check" << endl;
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

  // compress HangingKColAux array to HangingKCol by deleting all -1's
  // build the HangingRowPtr array
  N_=TestN_Hanging;
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
    } // endfor j
    HangingRowPtr[i]=oldindex;
    // cout << HangingRowPtr[i] << endl;
  } // endfor i
  HangingRowPtr[N_]=index;

  // free HangingKColAux
  delete [] HangingKColAux;

  // add the additional columns from hanging nodes to other nodes
  Offset = TestActiveBound;
  for(i=0;i<TestN_Hanging;i++)
  {
    // cout << "hanging node: " << i << endl;
    hn=TestHangingNodes[i];
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();

    for(j=0;j<n;j++) // loop over all nodes in coupling
    {
      k=DOF[j];
      // cout << "k= " << k << endl;
      if(k<TestActiveBound)
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
        } // endfor
      } // endif
    } // endfor j
  } // endfor i

  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=N_Rows;
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
    } // endfor j
    RowPtr[i]=oldindex;
    // cout << setw(4) << i << " RowPtr[i]: " << RowPtr[i] << endl;
  } // endfor i
  RowPtr[N_Rows]=N_Entries;

  Sort();
  
/*
  // print out the whole matrix structure
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

#ifdef __MORTAR__
#include <MortarBaseJoint.h>

/** generate the matrix structure, one space 1D and one 2D */
TStructure2D::TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace)
{
  TestSpace1D = testspace;
  AnsatzSpace2D = ansatzspace;
  TestSpace2D = NULL;
  AnsatzSpace1D = NULL;

  TBaseCell *CurrCell1D, *CurrCell2D;
  TCollection *coll1D = TestSpace1D->GetCollection();
  TCollection *coll2D = AnsatzSpace2D->GetCollection();
  const int *TmpEV;
  bool MortarSide;
  double X, Y, X0, Y0, X1, Y1, X2, Y2, X3, Y3;
  TVertex *CurrVertex;
  int **JointDOFs, N_JointDOFs;
  int N_, i, j, k, l, m, n, n2, index, oldindex, end, AN, TN, entry, l0, l1;
  int TestN_Inner, TestN_Dirichlet, TestN_Neumann;
  int AnsatzN_Inner, AnsatzN_Dirichlet, AnsatzN_Neumann, AnsatzN_Hanging;
  int *AuxPtr, *KColAux;
  int *TestGlobalNumbers, *AnsatzGlobalNumbers;
  int *TestBeginIndex, *AnsatzBeginIndex;
  int *TestNumbers, *AnsatzNumbers;
  FE2D CurrElementID;
  TFE2D *CurrElement;
  TFEDesc2D *CurrDesc;
  JointType type;

  TestN_Dirichlet = TestSpace1D->GetN_Dirichlet();
  TestN_Neumann = TestSpace1D->GetN_BoundaryNodes()[0];
  TestN_Inner = TestSpace1D->GetN_Inner();

  AnsatzN_Dirichlet = AnsatzSpace2D->GetN_Dirichlet();
  AnsatzN_Neumann = AnsatzSpace2D->GetN_BoundaryNodes()[0];
  AnsatzN_Inner = AnsatzSpace2D->GetN_Inner();
  AnsatzN_Hanging = AnsatzSpace2D->GetN_Hanging();

  N_Rows = TestN_Inner + TestN_Neumann + TestN_Dirichlet;
  N_Columns = AnsatzN_Inner + AnsatzN_Neumann + AnsatzN_Dirichlet +
              AnsatzN_Hanging;
  N_Entries = 0;

  // AuxPtr[i] will contain an upper bound for the number of
  // matrix entries in row i
  AuxPtr = new int[N_Rows+1];
  memset(AuxPtr, 0, (N_Rows+1)*SizeOfInt);

  TestGlobalNumbers = TestSpace1D->GetGlobalNumbers();
  AnsatzGlobalNumbers = AnsatzSpace2D->GetGlobalNumbers();

  TestBeginIndex = TestSpace1D->GetBeginIndex();
  AnsatzBeginIndex = AnsatzSpace2D->GetBeginIndex();

  // loop over all elements
  N_ = AnsatzSpace2D->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    CurrCell2D = coll2D->GetCell(i);
    k = CurrCell2D->GetN_Edges();
    for (j=0;j<k;j++)
    {
      type = CurrCell2D->GetJoint(j)->GetType();
      if (type == MortarJoint || type == MortarBaseJoint)
      {
  // !!! important
  // at the moment, this is only made for conforming elements

        CurrElementID = AnsatzSpace2D->GetFE2D(i, CurrCell2D);
        CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc2D();

        JointDOFs = CurrDesc->GetJointDOF();
        N_JointDOFs = CurrDesc->GetN_JointDOF();

        CurrCell2D->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV);
        CurrVertex = CurrCell2D->GetVertex(TmpEV[2*j]);
        X0 = CurrVertex->GetX();
        Y0 = CurrVertex->GetY();
        CurrVertex = CurrCell2D->GetVertex(TmpEV[2*j + 1]);
        X1 = CurrVertex->GetX();
        Y1 = CurrVertex->GetY();

        if (type == MortarJoint)
          l = ((TMortarJoint *) CurrCell2D->GetJoint(j))->GetMEdgeInColl();
        else
          l = ((TMortarBaseJoint *) CurrCell2D->GetJoint(j))->GetMEdgeInColl();
        if (l < 0)
        {
          MortarSide = true;
          l = -(l + 1);
        }
        else
          MortarSide = false;

        CurrCell1D = coll1D->GetCell(l);

        CurrVertex = CurrCell1D->GetVertex(0);
        X2 = CurrVertex->GetX();
        Y2 = CurrVertex->GetY();
        CurrVertex = CurrCell1D->GetVertex(1);
        X3 = CurrVertex->GetX();
        Y3 = CurrVertex->GetY();

        TestNumbers = TestGlobalNumbers + TestBeginIndex[l];
        n = TestBeginIndex[l+1] - TestBeginIndex[l];

        for (m=0;m<n;m++)
          AuxPtr[TestNumbers[m]] += N_JointDOFs;

        while ((CurrCell1D = CurrCell1D->GetJoint(1)->
                 GetNeighbour(CurrCell1D)) && MortarSide)
        {
          CurrVertex = CurrCell1D->GetVertex(0);
          X = CurrVertex->GetX();
          Y = CurrVertex->GetY();
          if ((X0-X)*(X0-X1) + (Y0-Y)*(Y0-Y1) < 1e-16) break;

          TestNumbers = TestGlobalNumbers + TestBeginIndex[++l];
          n = TestBeginIndex[l+1] - TestBeginIndex[l];

          for (m=0;m<n;m++)
            AuxPtr[TestNumbers[m]] += N_JointDOFs;

        }
      }
    }
  }

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  l = AuxPtr[0];
  AuxPtr[0] = 0;
  for(i=0;i<N_Rows;i++)
  {
    k = AuxPtr[i+1];
    AuxPtr[i+1] = AuxPtr[i] + l;
    l = k;
  }

  // cout << "Upper bound: " << AuxPtr[N_Rows] << endl;

  // get memory for KColAux array, initialize it with -1
  KColAux = new int[AuxPtr[N_Rows]];
  memset(KColAux, -1, AuxPtr[N_Rows]*SizeOfInt);

  N_ = AnsatzSpace2D->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    CurrCell2D = coll2D->GetCell(i);
    k = CurrCell2D->GetN_Edges();
    for (j=0;j<k;j++)
    {
      type = CurrCell2D->GetJoint(j)->GetType();
      if (type == MortarJoint || type == MortarBaseJoint)
      {
  // !!! important
  // at the moment, this is only made for conforming elements

        CurrElementID = AnsatzSpace2D->GetFE2D(i, CurrCell2D);
        CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc2D();

        JointDOFs = CurrDesc->GetJointDOF();
        N_JointDOFs = CurrDesc->GetN_JointDOF();

        CurrCell2D->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV);
        CurrVertex = CurrCell2D->GetVertex(TmpEV[2*j]);
        X0 = CurrVertex->GetX();
        Y0 = CurrVertex->GetY();
        CurrVertex = CurrCell2D->GetVertex(TmpEV[2*j + 1]);
        X1 = CurrVertex->GetX();
        Y1 = CurrVertex->GetY();

        if (type == MortarJoint)
          l0 = ((TMortarJoint *) CurrCell2D->GetJoint(j))->GetMEdgeInColl();
        else
          l0 = ((TMortarBaseJoint *) CurrCell2D->GetJoint(j))->GetMEdgeInColl();
        if (l0 < 0)
        {
          MortarSide = true;
          l0 = -(l0 + 1);
        }
        else
          MortarSide = false;

        l1 = l0;

        CurrCell1D = coll1D->GetCell(l1);

        CurrVertex = CurrCell1D->GetVertex(0);
        X2 = CurrVertex->GetX();
        Y2 = CurrVertex->GetY();
        CurrVertex = CurrCell1D->GetVertex(1);
        X3 = CurrVertex->GetX();
        Y3 = CurrVertex->GetY();

        while ((CurrCell1D = CurrCell1D->GetJoint(1)->
                 GetNeighbour(CurrCell1D)) && MortarSide)
        {
          CurrVertex = CurrCell1D->GetVertex(0);
          X = CurrVertex->GetX();
          Y = CurrVertex->GetY();
          if ((X0-X)*(X0-X1) + (Y0-Y)*(Y0-Y1) < 1e-16) break;

          l1++;
        }

        AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];
        for (l=l0;l<=l1;l++)
        {
          TestNumbers = TestGlobalNumbers + TestBeginIndex[l];

          n2 = TestBeginIndex[l+1] - TestBeginIndex[l];

          for (m=0;m<n2;m++)
          {
            TN = TestNumbers[m];
            for (n=0;n<N_JointDOFs;n++)
            {
              AN = AnsatzNumbers[JointDOFs[j][n]];
              index = AuxPtr[TN];
              entry = KColAux[index];
              // check whether this column is already in this row
              while (entry != -1 && entry != AN)
                entry = KColAux[++index];

              if (entry == -1)
              {
                // this is a new column for this row
                KColAux[index] = AN;
                N_Entries++;
              }
            }
          }
        }
      }
    }
  }

  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  KCol = new int[N_Entries];
  RowPtr = AuxPtr;
  
  index = 0;
  for (i=0;i<N_Rows;i++)
  {
    oldindex = index;
    m = AuxPtr[i+1];
    for (j=AuxPtr[i];j<m && KColAux[j] != -1;j++)
      KCol[index++] = KColAux[j];

    RowPtr[i] = oldindex;
  }

  RowPtr[i] = index;

  // free KColAux
  delete [] KColAux;
  
/*
  // print out the whole matrix structure (mortar)
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
  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
  #endif
  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
  cout << "Information on the stored matrix structure of B" << endl;
  cout << "Number of rows: " << N_Rows << endl;
  cout << "Number of columns: " << N_Columns << endl;
  cout << "Number of matrix entries: " << N_Entries << endl;
  }
}
#endif  // __MORTAR__

/** generate the matrix structure, both spaces are 2D */
/** both spaces are defined on different grids */
TStructure2D::TStructure2D(TFESpace2D *testspace, int test_level, 
                           TFESpace2D *ansatzspace, int ansatz_level)
{
  TCollection *coll_coarse, *coll_fine;
  TBaseCell *cell_coarse, *cell_fine;
  TBaseCell *cell_child_1, *cell_child_2, *cell_child_3, *cell_child_4, *cell_child_5;
  TBaseCell *cell_parent, *cell_tmp, *cell_child_6;
  int i,j,k,l,m,n, N_, n1, n2, index, oldindex;
  int TestN_BoundNodeTypes, AnsatzN_BoundNodeTypes;
  int TestN_Dirichlet, AnsatzN_Dirichlet;
  int TestN_Inner, AnsatzN_Inner;
  int *TestN_BoundNodes, *AnsatzN_BoundNodes;
  int TestSumBoundNodes, AnsatzSumBoundNodes;
  int *TestGlobalNumbers, *AnsatzGlobalNumbers;
  int *TestNumbers, *AnsatzNumbers;
  int *TestBeginIndex, *AnsatzBeginIndex;
  int *AuxPtr, *KColAux, end;
  int  N_child,N_child_1, N_child_2, N_child_3, N_child_4, N_child_5;
  int N_child_6;
  int ii, j1, j2, j3, j4, j5, j6, level_diff;
  FE2D CurrentElement; 

  level_diff = abs(test_level-ansatz_level);
  if (level_diff>6)
  {
     OutPut("level difference greater than 5 not implemented !!!" << endl);
     exit(4711);
  }

  if (test_level < ansatz_level)
  {
    // get collection of mesh cells
    coll_coarse = testspace->GetCollection();
    coll_fine = ansatzspace->GetCollection();
  }
  else
  {
    // get collection of mesh cells
    coll_fine = testspace->GetCollection();
    coll_coarse = ansatzspace->GetCollection();
  }

  // test space and ansatz space differ
  TestSpace2D = testspace;
  AnsatzSpace2D = ansatzspace;
  TestSpace1D = NULL;
  AnsatzSpace1D = NULL;

  // get information from the spaces
  AnsatzN_BoundNodeTypes = ansatzspace->GetN_DiffBoundaryNodeTypes();
  TestN_BoundNodeTypes = testspace->GetN_DiffBoundaryNodeTypes();

  AnsatzN_Dirichlet = ansatzspace->GetN_Dirichlet();
  TestN_Dirichlet = testspace->GetN_Dirichlet();

  AnsatzN_Inner = ansatzspace->GetN_Inner();
  TestN_Inner = testspace->GetN_Inner();

  AnsatzN_BoundNodes = ansatzspace->GetN_BoundaryNodes();
  TestN_BoundNodes = testspace->GetN_BoundaryNodes();

  AnsatzSumBoundNodes = 0;
  for(i=0;i<AnsatzN_BoundNodeTypes;i++)
    AnsatzSumBoundNodes += AnsatzN_BoundNodes[i];

  TestSumBoundNodes = 0;
  for(i=0;i<TestN_BoundNodeTypes;i++)
    TestSumBoundNodes += TestN_BoundNodes[i];

  AnsatzGlobalNumbers = ansatzspace->GetGlobalNumbers();
  TestGlobalNumbers = testspace->GetGlobalNumbers();

  AnsatzBeginIndex = ansatzspace->GetBeginIndex();
  TestBeginIndex = testspace->GetBeginIndex();

  N_Rows = TestSpace2D->GetN_DegreesOfFreedom();
  N_Columns = AnsatzSpace2D->GetN_DegreesOfFreedom();

  l = N_Rows + 1;
  AuxPtr = new int[l];
  memset(AuxPtr, 0, l*SizeOfInt);

  // set numeration of the cells
  N_ = coll_coarse->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell_coarse = coll_coarse->GetCell(i);
    cell_coarse->SetClipBoard(i);
  }
  N_ = coll_fine->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell_fine = coll_fine->GetCell(i);
    cell_fine->SetClipBoard(i);
  }

  if (test_level < ansatz_level)
  {
    // number of mesh cells
    N_ = coll_coarse->GetN_Cells();
    // loop over all mesh cells of the test space
    for(i=0;i<N_;i++)
    {
      // get cell i 
       cell_coarse = coll_coarse->GetCell(i);
       // get global number of d.o.f. of test and ansatz function 
       // which live on this mesh cell
       TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
       // get fe spaces on the mesh cell
       CurrentElement = TestSpace2D->GetFE2D(i, cell_coarse);
       n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
       
       // find all children of the coarse mesh cell 
       N_child = cell_coarse->GetN_Children();
       for (j1=0;j1<N_child;j1++)
       {
         // first refinement level 
         cell_child_1 =  cell_coarse->GetChild(j1);
         N_child_1 = cell_child_1->GetN_Children();
         // finest level reached 
         if (N_child_1==0)
         {
           ii = cell_child_1->GetClipBoard();
           AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

           CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_1);
           n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
           // update vector which stores information on the length of the rows
           for(j=0;j<n1;j++)
           {
              k = TestNumbers[j];
              AuxPtr[k] += n2;
           } // endfor j
         }
         else
         {
            // second refinement level 
           for (j2=0;j2<N_child_1;j2++)
           {
             cell_child_2 =  cell_child_1->GetChild(j2);
             N_child_2 = cell_child_2->GetN_Children();
             // finest level reached
             if (N_child_2==0)
             {
                ii = cell_child_2->GetClipBoard();
                AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_2);
                n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                // update vector which stores information on the length of the rows
                for(j=0;j<n1;j++)
                {
                   k = TestNumbers[j];
                   AuxPtr[k] += n2;
                } // endfor j
             }
             else
             {
               // third refinement level 
               for (j3=0;j3<N_child_2;j3++)
               {
                 cell_child_3 =  cell_child_2->GetChild(j3);
                 N_child_3 = cell_child_3->GetN_Children();
                 // finest level reached
                 if (N_child_3==0)
                 {
                   ii = cell_child_3->GetClipBoard();
                   AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                   CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_3);
                   n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                   // update vector which stores information on the length of the rows
                   for(j=0;j<n1;j++)
                   {
                      k = TestNumbers[j];
                      AuxPtr[k] += n2;
                   } // endfor j
                 }
                 else
                 {                
                   // fourth refinement level 
                   for (j4=0;j4<N_child_3;j4++)
                   {
                     cell_child_4 =  cell_child_3->GetChild(j4);
                     N_child_4 = cell_child_4->GetN_Children();
                     // finest level reached
                     if (N_child_4==0)
                     {
                       ii = cell_child_4->GetClipBoard();
                       AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                       CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_4);
                       n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                       // update vector which stores information on the length of the rows
                       for(j=0;j<n1;j++)
                       {
                         k = TestNumbers[j];
                         AuxPtr[k] += n2;
                       } // endfor j
                     }
                     else
                     {
                       // fifth refinement level 
                       for (j5=0;j5<N_child_4;j5++)
                       {
                          cell_child_5 =  cell_child_4->GetChild(j5);
                          N_child_5 = cell_child_5->GetN_Children();
                          // finest level reached
                          if (N_child_5==0)
                          {
                             ii = cell_child_5->GetClipBoard();
                             AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                             CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_5);
                             n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                             // update vector which stores information on the length of the rows
                             for(j=0;j<n1;j++)
                             {
                               k = TestNumbers[j];
                               AuxPtr[k] += n2;
                             } // endfor j
                          }
			  else 
			  {
			     // sixth refinement level 
			     for (j6=0;j6<N_child_5;j6++)
			     {
				cell_child_6 =  cell_child_5->GetChild(j6);
				N_child_6 = cell_child_6->GetN_Children();
				// finest level reached
				if (N_child_6==0)
				{
				    ii = cell_child_6->GetClipBoard();
				    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
				    
				    CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_6);
				    n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
				    // update vector which stores information on the length of the rows
				    for(j=0;j<n1;j++)
				    {
					k = TestNumbers[j];
					AuxPtr[k] += n2;
				    } // endfor j
				}
				else
				{
				    exit(4711);
				}
			     } // endfor j6
			  }  // endfor else (j5) 
                       } // endfor j5
                     } // endfor else (j4)
                   } // endfor j4
                 } // endfor else (j3)
               } // endfor j3
             } // endfor else (j2)
           } // endfor j2
         } //endfor else (j1)
       } // endfor j1   
    } // endfor i
  }
  else
  {
    // test_level >= ansatz_level
    // number of mesh cells
    N_ = coll_fine->GetN_Cells();
    // loop over all mesh cells of the test space
    for(i=0;i<N_;i++)
    {
      // get cell i 
      cell_fine = coll_fine->GetCell(i);
      // get global number of d.o.f. of test and ansatz function 
      // which live on this mesh cell
      TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
      // get fe spaces on the mesh cell
      CurrentElement = TestSpace2D->GetFE2D(i, cell_fine);
      n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // find parent cell on coarse grid
      cell_parent = cell_fine;
      for (j1=0;j1<level_diff;j1++)
      {
        cell_tmp = cell_parent->GetParent();
        cell_parent = cell_tmp;
      }
      ii = cell_parent->GetClipBoard();
      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

      CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_parent);
      n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // update vector which stores information on the length of the rows
      for(j=0;j<n1;j++)
      {
        k = TestNumbers[j];
        AuxPtr[k] += n2;
      } // endfor j
    }
  }

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=N_Rows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  // cout << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    // cout << "i: " << i << "  l: " << l << endl;
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    // cout << AuxPtr[i+1] << endl;
  }

  // cout << "Upper bound: " << AuxPtr[N_Rows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows]; // upper bound for number of matrix entries 
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  // compute column indices
  N_Entries=0;
  if (test_level < ansatz_level)
  {
    // number of mesh cells
    N_ = coll_coarse->GetN_Cells();
    for(i=0;i<N_;i++)
    {
      cell_coarse = coll_coarse->GetCell(i);
      TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
      CurrentElement = TestSpace2D->GetFE2D(i, cell_coarse);
      n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

      // find all children of the coarse mesh cell 
      N_child = cell_coarse->GetN_Children();
      // first level of refinement
      for (j1=0;j1<N_child;j1++)
      {
        cell_child_1 =  cell_coarse->GetChild(j1);
        N_child_1 = cell_child_1->GetN_Children();
        // finest level reached
        if (N_child_1==0)
        {
          ii = cell_child_1->GetClipBoard();
          AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
          CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_1);
          n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
          // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
          for(j=0;j<n1;j++)
          {
            for(k=0;k<n2;k++)
            {
              m=AnsatzNumbers[k];
              n=TestNumbers[j];
              index=AuxPtr[n];
              l=KColAux[index];
              // check whether this column is already in this row
              while(l!=-1 && l!=m && index<AuxPtr[n+1])
              {
                index++; 
                l=KColAux[index];
              }
              if(l==-1)
              {
                 // this is a new column for this row
                KColAux[index]=m;
                N_Entries++;
              }
            } // endfor k
          } // endfor j
        }
        else
        {
          // second level of refinement
          for (j2=0;j2<N_child_1;j2++)
          {
            cell_child_2 =  cell_child_1->GetChild(j2);
            N_child_2 = cell_child_2->GetN_Children();
            // finest level reached
            if (N_child_2==0)
            {
              ii = cell_child_2->GetClipBoard();
              AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
              CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_2);
              n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
              // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
              for(j=0;j<n1;j++)
              {
                for(k=0;k<n2;k++)
                {
                  m=AnsatzNumbers[k];
                  n=TestNumbers[j];
                  index=AuxPtr[n];
                  l=KColAux[index];
                  // check whether this column is already in this row
                  while(l!=-1 && l!=m && index<AuxPtr[n+1])
                  {
                    index++; 
                    l=KColAux[index];
                  }
                  if(l==-1)
                  {
                    // this is a new column for this row
                    KColAux[index]=m;
                    N_Entries++;
                  }
                } // endfor k
              } // endfor j
            }
            else
            {
              // third level of refinement
              for (j3=0;j3<N_child_2;j3++)
              {
                cell_child_3 =  cell_child_2->GetChild(j3);
                N_child_3 = cell_child_3->GetN_Children();
                // finest level reached
                if (N_child_3==0)
                {
                  ii = cell_child_3->GetClipBoard();
                  AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
                  CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_3);
                  n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                  // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
                  for(j=0;j<n1;j++)
                  {
                    for(k=0;k<n2;k++)
                    {
                      m=AnsatzNumbers[k];
                      n=TestNumbers[j];
                      index=AuxPtr[n];
                      l=KColAux[index];
                      // check whether this column is already in this row
                      while(l!=-1 && l!=m && index<AuxPtr[n+1])
                      {
                        index++; 
                        l=KColAux[index];
                      }
                      if(l==-1)
                      {
                        // this is a new column for this row
                        KColAux[index]=m;
                        N_Entries++;
                      }
                    } // endfor k
                  } // endfor j
                }
                else
                {
                  // fourth level of refinement
                  for (j4=0;j4<N_child_3;j4++)
                  {
                    cell_child_4 =  cell_child_3->GetChild(j4);
                    N_child_4 = cell_child_4->GetN_Children();
                    // finest level reached
                    if (N_child_4==0)
                    {
                      ii = cell_child_4->GetClipBoard();
                      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
                      CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_4);
                      n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                      // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
                      for(j=0;j<n1;j++)
                      {
                        for(k=0;k<n2;k++)
                        {
                          m=AnsatzNumbers[k];
                          n=TestNumbers[j];
                          index=AuxPtr[n];
                          l=KColAux[index];
                          // check whether this column is already in this row
                          while(l!=-1 && l!=m && index<AuxPtr[n+1])
                          {
                            index++; 
                            l=KColAux[index];
                          }
                          if(l==-1)
                          {
                            // this is a new column for this row
                            KColAux[index]=m;
                            N_Entries++;
                          }
                        } // endfor k
                      } // endfor j
                    }
                    else
                    {
                      // fifth level of refinement
                      for (j5=0;j5<N_child_4;j5++)
                      {
                        cell_child_5 =  cell_child_4->GetChild(j5);
                        N_child_5 = cell_child_5->GetN_Children();
                        // finest level reached
                        if (N_child_5==0)
                        {
                          ii = cell_child_5->GetClipBoard();
                          AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
                          CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_child_5);
                          n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                          // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
                          for(j=0;j<n1;j++)
                          {
                            for(k=0;k<n2;k++)
                            {
                              m=AnsatzNumbers[k];
                              n=TestNumbers[j];
                              index=AuxPtr[n];
                              l=KColAux[index];
                              // check whether this column is already in this row
                              while(l!=-1 && l!=m && index<AuxPtr[n+1])
                              {
                                index++; 
                                l=KColAux[index];
                              }
                              if(l==-1)
                              {
                                // this is a new column for this row
                                KColAux[index]=m;
                                N_Entries++;
                              }
                            } // endfor k
                          } // endfor j
                        }
                        else
                        {
                          exit(4711);
                        } // end else (j5)
                      } // endfor j5                                        
                    } // end else (j4)
                  } // endfor j4                  
                } // end else (j3)
              } // endfor j3
            } // end else (j2)
          } // endfor j2
        } // end else (j1)
      } // endfor j1       
    } // endfor i
  } // endfor test_level < ansatz_level
  else
  {
    // test_level >= ansatz_level
    // number of mesh cells
    N_ = coll_fine->GetN_Cells();
    // loop over all mesh cells of the test space
    for(i=0;i<N_;i++)
    {
      // get cell i 
      cell_fine = coll_fine->GetCell(i);
      // get global number of d.o.f. of test and ansatz function 
      // which live on this mesh cell
      TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
      // get fe spaces on the mesh cell
      CurrentElement = TestSpace2D->GetFE2D(i, cell_fine);
      n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // find parent cell on coarse grid
      cell_parent = cell_fine;
      for (j1=0;j1<level_diff;j1++)
      {
        cell_tmp = cell_parent->GetParent();
        cell_parent = cell_tmp;
      }
      ii = cell_parent->GetClipBoard();
      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

      ii = cell_parent->GetClipBoard();
      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
      CurrentElement = AnsatzSpace2D->GetFE2D(ii, cell_parent);
      n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
      for(j=0;j<n1;j++)
      {
        for(k=0;k<n2;k++)
        {
          m=AnsatzNumbers[k];
          n=TestNumbers[j];
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m && index<AuxPtr[n+1])
          {
            index++; 
            l=KColAux[index];
          }
          if(l==-1)
          {
                                // this is a new column for this row
            KColAux[index]=m;
            N_Entries++;
          }
        } // endfor k
      } // endfor j
    } // endfor i
  } // end else

/*
  // check
  cout << endl;
  cout << "check" << endl;
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
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  { 
  cout << "Number of matrix entries: ";
  cout << N_Entries << endl;
  cout << endl;
  }
  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=N_Rows;
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
    } // endfor j
    RowPtr[i]=oldindex;
    // cout << setw(4) << i << " RowPtr[i]: " << RowPtr[i] << endl;
  } // endfor i
  RowPtr[N_Rows]=N_Entries;
  
/*
  // print out the whole matrix structure
  OutPut(endl);
  N_=N_Rows;
  for(i=0;i<N_;i++)
  {
     OutPut( RowPtr[i] << "---" << RowPtr[i+1]-1 << endl);
     OutPut("Rows: " << setw(4) << i << ": ");
     end=RowPtr[i+1];
     for(j=RowPtr[i];j<end;j++)
        OutPut( setw(4) << KCol[j]);
     OutPut(endl);
  }
*/
  
  // free KColAux
  delete [] KColAux;
  
#ifdef _MPI
  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif
  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
  cout << "Information on the stored matrix structure" << endl;
  cout << "Number of rows: " << N_Rows << endl;
  cout << "Number of columns: " << N_Columns << endl;
  cout << "Number of matrix entries: " << N_Entries << endl;
  }
}


/** generate the matrix structure, one space 1D and one 2D */
TStructure2D::TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace, int **ansatzcelljoints)
{
  int i,j,k,l,n,N_, n1,n2, m, Cell2D_No, IJoint;
  int *GlobalNumbers, N_JointDOF, *JointDOF, N_BaseFunct, *DOF;
  int *BeginIndex, *GlobalNumbers2D, *BeginIndex2D;
  int *Numbers, *nieb_Numbers;
  int N_Dirichlet, end, N_Hanging, Offset;
  int N_Inner, NE, nieb_i, nieb_e, nieb_n1;
  int *AuxPtr, *AuxPtr_delete, *HangingAuxPtr;
  int *KColAux, *HangingKColAux, *ansatzcell, *ansatzjoint;
  int index, oldindex, EdgeIntegrals=testspace->IsDGSpace();

#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  TCollection *Coll, *Coll_2D;

  TBaseCell *cell, *neigh, *Me;
  FE1D CurrentElement, CurrentNeighbour;
  FE2D FEId;
  TFEDesc2D *FeDesc;

  TestSpace1D = testspace;
  AnsatzSpace2D = ansatzspace;
  TestSpace2D = NULL;
  AnsatzSpace1D = NULL;
  TestMortarSpaceGlobNo = NULL;

  ansatzcell  = ansatzcelljoints[0];
  ansatzjoint = ansatzcelljoints[1];

  // all dof are treated as unknowns !!!
  // no boundary description is used so far!!!
  N_Rows = TestSpace1D->GetN_DegreesOfFreedom();
  N_Columns =TestSpace1D->GetN_DegreesOfFreedom();
  N_Entries = 0;

  AnsatzMortarSpaceGlobNo = new int[N_Rows];
  for(i=0;i<N_Rows;i++)
   AnsatzMortarSpaceGlobNo[i] = -1;

  // AuxPtr[i] will contain an upper bound for the number of matrix entries in row i
  l=N_Rows+1;
  AuxPtr = new int[l];
  memset(AuxPtr, 0, l*sizeof(int));

  GlobalNumbers=TestSpace1D->GetGlobalNumbers();
  BeginIndex=TestSpace1D->GetBeginIndex();

  GlobalNumbers2D=AnsatzSpace2D->GetGlobalNumbers();
  BeginIndex2D=AnsatzSpace2D->GetBeginIndex();

  // loop over all elements
  N_=TestSpace1D->GetN_Cells();

  Coll = TestSpace1D->GetCollection();
  Coll_2D = AnsatzSpace2D->GetCollection();

  // associate each cell with her number in the collection
  for(i=0;i<N_;i++)
   {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
   }

  for(i=0;i<N_;i++)
   {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    //======================================================================================
    Cell2D_No = ansatzcell[i];
    Me = Coll_2D->GetCell(Cell2D_No);
    IJoint = ansatzjoint[i];
    FEId = AnsatzSpace2D->GetFE2D(Cell2D_No, Me);

    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF = FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers2D + BeginIndex2D[Cell2D_No];
    //======================================================================================

    CurrentElement = TestSpace1D->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();

    if(N_JointDOF != n1 )
     {
      cout<< " N_JointDOF != n1 !! Mortar space and 1DSpace must be same Fespace " <<endl;
      exit(0);
     }

    for(j=0;j<n1;j++)
     AnsatzMortarSpaceGlobNo[Numbers[j]] =  DOF[JointDOF[j]];

    for(j=0;j<n1;j++)
    {
      k=Numbers[j];
      AuxPtr[k]+=n1;

     //DG additional entries necessary if integrals over the edges(points in 1D)
     // appear in the discretization
     if(EdgeIntegrals)
     {
      l=cell->GetN_Edges();                   // # edges
      for(m=0;m<l;m++)                        // for all edges
      {
       neigh = (cell->GetJoint(m))->GetNeighbour(cell);

       if(neigh)
       {
        n = neigh->GetClipBoard();
        CurrentNeighbour = TestSpace1D->GetFE1D(n, neigh);
        n2 = TFEDatabase2D::GetFE1D(CurrentNeighbour)->GetN_DOF();
        AuxPtr[k]+=n2;
       }//if(neigh)
      }
     }//if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)

    }                                             // for(j=0;j<n1;
   }                                               //  for(i=0;i<N

  N_Entries = 0;

//   for(i=0;i<=N_;i++)
//     cout << i << "   " << AuxPtr[i] << endl;
//   cout << endl;
// exit(0);

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
//         cout << AuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "AuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << AuxPtr[i] << endl;
  cout << endl;
exit(0);
  */

  // if(TDatabase::ParamDB->SC_VERBOSE)
    // cout << "Upper bound: " << AuxPtr[N_Rows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows]; // upper bound for number of matrix entries
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  N_=TestSpace1D->GetN_Cells();
  Coll = TestSpace1D->GetCollection();

  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = TestSpace1D->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)                             // test space
    {
      n=Numbers[j];
      for(k=0;k<n1;k++)                           // ansatz space
      {
        m=Numbers[k];
        // this  node is a real node (inner or Neumann)
        index=AuxPtr[n];
        l=KColAux[index];
        // check whether this column is already in this row
        while(l!=-1 && l!=m)
        {
          index++;
          l=KColAux[index];
        }
        if(l==-1)
        {
          // this is a new column for this row
          KColAux[index]=m;
          N_Entries++;
        }
      }                                           // endfor k
      
     // DG part
     if(EdgeIntegrals)
     {
      NE=cell->GetN_Edges();                   // # edges
      for(nieb_e=0;nieb_e<NE;nieb_e++)                        // for all edges
      {
       neigh = (cell->GetJoint(nieb_e))->GetNeighbour(cell);

       if(neigh)
       {
        nieb_i = neigh->GetClipBoard();
        CurrentNeighbour = TestSpace1D->GetFE1D(nieb_i, neigh);
        nieb_Numbers=GlobalNumbers+BeginIndex[nieb_i];
        nieb_n1 = TFEDatabase2D::GetFE1D(CurrentNeighbour)->GetN_DOF();
 
        for(k=0;k<nieb_n1;k++)                           // ansatz space
         {
          m=nieb_Numbers[k];
          // this  node is a real node (inner or Neumann)
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
           index++;
           l=KColAux[index];
          }
         if(l==-1)
          {
           // this is a new column for this row
           KColAux[index]=m;
           N_Entries++;
          }
         } // for(k=0;k<nieb_n1;k++)	
       }  //if(neigh)
      } //  for(nieb_e=0;nieb_e<NE;nieb_e++)
     }//if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
    }                                             // endfor j
  }                                               // for(i=0;

  /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=N_Rows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index;j++) //  && KColAux[j]!=-1
  cout << setw(4) << KColAux[j];
  cout << endl;
  }
  exit(0);
  */

  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=N_Rows;
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
    //     cout << setw(4) << i << " RowPtr[i]: " << RowPtr[i] << endl;
  }

//  cout << "index: " << index << endl;
//   cout << "RowPtr[N_]: " << RowPtr[N_] << endl;
  Offset=index-RowPtr[N_];
  for(i=0,j=N_Rows;i<=N_Hanging;i++,j++)
   {
    RowPtr[j]+=Offset;
    // cout << setw(4) << j << " RowPtr[j]: " << RowPtr[j] << endl;
   }

  /*
    // print out the whole matrix structure
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

//    for(i=0;i<N_Rows;i++)
//     cout << i<< " AnsatzMortarSpaceGlobNo: " << AnsatzMortarSpaceGlobNo[i] << endl;

  delete [] KColAux;

#ifdef _MPI
  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>1)
#endif 
  {
  cout << "Information on the stored matrix structure" << endl;
  cout << "Number of rows: " << N_Rows << endl;
  cout << "Number of columns: " << N_Columns << endl;
  cout << "Number of matrix entries: " << N_Entries << endl;
  }

// exit(0);
}


TStructure2D::TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace,  
                           TNonMortarData *NonMortarFEData)
{
  int i,j,k,l,n, N_, n1,n2, m, Cell2D_No, IJoint;
  int *GlobalNumbers, N_JointDOF, *JointDOF, N_BaseFunct, *DOF;
  int *BeginIndex, *GlobalNumbers2D, *BeginIndex2D;
  int *Numbers, *nieb_Numbers, N_NonMotEdges;
  int N_Dirichlet, end, N_Hanging, Offset;
  int N_Inner, NE, nieb_i, nieb_e, nieb_n1;
  int *AuxPtr, *AuxPtr_delete, *HangingAuxPtr;
  int *KColAux, *HangingKColAux, N_ColAnsatzs, dof, ColIndex;
  int index, oldindex, EdgeIntegrals=testspace->IsDGSpace();
  int *N_DofPerEdge;

  TCollection *Coll;

  TBaseCell *cell, *neigh, *Me;
  FE1D CurrentElement, CurrentNeighbour;
  FE2D FEId;
  TFEDesc2D *FeDesc;

#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  TestSpace1D = testspace;
  AnsatzSpace2D = ansatzspace;
  TestSpace2D = NULL;
  AnsatzSpace1D = NULL;
  TestMortarSpaceGlobNo = NULL;

  // all dof are treated as unknowns !!!
  // no boundary description is used so far!!!
  N_Rows = TestSpace1D->GetN_DegreesOfFreedom();
  N_Columns = NonMortarFEData->N_NonMortaDofs;
  N_Entries = 0;

  AnsatzNonMortarSpaceGlobNo = new int[N_Columns];

  // AuxPtr[i] will contain an upper bound for the number of matrix entries in row i
  l=N_Rows+1;
  AuxPtr = new int[l];
  memset(AuxPtr, 0, l*sizeof(int));

  GlobalNumbers=TestSpace1D->GetGlobalNumbers();
  BeginIndex=TestSpace1D->GetBeginIndex();

  GlobalNumbers2D= NonMortarFEData->EdgeNonMotLocGlobalNo;
  BeginIndex2D= NonMortarFEData->EdgeNonMotBeginIndex;
  N_DofPerEdge = NonMortarFEData->N_DofPerEdge;

  // loop over all elements
  N_=TestSpace1D->GetN_Cells();

  Coll = TestSpace1D->GetCollection();

  for(i=0;i<N_;i++)
   {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = TestSpace1D->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();

    N_ColAnsatzs = N_DofPerEdge[i];

    for(j=0;j<n1;j++)
      AuxPtr[Numbers[j]]+=N_ColAnsatzs;

   }//  for(i=0;i<N_;i++)

  N_Entries = 0;

//   for(i=0;i<N_Rows;i++)
//     cout << i << "   " << AuxPtr[i] << endl;
//   cout << endl;
// // exit(0);

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
exit(0);
  */
  
  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows]; // upper bound for number of matrix entries
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  N_=TestSpace1D->GetN_Cells();
  Coll = TestSpace1D->GetCollection();

  for(i=0;i<N_;i++)
   {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = TestSpace1D->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();

    N_ColAnsatzs = N_DofPerEdge[i]; // ansatz space  
    DOF = GlobalNumbers2D + BeginIndex2D[i];
       
    for(j=0;j<n1;j++) // test space
    {
     n=Numbers[j];

     for(l=0;l<N_ColAnsatzs;l++) // insatz space  
      {          
       dof = DOF[l];
       
       index=AuxPtr[n];
       m=KColAux[index];    
       
        // check whether this column is already in this row
        while(m!=-1 && m!=dof)
         {
          index++;
          m=KColAux[index];
         }

        if(m==-1)
        {
          // this is a new column for this row
          KColAux[index]=dof;
          N_Entries++;
        }
      }
    } //  for(j=0;j<n1;j++
   }  // for(i=0;  

   
    /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=N_Rows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index;j++) //  && KColAux[j]!=-1
  cout << setw(4) << KColAux[j];
  cout << endl;
  }
//   exit(0);
  */ 
   
  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=N_Rows;
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
  }  
   
//  cout << "index: " << index << endl;
//   cout << "RowPtr[N_]: " << RowPtr[N_] << endl;
  Offset=index-RowPtr[N_];
  for(i=0,j=N_Rows;i<=N_Hanging;i++,j++)
   {
    RowPtr[j]+=Offset;
    // cout << setw(4) << j << " RowPtr[j]: " << RowPtr[j] << endl;
   }

  /*
    // print out the whole matrix structure
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

  //    for(i=0;i<N_Rows;i++)
//     cout << i<< " AnsatzMortarSpaceGlobNo: " << AnsatzMortarSpaceGlobNo[i] << endl;

  delete [] KColAux;

#ifdef _MPI
  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>1)
#endif 
  {
  cout << "Information on the stored matrix structure" << endl;
  cout << "Number of rows: " << N_Rows << endl;
  cout << "Number of columns: " << N_Columns << endl;
  cout << "Number of matrix entries: " << N_Entries << endl;
  }


} // TStructure2D


TStructure2D::TStructure2D(TFESpace2D *testspace, TFESpace1D *ansatzspace,  
                           TNonMortarData *NonMortarFEData)
{
  int i,j,k,l,n, N_, n1,n2, m, Cell2D_No, IJoint;
  int *GlobalNumbers, N_JointDOF, *JointDOF, N_BaseFunct, *DOF;
  int *BeginIndex, *GlobalNumbers2D, *BeginIndex2D;
  int *Numbers, *nieb_Numbers, N_NonMotEdges;
  int N_Dirichlet, end, N_Hanging, Offset;
  int N_Inner, NE, nieb_i, nieb_e, nieb_n1;
  int *AuxPtr, *AuxPtr_delete, *HangingAuxPtr;
  int *KColAux, *HangingKColAux, N_ColAnsatzs, dof, ColIndex;
  int index, oldindex, EdgeIntegrals=testspace->IsDGSpace();
  int *N_DofPerEdge;
  
  
  TCollection *Coll;

  TBaseCell *cell, *neigh, *Me;
  FE1D CurrentElement, CurrentNeighbour;
  FE2D FEId;
  TFEDesc2D *FeDesc;
  
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  TestSpace1D = NULL;
  AnsatzSpace2D = NULL;
  TestSpace2D = testspace;
  AnsatzSpace1D = ansatzspace;
  TestMortarSpaceGlobNo = NULL; 
  AnsatzNonMortarSpaceGlobNo = NULL;
  
  // all dof are treated as unknowns !!!
  // no boundary description is used so far!!!
  N_Rows = NonMortarFEData->N_NonMortaDofs ;
  N_Columns = AnsatzSpace1D->GetN_DegreesOfFreedom() ;
  N_Entries = 0;

  TestNonMortarSpaceGlobNo = new int[N_Rows];
  
  // AuxPtr[i] will contain an upper bound for the number of matrix entries in row i
  l=N_Rows+1;
  AuxPtr = new int[l];
  memset(AuxPtr, 0, l*sizeof(int));

  GlobalNumbers=AnsatzSpace1D->GetGlobalNumbers();
  BeginIndex=AnsatzSpace1D->GetBeginIndex();

  GlobalNumbers2D= NonMortarFEData->EdgeNonMotLocGlobalNo;
  BeginIndex2D= NonMortarFEData->EdgeNonMotBeginIndex;
  N_DofPerEdge = NonMortarFEData->N_DofPerEdge;

  // loop over all elements
  N_=AnsatzSpace1D->GetN_Cells();

  Coll = AnsatzSpace1D->GetCollection();

  for(i=0;i<N_;i++)
   {
    cell = Coll->GetCell(i);
    CurrentElement = AnsatzSpace1D->GetFE1D(i, cell);
    N_ColAnsatzs = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();

    DOF = GlobalNumbers2D+BeginIndex2D[i];
    n1 = N_DofPerEdge[i];

    for(j=0;j<n1;j++)
      AuxPtr[DOF[j]]+=N_ColAnsatzs;

   }//  for(i=0;i<N_;i++)  

  N_Entries = 0;

//   for(i=0;i<N_Rows;i++)
//     cout << i << "   " << AuxPtr[i] << endl;
//   cout << endl;
// // exit(0);

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
exit(0);
  */
  
  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows]; // upper bound for number of matrix entries
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  N_=AnsatzSpace1D->GetN_Cells();

  for(i=0;i<N_;i++)
   {
    cell = Coll->GetCell(i);
    CurrentElement = AnsatzSpace1D->GetFE1D(i, cell);
    
    N_ColAnsatzs = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[i];
 
    n1 = N_DofPerEdge[i]; // test space  
    Numbers = GlobalNumbers2D + BeginIndex2D[i];
       
    for(j=0;j<n1;j++) // test space
    {
     n=Numbers[j];

     for(l=0;l<N_ColAnsatzs;l++) // insatz space  
      {          
       dof = DOF[l];
       
       index=AuxPtr[n];
       m=KColAux[index];    
       
        // check whether this column is already in this row
        while(m!=-1 && m!=dof)
         {
          index++;
          m=KColAux[index];
         }

        if(m==-1)
         {
          // this is a new column for this row
          KColAux[index]=dof;
          N_Entries++;
        }
      }
    } //  for(j=0;j<n1;j++
   }  // for(i=0;    
  
     /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=N_Rows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index;j++) //  && KColAux[j]!=-1
  cout << setw(4) << KColAux[j];
  cout << endl;
  }
//   exit(0);
  */  

  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=N_Rows;
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
  }  
   
//  cout << "index: " << index << endl;
//   cout << "RowPtr[N_]: " << RowPtr[N_] << endl;
  Offset=index-RowPtr[N_];
  for(i=0,j=N_Rows;i<=N_Hanging;i++,j++)
   {
    RowPtr[j]+=Offset;
    // cout << setw(4) << j << " RowPtr[j]: " << RowPtr[j] << endl;
   }

  /*
    // print out the whole matrix structure
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

  //    for(i=0;i<N_Rows;i++)
//     cout << i<< " AnsatzMortarSpaceGlobNo: " << AnsatzMortarSpaceGlobNo[i] << endl;

  delete [] KColAux;

#ifdef _MPI
  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>1)
#endif 
  {
  cout << "Information on the stored matrix structure" << endl;
  cout << "Number of rows: " << N_Rows << endl;
  cout << "Number of columns: " << N_Columns << endl;
  cout << "Number of matrix entries: " << N_Entries << endl;
  }
     
     
} // TStructure2D




/** destructor */
TStructure2D::~TStructure2D()
{
}
