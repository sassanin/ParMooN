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
// @(#)Structure3D.C        1.10 11/24/99
// 
// Class:       TStructure3D
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

#include <FEDatabase3D.h>
#include <Structure3D.h>
#include <string.h>
#include <Database.h>

/** generate the matrix structure, both spaces are 3D */
TStructure3D::TStructure3D(TFESpace3D *testspace, TFESpace3D *ansatzspace)
 : TStructure()
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
  FE3D CurrentElement; 

  if(testspace->GetCollection() != ansatzspace->GetCollection())
  {
    return;
  }

  coll = testspace->GetCollection();

  // test space and ansatz space differ
  TestSpace3D = testspace;
  AnsatzSpace3D = ansatzspace;
  TestSpace2D = NULL;
  AnsatzSpace2D = NULL;

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

  N_Rows = TestSpace3D->GetN_DegreesOfFreedom();
  N_Columns = AnsatzSpace3D->GetN_DegreesOfFreedom();

  l = N_Rows + 1;
  AuxPtr = new int[l];
  memset(AuxPtr, 0, l*SizeOfInt);

  N_ = coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell = coll->GetCell(i);
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];

    CurrentElement = TestSpace3D->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    CurrentElement = AnsatzSpace3D->GetFE3D(i, cell);
    n2 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      k = TestNumbers[j];
      AuxPtr[k] += n2;
    } // endfor j
  } // endfor i

#ifdef __MORTAR__
#ifdef __ADD_LINK__
  // reserve space for additional link term
  TIt_Mortar *It1 = (TIt_Mortar *) TDatabase::IteratorDB[It_Mortar1];
  TIt_Mortar *It2 = (TIt_Mortar *) TDatabase::IteratorDB[It_Mortar2];
  TBaseCell *CellM, *CellNM;
  FE3D CurrElementID;
  TFE3D *CurrElement;
  TFEDesc3D *CurrDesc;
  int lev, infoNM, infoM;
  double lamM0, lamM1 = 0, lamNM0, lamNM1, startX, startY, endX, endY;
  double delX, delY;
  bool NewMortarSideEle;
  int **J_DOF_NM, N_J_DOF_NM, **J_DOF_M, N_J_DOF_M;
  int indexM = 0, indexNM;

  // set ClipBoard to number in collection
  N_ = TestSpace3D->GetN_Cells();
  coll = TestSpace3D->GetCollection();
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
        cerr << "Error in MatrixStructure3D: cell out of collection!!!"
             << endl;
        exit(-3);
      }

      indexNM = CellNM->GetClipBoard();
      CurrElementID = AnsatzSpace3D->GetFE3D(indexNM, CellNM);
      CurrElement = TFEDatabase3D::GetFE3D(CurrElementID);
      CurrDesc = CurrElement->GetFEDesc3D();

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
          cerr << "Error in MatrixStructure3D: cell out of collection!!!"
               << endl;
          exit(-3);
        }

        indexM = CellM->GetClipBoard();
        CurrElementID = AnsatzSpace3D->GetFE3D(indexM, CellM);
        CurrElement = TFEDatabase3D::GetFE3D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc3D();

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

  N_Entries=0;
  N_ = coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];

    cell = coll->GetCell(i);

    CurrentElement = TestSpace3D->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    CurrentElement = AnsatzSpace3D->GetFE3D(i, cell);
    n2 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

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
          index++; l=KColAux[index];
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
        cerr << "Error in MatrixStructure3D: cell out of collection!!!"
             << endl;
        exit(-3);
      }

      indexNM = CellNM->GetClipBoard();
      CurrElementID = AnsatzSpace3D->GetFE3D(indexNM, CellNM);
      CurrElement = TFEDatabase3D::GetFE3D(CurrElementID);
      CurrDesc = CurrElement->GetFEDesc3D();

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
          cerr << "Error in MatrixStructure3D: cell out of collection!!!"
               << endl;
          exit(-3);
        }

        indexM = CellM->GetClipBoard();
        CurrElementID = AnsatzSpace3D->GetFE3D(indexM, CellM);
        CurrElement = TFEDatabase3D::GetFE3D(CurrElementID);
        CurrDesc = CurrElement->GetFEDesc3D();

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
  
  // cout << "Number of matrix entries: ";
  // cout << N_Entries << endl;
  // cout << endl;

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
  delete KColAux;

#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif
  {
  OutPut("Information on the stored matrix structure" << endl);
  OutPut("Number of rows: " << N_Rows << endl);
  OutPut("Number of columns: " << N_Columns << endl);
  OutPut("Number of matrix entries: " << N_Entries << endl);
 }
}
