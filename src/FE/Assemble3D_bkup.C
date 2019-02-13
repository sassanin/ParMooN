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
// %W% %G%
// 
// Class:       TAssemble3D
// Purpose:     bilinear form (discretized and stabilized assemble)
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <DefineParams.h>
#include <MooNMD_Io.h>
#include <Assemble3D.h>
#include <Enumerations.h>
#include <Matrix3D.h>
#include <AuxParam3D.h>
#include <DiscreteForm3D.h>
#include <TNSE3D_FixPo.h>
#include <Joint.h>
#include <BoundFace.h>
#include <FEDatabase3D.h>
#include <SquareMatrix3D.h>
#include <HNDesc.h>
#include <HangingNode.h>
#include <Database.h>
#include <Convolution.h>
#include <Edge.h>
#include <string.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int compare_int(const void *a, const void *b)
{
  return (*(int *)a) - (*(int *)b);
}

void static PrintAllMat(int N_SqMatrices, TSquareMatrix3D **SqMatrices,
                        int N_Matrices, TMatrix3D **RecMatrices,
                        int N_Rhs, double **Rhs, TFESpace3D **FeRhs)
{
 int i, j, k, end, *RowPtr, *ColInd, N_Rows, ActiveBound;
 bool print = false;
 double *Entries, *RHS, temp1=0,temp2=0;
 
 TFESpace3D *fespace;  
  
  for(k=0;k<N_SqMatrices;k++)
  {
    OutPut(endl);
    cout << "sqmatrix: " << k << endl;
    RowPtr = SqMatrices[k]->GetRowPtr();
    Entries = SqMatrices[k]->GetEntries();
    ColInd = SqMatrices[k]->GetKCol();
    N_Rows = SqMatrices[k]->GetN_Rows();
    
    fespace = SqMatrices[k]->GetFESpace();
    ActiveBound = fespace->GetActiveBound();
    for(i=0;i<ActiveBound;i++)
    {
     end=RowPtr[i+1];
     
     if(print)
       OutPut("\033[1;31ma"<<k<<"("<< i+1 << ", " << "j):"<<"\033[0m");
       
     for(j=RowPtr[i];j<end;j++){
       temp1 += Entries[j];
       temp2 += Entries[j]*Entries[j];
       
       if(print)
	 OutPut("\033[1;31m" << ColInd[j]+1 <<"=\033[0m"<<Entries[j] << ", ");
     }
    
     if(print)
      OutPut(" "<< endl);
    }
    
    if(print)
      cout << endl;
    
    printf("sum of all entries            :: %lf\n",temp1);
    printf("sum of squares of all entries :: %lf\n",temp2);
    
    temp1=0;
    temp2=0;
  } // endfor k
  
  for(k=0;k<N_Matrices;k++)
  {
    cout << endl;
    cout << "matrix: " << k << endl;
    RowPtr = RecMatrices[k]->GetRowPtr();
    Entries = RecMatrices[k]->GetEntries();
    ColInd = RecMatrices[k]->GetKCol();
    N_Rows = RecMatrices[k]->GetN_Rows();
    for(i=0;i<N_Rows;i++)
    {
      end=RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
      {
        // cout << j << endl;
	if(print){
	  OutPut("b"<<k<<"("<< i+1 << "," << ColInd[j]+1 << " )=  ");
	  OutPut(Entries[j] << ";" << endl);
	}	
      }
    }
    cout << endl;
  } // endfor k

  for(k=0;k<N_Rhs;k++)
  {
    cout << "rhs: " << k << endl;
    N_Rows = FeRhs[k]->GetN_DegreesOfFreedom();
    RHS=Rhs[k];
    for(i=0;i<N_Rows;i++){
      if(print)
	cout << setw(5) << i << setw(20) << RHS[i] << endl;
    }
  }
  
  exit(0);
} //TAssembleMat3D::PrintAllMat()



void Assemble3D(int n_fespaces, TFESpace3D **fespaces,
                int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                int n_matrices, TMatrix3D **matrices,
                int n_rhs, double **rhs, TFESpace3D **ferhs,
                TDiscreteForm3D *DiscreteForm3D,
                BoundCondFunct3D **BoundaryConditions,
                BoundValueFunct3D **BoundaryValues,
                TAuxParam3D *Parameters)
{
  double hK;
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,k,l,l1,l2,l3,n,m, N_LocalUsedElements,ij,N_Vertex;
  int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints;
  int Used[N_FEs3D];
  int *N_BaseFunct;
  int N_Rows;
  BaseFunct3D *BaseFuncts;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  FE3D TestElement, AnsatzElement;
  QuadFormula2D FaceQuadFormula;
  TQuadFormula2D *qf2;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  TBoundFace *boundface;
  //  TIsoBoundEdge *isoboundedge;
  int **GlobalNumbers, **BeginIndex;
  int **RhsGlobalNumbers, **RhsBeginIndex;
  int **TestGlobalNumbers, **TestBeginIndex;
  int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
  TFE3D *ele;
  TFEDesc3D *FEDesc_Obj;
  double *weights, *xi, *eta, *zeta;
  double *t, *s;
  double x, y, z;
  double xf, yf, zf;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *Param[MaxN_QuadPoints_3D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts3D];
  BaseFunct3D LocBF[N_BaseFuncts3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, begin, end, last, middle;
  int *TestDOF, *AnsatzDOF;
  double *Entries;
  int *ColInd, *RowPtr;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp3D *BoundComp;
  double t0, t1, t2;
  int comp;
  BoundCond Cond0, Cond1;
  BoundCondFunct3D *BoundaryCondition;
  BoundValueFunct3D *BoundaryValue;
  TNodalFunctional3D *nf;
  RefTrans3D reftrans;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_BaseFunctions3D];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *FaceWeights, *theta;
  double x0, x1, y0, y1, z0, z1, hE;
  double **JointValues, *JointValue;
  // static bool *SecondDer = NULL;
  bool *SecondDer = NULL;
  double Param1[4], Param2[4];
  double LinComb[4];
 
  const int *TmpFV, *TmpLen, *EdgeVertex, *TmpFE, *ETmpLen;
  int MaxLen, EMaxLen, N_FaceEdges;
  double xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3;
  double nx, ny, nz;
  double t11,t22,time=0.0;

  double time1, time2, time_all, time_total;
  int *int_ptr, N_Edges, N_VertInCell, dof;

  TEdge *edge;
  TVertex *Vert; 
  
  bool InnerBoundary, OuterBoundary;

  time_total = GetTime();
  time_all = 0;

// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
  if(n_sqmatrices)
  {
    GlobalNumbers  = new int*    [n_sqmatrices];
    BeginIndex     = new int*    [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];
    
    for(i=0;i<n_sqmatrices;i++)
    {
      fespace          = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i]    = fespace->GetBeginIndex();

      j = sqmatrices[i]->GetHangingN_Entries();
//       cout << "N_Hanging: " << j << endl;
      HangingEntries[i] = new double [j];
      memset(HangingEntries[i], 0, SizeOfDouble*j);
    } // endfor
  } // endif n_sqmatrices

  if(n_matrices)
  {
    TestGlobalNumbers = new int* [n_matrices];
    AnsatzGlobalNumbers = new int* [n_matrices];
    TestBeginIndex = new int* [n_matrices];
    AnsatzBeginIndex = new int* [n_matrices];
    for(i=0;i<n_matrices;i++)
    {
      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetTestSpace();
      TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
      TestBeginIndex[i] = fespace->GetBeginIndex();
  
      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetAnsatzSpace();
      AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
      AnsatzBeginIndex[i] = fespace->GetBeginIndex();
    } // endfor
  } // endif n_matrices

  if(n_rhs)
  {
    HangingRhs       = new double* [n_rhs];
    RhsBeginIndex    = new int*    [n_rhs];
    RhsGlobalNumbers = new int*    [n_rhs];
    
    for(i=0;i<n_rhs;i++)
    {
      fespace             = ferhs[i];
      RhsBeginIndex[i]    = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
  
      j = fespace->GetN_Hanging();
//       cout << "N_Hanging: " << j << endl;
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    } // endfor

    LocRhs    = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions3D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions3D;

  } // endif n_rhs

  N_Parameters = Parameters->GetN_Parameters();
//   cout << "N_Parameters: " << N_Parameters << endl;
  if(N_Parameters)
  {
//     cout<<"asdasdasd"<<endl;
    aux = new double [MaxN_QuadPoints_3D*N_Parameters];
    for(j=0;j<MaxN_QuadPoints_3D;j++)
      Param[j] = aux + j*N_Parameters;
  }
  // 20 <= number of term in bilinear form
  // DUE NOTE CHANGE 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  if(N_AllMatrices)
  {
    aux      = new double  [N_AllMatrices*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions3D];
    
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions3D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions3D;

    LocMatrices = new double** [N_AllMatrices];
    
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions3D;
  } // endif N_AllMatrices

//   SecondDer = DiscreteForm3D->GetNeeds2ndDerivatives();
  
  if(SecondDer == NULL)
  {
    SecondDer = new bool[n_fespaces];
    for(i=0;i<n_fespaces;i++)
      SecondDer[i] = false;
  }
  
  // all spaces use same Coll
  Coll = fespaces[0]->GetCollection();
  N_Cells = Coll->GetN_Cells();
//   cout << "N_Cells: " << N_Cells << endl;

  // reset clipboards of all neighbours to -1
// #ifdef _OPENMP
//  omp_set_num_threads(TDatabase::ParamDB->OMPNUMTHREADS);  
// #pragma omp parallel default(shared) private(i,k,j,cell,neigh,N_Edges,N_VertInCell)
//   {
//     #pragma omp for schedule(guided) 
// #endif
    for(i=0;i<N_Cells;i++)
      {
	cell = Coll->GetCell(i);
	k = cell->GetN_Joints();
	for(j=0;j<k;j++)
	{
	  neigh=cell->GetJoint(j)->GetNeighbour(cell);
	  if(neigh) neigh->SetClipBoard(-1);
	}
	cell->SetClipBoard(-1);

#ifdef _MPI       
	  N_Edges=cell->GetN_Edges();
	  for(j=0;j<N_Edges;j++)
	  (cell->GetEdge(j))->SetClipBoard(-1);
	  
	  N_VertInCell = cell->GetN_Vertices();
	  for(j=0;j<N_VertInCell;j++)
	  (cell->GetVertex(j))->SetClipBoard(-1); 
#endif     
      } // endfor i
// #ifdef _OPENMP
//       #pragma omp for schedule(guided)
// #endif
      for(i=0;i<N_Cells;i++)
	Coll->GetCell(i)->SetClipBoard(i);
      
//       cout<<"asd-------------------"<<endl;
// #ifdef _OPENMP
//   }
// #endif  
// return;


// ########################################################################
// loop over all cells
// ########################################################################
//hybrid
 double temp=0;

#ifdef _OPENMP
 TBaseCell *temp_cell;     
 omp_set_num_threads(TDatabase::ParamDB->OMPNUMTHREADS);  
#pragma omp parallel default(shared) firstprivate(i,k,j,m,l,n,cell,hK,CurrentElement,LocalUsedElements,LocN_BF,\
                                             LocBF, reftrans, N_Points, xi, eta, zeta, X, Y, Z, Param,\
                                             ij, N_Vertex, weights, AbsDetjk, AuxArray, N_, \
                                             DOF, MatrixRow,edge,\
                                             l2,end,TestElement,AnsatzElement,TestDOF,AnsatzDOF,N_Test,\
                                             N_Ansatz,local_rhs,EdgeDOF,N_EdgeDOF,\
                                             ele, nf,int_ptr,\
                                             FEDesc_Obj,N_Joints,joint,InnerBoundary,\
                                             OuterBoundary,TmpFV, TmpLen, MaxLen,TmpFE, ETmpLen, EMaxLen,\
                                             xf,yf,zf,t0,l1,Cond0,t,s,LinComb,FunctionalValues,PointValues,\
					     FaceQuadFormula,qf2,FaceWeights,JointValues,\
					     xc1,xc2,xc3,yc1,yc2,yc3,zc1,zc2,zc3,nx,ny,nz,t2,JointValue,t1,l3,\
					     N_FaceEdges,dof,Vert,RhsBeginIndex,BeginIndex,\
					     N_Edges,N_VertInCell,N_AllMatrices,n_rhs,LocMatrices,LocRhs,fespace,\
					     RowPtr,ColInd,Entries,CurrentHangingEntries,RHS,CurrentHangingRhs,SecondDer,\
					     Matrix,HangingRowPtr,HangingColInd,ActiveBound,DirichletBound,aux,Matrices,\
					     EdgeVertex,Param1,Param2,x0, x1, y0, y1, z0, z1, hE,theta,N_LinePoints,\
					     EdgePoints,N_EdgePoints,BoundaryValue,BoundaryCondition,Cond1,comp,BoundComp,\
					     Coupling,v,HNDescr_Obj,HNDescr,hn,HangingNodes,HangingEntries,HangingRhs,begin,\
					     last,middle,righthand,x, y, z,AnsatzGlobalNumbers,AnsatzBeginIndex,\
					     TestGlobalNumbers,TestBeginIndex,RhsGlobalNumbers,GlobalNumbers,boundface,\
					     neigh,BaseFuncts,N_Rows,N_BaseFunct,Used,N_Cells,N_Parameters,N_Hanging,\
					     N_LocalUsedElements,DiscreteForm3D)

  {
//     cout<<"TDatabase::ParamDB->OMPNUMTHREADS"<<TDatabase::ParamDB->OMPNUMTHREADS<<endl;
     #pragma omp for schedule(guided)
#endif
      for(i=0;i<N_Cells;i++)
      {
	cell = Coll->GetCell(i);
	switch (TDatabase::ParamDB->CELL_MEASURE)
	{
	  case 0: // diameter
	    hK = cell->GetDiameter();
	    //OutPut("diameter " << hK <<endl);
	    break;

	  case 2: // shortest edge
	    hK = cell->GetShortestEdge();
	    //OutPut("shortest edge " << hK <<endl);
	    break;
	  case 1: // with reference map
	  case 3: // measure
	    hK = cell->GetMeasure();
	    hK = pow(hK,1.0/3.0);
	    //OutPut("measure " << hK <<endl);
	    //exit(1);
	    break;
	  default: // diameter
	    hK = cell->GetDiameter();
	    break;
	}//switch (TDatabase::ParamDB->CELL_MEASURE)
	
//       } // endfor i

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
      for(j=0;j<n_fespaces;j++)
      {
	CurrentElement       = fespaces[j]->GetFE3D(i, cell);
	LocalUsedElements[j] = CurrentElement;
	LocN_BF[j]           = N_BaseFunct[CurrentElement];
	LocBF[j]             = BaseFuncts[CurrentElement];
      }

      N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
#ifdef _OPENMP
#pragma omp critical
#endif
    //OutPut("CELL " << i << endl);
      reftrans=TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                                    Coll, cell, SecondDer,
                                    N_Points, xi, eta, zeta, weights,
                                    X, Y, Z, AbsDetjk);
#ifdef _OPENMP
#pragma omp critical
#endif      
      Parameters->GetParameters(N_Points, cell, i, xi, eta, zeta,
                              X, Y, Z, Param); 
			      
      if ((TDatabase::ParamDB->DISCTYPE == SDFEM)
        || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
        || (TDatabase::ParamDB->CELL_MEASURE == 4)
        || (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 105)
        || (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 108)
        || (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 109))
      {
        N_Vertex = cell->GetN_Vertices();
        for (ij=0;ij<N_Vertex;ij++)
        {
#ifdef _OPENMP
#pragma omp critical
#endif
	  {
	    TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
	    TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
	    TDatabase::ParamDB->INTERNAL_VERTEX_Z[ij] = cell->GetVertex(ij)->GetZ();
	  }
	}
        if (N_Vertex==4)
          TDatabase::ParamDB->INTERNAL_VERTEX_X[4] = -4711;
        TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
      }//if ((TDatabase::ParamDB->DISCTYP
#ifdef _OPENMP
#pragma omp critical
#endif
      if(DiscreteForm3D)
        DiscreteForm3D->GetLocalForms(N_Points, weights, AbsDetjk, 
                                    hK, X, Y, Z,
                                    LocN_BF, LocBF,
                                    Param, AuxArray,
                                    cell,
                                    N_AllMatrices, n_rhs,
                                    LocMatrices, LocRhs);

    // ####################################################################
    // add local matrices to global matrices (ansatz == test)
    // ####################################################################
    for(j=0;j<n_sqmatrices;j++)
    {
      // find space for this bilinear form
      fespace        = sqmatrices[j]->GetFESpace();
      CurrentElement = fespace->GetFE3D(i, cell);
      N_             = N_BaseFunct[CurrentElement];

      Matrix  = Matrices+j*MaxN_BaseFunctions3D;
      Entries = sqmatrices[j]->GetEntries();
      RowPtr  = sqmatrices[j]->GetRowPtr();
      ColInd  = sqmatrices[j]->GetKCol();

      CurrentHangingEntries = HangingEntries[j];
      HangingRowPtr         = sqmatrices[j]->GetHangingRowPtr();
      HangingColInd         = sqmatrices[j]->GetHangingKCol();

      ActiveBound    = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF            = GlobalNumbers[j] + BeginIndex[j][i];

      // add local matrix to global
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          // for all dof
          for(k=0;k<N_;k++)
          {
            // for all columns
            l2 = 0;
            end=RowPtr[l+1];
            for(n=RowPtr[l];n<end;n++)
            {
              // if column of current dof found
              if(DOF[k] == ColInd[n])
              {
                // cout << m << "   " << k << endl << n << endl;
#ifdef _OPENMP
                #pragma omp critical
#endif
		{
		  Entries[n] += MatrixRow[k];
		  temp += MatrixRow[k];
		  l2 = 1;
		}
		break;
              } // endif
            } // endfor n
            if(l2 == 0)
              cout << "(assemble3d.c) not found" << endl;
          } // endfor k
        } // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
            end = HangingRowPtr[l+1];
            for(n=HangingRowPtr[l];n<end;n++)
            {
              for(k=0;k<N_;k++)
              {
                if(DOF[k] == HangingColInd[n])
                {
#ifdef _OPENMP
		  #pragma omp critical
#endif
                  CurrentHangingEntries[n] += MatrixRow[k];
                  break;
                } // endif
              } // endfor k
            } // endfor n 
          }
          else
          {
            // Dirichlet node
            n=RowPtr[l];
            if(ColInd[n]==l)
            {
              Entries[n]=1.0;
            }
          }
        }
      } // endfor m
    } // endfor j


    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(j=0;j<n_matrices;j++)
    {
      TestElement = ((TFESpace3D *) matrices[j]->GetStructure()->
                    GetTestSpace())->GetFE3D(i, cell);
      AnsatzElement = ((TFESpace3D *) matrices[j]->GetStructure()->
                      GetAnsatzSpace())->GetFE3D(i, cell);

      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement << endl;
      // cout << "AnsatzElement: " << AnsatzElement << endl;

      N_Test = N_BaseFunct[TestElement];
      N_Ansatz = N_BaseFunct[AnsatzElement];

      Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions3D;

      Entries = matrices[j]->GetEntries();
      RowPtr = matrices[j]->GetRowPtr();
      ColInd = matrices[j]->GetKCol();

      TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
      AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];

      // add local matrix to global
      for(m=0;m<N_Test;m++)
      {
        l=TestDOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        for(k=0;k<N_Ansatz;k++)
        {
          end = RowPtr[l+1];
          for(n=RowPtr[l];n<end;n++)
          {
            if(AnsatzDOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
#ifdef _OPENMP
              #pragma omp critical
#endif
	      Entries[n] += MatrixRow[k];
              break;
            } // endif
          } // endfor n
        } // endfor k
      } // endfor m
    } // endfor j
//     time2 = GetTime();
//     time_all += time2-time1;

    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    //for(j=0;j<n_rhs;j++)
    for(j=0;j<-1;j++)
    {
      fespace = ferhs[j];
      ActiveBound = fespace->GetActiveBound();
      CurrentElement = fespace->GetFE3D(i, cell);

      N_ = N_BaseFunct[CurrentElement];

      local_rhs = righthand+j*MaxN_BaseFunctions3D; 
      RHS = rhs[j];
      CurrentHangingRhs = HangingRhs[j];
      // find space for this linear form
      
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

      // add local right-hand side to the global one
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
#ifdef _OPENMP
          #pragma omp critical
#endif
	  RHS[l] += local_rhs[m]; 
        } // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
#ifdef _OPENMP
	    #pragma omp critical
#endif
            CurrentHangingRhs[l] += local_rhs[m];
          }
        }
      } // endfor m

      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase3D::GetFE3D(CurrentElement);
      nf = ele->GetNodalFunctional3D();

      FEDesc_Obj = ele->GetFEDesc3D();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();

      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Faces();
      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        InnerBoundary = false;
        OuterBoundary = false;

        if(joint->GetType() == BoundaryFace ||
           joint->GetType() == IsoBoundFace)
          OuterBoundary = true;

        /*
        // check whether neighbour does not belong to Coll
        neigh = joint->GetNeighbour(cell);
        if(neigh)
        {
          // check for neighbour's clipboard
          if(neigh->GetClipBoard() == -1)
            InnerBoundary = true;
        } 
        */

        if(InnerBoundary || OuterBoundary)
        {          
          cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
          cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
          t0 = 1.0/TmpLen[m];
          
          xf = 0; yf = 0; zf = 0;
          for (l1=0;l1<TmpLen[m];l1++)
          {
            cell->GetVertex(TmpFV[m*MaxLen+l1])
                        ->GetCoords(X[l1], Y[l1], Z[l1]);
            //LinComb[l1] = t0;
            xf += t0*X[l1];
            yf += t0*Y[l1];
            zf += t0*Z[l1];
          }

          /*if(OuterBoundary)
          {
            boundface = (TBoundFace *)joint;
            BoundComp = boundface->GetBoundComp();
          
            boundface->GetParameters(Param1, Param2);
            comp=BoundComp->GetID();

            BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                   X, Y, Z, Param1, Param2,
                                   xf, yf, zf, t0, t1);
				   }*/
	  // the face gets the b.c. which is valid at its center
          BoundaryCondition(xf, yf, zf, Cond0);

          switch(Cond0)
          {
            case DIRICHLET:
              if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
              {
                nf->GetPointsForFace(N_Points, t, s);
  
                for(l1=0;l1<N_Points;l1++)
                {
                  switch(TmpLen[m])
                  {
                    case 4:
                      LinComb[0] = (1-t[l1])*(1-s[l1]);
                      LinComb[1] = t[l1]*(1-s[l1]);
                      LinComb[2] = t[l1]*s[l1];
                      LinComb[3] = (1-t[l1])*s[l1];
                      xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                          +LinComb[2]*X[2] + LinComb[3]*X[3];
                      yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                          +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                      zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                          +LinComb[2]*Z[2] + LinComb[3]*Z[3];
                    break;
                    case 3:
                      LinComb[0] = 1-t[l1]-s[l1];
                      LinComb[1] = t[l1];
                      LinComb[2] = s[l1];
                      xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                          +LinComb[2]*X[2];
                      yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                          +LinComb[2]*Y[2];
                      zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                          +LinComb[2]*Z[2];
                    break;
                  }
                  /*if(OuterBoundary)
                    BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                           X, Y, Z, Param1, Param2,
                                           xf, yf, zf, t0, t1);
		  */
                  //cout << l1 << ": " << t0 << " " << t1 << " <> ";
                  //cout << xf << " " << yf << " " << zf << endl;
  
                  BoundaryValue(xf, yf, zf, PointValues[l1]);
                  // cout << " PV: " << PointValues[l1] << endl;
                }
              }
              else
              {
                // no isoparametric
                nf->GetPointsForFace(m, N_Points, xi, eta, zeta);
  
                TFEDatabase3D::GetOrigFromRef(reftrans, N_Points,
                                              xi, eta, zeta,
                                              X, Y, Z, AbsDetjk);
  
                for(l1=0;l1<N_Points;l1++)
                {
                  BoundaryValue(X[l1], Y[l1], Z[l1], PointValues[l1]);
                  // PointValues[l1] = 0;
                  // cout << " PV: " << PointValues[l1] << endl;
                }
              }
              
              nf->GetFaceFunctionals(Coll, cell, m, PointValues,
                                     FunctionalValues);
              EdgeDOF = FEDesc_Obj->GetJointDOF(m);
              N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              for(l=0;l<N_EdgeDOF;l++)
              {
// #pragma omp critical
                RHS[DOF[EdgeDOF[l]]] += FunctionalValues[l];
                // cout << " DOF: " << DOF[EdgeDOF[l]];
                // cout << " FV: " << FunctionalValues[l] << endl;
              }
            break;

            case NEUMANN:
            case ROBIN:  // lhs of robin will be assembled in main program
              // cout << "Neumann condition in Assemble3D" << endl;
              l = TFEDatabase3D::GetPolynomialDegreeFromFE3D
                (CurrentElement);
              switch(TmpLen[m])
              {
                case 3:
                  // triangular face
                  FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*l);
                  FaceQuadFormula = Gauss3Tria;
                break;
                case 4:
                  // quadrilateral face
		    FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*l);
                break;
              }
              qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
              qf2->GetFormulaData(N_Points, FaceWeights, t, s);
              // generate data on reference mesh cell for the 2d face of 3d cell
              TFEDatabase3D::GetBaseFunct3DFromFE3D(CurrentElement)
                ->MakeRefElementData(FaceQuadFormula);
              // values of base functions in all quadrature points on face 
              JointValues = TFEDatabase3D::GetJointValues3D
                (BaseFuncts[CurrentElement], FaceQuadFormula, m);
              TFEDatabase3D::GetBaseFunct3D(BaseFuncts[CurrentElement])
                      ->ChangeBF(Coll, cell, N_Points, JointValues);

              switch(TmpLen[m])
              {
                case 3:
                  xc1 = X[1] - X[0];
                  xc2 = X[2] - X[0];

                  yc1 = Y[1] - Y[0];
                  yc2 = Y[2] - Y[0];

                  zc1 = Z[1] - Z[0];
                  zc2 = Z[2] - Z[0];

                  // normal vector
                  nx = yc1*zc2 - zc1*yc2;
                  ny = zc1*xc2 - xc1*zc2;
                  nz = xc1*yc2 - yc1*xc2;
                  // determinant of reference trafo
                  t2 = sqrt(nx*nx + ny*ny + nz*nz);
                  
                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    t0 = t[l];
                    t1 = s[l];
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = 1-t0-t1;
                    LinComb[1] = t0;
                    LinComb[2] = t1;

                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                        +LinComb[2]*X[2];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                        +LinComb[2]*Y[2];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                        +LinComb[2]*Z[2];

                    /*if(OuterBoundary)
                      BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                             X, Y, Z, Param1, Param2,
                                             xf, yf, zf, t0, t1);
		    */
                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    // cout << t1 << endl;
                    t0 *= FaceWeights[l]*t2;
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound){
#ifdef _OPENMP
			#pragma omp critical
#endif
                        RHS[l3] += t0*JointValue[k];
		      }
                  } // endfor l
                break;
                
                case 4:
                  xc1=(-X[0] + X[1] + X[2] - X[3]) * 0.25;
                  xc2=(-X[0] - X[1] + X[2] + X[3]) * 0.25;
                  xc3=( X[0] - X[1] + X[2] - X[3]) * 0.25;

                  yc1=(-Y[0] + Y[1] + Y[2] - Y[3]) * 0.25;
                  yc2=(-Y[0] - Y[1] + Y[2] + Y[3]) * 0.25;
                  yc3=( Y[0] - Y[1] + Y[2] - Y[3]) * 0.25;

                  zc1=(-Z[0] + Z[1] + Z[2] - Z[3]) * 0.25;
                  zc2=(-Z[0] - Z[1] + Z[2] + Z[3]) * 0.25;
                  zc3=( Z[0] - Z[1] + Z[2] - Z[3]) * 0.25;

                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    t0 = 0.5*(t[l]+1);
                    t1 = 0.5*(s[l]+1);
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = (1-t0)*(1-t1);
                    LinComb[1] = t0*(1-t1);
                    LinComb[2] = t0*t1;
                    LinComb[3] = (1-t0)*t1;

                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                        +LinComb[2]*X[2] + LinComb[3]*X[3];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                        +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                        +LinComb[2]*Z[2] + LinComb[3]*Z[3];

                    /*if(OuterBoundary)
                      BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                             X, Y, Z, Param1, Param2,
                                             xf, yf, zf, t0, t1);
		    */
                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    nx = (yc1+s[l]*yc3)*(zc2+t[l]*zc3)
                        -(zc1+s[l]*zc3)*(yc2+t[l]*yc3);
                    ny = (zc1+s[l]*zc3)*(xc2+t[l]*xc3)
                        -(xc1+s[l]*xc3)*(zc2+t[l]*zc3);
                    nz = (xc1+s[l]*xc3)*(yc2+t[l]*yc3)
                        -(yc1+s[l]*yc3)*(xc2+t[l]*xc3);
                    t1 = nx*nx+ny*ny+nz*nz;
                    // cout << t1 << endl;
                    t0 *= FaceWeights[l]*sqrt(t1);
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound){
#ifdef _OPENMP
			#pragma omp critical
#endif
                        RHS[l3] += t0*JointValue[k];
		      }
                  } // endfor l
                break;
              }
              TFEDatabase3D::GetBaseFunct3D(BaseFuncts[CurrentElement])
                      ->ChangeBF(Coll, cell, N_Points, JointValues);
            break;
          } // endswitch Cond0
          
#ifdef _MPI         
        /** edges on this face are already set, so no need of checking in BD edges on this face */
        N_FaceEdges = ETmpLen[m];
        for(l1=0;l1<N_FaceEdges;l1++)
         {
// 	   #pragma omp critical 
	   {
          edge = cell->GetEdge(TmpFE[m*EMaxLen+l1]);
          edge->SetClipBoard(i);
	   }
         } // fo 
         
         /** vertices on this face are already set, so no need of checking in BD vert on this face, 30.06.12, sashi */        
         for (l1=0;l1<TmpLen[m];l1++)
	 { 
          cell->GetVertex(TmpFV[m*MaxLen+l1])->SetClipBoard(i);
	 }
#endif               
        } // endif
      } // endfor m
      
#ifdef _MPI   
    if(N_EdgeDOF = FEDesc_Obj->GetN_EdgeDOF() > 0) //conforming FE
     {
      cell->GetShapeDesc()->GetEdgeVertex(EdgeVertex);   
      N_Edges=cell->GetN_Edges();
      
      for(m=0;m<N_Edges;m++)
       {
        edge = cell->GetEdge(m);  
    
        if( (edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D) && edge->GetClipBoard()==-1)
         {      
          // BD edge is not yet set 
// 	   #pragma omp critical 
          edge->SetClipBoard(i);  
          EdgeDOF = FEDesc_Obj->GetEdgeDOF(m);   
          
          for(l=0;l<N_EdgeDOF;l++)
           {
            dof = DOF[EdgeDOF[l]];
            fespace->GetDOFPosition(dof, xf, yf, zf);
            BoundaryCondition(xf, yf, zf, Cond0);
 
            if(Cond0==DIRICHLET) // nothing to do for nin-Dirichlet
             {
             BoundaryValue(xf, yf, zf, t0);
             RHS[dof] = t0;  
             }       
          }  
          /** vertices on this edge are already set, so no need of checking in BD vert on this edge, 30.06.12, sashi */
          cell->GetVertex(EdgeVertex[2*m])->SetClipBoard(i);
          cell->GetVertex(EdgeVertex[2*m+1])->SetClipBoard(i);      
         }
       }//   endfor m
     } 
    /** correct the Vert Dirichlet   */
    if( FEDesc_Obj->GetN_VertDOF() > 0) //conforming FE
     {
      N_VertInCell = cell->GetN_Vertices(); 
     
      for(m=0;m<N_VertInCell;m++)
       {  
        Vert = cell->GetVertex(m);
       
        if(Vert->IsBoundVert() && Vert->GetClipBoard()==-1 )
         {
          // BD vert is not yet set 
//           #pragma omp critical 
	  Vert->SetClipBoard(i);
          Vert->GetCoords(xf, yf, zf);   
          BoundaryCondition(xf, yf, zf, Cond0);  
          
          if(Cond0==DIRICHLET) // nothing to do for nin-Dirichlet
           {
            dof = DOF[FEDesc_Obj->GetVertDOF(m)];   
            BoundaryValue(xf, yf, zf, t0);
            RHS[dof] = t0;  	    
           }
     
         }
       }//   endfor m    
      } //   if( FEDesc_Obj->GetN_VertDOF() > 0) 
#endif      
      
    } // endfor j
    // cout << "end i:" << i << endl;
  } // endfor i
#ifdef _OPENMP
  }
#endif

// printf("MatrixRow :: %lf\n",temp);
//   // ####################################################################
//   // modify matrix according to coupling
//   // ####################################################################
//   for(j=0;j<n_sqmatrices;j++)
//   {
//     fespace = sqmatrices[j]->GetFESpace();
//     N_ = fespace->GetN_Hanging();
//     // there are no hanging nodes
//     if (N_ == 0)
//       continue;
//     HangingNodes = fespace->GetHangingNodes();
// 
//     Entries = sqmatrices[j]->GetEntries();
//     RowPtr = sqmatrices[j]->GetRowPtr();
//     ColInd = sqmatrices[j]->GetKCol();
// 
//     CurrentHangingEntries = HangingEntries[j];
//     HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
//     HangingColInd = sqmatrices[j]->GetHangingKCol();
// 
//     ActiveBound = fespace->GetActiveBound();
// 
//     for(i=0;i<N_;i++)
//     {
//       hn = HangingNodes[i];
//       HNDescr = hn->GetType();
//       HNDescr_Obj = TFEDatabase3D::GetHNDesc3D(HNDescr);
//       k = HNDescr_Obj->GetN_Nodes();
//       Coupling = HNDescr_Obj->GetCoeff();
//       DOF = hn->GetDOF();
// 
//       end = HangingRowPtr[i+1];
//       for(n=HangingRowPtr[i];n<end;n++)
//       {
//         v = CurrentHangingEntries[n];
//         m = HangingColInd[n];
//         for(l=0;l<k;l++)
//         {
//           l1 = DOF[l]; 
//           if(l1<ActiveBound)
//           {
//             last=RowPtr[l1+1];
//             for(l2=RowPtr[l1];l2<last;l2++)
//             {
//               if(ColInd[l2] == m)
//               {
//                 Entries[l2] += Coupling[l] * v;
//               }
//             } // endfor l2
//           } // endif
//         } // endfor l
//       } // endfor n
//     } // endfor i
// 
//   } // endfor j
// 
//   for(j=0;j<n_rhs;j++)
//   {
//     fespace = ferhs[j];
//     N_Hanging = fespace->GetN_Hanging();
//     // there are no hanging nodes
//     if (N_Hanging == 0)
//       continue;
//     HangingNodes = fespace->GetHangingNodes();
// 
//     RHS = rhs[j];
//     CurrentHangingRhs = HangingRhs[j];
// 
//     ActiveBound = fespace->GetActiveBound();
// 
//     for(i=0;i<N_Hanging;i++)
//     {
//       hn = HangingNodes[i];
//       HNDescr = hn->GetType();
//       HNDescr_Obj = TFEDatabase3D::GetHNDesc3D(HNDescr);
//       N_ = HNDescr_Obj->GetN_Nodes();
//       Coupling = HNDescr_Obj->GetCoeff();
//       DOF = hn->GetDOF();
// 
//       for(k=0;k<N_;k++)
//       {
//         l = DOF[k]; 
//         if(l<ActiveBound)
//         {
//           RHS[l] += Coupling[k] * CurrentHangingRhs[i];
//         }
//       } // endfor k
//     } // endfor i
//   } // endfor j
// 
//   // ####################################################################
//   // write coupling into matrix
//   // ####################################################################
//   for(j=0;j<n_sqmatrices;j++)
//   {
//     fespace = sqmatrices[j]->GetFESpace();
//     N_ = fespace->GetN_Hanging();
//     // there are no hanging nodes
//     if (N_ == 0)
//       continue;
//     HangingNodes = fespace->GetHangingNodes();
// 
//     Entries = sqmatrices[j]->GetEntries();
//     RowPtr = sqmatrices[j]->GetRowPtr();
//     ColInd = sqmatrices[j]->GetKCol();
// 
//     ActiveBound = fespace->GetActiveBound();
// 
//     n = RowPtr[ActiveBound];
// 
//     for(i=0;i<N_;i++)
//     {
//       hn = HangingNodes[i];
//       HNDescr = hn->GetType();
//       HNDescr_Obj = TFEDatabase3D::GetHNDesc3D(HNDescr);
//       k = HNDescr_Obj->GetN_Nodes();
//       Coupling = HNDescr_Obj->GetCoeff();
//       DOF = hn->GetDOF();
// 
//       Entries[n] = 1.0;
//       n++;
// 
//       for(l=0;l<k;l++)
//       {
//         Entries[n] = - Coupling[l];
//         n++;
//       } // endfor l
// 
//     } // endfor i
//   } // endfor j

  if(n_sqmatrices)
  {
    delete [] GlobalNumbers;
    delete [] BeginIndex;

    for(i=0;i<n_sqmatrices;i++)
      delete [] HangingEntries[i];
    delete [] HangingEntries;
  }

  if(n_matrices)
  {
    delete [] AnsatzGlobalNumbers;
    delete [] AnsatzBeginIndex;
    delete [] TestGlobalNumbers;
    delete [] TestBeginIndex;
  }

  if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete [] HangingRhs[i];
    delete [] HangingRhs;

    delete [] righthand;
    delete [] LocRhs;
    delete [] RhsBeginIndex;
    delete [] RhsGlobalNumbers;
  }

  if(N_Parameters)
  {
    delete [] Param[0];
  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }

  delete [] AuxArray[0];

  time_total = GetTime() - time_total;
  //OutPut("LocalToGlobal: " << time_all << " " << time_total << " " << time_all/time_total << endl);
  PrintAllMat(n_sqmatrices, sqmatrices, n_matrices, matrices, n_rhs, rhs, ferhs);
} // end of Assemble

// =======================================================================
//
// Assemble3DSlipBC
//
// some manipulations in matrices and the rhs are necessary
//
// =======================================================================
void Assemble3DSlipBC(int n_fespaces, TFESpace3D **fespaces,
                      int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                      int n_matrices, TMatrix3D **matrices,
                      int n_rhs, double **rhs, TFESpace3D **ferhs,
                      TDiscreteForm3D *DiscreteForm3D,
                      BoundCondFunct3D **BoundaryConditions,
                      BoundValueFunct3D **BoundaryValues,
                      TAuxParam3D *Parameters)
{
  double hK;
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,k,l,l1,l2,l3,n,m, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints;
  int Used[N_FEs3D];
  int *N_BaseFunct;
  int N_Rows;
  BaseFunct3D *BaseFuncts;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  FE3D TestElement, AnsatzElement;
  QuadFormula2D FaceQuadFormula;
  TQuadFormula2D *qf2;
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  TBoundFace *boundface;
  TIsoBoundEdge *isoboundedge;
  //  TIsoBoundEdge *isoboundedge;
  int **GlobalNumbers, **BeginIndex;
  int **RhsGlobalNumbers, **RhsBeginIndex;
  int **TestGlobalNumbers, **TestBeginIndex;
  int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
  TFE3D *ele;
  TFEDesc3D *FEDesc_Obj;
  double *weights, *xi, *eta, *zeta;
  double *t, *s;
  double x, y, z;
  double xf, yf, zf;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *Param[MaxN_QuadPoints_3D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts3D];
  BaseFunct3D LocBF[N_BaseFuncts3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries;
  int *ColInd, *RowPtr;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp3D *BoundComp;
  double t0, t1;
  int comp;
  BoundCond Cond0, Cond1;
  BoundCondFunct3D *BoundaryCondition;
  BoundValueFunct3D *BoundaryValue;
  TNodalFunctional3D *nf;
  RefTrans3D reftrans;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_BaseFunctions3D];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *FaceWeights, *theta;
  double x0, x1, y0, y1, z0, z1, hE;
  double **JointValues, *JointValue;
  double Param1[4], Param2[4];
  double LinComb[4];
  bool *SecondDer;
 
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3;
  double nn, nx, ny, nz, t11, t12, t13, t21, t22, t23;
  double x10, y10, z10, x20, y20, z20;
  double *EdgePoints1, *EdgePoints2;
  double *Entries1,*Entries2,*Entries3, *Entries4, *Entries5;
  double *Entries6,*Entries7;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  int *ColInd4, *RowPtr4, *ColInd5, *RowPtr5, *ColInd6, *RowPtr6;
  int *ColInd7, *RowPtr7;
  double penetration_penalty, friction_parameter;
  double friction_constant= TDatabase::ParamDB->FRICTION_CONSTANT;
  double friction_power = TDatabase::ParamDB->FRICTION_POWER;
  double penetration_constant = TDatabase::ParamDB->PENETRATION_CONSTANT;
  double penetration_power = TDatabase::ParamDB->PENETRATION_POWER;
  double integral[3], val, det;
  int ii, jj, dof_ii, dof_jj, ll, found;

// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();


  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];
    for(i=0;i<n_sqmatrices;i++)
    {
      fespace = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();

      j = sqmatrices[i]->GetHangingN_Entries();
      HangingEntries[i] = new double [j];
      memset(HangingEntries[i], 0, SizeOfDouble*j);
    } // endfor
  } // endif n_sqmatrices
  
  if(n_matrices)
  {
    TestGlobalNumbers = new int* [n_matrices];
    AnsatzGlobalNumbers = new int* [n_matrices];
    TestBeginIndex = new int* [n_matrices];
    AnsatzBeginIndex = new int* [n_matrices];
    for(i=0;i<n_matrices;i++)
    {
      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetTestSpace();
      TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
      TestBeginIndex[i] = fespace->GetBeginIndex();
  
      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetAnsatzSpace();
      AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
      AnsatzBeginIndex[i] = fespace->GetBeginIndex();
    } // endfor
  } // endif n_matrices

  if(n_rhs)
  {
    HangingRhs = new double* [n_rhs];
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
  
      j = fespace->GetN_Hanging();
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    } // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions3D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions3D;

  } // endif n_rhs

  N_Parameters = Parameters->GetN_Parameters();
  if(N_Parameters)
  {
    aux = new double [MaxN_QuadPoints_3D*N_Parameters];
    for(j=0;j<MaxN_QuadPoints_3D;j++)
      Param[j] = aux + j*N_Parameters;
  }

  // 20 <= number of term in bilinear form
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  if(N_AllMatrices)
  {
    aux = new double
            [N_AllMatrices*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions3D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions3D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions3D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions3D;
  } // endif N_AllMatrices

 SecondDer = new bool[n_fespaces];
 SecondDer[0] = false;
// ########################################################################
// loop over all cells
// ########################################################################
  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  // cout << "Anzahl Zellen: " << N_Cells << endl;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    //cout << "Zellenr.: " << i << endl;
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE3D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################

    // SecondDer not set !!!
    reftrans=TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                                    Coll, cell, SecondDer,
                                    N_Points, xi, eta, zeta, weights,
                                    X, Y, Z, AbsDetjk);

    // ####################################################################
    // manipulate global matrices
    // manipulate global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      // get fe space
      fespace = ferhs[j];
      // get current finite element
      CurrentElement = fespace->GetFE3D(i, cell);
      // get number of basis function of current finite element
      N_ = N_BaseFunct[CurrentElement];

      // get information on rhs 
      local_rhs = righthand+j*MaxN_BaseFunctions2D;
      RHS = rhs[j];
      CurrentHangingRhs = HangingRhs[j];
    
      // find bounds in fe space
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      // dof of the rhs nodes connected to this cell
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

      // only for faces on the boundary
      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase3D::GetFE3D(CurrentElement);
      nf = ele->GetNodalFunctional3D();
      nf->GetPointsForFace(N_EdgePoints, EdgePoints1, EdgePoints2);

      FEDesc_Obj = ele->GetFEDesc3D();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();

      // find slip type bc
      N_Joints = cell->GetN_Faces();
      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        if (joint->GetType() == BoundaryFace)
           //|| (joint->GetType() == IsoBoundaryFace))
        {
	    // !!! everything is commented which is not necessary for computing
	    // !!! (xf, yf, zf)
	    /*if(joint->GetType() == BoundaryFace)
          {
            boundface = (TBoundFace *)joint;
            BoundComp = boundface->GetBoundComp();
            boundface->GetParameters(Param1, Param2);
	    }*/
          /* else
          {
            isoboundface = (TIsoBoundFace *)joint;
            BoundComp = isoboundface->GetBoundComp();
            isoboundface->GetParameters(Param1, Param2);
            }*/
          // get id of the boundary component
          //comp=BoundComp->GetID();

          // compute point on the boundary face
          cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

          t0 = 1.0/TmpLen[m];
	  xf = 0; yf = 0; zf = 0;
          for (l1=0;l1<TmpLen[m];l1++)
          {
            cell->GetVertex(TmpFV[m*MaxLen+l1])
                        ->GetCoords(X[l1], Y[l1], Z[l1]);
            //LinComb[l1] = t0;
            xf += t0*X[l1];
            yf += t0*Y[l1];
            zf += t0*Z[l1];
          }
	  /*
	  OutPut(xf  << " " << yf << " " << zf << " ");

          t0 = 1.0/TmpLen[m];
          
          for (l1=0;l1<TmpLen[m];l1++)
          {
            // get coordinates of the vertices belonging to the 
            // boundary face
            cell->GetVertex(TmpFV[m*MaxLen+l1])
                        ->GetCoords(X[l1], Y[l1], Z[l1]);
            LinComb[l1] = t0;
          }
          // get coordinates (xf, yf, zf) and parameters (t0, t1) of the 
          // barycentre of  the boundary face
          BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                 X, Y, Z, Param1, Param2,
                                 xf, yf, zf, t0, t1);
	  OutPut(xf  << " " << yf << " " << zf << " "<<endl);
	  */
          // get boundary condition
          // it is implicitly assumed that only one bc is 
          // assigned to the face
          BoundaryCondition(xf, yf, zf, Cond0);

          switch(Cond0)
          {
            case DIRICHLET:
              break;
              
            case NEUMANN:
              break;
              
            case SLIP:
              exit(4711);
              break;
              
            case SLIP_FRICTION_PENETRATION_RESISTANCE:
              // get polynomial degree of finite element in current cell
              l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(CurrentElement);
              // definequadrature formulas
              switch(TmpLen[m])
              {
                case 3:
                  // triangular face
                  FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*l);
                  //FaceQuadFormula = Gauss3Tria;
                break;
                case 4:
                  // quadrilateral face
                  FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*l);
                break;
              }
              // get quadrature formulas for 2D face
              qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
              // get quad points and weights
              qf2->GetFormulaData(N_Points, FaceWeights, t, s);
              TFEDatabase3D::GetBaseFunct3DFromFE3D(CurrentElement)
                ->MakeRefElementData(FaceQuadFormula);
              // get values of test functions in all quadrature points
              // on joint m 
              JointValues = TFEDatabase3D::GetJointValues3D
                (BaseFuncts[CurrentElement], FaceQuadFormula, m);

              // compute unit normal vector of the face
              // the face is assumed to be a plane spanned by the first 
              // three vertices
              x10 = X[1] - X[0];
              y10 = Y[1] - Y[0];
              z10 = Z[1] - Z[0];
              x20 = X[2] - X[0];
              y20 = Y[2] - Y[0];
              z20 = Z[2] - Z[0];
              nx = y10*z20-z10*y20;
              ny = z10*x20-z20*x10;
              nz = x10*y20-x20*y10;
              hE= sqrt(nx*nx+ny*ny+nz*nz);
              nx /= hE;
              ny /= hE;
              nz /= hE;
              //OutPut(X[0] << " "<< Y[0] << " " << Z[0] << "    "); 
              //OutPut(X[1] << " "<< Y[1] << " " << Z[1] << "    "); 
              //OutPut(X[2] << " "<< Y[2] << " " << Z[2] << endl); 
              // compute two tangential vectors of the face
              if ( (fabs(nx)>=0.5) || (fabs(ny)>=0.5))
              {
                nn = sqrt(nx*nx+ny*ny);
                t11 = ny/nn;
                t12 = -nx/nn;
                t13 = 0;
                t21 = -t12*nz;
                t22 = t11*nz;
                t23 = t12*nx-t11*ny;
              }
              else
              {
                nn = sqrt(ny*ny+nz*nz);
                t11 = 0;
                t12 = -nz/nn;
                t13 = ny/nn;
                t21 = t13*ny-t12*nz;
                t22 = - t13*nx;
                t23 = t12*nx;
              }
              //OutPut("nx " << nx << " ny " << ny << " nz " << nz);
              //OutPut(" t11 " << t11 << " t12 " << t12 << " t13 " << t13);
              //OutPut(" t21 " << t21 << " t22 " << t22 << " t23 " << t23 << endl);
 
              switch(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY)
              {
                case 0:
                  // standard case, same value on all walls
                  // penalty value for weak imposition of no penetration bc
                  penetration_penalty = penetration_constant*pow(hK,penetration_power);              
                  // parameter for friction
                  //BoundaryValue(comp, (t0+t1)/2.0,friction_parameter);
                  friction_parameter = friction_constant * pow(hK,friction_power);
                  //OutPut(penetration_penalty << " " << friction_parameter << endl);
                  break;
                case 1:
                  // ChannelStepSlip3D.h
                  penetration_penalty = penetration_constant*pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * pow(hK,friction_power);
                  // free slip else
                  if ((fabs(Z[0])<1e-6) && (fabs(Z[1])<1e-6) && (fabs(Z[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((fabs(Z[0]-10)<1e-6) && (fabs(Z[1]-10)<1e-6) && (fabs(Z[2]-10)<1e-6))
                    friction_parameter = 0.0;
                  break;
                case 2:
                  // ChannelStepSlipFreeSlip.h
                  penetration_penalty = penetration_constant*pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * pow(hK,friction_power);
                  // free slip else
                  if ((fabs(Z[0])<1e-6) && (fabs(Z[1])<1e-6) && (fabs(Z[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((fabs(Z[0]-10)<1e-6) && (fabs(Z[1]-10)<1e-6) && (fabs(Z[2]-10)<1e-6))
                    friction_parameter = 0.0;
                  if ((fabs(Y[0]-10)<1e-6) && (fabs(Y[1]-10)<1e-6) && (fabs(Y[2]-10)<1e-6))
                    friction_parameter = 0.0;
                  break;
                case 3:
                  // CircularCylinder22000.h
                  penetration_penalty = penetration_constant*pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * pow(hK,friction_power);
                  //OutPut(penetration_penalty << " : " << friction_parameter << endl);
                  // free slip left and right
                  if ((fabs(Y[0])<1e-6) && (fabs(Y[1])<1e-6) && (fabs(Y[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((fabs(Y[0]-1.4)<1e-6) && (fabs(Y[1]-1.4)<1e-6) && (fabs(Y[2]-1.4)<1e-6))
                    friction_parameter = 0.0;
                  // free slip top and bottom
                  if ((fabs(Z[0])<1e-6) && (fabs(Z[1])<1e-6) && (fabs(Z[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((fabs(Z[0]-0.4)<1e-6) && (fabs(Z[1]-0.4)<1e-6) && (fabs(Z[2]-0.4)<1e-6))
                    friction_parameter = 0.0;
		  //if ( friction_parameter != 0.0)
		  // OutPut(penetration_penalty << " : " << friction_parameter << 
		  //	   " " << X[0] << " " << Y[0] << " " << Z[0] << endl);
                  break;
                case 4:
                  // WallMountedCube.h
                  penetration_penalty = penetration_constant*pow(hK,penetration_power);              
                  // friction for y = 0
                  friction_parameter = friction_constant * pow(hK,friction_power);
                  //OutPut(penetration_penalty << " : " << friction_parameter << endl);
                  // free slip left and right
                  if ((fabs(Y[0])<1e-6) && (fabs(Y[1])<1e-6) && (fabs(Y[2])<1e-6))
                    friction_parameter = 0.0;
                  if ((fabs(Y[0]-0.7)<1e-6) && (fabs(Y[1]-0.7)<1e-6) && (fabs(Y[2]-0.7)<1e-6))
                    friction_parameter = 0.0;
		  //if ( friction_parameter != 0.0)
		  //OutPut(penetration_penalty << " : " << friction_parameter << 
		  //  " " << X[0] << " " << Y[0] << " " << Z[0] << endl);
                  break;
                case 5:
                  // windtunnel_fine.h
                  // upper wall, symmetrie wall -> free slip, no penetration
                  if ((fabs(Z[0])<1e-6) && (fabs(Z[1])<1e-6) && (fabs(Z[2])<1e-6))
                  {
                    penetration_penalty = penetration_constant*pow(hK,penetration_power);
                  }
                  else
                    penetration_penalty = 0.0;
                  friction_parameter = 0.0;
                  break;
                case 6:
                  // windtunnel_m2.h
                  // free slip, no penetration
                  penetration_penalty = 1e12;
                  friction_parameter = 0.0;
                  break;
                default:
                  OutPut("INTERNAL_SLIP_WITH_FRICTION_IDENTITY not implemented !!!"<< endl);
                  exit(4711);
                   
              }  
              
              EdgeDOF = FEDesc_Obj->GetJointDOF(m);
              
              // compute additional matrix entries
              // for all velo dof in the mesh cell
              // ii - test function
              for (ii=0;ii<N_;ii++)
              {
                // look for 'ii'-th row in all matrices
                dof_ii = DOF[ii];
                // Dirichlet node
                if (dof_ii>=ActiveBound)
                  continue;
                
                // !!!!!! assumed that A_11 - A_33 are in sqmatrices[0] - [8]
                // ordered as in main program 
                // first velocity component -> matrices A_11, A_12, A13 (and M_11)
                if (j==0)
                {
                  // A11
                  Entries1 = sqmatrices[0]->GetEntries();
                  RowPtr1 = sqmatrices[0]->GetRowPtr();
                  ColInd1 = sqmatrices[0]->GetKCol();
                  
                  if (n_sqmatrices>3)
                  {
                    // A12
                    Entries2 = sqmatrices[3]->GetEntries();
                    RowPtr2 = sqmatrices[3]->GetRowPtr();
                    ColInd2 = sqmatrices[3]->GetKCol();
                     // A13
                    Entries3 = sqmatrices[4]->GetEntries();
                    RowPtr3 = sqmatrices[4]->GetRowPtr();
                    ColInd3 = sqmatrices[4]->GetKCol();
                  }
                  
                  // time dependent problem and NSTYPE 4
                  if (n_sqmatrices==18)
                  {
                    // M11
                    Entries4 = sqmatrices[9]->GetEntries();
                    RowPtr4 = sqmatrices[9]->GetRowPtr();
                    ColInd4 = sqmatrices[9]->GetKCol();
                    // M12
                    Entries5 = sqmatrices[12]->GetEntries();
                    RowPtr5 = sqmatrices[12]->GetRowPtr();
                    ColInd5 = sqmatrices[12]->GetKCol();                      
                    // M13
                    Entries6 = sqmatrices[13]->GetEntries();
                    RowPtr6 = sqmatrices[13]->GetRowPtr();
                    ColInd6 = sqmatrices[13]->GetKCol();                      
                  }
                  
                  if (n_matrices==3)
                  {
                    // B1T
                    Entries7 = matrices[0]->GetEntries();
                    RowPtr7 = matrices[0]->GetRowPtr();
                    ColInd7 = matrices[0]->GetKCol();                      
                  }
                }
                // second velocity component -> matrices A_21, A_22, A23
                if (j==1)
                {
                  if (n_sqmatrices>3)
                  {
                    // A21 
                    Entries1 = sqmatrices[5]->GetEntries();
                    RowPtr1 = sqmatrices[5]->GetRowPtr();
                    ColInd1 = sqmatrices[5]->GetKCol();
                    // A23 
                    Entries3 = sqmatrices[6]->GetEntries();
                    RowPtr3 = sqmatrices[6]->GetRowPtr();
                    ColInd3 = sqmatrices[6]->GetKCol();
                  }
                  // A22
                  Entries2 = sqmatrices[1]->GetEntries();
                  RowPtr2 = sqmatrices[1]->GetRowPtr();
                  ColInd2 = sqmatrices[1]->GetKCol();
                  
                  // time dependent problem and NSTYPE 4
                  if (n_sqmatrices==18)
                  {
                    // M22
                    Entries4 = sqmatrices[10]->GetEntries();
                    RowPtr4 = sqmatrices[10]->GetRowPtr();
                    ColInd4 = sqmatrices[10]->GetKCol();   
                    // M21
                    Entries5 = sqmatrices[14]->GetEntries();
                    RowPtr5 = sqmatrices[14]->GetRowPtr();
                    ColInd5 = sqmatrices[14]->GetKCol(); 
                    // M23
                    Entries6 = sqmatrices[15]->GetEntries();
                    RowPtr6 = sqmatrices[15]->GetRowPtr();
                    ColInd6 = sqmatrices[15]->GetKCol();                      
                  }

                  if (n_matrices==3)
                  {
                    // B2T
                    Entries7 = matrices[1]->GetEntries();
                    RowPtr7 = matrices[1]->GetRowPtr();
                    ColInd7 = matrices[1]->GetKCol();                      
                  }
                }
                // third velocity component -> matrices A_31, A_32, A33
                if (j==2)
                {
                  if (n_sqmatrices>3)
                  {
                    // A31 
                    Entries1 = sqmatrices[7]->GetEntries();
                    RowPtr1 = sqmatrices[7]->GetRowPtr();
                    ColInd1 = sqmatrices[7]->GetKCol();
                    // A32 
                    Entries2 = sqmatrices[8]->GetEntries();
                    RowPtr2 = sqmatrices[8]->GetRowPtr();
                    ColInd2 = sqmatrices[8]->GetKCol();
                  }
                  // A33
                  Entries3 = sqmatrices[2]->GetEntries();
                  RowPtr3 = sqmatrices[2]->GetRowPtr();
                  ColInd3 = sqmatrices[2]->GetKCol();
                  
                  // time dependent problem and NSTYPE 4
                  if (n_sqmatrices==18)
                  {
                    // M33
                    Entries4 = sqmatrices[11]->GetEntries();
                    RowPtr4 = sqmatrices[11]->GetRowPtr();
                    ColInd4 = sqmatrices[11]->GetKCol();   
                    // M31
                    Entries5 = sqmatrices[16]->GetEntries();
                    RowPtr5 = sqmatrices[16]->GetRowPtr();
                    ColInd5 = sqmatrices[16]->GetKCol(); 
                    // M32
                    Entries6 = sqmatrices[17]->GetEntries();
                    RowPtr6 = sqmatrices[17]->GetRowPtr();
                    ColInd6 = sqmatrices[17]->GetKCol();                      
                  }

                  if (n_matrices==3)
                  {
                    // B3T
                    Entries7 = matrices[2]->GetEntries();
                    RowPtr7 = matrices[2]->GetRowPtr();
                    ColInd7 = matrices[2]->GetKCol();                      
                  }
                }
                //OutPut("ii " << dof_ii << endl);
                
                // for all dof in the mesh cell
                // jj - ansatz function
                for (jj=0;jj<N_;jj++)
                {
                  dof_jj = DOF[jj];
                  // OutPut("jj " << dof_jj << endl);
                  // initialize the boundary integrals
                  for (l=0;l<3;l++)
                    integral[l] = 0;
                  
                  // compute boundary integrals
                  // first component of the velocity
                  if (j==0)
                  {
                    // for all quadrature points
                    for(l=0;l<N_Points;l++)
                    {
                      // values of test functions in this quadrature point
                      JointValue = JointValues[l];                      
                      //  weight times determinant of reference trafo
                      det = FaceWeights[l]*hE;
                      // (A_11)_{ii,jj}
                      val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*nx;
                      val += friction_parameter*JointValue[jj]*t11*JointValue[ii]*t11;
                      val += friction_parameter*JointValue[jj]*t21*JointValue[ii]*t21;
                      integral[0] += val*det;
                      // (A_12)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*nx;
                      val+= friction_parameter*JointValue[jj]*t11*JointValue[ii]*t12;
                      val+= friction_parameter*JointValue[jj]*t21*JointValue[ii]*t22;
                      integral[1] += val*det;
                      // (A_13)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*nz*JointValue[ii]*nx;
                      val+= friction_parameter*JointValue[jj]*t11*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t21*JointValue[ii]*t23;
                      integral[2] += val*det;
                    } // endfor l
                    
                    // update first matrix
                    found = 0;
                    for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                    {
                      if (ColInd1[ll] == dof_jj)
                      {
			  //OutPut("a11 " << integral[0] << endl);
                        Entries1[ll] += integral[0];
                        //OutPut("a11 " << integral[0] << " " << Entries1[ll] << endl);
                        found = 1;
                        break;
                      }
                    }
                    if (!found)
                    {
                      OutPut("ERROR A_11 " << endl);
                      exit(4711);
                    }

                    if (n_sqmatrices>3)
                    {
                      // update second matrix
                      found = 0;
                      for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                      {
                        if (ColInd2[ll] == dof_jj)
                        {
                          found = 1;
                          //OutPut("a12 " << integral[1] << endl);
                          Entries2[ll] += integral[1];
                          //OutPut("a12 " << integral[1] << " " << Entries2[ll] << endl);
                          break;
                        }
                      }
                      if (!found)
                      {
                        OutPut("ERROR A_12 "<< endl);
                        exit(4711);
                      }
                      // update third matrix
                      found = 0;
                      for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                      {
                        if (ColInd3[ll] == dof_jj)
                        {
                          found = 1;
                          //OutPut("a13 " << integral[2] << endl);
                          Entries3[ll] += integral[2];
                          //OutPut("a13 " << integral[2] << " " << Entries3[ll] << endl);
                          break;
                        }
                      }
                      if (!found)
                      {
                        OutPut("ERROR A_13 "<< endl);
                        exit(4711);
                      }
                    }
                    
                    /*if (n_sqmatrices==18)
                    { // M_11, set off diagonal to zero 
                      for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                      {
                        if (ColInd4[ll] != dof_ii)        
                          Entries4[ll] = 0;
                      }
                      // M_12, set row to zero 
                      for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                        Entries5[ll] = 0;
                      // M_13, set row to zero 
                      for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                        Entries6[ll] = 0;
                    }
                    
                    if (n_matrices==3)
                      for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                        Entries7[ll] = 0;
                    
                    // set rhs to zero
                    RHS[dof_ii] = 0;  */    
                  } // end first component (j==0)
                  
                  // compute boundary integrals
                  // second component of the velocity
                  if (j==1)
                  {
                    // for all quadrature points
                    for(l=0;l<N_Points;l++)
                    {
                      // values of test functions in this quadrature point
                      JointValue = JointValues[l];                      
                      //  weight times determinant of reference trafo
                      det = FaceWeights[l]*hE;
                      // (A_21)_{ii,jj}
                      val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*ny;
                      val += friction_parameter*JointValue[jj]*t11*JointValue[ii]*t12;
                      val += friction_parameter*JointValue[jj]*t21*JointValue[ii]*t22;
                      integral[0] += val*det;
                      // (A_22)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*ny;
                      val+= friction_parameter*JointValue[jj]*t12*JointValue[ii]*t12;
                      val+= friction_parameter*JointValue[jj]*t22*JointValue[ii]*t22;
                      integral[1] += val*det;
                      // (A_23)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*nz*JointValue[ii]*ny;
                      val+= friction_parameter*JointValue[jj]*t12*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t22*JointValue[ii]*t23;
                      integral[2] += val*det;
                    } // endfor l
                    
                    // update first matrix
                    found = 0;
                    for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                    {
                      if (ColInd2[ll] == dof_jj)
                      {
			  //OutPut("a22 " << integral[1] << endl);
                        Entries2[ll] += integral[1];
                        // if (fabs(ny)!>1e-6)
			//OutPut("a22 " << integral[1] << " " << Entries2[ll] << endl);
                        found = 1;
                        break;
                      }
                    }
                    if (!found)
                    {
                      OutPut("ERROR A_22 " << endl);
                      exit(4711);
                    }

                    if (n_sqmatrices>3)
                    {
                      // update second matrix
                      found = 0;
                      for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                      {
                        if (ColInd1[ll] == dof_jj)
                        {
                          found = 1;
                          //OutPut("a21 " << integral[0] << endl);
                          Entries1[ll] += integral[0];
                          //OutPut("a21 " << Entries1[ll] << endl);
                          break;
                        }
                      }
                      if (!found)
                      {
                        OutPut("ERROR A_21 "<< endl);
                        exit(4711);
                      }
                      // update third matrix
                      found = 0;
                      for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                      {
                        if (ColInd3[ll] == dof_jj)
                        {
                          found = 1;
                          //OutPut("a23 " << integral[2] << endl);
                          Entries3[ll] += integral[2];
                          //OutPut("a23 " <<  Entries3[ll] << endl);
                          break;
                        }
                      }
                      if (!found)
                      {
                        OutPut("ERROR A_23 "<< endl);
                        exit(4711);
                      }
                    }
                    
                    /* if (n_sqmatrices==18)
                    { // M_11, set off diagonal to zero 
                      for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                      {
                        if (ColInd4[ll] != dof_ii)        
                          Entries4[ll] = 0;
                      }
                      // M_12, set row to zero 
                      for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                        Entries5[ll] = 0;
                      // M_13, set row to zero 
                      for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                        Entries6[ll] = 0;
                    }
                    
                    if (n_matrices==3)
                      for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                        Entries7[ll] = 0;
                    
                    // set rhs to zero
                    RHS[dof_ii] = 0;   */   
                  } // end second component (j==1)

                  // compute boundary integrals
                  // third component of the velocity
                  if (j==2)
                  {
                    // for all quadrature points
                    for(l=0;l<N_Points;l++)
                    {
                      // values of test functions in this quadrature point
                      JointValue = JointValues[l];                      
                      //  weight times determinant of reference trafo
                      det = FaceWeights[l]*hE;
                      // (A_31)_{ii,jj}
                      val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*nz;
                      val += friction_parameter*JointValue[jj]*t11*JointValue[ii]*t13;
                      val += friction_parameter*JointValue[jj]*t21*JointValue[ii]*t23;
                      integral[0] += val*det;
                      // (A_23)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*nz;
                      val+= friction_parameter*JointValue[jj]*t12*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t22*JointValue[ii]*t23;
                      integral[1] += val*det;
                      // (A_33)_{ii,jj}
                      val =  penetration_penalty*JointValue[jj]*nz*JointValue[ii]*nz;
                      val+= friction_parameter*JointValue[jj]*t13*JointValue[ii]*t13;
                      val+= friction_parameter*JointValue[jj]*t23*JointValue[ii]*t23;
                      integral[2] += val*det;
                    } // endfor l
                    
                    // update first matrix
                    found = 0;
                    for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                    {
                      if (ColInd3[ll] == dof_jj)
                      {
                        Entries3[ll] += integral[2];
                        //OutPut("a33 " << Entries3[ll] << endl);
                        found = 1;
                        break;
                      }
                    }
                    if (!found)
                    {
                      OutPut("ERROR A_33 " << endl);
                      exit(4711);
                    }

                    if (n_sqmatrices>3)
                    {
                      // update second matrix
                      found = 0;
                      for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                      {
                        if (ColInd1[ll] == dof_jj)
                        {
                          found = 1;
                          //OutPut("a31 " << integral[0] << endl);
                          Entries1[ll] += integral[0];
                          // OutPut("a31 " << Entries1[ll] << endl);
                          break;
                        }
                      }
                      if (!found)
                      {
                        OutPut("ERROR A_31 "<< endl);
                        exit(4711);
                      }
                      // update third matrix
                      found = 0;
                      for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                      {
                        if (ColInd2[ll] == dof_jj)
                        {
                          found = 1;
                          //OutPut("a32 " << integral[1] << endl);
                          Entries2[ll] += integral[1];
                          //OutPut("a32 " << Entries2[ll] << endl);
                          break;
                        }
                      }
                      if (!found)
                      {
                        OutPut("ERROR A_32 "<< endl);
                        exit(4711);
                      }
                    }
                    
                    /*if (n_sqmatrices==18)
                    { // M_11, set off diagonal to zero 
                      for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                      {
                        if (ColInd4[ll] != dof_ii)        
                          Entries4[ll] = 0;
                      }
                      // M_12, set row to zero 
                      for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                        Entries5[ll] = 0;
                      // M_13, set row to zero 
                      for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                        Entries6[ll] = 0;
                    }
                    
                    if (n_matrices==3)
                      for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                        Entries7[ll] = 0;
                    
                    // set rhs to zero
                    RHS[dof_ii] = 0; */    
                  } // end third component (j==2)
                  
                }// endfor ansatz functions (jj)
              }// enfor test functions (ii)
         } // end switch (Cond0)
        } // endif (boundary joint)
      }  // endfor m (N_Joints)
    } // endfor j (n_rhs)
  }  // endfor i (N_Cells)

  if(n_sqmatrices)
  {
    delete [] GlobalNumbers;
    delete [] BeginIndex;

    for(i=0;i<n_sqmatrices;i++)
      delete [] HangingEntries[i];
    delete [] HangingEntries;
  }

  if(n_matrices)
  {
    delete [] AnsatzGlobalNumbers;
    delete [] AnsatzBeginIndex;
    delete [] TestGlobalNumbers;
    delete [] TestBeginIndex;
  }

  if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete [] HangingRhs[i];
    delete [] HangingRhs;

    delete [] righthand;
    delete [] LocRhs;
    delete [] RhsBeginIndex;
    delete [] RhsGlobalNumbers;
  }

  if(N_Parameters)
  {
    delete [] Param[0];
  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }

  delete [] AuxArray[0];
  delete [] SecondDer;

/*
  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
  for(k=0;k<n_sqmatrices;k++)
  {
    cout << endl;
    cout << "sqmatrix: " << k << endl;
    RowPtr = sqmatrices[k]->GetRowPtr();
    Entries = sqmatrices[k]->GetEntries();
    ColInd = sqmatrices[k]->GetKCol();
    N_Rows = sqmatrices[k]->GetN_Rows();
    for(i=0;i<N_Rows;i++)
    {
      end=RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
      {
        // cout << j << endl;
        cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
        cout << setw(10) << Entries[j] << endl;
      }
    }
    cout << endl;
  } // endfor k
  
  for(k=0;k<n_matrices;k++)
  {
    cout << endl;
    cout << "matrix: " << k << endl;
    RowPtr = matrices[k]->GetRowPtr();
    Entries = matrices[k]->GetEntries();
    ColInd = matrices[k]->GetKCol();
    N_Rows = matrices[k]->GetN_Rows();
    for(i=0;i<N_Rows;i++)
    {
      end=RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
      {
        // cout << j << endl;
        cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
        cout << setw(10) << Entries[j] << endl;
      }
    }
    cout << endl;
  } // endfor k

  for(k=0;k<n_rhs;k++)
  {
    cout << "rhs: " << k << endl;
    N_Rows = ferhs[k]->GetN_DegreesOfFreedom();
    RHS=rhs[k];
    for(i=0;i<N_Rows;i++)
      cout << setw(5) << i << setw(20) << RHS[i] << endl;
  }
*/

} // end of Assemble

/*************************************+******************************/
/*
/* Modification of Matrices for slip with friction bc for 
/* better condition
/*
/*******************************************************************/
void ModifyMatrixSlipBC(TSquareMatrix3D **sqmatrices, TMatrix3D **matrices,
    int N_U, double *rhs)
{
    int i, j, k, j0, j1, index, *ColInd, *RowPtr, row_off, col_off;
    double *setzero, bound, *Entries, maximal = -1;; 

    return;

    bound = TDatabase::ParamDB->PENETRATION_CONSTANT*1e-4;

    if (fabs(TDatabase::ParamDB->PENETRATION_POWER+2) >1e6)
	OutPut("ModifyMatrixSlipBC : recommended to set PENETRATION_POWER to -2 !!!" << endl);

    setzero = new double[3*N_U];
    memset(setzero, 0, SizeOfDouble*3*N_U);
    
    // find largest values of the matrices A_1k
    for (k=0;k<3;k++)
    {
	ColInd = sqmatrices[3*k]->GetKCol();
	RowPtr = sqmatrices[3*k]->GetRowPtr();
	Entries = sqmatrices[3*k]->GetEntries();
	for (i=0;i<N_U;i++)
	{
	    // i-th row of sqmatrix
	    j0 = RowPtr[i];
	    j1 = RowPtr[i+1];
	    for(j=j0;j<j1;j++)
	    {
		// column
		index = ColInd[j];
		if (fabs(Entries[j])>bound)
		    setzero[index] = 1;
		if (fabs(Entries[j])>maximal)
		    maximal = fabs(Entries[j]);		
	    }
	}
    }

    // find largest values of the matrices A_2k
    for (k=0;k<3;k++)
    {
	ColInd = sqmatrices[3*k+1]->GetKCol();
	RowPtr = sqmatrices[3*k+1]->GetRowPtr();
	Entries = sqmatrices[3*k+1]->GetEntries();
	for (i=0;i<N_U;i++)
	{
	    // i-th row of sqmatrix
	    j0 = RowPtr[i];
	    j1 = RowPtr[i+1];
	    for(j=j0;j<j1;j++)
	    {
		// column
		index = ColInd[j];
		if (fabs(Entries[j])>bound)
		    setzero[N_U+index] = 1;
		if (fabs(Entries[j])>maximal)
		    maximal = fabs(Entries[j]);		
	    }
	}
    }
    // find largest values of the matrices A_3k
    for (k=0;k<3;k++)
    {
	ColInd = sqmatrices[3*k+2]->GetKCol();
	RowPtr = sqmatrices[3*k+2]->GetRowPtr();
	Entries = sqmatrices[3*k+2]->GetEntries();
	for (i=0;i<N_U;i++)
	{
	    // i-th row of sqmatrix
	    j0 = RowPtr[i];
	    j1 = RowPtr[i+1];
	    for(j=j0;j<j1;j++)
	    {
		// column
		index = ColInd[j];
		if (fabs(Entries[j])>bound)
		    setzero[2*N_U+index] = 1;
		if (fabs(Entries[j])>maximal)
		    maximal = fabs(Entries[j]);		
	    }
	}
    }
    OutPut("maximal matrix entry " << maximal << endl);

    // set homogeneous Dirichlet bc for the marked nodes
    // loop over the matrix blocks A_ij
    for (k=0;k<9;k++)
    {
	ColInd = sqmatrices[k]->GetKCol();
	RowPtr = sqmatrices[k]->GetRowPtr();
	Entries = sqmatrices[k]->GetEntries();
	switch(k)
	{
	    case 0:
		row_off = 0;
		col_off = 0;
		break;
	    case 1:
		row_off = 0;
		col_off = N_U;
		break;
	    case 2:
		row_off = 0;
		col_off = 2*N_U;
		break;
	    case 3:
		row_off = N_U;
		col_off = 0;
		break;
	    case 4:
		row_off = N_U;
		col_off = N_U;
		break;
	    case 5:
		row_off = N_U;
		col_off = 2*N_U;
		break;
	    case 6:
		row_off = 2*N_U;
		col_off = 0;
		break;
	    case 7:
		row_off = 2*N_U;
		col_off = N_U;
		break;
	    case 8:
		row_off = 2*N_U;
		col_off = 2*N_U;
		break;
	}
	// loop over the dof
	for (i=0;i<N_U;i++)
	{
	    // i-th row of sqmatrix
	    j0 = RowPtr[i];
	    j1 = RowPtr[i+1];
	    // find diagonal
	    for(j=j0;j<j1;j++)
	    {
		index = ColInd[j];
		if (setzero[i+row_off])
		{
		    if ((index==i)&&(row_off==col_off))
			Entries[j] = 1.0;
		    else 
			Entries[j] = 0.0;
		    continue;
		}
		if (setzero[index+col_off])
		    Entries[j] = 0.0;
	    }
	}
    }
   
    // loop over the matrix blocks B_k
    for (k=0;k<3;k++)
    {
	ColInd = matrices[k]->GetKCol();
	RowPtr = matrices[k]->GetRowPtr();
	Entries = matrices[k]->GetEntries();
	switch(k)
	{
	    case 0:
		row_off = 0;
		break;
	    case 1:
		row_off = N_U;
		break;
	    case 2:
		row_off = 2*N_U;
		break;
	}
	for (i=0;i<N_U;i++)
	{
	    // i-th row of sqmatrix
	    j0 = RowPtr[i];
	    j1 = RowPtr[i+1];
	    for(j=j0;j<j1;j++)
	    {
		if (setzero[i+row_off])
		{
		    Entries[j] = 0.0;
		}
	    }
	}
    }
    // rhs
    for (i=0;i<3*N_U;i++)
    {
	if (setzero[i])
	    rhs[i] = 0;
    }

    delete  [] setzero;
}

/*
  Assemble3D_mixed:
    Assemble for vector finite elements (Raviart-Thomas)
    Need the global orientation of normal at each inner edge/face

    implementation: Alfonso (07.09.2010)
*/
void Assemble3D_mixed(int n_fespaces, TFESpace3D **fespaces,
int n_sqmatrices, TSquareMatrix3D **sqmatrices,
int n_matrices, TMatrix3D **matrices,
int n_rhs, double **rhs, TFESpace3D **ferhs,
TDiscreteForm3D *DiscreteForm3D,
BoundCondFunct3D **BoundaryConditions,
BoundValueFunct3D **BoundaryValues,
TAuxParam3D *Parameters)
{
  if(Parameters)
  {
    ErrMsg("input 'Parameters' of type 'TAuxParam3D*' is "<<
          "not set to NULL. This is usually done if you want values of FE "<<
          "functions during local assembling, for example in nonlinear "
          "problems. This is not yet supported. Exiting.");
    exit(1);
  }
  
  double hK;
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,k,l,l1,l3,n,m, N_LocalUsedElements;
  int N_Points, N_;
  int N_Test, N_Ansatz, N_Joints;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  FE3D TestElement, AnsatzElement;
  QuadFormula2D FaceQuadFormula;
  TQuadFormula2D *qf2;
  TJoint *joint;
  int **GlobalNumbers, **BeginIndex;
  int **RhsGlobalNumbers, **RhsBeginIndex;
  int **TestGlobalNumbers, **TestBeginIndex;
  int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta;
  double *t, *s;
  double xf, yf, zf;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *Param[MaxN_QuadPoints_3D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts3D];
  BaseFunct3D LocBF[N_BaseFuncts3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, begin, end, middle;
  int *TestDOF, *AnsatzDOF;
  double *Entries;
  int *ColInd, *RowPtr;
  double *RHS, *MatrixRow;
  double t0, t1, t2;
  BoundCond Cond0;
  BoundCondFunct3D *BoundaryCondition;
  BoundValueFunct3D *BoundaryValue;
  //TOutput3D *Output;
  TNodalFunctional3D *nf;
  RefTrans3D reftrans;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_BaseFunctions3D];
  int *EdgeDOF, N_EdgeDOF;
  double *FaceWeights;
  double **JointValues, *JointValue;

  // static bool *SecondDer = NULL;
  bool *SecondDer;
  double LinComb[4];

  const int *TmpFV, *TmpLen;
  int MaxLen;
  double xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3;
  double nx, ny, nz;
  
  double time1, time2;
  
  bool OuterBoundary;

  double time_total = GetTime();
  double time_all = 0;

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFunct3D *BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  int *N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    for(i=0;i<n_sqmatrices;i++)
    {
      fespace = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                               // endif n_sqmatrices

  if(n_matrices)
  {
    TestGlobalNumbers = new int* [n_matrices];
    AnsatzGlobalNumbers = new int* [n_matrices];
    TestBeginIndex = new int* [n_matrices];
    AnsatzBeginIndex = new int* [n_matrices];
    for(i=0;i<n_matrices;i++)
    {
      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetTestSpace();
      TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
      TestBeginIndex[i] = fespace->GetBeginIndex();

      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetAnsatzSpace();
      AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
      AnsatzBeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                               // endif n_matrices
  if(n_rhs)
  {
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions3D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions3D;

  }                                               // endif n_rhs
  
  // 20 <= number of term in bilinear form
  // DUE NOTE CHANGE 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
  aux = new double [MaxN_QuadPoints_3D*20];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;
  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions3D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions3D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions3D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions3D;
  }                                               // endif N_AllMatrices
  SecondDer = DiscreteForm3D->GetNeeds2ndDerivatives();
  

  // all spaces use same Coll
  TCollection *Coll = fespaces[0]->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  // ########################################################################
  // loop over all cells
  // ########################################################################
  for(i=0;i<N_Cells;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    
    int n_faces = cell->GetN_Faces(); // number of faces of this cell
    cell->SetNormalOrientation();
    for (int ijo = 0; ijo<n_faces; ijo++)
    {
      if (n_faces==4)
      {                                   // tetrahedra
        TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[ijo] = 
            cell->GetNormalOrientation(ijo);
      }
      else if (n_faces==6)
      {                                   // hexahedra
        TDatabase::ParamDB->NORMAL_ORIENTATION_HEXA[ijo] = 
            cell->GetNormalOrientation(ijo);
      }
    }
    
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
      // cases 4 and 5 are specially treated for special problems
      case 0:                                     // diameter
      case 4:                                     // diameter
      case 5:                                     // piecewise constant array
        hK = cell->GetDiameter();
        break;
      // case 3 is specially treated for special problem
      case 1:                                     // with reference map
      case 3:                                     // measure
        hK = cell->GetMeasure();
        hK = pow(hK,1.0/3.0);
        break;
      case 2:                                     // shortest edge
        hK = cell->GetShortestEdge();
        break;
      default:                                   // diameter
        hK = cell->GetDiameter();
        OutPut("CELL_MEASURE " << TDatabase::ParamDB->CELL_MEASURE <<
             " not available, set CELL_MEASURE: 0 !!!" << endl);
        TDatabase::ParamDB->CELL_MEASURE = 0;
        break;
    }
    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE3D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    //OutPut("CELL " << i << endl);
    reftrans = TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      Coll, cell, SecondDer,
      N_Points, xi, eta, zeta, weights,
      X, Y, Z, AbsDetjk);
    
    // this could provide values of FE functions during the local assemble
    // routine, not yet supported.
    //Parameters->GetParameters(N_Points, cell, i, xi, eta, zeta,
    //  X, Y, Z, Param);
    
    // use DiscreteForm to assemble a few matrices and
    // right-hand sides at once
    if(DiscreteForm3D)
      DiscreteForm3D->GetLocalForms(N_Points, weights, AbsDetjk,
        hK, X, Y, Z,
        LocN_BF, LocBF,
        Param, AuxArray,
        cell,
        N_AllMatrices, n_rhs,
        LocMatrices, LocRhs);
    else
    {
      ErrMsg("Assemble3D_mixed: no DiscreteForm3D given. Exit");
      exit(0);
    }
   

    time1 = GetTime();
    // ####################################################################
    // add local matrices to global matrices (ansatz == test)
    // ####################################################################

    for(j=0;j<n_sqmatrices;j++)
    {
      // find space for this bilinear form
      fespace = sqmatrices[j]->GetFESpace();
      CurrentElement = fespace->GetFE3D(i, cell);
      N_ = N_BaseFunct[CurrentElement];
      Matrix = Matrices+j*MaxN_BaseFunctions3D;
      Entries = sqmatrices[j]->GetEntries();
      RowPtr = sqmatrices[j]->GetRowPtr();
      ColInd = sqmatrices[j]->GetKCol();

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = GlobalNumbers[j] + BeginIndex[j][i];
      if (TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
      {
        ActiveBound = DirichletBound = sqmatrices[j]->GetN_Rows();
      }
      // add local matrix to global
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          // for all dof
          for(k=0;k<N_;k++)
          {

            begin = RowPtr[l];
            end = RowPtr[l+1]-1;
            l1 = DOF[k];
            if(ColInd[begin] == l1)
            {
              Entries[begin] += MatrixRow[k];
            }
            else
            {
              begin++;
              middle = (begin+end)/2;
              while(middle>begin)
              {
                if(l1 == ColInd[middle]) break;
                if(l1 < ColInd[middle]) end = (begin+end)/2;
                if(l1 > ColInd[middle]) begin = (begin+end)/2;
                middle = (begin+end)/2;
              }
              if(ColInd[middle] == l1)
              {
                Entries[middle] += MatrixRow[k];
              }
              else
              {
                if(ColInd[middle+1] == l1)
                {
                  Entries[middle+1] += MatrixRow[k];
                }
              }
            }
          }                                       // endfor k
        }                                         // endif l
        else
        {
          // Dirichlet node
          n=RowPtr[l];
          if(ColInd[n]==l)
          {
            Entries[n]=1.0;
          }
        }
      }                                           // endfor m
    }                                             // endfor j
    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(j=0;j<n_matrices;j++)
    {
      TestElement = ((TFESpace3D *) matrices[j]->GetStructure()->
        GetTestSpace())->GetFE3D(i, cell);
      AnsatzElement = ((TFESpace3D *) matrices[j]->GetStructure()->
        GetAnsatzSpace())->GetFE3D(i, cell);

      N_Test = N_BaseFunct[TestElement];
      N_Ansatz = N_BaseFunct[AnsatzElement];

      Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions3D;

      Entries = matrices[j]->GetEntries();
      RowPtr = matrices[j]->GetRowPtr();
      ColInd = matrices[j]->GetKCol();

      TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
      AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];

      int ActiveBound = ((TFESpace3D *) matrices[j]->GetStructure()->
        GetTestSpace())->GetActiveBound();
      
      // add local matrix to global
      for(m=0;m<N_Test;m++)
      {
        l=TestDOF[m];
        if(l >= ActiveBound)
          continue;
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        for(k=0;k<N_Ansatz;k++)
        {
          l1 = AnsatzDOF[k];
          begin = RowPtr[l];
          end = RowPtr[l+1]-1;
          middle = (begin+end)/2;
          while(middle>begin)
          {
            if(l1 == ColInd[middle]) break;
            if(l1 < ColInd[middle]) end = (begin+end)/2;
            if(l1 > ColInd[middle]) begin = (begin+end)/2;
            middle = (begin+end)/2;
          }
          if(ColInd[middle] == l1)
          {
            Entries[middle] += MatrixRow[k];
          }
          else
          {
            if(ColInd[middle+1] == l1)
            {
              Entries[middle+1] += MatrixRow[k];
            }
          }
        }                                         // endfor k
      }                                           // endfor m
    }                                             // endfor j
    time2 = GetTime();
    time_all += time2-time1;
    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      //OutPut("rhs " << j << endl);
      fespace = ferhs[j];
      ActiveBound = fespace->GetActiveBound();
      CurrentElement = fespace->GetFE3D(i, cell);

      N_ = N_BaseFunct[CurrentElement];

      local_rhs = righthand+j*MaxN_BaseFunctions3D;
      RHS = rhs[j];
      // find space for this linear form

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

      // add local right-hand side to the global one
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m];
        }                                         // endif l
      }                                           // endfor m

      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase3D::GetFE3D(CurrentElement);
      nf = ele->GetNodalFunctional3D();
      TFEDesc3D *FEDesc_Obj = ele->GetFEDesc3D();
      
      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Faces();
      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        OuterBoundary = false;

        if(joint->GetType() == BoundaryFace ||
          joint->GetType() == IsoBoundFace)
          OuterBoundary = true;

       
        if(OuterBoundary)
        {
          cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

          t0 = 1.0/TmpLen[m];

          // coordines of center are computed, this is where boundary 
          // conditions are taken
          xf = 0; yf = 0; zf = 0; 
          for (l1=0;l1<TmpLen[m];l1++)
          {
            cell->GetVertex(TmpFV[m*MaxLen+l1])
              ->GetCoords(X[l1], Y[l1], Z[l1]);
            //LinComb[l1] = t0;
            xf += t0*X[l1];
            yf += t0*Y[l1];
            zf += t0*Z[l1];
          }

          // the face gets the b.c. which is valid at its center
          BoundaryCondition(xf, yf, zf, Cond0);

          switch(Cond0)
          {
            case DIRICHLET:
              nf->GetPointsForFace(m, N_Points, xi, eta, zeta);
              TFEDatabase3D::GetOrigFromRef(reftrans, N_Points, xi, eta, zeta,
                                            X, Y, Z, AbsDetjk);
              
              for(l1=0;l1<N_Points;l1++)
              {
                BoundaryValue(X[l1], Y[l1], Z[l1], PointValues[l1]);
              }

              nf->GetFaceFunctionals(Coll, cell, m, PointValues,
                                     FunctionalValues);
              EdgeDOF = FEDesc_Obj->GetJointDOF(m);
              N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              
              for(l=0;l<N_EdgeDOF;l++)
              {
                RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
              }
              break;

            case NEUMANN:
              // cout << "Neumann condition in Assemble3D" << endl;
              l = TFEDatabase3D::GetPolynomialDegreeFromFE3D
                (CurrentElement);
              switch(TmpLen[m])
              {
                case 3:
                  // triangular face
                  FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*l);
                  FaceQuadFormula = Gauss3Tria;
                  break;
                case 4:
                  // quadrilateral face
                  FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*l);
                  break;
              }
              qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
              qf2->GetFormulaData(N_Points, FaceWeights, t, s);
              // generate data on reference mesh cell for the 2d face of 3d cell
              TFEDatabase3D::GetBaseFunct3DFromFE3D(CurrentElement)
                ->MakeRefElementData(FaceQuadFormula);
              // values of base functions in all quadrature points on face
              JointValues = TFEDatabase3D::GetJointValues3D
                (BaseFuncts[CurrentElement], FaceQuadFormula, m);
              TFEDatabase3D::GetBaseFunct3D(BaseFuncts[CurrentElement])
                ->ChangeBF(Coll, cell, N_Points, JointValues);

              switch(TmpLen[m])
              {
                case 3:
                  xc1 = X[1] - X[0];
                  xc2 = X[2] - X[0];

                  yc1 = Y[1] - Y[0];
                  yc2 = Y[2] - Y[0];

                  zc1 = Z[1] - Z[0];
                  zc2 = Z[2] - Z[0];

                  // normal vector
                  nx = yc1*zc2 - zc1*yc2;
                  ny = zc1*xc2 - xc1*zc2;
                  nz = xc1*yc2 - yc1*xc2;
                  // determinant of reference trafo
                  t2 = sqrt(nx*nx + ny*ny + nz*nz);

                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    t0 = t[l];
                    t1 = s[l];
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = 1-t0-t1;
                    LinComb[1] = t0;
                    LinComb[2] = t1;

                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                      +LinComb[2]*X[2];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                      +LinComb[2]*Y[2];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                      +LinComb[2]*Z[2];

                    /*if(OuterBoundary)
                      BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                             X, Y, Z, Param1, Param2,
                                             xf, yf, zf, t0, t1);
                    */
                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    // cout << t1 << endl;
                    t0 *= FaceWeights[l]*t2;
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += t0*JointValue[k];
                  }                               // endfor l
                  break;

                case 4:
                  xc1=(-X[0] + X[1] + X[2] - X[3]) * 0.25;
                  xc2=(-X[0] - X[1] + X[2] + X[3]) * 0.25;
                  xc3=( X[0] - X[1] + X[2] - X[3]) * 0.25;

                  yc1=(-Y[0] + Y[1] + Y[2] - Y[3]) * 0.25;
                  yc2=(-Y[0] - Y[1] + Y[2] + Y[3]) * 0.25;
                  yc3=( Y[0] - Y[1] + Y[2] - Y[3]) * 0.25;

                  zc1=(-Z[0] + Z[1] + Z[2] - Z[3]) * 0.25;
                  zc2=(-Z[0] - Z[1] + Z[2] + Z[3]) * 0.25;
                  zc3=( Z[0] - Z[1] + Z[2] - Z[3]) * 0.25;

                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    t0 = 0.5*(t[l]+1);
                    t1 = 0.5*(s[l]+1);
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = (1-t0)*(1-t1);
                    LinComb[1] = t0*(1-t1);
                    LinComb[2] = t0*t1;
                    LinComb[3] = (1-t0)*t1;

                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                      +LinComb[2]*X[2] + LinComb[3]*X[3];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                      +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                      +LinComb[2]*Z[2] + LinComb[3]*Z[3];

                    /*if(OuterBoundary)
                      BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                                             X, Y, Z, Param1, Param2,
                                             xf, yf, zf, t0, t1);
                    */
                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    nx = (yc1+s[l]*yc3)*(zc2+t[l]*zc3)
                      -(zc1+s[l]*zc3)*(yc2+t[l]*yc3);
                    ny = (zc1+s[l]*zc3)*(xc2+t[l]*xc3)
                      -(xc1+s[l]*xc3)*(zc2+t[l]*zc3);
                    nz = (xc1+s[l]*xc3)*(yc2+t[l]*yc3)
                      -(yc1+s[l]*yc3)*(xc2+t[l]*xc3);
                    t1 = nx*nx+ny*ny+nz*nz;
                    // cout << t1 << endl;
                    t0 *= FaceWeights[l]*sqrt(t1);
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += t0*JointValue[k];
                  }                               // endfor l
                  break;
              }
              TFEDatabase3D::GetBaseFunct3D(BaseFuncts[CurrentElement])
                ->ChangeBF(Coll, cell, N_Points, JointValues);
              break;
          }                                       // endswitch Cond0
        }                                         // endif
      }                                           // endfor m
    }                                             // endfor j
    // cout << "end i:" << i << endl;

  }                                               // endfor i (loop over all cells)
  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  
  // ####################################################################
  // write coupling into matrix
  // ####################################################################
  

  if(n_sqmatrices)
  {
    delete GlobalNumbers;
    delete BeginIndex;
  }

  if(n_matrices)
  {
    delete AnsatzGlobalNumbers;
    delete AnsatzBeginIndex;
    delete TestGlobalNumbers;
    delete TestBeginIndex;
  }

  if(n_rhs)
  {
    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }


  if(N_AllMatrices)
  {
    delete LocMatrices;
    delete Matrices[0];
    delete Matrices;
  }
  delete AuxArray[0];

  time_total = GetTime() - time_total;
}
