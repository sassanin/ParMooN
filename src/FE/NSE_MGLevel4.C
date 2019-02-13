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
// @(#)NSE3DMGLevel4.C        1.10 07/03/00
//
// Class:       TNSE3DMGLevel4
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system
//              type 2 (A, B1, B2, B1T, B2T)
//
// Author:      Volker John 25.08.1999
//
// History:     24.08.1999 start of implementation
//              25.08.1999 CellVanka
//              25.08.1999 step length control
//              25.08.1999 ExactSolve
//              26.08.1999 NodalVanka
//
// =======================================================================

#include <NSE_MGLevel4.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Solver.h>
#ifdef __2D__
  #include <FESpace2D.h>
  #include <FEDatabase2D.h>
#endif
#ifdef __3D__
  #include <FESpace3D.h>
  #include <FEDatabase3D.h>
  #include <Joint.h>
  #include <Edge.h>
  #include <BaseCell.h>  
  #include <Collection.h>
#endif

#include <stdlib.h>
#include <string.h>

#include <LinAlg.h>
#include <Solver.h>
#include <ItMethod.h>
#include <FgmresIte.h>
#include <DirectSolver.h>

#include <LocalProjection.h>

#ifdef _MPI 
#include <ParFECommunicator3D.h>
#include <ParDirectSolver.h>
#endif

/** constructor */
#ifdef __2D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix2D *a11, TSquareMatrix2D *a12,
                               TSquareMatrix2D *a21, TSquareMatrix2D *a22,
                               TMatrix2D *b1, TMatrix2D *b2,
                               TMatrix2D *b1t, TMatrix2D *b2t,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
			       int *dw)
#endif
#ifdef __3D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix3D *a11, TSquareMatrix3D *a12,
                               TSquareMatrix3D *a13, TSquareMatrix3D *a21,
                               TSquareMatrix3D *a22, TSquareMatrix3D *a23,
                               TSquareMatrix3D *a31, TSquareMatrix3D *a32,
                               TSquareMatrix3D *a33,
                               TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,
                               TMatrix3D *b1t, TMatrix3D *b2t, TMatrix3D *b3t,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
                               int *dw
			      #ifdef _MPI
 ,TParFECommunicator3D *parComm_U, TParFECommunicator3D *parComm_P
#endif
    
  )
#endif
  : TNSE_MGLevel(level, f1, u1, n_aux, al,
                 velocity_space, pressure_space, Coll)
{
  int i;
  double *aux;

  A11 = a11;
  A12 = a12;
  A21 = a21;
  A22 = a22;
  StructureA = A11->GetMatrixStructure();
  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
#ifdef __3D__
  A13 = a13;
  A23 = a23;
  A31 = a31;
  A32 = a32;
  A33 = a33;
  A13Entries = A13->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
#endif

  B1T = b1t;
  B2T = b2t;
  StructureBT = B1T->GetStructure();
  BTRowPtr = StructureBT->GetRowPtr();
  BTKCol = StructureBT->GetKCol();
  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  B1 = b1;
  B2 = b2;
  StructureB = B1->GetStructure();
  BRowPtr = StructureB->GetRowPtr();
  BKCol = StructureB->GetKCol();
  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

#ifdef __3D__
  B3T = b3t;
  B3TEntries = B3T->GetEntries();
  B3 = b3;
  B3Entries = B3->GetEntries();
#endif

  C = NULL;
  StructureC = NULL;
  CRowPtr = NULL;
  CKCol = NULL;
  CEntries = NULL;

  USpace = A11->GetFESpace();
#ifdef __2D__
  PSpace = (TFESpace2D *)StructureB->GetTestSpace();
#endif
#ifdef __3D__
  PSpace = (TFESpace3D *)StructureB->GetTestSpace();
#endif

  N_Active = USpace->GetActiveBound();
  HangingNodeBound = USpace->GetHangingBound();
  N_Dirichlet = USpace->GetN_Dirichlet();

  N_UDOF = USpace->GetN_DegreesOfFreedom();
  N_PDOF = PSpace->GetN_DegreesOfFreedom();

  N_DOF = GEO_DIM*N_UDOF+N_PDOF;

  U1 = u1;
  U2 = u1 + N_UDOF;
  P  = u1 + GEO_DIM*N_UDOF;

  Rhs1 = f1;
  Rhs2 = f1 + N_UDOF;
  RhsP = f1 + GEO_DIM*N_UDOF;

#ifdef __3D__
  U3 = u1 + 2*N_UDOF;
  Rhs3 = f1 + 2*N_UDOF;
#endif

  N_Aux = n_aux;
  Aux = new double* [N_Aux];
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Type = 4;

  alpha = al[0];
  downwind = dw;
  
#ifdef _MPI
  ParCommU =  parComm_U;
  ParCommP =  parComm_P;
#endif
}

/** constructor with matrix C */
#ifdef __2D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix2D *a11, TSquareMatrix2D *a12,
                               TSquareMatrix2D *a21, TSquareMatrix2D *a22,
                               TMatrix2D *b1, TMatrix2D *b2,
                               TMatrix2D *b1t, TMatrix2D *b2t,
                               TMatrix2D *c,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
                   int *dw)
#endif
#ifdef __3D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix3D *a11, TSquareMatrix3D *a12,
                               TSquareMatrix3D *a13, TSquareMatrix3D *a21,
                               TSquareMatrix3D *a22, TSquareMatrix3D *a23,
                               TSquareMatrix3D *a31, TSquareMatrix3D *a32,
                               TSquareMatrix3D *a33,
                               TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,
                               TMatrix3D *b1t, TMatrix3D *b2t, TMatrix3D *b3t,
                               TMatrix3D *c,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
                               int *dw
  )
#endif
  : TNSE_MGLevel(level, f1, u1, n_aux, al,
                 velocity_space, pressure_space, Coll)
{
  int i;
  double *aux;

  A11 = a11;
  A12 = a12;
  A21 = a21;
  A22 = a22;
  StructureA = A11->GetMatrixStructure();
  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
#ifdef __3D__
  A13 = a13;
  A23 = a23;
  A31 = a31;
  A32 = a32;
  A33 = a33;
  A13Entries = A13->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
#endif

  B1T = b1t;
  B2T = b2t;
  StructureBT = B1T->GetStructure();
  BTRowPtr = StructureBT->GetRowPtr();
  BTKCol = StructureBT->GetKCol();
  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  B1 = b1;
  B2 = b2;
  StructureB = B1->GetStructure();
  BRowPtr = StructureB->GetRowPtr();
  BKCol = StructureB->GetKCol();
  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

#ifdef __3D__
  B3T = b3t;
  B3TEntries = B3T->GetEntries();
  B3 = b3;
  B3Entries = B3->GetEntries();
#endif

  C = c;
  StructureC = C->GetStructure();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  USpace = A11->GetFESpace();
#ifdef __2D__
  PSpace = (TFESpace2D *)StructureB->GetTestSpace();
#endif
#ifdef __3D__
  PSpace = (TFESpace3D *)StructureB->GetTestSpace();
#endif

  N_Active = USpace->GetActiveBound();
  HangingNodeBound = USpace->GetHangingBound();
  N_Dirichlet = USpace->GetN_Dirichlet();

  N_UDOF = USpace->GetN_DegreesOfFreedom();
  N_PDOF = PSpace->GetN_DegreesOfFreedom();

  N_DOF = GEO_DIM*N_UDOF+N_PDOF;

  U1 = u1;
  U2 = u1 + N_UDOF;
  P  = u1 + GEO_DIM*N_UDOF;

  Rhs1 = f1;
  Rhs2 = f1 + N_UDOF;
  RhsP = f1 + GEO_DIM*N_UDOF;

#ifdef __3D__
  U3 = u1 + 2*N_UDOF;
  Rhs3 = f1 + 2*N_UDOF;
#endif

  N_Aux = n_aux;
  Aux = new double* [N_Aux];
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Type = 4;

  alpha = al[0];
  downwind = dw;

}

/** destructor */
TNSE_MGLevel4::~TNSE_MGLevel4()
{
} // ~TNSE3DMGLevel4

/** calculate defect */
void TNSE_MGLevel4::Defect(double *u1, double *f1,  double *d1, double &res)
{
  
#ifdef _MPI
  int i, rank, *MasterOfDofU,*MasterOfDofP, dof,numThreads= TDatabase::ParamDB->OMPNUMTHREADS;
  double res_global=0.;
  double res1 = 0.;
  
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
  MasterOfDofU = ParCommU->GetMaster();
  MasterOfDofP = ParCommP->GetMaster();

#endif

#ifdef __2D__

    CoupledDefect(A11,A12,A21,A22,B1,B2,B1T,B2T,u1,f1,d1);

    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(d1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif
#ifdef __3D__
//  u1 = x , f1 = b , d1 = r::
    
    CoupledDefect(A11,A12,A13,A21,A22,A23,A31,A32,A33,
                B1,B2,B3,B1T,B2T,B3T,u1,f1,d1);

#ifdef _MPI
    ParCommU->CommUpdate(d1);
    ParCommP->CommUpdate(d1 + GEO_DIM*N_UDOF);
#endif
    
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector3D(d1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif

  memset(d1+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
  memset(d1+N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#ifdef __3D__
  memset(d1+2*N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#endif
      
#ifdef _MPI
  for(i=0; i<N_UDOF; i++)
    if(MasterOfDofU[i] == rank)
  {
      res1 += d1[i	      ]*d1[i	  	];
      res1 += d1[i  +   N_UDOF]*d1[i  +   N_UDOF];
#ifdef __3D__
      res1 += d1[i  + 2*N_UDOF]*d1[i  + 2*N_UDOF];
#endif    
  }

  for(i=0;i<N_PDOF;i++)
    if(MasterOfDofP[i] == rank)
#ifdef __3D__
      res1 += d1[i + 3*N_UDOF]*d1[i + 3*N_UDOF];
#endif
#ifdef __2D__
      res1 += d1[i + 2*N_UDOF]*d1[i + 2*N_UDOF];
#endif
  
  TDatabase::ParamDB->time_communication_start = MPI_Wtime();
 
  MPI_Allreduce(&res1, &res_global, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  
  TDatabase::ParamDB->time_communication_end =  MPI_Wtime();
  TDatabase::ParamDB->time_communication += TDatabase::ParamDB->time_communication_end - TDatabase::ParamDB->time_communication_start;
  
  
  res = sqrt(res_global);
#else
  res = sqrt(Ddot(N_DOF,d1,d1));
#endif
}

/** correct Dirichlet and hanging nodes */
void TNSE_MGLevel4::CorrectNodes(double *u1)
{
  int i,j,k, index;
  double s, t, u, *u2, *u3;

  u2 = u1+N_UDOF;
#ifdef __3D__
  u3 = u2+N_UDOF;
#endif
  
#ifdef _MPI
  int *MasterOfDofU = ParCommU->GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

  // set Dirichlet nodes
  memset(u1+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
  memset(u2+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memset(u3+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
#endif

  // set hanging nodes
  j = ARowPtr[N_Active];
//   if(N_Active<HangingNodeBound)
//   {
    for(i=N_Active;i<HangingNodeBound;i++)
    {
      s = 0;
      t = 0;
  #ifdef __3D__
      u = 0;
  #endif
      k = ARowPtr[i+1];
      for(;j<k;j++)
      {
	index = AKCol[j];
	if(index != i)
	{
	  s -= A11Entries[j] * u1[index];
	  t -= A22Entries[j] * u2[index];
  #ifdef __3D__
	  u -= A22Entries[j] * u3[index];
  #endif
	}
      } // endfor j
//       #ifdef _MPI
//       if(MasterOfDofU[i]== rank)
//       #endif
      {
      u1[i] = s;
      u2[i] = t;
  #ifdef __3D__
      u3[i] = u;
  #endif
      }
    } // endfor i
  //}
}

/** cellwise Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel4::CellVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother,int N_Levels)
{
  

#ifdef __2D__
  const int RhsDim =  3*MaxN_BaseFunctions2D;
  TFE2D *UEle, *PEle;
  TSquareMatrix2D *sqmatrix[1];
#endif
#ifdef __3D__
  const int RhsDim =  4*MaxN_BaseFunctions3D;
  TFE3D *UEle, *PEle;
  TSquareMatrix3D *sqmatrix[1];
#endif
  int i,j,k,l,m, N_Cells, N_LocalDOFs, ii;
  int j1, j2, j3, j4, k1, k2, k3;
  double value, value1, value2, value3;
  double value11,value12,value13,value21,value22;
  double value23,value31,value32,value33;
  double *uold, *pold;
  TCollection *Coll;
  double System[RhsDim*RhsDim];
  double Rhs[RhsDim], sol[RhsDim];
  int *UGlobalNumbers, *UBeginIndex, *UDOFs, UDOF, N_U;
  int *PGlobalNumbers, *PBeginIndex, *PDOFs, PDOF, N_P;
  int N_LocalDOF, verbose;
  int begin, end, ActiveBound, begin1, end1;
  double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  TBaseCell *Cell;
  double *u2, *u3, *p, *rhs2, *rhs3, *rhsp;
  TItMethod *itmethod = NULL;
  int LargestDirectSolve = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  MatVecProc *MatVect=MatVectFull;
  DefectProc *Defect=DefectFull;
  TSquareMatrix **matrix= (TSquareMatrix **)sqmatrix;
 
  TDatabase::ParamDB->INTERNAL_LOCAL_DOF = -1;
#ifdef __2D__
  sqmatrix[0] = (TSquareMatrix2D *)System;
#endif
#ifdef __3D__
  sqmatrix[0] = (TSquareMatrix3D *)System;
#endif

  if(VankaColl)
    Coll = VankaColl;
  else
    Coll = USpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();
  ActiveBound = USpace->GetActiveBound();

  PGlobalNumbers = PSpace->GetGlobalNumbers();
  PBeginIndex = PSpace->GetBeginIndex();

  // set pointers
  u2 = u1 + N_UDOF;
#ifdef __3D__
  u3 = u2 + N_UDOF;
#endif
  p  = u1 + GEO_DIM*N_UDOF;

  rhs2 = rhs1 + N_UDOF;
#ifdef __3D__
  rhs3 = rhs2 + N_UDOF;
#endif
  rhsp = rhs1 + GEO_DIM*N_UDOF;

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#endif

  SetHangingNodes(u1);

  // old values
  uold = aux;
  pold  = uold+GEO_DIM*N_UDOF;

  // save current solution on 'old' vectors
  memcpy(uold, u1, N_DOF*SizeOfDouble);

  // start of cell loop
  for(ii=0;ii<N_Cells;ii++)
  {
    //ii = downwind[i];
    
    Cell = Coll->GetCell(ii);
    
#ifdef _MPI
    if(Cell->IsHaloCell()){
      cout << "this should" << endl;
      continue;
      cout << "this shouldnt" << endl;
    }   
#endif
//    OutPut(i << downwind[i] << endl);
#ifdef __2D__
    UEle = TFEDatabase2D::GetFE2D(USpace->GetFE2D(ii, Cell));
    PEle = TFEDatabase2D::GetFE2D(PSpace->GetFE2D(ii, Cell));
#endif
#ifdef __3D__
    UEle = TFEDatabase3D::GetFE3D(USpace->GetFE3D(ii, Cell));
    PEle = TFEDatabase3D::GetFE3D(PSpace->GetFE3D(ii, Cell));
#endif

    // get local number of dof
    N_U = UEle->GetN_DOF();
    N_P = PEle->GetN_DOF();
    N_LocalDOF = GEO_DIM*N_U+N_P;

    // reset local systems
    /*if (N_LocalDOF > RhsDim)
    {
      OutPut(
    "TNSE_MGLevel4::CellVanka - Not enough memory in array Rhs!!!"
    << endl << "available " << RhsDim << " needed " <<
    N_LocalDOF << endl);
      exit(4711);
      }*/
    // reset local systems
    memset(System, 0, SizeOfDouble*N_LocalDOF*N_LocalDOF);
    if (N_LocalDOF > LargestDirectSolve)
    {
      // size of local system has changed
      if (N_LocalDOF != TDatabase::ParamDB->INTERNAL_LOCAL_DOF)
      {
        // itmethod exists already
        if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
          delete itmethod;
        // allocate new itmethod
        itmethod = new TFgmresIte(MatVect, Defect, NULL, 0, N_LocalDOF, 1);
        TDatabase::ParamDB->INTERNAL_LOCAL_DOF = N_LocalDOF;
      }
    }

    UDOFs = UGlobalNumbers+UBeginIndex[ii];
    PDOFs = PGlobalNumbers+PBeginIndex[ii];

    // fill local matrix
    for(j=0;j<N_U;j++)
    {
      j1 = j;
      j2 = j+N_U;
#ifdef __3D__
      j3 = j2+N_U;
#endif
      UDOF = UDOFs[j];

      // A block
      begin = ARowPtr[UDOF];
      end = ARowPtr[UDOF+1];

      Rhs[j1] = rhs1[UDOF];
      Rhs[j2] = rhs2[UDOF];
#ifdef __3D__
      Rhs[j3] = rhs3[UDOF];
#endif

      for(k=begin;k<end;k++)
      {
        l = AKCol[k];

        value11 = A11Entries[k];
        value12 = A12Entries[k];
        value21 = A21Entries[k];
        value22 = A22Entries[k];
#ifdef __3D__
        value13 = A13Entries[k];
        value23 = A23Entries[k];
        value31 = A31Entries[k];
        value32 = A32Entries[k];
        value33 = A33Entries[k];
#endif

#ifdef __2D__
        if (UDOF>=ActiveBound) // Dirichlet node
          value12 = value21 = 0;

        Rhs[j1] -= value11*u1[l]+value12*u2[l];
        Rhs[j2] -= value21*u1[l]+value22*u2[l];
#endif
#ifdef __3D__
        if (UDOF>=ActiveBound) // Dirichlet node
          value12 = value13 = value21 = value23 = value31 = value32 = 0;

        Rhs[j1] -= value11*u1[l]+value12*u2[l]+value13*u3[l];
        Rhs[j2] -= value21*u1[l]+value22*u2[l]+value23*u3[l];
        Rhs[j3] -= value31*u1[l]+value32*u2[l]+value33*u3[l];
#endif

        for(m=0;m<N_U;m++)
          if(UDOFs[m]==l)
          {
            // column belongs to local system
            k1 = m*N_LocalDOF;
            k2 = (m+N_U)*N_LocalDOF;
            System[k1+j1] = value11;
            System[k2+j1] = value12;
            System[k1+j2] = value21;
            System[k2+j2] = value22;
#ifdef __3D__
            k3 = (m+2*N_U)*N_LocalDOF;
            System[k3+j1] = value13;
            System[k3+j2] = value23;
            System[k1+j3] = value31;
            System[k2+j3] = value32;
            System[k3+j3] = value33;
#endif
            break;
          }
      } // endfor k

      if(UDOF<ActiveBound)  // active dof
      {
        // transpose(B) block for non-Dirichlet nodes
        begin = BTRowPtr[UDOF];
        end = BTRowPtr[UDOF+1];

        for(k=begin;k<end;k++)
        {
          l = BTKCol[k];
          value1 = B1TEntries[k];
          value2 = B2TEntries[k];
#ifdef __3D__
          value3 = B3TEntries[k];
#endif
          value = p[l];
          Rhs[j1] -= value1*value;
          Rhs[j2] -= value2*value;
#ifdef __3D__
          Rhs[j3] -= value3*value;
#endif

          for(m=0;m<N_P;m++)
            if(PDOFs[m]==l)
            {
              // column belongs to local system
              k1 = (m+GEO_DIM*N_U)*N_LocalDOF;
              System[k1+j1] = value1;
              System[k1+j2] = value2;
#ifdef __3D__
              System[k1+j3] = value3;
#endif
              break;
            }

        } // endfor k
      } // endif UDOF<ActiveBound
    } // endfor j

    for(j=0;j<N_P;j++)
    {
      j1 = j+GEO_DIM*N_U;
      PDOF = PDOFs[j];
      begin = BRowPtr[PDOF];
      end = BRowPtr[PDOF+1];
      Rhs[j1] = rhsp[PDOF];

      for(k=begin;k<end;k++)
      {
        l=BKCol[k];
        value1 = B1Entries[k];
        value2 = B2Entries[k];
#ifdef __3D__
        value3 = B3Entries[k];
#endif
        Rhs[j1] -= value1*u1[l];
        Rhs[j1] -= value2*u2[l];
#ifdef __3D__
        Rhs[j1] -= value3*u3[l];
#endif

        for(m=0;m<N_U;m++)
          if(UDOFs[m]==l)
          {
            // column belongs to local system
            k1 = m;
            k2 = m+N_U;
#ifdef __3D__
            k3 = k2 + N_U;
#endif
            System[k1*N_LocalDOF+j1] = value1;
            System[k2*N_LocalDOF+j1] = value2;
#ifdef __3D__
            System[k3*N_LocalDOF+j1] = value3;
#endif
            break;
          }
      } // endfor k
    } // endfor j

    if(C) // handle matrix C if present
    {
      for(j=0;j<N_P;j++)
      {
        j1 = j+GEO_DIM*N_U;
        PDOF = PDOFs[j];
        begin = CRowPtr[PDOF];
        end = CRowPtr[PDOF+1];
        for(k=begin;k<end;k++)
        {
          l = CKCol[k];
          value = -CEntries[k]; // minus is right sign
          Rhs[j1] -= value*p[l];

          for(m=0;m<N_P;m++)
            if(PDOFs[m] == l)
            {
              // column belongs to local system
              System[(m+GEO_DIM*N_U)*N_LocalDOF + j1] = value;
              break;
            } // endif
        } // endfor k
      } // endfor j
    } // endif C


    // solve local system
    if (smoother==1 && !C) // no diagonal Vanka for matrix C
    {
      // diagonal Vanka
#ifdef __2D__
      SolveDiagonalVanka2D(System, Rhs, N_U, N_P, N_LocalDOF);
#endif
#ifdef __3D__
      SolveDiagonalVanka3D(System, Rhs, N_U, N_P, N_LocalDOF);
#endif
    }
    else
    {
      // full Vanka
//       cout << "full vanka dof :: " << N_LocalDOF << endl;
      
      if (N_LocalDOF > LargestDirectSolve)
      {
        memset(sol,0,N_LocalDOF*SizeOfDouble);
        verbose =  TDatabase::ParamDB->SC_VERBOSE;
        TDatabase::ParamDB->SC_VERBOSE = -1;
        itmethod->Iterate(matrix,NULL,sol,Rhs);
        TDatabase::ParamDB->SC_VERBOSE = verbose;
        memcpy(Rhs, sol, N_LocalDOF*SizeOfDouble);
      }
      else
      {
        SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
      }
    }
#ifdef __3D__
    j1 = 2*N_U;
#endif
    for(j=0;j<N_U;j++)
    {
      l = UDOFs[j];
      u1[l] += damp*Rhs[j];
      u2[l] += damp*Rhs[j+N_U];
#ifdef __3D__
      u3[l] += damp*Rhs[j+j1];
#endif  
    }

    j1 = GEO_DIM*N_U;
    for(j=0;j<N_P;j++)
    {
      l = PDOFs[j];
      p[l] += damp*Rhs[j+j1];
    }
  } // endfor loop over cells

  // apply damping
  if (fabs(1-alpha)>1e-3)
    for(j=0;j<N_DOF;j++)
       u1[j] = uold[j]+alpha*(u1[j]-uold[j]);

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#endif

  SetHangingNodes(u1);

  if(TDatabase::ParamDB->INTERNAL_LOCAL_DOF > 0)
  {
    delete itmethod;
  }

} // end Vanka

/** nodal Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel4::NodalVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother, int N_Levels)
{

  
#ifdef _MPI
  TDatabase::ParamDB->time_vanka_start = MPI_Wtime();
#else
  TDatabase::ParamDB->time_vanka_start = GetTime();
#endif
  
#ifdef _MPI
  int rank, *MasterOfDofU,*MasterOfDofP;
  
  char *DofmarkerP = ParCommP->Get_DofMarker();
  char *DofmarkerU = ParCommU->Get_DofMarker();
  
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
  MasterOfDofU = ParCommU->GetMaster();
  MasterOfDofP = ParCommP->GetMaster();    
#endif
  
#ifdef __2D__
  const int MaxN_LocalU = 2*MaxN_BaseFunctions2D;
  const int SystemRhs = 3*MaxN_BaseFunctions2D;
  TSquareMatrix2D *sqmatrix[1];
#endif
#ifdef __3D__
  const int MaxN_LocalU = 4*MaxN_BaseFunctions3D;
  const int SystemRhs = 8*MaxN_BaseFunctions3D;
  TSquareMatrix3D *sqmatrix[1];
#endif
  int i,j,k,l,m;
  int order;
  int j1, j2, j3, j4, k1, k2, k3;
  double value, value1, value2, value3;
  double value11,value12,value13,value21,value22;
  double value23,value31,value32,value33;
  double *uold, *pold;
  double System[SystemRhs*SystemRhs];
  double Rhs[SystemRhs], sol[SystemRhs];
  int N_LocalDOF;
  int begin, end, HangingBound, begin1, end1, verbose;
  int UDOFs[MaxN_LocalU], UDOF, N_U, N_U2, N_UGEO;
  double *u2, *u3, *p, *rhs2, *rhs3, *rhsp;
  TItMethod *itmethod = NULL;
  double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  int LargestDirectSolve = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  MatVecProc *MatVect=MatVectFull;
  DefectProc *Defect=DefectFull;
  TSquareMatrix **matrix= (TSquareMatrix **)sqmatrix;
  
  double gar;

  TDatabase::ParamDB->INTERNAL_LOCAL_DOF = -1;
  
#ifdef __2D__
  sqmatrix[0] = (TSquareMatrix2D *)System;
#endif
#ifdef __3D__
  sqmatrix[0] = (TSquareMatrix3D *)System;
#endif

  HangingBound = USpace->GetHangingBound();

  // set pointers
  u2 = u1 + N_UDOF;
#ifdef __3D__
  u3 = u2 + N_UDOF;
#endif
  p  = u1 + GEO_DIM*N_UDOF;

  rhs2 = rhs1 + N_UDOF;
#ifdef __3D__
  rhs3 = rhs2 + N_UDOF;
#endif
  rhsp = rhs1 + GEO_DIM*N_UDOF;

  // set Dirichlet values
  
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#endif

  // old values
  uold = aux;
  pold = uold+GEO_DIM*N_UDOF;

  // save current solution on 'old' vectors
  memcpy(uold, u1, N_DOF*SizeOfDouble);
  
	  #ifdef _MPI      
	    for(order=0;order<N_PDOF;order++)	      
	  #else
	    for(i=0;i<N_PDOF;i++)    
	  #endif
	      {
	  #ifdef _MPI
	      if(TDatabase::ParamDB->DOF_Reorder == 1)
		i = re_pdof[N_PDOF-1-order];    	      
	      else
		i = order;	    
	  #endif			
		N_U = 0;

	  #ifdef _MPI    
	      
	      if(TDatabase::ParamDB->DOF_Average){
		if(DofmarkerP[i] == 'h' || DofmarkerP[i] == 'H')
		  continue;
	      }     
	      else{
		if(DofmarkerP[i] == 'h' || DofmarkerP[i] == 'H'  ||  DofmarkerP[i] == 's')
		  continue;
	      }
	      
	  #endif
	      // go through row i of B1 and B2
	      begin = BRowPtr[i];
	      end = BRowPtr[i+1];
	      value = rhsp[i];      // rhs of this pressure value
	    for(k=begin;k<end;k++)
	      { 
		l=BKCol[k]; 
		value1 = B1Entries[k];
		value2 = B2Entries[k];
	  #ifdef __3D__
		value3 = B3Entries[k];
	  #endif
	  #ifdef __2D__
		value -= value1*u1[l]+value2*u2[l]; // update rhs
	  #endif      
	  #ifdef __3D__
		value -= value1*u1[l]+value2*u2[l]+value3*u3[l]; // update rhs
	  #endif    

	  #ifdef _MPI   
	      if(TDatabase::ParamDB->DOF_Average){
		if(DofmarkerU[l] == 'h' || DofmarkerU[l] == 'H')
		  continue;
	      }
	      
	      else{
		if(DofmarkerU[l] == 'h' || DofmarkerU[l] == 'H' || DofmarkerU[l] == 's')
		  continue;
	      }
	  #endif

	      UDOFs[N_U] = l;
	      j1 = GEO_DIM*N_U;
	      System[j1] = value1;  // save values for local B
	      System[j1+1] = value2;
	  #ifdef __3D__
	      System[j1+2] = value3;
	  #endif  
	      N_U++;           // count # velo dof connected to the pressure dof
	      }                    // row done

	      N_U2 = 2 * N_U;
	      N_UGEO = GEO_DIM * N_U;
	      N_LocalDOF = N_UGEO +1;

	      if (N_LocalDOF > LargestDirectSolve)
	      {
		if (N_LocalDOF != TDatabase::ParamDB->INTERNAL_LOCAL_DOF)
		{
		  if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
		  {
		    delete itmethod;
		  }
		  itmethod = new TFgmresIte(MatVect, Defect, NULL, 0, N_LocalDOF, 1);
		  TDatabase::ParamDB->INTERNAL_LOCAL_DOF = N_LocalDOF;
		}
	      }
	      
	      memset(System+N_UGEO, 0, SizeOfDouble*(N_LocalDOF*N_LocalDOF-N_UGEO));

	      Rhs[N_LocalDOF-1] = value;  // set rhs

	      for (k=0;k<N_U;k++)         // copy local B to the right place
	      {
		j4 = GEO_DIM*k;
		System[k*N_LocalDOF+N_UGEO]=System[j4];
		System[(k+N_U)*N_LocalDOF+N_UGEO]=System[j4+1];
		#ifdef __3D__
		System[(k+N_U2)*N_LocalDOF+N_UGEO]=System[j4+2];
		#endif
	      }
	      
	      memset(System, 0, SizeOfDouble*N_UGEO);
	      
	      for(j=0;j<N_U;j++)
	      {
		
		j1 = j;
		j2 = j+N_U;
		#ifdef __3D__
		j3 = j2+N_U;
		#endif
		
		UDOF = UDOFs[j];

		// A block
		begin = ARowPtr[UDOF];
		end = ARowPtr[UDOF+1];

		Rhs[j1] = rhs1[UDOF];
		Rhs[j2] = rhs2[UDOF];
		#ifdef __3D__
		Rhs[j3] = rhs3[UDOF];
		#endif

		for(k=begin;k<end;k++)
		{
		  l = AKCol[k];
		  value11 = A11Entries[k];
		  value12 = A12Entries[k];
		  value21 = A21Entries[k];
		  value22 = A22Entries[k];
		  #ifdef __3D__
		  value13 = A13Entries[k];
		  value23 = A23Entries[k];
		  value31 = A31Entries[k];
		  value32 = A32Entries[k];
		  value33 = A33Entries[k];
		  #endif

		  #ifdef __2D__
		  if (UDOF>=HangingBound) // Dirichlet node
		    value21 = value12 = 0;

		  Rhs[j1] -= value11*u1[l]+value12*u2[l];
		  Rhs[j2] -= value21*u1[l]+value22*u2[l];
		  #endif
		  
		  #ifdef __3D__
		  if (UDOF>=HangingBound) // Dirichlet node
		    value12 = value13 = value21 = value23 = value31 = value32 = 0;

// 		  #ifdef _MPI
// 		  if(DofmarkerU[UDOF] != 'h' || DofmarkerU[UDOF] != 'H' )
// 		  #endif
		  {
		    Rhs[j1] -= value11*u1[l]+value12*u2[l]+value13*u3[l];
		    Rhs[j2] -= value21*u1[l]+value22*u2[l]+value23*u3[l];
		    Rhs[j3] -= value31*u1[l]+value32*u2[l]+value33*u3[l];
		  }
		  #endif

		  for(m=0;m<N_U;m++)
		    if(UDOFs[m]==l)
		    {
		      k1 = m*N_LocalDOF;
		      k2 = (m+N_U)*N_LocalDOF;

		      System[k1+j1] = value11;
		      System[k2+j1] = value12;
		      System[k1+j2] = value21;
		      System[k2+j2] = value22;
		      #ifdef __3D__
		      k3 = (m + 2*N_U)*N_LocalDOF;
		      System[k3+j1] = value13;
		      System[k3+j2] = value23;
		      System[k1+j3] = value31;
		      System[k2+j3] = value32;
		      System[k3+j3] = value33;
		      #endif
		      break;
		    }
		} // endfor k

		if(UDOF<HangingBound)  // active dof
		{
		  // transpose(B) block for non-Dirichlet nodes
		  begin = BTRowPtr[UDOF];
		  end = BTRowPtr[UDOF+1];

		  for(k=begin;k<end;k++)
		  {
		    l = BTKCol[k];
		    value1 = B1TEntries[k];
		    value2 = B2TEntries[k];
		    #ifdef __3D__
		    value3 = B3TEntries[k];
		    #endif
		    value = p[l];
		 
		    {
		      Rhs[j1] -= value1*value;
		      Rhs[j2] -= value2*value;
		      #ifdef __3D__
		      Rhs[j3] -= value3*value;
		      #endif
		    }

		    if(i==l)
		    {
		      j4 = N_UGEO*N_LocalDOF;
		      System[j4+j1] = value1;
		      System[j4+j2] = value2;
		      #ifdef __3D__
		      System[j4+j3] = value3;
		      #endif
		    }
		  } // endfor k
		} // endif UDOF<HangingBound
	      } // endfor j

	      if(C)
	      {
		// fill C block if present
		begin = CRowPtr[i];
		end = CRowPtr[i+1];
		for(k=begin;k<end;k++)
		{
		  l = CKCol[k];
		  value = -CEntries[k]; // minus is right sign
		  Rhs[N_LocalDOF-1] -= value*p[l];
		  if(l==i) // main diagonal
		    System[N_LocalDOF*N_LocalDOF-1] = value;
		} // endfor k
	      } // endif C
	      
//=================================================Start Solving local system ===========================================

#ifdef _MPI
  TDatabase::ParamDB->time_vanka_solve_start = MPI_Wtime();
#else
  TDatabase::ParamDB->time_vanka_solve_start = GetTime();
#endif
	      if (smoother==3 && !C) // no diagonal Vanka for matrix C
	      {
	  #ifdef __2D__
		// diagonal Vanka
		SolveDiagonalVanka2D(System,  Rhs, N_U, 1, N_LocalDOF);
	  #endif
	  #ifdef __3D__
		// diagonal Vanka
		SolveDiagonalVanka3D(System,  Rhs, N_U, 1, N_LocalDOF);
	  #endif
	      }
	      else
	      {
		// full Vanka
		if (N_LocalDOF > LargestDirectSolve)
		{
		  int iter =0;
		  
		  double check2 = 0.0;
		  for(iter=0;iter<N_LocalDOF;iter++)
		    check2 += Rhs[iter]*Rhs[iter];
		  
		  memset(sol,0,N_LocalDOF*SizeOfDouble);
		  verbose =  TDatabase::ParamDB->SC_VERBOSE;
		  TDatabase::ParamDB->SC_VERBOSE = -1;
		  
		  int tem = TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER;
		  TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER = 0;
		  
		  if(check2)
		  itmethod->Iterate(matrix,NULL,sol,Rhs);
			  
		  TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER=tem;
		  TDatabase::ParamDB->SC_VERBOSE = verbose;
		  memcpy(Rhs, sol, N_LocalDOF*SizeOfDouble);
		}
		else
		{           
		  SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
		}
	      }
	      
#ifdef _MPI
  TDatabase::ParamDB->time_vanka_solve_end = MPI_Wtime();
#else
  TDatabase::ParamDB->time_vanka_solve_end = GetTime();
#endif	      
	      
  TDatabase::ParamDB->time_vanka_solve += TDatabase::ParamDB->time_vanka_solve_end - TDatabase::ParamDB->time_vanka_solve_start; 

//================================================= End Solving local system ===========================================

  
  // update dof
	      for(j=0;j<N_U;j++)
	      {
		l = UDOFs[j];
		u1[l] += damp*Rhs[j];
		u2[l] += damp*Rhs[j+N_U];
		#ifdef __3D__
		u3[l] += damp*Rhs[j+N_U2];
		#endif  		
	      }
	      p[i] += damp*Rhs[N_UGEO];
		//	----------------------------===================================-----------------------=====================================-------------------	      
	    } // endfor loop over pressure nodes

	    
	    // apply damping
	    for(j=0;j<N_DOF;j++)
	      u1[j] = uold[j]+alpha*(u1[j]-uold[j]);

	    // set Dirichlet values
	    memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
		  N_Dirichlet*SizeOfDouble);
	    memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
		  N_Dirichlet*SizeOfDouble);
	  #ifdef __3D__
	    memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
		  N_Dirichlet*SizeOfDouble);
	  #endif
	    
	      // itmethod exists
	    if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
	    {
	      TDatabase::ParamDB->INTERNAL_LOCAL_DOF=0;
	      delete itmethod;
	    }
	    
  if(TDatabase::ParamDB->DOF_Reorder)
  {  
	  u2 = u1 + N_UDOF;
	#ifdef __3D__
	  u3 = u2 + N_UDOF;
	#endif
	  p  = u1 + GEO_DIM*N_UDOF;

	  rhs2 = rhs1 + N_UDOF;
	#ifdef __3D__
	  rhs3 = rhs2 + N_UDOF;
	#endif
	  rhsp = rhs1 + GEO_DIM*N_UDOF;

	  // set Dirichlet values
	  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
		N_Dirichlet*SizeOfDouble);
	  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
		N_Dirichlet*SizeOfDouble);
	#ifdef __3D__
	  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
		N_Dirichlet*SizeOfDouble);
	#endif

	  // old values
	  uold = aux;
	  pold = uold+GEO_DIM*N_UDOF;

	  // save current solution on 'old' vectors
	  memcpy(uold, u1, N_DOF*SizeOfDouble);
		
	  /////////////////////////////////////////////////////////=====================================================
	#ifdef _MPI 
	    for(order=0;order<N_PDOF;order++)
	#else
	    for(i=0;i<N_PDOF;i++)
	#endif
	    {
	#ifdef _MPI
	      if(TDatabase::ParamDB->DOF_Reorder == 1)
		i = re_pdof[order];    
	      else
		i = order;	    
	  #endif			
		N_U = 0;

	  #ifdef _MPI    
	      
	      if(TDatabase::ParamDB->DOF_Average){
		if(DofmarkerP[i] == 'h' || DofmarkerP[i] == 'H')
		  continue;
	      }     
	      else{
		if(DofmarkerP[i] == 'h' || DofmarkerP[i] == 'H'  ||  DofmarkerP[i] == 's')
		  continue;
	      }
	      
	  #endif
	      // go through row i of B1 and B2
	      begin = BRowPtr[i];
	      end = BRowPtr[i+1];
	      value = rhsp[i];      // rhs of this pressure value
	    for(k=begin;k<end;k++)
	      { 
		l=BKCol[k]; 
		value1 = B1Entries[k];
		value2 = B2Entries[k];
	  #ifdef __3D__
		value3 = B3Entries[k];
	  #endif
	  #ifdef __2D__
		value -= value1*u1[l]+value2*u2[l]; // update rhs
	  #endif      
	  #ifdef __3D__
		value -= value1*u1[l]+value2*u2[l]+value3*u3[l]; // update rhs
	  #endif    

	  #ifdef _MPI   
	    if(TDatabase::ParamDB->DOF_Average != 3){
	      if(TDatabase::ParamDB->DOF_Average){
		if(DofmarkerU[l] == 'h' || DofmarkerU[l] == 'H')
		  continue;
	      }
	      
	      else{
		if(DofmarkerU[l] == 'h' || DofmarkerU[l] == 'H' || DofmarkerU[l] == 's')
		  continue;
	      }
	    }
	  #endif

	      UDOFs[N_U] = l;
	      j1 = GEO_DIM*N_U;
	      System[j1] = value1;  // save values for local B
	      System[j1+1] = value2;
	  #ifdef __3D__
	      System[j1+2] = value3;
	  #endif  
	      N_U++;           // count # velo dof connected to the pressure dof
	      }                    // row done

	      N_U2 = 2 * N_U;
	      N_UGEO = GEO_DIM * N_U;
	      N_LocalDOF = N_UGEO +1;

	      if (N_LocalDOF > LargestDirectSolve){
		if (N_LocalDOF != TDatabase::ParamDB->INTERNAL_LOCAL_DOF){
		  // itmethod exists already
		  if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0){
		    delete itmethod;
		  }
		  itmethod = new TFgmresIte(MatVect, Defect, NULL, 0, N_LocalDOF, 1);
		  TDatabase::ParamDB->INTERNAL_LOCAL_DOF = N_LocalDOF;
		}
	      }
	      
	      memset(System+N_UGEO, 0, SizeOfDouble*(N_LocalDOF*N_LocalDOF-N_UGEO));

	      Rhs[N_LocalDOF-1] = value;  // set rhs

	      for (k=0;k<N_U;k++)         // copy local B to the right place
	      {
		j4 = GEO_DIM*k;
		System[k*N_LocalDOF+N_UGEO]=System[j4];
		System[(k+N_U)*N_LocalDOF+N_UGEO]=System[j4+1];
	  #ifdef __3D__
		System[(k+N_U2)*N_LocalDOF+N_UGEO]=System[j4+2];
	  #endif
	      }
	      memset(System, 0, SizeOfDouble*N_UGEO);


	      for(j=0;j<N_U;j++){
		
		j1 = j;
		j2 = j+N_U;
		#ifdef __3D__
		j3 = j2+N_U;
		#endif
		
		UDOF = UDOFs[j];

		// A block
		begin = ARowPtr[UDOF];
		end = ARowPtr[UDOF+1];

		Rhs[j1] = rhs1[UDOF];
		Rhs[j2] = rhs2[UDOF];
		#ifdef __3D__
		Rhs[j3] = rhs3[UDOF];
		#endif

		for(k=begin;k<end;k++)
		{
		  l = AKCol[k];
		  value11 = A11Entries[k];
		  value12 = A12Entries[k];
		  value21 = A21Entries[k];
		  value22 = A22Entries[k];
		  #ifdef __3D__
		  value13 = A13Entries[k];
		  value23 = A23Entries[k];
		  value31 = A31Entries[k];
		  value32 = A32Entries[k];
		  value33 = A33Entries[k];
		  #endif

		  #ifdef __2D__
		  if (UDOF>=HangingBound) // Dirichlet node
		    value21 = value12 = 0;

		  Rhs[j1] -= value11*u1[l]+value12*u2[l];
		  Rhs[j2] -= value21*u1[l]+value22*u2[l];
		  #endif
		  
		  #ifdef __3D__
		  if (UDOF>=HangingBound) // Dirichlet node
		    value12 = value13 = value21 = value23 = value31 = value32 = 0;
		  
		  Rhs[j1] -= value11*u1[l]+value12*u2[l]+value13*u3[l];
		  Rhs[j2] -= value21*u1[l]+value22*u2[l]+value23*u3[l];
		  Rhs[j3] -= value31*u1[l]+value32*u2[l]+value33*u3[l];

		  #endif

		  for(m=0;m<N_U;m++)
		    if(UDOFs[m]==l)
		    {
		      k1 = m*N_LocalDOF;
		      k2 = (m+N_U)*N_LocalDOF;

		      System[k1+j1] = value11;
		      System[k2+j1] = value12;
		      System[k1+j2] = value21;
		      System[k2+j2] = value22;
		      #ifdef __3D__
		      k3 = (m + 2*N_U)*N_LocalDOF;
		      System[k3+j1] = value13;
		      System[k3+j2] = value23;
		      System[k1+j3] = value31;
		      System[k2+j3] = value32;
		      System[k3+j3] = value33;
		      #endif
		      break;
		    }
		} // endfor k

		if(UDOF<HangingBound)  // active dof
		{
		  // transpose(B) block for non-Dirichlet nodes
		  begin = BTRowPtr[UDOF];
		  end = BTRowPtr[UDOF+1];

		  for(k=begin;k<end;k++)
		  {
		    l = BTKCol[k];
		    value1 = B1TEntries[k];
		    value2 = B2TEntries[k];
		    #ifdef __3D__
		    value3 = B3TEntries[k];
		    #endif
		    value = p[l];
		 
		    {
		    Rhs[j1] -= value1*value;
		    Rhs[j2] -= value2*value;
		    #ifdef __3D__
		    Rhs[j3] -= value3*value;
		    #endif
		    }

		    if(i==l)
		    {
		      j4 = N_UGEO*N_LocalDOF;
		      System[j4+j1] = value1;
		      System[j4+j2] = value2;
		      #ifdef __3D__
		      System[j4+j3] = value3;
		      #endif
		    }
		  } // endfor k
		} // endif UDOF<HangingBound
	      } // endfor j

	      if(C)
	      {
		// fill C block if present
		begin = CRowPtr[i];
		end = CRowPtr[i+1];
		for(k=begin;k<end;k++)
		{
		  l = CKCol[k];
		  value = -CEntries[k]; // minus is right sign
		  Rhs[N_LocalDOF-1] -= value*p[l];
		  if(l==i) // main diagonal
		    System[N_LocalDOF*N_LocalDOF-1] = value;
		} // endfor k
	      } // endif C
	      
//=================================================Start Solving local system ===========================================

#ifdef _MPI
  TDatabase::ParamDB->time_vanka_solve_start = MPI_Wtime();
#else
  TDatabase::ParamDB->time_vanka_solve_start = GetTime();
#endif
	      // solve local system
	      if (smoother==3 && !C) // no diagonal Vanka for matrix C
	      {
	  #ifdef __2D__
		// diagonal Vanka
		SolveDiagonalVanka2D(System,  Rhs, N_U, 1, N_LocalDOF);
	  #endif
	  #ifdef __3D__
		// diagonal Vanka
		SolveDiagonalVanka3D(System,  Rhs, N_U, 1, N_LocalDOF);
	  #endif
	      }
	      else
	      {
		// full Vanka
		if (N_LocalDOF > LargestDirectSolve)
		{
		  int iter =0;
		  
		  double check2 = 0.0;
		  for(iter=0;iter<N_LocalDOF;iter++)
		    check2 += Rhs[iter]*Rhs[iter];
		  
		  memset(sol,0,N_LocalDOF*SizeOfDouble);
		  verbose =  TDatabase::ParamDB->SC_VERBOSE;
		  TDatabase::ParamDB->SC_VERBOSE = -1;
		  
		  int tem = TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER;
		  TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER = 0;
		  
		  if(check2)
		  itmethod->Iterate(matrix,NULL,sol,Rhs);
			  
		  TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER=tem;
		  TDatabase::ParamDB->SC_VERBOSE = verbose;
		  memcpy(Rhs, sol, N_LocalDOF*SizeOfDouble);
		}
		else
		{           
		  SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
		}
	      }
	      
	      	      
#ifdef _MPI
  TDatabase::ParamDB->time_vanka_solve_end = MPI_Wtime();
#else
  TDatabase::ParamDB->time_vanka_solve_end = GetTime();
#endif	      
	      
  TDatabase::ParamDB->time_vanka_solve += TDatabase::ParamDB->time_vanka_solve_end - TDatabase::ParamDB->time_vanka_solve_start; 

//================================================= End Solving local system ===========================================


	      // update dof
	      for(j=0;j<N_U;j++)
	      {
		l = UDOFs[j];
		u1[l] += damp*Rhs[j];
		u2[l] += damp*Rhs[j+N_U];
		#ifdef __3D__
		u3[l] += damp*Rhs[j+N_U2];
		#endif  		
	      }
	      p[i] += damp*Rhs[N_UGEO];

	  } // endfor loop over pressure nodes

		      // apply damping
		      for(j=0;j<N_DOF;j++)
			u1[j] = uold[j]+alpha*(u1[j]-uold[j]);

		      // set Dirichlet values
		      memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
			    N_Dirichlet*SizeOfDouble);
		      memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
			    N_Dirichlet*SizeOfDouble);
		    #ifdef __3D__
		      memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
			    N_Dirichlet*SizeOfDouble);
		    #endif

		      // itmethod exists
		      if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
		      {
			delete itmethod;
		      }       
      }
      
#ifdef _MPI      
   ParCommU->CommUpdate(u1);   
   ParCommP->CommUpdate(p);
#endif
   
#ifdef _MPI
  TDatabase::ParamDB->time_vanka_end = MPI_Wtime();
#else
  TDatabase::ParamDB->time_vanka_end = GetTime();
#endif

  TDatabase::ParamDB->time_vanka += TDatabase::ParamDB->time_vanka_end - TDatabase::ParamDB->time_vanka_start; 
} // end Vanka

/** solve exact on this level */
void TNSE_MGLevel4::SolveExact(double *u1, double *rhs1)
{
  double *a, *b;
  int i,j,k,l,index, N_DOF2 = N_DOF*N_DOF, end, begin, m;
  int N_UDOF2 = 2*N_UDOF, N_UDOFGEO = GEO_DIM *N_UDOF ;
  double value, value1, value2, value3;

  a = new double[N_DOF2];
  b = new double[N_DOF];

  memset(a, 0, N_DOF2*SizeOfDouble);

  j = ARowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      index = AKCol[j];
      a[index*N_DOF+i] =  A11Entries[j];
      a[(index+N_UDOF)*N_DOF+i] = A12Entries[j];
      a[index*N_DOF+i+N_UDOF] = A21Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOF] = A22Entries[j];
#ifdef __3D__
      a[(index+N_UDOF2)*N_DOF+i] = A13Entries[j];
      a[(index+N_UDOF2)*N_DOF+i+N_UDOF] = A23Entries[j];
      a[index*N_DOF+i+N_UDOF2] = A31Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOF2] = A32Entries[j];
      a[(index+N_UDOF2)*N_DOF+i+N_UDOF2] = A33Entries[j];
#endif
    }
  } // endfor i

  // Dirichlet and hanging nodes
  for(i=N_Active;i<N_UDOF;i++)
  {
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      a[i+N_DOF*index] = A11Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOF] = A22Entries[j];
#ifdef __3D__
      a[(index+N_UDOF2)*N_DOF+i+N_UDOF2] = A33Entries[j];
#endif
    }
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      l = (index+GEO_DIM*N_UDOF)*N_DOF;
      a[l+i] = B1TEntries[j];
      a[l+i+N_UDOF] = B2TEntries[j];
#ifdef __3D__
      a[l+i+2*N_UDOF] = B3TEntries[j];
#endif
    } // endfor j
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      a[(index)*N_DOF+i+N_UDOFGEO] = B1Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOFGEO] = B2Entries[j];
#ifdef __3D__
      a[(index+N_UDOF2)*N_DOF+i+N_UDOFGEO] = B3Entries[j];
#endif
    } // endfor j
  } // endfor i

  if(C)
  {
    j = CRowPtr[0];
    for(i=0;i<N_PDOF;i++)
    {
      k = CRowPtr[i+1];
      for(;j<k;j++)
      {
        index = CKCol[j];
        value = -CEntries[j]; // minus is right sign
        a[(index+GEO_DIM*N_UDOF)*N_DOF + i+GEO_DIM*N_UDOF] = value;
      }
    } // endfor i
  } // endif C

  // condition for pressure, fix first value
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    for (i=0;i<N_DOF;i++)
      a[i*N_DOF+N_UDOFGEO] = 0;
    a[(N_UDOFGEO)*N_DOF+N_UDOFGEO] = 1;
  }

  // copy into local data
  memcpy(b, rhs1, N_DOF*SizeOfDouble);
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    b[N_UDOFGEO] = 0;

  /* for (i=0;i<N_DOF;i++)
      for (j=0;j<N_DOF;j++)
      OutPut("a("<<i+1<<","<<j+1<<") = " << a[i*N_DOF+j] << endl);*/
 SolveLinearSystemLapack(a, b, N_DOF, N_DOF);
	  
  // copy from local data
  memcpy(u1, b, N_DOF*SizeOfDouble);

  delete a;
  delete b;
}
/** solve exact on this level */
void TNSE_MGLevel4::SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag)
{
 
#ifdef __3D__
    if (umfpack_flag==-1)
    {
	DirectSolver(A11, A12, A13, A21, A22, A23, A31, A32, A33,
		     B1T, B2T, B3T, B1, B2, B3, rhs1, u1, 4);
	umfpack_flag = 0;
    }
    DirectSolver(A11, A12, A13, A21, A22, A23, A31, A32, A33,
		 B1T, B2T, B3T, B1, B2, B3, rhs1, u1, umfpack_flag);
    umfpack_flag = 1;
#endif
#ifdef __2D__

    // free memory from previous time
    if (umfpack_flag==-1)
    {
	DirectSolver(A11, A12, A21, A22, B1T, B2T, B1, B2, rhs1, u1, 4);
	umfpack_flag = 0;
    }
	
    //OutPut("TNSE_MGLevel4::SolveExactUMFPACK: Are we here?" << endl);
    DirectSolver(A11, A12, A21, A22, B1T, B2T, B1, B2, rhs1, u1, umfpack_flag);
    umfpack_flag = 1;
#endif
}

/** step length control for smoother */
double TNSE_MGLevel4::StepLengthControl (double *u1, double *u1old,
                                         double *def1,
                                         int N_Parameters,
                                         double *Parameters)
{
  //==================================================
#ifdef _MPI
  int ii, rank, *MasterOfDofU,*MasterOfDofP, dof;
  double res_global=0.0;
  double res1 = 0.0;
  
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
  MasterOfDofU = ParCommU->GetMaster();
  MasterOfDofP = ParCommP->GetMaster();

#endif
  
  //====================================================
  double *x,*y,omega,numerator,nominator;
  int i,j;

  // allocate auxiliary array
  x = new double[2*N_DOF];
  y = x+N_DOF;

  // save current update in array x
  for (i=0;i<N_DOF;i++)
    x[i] = u1[i]-u1old[i];
  memset(y,0,N_DOF*SizeOfDouble);

  // compute matrix times update
#ifdef __2D__
//  if(C)
//    CoupledMatVect(A11,A12,A21,A22,B1,B2,B1T,B2T,C,x,y);
//  else
    CoupledMatVect(A11,A12,A21,A22,B1,B2,B1T,B2T,x,y);
#endif
#ifdef __3D__
  // if(C)
  //   CoupledMatVect(A11,A12,A13,A21,A22,A23,A31,A32,A33,
  //                  B1,B2,B3,B1T,B2T,B3T,C,x,y);
  // else
  CoupledMatVect(A11,A12,A13,A21,A22,A23,A31,A32,A33,
                 B1,B2,B3,B1T,B2T,B3T,x,y);
  #ifdef _MPI
    ParCommU->CommUpdate(y);
    ParCommP->CommUpdate(y + GEO_DIM*N_UDOF);
  #endif
  
#endif

#ifdef _MPI
  
  res1=0.0;
  
  for(i=0; i<N_UDOF; i++)
    if(MasterOfDofU[i] == rank)
    {      
      res1 += def1[i	       ]*y[i  	       ];      
      res1 += def1[i +   N_UDOF]*y[i +   N_UDOF];
#ifdef __3D__
      res1 += def1[i + 2*N_UDOF]*y[i + 2*N_UDOF];
#endif    
    }

  for(i=0;i<N_PDOF;i++)
    if(MasterOfDofP[i] == rank)
#ifdef __3D__
      res1 += def1[i + 3*N_UDOF]*y[i + 3*N_UDOF];
#endif
#ifdef __2D__
      res1 += def1[i + 2*N_UDOF]*y[i + 2*N_UDOF];
#endif
  
  res_global = 0.0;
  
  TDatabase::ParamDB->time_communication_start = MPI_Wtime();

  MPI_Allreduce(&res1, &res_global, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  
  TDatabase::ParamDB->time_communication_end =  MPI_Wtime();
  TDatabase::ParamDB->time_communication += TDatabase::ParamDB->time_communication_end - TDatabase::ParamDB->time_communication_start;
  
  
  numerator = (res_global);
#else
  numerator = Ddot(N_DOF,def1,y);
  
#endif
  
  
#ifdef _MPI  
  res1 = 0.0;  
  for(i=0; i<N_UDOF; i++)
    if(MasterOfDofU[i] == rank){
      
      res1 += y[i           ]*y[i           ];      
      res1 += y[i +   N_UDOF]*y[i +   N_UDOF];
#ifdef __3D__
      res1 += y[i + 2*N_UDOF]*y[i + 2*N_UDOF];
#endif    
  }

  for(i=0;i<N_PDOF;i++)
    if(MasterOfDofP[i] == rank)
#ifdef __3D__
      res1 += y[i + 3*N_UDOF]*y[i + 3*N_UDOF];
#endif
#ifdef __2D__
      res1 += y[i + 2*N_UDOF]*y[i + 2*N_UDOF];
#endif
      
  res_global = 0.0;
  
  TDatabase::ParamDB->time_communication_start = MPI_Wtime();
  
  MPI_Allreduce(&res1, &res_global, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);  
  
  TDatabase::ParamDB->time_communication_end =  MPI_Wtime();
  TDatabase::ParamDB->time_communication += TDatabase::ParamDB->time_communication_end - TDatabase::ParamDB->time_communication_start;
  
  
  
  nominator = (res_global);
  
#else  
  nominator = Ddot(N_DOF,y,y);
#endif
  
  if (nominator > 0)
    omega = numerator/nominator;
  else
    {
      if(N_Parameters>0)
        omega = alpha;
      else
        omega = 0.5;
      OutPut("MESSAGE : Step length control failed. Set omega = " << omega<< endl);
    }
  if (fabs(omega)<0.0001)
    {
       if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.9;
    }
  delete x;
  if (TDatabase::ParamDB->SC_VERBOSE>=2)
    OutPut("step length control " << omega << endl);
  return(omega);
  
}


/** Braess--Sarazin smoother  */
void TNSE_MGLevel4::BraessSarazin(double *u1, double *rhs1,
                                  double *aux, int N_Parameters,
                                  double *Parameters,int N_Levels)
{
  double *sol;
  int j;

  // array for update
  sol = aux;

  // initialize array for update
  memset(sol, 0, (N_DOF)*SizeOfDouble);

  // set pointer to rhs and solution
  j = 0;

  // call the algebraic solver
  // the last input parameter is only a dummy
#ifdef __2D__
  if(C)
  {
    OutPut("Braess-Sarazin smoother not implemented for C block !!!" << endl);
    exit(-1);
  }
  else
    Solver(A11, A12, A21, A22, B1, B2, B1T, B2T, rhs1, sol, j);
#endif
#ifdef __3D__
  OutPut("Braess-Sarazin smoother not implemented !!!" << endl);
  exit(4711);
#endif

  // update and apply damping
  for(j=0;j<N_DOF;j++)
    u1[j] += alpha*sol[j];

}
/** print all matrices and oth right hand sides */
void TNSE_MGLevel4::PrintAll()
{
}

#ifdef _MPI
void TNSE_MGLevel4::Par_Directsolve(double *sol, double * rhs)
{

  TSquareMatrix3D** SQMATRICES;
  SQMATRICES = new TSquareMatrix3D*[9];
  
  TMatrix3D** MATRICES;
  MATRICES = new TMatrix3D*[6];
  
    SQMATRICES[0] = A11;
    SQMATRICES[1] = A12;
    SQMATRICES[2] = A13;	  
    SQMATRICES[3] = A21;
    SQMATRICES[4] = A22;
    SQMATRICES[5] = A23; 
    SQMATRICES[6] = A31;
    SQMATRICES[7] = A32;
    SQMATRICES[8] = A33;  
    MATRICES[0] = B1;
    MATRICES[1] = B2;
    MATRICES[2] = B3;
    MATRICES[3] = B1T;
    MATRICES[4] = B2T;
    MATRICES[5] = B3T;

  int temp = TDatabase::ParamDB->SC_VERBOSE ;
  TDatabase::ParamDB->SC_VERBOSE =24 ;
  TParDirectSolver* DS = new TParDirectSolver(ParCommU,ParCommP,SQMATRICES,MATRICES);
  DS->Solve(sol,rhs,true);
  
  TDatabase::ParamDB->SC_VERBOSE =temp;
  
    SQMATRICES[0] = NULL;
    SQMATRICES[1] = NULL;
    SQMATRICES[2] = NULL;	  
    SQMATRICES[3] = NULL;
    SQMATRICES[4] = NULL;
    SQMATRICES[5] = NULL; 
    SQMATRICES[6] = NULL;
    SQMATRICES[7] = NULL;
    SQMATRICES[8] = NULL;  
    MATRICES[0] = NULL;
    MATRICES[1] = NULL;
    MATRICES[2] = NULL;
    MATRICES[3] = NULL;
    MATRICES[4] = NULL;
    MATRICES[5] = NULL;

    delete SQMATRICES;
    delete MATRICES;
  //DS->~TParDirectSolver();

}
void TNSE_MGLevel4::gmres_solve(double *sol, double * rhs)
{

  TSquareMatrix** SQMATRICES;
  SQMATRICES = new TSquareMatrix*[9];
  int dof = N_DOF;
  TMatrix** MATRICES;
  MATRICES = new TMatrix*[6];
  
    SQMATRICES[0] = A11;
    SQMATRICES[1] = A12;
    SQMATRICES[2] = A13;	  
    SQMATRICES[3] = A21;
    SQMATRICES[4] = A22;
    SQMATRICES[5] = A23; 
    SQMATRICES[6] = A31;
    SQMATRICES[7] = A32;
    SQMATRICES[8] = A33;  
    MATRICES[0] = B1;
    MATRICES[1] = B2;
    MATRICES[2] = B3;
    MATRICES[3] = B1T;
    MATRICES[4] = B2T;
    MATRICES[5] = B3T;
int temp = 0;
temp =  TDatabase::ParamDB->SC_VERBOSE;
 TDatabase::ParamDB->SC_VERBOSE = -24;
TFgmresIte* test = new TFgmresIte(MatVect_NSE4,Defect_NSE4,NULL,0,dof,0,ParCommU,ParCommP);
test->Iterate(SQMATRICES,MATRICES,sol,rhs);
   TDatabase::ParamDB->SC_VERBOSE = temp;
    SQMATRICES[0] = NULL;
    SQMATRICES[1] = NULL;
    SQMATRICES[2] = NULL;	  
    SQMATRICES[3] = NULL;
    SQMATRICES[4] = NULL;
    SQMATRICES[5] = NULL; 
    SQMATRICES[6] = NULL;
    SQMATRICES[7] = NULL;
    SQMATRICES[8] = NULL;  
    MATRICES[0] = NULL;
    MATRICES[1] = NULL;
    MATRICES[2] = NULL;
    MATRICES[3] = NULL;
    MATRICES[4] = NULL;
    MATRICES[5] = NULL;

    delete SQMATRICES;
    delete MATRICES;
  //DS->~TParDirectSolver();

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * 
 * 
 * 
 * 
 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TNSE_MGLevel4::UpdateHaloRhs(double * CurrentRhsU,double* CurrentDefectU){

  char *DofmarkerU = ParCommU->Get_DofMarker();  
  int i;
  
  for(i=0;i<N_UDOF;i++){
    
    if(DofmarkerU[i] == 'h' || DofmarkerU[i]=='H'){
      
      CurrentRhsU[i           ] = CurrentDefectU[i           ];
      CurrentRhsU[i +   N_UDOF] = CurrentDefectU[i +   N_UDOF];
      CurrentRhsU[i + 2*N_UDOF] = CurrentDefectU[i + 2*N_UDOF];
    }      
  }  
}

void TNSE_MGLevel4::ComputeOrder(){
    
  TCollection* coll = USpace->GetCollection();  
  TBaseCell * cell; 
  int iter ;

  if(TDatabase::ParamDB->DOF_Reorder == 1){
    
    cell_queue = createQ();
    u_clip = new int[N_UDOF];
    p_clip = new int[N_PDOF];
    cell_clip = new int[coll->GetN_Cells()];
      
    memset(u_clip,-1,sizeof(int)*N_UDOF);
    memset(p_clip,-1,sizeof(int)*N_PDOF);
    memset(cell_clip,-1,sizeof(int)*(coll->GetN_OwnCells()));
    memset(cell_clip + coll->GetN_OwnCells(),0,sizeof(int)*(coll->GetN_Cells() - coll->GetN_OwnCells()));        
    
    for(iter =coll->GetN_OwnCells();iter < coll->GetN_Cells(); iter++){
	  cell = coll->GetCell(iter);
	  Update(cell);
	  NeighUpdate(cell,iter);
    }
    
    
    while((iter = dequeue(cell_queue)) > -1)
      NeighUpdate(coll->GetCell(iter),iter);    
    delete cell_queue;
    
    Reorder();
  }
  else
  {    
     re_pdof = new int[N_PDOF];     
     for(iter=0;iter < N_PDOF; iter++)
       re_pdof[iter] = iter;
  }
}

void TNSE_MGLevel4::Update(TBaseCell* cell){

  int id = cell->GetCellIndex();
  
  int *u_dof,*p_dof;
  int min = cell_clip[id]; // if clip != min, then min+1, must match cell_clip ?? 
  int * begin;
  int *Glob;
  int i,temp;
  int count;
  
  begin = USpace->GetBeginIndex();
  Glob = USpace->GetGlobalNumbers();
  
  u_dof = Glob + begin[id];
  count = begin[id+1] - begin[id];
  for(i = 0;i < count; i++ ){
    if(cell->IsHaloCell()){
      u_clip[u_dof[i]] = 0;
    }
    else{
      temp = u_clip[u_dof[i]];
      if(temp != -1 && temp < min)
	min = temp;
    }
  }
    
  begin = PSpace->GetBeginIndex();
  Glob = PSpace->GetGlobalNumbers();
  
  p_dof = Glob + begin[id];
  count = begin[id+1] - begin[id];
  for(i = 0;i < count; i++ ){
    if(cell->IsHaloCell()){
      p_clip[p_dof[i]] = 0;
    }
    else{
      temp = p_clip[p_dof[i]];
      if(temp != -1 && temp < min)
	min = temp;
    }
  }
  
  
  begin = USpace->GetBeginIndex();
  Glob = USpace->GetGlobalNumbers();
  
  u_dof = Glob + begin[id];
  count = begin[id+1] - begin[id];
  for(i = 0;i < count; i++ ){
    temp = u_clip[u_dof[i]];
    if(temp == -1 || temp > min+1) // temp cannot be more check this::
      u_clip[u_dof[i]] = min+1;
  }
    
  begin = PSpace->GetBeginIndex();
  Glob = PSpace->GetGlobalNumbers();
  
  p_dof = Glob + begin[id];
  count = begin[id+1] - begin[id];
  for(i = 0;i < count; i++ ){
    temp = p_clip[p_dof[i]];
    if(temp == -1 || temp > min+1)
      p_clip[p_dof[i]] = min+1;
  }

}

void TNSE_MGLevel4::NeighUpdate(TBaseCell* cell, int iter){

  int n_joints;
  int n_edges;
  TJoint* joint;
  TBaseCell* neigh;
  TBaseCell** neighs;
  TEdge* edge;
  int e_neigh;
  
  int i,j,id;
  
      n_joints = cell->GetN_Joints();
      n_edges = cell->GetN_Edges();
      
    for(i=0;i<n_joints;i++){
      
      joint = cell->GetJoint(i);
      neigh = joint->GetNeighbour(cell);
      
      if(neigh){
	id = neigh->GetCellIndex();	
	if(cell_clip[id] == -1 ){
	  cell_clip[id] = cell_clip[iter]+1;
	  Update(neigh);
	  enqueue(id,cell_queue);
	}
      }
    }    
     
    for(i=0;i<n_edges;i++){
      edge = cell->GetEdge(i);
      edge->GetNeibs(e_neigh, neighs);
      if(neighs){
	for(j=0;j<e_neigh;j++){
	  neigh = neighs[j];
	  id = neigh->GetCellIndex();
	  if(cell_clip[id] == -1 ){
	    cell_clip[id] = cell_clip[iter]+1;
	    Update(neigh);
	    enqueue(id,cell_queue);
	  }	
	}
      }  
    }
}

void TNSE_MGLevel4::Reorder()
{

  int n=0,i;
  int *l;
  int *off;
  int type;
  int loc;

  for(i=0;i<N_PDOF;i++)
    if(n < p_clip[i])     
	n = p_clip[i];
  
  l = new int[n+1];
  off = new int[n+1];

  re_pdof = new int[N_PDOF];
  memset(l,0,sizeof(int)*(n+1));
  memset(off,0,sizeof(int)*(n+1));
  
  for(i=0;i<N_PDOF;i++)
    l[p_clip[i]+1]++;    
    
  for(i=1;i<n+1;i++){
    l[i] += l[i-1];
  }
    
  
  for(i=0;i<N_PDOF;i++){
    type = p_clip[i];
    loc = l[type]+off[type];
    if(loc < N_PDOF)
      re_pdof[loc] = i;
    else 
    {
      cout <<"This shouldnt happen" << endl;
      exit(0);
      
    }
    off[type]++;
  }
        
}



#endif
