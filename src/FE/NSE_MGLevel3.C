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
// @(#)NSE_MGLevel3.C        1.10 07/03/00
//
// Class:       TNSE_MGLevel3
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system
//              type 2 (A, B1, B2, B1T, B2T)
//
// Author:      Gunar Matthies 24.08.1999
// Author:      Volker John 25.08.1999
//
// History:     24.08.1999 start of implementation
//              25.08.1999 CellVanka 
//              25.08.1999 step length control
//              25.08.1999 ExactSolve
//              26.08.1999 NodalVanka
//
// =======================================================================

#include <NSE_MGLevel3.h>
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
#endif  

#include <stdlib.h>
#include <string.h>

#include <LinAlg.h>
#include <Solver.h>
#include <ItMethod.h>
#include <FgmresIte.h>

/** constructor */
#ifdef __2D__
  TNSE_MGLevel3::TNSE_MGLevel3(int level, 
                              TSquareMatrix2D *a11, TSquareMatrix2D *a12,
                              TSquareMatrix2D *a21, TSquareMatrix2D *a22,
                              TMatrix2D *b1, TMatrix2D *b2,  
                              TStructure2D *structureBT,
                              double *f1, double *u1, 
                              int n_aux, double *al, int velocity_space, 
                              int pressure_space, TCollection *Coll,
			      int *dw)
#endif    
#ifdef __3D__
  TNSE_MGLevel3::TNSE_MGLevel3(int level, 
                              TSquareMatrix3D *a11, TSquareMatrix3D *a12,
                              TSquareMatrix3D *a13, TSquareMatrix3D *a21,
                              TSquareMatrix3D *a22, TSquareMatrix3D *a23,
                              TSquareMatrix3D *a31, TSquareMatrix3D *a32,
                              TSquareMatrix3D *a33,
                              TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,  
                              TStructure3D *structureBT,
                              double *f1, double *u1, 
                              int n_aux, double *al, int velocity_space, 
                              int pressure_space, TCollection *Coll,
			      int *dw)
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

  B1 = b1;
  B2 = b2;
  StructureB = B1->GetStructure();
  BRowPtr = StructureB->GetRowPtr();
  BKCol = StructureB->GetKCol();
  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
#ifdef __3D__
  B3 = b3;
  B3Entries = B3->GetEntries();
#endif  

  StructureBT = structureBT;
  BTRowPtr = StructureBT->GetRowPtr();
  BTKCol = StructureBT->GetKCol();

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

  Type = 3;

  alpha = al[0];
  downwind = dw;
}

/** destructor */
TNSE_MGLevel3::~TNSE_MGLevel3()
{
} // ~TNSE_MGLevel3

/** calculate defect */
void TNSE_MGLevel3::Defect(double *u1, double *f1,  double *d1, double &res)
{

#ifdef __2D__
  // compute defect
  CoupledDefect(A11,A12,A21,A22,B1,B2,u1,f1,d1);
  // project defect into L20
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(d1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  
#ifdef __3D__
  // compute defect
  CoupledDefect(A11,A12,A13,A21,A22,A23,A31,A32,A33,B1,B2,B3,u1,f1,d1);
  // project defect into L20
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(d1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  

  memset(d1+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
  memset(d1+N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#ifdef __3D__
  memset(d1+2*N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#endif  

  // compute residual
  res = sqrt(Ddot(N_DOF,d1,d1));

}

/** correct Dirichlet and hanging nodes */
void TNSE_MGLevel3::CorrectNodes(double *u1)
{
  int i,j,k, index;
  double s, t, u, *u2, *u3;
  
  u2 = u1+N_UDOF;
#ifdef __3D__
  u3 = u2+N_UDOF;
#endif  
  
  i = N_Dirichlet*SizeOfDouble;
  // set Dirichlet nodes
  memset(u1+HangingNodeBound, 0, i);
  memset(u2+HangingNodeBound, 0, i);
#ifdef __3D__
  memset(u3+HangingNodeBound, 0, i);
#endif  
  
#ifdef __2D__
  // set hanging nodes 
  j = ARowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      if(index != i)
      {
        s -= A11Entries[j] * u1[index];
        t -= A22Entries[j] * u2[index];
      }
    } // endfor j
    u1[i] = s;
    u2[i] = t;
  } // endfor i
#endif
}

/** cellwise Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel3::CellVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother, int N_Levels)
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
  int i,j,k,l,m, N_Cells, N_LocalDOFs, ii, verbose;
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
  int N_LocalDOF;
  int begin, end, ActiveBound, begin1, end1;
  TBaseCell *Cell;
  double *u2, *u3, *p, *rhs2, *rhs3, *rhsp;
  TItMethod *itmethod = NULL;
  double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
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
	"TNSE_MGLevel3::CellVanka - Not enough memory in array Rhs!!!" 
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
          
          begin1 = BRowPtr[l];
          end1 =  BRowPtr[l+1];
          for (m=begin1;m<end1;m++)
            if (BKCol[m] == UDOF)
            {
              value1 = B1Entries[m];
              value2 = B2Entries[m];
#ifdef __3D__
              value3 = B3Entries[m];
#endif  
              break;
            }  
          value = p[l];
          Rhs[j1] -= value1*value;
          Rhs[j2] -= value2*value;
#ifdef __3D__
          Rhs[j3] -= value3*value;
#endif  
          
        } // endfor k
      } // endif UDOF<HangingBound
    } // endfor j

    for(j=0;j<N_P;j++)
    {
      j1 = j+GEO_DIM*N_U;       // local index
      PDOF = PDOFs[j];    // global index
      begin = BRowPtr[PDOF]; // PDOF-th row in B1, B2
      end = BRowPtr[PDOF+1];
      Rhs[j1] = rhsp[PDOF];  // set local rhs

      for(k=begin;k<end;k++) // go through row
      {
        l=BKCol[k];            // column in B, row in BT
        value1 = B1Entries[k]; // entries 
        value2 = B2Entries[k];
#ifdef __3D__
        value3 = B3Entries[k];
#endif  
        Rhs[j1] -= value1*u1[l]; // update rhs
        Rhs[j1] -= value2*u2[l];
#ifdef __3D__
        Rhs[j1] -= value3*u3[l];
#endif  
         
        for(m=0;m<N_U;m++)      // through local velo dofs
          if(UDOFs[m]==l)       // global velo dof coincides with l 
          {
            // column belongs to local system
            k1 = m;                 // local index
            k2 = m+N_U;
#ifdef __3D__
            k3 = k2 + N_U;
#endif  
            System[k1*N_LocalDOF+j1] = value1; // fill j1-th row
            System[k2*N_LocalDOF+j1] = value2;
#ifdef __3D__
            System[k3*N_LocalDOF+j1] = value3;
#endif  
            if (l<ActiveBound)                // fill j1-th column
            {
	      j4 = j1*N_LocalDOF;
              System[j4+k1] = value1; //transposed block
              System[j4+k2] = value2;
#ifdef __3D__
              System[j4+k3] = value3;
#endif  
            }
            break;
          }
      } // endfor k
    } // endfor j

   // solve local system
    if (smoother==1)       
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
         SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
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
void TNSE_MGLevel3::NodalVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother,int N_Levels)
{
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
  int LargestDirectSolve = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  MatVecProc *MatVect=MatVectFull;
  DefectProc *Defect=DefectFull;
  TSquareMatrix **matrix= (TSquareMatrix **)sqmatrix;  

  HangingBound = USpace->GetHangingBound();
  TDatabase::ParamDB->INTERNAL_LOCAL_DOF = -1;
#ifdef __2D__
  sqmatrix[0] = (TSquareMatrix2D *)System;
#endif  
#ifdef __3D__
  sqmatrix[0] = (TSquareMatrix3D *)System;
#endif  

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
  
  // start of loop over pressure dofs
  for(i=0;i<N_PDOF;i++)
  {
    N_U = 0;
    // go through row i of B1 and B2
    begin = BRowPtr[i];
    end = BRowPtr[i+1];
    value = rhsp[i];      // rhs of this pressure value
    for(k=begin;k<end;k++)
    {
      l=BKCol[k];
      UDOFs[N_U] = l;
      value1 = B1Entries[k];
      value2 = B2Entries[k];
#ifdef __3D__
      value3 = B3Entries[k];
#endif  
      j1 = GEO_DIM*N_U;
      System[j1] = value1;  // save values for local B
      System[j1+1] = value2;
#ifdef __2D__
      value -= value1*u1[l]+value2*u2[l]; // update rhs
#endif  
#ifdef __3D__
      System[j1+2] = value3;
      value -= value1*u1[l]+value2*u2[l]+value3*u3[l]; // update rhs
#endif  
      N_U++;             // count # velo dof connected to the pressure dof
    }                    // row done
    /*if (N_U>=MaxN_LocalU)
    {
      OutPut("TNSE_MGLevel3::NodalVanka - N_U too large !!! " << N_U 
	     << " " << MaxN_LocalU << __FILE__ << endl);
      exit(4711);
      }*/
    N_U2 = 2 * N_U; 
    N_UGEO = GEO_DIM * N_U; 
    N_LocalDOF = N_UGEO +1;
    /*if (N_LocalDOF > SystemRhs)
    {
      OutPut(
	     "TNSE_MGLevel3::NodalVanka - Not enough memory in array Rhs !!!" 
	     << endl << "available " << SystemRhs << " needed " <<
	     N_LocalDOF << endl);
      exit(4711);
      }*/
    if (N_LocalDOF > LargestDirectSolve)
    {
      // size of local system has changed
      if (N_LocalDOF != TDatabase::ParamDB->INTERNAL_LOCAL_DOF)
      {
        // itmethod exists already
        if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
        {
          //if ( TDatabase::ParamDB->P7 >0.0)
          //  delete prec;
          delete itmethod;
        }
        // allocate new itmethod
        //if ( TDatabase::ParamDB->P7 >0.0)
        //  prec = new TDiagVankaIte(MatVect, Defect, NULL, 0, N_LocalDOF);
        //else
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
    // reset first part of array System
    memset(System, 0, SizeOfDouble*N_UGEO);
        
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
        if (UDOF>=HangingBound) // Dirichlet node
          value12 = value21 = 0;

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

      if(UDOF<HangingBound)  // active dof
      {
       // transpose(B) block for non-Dirichlet nodes
        begin = BTRowPtr[UDOF];
        end = BTRowPtr[UDOF+1];

        for(k=begin;k<end;k++)
        {
          l = BTKCol[k];
          
          begin1 = BRowPtr[l];
          end1 =  BRowPtr[l+1];
          for (m=begin1;m<end1;m++)
            if (BKCol[m] == UDOF)
            {
              value1 = B1Entries[m];
              value2 = B2Entries[m];
#ifdef __3D__
              value3 = B3Entries[m];
#endif  
              break;
            }  
          value = p[l];
          Rhs[j1] -= value1*value;
          Rhs[j2] -= value2*value;
#ifdef __3D__
          Rhs[j3] -= value3*value;
#endif  
          if(i==l)
          {
            // column belongs to local system
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
  
    // solve local system
    if (smoother==3)       
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
        memset(sol,0,N_LocalDOF*SizeOfDouble);
        //OutPut(i << " local dof " << N_LocalDOF << endl);
         //for (ii=0;ii<10;ii++)
        //  cout << ii << " " << System[ii] << endl;
        verbose =  TDatabase::ParamDB->SC_VERBOSE; 
        TDatabase::ParamDB->SC_VERBOSE = -1;
        itmethod->Iterate(matrix,NULL,sol,Rhs);
        TDatabase::ParamDB->SC_VERBOSE = verbose;
        //SolveLinearSystem(System, Rhs, N_LocalDOF, N_LocalDOF);
        //norm =0 ;
        //for (ii=0;ii<N_LocalDOF;ii++)
        //  norm+=(sol[i]-Rhs[i])*(sol[i]-Rhs[i]);
        //OutPut(" difference " << sqrt(norm) << endl);       
        memcpy(Rhs, sol, N_LocalDOF*SizeOfDouble);
      }
      else
         SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
    }
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
} // end Vanka

/** solve exact on this level */
void TNSE_MGLevel3::SolveExact(double *u1, double *rhs1)
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
      begin = BRowPtr[index];
      end = BRowPtr[index+1];
      for (m=begin;m<end;m++)
        if (BKCol[m] == i)
        {
          value1 = B1Entries[m];
          value2 = B2Entries[m];
#ifdef __3D__
          value3 = B3Entries[m];
#endif  
          break;
        }  
      a[(index+N_UDOFGEO)*N_DOF+i] = value1;
      a[(index+N_UDOFGEO)*N_DOF+i+N_UDOF] = value2;
#ifdef __3D__
      a[(index+N_UDOFGEO)*N_DOF+i+N_UDOF2] = value3;
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

  SolveLinearSystemLapack(a, b, N_DOF, N_DOF);

  // copy from local data
  memcpy(u1, b, N_DOF*SizeOfDouble);

  delete a;
  delete b;
}

/** solve exact on this level */
void TNSE_MGLevel3::SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag)
{

  OutPut("TNSE_MGLevel3::SolveExactUMFPACK: Are we here?" << endl);
  OutPut("Not yet implemented!!!" << endl);
  exit(4711);
}

/** step length control for smoother */
double TNSE_MGLevel3::StepLengthControl (double *u1, double *u1old, 
                                         double *def1, 
                                         int N_Parameters, 
                                         double *Parameters)
{
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
  CoupledMatVect(A11,A12,A21,A22,B1,B2,x,y);
#endif  
#ifdef __3D__
  CoupledMatVect(A11,A12,A13,A21,A22,A23,A31,A32,A33,B1,B2,B3,x,y);
#endif  

  numerator = Ddot(N_DOF,def1,y);
  nominator = Ddot(N_DOF,y,y);

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
        omega = alpha;
      else
        omega = 0.9;
    }     
  delete x;
  if (TDatabase::ParamDB->SC_VERBOSE>=2)
    OutPut("step length control " << omega << endl);
  return(omega); 
}


/** Braess--Sarazin smoother  */
void TNSE_MGLevel3::BraessSarazin(double *u1, double *rhs1,
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
  Solver(A11, A12, A21, A22, B1, B2, rhs1, sol, j);
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
void TNSE_MGLevel3::PrintAll()
{
}
