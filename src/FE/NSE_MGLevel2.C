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
// @(#)NSE_MGLevel2.C        1.10 07/03/00
//
// Class:       TNSE_MGLevel2
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system
//              type 2 (A, B1, B2, B1T, B2T)
//
// Author:      Gunar Matthies 24.08.1999
// Author:      Volker John 24.08.1999
//
// History:     24.08.1999 start of implementation
//              24.08.1999 CellVanka 
//              24.08.1999 step length control
//              24.08.1999 ExactSolve
//              26.08.1999 NodalVanka
//
// =======================================================================

#include <NSE_MGLevel2.h>
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

#include <LocalProjection.h>

/** constructor */
#ifdef __2D__
TNSE_MGLevel2::TNSE_MGLevel2(int level, TSquareMatrix2D *a, 
                             TMatrix2D *b1, TMatrix2D *b2, 
                             TMatrix2D *b1t, TMatrix2D *b2t,
                             double *f1, double *u1,
                             int n_aux, double *al, int velocity_space, 
                             int pressure_space, TCollection *Coll,
                             int *dw)
#endif    
#ifdef __3D__
TNSE_MGLevel2::TNSE_MGLevel2(int level, TSquareMatrix3D *a, 
                             TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,
                             TMatrix3D *b1t, TMatrix3D *b2t, TMatrix3D *b3t,
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
  
  A = a;
  StructureA = A->GetMatrixStructure();
  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  B1T = b1t;
  B2T = b2t;
  StructureBT = B1T->GetStructure();
  BTRowPtr = StructureBT->GetRowPtr();
  BTKCol = StructureBT->GetKCol();
  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
#ifdef __3D__
  B3T = b3t;
  B3TEntries = B3T->GetEntries();
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

  C = NULL;
  StructureC = NULL;
  CRowPtr = NULL;
  CKCol = NULL;
  CEntries = NULL;

  USpace = A->GetFESpace();
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

  Type = 2;

  alpha = al[0];
  downwind = dw;
}

/** constructor with matrix C */
#ifdef __2D__
TNSE_MGLevel2::TNSE_MGLevel2(int level, TSquareMatrix2D *a, 
                             TMatrix2D *b1, TMatrix2D *b2, 
                             TMatrix2D *b1t, TMatrix2D *b2t,
                             TMatrix2D *c,
                             double *f1, double *u1,
                             int n_aux, double *al, int velocity_space, 
                             int pressure_space, TCollection *Coll,
                             int *dw)
#endif    
#ifdef __3D__
TNSE_MGLevel2::TNSE_MGLevel2(int level, TSquareMatrix3D *a, 
                             TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,
                             TMatrix3D *b1t, TMatrix3D *b2t, TMatrix3D *b3t,
                             TMatrix3D *c,
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
  
  A = a;
  StructureA = A->GetMatrixStructure();
  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  B1T = b1t;
  B2T = b2t;
  StructureBT = B1T->GetStructure();
  BTRowPtr = StructureBT->GetRowPtr();
  BTKCol = StructureBT->GetKCol();
  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
#ifdef __3D__
  B3T = b3t;
  B3TEntries = B3T->GetEntries();
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

  C = c;
  StructureC = C->GetStructure();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  USpace = A->GetFESpace();
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

  Type = 2;

  alpha = al[0];
  downwind = dw;
}

/** destructor */
TNSE_MGLevel2::~TNSE_MGLevel2()
{
} // ~TNSE_MGLevel2

/** calculate defect */
void TNSE_MGLevel2::Defect(double *u, double *f, double *d, double &res)
{
#ifdef __2D__
  // compute defect
  //if(C)
  //  CoupledDefect(A,B1,B2,B1T,B2T,C,u,f,d);
  //else
    CoupledDefect(A,B1,B2,B1T,B2T,u,f,d);

  // project d2 into L20
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(d+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  
#ifdef __3D__
  // compute defect
  // if(C)
  //   CoupledDefect(A,B1,B2,B3,B1T,B2T,B3T,C,u,f,d);
  // else
    CoupledDefect(A,B1,B2,B3,B1T,B2T,B3T,u,f,d);

  // project d3 into L20
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(d+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  

  memset(d+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
  memset(d+N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#ifdef __3D__
  memset(d+2*N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#endif  

  // compute residual
  res = sqrt(Ddot(N_DOF,d,d));
}

/** correct Dirichlet and hanging nodes */
void TNSE_MGLevel2::CorrectNodes(double *u1)
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
        s -= AEntries[j] * u1[index];
        t -= AEntries[j] * u2[index];
      }
    } // endfor j
    u1[i] = s;
    u2[i] = t;
  } // endfor i
#endif
}

/** cellwise Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel2::CellVanka(double *u1, double *rhs1, double *aux,
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
  int i,j,k,l,m, N_Cells, ii;
  int j1, j2, j3, k1, k2, k3;
  double value, value1, value2, value3;
  double *uold, *pold; 
  TCollection *Coll;
  double System[RhsDim*RhsDim];
  double Rhs[RhsDim], sol[RhsDim];
  int *UGlobalNumbers, *UBeginIndex, *UDOFs, UDOF, N_U;
  int *PGlobalNumbers, *PBeginIndex, *PDOFs, PDOF, N_P;
  int N_LocalDOF, verbose;
  int begin, end, HangingBound;
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
  HangingBound = USpace->GetHangingBound();

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

    UDOFs = UGlobalNumbers+UBeginIndex[ii]; // global dofs
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
        value = AEntries[k];

        Rhs[j1] -= value*u1[l];
        Rhs[j2] -= value*u2[l];
#ifdef __3D__
        Rhs[j3] -= value*u3[l];
#endif  

        for(m=0;m<N_U;m++)
          if(UDOFs[m]==l)
          {
            // column belongs to local system
            k1 = m;
            k2 = m+N_U;
#ifdef __3D__
            k3 = k2+N_U;
#endif  
            System[k1*N_LocalDOF+j1] = value;
            System[k2*N_LocalDOF+j2] = value;
#ifdef __3D__
            System[k3*N_LocalDOF+j3] = value;
#endif  
            break;
          }

      } // endfor k

      if(UDOF<N_Active)  // active dof
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
      } // endif UDOF<HangingBound
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
#ifdef __2D__
        Rhs[j1] -= value1*u1[l] + value2*u2[l];
#endif  
#ifdef __3D__
        value3 = B3Entries[k];
        Rhs[j1] -= value1*u1[l] + value2*u2[l] + value3 *u3[l];
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

    /*
    for(j=0;j<N_LocalDOF;j++)
    {
      for(k=0;k<N_LocalDOF;k++)
        cout << setw(6) << j << setw(6) << k << setw(20) << System[k*N_LocalDOF+j] << endl;
    }
    */

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
    
    /*
    for(j=0;j<N_LocalDOF;j++)
    {
      for(k=0;k<N_LocalDOF;k++)
        cout << setw(6) << j << setw(6) << k << setw(20) << System[k*N_LocalDOF+j] << endl;
    }
    */

    // solve local system
    if (smoother==1 && !C) // no diagonal Vanka for matrix C
    {
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
      {
        SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
        /*
        for(j=0;j<N_LocalDOF;j++)
          cout << setw(6) << j << setw(23) << Rhs[j] << endl;
        */
      }
    }

    // update dof
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

  if(TDatabase::ParamDB->INTERNAL_LOCAL_DOF > 0)
  {
    delete itmethod;
  }

  SetHangingNodes(u1);
} // end Vanka

/** nodal Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel2::NodalVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother, int N_Levels)
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
  double *uold, *pold; 
  double System[SystemRhs*SystemRhs];
  double Rhs[SystemRhs], sol[SystemRhs];
  int N_LocalDOF, verbose;
  int begin, end, HangingBound;
  double *u2, *u3, *p, *rhs2, *rhs3, *rhsp;
  int UDOFs[MaxN_LocalU], UDOF, N_U, N_U2, N_UGEO;
  TItMethod *itmethod = NULL, *prec;
  int LargestDirectSolve = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
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

  HangingBound = USpace->GetHangingBound();

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
    N_U2 = 2 * N_U; 
    N_UGEO = GEO_DIM * N_U; 
    N_LocalDOF = N_UGEO +1;

    if (N_LocalDOF > LargestDirectSolve)
    {
      // size of local system has changed
      if (N_LocalDOF != TDatabase::ParamDB->INTERNAL_LOCAL_DOF)
      {
        // itmethod exists already
        if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
          delete itmethod;

        // allocate new itmethod
        prec = NULL;
        itmethod = new TFgmresIte(MatVect, Defect, prec, 0, N_LocalDOF, 1);
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
        value = AEntries[k];

        Rhs[j1] -= value*u1[l];
        Rhs[j2] -= value*u2[l];
#ifdef __3D__
        Rhs[j3] -= value*u3[l];
#endif  

        for(m=0;m<N_U;m++)
          if(UDOFs[m]==l)
          {
            // column belongs to local system
            k1 = m;
            k2 = m+N_U;
#ifdef __3D__
            k3 = k2+N_U;
#endif   
            System[k1*N_LocalDOF+j1] = value;
            System[k2*N_LocalDOF+j2] = value;
#ifdef __3D__
            System[k3*N_LocalDOF+j3] = value;
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
          Rhs[j1] -= value1*value;
          Rhs[j2] -= value2*value;
#ifdef __3D__
          Rhs[j3] -= value3*value;
#endif  

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
void TNSE_MGLevel2::SolveExact(double *u1, double *rhs1)
{
  double *a, *b;
  int i,j,k,index, N_DOF2 = N_DOF*N_DOF;
  double value;

  a = new double[N_DOF2];
  b = new double[N_DOF];

  memset(a, 0, N_DOF2*SizeOfDouble);

  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      a[index*N_DOF+i] = value;
      a[(index+N_UDOF)*N_DOF+i+N_UDOF] = value;
#ifdef __3D__
      a[(index+2*N_UDOF)*N_DOF+i+2*N_UDOF] = value;
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
      a[(index+GEO_DIM*N_UDOF)*N_DOF+i] = B1TEntries[j];
      a[(index+GEO_DIM*N_UDOF)*N_DOF+i+N_UDOF] = B2TEntries[j];
#ifdef __3D__
      a[(index+GEO_DIM*N_UDOF)*N_DOF+i+2*N_UDOF] = B3TEntries[j];
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
      a[(i+GEO_DIM*N_UDOF)+N_DOF*index] = B1Entries[j];
      a[(i+GEO_DIM*N_UDOF)+N_DOF*(index+N_UDOF)] = B2Entries[j];
#ifdef __3D__
      a[(i+GEO_DIM*N_UDOF)+N_DOF*(index+2*N_UDOF)] = B3Entries[j];
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
        value = -CEntries[j]; // minux is right sign
        a[(index+GEO_DIM*N_UDOF)*N_DOF + i+GEO_DIM*N_UDOF] = value;
      }
    } // endfor i
  } // endif C

  // condition for pressure, fix first value
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  { 
    for (i=0;i<N_DOF;i++)
      a[i*N_DOF+GEO_DIM*N_UDOF] = 0;
    a[(GEO_DIM*N_UDOF)*N_DOF+GEO_DIM*N_UDOF] = 1;
  }

  // copy into local data
  memcpy(b, rhs1, N_DOF*SizeOfDouble);

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    b[GEO_DIM*N_UDOF] = 0;

  // matrix is stored column-wise !!
  SolveLinearSystemLapack(a, b, N_DOF, N_DOF);

  // copy from local data
  memcpy(u1, b, N_DOF*SizeOfDouble);

  delete a;
  delete b;
}

/** solve exact on this level */
void TNSE_MGLevel2::SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag)
{

  OutPut("TNSE_MGLevel2::SolveExactUMFPACK: Are we here?" << endl);
  OutPut("Not yet implemented!!!" << endl);
  exit(4711);
}

/** step length control for Vanka */
double TNSE_MGLevel2::StepLengthControl (double *u1, double *u1old, 
         double *def1,
         int N_Parameters, double *Parameters)
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

#ifdef __2D__
  // compute matrix times update
  //if(C)
  //  CoupledMatVect(A,B1,B2,B1T,B2T,C,x,y);
  //else
    CoupledMatVect(A,B1,B2,B1T,B2T,x,y);
#endif  
#ifdef __3D__
  // if(C)
  //   CoupledMatVect(A,B1,B2,B3,B1T,B2T,B3T,C,x,y);
  // else
    CoupledMatVect(A,B1,B2,B3,B1T,B2T,B3T,x,y);
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
      OutPut("MESSAGE : Step length control failed. Set omega = " << omega << endl);
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
  {
    OutPut("step length control " << omega << endl);
  }
  return(omega);
}

/** print all matrices and all right hand sides */
void TNSE_MGLevel2::PrintAll()
{
  int i,j,k,index;
  double value, value1, value2;

  cout << endl;
  cout << "A block" << endl;
  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      cout << setw(5) << i << setw(5) << index << setw(20) << value << endl;
    }
  } // endfor i

  cout << endl;
  cout << "B1T block" << endl;
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      cout << setw(5) << i << setw(5) << index << setw(20) << value1 << endl;
    } // endfor j
  } // endfor i

  cout << endl;
  cout << "B2T block" << endl;
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value2 = B2TEntries[j];
      cout << setw(5) << i << setw(5) << index << setw(20) << value2 << endl;
    } // endfor j
  } // endfor i

  cout << endl;
  cout << "B1 block" << endl;
  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      cout << setw(5) << i << setw(5) << index << setw(20) << value1 << endl;
    } // endfor j
  } // endfor i

  cout << endl;
  cout << "B2 block" << endl;
  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value2 = B2Entries[j];
      cout << setw(5) << i << setw(5) << index << setw(20) << value2 << endl;
    } // endfor j
  } // endfor i

  if(C)
  {
    cout << endl;
    cout << "C block" << endl;
    j = CRowPtr[0];
    for(i=0;i<N_UDOF;i++)
    {
      k = CRowPtr[i+1];
      for(;j<k;j++)
      {
        index = CKCol[j];
        value = -CEntries[j]; // minus is right sign
        cout << setw(5) << i << setw(5) << index << setw(20) << value << endl;
      }
    } // endfor i
  } // endif C

  cout << endl;
  cout << "rhs(u1)" << endl;
  for(i=0;i<N_UDOF;i++)
    cout << Rhs1[i] << endl;

  cout << endl;
  cout << "rhs(u2)" << endl;
  for(i=0;i<N_UDOF;i++)
    cout << Rhs2[i] << endl;

  cout << endl;
  cout << "rhs(p)" << endl;
  for(i=0;i<N_PDOF;i++)
    cout << RhsP[i] << endl;
}

/** Braess--Sarazin smoother  */
void TNSE_MGLevel2::BraessSarazin(double *u1, double *rhs,
    double *aux, int N_Parameters, double *Parameters,int N_Levels)
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
    Solver(A, B1T, B2T, B1, B2, rhs, sol,j);
#endif  
#ifdef __3D__
  OutPut("Braess-Sarazin smoother not implemented !!!" << endl);
  exit(4711);
#endif  

  // update and apply damping
  for(j=0;j<N_DOF;j++)
    u1[j] += alpha*sol[j];
}
