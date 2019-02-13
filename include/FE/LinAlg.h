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
// @(#)LinAlg.h        1.8 12/07/99
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          27.01.1999
//
//              Sashikumaar Ganesan     08.10.2009 (eigen values)
// =======================================================================

#ifdef _MPI
#  include "mpi.h"
// #include <ParVector3D.h>
#endif

#ifndef __LINALG__
#define __LINALG__

#include <SquareMatrix2D.h>
#include <Matrix2D.h>

#ifdef __3D__
  #include <SquareMatrix3D.h>
  #include <Matrix3D.h>
#endif

// =======================================================================
// BLAS 1
// =======================================================================

/** return inner product (x,y) */
double Ddot(int n, double *x, double *y);

 #ifdef _MPI
/** return inner product (x,y) */
// double ParDdot(double *x, double *y, TParFECommunicator2D *ParCommunicator);

/** return inner product (x,y) */
void ParDdotNSE2D(double *x, int N_U, int N_P, TParFECommunicator2D *ParCommunicatorU, TParFECommunicator2D *ParCommunicatorP,
                  double &res, double &imp_res);

/** return inner product (x,y) */
double ParDdotNSE3D(double *x, double *y, int N_U, int N_P, TParFECommunicator2D *ParCommunicatorU, TParFECommunicator2D *ParCommunicatorP);
 #endif

/** y := alpha*x + y */
void Daxpy(int n, double alpha, double *x, double *y);

/** z := alpha*x + beta*y */
void Dsum(int n, double alpha, double beta, double *x, double *y, double *z);

/** b := a */
void Dcopy(int n, double *x, double *y);

/** x := alpha*x */
void Dscal(int n, double alpha, double *x);

/** x := alpha+x */
void Dadd(int n, double alpha, double *x);

/** return Euclidian norm of x */
double Dnorm(int n, double *x);

// =======================================================================
// BLAS 2 
// =======================================================================

/** y := A * x (for ONE a block) */
void MatVect(TSquareMatrix *A, double *x, double *y);
void MatVect_Scalar(TSquareMatrix **A, TMatrix **B, double *x, double *y);
/** y := b - A *x */
void ScalarDefect(TSquareMatrix *A, double *sol, double *f, double *d,
                  double &res);

void Defect_Scalar(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** y := A * x (for ONE a block), active rows only */
void MatVectActive(TSquareMatrix *A, double *x, double *y);

/** C := C + alpha*D, active rows only */
void MatAdd(TSquareMatrix *C, TSquareMatrix *D, double alpha);
void MatAdd(TMatrix *C, TMatrix *D, double alpha);

/** C := alpha*C + D , active rows only */
void MatAdd2(TSquareMatrix *C, TSquareMatrix *D, double alpha);

/** y := B * x */
void MatVect(TMatrix *A, double *x, double *y);

/** y := B * x (tmatrix * vector) */
void MatVect1(TMatrix *A, double *x, double *y);

/** y := B' * x (transposed matrix * vector) */
void TransMatVect(TMatrix *A, double *x, double *y);

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *y);

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *b, double *r);

void MatVectFull(TSquareMatrix **A, TMatrix **B, double *x, double *y);
void DefectFull(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);


void MatVect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T,
        double *x, double *y);

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T,
        double *x, double *b, double *r);

void MatVect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

#ifdef _MPI
void ParDefect_NSE2(TSquareMatrix **, TMatrix **, TParVectorNSE3D  *, TParVectorNSE3D *, TParVectorNSE3D *);
#endif

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12,
                    TSquareMatrix *A21, TSquareMatrix *A22,
                    TMatrix *B1, TMatrix *B2,
                    double *x, double *y);

/** defect for coupled Stokes / Navier-Stokes system NSTYPE==3*/
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12,
                   TSquareMatrix *A21, TSquareMatrix *A22, 
                   TMatrix *B1, TMatrix *B2,
                   double *x, double *b, double *r);

void MatVect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** matrix * vector for coupled Stokes / Navier-Stokes system NSTYPE==4*/
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12,
                    TSquareMatrix *A21, TSquareMatrix *A22,
                    TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    double *x, double *y);

/** defect for coupled Stokes / Navier-Stokes system NSTYPE==4*/
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12,
                   TSquareMatrix *A21, TSquareMatrix *A22,
                   TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   double *x, double *b, double *r);

void MatVect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** matrix * vector for coupled Stokes / Navier-Stokes system NSTYPE==14*/
/** equal order interpolation */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12,
                    TSquareMatrix *A21, TSquareMatrix *A22,
		    TSquareMatrix *C,
                    TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    double *x, double *y);

/** defect for coupled Stokes / Navier-Stokes system NSTYPE==4*/
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12,
                   TSquareMatrix *A21, TSquareMatrix *A22,
		   TSquareMatrix *C,
                   TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   double *x, double *b, double *r);

void MatVect_EquOrd_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_EquOrd_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);


/**  Darcy using Raviart-Thomas (RT) or Brexxi-Douglas-Marini (DBM) elements
 ( A B' )
 ( B 0  )
*/
void CoupledMatVect(TSquareMatrix *A, TMatrix *B, double *x, double *y);
void CoupledDefect(TSquareMatrix *A, TMatrix *B, 
                   double *x, double *b, double *r);

// ========================================================================
// (Navier-)Stokes system with almost divergence-free basis, NSTYP==5
// ========================================================================

/** matrix * vector for coupled Stokes / Navier-Stokes system NSTYPE==5*/
void CoupledMatVect(TSquareMatrix *A, double *x, double *y);

/** defect for coupled Stokes / Navier-Stokes system NSTYPE==5*/
void CoupledDefect(TSquareMatrix *A, double *x, double *b, double *r);
                   
void MatVect_NSE5(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE5(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** Convection-diffusion problem with VMM */
void MatVectCD_VMM(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TSquareMatrix *D, 
                   double *x, double *y);
void MatVect_CD_VMM(TSquareMatrix **A, TMatrix **B, double *x, double *y);
void DefectCD_VMM(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TSquareMatrix *D,
                  double *x, double *b, double *r);
void Defect_CD_VMM(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** Convection-diffusion problem with VMM [KL02] */
void MatVectCD_VMM_KL02(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, 
                        TMatrix *C1, TMatrix *C2, TSquareMatrix *D, 
                        double *x, double *y);
void MatVect_CD_VMM_KL02(TSquareMatrix **A, TMatrix **B, double *x, double *y);
void DefectCD_VMM_KL02(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, 
                        TMatrix *C1, TMatrix *C2, TSquareMatrix *D, 
                        double *x, double *y, double *r);
void Defect_CD_VMM_KL02(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** matrix * vector for coupled Stokes / Navier-Stokes system NSTYPE==3*/
void CoupledMatVectLV96(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *y, double delta);

/** r := b - A * x */
void CoupledDefectLV96(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *b, double *r, double delta);

/** solve a system of linear equations */
void SolveLinearSystemLapack(double *a, double *b, int N_Eqn, int LDA);

/** solve a system of linear equations */
//void SolveLinearSystem(double *a, double *b, int N_Eqn, int LDA);

/** solve a system of linear equations with transposed matrix*/
void SolveLinearSystemTranspose(double *a, double *b, int N_Eqn, int LDA);

/* subroutine for solving a multiple systems of linear equations */
/* solution is transposed */
void SolveMultipleSystems(double *a, double *b, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs);

/* subroutine for solving a multiple systems of linear equations */
void SolveMultipleSystemsNew(double *a, double *b, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs);

/* subroutine for solving a multiple systems of linear equations 
   using the LAPACK routines*/

void SolveMultipleSystemsLapack(double *a, double *b, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs);

/* solve diagonal Vanka system */
void SolveDiagonalVanka2D(double *a, double *b, int N_U, int N_P, int LDA);

/* solve diagonal Vanka system */
void SolveDiagonalVanka3D(double *a, double *b, int N_U, int N_P, int LDA);

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVectMortar(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *matrix_mortar, 
        double *x, double *y);

/** calculate the eigenvalue of the system using Lapack routines*/
void FindEigenValues(double *ap, char &UPLO, int N_Eqn, char &COMPZ, double *d, double *z);


// =======================================================================
// Multi grid components
// =======================================================================

#ifdef __2D__
/** prolongate */
void Prolongate(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
        double *CoarseFunctions, double *FineFunctions, double *aux);

void Prolongate(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
        int N_Functions,
        double *CoarseFunctions, double *FineFunctions, double *aux);

/** defect restriction from level+1 to level */
void DefectRestriction(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
    double *CoarseFunctions, double *FineFunctions, double *aux);

/** defect restriction from level+1 to level */
void DefectRestriction(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
    int N_Functions,
    double *CoarseFunctions, double *FineFunctions, double *aux);

/** function restriction from level+1 to level */
void RestrictFunction(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
    double *CoarseFunction, double *FineFunction, double *aux);

/** function restriction from level+1 to level */
void RestrictFunction(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
    int N_Functions, 
    double *CoarseFunction, double *FineFunction, double *aux);

/** project vector v into L20 */
void IntoL20Vector2D(double *v, int Length, int order);

#ifdef _MPI
/** project parallel vector v into L20 */
void IntoL20Vector2D(double *v, int Length, int order, TParFECommunicator2D *ParCommunicator);
#endif

/** project fe function v into L20 */
void IntoL20FEFunction(double *v, int Length, TFESpace2D *FESpace,
                       int velocity_space, int pressure_space
                       #ifdef _MPI
                       , MPI_Comm comm
                       #endif
                      );

void VMS_ProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
                             TSquareMatrix2D **SQMATRICES, 
                             TMatrix2D **MATRICES);
#endif // __2D__

#ifdef __3D__

/** coupled matrix vector products */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
        double *x, double *y);

void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   double *x, double *b, double *r);

void MatVect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T, 
                    double *x, double *y);

void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r);

void MatVect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    double *x, double *y);

void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   double *x, double *b, double *r);

void MatVect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                    double *x, double *y);

void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r);

void MatVect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
		    TSquareMatrix *C,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                    double *x, double *y);

void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
		   TSquareMatrix *C,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r);

void MatVect_NSE14(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void Defect_NSE14(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

/** prolongate */
void Prolongate(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
        double *CoarseFunctions, double *FineFunctions,
        double *aux);

/** prolongate */
void Prolongate(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
        int N_Functions,
        double *CoarseFunctions, double *FineFunctions,
        double *aux);

/** defect restriction from level+1 to level */
void DefectRestriction(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
    double *CoarseFunctions, double *FineFunctions, double *aux);

/** defect restriction from level+1 to level */
void DefectRestriction(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
    int N_Functions,
    double *CoarseFunctions, double *FineFunctions, double *aux);

/** function restriction from level+1 to level */
void RestrictFunction(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
    double *CoarseFunction, double *FineFunction, double *aux);

/** function restriction from level+1 to level */
void RestrictFunction(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
    int N_Functions,
    double *CoarseFunction, double *FineFunction, double *aux);

/** project vector v into L20 */
void IntoL20Vector3D(double *v, int Length, int order);

/** project fe function into L20 */
void IntoL20FEFunction3D(double *v, int Length, TFESpace3D *FESpace);

void IntoL20FEFunction3D(double *v, int Length, TFESpace3D *FESpace,
                       int velocity_space, int pressure_space);

void VMS_ProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
                             TSquareMatrix3D **SQMATRICES, TMatrix3D **MATRICES);

void LumpMassMatrixToDiag(TSquareMatrix3D *M);
#endif // __3D__ 



#endif
