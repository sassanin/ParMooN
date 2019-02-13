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
   
// ======================================================================
// @(#)MovingNavierStokes.h        1.1 10/19/99
//
// Assembling routine for all matrices and right-hand sides
// ======================================================================

#ifndef __MOVINGNAVIERSTOKES__
#define __MOVINGNAVIERSTOKES__

// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one A block, 
//      M block from time discretization
//      B1, B2 (divergence blocks)
// ======================================================================

int MovingNSType1N_Terms = 4;
MultiIndex2D MovingNSType1Derivatives[4] = { D10, D01, D00, D00 };
int MovingNSType1SpaceNumbers[4] = { 0, 0, 0, 1 };
int MovingNSType1N_Matrices = 4;
int MovingNSType1RowSpace[4] = { 0, 0, 1, 1 };
int MovingNSType1ColumnSpace[4] = { 0, 0, 0, 0 };
int MovingNSType1N_Rhs = 2;
int MovingNSType1RhsSpace[2] = { 0, 0 };

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void MovingNSType1Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void MovingNSType1Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      M block from time discretization
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int MovingNSType2N_Terms = 4;
MultiIndex2D MovingNSType2Derivatives[4] = { D10, D01, D00, D00 };
int MovingNSType2SpaceNumbers[4] = { 0, 0, 0, 1 };
int MovingNSType2N_Matrices = 6;
int MovingNSType2RowSpace[6] = { 0, 0, 1, 1, 0, 0 };
int MovingNSType2ColumnSpace[6] = { 0, 0, 0, 0, 1, 1 };
int MovingNSType2N_Rhs = 2;
int MovingNSType2RhsSpace[2] = { 0, 0 };

// ======================================================================
// Type 2, Standard Galerkin
// ======================================================================
void MovingNSType2Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 2 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void MovingNSType2Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all four A blocks,
//      B1, B2 (divergence blocks), 
// ======================================================================

int MovingNSType3N_Terms = 4;
MultiIndex2D MovingNSType3Derivatives[4] = { D10, D01, D00, D00 };
int MovingNSType3SpaceNumbers[4] = { 0, 0, 0, 1 };
int MovingNSType3N_Matrices = 8;
int MovingNSType3RowSpace[8] = { 0, 0, 0, 0, 0, 0, 1, 1 };
int MovingNSType3ColumnSpace[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
int MovingNSType3N_Rhs = 2;
int MovingNSType3RhsSpace[2] = { 0, 0 };

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void MovingNSType3Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void MovingNSType3GalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void MovingNSType3Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void MovingNSType3UpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all four A blocks,
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int MovingNSType4N_Terms = 4;
MultiIndex2D MovingNSType4Derivatives[4] = { D10, D01, D00, D00 };
int MovingNSType4SpaceNumbers[4] = { 0, 0, 0, 1 };
int MovingNSType4N_Matrices = 10;
int MovingNSType4RowSpace[10] = { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
int MovingNSType4ColumnSpace[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 };
int MovingNSType4N_Rhs = 2;
int MovingNSType4RhsSpace[2] = { 0, 0 };

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void MovingNSType4Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void MovingNSType4GalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void MovingNSType4Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void MovingNSType4UpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Assembling routine for all nonlinear matrices
// ======================================================================

// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one nonlinear A block
//      WITHOUT right hand sides
// ======================================================================

int MovingNSType1NLN_Terms = 3;
MultiIndex2D MovingNSType1NLDerivatives[3] = { D10, D01, D00 };
int MovingNSType1NLSpaceNumbers[4] = { 0, 0, 0 };
int MovingNSType1NLN_Matrices = 1;
int MovingNSType1NLRowSpace[1] = { 0 };
int MovingNSType1NLColumnSpace[1] = { 0 };
int MovingNSType1NLN_Rhs = 0;
int *MovingNSType1NLRhsSpace = NULL;

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void MovingNSType1NLGalerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void MovingNSType1NLUpwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int MovingNSType2NLN_Terms = 3;
MultiIndex2D MovingNSType2NLDerivatives[3] = { D10, D01, D00 };
int MovingNSType2NLSpaceNumbers[3] = { 0, 0, 0 };
int MovingNSType2NLN_Matrices = 1;
int MovingNSType2NLRowSpace[1] = { 0 };
int MovingNSType2NLColumnSpace[1] = { 0 };
int MovingNSType2NLN_Rhs = 0;
int *MovingNSType2NLRhsSpace = NULL;

// ======================================================================
// Type 2, Standard Galerkin, only nonlinear part
// ======================================================================
void MovingNSType2NLGalerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void MovingNSType2NLUpwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      main diagonal blocks A11, A22
//      WITHOUT right hand sides
// ======================================================================

int MovingNSType3NLN_Terms = 3;
MultiIndex2D MovingNSType3NLDerivatives[3] = { D10, D01, D00 };
int MovingNSType3NLSpaceNumbers[3] = { 0, 0, 0 };
int MovingNSType3NLN_Matrices = 2;
int MovingNSType3NLRowSpace[2] = { 0, 0 };
int MovingNSType3NLColumnSpace[2] = { 0, 0 };
int MovingNSType3NLN_Rhs = 0;
int *MovingNSType3NLRhsSpace = NULL;

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void MovingNSType3NLGalerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void MovingNSType3NLGalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void MovingNSType3NLUpwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void MovingNSType3NLUpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      main diagonal blocks A11, A22
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int MovingNSType4NLN_Terms = 3;
MultiIndex2D MovingNSType4NLDerivatives[3] = { D10, D01, D00 };
int MovingNSType4NLSpaceNumbers[3] = { 0, 0, 0 };
int MovingNSType4NLN_Matrices = 2;
int MovingNSType4NLRowSpace[2] = { 0, 0 };
int MovingNSType4NLColumnSpace[2] = { 0, 0 };
int MovingNSType4NLN_Rhs = 0;
int *MovingNSType4NLRhsSpace = NULL;

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void MovingNSType4NLGalerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear part
// ======================================================================
void MovingNSType4NLGalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void MovingNSType4NLUpwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void MovingNSType4NLUpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Assembling routine for right-hand sides only
// ======================================================================

// ======================================================================
// declaration for all Navier-Stokes problems
//      ONLY right hand sides
// ======================================================================

int MovingNSRHSN_Terms = 1;
MultiIndex2D MovingNSRHSDerivatives[1] = { D00 };
int MovingNSRHSSpaceNumbers[1] = { 0 };
int MovingNSRHSN_Matrices = 0;
int *MovingNSRHSRowSpace = NULL;
int *MovingNSRHSColumnSpace = NULL;
int MovingNSRHSN_Rhs = 2;
int MovingNSRHSRhsSpace[2] = { 0, 0 };

// ======================================================================
// right-hand side ONLY
// ======================================================================
void MovingNSRHS(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);

// ======================================================================
// setting for error calculation for all types
// ======================================================================
MultiIndex2D MovingNSAllDerivatives[3] = { D00, D10, D01 };

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void MovingNSParams2(double *in, double *out);

int MovingNSN_FESpaces2 = 1;
int MovingNSN_Fct2 = 2;
int MovingNSN_ParamFct2 = 1;
int MovingNSN_FEValues2 = 2;
int MovingNSN_Params2 = 2;
int MovingNSFEFctIndex2[2] = { 0, 1 };
MultiIndex2D MovingNSFEMultiIndex2[2] = { D00, D00 };
ParamFct *MovingNSFct2[1] = { MovingNSParams2 };
int MovingNSBeginParam2[1] = { 0 };

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void GridParams4(double *in, double *out);

int GridN_FESpaces4 = 2;
int GridN_Fct4 = 4;
int GridN_ParamFct4 = 1;
int GridN_FEValues4 = 4;
int GridN_Params4 = 2;
int GridFEFctIndex4[4] = { 0, 1, 2, 3 };
MultiIndex2D GridFEMultiIndex4[4] = { D00, D00, D00, D00 };
ParamFct *GridFct4[1] = { GridParams4 };
int GridBeginParam4[1] = { 0 };

// ======================================================================
// declaration for grid moving matrix
//      one A block, 
// ======================================================================

int GridN_Terms = 2;
MultiIndex2D GridDerivatives[2] = { D10, D01 };
int GridSpaceNumbers[2] = { 0, 0 };
int GridN_Matrices = 4;
int GridRowSpace[4] = { 0,0,0,0 };
int GridColumnSpace[4] = { 0,0,0,0 };
int GridN_Rhs = 0;
int *GridRhsSpace = NULL;

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void GridAssemble(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// // kind of boundary condition (for FE space needed)
// void GridBoundCondition(int BdComp, double t, BoundCond &cond);
// 
// // value of boundary condition
// void GridBoundValue(int BdComp, double Param, double &value);
// 
// void GridCoeffs(int n_points, double *x, double *y,
//         double **parameters, double **coeffs);

#endif
