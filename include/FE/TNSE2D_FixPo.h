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
// @(#)TNSE2D_FixPo.h        1.2 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#ifndef __TNSE2D_FIXPO__
#define __TNSE2D_FIXPO__

// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one A block, 
//      M block from time discretization
//      B1, B2 (divergence blocks)
// ======================================================================

int TimeNSType1N_Terms = 4;
MultiIndex2D TimeNSType1Derivatives[4] = { D10, D01, D00, D00 };
int TimeNSType1SpaceNumbers[4] = { 0, 0, 0, 1 };
int TimeNSType1N_Matrices = 4;
int TimeNSType1RowSpace[4] = { 0, 0, 1, 1 };
int TimeNSType1ColumnSpace[4] = { 0, 0, 0, 0 };
int TimeNSType1N_Rhs = 2;
int TimeNSType1RhsSpace[2] = { 0, 0 };

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void TimeNSType1Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Layton 96
// ======================================================================
void TimeNSType1Layton96(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Smagorinsky
// ======================================================================
void TimeNSType1Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Coletti
// ======================================================================
void TimeNSType1Coletti(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, GL00Convolution
// ======================================================================
void TimeNSType1GL00Convolution(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, GL00AuxProblem
// ======================================================================
int TimeNSType1GL00AuxProblemN_Matrices = 5;
int TimeNSType1GL00AuxProblemRowSpace[5] = { 0, 0, 2, 1, 1 };
int TimeNSType1GL00AuxProblemColumnSpace[5] = { 0, 0, 2, 0, 0 };

void TimeNSType1GL00AuxProblem(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, group fem
// ======================================================================
void TimeNSType1GroupFEM(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs);


// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      M block from time discretization
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int TimeNSType2N_Terms = 4;
MultiIndex2D TimeNSType2Derivatives[4] = { D10, D01, D00, D00 };
int TimeNSType2SpaceNumbers[4] = { 0, 0, 0, 1 };
int TimeNSType2N_Matrices = 6;
int TimeNSType2RowSpace[6] = { 0, 0, 1, 1, 0, 0 };
int TimeNSType2ColumnSpace[6] = { 0, 0, 0, 0, 1, 1 };
int TimeNSType2N_Rhs = 2;
int TimeNSType2RhsSpace[2] = { 0, 0 };

// ======================================================================
// Type 2, Standard Galerkin
// ======================================================================
void TimeNSType2Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 2 SUPG
//      one A block, 
//      M block from time discretization
//      K block from supg discretization 
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================
int TimeNSType2N_TermsSUPG = 8;
MultiIndex2D TimeNSType2DerivativesSUPG[8] = { D10, D01, D00, D00 };
int TimeNSType2SpaceNumbersSUPG[8] = { 0, 0, 0, 1 };
int TimeNSType2N_MatricesSUPG = 7;
int TimeNSType2RowSpaceSUPG[7] = { 0, 0, 1, 1, 0, 0 };
int TimeNSType2ColumnSpaceSUPG[7] = { 0, 0, 0, 0, 1, 1 };
int TimeNSType2N_RhsSUPG = 2;
int TimeNSType2RhsSpaceSUPG[2] = { 0, 0 };

// ======================================================================
// Type 2, Standard SUPG
// ======================================================================
void TimeNSType2SUPG(double Mult, double *coeff, double *param, double hK, 
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);
// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, one K block from time discretization
//      B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================
int TimeNSType2NLN_TermsSUPG = 8;
MultiIndex2D TimeNSType2NLDerivativesSUPG[8] = { D10, D01, D00 };
int TimeNSType2NLSpaceNumbersSUPG[8] = { 0, 0, 0 };
int TimeNSType2NLN_MatricesSUPG = 4;
int TimeNSType2NLRowSpaceSUPG[4] = { 0 };
int TimeNSType2NLColumnSpaceSUPG[4] = { 0 };
int TimeNSType2NLN_RhsSUPG = 2;
int TimeNSType2NLRhsSpaceSUPG[2] = {0, 0};
// ======================================================================
// Type 2, Standard SUPG
// ======================================================================
void TimeNSType2NLSUPG(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs);
// ======================================================================


// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void TimeNSType2Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Layton96
// ======================================================================
void TimeNSType2Layton96(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Coletti
// ======================================================================
void TimeNSType2Coletti(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, GL00Convolution
// ======================================================================
void TimeNSType2GL00Convolution(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 2, GL00AuxProblem
// ======================================================================
int TimeNSType2GL00AuxProblemN_Matrices = 7;
int TimeNSType2GL00AuxProblemRowSpace[7] = { 0, 0, 2, 1, 1, 0, 0 };
int TimeNSType2GL00AuxProblemColumnSpace[7] = { 0, 0, 2, 0, 0, 1, 1 };

void TimeNSType2GL00AuxProblem(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all four A blocks,
//      B1, B2 (divergence blocks), 
// ======================================================================

int TimeNSType3N_Terms = 4;
MultiIndex2D TimeNSType3Derivatives[4] = { D10, D01, D00, D00 };
int TimeNSType3SpaceNumbers[4] = { 0, 0, 0, 1 };
int TimeNSType3N_Matrices = 8;
int TimeNSType3RowSpace[8] = { 0, 0, 0, 0, 0, 0, 1, 1 };
int TimeNSType3ColumnSpace[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
int TimeNSType3N_Rhs = 2;
int TimeNSType3RhsSpace[2] = { 0, 0 };

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void TimeNSType3Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3UpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Layton96, (grad u, grad v)
// ======================================================================
void TimeNSType3Layton96(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Layton96, D(u):D(v)
// ======================================================================
void TimeNSType3Layton96DD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Coletti, (grad u, grad v)
// ======================================================================
void TimeNSType3Coletti(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Coletti, D(u):D(v)
// ======================================================================
void TimeNSType3ColettiDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType3GL00Convolution(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType3GL00ConvolutionDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00AuxProblem, (grad u, grad v)
// ======================================================================
int TimeNSType3GL00AuxProblemN_Matrices = 9;
int TimeNSType3GL00AuxProblemRowSpace[9] =  { 0, 0, 0, 0, 0, 0, 2, 1, 1 };
int TimeNSType3GL00AuxProblemColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 2, 0, 0 };

void TimeNSType3GL00AuxProblem(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType3GL00AuxProblemDD(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, VMSProjection, D(u):D(v)
// ======================================================================
int TimeNSType3VMSProjectionN_Terms = 5;
MultiIndex2D TimeNSType3VMSProjectionDerivatives[5] = { D10, D01, D00, D00, D00 };
int TimeNSType3VMSProjectionSpaceNumbers[5] = { 0, 0, 0, 1, 3 };
int TimeNSType3VMSProjectionN_Matrices = 13;
int TimeNSType3VMSProjectionRowSpace[13] = { 0, 0, 0, 0, 0, 0, 
                                              3, 1, 1, 0, 0, 3, 3};
int TimeNSType3VMSProjectionColumnSpace[13] = { 0, 0, 0, 0, 0, 0, 
                                                 3, 0, 0, 3, 3, 0, 0};
                                           
void TimeNSType3VMSProjectionDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all four A blocks,
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

int TimeNSType4N_Terms = 4;
MultiIndex2D TimeNSType4Derivatives[4] = { D10, D01, D00, D00 };
int TimeNSType4SpaceNumbers[4] = { 0, 0, 0, 1 };
int TimeNSType4N_Matrices = 10;
int TimeNSType4RowSpace[10] = { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
int TimeNSType4ColumnSpace[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 };
int TimeNSType4N_Rhs = 2;
int TimeNSType4RhsSpace[2] = { 0, 0 };

// if convolution of velocity should be computed
//int TimeNSType4N_MatricesConvU = 11;
//int TimeNSType4RowSpaceConvU[11] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
//int TimeNSType4ColumnSpaceConvU[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 };

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void TimeNSType4Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void TimeNSType4UpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Layton96, (grad u, grad v)
// ======================================================================
void TimeNSType4Layton96(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Layton96, D(u):D(v)
// ======================================================================
void TimeNSType4Layton96DD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Coletti, (grad u, grad v)
// ======================================================================
void TimeNSType4Coletti(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Coletti, D(u):D(v)
// ======================================================================
void TimeNSType4ColettiDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4GL00Convolution(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GL00ConvolutionDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00AuxProblem, (grad u, grad v)
// ======================================================================
int TimeNSType4GL00AuxProblemN_Matrices = 11;
int TimeNSType4GL00AuxProblemRowSpace[11] =  
    { 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0 };
int TimeNSType4GL00AuxProblemColumnSpace[11] = 
    { 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1 };

void TimeNSType4GL00AuxProblem(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType4GL00AuxProblemDD(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 4, VMSProjection, D(u):D(v)
// ======================================================================
int TimeNSType4VMSProjectionN_Terms = 5;
MultiIndex2D TimeNSType4VMSProjectionDerivatives[5] = { D10, D01, D00, D00, D00 };
int TimeNSType4VMSProjectionSpaceNumbers[5] = { 0, 0, 0, 1, 3 };
int TimeNSType4VMSProjectionN_Matrices = 15;
int TimeNSType4VMSProjectionRowSpace[15] = { 0, 0, 0, 0, 0, 0, 
                                              3, 1, 1, 0, 0, 0, 0, 3, 3};
int TimeNSType4VMSProjectionColumnSpace[15] = { 0, 0, 0, 0, 0, 0, 
                                                 3, 0, 0, 1, 1, 3, 3, 0, 0};
                                           
void TimeNSType4VMSProjectionDD(double Mult, double *coeff, 
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

int TimeNSType1NLN_Terms = 3;
MultiIndex2D TimeNSType1NLDerivatives[3] = { D10, D01, D00 };
int TimeNSType1NLSpaceNumbers[3] = { 0, 0, 0 };
int TimeNSType1NLN_Matrices = 1;
int TimeNSType1NLRowSpace[1] = { 0 };
int TimeNSType1NLColumnSpace[1] = { 0 };
int TimeNSType1NLN_Rhs = 0;
int *TimeNSType1NLRhsSpace = NULL;

// ======================================================================
// declaration for all Navier-Stokes problems of type 2
//      one A block, 
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int TimeNSType2NLN_Terms = 3;
MultiIndex2D TimeNSType2NLDerivatives[3] = { D10, D01, D00 };
int TimeNSType2NLSpaceNumbers[3] = { 0, 0, 0 };
int TimeNSType2NLN_Matrices = 1;
int TimeNSType2NLRowSpace[1] = { 0 };
int TimeNSType2NLColumnSpace[1] = { 0 };
int TimeNSType2NLN_Rhs = 0;
int *TimeNSType2NLRhsSpace = NULL;

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// Type 2, Standard Galerkin, only nonlinear part
// Type 1, Coletti, only nonlinear part
// Type 2, Coletti, only nonlinear part
// Type 1, GL00Convolution, only nonlinear part
// Type 2, GL00Convolution, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// Type 2, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1_2NLUpwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Layto96, only nonlinear part
// Type 2, Layto96, only nonlinear part
// ======================================================================
void TimeNSType1_2NLLayton96(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      main diagonal blocks A11, A22
//      WITHOUT right hand sides
// ======================================================================

int TimeNSType3NLN_Terms = 3;
MultiIndex2D TimeNSType3NLDerivatives[3] = { D10, D01, D00 };
int TimeNSType3NLSpaceNumbers[3] = { 0, 0, 0 };
int TimeNSType3NLN_Matrices = 2;
int TimeNSType3NLRowSpace[2] = { 0, 0 };
int TimeNSType3NLColumnSpace[2] = { 0, 0 };
int TimeNSType3NLN_Rhs = 0;
int *TimeNSType3NLRhsSpace = NULL;

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      main diagonal blocks A11, A22
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int TimeNSType4NLN_Terms = 3;
MultiIndex2D TimeNSType4NLDerivatives[3] = { D10, D01, D00 };
int TimeNSType4NLSpaceNumbers[3] = { 0, 0, 0 };
int TimeNSType4NLN_Matrices = 2;
int TimeNSType4NLRowSpace[2] = { 0, 0 };
int TimeNSType4NLColumnSpace[2] = { 0, 0 };
int TimeNSType4NLN_Rhs = 0;
int *TimeNSType4NLRhsSpace = NULL;

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, Coletti, (grad u, grad v), only nonlinear part
// Type 4, Coletti, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 3, GL00Convolution, D(u):D(v), only nonlinear diagonal blocks
// Type 4, GL00Convolution, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// Type 4, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3_4NLUpwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// Type 4, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3_4NLUpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Layton96, (grad u, grad v), only nonlinear part
// Type 4, Standard Layton96, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLLayton96(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Layton96, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Layton96, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLLayton96DD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Standard Smagorinsky, (grad u, grad v), only nonlinear part
// ======================================================================
int TimeNSType3NLSmagorinskyN_Matrices = 4;
int TimeNSType3NLSmagorinskyRowSpace[4] = { 0, 0, 0, 0 };
int TimeNSType3NLSmagorinskyColumnSpace[4] = { 0, 0, 0, 0 };

void TimeNSType3_4NLSmagorinsky(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, VMSProjection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMSProjection, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
int TimeNSType3NLVMSProjectionN_Terms = 4;
MultiIndex2D TimeNSType3NLVMSProjectionDerivatives[4] = { D10, D01, D00, D00 };
int TimeNSType3NLVMSProjectionSpaceNumbers[4] = { 0, 0, 0, 3 };
int TimeNSType3NLVMSProjectionN_Matrices = 6;
int TimeNSType3NLVMSProjectionRowSpace[6] = { 0, 0, 0, 0, 0, 0};
int TimeNSType3NLVMSProjectionColumnSpace[6] = { 0, 0, 0, 0, 3, 3 };

void TimeNSType3_4NLVMSProjectionDD(double Mult, double *coeff, 
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

int TimeNSRHSN_Terms = 1;
MultiIndex2D TimeNSRHSDerivatives[1] = { D00 };
int TimeNSRHSSpaceNumbers[1] = { 0 };
int TimeNSRHSN_Matrices = 0;
int *TimeNSRHSRowSpace = NULL;
int *TimeNSRHSColumnSpace = NULL;
int TimeNSRHSN_Rhs = 2;
int TimeNSRHSRhsSpace[2] = { 0, 0 };

// ======================================================================
// right-hand side ONLY
// ======================================================================
void TimeNSRHS(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);
//=======================================================================
//=======================================================================
int TimeNSRHSN_TermsSUPG = 3;
MultiIndex2D TimeNSRHSDerivativesSUPG[3] = { D10, D01, D00 };
int TimeNSRHSSpaceNumbersSUPG[3] = { 0, 0, 0 };
int TimeNSRHSN_MatricesSUPG = 0;
int *TimeNSRHSRowSpaceSUPG = NULL;
int *TimeNSRHSColumnSpaceSUPG = NULL;
int TimeNSRHSN_RhsSUPG = 4;
int TimeNSRHSRhsSpaceSUPG[4] = { 0, 0, 0, 0 };

void TimeNSRHSSUPG(double Mult, double *coeff, double *param, double hK, 
                   double **OrigValues, int *N_BaseFuncts,
                   double ***LocMatrices, double **LocRhs);


void TimeNSRHSAuxProblemU(double Mult, double *coeff, 
                          double *param, double hK, 
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs);

int TimeNSRHSColN_Terms = 3;
MultiIndex2D TimeNSRHSColDerivatives[3] = { D10, D01, D00 };
int TimeNSRHSColSpaceNumbers[3] = { 0, 0, 0 };
int TimeNSRHSColN_Matrices = 0;
int *TimeNSRHSColRowSpace = NULL;
int *TimeNSRHSColColumnSpace = NULL;
int TimeNSRHSColN_Rhs = 2;
int TimeNSRHSColRhsSpace[2] = { 0, 0 };

// ======================================================================
// right-hand side ONLY, Coletti model
// ======================================================================
void TimeNSRHSColetti(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, Galdi-Layton model with convolution
// ======================================================================
void TimeNSRHSLESModel(double Mult, double *coeff, 
                   double *param, double hK, 
                   double **OrigValues, int *N_BaseFuncts,
                   double ***LocMatrices, double **LocRhs);


// ======================================================================
// right-hand side ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================

void TimeNSRHSGL00AuxProblemPaper2(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side for auxiliary problem 
// ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================

int TimeNSGL00AuxProblemRHSN_Terms = 1;
MultiIndex2D TimeNSGL00AuxProblemRHSDerivatives[1] = { D00 };
int TimeNSGL00AuxProblemRHSSpaceNumbers[1] = { 0 };
int TimeNSGL00AuxProblemRHSN_Matrices = 0;
int *TimeNSGL00AuxProblemRHSRowSpace = NULL;
int *TimeNSGL00AuxProblemRHSColumnSpace = NULL;
int TimeNSGL00AuxProblemRHSN_Rhs = 3;
int TimeNSGL00AuxProblemRHSRhsSpace[3] = { 0, 0, 0 };

void TimeNSGL00AuxProblemRHS(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);

void TimeNSGL00AuxProblemRHSPaper2(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, Smagorinsky Explicit
// ======================================================================
void TimeNSRHSSmagorinskyExplicit(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// + convection with higher order velocity
// ======================================================================
void TimeNSType3_4NLGalerkin_VMS_1_DD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side for additional terms in rhs of small scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_SmallRhs2D(double Mult, double *coeff, 
                           double *param, double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side for additional terms in rhs of large scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_Large_0_Rhs2D(double Mult, double *coeff, 
                              double *param, double hK, 
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side for additional terms in rhs of large scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_Large_1_Rhs2D(double Mult, double *coeff, 
                              double *param, double hK, 
                              double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

// ======================================================================
// Assembling routine for matrices for auxiliary problems
// ======================================================================

// ======================================================================
// declarations
// ======================================================================
/*
int TimeNSGL00AuxProblemN_Terms = 3;
MultiIndex2D TimeNSGL00AuxProblemDerivatives[3] = { D10, D01, D00 };
int TimeNSGL00AuxProblemSpaceNumbers[3] = { 0, 0, 0 };
int TimeNSGL00AuxProblemN_Matrices = 1;
int TimeNSGL00AuxProblemRowSpace[2] = { 0 };
int TimeNSGL00AuxProblemColumnSpace[2] = { 0 };
int TimeNSGL00AuxProblemN_Rhs = 0;
int *TimeNSGL00AuxProblemRhsSpace = NULL;

void TimeNSAuxMatrixGL00AuxProblem(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);
*/
// ======================================================================
// assemble matrix for auxiliary problem
// ======================================================================

int MatrixAuxiliaryProblemN_Terms = 3;
MultiIndex2D MatrixAuxiliaryProblemDerivatives[3] = { D10, D01, D00};
int MatrixAuxiliaryProblemSpaceNumbers[3] = { 0, 0, 0};
int MatrixAuxiliaryProblemN_Matrices = 1;
int MatrixAuxiliaryProblemRowSpace[1] = { 0 };
int MatrixAuxiliaryProblemColumnSpace[1] = { 0 };
int MatrixAuxiliaryProblemN_Rhs = 0;
int *MatrixAuxiliaryProblemRhsSpace = NULL;

void MatrixAuxiliaryProblem(double Mult, double *coeff, 
                            double *param, double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);


// ======================================================================
// declaration for VMS
// ONLY right hand sides
// ======================================================================

static int TimeNS_ho_RHSN_Terms = 3;
static MultiIndex2D TimeNS_ho_RHSDerivatives[4] = {  D10, D01, D00};
static int TimeNS_ho_RHSSpaceNumbers[1] = { 0 };
static int TimeNS_ho_RHSN_Matrices = 0;
static int *TimeNS_ho_RHSRowSpace = NULL;
static int *TimeNS_ho_RHSColumnSpace = NULL;
static int TimeNS_ho_RHSN_Rhs = 2;
static int TimeNS_ho_RHSRhsSpace[2] = { 0, 0 };


void TimeNSType1GalerkinRHS(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSType1GalerkinJ(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSType1GalerkinC(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


void TimeNSType3GalerkinJ(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);
                
// ======================================================================
// right-hand side ONLY, defect correction type 0, u2
// ======================================================================
void TimeNSRHSDefectCorrectionU2(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);
// ======================================================================
// right-hand side ONLY, defect correction type 1, u2
// ======================================================================
void TimeNSRHSDefectCorrectionU2_1(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);

// ======================================================================

void TimeNSType4GalerkinDD_Axial3D(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLGalerkinDD_Axial3D(double Mult, double *coeff,
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSType4GalerkinDD_2PhaseAxial3D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLGalerkinDD_2PhaseAxial3D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void GridAssemble4(double Mult, double *coeff,
                  double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs);

// ======================================================================

// D(u):D(v) TNSECST Galerkin
void Time_NSEType4_Galerkin_CST_Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// D(u):D(v) TNSECST Galerkin for two phase Axial3D flows 
void Time_NSEType4_Galerkin_CST_Galerkin_2PhaseAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


// D(u):D(v) TNSECST DEVSS/SUPG
void Time_NSEType4_Galerkin_CST_SUPG_DEVSS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// D(u):D(v) TNSECST DEVSS/SUPG non-linear terms only
void Time_NSEType4_Galerkin_CST_SUPG_DEVSS_NLTerms(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// D(u):D(v) TNSECST LPS non-linear terms only
void Time_NSEType4_Galerkin_CST_Galerkin_NLTerms(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// D(u):D(v) TNSECST LPS non-linear terms only for two phase axial flows
void Time_NSEType4_Galerkin_CST_Galerkin_NLTerms_2PhaseAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


// ======================================================================
// right-hand side ONLY, for TNSECST, DEVSS/SUPG
// ======================================================================
void Time_NSEType4_Galerkin_CST_SUPG_RhsOnly(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, for TNSECST, LPS
// ======================================================================
void Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, for TNSECST, LPS for 2Phase Axial3D flows
// ======================================================================
void Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly_2PhaseAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// D(u):D(v) TNSECST Galerkin for two phase Axial3D flows 
void Time_NSEType4_Galerkin_CST_Galerkin_ImpDropAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


// D(u):D(v) TNSECST LPS non-linear terms only for impinging droplet 3D-axial flows
void Time_NSEType4_Galerkin_CST_Galerkin_NLTerms_ImpDropAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, for TNSECST, LPS for Impinging droplet Axial3D flows
// ======================================================================
void Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly_ImpDropAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);






#endif  // __TNSE2D_FIXPO__
