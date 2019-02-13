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
// @(#)TNSE3D_FixPo.h        1.2 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#ifndef __TNSE3D_FIXPO__
#define __TNSE3D_FIXPO__

// ======================================================================
// declaration for all Navier-Stokes problems of type 1
//      one A block, 
//      M block from time discretization
//      B1, B2 (divergence blocks)
// ======================================================================

static int TimeNSType1N_Terms = 5;
static MultiIndex3D TimeNSType1Derivatives[5] = { D100, D010, D001, D000, D000 };
static int TimeNSType1SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
static int TimeNSType1N_Matrices = 5;
static int TimeNSType1RowSpace[5] = { 0, 0, 1, 1, 1 };
static int TimeNSType1ColumnSpace[5] = { 0, 0, 0, 0, 0 };
static int TimeNSType1N_Rhs = 3;
static int TimeNSType1RhsSpace[3] = { 0, 0, 0 };


double TurbulentViscosity3D(double delta, double* gradU, double* u, 
			    double* uConv, double* x, double* y,
			    double* z, double proj_space);

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void TimeNSType1Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Smagorinsky
// ======================================================================
void TimeNSType1Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, ClassicalLES
// ======================================================================
void TimeNSType1ClassicalLES3D(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, GL00Convolution
// ======================================================================
void TimeNSType1GL00Convolution3D(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, GL00AuxProblem
// ======================================================================
static int TimeNSType1GL00AuxProblemN_Matrices = 6;
static int TimeNSType1GL00AuxProblemRowSpace[6] = { 0, 0, 2, 1, 1, 1 };
static int TimeNSType1GL00AuxProblemColumnSpace[6] = { 0, 0, 2, 0, 0, 0 };

void TimeNSType1GL00AuxProblem3D(double Mult, double *coeff,
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

static int TimeNSType2N_Terms = 5;
static MultiIndex3D TimeNSType2Derivatives[5] = { D100, D010, D001, D000, D000 };
static int TimeNSType2SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
static int TimeNSType2N_Matrices = 8;
static int TimeNSType2RowSpace[8] = { 0, 0, 1, 1, 1, 0, 0, 0 };
static int TimeNSType2ColumnSpace[8] = { 0, 0, 0, 0, 0, 1, 1, 1 };
static int TimeNSType2N_Rhs = 3;
static int TimeNSType2RhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// Type 2, Standard Galerkin
// ======================================================================
void TimeNSType2Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 2 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void TimeNSType2Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, ClassicalLES
// ======================================================================
void TimeNSType2ClassicalLES3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, GL00Convolution
// ======================================================================
void TimeNSType2GL00Convolution3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 2, GL00AuxProblem
// ======================================================================
static int TimeNSType2GL00AuxProblemN_Matrices = 9;
static int TimeNSType2GL00AuxProblemRowSpace[9] = { 0, 0, 2, 1, 1, 1, 0, 0, 0 };
static int TimeNSType2GL00AuxProblemColumnSpace[9] = { 0, 0, 2, 0, 0, 0, 1, 1, 1 };

void TimeNSType2GL00AuxProblem3D(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all four A blocks,
//      B1, B2 (divergence blocks), 
// ======================================================================

static int TimeNSType3N_Terms = 5;
static MultiIndex3D TimeNSType3Derivatives[5] = { D100, D010, D001, D000, D000 };
static int TimeNSType3SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
static int TimeNSType3N_Matrices = 15;
static int TimeNSType3RowSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                1, 1, 1 };
static int TimeNSType3ColumnSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0 };
static int TimeNSType3N_Rhs = 3;
static int TimeNSType3RhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void TimeNSType3Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, ClassicalLES, (grad u, grad v)
// ======================================================================
void TimeNSType3ClassicalLES3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, ClassicalLES, D(u):D(v)
// ======================================================================
void TimeNSType3ClassicalLESDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType3GL00Convolution3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType3GL00ConvolutionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00AuxProblem, (grad u, grad v)
// ======================================================================
static int TimeNSType3GL00AuxProblemN_Matrices = 16;
static int TimeNSType3GL00AuxProblemRowSpace[16] =  { 0, 0, 0, 0, 0, 0, 
                                               0, 0, 0, 0, 0, 0, 
                                               2, 1, 1, 1 };
static int TimeNSType3GL00AuxProblemColumnSpace[16] = { 0, 0, 0, 0, 0, 0, 
                                                 0, 0, 0, 0, 0, 0, 
                                                 2, 0, 0, 0 };

void TimeNSType3GL00AuxProblem3D(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType3GL00AuxProblemDD3D(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      all four A blocks,
//      B1, B2 (divergence blocks), 
//      B1T, B2T (gradient blocks)
// ======================================================================

static int TimeNSType4N_Terms = 5;
static MultiIndex3D TimeNSType4Derivatives[5] = { D100, D010, D001, D000, D000 };
static int TimeNSType4SpaceNumbers[5] = { 0, 0, 0, 0, 1 };
static int TimeNSType4N_Matrices = 18;
static int TimeNSType4RowSpace[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                1, 1, 1, 0, 0, 0 };
static int TimeNSType4ColumnSpace[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 1, 1, 1 };
static int TimeNSType4N_Rhs = 3;
static int TimeNSType4RhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void TimeNSType4Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void TimeNSType4UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, ClassicalLES, (grad u, grad v)
// ======================================================================
void TimeNSType4ClassicalLES3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, ClassicalLES, D(u):D(v)
// ======================================================================
void TimeNSType4ClassicalLESDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4GL00Convolution3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GL00ConvolutionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00AuxProblem, (grad u, grad v)
// ======================================================================
static int TimeNSType4GL00AuxProblemN_Matrices = 19;
static int TimeNSType4GL00AuxProblemRowSpace[19] =  
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 0, 0, 0 };
static int TimeNSType4GL00AuxProblemColumnSpace[19] = 
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1 };

void TimeNSType4GL00AuxProblem3D(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType4GL00AuxProblemDD3D(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, VMS_Projection, D(u):D(v)
// ======================================================================
static int TimeNSType4VMS_ProjectionN_Terms = 6;
static MultiIndex3D TimeNSType4VMS_ProjectionDerivatives[6] = { D100, D010, D001, D000, D000, D000 };
static int TimeNSType4VMS_ProjectionSpaceNumbers[6] = { 0, 0, 0, 0, 1, 2 };
static int TimeNSType4VMS_ProjectionN_Matrices = 25;
static int TimeNSType4VMS_ProjectionRowSpace[25] = { 0, 0, 0, 0, 0, 0, 
                                              0, 0, 0, 0, 0, 0,
                                              2, 1, 1, 1, 0, 0, 0,
                                              0, 0, 0, 2, 2, 2};
                                            
static int TimeNSType4VMS_ProjectionColumnSpace[25] = { 0, 0, 0, 0, 0, 0, 
                                                 0, 0, 0, 0, 0, 0,
                                                 2, 0, 0, 0, 1, 1, 1,
                                                 2, 2, 2, 0, 0, 0};
                                                 
void TimeNSType4VMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 4, VMS_Projection with streamline formulation D(u):D(v)
// ======================================================================

void TimeNSType4VMS_ProjectionStreamlineDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
		double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for SUPG
// ONLY right hand sides
// ======================================================================

static int TimeNSVMS_Rhs_SUPGN_Terms = 7;
static MultiIndex3D TimeNSVMS_Rhs_SUPGDerivatives[7] = { D100, D010, D001, D000,
                                                         D100, D010, D001 };
static int TimeNSVMS_Rhs_SUPGSpaceNumbers[7] = { 0, 0, 0, 0, 1, 1, 1 };
static int TimeNSVMS_Rhs_SUPGN_Matrices = 0;
static int *TimeNSVMS_Rhs_SUPGRowSpace = NULL;
static int *TimeNSVMS_Rhs_SUPGColumnSpace = NULL;
static int TimeNSVMS_Rhs_SUPGN_Rhs = 7;
static int TimeNSVMS_Rhs_SUPGRhsSpace[7] = { 0, 0, 0, 0, 0, 0, 1};

// compute stabilization parameters
void SUPG_Param3D(double u1, double u2, double u3, double* coeff, double* params);

void TimeNSType14VMS_Rhs_SUPGDD3D(double Mult, double *coeff, double *param,
                                  double hK, double **OrigValues,
                                  int *N_BaseFuncts, double ***LocMatrices,
                                  double **LocRhs);


// ======================================================================
// Type 4, VMS_SUPG, D(u):D(v)
// ======================================================================
static int TimeNSType4VMS_SUPGN_Terms = 8;
static MultiIndex3D TimeNSType4VMS_SUPGDerivatives[8] = { D100, D010, D001, D000, 
                                                          D100, D010, D001, D000 };
static int TimeNSType4VMS_SUPGSpaceNumbers[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
static int TimeNSType4VMS_SUPGN_Matrices = 28;
static int TimeNSType4VMS_SUPGRowSpace[28] = { 0, 0, 0, 0, 0, 0, 
                                               0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 
					       1, 1, 1, 0, 0, 0};
                                            
static int TimeNSType4VMS_SUPGColumnSpace[28] = {0, 0, 0, 0, 0, 0, 
                                               0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 
					       0, 0, 0, 1, 1, 1};
                                                 
static int TimeNSType4VMS_SUPGN_Rhs = 4;
static int TimeNSType4VMS_SUPGRhsSpace[4] = { 0, 0, 0, 1};

void TimeNSType4VMS_SUPGDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

static int TimeNSType4NLVMS_SUPGN_Terms = 8;
static MultiIndex3D TimeNSType4NLVMS_SUPGDerivatives[8] = { D100, D010, D001, D000, 
							  D100, D010, D001, D000 };
static int TimeNSType4NLVMS_SUPGSpaceNumbers[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
static int TimeNSType4NLVMS_SUPGN_Matrices = 22;
static int TimeNSType4NLVMS_SUPGRowSpace[22] = { 0, 0, 0, 0, 0, 0, 
                                               0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0};
                                            
static int TimeNSType4NLVMS_SUPGColumnSpace[22] = {0, 0, 0, 0, 0, 0, 
                                               0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0,
						   0, 0, 0, 0};

static int TimeNSType4NLVMS_SUPGN_Rhs = 0;
static int *TimeNSType4NLVMS_SUPGRhsSpace = NULL;

void TimeNSType4NLVMS_SUPGDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for SUPG
// ONLY right hand sides
// ======================================================================

void TimeNSType4VMS_Rhs_SUPGDD3D(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 14, VMS_SUPG, D(u):D(v)
// ======================================================================
static int TimeNSType14VMS_SUPGN_Terms =8;
static MultiIndex3D TimeNSType14VMS_SUPGDerivatives[8] = { D100, D010, D001, D000, 
                             D100, D010, D001, D000};
static int TimeNSType14VMS_SUPGSpaceNumbers[11] = { 0, 0, 0, 0, 1, 1, 1, 1};
static int TimeNSType14VMS_SUPGN_Matrices = 37;
static int TimeNSType14VMS_SUPGRowSpace[37] = { 0, 0, 0, 0, 0, 0, 
                                                0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0,
                  1, 1, 1, 1, 0, 0, 0};
                                           
static int TimeNSType14VMS_SUPGColumnSpace[37] = {0, 0, 0, 0, 0, 0, 
                                                  0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0,
                    1, 0, 0, 0, 1, 1, 1};
                                                 
static int TimeNSType14VMS_SUPGN_Rhs = 7;
static int TimeNSType14VMS_SUPGRhsSpace[7] = { 0, 0, 0, 0, 0, 0, 1};

void TimeNSType14VMS_SUPGDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

static int TimeNSType14NLVMS_SUPGN_Terms = 8;
static MultiIndex3D TimeNSType14NLVMS_SUPGDerivatives[8] = { D100, D010, D001, D000, 
                   D100, D010, D001, D000};
static int TimeNSType14NLVMS_SUPGSpaceNumbers[8] = { 0, 0, 0, 0, 1, 1, 1, 1};
static int TimeNSType14NLVMS_SUPGN_Matrices = 34;
static int TimeNSType14NLVMS_SUPGRowSpace[34] = { 0, 0, 0, 0, 0, 0, 
                                                  0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 1, 1,
              1, 0, 0, 0};
                                            
static int TimeNSType14NLVMS_SUPGColumnSpace[34] = {0, 0, 0, 0, 0, 0, 
                                                    0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 1, 0, 0,
                0, 1, 1, 1};

static int TimeNSType14NLVMS_SUPGN_Rhs = 3;
static int TimeNSType14NLVMS_SUPGRhsSpace[3] = {0, 0, 0};

void TimeNSType14NLVMS_SUPGDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Leray-alpha  D(u):D(v)
// ======================================================================
void TimeNSType4LerayAlphaDD3D(double Mult, double *coeff,
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

static int TimeNSType1_2NLN_Terms = 4;
static MultiIndex3D TimeNSType1_2NLDerivatives[4] = { D100, D010, D001, D000 };
static int TimeNSType1_2NLSpaceNumbers[4] = { 0, 0, 0, 0 };
static int TimeNSType1_2NLN_Matrices = 1;
static int TimeNSType1_2NLRowSpace[1] = { 0 };
static int TimeNSType1_2NLColumnSpace[1] = { 0 };
static int TimeNSType1_2NLN_Rhs = 0;
static int *TimeNSType1_2NLRhsSpace = NULL;

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// Type 2, Standard Galerkin, only nonlinear part
// Type 1, ClassicalLES, only nonlinear part
// Type 2, ClassicalLES, only nonlinear part
// Type 1, GL00Convolution, only nonlinear part
// Type 2, GL00Convolution, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// Type 2, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1_2NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      main diagonal blocks A11, A22
//      WITHOUT right hand sides
// ======================================================================

static int TimeNSType3_4NLN_Terms = 4;
static MultiIndex3D TimeNSType3_4NLDerivatives[4] = { D100, D010, D001, D000 };
static int TimeNSType3_4NLSpaceNumbers[4] = { 0, 0, 0, 0 };
static int TimeNSType3_4NLN_Matrices = 3;
static int TimeNSType3_4NLRowSpace[3] = { 0, 0, 0 };
static int TimeNSType3_4NLColumnSpace[3] = { 0, 0, 0 };
static int TimeNSType3_4NLN_Rhs = 0;
static int *TimeNSType3_4NLRhsSpace = NULL;

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, ClassicalLES, (grad u, grad v), only nonlinear part
// Type 4, ClassicalLES, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, ClassicalLES, D(u):D(v), only nonlinear diagonal blocks
// Type 4, ClassicalLES, D(u):D(v), only nonlinear diagonal blocks
// Type 3, GL00Convolution, D(u):D(v), only nonlinear diagonal blocks
// Type 4, GL00Convolution, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// Type 4, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3_4NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// Type 4, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3_4NLUpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Standard Smagorinsky, (grad u, grad v), only nonlinear part
// ======================================================================
static int TimeNSType3_4NLSmagorinskyN_Matrices = 9;
static int TimeNSType3_4NLSmagorinskyRowSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static int TimeNSType3_4NLSmagorinskyColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };


void TimeNSType3_4NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
static int TimeNSType3_4NLVMS_ProjectionN_Terms = 5;
static MultiIndex3D TimeNSType3_4NLVMS_ProjectionDerivatives[5] = { D100, D010, D001, D000, D000 };
static int TimeNSType3_4NLVMS_ProjectionSpaceNumbers[5] = { 0, 0, 0, 0, 2 };
static int TimeNSType3_4NLVMS_ProjectionN_Matrices = 12;
static int TimeNSType3_4NLVMS_ProjectionRowSpace[12] = { 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0};
static int TimeNSType3_4NLVMS_ProjectionColumnSpace[12] = { 0, 0, 0, 0, 0, 0, 
                                                   0, 0, 0, 2, 2, 2};

void TimeNSType3_4NLVMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), adaptive coarse space
// Type 4, VMS_Projection, D(u):D(v), adaptive coarse space
// ======================================================================
static int TimeNSType3_4NL_Adap_VMS_ProjectionN_Terms = 5;
static MultiIndex3D TimeNSType3_4NL_Adap_VMS_ProjectionDerivatives[5] = { D100, D010, D001, D000, D000 };
static int TimeNSType3_4NL_Adap_VMS_ProjectionSpaceNumbers[5] = { 0, 0, 0, 0, 2 };
static int TimeNSType3_4NL_Adap_VMS_ProjectionN_Matrices = 16;
static int TimeNSType3_4NL_Adap_VMS_ProjectionRowSpace[16] = { 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 2, 0, 0, 0, 2, 2, 2};
static int TimeNSType3_4NL_Adap_VMS_ProjectionColumnSpace[16] = { 0, 0, 0, 0, 0, 0, 
                                                   0, 0, 0, 2, 2, 2, 2, 0, 0, 0};

void TimeNSType3_4NL_Adap_VMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// streamline projection
// ======================================================================
void TimeNSType3_4NLVMS_ProjectionStreamlineDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
		double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, VMS_Projection explicit, only Matrix_tilde_G??
// Type 4, VMS_Projection explicit, only Matrix_tilde_G??
// ======================================================================

static int TimeNSType3_4NLVMS_ProjectionExplN_Terms = 4;
static MultiIndex3D TimeNSType3_4NLVMS_ProjectionExplDerivatives[4] = { D100, D010, D001, D000};
static int TimeNSType3_4NLVMS_ProjectionExplSpaceNumbers[4] = { 0, 0, 0, 1 };
static int TimeNSType3_4NLVMS_ProjectionExplN_Matrices = 3;
static int TimeNSType3_4NLVMS_ProjectionExplRowSpace[3] = {0, 0, 0  };
static int TimeNSType3_4NLVMS_ProjectionExplColumnSpace[3] = {1, 1, 1 };

void TimeNSType3_4VMS_ProjectionExpl3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, div-div stabilization
// Type 4, div-div stabilization
// ======================================================================
void TimeNSType3_4NLDivDivDD3D(double Mult, double *coeff, 
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

static int TimeNSRHSN_Terms = 1;
static MultiIndex3D TimeNSRHSDerivatives[1] = { D000 };
static int TimeNSRHSSpaceNumbers[1] = { 0 };
static int TimeNSRHSN_Matrices = 0;
static int *TimeNSRHSRowSpace = NULL;
static int *TimeNSRHSColumnSpace = NULL;
static int TimeNSRHSN_Rhs = 3;
static int TimeNSRHSRhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// right-hand side ONLY
// ======================================================================
void TimeNSRHS3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);

static int TimeNSRHSLESN_Terms = 4;
static MultiIndex3D TimeNSRHSLESDerivatives[4] = { D100, D010, D001, D000 };
static int TimeNSRHSLESSpaceNumbers[4] = { 0, 0, 0, 0 };
static int TimeNSRHSLESN_Matrices = 0;
static int *TimeNSRHSLESRowSpace = NULL;
static int *TimeNSRHSLESColumnSpace = NULL;
static int TimeNSRHSLESN_Rhs = 3;
static int TimeNSRHSLESRhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// right-hand side ONLY, ClassicalLES model
// ======================================================================
void TimeNSRHSClassicalLES3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, Galdi-Layton model with convolution
// right-hand side ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================

void TimeNSRHSLESModel3D(double Mult, double *coeff, 
                         double *param, double hK, 
                         double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side for auxiliary problem 
// ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================

static int TimeNSGL00AuxProblemRHSN_Terms = 1;
static MultiIndex3D TimeNSGL00AuxProblemRHSDerivatives[1] = { D000 };
static int TimeNSGL00AuxProblemRHSSpaceNumbers[1] = { 0 };
static int TimeNSGL00AuxProblemRHSN_Matrices = 0;
static int *TimeNSGL00AuxProblemRHSRowSpace = NULL;
static int *TimeNSGL00AuxProblemRHSColumnSpace = NULL;
static int TimeNSGL00AuxProblemRHSN_Rhs = 6;
static int TimeNSGL00AuxProblemRHSRhsSpace[6] = { 0, 0, 0, 0, 0, 0 };

void TimeNSGL00AuxProblemRHS3D(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);


// ======================================================================
// right-hand side ONLY for auxiliary problem applied to velocity
// ======================================================================
void TimeNSRHSAuxProblemU(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for computation of rhs for RFB explicit
//    only rhs
// ======================================================================

static int TimeNSRFBExplRhsN_Terms = 1;
static MultiIndex3D TimeNSRFBExplRhsDerivatives[1] = { D000 };
static int TimeNSRFBExplRhsSpaceNumbers[1] = { 0 };
static int TimeNSRFBExplRhsN_Matrices = 0;
static int *TimeNSRFBExplRhsRowSpace = NULL;
static int *TimeNSRFBExplRhsColumnSpace = NULL;
static int TimeNSRFBExplRhsN_Rhs = 3;
static int TimeNSRFBExplRhsRhsSpace[3] = { 0, 0, 0 };

// ======================================================================
// Assembling routine for matrices for auxiliary problems
// ======================================================================

// ======================================================================
// declarations
// ======================================================================

//int TimeNSGL00AuxProblemN_Terms = 3;
//MultiIndex3D TimeNSGL00AuxProblemDerivatives[3] = { D10, D01, D00 };
//int TimeNSGL00AuxProblemSpaceNumbers[3] = { 0, 0, 0 };
//int TimeNSGL00AuxProblemN_Matrices = 1;
//int TimeNSGL00AuxProblemRowSpace[2] = { 0 };
//int TimeNSGL00AuxProblemColumnSpace[2] = { 0 };
//int TimeNSGL00AuxProblemN_Rhs = 0;
//int *TimeNSGL00AuxProblemRhsSpace = NULL;

//void TimeNSAuxMatrixGL00AuxProblem3D(double Mult, double *coeff,
//               double *param, double hK,
//               double **OrigValues, int *N_BaseFuncts,
//               double ***LocMatrices, double **LocRhs);


static int MatrixAuxiliaryProblemN_Terms = 4;
static MultiIndex3D MatrixAuxiliaryProblemDerivatives[4] = { D100, D010, D001, D000};
static int MatrixAuxiliaryProblemSpaceNumbers[4] = { 0, 0, 0, 0};
static int MatrixAuxiliaryProblemN_Matrices = 1;
static int MatrixAuxiliaryProblemRowSpace[1] = { 0 };
static int MatrixAuxiliaryProblemColumnSpace[1] = { 0 };
static int MatrixAuxiliaryProblemN_Rhs = 0;
static int *MatrixAuxiliaryProblemRhsSpace = NULL;

void MatrixAuxiliaryProblem(double Mult, double *coeff, 
                            double *param, double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);


// ======================================================================
// right-hand side for additional terms in rhs of small scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_SmallRhs3D(double Mult, double *coeff, 
                           double *param, double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for VMS
// ONLY right hand sides
// ======================================================================

static int TimeNS_ho_RHSN_Terms = 4;
static MultiIndex3D TimeNS_ho_RHSDerivatives[4] = {  D100, D010, D001, D000};
static int TimeNS_ho_RHSSpaceNumbers[1] = { 0 };
static int TimeNS_ho_RHSN_Matrices = 0;
static int *TimeNS_ho_RHSRowSpace = NULL;
static int *TimeNS_ho_RHSColumnSpace = NULL;
static int TimeNS_ho_RHSN_Rhs = 3;
static int TimeNS_ho_RHSRhsSpace[3] = { 0, 0, 0 };

void TimeNSType1GalerkinJ3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSGalerkinC3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


void TimeNSType3GalerkinJ3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);




#endif
