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
// TNSE2D_FixPo_SSMUM.h        
//
// declarations for SSMUM
//
// author: Volker John 08/05/22
//
// ======================================================================

// ======================================================================
//
// WITH ROTATING FRAME
//
// ======================================================================

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, Coletti, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
			       double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, Coletti, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, Coletti, (grad u, grad v), only nonlinear part
// Type 4, Coletti, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				   double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				     double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, for SSMUM
// ======================================================================
void TimeNSRHS_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
		     double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, for SSMUM
// ======================================================================
void TimeNS_REL_VELO_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
			   double ***LocMatrices, double **LocRhs);


// ======================================================================
//
// ALE
//
// ======================================================================

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, Coletti, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				   double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, Coletti, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				     double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, Coletti, (grad u, grad v), only nonlinear part
// Type 4, Coletti, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				       double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
					 double ***LocMatrices, double **LocRhs);

// ======================================================================
// right-hand side ONLY, for NSE
// ======================================================================
void TimeNSRHS_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
			 double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 4
//      main diagonal blocks A11, A22
//      WITHOUT B1T, B2T (gradient blocks)
//      WITHOUT right hand sides
// ======================================================================

int TimeNSType4_SSMUM_NLN_Terms = 3;
MultiIndex2D TimeNSType4_SSMUM_NLDerivatives[3] = { D10, D01, D00 };
int TimeNSType4_SSMUM_NLSpaceNumbers[3] = { 0, 0, 0 };
int TimeNSType4_SSMUM_NLN_Matrices = 4;
int TimeNSType4_SSMUM_NLRowSpace[4] = { 0, 0, 0, 0 };
int TimeNSType4_SSMUM_NLColumnSpace[4] = { 0, 0, 0, 0 };
int TimeNSType4_SSMUM_NLN_Rhs = 0;
int *TimeNSType4_SSMUM_NLRhsSpace = NULL;
