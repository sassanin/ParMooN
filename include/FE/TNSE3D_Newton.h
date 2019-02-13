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
// TNSE3D_Newton.h
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#ifndef __TNSE3D_NEWTON__
#define __TNSE3D_NEWTON__

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      all four A blocks,
//      B1, B2 (divergence blocks), 
// ======================================================================

static int TimeNSType3NewtonN_Terms = 5;
static MultiIndex3D TimeNSType3NewtonDerivatives[5] = { D100, D010, D001, D000, D000 };
static int TimeNSType3NewtonSpaceNumbers[5] = { 0, 0, 0, 0, 1 };
static int TimeNSType3NewtonN_Matrices = 15;
static int TimeNSType3NewtonRowSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                1, 1, 1 };
static int TimeNSType3NewtonColumnSpace[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0 };
static int TimeNSType3NewtonN_Rhs = 3;
static int TimeNSType3NewtonRhsSpace[3] = { 0, 0, 0 };

void TimeNSType3GalerkinNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);

void TimeNSType3UpwindNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);

void TimeNSType4VMS_ProjectionNewtonDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
					 double ***LocMatrices, double **LocRhs);

void TimeNSType4SmagorinskyNewtonDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// declaration for all Navier-Stokes problems of type 3
//      WITHOUT right hand sides
// ======================================================================

static int TimeNSType3_4NLNewtonN_Terms = 4;
static MultiIndex3D TimeNSType3_4NLNewtonDerivatives[4] = { D100, D010, D001, D000 };
static int TimeNSType3_4NLNewtonSpaceNumbers[4] = { 0, 0, 0, 0 };
static int TimeNSType3_4NLNewtonN_Matrices = 3;
static int TimeNSType3_4NLNewtonRowSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0 ,0 };
static int TimeNSType3_4NLNewtonColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0 ,0 };
static int TimeNSType3_4NLNewtonN_Rhs = 0;
static int *TimeNSType3_4NLNewtonRhsSpace = NULL;

void TimeNSType3_4NLGalerkinNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
				     double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLUpwindNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLSmagorinskyNewtonDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLVMS_ProjectionNewtonDD3D(double Mult, double *coeff, 
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

static int TimeNSRHSNewtonN_Terms = 1;
static MultiIndex3D TimeNSRHSNewtonDerivatives[1] = { D000 };
static int TimeNSRHSNewtonSpaceNumbers[1] = { 0 };
static int TimeNSRHSNewtonN_Matrices = 0;
static int *TimeNSRHSNewtonRowSpace = NULL;
static int *TimeNSRHSNewtonColumnSpace = NULL;
static int TimeNSRHSNewtonN_Rhs = 6;
static int TimeNSRHSNewtonRhsSpace[6] = { 0, 0, 0, 0, 0, 0 };

// ======================================================================
// right-hand side ONLY
// ======================================================================
void TimeNSRHSNewton3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);


void TimeNSRHSNewtonNL3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);


#endif
