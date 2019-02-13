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
// NSE2D_FixPoSkew.h        
//
// common declaration for all Navier-Stokes problems
// ======================================================================

#ifndef __NSE2DFIXPOSKEW__
#define __NSE2DFIXPOSKEW__

#include <Enumerations.h>


// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1GalerkinSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, (reduced) SDFEM
// ======================================================================
void NSType1SDFEMSkew(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, for upwind not available
// ======================================================================

// ======================================================================
// Type 1, Smagorinsky 
// ======================================================================
void NSType1SmagorinskySkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Standard Galerkin
// ======================================================================
void NSType2GalerkinSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2SDFEMSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Upwind  not available
// ======================================================================

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void NSType2SmagorinskySkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3GalerkinSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3GalerkinSkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind not available
// ======================================================================

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3SmagorinskySkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3SmagorinskySkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType4GalerkinSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4GalerkinSkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEMSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SDFEMSkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, not available
// ======================================================================

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4SmagorinskySkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType4SmagorinskySkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkinSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, (reduced) SDFEM or (simplified RFB)
// ======================================================================
void NSType1NLSDFEMSkew(double Mult, double *coeff,
			double *param, double hK,
			double **OrigValues, int *N_BaseFuncts,
			double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, upwind not available
// ======================================================================

// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// ======================================================================
void NSType1_2NLSmagorinskySkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2NLSDFEMSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkinSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinSkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 3, Upwind not available
// ======================================================================

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3_4NLSmagorinskySkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3_4NLSmagorinskySkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEMSkew(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMSkewDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

#endif
