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
// @(#)NSE2D_Newton.h        1.2 04/13/00
//
// common declaration for all Navier-Stokes problems
// ======================================================================

#ifndef __NSE2D_Newton__
#define __NSE2D_Newton__

#include <Enumerations.h>

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3GalerkinNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3GalerkinDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3UpwindNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3UpwindDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3SmagorinskyNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3SmagorinskyDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType4GalerkinNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4GalerkinDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEMNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SDFEMDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void NSType4UpwindNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void NSType4UpwindDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4SmagorinskyNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType4SmagorinskyDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3NLGalerkinNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3NLGalerkinDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3NLUpwindNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3NLUpwindDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3NLSmagorinskyNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3NLSmagorinskyDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType4NLGalerkinNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4NLGalerkinDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEMNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void NSType4NLUpwindNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void NSType4NLUpwindDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v), nonlinear
// ======================================================================
void NSType4NLSmagorinskyNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v), nonlinear
// ======================================================================
void NSType4NLSmagorinskyDDNewton(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


#endif
