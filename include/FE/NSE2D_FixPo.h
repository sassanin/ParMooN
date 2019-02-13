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
// @(#)NSE2D_FixPo.h        1.2 04/13/00
//
// common declaration for all Navier-Stokes problems
// ======================================================================

#ifndef __NSE2DFIXPO__
#define __NSE2DFIXPO__

#include <Enumerations.h>

// ======================================================================
// compute parameter for RFB stabilization
// Brezzi, Marini, Russo, CMAME 194 (2005) 127 - 148
// ======================================================================

double RFB_Parameter(double hK, double eps, double* b);

// ======================================================================
// compute parameter for SUPG stabilization
// ======================================================================
double SUPG_Parameter(double hK, double eps, double b1, double b2, double c);


void NSParamsVelo(double *in, double *out);

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, (reduced) SDFEM
// ======================================================================
void NSType1SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Smagorinsky
// ======================================================================
void NSType1Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, VMSProjection, nabla form
// ======================================================================
void NSType1VMSProjection(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Standard Galerkin + divergence term 
// ======================================================================
void NSType1GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Standard Galerkin
// ======================================================================
void NSType2Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void NSType2Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void NSType2Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Standard Galerkin + divergence term
// ======================================================================
void NSType2GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3GalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin +div term, (grad u, grad v)
// ======================================================================
void NSType3GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType4Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4GalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SDFEMDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void NSType4Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void NSType4UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType4SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, VMSProjection, D(u):D(v)
// ======================================================================
void NSType4VMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1,(reduced) SDFEM
// ======================================================================
void NSType1NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1_2NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// ======================================================================
void NSType1_2NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, VMSProjection, nabla form
// Type 2, VMSProjection, nabla form
// ======================================================================
void NSType1_2NLVMSProjection(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Standard Galerkin + div term, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3_4NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3_4NLUpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3_4NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3_4NLSmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin + div term, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// rhs for RFB stabilization
// ======================================================================
void NSRFBRhs(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
	      double ***LocMatrices, double **LocRhs);

// ======================================================================
// pressure separation
// ======================================================================
void NSPressSep(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// pressure separation with auxiliary problem
// ======================================================================
void NSPressSepAuxProb(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


// ======================================================================
// Standard Galerkin NSE-Type 4 and SUPG for CST 
// ======================================================================
void NSEType4_Galerkin_CST_SUPG_DEVSS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


// ======================================================================
// Standard Galerkin NSE-Type 4 and Galerkin for CST 
// ======================================================================
void NSEType4_Galerkin_CST_Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Standard Galerkin NSE-Type 4 and Galerkin for CST (Axial 3D case)
// ======================================================================
void NSEType4_Galerkin_CST_Galerkin_Axial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


// ======================================================================
// Standard Galerkin for CST
// ======================================================================
void CSTGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);
// ======================================================================
// SUPG for CST
// ======================================================================
void CST_SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);
// ======================================================================
// SUPG non-consistent for CST
// ======================================================================
void CST_SUFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Standard Galerkin for CST (Giesekus type)
// ======================================================================
void CST_GiesekusGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Standard Galerkin for DFT
// ======================================================================
void DFTGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);


// ======================================================================
// Rhs only for NSE  due to CST 
// ======================================================================
void NSCSTRhs(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Rhs only for NSE  due to CST 
// ======================================================================
void NSCSTRhs_DEVSS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ========================================================================
// routines for FJMT07
// ========================================================================

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1GalerkinFJMT07(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkinFJMT07(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

/*// ======================================================================
// auxiliary problem for differential filter
// ======================================================================
void Filter_Galerkin(double Mult, double *coeff,
                     double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);
*/
/*
// ======================================================================
// auxiliary problem
// ======================================================================
void NSAuxProblem(double Mult, double *coeff,
                   double *param, double hK,
                   double **OrigValues, int *N_BaseFuncts,
                   double ***LocMatrices, double **LocRhs);

int NSAuxProblemN_Terms = 3;
MultiIndex2D NSAuxProblemDerivatives[3] = { D10, D01, D00};
int NSAuxProblemSpaceNumbers[3] = { 0, 0, 0};
int NSAuxProblemN_Matrices = 1;
int NSAuxProblemRowSpace[1] = { 0 };
int NSAuxProblemColumnSpace[1] = { 0 };
int NSAuxProblemN_Rhs = 2;
int NSAuxProblemRhsSpace[2] = { 0, 0 };
*/
#endif
