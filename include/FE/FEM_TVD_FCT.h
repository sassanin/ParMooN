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
// FEM_TVD_FCT.h
//
// Purpose:     fem-tvd and fem-fct methods for scalar cdr equations
//
// Authors:     Volker John, 2008/02/07
//
// =======================================================================

#ifndef __FEM_TVD_FCT__
#define __FEM_TVD_FCT__

//**************************************************************
//compute the lumped matrix
//output is a vector
//**************************************************************
#ifdef __2D__
void LumpMassMatrixToVector(TSquareMatrix2D *M, double *lump_mass);
#endif    
#ifdef __3D__
void LumpMassMatrixToVector(TSquareMatrix3D *M, double *lump_mass);
#endif

/*******************************************************************************/
//
// FEM_TVD_ForConvDiff for steady-state cdr equations
// following D. Kuzmin (2007)
//
/*******************************************************************************/
#ifdef __2D__
void FEM_TVD_ForConvDiff(TSquareMatrix2D *sqmatrix, int N_U, int N_Active,
double *matrix_D_Entries, double *sol, double *rhs,
			 int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D);
#endif    
#ifdef __3D__
void FEM_TVD_ForConvDiff(TSquareMatrix3D *sqmatrix, int N_U, int N_Active,
double *matrix_D_Entries, double *sol, double *rhs,
int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D);
#endif

//**************************************************************
//compute system matrix for FEM-FCT
//output is a vector
//**************************************************************
#ifdef __2D__
void FEM_FCT_SystemMatrix(TSquareMatrix2D *M_C, TSquareMatrix2D *A, 
			  double *lump_mass,int N_U);
#endif    
#ifdef __3D__
void FEM_FCT_SystemMatrix(TSquareMatrix3D *M_C, TSquareMatrix3D *A, 
			  double *lump_mass,int N_U);
#endif

//**************************************************************
// MINMOD prelimiter
//**************************************************************
double MinMod(double a, double b);

/*******************************************************************************/
//
// FCT-FEM algorithm
// following D. Kuzmin, M. M"oller (2005) (nonlinear scheme)
//           D. Kuzmin (2008) (linear scheme)
//
/*******************************************************************************/
#ifdef __2D__
void FEM_FCT_ForConvDiff(TSquareMatrix2D *M_C,TSquareMatrix2D *A,
			 int N_U, int N_Active, 
			 double *lump_mass, double *matrix_D_Entries, 
			 double *sol, double *oldsol, 
			 double *B, double *rhs, double *rhs_old,
			 double *tilde_u,
			 int N_neum_to_diri, int *neum_to_diri,
			 int *neum_to_diri_bdry,
			 double *neum_to_diri_param, 
			 int compute_matrix_D, BoundValueFunct2D *BoundaryValue,
			 double *BoundaryValues);
#endif    
#ifdef __3D__
void FEM_FCT_ForConvDiff(TSquareMatrix3D *M_C,TSquareMatrix3D *A,
    int N_U, int N_Active, 
    double *lump_mass, double *matrix_D_Entries, 
    double *sol, double *oldsol, 
    double *B, double *rhs, double *rhs_old,
    double *tilde_u,
    int N_neum_to_diri, int *neum_to_diri,
    double *neum_to_diri_x, double *neum_to_diri_y, double *neum_to_diri_z, 
    int compute_matrix_D, BoundValueFunct3D *BoundaryValue,
    double *BoundaryValues);
#endif    

#endif
