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
// %W% %G%
// 
// Purpose:     assemble matrix and right-hand side
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __ASSEMBLE3D__
#define __ASSEMBLE3D__

#include <AllClasses.h>
#include <Constants.h>

/** a function from a finite element space */
void Assemble3D(int n_fespaces, TFESpace3D **fespaces,
                int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                int n_matrices, TMatrix3D **matrices,
                int n_rhs, double **rhs, TFESpace3D **ferhs,
                TDiscreteForm3D *DiscreteForm,
                BoundCondFunct3D **BoundaryConditions,
                BoundValueFunct3D **BoundaryValues,
                TAuxParam3D *parameters);


/** a function from a finite element space */
void Assemble3DSlipBC(int n_fespaces, TFESpace3D **fespaces,
                int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                int n_matrices, TMatrix3D **matrices,
                int n_rhs, double **rhs, TFESpace3D **ferhs,
                TDiscreteForm3D *DiscreteForm,
                BoundCondFunct3D **BoundaryConditions,
                BoundValueFunct3D **BoundaryValues,
                TAuxParam3D *parameters);

void ModifyMatrixSlipBC(TSquareMatrix3D **sqmatrices, TMatrix3D **matrices,
			int N_U, double *rhs);

/** assemble mixed finite elements such as Raviart-Thomas or
 * Brezzi-Douglas-Marini.
 */
void Assemble3D_mixed(int n_fespaces, TFESpace3D **fespaces,
int n_sqmatrices, TSquareMatrix3D **sqmatrices,
int n_matrices, TMatrix3D **matrices,
int n_rhs, double **rhs, TFESpace3D **ferhs,
TDiscreteForm3D *DiscreteForm3D,
BoundCondFunct3D **BoundaryConditions,
BoundValueFunct3D **BoundaryValues,
TAuxParam3D *Parameters);



#endif // __ASSEMBLE3D__
