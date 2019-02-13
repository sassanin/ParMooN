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
// @(#)Assemble2D.h        1.5 04/13/00
// 
// Purpose:     assemble matrix and right-hand side
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __ASSEMBLE2D__
#define __ASSEMBLE2D__

#include <AllClasses.h>
#include <Constants.h>
#include <FEDatabase2D.h>
#include <LocalAssembling2D.h>

#ifdef __3D__
  #include <Aux2D3D.h>
#endif

/** a function from a finite element space */
void Assemble2D(int n_fespaces, TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, TFESpace2D **ferhs,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                LocalAssembling2D& la
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
, int AssemblePhaseID = -1 
               );





/** a function from a finite element space */
void Assemble2D(int n_fespaces, TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, TFESpace2D **ferhs,
                TDiscreteForm2D *DiscreteForm,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                TAuxParam2D *parameters
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
, int AssemblePhaseID = -1 
               );

/** assembling of matrices multiplied by a factor */
void Assemble2D_FCT(int n_fespaces, TFESpace2D **fespaces, int n_sqmatrices, 
                    TSquareMatrix2D **sqmatrices, int n_matrices, 
                    TMatrix2D **matrices, int n_rhs, double **rhs, 
                    TFESpace2D **ferhs, TDiscreteForm2D *DiscreteForm,
                    BoundCondFunct2D **BoundaryConditions,
                    BoundValueFunct2D **BoundaryValues, TAuxParam2D *Parameters,
                    double factor
#ifdef __3D__
                    , TAux2D3D *Aux2D3D
#endif
                   );



/** assembling of slip type bc */
void Assemble2DSlipBC(int n_fespaces, TFESpace2D **fespaces,
                      int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                      int n_matrices, TMatrix2D **matrices,
                      int n_rhs, double **rhs, TFESpace2D **ferhs,
                      TDiscreteForm2D *DiscreteForm,
                      BoundCondFunct2D **BoundaryConditions,
                      BoundValueFunct2D **BoundaryValues,
                      TAuxParam2D *parameters,
                      TFEFunction2D *u1, TFEFunction2D *u2);

/** assembling for methods which need values on neighbour cells */
void Assemble2D_neigh(int n_fespaces, TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, TFESpace2D **ferhs,
                TDiscreteForm2D *DiscreteForm,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                TAuxParam2D *Parameters
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
                );

/** assembling for discontinuous Galerkin discretization */
void Assemble2D_DG(CoeffFct2D *Coeff, int n_fespaces, TFESpace2D **fespaces,
		   int n_sqmatrices, TSquareMatrix2D **sqmatrices,
		   int n_matrices, TMatrix2D **matrices,
		   int n_rhs, double **rhs, TFESpace2D **ferhs,
		   BoundCondFunct2D **BoundaryConditions,
		   BoundValueFunct2D **BoundaryValues,
		   TAuxParam2D *Parameters);

/** assembling for continuous interior penalty discretization */
void Assemble2D_CIP(CoeffFct2D *Coeff,int n_fespaces, TFESpace2D **fespaces,
		    int n_sqmatrices, TSquareMatrix2D **sqmatrices,
		    int n_matrices, TMatrix2D **matrices,
		    int n_rhs, double **rhs, TFESpace2D **ferhs,
		    BoundCondFunct2D **BoundaryConditions,
		    BoundValueFunct2D **BoundaryValues,
		    TAuxParam2D *Parameters);

/** assembling for vector finite elements (Raviart-Thomas (RT) and 
 * Brezzi-Douglas-Marini (BDM)) */
void Assemble2D_VectFE(int n_fespaces, TFESpace2D **fespaces,
           int n_sqmatrices, TSquareMatrix2D **sqmatrices,
           int n_matrices, TMatrix2D **matrices,
           int n_rhs, double **rhs, TFESpace2D **ferhs,
           TDiscreteForm2D *DiscreteForm,
           BoundCondFunct2D **BoundaryConditions,
           BoundValueFunct2D **BoundaryValues,
           TAuxParam2D *Parameters=NULL
           );

#ifdef __MORTAR__
void Assemble(TMatrix2D *matrix);

  /** add link term (SDFEM) */
  #ifdef __ADD_LINK_SDFEM__
  void AddLinkSDFEM(TFESpace2D *Space2D, TSquareMatrix2D *Matrix,
  #endif
#endif // __MORTAR__

void Assemble2D(int n_fespaces, TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, TFESpace2D **ferhs,
                TDiscreteForm2D *DiscreteForm,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                TAuxParam2D *Parameters,
		TAuxParam2D *ParametersBound, 
		TypeBoundSwitchFunct2D *TypeBoundSwitcher,
		int *CounterBoundaryParam
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
                );
#ifdef __2D__
void ComputeInterfaceConditionDirichlet(TSquareMatrix2D  **sqmatrices,
					TMatrix2D **matrices,
					TFEVectFunct2D *u_other,
					double *rhs,
					int *interface_dof_Array,
					int *interface_dof_other_Array,
					int interface_dof, double density, 
					double density_other);


void ComputeInterfaceConditionStress(TSquareMatrix2D  **sqmatrices,
				     TMatrix2D **matrices,
				     TFEVectFunct2D *u_other,
				     double *rhs,
				     int *interface_dof_Array,
				     int *interface_dof_other_Array,
				     int interface_dof, double density, 
				     double density_other, 
				     double *rhs_for_stress_other);
     
#endif // __2D__

#endif // __ASSEMBLE2D__
