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
// TimeDiscRout.C       07/12/18
//
// Volker John and Joachim Rang
//
// routines which are used for the different temporal discretizations
//
// ======================================================================
#ifndef __TIMEDISCROUT__
#define __TIMEDISCROUT__

/******************************************************************************/
/*                                                                            */
/* GENERAL                                                                    */
/*                                                                            */
/******************************************************************************/
void SetTimeDiscParameters(int increase_count = 0);

int GetN_SubSteps();

/******************************************************************************/
/*                                                                            */
/* ROSENBROCK                                                                 */
/*                                                                            */
/******************************************************************************/

void AllocateAuxiliaryVectorsRB(int &rb_order, int &RB_s,
				double* &RB_m, double* &RB_ms,
				double* &RB_alpha, double* &RB_gamma,
				double* &RB_sigma, double* &RB_RHS_YN,
				double* &old_sol_rbU, double* &rb_mein,
				double* &rb_diff, double* &sol_tilde,				
				double* &B1, double* &B2, double* &dq,				
				double* RB_A, double* RB_C, double *RB_S,
				int N_Unknowns);

void AssembleRHS_RB_DIRK(TFESpace2D **fesp, TFEFunction2D **fefct, TFESpace2D **ferhs,
TFESpace2D **USpaces, TFEFunction2D **U1Array, TFEFunction2D **U2Array,
TDiscreteForm2D *DiscreteForm, TDiscreteForm2D *DiscreteFormRHS,
double **RHSs, double *rhs, TAuxParam2D *aux,
BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **my_BoundValues,
			 int mg_type, int mg_level, int N_U, int N_Unknowns);


/******************************************************************************/
/*                                                                            */
/* DIRK                                                                       */
/*                                                                            */
/******************************************************************************/

void AllocateAuxiliaryVectorsDIRK(int &rb_order, int &RB_s,
				  int &stiff_acc1, int &stiff_acc2,
				  double* &RB_m, double* &RB_ms,
				  double* &RB_alpha, 
				  double* &RB_RHS_YN,
				  double* &old_sol_rbU, double* &rb_mein,
				  double* &old_sol_rbK, double* &old_rhs,				  
				  double* &rb_diff, double* &sol_tilde,				
				  double* &B1, double* &B2,				
				  double* RB_A,
				  int N_Unknowns);

#endif // __TIMEDISCROUT__
