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
   
#ifdef __2D__
#endif

#ifdef __3D__
void ComputeConvolutionOfNabla_uNabla_uTrans3D(TNSE_MultiGrid *MG, 
					       TFEVectFunct3D **UArray,
					       TFEVectFunct3D *duConv,
					       TFESpace3D **duConvSpaces,
					       TFEFunction3D **du11ConvArray,
					       TFEFunction3D **du12ConvArray, 
					       TFEFunction3D **du13ConvArray,
					       TFEFunction3D **du22ConvArray, 
					       TFEFunction3D **du23ConvArray, 
					       TFEFunction3D **du33ConvArray,
					       int mg_level,
					       int N_Unknowns);

void ComputeConvolutionForTurbVisType4(TNSE_MultiGrid *MG, 
				       TFESpace3D **USpaces,
				       TFEFunction3D **U1Array,
				       TFEFunction3D **U2Array,
				       TFEFunction3D **U3Array,
				       TFEVectFunct3D **uConvArray,
				       TFEFunction3D **u1ConvArray,
				       TFEFunction3D **u2ConvArray,
				       TFEFunction3D **u3ConvArray,
				       TDiscreteForm3D *DiscreteForm,
				       TSquareMatrix3D *sqmatrixGL00AuxProblem,
				       int mg_level, int N_U);

void PrepareRHSLES(TFESpace3D **USpaces,
		   TFEFunction3D **U1Array,
		   TFEFunction3D **U2Array,
		   TFEFunction3D **U3Array,
		   TFESpace3D **uConvSpaces,
		   TFEFunction3D **u1ConvArray,
		   TFEFunction3D **u2ConvArray,
		   TFEFunction3D **u3ConvArray,
		   TFESpace3D **duConvSpaces,
		   TFEFunction3D **du11ConvArray,
		   TFEFunction3D **du12ConvArray,
		   TFEFunction3D **du13ConvArray,
		   TFEFunction3D **du22ConvArray,
		   TFEFunction3D **du23ConvArray,
		   TFEFunction3D **du33ConvArray,
		   TFEFunction3D **GL00AuxProblemSol11Array,
		   TFEFunction3D **GL00AuxProblemSol12Array,
		   TFEFunction3D **GL00AuxProblemSol13Array,
		   TFEFunction3D **GL00AuxProblemSol22Array,
		   TFEFunction3D **GL00AuxProblemSol23Array,
		   TFEFunction3D **GL00AuxProblemSol33Array,
		   TDiscreteForm3D *DiscreteFormRHSClassicalLES,
		   TDiscreteForm3D *DiscreteFormRHSLESModel,
		   TDiscreteForm3D *DiscreteFormGL00AuxProblemRHS,
		   BoundCondFunct3D **BoundaryConditions,
		   BoundValueFunct3D **BoundValues,
		   TSquareMatrix3D *sqmatrixGL00AuxProblem,
		   double *rhs,
		   double *solGL00AuxProblem,
		   double *LESModelRhs,
		   int mg_level, int N_U, int N_P);

void ConvolveSolution(TNSE_MultiGrid *MG,
		      TFESpace3D **USpaces,
		      TFEFunction3D **U1Array,
		      TFEFunction3D **U2Array,
		      TFEFunction3D **U3Array,
		      TDiscreteForm3D *DiscreteForm,
		      TSquareMatrix3D *sqmatrixGL00AuxProblem,
		      double *rhsGL00AuxProblem,
		      double *u_uConv,
		      int mg_level, int N_U);

void ApplyDifferentialFilterToVelocity(TFESpace3D **USpaces,
                                       TFEVectFunct3D **UArray,
                                       TSquareMatrix3D *sqmatrixGL00AuxProblem,
                                       TDiscreteForm3D *DiscreteFormRHSAuxProblemU,
                                       double *solGL00AuxProblem,
                                       BoundCondFunct3D **BoundaryConditionsAuxProblem, 
                                       BoundValueFunct3D **BoundValuesAuxProblem,
                                       int mg_level);


void BoundConditionAuxProblem(int CompID, double x, double y, double z,
                              BoundCond &cond);

void BoundValueAuxProblem(int CompID, double x, double y, double z,
                          double &value);
void BoundValueAuxProblemU1(int CompID, double x, double y, double z,
                          double &value);


#endif // __3D__
