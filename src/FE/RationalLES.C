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
// RationalLES.C
//
// Purpose:     routines for the rational LES model
//
// Author:      Volker John       start 2005/12/16
//
// =======================================================================

#include <Assemble3D.h>
#include <Constants.h>
#include <Convolution.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>
#include <RationalLES.h>
#include <Solver.h>
#include <TNSE3D_LESParamRout.h>

#include <MainUtilities.h>

#ifdef __3D__
/*******************************************************************************/
//
// ComputeConvolutionOfNabla_uNabla_uTrans3D
//
// result is on duConv->GetValues()
//
/*******************************************************************************/

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
					       int N_Unknowns)
{
    int ii, level_down;
    double *auxConv;
    
    OutPut("ComputeConvolutionOfNabla_uNabla_uTrans3D not tested after removing from main"<<endl);
    level_down =0;
    // compute convolution
    MG->RestrictToAllGrids();
    
    ConvoluteSymmetricTensor3D(UArray[mg_level-1-level_down], duConv);

    auxConv = new double[N_Unknowns];
    // prolongate convolved function
    for (ii=mg_level-1-level_down ; ii< mg_level-1; ii++)
    {
	Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
		   du11ConvArray[ii]->GetValues(),
		   du11ConvArray[ii+1]->GetValues(),
		   auxConv);
	Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
		   du12ConvArray[ii]->GetValues(),
		   du12ConvArray[ii+1]->GetValues(),
		   auxConv);
	Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
		   du13ConvArray[ii]->GetValues(),
		   du13ConvArray[ii+1]->GetValues(),
		   auxConv);
	Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
		   du22ConvArray[ii]->GetValues(),
		   du22ConvArray[ii+1]->GetValues(),
		   auxConv);
	Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
		   du23ConvArray[ii]->GetValues(),
		   du23ConvArray[ii+1]->GetValues(),
		   auxConv);
	Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
		   du33ConvArray[ii]->GetValues(),
		   du33ConvArray[ii+1]->GetValues(),
		   auxConv);
    }
    delete auxConv ;
    return;
}

/*******************************************************************************/
//
//  ComputeConvolutionForTurbVisType4
//
//  result is in uConvArray[mg_level-1]->GetValues()
//
/*******************************************************************************/

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
				       int mg_level, int N_U)
{
    int ii, N_FESpaces, N_Rhs, N_SquareMatrices, N_RectMatrices;
    double *rhsGL00AuxProblem, *RHSs[3], t1, t2, *u_uConv;
    TFESpace3D *fesp[1], *ferhs[3];
    TFEFunction3D *fefct[3];
    BoundCondFunct3D *BoundaryConditionsAuxProblem[3];
    BoundValueFunct3D *BoundValuesAuxProblem[3];
    TAuxParam3D *aux;

    OutPut("ComputeConvolutionForTurbVisType4 not tested after removing from main"<<endl);
    BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;

    BoundValuesAuxProblem[0] = BoundValueAuxProblem;
    BoundValuesAuxProblem[1] = BoundValueAuxProblem;
    BoundValuesAuxProblem[2] = BoundValueAuxProblem;

    // restrict solution to all grids
    // prepare auxiliary problem
    MG->RestrictToAllGrids();
       
    fesp[0] = USpaces[mg_level-1];
    
    fefct[0] = U1Array[mg_level-1];
    fefct[1] = U2Array[mg_level-1];
    fefct[2] = U3Array[mg_level-1];

    aux =  new TAuxParam3D(TimeNSN_FESpacesVeloLES, TimeNSN_FctVeloLES,
			   TimeNSN_ParamFctVeloLES,
			   TimeNSN_FEValuesVeloLES,
			   fesp, fefct,
			   TimeNSFctVeloLES,
			   TimeNSFEFctIndexVeloLES, TimeNSFEMultiIndexVeloLES,
			   TimeNSN_ParamsVeloLES, TimeNSBeginParamVeloLES);

    ferhs[0] = USpaces[mg_level-1];
    ferhs[1] = USpaces[mg_level-1];
    ferhs[2] = USpaces[mg_level-1];    
   
    // assemble rhs
    N_FESpaces = 1;
    N_Rhs = 3;
    N_SquareMatrices = 0;
    N_RectMatrices = 0;

    rhsGL00AuxProblem = new double[3*N_U*SizeOfDouble];
    memset(rhsGL00AuxProblem,0,3*N_U*SizeOfDouble);
    RHSs[0] = rhsGL00AuxProblem;
    RHSs[1] = rhsGL00AuxProblem + N_U;
    RHSs[2] = rhsGL00AuxProblem + 2*N_U;

    Assemble3D(N_FESpaces, fesp,
	       N_SquareMatrices, NULL,
	       N_RectMatrices, NULL,
	       N_Rhs, RHSs, ferhs,
	       DiscreteForm,
	       BoundaryConditionsAuxProblem,
	       BoundValuesAuxProblem,
	       aux);
    delete aux;

    u_uConv = uConvArray[mg_level-1]->GetValues();
    // solve auxiliary problem
    t1 = GetTime();
    // result is on u_uConv
    Solver(sqmatrixGL00AuxProblem, RHSs[0], u_uConv,3);
    t2 = GetTime();
    OutPut( "time for AMG solving: " << t2-t1 << endl);
    
    // restrict computed solution to coarser levels to use it for assembling
    for(ii = mg_level-1 ; ii > 0;ii--)
    {
	RestrictFunction(USpaces[ii-1], USpaces[ii],
			 u1ConvArray[ii-1]->GetValues(),
			 u1ConvArray[ii]->GetValues(),
			 MG->GetLevel(ii-1)->GetAuxVector(0));
	RestrictFunction(USpaces[ii-1], USpaces[ii],
			 u2ConvArray[ii-1]->GetValues(),
			 u2ConvArray[ii]->GetValues(),
			 MG->GetLevel(ii-1)->GetAuxVector(0));
	RestrictFunction(USpaces[ii-1], USpaces[ii],
			 u3ConvArray[ii-1]->GetValues(),
			 u3ConvArray[ii]->GetValues(),
			 MG->GetLevel(ii-1)->GetAuxVector(0));
    }
    delete rhsGL00AuxProblem;
}

/*******************************************************************************/
//
//  PrepareRHSLES 
//
//  result is in LESModelRhs
//
/*******************************************************************************/

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
		   int mg_level, int N_U, int N_P)
{
    int N_FESpaces, N_Rhs;
    double *rhsGL00AuxProblem, *RHSs[6], t1, t2;
    TFESpace3D *fesp[4], *ferhs[6];
    TFEFunction3D *fefct[12];
    TAuxParam3D *aux;
    BoundCondFunct3D *BoundaryConditionsAuxProblem[6];
    BoundValueFunct3D *BoundValuesAuxProblem[6]; 

    // classical LES and rational LES with aux problem somewhat checked
    OutPut("PrepareRHSLES not tested after removing from main"<<endl);
    
    BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[3] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[4] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[5] = BoundConditionAuxProblem;
    BoundValuesAuxProblem[0] = BoundValueAuxProblem;
    BoundValuesAuxProblem[1] = BoundValueAuxProblem;
    BoundValuesAuxProblem[2] = BoundValueAuxProblem;
    BoundValuesAuxProblem[3] = BoundValueAuxProblem;
    BoundValuesAuxProblem[4] = BoundValueAuxProblem;
    BoundValuesAuxProblem[5] = BoundValueAuxProblem;
    
    N_FESpaces = 1;
    N_Rhs = 3;
    RHSs[0] = rhs;
    RHSs[1] = rhs+N_U;
    RHSs[2] = rhs+2*N_U;
    
    switch(TDatabase::ParamDB->DISCTYPE)
    {
	case CLASSICAL_LES :
            fesp[0] = USpaces[mg_level-1];
	    
            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];
            fefct[2] = U3Array[mg_level-1];

            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];
            ferhs[2] = USpaces[mg_level-1];

      if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE<4
        ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE>=100))
            {
		aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVeloLES,
				       TimeNSN_FctVelo_GradVeloLES,
				       TimeNSN_ParamFctVelo_GradVeloLES,
				       TimeNSN_FEValuesVelo_GradVeloLES,
				       fesp, fefct,
				       TimeNSFctVelo_GradVeloLES,
				       TimeNSFEFctIndexVelo_GradVeloLES,
				       TimeNSFEMultiIndexVelo_GradVeloLES,
				       TimeNSN_ParamsVelo_GradVeloLES,
				       TimeNSBeginParamVelo_GradVeloLES);
            }
	    
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 4)
            {
		// convolve velocity
		fesp[1] = uConvSpaces[mg_level-1];
		N_FESpaces++;
		
		fefct[3] = u1ConvArray[mg_level-1];
		fefct[4] = u2ConvArray[mg_level-1];
		fefct[5] = u3ConvArray[mg_level-1];
		
		aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVeloNuT4,
				       TimeNSN_FctVelo_GradVeloNuT4,
				       TimeNSN_ParamFctVelo_GradVeloNuT4,
				       TimeNSN_FEValuesVelo_GradVeloNuT4,
				       fesp, fefct,
				       TimeNSFctVelo_GradVeloNuT4,
				       TimeNSFEFctIndexVelo_GradVeloNuT4,
				       TimeNSFEMultiIndexVelo_GradVeloNuT4,
				       TimeNSN_ParamsVelo_GradVeloNuT4,
				       TimeNSBeginParamVelo_GradVeloNuT4);
            }
	    
            // initialize array
            memset(RHSs[0], 0, (3*N_U+N_P)*SizeOfDouble);
	    //  (3*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);

            Assemble3D(N_FESpaces, fesp,
		       0,NULL,
		       0,NULL,
		       N_Rhs, RHSs, ferhs,
		       DiscreteFormRHSClassicalLES,
		       BoundaryConditions,
		       BoundValues,
		       aux);

	    delete aux;

            memcpy(LESModelRhs, RHSs[0], 3*N_U*SizeOfDouble);
	    break;
  
         case GL00_CONVOLUTION :
            // current solution
            fesp[0] = USpaces[mg_level-1];

            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];
            fefct[2] = U3Array[mg_level-1];

            // current rhs
            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];
            ferhs[2] = USpaces[mg_level-1];

            // convolution of grad u grad u^T
            fesp[1] = duConvSpaces[mg_level-1];
            N_FESpaces++;

            fefct[3] = du11ConvArray[mg_level-1];
            fefct[4] = du12ConvArray[mg_level-1];
            fefct[5] = du13ConvArray[mg_level-1];
            fefct[6] = du22ConvArray[mg_level-1];
            fefct[7] = du23ConvArray[mg_level-1];
            fefct[8] = du33ConvArray[mg_level-1];

            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE<=4
              ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE>=100))
            {
              aux =  new TAuxParam3D(TimeNSN_FESpacesRHSGL00Convolution,
                TimeNSN_FctRHSGL00Convolution,
                TimeNSN_ParamFctRHSGL00Convolution,
                TimeNSN_FEValuesRHSGL00Convolution,
                fesp, fefct,
                TimeNSFctRHSGL00Convolution,
                TimeNSFEFctIndexRHSGL00Convolution,
                TimeNSFEMultiIndexRHSGL00Convolution,
                TimeNSN_ParamsRHSGL00Convolution,
                TimeNSBeginParamRHSGL00Convolution);
            }
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE > 4
              ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE>=100))
            {
              // convoluted velocity
              fesp[2] = uConvSpaces[mg_level-1];
              N_FESpaces++;

              fefct[9] = u1ConvArray[mg_level-1];
              fefct[10] = u2ConvArray[mg_level-1];
              fefct[11] = u3ConvArray[mg_level-1];

              aux =  new TAuxParam3D(TimeNSN_FESpacesRHSGL00ConvolutionNuT4,
                TimeNSN_FctRHSGL00ConvolutionNuT4,
                TimeNSN_ParamFctRHSGL00ConvolutionNuT4,
                TimeNSN_FEValuesRHSGL00ConvolutionNuT4,
                fesp, fefct,
                TimeNSFctRHSGL00ConvolutionNuT4,
                TimeNSFEFctIndexRHSGL00ConvolutionNuT4,
                TimeNSFEMultiIndexRHSGL00ConvolutionNuT4,
                TimeNSN_ParamsRHSGL00ConvolutionNuT4,
                TimeNSBeginParamRHSGL00ConvolutionNuT4);
            }
	   
            Assemble3D(N_FESpaces, fesp,
		       0, NULL,
		       0, NULL,
		       N_Rhs, RHSs, ferhs,
		       DiscreteFormRHSLESModel,
		       BoundaryConditions,
		       BoundValues,
		       aux);
	    delete aux;

            memcpy(LESModelRhs, RHSs[0], 3*N_U*SizeOfDouble);
	    break;
    
	case GL00_AUX_PROBLEM :
	    fesp[0] = USpaces[mg_level-1];
	   
            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];
            fefct[2] = U3Array[mg_level-1];

            aux =  new TAuxParam3D(TimeNSN_FESpacesGradVelo,
              TimeNSN_FctGradVelo,
              TimeNSN_ParamFctGradVelo,
              TimeNSN_FEValuesGradVelo,
              fesp, fefct,
              TimeNSFctGradVelo,
              TimeNSFEFctIndexGradVelo,
              TimeNSFEMultiIndexGradVelo,
              TimeNSN_ParamsGradVelo,
              TimeNSBeginParamGradVelo);


            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];
            ferhs[2] = USpaces[mg_level-1];
            ferhs[3] = USpaces[mg_level-1];
            ferhs[4] = USpaces[mg_level-1];
            ferhs[5] = USpaces[mg_level-1];

            N_FESpaces = 1;
            N_Rhs = 6;
	       
	    rhsGL00AuxProblem = new double[6*N_U*SizeOfDouble];	    
            memset(rhsGL00AuxProblem,0,6*N_U*SizeOfDouble);
            RHSs[0] = rhsGL00AuxProblem;
            RHSs[1] = rhsGL00AuxProblem+ N_U;
            RHSs[2] = rhsGL00AuxProblem+ 2*N_U;
            RHSs[3] = rhsGL00AuxProblem+ 3*N_U;
            RHSs[4] = rhsGL00AuxProblem+ 4*N_U;
            RHSs[5] = rhsGL00AuxProblem+ 5*N_U;

	    Assemble3D(N_FESpaces, fesp,
		       0, NULL,
		       0, NULL,
		       N_Rhs, RHSs, ferhs,
		       DiscreteFormGL00AuxProblemRHS,
		       BoundaryConditionsAuxProblem,
		       BoundValuesAuxProblem,
		       aux);
	    delete aux;

           // solve auxiliary problem
            t1 = GetTime();
	    // solution on solGL00AuxProblem
            Solver(sqmatrixGL00AuxProblem, RHSs[0], solGL00AuxProblem, 6);
            t2 = GetTime();
            OutPut( "time for AMG solving: " << t2-t1 << endl);
	    delete rhsGL00AuxProblem;
	    
            // assemble rhs for Galdi/Layton model
            fefct[3] = GL00AuxProblemSol11Array[mg_level-1];
            fefct[4] = GL00AuxProblemSol12Array[mg_level-1];
            fefct[5] = GL00AuxProblemSol13Array[mg_level-1];
            fefct[6] = GL00AuxProblemSol22Array[mg_level-1];
            fefct[7] = GL00AuxProblemSol23Array[mg_level-1];
            fefct[8] = GL00AuxProblemSol33Array[mg_level-1];

            N_Rhs = 3;
          
            RHSs[0] = rhs;
            RHSs[1] = rhs+N_U;
            RHSs[2] = rhs+2*N_U;

	    if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE<=4
        ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE>=100))
            {
              aux =  new TAuxParam3D(
                TimeNSN_FESpacesGL00AuxProblem,
                TimeNSN_FctGL00AuxProblem,
                TimeNSN_ParamFctGL00AuxProblem,
                TimeNSN_FEValuesGL00AuxProblem,
                fesp, fefct,
                TimeNSFctGL00AuxProblem,
                TimeNSFEFctIndexGL00AuxProblem,
                TimeNSFEMultiIndexGL00AuxProblem,
                TimeNSN_ParamsGL00AuxProblem,
                TimeNSBeginParamGL00AuxProblem);
            }

            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE > 4
              ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE>=100))
            {
              // convoluted velocity
              fesp[1] = uConvSpaces[mg_level-1];
              N_FESpaces++;

              fefct[9] = u1ConvArray[mg_level-1];
              fefct[10] = u2ConvArray[mg_level-1];
              fefct[11] = u3ConvArray[mg_level-1];

              aux =  new TAuxParam3D(TimeNSN_FESpacesGL00AuxProblemNuT4,
                TimeNSN_FctGL00AuxProblemNuT4,
                TimeNSN_ParamFctGL00AuxProblemNuT4,
                TimeNSN_FEValuesGL00AuxProblemNuT4,
                fesp, fefct,
                TimeNSFctGL00AuxProblemNuT4,
                TimeNSFEFctIndexGL00AuxProblemNuT4,
                TimeNSFEMultiIndexGL00AuxProblemNuT4,
                TimeNSN_ParamsGL00AuxProblemNuT4,
                TimeNSBeginParamGL00AuxProblemNuT4);
            }

            Assemble3D(N_FESpaces, fesp,
		       0, NULL,
		       0, NULL,
		       N_Rhs, RHSs, ferhs,
		       DiscreteFormRHSLESModel,
		       BoundaryConditions,
		       BoundValues,
		       aux);
	    delete aux;

            memcpy(LESModelRhs, RHSs[0], 3*N_U*SizeOfDouble);
	    break;
    }
    return;
}
/*******************************************************************************/
//
//  ConvolveSolution
//
//  result is in u_uConv
//
/*******************************************************************************/

void ConvolveSolution(TNSE_MultiGrid *MG,
		      TFESpace3D **USpaces,
		      TFEFunction3D **U1Array,
		      TFEFunction3D **U2Array,
		      TFEFunction3D **U3Array,
		      TDiscreteForm3D *DiscreteForm,
		      TSquareMatrix3D *sqmatrixGL00AuxProblem,
		      double *rhsGL00AuxProblem,
		      double *u_uConv,
		      int mg_level, int N_U)
{
    int N_FESpaces, N_Rhs;
    double *RHSs[3];
    TFESpace3D *fesp[1], *ferhs[3];
    TFEFunction3D *fefct[3];
    TAuxParam3D *aux;
    BoundCondFunct3D *BoundaryConditionsAuxProblem[3];
    BoundValueFunct3D *BoundValuesAuxProblem[3];

    OutPut("ConvolveSolution not tested after removing from main"<<endl);
    BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;

    BoundValuesAuxProblem[0] = BoundValueAuxProblem;
    BoundValuesAuxProblem[1] = BoundValueAuxProblem;
    BoundValuesAuxProblem[2] = BoundValueAuxProblem;

      // prepare auxiliary problem
    MG->RestrictToAllGrids();
    
    fesp[0] = USpaces[mg_level-1];

    fefct[0] = U1Array[mg_level-1];
    fefct[1] = U2Array[mg_level-1];
    fefct[2] = U3Array[mg_level-1];
    
    ferhs[0] = USpaces[mg_level-1];
    ferhs[1] = USpaces[mg_level-1];
    ferhs[2] = USpaces[mg_level-1];

    // assemble rhs
    N_FESpaces = 1;
    N_Rhs = 3;
    
    memset(rhsGL00AuxProblem,0,3*N_U*SizeOfDouble);
    RHSs[0] = rhsGL00AuxProblem;
    RHSs[1] = rhsGL00AuxProblem + N_U;
    RHSs[2] = rhsGL00AuxProblem + 2*N_U;
    
    aux =  new TAuxParam3D(TimeNSN_FESpacesVeloLES, TimeNSN_FctVeloLES,
			   TimeNSN_ParamFctVeloLES,
			   TimeNSN_FEValuesVeloLES,
			   fesp, fefct,
			   TimeNSFctVeloLES,
			   TimeNSFEFctIndexVeloLES, TimeNSFEMultiIndexVeloLES,
			   TimeNSN_ParamsVeloLES, TimeNSBeginParamVeloLES);
    
    Assemble3D(N_FESpaces, fesp,
	       0, NULL,
	       0, NULL,
	       N_Rhs, RHSs, ferhs,
	       DiscreteForm,
	       BoundaryConditionsAuxProblem,
	       BoundValuesAuxProblem,
	       aux);
    delete aux;

    // solve auxiliary problem
    Solver(sqmatrixGL00AuxProblem, RHSs[0], u_uConv,3);
 
    return;
}

// ========================================================================
// applies differential filter to the velocity
// ========================================================================

void ApplyDifferentialFilterToVelocity(TFESpace3D **USpaces,TFEVectFunct3D **UArray,
       TSquareMatrix3D *sqmatrixGL00AuxProblem,
       TDiscreteForm3D *DiscreteFormRHSAuxProblemU,
       double *solGL00AuxProblem, BoundCondFunct3D **BoundaryConditionsAuxProblem, 
       BoundValueFunct3D **BoundValuesAuxProblem,
               int mg_level)

{
  int N_U;
  double *rhsGL00AuxProblem, *RHSs[3], t1, t2;
  TFESpace3D *fesp[1], *ferhs[3];
  TFEFunction3D *fefct[3];
  TAuxParam3D *aux;
    //BoundCondFunct3D *BoundaryConditionsAuxProblem[3];
    //BoundValueFunct3D *BoundValuesAuxProblem[3]; 
  
  /*   BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
    BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;
    BoundValuesAuxProblem[0] = BoundValueAuxProblem;
    BoundValuesAuxProblem[1] = BoundValueAuxProblem;
    BoundValuesAuxProblem[2] = BoundValueAuxProblem;*/
 
    // prepare assembling
  fesp[0] = USpaces[mg_level-1];
     
  fefct[0] = UArray[mg_level-1]->GetComponent(0);
  fefct[1] = UArray[mg_level-1]->GetComponent(1);
  fefct[2] = UArray[mg_level-1]->GetComponent(2);

          aux =  new TAuxParam3D(TimeNSN_FESpacesVelo,
            TimeNSN_FctVelo,
            TimeNSN_ParamFctVelo,
            TimeNSN_FEValuesVelo,
            fesp, fefct,
            TimeNSFctVelo,
            TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
            TimeNSN_ParamsVelo, TimeNSBeginParamVelo);

            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];
            ferhs[2] = USpaces[mg_level-1];
    
      N_U = fesp[0]->GetN_DegreesOfFreedom();
      rhsGL00AuxProblem = new double[3*N_U*SizeOfDouble];     
            memset(rhsGL00AuxProblem,0,3*N_U*SizeOfDouble);
            RHSs[0] = rhsGL00AuxProblem;
            RHSs[1] = rhsGL00AuxProblem+ N_U;
            RHSs[2] = rhsGL00AuxProblem+ 2*N_U;

      Assemble3D(1, fesp,
           0, NULL,
           0, NULL,
           3, RHSs, ferhs,
           DiscreteFormRHSAuxProblemU,
           BoundaryConditionsAuxProblem,
           BoundValuesAuxProblem,
           aux);
      delete aux;

           // solve auxiliary problem
            t1 = GetTime();
      // solution on solGL00AuxProblem
            Solver(sqmatrixGL00AuxProblem, RHSs[0], solGL00AuxProblem, 3);
            t2 = GetTime();
            OutPut( "differential filter: time for AMG solving: " << t2-t1 << endl);
      delete rhsGL00AuxProblem;
}

// ========================================================================
// boundary values for auxiliary problem in Galdi/Layton model
// ========================================================================

void BoundConditionAuxProblem(int CompID, double x, double y, double z,
                              BoundCond &cond)
{
  double eps = 1e-6;
 
  /** for Channel180.h */
  /*  if ((fabs(z) < eps) || (fabs(z-2.0) < eps))
    cond = DIRICHLET;
  else
    cond = NEUMANN;*/
  /** for  CylinderSquare22000.h */
  
  cond = NEUMANN;
  cond = DIRICHLET;
 
 
    // inflow boundary condition
    if (fabs(x)<1e-6)
      cond  = DIRICHLET;

  if ((fabs(x-0.45)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
  {
      cond  = DIRICHLET;
      //cout << "left ";
  }
  if ((fabs(x-0.55)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
  {
      cond  = DIRICHLET;
      //cout << "right ";
  }
  if ((fabs(y-0.65)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
  {
      cond  = DIRICHLET;
      //cout << "lower ";
  }
  if ((fabs(y-0.75)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
  {
      cond  = DIRICHLET;
      //cout << "upper ";
  }
    
  if (fabs(x-2.5)<1e-6)
    cond = NEUMANN;
}

void BoundValueAuxProblem(int CompID, double x, double y, double z, double &value)
{
  /** for Channel180.h */
  // value = 0;
  /** for  CylinderSquare22000.h */
  value = 0;
}

void BoundValueAuxProblemU1(int CompID,double x, double y, double z, double &value)
{
  double val[4];
  /** for Channel180.h */
  // value = 0;
  /** for  CylinderSquare22000.h */
  value = 1.5;
   if (fabs(x)<1e-6)
      value = 1;

  if ((fabs(x-0.45)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
  {
     value = 2;
      //cout << "left ";
  }
  if ((fabs(x-0.55)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
  {
      value = 3;
      //cout << "right ";
  }
  if ((fabs(y-0.65)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
  {
      value = 4;
      //cout << "lower ";
  }
  if ((fabs(y-0.75)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
  {
       value = 5; 
      //cout << "upper ";
  }
   
  if (fabs(x-2.5)<1e-6)
    value = 0;  
  OutPut(x << " " << y << " " << z << " " << value << " " << endl);
 }


#endif
