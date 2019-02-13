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
// @(#)TNSE2D_ParamRout.h        1.2 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#ifndef __TNSE2D_PARAMROUT__
#define __TNSE2D_PARAMROUT__

// ======================================================================
// setting for error calculation for all types
// ======================================================================
MultiIndex2D TimeNSAllDerivatives[3] = { D00, D10, D01 };

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================

void TimeNSParams2(double *in, double *out);

int TimeNSN_FESpaces2 = 1;
int TimeNSN_Fct2 = 2;
int TimeNSN_ParamFct2 = 1;
int TimeNSN_FEValues2 = 2;
int TimeNSN_Params2 = 2;
int TimeNSFEFctIndex2[2] = { 0, 1 };
MultiIndex2D TimeNSFEMultiIndex2[2] = { D00, D00 };
ParamFct *TimeNSFct2[1] = { TimeNSParams2 };
int TimeNSBeginParam2[1] = { 0 };



// ========================================================================
// Rosenbrock Methods
// ========================================================================

/*
void TimeNSParams2RB(double *in, double *out);

int TimeNSN_FESpaces2 = 1;
int TimeNSN_Fct2 = 2;
int TimeNSN_ParamFct2 = 1;
int TimeNSN_FEValues2 = 6;
int TimeNSN_Params2 = 6;
int TimeNSFEFctIndex2[6] = { 0, 1, 0, 1, 0, 1 };
MultiIndex2D TimeNSFEMultiIndex2[6] = { D00, D00, D10, D10, D01, D01 };
ParamFct *TimeNSFct2[1] = { TimeNSParams2RB };
int TimeNSBeginParam2[1] = { 0 };
 */

// ========================================================================
// coletti, without g_\delta \ast u
// ========================================================================

void TimeNSParamsVelo_GradVelo(double *in, double *out);
int TimeNSN_FESpacesVelo_GradVelo = 1;
int TimeNSN_FctVelo_GradVelo = 2;
int TimeNSN_ParamFctVelo_GradVelo = 1;
int TimeNSN_FEValuesVelo_GradVelo = 6;
int TimeNSN_ParamsVelo_GradVelo = 6;
int TimeNSFEFctIndexVelo_GradVelo[6] = { 0, 1, 0, 1, 0, 1 };
MultiIndex2D TimeNSFEMultiIndexVelo_GradVelo[6] = { D00, D00, D10, D10, D01, D01 };
ParamFct *TimeNSFctVelo_GradVelo[1] = { TimeNSParamsVelo_GradVelo };
int TimeNSBeginParamVelo_GradVelo[1] = { 0 };

// ========================================================================
// coletti, without g_\delta \ast u, in ALE 
// ========================================================================
void TimeNSParamsVelo_GradVelo_ALE(double *in, double *out);
int TimeNSN_FESpacesVelo_GradVelo_ALE = 2;
int TimeNSN_FctVelo_GradVelo_ALE = 4;
int TimeNSN_ParamFctVelo_GradVelo_ALE  = 1;
int TimeNSN_FEValuesVelo_GradVelo_ALE  = 8;
int TimeNSN_ParamsVelo_GradVelo_ALE  = 8;
int TimeNSFEFctIndexVelo_GradVelo_ALE[8] = { 0, 1, 0, 1, 0, 1, 2, 3 };
MultiIndex2D TimeNSFEMultiIndexVelo_GradVelo_ALE[8] = { D00, D00, D10, D10, D01, D01, D00, D00};
ParamFct *TimeNSFctVelo_GradVelo_ALE[1] = { TimeNSParamsVelo_GradVelo_ALE };
int TimeNSBeginParamVelo_GradVelo_ALE[1] = { 0 };

// ========================================================================
// velocity, gradient and convolution of velocity
// ========================================================================

void TimeNSParamsVelo_GradVelo_ConvVelo(double *in, double *out);
int TimeNSN_FESpacesVelo_GradVelo_ConvVelo = 2;
int TimeNSN_FctVelo_GradVelo_ConvVelo = 4;
int TimeNSN_ParamFctVelo_GradVelo_ConvVelo = 1;
int TimeNSN_FEValuesVelo_GradVelo_ConvVelo = 8;
int TimeNSN_ParamsVelo_GradVelo_ConvVelo = 8;
int TimeNSFEFctIndexVelo_GradVelo_ConvVelo[8] = { 0, 1, 0, 1, 0, 1, 2, 3 };
MultiIndex2D TimeNSFEMultiIndexVelo_GradVelo_ConvVelo[8] = { D00, D00, D10, 
                                                             D10, D01, D01,
                                                             D00, D00 };
ParamFct *TimeNSFctVelo_GradVelo_ConvVelo[1] = { TimeNSParamsVelo_GradVelo_ConvVelo };
int TimeNSBeginParamVelo_GradVelo_ConvVelo[1] = { 0 };

// ========================================================================
// Coletti,
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================
void TimeNSParamsVelo_GradVeloNuT4(double *in, double *out);
int TimeNSN_FESpacesVelo_GradVeloNuT4 = 2;
int TimeNSN_FctVelo_GradVeloNuT4 = 4;
int TimeNSN_ParamFctVelo_GradVeloNuT4 = 1;
int TimeNSN_FEValuesVelo_GradVeloNuT4 = 8;
int TimeNSN_ParamsVelo_GradVeloNuT4 = 8;
int TimeNSFEFctIndexVelo_GradVeloNuT4[8] = { 0, 1, 0, 1, 0, 1, 2, 3 };
MultiIndex2D TimeNSFEMultiIndexVelo_GradVeloNuT4[8] = { D00, D00, D10, D10, D01, D01, D00, D00 };
ParamFct *TimeNSFctVelo_GradVeloNuT4[1] = { TimeNSParamsVelo_GradVeloNuT4 };
int TimeNSBeginParamVelo_GradVeloNuT4[1] = { 0 };

// ========================================================================
// Galdi/Layton with convolution, without g_\delta \ast u 
// ========================================================================
void TimeNSParamsGL00Convolution(double *in, double *out);
int TimeNSN_FESpacesGL00Convolution = 2;
int TimeNSN_FctGL00Convolution = 5;
int TimeNSN_ParamFctGL00Convolution = 1;
int TimeNSN_FEValuesGL00Convolution = 9;
int TimeNSN_ParamsGL00Convolution = 9;
int TimeNSFEFctIndexGL00Convolution[9] = { 0, 1, 0, 1, 0, 1, 2, 3, 4 };
MultiIndex2D TimeNSFEMultiIndexGL00Convolution[9] = { D00, D00, D10, D10, D01, D01,
                                             D00, D00, D00};
ParamFct *TimeNSFctGL00Convolution[1] = { TimeNSParamsGL00Convolution };
int TimeNSBeginParamGL00Convolution[1] = { 0 };

// ========================================================================
// Galdi/Layton with convolution, rhs assembling, without g_\delta \ast u
// ========================================================================
/*void TimeNSParamsRHSGL00Convolution(double *in, double *out);
int TimeNSN_FESpacesRHSGL00Convolution = 2;
int TimeNSN_FctRHSGL00Convolution = 5;
int TimeNSN_ParamFctRHSGL00Convolution = 1;
int TimeNSN_FEValuesRHSGL00Convolution = 9;
int TimeNSN_ParamsRHSGL00Convolution = 7;
int TimeNSFEFctIndexRHSGL00Convolution[9] = { 0, 1, 0, 1, 0, 1, 2, 3, 4 };
MultiIndex2D TimeNSFEMultiIndexRHSGL00Convolution[9] = 
{ D00, D00, D10, D10, D01, D01, D00, D00, D00};
ParamFct *TimeNSFctRHSGL00Convolution[1] = { TimeNSParamsRHSGL00Convolution };
int TimeNSBeginParamRHSGL00Convolution[1] = { 0 }; */

void TimeNSParamsRHSLES(double *in, double *out);
int TimeNSN_FESpacesRHSGL00Convolution = 2;
int TimeNSN_FctRHSGL00Convolution = 5;
int TimeNSN_ParamFctRHSGL00Convolution = 1;
int TimeNSN_FEValuesRHSGL00Convolution = 3;
int TimeNSN_ParamsRHSGL00Convolution = 3;
int TimeNSFEFctIndexRHSGL00Convolution[3] = { 2, 3, 4 };
MultiIndex2D TimeNSFEMultiIndexRHSGL00Convolution[3] = 
{ D00, D00, D00};
ParamFct *TimeNSFctRHSGL00Convolution[1] = { TimeNSParamsRHSLES };
int TimeNSBeginParamRHSGL00Convolution[1] = { 0 };

// ========================================================================
// Galdi/Layton with convolution, rhs assembling, 
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================
void TimeNSParamsRHSGL00ConvolutionNuT4(double *in, double *out);
int TimeNSN_FESpacesRHSGL00ConvolutionNuT4 = 3;
int TimeNSN_FctRHSGL00ConvolutionNuT4 = 7;
int TimeNSN_ParamFctRHSGL00ConvolutionNuT4 = 1;
int TimeNSN_FEValuesRHSGL00ConvolutionNuT4 = 11;
int TimeNSN_ParamsRHSGL00ConvolutionNuT4 = 11;
int TimeNSFEFctIndexRHSGL00ConvolutionNuT4[11] = { 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6 };
MultiIndex2D TimeNSFEMultiIndexRHSGL00ConvolutionNuT4[11] = 
{ D00, D00, D10, D10, D01, D01, D00, D00, D00, D00, D00};
ParamFct *TimeNSFctRHSGL00ConvolutionNuT4[1] = { TimeNSParamsRHSGL00ConvolutionNuT4 };
int TimeNSBeginParamRHSGL00ConvolutionNuT4[1] = { 0 };

// ========================================================================
// parameters for Galdi/Layton model with auxiliary problem
// ========================================================================

int TimeNSN_FESpacesGL00AuxProblem = 1;
int TimeNSN_FctGL00AuxProblem = 5;
int TimeNSN_ParamFctGL00AuxProblem = 1;
int TimeNSN_FEValuesGL00AuxProblem = 3;
int TimeNSN_ParamsGL00AuxProblem = 3;
int TimeNSFEFctIndexGL00AuxProblem[3] = { 2, 3, 4};
MultiIndex2D TimeNSFEMultiIndexGL00AuxProblem[3] = { D00, D00, D00};
ParamFct *TimeNSFctGL00AuxProblem[1] = { TimeNSParamsRHSLES };
int TimeNSBeginParamGL00AuxProblem[1] = { 0 };

// ========================================================================
// parameters for Galdi/Layton model with auxiliary problem
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================

void TimeNSParamsGL00AuxProblemNuT4(double *in, double *out);
int TimeNSN_FESpacesGL00AuxProblemNuT4 = 2;
int TimeNSN_FctGL00AuxProblemNuT4 = 7;
int TimeNSN_ParamFctGL00AuxProblemNuT4 = 1;
int TimeNSN_FEValuesGL00AuxProblemNuT4 = 11;
int TimeNSN_ParamsGL00AuxProblemNuT4 = 11;
int TimeNSFEFctIndexGL00AuxProblemNuT4[11] = { 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6 };
MultiIndex2D TimeNSFEMultiIndexGL00AuxProblemNuT4[11] = { D00, D00, D10, D10, D01, D01,
                                             D00, D00, D00, D00, D00 };
ParamFct *TimeNSFctGL00AuxProblemNuT4[1] = { TimeNSParamsGL00AuxProblemNuT4 };
int TimeNSBeginParamGL00AuxProblemNuT4[1] = { 0 };

void TimeNSParamsGL00AuxProblemPaper2(double *in, double *out);
int TimeNSN_FESpacesGL00AuxProblemPaper2 = 1;
int TimeNSN_FctGL00AuxProblemPaper2 = 8;
int TimeNSN_ParamFctGL00AuxProblemPaper2 = 1;
int TimeNSN_FEValuesGL00AuxProblemPaper2 = 13;
int TimeNSN_ParamsGL00AuxProblemPaper2 = 11;
int TimeNSFEFctIndexGL00AuxProblemPaper2[13] = { 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6, 6, 7};
MultiIndex2D TimeNSFEMultiIndexGL00AuxProblemPaper2[13] = { D00, D00, D10, D10, D01, D01,
                                             D00, D00, D00, D10, D10, D01, D01};
ParamFct *TimeNSFctGL00AuxProblemPaper2[1] = { TimeNSParamsGL00AuxProblemPaper2 };
int TimeNSBeginParamGL00AuxProblemPaper2[1] = { 0 };

// ========================================================================
// parameters: gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsGrad(double *in, double *out);

int TimeNSN_FESpacesGrad = 1;
int TimeNSN_FctGrad = 2;
int TimeNSN_ParamFctGrad = 1;
int TimeNSN_FEValuesGrad = 6;
int TimeNSN_ParamsGrad = 4;
int TimeNSFEFctIndexGrad[6] = { 0, 1, 0, 1, 0, 1 };
MultiIndex2D TimeNSFEMultiIndexGrad[6] = { D00, D00, D10, D10, D01, D01 };
ParamFct *TimeNSFctGrad[1] = { TimeNSParamsGrad };
int TimeNSBeginParamGrad[1] = { 0 };

// ========================================================================
// parameters for VMS
// ========================================================================
void TimeNSParams_VMS_SmallRhs2D(double *in, double *out);

int TimeNSN_FESpaces_VMS_SmallRhs = 3;
int TimeNSN_Fct_VMS_SmallRhs  = 5;
int TimeNSN_ParamFct_VMS_SmallRhs  = 1;
int TimeNSN_FEValues_VMS_SmallRhs = 13;
int TimeNSN_Params_VMS_SmallRhs = 13;
int TimeNSFEFctIndex_VMS_SmallRhs[13] = { 0, 1, 0, 1, 0, 1, 2, 3, 
                                          2, 3, 2, 3, 4};
MultiIndex2D TimeNSFEMultiIndex_VMS_SmallRhs[13] = { D00, D00,
                                                     D10, D10,
                                                     D01, D01,
                                                     D00, D00,
                                                     D10, D10,
                                                     D01, D01,
                                                     D00 };

ParamFct *TimeNSFct_VMS_SmallRhs[1] = {  TimeNSParams_VMS_SmallRhs2D };
int TimeNSBeginParam_VMS_SmallRhs[1] = { 0 };

void TimeNSParams_VMS_LargeRhs2D(double *in, double *out);

int TimeNSN_FESpaces_VMS_LargeRhs = 3;
int TimeNSN_Fct_VMS_LargeRhs  = 5;
int TimeNSN_ParamFct_VMS_LargeRhs  = 1;
int TimeNSN_FEValues_VMS_LargeRhs = 13;
int TimeNSN_Params_VMS_LargeRhs = 13;
int TimeNSFEFctIndex_VMS_LargeRhs[13] = { 0, 1, 0, 1, 0, 1, 2, 3, 
                                          2, 3, 2, 3, 4};
MultiIndex2D TimeNSFEMultiIndex_VMS_LargeRhs[13] = { D00, D00,
                                                     D10, D10,
                                                     D01, D01,
                                                     D00, D00,
                                                     D10, D10,
                                                     D01, D01,
                                                     D00 };

ParamFct *TimeNSFct_VMS_LargeRhs[1] = {  TimeNSParams_VMS_LargeRhs2D };
int TimeNSBeginParam_VMS_LargeRhs[1] = { 0 };

// ========================================================================
// parameters: low order : u1old, u2old, higher order  : u1old, u2old
// ========================================================================
void TimeNSParams_NLGalerkin_VMS_1_2D(double *in, double *out);

int TimeNSN_FESpaces_NLGalerkin_VMS_1 = 2;
int TimeNSN_Fct_NLGalerkin_VMS_1 = 4;
int TimeNSN_ParamFct_NLGalerkin_VMS_1 = 1;
int TimeNSN_FEValues_NLGalerkin_VMS_1 = 4;
int TimeNSN_Params_NLGalerkin_VMS_1 = 4;
int TimeNSFEFctIndex_NLGalerkin_VMS_1[4] = { 0, 1, 2, 3 };
MultiIndex2D TimeNSFEMultiIndex_NLGalerkin_VMS_1[4] = { D00, D00, D00, D00 };
ParamFct *TimeNSFct_NLGalerkin_VMS_1[1] = { TimeNSParams_NLGalerkin_VMS_1_2D };
int TimeNSBeginParam_NLGalerkin_VMS_1[1] = { 0 };

// ========================================================================
// right-hand side ONLY, defect correction type 1, u2
// ========================================================================

void TimeNSParamsVelo_GradVeloOld2(double *in, double *out);
int TimeNSN_FESpacesVelo_GradVeloOld2 = 1;
int TimeNSN_FctVelo_GradVeloOld2 = 6;
int TimeNSN_ParamFctVelo_GradVeloOld2 = 1;
int TimeNSN_FEValuesVelo_GradVeloOld2 = 10;
int TimeNSN_ParamsVelo_GradVeloOld2 = 10;
int TimeNSFEFctIndexVelo_GradVeloOld2[10] = { 0, 1, 0, 1, 0, 1, 2, 3, 4, 5 };
MultiIndex2D TimeNSFEMultiIndexVelo_GradVeloOld2[10] = { D00, D00, D10, D10, D01, D01,
                                                        D00, D00, D00, D00};
ParamFct *TimeNSFctVelo_GradVeloOld2[1] = { TimeNSParamsVelo_GradVeloOld2 };
int TimeNSBeginParamVelo_GradVeloOld2[1] = { 0 };

// ========================================================================
// parameters: x, y, u1old, u2old
// ========================================================================

void TimeNSParamsVeloPos(double *in, double *out);

int TimeNSN_FESpacesVeloPos = 1;
int TimeNSN_FctVeloPos = 2;
int TimeNSN_ParamFctVeloPos = 1;
int TimeNSN_FEValuesVeloPos = 2;
int TimeNSN_ParamsVeloPos = 4;
int TimeNSFEFctIndexVeloPos[2] = { 0, 1 };
MultiIndex2D TimeNSFEMultiIndexVeloPos[2] = { D00, D00 };
ParamFct *TimeNSFctVeloPos[1] = { TimeNSParamsVeloPos };
int TimeNSBeginParamVeloPos[1] = { 0 };

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams(double *in, double *out);

int MovingTNSN_FESpaces = 2;
int MovingTNSN_Fct = 4;
int MovingTNSN_ParamFct = 1;
int MovingTNSN_FEValues = 4;
int MovingTNSN_Params = 3;
int MovingTNSFEFctIndex[4] = { 0, 1, 2, 3 };
MultiIndex2D MovingTNSFEMultiIndex[4] = { D00, D00, D00, D00 };
ParamFct *MovingTNSFct[1] = { MovingTNSParams };
int MovingTNSBeginParam[1] = { 0 };



// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
//  for axial symmetric case
// ======================================================================
void MovingTNSParams_Axial3D(double *in, double *out);

int MovingTNSN_FESpaces_Axial3D = 2;
int MovingTNSN_Fct_Axial3D = 4;
int MovingTNSN_ParamFct_Axial3D = 1;
int MovingTNSN_FEValues_Axial3D = 4;
int MovingTNSN_Params_Axial3D = 4;
int MovingTNSFEFctIndex_Axial3D[4] = { 0, 1, 2, 3 };
MultiIndex2D MovingTNSFEMultiIndex_Axial3D[4] = { D00, D00, D00, D00 };
ParamFct *MovingTNSFct_Axial3D[1] = { MovingTNSParams_Axial3D };
int MovingTNSBeginParam_Axial3D[1] = { 0 };


// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
//  for axial symmetric case
// ======================================================================
void MovingTNSParams_Axial3D_HeatLine(double *in, double *out);

int MovingTNSN_FESpaces_Axial3D_HeatLine = 3;
int MovingTNSN_Fct_Axial3D_HeatLine = 7;
int MovingTNSN_ParamFct_Axial3D_HeatLine = 1;
int MovingTNSN_FEValues_Axial3D_HeatLine = 11;
int MovingTNSN_Params_Axial3D_HeatLine = 8;
int MovingTNSFEFctIndex_Axial3D_HeatLine[11] = { 0, 0, 0, 1, 2, 3, 4, 1, 2, 3, 4};
MultiIndex2D MovingTNSFEMultiIndex_Axial3D_HeatLine[11] = { D00, D10, D01, D00, D00, D00, D00, D01, D10, D01, D10};
ParamFct *MovingTNSFct_Axial3D_HeatLine[1] = { MovingTNSParams_Axial3D_HeatLine };
int MovingTNSBeginParam_Axial3D_HeatLine[1] = { 0 };


// ======================================================================
// parameters for heatline   
// ======================================================================
void  ParamsFct_HeatLine(double *in, double *out);

int N_FESpaces_HeatLine = 2;
int N_FEFct_HeatLine = 5;
int N_ParamFct_HeatLine = 1;
int N_FEValues_HeatLine = 7;
int N_Parameters_HeatLine = 7;
int FEFctIndex_HeatLine[7] = { 0, 0, 0, 1, 2, 1, 2};
MultiIndex2D FEValueMultiIndex_HeatLine[7] = {D00, D10, D01, D00, D00, D01, D10};
ParamFct *ParamFct_HeatLineAll[1] = { ParamsFct_HeatLine };
int BeginParam_HeatLine[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, u1 previous, u2 previous
// used for : SUPG
// ========================================================================
void TimeNSParams4(double *in, double *out);
int TimeNSN_FESpaces4 = 1; 
int TimeNSN_Fct4 = 4;
int TimeNSN_FEValues4 = 4;
int TimeNSN_Params4 = 4; 
int TimeNSFEFctIndex4[4] = { 0, 1, 2, 3 }; // size: TimeNSN_FEValues2
MultiIndex2D TimeNSFEMultiIndex4[4] = { D00, D00, D00, D00 }; // size:TimeNSN_FEValues2

int TimeNSN_ParamFct4 = 1;
ParamFct *TimeNSFct4[1] = { TimeNSParams4 }; // size: TimeNSN_ParamFct2
int TimeNSBeginParam4[1] = { 0 }; // size: TimeNSN_ParamFct2
 
 
 
// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, gradient(tau1), gradient(tau2), gradient(tau3)
// ========================================================================
void TimeNSParamsVelo_GradVelo_CST(double *in, double *out);

int TimeNSN_FESpacesVelo_GradVelo_CST = 3;
int TimeNSN_FctVelo_GradVelo_CST = 6;
int TimeNSN_ParamFctVelo_GradVelo_CST = 1;
int TimeNSN_FEValuesVelo_GradVelo_CST = 15;
int TimeNSN_ParamsVelo_GradVelo_CST = 15;
int TimeNSFEFctIndexVelo_GradVelo_CST[15] = { 0, 1, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5 };
MultiIndex2D TimeNSFEMultiIndexVelo_GradVelo_CST[15] = { D00, D00, D10, D10, D01, D01, D00, D00, D00, D10, D10, D10, D01, D01, D01 };
ParamFct *TimeNSFctVelo_GradVelo_CST[1] = { TimeNSParamsVelo_GradVelo_CST };
int TimeNSBeginParamVelo_GradVelo_CST[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), grad(tau1, tau2, tau3) , D1, D2, D3
// ========================================================================
void TimeNSParamsVelo_GradVelo_CST_DEVSS(double *in, double *out);

int TimeNSN_FESpacesVelo_GradVelo_CST_DEVSS = 4;
int TimeNSN_FctVelo_GradVelo_CST_DEVSS = 9;
int TimeNSN_ParamFctVelo_GradVelo_CST_DEVSS = 1;
int TimeNSN_FEValuesVelo_GradVelo_CST_DEVSS = 18;
int TimeNSN_ParamsVelo_GradVelo_CST_DEVSS = 18;
int TimeNSFEFctIndexVelo_GradVelo_CST_DEVSS[18] = { 0, 1, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5, 6, 7, 8};
MultiIndex2D TimeNSFEMultiIndexVelo_GradVelo_CST_DEVSS[18] = { D00, D00, D10, D10, D01, D01, D00, D00, D00, D10, D10, D10, D01, D01, D01, D00, D00, D00};
ParamFct *TimeNSFctVelo_GradVelo_CST_DEVSS[1] = { TimeNSParamsVelo_GradVelo_CST_DEVSS };
int TimeNSBeginParamVelo_GradVelo_CST_DEVSS[1] = { 0 };
 
 
// ========================================================================
// parameters: tau1old, tau2old, tau3old, gradient(tau1), gradient(tau2), gradient(tau3)
// ========================================================================
void TimeNSParams_CST(double *in, double *out);
int TimeNSN_FESpaces_CST = 1;
int TimeNSN_Fct_CST = 3;
int TimeNSN_ParamFct_CST = 1;
int TimeNSN_FEValues_CST = 9;
int TimeNSN_ParamsVelo_CST = 9;
int TimeNSFEFctIndex_CST[9] = { 0, 1, 2, 0, 1, 2, 0, 1, 2 };
MultiIndex2D TimeNSFEMultiIndex_CST[9] = { D00, D00, D00, D10, D10, D10, D01, D01, D01 };
ParamFct *TimeNSFct_CST[1] = { TimeNSParams_CST };
int TimeNSBeginParam_CST[1] = { 0 };
 
 
// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, gradient(tau1), gradient(tau2), gradient(tau3),  gridv_x, gridv_y, gridv_x, gridv_y
// ========================================================================
void MovingNSParamsVelo_GradVelo_CST_Axial3D(double *in, double *out);

int MovingNSN_FESpacesVelo_GradVelo_CST_Axial3D = 4;
int MovingNSN_FctVelo_GradVelo_CST_Axial3D = 8;
int MovingNSN_ParamFctVelo_GradVelo_CST_Axial3D = 1;
int MovingNSN_FEValuesVelo_GradVelo_CST_Axial3D = 19;
int MovingNSN_ParamsVelo_GradVelo_CST_Axial3D = 19;
int MovingNSFEFctIndexVelo_GradVelo_CST_Axial3D[19] = { 0, 1, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5, 0, 1, 6, 7 };
MultiIndex2D MovingNSFEMultiIndexVelo_GradVelo_CST_Axial3D[19] = { D00, D00, D10, D10, D01, D01, D00, D00, D00, D10, D10, D10, D01, D01, D01, D00, D00, D00, D00};
ParamFct *MovingNSFctVelo_GradVelo_CST_Axial3D[1] = { MovingNSParamsVelo_GradVelo_CST_Axial3D };
int MovingNSBeginParamVelo_GradVelo_CST_Axial3D[1] = { 0 };
 
 
 
 
 
// ========================================================================
// boundary values for auxiliary problem in Galdi/Layton model
// ========================================================================

void BoundConditionAuxProblem(int i, double t, BoundCond &cond)
{
  cond = NEUMANN;
  //cond = DIRICHLET;
}

void BoundValueAuxProblem(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// boundary values for higher order fe in VMS
// ========================================================================

void ho_BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void ho_BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}
#endif
