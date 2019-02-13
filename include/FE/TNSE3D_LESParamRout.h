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
// @(#)TNSE3D_LESParamRout.h        1.2 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#ifndef __TNSE3D_LESPARAMROUT__
#define __TNSE3D_LESPARAMROUT__

// ========================================================================
// parameters: u1old, u2old, u3old
// ========================================================================
void TimeNSParamsVelo3D(double *in, double *out);

int TimeNSN_FESpacesVeloLES = 1;
int TimeNSN_FctVeloLES = 3;
int TimeNSN_ParamFctVeloLES = 1;
int TimeNSN_FEValuesVeloLES = 3;
int TimeNSN_ParamsVeloLES = 3;
int TimeNSFEFctIndexVeloLES[3] = { 0, 1, 2 };
MultiIndex3D TimeNSFEMultiIndexVeloLES[3] = { D000, D000, D000 };
ParamFct *TimeNSFctVeloLES[1] = { TimeNSParamsVelo3D };
int TimeNSBeginParamVeloLES[1] = { 0 };


// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsVelo_GradVelo3D(double *in, double *out);

int TimeNSN_FESpacesVelo_GradVeloLES = 1;
int TimeNSN_FctVelo_GradVeloLES = 3;
int TimeNSN_ParamFctVelo_GradVeloLES = 1;
int TimeNSN_FEValuesVelo_GradVeloLES = 12;
int TimeNSN_ParamsVelo_GradVeloLES = 12;
int TimeNSFEFctIndexVelo_GradVeloLES[12] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
MultiIndex3D TimeNSFEMultiIndexVelo_GradVeloLES[12] = { D000, D000, D000, 
                                                 D100, D100, D100,
                                                 D010, D010, D010,
                                                 D001, D001, D001};
ParamFct *TimeNSFctVelo_GradVeloLES[1] = { TimeNSParamsVelo_GradVelo3D };
int TimeNSBeginParamVelo_GradVeloLES[1] = { 0 };

// ========================================================================
// parameters: gradient(u1), gradient(u2), gradient(u3)
// ========================================================================
void TimeNSParamsGradVelo3D(double *in, double *out);

int TimeNSN_FESpacesGradVelo = 1;
int TimeNSN_FctGradVelo = 3;
int TimeNSN_ParamFctGradVelo = 1;
int TimeNSN_FEValuesGradVelo = 12;
int TimeNSN_ParamsGradVelo = 9;
int TimeNSFEFctIndexGradVelo[12] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2}; 
MultiIndex3D TimeNSFEMultiIndexGradVelo[12] = { D000, D000, D000, 
                                            D100, D100, D100,
                                            D010, D010, D010,
                                            D001, D001, D001};
ParamFct *TimeNSFctGradVelo[1] = { TimeNSParamsGradVelo3D };
int TimeNSBeginParamGradVelo[1] = { 0 };

// ========================================================================
// Classical LES
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================
void TimeNSParamsVelo_GradVeloNuT4_3D(double *in, double *out);

int TimeNSN_FESpacesVelo_GradVeloNuT4 = 2;
int TimeNSN_FctVelo_GradVeloNuT4 = 6;
int TimeNSN_ParamFctVelo_GradVeloNuT4 = 1;
int TimeNSN_FEValuesVelo_GradVeloNuT4 = 12;
int TimeNSN_ParamsVelo_GradVeloNuT4 = 15;
int TimeNSFEFctIndexVelo_GradVeloNuT4[15] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 
                                          3, 4, 5};
MultiIndex3D TimeNSFEMultiIndexVelo_GradVeloNuT4[15] = { D000, D000, D000, 
                                                     D100, D100, D100,
                                                     D010, D010, D010,
                                                     D001, D001, D001,
                                                     D000, D000, D000};
ParamFct *TimeNSFctVelo_GradVeloNuT4[1] = { TimeNSParamsVelo_GradVeloNuT4_3D };
int TimeNSBeginParamVelo_GradVeloNuT4[1] = { 0 };

// ========================================================================
// Galdi/Layton with convolution, rhs assembling, without g_\delta \ast u
// ========================================================================
void TimeNSParamsRHSLES3D(double *in, double *out);

int TimeNSN_FESpacesRHSGL00Convolution = 2;
int TimeNSN_FctRHSGL00Convolution = 9;
int TimeNSN_ParamFctRHSGL00Convolution = 1;
int TimeNSN_FEValuesRHSGL00Convolution = 6;
int TimeNSN_ParamsRHSGL00Convolution = 6;
int TimeNSFEFctIndexRHSGL00Convolution[6] = { 3, 4, 5, 6, 7, 8};
MultiIndex3D TimeNSFEMultiIndexRHSGL00Convolution[6] =  { D000, D000, D000,
                                                          D000, D000, D000};
ParamFct *TimeNSFctRHSGL00Convolution[1] = { TimeNSParamsRHSLES3D };
int TimeNSBeginParamRHSGL00Convolution[1] = { 0 };

// ========================================================================
// Galdi/Layton with convolution, rhs assembling, 
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================
void TimeNSParamsRHSGL00ConvolutionNuT4_3D(double *in, double *out);
int TimeNSN_FESpacesRHSGL00ConvolutionNuT4 = 3;
int TimeNSN_FctRHSGL00ConvolutionNuT4 = 12;
int TimeNSN_ParamFctRHSGL00ConvolutionNuT4 = 1;
int TimeNSN_FEValuesRHSGL00ConvolutionNuT4 = 24;
int TimeNSN_ParamsRHSGL00ConvolutionNuT4 = 21;
int TimeNSFEFctIndexRHSGL00ConvolutionNuT4[21] = { 0, 1, 2, 0, 1, 2,
                                               0, 1, 2, 0, 1, 2,
                                               3, 4, 5, 6, 7, 8,
                                               9, 10, 11};
MultiIndex3D TimeNSFEMultiIndexRHSGL00ConvolutionNuT4[21] =  { D000, D000, D000,
                                                          D100, D100, D100,
                                                          D010, D010, D010,
                                                          D001, D001, D001,
                                                          D000, D000, D000,
                                                          D000, D000, D000,
                                                          D000, D000, D000};
ParamFct *TimeNSFctRHSGL00ConvolutionNuT4[1] = { TimeNSParamsRHSGL00ConvolutionNuT4_3D };
int TimeNSBeginParamRHSGL00ConvolutionNuT4[1] = { 0 };


// ========================================================================
// parameters for Galdi/Layton model with auxiliary problem
// ========================================================================

int TimeNSN_FESpacesGL00AuxProblem = 1;
int TimeNSN_FctGL00AuxProblem = 9;
int TimeNSN_ParamFctGL00AuxProblem = 1;
int TimeNSN_FEValuesGL00AuxProblem = 6;
int TimeNSN_ParamsGL00AuxProblem = 6;
int TimeNSFEFctIndexGL00AuxProblem[6] = { 3, 4, 5, 6, 7, 8};
MultiIndex3D TimeNSFEMultiIndexGL00AuxProblem[6] = { D000, D000, D000,
                                                     D000, D000, D000};
ParamFct *TimeNSFctGL00AuxProblem[1] = { TimeNSParamsRHSLES3D };
int TimeNSBeginParamGL00AuxProblem[1] = { 0 };


// ========================================================================
// parameters for Galdi/Layton model with auxiliary problem
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================

void TimeNSParamsGL00AuxProblemNuT4_3D(double *in, double *out);

int TimeNSN_FESpacesGL00AuxProblemNuT4 = 2;
int TimeNSN_FctGL00AuxProblemNuT4 = 12;
int TimeNSN_ParamFctGL00AuxProblemNuT4 = 1;
int TimeNSN_FEValuesGL00AuxProblemNuT4 = 24;
int TimeNSN_ParamsGL00AuxProblemNuT4 = 21;
int TimeNSFEFctIndexGL00AuxProblemNuT4[21] = { 0, 1, 2, 0, 1, 2, 
                                               0, 1, 2, 0, 1, 2,
                                               3, 4, 5, 6, 7, 8,
                                               9, 10, 11};
MultiIndex3D TimeNSFEMultiIndexGL00AuxProblemNuT4[21] = { D000, D000, D000,
                                                          D100, D100, D100,
                                                          D010, D010, D010,
                                                          D001, D001, D001, 
                                                          D000, D000, D000,
                                                          D000, D000, D000,
                                                          D000, D000, D000};
ParamFct *TimeNSFctGL00AuxProblemNuT4[1] = { TimeNSParamsRHSGL00ConvolutionNuT4_3D };
int TimeNSBeginParamGL00AuxProblemNuT4[1] = { 0 };


int TimeNSN_FESpacesVelo;
int TimeNSN_FctVelo;
int TimeNSN_ParamFctVelo;
int TimeNSN_FEValuesVelo;
ParamFct *TimeNSFctVelo[1];
int TimeNSFEFctIndexVelo[3];
MultiIndex3D TimeNSFEMultiIndexVelo[3];
int TimeNSN_ParamsVelo;
int TimeNSBeginParamVelo[1];


#endif
