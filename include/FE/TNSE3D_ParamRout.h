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
// @(#)TNSE3D_ParamRout.h        1.2 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#ifndef __TNSE3D_PARAMROUT__
#define __TNSE3D_PARAMROUT__

// ======================================================================
// setting for error calculation for all types
// ======================================================================
MultiIndex3D TimeNSAllDerivatives[4] = { D000, D100, D010, D001 };

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old, u3old
// ========================================================================
void TimeNSParamsVelo3D(double *in, double *out);

int TimeNSN_FESpacesVelo = 1;
int TimeNSN_FctVelo = 3;
int TimeNSN_ParamFctVelo = 1;
int TimeNSN_FEValuesVelo = 3;
int TimeNSN_ParamsVelo = 3;
int TimeNSFEFctIndexVelo[3] = { 0, 1, 2 };
MultiIndex3D TimeNSFEMultiIndexVelo[3] = { D000, D000, D000 };
ParamFct *TimeNSFctVelo[1] = { TimeNSParamsVelo3D };
int TimeNSBeginParamVelo[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsVelo_GradVelo3D(double *in, double *out);

int TimeNSN_FESpacesVelo_GradVelo = 1;
int TimeNSN_FctVelo_GradVelo = 3;
int TimeNSN_ParamFctVelo_GradVelo = 1;
int TimeNSN_FEValuesVelo_GradVelo = 12;
int TimeNSN_ParamsVelo_GradVelo = 15;
int TimeNSFEFctIndexVelo_GradVelo[12] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
MultiIndex3D TimeNSFEMultiIndexVelo_GradVelo[12] = { D000, D000, D000, 
                                                 D100, D100, D100,
                                                 D010, D010, D010,
                                                 D001, D001, D001};
ParamFct *TimeNSFctVelo_GradVelo[1] = { TimeNSParamsVelo_GradVelo3D };
int TimeNSBeginParamVelo_GradVelo[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, u3old
// all partial derivatives
// convolution of u1old, u2old, u3old
// ========================================================================
void TimeNSParamsVelo_GradVelo_ConvVelo3D(double *in, double *out);

int TimeNSN_FESpacesVelo_GradVelo_ConvVelo = 2;
int TimeNSN_FctVelo_GradVelo_ConvVelo = 6;
int TimeNSN_ParamFctVelo_GradVelo_ConvVelo = 1;
int TimeNSN_FEValuesVelo_GradVelo_ConvVelo = 15;
int TimeNSN_ParamsVelo_GradVelo_ConvVelo = 16;
int TimeNSFEFctIndexVelo_GradVelo_ConvVelo[15] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
                                                   3, 4, 5};
MultiIndex3D TimeNSFEMultiIndexVelo_GradVelo_ConvVelo[15] = { D000, D000, D000, 
                                                              D100, D100, D100,
                                                              D010, D010, D010,
                                                              D001, D001, D001,
                                                              D000, D000, D000};
ParamFct *TimeNSFctVelo_GradVelo_ConvVelo[1] = { TimeNSParamsVelo_GradVelo_ConvVelo3D };
int TimeNSBeginParamVelo_GradVelo_ConvVelo[1] = { 0 };



/*// ========================================================================
// Galdi/Layton with convolution, without g_\delta \ast u 
// ========================================================================
void TimeNSParamsGL00Convolution(double *in, double *out);
int TimeNSN_FESpacesGL00Convolution = 2;
int TimeNSN_FctGL00Convolution = 5;
int TimeNSN_ParamFctGL00Convolution = 1;
int TimeNSN_FEValuesGL00Convolution = 9;
int TimeNSN_ParamsGL00Convolution = 9;
int TimeNSFEFctIndexGL00Convolution[9] = { 0, 1, 0, 1, 0, 1, 2, 3, 4 };
MultiIndex3D TimeNSFEMultiIndexGL00Convolution[9] = { D000, D000, D100, D100, D010, D010,
                                             D000, D000, D000};
ParamFct *TimeNSFctGL00Convolution[1] = { TimeNSParamsGL00Convolution3D };
int TimeNSBeginParamGL00Convolution[1] = { 0 };
*/
/*void TimeNSParamsGL00AuxProblemPaper2(double *in, double *out);
int TimeNSN_FESpacesGL00AuxProblemPaper2 = 1;
int TimeNSN_FctGL00AuxProblemPaper2 = 8;
int TimeNSN_ParamFctGL00AuxProblemPaper2 = 1;
int TimeNSN_FEValuesGL00AuxProblemPaper2 = 13;
int TimeNSN_ParamsGL00AuxProblemPaper2 = 11;
int TimeNSFEFctIndexGL00AuxProblemPaper2[13] = { 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6, 6, 7};
MultiIndex3D TimeNSFEMultiIndexGL00AuxProblemPaper2[13] = { D00, D00, D10, D10, D01, D01,
                                             D00, D00, D00, D10, D10, D01, D01};
ParamFct *TimeNSFctGL00AuxProblemPaper2[1] = { TimeNSParamsGL00AuxProblemPaper2 };
int TimeNSBeginParamGL00AuxProblemPaper2[1] = { 0 };
*/

// ========================================================================
// parameters: u, grad u, G^H
// ========================================================================
void TimeNSParamsVelo_GradVelo_LargeScale3D(double *in, double *out);

int TimeNSN_FESpacesVelo_GradVelo_LargeScale = 4;
int TimeNSN_FctVelo_GradVelo_LargeScale = 10;
int TimeNSN_ParamFctVelo_GradVelo_LargeScale = 1;
int TimeNSN_FEValuesVelo_GradVelo_LargeScale = 19;
int TimeNSN_ParamsVelo_GradVelo_LargeScale = 22;
int TimeNSFEFctIndexVelo_GradVelo_LargeScale[19] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
						     3, 4, 5, 6, 7, 8, 9};
MultiIndex3D TimeNSFEMultiIndexVelo_GradVelo_LargeScale[19] = { D000, D000, D000, 
								D100, D100, D100,
								D010, D010, D010,
								D001, D001, D001, 
								D000, D000, D000,
								D000, D000, D000,		
								D000};
ParamFct *TimeNSFctVelo_GradVelo_LargeScale[1] = { TimeNSParamsVelo_GradVelo_LargeScale3D };
int TimeNSBeginParamVelo_GradVelo_LargeScale[1] = { 0 };

//================used for VMS3D=============================
void TimeNSParamsVelo_GradVelo_VMS3D(double *in, double *out);
int TimeNSN_FESpacesVelo_GradVelo_VMS3D = 1; //number of fe spaces passed
int TimeNSN_FctVelo_GradVelo_VMS3D = 3;      // number of fe functions
int TimeNSN_ParamFctVelo_GradVelo_VMS3D  = 1; // always starts with 1
int TimeNSN_FEValuesVelo_GradVelo_VMS3D  = 12; //  u_1, u_2, u_3,u_1_x,u_2_x,u_3_x,u_1_y,u_2_y,u_3_y,u_1_z,u_2_z,u_3_z
int TimeNSN_ParamsVelo_GradVelo_VMS3D  = 12;  // same as above  
int TimeNSFEFctIndexVelo_GradVelo_VMS3D[12] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 }; // componets of velo 
MultiIndex3D TimeNSFEMultiIndexVelo_GradVelo_VMS3D[12] = { D000, D000, D000, D100, D100, D100, D010, D010, D010, D001, D001, D001 }; 
ParamFct* TimeNSFctVelo_GradVelo_VMS3D[1] = { TimeNSParamsVelo_GradVelo_VMS3D };
int TimeNSBeginParamVelo_GradVelo_VMS3D[1] = { 0 };
//=========================================================

//====================used for ALEVMS3D============================

void TimeNSParamsVelo_GradVelo_VMS3D_ALE(double *in, double *out);
int TimeNSN_FESpacesVelo_GradVelo_VMS3D_ALE = 2; //velocity and grid fespace
int TimeNSN_FctVelo_GradVelo_VMS3D_ALE = 6;      // u1,u2,u3,mesh1,mesh2,mesh3
int TimeNSN_ParamFctVelo_GradVelo_VMS3D_ALE  = 1; // always starts with 1
int TimeNSN_FEValuesVelo_GradVelo_VMS3D_ALE  = 15; //  u_1, u_2, u_3,
                                                   //u_1_x,u_2_x,u_3_x,
                                                   //u_1_y,u_2_y,u_3_y,
                                                   //u_1_z,u_2_z,u_3_z
                                                   //  w_1,w_2,w_3

int TimeNSN_ParamsVelo_GradVelo_VMS3D_ALE  = 15;  // same as above  
int TimeNSFEFctIndexVelo_GradVelo_VMS3D_ALE[15] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 , 3, 4, 5}; // componets of velo 
MultiIndex3D TimeNSFEMultiIndexVelo_GradVelo_VMS3D_ALE[15] = { D000, D000, D000,
                                                        D100, D100, D100,
							D010, D010, D010,
							D001, D001, D001,
                                                        D000, D000, D000}; 
							//  u_1, u_2, u_3,
							//  u_1_x,u_2_x,u_3_x,
							//  u_1_y,u_2_y,u_3_y,
							//  u_1_z,u_2_z,u_3_z
							//  w_1,w_2,w_3
							
ParamFct* TimeNSFctVelo_GradVelo_VMS3D_ALE[1] = { TimeNSParamsVelo_GradVelo_VMS3D_ALE };
int TimeNSBeginParamVelo_GradVelo_VMS3D_ALE[1] = { 0 };
//=======================================================================

// ========================================================================
// parameters for VMS
// ========================================================================
void TimeNSParams_VMS_SmallRhs3D(double *in, double *out);

int TimeNSN_FESpaces_VMS_SmallRhs = 4;
int TimeNSN_Fct_VMS_SmallRhs  = 7;
int TimeNSN_ParamFct_VMS_SmallRhs  = 1;
int TimeNSN_FEValues_VMS_SmallRhs = 25;
int TimeNSN_Params_VMS_SmallRhs = 25;
int TimeNSFEFctIndex_VMS_SmallRhs[25] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
                                          3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5,
                                          6};
MultiIndex3D TimeNSFEMultiIndex_VMS_SmallRhs[25] = { D000, D000, D000,
                                                     D100, D100, D100,
                                                     D010, D010, D010,
                                                     D001, D001, D001,
                                                     D000, D000, D000,
                                                     D100, D100, D100,
                                                     D010, D010, D010,
                                                     D001, D001, D001,
                                                     D000 };
ParamFct *TimeNSFct_VMS_SmallRhs[1] = {  TimeNSParams_VMS_SmallRhs3D };
int TimeNSBeginParam_VMS_SmallRhs[1] = { 0 };

// ========================================================================
// boundary values for higher order fe in VMS
// ========================================================================

void ho_BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

void ho_BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

#endif
