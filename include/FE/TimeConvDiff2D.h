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
// @(#)TimeConvDiff2D.h        1.6 04/13/00
//
// common declaration for time dependent convection diffusion problems
// ======================================================================
#include <Enumerations.h>
#include <FESpace2D.h>

#ifndef __TIMECONVDIFF2D__
#define __TIMECONVDIFF2D__

// part for standard Galerkin
// int N_Terms = 3;
// MultiIndex2D Derivatives[3] = { D10, D01, D00 };
// int SpacesNumbers[3] = { 0, 0, 0 };

// part for SDFEM
// int N_Terms_SD = 5;
// MultiIndex2D Derivatives_SD[5] = { D10, D01, D00, D20, D02 };
// int SpacesNumbers_SD[5] = { 0, 0, 0, 0 };

// part for all
// int N_Matrices = 1;
// int RowSpace[2] = { 0 };
// int ColumnSpace[2] = { 0 };
// int N_Rhs = 1;
// int RhsSpace[1] = { 0 };

// part for mass matrix
int N_Matrices_Mass = 1;
int RowSpace_Mass[1] = { 0 };
int ColumnSpace_Mass[1] = { 0 };
int N_Rhs_Mass= 0;
int N_Terms_Mass = 1;
MultiIndex2D Derivatives_Mass[1] = { D00 };
int SpacesNumbers_Mass[1] = { 0 };

// MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };

void TimeBilinearAssemble(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeMassAssemble(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeBilinearAssemble_SD(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeBilinearAssembleRB1(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeBilinearAssembleRB(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeRhsAssembleRB(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

// ======================================================================
// parameter routine settings
// convection is finite element velocity field 
// ======================================================================
void TimeCDParamsVeloField(double *in, double *out);

int TimeCDParamsVeloFieldN_FESpaces = 1;
int TimeCDParamsVeloFieldN_Fct = 2;
int TimeCDParamsVeloFieldN_ParamFct = 1;
int TimeCDParamsVeloFieldN_FEValues = 2;
int TimeCDParamsVeloFieldN_Params = 3;
int TimeCDParamsVeloFieldFEFctIndex[2] = { 0, 1 };
MultiIndex2D TimeCDParamsVeloFieldFEMultiIndex[2] = { D00, D00 };
ParamFct *TimeCDParamsVeloFieldFct[1] = { TimeCDParamsVeloField };
int TimeCDParamsVeloFieldBeginParam[1] = { 0 };


void TimeCDParamsVeloField_ALE(double *in, double *out);

int TimeCDParamsVeloFieldN_Params_ALE = 4;
int TimeCDParamsVeloFieldN_FEValues_ALE = 4;
int TimeCDParamsVeloFieldFEFctIndex_ALE[4] = { 0, 1, 2, 3 };
MultiIndex2D TimeCDParamsVeloFieldFEMultiIndex_ALE[4] = { D00, D00,  D10, D01 };
ParamFct *TimeCDParamsVeloFieldFct_ALE[1] = { TimeCDParamsVeloField_ALE };
#endif


void MatrixMARhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs);

void MatrixMARhsAssemble(double Mult, double *coeff, double *param,
                         double hK, double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs);

void MatrixMARhsALEAssemble_SUPG(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs);

void MatrixMRhsALEAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);

void MatrixARhsAssembleHeatLine(double Mult, double *coeff, double *param, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void MatrixARhsAssembleHeatLine_Axial3D(double Mult, double *coeff, double *param, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);
