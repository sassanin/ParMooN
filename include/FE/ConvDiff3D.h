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
// %W% %G%
//
// common declaration for all 3D convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF3D__
#define __CONVDIFF3D__

// part for standard Galerkin
static int N_Terms = 4;
static MultiIndex3D Derivatives[4] = { D100, D010, D001, D000 };
static int SpacesNumbers[4] = { 0, 0, 0, 0 };

// part for SDFEM (without 2nd derivatives)
static int N_Terms_SD = 4;
static MultiIndex3D Derivatives_SD[4] = { D100, D010, D001, D000};
static int SpacesNumbers_SD[4] = { 0, 0, 0, 0 };

/*
// part for UPWIND with lumping of reaction term and rhs
int N_Terms_UPW1 = 2;
MultiIndex3D Derivatives_UPW1[2] = { D10, D01 };
int SpacesNumbers_UPW1[2] = { 0, 0 };
*/

// part for UPWIND without lumping of reaction term and rhs
static int N_Terms_UPW2 = 4;
static MultiIndex3D Derivatives_UPW2[4] = { D100, D010, D001, D000 };
static int SpacesNumbers_UPW2[4] = { 0, 0, 0, 0 };

// part for all
static int N_Matrices = 1;
static int RowSpace[1] = { 0 };
static int ColumnSpace[1] = { 0 };
static int N_Rhs = 1;
static int RhsSpace[1] = { 0 };

static MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };

void BilinearAssemble(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SD(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

/*
void BilinearAssemble_UPW1(double Mult, double *coeff, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
*/

void BilinearAssemble_UPW2(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD_Orthogonal(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// parameters for DC/CD shock capturing scheme
void DC_CD_Params(double *in, double *out);

int DC_CD_N_FESpaces = 1;
int DC_CD_N_Fct = 2;
int DC_CD_N_ParamFct = 1;
int DC_CD_N_FEValues = 2;
int DC_CD_N_Params = 2;
int DC_CD_FEFctIndex[2] = { 0, 1 };
MultiIndex3D DC_CD_FEMultiIndex[2] = { D000, D000 };
ParamFct *DC_CD_Fct[1] = { DC_CD_Params };
int DC_CD_BeginParam[1] = { 0 };

// parameters for MBE shock capturing scheme
void MBE_Params(double *in, double *out);

int MBE_N_FESpaces = 1;
int MBE_N_Fct = 1;
int MBE_N_ParamFct = 1;
int MBE_N_FEValues = 4;
int MBE_N_Params = 4;
int MBE_FEFctIndex[4] = { 0, 0, 0, 0  };
MultiIndex3D MBE_FEMultiIndex[4] = { D000, D100, D010, D001 };
ParamFct *MBE_Fct[1] = { MBE_Params };
int MBE_BeginParam[1] = { 0 };

// parameters for SC_2 shock capturing scheme
void SC_2_Params(double *in, double *out);

int SC_2_N_FESpaces = 2;
int SC_2_N_Fct = 3;
int SC_2_N_ParamFct = 1;
int SC_2_N_FEValues = 5;
int SC_2_N_Params = 5;
int SC_2_FEFctIndex[5] = { 0, 1, 2, 2, 2 };
MultiIndex3D SC_2_FEMultiIndex[5] = { D000, D000, D100, D010, D001 };
ParamFct *SC_2_Fct[1] = { SC_2_Params };
int SC_2_BeginParam[1] = { 0 };

// parameters for SOLD schemes
void SOLD_Params(double *in, double *out);

int SOLD_N_FESpaces = 2;
int SOLD_N_Fct = 3;
int SOLD_N_ParamFct = 1;
int SOLD_N_FEValues = 6;
int SOLD_N_Params = 6;
int SOLD_FEFctIndex[6] = { 0, 0, 0, 0, 1, 2 };
MultiIndex3D SOLD_FEMultiIndex[6] = { D000, D100, D010, D001, D000, D000 };
ParamFct *SOLD_Fct[1] = { SOLD_Params };
int SOLD_BeginParam[1] = { 0 };


#endif // __CONVDIFF3D__
