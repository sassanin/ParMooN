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
   
#ifndef __NSE2D_PARAMROUT__
#define __NSE2D_PARAMROUT__
#include <NSE2D_FixPo.h>

// ======================================================================
// setting for error calculation for all types
// ======================================================================
MultiIndex2D NSAllDerivatives[3] = { D00, D10, D01 };
MultiIndex2D NSErrorEstiamtorU_Derivatives[8] = { D10, D01, D00, D20, D02, 
                                         D10, D01, D00  };
MultiIndex2D NSErrorEstiamtorU_DerivativesEstimator[5] = { D10, D01, D00, D20, D02};
MultiIndex2D NSErrorEstiamtorP_Derivatives[3] = { D10, D01, D00 };
MultiIndex2D NSErrorEstiamtorP_DerivativesEstimator[2] = { D10, D01 };

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
int NSN_FESpacesVelo = 1;
int NSN_FctVelo = 2;
int NSN_ParamFctVelo = 1;
int NSN_FEValuesVelo = 2;
int NSN_ParamsVelo = 2;
int NSFEFctIndexVelo[2] = { 0, 1 };
MultiIndex2D NSFEMultiIndexVelo[2] = { D00, D00 };
ParamFct *NSFctVelo[1] = { NSParamsVelo };
int NSBeginParamVelo[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVelo(double *in, double *out);

int NSN_FESpacesVelo_GradVelo = 1;
int NSN_FctVelo_GradVelo = 2;
int NSN_ParamFctVelo_GradVelo = 1;
int NSN_FEValuesVelo_GradVelo = 6;
int NSN_ParamsVelo_GradVelo = 6;
int NSFEFctIndexVelo_GradVelo[6] = { 0, 1, 0, 1, 0, 1 };
MultiIndex2D NSFEMultiIndexVelo_GradVelo[6] = { D00, D00, D10, D10, D01, D01 };
ParamFct *NSFctVelo_GradVelo[1] = { NSParamsVelo_GradVelo };
int NSBeginParamVelo_GradVelo[1] = { 0 };

// ======================================================================
// auxiliary problem
// ======================================================================
void NSAuxProblem(double Mult, double *coeff, 
                   double *param, double hK, 
                   double **OrigValues, int *N_BaseFuncts,
                   double ***LocMatrices, double **LocRhs);


int NSAuxProblemN_Terms = 3;
MultiIndex2D NSAuxProblemDerivatives[3] = { D10, D01, D00};
int NSAuxProblemSpaceNumbers[3] = { 0, 0, 0};
int NSAuxProblemN_Matrices = 1;
int NSAuxProblemRowSpace[1] = { 0 };
int NSAuxProblemColumnSpace[1] = { 0 };
int NSAuxProblemN_Rhs = 2;
int NSAuxProblemRhsSpace[2] = { 0, 0 };

void NSParamsVeloExact(double *in, double *out);
ParamFct *NSFctVeloExact[1] = { NSParamsVeloExact };

// ======================================================================
// auxiliary problem for differential filter
// ======================================================================
void Filter_Galerkin(double Mult, double *coeff, 
                     double *param, double hK, 
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);
// ======================================================================
//  declarations for auxiliary problem for differential filter
//      one matrix 
//      two rhs
// ======================================================================

int Filter_N_Terms = 3;
MultiIndex2D Filter_Derivatives[3] = { D10, D01, D00};
int Filter_SpaceNumbers[3] = { 0, 0, 0};
int Filter_N_Matrices = 1;
int Filter_RowSpace[1] = { 0 };
int Filter_ColumnSpace[1] = { 0 };
int Filter_N_Rhs = 2;
int Filter_RhsSpace[2] = { 0, 0 };

// ========================================================================
// parameters: separated pressure
// ========================================================================
void NSParamsPressSep(double *in, double *out);

int NSN_FESpacesPressSep = 1;
int NSN_FctPressSep = 1;
int NSN_ParamFctPressSep = 1;
int NSN_FEValuesPressSep = 2;
int NSN_ParamsPressSep = 2;
int NSFEFctIndexPressSep[2] = { 0, 0 };
MultiIndex2D NSFEMultiIndexPressSep[2] = { D10, D01 };
ParamFct *NSFctPressSep[1] = { NSParamsPressSep };
int NSBeginParamPressSep[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void NSParamsVeloAxialSymm3D(double *in, double *out);

int NSN_FESpacesVeloAxialSymm3D = 1;
int NSN_FctVeloAxialSymm3D = 2;
int NSN_ParamFctVeloAxialSymm3D = 1;
int NSN_FEValuesVeloAxialSymm3D = 2;
int NSN_ParamsVeloAxialSymm3D = 3;
int NSFEFctIndexVeloAxialSymm3D[2] = { 0, 1 };
MultiIndex2D NSFEMultiIndexVeloAxialSymm3D[2] = { D00, D00 };
ParamFct *NSFctVeloAxialSymm3D[1] = { NSParamsVeloAxialSymm3D };
int NSBeginParamVeloAxialSymm3D[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVeloAxialSymm3D(double *in, double *out);

int NSN_FESpacesVelo_GradVeloAxialSymm3D = 1;
int NSN_FctVelo_GradVeloAxialSymm3D = 2;
int NSN_ParamFctVelo_GradVeloAxialSymm3D = 1;
int NSN_FEValuesVelo_GradVeloAxialSymm3D = 6;
int NSN_ParamsVelo_GradVeloAxialSymm3D = 7;
int NSFEFctIndexVelo_GradVeloAxialSymm3D[6] = { 0, 1, 0, 1, 0, 1 };
MultiIndex2D NSFEMultiIndexVelo_GradVeloAxialSymm3D[6] = { D00, D00, D10, D10, D01, D01 };
ParamFct *NSFctVelo_GradVeloAxialSymm3D[1] = { NSParamsVelo_GradVeloAxialSymm3D };
int NSBeginParamVelo_GradVeloAxialSymm3D[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_Press(double *in, double *out);

int NSN_FESpacesVelo_Press = 2;
int NSN_FctVelo_Press = 3;
int NSN_ParamFctVelo_Press = 1;
int NSN_FEValuesVelo_Press = 4;
int NSN_ParamsVelo_Press = 4;
int NSFEFctIndexVelo_Press[4] = { 0, 1, 2, 2 };
MultiIndex2D NSFEMultiIndexVelo_Press[6] = { D00, D00, D10, D10};
ParamFct *NSFctVelo_Press[1] = { NSParamsVelo_Press};
int NSBeginParamVelo_Press[1] = { 0 };


// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, gradient(tau1), gradient(tau2), gradient(tau3)
// ========================================================================
void NSParamsVelo_GradVelo_CST(double *in, double *out);

int NSN_FESpacesVelo_GradVelo_CST = 3;
int NSN_FctVelo_GradVelo_CST = 6;
int NSN_ParamFctVelo_GradVelo_CST = 1;
int NSN_FEValuesVelo_GradVelo_CST = 15;
int NSN_ParamsVelo_GradVelo_CST = 15;
int NSFEFctIndexVelo_GradVelo_CST[15] = { 0, 1, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5 };
MultiIndex2D NSFEMultiIndexVelo_GradVelo_CST[15] = { D00, D00, D10, D10, D01, D01, D00, D00, D00, D10, D10, D10, D01, D01, D01 };
ParamFct *NSFctVelo_GradVelo_CST[1] = { NSParamsVelo_GradVelo_CST };
int NSBeginParamVelo_GradVelo_CST[1] = { 0 };



// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, gradient(tau1), gradient(tau2), gradient(tau3) (Axial3D)
// ========================================================================
void NSParamsVelo_GradVelo_CST_Axial3D(double *in, double *out);

int NSN_FESpacesVelo_GradVelo_CST_Axial3D = 3;
int NSN_FctVelo_GradVelo_CST_Axial3D = 6;
int NSN_ParamFctVelo_GradVelo_CST_Axial3D = 1;
int NSN_FEValuesVelo_GradVelo_CST_Axial3D = 16;
int NSN_ParamsVelo_GradVelo_CST_Axial3D = 16;
int NSFEFctIndexVelo_GradVelo_CST_Axial3D[16] = { 0, 1, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5, 2};
MultiIndex2D NSFEMultiIndexVelo_GradVelo_CST_Axial3D[16] = { D00, D00, D10, D10, D01, D01, D00, D00, D00, D10, D10, D10, D01, D01, D01, D00};
ParamFct *NSFctVelo_GradVelo_CST_Axial3D[1] = { NSParamsVelo_GradVelo_CST_Axial3D };
int NSBeginParamVelo_GradVelo_CST_Axial3D[1] = { 0 };





// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, D1, D2, D3, p
// ========================================================================
void NSParamsVelo_GradVelo_DEVSS(double *in, double *out);

int NSN_FESpacesVelo_GradVelo_DEVSS = 4;
int NSN_FctVelo_GradVelo_DEVSS = 9;
int NSN_ParamFctVelo_GradVelo_DEVSS = 1;
int NSN_FEValuesVelo_GradVelo_DEVSS = 12;
int NSN_ParamsVelo_GradVelo_DEVSS = 12;
int NSFEFctIndexVelo_GradVelo_DEVSS[12] = { 0, 1, 0, 1, 0, 1, 3, 4, 5, 6, 7, 8};
MultiIndex2D NSFEMultiIndexVelo_GradVelo_DEVSS[12] = { D00, D00, D10, D10, D01, D01, D00, D00, D00, D00, D00, D00};
ParamFct *NSFctVelo_GradVelo_DEVSS[1] = { NSParamsVelo_GradVelo_DEVSS };
int NSBeginParamVelo_GradVelo_DEVSS[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), grad(tau1, tau2, tau3) , D1, D2, D3
// ========================================================================
void NSParamsVelo_GradVelo_CST_DEVSS(double *in, double *out);

int NSN_FESpacesVelo_GradVelo_CST_DEVSS = 4;
int NSN_FctVelo_GradVelo_CST_DEVSS = 9;
int NSN_ParamFctVelo_GradVelo_CST_DEVSS = 1;
int NSN_FEValuesVelo_GradVelo_CST_DEVSS = 18;
int NSN_ParamsVelo_GradVelo_CST_DEVSS = 18;
int NSFEFctIndexVelo_GradVelo_CST_DEVSS[18] = { 0, 1, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 3, 4, 5, 6, 7, 8};
MultiIndex2D NSFEMultiIndexVelo_GradVelo_CST_DEVSS[18] = { D00, D00, D10, D10, D01, D01, D00, D00, D00, D10, D10, D10, D01, D01, D01, D00, D00, D00};
ParamFct *NSFctVelo_GradVelo_CST_DEVSS[1] = { NSParamsVelo_GradVelo_CST_DEVSS };
int NSBeginParamVelo_GradVelo_CST_DEVSS[1] = { 0 };

// ========================================================================
// parameters: tau1old, tau2old, tau3old, gradient(tau1), gradient(tau2), gradient(tau3)
// ========================================================================
void NSParams_CST(double *in, double *out);
int NSN_FESpaces_CST = 1;
int NSN_Fct_CST = 3;
int NSN_ParamFct_CST = 1;
int NSN_FEValues_CST = 9;
int NSN_ParamsVelo_CST = 9;
int NSFEFctIndex_CST[9] = { 0, 1, 2, 0, 1, 2, 0, 1, 2 };
MultiIndex2D NSFEMultiIndex_CST[9] = { D00, D00, D00, D10, D10, D10, D01, D01, D01 };
ParamFct *NSFct_CST[1] = { NSParams_CST };
int NSBeginParam_CST[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, temperature
// ========================================================================
/*
void NSParamsVeloTemp(double *in, double *out);

int NSN_FESpacesVeloTemp = 3;
int NSN_FctVeloTemp = 3;
int NSN_ParamFctVeloTemp = 1;
int NSN_FEValuesVeloTemp = 3;
int NSN_ParamsVeloTemp = 3;
int NSFEFctIndexVeloTemp[3] = { 0, 1, 2 };
MultiIndex2D NSFEMultiIndexVeloTemp[3] = { D00, D00, D00 };
ParamFct *NSFctVeloTemp[1] = { NSParamsVeloTemp };
int NSBeginParamVeloTemp[1] = { 0 };
*/

#endif // __NSE2D_PARAMROUT__
