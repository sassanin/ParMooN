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
// setting for error calculation for all types
// ======================================================================
MultiIndex3D NSAllDerivatives[4] = { D000, D100, D010, D001 };
MultiIndex3D NSErrorEstiamtorU_Derivatives[8] = { D100, D010, D001, D000,  
                                                  D100, D010, D001, D000  };
MultiIndex3D NSErrorEstiamtorP_Derivatives[4] = { D100, D010, D000 };

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old, u3old
// ========================================================================
void NSParamsVelo3D(double *in, double *out);

int NSN_FESpacesVelo = 1;
int NSN_FctVelo = 3;
int NSN_ParamFctVelo = 1;
int NSN_FEValuesVelo = 3;
int NSN_ParamsVelo = 3;
int NSFEFctIndexVelo[3] = { 0, 1, 2 };
MultiIndex3D NSFEMultiIndexVelo[3] = { D000, D000, D000 };
ParamFct *NSFctVelo[1] = { NSParamsVelo3D };
int NSBeginParamVelo[1] = { 0 };

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVelo3D(double *in, double *out);

int NSN_FESpacesVelo_GradVelo = 1;
int NSN_FctVelo_GradVelo = 3;
int NSN_ParamFctVelo_GradVelo = 1;
int NSN_FEValuesVelo_GradVelo = 12;
int NSN_ParamsVelo_GradVelo = 12;
int NSFEFctIndexVelo_GradVelo[12] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
MultiIndex3D NSFEMultiIndexVelo_GradVelo[12] = { D000, D000, D000, 
                                                 D100, D100, D100,
                                                 D010, D010, D010,
                                                 D001, D001, D001};
ParamFct *NSFctVelo_GradVelo[1] = { NSParamsVelo_GradVelo3D };
int NSBeginParamVelo_GradVelo[1] = { 0 };

// ======================================================================
// auxiliary problem
// ======================================================================
void NSAuxProblem(double Mult, double *coeff, 
                  double *param, double hK, 
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs);


int NSAuxProblemN_Terms = 4;
MultiIndex3D NSAuxProblemDerivatives[4] = { D100, D010, D001, D000};
int NSAuxProblemSpaceNumbers[4] = { 0, 0, 0, 0};
int NSAuxProblemN_Matrices = 1;
int NSAuxProblemRowSpace[1] = { 0 };
int NSAuxProblemColumnSpace[1] = { 0 };
int NSAuxProblemN_Rhs = 3;
int NSAuxProblemRhsSpace[3] = { 0, 0, 0 };

void NSParamsVeloExact(double *in, double *out);
ParamFct *NSFctVeloExact[1] = { NSParamsVeloExact };

// ========================================================================
// parameters: separated pressure
// ========================================================================
void NSParamsPressSep(double *in, double *out);

int NSN_FESpacesPressSep = 1;
int NSN_FctPressSep = 1;
int NSN_ParamFctPressSep = 1;
int NSN_FEValuesPressSep = 3;
int NSN_ParamsPressSep = 3;
int NSFEFctIndexPressSep[3] = { 0, 0, 0 };
MultiIndex3D NSFEMultiIndexPressSep[3] = { D100, D010, D001 };
ParamFct *NSFctPressSep[1] = { NSParamsPressSep };
int NSBeginParamPressSep[1] = { 0 };
