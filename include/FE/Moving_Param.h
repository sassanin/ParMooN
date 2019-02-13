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
   
#include <FEDatabase3D.h>

// ========================================================================
// parameter routine
// ========================================================================
void MovingParamsVelo3D(double *in, double *out);

// ========================================================================
// settings for parameter routine
// ========================================================================
int MovingN_FESpacesVelo = 3;
int MovingN_FctVelo = 6;
int MovingN_ParamFctVelo = 1;
int MovingN_FEValuesVelo = 6;
int MovingN_ParamsVelo = 3;
int MovingFEFctIndexVelo[6] = { 0, 1, 2, 3, 4, 5 };
MultiIndex3D MovingFEMultiIndexVelo[6] = { D000, D000, D000,
                                           D000, D000, D000 };
ParamFct *MovingFctVelo[1] = { MovingParamsVelo3D };
int MovingBeginParamVelo[1] = { 0 };
