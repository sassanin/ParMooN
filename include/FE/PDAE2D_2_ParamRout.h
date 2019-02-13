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
// common declaration for time dependent semiconductor device equation
// ======================================================================

#ifndef __PDAE2D_2_PARAMROUT__
#define __PDAE2D_2_PARAMROUT__

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
void PDAE2D_2_Params2(double *in, double *out);

int PDAE_N_FESpaces2 = 1;
int PDAE_N_Fct2 = 2;
int PDAE_N_ParamFct2 = 1; // richtig
int PDAE_N_FEValues2 = 6;
int PDAE_N_Params2 = 6;
int PDAE_FEFctIndex2[6] = { 0, 0, 0, 1, 1, 1 };
MultiIndex2D PDAE_FEMultiIndex2[6] = { D10, D01, D00, D10, D01, D00 };
ParamFct *PDAE_Fct2[1] = { PDAE2D_2_Params2 };
int PDAE_BeginParam2[1] = { 0 };

#endif
