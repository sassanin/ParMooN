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
   
#ifndef MPIMIS_brownANTS
#define MPIMIS_brownANTS

#include <math.h>

// General constants
// #define PI 3.1415926535897932384626433832795
#define SQRT_12 3.464101615137754587054892683011744733885610507620761256111613958903866
#define ONE_THIRD 0.333333333333333333333333333333333333333333333333333333333333333333333
#define TWO_THIRD 0.666666666666666666666666666666666666666666666666666666666666666666666

const int eins = 1;
//const double deins = 1.0;
//const double meins = -1.0;

// Breakage constants
/*
const double H_V = 1.37e8;
const double mu = 3.08e9;
const double Gamma_Kr = 6.7;
const double k_attr = 0.01;
const double n_attr = 1.0;
const double ro_d = 1323.0;
const double k_V = PI / 6.0;

const double C = pow(pow(H_V, 2.0/3.0)/(mu * Gamma_Kr), 0.25);
const double L_min = 32 * mu * Gamma_Kr / (3 * ipow(H_V, 2));
*/

// Aggregation constants
const double k_B = 1.3806504e-23;



#endif // MPIMIS_brownANTS
