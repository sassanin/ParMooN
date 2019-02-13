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
   
#ifndef MPIMIS_TOOLS
#define MPIMIS_TOOLS

#include "basictools.h"

#include <math.h>
#include <stdlib.h>

inline bool isZero(double p_dValue);
inline double modulus(double p_dA, double p_dB);

inline bool isZero(double p_dValue) { // Tells if a number is zero considering a threshold 
	return (abs(p_dValue) < ZERO_THRESHOLD ? true : false);
}
inline double modulus(double p_dA, double p_dB) { // Gives the modulus of a complex number: p_dA + p_dB * i
	return sqrt(sqr(p_dA) + sqr(p_dB));
}

#endif // MPIMIS_TOOLS
