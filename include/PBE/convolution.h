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
   
#ifndef MPIMIS_CONVOLUTION
#define MPIMIS_CONVOLUTION

#include "basictools.h"

// Computes discrete version of w(x) = \int u(x-y) * v(y) dy on h-refined grid (a, b]
// Here is assumed that a and b are pieacewise linear, and c is also projected onto space 
// of piecewise linear functions
void linear_convolution(double* u, double* v, double* w, int start_n, int L, double a, double b);

// Gets Legendre coefficients for linear functions
void get_Legendre_coefficients(double* fi, double* fiplus1, double* xi, double* xiplus1, double* coeff0, double* coeff1);

#endif // MPIMIS_CONVOLUTION
