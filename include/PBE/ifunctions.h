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
   
#ifndef IFUNCTIONS
#define IFUNCTIONS

#include "constants.h"
#include "basictools.h"

double nu_frag(double x, double* r, function_3D3D velocity, double t);
double V_attr(double x, double* r, function_3D3D velocity, double t);
double E_kin(double x, double* r, function_3D3D velocity, double t);
const double L_min = 32 * mu * Gamma_Kr / (3 * H_V * H_V);
double L_max(double x, double* r, function_3D3D velocity, double t);
double L_star(double x, double* r, function_3D3D velocity, double t);
double b(double* r, function_3D3D velocity, double t);

#endif

