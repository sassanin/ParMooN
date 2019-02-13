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
   
#include "ifunctions.h"
#include "constants.h"
#include "attrition.h"
#include "basictools.h"

double nu_frag(double x, double* r, function_3D3D velocity, double t)
{
	return V_attr(x, r, velocity, t) * (pow(L_min, -2.25) - pow(L_max(x, r, velocity, t), -2.25)) /
		   (3 * k_V * (pow(L_max(x, r, velocity, t), 0.75) - pow(L_min, 0.75)));
}

double V_attr(double x, double* r, function_3D3D velocity, double t)
{
	return 2 * pow(H_V, 2/3) * E_kin(x, r, velocity, t) / (3 * mu * Gamma_Kr);
}

double E_kin(double x, double* r, function_3D3D velocity, double t)
{
	return 0.5 * ro_d * k_V * x * x * x * pow(spectral_norm(velocity(r, t), 2), 2);
}

double L_max(double x, double* r, function_3D3D velocity, double t)
{
	return 0.5 * pow(H_V, 2/9) * E_kin(x, r, velocity, t) / (pow(mu,1/3) * pow(Gamma_Kr, 1/3));
}

double L_star(double x, double* r, function_3D3D velocity, double t)
{
	return pow(x * x * x - V_attr(x, r, velocity, t) / k_V, 1/3);
}

double b(double* r, function_3D3D velocity, double t)
{
	if(is_on_boundary(r))
		return 0;
	return 1; // TODO
}
