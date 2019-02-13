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
   
#include "breakage.h"
#include "aggregation.h"
#include "BrAgg.h"
#include <math.h>
#include <stdlib.h>

int main()
{
	int n = 9;
	int nx = 1, ny = 1, nz = 1, na = ((n - 1) / 2.0) * (cLevel + 1) + 1;

	srand(time(NULL));

	double* v;
	v = allocate_vector(nx * ny * nz * 3);
	fill_vector(v, nx * ny * nz * 3, 1);

	double* grad_v;
	grad_v = allocate_vector(nx * ny * nz * 9);
	fill_vector(grad_v, nx * ny * nz * 9, 1);

	double* input; input = allocate_vector(nx * ny * nz * na);
	fill_vector(input, nx * ny * nz * na, 1);

	double* output; output = allocate_vector(nx * ny * nz * na);
	clear_vector(output, nx * ny * nz * na);

	double* temp; temp = allocate_vector(nx * ny * nz);
	fill_vector(temp, nx * ny * nz, 1);

	double L_max = 1e-3, f_max = 1e3;
	Breakage_agglomeration(nx, ny, nz, n, input, v, grad_v, temp, output, L_max, f_max);

	return 0;
}
