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
   
#include "aggregation.h"
//#include "timer.h"
#include "constants.h"
#include "MooNMD_Io.h"

int apply_integral_operator(int nx, int ny, int nz, int na, double* input, double* grad_v, double* temp, double* output, double* grid, double L_max, double* params) {
	static paggregation pa = new aggregation(0, L_max, na, nx, ny, nz, grid, params);
//	print_vector(pa->psd.grid, pa->psd.n + 1);

	clear_vector(output, nx * ny * nz * na);

//	timer t;
//	t.start_timer();
//	scale_vector(input, nx * ny * nz * na, 1.0 / f_max);
//	printf("Input norm is %g\n", euclidean_norm(input, na * nx * ny * nz));
	pa->apply_aggregation(input, output, grad_v, temp);
//	t.stop_timer();
//	printf("Total time for aggregation was %g sec.\n", t.get_time());
//	scale_vector(output, nx * ny * nz * na, 1.0 / f_max);
//	print_vector(input, na);
//	printf("\n\n");
//	print_vector(output, na);
//	printf("Output norm is %g\n\n", euclidean_norm(output, na * nx * ny * nz));
//	pa->~Aggregation();

	//delete pa;

	return 0;
}
