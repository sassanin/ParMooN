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
#include <math.h>
#include <stdlib.h>

int main()
{
	int nxInit = 32, L = 5, nr1 = 17, nr2 = 17, nr3 = 17, nx = (L + 1) * nxInit / 2.0 + 1;
	double Lmin = L_min, Lmax = Lmin * pow(nxInit * pow(2, L - 1) + 1, 1.0/3.0);
	double l1 = 0.01, l2 = 0.01, l3 = 0.5;

	// TO BE CHANGED WITH REAL VELOCITY / SHEAR FLOW
	double* v = allocate_vector(nr1 * nr2 * nr3 * 3);
	fill_vector(v, nr1 * nr2 * nr3 * 3, fabs((double)rand() / (double)20000));
	double* grad_v = allocate_vector(nr1 * nr2 * nr3 * 9);
	fill_vector(grad_v, nr1 * nr2 * nr3 * 9, fabs((double)rand() / (double)10000));

	// TO BE CHANGED WITH REAL INPUT / OUTPUT
	double* input = allocate_vector(nr1 * nr2 * nr3 * nx);
	fill_vector(input, nr1 * nr2 * nr3 * nx, (double)rand() / (double)1000);
	double* output = allocate_vector(nr1 * nr2 * nr3 * nx);
	clear_vector(output, nr1 * nr2 * nr3 * nx);

	pbreakage pb = new breakage(Lmin, Lmax, nxInit, L, nr1, nr2, nr3, l1, l2, l3);
	pb->apply_breakage(input, output, (void*)v);

	paggregation pa = new aggregation(Lmin, Lmax, nxInit, L, nr1, nr2, nr3, l1, l2, l3);
	pa->apply_aggregation(input, output, (void*)grad_v);

	return 0;
}
