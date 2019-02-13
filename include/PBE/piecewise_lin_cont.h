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
   
#include "basictools.h"
#include "fepc_easy.h"
#include "faltung.h"

class piecewise_constant
{
public:
	int n;
	double* grid; // array of size n + 1 // grid[0] = a, grid[n] = b
	double* value; // array of size n

	double integral();
};

class piecewise_linear
{
public:
	int n;
	double* grid; // array of size n + 1 // grid[0] = a, grid[n] = b
	double* start_value; // array of size n
	double* end_value; // array of size n

	double integral();
};

class piecewise_linear_continuous
{
public:
	piecewise_linear_continuous() {};

	int n;
	double* grid; // array of size n + 1 // grid[0] = a, grid[n] = b
	double* coeff; // array of size n + 1 // coefficiets of the basis functions (hat functions). The same as the values at gridpoints.

	double L1();
	double integral();
	double mass();

	piecewise_linear_continuous& operator = (const piecewise_linear_continuous& f);

	int multiplication_with_function (double* f, piecewise_constant& w);
	int multiplication_with_function (double* f, piecewise_linear& w);

	int ProjectToMassConservedExact();
	int ProjectToMassConservedShift();
	int ProjectToMassConservedPSD();

	int IsContinuousProjectionOf(piecewise_constant& f, double* diag, double* subdiag);
	int IsContinuousProjectionOf(piecewise_linear& f, double* diag, double* subdiag);

	piecewise_constant DiscontinuosRepresentation();
};

int fill_diagonal_part_of_Gramm_matrix(int n, double* grid, double* diag);
int fill_subdiagonal_part_of_Gramm_matrix(int n, double* grid, double* subdiag);
