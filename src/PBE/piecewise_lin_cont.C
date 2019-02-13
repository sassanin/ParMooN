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
   
#include "piecewise_lin_cont.h"
#include "basictools.h"
//#include "timer.h"

double piecewise_linear_continuous::integral()
{
	double integral = 0;

	for(int i = 0; i < n; i++)
		integral += (coeff[i] + coeff[i + 1]) * (grid[i + 1] - grid[i]) / 2.0;

	return integral;
}

double piecewise_linear::integral()
{
	double integral = 0;

	for(int i = 0; i < n; i++)
		integral += (end_value[i] + start_value[i]) * (grid[i + 1] - grid[i]) / 2.0;

	return integral;
}

double piecewise_constant::integral()
{
	double integral = 0;

	for(int i = 0; i < n; i++)
		integral += value[i] * (grid[i + 1] - grid[i]);

	return integral;
}

double piecewise_linear_continuous::L1()
{
	double integral = 0;

	for(int i = 0; i < n; i++)
		integral += (fabs(coeff[i]) + fabs(coeff[i + 1])) * (grid[i + 1] - grid[i]) / 2.0;

	return integral;
}

double piecewise_linear_continuous::mass()
{
	double mass = 0;

	for(int i = 0; i < n; i++)
	{
		double a = (coeff[i + 1] - coeff[i]) / (grid[i + 1] - grid[i]);
		double b = (coeff[i] * grid[i + 1] - coeff[i + 1] * grid[i]) / (grid[i + 1] - grid[i]);
		mass += a * (cube(grid[i + 1]) - cube(grid[i])) / 3.0 +
				b * (sqr(grid[i + 1]) - sqr(grid[i])) / 2.0;
	}

	return mass;
}

int piecewise_linear_continuous::ProjectToMassConservedExact()
{
	// Exact projection
	double m3 = 3 * (*this).mass();

	static double faktor = 1. / (cube(grid[n]) - cube(grid[0]));

	for(int i = 0; i <= n; i++)
		coeff[i] -= m3 * grid[i] * faktor;

/*	(*this).coeff[n] = 0;
	double m = -(*this).mass(), local_mass = 0.0, c, a, b;

	int i, j;

	for(i = n; i > 1 && coeff[i-1] >= 0; i--);

	for(int j = i; j < n; j++) {
		a = -1 / (grid[j + 1] - grid[j]);
		b = ((n-j) * grid[j + 1] - (n-j-1) * grid[j]) / (grid[j + 1] - grid[j]);
		local_mass += a * (cube(grid[j + 1]) - cube(grid[j])) / 3.0 +
					  b * (sqr(grid[j + 1]) - sqr(grid[j])) / 2.0;
	}
	a = (n-i) / (grid[i] - grid[i-1]);
	b = -(n-i) * grid[i-1] / (grid[i] - grid[i - 1]);
	local_mass += a * (cube(grid[i]) - cube(grid[i-1])) / 3.0 +
				  b * (sqr(grid[i]) - sqr(grid[i-1])) / 2.0;

	c = m / local_mass;
	
	for(int j = i; j < n; j++) {
		coeff[j] += (n - j) * c;
	} */
	
	return 0;
}

int piecewise_linear_continuous::ProjectToMassConservedShift()
{
//	 (*this).coeff[n] = 0;
	double m2 = 2 * (*this).mass();

	static double faktor = 1. / (sqr(grid[n]) - sqr(grid[0]));

	for(int i = 0; i <= n; i++)
		coeff[i] -= m2 * faktor;

	return 0;
}

int fill_diagonal_part_of_Gramm_matrix(int n, double* grid, double* diag)
{
	for(int i = 0; i <= n; i++)
	{
		diag[i] = 0;
		if(i < n)
			diag[i] += (grid[i + 1] - grid[i]) / 3.0;

		if(i > 0)
			diag[i] += (grid[i] - grid[i - 1]) / 3.0;
	}
	return 0;
}

int fill_subdiagonal_part_of_Gramm_matrix(int n, double* grid, double* subdiag)
{
	for(int i = 0; i < n; i++)
	{
		subdiag[i] = (grid[i + 1] - grid[i]) / 6.0;
	}
	return 0;
}

int piecewise_linear_continuous::IsContinuousProjectionOf(piecewise_constant& f, double* diag, double* subdiag)
{
	int n = f.n;
/*	static double* rhs = new double[n + 1];

	rhs[0] = 3 * f.value[0] * subdiag[0];
	rhs[n] = 3 * f.value[n-1] * subdiag[n-1];
	for(int i = 1; i < n; i++) {
		rhs[i] = 3 * f.value[i] * subdiag[i] +
				 3 * f.value[i-1] * subdiag[i-1];
	}

	solve_tridiagonal_system(diag, subdiag, subdiag, rhs, (*this).coeff, n + 1); */

	(*this).coeff[0] = 0.5 * f.value[0];
	(*this).coeff[n] = 0.5 * f.value[n-1];
	for(int i = 1; i < n; i++) {
		(*this).coeff[i] = 0.5 * (f.value[i] + f.value[i-1]);
	}

	return 0;
}

int piecewise_linear_continuous::IsContinuousProjectionOf(piecewise_linear& f, double* diag, double* subdiag)
{
	int n = f.n;
	static double* rhs = new double[n + 1];

/*	rhs[0] = (2 * f.start_value[0] + f.end_value[0]) * subdiag[0];
	rhs[n] = (2 * f.end_value[n-1] + f.start_value[n-1]) * subdiag[n-1];
	for(int i = 1; i < n; i++) {
		rhs[i] = (2 * f.start_value[i] + f.end_value[i]) * subdiag[i] +
				 (2 * f.end_value[i-1] + f.start_value[i-1]) * subdiag[i-1];
	}

	solve_tridiagonal_system(diag, subdiag, subdiag, rhs, (*this).coeff, n + 1); */

	(*this).coeff[0] = 0.5 * f.start_value[0];
	(*this).coeff[n] = 0.5 * f.end_value[n-1];
	for(int i = 1; i < n; i++) {
		(*this).coeff[i] = 0.5 * (f.start_value[i] + f.end_value[i-1]);
	}

	return 0;
}

int piecewise_linear_continuous::multiplication_with_function (double* f, piecewise_constant& w)
{
	for(int i = 0; i < n; i++)
		w.value[i] = 0.5 * (coeff[i] + coeff[i+1]) * f[i];
	return 0;
}

int piecewise_linear_continuous::multiplication_with_function (double* f, piecewise_linear& w)
{
	for(int i = 0; i < n; i++)
	{
		w.start_value[i] = coeff[i] * f[i];
		w.end_value[i] = coeff[i+1] * f[i];
	}
	return 0;
}
