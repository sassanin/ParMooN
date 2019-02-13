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

static double C_r;

/***************************************************************************/
/***************************************************************************/

int Breakage::is_on_boundary(int ir, double* vr) {
	if(euclidean_norm(vr, 3) > ZERO_THRESHOLD)
		return 0;

	if(ir % nx == 0)
		return 1;
	ir -= nx * (int)(ir / nx);
	if(ir % ny == 0)
		return 1;
	ir -= ny * (int)(ir / ny);
	return (ir % nz == 0);
}

double Breakage::get_value_at_x(double* f, double x0) {
	int i = get_uniform_index(pow(xmax, 3) / na, pow(x0, 3));

	if(i == na - 1)	return f[na - 1];
	if(i >= na || i < 0) return 0;

	double alpha = (f[i + 1] - f[i]) / (a[i + 1] - a[i]);
	double beta = f[i] - alpha * a[i];
	return alpha * x0 + beta;
}

void Breakage::fill_g(double& CE) {
	for(int i = 0; i < na; i++)
		if(CE != 0)
			g[i] = pow(2 * a[i], 0.75) / CE;
		else
			g[i] = 1;
}
void Breakage::fill_lg(double& CE) {
	int i, c = 0;

	for(i = 0; i < na; i++) {
		c = get_uniform_index(xmin, g[i]);
		if(c >= na - 1) 
			lg[i] = na - 1;
		else if(c < 1)
			lg[i] = 1;
		else
			lg[i] = c + 1;
	}

	lg[na] = na - 1;
}

void Breakage::fill_freq_integral_coeff() {
	int i;
	freq_integral_coeff = allocate_vector(na);
	clear_vector(freq_integral_coeff, na);

	#define A a[i-1]
	#define B a[i]

	for(i = 1; i < na; i++) {
		freq_integral_coeff[i-1] -= 0.2 * (A*A*A*A + A*A*A*B + A*A*B*B + A*B*B*B + B*B*B*B) -
									0.25 * B * (A*A*A + A*A*B + A*B*B + B*B*B);
		freq_integral_coeff[i] += 0.2 * (A*A*A*A + A*A*A*B + A*A*B*B + A*B*B*B + B*B*B*B) -
									0.25 * A * (A*A*A + A*A*B + A*B*B + B*B*B);
	}
	#undef A
	#undef B
}

// Definite integral of (x^4 * f(x))/(x-c) on [a,b], where f is linear, f(a)=fa, f(b)=fb
double Breakage::definite_integral(double& fa, double& fb, double& a, double& b, double c, double& x1, double& x2) {
	if(a == b) return 0;

	double alpha = (fb - fa) / (b - a);
	double beta = alpha * a - fa;

	double d4 = ipow(x2, 4) - ipow(x1, 4);
	double d3 = ipow(x2, 3) - ipow(x1, 3);
	double d2 = ipow(x2, 2) - ipow(x1, 2);

	return alpha * d4 / 5 + (alpha * c - beta) * 
		   (d4 / 4 + c * d3 / 3 + ipow(c, 2) * d2 / 2 + ipow(c, 3) * (x2 - x1) + ipow(c, 4) * (log(fabs(x2 - c)) - log(fabs(x1 - c))));
}

double Breakage::frequency(double* input) {
	double s = 0;
	for(int i = 0; i < na; i++)
		s += input[i] * freq_integral_coeff[i];
	return k_V * s;
}

int Breakage::solvequarticequation(double* roots, int& cRealRoots, double C_r, double cube_x) {
	double a = -1 / C_r, b = 0, c = 0, d = cube_x / C_r;
	int i;
	
	double quarter_of_a = 0.25 * a;
	double sqr_quarter_of_a = sqr(quarter_of_a);
	double e = -6 * sqr_quarter_of_a;
	double f = cube(a) / 8;
	double g = d - 3 * sqr(sqr_quarter_of_a);
	
	if(g == 0) {
		roots[0] = -quarter_of_a ;
		
		double cres[3];
		int croots;
		solveCubicEquationReal1(cres, croots, 0, e, f);
		
		for(i = 0; i < croots; i++)
				roots[i + 1] = cres[i] - quarter_of_a;

		cRealRoots = croots + 1;
	} else {
		double cres[3]; int croots;
		solveCubicEquationReal1(cres, croots, 2 * e, sqr(e) - 4 * g, -sqr(f));
		
		double h = sqrt(cres[0]);
		double j = (e + cres[0] - f/h) / 2;

		double eres[2];
		solveQuadraticEquationReal1(eres, croots, h, j);
		cRealRoots = 0;
		for(i = 0; i < croots; i++)
			roots[cRealRoots++] = eres[i] - quarter_of_a;

		solveQuadraticEquationReal1(eres, croots, -h, g/j);
		for(i = 0; i < croots; i++)
			roots[cRealRoots++] = eres[i] - quarter_of_a;
	}

	return 0;
}

void Breakage::breakage_plus(double* input, double* output, double* v) {
	int i, j; 
	double b = 0;
    double E = pow(0.5 * ro_d * k_V * ipow(euclidean_norm(v, 3), 2), 1.0/3.0);
	double CE = C * E;

	if(CE < ZERO_THRESHOLD)
		return;

	double factor = pow(0.5, 0.25) * ipow(CE, 3) / k_V;
	double y_c = pow(2 * L_min, 0.75) / CE;

	double int1D = 0, si = 0;
	fill_g(CE);
	fill_lg(CE);
	for(i = na - 1; i >= 0; i--) {
		for(j = lg[i]; j < lg[i + 1]; j++)
			int1D += definite_integral(input[lg[j]], input[lg[j] + 1], a[lg[j]], a[lg[j] + 1], y_c, a[lg[j]], a[lg[j] + 1]);

		si = definite_integral(input[lg[i] - 1], input[lg[i]], a[lg[i] - 1], a[lg[i]], y_c, g[i], a[lg[i]]);
		output[i] += (int1D + si) * factor * pow(a[i], a[i] == 0 ? 1 : -3.25);
	}
	
	// Second part of P(x,y) (Dirac delta)
	#define yj roots[j]

	C_r = 2 * ipow(CE, 4) / (3 * k_V);
	static double roots[4]; int croots = 0;

	for(i = 0; i < na; i++) {
		solvequarticequation(roots, croots, C_r, ipow(a[i], 3));
		for(int j = 0; j < croots; j++)
			if(yj > g[i] && yj < xmax) {
				output[i] += get_value_at_x(input, yj) * 3 * pow(sqr(C_r * yj - 1),	1.0/3.0)
							 / (4 * C_r * yj - 3);
			}
	}

	#undef yj
}

void Breakage::breakage_minus(double* input, double* output) {
	for(int i = 0; i < na; i++)
		output[i] -= input[i];
}

void Breakage::apply_breakage(double* input, double* output, double* v) {
	double* f = allocate_vector(na);
	double* rhs = allocate_vector(na);
	clear_vector(rhs, na);
	int Nr = nx * ny * nz, i, j; double b;

	for(i = 0; i < Nr; i++) {
		for(j = 0; j < na; j++)
			f[j] = input[i + j * Nr];

		if(!is_on_boundary(i, v + 3 * i)) {
			b = frequency(f);
			breakage_plus(f, rhs, v + 3 * i);
			breakage_minus(f, rhs);
			scale_vector(rhs, na, b);
		} else
			clear_vector(rhs, na);

		for(j = 0; j < na; j++)
			output[i + j * Nr] += rhs[j];
	}

	free_vector(f);
	free_vector(rhs);
}
