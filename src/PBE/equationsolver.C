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
   
#include "equationsolver.h"

int solveQuadraticEquation(idouble*  res, double p_dA, double p_dB, double p_dC)
{
	double a0 = p_dC / p_dA;
	double a1 = p_dB / p_dA;
	
	double delta = sqr(a1) - 4 * a0; 
	
	if (delta >= 0) { // Two real solutions
		res[0].m_dReal = (-a1 + sqrt(delta)) / 2;
		res[1].m_dReal = (-a1 - sqrt(delta)) / 2;
		
		res[0].m_dImm = 0;
		res[1].m_dImm = 0;
	} else { // Two complex solutions
		delta *= -1;
		
		res[0].m_dReal = -a1 / 2;
		res[1].m_dReal = res[0].m_dReal;

		res[0].m_dImm = sqrt(delta) / 2;
		res[1].m_dImm = -res[0].m_dImm;
	}
	
	return 0;
}

int solveQuadraticEquationReal(double*  res, int& cRealRoots, double p_dA, double p_dB, double p_dC)
{
	double a1 = p_dB / p_dA;
	double delta = sqr(a1) - 4 * p_dC / p_dA; 
	
	cRealRoots = 0;
	if (delta >= 0) {
		res[0] = (-a1 + sqrt(delta)) / 2;
		res[1] = (-a1 - sqrt(delta)) / 2;

		cRealRoots = 2;
	}
	return 0;
}

int solveQuadraticEquationReal1(double*  res, int& cRealRoots, double p_dB, double p_dC)
{
	double delta = sqr(p_dB) - 4 * p_dC; 
	
	if (delta >= 0) {
		res[0] = 0.5 * (-p_dB + sqrt(delta));
		res[1] = 0.5 * (-p_dB - sqrt(delta)) ;

		cRealRoots = 2;
	} else
		cRealRoots = 0;

	return 0;
}

int solveCubicEquation(idouble* res, double p_dA, double p_dB, double p_dC, double p_dD)
{
	#define a p_dA
	#define b p_dB
	#define c p_dC
	#define d p_dD

	a *= 3;
	double e = 1 / a * (c - sqr(b) / a);
	double f = (3 / a) * (d + (2 * cube(b)) / (3 * sqr(a)) - (b * c) / a);
	
	idouble wres[2]; 
	solveQuadraticEquation(wres, 1, f, -cube(e));
	
	idouble wr;
	for(int i = 0; i < 2; i++)
		if (!isZero(wres[i])) {
			wr = wres[i];
			break;
		}
	
	double r = modulus(wr.m_dReal, wr.m_dImm);
	double theta = atan2(wr.m_dImm, wr.m_dReal);
	
	idouble zres[3];
	
	double dCR = cubert(r);
	
	zres[0].m_dReal = dCR * cos(theta/3);
	zres[0].m_dImm = dCR * sin(theta/3);
	
	zres[1].m_dReal = dCR * cos((theta + TWOPI)/3);
	zres[1].m_dImm = dCR * sin((theta + TWOPI)/3);
	
	zres[2].m_dReal = dCR * cos((theta + FOURPI)/3);
	zres[2].m_dImm = dCR * sin((theta + FOURPI)/3);
	
	idouble yres[3];
	
	for(int i = 0; i < 3; i++) {
		if (!isZero(zres[i])) {
			double ro = sqr(zres[i].m_dReal) + sqr(zres[i].m_dImm);
			yres[i].m_dReal = zres[i].m_dReal - (zres[i].m_dReal*e) / ro;
			yres[i].m_dImm = zres[i].m_dImm  + (zres[i].m_dImm*e)  / ro;
		}
	}
	
	double T = b/a;
	
	for(int i = 0; i < 3; i++) {
		res[i].m_dReal = yres[i].m_dReal - T;
		res[i].m_dImm = yres[i].m_dImm;
	}
	
	return 0;
	
	#undef a
	#undef b
	#undef c
	#undef d
}

int solveCubicEquationReal(double* res, int& cRealRoots, double p_dA, double p_dB, double p_dC, double p_dD)
{
	#define a p_dA
	#define b p_dB
	#define c p_dC
	#define d p_dD
	
	a *= 3;
	double e = 1 / a * (c - sqr(b) / a);
	double f = (3 / a) * (d + (2 * cube(b)) / (3 * sqr(a)) - (b * c) / a);
	
	idouble wres[2];
	solveQuadraticEquation(wres, 1, f, -cube(e));
	
	idouble wr;
	for(int i = 0; i < 2; i++)
		if (!isZero(wres[i])) {
			wr = wres[i];
			break;
		}
	
	double r = modulus(wr.m_dReal, wr.m_dImm);
	double theta = atan2(wr.m_dImm, wr.m_dReal);
	
	idouble zres[3];
	
	double dCR = cubert(r);
	
	zres[0].m_dReal = dCR * cos(theta/3);
	zres[0].m_dImm = dCR * sin(theta/3);
	
	zres[1].m_dReal = dCR * cos((theta + TWOPI)/3);
	zres[1].m_dImm = dCR * sin((theta + TWOPI)/3);
	
	zres[2].m_dReal = dCR * cos((theta + FOURPI)/3);
	zres[2].m_dImm = dCR * sin((theta + FOURPI)/3);
	
	idouble yres[3];
	double S = e;
	
	for(int i = 0; i < 3; i++) {
		if (!isZero(zres[i])) {
			double ro = sqr(zres[i].m_dReal) + sqr(zres[i].m_dImm);
			yres[i].m_dReal = zres[i].m_dReal - (zres[i].m_dReal*e) / ro;
			yres[i].m_dImm = zres[i].m_dImm  + (zres[i].m_dImm*e)  / ro;
		}
	}
	
	double T = b/a;
	
	cRealRoots = 0;
	for(int i = 0; i < 3; i++)
		if(isReal(yres[i]))
			res[cRealRoots++] = yres[i].m_dReal - T;
	
	return 0;
	
	#undef a
	#undef b
	#undef c
	#undef d
}

int solveCubicEquationReal1(double* res, int& cRealRoots, double p_dB, double p_dC, double p_dD)
{
	#define one_over_a 0.33333333333333 // In any case 'a' should be initialised as 3
	#define b p_dB
	#define c p_dC
	#define d p_dD
	
	double e = one_over_a * (c - one_over_a * sqr(b));
	double f = d + (2 * cube(b)) / 27 - one_over_a * b * c;
	
	idouble wres[2];
	solveQuadraticEquation(wres, 1, f, -cube(e));
	
	idouble wr;
	for(int i = 0; i < 2; i++)
		if (!isZero(wres[i])) {
			wr = wres[i];
			break;
		}
	
	double r = modulus(wr.m_dReal, wr.m_dImm);
	double theta = atan2(wr.m_dImm, wr.m_dReal);
	
	idouble zres[3];
	
	double dCR = cubert(r);
	
	zres[0].m_dReal = dCR * cos(theta/3);
	zres[0].m_dImm = dCR * sin(theta/3);
	
	zres[1].m_dReal = dCR * cos((theta + TWOPI)/3);
	zres[1].m_dImm = dCR * sin((theta + TWOPI)/3);
	
	zres[2].m_dReal = dCR * cos((theta + FOURPI)/3);
	zres[2].m_dImm = dCR * sin((theta + FOURPI)/3);
	
	idouble yres[3];
	double S = e;
	
	for(int i = 0; i < 3; i++) {
		if (!isZero(zres[i])) {
			double ro = sqr(zres[i].m_dReal) + sqr(zres[i].m_dImm);
			yres[i].m_dReal = zres[i].m_dReal - (zres[i].m_dReal*e) / ro;
			yres[i].m_dImm = zres[i].m_dImm  + (zres[i].m_dImm*e)  / ro;
		}
	}
	
	double T = one_over_a * b;
	
	cRealRoots = 0;
	for(int i = 0; i < 3; i++)
		if(isReal(yres[i]))
			res[cRealRoots++] = yres[i].m_dReal - T;
	
	return 0;
	
	#undef a
	#undef b
	#undef c
	#undef d
}

// Reference: http://mathforum.org/dr.math/faq/faq.cubic.equations.html
int solveQuarticEquation(idouble* res, double p_dA, double p_dB, double p_dC, double p_dD, double p_dE)
{
	double a = p_dB / p_dA;
	double b = p_dC / p_dA;
	double c = p_dD / p_dA;
	double d = p_dE / p_dA;
	
	double quarter_of_a = 0.25 * a;
	double e = b - 6 * sqr(quarter_of_a);
	double f = c + cube(a) / 8 - (a * b) / 2;
	double g = d - 3 * sqr(sqr(quarter_of_a)) + sqr(quarter_of_a) * b - quarter_of_a * c;
	
	if (f == 0 && g == 0 && e == 0) {
		res[0].m_dReal = 0; res[0].m_dImm = 0;
		res[1] = res[2] = res[3] = res[0];
	} else if(g == 0) {
		res[0].m_dReal = -quarter_of_a ;
		res[0].m_dImm = 0;
		
		idouble cres[3];
		solveCubicEquation(cres, 1, 0, e, f);
		
		res[1].m_dReal = cres[0].m_dReal - quarter_of_a;
		res[1].m_dImm = cres[0].m_dImm;
		res[2].m_dReal = cres[1].m_dReal - quarter_of_a;
		res[2].m_dImm = cres[1].m_dImm;
		res[3].m_dReal = cres[2].m_dReal - quarter_of_a;
		res[3].m_dImm = cres[2].m_dImm;
	} else if (f == 0) {
		idouble qres[2];
		solveQuadraticEquation(qres, 1, e, g);
		
		idouble z1 = isqrt(qres[0]);
		idouble z2 = isqrt(qres[1]);

		res[0].m_dReal = z1.m_dReal - quarter_of_a;
		res[0].m_dImm = z1.m_dImm;
		res[1].m_dReal = -z1.m_dReal - quarter_of_a;
		res[1].m_dImm = -z1.m_dImm;
		res[2].m_dReal = z2.m_dReal - quarter_of_a;
		res[2].m_dImm = z2.m_dImm;
		res[3].m_dReal = -z2.m_dReal - quarter_of_a;
		res[3].m_dImm = -z2.m_dImm;
	} else {
		double m = 2*e;
		double n = sqr(e) - 4*g;
		double o = -sqr(f);
		
		idouble cres[3];
		solveCubicEquation(cres, 1, m, n, o);
		
		double ires;
		
		for(int i = 0; i < 3; i++) {
			if (isReal(cres[i])) {
				ires = abs(cres[i].m_dReal);
				break;
			}
		}

		double h = sqrt(ires);
		double j = (e + ires - f/h) / 2;

		idouble e1res[2];
		solveQuadraticEquation(e1res, 1, h, j);
		idouble e2res[2];
		solveQuadraticEquation(e2res, 1, -h, g/j);
		
		res[0].m_dReal = e1res[0].m_dReal - quarter_of_a;
		res[0].m_dImm = e1res[0].m_dImm;
		res[1].m_dReal = e1res[1].m_dReal - quarter_of_a;
		res[1].m_dImm = e1res[1].m_dImm;
		res[2].m_dReal = e2res[0].m_dReal - quarter_of_a;
		res[2].m_dImm = e2res[0].m_dImm;
		res[3].m_dReal = e2res[1].m_dReal - quarter_of_a;
		res[3].m_dImm = e2res[1].m_dImm;
	}
	
	return 0;
}

int solveQuarticEquationReal(double* res, int& cRealRoots, double p_dA, double p_dB, double p_dC, double p_dD, double p_dE)
{
	double a = p_dB / p_dA, b = p_dC / p_dA, c = p_dD / p_dA, d = p_dE / p_dA;
	int i;
	
	double quarter_of_a = 0.25 * a;
	double e = b - 6 * sqr(quarter_of_a);
	double f = c + cube(a) / 8 - (a * b) / 2;
	double g = d - 3 * sqr(sqr(quarter_of_a)) + sqr(quarter_of_a) * b - quarter_of_a * c;
	
	if (f == 0 && g == 0 && e == 0) {
		res[0] = 0;
		res[1] = res[2] = res[3] = res[0];
		cRealRoots = 4;
	} else if(g == 0) {
		res[0] = -quarter_of_a ;
		
		double cres[3];
		int croots;
		solveCubicEquationReal(cres, croots, 1, 0, e, f);
		
		for(i = 0; i < croots; i++)
				res[i + 1] = cres[i] - quarter_of_a;

		cRealRoots = croots + 1;
	} else if (f == 0) {
		double qres[2];
		int croots;
		solveQuadraticEquationReal(qres, croots, 1, e, g);
		cRealRoots = 0;
		
		if(croots > 0) {
			if(qres[0] >= 0) {
				double z1 = sqrt(qres[0]);
				res[cRealRoots++] = z1 - quarter_of_a;
				res[cRealRoots++] = -z1 - quarter_of_a;
			}
			if(qres[1] >= 0) {
				double z2 = sqrt(qres[1]);
				res[cRealRoots++] = z2 - quarter_of_a;
				res[cRealRoots++] = -z2 - quarter_of_a;
			}
		}
	} else {
		double m = 2*e;
		double n = sqr(e) - 4*g;
		double o = -sqr(f);
		
		double cres[3]; int croots;
		solveCubicEquationReal(cres, croots, 1, m, n, o);
		
		double ires = abs(cres[0]);
		
		double h = sqrt(ires);
		double j = (e + ires - f/h) / 2;

		solveQuadraticEquationReal(cres, croots, 1, h, j);
		cRealRoots = 0;
		for(i = 0; i < croots; i++)
			res[cRealRoots++] = cres[i] - quarter_of_a;

		solveQuadraticEquationReal(cres, croots, 1, -h, g/j);
		for(i = 0; i < croots; i++)
			res[cRealRoots++] = cres[i] - quarter_of_a;
	}
	
	return 0;
}
