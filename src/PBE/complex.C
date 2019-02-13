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
   
#include "complex.h"

idouble add(idouble a, idouble b) {
	idouble c; 
	c.m_dReal = a.m_dReal + b.m_dReal; 
	c.m_dImm = a.m_dImm + b.m_dImm; 
	return c; 
}

idouble subtract(idouble a, idouble b) {
	idouble c; 
	c.m_dReal = a.m_dReal - b.m_dReal; 
	c.m_dImm = a.m_dImm - b.m_dImm;
	return c; 
}

bool isZero(Complex p_dValue) {
	return (isZero(p_dValue.m_dReal) && isZero(p_dValue.m_dImm));
}
bool isReal(idouble a) {
	return isZero(a.m_dImm); 
}

idouble isqrt(Complex p_pIValue) {
	#define a p_pIValue.m_dReal
	#define b p_pIValue.m_dImm
	
	Complex res;
	
	double r = modulus(a, b);
	
	int eps;
	if (b != 0)
		eps = (int)sgn(b);
	else
		eps = 1;
	
	double alpha = (1/SQRT2);
	double dSqrtA = sqrt(r + a);
	double dSqrtB = sqrt(r - a);
	
	res.m_dReal = alpha * eps * dSqrtA;
	res.m_dImm = alpha * dSqrtB;
	
	return res;
	
	#undef a
	#undef b
}
