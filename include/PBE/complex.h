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
   
#ifndef MPIMIS_COMPLEX
#define MPIMIS_COMPLEX

#include "complextools.h"
#include "basictools.h"

typedef class Complex complex;
typedef complex* pcomplex;
typedef complex idouble;
typedef idouble* pidouble;

class Complex
{
public:
	Complex() { m_dReal = 0; m_dImm = 0; };
	~Complex() {};

	double m_dReal;
	double m_dImm;
};

idouble add(idouble a, idouble b);
idouble subtract(idouble a, idouble b);  

bool isZero(Complex p_dValue);
bool isReal(idouble a);

Complex isqrt(Complex p_pIValue);

#endif // MPIMIS_COMPLEX
