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
   
#ifndef PBE
#define PBE

#include "hmatrix.h"
#include "sparsematrix.h"

typedef enum {
	kt_rkmatrix = 1,
	kt_convolution = 2,
	kt_lshape = 4,
	kt_hmatrix = 8
} KernelType;

typedef struct _rkdata rkdata;
typedef rkdata* prkdata;

struct _rkdata {
	int r;
	function_1D** fx;
	function_1D** fy;
};

typedef struct _integraloperator integraloperator;
typedef integraloperator* pintegraloperator;

struct _integraloperator {
	int nx;
	double* xgrid;

	KernelType kt;
	prkmatrix rk;
	phmatrix hmatr;
	
	int nr;
	psparsematrix s;
};

pintegraloperator new_integraloperator(int nx, int nr);

void del_integraloperator(pintegraloperator m);

int load_integraloperator(pintegraloperator mo, double* internal_grid, KernelType kt, prkdata lowr, function_2D kl, function_2D kh);

int apply_integraloperator(pintegraloperator mo, double* input, double* output);

int load_hmatrix(pintegraloperator mo, function_2D k);

void generate_mass_matrix(psparsematrix M);

#endif
