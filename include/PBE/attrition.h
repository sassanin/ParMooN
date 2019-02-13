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
   
#ifndef ATTRITION
#define ATTRITION

#include "sparsematrix.h"
//#include "projection.h"

typedef struct _attrition attrition;
typedef attrition* pattrition;

struct _attrition {
	int nx;

	int nr;
	psparsematrix s;

	function_3D3D* velocity;
};

#define CMAX 9

pattrition new_attrition(int nx, int nr, function_3D3D velocity);

void del_attrition(pattrition m);

//int load_attrition(pattrition mo, double* internal_grid, KernelType kt, prkdata lowr, function_2D kl, function_2D kh);

int apply_attrition(pattrition mo, double* input, double* output, double t);

void generate_mass_matrix(psparsematrix M);

double get_x(int ix, int nx);

double* get_r(int ir, int nr);

int is_on_boundary(double* r);

#endif
