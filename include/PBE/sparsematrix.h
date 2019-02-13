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
   
#ifndef SPARSEMATRIX
#define SPARSEMATRIX

#include "basictools.h"

typedef struct _sparsematrix sparsematrix;
typedef sparsematrix *psparsematrix;
typedef const sparsematrix *pcsparsematrix;

struct _sparsematrix {
	int rows;
	int cols;
	int cmax; // maximum possible number of entries for each column

	int* nzindices; // indices of rows of nonzero entries for each column 
					// (cmax * cols dimensional array). If -1 then 
	double* data; // values of nonzero entries (cmax * cols dimensional array)
};

psparsematrix new_sparsematrix(int r, int c, int cmax);
void del_sparsematrix(psparsematrix s);

void matrix_times_sparsematrix(psparsematrix s, double* A, int A_rows, double* B);

#endif
