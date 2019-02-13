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
   
#include "sparsematrix.h"
#include <malloc.h>

psparsematrix new_sparsematrix(int r, int c, int cmax)
{
	psparsematrix s;
	s = (psparsematrix)malloc(sizeof(sparsematrix));

	s->rows = r;
	s->cols = c;
	s->cmax = cmax;

	s->nzindices = (int*)malloc(cmax * c * sizeof(int));
	s->data = allocate_vector(cmax * c);

	return s;
}

void del_sparsematrix(psparsematrix s)
{
	s->rows = 0;
	s->cols = 0;
	s->cmax = 0;
	free((void*)s->nzindices);
	free_vector(s->data);
	s = 0;
}

void matrix_times_sparsematrix(psparsematrix s, double* A, int A_rows, double* B)
{
	int i, j, h;
	for(i = 0; i < s->cols; i++)
		for(j = 0; j < A_rows; j++)
		{
			B[i * A_rows + j]  = 0;
			for(h = 0; h < s->cmax; h++)
			{
				if(s->nzindices[s->cmax * i + h] >= 0)
					B[i * A_rows + j] += s->data[s->cmax * i + h] * A[s->nzindices[s->cmax * i + h] * A_rows + j];
			}
		}
}

