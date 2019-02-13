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
   
#include "fullmatrix.h"
#include <mkl.h>
#include <memory.h>
#include <malloc.h>

pfullmatrix new_fullmatrix(int r, int c)
{
	pfullmatrix f;

	f = (pfullmatrix)malloc(sizeof(fullmatrix));
	f->rows = r;
	f->cols = c;
	f->data = allocate_vector(r * c);

	return f;
}

void del_fullmatrix(pfullmatrix f)
{
	free_vector(f->data);
	f = 0;
}

void fullmatrix_times_matrix(pfullmatrix f, EvalMode mode, double* B, int n, double* C, double alpha, double beta)
{
	matrix_times_matrix(f->data, "n", f->rows, f->cols, B, "n", f->cols, n, C, alpha, beta);
}

void fullmatrix_times_vector(pfullmatrix f, EvalMode mode, double* b, double* c, double alpha, double beta)
{
	double* Ab = NULL;
	int eins = 1;
	switch(mode)
	{
	case EM_DEFAULT:
		matrix_times_vector(f->data, "n", f->rows, f->cols, b, c, alpha, beta);
		break;
	case EM_TRANSPOSED:
		matrix_times_vector(f->data, "T", f->rows, f->cols, b, c, alpha, beta);
		break;
	case EM_LOWER_TRIANGULAR:
		Ab = allocate_vector(f->rows);
		memcpy(Ab, b, f->rows * sizeof(double));
		triangular_matrix_times_vector(f->data, "L", f->rows, Ab);
		scale_vector(c, f->rows, beta);
		daxpy(&f->rows, &alpha, Ab, &eins , c, &eins);
		break;
	case EM_UPPER_TRIANGULAR:
		Ab = allocate_vector(f->rows);
		memcpy(Ab, b, f->rows * sizeof(double));
		triangular_matrix_times_vector(f->data, "U", f->rows, Ab);
		scale_vector(c, f->rows, beta);
		daxpy(&f->rows, &alpha, Ab, &eins , c, &eins);
		break;
	}
	if(Ab)
		free_vector(Ab);
}

