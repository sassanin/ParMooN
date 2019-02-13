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
   
#include "hmatrix.h"
#include "aca.h"
#include <assert.h>
#include <malloc.h>
#include <memory.h>

phmatrix new_hmatrix(int row, int col)
{
	phmatrix h;

	h = (phmatrix)malloc(sizeof(hmatrix));
	h->r = 0x0;
	h->f = 0x0;
	h->sons = 0x0;

	h->rows = row;
	h->cols = col;
	h->l = 5;

	return h;
}

void del_hmatrix(phmatrix h)
{
	int i;
	if(h->r != 0)
		del_rkmatrix(h->r);
	if(h->f != 0)
		del_fullmatrix(h->f);
	if(h->sons != 0)
		for(i = 0; i < h->block_cols * h->block_rows; i++)
		{
			del_hmatrix(h->sons[i]);
			h->sons[i] = 0;
		}
	h = 0;
}

int hmvm(phmatrix h, double* b, double* c)
{
	int i, j;
	int bIndex;
	int cIndex;
	hmatrix** s_el = h->sons;

	if(h->rows == 0 || h->cols == 0)
		return -1;

	if(h->sons)
	{
		cIndex = 0;
		for(i = 0; i < h->block_rows; i++)
		{
			bIndex = 0;
			for(j = 0; j < h->block_cols; j++)
			{
				hmvm(s_el[i + j * h->block_rows], b + bIndex, c + cIndex);

				bIndex += s_el[i + j * h->block_rows]->cols;
			}

			assert(bIndex == h->cols);
			cIndex += s_el[i]->rows;
		}
		assert(cIndex == h->rows);
	}

	if(h->r)
		rkmatrix_times_vector(h->r, EM_DEFAULT, b, c, 1.0, 1.0);

	if(h->f)
		fullmatrix_times_vector(h->f, EM_DEFAULT, b, c, 1.0, 1.0);

	return 0;
}

int hmmm(phmatrix h, double* B, int B_cols, double* C)
{
	int i;
	for(i = 0; i < B_cols; i++)
		hmvm(h, B + i * h->cols, C + i * h->rows);

	return 0;
}

int trihmvm(phmatrix h, char* uplo, double* b, double* c)
{
	int i, j;
	int bCompute;
	int bIndex;
	int cIndex;
	hmatrix** s_el = h->sons;
	EvalMode em;

	bCompute = 0;
	if(h->rows == 0 || h->cols == 0 || h->rows != h->cols)
		return -1;

	if(h->sons)
	{
		assert(h->block_rows == h->block_cols);
		cIndex = 0;
		for(i = 0; i < h->block_rows; i++)
		{
			bIndex = 0;
			for(j = 0; j < h->block_cols; j++)
			{
				if(uplo == "L")
					bCompute = i > j ? 1 : 0;
				if(uplo == "U")
					bCompute = i < j ? 1 : 0;

				if(bCompute)
					hmvm(s_el[i + j * h->block_rows], b + bIndex, c + cIndex);
				if(i == j)
					trihmvm(s_el[i + j * h->block_rows], uplo, b + bIndex, c + cIndex);

				bIndex += s_el[i + j * h->block_rows]->cols;
			}

			assert(bIndex == h->cols);
			cIndex += s_el[i]->rows;
		}
		assert(cIndex == h->rows);
	}

	em = uplo == "U" ? EM_UPPER_TRIANGULAR : EM_LOWER_TRIANGULAR;
	if(h->r)
		rkmatrix_times_vector(h->r, em, b, c, 1.0, 1.0);

	if(h->f)
		fullmatrix_times_vector(h->f, em, b, c, 1.0, 1.0);

	return 0;
}

int trihmmm(phmatrix h, char* uplo, double* B,int B_cols, double* C)
{
	int i;
	for(i = 0; i < B_cols; i++)
		trihmvm(h, uplo, B + i * h->cols, C + i * h->rows);

	return 0;
}

int create_hmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start)
{
	int ret = 0;

	if(h->l == 0)
		h->l = 5;

	if(h->rows < h->l || h->cols < h->l)
	{
		ret = create_fullmatrix(h, entry, data, row_start, col_start);
		return ret;
	}

	h->block_rows = 2;
	h->block_cols = 2;
	h->sons = (hmatrix**)malloc(4 * sizeof(hmatrix*));
	h->sons[0] = new_hmatrix(h->rows/2, h->cols/2);
	ret += create_hmatrix(h->sons[0], entry, data, row_start, col_start);

	h->sons[1] = new_hmatrix(h->rows - h->rows/2, h->cols/2);
	ret += create_rkmatrix(h->sons[1], entry, data, row_start + h->rows/2, col_start);

	h->sons[2] = new_hmatrix(h->rows/2, h->cols - h->cols/2);
	ret += create_rkmatrix(h->sons[2], entry, data, row_start, col_start + h->cols/2);

	h->sons[3] = new_hmatrix(h->rows - h->rows/2, h->cols - h->cols/2);
	ret += create_hmatrix(h->sons[3], entry, data, row_start + h->rows/2, col_start + h->cols/2);

	return ret;
}

int create_rkmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start)
{
	int kmax, k;
	double eps;
	double *A, *B;

	eps = 1e-12;
	kmax = 20;

	A = allocate_vector(kmax * h->rows);
	B = allocate_vector(kmax * h->cols);
	assert(A != 0);
	assert(B != 0);

	k = aca_fill_block(A, B, h->rows, h->cols, row_start, col_start, entry, data, kmax, eps, HLIB_ACA_DEFAULT); 
	assert(k < kmax);
//	k = newaca_fill_block(A, B, h->rows, h->cols, row_start, col_start, entry, data, kmax, eps, HLIB_ACA_DEFAULT);
	h->r = new_rkmatrix(k, h->rows, h->cols);
	assert(h->r != 0);
	memcpy(h->r->a, A, k * h->rows * sizeof(double));
	memcpy(h->r->b, B, k * h->cols * sizeof(double));

	free_vector(A);
	free_vector(B);

	return 0;
}

int create_fullmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start)
{
	int i, j;
	h->f = new_fullmatrix(h->rows, h->cols);
	for(i = 0; i < h->cols; i++)
		for(j = 0; j < h->rows; j++)
			h->f->data[h->rows * i + j] = entry(row_start + j, col_start + i, data);

	return 0;
}

