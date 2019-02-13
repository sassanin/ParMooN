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
   
#include "rkmatrix.h"

// C = alpha * op(R) * B + beta * C
void rkmatrix::rkmatrix_times_matrix(EvalMode mode, double* B, double* C, int n, double alpha, double beta) {
	double* bTB = allocate_vector(rk * n);

	switch(mode) {
	case EM_DEFAULT:
		matrix_times_matrix(b, "T", cols, rk, B, "N", cols, n, bTB, alpha, 0.0);
		matrix_times_matrix(a, "N", rows, rk, bTB, "N", rk, n, C, alpha, beta);
		break;
	case EM_TRANSPOSED:
		matrix_times_matrix(a, "T", rows, rk, B, "N", rows, n, bTB, alpha, 0.0);
		matrix_times_matrix(b, "N", cols, rk, bTB, "N", rk, n, C, alpha, beta);
		break;
	case EM_UPPER_TRIANGULAR:
		break;
	case EM_LOWER_TRIANGULAR:
		break;
	}
	free_vector(bTB);
}

// y = alpha * op(R) * x + beta * y
void rkmatrix::rkmatrix_times_vector(EvalMode mode, double* x, double* y, double alpha, double beta) {
	int i, j;
	int eins = 1;
	double tmp = 0;

	switch(mode) {
	case EM_DEFAULT:
		if(beta != 1.0)	{
			if(beta == 0.0) clear_vector(y, rows);
			else scale_vector(y, rows, beta);
		}

		for(i = 0; i < rk; i++)	{
			tmp = alpha * ddot_(&cols, x, &eins, &(b[i * cols]), &eins);
			daxpy_(&rows, &tmp, &(a[i * rows]), &eins, y, &eins);
		}
		break;

	case EM_TRANSPOSED:
		if(beta != 1.0)	{
			if(beta == 0.0) clear_vector(y, cols);
			else scale_vector(y, cols, beta);
		}

		for(i = 0; i < rk; i++) {
			tmp = alpha * ddot_(&rows, x, &eins, &(a[i * cols]), &eins);
			daxpy_(&cols, &tmp, &(b[i * rows]), &eins , y, &eins);
		}
		break;

	case EM_UPPER_TRIANGULAR:
		assert(rows == cols);
		if(beta != 1.0)	{
			if(beta == 0.0)	clear_vector(y, rows);
			else scale_vector(y, rows, beta);
		}

		for(i = 0; i < rk; i++)	{
			tmp = 0;
			for(j = cols - 1; j >= 0; j--) {
				tmp += b[i * cols + j] * x[j];
				y[j] += a[i * rows + j] * tmp;
			}
		}
		break;

	case EM_LOWER_TRIANGULAR:
		assert(rows == cols);
		if(beta != 1.0)	{
			if(beta == 0.0)	clear_vector(y, rows);
			else scale_vector(y, rows, beta);
		}

		for(i = 0; i < rk; i++) {
			tmp = 0;
			for(j = 0; j < cols; j++) {
				tmp += b[i * cols + j] * x[j];
				y[j] += a[i * rows + j] * tmp;
			}
		}

		break;
	}
}

// Sinfular value decomposition of low-rank matrix
// via QR decompositions of its components
int rkmatrix::rk_svd(double* u, double* sifma, double* v) {
	double *atmp = allocate_vector(rows * rk);
	double *btmp = allocate_vector(cols * rk);
	int i, j, h = 0;
	double* taua = allocate_vector(rows + cols);
	double* taub = allocate_vector(rows + cols);
	double* rarb = allocate_vector(rk * rk);
	int cwork = rows + cols + 10;
	double* work = allocate_vector(cwork);
	int info;
	double* ru = allocate_vector(rk * rk);
	double* rv = allocate_vector(rk * rk);

	if(rk == 0 || rows == 0 || cols == 0) {
		printf("(function rk_svd)\nRank k matrix has invalid data.\nOne of the parameters is null.\n");
		return -1;
	}
	if(rk > rows || rk > cols) {
		printf("(function rk_svd)\nRank k matrix has invalid data.\nRank is more than one of the dimensions.\n");
		return -1;
	}

	memcpy(atmp, a, rows * rk * sizeof(double));
	memcpy(btmp, b, cols * rk * sizeof(double));

	dgeqrf_(&rows, &rk, atmp, &rows, taua, work, &cwork, &info);
	if(info != 0) {
		printf("(function rk_svd)\nError %d on QR decomposition of r->a\n", info);
		return -1;
	}

	dgeqrf_(&cols, &rk, btmp, &cols, taub, work, &cwork, &info);
	if(info != 0) {
		printf("(function rk_svd)\nError %d on QR decomposition of r->b\n", info);
		return -1;
	}

	for(i = 0; i < rk; i++)
		for(j = 0; j < rk; j++)	{
			rarb[i + rk * j] = 0;
			for(h = max(i, j); h < rk; h++)
				rarb[i + rk * j] += atmp[i + h * rows] * btmp[j + h * cols];
		}

	cwork = rows + 10;
	free_vector(work);
	work = allocate_vector(cwork);
	dorgqr_(&rows, &rk, &rk, atmp, &rows, taua, work, &cwork, &info);
	if(info != 0) {
		printf("(function rk_svd)\nError %d on gettinf Q from QR decomposition of r->a\n", info);
		return -1;
	}

	cwork = cols + 10;
	free_vector(work);
	work = allocate_vector(cwork);
	dorgqr_(&cols, &rk, &rk, btmp, &cols, taub, work, &cwork, &info);
	if(info != 0) {
		printf("(function rk_svd)\nError %d on gettinf Q from QR decomposition of r->b\n", info);
		return -1;
	}

	cwork = 10 * rk * rk;
	free_vector(work);
	work = allocate_vector(cwork);

	dgesvd_("A", "A", &rk, &rk, rarb, &rk, sifma, ru, &rk, rv, &rk, work, &cwork, &info);

	if(info != 0) {
		printf("(function rk_svd)\nError %d on SVD of ra * rb\n", info);
		return -1;
	}

	if(u != NULL)
		matrix_times_matrix(atmp, "N", rows, rk, ru, "N", rk, rk, u, 1, 0);
	if(v != NULL)
		matrix_times_matrix(btmp, "N", rows, rk, rv, "T", rk, rk, v, 1, 0);

	free_vector(work);
	free_vector(rarb);
	free_vector(btmp);
	free_vector(atmp);
	free_vector(taub);
	free_vector(taua);
	free_vector(ru);
	free_vector(rv);

	return 0;
}

double rkmatrix::get_frobenius_norm() {
	int i; 
	double s = 0;
	double* ata = allocate_vector(rk * rk);
	double* btb = allocate_vector(rk * rk);

	matrix_times_matrix(a, "T", rows, rk, a, "N", rows, rk, ata, 1.0, 0.0);
	matrix_times_matrix(b, "T", cols, rk, b, "N", cols, rk, btb, 1.0, 0.0);
	for(i = 0; i < rk * rk; i++)
		s += ata[i] * btb[i];

	return sqrt(s);
}
