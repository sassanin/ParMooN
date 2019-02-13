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
   
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "basictools.h"
#include "exported.h"

/**************************************************************************
****************** Norms, scalar product, orthonormalisation **************
**************************************************************************/
double scalar_product(double* a, double* b, int n) {
	int i;
	double ret = 0;

	for(i = 0; i < n; i++)
		ret += a[i] * b[i];

	return ret;
}

double euclidean_norm(double* a, int n) {
	return sqrt(scalar_product(a, a, n));
}

double frobenius_norm(double* A, int rows, int cols) {
	int i, j;
	double ret = 0;

	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			ret += A[i * cols + j] * A[i * cols + j];

	return sqrt(ret);
}

double normalize_vector(double* a, int n) {
	int i; 
	double norm = euclidean_norm(a, n);

	if(norm == 0) return 0;

	for(i = 0; i < n; i++)
		a[i] /= norm;

	return norm;
}

// Distance of the vector 'b' from the line spanned by 'a'.
// The vector 'a' must have norm 1
double vector_distance(double* a, double* b, int n) {
	double l = scalar_product(a, b, n), ret = 0;
	for(int i = 0; i < n; i++)
		ret += pow(b[i] - l * a[i], 2);

	return sqrt(ret);
}

void make_orthonormal_to_system(double* a, double* A, int k, int n) {
	int i, j;
	double* scalars = allocate_vector(k);

	for(i = 0; i < k; i++)
		scalars[i] = scalar_product(a, A + i * n, n);

	for(i = 0; i < k; i++)
		for(j = 0; j < n; j++)
			a[j] -= scalars[i] * A[i * n + j];

	free_vector(scalars);
	normalize_vector(a, n);
}

inline double sum_vector(double* v, int n) {
	int i; 
	double ret = 0;

	for(i = 0; i < n; i++)
		ret += v[i];

	return ret;
}

inline double average_vector(double* v, int n) {
	if(n == 0) return 0;

	return sum_vector(v, n) / n;
}

double min_vector(double* v, int n) {
	int i;
	double ret = v[0];

	for(i = 1; i < n; i++)
		if(ret > v[i]) ret = v[i];

	return ret;
}

double max_vector(double* v, int n) {
	int i;
	double ret = v[0];

	for(i = 1; i < n; i++)
		if(ret < v[i]) ret = v[i];

	return ret;
}

double min_abs_vector(double* v, int n) {
	int i;
	double ret = fabs(v[0]);

	for(i = 1; i < n; i++)
		if(ret > fabs(v[i])) ret = fabs(v[i]);

	return ret;
}

double max_abs_vector(double* v, int n) {
	int i;
	double ret = fabs(v[0]);

	for(i = 1; i < n; i++)
		if(ret < fabs(v[i])) ret = fabs(v[i]);

	return ret;
}

int min_vector(int* v, int n) {
	int i;
	int ret = v[0];

	for(i = 1; i < n; i++)
		if(ret > v[i]) ret = v[i];

	return ret;
}

int max_vector(int* v, int n) {
	int i;
	int ret = v[0];

	for(i = 1; i < n; i++)
		if(ret < v[i]) ret = v[i];

	return ret;
}

int min_abs_vector(int* v, int n) {
	int i;
	int ret = abs(v[0]);

	for(i = 1; i < n; i++)
		if(ret > abs(v[i])) ret = abs(v[i]);

	return ret;
}

int max_abs_vector(int* v, int n) {
	int i;
	int ret = abs(v[0]);

	for(i = 1; i < n; i++)
		if(ret < abs(v[i])) ret = abs(v[i]);

	return ret;
}

/**************************************************************************
*********************** Standard algebraic operations *********************
**************************************************************************/

// C += alpha * op(C) * op(B) + beta * C
void matrix_times_matrix(double* A, char* A_transposed, int A_rows, int A_cols, 
						 double* B, char* B_transposed, int B_rows, int B_cols, 
						 double* C, double alpha, double beta) {
	int m, n, k, lda, ldb, ldc;

	m = (A_transposed == "N" ? A_rows : A_cols);
	n = (B_transposed == "N" ? B_cols : B_rows);
	assert(A_rows + A_cols - m == B_rows + B_cols - n);
	k = A_rows + A_cols - m;

	lda = (A_transposed == "N" ? m : k);
	ldb = (B_transposed == "N" ? k : n);
	ldc = m;

	dgemm_(A_transposed, B_transposed, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

// c = alpha * op(A) * b + beta * c
void matrix_times_vector(double* A, char* A_transposed, int A_rows, int A_cols, double* b, double* c, double alpha, double beta) {
	int incx = 1, incy = 1;

	dgemv_(A_transposed, &A_rows, &A_cols, &alpha, A, &A_rows, b, &incx, &beta, c, &incy);
}

// X = A * X
void triangular_matrix_times_matrix(double* A, char* uplo, double* X, int rows, int cols) {
	double alpha = 1.0;
	int lda = rows, ldb = rows;

	dtrmm_("L", uplo, "n", "n", &rows, &cols, &alpha, A, &lda, X, &ldb);
}

// x = A * x
void triangular_matrix_times_vector(double* A, char* uplo, int n, double* x) {
	int incx = 1, lda = n;

	dtrmv_(uplo, "n", "n", &n ,A, &lda, x, &incx);
}

// SLOW!!! Can be faster
void transpose_matrix(double* A, int m, int n) {
	int i, j;
	double* B = allocate_vector(m * n);

	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			B[i * n + j] = A[j * m + i];

	memcpy(A, B, m * n * sizeof(double));
	free_vector(B);
}

void get_matrix_block(double* A, int rows, int cols, int block_row_start, int block_col_start, int block_rows, int block_cols, double* B) {
	for(int i = 0; i < block_cols; i++)
		for(int j = 0; j < block_rows; j++)
			B[i * block_rows + j] = A[(block_row_start + j) + rows * (block_col_start + i)];
}

int invert_matrix(double* A, int n) {
	int info, cwork = n * n;
	int* ipiv = (int*)malloc(n * sizeof(int));

	dgetrf_(&n, &n, A, &n, ipiv, &info);
	if(info != 0) {
		free((void*)ipiv);
		return info;
	}
	double* work = allocate_vector(cwork);
	dgetri_(&n, A, &n, ipiv, work ,&cwork, &info);
	free_vector(work);

	return info;
}

/**************************************************************************
**************** Printing structures to the screen or files ***************
**************************************************************************/

void print_vector(double* v, int n) {
	for(int i = 0; i < n; i++) printf("%.16g \t", v[i]);
	printf("\n");
}

void print_matrix(double* A, int m, int n) {
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++)
			printf("%g \t", A[i * n + j]);
		printf("\n");
	}
}

void fprint_vector(char* filename, double* v, int n) {
	FILE* file = fopen(filename, "w");
	for(int i = 0; i < n; i++)
		fprintf(file, "%g\n", v[i]);
	fclose(file);
}

void fscan_vector(char* filename, double* v, int n) {
	FILE* file = fopen(filename, "r");
	for(int i = 0; i < n; i++)
		fscanf(file, "%lf", v + i);
	fclose(file);
}

void fprint_matrix(char* filename, double* A, int rows, int cols)
{
	FILE* file = fopen(filename, "w");
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++)
			fprintf(file, "%g\t", A[i + j * rows]);
		fprintf(file, "\n");
	}
	fclose(file);
}

void fill_matrix(double* A, double (*f)(int, int, void*), int rows, int cols, void* data) {
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			A[cols * i + j] = f(i, j, data);
}

void fill_matrix(double* A, function_2D* f, int rows, int cols, double* grdpoints) {
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			A[cols * i + j] = f(grdpoints[i], grdpoints[j]);
}

/**************************************************************************
**************************** Fast calculations ****************************
**************************************************************************/

// function pow for integer degree (y)
double ipow(double x, int y) {
	double r = 1.0;

	while(y) {
		if(y % 2) r *= x;
		x *= x; y /= 2;
	}
	return r;
}

double log(double x, double a) {
	if(a <= 0) return -1;
	
	return log(x) / log(a);
}

// function pow for integer number (x) and positive integer degree (y)
int iipow(int x, int y) {
	int r = 1;

	while(y) {
		if(y % 2) r *= x;

		x *= x; y >>= 1;
	}
	return r;
}

double sgn(double x) { // Gives a sign of a number
	return ( x >= 0 ? 1 : -1);
}
double sqr(double x) {
	return x * x;
}
double cube(double x) { // Gives the cube number
	return (x * x * x);
}
double cubert(double x) { // Gives the cube root number ( cubert(-8) = -2)
	return pow(fabs(x), 1.0/3.0) * sgn(x);
}

int factorial(int n) {
	if(n == 0)
		return 1;

	return n * factorial(n - 1);
}

// Fill data
void fill_hrefined_grid(double a, double b, int cLevels, int n, double* grid) {
	int c = 1, l, i;
	double h = (b - a) / (double)(n * pow(2.0, cLevels - 1));
	grid[0] = a;

	for(i = 0; i < n; i++) {
		a += h;	grid[++c] = a; 
	}
	
	for(l = 1; l < cLevels; l++) {
		a += h;	h *= 2;
		for(i = 0; i < n / 2; i++) {
			grid[++c] = a; a += h;
		}
	}
}

void fill_uniform_grid(double a, double b, int n, double* grid) {
	int i;
	double h = (b - a) / (double)(n - 1);
	for(i = 0; i < n; i++)
		grid[i] = a + i * h;
}

void fill_exponential_data(function_1D f, int cLevels, int n, double* d) {
	int c = 0, l, i;
	double x = 0, h = 1.0 / (double)(n * pow(2.0, cLevels - 1));

	for(i = 0; i < n; i++) {
		d[c++] = (f(x) + f(x + h)) / 2; x += h;
	}

	h *= 2; 
	for(l = 1; l < cLevels; l++) {
		for(i = 0; i < n / 2; i++) {
			d[c++] = (f(x) + f(x + h)) / 2;	x += h;
		}
		h *= 2;
	}
}

void fill_uniform_data(function_1D f, int n, double* d) {
	for(int i = 0; i < n; i++)
		d[i] = (f((double)i / (double)n) + f((double)(i + 1) / (double)n)) / 2;
}

void fill_vector_from_function(function_1D f, double* grid, int n, double* v, Discretization_Scheme ds) {
	switch(ds) {
		case DS_POINTWISE:
			for(int i = 0; i < n; i++)
				v[i] = f(grid[i]);
			break;
		case DS_MEANVALUE:
			for(int i = 0; i < n - 1; i++)
				v[i] = 0.5 * (f(grid[i]) + f(grid[i + 1]));
			break;
		case DS_LINEAR_FE:
			for(int i = 1; i < n - 1; i++)
				v[i] = 0.5 * f(grid[i]) * (grid[i + 1] - grid[i - 1]);
			v[0] = 0.5 * f(grid[0]) * (grid[1] - grid[0]);
			v[n - 1] = 0.5 * f(grid[n - 1]) * (grid[n - 1] - grid[n - 2]);
			break;
		case DS_L2_CONSTANT:
			static int cQuadraturePoints = 40;
			int i, c;
			double h_quad;
			for(i = 0; i < n; i++)
			{
				v[i] = 0;
				h_quad = (grid[i + 1] - grid[i]) / cQuadraturePoints;
				for(c = 0; c < cQuadraturePoints; c ++)
					v[i] += f(grid[i] + (c + 0.5) * h_quad) / cQuadraturePoints;
			}
			break;			

	}
}
void transform_grid_m2l(double& xmin, double& xmax, int n, double* x) {
	xmin = pow(xmin, 1.0 / 3.0);
	xmax = pow(xmax, 1.0 / 3.0);
	for(int i = 0; i < n; i++)
		x[i] = pow(x[i], 1.0/3.0);
}
void transform_grid_l2m(double& xmin, double& xmax, int n, double* x) {
	xmin = ipow(xmin, 3);
	xmax = ipow(xmax, 3);
	for(int i = 0; i < n; i++)
		x[i] = ipow(x[i], 3);
}

// Finds the lowest index for h-refined grid of [0, +inf) for given x
int get_hrefined_index(double h, int n, double x) {
	if(x <= n * h)
		return get_uniform_index(h, x);

	x /= h;
	int k = (int)log(x/n, 2.0);

	return iround((k + 1) * n / 2.0 + 1 + (int)((x - iipow(2, k) * n) / iipow(2, k + 1)));
}

// Finds the lowest index for uniform grid on [0, +inf] for given x
int get_uniform_index(double h, double x) {
	return (int)(x / h);
}

// Finds the lowest index for arbitrary grid for given x
int get_index(double* x, int n, double x0) {
	// TODO: via bisection, in O(log(n)) flops
	int i = 0;
	return i;
}

// Functions for testing
void generate_orthonormal_matrix(double* Q, int k, int n) {
	srand((unsigned)time(NULL));

	for(int i = 0; i < k; i++) {
		for(int j = 0; j < n; j++)
			Q[i * n + j] = (double)rand() / 16000.0;
		make_orthonormal_to_system(Q + i * n, Q, i, n);
	}
}

void generate_normalized_vector(double* v, int n) {
	generate_vector(v, n);
	normalize_vector(v, n);
}

void generate_vector(double* v, int n) {
	srand((unsigned)time(NULL));

	for(int i = 0; i < n; i++)
		v[i] = ((double)rand() - 16000.0) / 16000.0;
}

/*
Solves a general tridiagonal system of equations usinf LU
factorisation.

PARAMETERS
diag - array of the diagonal elements
ldiag - array of the lower diagonal elements
udiag - array of the upper diagonal elements
b - right-hand side of the equation
x - (output) solution vector
*/
int solve_tridiagonal_system(double* diag, double* ldiag, 
				   double* udiag, double* b, double* x, int n)
{
	static double* d = new double[n];
	static double* dl = new double[n - 1];
	static double* du = new double[n - 1];
	static double* du2 = new double[n - 2];
	static int* ipiv = new int[n];
	int i, info, c = 1;

	memcpy(d, diag, n * sizeof(double));
	memcpy(dl, ldiag, (n - 1) * sizeof(double));
	memcpy(du, udiag, (n - 1) * sizeof(double));

	for(i = 0; i < n; i++)
		ipiv[i] = i;

	dgttrf_(&n, dl, d, du, du2, ipiv, &info);

	if(info != 0)
	{
		printf("Error on LU factorisation.\n Error number is %d", info);

		delete [] d;
		delete [] dl;
		delete [] du;
		delete [] du2;
		delete [] ipiv;

		return info;
	}

	memcpy(x, b, n * sizeof(double));

	dgttrs_("N", &n, &c, dl, d, du, du2, ipiv, x, &n, &info);
	if(info != 0)
	{
		printf("Error on solving the tridiagonal system.\n Error number is %d", info);

		delete [] d;
		delete [] dl;
		delete [] du;
		delete [] du2;
		delete [] ipiv;

		return info;
	}

//	delete [] d;
//	delete [] dl;
//	delete [] du;
//	delete [] du2;
//	delete [] ipiv;

	return 0;
}
