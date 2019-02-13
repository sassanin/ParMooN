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
   
#ifndef MPIMIS_BASICTOOLS
#define MPIMIS_BASICTOOLS

#include <math.h>
#ifdef __MAC64__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <assert.h>
#include "exported.h"

typedef bool boolean_function(int, int, void*, int n); // For ns1k_scf used only
typedef double function_1D(double);
typedef double function_2D(double, double);

typedef enum {
	EM_DEFAULT,
	EM_TRANSPOSED,
	EM_LOWER_TRIANGULAR,
	EM_UPPER_TRIANGULAR
} EvalMode;

typedef enum {
	DS_POINTWISE,
	DS_MEANVALUE,
	DS_LINEAR_FE,
	DS_L2_CONSTANT
} Discretization_Scheme;

#define ZERO_THRESHOLD	0.00000000000001	/* If |x| < ZERO_THRESHOLD then x is considered zero */
#define SQRT2			1.414213562373095049
#define TWOPI			6.283185307179586476
#define FOURPI			12.566370614359172953

/**************************************************************************
*********************** Functions working with vectors ********************
**************************************************************************/

/**************************************************************************
****************** Norms, scalar product, orthonormalisation **************
**************************************************************************/

double scalar_product(double* a, double* b, int n);
double euclidean_norm(double* a, int n);
double frobenius_norm(double* A, int rows, int cols); // TODO: Make inline (?)
double normalize_vector(double* a, int n);

// returns the space distance (angle) between to directions
double vector_distance(double* a, double* b, int n);

void make_orthonormal_to_system(double* a, double* A, int k, int n);

/**************************************************************************
*********************** Standard arithmetic operations ********************
**************************************************************************/

// Initialization
inline double* allocate_vector(int size) {
	return (double*)malloc(size * sizeof(double)); 
}
inline void free_vector(double* v) {
	if(v != NULL) free((void*)v);
}
inline void fill_vector(double* v, int n, double value) {
	for(int i = 0; i < n; i++) v[i] = value;
}
inline void fill_vector_rand(double* v, int n) {
	for(int i = 0; i < n; i++) v[i] = rand();
}
inline void clear_vector(double* v, int n) {
	fill_vector(v, n, 0);
}

// Standard arithmetics
inline void scale_vector(double* v, int n, double c) {
	for(int i = 0; i < n; i++) v[i] *= c; 
}
inline double sum_vector(double* v, int n);
inline double average_vector(double* v, int n);

// min-max functions
inline double max(double a, double b){return a > b ? a : b;};
inline double min(double a, double b){return a < b ? a : b;};
inline int max(int a, int b){return a > b ? a : b;};
inline int min(int a, int b){return a < b ? a : b;};

double min_vector(double* v, int n);
double max_vector(double* v, int n);
double min_abs_vector(double* v, int n);
double max_abs_vector(double* v, int n);

int min_vector(int* v, int n);
int max_vector(int* v, int n);
int min_abs_vector(int* v, int n);
int max_abs_vector(int* v, int n);

// Algebraic operations
// C += alpha * op(C) * op(B) + beta * C
void matrix_times_matrix(double* A, char* A_transposed, int A_rows, int A_cols, double* B, char* B_transposed, int B_rows, int B_cols, double* C, double alpha, double beta);

// c = alpha * op(A) * b + beta * c
void matrix_times_vector(double* A, char* A_transposed, int A_rows, int A_cols, double* b, double* c, double alpha, double beta);

// X = A * X
void triangular_matrix_times_matrix(double* A, char* uplo, double* X, int rows, int cols);

// x = A * x
void triangular_matrix_times_vector(double* A, char* uplo, int n, double* x);

// A is m x n fullmatrix, on output is itself's transposed
void transpose_matrix(double* A, int m, int n);

// B is the corresponding (full) block of matrix A
void get_matrix_block(double* A, int rows, int cols, int block_row_start, int block_col_start, int block_rows, int block_cols, double* B);

// A is n x n square fullmatrix, on output is itself's inverse
// return 0, if everything is ok, otherwise error number
int invert_matrix(double* A, int n);

// Printing
void print_vector(double* v, int n);
void print_matrix(double* A, int m, int n);
void fprint_vector(char* filename, double* v, int n);
void fscan_vector(char* filename, double* v, int n);
void fprint_matrix(char* filename, double* A, int rows, int cols);

// Fills the entries of the matrix via given 2D function and mesh
void fill_matrix(double* A, double (*)(int, int, void*), int rows, int cols, void* data);
void fill_matrix(double* A, function_2D* k, int rows, int cols, double* grdpoints);

/**************************************************************************
**************************** Fast calculations ****************************
**************************************************************************/

// return x^y
double ipow(double x, int y);

// return log_a(x) (a is the base)
double log(double x, double a);

// return x^y
int iipow(int x, int y);

double sgn(double x);
double sqr(double x);
double cube(double x);
double cubert(double x);
int factorial(int n);

// Fill data
void fill_hrefined_grid(double a, double b, int cLevels, int n, double* grid);
void fill_uniform_grid(double a, double b, int n, double* grid);

void fill_exponential_data(function_1D f, int cLevels, int n, double* d);
void fill_uniform_data(function_1D f, int n, double* d);

void transform_grid_m2l(double& xmin, double& xmax, int n, double* x);
void transform_grid_l2m(double& xmin, double& xmax, int n, double* x);

void fill_vector_from_function(function_1D f, double* grid, int n, double* v, Discretization_Scheme ds = DS_POINTWISE);

int get_hrefined_index(double h, int n, double x);
int get_uniform_index(double h, double x);
int get_index(double* x, int n, double x0);

/**************************************************************************
****************************** Testing tools ******************************
**************************************************************************/
void generate_orthonormal_matrix(double* Q, int k, int n);
void generate_vector(double* v, int n);
void generate_normalized_vector(double* v, int n);

/*
Solves a feneral tridiafonal system of equations usinf LU
factorisation.

PARAMETERS
diaf - array of the diafonal elements
ldiaf - array of the lower diafonal elements
udiaf - array of the upper diafonal elements
b - rifht-hand side of the equation
x - (output) solution vector
*/
int solve_tridiagonal_system(double* diaf, double* ldiaf, double* udiaf, double* b, double* x, int n);

inline int iround(double x) { return (int)(x + 0.5); }

#endif // MPIMIS_BASICTOOLS
