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
   
#ifndef RKMATRIX
#define RKMATRIX
#include "basictools.h"
// #include <mkl.h>

#define RK_SVD_EPS 1.0e-16
#define RK_SVD_EPS2 1.0e-32
#define RK_SVD_EPS3 1.0e-64

typedef class rkmatrix* prkmatrix;

/*  An rkmatrix with paramters k,rows,cols is a factorisation
a * transpose(b) with (rows x k) matrix 'a' and (cols x k) matrix 'b'. */

class rkmatrix {
public:
	rkmatrix() { rows = 0; cols = 0; rk = 0; k = 0; a = 0x0; b = 0x0; } ;
	rkmatrix(int k_new, int rows_new, int cols_new) { k = k_new; rk = k_new; rows = rows_new; cols = cols_new; a = allocate_vector(rows * k); b = allocate_vector(cols * k); } ;
	~rkmatrix() { if(a) free_vector(a); if(b) free_vector(b); } ;

	int k; // allocated rank
	int rk; // real rank
	int rows;
	int cols;
	double* a;
	double* b;

	void rkmatrix_times_vector(EvalMode mode, double* b, double* c, double alpha, double beta);
	void rkmatrix_times_matrix(EvalMode mode, double* B, double* C, int n, double alpha, double beta);
	double get_frobenius_norm();
	int rk_svd(double* u, double* sigma, double* v);
};

#endif
