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
   
#ifndef HMATRIX
#define HMATRIX

#include "rkmatrix.h"
#include "fullmatrix.h"

typedef enum {
	hs_strong, 
	hs_weak
} HStrategy;

typedef struct _hmatrix hmatrix;
typedef hmatrix *phmatrix;

struct _hmatrix {
	int rows;
	int cols;
	int block_rows;
	int block_cols;
	int l;
	prkmatrix r;
	pfullmatrix f;
	phmatrix* sons;
};

phmatrix new_hmatrix(int r, int c);
void del_hmatrix(phmatrix h);

int create_hmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start);
int create_rkmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start);
int create_fullmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start);

int hmvm(phmatrix h, double* b, double* c);
int hmmm(phmatrix h, double* B, int B_cols, double* C);
int trihmvm(phmatrix h, char* uplo, double* b,double* c);
int trihmmm(phmatrix h, char* uplo, double* B,int B_cols, double* C);

#endif
