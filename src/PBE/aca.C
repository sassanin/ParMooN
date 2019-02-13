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
   
#include "aca.h"
#include "basic.h"
#include <malloc.h>
#include <assert.h>
#include <math.h>

static void
compute_row(const double *A, const double *B, int rows, int cols, 
	    int row_start, int col_start, 
	    int k, int i, const int *skip, double *row,
	    double (*entry)(int row, int col, void* data),
	    void *data)
{
  int j;
  
  assert(i >= 0);
  assert(i < rows);
  assert(entry != NULL);

  for(j = 0; j < cols; j++)
	  row[j] = entry(row_start + i, col_start + j, data);
//  fillrow(row_start+i, col_start, cols, skip, row, data);

  if(k > 0)
    dgemv("Not Transposed",
	   &cols, &k,
	   &meins,
	   B, &cols,
	   A+i, &rows,
	   &deins,
	   row, &eins);
  
  if(skip)
    for(j=0; j<cols; j++)
      if(skip[j])
	row[j] = 0.0;
}

static void
compute_column(const double *A, const double *B, int rows, int cols, 
	       int row_start, int col_start, 
	       int k, int j, const int *skip, double *col,
	       double (*entry)(int row, int col, void* data),
	       void *data)
{
  int i;
  
  assert(j >= 0);
  assert(j < cols);
  assert(entry != NULL);

  for(i = 0; i < rows; i++)
	  col[i] = entry(row_start + i, col_start + j, data);
//  fillcol(col_start+j, row_start, rows, skip, col, data);

  if(k > 0)
    dgemv("Not Transposed",
	   &rows, &k,
	   &meins,
	   A, &rows,
	   B+j, &cols,
	   &deins,
	   col, &eins);

  if(skip)
    for(i=0; i<rows; i++)
      if(skip[i])
	col[i] = 0.0;
}

int
aca_fill_block(double *A, double *B, int rows, int cols,
	       int row_start, int col_start,
	       double (*entry)(int row, int col, void *data),
	       void *data, int kmax, double eps, ACAStrategy strategy)
{
  double *col;
  double *row;
  double maxval, invmax, blknorm, apxnorm;
  int i, itmp, j, jtmp, jold;
  int tries, stop_aca;
  int k;

  k = 0;
  tries = cols * rows + 1;
  j = 0;
  apxnorm = 0.0;
  stop_aca = 0;
  while(k < kmax && tries > 0) {
    tries--;

    col = A + k*rows;
    row = B + k*cols;

    compute_column(A, B, rows, cols, row_start, col_start, k, j, NULL, col, entry, data);

    maxval = fabs(col[0]);
    i = 0;
    for(itmp = 1; itmp < rows; itmp++)
      if(fabs(col[itmp]) > maxval) {
	maxval = fabs(col[itmp]);
	i = itmp;
      }

    if(maxval < ACA_EPS) {
/*      if(k==0) {
	compute_row(A, B, rows, cols, row_start, col_start, k, 0, NULL, row, entry, data);
	maxval = fabs(row[0]);
	j = 0;
	for(jtmp = 1; jtmp < cols; jtmp++)
	  if(fabs(row[jtmp]) > maxval) {
	    maxval = fabs(row[jtmp]);
	    j = jtmp;
	  }
      } else { */
	j = (j+1) % cols;
//      }
      break;
    }

    compute_row(A, B, rows, cols, row_start, col_start, k, i, NULL, row, entry, data);

    invmax = 1.0 / col[i];
    dscal(&rows, &invmax, col, &eins);

    maxval = fabs(row[0]);
    jold = j;
    j = 0; 
    if(jold==0) j=1;
    for(jtmp = 1; jtmp < cols; jtmp++)
      if(fabs(row[jtmp]) > maxval && jtmp!=jold) {
	maxval = fabs(row[jtmp]);
	j = jtmp;
      }
    
    switch(strategy) {
      default:
      case HLIB_ACA_DEFAULT:
      case HLIB_ACA_RELATIVE2:
	blknorm = (dnrm2(&rows, col, &eins) *
		   dnrm2(&cols, row, &eins));
	if(k == 0) {
	  apxnorm = blknorm;
	  stop_aca = 0;
	}
	else
	  stop_aca = (blknorm < apxnorm * eps);
	break;

      case HLIB_ACA_ABSOLUTE2:
	blknorm = (dnrm2(&rows, col, &eins) *
		   dnrm2(&cols, row, &eins));
	if(k == 0) {
	  apxnorm = blknorm;
	  stop_aca = 0;
	}
	else
	  stop_aca = (blknorm < eps);
	break;

      case HLIB_ACA_ABSOLUTE:
	blknorm = (dnrm2(&rows, col, &eins) *
		   dnrm2(&cols, row, &eins));
	stop_aca = (blknorm < eps);
	break;

      case HLIB_ACA_RELATIVE:
	for(jtmp = 0; jtmp < k; jtmp++)
	  apxnorm += (ddot(&rows, col, &eins, A+jtmp*rows, &eins) *
		      ddot(&cols, row, &eins, B+jtmp*cols, &eins));
	blknorm = (ddot(&rows, col, &eins, col, &eins) *
		   ddot(&cols, row, &eins, row, &eins));
	stop_aca = (blknorm < apxnorm * eps * eps);
	if(!stop_aca)
	  apxnorm += blknorm;
	break;
    }
    
    if(stop_aca)
      break;

    k++;
  }

  return k;
}






static void
newcompute_row(const double *A, const double *B, int rows, int cols, int k,
	       int i, double *row,
	       int row_off, int col_off,
	       double (*entry)(int row, int col, void *data), void *data,
	       int *col_done)
{
  int j;

  assert(i >= 0);
  assert(i < rows);

  for(j=0; j<cols; j++)
    if(col_done == NULL || col_done[j] == 0)
      row[j] = entry(row_off+i, col_off+j, data);

  if(k > 0)
    dgemv("Not Transposed",
	   &cols, &k,
	   &meins,
	   B, &cols,
	   A+i, &rows,
	   &deins,
	   row, &eins);
 
  for(j=0; j<cols; j++)
    if(col_done != NULL && col_done[j] != 0)
      row[j] = 0.0;
 
}

static void
newcompute_column(const double *A, const double *B, int rows, int cols, int k,
		  int j, double *col,
		  int row_off, int col_off,
		  double (*entry)(int row, int col, void *data), void *data,
		  int *row_done)
{
  int i;

  assert(j >= 0);
  assert(j < cols);

  for(i=0; i<rows; i++)
    if(row_done == NULL || row_done[i] == 0)
      col[i] = entry(row_off+i, col_off+j, data);

  if(k > 0)
    dgemv("Not Transposed",
	   &rows, &k,
	   &meins,
	   A, &rows,
	   B+j, &cols,
	   &deins,
	   col, &eins);

  for(i=0; i<rows; i++)
    if(row_done != NULL && row_done[i] != 0)
      col[i] = 0.0;
}

int
newaca_fill_block(double *A, double *B, int rows, int cols,
		  int row_off, int col_off,
		  double (*entry)(int row, int col, void *data), void *data,
		  int kmax, double eps, ACAStrategy strategy){
  
  double *col, *col_ref=0x0;
  double *row, *row_ref=0x0;
  double row_maxval, col_maxval, invmax;
  double blknorm, apxnorm, colnorm,rownorm;
  int i, itmp, j, jtmp, i_ref=0, j_ref=0;
  int tries, tries_col, tries_row, stop_aca;
  int k;
  int *row_done, *col_done;
  int *row_illicit, *col_illicit;
  int col_zeros=0;

  if(kmax<1 || rows<1 || cols<1) return 0;     /* rank 0 approximation */
  tries = kmax+2;
  tries_row = 2;
  tries_col = 2;

  row_done = (int *) malloc((size_t) sizeof(int) * rows);
  assert(row_done != NULL);
  col_done = (int *) malloc((size_t) sizeof(int) * cols);
  assert(col_done != NULL);
  row_illicit = (int *) malloc((size_t) sizeof(int) * rows);
  assert(row_illicit != NULL);
  col_illicit = (int *) malloc((size_t) sizeof(int) * cols);
  assert(col_illicit != NULL);

  k = 0;
  for(i=0; i<rows; i++){
    row_done[i] = 0;
    row_illicit[i] = 0;
  }
  for(j=0; j<cols; j++){
    col_done[j] = 0;
    col_illicit[j] = 0;
  }


  j_ref = 0;        /* reference column and row */
  col_ref = A + (kmax-1)*rows;
  newcompute_column(A, B, rows, cols, k, j_ref, col_ref, row_off, col_off,
		    entry, data, NULL);
  colnorm = dnrm2(&rows, col_ref, &eins);
  /*
  if(colnorm<ACA_EPS){
    for(j=j_ref; j!=(j_ref+cols-1)%cols && colnorm<ACA_EPS &&tries_col>0; 
	j=(j+1)%cols,tries_col--){
      newcompute_column(A, B, rows, cols, k, j, col_ref, row_off, col_off,
			entry, data, NULL);
      colnorm = dnrm2(&rows, col_ref, &eins);
      if(colnorm<ACA_EPS) col_done[j] = 1;
    }
    j_ref = (j+cols-1)%cols;
  }
  */
  /* special case: some entries of j_ref are very small */
  rownorm = 0.0;
  if(colnorm > ACA_EPS){
    i_ref = 0;
    rownorm = fabs(col_ref[0]);
    for(i=1; i<rows; i++){
      if(fabs(col_ref[i])<rownorm){
	rownorm = fabs(col_ref[i]);
	i_ref = i;
      }
    }
    row_ref = B + (kmax-1)*cols;
    newcompute_row(A, B, rows, cols, k, i_ref, row_ref, row_off, col_off,
		   entry, data, NULL);
    rownorm = dnrm2(&cols, row_ref, &eins);
    if(fabs(col_ref[i_ref]) < ACA_EPS) col_zeros = 1;
  }else{
    row_ref = B + (kmax-1)*cols;
    i_ref = 0; 
    newcompute_row(A, B, rows, cols, k, i_ref, row_ref, row_off, col_off,
		   entry, data, NULL);
    rownorm = dnrm2(&cols, row_ref, &eins);
    
    if(rownorm<ACA_EPS){
      for(i=i_ref; i!=(i_ref+rows-1)%rows && rownorm<ACA_EPS && tries_row>0; 
	  i=(i+1)%rows,tries_row--){
	newcompute_row(A, B, rows, cols, k, i, row_ref, row_off, col_off,
		       entry, data, NULL);
	rownorm = dnrm2(&cols, row_ref, &eins);
	if(rownorm<ACA_EPS) row_done[i] = 1;
      }
      i_ref = (i+rows-1)%rows;
    }
  }
  
  /*
  for(j=0; j<cols; j++){
    if(fabs(row_ref[j]) < ACA_EPS) row_zeros = 1;
  }
  */

  if(rownorm<ACA_EPS && colnorm<ACA_EPS){
    free(col_done);
    free(row_done);
    free(col_illicit);
    free(row_illicit);
    return 0;
  }
  
  if(col_zeros!=0){
    if(rownorm > ACA_EPS){
      row = B;
      for(i=0; i<rows && fabs(col_ref[i])<ACA_EPS; i++);
      if(i<rows){
	newcompute_row(A, B, rows, cols, 0, i, row, row_off, col_off, entry, data,col_done);
	i = 0;
	for(j=0; j<cols && i==0; j++){
	  if(fabs(row_ref[j])>ACA_EPS && fabs(row[j])>ACA_EPS) i=1;
	}
	if(i==0){
	  for(i=0; i<rows; i++)
	    if(row_done[i]==1 || fabs(col_ref[i])>ACA_EPS)
	      row_illicit[i] = 1;
	  for(j=0; j<cols; j++)
	    if(col_done[j]==1 || fabs(row_ref[j])>ACA_EPS)
	      col_illicit[j] = 1;
	}
      }
    }
  }
  apxnorm = 0.0;
  stop_aca = 0;
  
  while(k < kmax-1 && (tries_row > 0 || tries_col>0) && tries>0) {
    tries--;

    col = A + k*rows;
    row = B + k*cols;
    
    col_maxval = 0.0;
    if(colnorm>ACA_EPS){
      i = 0; col_maxval = fabs(col_ref[0]);
      for(itmp = 1; itmp < rows; itmp++){
	if(fabs(col_ref[itmp]) > col_maxval) {
	  col_maxval = fabs(col_ref[itmp]);
	  i = itmp;
	}
      }
    }
    row_maxval = 0.0;
    if(rownorm>ACA_EPS){
      j = 0; row_maxval = fabs(row_ref[0]);
      for(jtmp = 1; jtmp < cols; jtmp++){
	if(fabs(row_ref[jtmp]) > row_maxval) {
	  row_maxval = fabs(row_ref[jtmp]);
	  j = jtmp;
	}
      }
    }
    if(row_maxval>col_maxval){
      if(j!=j_ref){
	newcompute_column(A, B, rows, cols, k, j, col, row_off, col_off,
			  entry, data, row_done);
      }else{
	dcopy(&rows,col_ref,&eins,col,&eins);
      }
      i = 0; col_maxval = fabs(col[0]);
      for(itmp = 1; itmp < rows; itmp++){
	if(fabs(col[itmp]) > col_maxval) {
	  col_maxval = fabs(col[itmp]);
	  i = itmp;
	}
      }
	
      if(col_maxval < ACA_EPS) {
	stop_aca = 1;
      }else{
	newcompute_row(A, B, rows, cols, k, i, row, row_off, col_off,
		       entry, data, col_done);
	invmax = 1.0 / col[i];
	dscal(&rows, &invmax, col, &eins);
	
      }
    }else{
      if(i!=i_ref){
	newcompute_row(A, B, rows, cols, k, i, row, row_off, col_off,
		       entry, data, col_done);
      }else{
	dcopy(&cols,row_ref,&eins,row,&eins);
      }
      j = 0; row_maxval = fabs(row[0]);
      for(jtmp = 1; jtmp < cols; jtmp++){
	if(fabs(row[jtmp]) > row_maxval) {
	  row_maxval = fabs(row[jtmp]);
	  j = jtmp;
	}
      }
      
      if(row_maxval < ACA_EPS) {
	stop_aca = 1;
      }else{
	newcompute_column(A, B, rows, cols, k, j, col, row_off, col_off,
			  entry, data, row_done);
	invmax = 1.0 / row[j];
	dscal(&cols, &invmax, row, &eins);
      }
    }
    
    
    row_done[i] = 1;
    col_done[j] = 1;
    row_illicit[i] = 1;
    col_illicit[j] = 1;
    
    if(i!=i_ref){
      invmax = -col[i_ref];
      daxpy(&cols,&invmax,row,&eins,row_ref,&eins);
      rownorm = dnrm2(&cols, row_ref, &eins);
    }
    if(i==i_ref || rownorm<ACA_EPS){
      if(i==i_ref) tries_row++;
      if(tries_row>0){
	/* new reference row */
	rownorm = 0.0;
	for(i=i_ref; i!=(i_ref+rows-1)%rows && rownorm<ACA_EPS && tries_row>0; 
	    i=(i+1)%rows){
	  if(row_illicit[i]==0){
	    newcompute_row(A, B, rows, cols, k+1, i, row_ref, row_off, col_off,
			   entry, data, col_done);
	    rownorm = dnrm2(&cols, row_ref, &eins);
	    if(rownorm<ACA_EPS){
	      row_done[i] = 1;
	      row_illicit[i] = 1;
	    }
	    tries_row--;
	  }else{
	    rownorm = 0.0;
	  }
	}
	i_ref = (i+rows-1)%rows;
      }
    }
    if(j!=j_ref){
      /* update reference column */
      invmax = -row[j_ref];
      daxpy(&rows,&invmax,col,&eins,col_ref,&eins);
      colnorm = dnrm2(&rows, col_ref, &eins);
    }
    if(j==j_ref || colnorm<ACA_EPS){
      if(j==j_ref) tries_col++;
      if(tries_col>0){
	/* new reference column */
	colnorm = 0.0;
	for(j=j_ref; j!=(j_ref+cols-1)%cols && colnorm<ACA_EPS && tries_col>0; 
	    j=(j+1)%cols){
	  if(col_illicit[j]==0){
	    newcompute_column(A, B, rows, cols, k+1, j, col_ref, row_off, col_off,
			      entry, data, row_done);
	    colnorm = dnrm2(&rows, col_ref, &eins);
	    if(colnorm<ACA_EPS){
	      col_done[j] = 1;
	      col_illicit[j] = 1;
	    }
	    tries_col--;
	  }else{
	    colnorm = 0.0;
	  }
	}
	j_ref = (j+cols-1)%cols;
      }
    }
    if(colnorm<ACA_EPS && rownorm<ACA_EPS){
      stop_aca = 1;
      k++;
    }

    if(stop_aca==0){
      switch(strategy) {
      default:
      case HLIB_ACA_DEFAULT:
	
      case HLIB_ACA_RELATIVE2:
	blknorm = (dnrm2(&rows, col, &eins) *
		   dnrm2(&cols, row, &eins));
	if(k == 0) {
	  apxnorm = blknorm;
	  stop_aca = 0;
	}
	else{
	  if(blknorm < apxnorm * eps
	     && rownorm < apxnorm * eps 
	     && colnorm < apxnorm * eps) stop_aca = 1;  
	}
	break;
	
      case HLIB_ACA_ABSOLUTE2:
	blknorm = (dnrm2(&rows, col, &eins) *
		   dnrm2(&cols, row, &eins));
	if(k == 0) {
	  apxnorm = blknorm;
	  stop_aca = 0;
	}
	else
	  stop_aca = (blknorm < eps);
	break;
	
      case HLIB_ACA_ABSOLUTE:
	blknorm = (dnrm2(&rows, col, &eins) *
		   dnrm2(&cols, row, &eins));
	stop_aca = (blknorm < eps);
	break;
	
      case HLIB_ACA_RELATIVE:
	for(jtmp = 0; jtmp < k; jtmp++)
	  apxnorm += (ddot(&rows, col, &eins, A+jtmp*rows, &eins) *
		      ddot(&cols, row, &eins, B+jtmp*cols, &eins));
	blknorm = (ddot(&rows, col, &eins, col, &eins) *
		   ddot(&cols, row, &eins, row, &eins));
	stop_aca = (blknorm < apxnorm * eps * eps);
	if(!stop_aca)
	  apxnorm += blknorm;
	break;
      } 
    }
    if(stop_aca) break;
    k++;
  }
  
  free(col_done);
  free(row_done);
  free(col_illicit);
  free(row_illicit);
  
  return k;
}

