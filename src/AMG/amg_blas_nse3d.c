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
   
/****************************************************************************/
/*                                                                          */
/* File:          amg_blas_nse3d.c                                          */
/*                                                                          */
/* Purpose:   matrix vector operations for 3d Navier-Stokes problems        */
/*                                                                          */
/* History:   03 AUG 2000 file created                                      */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/* Author:    Volker John                                                   */
/*            Otto--von--Guericke Universitaet Magdeburg                    */
/*            Institut fuer Analysis und Numerik                            */
/*            Postfach 4120                                                 */
/*            39016 Magdeburg                                               */
/*            email : volker.john@mathematik.uni-magdeburg.de               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*    system include files                                                  */
/*    application include files                                             */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <amg_header.h>
#include <amg_low.h>
#include <amg_sp.h>
#include <amg_blas.h>
#include <amg_solve_main.h>
#include <amg_1d_prec.h>

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

extern AMG_VECTOR *velo_prolong[AMG_MAX_LEVELS],*pres_prolong[AMG_MAX_LEVELS];
extern AMG_VECTOR *velo_result;
extern AMG_GRAPH *G_schur[AMG_MAX_LEVELS],*G[AMG_MAX_LEVELS];
extern AMG_SolverContext *global_sc;
extern AMG_MATRIX *M[AMG_MAX_LEVELS];

/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR 2D SADDLE POINT PROBLEMS OF TYPE 1          */
/*                                                                          */
/* ( A    0     B_0)                                                        */
/* ( 0    A     B_1) x = y                                                  */
/* (B_0^T B_1^T 0  )                                                        */
/*                                                                          */
/* Input : A , B_0^T, B_1^T, x, y                                           */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
int AMG_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l,ll;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_1 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
        printf("A(%d,%d) = %g; \n",i+1,i+1,a[start]);
        printf("A(%d,%d) = %g; \n",n_a+i+1,n_a+i+1,a[start]);
      for (k=start+1; k<end; k++) 
{
        printf("A(%d,%d) = %g; \n",i+1,ja[k]+1,a[k]);
        printf("A(%d,%d) = %g; \n",n_a+i+1,n_a+ja[k]+1,a[k]);
}
    }
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);

 for (i=0; i<m_b0; i++)
    {
      start = ra[i];end = ra[i+1]-1;
      for (k=start; k<=end; k++) 
        {        
          printf("A(%d,%d) = %g; \n",2*n_a+i+1,ja[k]+1,a[k]);
          printf("A(%d,%d) = %g; \n",ja[k]+1,2*n_a+i+1,a[k]);
        }
    }
  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

 for (i=0; i<m_b1; i++)
    {
      start = ra[i];end = ra[i+1]-1;
      for (k=start; k<=end; k++) {
        printf("A(%d,%d) = %g; \n",i+1+2*n_a,ja[k]+n_a+1,a[k]);
        printf("A(%d,%d) = %g; \n",ja[k]+n_a+1,i+1+2*n_a,a[k]);}
    }
 exit(1);
 for (i=0;i<n_y;i++){

    y[i] = 1.0*i;
    printf("y(%d) = %g\n",i+1,y[i]);}
*/

  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[i] = s0;
      x[i+n_a] = s1;
    }

 /* multiplication with block B_0 and B_0^T */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
            x[l] += a[k]*y[2*n_a+i];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] = s0;     
    }
 /* multiplication with block B_1 and B_1^T*/

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
            x[l+n_a] += a[k]*y[2*n_a+i];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] += s0;
    }

  return(AMG_OK);
}

int AMG_dmatminus_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l,ll;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }

  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          l =ja[k]; 
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[i] -= s0;
      x[i+n_a] -= s1;
    }

 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
            x[l] -= a[k]*y[2*n_a+i];
          s0 += a[k]*y[l];
        } 
      x[2*n_a+i] -= s0;
    }

 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
            x[l+n_a] -= a[k]*y[2*n_a+i];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] -= s0;
    }
  return(AMG_OK);
}

int AMG_A_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s0,s1;
  register int b,bb;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[i] = s0;
      x[i+n_a] = s1;
    }

  return(AMG_OK);
}

int AMG_A_dmatminus_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s0,s1;
  register int b,bb;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[i] -= s0;
      x[i+n_a] -= s1;
    }
  return(AMG_OK);
}

int AMG_B_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                   AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l,ll;
  register double *x, *y, *a;
  register int *ra, *ja, *aja , *ara;
  register int n_b2,m_b2;
  
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if ((m_b0!=m_b1)||(m_b0!=m_b2))
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_3_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=n_b0+n_b1+n_b2)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_3_TYPE_1 - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_y!=m_b0)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_3_TYPE_1 - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  for (i=0; i<n_x; i++) x[i] = 0.0; 

 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll = ara[l];
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l] += a[k]*y[i];
    }
  }

 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll = ara[l];
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l+n_b0] += a[k]*y[i];
    }
  }

 /* multiplication with block B_2 */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  n_b2 = n_b0+n_b1;
 
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll = ara[l];
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l+n_b2] += a[k]*y[i];
    }
  }

  return(AMG_OK);
}

int AMG_B_Trans_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                         AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,n_b2,m_b2,l;
  register double *x, *y, *a;
  register int *ra, *ja;
  register double s0;
  
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if ((m_b0!=m_b1)||(m_b0!=m_b2))
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_3_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_y!=n_b0+n_b1+n_b2)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_3_TYPE_1 - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_x!=m_b0)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_3_TYPE_1 - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);

 /* multiplication with block B_0^T */
  
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
 
  for (i=0; i<m_b0; i++)
  {
    start = ra[i]; 
    end= ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      s0 += a[k]*y[l];        
    }
    x[i] = s0;
    }
  
 /* multiplication with block B_1^T */
  
  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);
 
  for (i=0; i<m_b1; i++)
  {
    start = ra[i]; 
    end= ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++) 
    {
      l = ja[k]+n_b0;
      s0 += a[k]*y[l];
    }        
    x[i] += s0;
  }

 /* multiplication with block B_2^T */
  
  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  n_b2 =n_b0 + n_b1;
  for (i=0; i<m_b2; i++)
  {
    start = ra[i]; 
    end= ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++) 
    {
      l = ja[k]+n_b2;
      s0 += a[k]*y[l];
    }        
    x[i] += s0;
  }
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR 2D SADDLE POINT PROBLEMS OF TYPE 2          */
/*                                                                          */
/* ( A    0     B_0)                                                        */
/* ( 0    A     B_1) x = y                                                  */
/* (B_2  B_3     0  )                                                       */
/*                                                                          */
/* Input : A , B_0, B_1, B_2, B_3, x, y                                     */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
int AMG_dmatmul_SADDLE_3_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_2 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_2 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[i] = s0;
      x[i+n_a] = s1;
    }

 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[2*n_a+l];
        }
      x[i] += s0;
   }
 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[2*n_a+l];
        }
      x[i+n_a] += s0;
   }

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] = s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] += s0;
    }
  return(AMG_OK);
}

int AMG_dmatminus_SADDLE_3_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_2 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_2 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[i] -= s0;
      x[i+n_a] -= s1;
    }

 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[2*n_a+l];
        }
      x[i] -= s0;
   }
 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[2*n_a+l];
        }
      x[i+n_a] -= s0;
   }

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] -= s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] -= s0;
    }
  return(AMG_OK);
}

int AMG_B_dmatmul_SADDLE_3_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l,ll;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *aja , *ara;
  register double s0,s1;
  register int b,bb;
  
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_b0!=n_b1)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_3_TYPE_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=m_b0+m_b1)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_3_TYPE_2 - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_y!=n_b0)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_3_TYPE_2 - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  for (i=0; i<n_x; i++) x[i] = 0.0; 

 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[i] = s0;
   }
 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[l];
        }
      x[i+m_b0] = s0;
   }

  return(AMG_OK);
}

int AMG_B_Trans_dmatmul_SADDLE_3_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s0,s1;
  register int b,bb;
  
  n_b0 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_3_TYPE_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_y!=n_b0+n_b1)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_3_TYPE_2 - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_x!=m_b0)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_3_TYPE_2 - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[i] = s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k]+n_b0;
          s0 += a[k]*y[l];
        }
      x[i] += s0;
    }

  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR 2D SADDLE POINT PROBLEMS OF TYPE 3          */
/*                                                                          */
/* ( A     B_2     B_0  )                                                   */
/* (B_3    B_4     B_1  ) x = y                                             */
/* (B_0^T  B_1^T     0  )                                                   */
/*                                                                          */
/* Input : A , B_0^T , B_1^T , B_2, B_3, B_4 x, y                           */
/*                                                                          */
/****************************************************************************/

int AMG_dmatmul_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a11,*a12,*a21,*a22;
  register double s0,s1;
  register int b,bb;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_3 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_3 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_3 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a11 = AMG_MATRIX_A(A);
  a12 = AMG_MATRIX_A(B[2]);
  a21 = AMG_MATRIX_A(B[3]);
  a22 = AMG_MATRIX_A(B[4]);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  for (i=0; i<n_a; i++)
  {
    start = ra[i]; 
    end = start+ja[start];
    s0 = a11[start]*y[i];
    s1 = a22[start]*y[i+n_a];
    for (k=start+1; k<end; k++) 
    {
      l = ja[k];
      s0 += a11[k]*y[l];            /* A11 */
      if (i<B[2]->active) /* no active node */
         s0 += a12[k]*y[l+n_a];      /* A12 */
      s1 += a22[k]*y[l+n_a];         /* A22 */ 
      if (i<B[3]->active) /* no active node */
        s1 += a21[k]*y[l];           /* A12 */
    }
    x[i] = s0;
    x[i+n_a] = s1;    
  }
   
 /* multiplication with block B_0 and B_0^T */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
            x[l] += a[k]*y[2*n_a+i];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] = s0;     
    }

  /* multiplication with block B_1 and B_1^T*/

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
            x[l+n_a] += a[k]*y[2*n_a+i];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] += s0;
    }

  return(AMG_OK);
}

int AMG_dmatminus_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register double *a11,*a12,*a21,*a22;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatminus SADDLE_3_TYPE_1 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a11 = AMG_MATRIX_A(A);
  a12 = AMG_MATRIX_A(B[2]);
  a21 = AMG_MATRIX_A(B[3]);
  a22 = AMG_MATRIX_A(B[4]);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  for (i=0; i<n_a; i++)
  {
    start = ra[i]; 
    end = start+ja[start];
    s0 = a11[start]*y[i];
    s1 = a22[start]*y[i+n_a];
    for (k=start+1; k<end; k++) 
    {
      l = ja[k];
      s0 += a11[k]*y[l];            /* A11 */
      if (i<B[2]->active) /* no active node */
         s0 += a12[k]*y[l+n_a];      /* A12 */
      s1 += a22[k]*y[l+n_a];         /* A22 */ 
      if (i<B[3]->active) /* no active node */
        s1 += a21[k]*y[l];           /* A12 */
    }
    x[i] -= s0;
    x[i+n_a] -= s1;    
  }
 
 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
            x[l] -= a[k]*y[2*n_a+i];
          s0 += a[k]*y[l];
        } 
      x[2*n_a+i] -= s0;
    }

 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
            x[l+n_a] -= a[k]*y[2*n_a+i];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] -= s0;
    }

  return(AMG_OK);
}

int AMG_A_dmatmul_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                   AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a11,*a12,*a21,*a22;
  register double s0,s1;
  register int b,bb;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  level = A->level;
  if (level==0) /* finest level */
  {
    a11 = AMG_MATRIX_A(A);
    a12 = AMG_MATRIX_A(B[2]);
    a21 = AMG_MATRIX_A(B[3]);
    a22 = AMG_MATRIX_A(B[4]);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a22[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];            /* A11 */
        if (i<B[2]->active) /* no active node */
          s0 += a12[k]*y[l+n_a];      /* A12 */
        s1 += a22[k]*y[l+n_a];         /* A22 */ 
        if (i<B[3]->active) /* no active node */
          s1 += a21[k]*y[l];           /* A12 */
      }
      x[i] = s0;
      x[i+n_a] = s1;    
    }
  }
  else /* coarser levels, there is only A */
  {
    a11 = AMG_MATRIX_A(A);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a11[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];
        s1 += a11[k]*y[l+n_a];
      }
      x[i] = s0;
      x[i+n_a] = s1;
    }
  }

  return(AMG_OK);
}
int AMG_A_dmatminus_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                     AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a11,*a12,*a21,*a22;
  register double s0,s1;
  register int b,bb;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  level = A->level;
  
  if (level==0) /* finest level */
  {
    a11 = AMG_MATRIX_A(A);
    a12 = AMG_MATRIX_A(B[2]);
    a21 = AMG_MATRIX_A(B[3]);
    a22 = AMG_MATRIX_A(B[4]);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a22[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];            /* A11 */
        if (i<B[2]->active) /* no active node */
          s0 += a12[k]*y[l+n_a];      /* A12 */
        s1 += a22[k]*y[l+n_a];         /* A22 */ 
        if (i<B[3]->active) /* no active node */
          s1 += a21[k]*y[l];           /* A12 */
      }
      x[i] -= s0;
      x[i+n_a] -= s1;    
    }
  }
  else /* coarser levels, there is only A */
  {
    a11 = AMG_MATRIX_A(A);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a11[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];
        s1 += a11[k]*y[l+n_a];
      }
      x[i] -= s0;
      x[i+n_a] -= s1;
    }
  }
  
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR 2D SADDLE POINT PROBLEMS OF TYPE 4          */
/*                                                                          */
/* ( A   B_4     B_0)                                                       */
/* (B_5  B_6     B_1) x = y                                                 */
/* (B_2  B_3     0  )                                                       */
/*                                                                          */
/* Input : A , B_0, B_1, B_2, B_3, B_4, B_5, B_6, x, y                      */
/*                                                                          */
/****************************************************************************/

int AMG_dmatmul_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a11,*a12,*a21,*a22;
  register double s0,s1;
  register int b,bb;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a11 = AMG_MATRIX_A(A);
  a12 = AMG_MATRIX_A(B[4]);
  a21 = AMG_MATRIX_A(B[5]);
  a22 = AMG_MATRIX_A(B[6]);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  for (i=0; i<n_a; i++)
  {
    start = ra[i]; 
    end = start+ja[start];
    s0 = a11[start]*y[i];
    s1 = a22[start]*y[i+n_a];
    for (k=start+1; k<end; k++) 
    {
      l = ja[k];
      s0 += a11[k]*y[l];            /* A11 */
      if (i<B[4]->active) /* no active node */
         s0 += a12[k]*y[l+n_a];      /* A12 */
      s1 += a22[k]*y[l+n_a];         /* A22 */ 
      if (i<B[5]->active) /* no active node */
        s1 += a21[k]*y[l];           /* A12 */
    }
    x[i] = s0;
    x[i+n_a] = s1;    
  }
   
  /* multiplication with block B_0 */
  
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);
  
  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[2*n_a+l];
        }
      x[i] += s0;
   }
 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[2*n_a+l];
        }
      x[i+n_a] += s0;
   }

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] = s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] += s0;
    }
  return(AMG_OK);
}

int AMG_dmatminus_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register double *a11,*a12,*a21,*a22;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatmul SADDLE_3_TYPE_4 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a11 = AMG_MATRIX_A(A);
  a12 = AMG_MATRIX_A(B[4]);
  a21 = AMG_MATRIX_A(B[5]);
  a22 = AMG_MATRIX_A(B[6]);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  for (i=0; i<n_a; i++)
  {
    start = ra[i]; 
    end = start+ja[start];
    s0 = a11[start]*y[i];
    s1 = a22[start]*y[i+n_a];
    for (k=start+1; k<end; k++) 
    {
      l = ja[k];
      s0 += a11[k]*y[l];            /* A11 */
      if (i<B[4]->active) /* no active node */
         s0 += a12[k]*y[l+n_a];      /* A12 */
      s1 += a22[k]*y[l+n_a];         /* A22 */ 
      if (i<B[5]->active) /* no active node */
        s1 += a21[k]*y[l];           /* A12 */
    }
    x[i] -= s0;
    x[i+n_a] -= s1;    
  }
 
 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[2*n_a+l];
        }
      x[i] -= s0;
   }
 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[2*n_a+l];
        }
      x[i+n_a] -= s0;
   }

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] -= s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] -= s0;
    }
  return(AMG_OK);
}

int AMG_A_dmatmul_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                   AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a11,*a12,*a21,*a22;
  register double s0,s1;
  register int b,bb;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  level = A->level;
  if (level==0) /* finest level */
  {
    a11 = AMG_MATRIX_A(A);
    a12 = AMG_MATRIX_A(B[4]);
    a21 = AMG_MATRIX_A(B[5]);
    a22 = AMG_MATRIX_A(B[6]);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a22[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];            /* A11 */
        if (i<B[4]->active) /* no active node */
          s0 += a12[k]*y[l+n_a];      /* A12 */
        s1 += a22[k]*y[l+n_a];         /* A22 */ 
        if (i<B[5]->active) /* no active node */
          s1 += a21[k]*y[l];           /* A12 */
      }
      x[i] = s0;
      x[i+n_a] = s1;    
    }
  }
  else /* coarser levels, there is only A */
  {
    a11 = AMG_MATRIX_A(A);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a11[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];
        s1 += a11[k]*y[l+n_a];
      }
      x[i] = s0;
      x[i+n_a] = s1;
    }
  }

  return(AMG_OK);
}
int AMG_A_dmatminus_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                   AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a11,*a12,*a21,*a22;
  register double s0,s1;
  register int b,bb;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  level = A->level;
  
  if (level==0) /* finest level */
  {
    a11 = AMG_MATRIX_A(A);
    a12 = AMG_MATRIX_A(B[4]);
    a21 = AMG_MATRIX_A(B[5]);
    a22 = AMG_MATRIX_A(B[6]);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a22[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];            /* A11 */
        if (i<B[4]->active) /* no active node */
          s0 += a12[k]*y[l+n_a];      /* A12 */
        s1 += a22[k]*y[l+n_a];         /* A22 */ 
        if (i<B[5]->active) /* no active node */
          s1 += a21[k]*y[l];           /* A12 */
      }
      x[i] -= s0;
      x[i+n_a] -= s1;    
    }
  }
  else /* coarser levels, there is only A */
  {
    a11 = AMG_MATRIX_A(A);
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a11[start]*y[i];
      s1 = a11[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a11[k]*y[l];
        s1 += a11[k]*y[l+n_a];
      }
      x[i] -= s0;
      x[i+n_a] -= s1;
    }
  }
  
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR BRAESS-SARAZIN SMOOTHERS FOR 3D SADDLE      */
/* POINT PROBLEMS OF TYPE 1                                                 */
/*                                                                          */
/* ( f(A)    0       0    B_1)                                              */
/* ( 0      f(A)     0    B_2) x = y                                        */
/* ( 0       0      f(A)  B_3)                                              */
/* (B_1^T   B_2^T  B_3^T   0 )                                              */
/*                                                                          */
/* Input : A , B_1^T, B_2^T, B_3^T, x, y                                    */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,n_b2,m_b2,l,ll;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1,s2;
  register int b,bb, n_a2, n_a3;
  double alpha;
  int level;

  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b0!=m_b1)||((m_b0!=m_b2)))
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=3*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (3*n_a!=n_b0+n_b1+n_b2)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_1 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  n_a2 = 2*n_a;
  n_a3 = 3*n_a;
  /* multiplication with block A */
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    for (i=0; i<n_a; i++)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = a[start]*y[i];
        s1 = a[start]*y[i+n_a];
        s2 = a[start]*y[i+n_a2];
        for (k=start+1; k<end; k++) 
          {
            l = ja[k];
            s0 += a[k]*y[l];
            s1 += a[k]*y[l+n_a];
            s2 += a[k]*y[l+n_a2];
          }
        x[i] = s0;
        x[i+n_a] = s1;
        x[i+n_a2] = s2;
      }
    break;
  case 1: /* idendity */
    for (i=0; i<n_a3; i++)
      x[i] = y[i];
    break;
  case 2: /* diagonal of original matrix */
    for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i] = a[start]*y[i];
        x[i+n_a] = a[start]*y[i+n_a];
        x[i+n_a2] = a[start]*y[i+n_a2];
      }
    break;
  case 3: /* ILU of original matrix  */
  case 4: /* ILUT of original matrix  */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y[i];       /* diagonal entry */
      s1 = a[start]*y[i+n_a];
      s2 = a[start]*y[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
          s2 += a[k]*y[l+n_a2];
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
      x[i+n_a2] = s2;
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = x[i];                /* diagonal entry */
      s1 = x[i+n_a];
      s2 = x[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s0 += a[k]*x[l];    /* these x values are still unchanged */
          s1 += a[k]*x[l+n_a];/* because of l < i */
          s2 += a[k]*x[l+n_a2];
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
      x[i+n_a2] = s2;
    }
    break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

  /* scala with damping factor if necessary */
  alpha = global_sc->braess_sarazin_alpha;
  if (alpha!=1.0)
    {
      for (i=0; i<n_a3; i++)
        x[i] *= alpha;
    }
      

 /* multiplication with block B_0 and B_0^T */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
            x[l] += a[k]*y[n_a3+i];
          s0 += a[k]*y[l];
        }
      x[n_a3+i] = s0;     
    }
 /* multiplication with block B_1 and B_1^T*/

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
            x[l+n_a] += a[k]*y[n_a3+i];
          s0 += a[k]*y[l+n_a];
        }
      x[n_a3+i] += s0;
    }

 /* multiplication with block B_2 and B_2^T*/

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);

  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          ll =ara[l]; 
          if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
            x[l+n_a2] += a[k]*y[n_a3+i];
          s0 += a[k]*y[l+n_a2];
        }
      x[n_a3+i] += s0;
    }

  return(AMG_OK);
}

int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,n_b2,m_b2,l,ll;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1,s2;
  register int b,bb,n_a2,n_a3;
  double alpha,*s;
  int level;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b0!=m_b1)||(m_b0!=m_b2))
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=3*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (3*n_a!=n_b0+n_b1+n_b2)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_1 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  alpha = global_sc->braess_sarazin_alpha;  
  n_a2 = 2*n_a;
  n_a3 = 3*n_a;
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    for (i=0; i<n_a; i++)
    {
      start = ra[i];
      end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      s2 = a[start]*y[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l =ja[k]; 
        s0 += a[k]*y[l];
        s1 += a[k]*y[l+n_a];
        s2 += a[k]*y[l+n_a2];
      }
      x[i] -= alpha*s0;
      x[i+n_a] -= alpha*s1;
      x[i+n_a2] -= alpha*s2;
    }
    break;
  case 1: /* idendity */
    for (i=0; i<n_a3; i++)
      x[i] -= alpha*y[i];
    break;
  case 2: /* diagonal of original matrix */
    for (i=0; i<n_a; i++)
    {
      start = ra[i];
      x[i] -=alpha*a[start]*y[i];
      x[i+n_a] -=alpha*a[start]*y[i+n_a];
      x[i+n_a2] -=alpha*a[start]*y[i+n_a2];
    }
    break;
  case 3: /* ILU of original matrix  */
  case 4: /* ILUT of original matrix  */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* allocate memory (this is not nicely coded) */
    s = malloc(n_x*sizeof(double)); 
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s[i] = a[start]*y[i];       /* diagonal entry */
      s[i+n_a] = a[start]*y[i+n_a];
      s[i+n_a2] = a[start]*y[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s[i] += a[k]*y[l];
          s[i+n_a] += a[k]*y[l+n_a];
          s[i+n_a2] += a[k]*y[l+n_a2];
        }
      }
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s[i] += a[k]*s[l];        /* these s values are still unchanged */
          s[i+n_a] += a[k]*s[l+n_a];/* because of l < i */
          s[i+n_a2] += a[k]*s[l+n_a2];
        }
      }
      x[i] -= alpha*s[i];
      x[i+n_a] -= alpha*s[i+n_a];
      x[i+n_a2] -= alpha*s[i+n_a2];
    }
    free(s);
    break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }
 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll =ara[l]; 
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l] -= a[k]*y[n_a3+i];
      s0 += a[k]*y[l];
    } 
    x[n_a3+i] -= s0;
  }

 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll =ara[l]; 
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l+n_a] -= a[k]*y[n_a3+i];
      s0 += a[k]*y[l+n_a];
    }
    x[n_a3+i] -= s0;
  }
 /* multiplication with block B_2 */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);

  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll =ara[l]; 
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l+n_a2] -= a[k]*y[n_a3+i];
      s0 += a[k]*y[l+n_a2];
    }
    x[n_a3+i] -= s0;
  }

  return(AMG_OK);
}

int AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s0,s1,s2;
  register int b,bb,n_a2;
  double alpha;
  int level;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=3*n_a)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  n_a2 = 2*n_a;
  /* multiplication with block A */
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      s2 = a[start]*y[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a[k]*y[l];
        s1 += a[k]*y[l+n_a];
        s2 += a[k]*y[l+n_a2];
      }
      x[i] = s0;
      x[i+n_a] = s1;
      x[i+n_a2] = s2;
    }
    break;
  case 1: /* idendity */
    k = 3*n_a;
    for (i=0; i<k; i++)
      x[i] = y[i];
    break;
  case 2: /* diagonal of original matrix */
    for (i=0; i<n_a; i++)
    {
      start = ra[i];
      x[i] = a[start]*y[i];
      x[i+n_a] = a[start]*y[i+n_a];
      x[i+n_a2] = a[start]*y[i+n_a2];
    }
    break;
  case 3: /* ILU of original matrix  */
  case 4: /* ILUT of original matrix  */
   level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y[i];       /* diagonal entry */
      s1 = a[start]*y[i+n_a];
      s2 = a[start]*y[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
          s2 += a[k]*y[l+n_a2];
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
      x[i+n_a2] = s2;
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = x[i];                /* diagonal entry */
      s1 = x[i+n_a];
      s2 = x[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s0 += a[k]*x[l];    /* these x values are still unchanged */
          s1 += a[k]*x[l+n_a];/* because of l < i */
          s2 += a[k]*x[l+n_a2];/* because of l < i */
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
      x[i+n_a2] = s2;
    }
    break;
 default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

  /* scala with damping factor if necessary */
  alpha = global_sc->braess_sarazin_alpha;
  if (alpha!=1.0)
    {
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] *= alpha;
    }

  return(AMG_OK);
}

int AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                    AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,n_a2;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s0,s1,s2;
  register int b,bb;
  double alpha,*s;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=3*n_a)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  n_a2 = 2*n_a;
  /* multiplication with block A */
  
  alpha = global_sc->braess_sarazin_alpha;  
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      s2 = a[start]*y[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l =ja[k]; 
        s0 += a[k]*y[l];
        s1 += a[k]*y[l+n_a];
        s2 += a[k]*y[l+n_a2];
      }
      x[i] -= alpha*s0;
      x[i+n_a] -= alpha*s1;
      x[i+n_a2] -= alpha*s2;
    }
    break;
  case 1: /* idendity */
    k = 3*n_a;
    for (i=0; i<k; i++)
      x[i] -= alpha*y[i];
    break;
  case 2: /* diagonal of original matrix */
    for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i] -=alpha*a[start]*y[i];
        x[i+n_a] -=alpha*a[start]*y[i+n_a];
        x[i+n_a2] -=alpha*a[start]*y[i+n_a2];
      }
    break;
  case 3: /* ILU of original matrix  */
  case 4: /* ILUT of original matrix  */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* allocate memory (this is not nicely coded) */
    s = malloc(n_x*sizeof(double)); 
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s[i] = a[start]*y[i];       /* diagonal entry */
      s[i+n_a] = a[start]*y[i+n_a];
      s[i+n_a2] = a[start]*y[i+n_a2];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s[i] += a[k]*y[l];
          s[i+n_a] += a[k]*y[l+n_a];
          s[i+n_a2] += a[k]*y[l+n_a2];
        }
      }
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s[i] += a[k]*s[l];        /* these s values are still unchanged */
          s[i+n_a] += a[k]*s[l+n_a];/* because of l < i */
          s[i+n_a2] += a[k]*s[l+n_a2];
        }
      }
      x[i] -= alpha*s[i];
      x[i+n_a] -= alpha*s[i+n_a];
      x[i+n_a2] -= alpha*s[i+n_a2];
    }
    free(s);
    break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR BRAESS-SARAZIN SMOOTHER FOR 2D SADDLE POINT */
/* PROBLEMS OF TYPE 2                                                       */
/*                                                                          */
/* ( al C(A)    0         B_0)                                              */
/* (     0    al C(A)     B_1) x = y                                        */
/* (B_2  B_3     0  )                                                       */
/*                                                                          */
/* Input : A , B_0, B_1, B_2, B_3, x, y                                     */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  double alpha;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_2 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_2 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with block A */
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    for (i=0; i<n_a; i++)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = a[start]*y[i];
        s1 = a[start]*y[i+n_a];
        for (k=start+1; k<end; k++) 
          {
            l = ja[k];
            s0 += a[k]*y[l];
            s1 += a[k]*y[l+n_a];
          }
        x[i] = s0;
        x[i+n_a] = s1;
      }
    break;
  case 1: /* idendity */
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] = y[i];
    break;
  case 2: /* diagonal of original matrix */
    for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i] = a[start]*y[i];
        x[i+n_a] = a[start]*y[i+n_a];
      }
    break;
  case 3: /* ILU of original matrix  */
  case 4: /* ILUT of original matrix  */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y[i];       /* diagonal entry */
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = x[i];                /* diagonal entry */
      s1 = x[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s0 += a[k]*x[l];    /* these x values are still unchanged */
          s1 += a[k]*x[l+n_a];/* because of l < i */
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
    }
    break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

  /* scala with damping factor if necessary */
  alpha = global_sc->braess_sarazin_alpha;
  if (alpha!=1.0)
    {
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] *= alpha;
    }
      
 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[2*n_a+l];
        }
      x[i] += s0;
   }
 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[2*n_a+l];
        }
      x[i+n_a] += s0;
   }

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] = s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] += s0;
    }
  return(AMG_OK);
}

int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A,
                                                  AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  double alpha,*s;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_2 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_2 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with block A */
  
  alpha = global_sc->braess_sarazin_alpha;  
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l =ja[k]; 
        s0 += a[k]*y[l];
        s1 += a[k]*y[l+n_a];
      }
      x[i] -= alpha*s0;
      x[i+n_a] -= alpha*s1;
    }
    break;
  case 1: /* idendity */
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] -= alpha*y[i];
    break;
  case 2: /* diagonal of original matrix */
    for (i=0; i<n_a; i++)
    {
      start = ra[i];
      x[i] -=alpha*a[start]*y[i];
      x[i+n_a] -=alpha*a[start]*y[i+n_a];
    }
    break;
  case 3: /* ILU of original matrix  */
  case 4: /* ILUT of original matrix  */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* allocate memory (this is not nicely coded) */
    s = malloc(n_x*sizeof(double)); 
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s[i] = a[start]*y[i];       /* diagonal entry */
      s[i+n_a] = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s[i] += a[k]*y[l];
          s[i+n_a] += a[k]*y[l+n_a];
        }
      }
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s[i] += a[k]*s[l];        /* these s values are still unchanged */
          s[i+n_a] += a[k]*s[l+n_a];/* because of l < i */
        }
      }
      x[i] -= alpha*s[i];
      x[i+n_a] -= alpha*s[i+n_a];
    }
    free(s);
    break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[2*n_a+l];
        }
      x[i] -= s0;
   }
 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[2*n_a+l];
        }
      x[i+n_a] -= s0;
   }

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] -= s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] -= s0;
    }
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR BRAESS-SARAZIN SMOOTHERS FOR 2D SADDLE      */
/* POINT PROBLEMS OF TYPE 3                                                 */
/*                                                                          */
/* ( f(A)  f(B2)  B0)                                                       */
/* ( f(B3) f(B4)  B1) x = y                                                 */
/* ( B0^T  B1^T   0  )                                                      */
/*                                                                          */
/* Input : A , B0^T , B1^T , B2, B3, B4 x, y                                */
/*                                                                          */
/****************************************************************************/
int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l,ll;
  register double *x, *y, *a, *xx, *aa, *yy, *a22;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  double alpha;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_3 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_3 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
   /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with block A */
  
  switch(global_sc->braess_sarazin_matrix)
  {
    case 0: /* original matrix (with alpha = 1) */
      AMG_dmatmul_SADDLE_3_TYPE_3(x_,A,B,y_);
      return(AMG_OK);
    break;
    case 1: /* idendity */
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] = y[i];
      break;
    case 2: /* diagonal of original matrix */
      level = A->level;
      if (level==0)  /* finest level */
      {
        a22 = AMG_MATRIX_A(B[4]);
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] = a[start]*y[i];
          x[i+n_a] = a22[start]*y[i+n_a];
        }
      }
      else /* coarser levels */
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] = a[start]*y[i];
          x[i+n_a] = a[start]*y[i+n_a];
        }
      break;
    case 3: /* ILU of original matrix (left upper block) */
    case 4: /* ILUT of original matrix (left upper block) */
      level = A->level;
      a = AMG_MATRIX_A(M[level]);
      ra = AMG_MATRIX_RA(M[level]);
      ja = AMG_MATRIX_JA(M[level]);
      /* x = Uy */
      for (i=0; i<n_a; i++)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = a[start]*y[i];       /* diagonal entry */
        s1 = a[start]*y[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l>i)              /* column > row */
          {
            s0 += a[k]*y[l];
            s1 += a[k]*y[l+n_a];
          }
        }
        x[i] = s0;
        x[i+n_a] = s1;
      }
      /* x = Ly */
      for (i=n_a-1; i>=0; i--)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = x[i];                /* diagonal entry */
        s1 = x[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i)              /* column < row */
          {
            s0 += a[k]*x[l];    /* these x values are still unchanged */
            s1 += a[k]*x[l+n_a];/* because of l < i */
          }
        }
        x[i] = s0;
        x[i+n_a] = s1;
      }
      break;
    default:
      AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
      exit(4711);
      break;
  }

  /* scala with damping factor if necessary */
  alpha = global_sc->braess_sarazin_alpha;
  if (alpha!=1.0)
  {
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] *= alpha;
  }
  
  
 /* multiplication with block B_0 and B_0^T */
  
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll =ara[l]; 
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
        x[l] += a[k]*y[2*n_a+i];
      s0 += a[k]*y[l];
    }
    x[2*n_a+i] = s0;     
  }
 /* multiplication with block B_1 and B_1^T*/

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll =ara[l]; 
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Dirichlet */
        x[l+n_a] += a[k]*y[2*n_a+i];
      s0 += a[k]*y[l+n_a];
    }
    x[2*n_a+i] += s0;
  }
  
  return(AMG_OK);
}

int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l,ll;
  register double *x, *y, *a, *xx, *aa, *yy, *a22;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  double alpha,*s;
  int level;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_3 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b0)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_3 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  alpha = global_sc->braess_sarazin_alpha;  
  
  switch(global_sc->braess_sarazin_matrix)
  {
    case 0: /* original matrix (with alpha = 1.0) */
      AMG_dmatminus_SADDLE_3_TYPE_3(x_,A,B,y_);
      return(AMG_OK);
      break;
    case 1: /* idendity */
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] -= alpha*y[i];
      break;
    case 2: /* diagonal of original matrix */
      level = A->level;
      if (level==0) /* finest level */
      {
        a22 = AMG_MATRIX_A(B[4]);
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i]  -=alpha*a[start]*y[i];
          x[i+n_a] -=alpha*a22[start]*y[i];
        }
      }
      else /* coarser levels */
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] -=alpha*a[start]*y[i];
          x[i+n_a] -=alpha*a[start]*y[i+n_a];
        }
      break;
    case 3: /* ILU of original matrix (left upper block) */
    case 4: /* ILUT of original matrix (left upper block) */
      level = A->level;
      a = AMG_MATRIX_A(M[level]);
      ra = AMG_MATRIX_RA(M[level]);
      ja = AMG_MATRIX_JA(M[level]);
      /* allocate memory (this is not nicely coded) */
      s = malloc(n_x*sizeof(double)); 
      /* x = Uy */
      for (i=0; i<n_a; i++)
      {
        start = ra[i]; 
        end = start+ja[start];
        s[i] = a[start]*y[i];       /* diagonal entry */
        s[i+n_a] = a[start]*y[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l>i)              /* column > row */
          {
            s[i] += a[k]*y[l];
            s[i+n_a] += a[k]*y[l+n_a];
          }
        }
      }
      /* x = Ly */
      for (i=n_a-1; i>=0; i--)
      {
        start = ra[i]; 
        end = start+ja[start];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i)              /* column < row */
          {
            s[i] += a[k]*s[l];        /* these s values are still unchanged */
            s[i+n_a] += a[k]*s[l+n_a];/* because of l < i */
          }
        }
        x[i] -= alpha*s[i];
        x[i+n_a] -= alpha*s[i+n_a];
      }
      free(s);
      break;
    default:
      AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
      exit(4711);
      break;
  }
  /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll =ara[l]; 
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l] -= a[k]*y[2*n_a+i];
      s0 += a[k]*y[l];
    } 
    x[2*n_a+i] -= s0;
  }
  
  /* multiplication with block B_1 */
  
  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);
  
  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      ll =ara[l]; 
      if (aja[ll] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
        x[l+n_a] -= a[k]*y[2*n_a+i];
      s0 += a[k]*y[l+n_a];
    }
    x[2*n_a+i] -= s0;
  }
  return(AMG_OK);
}

int AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l;
  register double *x, *y, *a, *xx, *aa, *yy, *a22;
  register int *ra, *ja;
  register double s0,s1;
  register int b,bb;
  double alpha;
  int level;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* multiplication with block A */
  
  switch(global_sc->braess_sarazin_matrix)
  {
    case 0: /* original matrix (with alpha = 1) */
      AMG_A_dmatmul_SADDLE_3_TYPE_3(x_,A,B,y_);
      return(AMG_OK);
      break;
    case 1: /* idendity */
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] = y[i];
      break;
    case 2: /* diagonal of original matrix */
      level = A->level;
      if (level==0)
      {
        a22 = AMG_MATRIX_A(B[4]);
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] = a[start]*y[i];
          x[i+n_a] = a22[start]*y[i+n_a];
        }
      }
      else
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] = a[start]*y[i];
          x[i+n_a] = a[start]*y[i+n_a];
        }
      break;
    case 3: /* ILU of original matrix (left upper block) */
    case 4: /* ILUT of original matrix (left upper block) */
      level = A->level;
      a = AMG_MATRIX_A(M[level]);
      ra = AMG_MATRIX_RA(M[level]);
      ja = AMG_MATRIX_JA(M[level]);
      /* x = Uy */
      for (i=0; i<n_a; i++)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = a[start]*y[i];       /* diagonal entry */
        s1 = a[start]*y[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l>i)              /* column > row */
          {
            s0 += a[k]*y[l];
            s1 += a[k]*y[l+n_a];
          }
        }
        x[i] = s0;
        x[i+n_a] = s1;
      }
      /* x = Ly */
      for (i=n_a-1; i>=0; i--)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = x[i];                /* diagonal entry */
        s1 = x[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i)              /* column < row */
          {
            s0 += a[k]*x[l];    /* these x values are still unchanged */
            s1 += a[k]*x[l+n_a];/* because of l < i */
          }
        }
        x[i] = s0;
        x[i+n_a] = s1;
      }
      break;
    default:
      AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
      exit(4711);
      break;
  }
  
  /* scala with damping factor if necessary */
  alpha = global_sc->braess_sarazin_alpha;
  if (alpha!=1.0)
  {
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] *= alpha;
  }
  
  return(AMG_OK);
}

int AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                    AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,l;
  register double *x, *y, *a, *xx, *aa, *yy, *a22;
  register int *ra, *ja;
  register double s0,s1;
  register int b,bb;
  double alpha,*s;
  int level;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */
  
  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatminus BRAESS_SARAZIN_SADDLE_3_TYPE_3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* multiplication with block A */
  
  alpha = global_sc->braess_sarazin_alpha;  
  
  switch(global_sc->braess_sarazin_matrix)
  {
    case 0: /* original matrix (with alpha = 1) */
      AMG_A_dmatminus_SADDLE_3_TYPE_3(x_,A,B,y_);
      break;
    case 1: /* idendity */
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] -= alpha*y[i];
      break;
    case 2: /* diagonal of original matrix */
      level = A->level;
      if (level==0)
      {
        a22 = AMG_MATRIX_A(B[4]);
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i]  -=alpha*a[start]*y[i];
          x[i+n_a] -=alpha*a22[start]*y[i];
        }
      }
      else
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] -=alpha*a[start]*y[i];
          x[i+n_a] -=alpha*a[start]*y[i+n_a];
        }
      break;
    case 3: /* ILU of original matrix  */
    case 4: /* ILUT of original matrix  */
      level = A->level;
      a = AMG_MATRIX_A(M[level]);
      ra = AMG_MATRIX_RA(M[level]);
      ja = AMG_MATRIX_JA(M[level]);
      /* allocate memory (this is not nicely coded) */
      s = malloc(n_x*sizeof(double)); 
      /* x = Uy */
      for (i=0; i<n_a; i++)
      {
        start = ra[i]; 
        end = start+ja[start];
        s[i] = a[start]*y[i];       /* diagonal entry */
        s[i+n_a] = a[start]*y[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l>i)              /* column > row */
          {
            s[i] += a[k]*y[l];
            s[i+n_a] += a[k]*y[l+n_a];
          }
        }
      }
      /* x = Ly */
      for (i=n_a-1; i>=0; i--)
      {
        start = ra[i]; 
        end = start+ja[start];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i)              /* column < row */
          {
            s[i] += a[k]*s[l];        /* these s values are still unchanged */
            s[i+n_a] += a[k]*s[l+n_a];/* because of l < i */
          }
        }
        x[i] -= alpha*s[i];
        x[i+n_a] -= alpha*s[i+n_a];
      }
      free(s);
      break;
    default:
      AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
      exit(4711);
      break;
  }
  
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR BRAESS-SARAZIN SMOOTHERS FOR 2D SADDLE      */
/* POINT PROBLEMS OF TYPE 4                                                 */
/*                                                                          */
/* ( f(A)  f(B4)     B_0)                                                   */
/* ( f(B5) f(B6)     B_1) x = y                                             */
/* ( B_2    B_3      0  )                                                   */
/*                                                                          */
/* Input : A , B_0, B_1, B_2, B_3, B_4, B_5, B_6, x, y                      */
/*                                                                          */
/****************************************************************************/
int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a22;
  register double s0,s1;
  register int b,bb;
  double alpha;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with block A */
  
  switch(global_sc->braess_sarazin_matrix)
  {
    case 0: /* original matrix */
      AMG_dmatmul_SADDLE_3_TYPE_4(x_,A,B,y_);
      return(AMG_OK);
    break;
    case 1: /* idendity */
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] = y[i];
      break;
    case 2: /* diagonal of original matrix */
      level = A->level;
      if (level==0)
      {
        a22 = AMG_MATRIX_A(B[6]);
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] = a[start]*y[i];
          x[i+n_a] = a22[start]*y[i+n_a];
        }
      }
      else
        for (i=0; i<n_a; i++)
        {
          start = ra[i];
          x[i] = a[start]*y[i];
          x[i+n_a] = a[start]*y[i+n_a];
        }
      break;
    case 3: /* ILU of original matrix (left upper block) */
    case 4: /* ILUT of original matrix (left upper block) */
      level = A->level;
      a = AMG_MATRIX_A(M[level]);
      ra = AMG_MATRIX_RA(M[level]);
      ja = AMG_MATRIX_JA(M[level]);
      /* x = Uy */
      for (i=0; i<n_a; i++)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = a[start]*y[i];       /* diagonal entry */
        s1 = a[start]*y[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l>i)              /* column > row */
          {
            s0 += a[k]*y[l];
            s1 += a[k]*y[l+n_a];
          }
        }
        x[i] = s0;
        x[i+n_a] = s1;
      }
            
      /* x = Ly */
      for (i=n_a-1; i>=0; i--)
      {
        start = ra[i]; 
        end = start+ja[start];
        s0 = x[i];                /* diagonal entry */
        s1 = x[i+n_a];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i)              /* column < row */
          {
            s0 += a[k]*x[l];    /* these x values are still unchanged */
            s1 += a[k]*x[l+n_a];/* because of l < i */
          }
        }
        x[i] = s0;
        x[i+n_a] = s1;
      }
      break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

  /* scala with damping factor if necessary */
  alpha = global_sc->braess_sarazin_alpha;
  if (alpha!=1.0)
  {
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] *= alpha;
  }
 
  /* multiplication with block B_0 */
  
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);
  
  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
  {
    /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
    if (i>=A->active) /* no active node */
      continue;
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;      
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      s0 += a[k]*y[2*n_a+l];
    }
    x[i] += s0;
  }
  /* multiplication with block B_1 */
  
  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);
  
  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
  {
    /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
    if (i>=A->active) /* no active node */
      continue;
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      s0+=a[k]*y[2*n_a+l];
    }
    x[i+n_a] += s0;
  }
  
 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      s0 += a[k]*y[l];
    }
    x[2*n_a+i] = s0;     
  }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
  {
    start=ra[i]; 
    end = ra[i+1]-1;
    s0 = 0;
    for (k=start; k<=end; k++)
    {
      l = ja[k];
      s0 += a[k]*y[l+n_a];
    }
    x[2*n_a+i] += s0;
  }
  return(AMG_OK);
}

int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *a, *xx, *aa, *yy;
  register double *a11,*a12,*a21,*a22;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  int level;
  double alpha, *s;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2)
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  alpha = global_sc->braess_sarazin_alpha;
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    AMG_dmatminus_SADDLE_3_TYPE_4(x_,A,B,y_);
    return(AMG_OK);
    break;
  case 1: /* idendity */
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] -= alpha*y[i];
    break;
  case 2: /* diagonal of original matrix */
    level = A->level;
    if (level==0)
    {
      a22 = AMG_MATRIX_A(B[4]);
      for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i]  -=alpha*a[start]*y[i];
        x[i+n_a] -=alpha*a22[start]*y[i];
      }
    }
    else
      for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i] -=alpha*a[start]*y[i];
        x[i+n_a] -=alpha*a[start]*y[i+n_a];
      }
    break;
  case 3: /* ILU of original matrix (left upper block) */
  case 4: /* ILUT of original matrix (left upper block) */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* allocate memory (this is not nicely coded) */
    s = malloc(n_x*sizeof(double)); 
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s[i] = a[start]*y[i];       /* diagonal entry */
      s[i+n_a] = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s[i] += a[k]*y[l];
          s[i+n_a] += a[k]*y[l+n_a];
        }
      }
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s[i] += a[k]*s[l];        /* these s values are still unchanged */
          s[i+n_a] += a[k]*s[l+n_a];/* because of l < i */
        }
      }
      x[i] -= alpha*s[i];
      x[i+n_a] -= alpha*s[i+n_a];
    }
    free(s);
    break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

  /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;      
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[2*n_a+l];
        }
      x[i] -= s0;
   }
  /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[1]);
  ra = AMG_MATRIX_RA(B[1]);
  ja = AMG_MATRIX_JA(B[1]);

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      /*if (aja[ara[i]] == 1)*/ /* row in A has only one entry -> Dirichlet */
      if (i>=A->active) /* no active node */
        continue;
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0+=a[k]*y[2*n_a+l];
        }
      x[i+n_a] -= s0;
   }

 /* multiplication with block B_2^T */

  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  for (i=0; i<m_b2; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
        }
      x[2*n_a+i] -= s0;     
    }

 /* multiplication with block B_3^T */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);

  for (i=0; i<m_b3; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l+n_a];
        }
      x[2*n_a+i] -= s0;
    }
  return(AMG_OK);
}

int AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                   AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a,*a22;
  register double s0,s1;
  register int b,bb;
  int level;
  double alpha;


  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatmul BRAESS_SARAZIN_SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  
  /* multiplication with block A */
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    AMG_A_dmatmul_SADDLE_3_TYPE_4(x_,A,B,y_);
    return(AMG_OK);
    break;
  case 1: /* idendity */
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] = y[i];
    break;
  case 2: /* diagonal of original matrix */
    level = A->level;
    if (level==0)
    {
      a22 = AMG_MATRIX_A(B[6]);
      for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i] = a[start]*y[i];
        x[i+n_a] = a22[start]*y[i+n_a];
      }
    }
    else
      for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i] = a[start]*y[i];
        x[i+n_a] = a[start]*y[i+n_a];
      }
    break;
  case 3: /* ILU of original matrix (left upper block) */
  case 4: /* ILUT of original matrix (left upper block) */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y[i];       /* diagonal entry */
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = x[i];                /* diagonal entry */
      s1 = x[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s0 += a[k]*x[l];    /* these x values are still unchanged */
          s1 += a[k]*x[l+n_a];/* because of l < i */
        }
      }
      x[i] = s0;
      x[i+n_a] = s1;
    }
    break;
    default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }

  /* scala with damping factor if necessary */
  alpha = global_sc->braess_sarazin_alpha;
  if (alpha!=1.0)
    {
      k = 2*n_a;
      for (i=0; i<k; i++)
        x[i] *= alpha;
    }

  return(AMG_OK);
}
int AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                   AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l,ll;
  register int n_b0,n_b1,n_b2,n_b3,m_b0,m_b1,m_b2,m_b3;
  register double *x, *y, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double *a, *a22;
  register double s0,s1;
  register int b,bb;
  double alpha, *s;
  int level;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_4 - vector lengths mismatched !!\n");
        exit(4711);
     }
     
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: A_dmatminus SADDLE_3_TYPE_4 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);

  /* multiplication with A */
  alpha = global_sc->braess_sarazin_alpha;  
  
  switch(global_sc->braess_sarazin_matrix)
  {
  case 0: /* original matrix */
    AMG_A_dmatminus_SADDLE_3_TYPE_4(x_,A,B,y_);
    break;
  case 1: /* idendity */
    k = 2*n_a;
    for (i=0; i<k; i++)
      x[i] -= alpha*y[i];
    break;
  case 2: /* diagonal of original matrix */
    level = A->level;
    if (level==0)
    {
      a22 = AMG_MATRIX_A(B[6]);
      for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i]  -=alpha*a[start]*y[i];
        x[i+n_a] -=alpha*a22[start]*y[i];
      }
    }
    else
      for (i=0; i<n_a; i++)
      {
        start = ra[i];
        x[i] -=alpha*a[start]*y[i];
        x[i+n_a] -=alpha*a[start]*y[i+n_a];
      }
    break;
  case 3: /* ILU of original matrix  */
  case 4: /* ILUT of original matrix  */
    level = A->level;
    a = AMG_MATRIX_A(M[level]);
    ra = AMG_MATRIX_RA(M[level]);
    ja = AMG_MATRIX_JA(M[level]);
    /* allocate memory (this is not nicely coded) */
    s = malloc(n_x*sizeof(double)); 
    /* x = Uy */
    for (i=0; i<n_a; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s[i] = a[start]*y[i];       /* diagonal entry */
      s[i+n_a] = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l>i)              /* column > row */
        {
          s[i] += a[k]*y[l];
          s[i+n_a] += a[k]*y[l+n_a];
        }
      }
    }
    /* x = Ly */
    for (i=n_a-1; i>=0; i--)
    {
      start = ra[i]; 
      end = start+ja[start];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        if (l<i)              /* column < row */
        {
          s[i] += a[k]*s[l];        /* these s values are still unchanged */
          s[i+n_a] += a[k]*s[l+n_a];/* because of l < i */
        }
      }
      x[i] -= alpha*s[i];
      x[i+n_a] -= alpha*s[i+n_a];
    }
    free(s);
    break;
  default:
    AMG_Print("Matrix in Braess-Sarazin smoother UNKNOWN !!!\n");
    exit(4711);
    break;
  }
  return(AMG_OK);
}
