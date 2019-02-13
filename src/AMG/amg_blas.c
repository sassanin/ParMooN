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
/* File:          amg_blas.c                                                */
/*                                                                          */
/* Purpose:   BLAS (Basic Linear Algebra Subroutines) for amg package       */
/*                        Contains Level 1 and 2                            */
/*                                                                          */
/* Author:        Peter Bastian                                             */
/*                        Institut fuer Computeranwendungen III             */
/*                        Universitaet Stuttgart                            */
/*                        Pfaffenwaldring 27                                */
/*                        70550 Stuttgart                                   */
/*                        email: peter@ica3.uni-stuttgart.de                */
/*                        phone: 0049-(0)711-685-7003                       */
/*                        fax  : 0049-(0)711-685-7000                       */
/*                                                                          */
/* History:   05 FEB 1996 Begin                                             */
/*            02 OKT 1997 redesign                                          */
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
/* History:   1998/02/19 start using this library for MooN_MD               */
/*                                                                          */
/* Remarks:   1998/05/27 start implementing matrix vector routines          */
/*                       for saddle point problems                          */
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
/* defines in the following order                                           */
/*                                                                          */
/*    compile time constants defining static data size (i.e. arrays)        */
/*    other constants                                                       */
/*    macros                                                                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/* in the corresponding include file!)                                      */
/*                                                                          */
/****************************************************************************/

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
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*D
   blas - BLAS routines based on the amg data structure 
   
   SYNOPSIS:
   int    AMG_dset      (AMG_VECTOR *x, double a );
   int    AMG_randomize (AMG_VECTOR *x );
   int    AMG_dcopy     (AMG_VECTOR *x, AMG_VECTOR *y );
   int    AMG_dscale    (AMG_VECTOR *x, double a );
   int    AMG_daxpy     (AMG_VECTOR *x, double a, AMG_VECTOR *y);
   double AMG_ddot      (AMG_VECTOR *x, AMG_VECTOR *y);
   int    AMG_dmatset   (AMG_MATRIX *A, double a);
   int    AMG_dmatcopy  (AMG_MATRIX *A, AMG_MATRIX *B );
   int    AMG_dmatmul   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_VECTOR *y);
   int    AMG_dmatminus (AMG_VECTOR *x, AMG_MATRIX *A, AMG_VECTOR *y);

   
   PARAMETERS:
.  x,y - vectors
.  A,B - sparse matrices
.  a - scalar
   
   DESCRIPTION:
.n AMG_dset sets all values of vector x to a.
.n AMG_randomize sets random values in vector x.
.n AMG_dcopy copies vector x to vector y.
.n AMG_dscale multiplies each entry of x with a.
.n AMG_daxpy computes x = x + ay
.n AMG_ddot computes euclidean scalar product of two vectors
.n AMG_dmatset sets al entries of matrix A to a.
.n AMG_dmatcopy copies matrix B to A (same structure assumed)
.n AMG_dmatmul computes x = Ay
.n AMG_dmatminus computes x = x - Ay

   All routines check the compatibility of their arguments and return
   errors if necessary. Where appropriate the routines have been 
   implemented in a scalar and a block version to enhance efficiency
   in the scalar case.
   
   RETURN VALUE:
.n AMG_OK
.n AMG_FATAL 
   
D*/
/****************************************************************************/

int AMG_dset (AMG_VECTOR *x, double a)
{
        register int i,n;
        register double *values;
        
        n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
        values = AMG_VECTOR_X(x);
        for (i=0; i<n; i++) *values++ = a;
        
        return(AMG_OK);
}

int AMG_randomize (AMG_VECTOR *x)
{
        register int i,n;
        register double *values;
        
        n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
        values = AMG_VECTOR_X(x);
        for (i=0; i<n; i++) *values++ = rand();
        
        return(AMG_OK);
}

int AMG_dcopy (AMG_VECTOR *x, AMG_VECTOR *y)
{
  register int i,n;
  register double *values_x, *values_y;
  
  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);

  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(y)) 
     {
        printf("x %d y %d \n",AMG_VECTOR_N(x),AMG_VECTOR_N(y));
        AMG_Print("ERROR: AMG_dcopy - vector lengths mismatched !!\n");
        exit(4711);
     }
  
     if (AMG_VECTOR_B(x)!=AMG_VECTOR_B(y)) return(AMG_FATAL);
  }

  values_x = AMG_VECTOR_X(x);
  values_y = AMG_VECTOR_X(y);
  for (i=0; i<n; i++) *values_x++ = *values_y++;
  
  return(AMG_OK);
}

int AMG_dscale (AMG_VECTOR *x, double a)
{
        register int i,n;
        register double *values;
        
        n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
        values = AMG_VECTOR_X(x);
        for (i=0; i<n; i++) *values++ *= a;
        
        return(AMG_OK);
}

int AMG_daxpy (AMG_VECTOR *x, double a, AMG_VECTOR *y)
{
  register int i,n;
  register double *values_x, *values_y;
  
  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(y)) return(AMG_FATAL);
     if (AMG_VECTOR_B(x)!=AMG_VECTOR_B(y)) return(AMG_FATAL);
  }
  values_x = AMG_VECTOR_X(x);
  values_y = AMG_VECTOR_X(y);
  if (a==1.0)
    {
      for (i=0; i<n; i++) *values_x++ += (*values_y++);
      return(AMG_OK);
    }
  if (a==-1.0)
    {
      for (i=0; i<n; i++) *values_x++ -= (*values_y++);
      return(AMG_OK);
    }

  for (i=0; i<n; i++) *values_x++ += a* (*values_y++);
  return(AMG_OK);
}

int AMG_daxpby (AMG_VECTOR *x, double a, AMG_VECTOR *y, double b, AMG_VECTOR *z)
{
  register int i,n;
  register double *values_x, *values_y, *values_z;

  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(y)) return(AMG_FATAL);
     if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(z)) return(AMG_FATAL);
     if (AMG_VECTOR_B(x)!=AMG_VECTOR_B(y)) return(AMG_FATAL);
  }
  values_x = AMG_VECTOR_X(x);
  values_y = AMG_VECTOR_X(y);
  values_z = AMG_VECTOR_X(z);
 
  if (a==1.0)
    {
      if (b==1.0)
        {
          for (i=0; i<n; i++) *values_x++ = (*values_y++) + (*values_z++);
          return(AMG_OK);
        }
      if (b==-1.0)
        {
          for (i=0; i<n; i++) *values_x++ = (*values_y++) - (*values_z++);
          return(AMG_OK);
        }
      for (i=0; i<n; i++) *values_x++ = (*values_y++) + b * (*values_z++);
      return(AMG_OK);
    }
  if (b==1.0)
    {
      if (a==-1.0)
        {
          for (i=0; i<n; i++) *values_x++ = - (*values_y++) + (*values_z++);
          return(AMG_OK);
        }
      for (i=0; i<n; i++) *values_x++ = a * (*values_y++) + (*values_z++);
      return(AMG_OK);
    }

  for (i=0; i<n; i++) *values_x++ = a * (*values_y++) + b * (*values_z++);
  
  return(AMG_OK);
}

double AMG_ddot (AMG_VECTOR *x, AMG_VECTOR *y)
{
  register int i,n;
  register double *values_x, *values_y;
  register double s=0.0;
  
  n = AMG_VECTOR_N(x)*AMG_VECTOR_B(x);
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x)!=AMG_VECTOR_N(y)) 
     {
        AMG_Print("vector length mismatched in AMG_ddot!\n");
        exit(4711);
     }
  }
  values_x = AMG_VECTOR_X(x);
  values_y = AMG_VECTOR_X(y);
  for (i=0; i<n; i++) s += (*values_x++) * (*values_y++);
        
  return(s);
}

int AMG_dmatset (AMG_MATRIX *A, double a)
{
  int size,i;
  double *values;
  
  size = AMG_MATRIX_N(A)*AMG_MATRIX_BB(A);
  values = AMG_MATRIX_A(A);
  for (i=0; i<size; i++) *values++ = a;
  
  return(AMG_OK);
}


int AMG_dmatcopy (AMG_MATRIX *A, AMG_MATRIX *B)
{
  int size_a,size_b,i;
  double *values_a, *values_b;
  
  size_a = AMG_MATRIX_N(A)*AMG_MATRIX_BB(A);
  size_b = AMG_MATRIX_N(B)*AMG_MATRIX_BB(B);
  if (global_sc->verbose>1)
     if (size_a!=size_b) return(AMG_FATAL);
  
  values_a = AMG_MATRIX_A(A);
  values_b = AMG_MATRIX_A(B);
  for (i=0; i<size_a; i++) *values_a++ = *values_b++;
  
  return(AMG_OK);
}

#define ZERO2(x)    x[0] = x[1] = 0.0;
#define ZERO3(x)    x[0] = x[1] = x[2] = 0.0;
#define ZERO4(x)    x[0] = x[1] = x[2] = x[3] = 0.0;


#define XAY2(x,a,y)        {x[0] += a[0]*y[0] + a[1]*y[1]; \
                    x[1] += a[2]*y[0] + a[3]*y[1];}

#define XAY3(x,a,y)        {x[0] += a[0]*y[0] + a[1]*y[1] + a[2]*y[2]; \
                    x[1] += a[3]*y[0] + a[4]*y[1] + a[5]*y[2]; \
                    x[2] += a[6]*y[0] + a[7]*y[1] + a[8]*y[2];}

#define XAY4(x,a,y)        {x[0] += a[0]*y[0] + a[1]*y[1] + a[2]*y[2] + a[3]*y[3]; \
                    x[1] += a[4]*y[0] + a[5]*y[1] + a[6]*y[2] + a[7]*y[3]; \
                    x[2] += a[8]*y[0] + a[9]*y[1] + a[10]*y[2] + a[11]*y[3]; \
                    x[3] += a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];} \

/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR SCALAR SYSTEMS                              */
/*                                                                          */
/****************************************************************************/

int AMG_dmatmul (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,start,end;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, l;
  register double s;
  register int b,bb;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_B(x_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
     if (AMG_VECTOR_B(y_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
     if (AMG_VECTOR_N(x_)!=AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatmul SCALAR - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  n = AMG_VECTOR_N(x_);
  b = AMG_VECTOR_B(x_);
  bb = AMG_MATRIX_BB(A);
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (b)
    {
    case 1:
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          s = a[start]*y[i];
          for (k=start+1; k<end; k++) 
          {
             l = ja[k];
             s+= a[k]*y[l];
          }
          x[i] = s;
        }
      break;
      
    case 2:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO2(xx);
          yy = y+(i*b); XAY2(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAY2(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;                
      
    case 3:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO3(xx);
          yy = y+(i*b); XAY3(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAY3(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;                

    case 4:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO4(xx);
          yy = y+(i*b); XAY4(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAY4(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;
      
    default:
      AMG_Print("dmatmul: blocksize>4 not implemented yet\n");
      break;        
    } 
  return(AMG_OK);
}

#define XAYM2(x,a,y)   {x[0] -= a[0]*y[0] + a[1]*y[1]; \
                            x[1] -= a[2]*y[0] + a[3]*y[1];}

#define XAYM3(x,a,y)   {x[0] -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2]; \
                            x[1] -= a[3]*y[0] + a[4]*y[1] + a[5]*y[2]; \
                            x[2] -= a[6]*y[0] + a[7]*y[1] + a[8]*y[2];}

#define XAYM4(x,a,y)   {x[0] -= a[0]*y[0] + a[1]*y[1] + a[2]*y[2] + a[3]*y[3]; \
                            x[1] -= a[4]*y[0] + a[5]*y[1] + a[6]*y[2] + a[7]*y[3]; \
                            x[2] -= a[8]*y[0] + a[9]*y[1] + a[10]*y[2] + a[11]*y[3]; \
                            x[3] -= a[12]*y[0] + a[13]*y[1] + a[14]*y[2] + a[15]*y[3];} \


int AMG_dmatminus (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,start,end, l;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x_)!=AMG_MATRIX_N(A))
     {
        printf("x %d A %d \n",AMG_VECTOR_N(x_),AMG_MATRIX_N(A));
        AMG_Print("ERROR: dmatminus SCALAR - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=AMG_MATRIX_N(A))
     {
        printf("y %d A %d \n",AMG_VECTOR_N(y_),AMG_MATRIX_N(A));
        AMG_Print("ERROR: dmatminus SCALAR - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatminus SCALAR - vector lengths mismatched !!\n");
        exit(4711);
     }

     if (AMG_VECTOR_B(x_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
     if (AMG_VECTOR_B(y_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
  }

  /* prepare data */
  n = AMG_VECTOR_N(x_);
  b = AMG_VECTOR_B(x_);
  bb = AMG_MATRIX_BB(A);
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (b)
    {
    case 1:
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          s = a[start]*y[i];
          for (k=start+1; k<end; k++)
          {
             l =ja[k];
             s += a[k]*y[ja[k]];
          }
          x[i] -= s;
        }
      break;
      
    case 2:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO2(xx);
          yy = y+(i*b); XAYM2(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAYM2(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;                
      
    case 3:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO3(xx);
          yy = y+(i*b); XAYM3(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAYM3(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;                
      
    case 4:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO4(xx);
          yy = y+(i*b); XAYM4(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAYM4(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;
      
    default:
	AMG_Print("dmatmul: blocksize>4 not implemented yet\n");
	break;        
    }
  
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR SCALAR SYSTEMS WITH TWO RIGHT HAND SIDES    */
/*                                                                          */
/****************************************************************************/

int AMG_dmatmul_SCALAR2 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                         AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,l,start,end;
  register double *x0, *y0, *x1, *y1, *a;
  register int *ra, *ja;
  register double s0, s1;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x_)!=2*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=2*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatmul SCALAR2 - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  n = AMG_MATRIX_N(A);
  x0 = AMG_VECTOR_X(x_);
  x1 = x0 + n;
  y0 = AMG_VECTOR_X(y_);
  y1 = y0 + n;
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<n; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y0[i];
      s1 = a[start]*y1[i];
      /*printf("a(%d,%d) = %g;\n", i+1,i+1,a[start]);*/
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a[k]*y0[l];
        s1 += a[k]*y1[l];
        /*printf("a(%d,%d) = %g;\n", i+1,l+1,a[k]);*/
      }
      x0[i] = s0;
      x1[i] = s1;
    }
 
  return(AMG_OK);
}

int AMG_dmatminus_SCALAR2 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                           AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,l,start,end;
  register double *x0, *y0, *x1, *y1, *a;
  register int *ra, *ja;
  register double s0, s1 ;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x_)!=2*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatminus SCALAR2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=2*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatminus SCALAR2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatminus SCALAR2 - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  n = AMG_MATRIX_N(A);
  x0 = AMG_VECTOR_X(x_);
  x1 = x0 + n;
  y0 = AMG_VECTOR_X(y_);
  y1 = y0 + n;
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<n; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y0[i];
      s1 = a[start]*y1[i];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a[k]*y0[l];
        s1 += a[k]*y1[l];
      }
      x0[i] -= s0;
      x1[i] -= s1;
    }

  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR SCALAR SYSTEMS WITH THREE RIGHT HAND SIDES  */
/*                                                                          */
/****************************************************************************/

int AMG_dmatmul_SCALAR3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                         AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,l,start,end;
  register double *x0, *y0, *x1, *y1, *x2, *y2, *a;
  register int *ra, *ja;
  register double s0, s1, s2;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x_)!=3*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=3*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatmul SCALAR3 - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  n = AMG_MATRIX_N(A);
  x0 = AMG_VECTOR_X(x_);
  x1 = x0 + n;
  x2 = x1 + n;
  y0 = AMG_VECTOR_X(y_);
  y1 = y0 + n;
  y2 = y1 + n;
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<n; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y0[i];
      s1 = a[start]*y1[i];
      s2 = a[start]*y2[i];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a[k]*y0[l];
        s1 += a[k]*y1[l];
        s2 += a[k]*y2[l];
      }
      x0[i] = s0;
      x1[i] = s1;
      x2[i] = s2;
    }

  return(AMG_OK);
}

int AMG_dmatminus_SCALAR3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                           AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,l,start,end;
  register double *x0, *y0, *x1, *y1, *x2, *y2, *a;
  register int *ra, *ja;
  register double s0, s1, s2;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x_)!=3*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatminus SCALAR3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=3*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatminus SCALAR3 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatminus SCALAR3 - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  n = AMG_MATRIX_N(A);
  x0 = AMG_VECTOR_X(x_);
  x1 = x0 + n;
  x2 = x1 + n;
  y0 = AMG_VECTOR_X(y_);
  y1 = y0 + n;
  y2 = y1 + n;
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<n; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y0[i];
      s1 = a[start]*y1[i];
      s2 = a[start]*y2[i];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a[k]*y0[l];
        s1 += a[k]*y1[l];
        s2 += a[k]*y2[l];
      }
      x0[i] -= s0;
      x1[i] -= s1;
      x2[i] -= s2;
    }

  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR SCALAR SYSTEMS WITH SIX RIGHT HAND SIDES    */
/*                                                                          */
/****************************************************************************/

int AMG_dmatmul_SCALAR6 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                         AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,l,start,end;
  register double *x0, *y0, *x1, *y1, *x2, *y2, *a;
  register double *x3, *y3, *x4, *y4, *x5, *y5;
  register int *ra, *ja;
  register double s0, s1, s2, s3, s4, s5;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x_)!=6*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR6 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=6*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatmul SCALAR6 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatmul SCALAR6 - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  n = AMG_MATRIX_N(A);
  x0 = AMG_VECTOR_X(x_);
  x1 = x0 + n;
  x2 = x1 + n;
  x3 = x2 + n;
  x4 = x3 + n;
  x5 = x4 + n;
  y0 = AMG_VECTOR_X(y_);
  y1 = y0 + n;
  y2 = y1 + n;
  y3 = y2 + n;
  y4 = y3 + n;
  y5 = y4 + n;
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<n; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y0[i];
      s1 = a[start]*y1[i];
      s2 = a[start]*y2[i];
      s3 = a[start]*y3[i];
      s4 = a[start]*y4[i];
      s5 = a[start]*y5[i];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a[k]*y0[l];
        s1 += a[k]*y1[l];
        s2 += a[k]*y2[l];
        s3 += a[k]*y3[l];
        s4 += a[k]*y4[l];
        s5 += a[k]*y5[l];
      }
      x0[i] = s0;
      x1[i] = s1;
      x2[i] = s2;
      x3[i] = s3;
      x4[i] = s4;
      x5[i] = s5;
    }

  return(AMG_OK);
}

int AMG_dmatminus_SCALAR6 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                           AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,l,start,end;
  register double *x0, *y0, *x1, *y1, *x2, *y2, *a;
  register double *x3, *y3, *x4, *y4, *x5, *y5;
  register int *ra, *ja;
  register double s0, s1, s2, s3, s4, s5;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(x_)!=6*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatminus SCALAR6 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=6*AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmatminus SCALAR6 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmatminus SCALAR6 - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  n = AMG_MATRIX_N(A);
  x0 = AMG_VECTOR_X(x_);
  x1 = x0 + n;
  x2 = x1 + n;
  x3 = x2 + n;
  x4 = x3 + n;
  x5 = x4 + n;
  y0 = AMG_VECTOR_X(y_);
  y1 = y0 + n;
  y2 = y1 + n;
  y3 = y2 + n;
  y4 = y3 + n;
  y5 = y4 + n;
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<n; i++)
    {
      start = ra[i]; 
      end = start+ja[start];
      s0 = a[start]*y0[i];
      s1 = a[start]*y1[i];
      s2 = a[start]*y2[i];
      s3 = a[start]*y3[i];
      s4 = a[start]*y4[i];
      s5 = a[start]*y5[i];
      for (k=start+1; k<end; k++) 
      {
        l = ja[k];
        s0 += a[k]*y0[l];
        s1 += a[k]*y1[l];
        s2 += a[k]*y2[l];
        s3 += a[k]*y3[l];
        s4 += a[k]*y4[l];
        s5 += a[k]*y5[l];
      }
      x0[i] -= s0;
      x1[i] -= s1;
      x2[i] -= s2;
      x3[i] -= s3;
      x4[i] -= s4;
      x5[i] -= s5;
    }

  return(AMG_OK);
}


/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR 1D SADDLE POINT SYSTEMS                     */
/*                                                                          */
/* ( A  B) y = x                                                            */
/* (B^T 0)                                                                  */
/*                                                                          */
/****************************************************************************/

int AMG_dmatmul_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **List_B, AMG_VECTOR *y_)
{
  register int n_a,n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;
  register AMG_MATRIX *B=List_B[0];

  
  n_a = AMG_MATRIX_N(A);          /* columns in A */
  n_b = AMG_MATRIX_N(B);          /* columns in B */
  m_b = AMG_MATRIX_M(B);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);          /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);         /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul SADDLE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=n_a+n_b)
     {
        AMG_Print("ERROR: dmatmul SADDLE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (n_a!=m_b)
     {
        AMG_Print("ERROR: dmatmul SADDLE_1 - matrix sizes mismatched !!\n");
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
      s = a[start]*y[i];
      for (k=start+1; k<end; k++) s += a[k]*y[ja[k]];
      x[i] = s;
    }

  a = AMG_MATRIX_A(B);
  ra = AMG_MATRIX_RA(B);
  ja = AMG_MATRIX_JA(B);

 /* multiplication with block B */
  
  for (i=0; i<n_a; i++)  /* n_a == m_b */
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      /*s = a[start]*y[i+n_a];*/
      s = 0;
      for (k=start; k<=end; k++)         
        s += a[k]*y[ja[k]+n_a];        
      x[i] += s;
    }

 /* multiplication with block B^T */

  for (i=0; i<n_b; i++) x[i+n_a] = 0.0; 

  for (i=0; i<n_a; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      /*x[n_a+i] += a[start]*y[i];*/
      for (k=start; k<=end; k++)
        x[n_a+ja[k]] += a[k]*y[i];
    }

  return(AMG_OK);
}

int AMG_dmatminus_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **List_B, AMG_VECTOR *y_)
{
  register int n_a,n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;
  register AMG_MATRIX *B=List_B[0];
  
  n_a = AMG_MATRIX_N(A);          /* columns in A */
  n_b = AMG_MATRIX_N(B);          /* columns in B */
  m_b = AMG_MATRIX_M(B);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);          /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);         /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus SADDLE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=n_a+n_b)
     {
        printf("x %d a %d b %d\n",n_x,n_a,n_b);
        AMG_Print("ERROR: dmatminus SADDLE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (n_a!=m_b)
     {
        printf("y %d a %d b %d\n",n_x,n_a,m_b);
        AMG_Print("ERROR: dmatminus SADDLE_1 - matrix sizes mismatched !!\n");
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
      s = a[start]*y[i];
      for (k=start+1; k<end; k++) s += a[k]*y[ja[k]];
      x[i] -= s;
    }

  a = AMG_MATRIX_A(B);
  ra = AMG_MATRIX_RA(B);
  ja = AMG_MATRIX_JA(B);

 /* multiplication with block B */
  
  for (i=0; i<n_a; i++)  /* n_a == m_b */
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      /*s = a[start]*y[i+n_a];*/
      s = 0;
      for (k=start; k<=end; k++)         
        s += a[k]*y[ja[k]+n_a];        
      x[i] -= s;
    }

 /* multiplication with block B^T */

  for (i=0; i<n_a; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      /*x[n_a+i] += a[start]*y[i];*/
      for (k=start; k<=end; k++)
        x[n_a+ja[k]] -= a[k]*y[i];
    }

  return(AMG_OK);
}

int AMG_A_dmatmul_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **List_B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;
  
  n_a = AMG_MATRIX_N(A);          /* columns in A */
  n_x = AMG_VECTOR_N(x_);          /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);         /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=n_a)
     {
        AMG_Print("ERROR: A_dmatmul SADDLE_1 - vector and matrix lengths mismatched !!\n");
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
      s = a[start]*y[i];
      for (k=start+1; k<end; k++) s += a[k]*y[ja[k]];
      x[i] = s;
    }

  return(AMG_OK);
}
int AMG_B_dmatmul_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **List_B, AMG_VECTOR *y_)
{
  register int n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;
  register AMG_MATRIX *B=List_B[0];
  
  n_b = AMG_MATRIX_N(B);          /* columns in B */
  m_b = AMG_MATRIX_M(B);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);          /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);         /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=m_b)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (n_y!=n_b)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(B);
  ra = AMG_MATRIX_RA(B);
  ja = AMG_MATRIX_JA(B);

 /* multiplication with block B */
  
  for (i=0; i<m_b; i++)  /* n_a == m_b */
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      /*s = a[start]*y[i+n_a];*/
      s = 0;
      for (k=start; k<=end; k++)         
        s += a[k]*y[ja[k]];        
      x[i] = s;
    }

  return(AMG_OK);
}
int AMG_B_Trans_dmatmul_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **List_B, AMG_VECTOR *y_)
{
  register int n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s;
  register int b,bb;
  register AMG_MATRIX *B=List_B[0];
  
  n_b = AMG_MATRIX_N(B);          /* columns in B */
  m_b = AMG_MATRIX_M(B);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);          /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);         /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_b)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (n_y!=m_b)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(B);
  ra = AMG_MATRIX_RA(B);
  ja = AMG_MATRIX_JA(B);

 /* multiplication with block B^T */

  for (i=0; i<n_b; i++) x[i] = 0.0; 

  for (i=0; i<m_b; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      /*x[n_a+i] += a[start]*y[i];*/
      for (k=start; k<=end; k++)
        x[ja[k]] += a[k]*y[i];
    }

  return(AMG_OK);
}
/****************************************************************************/
/*                                                                          */
/* MATRIX-VECTOR OPERATIONS FOR 1D SADDLE POINT PROBLEMS OF TYPE 1          */
/*                                                                          */
/* ( A    B)                                                                */
/* (       ) y = x                                                          */
/* (B^T   0)                                                                */
/*                                                                          */
/* Input : A , B^T, x, y                                                    */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
int AMG_dmatmul_SADDLE_1_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  
  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        printf("n_x %d n_y %d\n",n_x,n_y); 
        AMG_Print("ERROR: dmatmul SADDLE_1_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=n_a+m_b)
     {
        AMG_Print("ERROR: dmatmul SADDLE_1_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (n_a!=n_b)
     {
        AMG_Print("ERROR: dmatmul SADDLE_1_TYPE_1 - matrix sizes mismatched !!\n");
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
      for (k=start+1; k<end; k++) 
        {
          s0 += a[k]*y[ja[k]];
        }
      x[i] = s0;
    }

 /* multiplication with block B */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        {
          if (aja[ara[ja[k]]] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
            x[ja[k]] += a[k]*y[n_a+i];
        }
    }

 /* multiplication with block B^T */
   
  for (i=0; i<m_b; i++)
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++) 
        s0 += a[k]*y[ja[k]];
      x[n_a+i] = s0;
    }

  return(AMG_OK);
}

int AMG_dmatminus_SADDLE_1_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus SADDLE_1_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (n_x!=n_a+m_b)
     {
        AMG_Print("ERROR: dmatminus SADDLE_1_TYPE_1 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (n_a!=n_b)
     {
        AMG_Print("ERROR: dmatminus SADDLE_1_TYPE_1 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      for (k=start+1; k<end; k++) 
        {
          s0 += a[k]*y[ja[k]];
        }
      x[i] -= s0;
    }

 /* multiplication with block B */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);

  for (i=0; i<m_b; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        if (aja[ara[ja[k]]] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
          x[ja[k]] -= a[k]*y[n_a+i];
    }

 /* multiplication with block B^T */
   
  for (i=0; i<m_b; i++)
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++) 
        s0 += a[k]*y[ja[k]];
      x[n_a+i] -= s0;
    }

  return(AMG_OK);
}

int AMG_B_dmatmul_SADDLE_1_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  
  n_b = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_b)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_1_TYPE_1 - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_y!=m_b)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_1_TYPE_1 - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  for (i=0; i<n_x; i++) x[i] = 0.0; 

 /* multiplication with block B */

  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);


  for (i=0; i<m_b; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        if (aja[ara[ja[k]]] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
          x[ja[k]] += a[k]*y[i];    /* skip it */
    }

  return(AMG_OK);
}

int AMG_B_Trans_dmatmul_SADDLE_1_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b,n_x,n_y,i,k,start,end,m_b;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s0,s1;
  register int b,bb;
  
  n_b = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_y!=n_b)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_1_TYPE_1 - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_x!=m_b)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_1_TYPE_1 - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);

 /* multiplication with block B^T */
  
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
 
  for (i=0; i<m_b; i++)
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)         
        s0 += a[k]*y[ja[k]];        
      x[i] = s0;
    }

  return(AMG_OK);
}



int AMG_dmatminus_gal_InterpolateVectorToLevel0(int k)
{
  int i;
  
  for (i=k-1;i>=0;i--)                /* prolongate velo to finest level */
    {
      AMG_dset(velo_prolong[i],0.0);  
      pc_prolongate_auto_2d(G[i],G[i],velo_prolong[i],velo_prolong[i+1],global_sc->omega_p); /* result velo_prolong[0] */
    }
  for (i=k-1;i>=0;i--)                /* prolongate pressure to finest level */
    {
      AMG_dset(pres_prolong[i],0.0);  
      pc_prolongate_auto(G_schur[i],G_schur[i],pres_prolong[i],pres_prolong[i+1],global_sc->omega_p); /* result pres_prolong[0] */
    }
  AMG_dset(velo_result,0.0);          /* initialize result vectors */
  return(0);
}

int AMG_dmatminus_gal_RestrictResultFromLevel0(int k)
{
  int i;
  
  AMG_dcopy(velo_prolong[0],velo_result);          /* initialize result vectors */

  for (i=0;i<k;i++)                /* restrict velo to level k */
    {
      pc_restrict_2d(G[i],G[i],velo_prolong[i],velo_prolong[i+1]);
    }
  for (i=0;i<k;i++)                /* restrict pres to level k */
    {
      pc_restrict(G_schur[i],G_schur[i],pres_prolong[i],pres_prolong[i+1]);
    }
  return(0);
}

int AMG_dmatminus_gal_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1,n_b,n_a0;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
  AMG_VECTOR inter_y;
  int level;
   
  level = A->level;
  if (level==0)                      /* finest level, use normal routine */                
    {
      /*AMG_dmatminus_SADDLE_2_TYPE_1 (x_,A,B,y_);
        return(0)*/;
    }

  n_a = AMG_MATRIX_N(A);              /* columns in A */
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B on finest level */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B on finest level */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B on finest level */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B on finest level */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */
  n_b = n_x - 2*n_a;                  /* length of pressure part in vector */
 
  if (global_sc->verbose>1)
  {
     if (m_b0!=m_b1)                     /* # rows in B matrices must be equal */      
     {
        AMG_Print("ERROR: dmatminus_inter SADDLE_2_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=n_y)                       /* length of solution and rhs vector must be equal */
     {
        AMG_Print("ERROR: dmatminus_inter SADDLE_2_TYPE_1 - vector lengths mismatched !!\n");
        exit(4711);
     }
  }
  for (i=0;i<2*n_a;i++)                 /* prepare velocity for prolongation on level 0 */
    {velo_prolong[level]->x[i] = y_->x[i];}
  for (i=0;i<n_b;i++)                 /* prepare pressure for prolongation on level 0 */
    {pres_prolong[level]->x[i] = y_->x[i+2*n_a];}

  AMG_dmatminus_gal_InterpolateVectorToLevel0(level);  /* interpolate rhs to level 0 in order to multiply with B, B^T */
  
  /* prepare data */
  x = AMG_VECTOR_X(x_);               /* solution */
  y = AMG_VECTOR_X(y_);               /* rhs */
  a = AMG_MATRIX_A(A);                /* matrix A */
  ra = AMG_MATRIX_RA(A);              /* rows in A */
  ja = AMG_MATRIX_JA(A);              /* columns in A */
  
  /* multiplication with block A, this is standard since A is known on level k */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          s0 += a[k]*y[ja[k]];
          s1 += a[k]*y[ja[k]+n_a];
          }
      x[i] -= s0;
      x[i+n_a] -= s1;
    }

 /* multiplication with block B_0 */

  a = AMG_MATRIX_A(B[2]);         /* matrix B[0] */
  ra = AMG_MATRIX_RA(B[2]);       /* rows in B[0] */ 
  ja = AMG_MATRIX_JA(B[2]);       /* columns in B[0] */

  for (i=0; i<m_b0; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        velo_result->x[ja[k]] += a[k]*pres_prolong[0]->x[i];
    }

 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);
  n_a0 =velo_prolong[0]->n/2; 

  for (i=0; i<m_b1; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
          velo_result->x[ja[k]+n_a0] += a[k]*pres_prolong[0]->x[i];
    }

 /* multiplication with block B_0^T */
  
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
  AMG_dset(pres_prolong[0],0.0);
 
  for (i=0; i<m_b0; i++)
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)         
        s0 += a[k]*velo_prolong[0]->x[ja[k]];        
      pres_prolong[0]->x[i] = s0;
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
        s0 += a[k]*velo_prolong[0]->x[ja[k]+n_a0];        
      pres_prolong[0]->x[i] += s0;
    }

  AMG_dmatminus_gal_RestrictResultFromLevel0(level);   /* restrict result from level 0 */

  for (i=0;i<2*n_a;i++)              /* update velocity solution */
    {x_->x[i] -= velo_prolong[level]->x[i];}
  for (i=0;i<n_x-2*n_a;i++)          /* update pressure solution */
    {x_->x[i+2*n_a] -= pres_prolong[level]->x[i];} 
 
  return(AMG_OK);
}

int AMG_B_dmatmul_InterpolateVectorToLevel0(int k)
{
  int i;
  
  for (i=k-1;i>=0;i--)                /* prolongate pressure to finest level */
    {
      AMG_dset(pres_prolong[i],0.0);  
      pc_prolongate_auto(G_schur[i],G_schur[i],pres_prolong[i],pres_prolong[i+1],global_sc->omega_p); /* result pres_prolong[0] */
    }
  AMG_dset(velo_prolong[0],0.0);          /* initialize result vectors */
  return(0);
}

int AMG_B_dmatmul_gal_RestrictResultFromLevel0(int k)
{
  int i;
  
  for (i=0;i<k;i++)                /* restrict velo to level k */
    {
      pc_restrict_2d(G[i],G[i],velo_prolong[i],velo_prolong[i+1]);
    }
return(0);
}

int AMG_B_dmatmul_gal_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_) 
{
  register int n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *aja , *ara;
  register double s0,s1;
  register int b,bb;
  int level;
  
  n_b0 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  level = A->level;
  if (level==0)                       /* finest level */                   
    {
      /* AMG_B_dmatmul_SADDLE_2_TYPE_1 (x_,A,B,y_);
         return(0)*/;
    } 
  if (global_sc->verbose>1)
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: B_dmatmul_inter SADDLE_2_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
  for (i=0;i<n_y;i++)               /* copy rhs for prolongation */
    pres_prolong[level]->x[i] = y_->x[i];
  AMG_B_dmatmul_InterpolateVectorToLevel0(level); /* prolongate rhs */

 /* multiplication with block B_0 */

  x = AMG_VECTOR_X(velo_prolong[0]);
  y = AMG_VECTOR_X(pres_prolong[0]);
  a = AMG_MATRIX_A(B[2]);
  ra = AMG_MATRIX_RA(B[2]);
  ja = AMG_MATRIX_JA(B[2]);
  ara = AMG_MATRIX_RA(A);
  aja = AMG_MATRIX_JA(A);
  for  (i=0; i<m_b0; i++)
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        /*printf("b00(%d,%d) = %g;\n",i+1,ja[k]+1,a[k])*/;
    }

  for (i=0; i<m_b0; i++)           /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        x[ja[k]] += a[k]*y[i];
    }

 /* multiplication with block B_1 */

  a = AMG_MATRIX_A(B[3]);
  ra = AMG_MATRIX_RA(B[3]);
  ja = AMG_MATRIX_JA(B[3]);
  for  (i=0; i<m_b0; i++)
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        /*printf("b11(%d,%d) = %g;\n",i+1,ja[k]+1,a[k])*/;
    }
  for (i=0; i<m_b1; i++)           /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        x[ja[k]+n_b0] += a[k]*y[i];
    }
  AMG_B_dmatmul_gal_RestrictResultFromLevel0(level);   /* restrict result from level 0 */
  for (i=0;i<n_x;i++)              /* update result */
    x_->x[i] = velo_prolong[level]->x[i];

  return(AMG_OK);
}

int AMG_B_Trans_dmatmul_gal_InterpolateVectorToLevel0(int k)
{
  int i;
  
  for (i=k-1;i>=0;i--)                /* prolongate velo to finest level */
    {
      AMG_dset(velo_prolong[i],0.0);  
      pc_prolongate_auto_2d(G[i],G[i],velo_prolong[i],velo_prolong[i+1],global_sc->omega_p); /* result velo_prolong[0] */
    }
  return(0);
}

int AMG_B_Trans_dmatmul_gal_RestrictResultFromLevel0(int k)
{
  int i;
  
  for (i=0;i<k;i++)                /* restrict pres to level k */
    {
      pc_restrict(G_schur[i],G_schur[i],pres_prolong[i],pres_prolong[i+1]);
    }
  return(0);
}

int AMG_B_Trans_dmatmul_gal_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register double s0,s1;
  register int b,bb;
  int level;
  
  level = A->level;
  if (level==0)                       /* finest level */            
    {
      AMG_B_Trans_dmatmul_SADDLE_2_TYPE_1 (x_,A,B,y_);
      return(0);
    }

  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_2_TYPE_1 - B matrix sizes mismatched !!\n");
        exit(4711);
     }

  for (i=0;i<n_y;i++)                 /* copy rhs for prolongation */
    velo_prolong[level]->x[i] = y_->x[i];
  AMG_B_Trans_dmatmul_gal_InterpolateVectorToLevel0(level); /* prolongate rhs */

  x = AMG_VECTOR_X(pres_prolong[0]);
  y = AMG_VECTOR_X(velo_prolong[0]);
  a = AMG_MATRIX_A(B[0]);
  ra = AMG_MATRIX_RA(B[0]);
  ja = AMG_MATRIX_JA(B[0]);
 
 /* multiplication with block B_0^T */

  for (i=0; i<m_b0; i++)
    {
      start = ra[i]; 
      end= ra[i+1]-1;
      s0 = 0;
      for (k=start; k<=end; k++)         
        s0 += a[k]*y[ja[k]];        
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
        s0 += a[k]*y[ja[k]+n_b0];        
      x[i] += s0;
    }

  AMG_B_Trans_dmatmul_gal_RestrictResultFromLevel0(level); /* restrict solution */
  for (i=0;i<n_x;i++)          /* update result */ 
    x_->x[i] = pres_prolong[level]->x[i]; 

  return(AMG_OK);
}
/****************************************************************************/
/*                                                                            */
/* MATRIX-VECTOR OPERATIONS FOR SCALAR SYSTEMS OF VASSILEVSKI/LAZAROV            */
/* PRECONDITIONER TYPE                                                            */
/*                                                                            */
/* ( A    0) + 1/delta         B B^T  x = y                                            */
/* ( 0    A)                                                                */
/*                                                                            */
/* b^T = (B_1^T,B_2^T)                                                            */
/*                                                                            */
/* Input : A , B_1^T, B_2^T, x, y                                            */
/*                                                                            */
/****************************************************************************/
int AMG_dmatmul_SCALAR_VAL_LAZ_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1;
  register double *x, *y, *a, *xx, *aa, *yy, *tmp;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1,delta;
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
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }

  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  tmp = malloc(m_b0*sizeof(double)); /* to store results of multiplication with B^T */

  delta = global_sc->vas_laz_delta;

  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          s0 += a[k]*y[ja[k]];
          s1 += a[k]*y[ja[k]+n_a];
        }
      x[i] = s0;
      x[i+n_a] = s1;
    }

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
        s0 += a[k]*y[ja[k]];
      tmp[i] = s0;
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
        s0 += a[k]*y[ja[k]+n_a];        
      tmp[i] += s0;
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
      for (k=start; k<=end; k++)
        if (aja[ara[ja[k]]] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
          x[ja[k]] += a[k]*tmp[i]/delta;
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
        if (aja[ara[ja[k]]] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
          x[ja[k]+n_a] += a[k]*tmp[i]/delta;
    }

  free(tmp);
  return(AMG_OK);
}

int AMG_dmatminus_SCALAR_VAL_LAZ_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_b0,n_x,n_y,i,k,start,end,m_b0,n_b1,m_b1;
  register double *x, *y, *a, *xx, *aa, *yy, *tmp;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1,delta;
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
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - vector lengths mismatched !!\n");
        exit(4711);
     }
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a)
     {
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (2*n_a!=n_b0+n_b1)
     {
        AMG_Print("ERROR: dmatmul SCALAR_VAL_LAZ_2 - matrix sizes mismatched !!\n");
        exit(4711);
     }
  }

  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  tmp = malloc(m_b0*sizeof(double)); /* to store results of multiplication with B^T */

  delta = global_sc->vas_laz_delta;

  /* multiplication with block A */
  
  for (i=0; i<n_a; i++)
    {
      start = ra[i]; end = start+ja[start];
      s0 = a[start]*y[i];
      s1 = a[start]*y[i+n_a];
      for (k=start+1; k<end; k++) 
        {
          s0 += a[k]*y[ja[k]];
          s1 += a[k]*y[ja[k]+n_a];
        }
      x[i] -= s0;
      x[i+n_a] -= s1;
    }

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
        s0 += a[k]*y[ja[k]];
      tmp[i] = s0;
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
        s0 += a[k]*y[ja[k]+n_a];        
      tmp[i] += s0;
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
      for (k=start; k<=end; k++)
        if (aja[ara[ja[k]]] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
          x[ja[k]] -= a[k]*tmp[i]/delta;
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
        if (aja[ara[ja[k]]] > 1)    /* row ja[k] in A has only one entry -> Diriclet */
          x[ja[k]+n_a] -= a[k]*tmp[i]/delta;
    }
  free(tmp);

  return(AMG_OK);
}

/****************************************************************************/
/*                                                                            */
/* MATRIX-VECTOR OPERATIONS FOR 2D SADDLE POINT PROBLEMS OF TYPE 2            */
/* WITH MORTAR                                                                    */
/*                                                                            */
/* ( A    0     B_0 B_4 0  )                                                      */
/* ( 0    A     B_1 0  B_4 ) x = y                                            */
/* (B_2^T B_3^T 0   0  0   )                                                     */
/* (B_4^T  0    0   0  0   )                                                    */
/* ( 0    B_4^T 0   0  0   )                                                    */
/*                                                                            */
/* Input : A , B_0, B_1, B_2^T, B_3^T, B_4^T, x, y                            */
/*                                                                            */
/*                                                                            */
/****************************************************************************/
int AMG_dmatmul_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l;
  register int n_b0,n_b1,n_b2,n_b3,n_b4,m_b0,m_b1,m_b2,m_b3,m_b4;
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
  n_b4 = AMG_MATRIX_N(B[4]);          /* columns in B */
  m_b4 = AMG_MATRIX_M(B[4]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatmul SADDLE_2_TYPE_2_MORTAR  - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatmul SADDLE_2_TYPE_2_MORTAR  - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2+2*m_b4)
     {
        AMG_Print("ERROR: dmatmul SADDLE_2_TYPE_2_MORTAR  - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatmul SADDLE_2_TYPE_2_MORTAR  - matrix sizes mismatched !!\n");
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

  /* multiplication with block B_4^T */

  a = AMG_MATRIX_A(B[4]);
  ra = AMG_MATRIX_RA(B[4]);
  ja = AMG_MATRIX_JA(B[4]);
  for (i=0; i<m_b4; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = s1 =0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[2*n_a+m_b2+i] = s0;     
      x[2*n_a+m_b2+m_b4+i] = s1;
    }

  /* multiplication with block B_4  */
  for (i=0; i<m_b4; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          x[l] += a[k]*y[2*n_a+m_b2+i];
          x[l+n_a] += a[k]*y[2*n_a+m_b2+m_b4+i];
        }
    }
  return(AMG_OK);
}

int AMG_dmatminus_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_a,n_x,n_y,i,k,start,end,l;
  register int n_b0,n_b1,n_b2,n_b3,n_b4,m_b0,m_b1,m_b2,m_b3,m_b4;
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
  n_b4 = AMG_MATRIX_N(B[4]);          /* columns in B */
  m_b4 = AMG_MATRIX_M(B[4]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (n_x!=n_y)
     {
        AMG_Print("ERROR: dmatminus SADDLE_2_TYPE_2_MORTAR  - vector lengths mismatched !!\n");
        exit(4711);
     }
     if ((m_b2!=m_b3)||(n_b2!=n_b3)|| (m_b0!=m_b1)||(n_b0!=n_b1)) 
     {
        AMG_Print("ERROR: dmatminus SADDLE_2_TYPE_2_MORTAR  - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_x!=2*n_a+m_b2+2*m_b4)
     {
        AMG_Print("ERROR: dmatminus SADDLE_2_TYPE_2_MORTAR  - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if ((2*n_a!=m_b0+m_b1)||(2*n_a!=n_b2+n_b3))
     {
        AMG_Print("ERROR: dmatminus SADDLE_2_TYPE_2_MORTAR  - matrix sizes mismatched !!\n");
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

   /* multiplication with block B_4^T */

  a = AMG_MATRIX_A(B[4]);
  ra = AMG_MATRIX_RA(B[4]);
  ja = AMG_MATRIX_JA(B[4]);
  for (i=0; i<m_b4; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = s1 =0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_a];
        }
      x[2*n_a+m_b2+i] = -s0;     
      x[2*n_a+m_b2+m_b4+i] = -s1;
    }

  /* multiplication with block B_4  */
  for (i=0; i<m_b4; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          x[l] -= a[k]*y[2*n_a+m_b2+i];
          x[l+n_a] -= a[k]*y[2*n_a+m_b2+m_b4+i];
        }
    }
 return(AMG_OK);
}

int AMG_B_dmatmul_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_x,n_y,i,k,start,end,l;
  register int n_b0,n_b1,n_b4,m_b0,m_b1,m_b4;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;
 
  n_b0 = AMG_MATRIX_N(B[0]);          /* columns in B */
  m_b0 = AMG_MATRIX_M(B[0]);          /* rows in B */
  n_b1 = AMG_MATRIX_N(B[1]);          /* columns in B */
  m_b1 = AMG_MATRIX_M(B[1]);          /* rows in B */
  n_b4 = AMG_MATRIX_N(B[4]);          /* columns in B */
  m_b4 = AMG_MATRIX_M(B[4]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (m_b0!=m_b1)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_2_TYPE_2_MORTAR  - B matrix sizes mismatched (1)!!\n");
        exit(4711);
     }
     
     if (m_b0!=n_b4)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_2_TYPE_2_MORTAR  - B matrix sizes mismatched (2)!!\n");       
        exit(4711);
     }
     if (n_b0!=n_b1)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_2_TYPE_2_MORTAR  - B matrix sizes mismatched (3)!!\n");
        exit(4711);
     }
     if (n_x!=m_b0+m_b1)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_2_TYPE_2_MORTAR  - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_y!=n_b0+2*m_b4)
     {
        AMG_Print("ERROR: B_dmatmul SADDLE_2_TYPE_2_MORTAR  - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }


  /* prepare data */
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

  /* multiplication with block B_4  */
  a = AMG_MATRIX_A(B[4]);
  ra = AMG_MATRIX_RA(B[4]);
  ja = AMG_MATRIX_JA(B[4]);
  for (i=0; i<m_b4; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          x[l] += a[k]*y[n_b0+i];
          x[l+m_b0] += a[k]*y[n_b0+m_b4+i];
        }
    }
  return(AMG_OK);
}
int AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n_x,n_y,i,k,start,end,l;
  register int n_b2,n_b3,n_b4,m_b2,m_b3,m_b4;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ara, *aja;
  register double s0,s1;
  register int b,bb;

  n_b2 = AMG_MATRIX_N(B[2]);          /* columns in B */
  m_b2 = AMG_MATRIX_M(B[2]);          /* rows in B */
  n_b3 = AMG_MATRIX_N(B[3]);          /* columns in B */
  m_b3 = AMG_MATRIX_M(B[3]);          /* rows in B */
  n_b4 = AMG_MATRIX_N(B[4]);          /* columns in B */
  m_b4 = AMG_MATRIX_M(B[4]);          /* rows in B */
  n_x = AMG_VECTOR_N(x_);             /* length of solution vector */
  n_y = AMG_VECTOR_N(y_);             /* length of rhs vector */

  if (global_sc->verbose>1)
  {
     if (m_b2!=m_b3)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_2_TYPE_2_MORTAR  - B matrix sizes mismatched !!\n");
        exit(4711);
     }
     if (n_y!=n_b2+n_b3)
     {
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_2_TYPE_2_MORTAR  - vector and matrix lengths mismatched(1) !!\n");
        exit(4711);
     }
     if (n_x!=m_b2+2*m_b4)
     {
        printf("n_x %d m_b2 %d n_b4 %d m_b4 %d\n",n_x,m_b2,n_b4,m_b4);
        AMG_Print("ERROR: B_Trans_dmatmul SADDLE_2_TYPE_2_MORTAR  - vector and matrix lengths mismatched(2) !!\n");
        exit(4711);
     }
  }
  /* prepare data */
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);

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
      x[i] = s0;     
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
          l = ja[k]+n_b2;
          s0 += a[k]*y[l];
        }
      x[i] += s0;
    }
  /* multiplication with block B_4^T */

  a = AMG_MATRIX_A(B[4]);
  ra = AMG_MATRIX_RA(B[4]);
  ja = AMG_MATRIX_JA(B[4]);
  for (i=0; i<m_b4; i++)     /* multiply columnwise */ 
    {
      start=ra[i]; 
      end = ra[i+1]-1;
      s0 = s1 =0;
      for (k=start; k<=end; k++)
        {
          l = ja[k];
          s0 += a[k]*y[l];
          s1 += a[k]*y[l+n_b2];
        }
      x[m_b2+i] = s0;     
      x[m_b2+m_b4+i] = s1;
    }

  return(AMG_OK);
}

int AMG_dmattransmul (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_)
{
  register int n,i,k,start,end,col;
  register double *x, *y, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int b,bb;
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_B(x_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
     if (AMG_VECTOR_B(y_)!=AMG_MATRIX_B(A)) return(AMG_FATAL);
     if (AMG_VECTOR_N(x_)!=AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmattransmul SCALAR - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(y_)!=AMG_MATRIX_N(A))
     {
        AMG_Print("ERROR: dmattransmul SCALAR - vector and matrix lengths mismatched !!\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(x_)!=AMG_VECTOR_N(y_))
     {
        AMG_Print("ERROR: dmattransmul SCALAR - vector lengths mismatched !!\n");
        exit(4711);
     }
  }

  /* prepare data */
  n = AMG_VECTOR_N(x_);
  b = AMG_VECTOR_B(x_);
  bb = AMG_MATRIX_BB(A);
  x = AMG_VECTOR_X(x_);
  y = AMG_VECTOR_X(y_);
  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  AMG_dset(x_,0.0);

  /* loop */
  switch (b)
    {
    case 1:
      for (i=0; i<n; i++)
        {
          start = ra[i]; 
	  end = start+ja[start];
	  x[i] += a[start]*y[i];         /* diagonal entry */
	  for (k=start+1; k<end; k++)
	  {  
	    col = ja[k];
	    x[col] += a[k] * y[i];
	  }
	}
      break;
      
      /*case 2:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO2(xx);
          yy = y+(i*b); XAY2(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAY2(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;                
      
    case 3:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO3(xx);
          yy = y+(i*b); XAY3(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAY3(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;                

    case 4:
      xx = x; aa = a;
      for (i=0; i<n; i++)
        {
          start = ra[i]; end = start+ja[start];
          ZERO4(xx);
          yy = y+(i*b); XAY4(xx,aa,yy); aa+=bb;
          for (k=start+1; k<end; k++)
            {
              yy = y+(ja[k]*b); XAY4(xx,aa,yy); aa+=bb;
            }
          xx+=b;
        }
      break;
      */
    default:
      AMG_Print("dmattransmul: blocksize>1 not implemented yet\n");
      break;        
    } 
  return(AMG_OK);
}
