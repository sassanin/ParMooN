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
/* File:          amg_iter.c                                                */
/*                                                                          */
/* Purpose:   Simple iterative schemes for amg package                      */
/*                                                                          */
/* Author:          Peter Bastian                                           */
/*                          Institut fuer Computeranwendungen III           */
/*                          Universitaet Stuttgart                          */
/*                          Pfaffenwaldring 27                              */
/*                          70550 Stuttgart                                 */
/*                          email: peter@ica3.uni-stuttgart.de              */
/*                          phone: 0049-(0)711-685-7003                     */
/*                          fax  : 0049-(0)711-685-7000                     */
/*                                                                          */
/* Author:    Volker John                                                   */
/*            Otto--von--Guericke Universitaet Magdeburg                    */
/*            Institut fuer Analysis und Numerik                            */
/*            Postfach 4120                                                 */
/*            39016 Magdeburg                                               */
/*            email : volker.john@mathematik.uni-magdeburg.de               */
/*                                                                          */
/* History:   04 FEB 1996 Begin                                             */
/*                          01 OKT 1997 redesign                            */
/*                          21 OKT 1997 EX code due to Klaus Johannsen      */
/*              1998/02/24 ILU backward and forward                         */
/*              1998/06/14 2d backwards and forwards                        */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*                          system include files                            */
/*                          application include files                       */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <amg_header.h>
#include <amg_low.h>
#include <amg_sp.h>
#include <amg_blas.h>
#include <amg_solve_main.h>
#include <amg_iter.h>

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*                  in the corresponding include file!)                     */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/
static char buffer[128];
extern AMG_SolverContext *global_sc;
/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*D
   AMG_jac - jacobi kernel 
   
   SYNOPSIS:
   int AMG_jac (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d)

   PARAMETERS:
.  A - matrix
.  d - defect
.  v - correction
   
   DESCRIPTION:
   Solves the system Dv=d where D=diag(A).
   
   RETURN VALUE:
.n AMG_OK
.n AMG_FATAL 
   
D*/
/****************************************************************************/

int AMG_jac (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end;
  register double *v, *d, *a, *xx, *aa, *yy, diag;
  register int *ra;
  register int b,blocks;
  double om;
  char buf[60];
  
  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        AMG_Print("AMG_jac : vector and matrix size mismatched !\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_jac : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  
  /* loop */
  switch (blocks)
    {
    case 1 :
      om = omega[0];
      for (i=0; i<n; i++) 
      {
         diag = a[ra[i]];
         v[i] = om*d[i]/diag;
      /*printf("%d %g %g\n",i,d[i],a[ra[i]]);*/
      }
      break;
    case 2 :
      om = omega[0];
      for (i=0; i<n; i++) 
        {
          diag = a[ra[i]];
          v[i] = om*d[i]/diag;
          v[i+n] = om*d[i+n]/diag;
        }
      break;      
    case 3 :
      om = omega[0];
      for (i=0; i<n; i++) 
        {
          diag = a[ra[i]];
          v[i] = om*d[i]/diag;
          v[i+n] = om*d[i+n]/diag;
          v[i+2*n] = om*d[i+2*n]/diag;
        }
      break;      
    case 6 :
      om = omega[0];
      for (i=0; i<n; i++) 
        {
          diag = a[ra[i]];
          v[i] = om*d[i]/diag;
          v[i+n] = om*d[i+n]/diag;
          v[i+2*n] = om*d[i+2*n]/diag;
          v[i+3*n] = om*d[i+3*n]/diag;
          v[i+4*n] = om*d[i+4*n]/diag;
          v[i+5*n] = om*d[i+5*n]/diag;
        }
      break;      
    default:
      sprintf(buf,"jac: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }        
  return(AMG_OK);

}

/****************************************************************************/
/*D
   AMG_sorf - SOR kernels forward and backward
   
   SYNOPSIS:
   int AMG_sorf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
   int AMG_sorb (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)

   PARAMETERS:
.  A - matrix
.  d - defect
.  v - correction
.  omega - damping factor
   
   DESCRIPTION:
   Solve the system Lv=d or Uv=d with SOR damping where L is the
   lower triangle of A and U the upper triangle of A. The damping factor
   is an array containing one value per component.
   
   RETURN VALUE:
.n AMG_OK
.n AMG_FATAL 
   
D*/
/****************************************************************************/

int AMG_sorf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end,l;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int blocks;
  register double s,s1,s2,s3,s4,s5,om;
  char buf[60];

  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        sprintf(buf,"AMG_sorf : vector and matrix size mismatched (%d,%d) !\n",
                AMG_VECTOR_N(v_),blocks*n);
        AMG_Print(buf);
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_sorf : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (blocks)
  {
    case 1 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i) s += a[k]*v[l];
        }
        v[i] = om*(d[i]-s)/a[start];
      }
      break;
    case 2 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        }
      break;
    case 3 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
      }
      break;
      
    case 6 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = s3 = s4 = s5 = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
            s3 += a[k]*v[l+3*n];
            s4 += a[k]*v[l+4*n];
            s5 += a[k]*v[l+5*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
        v[i+3*n] = om*(d[i+3*n]-s3)/a[start];
        v[i+4*n] = om*(d[i+4*n]-s4)/a[start];
        v[i+5*n] = om*(d[i+5*n]-s5)/a[start];
     }
      break;
    default:
      sprintf(buf,"sorf: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
  return(AMG_OK);
}

int AMG_sorb (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end,l;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int blocks;
  register double s,s1,s2,s3,s4,s5,om;
  char buf[60];
  
  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        AMG_Print("AMG_sorb : vector and matrix size mismatched !\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_sorb : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (blocks)
  {
    case 1 :
      om=omega[0];
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l>i) s += a[k]*v[l];
        }
        v[i] = om*(d[i]-s)/a[start];
      }
      break;
    case 2 :
      om=omega[0];
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = 0.0;
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
          if (l>i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
      }
      break;
      
    case 3 :
      om=omega[0];
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = 0.0;
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
          if (l>i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
      }
      break;

    case 6 :
      om=omega[0];
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = s3 = s4 = s5 = 0.0;
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
          if (l>i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
            s3 += a[k]*v[l+3*n];
            s4 += a[k]*v[l+4*n];
            s5 += a[k]*v[l+5*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
        v[i+3*n] = om*(d[i+3*n]-s3)/a[start];
        v[i+4*n] = om*(d[i+4*n]-s4)/a[start];
        v[i+5*n] = om*(d[i+5*n]-s5)/a[start];
      }
      break;

    default:
      sprintf(buf,"sorb: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
   
  return(AMG_OK);
}
int AMG_sorfb (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end,l;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int blocks, n2,n3,n4,n5;
  register double s,s1,s2,s3,s4,s5,om;
  char buf[60];

  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        sprintf(buf,"AMG_sorfb : vector and matrix size mismatched (%d,%d) !\n",
                AMG_VECTOR_N(v_),blocks*n);
        AMG_Print(buf);
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_sorfb : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (blocks)
  {
    case 1 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s += a[k]*v[l];
        }
        v[i] += om*(d[i]-s)/a[start];
      }
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s += a[k]*v[l];
        }
        v[i] += om*(d[i]-s)/a[start];
      }
      break;
    case 2 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];          
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        }
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	}
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
      }
      break;
    case 3 :
      om=omega[0];
      n2 = 2*n;
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s += a[k]*v[l];
          s1 += a[k]*v[l+n];
          s2 += a[k]*v[l+n2];
	}
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
      }
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	  s2 += a[k]*v[l+n2];
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
      }
      break;
      
    case 6 :
      om=omega[0];
      n2 = 2*n;
      n3 = 3*n;
      n4 = 4*n;
      n5 = 5*n;
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        s3 = a[start]*v[i+n3];
	s4 = a[start]*v[i+n4];
	s5 = a[start]*v[i+n5];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	  s2 += a[k]*v[l+n2];
	  s3 += a[k]*v[l+n3];
	  s4 += a[k]*v[l+n4];
	  s5 += a[k]*v[l+n5];
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
        v[i+n3] += om*(d[i+n3]-s3)/a[start];
        v[i+n4] += om*(d[i+n4]-s4)/a[start];
        v[i+n5] += om*(d[i+n5]-s5)/a[start];
      }
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        s3 = a[start]*v[i+n3];
	s4 = a[start]*v[i+n4];
	s5 = a[start]*v[i+n5];
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	  s2 += a[k]*v[l+n2];
	  s3 += a[k]*v[l+n3];
	  s4 += a[k]*v[l+n4];
	  s5 += a[k]*v[l+n5];
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
        v[i+n3] += om*(d[i+n3]-s3)/a[start];
        v[i+n4] += om*(d[i+n4]-s4)/a[start];
        v[i+n5] += om*(d[i+n5]-s5)/a[start];
      }
      break;
    default:
      sprintf(buf,"sorf: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_iluf - ILU kernels forward and backward
   
   SYNOPSIS:
   int AMG_iluf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)
   int AMG_ilub (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)

   PARAMETERS:
.  A - matrix
.  d - defect
.  v - correction
.  omega - damping factor
   
   DESCRIPTION:
   Solve the system Lv=d or Uv=d with ILU  where L is the
   lower triangle of A and U the upper triangle of A.
   
   RETURN VALUE:
.n AMG_OK
.n AMG_FATAL 
   
D*/
/****************************************************************************/

int AMG_iluf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)
{
  register int n,i,k,start,end,l;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int blocks;
  register double s,s1,s2,s3,s4,s5;
  char buf[60];
  
  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        sprintf(buf,"AMG_iluf : vector and matrix size mismatched (%d,%d) !\n",
                AMG_VECTOR_N(v_),blocks*n);
        AMG_Print(buf);
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_iluf : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (blocks)
  {
    case 1 :
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];       /* row i */
        s = 0.0;
        for (k=start+1; k<end; k++)                 /* columns in row i */
        {
          l = ja[k];
          if (l<i) 
            s += a[k]*d[l];                     /* if column index lower than i, accumulate s */
        }
        v[i] = d[i]-s;                              /* diagonal is one */
      }
      break;
    case 2 :
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];       /* row i */
        s = s1 = 0.0;
        for (k=start+1; k<end; k++)                 /* columns in row i */
        {
          l = ja[k];
          if (l<i)                              /* if column index lower than i, accumulate s */
          {
            s += a[k]*d[l];
            s1 += a[k]*d[l+n];
          }      
        }   
        v[i] = d[i]-s;                              /* diagonal is one */
        v[i+n] = d[i+n]-s1;                         /* diagonal is one */
      }
      break;
    case 3 :
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];       /* row i */
        s = s1 = s2 = 0.0;
        for (k=start+1; k<end; k++)                 /* columns in row i */
        {
          l = ja[k];
          if (l<i)                              /* if column index lower than i, accumulate s */
          {
            s += a[k]*d[l];
            s1 += a[k]*d[l+n];
            s2 += a[k]*d[l+2*n];
          }      
        }   
        v[i] = d[i]-s;                              /* diagonal is one */
        v[i+n] = d[i+n]-s1;                         /* diagonal is one */
        v[i+2*n] = d[i+2*n]-s2;                         /* diagonal is one */
      }
      break;

    case 6 :
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];       /* row i */
        s = s1 = s2 = s3 = s4 = s5 = 0.0;
        for (k=start+1; k<end; k++)                 /* columns in row i */
        {
          l = ja[k];
          if (l<i)                              /* if column index lower than i, accumulate s */
          {
            s += a[k]*d[l];
            s1 += a[k]*d[l+n];
            s2 += a[k]*d[l+2*n];
            s3 += a[k]*d[l+3*n];
            s4 += a[k]*d[l+4*n];
            s5 += a[k]*d[l+5*n];
         }      
        }   
        v[i] = d[i]-s;                              /* diagonal is one */
        v[i+n] = d[i+n]-s1;                         /* diagonal is one */
        v[i+2*n] = d[i+2*n]-s2;                         /* diagonal is one */
        v[i+3*n] = d[i+3*n]-s3;                              /* diagonal is one */
        v[i+4*n] = d[i+4*n]-s4;                         /* diagonal is one */
        v[i+5*n] = d[i+5*n]-s5;                         /* diagonal is one */
      }
      break;
      
    default:
      sprintf(buf,"iluf: blocksize %d not implemented yet\n", blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
  
  
  return(AMG_OK);
}

int AMG_ilub (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)
{
  register int n,i,k,start,end,l;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int blocks;
  register double s,s1,s2,s3,s4,s5;
  char buf[60];
  
  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        AMG_Print("AMG_iluf : vector and matrix size mismatched !\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_iluf : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (blocks)
  {
    case 1 :
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];     /* row i */
        s = 0.0;
        for (k=start+1; k<end; k++)               /* columns of row i */
        {
          l = ja[k];
          if (l>i) s += a[k]*d[l];        /* if column index greater than i, accumulate s */
        }
        v[i] = (d[i]-s)/a[start];                 /* divide by the diagonal */
      }
      break;
    case 2 :
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];     /* row i */
        s = s1 = 0.0;
        for (k=start+1; k<end; k++)               /* columns of row i */
        {
          l = ja[k];
          if (l>i)                            /* if column index greater than i, accumulate s,s1 */
          {
            s += a[k]*d[l];        
            s1 += a[k]*d[l+n]; 
          }      
        } 
        v[i] = (d[i]-s)/a[start];                 /* divide by the diagonal */
        v[i+n] = (d[i+n]-s1)/a[start];            /* divide by the diagonal */
      }
      break;
      
    case 3 :
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];     /* row i */
        s = s1 = s2 = 0.0;
        for (k=start+1; k<end; k++)               /* columns of row i */
        {
          l = ja[k];
          if (l>i)                            /* if column index greater than i, accumulate s,s1 */
          {
            s += a[k]*d[l];        
            s1 += a[k]*d[l+n]; 
            s2 += a[k]*d[l+2*n]; 
          }      
        } 
        v[i] = (d[i]-s)/a[start];                 /* divide by the diagonal */
        v[i+n] = (d[i+n]-s1)/a[start];            /* divide by the diagonal */
        v[i+2*n] = (d[i+2*n]-s2)/a[start];            /* divide by the diagonal */
      }
      break;

    case 6 :
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];     /* row i */
        s = s1 = s2 = s3 = s4 = s5 = 0.0;
        for (k=start+1; k<end; k++)               /* columns of row i */
        {
          l = ja[k];
          if (l>i)                            /* if column index greater than i, accumulate s,s1 */
          {
            s += a[k]*d[l];        
            s1 += a[k]*d[l+n]; 
            s2 += a[k]*d[l+2*n]; 
            s3 += a[k]*d[l+3*n];        
            s4 += a[k]*d[l+4*n]; 
            s5 += a[k]*d[l+5*n]; 
          }      
        } 
        v[i] = (d[i]-s)/a[start];                 /* divide by the diagonal */
        v[i+n] = (d[i+n]-s1)/a[start];            /* divide by the diagonal */
        v[i+2*n] = (d[i+2*n]-s2)/a[start];            /* divide by the diagonal */
        v[i+3*n] = (d[i+3*n]-s3)/a[start];                 /* divide by the diagonal */
        v[i+4*n] = (d[i+4*n]-s4)/a[start];            /* divide by the diagonal */
        v[i+5*n] = (d[i+5*n]-s5)/a[start];            /* divide by the diagonal */
      }
      break;
      
    default:
      sprintf(buf,"ilub: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
        
  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_EXDecomposeMatrixdouble - LU decompose a band matrix (double numbers)

   SYNOPSIS:
   int AMG_EXDecomposeMatrixdouble (double *Mat, int bw, int n);

   PARAMETERS:
.  Mat - pointer to double array containing the bandmatrix
.  bw - bandwidth
.  n - number of rows (==columns) of the matrix

   DESCRIPTION:
   This function calculates the Gauss decomposition of the given band matrix; 
   the L and U factors are stored instead of the band matrix.

   RETURN VALUE:
   int  0: o.k.
        1: main diagonal element to small (0.0); can not divide
                
   SEE ALSO:
   EXDecomposeMatrixFLOAT, EXCopyMatrixdouble, EXApplyLUdouble
D*/                                                                        
/****************************************************************************/

int AMG_EXDecomposeMatrix (double *Mat, int bw, int n)
{
  int i,j,k;
  double f,d;
  
  for (i=0; i<n-1; i++)
    {
      d = AMG_EX_MAT(Mat,bw,i,i);
      if (AMG_ABS(d)<=1.0E-80) return (1);
      for (j=i+1; j<=AMG_MIN(i+bw,n-1); j++)
        {
          f = AMG_EX_MAT(Mat,bw,j,i)/d;
          AMG_EX_MAT(Mat,bw,j,i) = f;
          for (k=i+1; k<=AMG_MIN(i+bw,n-1); k++)
            AMG_EX_MAT(Mat,bw,j,k) -= f*AMG_EX_MAT(Mat,bw,i,k);
        }
    }
  return(0);
}

/****************************************************************************/
/*D
   AMG_EXApplyLUdouble - applies a LU decomposed band matrix (double numbers)

   SYNOPSIS:
   int AMG_EXApplyLUdouble (double *Mat, int bw, int n, double *Vec);

   PARAMETERS:
.  Mat - pointer to double array containing the bandmatrix
.  bw - bandwidth
.  n - number of rows (==columns) of the matrix
.  Vec - pointer to double array containing the vector

   DESCRIPTION:
   This function solves for the LU decomposed band matrix 'Mat' the equation
   L*U x = f. 
   f is provided in 'Vec' and the result x is returned again in 'Vec'.
   
   RETURN VALUE:
   int  0: o.k.
                
   SEE ALSO:
   EXApplyLUFLOAT, EXCopyMatrixdouble, EXDecomposeMatrixdouble
D*/                                                                        
/****************************************************************************/

int AMG_EXApplyLU (double *Mat, int bw, int n, int blocks, double *Vec)
{
  int i,j;
  double s;
    
  switch(blocks)
    {
    case 1 :
      /* invert lower */
      for (i=1; i<n; i++)
        for (j=AMG_MAX(i-bw,0); j<i; j++)
          Vec[i] -= AMG_EX_MAT(Mat,bw,i,j)*Vec[j];
      
      /* invert upper */
      for (i=n-1; i>=0; i--)
        {
          for (j=i+1; j<=AMG_MIN(i+bw,n-1); j++)
            Vec[i] -= AMG_EX_MAT(Mat,bw,i,j)*Vec[j];
          Vec[i] /= AMG_EX_MAT(Mat,bw,i,i);
        }
      break;
    case 2 :
      /* invert lower */
       for (i=1; i<n; i++)
        for (j=AMG_MAX(i-bw,0); j<i; j++)
          {
            s = AMG_EX_MAT(Mat,bw,i,j);
            Vec[i] -= s*Vec[j];
            Vec[i+n] -= s*Vec[j+n];
          }
      
      /* invert upper */
      for (i=n-1; i>=0; i--)
        {
          for (j=i+1; j<=AMG_MIN(i+bw,n-1); j++)
            {
              s = AMG_EX_MAT(Mat,bw,i,j);
              Vec[i] -= s*Vec[j];
              Vec[i+n] -= s*Vec[j+n];
            }
          s = AMG_EX_MAT(Mat,bw,i,i);
          Vec[i] /= s;
          Vec[i+n] /= s;
        }
      break;
    default:
      AMG_Print("EXApplyLU: blocksize>1 not implemented yet\n");
      exit(4711);
      break;
    }        
  return (0);
}

/****************************************************************************/
/*                                                                          */
/*                                                                          */
/* PRECONDITIONERS FOR TRANSPOSED MATRIX                                    */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
/*D
   AMG_sorf - SOR kernels forward and backward
   
   SYNOPSIS:
   int AMG_sorf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
   int AMG_sorb (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)

   PARAMETERS:
.  A - matrix
.  d - defect
.  v - correction
.  omega - damping factor
   
   DESCRIPTION:
   Solve the system Lv=d or Uv=d with SOR damping where L is the
   lower triangle of A and U the upper triangle of A. The damping factor
   is an array containing one value per component.
   
   RETURN VALUE:
.n AMG_OK
.n AMG_FATAL 
   
D*/
/****************************************************************************/

int AMG_sorf_trans(AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end,l,j;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ratr, *jatr, *postr;
  register int blocks;
  register double s,s1,s2,s3,s4,s5,om;
  char buf[60];

  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);
    
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        sprintf(buf,"AMG_sorf_trans : vector and matrix size mismatched (%d,%d) !\n",
                AMG_VECTOR_N(v_),blocks*n);
        AMG_Print(buf);
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_sorf_trans : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }
  
  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  ratr = A->ratr;
  jatr = A->jatr;
  postr = A->postr;
  
  /* loop */
  switch (blocks)
  {
    case 0 :
      om=omega[0];
      for (i=n-1; i>=0; i--)       /* component which will be updated */
      {
	s = 0;  
	for (j=i+1;j<n;j++)        /* look for components which are already updated */
	{
	  start = ra[j]; 
	  end = start+ja[start];
	  for (k=start+1; k<end; k++) 
	  {
	      l = ja[k];
	      if (l==i) 
	      {
		  /*  printf("%d %d %g %g\n",i,j,a[k],v[j]); */ 
		s += a[k]*v[j];
		break;
	      }
	  }
	}
	v[i] = om*(d[i]-s)/a[ra[i]];
      }
      break;
    case 1 :
      om=omega[0];
      for (i=n-1; i>=0; i--)       /* component which will be updated */
      {                            /* i: row number of transposed matrix */
	s = 0;
	start = ratr[i];
	end = start + jatr[start];
	for (k=start+1; k<end; k++) 
	{
	    l = jatr[k];           /* l: column number of transposed matrix */
	    if (i<l)
	    {
	      j = postr[k];
	      /* printf("k %d l %d j %d\n",k,l,j);*/
	      s += a[j]*v[l];
	    }
	}
	v[i] = om*(d[i]-s)/a[ra[i]];
      }
      break;
      /*  case 2 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        }
      break;
    case 3 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
      }
      break;
      
    case 6 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = s3 = s4 = s5 = 0.0;
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          if (l<i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
            s3 += a[k]*v[l+3*n];
            s4 += a[k]*v[l+4*n];
            s5 += a[k]*v[l+5*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
        v[i+3*n] = om*(d[i+3*n]-s3)/a[start];
        v[i+4*n] = om*(d[i+4*n]-s4)/a[start];
        v[i+5*n] = om*(d[i+5*n]-s5)/a[start];
     }
      break;
      */
    default:
      sprintf(buf,"sorf_trans: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
  return(AMG_OK);
}

int AMG_sorb_trans (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end,l,j;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int blocks;
  register double s,s1,s2,s3,s4,s5,om;
  char buf[60];
  
  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);
  
  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        AMG_Print("AMG_sorb : vector and matrix size mismatched !\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_sorb : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (blocks)
  {
    case 1 :
      om=omega[0];
      for (i=0; i<n; i++)       /* component which will be updated */
      {
	s = 0;  
	for (j=0;j<i;j++)        /* look for components which are already updated */
	{
	  start = ra[j]; 
	  end = start+ja[start];
	  for (k=start+1; k<end; k++) 
	  {
	      l = ja[k];
	      if (l==i) 
	      {
		  /*  printf("%d %d %g %g\n",i,j,a[k],v[j]); */ 
		s += a[k]*v[j];
		break;
	      }
	  }
	}
	v[i] = om*(d[i]-s)/a[ra[i]];
      }
      break;
      /*  case 2 :
      om=omega[0];
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = 0.0;
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
          if (l>i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
      }
      break;
      
    case 3 :
      om=omega[0];
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = 0.0;
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
          if (l>i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
      }
      break;

    case 6 :
      om=omega[0];
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = s1 = s2 = s3 = s4 = s5 = 0.0;
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
          if (l>i) 
          {
            s += a[k]*v[l];
            s1 += a[k]*v[l+n];
            s2 += a[k]*v[l+2*n];
            s3 += a[k]*v[l+3*n];
            s4 += a[k]*v[l+4*n];
            s5 += a[k]*v[l+5*n];
          }
        }
        v[i] = om*(d[i]-s)/a[start];
        v[i+n] = om*(d[i+n]-s1)/a[start];
        v[i+2*n] = om*(d[i+2*n]-s2)/a[start];
        v[i+3*n] = om*(d[i+3*n]-s3)/a[start];
        v[i+4*n] = om*(d[i+4*n]-s4)/a[start];
        v[i+5*n] = om*(d[i+5*n]-s5)/a[start];
      }
      break;
      */
    default:
      sprintf(buf,"sorb_trans: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
   
  return(AMG_OK);
}
int AMG_sorfb_trans (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega)
{
  register int n,i,k,start,end,l,j;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja;
  register int blocks, n2,n3,n4,n5;
  register double s,s1,s2,s3,s4,s5,om;
  char buf[60];

  printf("ssor for transposed matrix not implemented !\n");
  exit(4711);

  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        sprintf(buf,"AMG_sorfb_trans : vector and matrix size mismatched (%d,%d) !\n",
                AMG_VECTOR_N(v_),blocks*n);
        AMG_Print(buf);
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_sorfb_trans : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  /* loop */
  switch (blocks)
  {
    case 1 :
      om=omega[0];
      for (i=0; i<n; i++)       /* component which will be updated */
      {
	s = 0;  
	for (j=0;j<i;j++)        /* look for components which are already updated */
	{
	  start = ra[j]; 
	  end = start+ja[start];
	  for (k=start+1; k<end; k++) 
	  {
	      l = ja[k];
	      if (l==i) 
	      {
		  /*  printf("%d %d %g %g\n",i,j,a[k],v[j]); */ 
		s += a[k]*v[j];
		break;
	      }
	  }
	}
	v[i] += om*(d[i]-s)/a[ra[i]];
      }
      for (i=n-1; i>=0; i--)       /* component which will be updated */
      {
	s = 0;  
	for (j=i+1;j<n;j++)        /* look for components which are already updated */
	{
	  start = ra[j]; 
	  end = start+ja[start];
	  for (k=start+1; k<end; k++) 
	  {
	      l = ja[k];
	      if (l==i) 
	      {
		  /*  printf("%d %d %g %g\n",i,j,a[k],v[j]); */ 
		s += a[k]*v[j];
		break;
	      }
	  }
	}
	v[i] += om*(d[i]-s)/a[ra[i]];
      }
      break;
      /*
    case 2 :
      om=omega[0];
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];          
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        }
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	}
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
      }
      break;
    case 3 :
      om=omega[0];
      n2 = 2*n;
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
          s += a[k]*v[l];
          s1 += a[k]*v[l+n];
          s2 += a[k]*v[l+n2];
	}
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
      }
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	  s2 += a[k]*v[l+n2];
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
      }
      break;
      
    case 6 :
      om=omega[0];
      n2 = 2*n;
      n3 = 3*n;
      n4 = 4*n;
      n5 = 5*n;
      for (i=0; i<n; i++)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        s3 = a[start]*v[i+n3];
	s4 = a[start]*v[i+n4];
	s5 = a[start]*v[i+n5];
        for (k=start+1; k<end; k++) 
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	  s2 += a[k]*v[l+n2];
	  s3 += a[k]*v[l+n3];
	  s4 += a[k]*v[l+n4];
	  s5 += a[k]*v[l+n5];
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
        v[i+n3] += om*(d[i+n3]-s3)/a[start];
        v[i+n4] += om*(d[i+n4]-s4)/a[start];
        v[i+n5] += om*(d[i+n5]-s5)/a[start];
      }
      for (i=n-1; i>=0; i--)
      {
        start = ra[i]; end = start+ja[start];
        s = a[start]*v[i];
	s1 = a[start]*v[i+n];
	s2 = a[start]*v[i+n2];
        s3 = a[start]*v[i+n3];
	s4 = a[start]*v[i+n4];
	s5 = a[start]*v[i+n5];
        for (k=start+1; k<end; k++)
        {
          l = ja[k];
	  s += a[k]*v[l];
	  s1 += a[k]*v[l+n];
	  s2 += a[k]*v[l+n2];
	  s3 += a[k]*v[l+n3];
	  s4 += a[k]*v[l+n4];
	  s5 += a[k]*v[l+n5];
        }
        v[i] += om*(d[i]-s)/a[start];
        v[i+n] += om*(d[i+n]-s1)/a[start];
        v[i+n2] += om*(d[i+n2]-s2)/a[start];
        v[i+n3] += om*(d[i+n3]-s3)/a[start];
        v[i+n4] += om*(d[i+n4]-s4)/a[start];
        v[i+n5] += om*(d[i+n5]-s5)/a[start];
      }
      break;
      */
    default:
      sprintf(buf,"sorfb_trans: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_iluf - ILU kernels forward and backward
   
   SYNOPSIS:
   int AMG_iluf (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)
   int AMG_ilub (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)

   PARAMETERS:
.  A - matrix
.  d - defect
.  v - correction
.  omega - damping factor
   
   DESCRIPTION:
   Solve the system Lv=d or Uv=d with ILU  where L is the
   lower triangle of A and U the upper triangle of A.
   
   RETURN VALUE:
.n AMG_OK
.n AMG_FATAL 
   
D*/
/****************************************************************************/

int AMG_iluf_trans (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)
{
  register int n,i,k,start,end,l,j;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ratr, *jatr, *postr;
  register int blocks;
  register double s,s1,s2,s3,s4,s5;
  char buf[60];
  
  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        sprintf(buf,"AMG_iluf : vector and matrix size mismatched (%d,%d) !\n",
                AMG_VECTOR_N(v_),blocks*n);
        AMG_Print(buf);
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_iluf : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }
  
  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  ratr = A->ratr;
  jatr = A->jatr;
  postr = A->postr;
  
  /* loop */
  switch (blocks)
  {
    case 0 :
      for (i=n-1; i>=0; i--)       /* component which will be updated */
      {
	s = 0;  
	for (j=i+1;j<n;j++)        /* look for components which are already updated */
	{
	  start = ra[j]; 
	  end = start+ja[start];
	  for (k=start+1; k<end; k++) 
	  {
	      l = ja[k];
	      if (l==i) 
	      {
		  /*  printf("%d %d %g %g\n",i,j,a[k],v[j]); */ 
		s += a[k]*v[j];
		break;
	      }
	  }
	}
	v[i] = d[i]-s;
      }
      break;
    case 1 :
      for (i=n-1; i>=0; i--)       /* component which will be updated */
      {                            /* i: row number of transposed matrix */
	s = 0;
	start = ratr[i];
	end = start + jatr[start];
	for (k=start+1; k<end; k++) 
	{
	    l = jatr[k];           /* l: column number of transposed matrix */
	    if (i<l)
	    {
	      j = postr[k];
	      /* printf("k %d l %d j %d\n",k,l,j);*/
	      s += a[j]*v[l];
	    }
	}
	v[i] = d[i]-s;
      }
      break;
    default:
      sprintf(buf,"iluf_trans: blocksize %d not implemented yet\n", blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
  
  
  return(AMG_OK);
}

int AMG_ilub_trans (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_)
{
  register int n,i,k,start,end,l,j;
  register double *v, *d, *a, *xx, *aa, *yy;
  register int *ra, *ja, *ratr, *jatr, *postr;
  register int blocks;
  register double s,s1,s2,s3,s4,s5;
  char buf[60];
  
  /* prepare data */
  n  = AMG_MATRIX_N(A);
  blocks  = AMG_MATRIX_BLDI(A);
  v  = AMG_VECTOR_X(v_);
  d  = AMG_VECTOR_X(d_);

  /* plausi */
  if (global_sc->verbose>1)
  {
     if (AMG_VECTOR_N(v_)!= blocks*n) 
     {
        AMG_Print("AMG_ilub_trans : vector and matrix size mismatched !\n");
        exit(4711);
     }
     if (AMG_VECTOR_N(d_)!=blocks*n) 
     {
        AMG_Print("AMG_ilub_trans : vector and matrix size mismatched !\n");
        exit(4711);
     }
  }

  /* prepare data */
  a  = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  ratr = A->ratr;
  jatr = A->jatr;
  postr = A->postr;
    
  /* loop */
  switch (blocks)
  {
    case 0 :
      for (i=0; i<n; i++)       /* component which will be updated */
      {
	s = 0;  
	for (j=0;j<i;j++)        /* look for components which are already updated */
	{
	  start = ra[j]; 
	  end = start+ja[start];
	  for (k=start+1; k<end; k++) 
	  {
	      l = ja[k];
	      if (l==i) 
	      {
		  /*  printf("%d %d %g %g\n",i,j,a[k],v[j]); */ 
		s += a[k]*v[j];
		break;
	      }
	  }
	}
	v[i] = (d[i]-s)/a[ra[i]];
      }
      break;
    case 1 :
      for (i=0; i<n; i++)       /* component which will be updated */
      {                         /* i: row number of transposed matrix */
	s = 0;  
	start = ratr[i]; 
	end = start+jatr[start];
	for (k=start+1; k<end; k++) 
	{
	  l = jatr[k];          /* l: column number of transposed matrix */
	  if (i>l) 
	  {
	    j = postr[k];
	    s += a[j]*v[l];
	  }	
	}
	v[i] = (d[i]-s)/a[ra[i]];
      }
      break;       
    default:
      sprintf(buf,"ilub_trans: blocksize %d not implemented yet\n",blocks);
      AMG_Print(buf);
      exit(4711);
      break;        
    }
        
  return(AMG_OK);
}
