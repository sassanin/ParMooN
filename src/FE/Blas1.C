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
   
// =======================================================================
// Blas1.C
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          17.01.2003
//
// =======================================================================

#include <string.h>
#include <math.h>

/** return inner product (x,y) */
double Ddot(int n, double *x, double *y)
{
  register double r;
  register int i;
  register double *a, *b;

  r = 0.0;
  a = x;
  b = y;
  for(i=0;i<n;i++)
  {
    r += *a * *b;
    a++;
    b++;
  }

  return r;
}

/** y := alpha*x + y */
void Daxpy(int n, double alpha, double *x, double *y)
{
  register int i;
  register double *a, *b;
  register double scal;

  a = x;
  b = y;
  scal = alpha;
  for(i=0;i<n;i++)
  {
    *b += scal * *a;
    a++;
    b++;
  }

}

/** z := alpha*x + beta*y */
void Dsum(int n, double alpha, double beta, double *x, double *y, double *z)
{
  register int i;
  register double *a, *b, *c;
  register double scal1, scal2;

  a = x;
  b = y;
  c = z;
  scal1 = alpha;
  scal2 = beta;
  for(i=0;i<n;i++)
  {
    *c = scal1 * *a + scal2 * *b;
    a++;
    b++;
    c++;
  }

}

/** b := a */
void Dcopy(int n, double *a, double *b)
{
  memcpy(b, a, n*sizeof(double));
}

/** x := alpha*x */
void Dscal(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a *= scal;
    a++;
  }
}

/** x := alpha+x */
void Dadd(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a += scal;
    a++;
  }
}


/** return Euclidian norm of x */
double Dnorm(int n, double *x)
{
  register double r;
  register int i;
  register double *a;
  register double z;

  a = x;
  r = 0.0;
  for(i=0; i<n; i++)
  {
    z = *a; 
    r += z*z;
    a++;
  }

  return sqrt(r);
}
