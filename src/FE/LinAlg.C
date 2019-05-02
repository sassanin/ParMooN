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
// @(#)LinAlg.C        1.18 07/03/00
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          27.01.1999
//              Volker John             27.10.1999  
//              Sashikumaar Ganesan     08.10.2009 (eigen values)
// =======================================================================

#include <Database.h>
#include <MooNMD_Io.h>
#ifdef __2D__
  #include <FEDatabase2D.h>
#else
  #include <FEDatabase3D.h>
#endif  
#include <LinAlg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

extern "C" {

void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);

void dgetrs_(char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv,
             double *B, int *ldb, int *info);

void dsptrd_(char *UPLO, int *n,  double *AP, double *D, double *E, 
             double *TAU, int *info);

void dopgtr_(char *UPLO, int *n,  double *AP, double *TAU, double *Q, int *LDQ, double *work,int *info);

void  dsteqr_(char *compz, int *N, double *D,  double *E, double *Z, int *LDZ, double *work, int *info);
 
void dgetri_(int *m, double *A, int *lda, int *ipiv, double *WORK, int *LWORK, int *info);

void dgemm_(char *TRANSA, char *TRANSB, int *m, int *n, int *k, double *ALPHA, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

}

// ==========================================================================
// general routines independent of dimension
// ==========================================================================
#define AT(i,j) (a[j*LDA+i])
#define A(i,j) (a[i*LDA+j])

void MatVectFull(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  int dof =  TDatabase::ParamDB->INTERNAL_LOCAL_DOF; 
  double *a = (double *)A[0];
  int i,j;

  memset(y,0,dof*SizeOfDouble);
  // matrix is stored as transposed;
  for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
      y[j] += a[i*dof+j] * x[i];


  return;

}
 
void DefectFull(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int dof =  TDatabase::ParamDB->INTERNAL_LOCAL_DOF; 
  double *a = (double *)A[0];
  int i,j;

  memcpy(r,b,dof*SizeOfDouble);
  // matrix is stored as transposed;
  for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
      r[j] -= a[i*dof+j] * x[i];

  return;
}
void SolveLinearSystemLapack(double *a, double *b, int N_Eqn, int LDA)
{
// Arguments:
//    a         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
  int m, n, nrhs, lda, ldb;
  int *ipivot, info;
  char t='n';

  m = N_Eqn;
  n = N_Eqn;
  lda = N_Eqn;
  ldb = N_Eqn;
  nrhs = 1;

  ipivot = new int[n];
  
  dgetrf_(&m, &n, a, &lda, ipivot, &info);
  dgetrs_(&t, &n, &nrhs, a, &lda, ipivot, b, &ldb, &info);

  delete ipivot;
}
void SolveLinearSystemTranspose(double *a, double *b, int N_Eqn, int LDA)
{
// Arguments:
//    a         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
  int m, n, nrhs, lda, ldb;
  int *ipivot, info;
  char t='t';

  m = N_Eqn;
  n = N_Eqn;
  lda = N_Eqn;
  ldb = N_Eqn;
  nrhs = 1;

  ipivot = new int[n];
  
  dgetrf_(&m, &n, a, &lda, ipivot, &info);
  dgetrs_(&t, &n, &nrhs, a, &lda, ipivot, b, &ldb, &info);

  delete [] ipivot;
}

/* subroutine for solving a system of linear equations */
void SolveLinearSystem(double *a, double *b, int N_Eqn, int LDA)
// Arguments:
//    a         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
{
  int i,j,k,l,m, row;
  double pivot, tmp, eps=1e-12;


  OutPut("Use SolveLinearSystemNew !!!"<< endl);
  exit(4711);
  int ii, jj;
/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    for(jj=0;jj<N_Eqn;jj++)
    {
      cout << "a(" << ii+1 << "," << jj+1 << ") = " << setw(8);
      cout << AT(ii,jj) << ";" << endl;
    }
  }
  cout << endl;
*/
  for(i=0;i<N_Eqn-1;i++)
  {
    pivot = 0;
    row = i;
    for(l=i;l<N_Eqn;l++)
    {
      tmp = fabs(AT(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    //cout << i << " " << pivot << endl;
    if(pivot < eps) 
    {  
      OutPut("pivot " << pivot << " < " << eps);
      OutPut(" Error in solving linear System " << __FILE__ << endl);
      OutPut("equation: "<< i << endl);
      exit(4711);
    }
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
        pivot = AT(i,l);
        AT(i,l) = AT(row,l);
        AT(row, l) = pivot;
      }
      tmp = b[i];
      b[i] = b[row];
      b[row] = tmp;
    } // endif

    tmp = AT(i,i);
    for(l=i+1;l<N_Eqn;l++)
      AT(l,i) /= tmp;

    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = AT(i,l);
      for(m=i+1;m<N_Eqn;m++)
        AT(m,l) -= AT(m,i) * tmp;
    }

  } // endfor i

  for(i=0;i<N_Eqn;i++)
  {
    tmp = b[i];
    for(l=i+1;l<N_Eqn;l++)
      b[l] -= AT(l,i)*tmp;
  }

  for(i=N_Eqn-1;i>=0;i--)
  {
    b[i] /= AT(i,i);
    tmp = b[i];
    for(l=0;l<i;l++)
      b[l] -= AT(l,i)*tmp;
  }
}

void SolveMultipleSystemsLapack(double *A, double *B, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    A         double array which contains the matrix row wise
//    B         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
//
// This works for a -- stored row wise
//                b -- stored row wise
// The LAPACK Routines for the transposed matrices are used               
//
// The result will be stored row wise
{
  int info, *ipiv;
  char t='t';

  ipiv = new int[N_Eqn];

  dgetrf_(&N_Eqn, &N_Eqn, A, &LDA, ipiv, &info);
  dgetrs_(&t, &N_Eqn, &N_Rhs, A, &LDA, ipiv, B, &LDB, &info);

  delete ipiv;
}

void SolveMultipleSystemsLapack_old(double *A, double *B, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    A         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    B         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
// This works for a -- stored row wise
//                b -- stored row wise
//                     since the LAPACK routines work for
//                     column wise stored right hand sides,
//                     the matrix b will be transposed and
//                     the output matrix containing the solution
//                     will be transposed back
//                LDB = N_Rhs
// The result will be stored column wise
{
  int info, *ipiv;
  int i, j;
  double *transp;
  char t='t';
  transp = new double[N_Eqn*N_Rhs];
  ipiv = new int[N_Eqn];

  dgetrf_(&N_Eqn, &N_Eqn, A, &LDA, ipiv, &info);

  /*for(i=0;i<N_Eqn;i++)
    for(j=0;j<N_Rhs;j++)
      transp[j+i*N_Rhs]=B[j*N_Eqn+i];
*/
 /*for(i=0;i<N_Eqn*N_Rhs;i++)
      OutPut("B1("<<i<<")= "<<B[i]<<"         T1("<<i<<")= "<<transp[i]<<endl);*/

  dgetrs_(&t, &N_Eqn, &N_Rhs, A, &LDA, ipiv, B, &LDB, &info);

 /*for(i=0;i<N_Eqn*N_Rhs;i++)
      OutPut("B2("<<i<<")= "<<B[i]<<"         T2("<<i<<")= "<<transp[i]<<endl);*/


/*  for(j=0;j<N_Rhs;j++)
    for(i=0;i<N_Eqn;i++)
      B[j*N_Eqn+i]=transp[j+i*N_Rhs];
*/
  /*for(i=0;i<N_Eqn*N_Rhs;i++)
      OutPut("B3("<<i<<")= "<<B[i]<<"         T3("<<i<<")= "<<transp[i]<<endl);*/


//exit(1);

  delete transp;
  delete ipiv;
}


/* subroutine for solving a multiple systems of linear equations */
void SolveMultipleSystemsNew(double *a, double *b, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    a         double array which contains the matrix row wise
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
// This works for a -- stored row wise
//                b -- stored column wise, b = (b1_1,..., b_NRhs_1, b1_2, ...)
//                LDB = N_Rhs
// The result will be stored column wise
{
  int i,j,k,l,m, row,info, *ipiv;
  double pivot, tmp, *f, *frow,eps=0.0;

  int ii, jj;
/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    for(jj=0;jj<N_Eqn;jj++)
    {
      cout << "a(" << ii+1 << "," << jj+1 << ") = " << setw(8);
      cout << A(ii,jj) << ";" << endl;
    }
  }
  cout << endl;
  for(jj=0;jj<N_Rhs;jj++)
  {
    cout << jj+1 << ": ";
    for(ii=0;ii<N_Eqn;ii++)
      cout << setw(15) << b[ii*LDB+jj];
    cout << endl;
  }
*/

  // LU decomposition of matrix A with pivot search
  for(i=0;i<N_Eqn-1;i++)
  {
    pivot = 0;
    row = i;
    // find pivot
    for(l=i;l<N_Eqn;l++)
    {
      tmp = fabs(A(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    if(pivot <= eps)
    { 
      OutPut("Error in solving multiple Systems " << __FILE__ << endl);
      OutPut("equation: " << i << endl);
      Error("Error in solving multiple Systems " << __FILE__ << endl);
      exit(4711);
    }
    // change rows if necessary
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
        pivot = A(i,l);
        A(i,l) = A(row,l);
        A(row, l) = pivot;
      }
      // row-th row of rhs
      frow = b+row*LDB;
      // i-th row of rhs
      f = b+i*LDB;          
      for(j=0;j<N_Rhs;j++)
      {
        tmp = f[j];
        f[j] = frow[j];
        frow[j] = tmp;
      }
    } // endif

    tmp = A(i,i);
    for(l=i+1;l<N_Eqn;l++)
      A(l,i) /= tmp;

    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = A(i,l);
      for(m=i+1;m<N_Eqn;m++)
        A(m,l) -= A(m,i) * tmp;
    }

  } // endfor i

  frow = b;
   // solve left lower system
  for(i=0;i<N_Eqn;i++)
  { 
    // i-th row of rhs
    f = b+i*LDB;
    // for all rhs
    for(j=0;j<N_Rhs;j++)
    {    
      tmp = f[j];
      // for remaining rows
      for(l=i+1;l<N_Eqn;l++)
        {
          frow[l*LDB+j] -= A(l,i)*tmp; 
        }
    }
  }
  // solve right upper system
  for(i=N_Eqn-1;i>=0;i--)
  {
    // i-th row of rhs
    f = b+i*LDB;    
    // for all rhs
    for(j=0;j<N_Rhs;j++)
     {
       f[j] /= A(i,i);
       tmp = f[j];
       for(l=0;l<i;l++)
         {
           frow[l*LDB+j] -= A(l,i)*tmp;
         }          
    }
  }
}

void SolveMultipleSystems(double *a, double *b, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    a         double array which contains the matrix row wise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
//
// This works for a -- stored row wise
//                b -- stored row wise, b = (b1,b2,b3, ...)
//                LDA = LDB
// The result will be stored row wise
{
  int i,j,k,l,m, row;
  double pivot, tmp, *f;

  int ii, jj;

/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    cout << ii;
    for(jj=0;jj<N_Eqn;jj++)
      cout << setw(8) << A(ii,jj);
    cout << endl;
  }
  cout << endl;

  cout << "number of Rhs: " << N_Rhs << endl;
*/

  // for all columns
  for(i=0;i<N_Eqn-1;i++)
  {
    // compute pivot element
    pivot = 0;
    row = i;
    // check current column from diagonal element
    for(l=i;l<N_Eqn;l++)
    {
      // a_li 
      tmp = fabs(A(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    if(pivot == 0.0)
    { 
      OutPut("Error in solving multiple Systems" << __FILE__ << endl);
      Error("Error in solving multiple Systems" << __FILE__ << endl);
      exit(4711);
    }
    // change rows i and 'row' if necessary
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
	// N_Eqn subsequent entries since a is stored row wise
        pivot = A(i,l);
        A(i,l) = A(row,l);
        A(row, l) = pivot;
      }
      // same for rhs
      for(j=0;j<N_Rhs;j++)
      { 
	// since rhs is stored row wise, find corresponding index in each rhs
        f = b+j*LDB;
        tmp = f[i];
        f[i] = f[row];
        f[row] = tmp;
      }
    } // endif

    // apply pivoting
    tmp = A(i,i);
    // current column
    for(l=i+1;l<N_Eqn;l++)
      A(l,i) /= tmp;
    // remainder of the matrix
    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = A(i,l);
      for(m=i+1;m<N_Eqn;m++)
        A(m,l) -= A(m,i) * tmp;
    }
  } // endfor i

  for(i=0;i<N_Eqn;i++)
  {
    for(j=0;j<N_Rhs;j++)
    {
      f = b+j*LDB;
      tmp = f[i];
      for(l=i+1;l<N_Eqn;l++)
        f[l] -= A(l,i)*tmp;
    }
  }

  for(i=N_Eqn-1;i>=0;i--)
  {
    for(j=0;j<N_Rhs;j++)
    {
      f = b+j*LDB;
      f[i] /= A(i,i);
      tmp = f[i];
      for(l=0;l<i;l++)
        f[l] -= A(l,i)*tmp;
    }
  }
}

/** calculate the eigenvalue of the system using Lapack routines*/
void FindEigenValues(double *ap, char &UPLO, int N_Eqn, char &COMPZ, double *d, double *z)
{
// Arguments:
//  ap         double precision array which contains the packed upper triangular matrix column wise
//              a[i,j] = a[i +(j-1)*j/2]
//  UPLO    (input) CHARACTER*1
//          = 'U':  Upper triangle of A is packed;
//          = 'L':  Lower triangle of A is packed./**/
//  N_Eqn    : order of the matrix ap 
//  COMPZ   (input) CHARACTER*1
//           = 'N':  Compute eigenvalues only.
//           = 'V':  Compute eigenvalues and eigenvectors of the original
//                   symmetric matrix.  On entry, Z must contain the
//                   orthogonal matrix used to reduce the original matrix
//                   to tridiagonal form.
//           = 'I':  Compute eigenvalues and eigenvectors of the
//                   tridiagonal matrix.  Z is initialized to the identity
//                   matrix.
//  d       (output) DOUBLE PRECISION array, dimension (LDZ, N_Eqn) eigen values
//  z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
//           COMPZ = 'V', Z contains the
//           orthonormal eigenvectors of the original symmetric matrix,
//           and if COMPZ = 'I', Z contains the orthonormal eigenvectors
//           of the symmetric tridiagonal matrix.
//           If COMPZ = 'N', then Z is not referenced.


 int info, i;

//  double *d = new double[N_Eqn];
 double *e = new double[N_Eqn-1];
 double *tau = new double[N_Eqn-1];
//  double *q = new double[N_Eqn*N_Eqn];
 double *work =  new double[2*N_Eqn-1];
 
 dsptrd_(&UPLO, &N_Eqn,  ap, d, e, tau, &info);
  
 dopgtr_(&UPLO, &N_Eqn,  ap, tau, z, &N_Eqn, work, &info);

 dsteqr_(&COMPZ, &N_Eqn, d, e, z, &N_Eqn, work, &info);

 
//   for(i=0; i<N_Eqn; i++)
//    cout<< " Eigen values " << d[i] << endl;
//   
//   for (i=0;i<N_Eqn*N_Eqn;i++)
//     cout<< "z["<<i<<"]="<<z[i]<<endl;

  delete [] e;
  delete [] tau;
  delete [] work;
  
}



/** calculate the determinant of the matrix using Lapack routines*/
double MatrixDeterminant(double *a, int dimension)
{
  
  int m, n, lda;
  int *ipivot, info;
  char t='n';

  m = dimension;
  n = dimension;
  lda = dimension;
  ipivot = new int[n];
  double det =1;
  dgetrf_(&m, &n, a, &lda, ipivot, &info);
  //Determinant Evaluation
  for(int i=0; i<dimension ; i++) det *= a[dimension*i +i];

  int j;
  double detp=1.;
  for( j=0;j <n;j++){
    if(j+1!=ipivot[j]){
        // j+1 : following feedback of ead : ipiv is from Fortran, hence starts at 1.
        // hey ! This is a transpose !
        detp=-detp;
    }
  }
  det = det*detp;
 
  delete ipivot;
  return det;
}



/** calculate the inverse of a 3x3 matrix using Lapack routines*/
void MatrixInverse(double *a, int dimension, int flag)
{
  
  int m, n, lda;
  int *ipivot, info, LWORK;
  double *WORK = new double[3]; 
  char t='n';

  m = 3;
  n = 3;
  lda = 3;
  LWORK = 3;
  ipivot = new int[n];
 
  dgetrf_(&m, &n, a, &lda, ipivot, &info);
  dgetri_(&m, a, &lda, ipivot, WORK, &LWORK, &info );
  //cout<<"Inverse of passed matrix is :"<<endl;
  //for(int i =0; i<= 8; i++) cout<<a[i]<<" ";
  delete ipivot;
  //return a;
}


/** calculate a 3x3  Matrix-Matrix multiplication*/
void MatrixMult(double *a, double *b, double *c, char TRANSA, char TRANSB)
{  
  
  int m = 3, n = 3, k = 3;
  int lda =3, ldb =3, ldc =3;
  double beta = 1;
  double ALPHA = 1;
  //double c[9] = {0,0,0,0,0,0,0,0,0};
  //cout<<endl<<"Matrix-matrix multiplication is :"<<endl;
  dgemm_(&TRANSA, &TRANSB, &m, &n, &k, &ALPHA, a, &lda, b, &ldb, &beta, c, &ldc);
  //for(int i =0; i<= 8; i++) cout<<c[i]<<" ";
}


