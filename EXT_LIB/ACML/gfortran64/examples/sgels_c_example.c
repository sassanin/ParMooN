
/* sgels Example Program Text */
/*
 * ACML version 1.0 Copyright AMD,NAG 2003
 */

#include <acml.h>
#include <stdio.h>

int main(void)
{
#define MMAX 8
#define NMAX 8
#define NRHMAX 8
  const int lda=MMAX, ldb=NMAX;
  int i, info, j, m, n, nrhs;
  float a[MMAX*NMAX], b[NMAX*NRHMAX];

  /* These macros allow access to 1-d arrays as though
     they are 2-d arrays stored in column-major order,
     as required by ACML C routines. */
#define A(I,J) a[((J)-1)*lda+(I)-1]
#define B(I,J) b[((J)-1)*ldb+(I)-1]

  printf("ACML example: least-squares solution of overdetermined system of\n");
  printf("              linear equations using sgels\n");
  printf("----------------------------------------------------------------\n");
  printf("\n");

  /* Initialize matrix A */
  m = 6;
  n = 4;
  A(1,1) = -0.57;
  A(1,2) = -1.28;
  A(1,3) = -0.39;
  A(1,4) = 0.25;
  A(2,1) = -1.93;
  A(2,2) = 1.08;
  A(2,3) = -0.31;
  A(2,4) = -2.14;
  A(3,1) = 2.30;
  A(3,2) = 0.24;
  A(3,3) = 0.40;
  A(3,4) = -0.35;
  A(4,1) = -1.93;
  A(4,2) = 0.64;
  A(4,3) = -0.66;
  A(4,4) = 0.08;
  A(5,1) = 0.15;
  A(5,2) = 0.30;
  A(5,3) = 0.15;
  A(5,4) = -2.13;
  A(6,1) = -0.02;
  A(6,2) = 1.03;
  A(6,3) = -1.43;
  A(6,4) = 0.50;

  /* Initialize right-hand-side matrix B */
  nrhs = 2;
  B(1,1) = -3.15;
  B(1,2) = 2.19;
  B(2,1) = -0.11;
  B(2,2) = -3.64;
  B(3,1) = 1.99;
  B(3,2) = 0.57;
  B(4,1) = -2.70;
  B(4,2) = 8.23;
  B(5,1) = 0.26;
  B(5,2) = -6.35;
  B(6,1) = 4.50;
  B(6,2) = -1.48;

  printf("Matrix A:\n");
  for (i = 1; i <= m; i++)
    {
      for (j = 1; j <= n; j++)
        printf("%8.4f ", A(i,j));
      printf("\n");
    }

  printf("\n");
  printf("Right-hand-side matrix B:\n");
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= nrhs; j++)
        printf("%8.4f ", B(i,j));
      printf("\n");
    }

  /* Compute least-squares solution of AX = B */
  sgels('n',m,n,nrhs,a,lda,b,ldb,&info);

  /* Print solution */
  printf("\n");
  printf("Least-squares solution matrix X of equations A*X = B:\n");
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= nrhs; j++)
        printf("%8.4f ", B(i,j));
      printf("\n");
    }

  return 0;
}
