
/* sgetrf Example Program Text */
/*
 * ACML version 1.0 Copyright AMD,NAG 2003
 */

#include <acml.h>
#include <stdio.h>

int main(void)
{
#define NMAX 8
#define NRHMAX 8
  const int lda=NMAX, ldb=NMAX;
  int i, info, j, n, nrhs;
  float a[NMAX*NMAX], b[NMAX*NRHMAX];
  int ipiv[NMAX];

  /* These macros allow access to 1-d arrays as though
     they are 2-d arrays stored in column-major order,
     as required by ACML C routines. */
#define A(I,J) a[((J)-1)*lda+(I)-1]
#define B(I,J) b[((J)-1)*ldb+(I)-1]

  printf("ACML example: solution of linear equations using sgetrf/sgetrs\n");
  printf("--------------------------------------------------------------\n");
  printf("\n");

  /* Initialize matrix A */
  n = 4;
  A(1,1) = 1.80;
  A(1,2) = 2.88;
  A(1,3) = 2.05;
  A(1,4) = -0.89;
  A(2,1) = 5.25;
  A(2,2) = -2.95;
  A(2,3) = -0.95;
  A(2,4) = -3.80;
  A(3,1) = 1.58;
  A(3,2) = -2.69;
  A(3,3) = -2.90;
  A(3,4) = -1.04;
  A(4,1) = -1.11;
  A(4,2) = -0.66;
  A(4,3) = -0.59;
  A(4,4) = 0.80;

  /* Initialize right-hand-side matrix B */
  nrhs = 2;
  B(1,1) = 9.52;
  B(1,2) = 18.47;
  B(2,1) = 24.35;
  B(2,2) = 2.25;
  B(3,1) = 0.77;
  B(3,2) = -13.28;
  B(4,1) = -6.22;
  B(4,2) = -6.21;

  printf("Matrix A:\n");
  for (i = 1; i <= n; i++)
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

  /* Factorize A */
  sgetrf(n,n,a,lda,ipiv,&info);

  printf("\n");
  if (info == 0)
    {
      /* Compute solution */
      sgetrs('N',n,nrhs,a,lda,ipiv,b,ldb,&info);
      /* Print solution */
      printf("Solution matrix X of equations A*X = B:\n");
      for (i = 1; i <= n; i++)
        {
          for (j = 1; j <= nrhs; j++)
            printf("%8.4f ", B(i,j));
          printf("\n");
        }
    }
  else
    printf("The factor U of matrix A is singular\n");

  return 0;
}
