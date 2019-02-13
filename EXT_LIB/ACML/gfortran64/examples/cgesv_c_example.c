
/* cgesv Example Program Text */
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
  complex a[NMAX*NMAX], b[NMAX*NRHMAX];
  int ipiv[NMAX];

  /* These macros allow access to 1-d arrays as though
     they are 2-d arrays stored in column-major order,
     as required by ACML C routines. */
#define A(I,J) a[((J)-1)*lda+(I)-1]
#define B(I,J) b[((J)-1)*ldb+(I)-1]

  printf("\n");
  printf("ACML example: solution of linear equations using cgesv\n");
  printf("------------------------------------------------------\n");
  printf("\n");

  /* Initialize matrix A */
  n = 4;
  A(1,1) = compose_complex(-1.34, 2.55);
  A(1,2) = compose_complex( 0.28, 3.17);
  A(1,3) = compose_complex(-6.39,-2.20);
  A(1,4) = compose_complex( 0.72,-0.92);
  A(2,1) = compose_complex(-0.17,-1.41);
  A(2,2) = compose_complex( 3.31,-0.15);
  A(2,3) = compose_complex(-0.15, 1.34);
  A(2,4) = compose_complex( 1.29, 1.38);
  A(3,1) = compose_complex(-3.29,-2.39);
  A(3,2) = compose_complex(-1.91, 4.42);
  A(3,3) = compose_complex(-0.14,-1.35);
  A(3,4) = compose_complex( 1.72, 1.35);
  A(4,1) = compose_complex( 2.41, 0.39);
  A(4,2) = compose_complex(-0.56, 1.47);
  A(4,3) = compose_complex(-0.83,-0.69);
  A(4,4) = compose_complex(-1.96, 0.67);

  /* Initialize right-hand-side matrix B */
  nrhs = 2;
  B(1,1) = compose_complex(26.98, 50.86);
  B(1,2) = compose_complex(31.32, -6.70);
  B(2,1) = compose_complex( 7.72, -7.30);
  B(2,2) = compose_complex(15.86, -1.42);
  B(3,1) = compose_complex(-4.03, 26.66);
  B(3,2) = compose_complex(-2.15, 30.19);
  B(4,1) = compose_complex(-0.80,  3.24);
  B(4,2) = compose_complex(-2.56,  7.55);

  printf("Matrix A:\n");
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= n; j++)
        printf(" (%8.4f,%8.4f)", A(i,j).real, A(i,j).imag);
      printf("\n");
    }

  printf("\n");
  printf("Right-hand-side matrix B:\n");
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= nrhs; j++)
        printf(" (%8.4f,%8.4f)", B(i,j).real, B(i,j).imag);
      printf("\n");
    }

  /* Factorize A */
  cgesv(n,nrhs,a,lda,ipiv,b,ldb,&info);

  printf("\n");
  if (info == 0)
    {
      /* Print solution */
      printf("Solution matrix X of equations A*X = B:\n");
      for (i = 1; i <= n; i++)
        {
          for (j = 1; j <= nrhs; j++)
            printf(" (%8.4f,%8.4f)", B(i,j).real, B(i,j).imag);
          printf("\n");
        }
    }
  else
    printf("The matrix A is singular\n");

  return 0;
}
