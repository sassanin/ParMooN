
/* cgels Example Program Text */
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
  complex a[MMAX*NMAX], b[NMAX*NRHMAX];

  /* These macros allow access to 1-d arrays as though
     they are 2-d arrays stored in column-major order,
     as required by ACML C routines. */
#define A(I,J) a[((J)-1)*lda+(I)-1]
#define B(I,J) b[((J)-1)*ldb+(I)-1]

  printf("ACML example: least-squares solution of overdetermined system of\n");
  printf("              linear equations using cgels\n");
  printf("----------------------------------------------------------------\n");
  printf("\n");

  /* Initialize matrix A */
  m = 6;
  n = 4;
  A(1,1) = compose_complex( 0.96,-0.81);
  A(1,2) = compose_complex(-0.03, 0.96);
  A(1,3) = compose_complex(-0.91, 2.06);
  A(1,4) = compose_complex(-0.05, 0.41);
  A(2,1) = compose_complex(-0.98, 1.98);
  A(2,2) = compose_complex(-1.20, 0.19);
  A(2,3) = compose_complex(-0.66, 0.42);
  A(2,4) = compose_complex(-0.81, 0.56);
  A(3,1) = compose_complex( 0.62,-0.46);
  A(3,2) = compose_complex( 1.01, 0.02);
  A(3,3) = compose_complex( 0.63,-0.17);
  A(3,4) = compose_complex(-1.11, 0.60);
  A(4,1) = compose_complex(-0.37, 0.38);
  A(4,2) = compose_complex( 0.19,-0.54);
  A(4,3) = compose_complex(-0.98,-0.36);
  A(4,4) = compose_complex( 0.22,-0.20);
  A(5,1) = compose_complex( 0.83, 0.51);
  A(5,2) = compose_complex( 0.20, 0.01);
  A(5,3) = compose_complex(-0.17,-0.46);
  A(5,4) = compose_complex( 1.47, 1.59);
  A(6,1) = compose_complex( 1.08,-0.28);
  A(6,2) = compose_complex( 0.20,-0.12);
  A(6,3) = compose_complex(-0.07, 1.23);
  A(6,4) = compose_complex( 0.26, 0.35);

  /* Initialize right-hand-side matrix B */
  nrhs = 2;
  B(1,1) = compose_complex(-1.54, 0.76);
  B(1,2) = compose_complex( 3.17,-2.09);
  B(2,1) = compose_complex( 0.12,-1.92);
  B(2,2) = compose_complex(-6.53, 4.18);
  B(3,1) = compose_complex(-9.08,-4.31);
  B(3,2) = compose_complex( 7.28, 0.73);
  B(4,1) = compose_complex( 7.49, 3.65);
  B(4,2) = compose_complex( 0.91,-3.97);
  B(5,1) = compose_complex(-5.63,-2.12);
  B(5,2) = compose_complex(-5.46,-1.64);
  B(6,1) = compose_complex( 2.37, 8.03);
  B(6,2) = compose_complex(-2.84,-5.86);

  printf("Matrix A:\n");
  for (i = 1; i <= m; i++)
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

  /* Compute least-squares solution of AX = B */
  cgels('n',m,n,nrhs,a,lda,b,ldb,&info);

  /* Print solution */
  printf("\n");
  printf("Least-squares solution matrix X of equations A*X = B:\n");
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= nrhs; j++)
        printf(" (%8.4f,%8.4f)", B(i,j).real, B(i,j).imag);
      printf("\n");
    }

  return 0;
}
