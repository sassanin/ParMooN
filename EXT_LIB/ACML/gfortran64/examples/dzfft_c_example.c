
/* dzfft Example Program Text */
/*
 * ACML version 2.0 Copyright AMD,NAG 2004
 */

#include <acml.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
#define NMAX 20
  int i, info, j, n;
  double x[NMAX], xx[NMAX];
  double *comm;

  printf("ACML example: FFT of a real sequence using dzfft\n");
  printf("------------------------------------------------\n");
  printf("\n");

  /* The sequence of double data */
  n = 7;
  x[0] = 0.34907;
  x[1] = 0.54890;
  x[2] = 0.74776;
  x[3] = 0.94459;
  x[4] = 1.13850;
  x[5] = 1.32850;
  x[6] = 1.51370;
  for (i = 0; i < n; i++)
    xx[i] = x[i];

  /* Allocate communication work array */
  comm = (double *)malloc((3*n+100)*sizeof(double));

  /* Initialize communication work array */
  dzfft(0,n,x,comm,&info);

  /* Compute a real --> Hermitian transform */
  dzfft(1,n,x,comm,&info);

  printf("Components of discrete Fourier transform:\n");
  printf("\n");
  for (j = 0; j < n; j++)
    printf("%4d   %7.4f\n", j, x[j]);

  /* Conjugate the Vector X to simulate inverse transform */
  for (j = n/2+1; j < n; j++)
    x[j] = -x[j];

  /* Compute the complex Hermitian --> real transform */
  zdfft(2,n,x,comm,&info);

  printf("\n");
  printf("Original sequence as restored by inverse transform:\n");
  printf("\n");
  printf("       Original  Restored\n");
  for (j = 0; j < n; j++)
    printf("%4d   %7.4f   %7.4f\n", j, xx[j], x[j]);

  free(comm);
  return 0;
}
