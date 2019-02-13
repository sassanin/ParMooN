
/* zfft1d Example Program Text */
/*
 * ACML version 1.0 Copyright AMD,NAG 2003
 */

#include <acml.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
#define NMAX 20
  int i, info, j, n;
  doublecomplex x[NMAX], xx[NMAX];
  doublecomplex *comm;

  printf("ACML example: FFT of a complex sequence using zfft1d\n");
  printf("----------------------------------------------------\n");
  printf("\n");

  /* The sequence of doublecomplex data */
  n = 7;
  x[0] = compose_doublecomplex(0.34907,-0.37168);
  x[1] = compose_doublecomplex(0.54890,-0.35669);
  x[2] = compose_doublecomplex(0.74776,-0.31174);
  x[3] = compose_doublecomplex(0.94459,-0.23702);
  x[4] = compose_doublecomplex(1.13850,-0.13274);
  x[5] = compose_doublecomplex(1.32850,0.00074);
  x[6] = compose_doublecomplex(1.51370,0.16298);
  for (i = 0; i < n; i++)
    xx[i] = x[i];

  /* Allocate communication work array */
  comm = (doublecomplex *)malloc((3*n+100)*sizeof(doublecomplex));

  /* Initialize communication work array */
  zfft1d(0,n,x,comm,&info);

  /* Compute a forward transform */
  zfft1d(-1,n,x,comm,&info);

  printf("Components of discrete Fourier transform:\n");
  printf("\n");
  printf("          Real    Imag\n");
  for (j = 0; j < n; j++)
    printf("%4d   (%7.4f,%7.4f)\n", j, x[j].real, x[j].imag);

  /* Compute the reverse transform */
  zfft1d(1,n,x,comm,&info);

  printf("\n");
  printf("Original sequence as restored by inverse transform:\n");
  printf("\n");
  printf("            Original            Restored\n");
  printf("          Real    Imag        Real    Imag\n");
  for (j = 0; j < n; j++)
    printf("%4d   (%7.4f,%7.4f)   (%7.4f,%7.4f)\n", j, xx[j].real, xx[j].imag, x[j].real, x[j].imag);

  free(comm);
  return 0;
}
