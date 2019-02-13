
/* cdotu Example Program Text */
/*
 * ACML version 1.0 Copyright AMD,NAG 2003
 */

#include <acml.h>
#include <stdio.h>

int main(void)
{
  complex r;
  complex x[3], y[3];
  int n = 3;
  int i;

  printf("ACML example: dot product of two complex vectors using cdotu\n");
  printf("------------------------------------------------------------\n");
  printf("\n");

  x[0] = compose_complex(1.0, 2.0);
  x[1] = compose_complex(2.0, 1.0);
  x[2] = compose_complex(1.0, 3.0);
  y[0] = compose_complex(3.0, 1.0);
  y[1] = compose_complex(1.0, 4.0);
  y[2] = compose_complex(1.0, 2.0);

  printf("Vector x: ");
  for (i = 0; i < n; i++)
    printf(" (%7.4f,%7.4f)\n", x[i].real, x[i].imag);
  printf("Vector y: ");
  for (i = 0; i < n; i++)
    printf(" (%7.4f,%7.4f)\n", y[i].real, y[i].imag);

  r = cdotu(n, x, 1, y, 1);

  printf("r = x.y = (%12.3f,%12.3f)\n", r.real, r.imag);
  return 0;
}
