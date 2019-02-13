/*
** A program to test the C interfaces
** supplied for the fast scalar and array libm routines.
*/
#include <stdio.h>
#include <acml_mv.h>

int main()
{
  double A,B,C[5],D[5],E[5],S;
  float X,Y,U[8],V[8],W[8],SF,CF[1];

  /* Test double precision routines */

  printf ("ACML vector math library example C program\n");
  printf ("------------------------------------------\n");


  printf ("\nTesting double precision scalar functions\n");

  A=3.0;
  printf ("   fastexp(%9.3f) = %9.3f\n",A,fastexp(A));
  A=10.0;
  printf ("   fastlog(%9.3f) = %9.3f\n",A,fastlog(A));
  printf (" fastlog10(%9.3f) = %9.3f\n",A,fastlog10(A));
  printf ("  fastlog2(%9.3f) = %9.3f\n",A,fastlog2(A));
  A=3.0;
  B=2.5;
  printf ("   fastpow(%9.3f,%9.3f) = %9.3f\n",A,B,fastpow(A,B));
  A=3.14159265358979;
  printf ("   fastsin(%9.3f) = %9.3f\n",A,fastsin(A));
  printf ("   fastcos(%9.3f) = %9.3f\n",A,fastcos(A));
  fastsincos(A,&S,C);
  printf ("fastsincos(%g,*s,*c) = %9.3f,%9.3f\n",A,S,C[0]);
  A=1000000.0;
  printf ("   fastsin(%7.1f) = %9.3f\n",A,fastsin(A));

  printf ("\nTesting double precision array functions\n");
  C[0] = 2.0;
  C[1] = 6.0;
  C[2] = 8.0;
  C[3] = 9.0;
  C[4] = 9.5;
  vrda_exp(4,C,D);
  printf ("   vrda_exp(%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",D[0],D[1],D[2],D[3]);
  vrda_exp(5,C,D);
  printf ("   vrda_exp(%9.3f,%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3],C[4]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f, %9.3f\n",D[0],D[1],D[2],D[3],D[4]);
  vrda_log(4,C,D);
  printf ("   vrda_log(%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",D[0],D[1],D[2],D[3]);
  vrda_log10(4,C,D);
  printf (" vrda_log10(%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",D[0],D[1],D[2],D[3]);
  vrda_log2(4,C,D);
  printf ("  vrda_log2(%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",D[0],D[1],D[2],D[3]);
  vrda_sin(4,C,D);
  printf ("   vrda_sin(%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",D[0],D[1],D[2],D[3]);
  vrda_cos(4,C,D);
  printf ("   vrda_cos(%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",D[0],D[1],D[2],D[3]);
  vrda_sincos(4,C,D,E);
  printf ("vrda_sincos(%9.3f,%9.3f,%9.3f,%9.3f) = \n",C[0],C[1],C[2],C[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",D[0],D[1],D[2],D[3]);
  printf ("            %9.3f,%9.3f,%9.3f,%9.3f\n",E[0],E[1],E[2],E[3]);


  /* Test single precision routines */

  printf ("\nTesting single precision scalar functions\n");


  X=3.0;
  printf ("   fastexpf(%9.3f) = %9.3f\n",X,fastexpf(X));
  X=10.0;
  printf ("   fastlogf(%9.3f) = %9.3f\n",X,fastlogf(X));
  printf (" fastlog10f(%9.3f) = %9.3f\n",X,fastlog10f(X));
  printf ("  fastlog2f(%9.3f) = %9.3f\n",X,fastlog2f(X));
  X=3.0;
  Y=2.5;
  printf ("   fastpowf(%9.3f,%9.3f) = %9.3f\n",X,Y,fastpowf(X,Y));


  X=3.0;
  printf ("   fastsinf(%9.3f) = %9.3f\n",X,fastsinf(X));
  X=3.0;
  printf ("   fastcosf(%9.3f) = %9.3f\n",X,fastcosf(X));
  X=3.0;
  fastsincosf(X,&SF,CF);
  printf ("fastsincosf(%9.3f) = %9.3f,%9.3f\n",X,SF,CF[0]);


  printf ("\nTesting single precision array functions\n");
  U[0] = 1.0;
  U[1] = 2.0;
  U[2] = 3.0;
  U[3] = 4.0;
  U[4] = 6.0;
  U[5] = 7.0;
  U[6] = 8.0;
  U[7] = 9.0;


  vrsa_expf(8,U,V);
  printf ("   vrsa_expf(%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f) =\n",U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7]);
  printf ("             %15.3f,%15.3f,%15.3f,%15.3f,\n"
          "             %15.3f,%15.3f,%15.3f,%15.3f\n\n",V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7]);
  vrsa_logf(8,U,V);
  printf ("   vrsa_logf(%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f) =\n",U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7]);
  printf ("             %7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n\n",V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7]);

  vrsa_log10f(8,U,V);
  printf (" vrsa_log10f(%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f) =\n",U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7]);
  printf ("             %7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n\n",V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7]);

  vrsa_log2f(8,U,V);
  printf ("  vrsa_log2f(%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f) =\n",U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7]);
  printf ("             %7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f\n\n",V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7]);


  vrsa_sinf(8,U,V);
  printf ("   vrsa_sinf(%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f) =\n",U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7]);
  printf ("          %15.3f,%15.3f,%15.3f,%15.3f,\n"
          "          %15.3f,%15.3f,%15.3f,%15.3f\n\n",V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7]);

  vrsa_cosf(8,U,V);
  printf ("   vrsa_cosf(%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f) =\n",U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7]);
  printf ("          %15.3f,%15.3f,%15.3f,%15.3f,\n"
          "          %15.3f,%15.3f,%15.3f,%15.3f\n\n",V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7]);

  vrsa_sincosf(8,U,V,W);
  printf ("vrsa_sincosf(%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f) =\n",U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7]);
  printf ("          %15.3f,%15.3f,%15.3f,%15.3f,\n"
          "          %15.3f,%15.3f,%15.3f,%15.3f\n",V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7]);
  printf ("          %15.3f,%15.3f,%15.3f,%15.3f,\n"
          "          %15.3f,%15.3f,%15.3f,%15.3f\n\n",W[0],W[1],W[2],W[3],W[4],W[5],W[6],W[7]);

  vrsa_powf(8,U,V,W);
  printf ("   vrsa_powf(%7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f,\n",U[0],V[0],U[1],V[1],U[2],V[2],U[3],V[3]);
  printf ("             %7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f) =\n",U[4],V[4],U[5],V[5],U[6],V[6],U[7],V[7]);
  printf ("          %15.3f,%15.3f,%15.3f,%15.3f,\n"
          "          %15.3f,%15.3f,%15.3f,%15.3f\n\n",W[0],W[1],W[2],W[3],W[4],W[5],W[6],W[7]);

  vrsa_powxf(8,U,V[0],W);
  printf ("  vrsa_powxf(%7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f,\n",U[0],V[0],U[1],V[0],U[2],V[0],U[3],V[0]);
  printf ("             %7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f,%7.3f^%7.3f) =\n",U[4],V[0],U[5],V[0],U[6],V[0],U[7],V[0]);
  printf ("          %15.3f,%15.3f,%15.3f,%15.3f,\n"
          "          %15.3f,%15.3f,%15.3f,%15.3f\n\n",W[0],W[1],W[2],W[3],W[4],W[5],W[6],W[7]);

  return(0);
}
