
/* drandleapfrog and drandskipahead Example Program Text */
/*
 * ACML version 3.0 Copyright AMD,NAG 2005
 */

#include <acml.h>
#include <stdio.h>

int main(void)
{
#define MSTATE 20
#define MSEED 4
#define MN 12
 const int mn=12;
  double a, b;
  int genid, i, info, lseed, lstate, n, nskip, subid;
  double x1[MN], x2[MN], x3[MN], x4[MN];
  int seed[MSEED], state1[MSTATE], state2[MSTATE],
    state3[MSTATE], state4[MSTATE];

  printf("ACML example: Multiple streams using"
         " drandskipahead and drandleapfrog\n");
  printf("------------------------------------"
         "---------------------------------\n");
  printf("\n");

  /* Initialize the number of variates required */
  n = mn;

  /* Use the basic generator as the base generator */
  genid = 1;
  subid = 1;

  /* Populate the seed array, basic generator needs one seed,
     and a state array of length 16 */
  lstate = 16;
  lseed = 1;
  seed[0] = 122421;

  /* Generate values from a uniform U(0,1) distribution */
  a = 0.0;
  b = 1.0;

  /* Creating multiple streams using SKIPAHEAD */
  printf("Generating multiple streams using SKIPAHEAD:\n");

  /* Initialize the base generator */
  drandinitialize(genid,subid,seed,&lseed,state1,&lstate,&info);

  /* Going to be using four streams to generate N values,
     so initialize all four to the same values */
  for (i=0; i<lstate; i++)
    {
      state2[i] = state1[i];
      state3[i] = state1[i];
      state4[i] = state1[i];
    }

  /* Calculate the number of values from each stream */
  nskip = n/4;
  if (4*nskip < n)
    nskip++;

  /* Advance each stream the required amount */
  drandskipahead(nskip,state2,&info);
  drandskipahead(2*nskip,state3,&info);
  drandskipahead(3*nskip,state4,&info);

  /* Generate a sequence from all four streams */
  dranduniform(nskip,a,b,state1,x1,&info);
  dranduniform(nskip,a,b,state2,x2,&info);
  dranduniform(nskip,a,b,state3,x3,&info);
  dranduniform(nskip,a,b,state4,x4,&info);

  /* Report the values */
  for (i=0; i<nskip; i++)
    printf("%10.4f%10.4f%10.4f%10.4f\n", x1[i], x2[i], x3[i], x4[i]);

  /* Creating multiple streams using LEAPFROG */
  printf("\n");
  printf("Generating multiple streams using LEAPFROG:\n");

  /* Initialize the base generator
     In practice the generators should only be initialized once per
     application. They are reinitialized here so that the results from
     SKIPAHEAD and LEAPFROG can be compared */
  drandinitialize(genid,subid,seed,&lseed,state1,&lstate,&info);

  /* Going to be using four streams to generate N values,
     so initialize all four to the same values */
  for (i=0; i < lstate; i++)
    {
      state2[i] = state1[i];
      state3[i] = state1[i];
      state4[i] = state1[i];
    }

  /* Calculate the number of values from each stream */
  nskip = n/4;
  if (4*nskip < n)
    nskip++;

  /* Advance each stream the required amount */
  drandleapfrog(4,1,state1,&info);
  drandleapfrog(4,2,state2,&info);
  drandleapfrog(4,3,state3,&info);
  drandleapfrog(4,4,state4,&info);

  /* Generate a sequence from all four streams */
  dranduniform(nskip,a,b,state1,x1,&info);
  dranduniform(nskip,a,b,state2,x2,&info);
  dranduniform(nskip,a,b,state3,x3,&info);
  dranduniform(nskip,a,b,state4,x4,&info);

  /* Report the values */
  for (i=0; i<nskip; i++)
    printf("%10.4f%10.4f%10.4f%10.4f\n", x1[i], x2[i], x3[i], x4[i]);

  return 0;
}
