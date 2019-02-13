
/* drandinitializeuser Example Program Text */
/*
 * ACML version 3.0 Copyright AMD,NAG 2005
 */

#include <acml.h>
#include <stdio.h>

/* Declare prototypes for user-supplied initialization
   and generator routines */
void whini(int *genid, int *subid, int *seed, int *lseed,
           int *state, int *lstate, int *info);
void whgen(int *n, int *state, double *x, int *info);

int main(void)
{
#define MSTATE 28
#define MSEED 10
#define MN 12
  const int mn=MN;
  double a, b;
  int i, info, lseed, lstate, n;
  double x[MN];
  int seed[MSEED], state[MSTATE];

  printf("ACML example: drandinitializeuser\n");
  printf("---------------------------------\n");
  printf("Use of a user-supplied base generator for random numbers\n");
  printf("\n");

  /* Initialize the number of variates required */
  n = mn;

  /* Populate the seed array, and give the length of the state array.
     Note that the state array must be declared 3 elements larger
     than required by the user-supplied whini routine. The extra
     three elements are required by drandinitializeuser itself
     rather than by whini. */
  lstate = 28;
  lseed = 4;
  seed[0] = 10;
  seed[1] = 7;
  seed[2] = 19;
  seed[3] = 58;

  /* Initialize the base generator. */
  drandinitializeuser(whini,whgen,1,1,seed,&lseed,state,&lstate,&info);

  /* Generate a sequence from a uniform U(0,1) distribution */
  a = 0.0;
  b = 1.0;
  dranduniform(n,a,b,state,x,&info);

  /* Print the sequence */
  printf("Numbers from a uniform U(0,1) distribution:\n");
  for (i=0; i<n; i++)
    {
      printf("%10.4f", x[i]);
      if ((i+1)%4 == 0)
        printf("\n");
    }
  printf("\n");
  return 0;
}

void whini(int *genid, int *subid, int *seed, int *lseed,
           int *state, int *lstate, int *info)
{
  /*
     whini - Initialisation routine for the user-supplied
             base generator.

     Reference: the base generator implemented here is described in:
           "Generating good pseudo-random numbers", B.A.Wichmann
           and I.D.Hill (preprint version 0.92, April 2005)

     Input:
     genid  - Generator ID (irrelevant in this example)
     subid  - Sub generator ID (irrelevant in this example)
     seed   - Array of seeds, seed[i] > 0 for i = 0 to lseed-1
     lseed  - Length of array seed. In this example, lseed must be 4.
     lstate - Length of array state, lstate >= 25

     Output:
     state  - Array holding information about the generator:
                    0 = Amount of state used = 25 (i.e. elements 0 to 24)
                    1 = 'Magic' number
                    2 = genid
                    3 = subid
              4 to  7 = Current state of the generator
              8 to 23 = Parameters of constituent parts of the generator.
                   24 = 'Magic' number again
     info   - Failure code:
              0        = Everything OK
              -1 to -6 = Error in parameter abs(info)
  */
  const int eseed=4, estate=25;
  const int magic=54321;
  int i;
  /* Basic parameter checking */
  if (*lseed < eseed)
    {
      *info = -4;
      return;
    }
  else if (*lstate < estate)
    {
      *info = -6;
      return;
    }
  /* End of basic parameter checking */

  /* Check and copy the supplied seed into the state array */
  for (i=0; i<eseed; i++)
    {
      if (seed[i] < 1)
        {
          *info = -3;
          return;
        }
      else
        state[4+i] = seed[i];
    }

  /* Populate the rest of the STATE array */
  state[0] = estate;
  state[1] = magic;
  state[2] = *genid;
  state[3] = *subid;

  /* These define the parameters of the four constituent generators */
  state[8] = 11600;
  state[9] = 185127;
  state[10] = 10379;
  state[11] = 2147483579;
  state[12] = 47003;
  state[13] = 45688;
  state[14] = 10479;
  state[15] = 2147483543;
  state[16] = 23000;
  state[17] = 93368;
  state[18] = 19423;
  state[19] = 2147483423;
  state[20] = 33000;
  state[21] = 65075;
  state[22] = 8123;
  state[23] = 2147483123;

  /* Finally, repeat the magic number */
  state[24] = magic;

  *info = 0;
  return;
}

void whgen(int *n, int *state, double *x, int *info)
{
  /* whgen: Generate a vector of pseudo random numbers
            using the WH generator.

     See whini for reference.

     Input:
     n  - Number of values being generated, n >= 0

     Input / Output:
     state - On Entry: The current state of the generator
             On Exit : The state of the generator after
                       generating the n values

     Output:
     x     - Vector holding the generated values, x[n]
     info  - Error code:
             0        = Everything OK
             -1 to -2 = Error in parameter abs(info)
  */
  const int magic=54321;
  double w;
  int i, imagt, imagx, imagy, imagz, imodt, imodx,
    imody, imodz, imult, imulx, imuly, imulz, isubt,
    isubx, isuby, isubz, it, ix, iy, iz;

  /* Parameter checks */
  if (*n < 0)
    {
      *info = -1;
      return;
    }

  /* Check the magic number indicating that the initialization
     routine has previously been called */
  if (state[1] != magic || state[state[0]-1] != magic)
    {
      *info = -2;
      return;
    }
  /* End of parameter checks */

  /* Extract current state */
  ix = state[4];
  iy = state[5];
  iz = state[6];
  it = state[7];
  imulx = state[8];
  imodx = state[9];
  isubx = state[10];
  imagx = state[11];
  imuly = state[12];
  imody = state[13];
  isuby = state[14];
  imagy = state[15];
  imulz = state[16];
  imodz = state[17];
  isubz = state[18];
  imagz = state[19];
  imult = state[20];
  imodt = state[21];
  isubt = state[22];
  imagt = state[23];

  /* Generate the N values */
  for (i=0; i<*n; i++)
    {
      ix = imulx*(ix%imodx) - isubx*(ix/imodx);
      iy = imuly*(iy%imody) - isuby*(iy/imody);
      iz = imulz*(iz%imodz) - isubz*(iz/imodz);
      it = imult*(it%imodt) - isubt*(it/imodt);
      if (ix < 0)
        ix = ix + imagx;
      if (iy < 0)
        iy = iy + imagy;
      if (iz < 0)
        iz = iz + imagz;
      if (it < 0)
        it = it + imagt;

      /* Convert the current state into a value between 0 and 1 */
      w = ix/(double)imagx + iy/(double)imagy + iz/(double)imagz +
        it/(double)imagt;
      x[i] = w - (int)w;
    }

  /* Put the current generator state into the STATE array */
  state[4] = ix;
  state[5] = iy;
  state[6] = iz;
  state[7] = it;

  *info = 0;
  return;
}
