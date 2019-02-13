// ======================================================================
// sharp characteristic interior layer
// Knopp, Lube, Rapin, CMAME 2002
// ======================================================================
#define __TWO_INTERIOR_LAYERS__

#include <MooNMD_Io.h>
#include <Constants.h>

void ExampleFile()
{
  OutPut("Example: TwoInteriorLayers.h" << endl) ;
}
// exact solution
void Exact(double x, double y, double *values)
{

  values[0]= 0;
  values[1]= 0;
  values[2]= 0;
  values[3]= 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
    if (i==3)
	cond = NEUMANN;
    else
	cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
   if (BdComp==0)
   {
      if ((Param>1.0/3.0)&& (Param<2.0/3.0))
         value = 1;
      else
         value = 0;
   }
   else
      value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;
  int i;
  double *coeff;
  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = -y;
    coeff[2] = x;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = sqrt(x*x + y*y);
  }
}


 