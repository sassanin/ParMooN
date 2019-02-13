// ======================================================================
// Sine problem 3D
// ======================================================================
// #include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: Laplace.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] = Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] = Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] = Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdID, double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
//       cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdID, double x, double y, double z, double &value)
{
  value = sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = (3*Pi*Pi*eps)*sin(Pi*x[i])*sin(Pi*y[i])*sin(Pi*z[i]);
  }
}

