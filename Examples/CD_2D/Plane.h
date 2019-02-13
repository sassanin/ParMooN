// ======================================================================
// Plane problem 3D
// ======================================================================
// #include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: Plane.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 1+2*x+3*y+4*z;
  values[1] = 2;
  values[2] = 3;
  values[3] = 4;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = 1+2*x+3*y+4*z;
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
    coeff[1] = 0; // 1;
    coeff[2] = 0; // 2;
    coeff[3] = 0; // 3;
    coeff[4] = 0; // 0;
    coeff[5] = 0; // 1*2+2*3+3*4;
  }
}

