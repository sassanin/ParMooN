// Navier-Stokes problem, solution in ansatz space
// velocity pw quadratic, pressure linear
// 
// u(x,y) = (x^2+y^2+z^2, x^2+2xy+13, -2xz+5y^2)^T
// p(x,y) = 3x-2y+7z-4

void ExampleFile()
{
  OutPut("Example: AnsatzQuadLin.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = x*x+y*y+z*z;
  values[1] = 2*x;
  values[2] = 2*y;
  values[3] = 2*z;
  values[4] = 6;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = x*x+2*x*z+13;
  values[1] = 2*x+2*z;
  values[2] = 0;
  values[3] = 2*x;
  values[4] = 2;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = -2*x*z+5*y*y;
  values[1] = -2*z;
  values[2] = 10*y;
  values[3] = -2*x;
  values[4] = 10;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-2*y+7*z-4;
  values[1] = 3;
  values[2] = -2;
  values[3] = 7;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = x*x+y*y+z*z;
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = x*x+2*x*z+13;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = -2*x*z+5*y*y;
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = x*x+y*y+z*z;
  values[1] = 2*x;
  values[2] = 2*y;
  values[3] = 2*z;
  values[4] = 6;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] = x*x+2*x*z+13;
  values[1] = 2*x+2*z;
  values[2] = 0;
  values[3] = 2*x;
  values[4] = 2;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = -2*x*z+5*y*y;
  values[1] = -2*z;
  values[2] = 10*y;
  values[3] = -2*x;
  values[4] = 10;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-2*y+7*z-4;
  values[1] = 3;
  values[2] = -2;
  values[3] = 7;
  values[4] = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

  if (TDatabase::ParamDB->STOKES_PROBLEM)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      coeff[0] = eps;
      coeff[1] = -6*eps + 3; // f1
      coeff[2] = -2*eps - 2; // f2
      coeff[3] = -10*eps + 7; // f3
    }
  }
  else
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      x = X[i];
      y = Y[i];
      z = Z[i];
      coeff[0] = eps;
      coeff[1] = -6*eps + 3 + 2*x*x*x + 2*x*y*y-2*x*z*z +2*y*x*x+ 4*x*y*z
        +26*y + 10 *y*y*z; // f1
      coeff[2] = -2*eps - 2 +2*x*x*x -2 *x*x*z + 12*x*y*y + 2*z*y*y +2*x*z*z
        + 2*z*z*z; // f2
      coeff[3] = -10*eps + 7 + 2*x*x*z -2*z*y*y-2*z*z*z + 10 *y*x*x 
        + 20*x*y*z + 130*y -10 *x*y*y ; // f3
    }
  }
    
}
