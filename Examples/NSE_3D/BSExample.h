// Navier-Stokes problem, solution in ansatz space
// velocity pw quadratic, pressure linear
//
//
void ExampleFile()
{
  OutPut("Example: BSExample.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
  values[1] = cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)+4*x*x*x*cos(Pi*y);
  values[2] = sin(Pi*x)*cos(Pi*y)*Pi*sin(Pi*z)-x*x*x*x*sin(Pi*y)*Pi;
  values[3] = sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi;
  values[4] = -3*sin(Pi*x)*Pi*Pi*sin(Pi*y)*sin(Pi*z)
                       +12*x*x*cos(Pi*y)-x*x*x*x*cos(Pi*y)*Pi*Pi;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
  values[1] = -Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z)-9*y*y*z;
  values[3] = -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z)-3*y*y*y;
  values[4] = -3*cos(Pi*x)*Pi*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
  values[1] = -sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi
    -12*x*x*cos(Pi*y)*z-sin(Pi*z)*sin(Pi*x)*Pi*sin(Pi*y);
  values[2] = cos(Pi*z)*cos(Pi*x)*cos(Pi*y)*Pi
    +4*x*x*x*sin(Pi*y)*Pi*z+cos(Pi*x)*cos(Pi*y)*sin(Pi*z)*Pi+9*y*z*z;
  values[3] = -cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)+cos(Pi*x)*sin(Pi*y)*Pi*cos(Pi*z)+9*y*y*z;
  values[4] = -3*cos(Pi*x)*Pi*Pi*sin(Pi*y)*cos(Pi*z)
                       -24*x*cos(Pi*y)*z
                       -3*sin(Pi*z)*cos(Pi*x)*Pi*Pi*sin(Pi*y)
                       +4*x*x*x*cos(Pi*y)*Pi*Pi*z
                       +9*z*z+9*y*y;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-sin(y+4*z)-1.5-sin(1.0)*cos(1.0)
    -4*cos(1.0)*cos(1.0)*cos(1.0)*cos(1.0)*sin(1.0)
    +sin(1.0)*cos(1.0)*cos(1.0)-2*sin(1.0)*sin(1.0)*sin(1.0)
    +2*sin(1.0)*cos(1.0)*cos(1.0)*cos(1.0)+2*sin(1.0);
  values[1] = 3;
  values[2] = -cos(y+4*z);
  values[3] = -4*cos(y+4*z);
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
//   cond = NEUMANN;
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
}

// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
  values[1] = cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)+4*x*x*x*cos(Pi*y);
  values[2] = sin(Pi*x)*cos(Pi*y)*Pi*sin(Pi*z)-x*x*x*x*sin(Pi*y)*Pi;
  values[3] = sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi;
  values[4] = -3*sin(Pi*x)*Pi*Pi*sin(Pi*y)*sin(Pi*z)
                       +12*x*x*cos(Pi*y)-x*x*x*x*cos(Pi*y)*Pi*Pi;
  values[0] = 0;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
  values[1] = -Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z)-9*y*y*z;
  values[3] = -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z)-3*y*y*y;
  values[4] = -3*cos(Pi*x)*Pi*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z;
  values[0] = 0;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
  values[1] = -sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi
    -12*x*x*cos(Pi*y)*z-sin(Pi*z)*sin(Pi*x)*Pi*sin(Pi*y);
  values[2] = cos(Pi*z)*cos(Pi*x)*cos(Pi*y)*Pi
    +4*x*x*x*sin(Pi*y)*Pi*z+cos(Pi*x)*cos(Pi*y)*sin(Pi*z)*Pi+9*y*z*z;
  values[3] = -cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)+cos(Pi*x)*sin(Pi*y)*Pi*cos(Pi*z)+9*y*y*z;
  values[4] = -3*cos(Pi*x)*Pi*Pi*sin(Pi*y)*cos(Pi*z)
                       -24*x*cos(Pi*y)*z
                       -3*sin(Pi*z)*cos(Pi*x)*Pi*Pi*sin(Pi*y)
                       +4*x*x*x*cos(Pi*y)*Pi*Pi*z
                       +9*z*z+9*y*y;
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-sin(y+4*z)-1.5-sin(1.0)*cos(1.0)
    -4*cos(1.0)*cos(1.0)*cos(1.0)*cos(1.0)*sin(1.0)
    +sin(1.0)*cos(1.0)*cos(1.0)-2*sin(1.0)*sin(1.0)*sin(1.0)
    +2*sin(1.0)*cos(1.0)*cos(1.0)*cos(1.0)+2*sin(1.0);
  values[1] = 3;
  values[2] = -cos(y+4*z);
  values[3] = -4*cos(y+4*z);
  values[4] = 0;
  values[0] = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z, u1, u2, u3, ux, uy, uz;

  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];

      x = X[i];
      y = Y[i];
      z = Z[i];
      coeff[0] = eps;
      coeff[1] = -eps*(-3*sin(Pi*x)*Pi*Pi*sin(Pi*y)*sin(Pi*z)
                       +12*x*x*cos(Pi*y)-x*x*x*x*cos(Pi*y)*Pi*Pi) + 3; // f1
      coeff[2] = -eps*(-3*cos(Pi*x)*Pi*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z) 
        - cos(y+4*z); // f2
      coeff[3] = -eps*(-3*cos(Pi*x)*Pi*Pi*sin(Pi*y)*cos(Pi*z)
                       -24*x*cos(Pi*y)*z
                       -3*sin(Pi*z)*cos(Pi*x)*Pi*Pi*sin(Pi*y)
                       +4*x*x*x*cos(Pi*y)*Pi*Pi*z
                       +9*z*z+9*y*y) 
        - 4*cos(y+4*z);   // f3
    }
  }
  else if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      x = X[i];
      y = Y[i];
      z = Z[i];
      u1 = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
      u2 =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
      u3 = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
      ux =  cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)+4*x*x*x*cos(Pi*y);
      uy = sin(Pi*x)*cos(Pi*y)*Pi*sin(Pi*z)-x*x*x*x*sin(Pi*y)*Pi;
      uz = sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi;
      coeff[1] = -eps*(-3*sin(Pi*x)*Pi*Pi*sin(Pi*y)*sin(Pi*z)
                       +12*x*x*cos(Pi*y)-x*x*x*x*cos(Pi*y)*Pi*Pi) 
        + u1*ux+u2*uy+u3*uz+ 3; // f1

      ux = -Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
      uy = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z)-9*y*y*z;
      uz= -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z)-3*y*y*y;      
      coeff[2] = -eps*(-3*cos(Pi*x)*Pi*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z) 
        + u1*ux+u2*uy+u3*uz - cos(y+4*z); // f2

      ux = -sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi
        -12*x*x*cos(Pi*y)*z-sin(Pi*z)*sin(Pi*x)*Pi*sin(Pi*y);
      uy = cos(Pi*z)*cos(Pi*x)*cos(Pi*y)*Pi
        +4*x*x*x*sin(Pi*y)*Pi*z+cos(Pi*x)*cos(Pi*y)*sin(Pi*z)*Pi+9*y*z*z;
      uz = -cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)
        -4*x*x*x*cos(Pi*y)+cos(Pi*x)*sin(Pi*y)*Pi*cos(Pi*z)+9*y*y*z;
      coeff[3] = -eps*(-3*cos(Pi*x)*Pi*Pi*sin(Pi*y)*cos(Pi*z)
                       -24*x*cos(Pi*y)*z
                       -3*sin(Pi*z)*cos(Pi*x)*Pi*Pi*sin(Pi*y)
                       +4*x*x*x*cos(Pi*y)*Pi*Pi*z
                       +9*z*z+9*y*y) 
        + u1*ux+u2*uy+u3*uz - 4*cos(y+4*z);   // f3
      coeff[0] = eps;
    }
  }
  else
  {  OutPut("Oseen type not included.... Change FLOW_PROBLEM_TYPE"<<endl);
     exit(0);
    
  } 
}
