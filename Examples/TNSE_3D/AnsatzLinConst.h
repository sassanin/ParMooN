// time dependent Navier-Stokes problem 3D, ansatz
// 
// u(x,y) = t*(y+z, x-z, 2*x+y)
// p(x,y) = 0

void ExampleFile()
{
  OutPut("Example: AnsatzLinConst.h" << endl);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(y+z);
}

void InitialU2(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(x-z);
}

void InitialU3(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(2*x+y);

}

void InitialP(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(y+z);
  values[1] = 0;
  values[2] = t;
  values[3] = t;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(x-z);
  values[1] = t;
  values[2] = 0;
  values[3] = -t;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(2*x+y);
  values[1] = 2*t;
  values[2] = t;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;

}

// ========================================================================
// boundary conditions
// ========================================================================
// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*(y+z);
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*(x-z);
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*(2*x+y);    
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double t=TDatabase::TimeDB->CURRENTTIME;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = eps;
    coeff[1] = (y+z) + t*t*(x-z) + t*t*(2*x+y); // f1
    coeff[2] = (x-z) + t*t*(y+z) - t*t*(2*x+y); // f2
    coeff[3] =  (2*x+y) + 2*t*t*(y+z) + t*t*(x-z);// f3
  }  
}


