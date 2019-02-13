// ======================================================================
// instationary problem
// ======================================================================

#define __SINCOS1__

/// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: SinCos1.h" << endl); 
  TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 0;
}
// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = sin(t*t)*cos(x*y*y);
  values[1] = -sin(t*t)*sin(x*y*y)*y*y;
  values[2] = -sin(t*t)*sin(x*y*y)*2*x*y;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  switch(BdComp)
  {
    case 0: value = sin(t*t);
            break;
    case 1: value = sin(t*t)*cos(Param*Param);
            break;
    case 2: value = sin(t*t)*cos(1-Param);
            break;
    case 3: value = sin(t*t);
            break;
  }
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = sin(t*t)*cos(x*y*y);
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->PE_NR;
  double b1=2., b2=-1., c=1.;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  // previous discrete time
  tau = t-tau;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;

    coeff[4] = cos(t*t)*2*t*cos(x*y*y)
       - eps* sin(t*t)*(-cos(x*y*y)*(y*y*y*y+4*x*x*y*y) - sin(x*y*y)*2*x)  
       - b1*sin(t*t)*sin(x*y*y)*y*y - b2*sin(t*t)*sin(x*y*y)*2*x*y 
       + c*sin(t*t)*cos(x*y*y);
    // rhs from previous time step
    coeff[5] = cos(tau*tau)*2*tau*cos(x*y*y)
       - eps* sin(tau*tau)*(-cos(x*y*y)*(y*y*y*y+4*x*x*y*y) - sin(x*y*y)*2*x)  
       - b1*sin(tau*tau)*sin(x*y*y)*y*y - b2*sin(tau*tau)*sin(x*y*y)*2*x*y 
       + c*sin(tau*tau)*cos(x*y*y);
  }
}
