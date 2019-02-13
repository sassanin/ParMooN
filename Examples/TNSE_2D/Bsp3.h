// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (t*x, -t*y)
// p(x,y) = 0

void ExampleFile()
{
  OutPut("Example: Bsp3.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*x;
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -y*t;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*x;
  values[1] = t;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -y*t;
  values[1] = 0;
  values[2] = -t;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = t*Param;
    break;
    case 1:
      value = t*1;
    break;
    case 2:
      value = t*(1-Param);
    break;
    case 3:
      value = 0;
    break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 0;
    break;
    case 1:
      value = -t*Param;
    break;
    case 2:
      value = -t;
    break;
    case 3:
      value = -t*(1-Param);
    break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double t1, t2, t3, t5, t6, t8, t9, t10, t14, t16;
  double t17, t22, t23, t35, t38, t40;
  int i;
  double *coeff, x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    coeff[1] = x*(1+t*t);
    coeff[2] = y*(t*t-1);
  }
}


