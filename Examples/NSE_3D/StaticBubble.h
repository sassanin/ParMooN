// Stokes problem, two-phase Static Bubble

void ExampleFile()
{
  OutPut("Example: StaticBubble.h" << endl) ;
  
  TDatabase::ParamDB->INTERFACE_FLOW=TRUE;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


// void ExactNull(double x, double y, double z, double *values)
// {
//   values[0] =0;
//   values[1] =0;
//   values[2] =0;
//   values[3] =0;
//   values[4] =0;
// }

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{
  
  if(CompID!=6)
   { cond = DIRICHLET; }
  else
   { cond = INTERFACE; }
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
    value = 0 ;     
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

 if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      coeff[0] = 1.;
      coeff[1] = 0; // f1
      coeff[2] = 0; // f2
      coeff[3] = 0; // f3
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
      coeff[1] = 0; // f1
      coeff[2] = 0; // f2
      coeff[3] = 0; // f3
    }
  }
    
}
