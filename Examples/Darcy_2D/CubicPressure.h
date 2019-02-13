// Bechmark for Darcy problem, exact solution is
// 

void ExampleFile()
{
  OutPut("Example: CubicPressure.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = -x*x*y + y*y*y/3.0;
  values[1] = -2*x*y;
  values[2] = -x*x + y*y;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -x*x*x/3.0 + x*y*y;
  values[1] = -x*x + y*y;
  values[2] = 2*x*y;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = (x*x*x*y - x*y*y*y)/3.0;
  values[1] = x*x*y - y*y*y/3.0;
  values[2] = x*x*x/3.0 - x*y*y;
  values[3] = 0.0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int bdComp, double t, BoundCond &cond)  
{
  cond = (bdComp == 3) ? DIRICHLET : NEUMANN;
}

// u \cdot n
void FluxBoundValue(int bdComp, double t, double &value)
{
  
  switch(bdComp)
  {
    case 0: 
    {
      value = 0.0; // Neumann
      //value = t*t*t/3.0; // Dirichlet
      break;
    }
    case 1: 
    {
      value = (t - t*t*t) / 3.0; // Neumann
      //value = t*t*t/3.0 - t; // Dirichlet
      break;
    }
    case 2: 
    {
      double x = 1-t;
      value = (x*x*x - x) / 3.0; // Neumann
      //value = x - x*x*x/3.0; // Dirichlet
      break;
    }
    case 3:
    {
      double y = 1-t;
      //value = 0.0; // Neumann
      value = -y*y*y/3.0; // Dirichlet
      break;
    }
  default:
    ErrMsg("wrong boundary part number");
    exit(0);
  break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double eps = 1.0 / TDatabase::ParamDB->SIGMA_PERM;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    // RHS for exact solution
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // g
  }
}


