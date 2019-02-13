// Bechmark for Darcy problem, exact solution is
// 

void ExampleFile()
{
  OutPut("Example: SinSinSolution.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = -2*Pi*cos(2*Pi*x)*sin(2*Pi*y);
  values[1] = 4*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y);
  values[2] = -4*Pi*Pi*cos(2*Pi*x)*cos(2*Pi*y);
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -2*Pi*sin(2*Pi*x)*cos(2*Pi*y);
  values[1] = -4*Pi*Pi*cos(2*Pi*x)*cos(2*Pi*y);
  values[2] = 4*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y);
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = sin(2*Pi*x)*sin(2*Pi*y);
  values[1] = 2*Pi*cos(2*Pi*x)*sin(2*Pi*y);
  values[2] = 2*Pi*sin(2*Pi*x)*cos(2*Pi*y);
  values[3] = 0.;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int bdComp, double t, BoundCond &cond)  
{
  cond = (bdComp == 0 || bdComp == 3) ? NEUMANN : DIRICHLET;
}

// u \cdot n
void FluxBoundValue(int bdComp, double t, double &value)
{
  switch(bdComp)
  {
    case 0: 
    {
      value = 0; // Neumann
      //value = 2*Pi*sin(2*Pi*t); // Dirichlet
      break;
    }
    case 1: 
    {
      //value = 0; // Neumann
      value = -2*Pi*sin(2*Pi*t); // Drichlet
      break;
    }
    case 2: 
    {
      //value = 0; // Neumann
      value = -2*Pi*sin(2*Pi*(1-t)); // Dirichlet
      break;
    }
    case 3:
    {
      value = 0; // Neumann
      //value = 2*Pi*sin(2*Pi*(1-t)); // Dirichlet
      break;
    }
    default: cout << "wrong boundary part number" << endl;
    
    break;
  }
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double eps = 1.0/TDatabase::ParamDB->SIGMA_PERM;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    // RHS for exact solution
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 8*Pi*Pi*sin(2*Pi*X[i])*sin(2*Pi*Y[i]); // g
  }
}


