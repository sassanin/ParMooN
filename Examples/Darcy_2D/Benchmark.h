// Bechmark for Darcy problem, exact solution is
// 
// u(x,y) = (y,x)
// p(x,y) = 0

void ExampleFile()
{
  OutPut("Example: Benchmark.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = y;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = x;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int bdComp, double t, BoundCond &cond)  
{
  cond = (bdComp == 0) ? NEUMANN : DIRICHLET;
}

// u \cdot n
void UNBoundValue(int BdComp, double t, double &value)
{
  switch(BdComp)
  {
    case 0: 
      //value = 0; // Dirichlet
      value = 0; // Neumann
      break;
    case 1: 
      value = t; // Dirichlet
      //value = 0; // Neumann
      break;
    case 2: 
      value = 1 - t; // Dirichlet
      //value = 0; // Neumann
      break;
    case 3: 
      value = t - 1; // Dirichlet
      //value = 0; // Neumann
      break;
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
  const double eps = 1./TDatabase::ParamDB->SIGMA_PERM;
  
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    // RHS for exact solution
    coeffs[i][1] = Y[i]; // f1
    coeffs[i][2] = X[i]; // f2
    coeffs[i][3] = 0.0;  // g
  }
}


