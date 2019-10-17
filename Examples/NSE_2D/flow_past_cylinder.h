// Navier-Stokes problem, Driven cavity
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  OutPut("Example: flow_Cylinder.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
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
void BoundCondition(int i, double t, BoundCond &cond)
{
   // INLET = 0
   // OUTLET = 1
   // walls = 2
   // cylinder = 3
    
    if (i == 1)
        cond = NEUMANN;
    else
        cond = DIRICHLET;
   
}

void U1BoundValue(int BdComp, double Param, double &value)
{
       // INLET = 0 // OUTLET = 1 // walls = 2 // cylinder = 3
    if (Param > 0 ) cout << "Param Exists" << endl;
  switch(BdComp)
  {
    case 0: value = 1 ; //value = -(Param)*(Param-1)*4 ;
            break;
    case 1: value = 0;
            break;
    case 2: value =0;
            break;
    case 3: value = 0;
  
            break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>3) cout << "wrong boundary part number" << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}

