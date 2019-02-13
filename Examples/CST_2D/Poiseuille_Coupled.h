 
// u(x,y) = (4*y*(1-y), 0)
// p(x,y) = x-1/2

void ExampleFile()
{
  OutPut("Example: Poiseuille_Coupled.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 4*y*(1-y);
  values[1] = 0;
  values[2] = 4-8*y;
  values[3] = -8;
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
  values[0] = x-0.5;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
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
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=4*Param*(1-Param);
            break;
    case 2: value=0;
            break;
    case 3: value=4*Param*(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}



void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 1+8*eps; // f1
    coeff[2] = 0; // f2
  }
}

// ========================================================================
// exact solution
// ========================================================================
void ExactS1(double x, double y, double *values)
{
    values[0] = x+3.0;
    values[1] = 1;
    values[2] = 0;
    values[3] = 0;
}



void ExactS2(double x, double y, double *values)
{
    values[0] = y;
    values[1] = 0;
    values[2] = 1;
    values[3] = 0;
}



void ExactS3(double x, double y, double *values)
{
  values[0] = (x*y) + 5.0;
  values[1] = y;
  values[2] = x;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_CST(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;

}

void S1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=Param + 3.0;
            break;
    case 1: value=1 + 3.0;
            break;
    case 2: value=(1-Param) + 3.0;
            break;
    case 3: value=0 + 3.0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}



void S2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=Param;
            break;
    case 2: value=1;
            break;
    case 3: value=(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}




void S3BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0 + 5.0;
            break;
    case 1: value=Param + 5.0;
            break;
    case 2: value=(1-Param) + 5.0;
            break;
    case 3: value=0 + 5.0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs_CST(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double nu=1./TDatabase::ParamDB->RE_NR, alpha = TDatabase::ParamDB->P2;
  int i;
  double *coeff, *param, x, y;
  double u1, u1x, u1y, u2, u2x, u2y, tau1, tau2, tau3, tau1x, tau1y, tau2x, tau2y, tau3x, tau3y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    param = parameters[i];
    
  u1 = param[0];                 
  u2 = param[1];                 
  u1x = param[2];
  u2x = param[3];
  u1y = param[4];
  u2y = param[5];
    
    tau1 =  x + 3.0;
    tau1x = 1;
    tau1y = 0;   
    tau2 =  y;
    tau2x =  0;
    tau2y = 1;
    tau3 =   (x*y) + 5.0;
    tau3x =  y;
    tau3y = x;
    
    coeff[0] = nu;
    coeff[1] = (u1*tau1x) + (u2*tau1y)  + ((1-(2.0*alpha))*nu*(tau1)) - (2*(tau1*u1x + tau2*u2x)) + (alpha*nu*((tau1*tau1) + (tau2*tau2)));
      
    coeff[2] =  (u1*tau2x) + (u2*tau2y)  + ((1-(2.0*alpha))*nu*(tau2))  - (tau1*u1y + tau3*u2x) + (alpha*nu*((tau1*tau2) + (tau2*tau3)));
    
    coeff[3] =  (u1*tau3x) + (u2*tau3y)  + ((1-(2.0*alpha))*nu*(tau3)) - (2*(tau2*u1y + tau3*u2y)) + (alpha*nu*((tau2*tau2) + (tau3*tau3)));
    
    coeff[4] = alpha;
    

  }
}