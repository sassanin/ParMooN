// Stationary Conformation-Stress tensor problem with sine and cosine functions
// 

void ExampleFile()
{
  
  OutPut("Example: SinCos_CST_Giesekus.h ") ;
}



// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = sin(Pi*x);
    values[1] = Pi*cos(Pi*x);
    values[2] = 0;
    values[3] = -Pi*Pi*sin(Pi*x);
}

void ExactU2(double x, double y, double *values)
{
    values[0] = -Pi*y*cos(Pi*x);
    values[1] = Pi*Pi*y*sin(Pi*x);
    values[2] = -Pi*cos(Pi*x);
    values[3] = Pi*Pi*Pi*y*cos(Pi*x);
}

void ExactP(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*cos(Pi*y);
  values[1] = Pi*cos(Pi*x)*cos(Pi*y);
  values[2] = -Pi*sin(Pi*x)*sin(Pi*y);
  values[3] = -Pi*Pi*sin(Pi*x)*cos(Pi*y)-Pi*Pi*sin(Pi*x)*cos(Pi*y);
}

void InitialU1(double x, double y, double *values)
{
  values[0] = sin(Pi*x);
}

void InitialU2(double x, double y, double *values)
{
  values[0] = -Pi*y*cos(Pi*x);
}

void InitialP(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*cos(Pi*y);
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=sin(Pi*Param);
            break;
    case 1: value=sin(Pi);
            break;
    case 2: value=sin(Pi*(1-Param));
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=Pi*Param;
            break;
    case 2: value=-Pi*cos(Pi*(1-Param));
            break;
    case 3: value=-Pi*(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double nu=1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;
  double u1, u1x, u1y, u1lap, u2, u2x, u2y, u2lap, px, py;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    coeff[0] = nu;
    // prescribed solution
    u1 =  sin(Pi*x);
    u1x = Pi*cos(Pi*x);
    u1y = 0;
    u1lap =  -Pi*Pi*sin(Pi*x);
    u2 =  -Pi*y*cos(Pi*x);
    u2x =  Pi*Pi*y*sin(Pi*x);
    u2y = -Pi*cos(Pi*x);
    u2lap = Pi*Pi*Pi*y*cos(Pi*x);
    px =  Pi*cos(Pi*x)*cos(Pi*y);
    py = -Pi*sin(Pi*x)*sin(Pi*y);

    coeff[1] = -nu*u1lap+u1*u1x+u2*u1y+px;
    coeff[2] = -nu*u2lap+u1*u2x+u2*u2y+py;

    coeff[3] = u1;
    coeff[4] = u2;
    /*coeff[5] = px;
    coeff[6] = py;
    coeff[7] = u1y;
    coeff[8] = u2y;*/
  }
}

// ========================================================================
// exact solution
// ========================================================================
void ExactS1(double x, double y, double *values)
{
    values[0] = sin(Pi*x);
    values[1] = Pi*cos(Pi*x);
    values[2] = 0;
    values[3] = -Pi*Pi*sin(Pi*x);
}

void ExactS2(double x, double y, double *values)
{
    values[0] = -Pi*y*cos(Pi*x);
    values[1] = Pi*Pi*y*sin(Pi*x);
    values[2] = -Pi*cos(Pi*x);
    values[3] = Pi*Pi*Pi*y*cos(Pi*x);
}



void ExactS3(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*cos(Pi*y);
  values[1] = Pi*cos(Pi*x)*cos(Pi*y);
  values[2] = -Pi*sin(Pi*x)*sin(Pi*y);
  values[3] = -Pi*Pi*sin(Pi*x)*cos(Pi*y)-Pi*Pi*sin(Pi*x)*cos(Pi*y);
}

void InitialS1(double x, double y, double *values)
{
  values[0] = sin(Pi*x);
}

void InitialS2(double x, double y, double *values)
{
  values[0] = -Pi*y*cos(Pi*x);
}

void InitialS3(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*cos(Pi*y);
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
    case 0: value=sin(Pi*Param);
            break;
    case 1: value=0;
            break;
    case 2: value=sin(Pi*(1-Param));
            break;
    case 3: value=0;
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
    case 1: value=Pi*Param;
            break;
    case 2: value=-Pi*cos(Pi*(1-Param));
            break;
    case 3: value=-Pi*(1-Param);
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
    case 0: value=sin(Pi*Param);
            break;
    case 1: value=0;
            break;
    case 2: value=-sin(Pi*(1-Param));
            break;
    case 3: value=0;
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
  double nu=1./TDatabase::ParamDB->WEI_NR, alpha = TDatabase::ParamDB->P2;
  int i;
  double *coeff, x, y;
  double u1, u1x, u1y, u2, u2x, u2y, tau1, tau2, tau3, tau1x, tau1y, tau2x, tau2y, tau3x, tau3y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];


    // prescribed solution
    u1 =  sin(Pi*x);
    u1x = Pi*cos(Pi*x);
    u1y = 0;
    u2 =  -Pi*y*cos(Pi*x);
    u2x =  Pi*Pi*y*sin(Pi*x);
    u2y = -Pi*cos(Pi*x);
    tau1 =  sin(Pi*x);
    tau1x = Pi*cos(Pi*x);
    tau1y = 0;
    tau2 =  -Pi*y*cos(Pi*x);
    tau2x =  Pi*Pi*y*sin(Pi*x);
    tau2y = -Pi*cos(Pi*x);

    tau3 =   sin(Pi*x)*cos(Pi*y);
    tau3x =  Pi*cos(Pi*x)*cos(Pi*y);
    tau3y = -Pi*sin(Pi*x)*sin(Pi*y);
    
    coeff[0] = nu;
    coeff[1] = (u1*tau1x) + (u2*tau1y)  + ((1-(2.0*alpha))*nu*(tau1)) - (2*(tau1*u1x + tau2*u1y)) + (alpha*nu*((tau1*tau1) + (tau2*tau2)));
      
    coeff[2] =  (u1*tau2x) + (u2*tau2y)  + ((1-(2.0*alpha))*nu*(tau2))  - (tau1*u2x + tau3*u1y) + (alpha*nu*((tau1*tau2) + (tau2*tau3)));
    
    coeff[3] =  (u1*tau3x) + (u2*tau3y)  + ((1-(2.0*alpha))*nu*(tau3)) - (2*(tau2*u2x + tau3*u2y)) + (alpha*nu*((tau2*tau2) + (tau3*tau3)));
    
    coeff[4] = alpha;

    

  }
}
