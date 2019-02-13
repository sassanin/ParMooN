// Deformation tensor problem in DEVSS with sine and cosine functions for velocity
// 

void ExampleFile()
{
  
  OutPut("Example: SinCos_DFT.h ") ;
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
void ExactD1(double x, double y, double *values)
{
    values[0] = Pi*cos(Pi*x);
    values[1] = -Pi*Pi*sin(Pi*x);
    values[2] = 0;
    values[3] = -Pi*Pi*Pi*cos(Pi*x);
}

void ExactD2(double x, double y, double *values)
{
    values[0] = 0.5*Pi*Pi*y*sin(Pi*x);
    values[1] = 0.5*Pi*Pi*Pi*y*cos(Pi*x);
    values[2] = 0.5*Pi*Pi*sin(Pi*x);
    values[3] = -0.5*Pi*Pi*Pi*Pi*y*sin(Pi*x);
}

void ExactD3(double x, double y, double *values)
{
    values[0] = -Pi*cos(Pi*x);
    values[1] = Pi*Pi*sin(Pi*x);
    values[2] = 0;
    values[3] = Pi*Pi*Pi*cos(Pi*x);
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_DFT(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;

}

void D1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=Pi*cos(Pi*Param);
            break;
    case 1: value=Pi*cos(Pi);
            break;
    case 2: value=Pi*cos(Pi*(1-Param));
            break;
    case 3: value=Pi;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}

void D2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0.5*Pi*Pi*Param*sin(Pi);
            break;
    case 2: value=0.5*Pi*Pi*sin(Pi*(1-Param));
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}



void D3BoundValue(int BdComp, double Param, double &value)
{
switch(BdComp)
  {
    case 0: value=-Pi*cos(Pi*Param);
            break;
    case 1: value=-Pi*cos(Pi);
            break;
    case 2: value=-Pi*cos(Pi*(1-Param));
            break;
    case 3: value=-Pi;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}


void LinCoeffs_DFT(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y, *param;
  double u1x, u1y, u2x, u2y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    coeff[0] = 1;

    
  }
}
