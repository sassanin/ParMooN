void ExampleFile()
{
  OutPut("Example: ShearFlow_Oldroyd_DEVSS.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 1-y*y;
  values[1] = 0;
  values[2] = -2*y;
  values[3] = -2.0;
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
  values[0] = -x;
  values[1] = -1.0;
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
    case 0: value=1;
            break;
    case 1: value=1-pow(Param,2);
            break;
    case 2: value=0;
            break;
    case 3: value=1-pow((1-Param),2);
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
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
   double eps = 1/TDatabase::ParamDB->RE_NR, beta = TDatabase::ParamDB->P3, nu = 1/TDatabase::ParamDB->WEI_NR;
  int i;
  double *coeff, y, x, nondim;
  
  if (TDatabase::ParamDB->TENSOR_TYPE == 1)
  nondim = beta*eps;
  else if (TDatabase::ParamDB->TENSOR_TYPE == 2)
  nondim = eps;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    y = Y[i];
    x = X[i];
    coeff[0] = nondim;
    coeff[1] = 2*eps - 1.0; // f1
    coeff[2] = 0; // f2
  }
}

// ========================================================================
// exact solution
// ========================================================================
void ExactS1(double x, double y, double *values)
{
    values[0] = 1.0 + 8*TDatabase::ParamDB->WEI_NR*TDatabase::ParamDB->WEI_NR*pow(y,2);
    values[1] = 0;
    values[2] = 16*TDatabase::ParamDB->WEI_NR*TDatabase::ParamDB->WEI_NR*y;
    values[3] = 16*TDatabase::ParamDB->WEI_NR*TDatabase::ParamDB->WEI_NR;
}

void ExactS2(double x, double y, double *values)
{
    values[0] = -2*y*TDatabase::ParamDB->WEI_NR;
    values[1] = 0;
    values[2] = -2*TDatabase::ParamDB->WEI_NR;
    values[3] = 0;
}


void ExactS3(double x, double y, double *values)
{
  values[0] = 1.0;
  values[1] = 0;
  values[2] = 0;
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
  double Wei = TDatabase::ParamDB->WEI_NR;
  switch(BdComp)
  {
    case 0: value=1.0;
            break;
    case 1: value=1.0 + 8*Wei*Wei*pow(Param,2);
            break;
    case 2: value=1.0 + 8*Wei*Wei;
            break;
    case 3: value=1.0 + 8*Wei*Wei*pow(1-Param,2);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}



void S2BoundValue(int BdComp, double Param, double &value)
{
  double Wei = TDatabase::ParamDB->WEI_NR;
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=-2*Wei*Param;
            break;
    case 2: value=-2*Wei;
            break;
    case 3: value=-2*Wei*(1-Param);
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
    case 0: value=1.0;
            break;
    case 1: value=1.0;
            break;
    case 2: value=1.0;
            break;
    case 3: value=1.0;
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
  double nu=1./TDatabase::ParamDB->WEI_NR;
  int i;
  double *coeff, *param, x, y;
  double u1, u1x, u1y, u2, u2x, u2y, tau1, tau2, tau3, tau1x, tau1y, tau2x, tau2y, tau3x, tau3y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
  
    coeff[0] = nu;
    coeff[1] = nu;          //  f1
    coeff[2] =  0;          //  f2
    coeff[3] =  nu;         //  f3
    
    
    

  }
}

// ========================================================================
// exact solution
// ========================================================================
void ExactD1(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactD2(double x, double y, double *values)
{
    values[0] = -y;
    values[1] = 0;
    values[2] = -1.0;
    values[3] = 0;
}

void ExactD3(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
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
  return;
}

void D2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=-Param;
            break;
    case 2: value=-1.0;
            break;
    case 3: value=-(1-Param);
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

//     param = parameters[i];
//       
//     u1x = param[2];
//     u1y = param[4];
//     u2x =  param[3];
//     u2y = param[5];

    coeff[0] = 1;
//     coeff[1] = u1x;
//     coeff[2] = (u1y+u2x)*0.5;
//     coeff[3] = u2y;
    
  }
}