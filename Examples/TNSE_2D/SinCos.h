// Time-dependent NSE problem with sine and cosine functions
// 

void ExampleFile()
{
  
  OutPut("Example: SinCos.h "<<endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = exp(-t)*sin(Pi*x);

}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = exp(-t)*-Pi*y*cos(Pi*x);

}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = exp(-t)*sin(Pi*x)*cos(Pi*y);

}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    double t=TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);

    values[0] = sin(Pi*x)*t1;
    values[1] = Pi*cos(Pi*x)*t1;
    values[2] = 0;
    values[3] = -Pi*Pi*sin(Pi*x)*t1;
}

void ExactU2(double x, double y, double *values)
{
    double t=TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);

    values[0] = -Pi*y*cos(Pi*x)*t1;
    values[1] = Pi*Pi*y*sin(Pi*x)*t1;
    values[2] = -Pi*cos(Pi*x)*t1;
    values[3] = Pi*Pi*Pi*y*cos(Pi*x)*t1;
}

void ExactP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);

  values[0] = sin(Pi*x)*cos(Pi*y)*t1;
  values[1] = Pi*cos(Pi*x)*cos(Pi*y)*t1;
  values[2] = -Pi*sin(Pi*x)*sin(Pi*y)*t1;
  values[3] = (-Pi*Pi*sin(Pi*x)*cos(Pi*y)-Pi*Pi*sin(Pi*x)*cos(Pi*y))*t1;
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
  double t=TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);

  switch(BdComp)
  {
    case 0: value=sin(Pi*Param)*t1;
            break;
    case 1: value=0;
            break;
    case 2: value=sin(Pi*(1-Param))*t1;
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
  double t=TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);

  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=Pi*Param*t1;
            break;
    case 2: value=-Pi*cos(Pi*(1-Param))*t1;
            break;
    case 3: value=-Pi*(1-Param)*t1;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}




void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double nu = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;
  double u1, u1x, u1y, u2, u2x, u2y, u1t, u2t; 
  double u1lap, u2lap, px, py;
  double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t), t2 = -exp(-t);
   
   
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    // prescribed solution
    u1 =  sin(Pi*x)*t1;
    u1t = sin(Pi*x)*t2;
    u1x = Pi*cos(Pi*x)*t1;
    u1y = 0;
    u1lap =  -Pi*Pi*sin(Pi*x)*t1;
    
    u2 =  -Pi*y*cos(Pi*x)*t1;
    u2t = -Pi*y*cos(Pi*x)*t2;
    u2x =  Pi*Pi*y*sin(Pi*x)*t1;
    u2y = -Pi*cos(Pi*x)*t1;
    u2lap = Pi*Pi*Pi*y*cos(Pi*x)*t1;
    
    px =  Pi*cos(Pi*x)*cos(Pi*y)*t1;
    py = -Pi*sin(Pi*x)*sin(Pi*y)*t1;
    
    
    coeff[0] = nu; // D(u):D(v) term 
    coeff[1] = -nu*u1lap + px + u1t + u1*u1x + u2*u1y;
    coeff[2] = -nu*u2lap + py + u2t + u1*u2x + u2*u2y;

    
    
  }
}

