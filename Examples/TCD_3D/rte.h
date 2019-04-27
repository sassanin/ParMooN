// #define __UREA__
// #define __SIMPATURS__
//   
// #include <Urea_3d4d.h>
// #include <MacroCell.h>

void ExampleFile()
{
  
  OutPut("Example: amc.h" << endl); //OutPut("GRID_TYPE set to " << TDatabase::ParamDB->GRID_TYPE << endl);
}

// ========================================================================
// definitions for the temperature
// ========================================================================

void Exact( double x, double y, double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
  
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{ 
 double t = 0;
 double k = 0.1;
  
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
 /* values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);*/  
}

void BoundCondition(int dummy,double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int dummy,double x, double y, double z, double &value)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
  
  value = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);  
}

// ========================================================================
// BilinearCoeffs for Heat 
// ========================================================================
void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
//   double eps= 1.0/TDatabase::ParamDB->RE_NR;
  double eps= 1.0;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  b[0] = 0;
  b[1] = 0;
  b[2] = 0;  
  c = 0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = b[0];
    // convection in y direction
    coeff[2] = b[1];
    // convection in z direction
    coeff[3] = b[2];
    // reaction
    coeff[4] = c;
     // rhs
    coeff[5] = (3.*eps*Pi*Pi - 0.1)*(exp(-0.1*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
    coeff[6] = 0;

  }
}
