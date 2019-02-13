// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================
// exact solution in unit cube
void ExampleFile()
{
  OutPut("Example: ConstTSmooth.h" << endl);
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  double t, PI =3.14159265;   
  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = (exp(-0.1*t))*sin(PI*x)*cos(PI*y)*cos(PI*z);
  values[1] = PI*(exp(-0.1*t))*cos(PI*x)*cos(PI*y)*cos(PI*z);
  values[2] = -PI*(exp(-0.1*t))*cos(PI*x)*sin(PI*y)*cos(PI*z);
  values[3] = -PI*(exp(-0.1*t))*cos(PI*x)*cos(PI*y)*sin(PI*z);
  values[4] = -3.*PI*PI*(exp(-0.1*t))*sin(PI*x)*cos(PI*y)*cos(PI*z); ;  
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
  double t, PI =3.14159265;

  t =0.;  
  values[0] = (exp(-0.1*t))*sin(PI*x)*cos(PI*y)*cos(PI*z);
}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  double t, PI =3.14159265;
  t = TDatabase::TimeDB->CURRENTTIME;
  
  value = (exp(-0.1*t))*sin(PI*x)*cos(PI*y)*cos(PI*z);
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, PI =3.14159265;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  
  if(TDatabase::ParamDB->RE_NR!=0)
   eps = 1.0/TDatabase::ParamDB->RE_NR;
  else
   eps = 0;
  
//   cout << "eps  eps = " << eps << endl;
  b[0] = 0;
  b[1] = 0;
  b[2] = 0;  
  c = 0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

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
    coeff[5] = (3.*eps*PI*PI - 0.1)*(exp(-0.1*t))*sin(PI*x)*cos(PI*y)*cos(PI*z); // f
    coeff[6] = 0;
    
  }
}

