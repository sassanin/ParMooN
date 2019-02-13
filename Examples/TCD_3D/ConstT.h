// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================
// Chang et al., Numerical Heat Transfer, Part B vol 19, pp. 69-84, 1991
// exact solution for constant wall temperature for parallelepiped
void ExampleFile()
{
  OutPut("Example: ConstT.h" << endl);
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  double t, T_int;
  double PI =3.14159265;
  
//   x = 2.*x-1.;
//   y = 2.*y-1.;
//   z = 2.*z-1.;
  
  T_int =1.;    
  t = TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = (64.*T_int/(PI*PI*PI))*(exp(-3.*PI*PI*t/4.))*cos(PI*x/2.)*cos(PI*y/2.)*cos(PI*z/2.);
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  
}

// void Exact(double x, double y, double z, double *values)
// {
//  int i, j, k, N;
//  double t, a, s, T, T_int, T_w, temp;
//  double m, n, l;
//  double PI =3.14159265;
// 
//   x = 2.*x-1.;
//   y = 2.*y-1.;
//   z = 2.*z-1.;
// 
//  N= 20;
//  t = TDatabase::TimeDB->CURRENTTIME;
// 
//   a = 0.;
//   T_int =1.;
//   T_w = 2.;
// 
//   x = 0.5;
//   y = 0;
//   z = 0;
//   s = T_w;
//   for(i=1;i<N;i++)
//    for(j=1;j<N;j++)
//     for(k=1;k<N;k++)
//      {
//       m = (double)i;
//       n = (double)j;
//       l = (double)k;
//       a =(64.*(T_int - T_w)/(PI*PI*PI*(2.*m-1.)*(2.*n-1.)*(2.*l-1.)))*sin((2.*m-1.)*PI/2.)*sin((2.*n-1.)*PI/2.)*sin((2.*l-1.)*PI/2.);
//       temp = ((2.*m-1.)*PI/2.)*((2.*m-1.)*PI/2.) + ((2.*n-1.)*PI/2.)*((2.*n-1.)*PI/2.) + ((2.*l-1.)*PI/2.)*((2.*l-1.)*PI/2.);
//       s += a*(exp(-temp*t))*cos((2.*m-1.)*PI*x/2.)*cos((2.*n-1.)*PI*y/2.)*cos((2.*l-1.)*PI*z/2.) ;
//      }
// 
//   values[0] = s;
//   values[1] = 0;
//   values[2] = 0;
//   values[3] = 0;
//   values[4] = 0;
// }

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
  double t, T_int;
  double PI =3.14159265;
  
//   x = 2.*x-1.;
//   y = 2.*y-1.;
//   z = 2.*z-1.;

  T_int =1.;    
  t =0.;
  
  values[0] = (64.*T_int/(PI*PI*PI))*(exp(-3.*PI*PI*t/4.))*cos(PI*x/2.)*cos(PI*y/2.)*cos(PI*z/2.);
}

// void InitialCondition(double x, double y, double z, double *values)
// {
//  int i, j, k, N;
//  double t, a, s, T, T_int, T_w, temp;
//  double m, n, l;
//  double PI =3.14159265;
// 
//   x = 2.*x-1.;
//   y = 2.*y-1.;
//   z = 2.*z-1.;
// 
//  N= 20;
//  t = 0.;
// 
//   a = 0.;
//   T_int =1.;
//   T_w = 2.;
// 
//   x = 0.5;
//   y = 0;
//   z = 0;
//   s = T_w;
//   for(i=1;i<N;i++)
//    for(j=1;j<N;j++)
//     for(k=1;k<N;k++)
//      {
//       m = (double)i;
//       n = (double)j;
//       l = (double)k;
//       a =(64.*(T_int - T_w)/(PI*PI*PI*(2.*m-1.)*(2.*n-1.)*(2.*l-1.)))*sin((2.*m-1.)*PI/2.)*sin((2.*n-1.)*PI/2.)*sin((2.*l-1.)*PI/2.);
//       temp = ((2.*m-1.)*PI/2.)*((2.*m-1.)*PI/2.) + ((2.*n-1.)*PI/2.)*((2.*n-1.)*PI/2.) + ((2.*l-1.)*PI/2.)*((2.*l-1.)*PI/2.);
//       s += a*(exp(-temp*t))*cos((2.*m-1.)*PI*x/2.)*cos((2.*n-1.)*PI*y/2.)*cos((2.*l-1.)*PI*z/2.) ;
//      }
// 
//   values[0] = s;
// }
// 
// 


// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
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
    coeff[5] = 0;
    coeff[6] = 0;
  }
}

