// ======================================================================
// Sine problem 3D
// ======================================================================
// #include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: HeatChanel.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  if (fabs(z)<1e-8)
  { cond = DIRICHLET; }
  else
  { cond  = NEUMANN; }
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
   double r2;
   
   value = 0;
       
   r2 = (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
   if (fabs(z)<1e-8 && r2<0.25)
   {
      value =  -4.*(0.25-r2);
   }
   
    r2 = (x+0.5)*(x+0.5)+(y+0.5)*(y+0.5);
   if (fabs(z)<1e-8 && r2<0.25)
   {
      value =  -4.*(0.25-r2);
   }
   
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  static double eps=1.9e-5, u=TDatabase::ParamDB->P0;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = u;
    coeff[4] = 0;
    coeff[5] = 0.;
  }
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{ 
 double t = 0;
 double k = 0.1;
  
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;  
}
