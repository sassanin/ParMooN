// ======================================================================
// Sine problem
// ======================================================================

#define __HEATLINE__

void ExampleFile()
{
  OutPut("Example: SineLaplace_hl.h" << endl) ;
}

void GetVelo(double xi, double yi, double*val)
{
  double theta = atan2(yi,xi);
  double r = sqrt(xi*xi + yi*yi); 
 
   if(r<1.)
    {
     val[0] = (1.-r)*sin(theta);
     val[1] = -(1.-r)*cos(theta);
    }
    else
    {
     val[0] = 0.;
     val[1] = 0.;       
    }
     
}

// exact solution
void ExactU1(double x, double y, double *values)
{
  double xi = x - 0.5;
  double yi = y - 0.5;
  //convective velo
    GetVelo(xi, yi, values); 

  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

// exact solution
void ExactU2(double x, double y, double *values)
{
  double xi = x - 0.5;
  double yi = y - 0.5;
  //convective velo
    GetVelo(xi, yi, values);  
    
  values[0] = values[1] ;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}


// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if(BdComp==0 || BdComp==2)
  { cond = NEUMANN;}
  else  
  { cond = DIRICHLET;}
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;

  if(BdComp==3)
    value = 1.;
   else  
    value = 0;

}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;
  int i;
  double *coeff, *param, theta, r, xi, yi;
//   int xmin,xmax,ymin,ymax,x0,y0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    
    xi = x[i] - 0.5;
    yi = y[i] - 0.5;
    
    coeff[0] = eps;
    
    //convective velo
    GetVelo(xi, yi, coeff+1);
    
    coeff[3] = 0;
    coeff[4] = 0;
    
//   cout << xi << " , " << yi << " eps " << eps << " : u " <<   coeff[1] << " v " << coeff[2] << endl;
  }
}

// kind of boundary condition (for FE space needed)
void HeatFuncBoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void HeatFuncBoundValue(int BdComp, double Param, double &value)
{
  static double eps=1./TDatabase::ParamDB->PE_NR;

  if(BdComp==0)
    value = 0;
  else if(BdComp==1)
    value = Param;
  else if(BdComp==2)
    value = 1;
  else
    value = 1.-Param;
}

void HeatfuncCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double r2;
 
 for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = 1.0;
    coeff[1] = 0; 
    coeff[2] = 0; 
    coeff[3] = 0;
    coeff[4] = 0;
  }
}