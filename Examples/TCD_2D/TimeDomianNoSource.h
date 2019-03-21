// ==========================================================================
// instationary problem
// ==========================================================================
// example file
// =========================================================================

//==========================================================================
//a zero field with non-zero Direchlet boundary condition
//of the unit square (0,1)X(0,1)
//==========================================================================
#define __MOVINGMESH__ 
// #define __TIMECONSISTSUPG__
// #define __CONSERVATIVEALE__ 
// #define __CONSISTANTSUPG__


void ExampleFile()
{
  OutPut("Example: TimeDomianNoSource.h" << endl); 
  OutPut("__MOVINGMESH__" << endl); 
  
#ifdef  __CONSERVATIVEALE__  
  TDatabase::ParamDB->P6 = 1;// con-ALE
  OutPut("Conservative ALE form with, that is, with - div w" << endl);   
#else
  TDatabase::ParamDB->P6=0;// non-conservative ALE
  OutPut("Non-Conservative ALE form, that is, witout - div w term" << endl);  
#endif
  
#ifdef  __CONSISTANTSUPG__  
  TDatabase::ParamDB->REACTOR_P28=1;// consist-SUPG
  OutPut("Consistant SUPG ALE form" << endl);  
#else
  TDatabase::ParamDB->REACTOR_P28=0;// in consist-SUPG
  OutPut("IN-Consistant SUPG ALE form" << endl);    
#endif
  
  
  
//   if(TDatabase::ParamDB->P6==1)  
//    {
// 
//    }
//   else 
//    {
// 
//    }
   
    if(TDatabase::ParamDB->P6==1) // con-ALE
   {
    OutPut("Conservative ALE form with, that is, with - div w" << endl); 
   }
  else// non-conservative ALE
   {
    OutPut("Non-Conservative ALE form, that is, witout - div w term" << endl); 
   }
   
   
}
// exact solution
void Exact(double x, double y, double *values)
{    
  values[0]= 0.;
  values[1]= 0; 
  values[2]= 0;
  values[3]= 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
//  cond = DIRICHLET;
  cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{  
 value = 0;
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
 values[0] = 1600.*x*(1.-x)*y*(1.-y);
//   values[0] = 0.;
//  values[0] = 1. - sqrt(x*x + y*y);
}

void ModifyCoords(double x, double y, double &X, double &Y, double t, double &W1, double &W2)
{
//  double t = TDatabase::TimeDB->CURRENTTIME; 
//  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
 double r,d, norm, r1, n1, n2;
 
  d = 2. - cos(20.*Pi*t);
 
 X = x*d;
 Y = y*d;
   
 
 // mesh velocity
 W1 = (  20.*Pi*x*sin(20.*Pi*t))/d;
 W2 = (  20.*Pi*y*sin(20.*Pi*t))/d;
 
 // - sign due to -w\cdot\nabla C in the equation    
//cout << x << " x, y " << y  << endl; 
//  r = sqrt( x*x + y*y);
//  
//  //cout << "r" << r <<  endl; 
//   d = 0.5*r*(1.-r)*sin(8*Pi*t);
//  // the = atan(y/x);
// //cout << " d " << d   <<  endl; 
//   
//   if(r>1e-8)
//   {
//    n1 = x/r;
//    n2 = y/r;
//    
//    norm = sqrt( n1*n1 + n2*n2);
//    n1 /=norm;
//    n2 /=norm;
//    
//     //cout << n1 << " n1, n2 " << n2  << endl;
//    //cout << X << " x    , y " << Y  << endl;
//    X =  x + d * n1;
//    Y =  y + d * n2;
// //   r1 = sqrt( X*X + Y*Y);  
// 
//                
//   }
/*   X =  x ;
   Y =  y ;*/  
 
      
   cout << X << " x, y " << Y  << endl;
  }

void ModifyCoords(double x, double y, double &X, double &Y, double t)
{
//  double t = TDatabase::TimeDB->CURRENTTIME; 
//  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
 double r,d, norm, r1, n1, n2;
 
  d = 2. - cos(20.*Pi*t);
 
 X = x*d;
 Y = y*d;
   
 
 // mesh velocity
//  W1 = (-1.)*20.*Pi*x*sin(20.*Pi*t);
//  W2 = (-1.)*20.*Pi*x*sin(20.*Pi*t);
 
 // - sign due to -w\cdot\nabla C in the equation    
// cout << x << " x, y " << y  << endl; 
//  r = sqrt( x*x + y*y);
//  
//  //cout << "r" << r <<  endl; 
//   d = 0.5*r*(1.-r)*sin(8*Pi*t);
//  // the = atan(y/x);
// //cout << " d " << d   <<  endl; 
//   
//   if(r>1e-8)
//   {
//    n1 = x/r;
//    n2 = y/r;
//    
//    norm = sqrt( n1*n1 + n2*n2);
//    n1 /=norm;
//    n2 /=norm;
//    
//     //cout << n1 << " n1, n2 " << n2  << endl;
//    //cout << X << " x    , y " << Y  << endl;
//    X =  x + d * n1;
//    Y =  y + d * n2;
// //   r1 = sqrt( X*X + Y*Y);  
// 
//                
//   }
/*   X =  x ;
   Y =  y ;*/  
 
      
//    cout << X << " x, y " << Y  << endl;
  }

  
  
void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->PE_NR;
  int i;
  double  *coeff;// *param;
//   double divW, t = TDatabase::TimeDB->CURRENTTIME;
 
//   divW=40*Pi*sin(20*Pi*t)/(2-cos(20*Pi*t));

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    //TCD2D.C, MatrixMARhsAssemble()
    coeff[0] = eps; // diffusion
    coeff[1] = 0.0; //  b1
    coeff[2] = 0.0; //  b2
    coeff[3] = 0.0;  // c 
    coeff[4] = 0.0;  // rhs  
    
  }
}

// kind of boundary condition (for FE space needed)
void GridBoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void GridBoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void GridCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}
 
