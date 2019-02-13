// Bechmark for Darcy problem, exact solution is known
// see Arif Masud, Thomas J.R. Hughes: "A stabilized mixed finite element
// method for Darcy flow", chapter 4.2
//
// Problem: The source term consists of two point sources at (0,0) and (1,1). 
// These are Dirac delta functions. To model these we use equivalent 
// prescribed flows at the first edges touching the corner points (0,0) and 
// (1,1). However for this to work we need to know the edge length and the  
// order of the used finite elements. That's why one has to modify the 
// function "FluxBoundValue" every time the grid is refined..
//
// Problem 2: On all other parts of the boundary we prescribe u.n=0 for 
// symmetry reasons. However the solution given in ExactP, ExactU1 and 
// ExactU2 does not fulfill this. To achieve this, we would have to continue 
// the source/sink pattern until infinity. 

void ExampleFile()
{
  OutPut("Example: 5SpotProblem.h\n");
  OutPut(" This works so far only on one level with a regular grid\n");
}

// ========================================================================
// exact solution
// ========================================================================
inline double r(double x,double y,double x0,double y0)
{  return (x-x0)*(x-x0) + (y-y0)*(y-y0); }
void ExactU1(double x, double y, double *values)
{
  double r0,r1,r2,r3,r4;
  r0 = r(x,y, 0, 0);
  r1 = r(x,y, 1, 1);
  r2 = r(x,y,-1, 1);
  r3 = r(x,y,-1,-1);
  r4 = r(x,y, 1,-1);
  values[0] = (x/r0 - (x-1)/r1 - (x+1)/r2 - (x+1)/r3 - (x-1)/r4)/4;
  values[1] = 0;
  values[2] = -(x*y/(r0*r0) - (x-1)*(y-1)/(r1*r1) - (x+1)*(y-1)/(r2*r2) 
                           - (x+1)*(y+1)/(r3*r3) - (x-1)*(y+1)/(r4*r4) )/2;
  values[3] = 0;
}
void ExactU2(double x, double y, double *values)
{
  double r0,r1,r2,r3,r4;
  r0 = r(x,y, 0, 0);
  r1 = r(x,y, 1, 1);
  r2 = r(x,y,-1, 1);
  r3 = r(x,y,-1,-1);
  r4 = r(x,y, 1,-1);
  values[0] = (y/r0 - (y-1)/r1 - (y-1)/r2 - (y+1)/r3 - (y+1)/r4)/4;;
  values[1] = -(x*y/(r0*r0) - (x-1)*(y-1)/(r1*r1) - (x+1)*(y-1)/(r2*r2) 
                           - (x+1)*(y+1)/(r3*r3) - (x-1)*(y+1)/(r4*r4) )/2;
  values[2] = 0;
  values[3] = 0;
}
void ExactP(double x, double y, double *values)
{
  
  int n_plus = 1;
  int n_minus = 4;
  double plus[] = {0,0, 
                   2,0, 2,2, 0,2, -2,2, -2,0, -2,-2, 0,-2, 2,-2,
                   4,0, 4,2, 4,4, 2,4, 0,4, -2,4, -4,4, -4,2, -4,0, -4,-2, -4,-4, -2,-4, 0,-4, 2,-4, 4,-4, 4,-2};
  double minus[] = {1,1, -1,1, -1,-1, 1,-1,
                    3,1, 3,3, 1,3, -1,3, -3,3, -3,1, -3,-1, -3,-3, -1,-3, 1,-3, 3,-3, 3,-1,
                    5,1, 5,3, 5,5, 3,5, 1,5, -1,5, -3,5, -5,5, -5,3, -5,1, -5,-1, -5,-3, -5,-5, -3,-5, -1,-5, 1,-5, 3,-5, 5,-5, 5,-3, 5,-1};
  values[0] = 0;
  for(int i=0; i<n_plus; i++)
    values[0] -= log(sqrt(r(x,y,plus[2*i],plus[2*i+1])));
  for(int i=0; i<n_minus; i++)
    values[0] += log(sqrt(r(x,y,minus[2*i],minus[2*i+1])));
  values[0] *= 0.25;
  
  values[1] = 0;
  for(int i=0; i<n_plus; i++)
    values[1] -= (x-plus[2*i])/r(x,y,plus[2*i],plus[2*i+1]);
  for(int i=0; i<n_minus; i++)
    values[1] += (x-minus[2*i])/r(x,y,minus[2*i],minus[2*i+1]);
  values[1] *= 0.25;
  
  values[2] = 0;
  for(int i=0; i<n_plus; i++)
    values[2] -= (y-plus[2*i+1])/r(x,y,plus[2*i],plus[2*i+1]);
  for(int i=0; i<n_minus; i++)
    values[2] += (y-minus[2*i+1])/r(x,y,minus[2*i],minus[2*i+1]);
  values[2] *= 0.25;
  /*
  double r0,r1,r2,r3,r4;
  r0 = r(x,y, 0, 0);
  r1 = r(x,y, 1, 1);
  r2 = r(x,y,-1, 1);
  r3 = r(x,y,-1,-1);
  r4 = r(x,y, 1,-1);
  values[0] = -(log(sqrt(r0))-log(sqrt(r1))-log(sqrt(r2))
                            -log(sqrt(r3))-log(sqrt(r4)) )/4;
  values[1] = -(x/r0 - (x-1)/r1 - (x+1)/r2 - (x+1)/r3 - (x-1)/r4)/4;
  values[2] = -(y/r0 - (y-1)/r1 - (y-1)/r2 - (y+1)/r3 - (y+1)/r4)/4;
  */
  values[3] = 0;
}



// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int bdComp, double t, BoundCond &cond)  
{ cond = DIRICHLET; }
// u \cdot n
void FluxBoundValue(int bdComp, double Param, double &value)
{
  value = 0;
  double h = 1.0/32.0;
  int order = 1;
  switch(bdComp)
  {
    case 0: 
      if (Param<h)
      {
        value=-(1-Param/h)/(4*h);
        if( order == 2)
          value = -(3/(8*h))*(1-2*Param/h + Param*Param/(h*h));
      }
      break;
    case 1:
      if(Param>1-h)
      {
        value=(1-(1-Param)/h)/(4*h);
        if( order == 2)
          value = (3/(8*h))*(1-2*(1-Param)/h + (1-Param)*(1-Param)/(h*h));
      }
      break;
    case 2:
      if (Param<h)
      {
        value=(1-Param/h)/(4*h);
        if( order == 2)
          value = (3/(8*h))*(1-2*Param/h + Param*Param/(h*h));
      }
      break;
    case 3:
      if(Param>1-h)
      {
        value=-(1-(1-Param)/h)/(4*h);
        if( order == 2)
          value = -(3/(8*h))*(1-2*(1-Param)/h + (1-Param)*(1-Param)/(h*h));
      }
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
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    if ( (X[i] < 0.5 && Y[i]<0.5) || (X[i] > 0.5 && Y[i]>0.5) )
      coeff[0] = 2.0;
    else 
      coeff[0] = 200.0;
    // RHS for exact solution
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // g
  }
}


