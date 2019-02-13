// Bechmark for Darcy problem, exact solution of the flux is linear 
// 

void ExampleFile()
{
  OutPut("Example: Obstacle.h.\nThe obstacle is 10 times less conductive " <<
         "and has the shape of\n ");
  switch((int)TDatabase::ParamDB->P0)
  {
    case 0:
      OutPut("a circle in the center\n");
      break;
    case 1:
      OutPut("a square in the center\n");
      break;
    case 2:
      OutPut("a diamond in the center\n");
      break;
    case 3:
      OutPut("two boxes on one vertical line, one at the top, the other "<<             
             "at the bottom\n");
      break;
    case 4:
      OutPut("two half circles on one vertical line, one at the top, " <<
             "one at the bottom\n");
      break;
    case 5:
      OutPut("two boxes with a vertical offset, on at the top, one at " <<
             "the bottom\n");
      break;
    case 6:
      OutPut("two half ellipsoids with a vertical offset, on at the top, " <<
             "one at the bottom\n");
      break;
    case 7:
      OutPut("hump in the center (continuous change in permeability)\n");
      break;
    default:
      ErrMsg("Unknown obstacle. Set P0 in the input file to 0,...,7\n");
      exit(0);
      break;
  }
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

// ========================================================================
// multi indices used for various things
// ========================================================================
// no exact solution known
void ExactU1(double x, double y, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int bdComp, double t, BoundCond &cond)  
{
  cond = DIRICHLET;
}

// u \cdot n
void FluxBoundValue(int bdComp, double Param, double &value)
{
  switch(bdComp)
  {
    case 0: 
      value = 0.0;
      break;
    case 1: 
      value = 1.0;
      break;
    case 2: 
      value = 0.0;
      break;
    case 3:
      value = -1.0;
      break;
    default: cout << "wrong boundary part number" << endl;
      break;
  }
}

// ========================================================================
// coefficient and right hand side
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double eps = 1.0/TDatabase::ParamDB->SIGMA_PERM;
  double *coeff;
  double x,y;
  double r;
  
  // choose different obstacles:
  // 0 - circle in center
  // 1 - square in center
  // 2 - diamond in center
  // 3 - two boxes, one at the top, one at the bottom
  // 4 - two half circles, one at the top, one at the bottom
  // 5 - two boxes, one at the top, one at the bottom, boxes have an offset
  // 6 - two half ellipsoids with an offset, one at the top, one at the bottom,
  // 7 - continuous change in permeability in form of a hump
  int obstacle = TDatabase::ParamDB->P0;
  // inside the obstacle the porousity sigma is multiplied by factor
  double factor = 10;

  for(int i = 0; i < n_points; i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    
    coeff[0] = eps;
    switch (obstacle)
    {
      case 0:
        r = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
        if( r < 0.2 )
          coeff[0] *= factor;
        break;
      case 1:
        if( x>0.45 && x<0.55 && y>0.35 && y<0.65)
          coeff[0] *= factor;
        break;
      case 2:
        if (x+y>0.8 && x+y<1.2 && x-y<0.2 && x-y>-0.2)
          coeff[0] *= factor;
        break;
      case 3:
        if(x>0.35 && x<0.65 && (y<0.2 || y>0.8))
          coeff[0] *= factor;
        break;
      case 4:
        r = sqrt((x-0.5)*(x-0.5)+y*y);
        if( r < 0.2 )
          coeff[0] *= factor;
        r = sqrt((x-0.5)*(x-0.5)+(y-1)*(y-1));
        if( r < 0.2 )
          coeff[0] *= factor;
        break;
      case 5:
        if((x>0.65 && x<0.75 && y<0.6) || (x>0.25 && x<0.35 && y>0.4))
          coeff[0] *= factor;
        break;
      case 6:
        r = sqrt((x-0.75)*(x-0.75)+0.1*y*y);
        if( r < 0.2 )
          coeff[0] *= factor;
        r = sqrt((x-0.25)*(x-0.25)+0.1*(y-1)*(y-1));
        if( r < 0.2 )
          coeff[0] *= factor;
        break;
      case 7:
        r = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
          if( r < 0.2 )
            coeff[0] *= 1+0.5*factor*(cos(5*Pi*r)+1);
        break;
    }
    
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // g
  }
}
