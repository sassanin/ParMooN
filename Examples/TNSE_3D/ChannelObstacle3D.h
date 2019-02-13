// Navier-Stokes problem, 3D Channel flow
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
#ifdef _MPI
 int rank;
 MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

 if(rank==TDatabase::ParamDB->Par_P0)
#endif
  { 
   OutPut("Example: ChannelObstacle3D.h " << endl);
  }

}

void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{ 
   TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  
  if(TDatabase::ParamDB->MESH_TYPE==0)
  {
    if(x==1)
       cond = NEUMANN;
    else if
      (x==0)
      cond = DIRICHLET;
    else if
      (y==0)
      cond = DIRICHLET;
    else if
      (y==1)
      cond = DIRICHLET;
    else if
      (z==0)
      cond = DIRICHLET;
    else if
      (z==1)
      cond = DIRICHLET;
    
 }
 else
 {
     switch(CompID)
  {
      // inflow through the face (x==0) and outflow at (x==1) 
      // in ChannelObstacle3D.mesh (x==0)is the face tagged 1000 and 
      // in the face (x==1) tagged as 1002
    
    case 0: case 1: case 4: case 3: case 5:case 6:case 7:case 8:case 9:
      cond = DIRICHLET;
       break;
      
    case 2:  
     cond = NEUMANN;
      break;
      
    default: 
      cout << "wrong boundary part number" << endl;
            break;     
   
  } 
 }
 
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{


 if(TDatabase::ParamDB->MESH_TYPE==0)
 {
    if(x==0)
      value=y*(1-y)*z*(1-z);
//       value=1;
    else if
      (x==1)
      value=0;
    else if
      (y==0)
      value=0;
    else if
      (y==1)
      value=0;
    else if
      (z==0)
      value=0;
    else if
      (z==1)
     value=0;
 }
 else
 {
//    cout<<"here in drivencavity.h"<<endl;
   switch(CompID)
  {
            
    case 0: 
            value=((6-y)*(6+y)*z*(10-z))/900;
           //value=y*(1-y)*z*(1-z);
            break;
    case 1: 
            value=0;
            break;
    case 2: 
            value=0;
            break;
    case 3: 
            value=0;
            break;
    case 4: 
            value=0;
            break;
    case 5: 
            value=0;
            break;
    case 6: 
            value=0;
            break;
    case 7: 
            value=0;
            break;
    case 8: 
            value=0;
            break;
    case 9: 
            value=0;
            break;
    default: 
      cout << "wrong boundary part number" << endl;
            break;
  }  
}
 
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // f3
  }
}
