// Navier-Stokes problem, Driven cavity
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
   OutPut("Example: DrivenCavity3D.h " << endl);
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
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
  double eps = 1e-8;
     
    if(TDatabase::ParamDB->MESH_TYPE==1){
      switch(CompID)
  {
            
    case 0: 
            value=0;
            break;
    case 5: 
      
    if (fabs(y)==1)
	{
	  if ((fabs(x)>eps)&&(fabs(1-x)>eps)&&(fabs(z)>eps)&&(fabs(1-z)>eps))
	    value = 1.0;
	  else
	    value = 0.0;
	}
	else
    value =0.0 ;
//        value =1.0 ;
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
    case 1: 
            value=0;
            break;
    default: 
      cout << "wrong boundary part number" << endl;
            break;
  }  
    }
    else{
        if (fabs(y)==1)
	{
	  if ((fabs(x)>eps)&&(fabs(1-x)>eps)&&(fabs(z)>eps)&&(fabs(1-z)>eps))
	    value = 1.0;
	  else
	    value = 0.0;
	}
	else
    value =0.0 ;}  
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
