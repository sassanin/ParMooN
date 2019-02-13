// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================
#define __SIN4__

void ExampleFile()
{
  OutPut("Example: Sin4.h" << endl);
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  double t;

    t = TDatabase::TimeDB->CURRENTTIME;
    values[0] = sin(t)*(sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)+1);
    values[1] = sin(t)*2*Pi*(cos(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z));
    values[2] = sin(t)*2*Pi*(sin(2*Pi*x)*cos(2*Pi*y)*sin(2*Pi*z));
    values[3] = sin(t)*2*Pi*(sin(2*Pi*x)*sin(2*Pi*y)*cos(2*Pi*z));
    values[4] = sin(t)*4*Pi*Pi*(-sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)
				-sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)
				-sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z));
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
    values[0] = 0;
}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  value = sin(t);
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  b[0] = 1;
  b[1] = 1;
  b[2] = 1;  
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
    coeff[5] = cos(t)*(sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)+1)
	-coeff[0] * (sin(t)*4*Pi*Pi*(-sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)
			-sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)
			-sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)))
	+coeff[1] * sin(t)*2*Pi*(cos(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z))
	+coeff[2] * sin(t)*2*Pi*(sin(2*Pi*x)*cos(2*Pi*y)*sin(2*Pi*z))
	+coeff[3] * sin(t)*2*Pi*(sin(2*Pi*x)*sin(2*Pi*y)*cos(2*Pi*z))
	+coeff[4] * sin(t)*(sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z)+1);
    // rhs from previous time step
    coeff[6] = 0;
  }
}

/****************************************************************/
//
// for FEM_TVD
//
/****************************************************************/
void CheckWrongNeumannNodes(TCollection *Coll, TFESpace3D *fespace,
int &N_neum_to_diri, int* &neum_to_diri,
			    double* &neum_to_diri_x, double* &neum_to_diri_y, double* &neum_to_diri_z)
{
   const int max_entries = 50000;  
  int i, j, N_, min_val;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[8], tmp_diri[max_entries];
  double x[8], y[8], z[8], eps = 1e-6, tmp_x[max_entries], tmp_y[max_entries], tmp_z[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  
  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();

  diri_counter = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    found = 0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      boundary_vertices[j] = 0;
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      // vertex on the upper lid
      if ((fabs(x[j])<eps)||(fabs(y[j])<eps)||(fabs(x[j]-1)<eps)||(fabs(y[j]-1)<eps)
	  ||(fabs(z[j])<eps)||(fabs(z[j]-1)<eps))
      {
	  // Dirichlet boundary
	  boundary_vertices[j] = 1;
	  found++;
      }
    }
    
    // no cell with face with vertex on the boundary
    if (found<3) 
	continue;
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase3D::GetN_BaseFunctFromFE3D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_3D_T_A:
	case C_Q1_3D_H_A:
	case C_Q1_3D_H_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_3D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if ((j==0)||(j==1)||(j==4)||(j==5))
			{
			    tmp_diri[diri_counter] = dof[j];
			}
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    if (j==3)
				tmp_diri[diri_counter] = dof[2];
			    if (j==6)
				tmp_diri[diri_counter] = dof[7];
			    if (j==7)
				tmp_diri[diri_counter] = dof[6];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		    }
		    tmp_x[diri_counter] = x[j];
		    tmp_y[diri_counter] = y[j];
		    tmp_z[diri_counter] = z[j];
		    diri_counter++;
		      OutPut(j<< " " << tmp_x[diri_counter-1] << " "  << x[j] << endl);
		}
	    }
	    break;
	default:
	    OutPut("CheckNeumannNodesForVelocity not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckNeumannNodesForVelocity, just delete the exit" << endl);
	    exit(4711);
    }	    
  }
    
  // condense
  for (i=0;i<diri_counter;i++)
  {
      if (tmp_diri[i] == -1)
	  continue;
      diri_counter_1++;
      for (j=i+1;j<diri_counter;j++)
      {
	  if (tmp_diri[i] == tmp_diri[j])
	  {
	      tmp_diri[j] = -1;
	  }
      }
  }
  
  OutPut("CheckNeumannNodesForVelocity: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary coordinates
  neum_to_diri_x = new double[diri_counter_1];
  neum_to_diri_y = new double[diri_counter_1];
  neum_to_diri_z = new double[diri_counter_1];
  // fill array and sort
  for (i=0;i<diri_counter_1;i++)
  {
      min_val = tmp_diri[0];
      found = 0;
      for (j=1;j<diri_counter;j++)
      {
	  if ((tmp_diri[j]>0) && ((tmp_diri[j] < min_val) || 
				  (min_val == -1)))
	  {
	      min_val =  tmp_diri[j];
	      found = j;
	  }
      }
      neum_to_diri[i] = tmp_diri[found];
      neum_to_diri_x[i] = tmp_x[found];
      neum_to_diri_y[i] = tmp_y[found];
      neum_to_diri_z[i] = tmp_z[found];
      tmp_diri[found] = -1;
  }
  
  for (i=0;i<diri_counter_1;i++)
  {
      OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_x[i] <<
	     " " << neum_to_diri_y[i] << " " << neum_to_diri_z[i] <<  endl);
  }
}

  
