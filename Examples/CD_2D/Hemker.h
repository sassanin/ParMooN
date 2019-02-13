// ======================================================================
// Example from P.W. Hemker
// ======================================================================
#define __HEMKER__

void ExampleFile()
{
  OutPut("Example: Hemker.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void BoundCondition(int BdComp, double t, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 1:
	case 2:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 1:
      value = 0;
      break;
    case 4:
      value = 1;
      break;
    default:
      value = 0;
  }
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  values[0] = 0;
}

void BoundConditionAdjoint(int BdComp, double t, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 2:
	case 3:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValueAdjoint(int BdComp, double Param, double &value)
{
    value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double angle = 0, v1, v2;
  int i;
  double *coeff;

  v1 = cos(angle);
  v2 = sin(angle);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = v1;
    coeff[2] = v2;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

void ComputeExtremalValues(int N, double *sol, double  *values)
{
   int i;
   double max, min;

   min = 1e10;
   max = -1e10;
   
   for(i=0;i<N;i++)
   {
      if(sol[i]-1 > max)
         max = sol[i]-1;
      if(sol[i] < min)
         min = sol[i];
   }

   values[0] = min;
   values[1] = max;
}
 
/** compute curve of the outflow boundary */
void ComputeOutflowBoundary(int level, TFEFunction2D *ufct)
{
  double h, x=4,values[3],y;
  int i, bound_points = 401;
  h = 6.0/(bound_points-1);
  for (i=0;i<bound_points; i++)
  {
      y = -3+i*h;
      ufct->FindGradient(x,y,values);
      OutPut("cutline " << x << " " <<  y << 
            " " <<  values[0] << endl);
  }
}

// computation of some global errors, only for P1 or Q1 !!!
//
// values[0] : absolute value of largest negative undershoot in 
//             a circle around the cylinder
// values[1] : difference of largest positive value in a circle
//             around the cylinder and 1
// values[2] : absolute value of largest negative undershoot 
//             for x > 2
// values[3] : difference of largest positive value for x>2 and 1
// 
void ComputeLocalExtrema(TFEFunction2D *ufct, double *values)
{
  TBaseCell *cell;
  TCollection *Coll;
  TFESpace2D *FESpace2D;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  int N_BaseFunct;
  double xi, eta;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, u;
  double *Values, x, y, val;
  int *GlobalNumbers, *BeginIndex, N_Cells, N_Edges;
  int i, j, k, *Numbers;
  double extr[4];

  extr[0] = -1;
  extr[1] = -1;
  extr[2] = -1;
  extr[3] = 0;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // loop over all edges
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    //double diam = cell->GetDiameter();
    
    FE_ID = FESpace2D->GetFE2D(i, cell);
    FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
    RefTrans = FE_Obj->GetRefTransID();

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    N_BaseFunct = bf->GetDimension();
    
    uorig = new double[N_BaseFunct];
    uxorig = new double[N_BaseFunct];
    uyorig = new double[N_BaseFunct];
    
    uref = new double[N_BaseFunct];
    uxiref = new double[N_BaseFunct];
    uetaref = new double[N_BaseFunct];
    
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    for (j=0;j<N_Edges;j++)
    {
      // compute coordinates
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      if (x<-1.5)
	  continue;
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
                uref, uxiref, uetaref, uorig, uxorig, uyorig);
      // compute value of fe function at (x,y)
      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      //OutPut(x << " " << y << " " << u << endl);
      // strip with the circle
      if ((x>-1.5)&&(x<1.5))
      {
	  if ((u<=0)&&(fabs(u)>extr[0]))
	      extr[0] = fabs(u);
	  if ((u>=1)&&(fabs(u-1)>extr[1]))
	      extr[1] = fabs(u-1);
      }
      if (x>2)
      {
	  if ((u<=0)&&(fabs(u)>extr[2]))
	      extr[2] = fabs(u);
	  if ((u>=1)&&(fabs(u-1)>extr[3]))
	      extr[3] = fabs(u-1);
      }
    } // endfor (j) N_Edges
  
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
    
  values[0] = extr[0];
  values[1] = extr[1];
  values[2] = extr[2];
  values[3] = extr[3];
 }

/****************************************************************/
//
// for FEM_TVD
//
/****************************************************************/

void CheckWrongNeumannNodes(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
  const int max_entries = 50000;  
  int i, j, min_val, type;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-6, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;

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
      vertex->GetCoords(x[j], y[j]);
      if ((fabs(x[j]+3)<eps)//||(fabs(y[j]+3)<eps)||(fabs(y[j]-3)<eps)
	  || (fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps))
      {
	   boundary_vertices[j] = 1;
	   found++;
      }
    }
    // no cell with edge with vertex on the boundary
    if (found<2) 
	continue;
       // finite element on the mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    //int N_ = TFEDatabase2D::GetN_BaseFunctFromFE2D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_2D_T_A:
	case C_Q1_2D_Q_A:
	case C_Q1_2D_Q_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_2D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if (j<2){
			    tmp_diri[diri_counter] = dof[j];
			}
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    else
				tmp_diri[diri_counter] = dof[2];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		    }
		    // inflow x = -3
		    if (fabs(x[j]+3)<eps) 
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    }
		    // circle
		    if (fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) 
		    {
			tmp_bdry[diri_counter] = 4;
			// parameter does not matter, since b.c. equal to 1
			tmp_param[diri_counter] =  0;
		    }
		    diri_counter++;
		}
	    }
	    break;
	// P_2, Q_2
	case C_P2_2D_T_A:
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[4];
                       tmp_diri[diri_counter+2] = dof[5];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[5];
                       tmp_diri[diri_counter+2] = dof[8];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[5];
                       tmp_diri[diri_counter+1] = dof[3];
                       tmp_diri[diri_counter+2] = dof[0];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[8];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[6];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[6];
                     tmp_diri[diri_counter+1] = dof[3];
                     tmp_diri[diri_counter+2] = dof[0];
                   break;

                }
              
		if (diri_counter+2 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}

		// inflow x = -3
		if ((fabs(x[j]+3)<eps)&&(fabs(x[(j+1)%N_V]+3)<eps)) 
		{
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    tmp_param[diri_counter+2] = (-y[(j+1)%N_V]+3)/6.0;
		    tmp_param[diri_counter+1] = (tmp_param[diri_counter] +  tmp_param[diri_counter+2])/2.0;
		}
		// circle
		if ((fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) && 
		    (fabs(sqrt(x[(j+1)%N_V]*x[(j+1)%N_V]+y[(j+1)%N_V]*y[(j+1)%N_V])-1)<eps))
		{
		    tmp_bdry[diri_counter] = 4;
		    tmp_bdry[diri_counter+1] = 4;
		    tmp_bdry[diri_counter+2] = 4;
		    // parameter does not matter, since b.c. equal to 1
		    tmp_param[diri_counter] = 0;
		    tmp_param[diri_counter+1] = 0;
		    tmp_param[diri_counter+2] = 0;
		}
		diri_counter +=3;
	      }
	    }
	    break;
	// P_3, Q_3
	case C_P3_2D_T_A:
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;

               // P3: local dof 0, 1, 2, 3 are on the boundary
               // Q3: local dof 0, 1, 2, 3 are on the boundary
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
		     tmp_diri[diri_counter+3] = dof[3];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[6];
                       tmp_diri[diri_counter+2] = dof[8];
		       tmp_diri[diri_counter+3] = dof[9];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[11];
		       tmp_diri[diri_counter+3] = dof[15];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[9];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[4];
                       tmp_diri[diri_counter+3] = dof[0];
		     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[15];
                       tmp_diri[diri_counter+1] = dof[14];
                       tmp_diri[diri_counter+2] = dof[13];
			tmp_diri[diri_counter+3] = dof[12];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[12];
                     tmp_diri[diri_counter+1] = dof[8];
                     tmp_diri[diri_counter+2] = dof[4];
		     tmp_diri[diri_counter+3] = dof[0];
                   break;
                }
              
		if (diri_counter+3 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}

		// inflow x = -3
		if ((fabs(x[j]+3)<eps)&&(fabs(x[(j+1)%N_V]+3)<eps)) 
		{
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_bdry[diri_counter+3] = 3;
		    tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    tmp_param[diri_counter+3] = (-y[(j+1)%N_V]+3)/6.0;
		    tmp_param[diri_counter+1] = (2*tmp_param[diri_counter] +  tmp_param[diri_counter+3])/3.0;
		    tmp_param[diri_counter+2] = (tmp_param[diri_counter] +  2*tmp_param[diri_counter+3])/2.0;
		}
		// circle
		if ((fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) && 
		    (fabs(sqrt(x[(j+1)%N_V]*x[(j+1)%N_V]+y[(j+1)%N_V]*y[(j+1)%N_V])-1)<eps))
		{
		    tmp_bdry[diri_counter] = 4;
		    tmp_bdry[diri_counter+1] = 4;
		    tmp_bdry[diri_counter+2] = 4;
		    tmp_bdry[diri_counter+3] = 4;
		    // parameter does not matter, since b.c. equal to 1
		    tmp_param[diri_counter] = 0;
		    tmp_param[diri_counter+1] = 0;
		    tmp_param[diri_counter+2] = 0;
		    tmp_param[diri_counter+3] = 0;
		}
		diri_counter +=4;
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

  //OutPut("CheckNeumannNodesForVelocity: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary numbers
  neum_to_diri_bdry = new int[diri_counter_1];
  // allocate array for the corresponding boundary parameters
  neum_to_diri_param = new double[diri_counter_1];
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
      neum_to_diri_bdry[i] = tmp_bdry[found];
      neum_to_diri_param[i] = tmp_param[found];
      tmp_diri[found] = -1;
  }
/*
  for (i=0;i<diri_counter_1;i++)
  {
      OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_bdry[i] <<
	     " " << neum_to_diri_param[i] <<  endl);
  }
*/
}

void SetDirichletNodesFromNeumannNodes(TSquareMatrix2D **SQMATRICES, 
				       double *rhs,
				       int N_neum_to_diri,
				       int *neum_to_diri,
				       int *neum_to_diri_bdry,
				       double *neum_to_diri_param)
{
    TSquareMatrix2D *MatrixA;
    double *Entries_A;
    int i, l, l0, l1, index, *RowPtr_A, *KCol_A;
    
    MatrixA = SQMATRICES[0];
    RowPtr_A      = MatrixA->GetRowPtr();
    KCol_A        = MatrixA->GetKCol();
    Entries_A     = MatrixA->GetEntries();
    // loop over dof to change
    for (i=0;i<N_neum_to_diri;i++)
    {
	index = neum_to_diri[i];
	l0 = RowPtr_A[index];
	l1 = RowPtr_A[index+1];
	for (l=l0;l<l1;l++)
	{
	    // diagonal entry
	    if (KCol_A[l]==index)
		Entries_A[l] = 1;  
	    else
		Entries_A[l] = 0;
	}
	// set boundary condition
	BoundValue(neum_to_diri_bdry[i], neum_to_diri_param[i], rhs[index]);
    }
}

void ComputeDifferenceToCoarseLevel(TCollection *Coll_fine,
				    TCollection *Coll_coarse,
				    TFEFunction2D *u_fine, 
				    TFEFunction2D *u_coarse)
{
    int i, j, k, N_Cells, N_Edges, coarse_no;
    double x, y, x_c, y_c, val_fine[4], val_coarse[4], c1err = -1, c1err_coarse = -1;
    double x_err, y_err, x_err_c, y_err_c;
    TBaseCell *cell, *parent;
    
    // number of cells
    N_Cells = Coll_fine->GetN_Cells();
    
    // loop over all edges
    for(i=0;i<N_Cells;i++)
    {
	// cell
	cell = Coll_fine->GetCell(i);
	// parent cell
	parent = cell->GetParent();
	coarse_no = Coll_coarse->GetIndex(parent);
	//OutPut(coarse_no << " ");
	// number of edges
	N_Edges=cell->GetN_Edges();
	for (j=0;j<N_Edges;j++)
	{
	    cell->GetVertex(j)->GetCoords(x, y);
	    u_fine->FindGradientLocal(cell, i, x, y, val_fine);
	    u_coarse->FindGradientLocal(parent, coarse_no, x, y, val_coarse);
	    if (fabs(val_fine[0] - val_coarse[0]) > c1err)
	    {
		c1err = fabs(val_fine[0] - val_coarse[0]);
		x_err = x;
		y_err = y;
	    }
	    for (k=0;k<N_Edges;k++)
	    {
		parent->GetVertex(k)->GetCoords(x_c, y_c);
		if ((fabs(x_c -x ) < 1e-6) && (fabs(y_c -y ) < 1e-6))
		{
		    if (fabs(val_fine[0] - val_coarse[0]) > c1err_coarse)
		    {
			c1err_coarse = fabs(val_fine[0] - val_coarse[0]);
			x_err_c = x;
			y_err_c = y;
		    }
		}
	    }
	}
    } 

    //OutPut("C1 error f " << c1err  << " \\& " << x_err <<"," << y_err << endl);
    //OutPut(" C1 error c " << c1err_coarse  << " \\& " << x_err_c <<"," << y_err_c << endl);

    OutPut("C1 error f "<<" & " << c1err  << " &  ( " << x_err <<"," << y_err << ")"<< "\\\\\\hline" << endl);
    OutPut("C1 error c "<< " & "<< c1err_coarse  << " &  ( " << x_err_c <<"," << y_err_c << ")" << "\\\\\\hline"<< endl);
}
void ComputeCutLines_X(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues, y0=0;
  int i, j, N_Cells, bound_points = 10001;
    TBaseCell *cell;
  
  h = 3.0/(bound_points-1);

  cutvalues = new double[6*bound_points];
  memset(cutvalues , 0 , 6*bound_points*SizeOfDouble);

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    x = -1;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	cutvalues[j] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    x = 0;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
    x = 1;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+3*bound_points] = val[0];
	}
    }
    x = 2;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+4*bound_points] = val[0];
	}
    }
    x = 4;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+5*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cutx"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points]  << " " << cutvalues[j+3*bound_points]  << " " 
	     << cutvalues[j+4*bound_points]  << " " << cutvalues[j+5*bound_points] << endl);
  }
}

void ComputeCutLines_Y(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues;
  int i, j, N_Cells, bound_points = 20001;
    TBaseCell *cell;
  
  h = 10.0/(bound_points-1);

  cutvalues = new double[3*bound_points];
  memset(cutvalues , 0 , 3*bound_points*SizeOfDouble);
  
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    y = 0;
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	cutvalues[j] = x;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    y = 1;
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    /*if (cutvalues[j+2*bound_points]!=4711)
	    {
		OutPut("belegt"<<endl);
		exit(1);
	    }
	    else*/
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cuty"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points] << endl);
  }
}

void ComputeCutLines_epsY(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues;
  int i, j, N_Cells, bound_points = 20001;
    TBaseCell *cell;
    double eps = 1.0/TDatabase::ParamDB->RE_NR;
  
  h = 10.0/(bound_points-1);

  cutvalues = new double[3*bound_points];
  memset(cutvalues , 0 , 3*bound_points*SizeOfDouble);

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    y = 1-sqrt(eps);
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	cutvalues[j] = x;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    y = 1+sqrt(eps);
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cutyeps"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points] << endl);
  }
}
void ComputeCutLines_eps_radial(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues, tmp, r;
  int i, j, N_Cells, bound_points = 10001;
    TBaseCell *cell;
    double eps = 1.0/TDatabase::ParamDB->RE_NR;
  
  h = 2*Pi/(bound_points-1);

  cutvalues = new double[10*bound_points];
  memset(cutvalues , 0 , 10*bound_points*SizeOfDouble);

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    r = (1+eps)*(1+eps);
    for (j=0;j<bound_points;j++)
    {
	x = r*cos(h * j-Pi);
	y = r*sin(h*j-Pi);
	cutvalues[j] = h*j -Pi;
	cutvalues[j+bound_points] = x;
	cutvalues[j+2*bound_points] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+3*bound_points] = val[0];
	}
    }
    tmp=pow(eps,2.0/3.0);
    r = (1+tmp)*(1+tmp);
    for (j=0;j<bound_points;j++)
    {
	x = r*cos(h * j-Pi);
	y = r*sin(h*j-Pi);
	cutvalues[j+4*bound_points] = x;
	cutvalues[j+5*bound_points] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+6*bound_points] = val[0];
	}
    }
    tmp=sqrt(eps);
    r = (1+tmp)*(1+tmp);
    for (j=0;j<bound_points;j++)
    {
	x = r*cos(h * j-Pi);
	y = r*sin(h*j-Pi);
	cutvalues[j+7*bound_points] = x;
	cutvalues[j+8*bound_points] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+9*bound_points] = val[0];
	}
    }

  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cutradeps"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points] << " " << cutvalues[j+3*bound_points] 
	     << " " << cutvalues[j+4*bound_points] << " " << cutvalues[j+5*bound_points] 
	     << " " << cutvalues[j+6*bound_points] << " " << cutvalues[j+7*bound_points] 
	     << " " << cutvalues[j+8*bound_points] << " " << cutvalues[j+9*bound_points]  << endl);
  }
}
