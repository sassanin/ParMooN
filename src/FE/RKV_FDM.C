/** ==========================================================================
#    This file is part of the finite element software ParMooN.
# 
#    ParMooN (cmg.cds.iisc.ac.in/parmoon) is a free finite element software  
#    developed by the research groups of Prof. Sashikumaar Ganesan (IISc, Bangalore),
#    Prof. Volker John (WIAS Berlin) and Prof. Gunar Matthies (TU-Dresden):
#
#    ParMooN is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    ParMooN is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
#
#    If your company is selling a software using ParMooN, please consider 
#    the option to obtain a commercial license for a fee. Please send 
#    corresponding requests to sashi@iisc.ac.in

# =========================================================================*/ 
   
// =======================================================================
// @(#)RKV_FDM.C
//
// Purpose: contains routines for the RKV temporal discretization in
//          connection with the FDM spatial discretization
//
// Author: Volker John
//
// History: start of implementation 18.05.2010
//
// =======================================================================
#include <Database.h>
#include <MainUtilities.h>
#include <LinAlg.h>
#include <RKV_FDM.h>
#include <FEFunction2D.h>

#include <stdlib.h>
// #include <malloc.h>
#include <string.h>
#include <sstream>

/******************************************************************************/
// Compute_DOF_Conversion_Q1_FDM_UnitSquare2D
// computes the conversion of the dof from Q1 finite element space
// to a FD grid on the unit square
// dof_conversion[k] constains the corresponding fem dof
// works only for equi-distant grids
/******************************************************************************/
#ifdef __2D__
void Compute_DOF_Conversion_Q1_FDM_UnitSquare2D(TCollection *coll,
TFEFunction2D *fesol, int N_x, int N_y, double h,
int* &dof_conversion, double* &x_coord, double* &y_coord)
{
  int i, j, k, dof, N_Cells, N_U, N_V, N2;
  int *begin_index, *global_numbers, *globdof;
  double eps=1e-6, x[4], y[4];
  TFESpace2D *fespace;
  TBaseCell *cell;
  TVertex *vertex;

  N2 = (N_x+1)*(N_y+1);
  // uniform grid
  if(TDatabase::ParamDB->P9!=4711.0)
  {
    x_coord = new double[N2];
    y_coord = new double[N2];
  }
  else
  {
    // non-uniform grid
    i = (int)(TDatabase::ParamDB->INTERNAL_VERTEX_X[0]+1e-6);
    x_coord = new double[i];
    y_coord = new double[i];
  }

  for (i=0 ; i<N_x+1 ; i++ )
  {
    for (j=0; j<N_y+1; j++)
    {
      // compute x and y
      x_coord[i*(N_x+1)+j] = j*h;
      y_coord[i*(N_x+1)+j] = i*h;
    }
  }
  
  // extract information about the fe space and the global d.o.f.
  fespace = fesol->GetFESpace2D();
  // number of dof
  N_U = fespace->GetN_DegreesOfFreedom();
  // allocate array that stores information on conversion
  dof_conversion = new int[N_U];
  // change enumeration of coordinates
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // number of mesh cells
  N_Cells = coll->GetN_Cells();

  // loop over the mesh cells of the global grid
  // find conversion of numbering
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 4)
    {
      OutPut("FDM only for quadrilateral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    globdof = global_numbers+begin_index[i];
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      // compute global dof at this vertex
      switch(j)
      {
        case 0:
          dof = globdof[0];
          break;
        case 1:
          dof = globdof[1];
          break;
        case 2:
          dof = globdof[3];
          break;
        case 3:
          dof = globdof[2];
          break;
      }
      // this is not elegant
      for (k=0;k<N2;k++)
      {
        if ((fabs(x[j]-x_coord[k])<eps)&&(fabs(y[j]-y_coord[k])<eps))
        {
          dof_conversion[k] = dof;
          break;
        }
      }
    }
  }
}

void Compute_DOF_Conversion_Q1_FDM_Rectangle2D(TCollection *coll,
TFEFunction2D *fesol, int &N_x, int &N_y,      
int* &dof_conversion, double* &x_coord, double* &y_coord)
{
  int i, j, k, dof, N_Cells, N_U, N_V, N2, N_1, N_2, found;
  int *begin_index, *global_numbers, *globdof, max_lines = 10000;
  double eps=1e-6, x[4], y[4], *x_tmp, *y_tmp, val;
  TFESpace2D *fespace;
  TBaseCell *cell;
  TVertex *vertex;

  // compute coordinates of the lines
  x_tmp = new double[max_lines];
  y_tmp = new double[max_lines];
  // initialize
  for (i=0;i<max_lines;i++)
    x_tmp[i] =  y_tmp[i] = -4711.0;
  
  N_1 = N_2 = 0;
  // number of mesh cells
  N_Cells = coll->GetN_Cells();
  // loop over the meshcells
  for(i=0;i<N_Cells;i++)
  {
    // get cell
    cell = coll->GetCell(i);
    // number of vertices
    N_V = cell->GetN_Vertices();
    if (N_V < 4)
    {
      OutPut("FDM only for quadrilateral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    // loop over the vertices
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      // check x-coordinate 
      found = 0;
      for (k=0;k<N_1;k++)
      {
	// coordinate already found
	if (fabs(x[j] - x_tmp[k]) < eps)
	{
	  found = 1;
	  break;
	}
      }
      // new coordinate
      if (!found)
      {
	// extend in list
	x_tmp[N_1] = x[j];
	N_1++;
      }
      if (N_1 > max_lines)
      {
	OutPut("Too many lines, increase max_lines !!!" << endl);
	exit(4711);
      }
      // check y-coordinate 
      found = 0;
      for (k=0;k<N_2;k++)
      {
	// coordinate already found
	if (fabs(y[j] - y_tmp[k]) < eps)
	{
	  found = 1;
	  break;
	}
      }
      // new coordinate
      if (!found)
      {
	// extend in list
	y_tmp[N_2] = y[j];
	N_2++;
      }      
      if (N_2 > max_lines)
      {
	OutPut("Too many lines, increase max_lines !!!" << endl);
	exit(4711);
      }
    }
  }
  // all coordinates found, sort from lowest to largest
  for(i=0;i<N_1-1;i++)
    for(j=i+1;j<N_1;j++)
      if(x_tmp[i]>x_tmp[j])
      {
        val = x_tmp[i];
        x_tmp[i] = x_tmp[j];
        x_tmp[j] = val;
      }
  for(i=0;i<N_2-1;i++)
    for(j=i+1;j<N_2;j++)
      if(y_tmp[i]>y_tmp[j])
      {
        val = y_tmp[i];
        y_tmp[i] = y_tmp[j];
        y_tmp[j] = val;
      }

  N_x = N_1 - 1; 
  N_y = N_2 - 1;   
  OutPut("rectangle ["<< x_tmp[0] <<"," << x_tmp[N_1-1] << "] x [" 
         << y_tmp[0] <<"," << y_tmp[N_2-1] << "] with " 
         << N_x << " x " << N_y << " intervals"<<endl);

  N2 = (N_x+1)*(N_y+1);
  if(TDatabase::ParamDB->P9!=4711.0)
  {
    x_coord = new double[N2];
    y_coord = new double[N2];
  }
  else
  {
    i = (int)(TDatabase::ParamDB->INTERNAL_VERTEX_X[0]+1e-6);
    x_coord = new double[i];
    y_coord = new double[i];
  }

  for (i=0 ; i<N_x+1 ; i++ )
  {
    for (j=0; j<N_y+1; j++)
    {
      // compute x and y
      x_coord[i*(N_x+1)+j] = x_tmp[j];
      y_coord[i*(N_x+1)+j] = y_tmp[i];
    }
  }
  // free memory
  delete[] x_tmp;
  delete[] y_tmp;
  
  // extract information about the fe space and the global d.o.f.
  fespace = fesol->GetFESpace2D();
  // number of dof
  N_U = fespace->GetN_DegreesOfFreedom();
  // allocate array that stores information on conversion
  dof_conversion = new int[N_U];
  // change enumeration of coordinates
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();

  // loop over the mesh cells of the global grid
  // find conversion of numbering
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    globdof = global_numbers+begin_index[i];
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      // compute global dof at this vertex
      switch(j)
      {
        case 0:
          dof = globdof[0];
          break;
        case 1:
          dof = globdof[1];
          break;
        case 2:
          dof = globdof[3];
          break;
        case 3:
          dof = globdof[2];
          break;
      }
      // this is not elegant
      for (k=0;k<N2;k++)
      {
        if ((fabs(x[j]-x_coord[k])<eps)&&(fabs(y[j]-y_coord[k])<eps))
        {
          dof_conversion[k] = dof;
          break;
        }
      }
    }
  }
}

/******************************************************************************/
// CheckAdaptiveCriterion
//
// adaptive algorithm for FDM-ENO scheme, the grid will be adaptively changed
//
// input:  sol       - solution on the current grid
//         N_x       - number of intervals in x direction of current grid
//         N_y       - number of intervals in y direction of current grid
//         N_Unknowns- number of unknowns, current grid
//         x_coord   - x coordinates of current grid
//         y_coord   - y coordinates of current grid
//
// output: sol       - solution on the new grid
//         N_x       - number of intervals in x direction of new grid
//         N_y       - number of intervals in y direction of new grid
//         N_Unknowns- number of unknowns, new grid
//         x_coord   - x coordinates of new grid
//         y_coord   - y coordinates of new grid
/******************************************************************************/
void CheckAdaptiveCriterion(double *sol, int &N_x, int &N_y,  int &N_Unknowns,
double *x_coord, double *y_coord, double *aux_adaptive, double mesh_size)
{
  const int max_n = 2048;
  int i, j, index, index1, index2, index3, x_steps, y_steps, steps, size, keep;
  int Nx1, Ny1;
  double val, maximal, x[max_n], y[max_n], h[max_n];
  double  max_interval = 0.1, min_interval = 10e-4, lower = 1e-3, upper = 1e-1;
  double *sol_curr, *sol_curr_x, *x_curr_x, *y_curr_x;

  /*val = TDatabase::TimeDB->CURRENTTIME/TDatabase::TimeDB->TIMESTEPLENGTH;
  if (mesh_size/(1+val) > min_interval)
      min_interval = mesh_size/(1+val);
  if ((1+val)*mesh_size < max_interval)
  max_interval = (1+val)*mesh_size;*/

  //if (TDatabase::TimeDB->CURRENTTIME < 0.1)
  // min_interval = 1e-3;

  // OutPut(mesh_size/val << " min interval " << 1.0/min_interval << " max interval " << 1.0/max_interval << endl);
  Nx1 = N_x+1;
  Ny1 = N_y+1;
  size = (int)(TDatabase::ParamDB->INTERNAL_VERTEX_X[0]);
  sol_curr = aux_adaptive;
  memcpy(sol_curr, sol, Nx1*Ny1*SizeOfDouble);
  sol_curr_x = aux_adaptive + size;
  x_curr_x = sol_curr_x + size;
  y_curr_x = x_curr_x + size;

  // number of already taken lines
  x_steps = -1;
  // for x=0
  keep = 1;
  // loop over the x - values
  for(i=1;i<=N_x;i++)
  {
    // compute maximal difference of function values in the interval
    maximal = 0;
    for (j=0;j<=N_y;j++)
    {
      index = j*Nx1 + i;
      val = fabs(sol_curr[index]-sol_curr[index-1]);
      if (val > maximal)
        maximal = val;
    }
    // remove coordinate i-1
    // it is keep==0
    if ((maximal < lower)&&( x_coord[i] - x[x_steps] <  max_interval)&&(i>1)&&(!keep))
    {
      continue;
    }
    // insert coordinate
    if ((maximal >= upper)&&((x_coord[i] -  x_coord[i-1])/2.0 > min_interval))
    {
      x_steps++;
      x[x_steps] = x_coord[i-1];
      x_steps++;
      x[x_steps] = val = (x_coord[i] +  x_coord[i-1])/2.0;
      keep = 1;
      continue;
    }
    // keep coordinate
    x_steps++;
    x[x_steps] = x_coord[i-1];
    keep = 0;
  }
  // last coordinate
  x[x_steps+1] = x_coord[N_x];
  x_steps++;
  //for (i=0;i<=x_steps;i++)
  //OutPut(i << " " << x[i] << endl);

  GradeMesh(x_steps, x, h);
  //for (i=0;i<=x_steps;i++)
  //OutPut(i << " " << x[i] << endl);

  steps = 0;
  // solution on new grid
  // loop over new x coordinates
  for (i=0;i<x_steps+1;i++)
  {
    while(1)
    {
      // take values
      // same x coordinate
      if (fabs(x[i]-x_coord[steps])<1e-6)
      {
        for (j=0;j<=N_y;j++)
        {
          index = j*Nx1 + steps;
          index1 = j*(x_steps+1) + i;
          sol_curr_x[index1] = sol_curr[index];
          x_curr_x[index1] = x_coord[index];
          y_curr_x[index1] = y_coord[index];
        }
        // increase counter for old x coordinates
        steps++;
        break;
      }
      // interpolate
      if (x[i] < x_coord[steps])
      {
        for (j=0;j<=N_y;j++)
        {
          index = j* Nx1 + steps-1;
          // left point
          index1 = j*(x_steps+1) + i;
          // right point, y-coordinate stays the same
          index2 = index++;
          sol_curr_x[index1] = (sol_curr[index]+sol_curr[index2])/2.0;
          x_curr_x[index1] = (x_coord[index]+ x_coord[index2])/2.0;
          y_curr_x[index1] = y_coord[index];
        }
        break;
      }
      // check next old x coordinate
      steps++;
    }
    // otherwise, coordinate was deleted
  }
  // set new number of intervals in x direction
  N_x = x_steps;

  // number of already taken lines
  y_steps = -1;
  // keep y=0
  keep = 1;
  // loop over the y - values
  for(i=1;i<=N_y;i++)
  {
    // compute maximal difference of function values in the interval
    maximal = 0;
    for (j=0;j<=x_steps;j++)
    {
      index = i*(x_steps+1) + j;
      val = fabs(sol_curr_x[index]-sol_curr_x[index-x_steps-1]);
      if (val > maximal)
        maximal = val;
    }
    //OutPut("y " << y_curr_x[i*(x_steps+1)] << " max " << maximal << endl);
    // remove coordinate i-1
    index1 = i*(x_steps+1);
    if ((maximal < lower)&&( y_curr_x[index1] - y[y_steps] <  max_interval)&&(i>1)&&(!keep))
    {
      continue;
    }
    // insert coordinate
    index = index1 - (x_steps+1);
    if ((maximal >= upper)&& ((y_curr_x[index1] - y_curr_x[index])/2.0> min_interval))
    {
      y_steps++;
      y[y_steps] = y_curr_x[index];
      y_steps++;
      y[y_steps] = (y_curr_x[index1] + y_curr_x[index])/2.0;
      keep = 1;
      continue;
    }
    // keep coordinate
    y_steps++;
    y[y_steps] = y_curr_x[index];
    keep = 0;
  }
  // last coordinate
  y[y_steps+1] = y_curr_x[(x_steps+1)*(N_y+1)-1];
  y_steps++;
  //  for (i=0;i<=y_steps;i++)
  //    OutPut(i << " " << y[i] << endl);
  GradeMesh(y_steps, y, h);
  //  for (i=0;i<=y_steps;i++)
  //    OutPut(i << " " << y[i] << endl);
  //exit(1);
  // new number of y-intervals
  N_y = y_steps;
  // new number of unknowns
  N_Unknowns = (N_x+1)*(N_y+1);
  OutPut("new grid "<< N_x+1 << " x " << N_y+1 << " entries " <<   N_Unknowns  << endl);
  if (N_Unknowns > TDatabase::ParamDB->INTERNAL_VERTEX_X[0])
  {
    OutPut("grid too fine, not enough memory !!!" << endl);
    exit(4711);
  }
  steps = 0;
  // solution on new grid
  // loop over the new y coordinates
  for (i=0;i<y_steps+1;i++)
  {
    while(1)
    {
      // left coordinate of y level on old grid
      index2 = steps*(x_steps+1);
      index3 = i*(x_steps+1);
      // same coordinate take values
      if (fabs(y[i]-y_curr_x[index2])<1e-6)
      {
        for (j=0;j<x_steps+1;j++)
        {
          sol_curr[index3+j] = sol_curr_x[index2+j];
        }
        steps++;
        break;
      }
      // interpolate
      if (y[i] < y_curr_x[index2])
      {
        for (j=0;j<x_steps+1;j++)
        {
          index = index2 + j;
          index1 = index3 + j;
          sol_curr[index1] = (sol_curr_x[index]+sol_curr_x[index-x_steps-1])/2.0;
        }
        break;
      }
      steps++;
    }
    // otherwise, coordinate was deleted
  }

  // copy solution to array sol
  memcpy(sol, sol_curr, N_Unknowns*SizeOfDouble);
  // fill arrays with new coordinates
  index = 0;
  for (i=0;i<=N_y;i++)
  {
    for (j=0;j<=N_x;j++)
    {
      x_coord[index] = x[j];
      y_coord[index] = y[i];
      index++;
    }
  }
}
#endif

/******************************************************************************/
// GradeMesh
//
// change mesh such that ratio of two successive intervals (in both directions)
// is not more than 2
//
/******************************************************************************/
void GradeMesh(int &N_x, double *x, double *h)
{
  int i, steps, changed = 1;

  while (changed)
  {
    changed = 0;
    h[0] = x[0];
    h[1] = x[1];
    steps = 1;
    for (i=1;i<N_x;i++)
    {
      // next interval too large
      if ((x[i+1]-x[i])/(x[i]-x[i-1]) > 2.01)
      {
        // insert new coordinate
        steps++;
        h[steps] = (x[i]+x[i+1])/2.0;
        steps++;
        h[steps] = x[i+1];
        changed = 1;
      }
      else
      {
        // take next coordinate
        steps++;
        h[steps] = x[i+1];
      }
    }
    N_x = steps;
    memcpy(x,h,(N_x+1)*SizeOfDouble);

    // back direction
    h[0] = x[N_x];
    h[1] = x[N_x-1];
    steps = 1;
    for (i=N_x-1;i>=1;i--)
    {
      // previous interval too large
      if ((x[i]-x[i-1])/(x[i+1]-x[i]) > 2.01)
      {
        // insert new coordinate
        steps++;
        h[steps] = (x[i]+x[i-1])/2.0;
        steps++;
        h[steps] = x[i-1];
        changed = 1;
      }
      else
      {
        // take next coordinate
        steps++;
        h[steps] = x[i-1];
      }
    }
    N_x = steps;
    for (i=0;i<=N_x;i++)
    {
      x[i] = h[N_x-i];
    }
  }
}

void ComputeMatrixForPressureLaplacian(TSquareMatrix2D *SqmatrixPres, 
BoundCondFunct2D *BoundConditionPressureLaplace,
BoundValueFunct2D *PressureBoundValueLaplace,
int N_x, int N_y, double *x_coord, double *y_coord)
{
  int i, j, N_Unknowns, N_Entries, N2, N_x1, begin, end, bdr_number;
  int *Row, *KCol;
  double h_i, h_i1, *values;
  
  // get data of the sparse matrix
  N_Unknowns = SqmatrixPres->GetN_Columns();
  Row = SqmatrixPres->GetRowPtr();
  KCol = SqmatrixPres->GetKCol();
  N_Entries = SqmatrixPres->GetN_Entries();
  values = SqmatrixPres->GetEntries();
  // initialize
  memset(values, 0, N_Entries*SizeOfDouble);

  N_x1 = N_x + 1;
  N2 = N_x1 * (N_y + 1);
  
  // loop over all nodes of the FDM grid
  for ( i=0 ; i<N2 ; i++ )
  {
    // interior nodes
    // HAS TO BE EXTENDED FOR NON-RECTANGULAR DOMAINS
    if (( i%(N_x1)!=0 )&&( (i+1)%(N_x1)!=0 )&&(i<= N_x1*N_y)&&(i>N_x))
    {      
      // viscous term wrt x
      h_i = x_coord[i] - x_coord[i-1];
      h_i1 = x_coord[i+1] - x_coord[i];
      // find entries
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
	if ( KCol[j] == i-1)
	  values[j] -= 2.0/(h_i * (h_i + h_i1));
	if ( KCol[j] == i)
	  values[j] += 2.0/(h_i * (h_i + h_i1)) + 2.0/(h_i1 * (h_i + h_i1));
	if ( KCol[j] == i+1)
	  values[j] -= 2.0/(h_i1 * (h_i + h_i1));	
      }
      // viscous term wrt y
      h_i = y_coord[i] - y_coord[i-N_x1];
      h_i1 = y_coord[i+N_x1] - y_coord[i];
      // find entries
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
	if ( KCol[j] == i - N_x1)
	  values[j] -= 2.0/(h_i * (h_i + h_i1));
	if ( KCol[j] == i)
	  values[j] += 2.0/(h_i * (h_i + h_i1)) + 2.0/(h_i1 * (h_i + h_i1));
	if ( KCol[j] == i + N_x1)
	  values[j] -= 2.0/(h_i1 * (h_i + h_i1));	
      }
    }
    else
    {
      // boundary nodes
    }
  }
}

#ifdef __3D__
void Compute_DOF_Conversion_Q1_FDM_UnitCube3D(TCollection *coll,
TFEFunction3D *fesol, int N_x, int N_y, int N_z, double h,
int* &dof_conversion, double* &x_coord, double * &y_coord, double* &z_coord)
{
  int i, j, k, index, dof, N_Cells, N_U, N_V, N3, N2;
  int *begin_index, *global_numbers, *globdof;
  double eps=1e-6, x[8], y[8], z[8];
  TFESpace3D *fespace;
  TBaseCell *cell;
  TVertex *vertex;

  N2 = (N_x+1)*(N_y+1);
  N3 = (N_x+1)*(N_y+1)*(N_z+1);
  if(TDatabase::ParamDB->P9!=4711.0)
  {
    x_coord = new double[N3];
    y_coord = new double[N3];
    z_coord = new double[N3];
  }
  else
  {
    i = (int)(TDatabase::ParamDB->INTERNAL_VERTEX_X[0]+1e-6);
    x_coord = new double[i];
    y_coord = new double[i];
    z_coord = new double[i];
  }

  for (i=0 ; i<N_z+1 ; i++ )
  {
    for (j=0; j<N_y+1; j++)
    {
      for (k=0; k<N_x+1; k++)
      {
        index = i*N2 + j * (N_x+1) + k;
        // compute x and y
        x_coord[index] = k*h;
        y_coord[index] = j*h;
        z_coord[index] = i*h;
        //OutPut(index << " " << x_coord[index] << " " << y_coord[index] << " " << z_coord[index] << endl);
      }
    }
  }

  // extract information about the fe space and the global d.o.f.
  fespace = fesol->GetFESpace3D();
  // number of dof
  N_U = fespace->GetN_DegreesOfFreedom();
  // allocate array that stores information on conversion
  dof_conversion = new int[N_U];
  // change enumeration of coordinates
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // number of mesh cells
  N_Cells = coll->GetN_Cells();

  // loop over the mesh cells of the global grid
  // find conversion of numbering
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 8)
    {
      OutPut("FDM only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    globdof = global_numbers+begin_index[i];
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j],z[j]);
      // compute global dof at this vertex
      switch(j)
      {
        case 0:
          dof = globdof[0];
          break;
        case 1:
          dof = globdof[1];
          break;
        case 2:
          dof = globdof[3];
          break;
        case 3:
          dof = globdof[2];
          break;
        case 4:
          dof = globdof[4];
          break;
        case 5:
          dof = globdof[5];
          break;
        case 6:
          dof = globdof[7];
          break;
        case 7:
          dof = globdof[6];
          break;
      }
      // this is not elegant
      for (k=0;k<N3;k++)
      {
        if ((fabs(x[j]-x_coord[k])<eps)&&(fabs(y[j]-y_coord[k])<eps)&&(fabs(z[j]-z_coord[k])<eps))
        {
          dof_conversion[k] = dof;
          break;
        }
      }
    }
  }
}


/******************************************************************************/
// CheckAdaptiveCriterion
//
// adaptive algorithm for FDM-ENO scheme, the grid will be adaptively changed
//
// input:  sol       - solution on the current grid
//         N_x       - number of intervals in x direction of current grid
//         N_y       - number of intervals in y direction of current grid
//         N_Unknowns- number of unknowns, current grid
//         x_coord   - x coordinates of current grid
//         y_coord   - y coordinates of current grid
//
// output: sol       - solution on the new grid
//         N_x       - number of intervals in x direction of new grid
//         N_y       - number of intervals in y direction of new grid
//         N_Unknowns- number of unknowns, new grid
//         x_coord   - x coordinates of new grid
//         y_coord   - y coordinates of new grid
/******************************************************************************/
void CheckAdaptiveCriterion(double *sol, int &N_x, int &N_y,  int &N_z, int &N_Unknowns,
double *x_coord, double *y_coord, double *z_coord,
double *aux_adaptive)
{
  const int max_n = 2048;
  int i, j, k, index, index1, index2, index3, x_steps, y_steps, z_steps;
  int steps, size, keep;
  int Nx1, Ny1, Nz1, N2;
  double val, maximal, x[max_n], y[max_n], z[max_n], h[max_n];
  double  max_interval = 0.1, min_interval = 1e-3, lower = 1e-3, upper = 0.1;
  double *sol_curr, *sol_curr_x, *x_curr_x, *y_curr_x, *z_curr_x;

  Nx1 = N_x+1;
  Ny1 = N_y+1;
  N2 = Nx1*Ny1;
  Nz1 = N_z+1;
  size = (int)(TDatabase::ParamDB->INTERNAL_VERTEX_X[0]);
  sol_curr = aux_adaptive;
  index = 0;
  for (k=0;k<Nz1;k++)
    for (j=0;j<Ny1;j++)
      for (i=0;i<Nx1;i++)
      {
        sol[index] = x_coord[index] +1/(y_coord[index]*  y_coord[index]+1e-6) - exp(z_coord[index]*  z_coord[index] *  z_coord[index]);
        index++;
      }
  memcpy(sol_curr, sol, Nx1*Ny1*Nz1*SizeOfDouble);
  sol_curr_x = aux_adaptive + size;
  x_curr_x = sol_curr_x + size;
  y_curr_x = x_curr_x + size;
  z_curr_x = y_curr_x + size;

  // number of already taken lines
  x_steps = -1;
  // for x=0
  keep = 1;
  // loop over the x - values
  for(i=1;i<=N_x;i++)
  {
    // compute maximal difference of function values in the interval
    maximal = 0;
    for (j=0;j<=N_y;j++)
    {
      for (k=0;k<=N_z;k++)
      {
        index = k*N2 + j*Nx1 + i;
        val = fabs(sol_curr[index]-sol_curr[index-1]);
        if (val > maximal)
          maximal = val;
      }
    }
    // remove coordinate i-1
    // it is keep==0
    if ((maximal < lower)&&( x_coord[i] - x[x_steps] <  max_interval)&&(i>1)&&(!keep))
    {
      continue;
    }
    // insert coordinate
    if ((maximal >= upper)&&((x_coord[i] -  x_coord[i-1])/2.0 > min_interval))
    {
      x_steps++;
      x[x_steps] = x_coord[i-1];
      x_steps++;
      x[x_steps] = val = (x_coord[i] +  x_coord[i-1])/2.0;
      keep = 1;
      continue;
    }
    // keep coordinate
    x_steps++;
    x[x_steps] = x_coord[i-1];
    keep = 0;
  }
  // last coordinate
  x[x_steps+1] = x_coord[N_x];
  x_steps++;
  //for (i=0;i<=x_steps;i++)
  //OutPut(i << " " << x[i] << endl);

  GradeMesh(x_steps, x, h);
  for (i=0;i<=x_steps;i++)
    OutPut(i << " " << x[i] << endl);
  steps = 0;
  // solution on new grid
  // loop over new x coordinates
  for (i=0;i<x_steps+1;i++)
  {
    while(1)
    {
      // take values
      // same x coordinate
      if (fabs(x[i]-x_coord[steps])<1e-6)
      {
        for (j=0;j<N_y+1;j++)
        {
          for (k=0;k<=N_z;k++)
          {
            index =  k*N2 + j* Nx1 + steps;
            index1 = k*(x_steps+1)*Ny1 + j*(x_steps+1) + i;
            sol_curr_x[index1] = sol_curr[index];
            x_curr_x[index1] = x_coord[index];
            y_curr_x[index1] = y_coord[index];
            z_curr_x[index1] = z_coord[index];
          }
        }
        // increase counter for old x coordinates
        steps++;
        break;
      }
      // interpolate
      if (x[i] < x_coord[steps])
      {
        for (j=0;j<N_y+1;j++)
        {
          for (k=0;k<=N_z;k++)
          {
            index = k * N2 + j*Nx1 + steps-1;
            // left point
            index1 = k * (x_steps+1)*Ny1 + j*(x_steps+1) + i;
            // right point, y-coordinate stays the same
            index2 = index++;
            sol_curr_x[index1] = (sol_curr[index]+sol_curr[index2])/2.0;
            x_curr_x[index1] = (x_coord[index]+ x_coord[index2])/2.0;
            y_curr_x[index1] = y_coord[index];
            z_curr_x[index1] = z_coord[index];
          }
        }
        break;
      }
      // check next old x coordinate
      steps++;
    }
    // otherwise, coordinate was deleted
  }
  // set new number of intervals in x direction
  N_x = x_steps;
  N2 = (N_x+1)*(N_y+1);

  // number of already taken lines
  y_steps = -1;
  // keep y=0
  keep = 1;
  // loop over the y - values
  for(i=1;i<=N_y;i++)
  {
    // compute maximal difference of function values in the interval
    maximal = 0;
    for (j=0;j<=x_steps;j++)
    {
      for (k=0;k<=N_z;k++)
      {
        index = k*N2 + i*(x_steps+1) + j;
        val = fabs(sol_curr_x[index]-sol_curr_x[index-x_steps-1]);
        if (val > maximal)
          maximal = val;
      }
    }
    //OutPut("y " << y_curr_x[i*(x_steps+1)] << " max " << maximal << endl);
    // remove coordinate i-1
    // coordinate on bottom level
    index1 = i*(x_steps+1);
    if ((maximal < lower)&&( y_curr_x[index1] - y[y_steps] <  max_interval)&&(i>1)&&(!keep))
    {
      continue;
    }
    // insert coordinate
    index = index1 - (x_steps+1);
    if ((maximal >= upper)&& ((y_curr_x[index1] - y_curr_x[index])/2.0> min_interval))
    {
      y_steps++;
      y[y_steps] = y_curr_x[index];
      y_steps++;
      y[y_steps] = (y_curr_x[index1] + y_curr_x[index])/2.0;
      keep = 1;
      continue;
    }
    // keep coordinate
    y_steps++;
    y[y_steps] = y_curr_x[index];
    keep = 0;
  }
  // last coordinate
  y[y_steps+1] = y_curr_x[(x_steps+1)*(N_y+1)-1];
  y_steps++;
  //  for (i=0;i<=y_steps;i++)
  //    OutPut(i << " " << y[i] << endl);
  GradeMesh(y_steps, y, h);
  for (i=0;i<=y_steps;i++)
    OutPut(i << " " << y[i] << endl);

  // new number of y-intervals
  N_y = y_steps;

  steps = 0;
  // solution on new grid
  // loop over the new y coordinates
  for (i=0;i<y_steps+1;i++)
  {
    while(1)
    {
      // left coordinate of y level on old grid
      index2 = steps*(x_steps+1);
      index3 = i*(x_steps+1);
      // same coordinate take values
      if (fabs(y[i]-y_curr_x[index2])<1e-6)
      {
        for (j=0;j<x_steps+1;j++)
        {
          for (k=0;k<=N_z;k++)
          {
            sol_curr[k*(N_x+1)*(N_y+1)+index3+j] = sol_curr_x[k*N2+index2+j];
            x_coord[k*(N_x+1)*(N_y+1)+index3+j] = x_curr_x[k*N2+index2+j];
            y_coord[k*(N_x+1)*(N_y+1)+index3+j] = y_curr_x[k*N2+index2+j];
            z_coord[k*(N_x+1)*(N_y+1)+index3+j] = z_curr_x[k*N2+index2+j];
            //OutPut(k*N2+index2+j << " " <<  z_curr_x[k*N2+index2+j] << endl);
            //OutPut(k*(N_x+1)*(N_y+1)+index3+j << " " << z_coord[k*(N_x+1)*(N_y+1)+index3+j] << endl);
          }
        }
        steps++;
        break;
      }
      // interpolate
      if (y[i] < y_curr_x[index2])
      {
        for (j=0;j<x_steps+1;j++)
        {
          for (k=0;k<=N_z;k++)
          {
            index = k*N2+index2 + j;
            index1 = k*(N_x+1)*(N_y+1)+index3 + j;
            sol_curr[index1] = (sol_curr_x[index]+sol_curr_x[index-x_steps-1])/2.0;
            x_coord[index1] = x_curr_x[index];
            y_coord[index1] = (y_curr_x[index]+y_curr_x[index-x_steps-1])/2.0;
            z_coord[index1] = z_curr_x[index];
          }
        }
        break;
      }
      steps++;
    }
    // otherwise, coordinate was deleted
  }
  /*  index = 0;
    for (k=0;k<=N_z;k++)
    {
        for (j=0;j<=N_y;j++)
        {
      for (i=0;i<=N_x;i++)
      {
          OutPut( index << " " << x_coord[index] << " " << y_coord[index]  << " " << z_coord[index]
            << " " << sol_curr[index] << " " << x_coord[index] + 1/(y_coord[index]*  y_coord[index]+1e-6) -
            exp(z_coord[index]*  z_coord[index] *  z_coord[index]) << endl);
          if (fabs(sol_curr[index]-(x_coord[index] + 1/(y_coord[index]*  y_coord[index]+1e-6) -
  exp(z_coord[index]*  z_coord[index] *  z_coord[index])))>1e-6)
  {
  OutPut("error "<< endl);
  exit(1);
  }
  index++;
  }
  }
  }
  exit(1);*/
  N2 = (N_x+1)*(N_y+1);
  // number of already taken lines
  z_steps = -1;
  // keep y=0
  keep = 1;
  // loop over the y - values
  for(i=1;i<=N_z;i++)
  {
    // compute maximal difference of function values in the interval
    maximal = 0;
    for (j=0;j<=x_steps;j++)
    {
      for (k=0;k<=y_steps;k++)
      {
        index = i*N2 + k*(x_steps+1) + j;
        val = fabs(sol_curr_x[index]-sol_curr_x[index-N2]);
        if (val > maximal)
          maximal = val;
      }
    }
    //OutPut("y " << y_curr_x[i*(x_steps+1)] << " max " << maximal << endl);
    // remove coordinate i-1
    index1 = i*N2;
    if ((maximal < lower)&&(z_coord[index1] - z[z_steps] <  max_interval)&&(i>1)&&(!keep))
    {
      continue;
    }
    // insert coordinate
    index = index1 - N2;
    if ((maximal >= upper)&& ((z_coord[index1] - z_coord[index])/2.0> min_interval))
    {
      z_steps++;
      z[z_steps] = z_coord[index];
      z_steps++;
      z[z_steps] = (z_coord[index1] + z_coord[index])/2.0;
      keep = 1;
      continue;
    }
    // keep coordinate
    z_steps++;
    z[z_steps] = z_coord[index];
    //OutPut(i << " " << N2 << " " << index << "  "  << z_coord[index] << " " << endl);
    keep = 0;
  }
  // last coordinate
  z[z_steps+1] = z_coord[N2*(N_z+1)-1];
  z_steps++;
  //  for (i=0;i<=y_steps;i++)
  //    OutPut(i << " " << y[i] << endl);
  GradeMesh(z_steps, z, h);
  for (i=0;i<=z_steps;i++)
    OutPut(i << " " << z[i] << endl);
  // new number of unknowns
  N_z = z_steps;
  N_Unknowns = (N_x+1)*(N_y+1)*(N_z+1);
  OutPut("new grid "<< N_x+1 << " x " << N_y+1 << " x " << N_z+1 << " entries " <<   N_Unknowns  << endl);
  if (N_Unknowns > TDatabase::ParamDB->INTERNAL_VERTEX_X[0])
  {
    OutPut("grid too fine, not enough memory !!!" << endl);
    exit(4711);
  }
  steps = 0;
  // solution on new grid
  // loop over the new y coordinates
  for (i=0;i<z_steps+1;i++)
  {
    while(1)
    {
      // left coordinate of y level on old grid
      index2 = steps*N2;
      index3 = i*N2;
      //OutPut(i << " " << z[i] << " " << index2 << " " << index3 << endl);
      // same coordinate take values
      if (fabs(z[i]-z_coord[index2])<1e-6)
      {
        for (j=0;j<N_x+1;j++)
        {
          for (k=0;k<N_y+1;k++)
          {
            //OutPut(index3+k*(N_x+1)+j << " " << index2+k*(N_x+1)+j << endl);
            sol_curr_x[index3+k*(N_x+1)+j] = sol_curr[index2+k*(N_x+1)+j];
            //OutPut(index3+k*(N_x+1)+j << " done " << index2+k*(N_x+1)+j << endl);
            //OutPut(k*N2+index2+j << " " <<  z_curr_x[k*N2+index2+j] << endl);
            //OutPut(k*(N_x+1)*(N_y+1)+index3+j << " " << z_coord[k*(N_x+1)*(N_y+1)+index3+j] << endl);
          }
        }
        steps++;
        break;
      }
      // interpolate
      if (z[i] < z_coord[index2])
      {
        for (j=0;j<N_x+1;j++)
        {
          for (k=0;k<N_y+1;k++)
          {
            index = steps*N2+ k*(N_x+1)+j;
            index1 = i*N2+k*(N_x+1)+j;
            sol_curr_x[index1] = (sol_curr[index]+sol_curr[index-N2])/2.0;
          }
        }
        break;
      }
      steps++;
    }
    // otherwise, coordinate was deleted
  }
  // copy solution to array sol
  memcpy(sol, sol_curr_x, N_Unknowns*SizeOfDouble);
  OutPut("copy done" << endl);
  // fill arrays with new coordinates
  index = 0;
  for (k=0;k<=N_z;k++)
  {
    for (i=0;i<=N_y;i++)
    {
      for (j=0;j<=N_x;j++)
      {
        x_coord[index] = x[j];
        y_coord[index] = y[i];
        z_coord[index] = z[k];
        index++;
      }
    }
  }
  OutPut("coord done" << endl);
  index = 0;
  for (k=0;k<=N_z;k++)
  {
    for (j=0;j<=N_y;j++)
    {
      for (i=0;i<=N_x;i++)
      {
        OutPut( index << " " << x_coord[index] << " " << y_coord[index]  << " " << z_coord[index]
          << " " << sol[index] << " " <<
          x_coord[index] + 1/(y_coord[index]*  y_coord[index]+1e-6) -
          exp(z_coord[index]*  z_coord[index] *  z_coord[index]) << endl);
        if (fabs(sol[index]-(x_coord[index] + 1/(y_coord[index]*  y_coord[index]+1e-6) -
          exp(z_coord[index]*  z_coord[index] *  z_coord[index]))) > 1e-6)
        {
          OutPut("error "<< fabs(sol[index]-x_coord[index] + 1/(y_coord[index]*  y_coord[index]+1e-6) -
            exp(z_coord[index]*  z_coord[index] *  z_coord[index])) <<endl);
          //exit(1);
        }
        index++;
      }
    }
  }
  WriteVTKFDM("test_vtk.vtk", sol, x_coord,
    y_coord, z_coord, N_x, N_y, N_z);

  exit(1);
}
#endif

/******************************************************************************/
// DividedDifferences
// compute divided differences
/******************************************************************************/
double DividedDifferences(int n, double *val, double *coord)
{
  double value, value1;

  switch(n)
  {
    case 0: value = val[0];
    break;
    case 1: value = (val[1]-val[0])/(coord[1]-coord[0]);
    break;
    case 2: value = (val[1]-val[0])/(coord[1]-coord[0]);
    value1 = (val[2]-val[1])/(coord[2]-coord[1]);
    value = (value1 -  value)/(coord[2]-coord[0]);
    break;
    case 3: value = DividedDifferences(2, val, coord);
    value1 = DividedDifferences(2, val+1, coord+1);
    value = (value1 -  value)/(coord[3]-coord[0]);
    break;
    default:
      OutPut("Divided differences not implemented !!! "<< endl);
      exit(4711);
  }
  return(value);
}


void DividedDifferences_3_2(double *val, double *coord, double *result)
{

  double value, value1, value2;

  value = DividedDifferences(2, val, coord);
  value1 = DividedDifferences(2, val+1, coord+1);
  value2 = DividedDifferences(2, val+2, coord+2);
  result[0] = (value1 -  value)/(coord[3]-coord[0]);
  result[1] = (value2 -  value1)/(coord[4]-coord[1]);
}


double DividedDifferences2(double val0, double val1, double val2,
double coord0, double coord1, double coord2)
{
  double value, value1;

  value = (val1-val0)/(coord1-coord0);
  value1 = (val2-val1)/(coord2-coord1);
  value = (value1 -  value)/(coord2-coord0);
  return(value);
}


/******************************************************************************/
// Conversion of arrays
/******************************************************************************/

void FEM2FDM(int N, double *fem_array, double *fdm_array, int *dof_conversion)
{
  int i, j;

  for (i=0;i<N;i++)
  {
    j =  dof_conversion[i];
    fdm_array[i] = fem_array[j];
  }
}


void FDM2FEM(int N, double *fdm_array, double *fem_array, int *dof_conversion)
{
  int i, j;

  for (i=0;i<N;i++)
  {
    j =  dof_conversion[i];
    fem_array[j] = fdm_array[i];
  }
}


/******************************************************************************/
// ComputeStages_FDM2D
// computes the stages in explicit RK methods combined with a FDM
/******************************************************************************/

#ifdef __2D__
void ComputeStages_FDM2D(int dim, CoeffFct2D *Coeffs, BoundCondFunct2D *bound_cond,
BoundValueFunct2D *bound_val,
double *sol, double **stages, int current_stage, int N_x, int N_y,
int *dof_conversion, double *x_coord, double *y_coord)
{
  int i, j, N_x1, N_y1, N2, N1_[2];
  double tau, hx_i, hx_i1, hy_i, hy_i1, valx, valy, valz, val[6], coord[6], value;
  double *sol_curr, coeff_array[20], *coeff, *current_stage_fdm, *sol_help;
  double d1, d2, d3, d4, uhx[3], omega[3], alpha[3], d[3], beta[3], av[5], c_e = 1e-6;
  double *coordinates[2];
  int *offset_ = NULL, *offset1_;
  BoundCond cond;

  // set time of the stage, in main already new time level
  tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  TDatabase::TimeDB->CURRENTTIME += (TDatabase::TimeDB->RK_c[current_stage] - 1)* tau;
  N_x1 = N_x+1;
  N_y1 = N_y+1;
  N2 = (N_x1)*(N_y1);
  // allocate arrays
  sol_curr = stages[5];
  coeff = coeff_array;

  // compute linear combination of current stages
  // numeration of dof corresponding to fe space
  current_stage_fdm = stages[current_stage];
  memset(current_stage_fdm, 0, N2 * SizeOfDouble);
  // initialize current stage
  memcpy(sol_curr,sol,N2*SizeOfDouble);
  // add previous stages
  for (i=0;i<current_stage;i++)
    Daxpy(N2, tau*TDatabase::TimeDB->RK_A[current_stage][i], stages[i], sol_curr);

  //N1_[0] = N_x1;
  //N1_[1] = N_y1;
  //coordinates[0] = x_coord;
  //coordinates[1] = y_coord;
  //InitializeConvectiveTermFDM(2, offset_, offset1_, N1_);

  // discretization of pde with FDM
  // loop over all nodes of the FDM grid
  for ( i=0 ; i<N2 ; i++ )
  {
    // coefficients of the equation
    Coeffs(1, &x_coord[i], &y_coord[i], NULL, &coeff);
    // compute the diffusive term
    if (( i%(N_x1)!=0 )&&( (i+1)%(N_x1)!=0 )&&(i<= N_x1*N_y)&&(i>N_x))
    {
      // diffusive term wrt x
      val[0] = sol_curr[i-1];
      val[1] = sol_curr[i];
      val[2] = sol_curr[i+1];
      coord[0] = x_coord[i-1];
      coord[1] = x_coord[i];
      coord[2] = x_coord[i+1];
      valx = 2*DividedDifferences(2, val, coord);
      // diffusive term wrt y
      val[0] = sol_curr[i-N_x1];
      val[2] = sol_curr[i+N_x1];
      coord[0] = y_coord[i-N_x1];;
      coord[1] = y_coord[i];
      coord[2] = y_coord[i+N_x1];
      valy = 2*DividedDifferences(2, val, coord);
      // diffusive term wrt x
      // add also reactive term and rhs
      current_stage_fdm[i] = coeff[0]*(valx + valy) - coeff[3]*sol_curr[i] + coeff[4];
      //OutPut(coeff[0] << " " << coeff[3] << " " << coeff[4] << endl);
    }

    /* ConvectiveTermFDM(2, i,
          coeff+1, sol_curr, current_stage_fdm, coordinates,
          offset_, offset1_);
          continue;*/
    switch(TDatabase::ParamDB->DISCTYPE)
    {
      case GALERKIN:
        exit(1);
        break;
        // simple upwind scheme
      case UPWIND:
        // compute the term A_x
        if (coeff[1] >= 0)
        {
          // not on boundary x = x_min -> in the "x-sense" exists a left neighbour
          if ( i%(N_x1)!=0 )
            current_stage_fdm[i] -= coeff[1]*(sol_curr[i]-sol_curr[i-1])/(x_coord[i]-x_coord[i-1]);
        }
        else
        {
          // not on boundary x = x_max -> in the "x-sense" exists a right neighbour
          if ( (i+1)%(N_x1)!=0 )
            current_stage_fdm[i] -= coeff[1]*(sol_curr[i+1]-sol_curr[i])/(x_coord[i+1]-x_coord[i]);
        }

        // compute the term A_y
        if (coeff[2] >= 0)
        {
          // not on boundary y = y_min -> in the "y-sense" exists a left neighbour
          if ((i%((N_x1)*(N_y1)))>N_x )
            current_stage_fdm[i] -= coeff[2]*(sol_curr[i]-sol_curr[i-(N_x1)])/(y_coord[i]-y_coord[i-(N_x1)]);
        }
        else
        {
          // not on boundary y = y_max -> in the "y-sense" exists a right neighbour
          if ((i%((N_x1)*(N_y1)))<((N_x1)*N_y) )
            current_stage_fdm[i] -= coeff[2]*(sol_curr[i+(N_x1)]-sol_curr[i])/(y_coord[i+(N_x1)]-y_coord[i]);
        }
        break;
      case ENO_3:
        // compute the term A_x
        coeff[1] = - coeff[1];
        coeff[2] = - coeff[2];
        if (coeff[1] <=  0)
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+1];
            val[5] = sol_curr[i+2];
            coord[0] = 2 * x_coord[i]-x_coord[i+3];
            coord[1] = 2 * x_coord[i]-x_coord[i+2];
            coord[2] = 2 * x_coord[i]-x_coord[i+1];
            coord[3] = x_coord[i];
            coord[4] = x_coord[i+1];
            coord[5] = x_coord[i+2];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )
            {
              val[0] = val[1] = val[2] = sol_curr[i-1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+1];
              val[5] = sol_curr[i+2];
              coord[0] = 2 * x_coord[i-1]-x_coord[i+1];
              coord[1] = 2 * x_coord[i-1]-x_coord[i];
              coord[2] = x_coord[i-1];
              coord[3] = x_coord[i];
              coord[4] = x_coord[i+1];
              coord[5] = x_coord[i+2];
            }
            else
            {
              // next layer on left boundary
              if ( (i-2)%(N_x1)==0 )
              {
                val[0] = val[1] = sol_curr[i-2];
                val[2] = sol_curr[i-1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+1];
                val[5] = sol_curr[i+2];
                coord[0] = 2 * x_coord[i-2]-x_coord[i-1];
                coord[1] = x_coord[i-2];
                coord[2] = x_coord[i-1];
                coord[3] = x_coord[i];
                coord[4] = x_coord[i+1];
                coord[5] = x_coord[i+2];
              }
              else
              {
                // point on right boundary
                if ( (i+1)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-3];
                  val[1] = sol_curr[i-2];
                  val[2] = sol_curr[i-1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = x_coord[i-3];
                  coord[1] = x_coord[i-2];
                  coord[2] = x_coord[i-1];
                  coord[3] = x_coord[i];
                  coord[4] = 2 *  x_coord[i] - x_coord[i-1];
                  coord[5] = 2 *  x_coord[i] - x_coord[i-2];
                }
                else
                {
                  // point next to right boundary
                  if ( (i+2)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+1];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = 2 *  x_coord[i+1] - x_coord[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+1];
                    val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = x_coord[i+2];
                    /*sol_help = sol_curr+i-3;
                                val[0] = *sol_help;
                                val[1] = *(sol_help+1);
                                val[2] = *(sol_help+2);
                                val[3] = *(sol_help+3);
                                val[4] = *(sol_help+4);
                                val[5] = *(sol_help+5);
                                sol_help = x_coord+i-3;
                                coord[0] = *sol_help;
                                coord[1] = *(sol_help+1);
                                coord[2] = *(sol_help+2);
                    coord[3] = *(sol_help+3);
                    coord[4] = *(sol_help+4);
                    coord[5] = *(sol_help+5);*/
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+2,coord+2);
          d2 = DividedDifferences(2,val+1,coord+1);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              //value = DividedDifferences(1,val+1,coord+1);
              value =  (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d3 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3a ");
            }
            else
            {
              //value = DividedDifferences(1,val+2,coord+2);
              valx = coord[3]-coord[2];
              value = (val[3]-val[2])/valx;
              value += d1 * valx;
              value += d4 * valx*(coord[3]-coord[4]);
              //OutPut("a2 ");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              //value = DividedDifferences(1,val,coord);
              value = (val[1]-val[0])/(coord[1]-coord[0]);
              value += DividedDifferences(2,val,coord) * (coord[3]-coord[0] + coord[3]-coord[1]);
              valx = coord[3]-coord[0];
              valy = coord[3]-coord[1];
              valz = coord[3]-coord[2];
              //value += d3 * ((coord[3]-coord[1])*(coord[3]-coord[2])
              //  + (coord[3]-coord[0])*(coord[3]-coord[2]) + (coord[3]-coord[0])*(coord[3]-coord[1]));
              value += d3 * (valy * valz + valx * valz + valx * valy);
              //OutPut("a4 ");
            }
            else
            {
              //value = DividedDifferences(1,val+1,coord+1);
              value = (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3b ");
              // OutPut("a3 " << value << "  " << DividedDifferences(1,val+3,coord+3)<<":");
            }
          }
        }
        else
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+1];
            val[4] = sol_curr[i+2];
            val[5] = sol_curr[i+3];
            coord[0] = 2 * x_coord[i]-x_coord[i+2];
            coord[1] = 2 * x_coord[i]-x_coord[i+1];
            coord[2] = x_coord[i];
            coord[3] = x_coord[i+1];
            coord[4] = x_coord[i+2];
            coord[5] = x_coord[i+3];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )
            {
              val[0] = val[1] = sol_curr[i-1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+1];
              val[4] = sol_curr[i+2];
              val[5] = sol_curr[i+3];
              coord[0] = 2 * x_coord[i-1]-x_coord[i];
              coord[1] = x_coord[i-1];
              coord[2] = x_coord[i];
              coord[3] = x_coord[i+1];
              coord[4] = x_coord[i+2];
              coord[5] = x_coord[i+3];
            }
            else
            {
              // point on right boundary
              if ( (i+1)%(N_x1)==0 )
              {
                val[0] = sol_curr[i-2];
                val[1] = sol_curr[i-1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = x_coord[i-2];
                coord[1] = x_coord[i-1];
                coord[2] = x_coord[i];
                coord[3] = 2 *  x_coord[i] - x_coord[i-1];
                coord[4] = 2 *  x_coord[i] - x_coord[i-2];
                coord[5] = 2 *  x_coord[i] - x_coord[i-3];
              }
              else
              {
                // point next to right boundary
                if ( (i+2)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-2];
                  val[1] = sol_curr[i-1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+1];
                  coord[0] = x_coord[i-2];
                  coord[1] = x_coord[i-1];
                  coord[2] = x_coord[i];
                  coord[3] = x_coord[i+1];
                  coord[4] = 2 *  x_coord[i+1] - x_coord[i];
                  coord[5] = 2 *  x_coord[i+1] - x_coord[i-1];
                }
                else
                {
                  // point in next layer to right boundary
                  if ( (i+3)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = 2 *  x_coord[i+2] - x_coord[i+1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = sol_curr[i+2];
                    val[5] = sol_curr[i+3];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = x_coord[i+3];
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+1,coord+1);
          d2 = DividedDifferences(2,val+2,coord+2);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[2]-coord[0] + coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[0])*(coord[2]-coord[1]);
              //OutPut("A3 ");
            }
            else
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d4 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              //OutPut("A2a ");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              valx = coord[2]-coord[1];
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * valx;
              value += d3 * valx*(coord[2]-coord[3]);
              //OutPut("A2b ");
            }
            else
            {
              valx = coord[2]-coord[3];
              value = DividedDifferences(1,val+2,coord+2);
              value += d2 * valx;
              value += d4 * valx*(coord[2]-coord[4]);
              //comparison with backward difference in node i
              //OutPut("A1b ");
              //OutPut("A1 " << value <<"  " <<   DividedDifferences(1,val+1,coord+1) << ":");
            }
          }
        }                                         // end coeff[1] <=0
        //OutPut(value << " ");
        current_stage_fdm[i] += coeff[1] * value;
        // compute the term A_y
        if (coeff[2]  <= 0)
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<=N_x )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+N_x1];
            val[5] = sol_curr[i+2 * N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+3*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[2] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[3] = y_coord[i];
            coord[4] = y_coord[i+N_x1];
            coord[5] = y_coord[i+2 * N_x1];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>=N_x1)&&(i<2*N_x1) )
            {
              val[0] = val[1] = val[2] = sol_curr[i-N_x1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+N_x1];
              val[5] = sol_curr[i+2*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i+N_x1];
              coord[1] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[2] = y_coord[i-N_x1];
              coord[3] = y_coord[i];
              coord[4] = y_coord[i+N_x1];
              coord[5] = y_coord[i+2*N_x1];
            }
            else
            {
              // next layer on bottom boundary
              if ( (i>=2*N_x1)&&(i<3*N_x1) )
              {
                val[0] = val[1] = sol_curr[i-2*N_x1];
                val[2] = sol_curr[i-N_x1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+N_x1];
                val[5] = sol_curr[i+2*N_x1];
                coord[0] = 2 * y_coord[i-2*N_x1]-y_coord[i-N_x1];
                coord[1] = y_coord[i-2*N_x1];
                coord[2] = y_coord[i-N_x1];
                coord[3] = y_coord[i];
                coord[4] = y_coord[i+N_x1];
                coord[5] = y_coord[i+2*N_x1];
              }
              else
              {
                // point on top boundary
                if ( i>= N_y*N_x1)
                {
                  val[0] = sol_curr[i-3*N_x1];
                  val[1] = sol_curr[i-2*N_x1];
                  val[2] = sol_curr[i-N_x1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = y_coord[i-3*N_x1];
                  coord[1] = y_coord[i-2*N_x1];
                  coord[2] = y_coord[i-N_x1];
                  coord[3] = y_coord[i];
                  coord[4] = 2 *  y_coord[i] - y_coord[i-N_x1];
                  coord[5] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                }
                else
                {
                  // point next to top boundary
                  if (( i>= (N_y-1)*N_x1)&&(i<N_y*N_x1))
                  {
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i];
                  }
                  // inner point
                  else
                  {
                    // if (i==14382)
                    //  OutPut("herea "<< coeff[2] << endl);
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+N_x1];
                    val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = y_coord[i+2*N_x1];
                    // if (i==14382)
                    //	OutPut("val " << val[0] -1 << endl);
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+2,coord+2);
          d2 = DividedDifferences(2,val+1,coord+1);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d3 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3a ");
              //OutPut("a3 " << value << " ");
            }
            else
            {
              value = DividedDifferences(1,val+2,coord+2);
              value += d1 * (coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[2])*(coord[3]-coord[4]);
              //OutPut("a2 ");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[3]-coord[0] + coord[3]-coord[1]);
              value += d3 * ((coord[3]-coord[1])*(coord[3]-coord[2])
                + (coord[3]-coord[0])*(coord[3]-coord[2]) + (coord[3]-coord[0])*(coord[3]-coord[1]));
              // OutPut("a4 ");
              //OutPut("a4 " << coeff[2]  << "  " << DividedDifferences(2,val+2,coord+2) << "  "
              // << DividedDifferences(2,val+1,coord+1) << "  " << DividedDifferences(3,val,coord)
              // << "  " << DividedDifferences(3,val+1,coord+1)
              //<< endl);
            }
            else
            {
              // if (i==14382)
              //  OutPut("case2a "<< d3 << " " << d4 << endl);

              value = DividedDifferences(1,val+1,coord+1);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              // comparison with backward difference in node i
              // OutPut("a3b ");
              //OutPut("a3 " << value << "  " << DividedDifferences(1,val+3,coord+3)<<":");
            }
          }
        }
        else
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<=N_x )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+N_x1];
            val[4] = sol_curr[i+2*N_x1];
            val[5] = sol_curr[i+3*N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[2] = y_coord[i];
            coord[3] = y_coord[i+N_x1];
            coord[4] = y_coord[i+2*N_x1];
            coord[5] = y_coord[i+3*N_x1];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>=N_x1)&&(i<2*N_x1) )
            {
              val[0] = val[1] = sol_curr[i-N_x1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+N_x1];
              val[4] = sol_curr[i+2*N_x1];
              val[5] = sol_curr[i+3*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[1] = y_coord[i-N_x1];
              coord[2] = y_coord[i];
              coord[3] = y_coord[i+N_x1];
              coord[4] = y_coord[i+2*N_x1];
              coord[5] = y_coord[i+3*N_x1];
            }
            else
            {
              // point on top boundary
              if ( i>= N_y*N_x1)
              {
                val[0] = sol_curr[i-2*N_x1];
                val[1] = sol_curr[i-N_x1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = y_coord[i-2*N_x1];
                coord[1] = y_coord[i-N_x1];
                coord[2] = y_coord[i];
                coord[3] = 2 *  y_coord[i] - y_coord[i-N_x1];
                coord[4] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                coord[5] = 2 *  y_coord[i] - y_coord[i-3*N_x1];
              }
              else
              {
                // point next to top boundary
                if (( i>= (N_y-1)*N_x1)&&(i<N_y*N_x1))
                {
                  val[0] = sol_curr[i-2*N_x1];
                  val[1] = sol_curr[i-N_x1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+N_x1];
                  coord[0] = y_coord[i-2*N_x1];
                  coord[1] = y_coord[i-N_x1];
                  coord[2] = y_coord[i];
                  coord[3] = y_coord[i+N_x1];
                  coord[4] = 2 *  y_coord[i+N_x1] - y_coord[i];
                  coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i-N_x1];
                }
                else
                {
                  // point in next layer to top boundary
                  if (( i>= (N_y-2)*N_x1)&&(i<(N_y-1)*N_x1))
                  {
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = 2 *  y_coord[i+2*N_x1] - y_coord[i+N_x1];
                  }
                  // inner point
                  else
                  {
                    // if (i==14382)
                    //  OutPut("here "<< coeff[2] << endl);
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = sol_curr[i+2*N_x1];
                    val[5] = sol_curr[i+3*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = y_coord[i+3*N_x1];
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+1,coord+1);
          d2 = DividedDifferences(2,val+2,coord+2);
          // if (i==14382)
          // OutPut("d1 "<< d1 << " d2 " << d2 << endl);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[2]-coord[0] + coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[0])*(coord[2]-coord[1]);
              // OutPut("A3 ");
            }
            else
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d4 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              // OutPut("A2a ");
              //OutPut("A1 " << value << "  " << DividedDifferences(1,val+1,coord+1)<<":");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            // if (i==14382)
            //OutPut("d3 "<< d3 << " d4 " << d4 << endl);
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              // if (i==14382)
              //	OutPut("A2b " << value << endl);
            }
            else
            {
              value = DividedDifferences(1,val+2,coord+2);
              value += d2 * (coord[2]-coord[3]);
              value += d4 * (coord[2]-coord[3])*(coord[2]-coord[4]);
              // comparison with backward difference in node i
              //if (i==14382)
              //	OutPut("A1b " <<  value << endl);
              //OutPut("A1 " << value << "  " << DividedDifferences(1,val+1,coord+1)<<":");
            }
          }
        }                                         // end coeff[2] <=0
        //OutPut(value << endl);
        current_stage_fdm[i] += coeff[2] * value;
        break;
      case WENO_5:
        d[0] = 0.3;
        d[1] = 0.6;
        d[2] = 0.1;
        // compute the term A_x
        coeff[1] = - coeff[1];
        coeff[2] = - coeff[2];
        if (coeff[1] <=  0)
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+1];
            val[5] = sol_curr[i+2];
            coord[0] = 2 * x_coord[i]-x_coord[i+3];
            coord[1] = 2 * x_coord[i]-x_coord[i+2];
            coord[2] = 2 * x_coord[i]-x_coord[i+1];
            coord[3] = x_coord[i];
            coord[4] = x_coord[i+1];
            coord[5] = x_coord[i+2];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )
            {
              val[0] = val[1] = val[2] = sol_curr[i-1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+1];
              val[5] = sol_curr[i+2];
              coord[0] = 2 * x_coord[i-1]-x_coord[i+1];
              coord[1] = 2 * x_coord[i-1]-x_coord[i];
              coord[2] = x_coord[i-1];
              coord[3] = x_coord[i];
              coord[4] = x_coord[i+1];
              coord[5] = x_coord[i+2];
            }
            else
            {
              // next layer on left boundary
              if ( (i-2)%(N_x1)==0 )
              {
                val[0] = val[1] = sol_curr[i-2];
                val[2] = sol_curr[i-1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+1];
                val[5] = sol_curr[i+2];
                coord[0] = 2 * x_coord[i-2]-x_coord[i-1];
                coord[1] = x_coord[i-2];
                coord[2] = x_coord[i-1];
                coord[3] = x_coord[i];
                coord[4] = x_coord[i+1];
                coord[5] = x_coord[i+2];
              }
              else
              {
                // point on right boundary
                if ( (i+1)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-3];
                  val[1] = sol_curr[i-2];
                  val[2] = sol_curr[i-1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = x_coord[i-3];
                  coord[1] = x_coord[i-2];
                  coord[2] = x_coord[i-1];
                  coord[3] = x_coord[i];
                  coord[4] = 2 *  x_coord[i] - x_coord[i-1];
                  coord[5] = 2 *  x_coord[i] - x_coord[i-2];
                }
                else
                {
                  // point next to right boundary
                  if ( (i+2)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+1];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = 2 *  x_coord[i+1] - x_coord[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+1];
                    val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = x_coord[i+2];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/3.0-val[3]/2.0 + val[4] - val[5]/6.0)/(coord[3]-coord[2]);
          uhx[1] = (val[1]/6.0-val[2]+val[3]/2.0 + val[4]/3.0)/(coord[2]-coord[1]);
          uhx[2] = (-val[0]/3.0+1.5*val[1]-3*val[2] + 11.0*val[3]/6.0)/(coord[1]-coord[0]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[3] + av[4];
          valy = 3*av[2] - 4*av[3] + av[4];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[1] - 2*av[2] + av[3];
          valy = av[1] - av[3];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[0] - 2*av[1] +av[2];
          valy = av[0] - 4*av[1] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }
          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }
        else                                      // - coeff[1] > 0
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+1];
            val[4] = sol_curr[i+2];
            val[5] = sol_curr[i+3];
            coord[0] = 2 * x_coord[i]-x_coord[i+2];
            coord[1] = 2 * x_coord[i]-x_coord[i+1];
            coord[2] = x_coord[i];
            coord[3] = x_coord[i+1];
            coord[4] = x_coord[i+2];
            coord[5] = x_coord[i+3];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )
            {
              val[0] = val[1] = sol_curr[i-1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+1];
              val[4] = sol_curr[i+2];
              val[5] = sol_curr[i+3];
              coord[0] = 2 * x_coord[i-1]-x_coord[i];
              coord[1] = x_coord[i-1];
              coord[2] = x_coord[i];
              coord[3] = x_coord[i+1];
              coord[4] = x_coord[i+2];
              coord[5] = x_coord[i+3];
            }
            else
            {
              // point on right boundary
              if ( (i+1)%(N_x1)==0 )
              {
                val[0] = sol_curr[i-2];
                val[1] = sol_curr[i-1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = x_coord[i-2];
                coord[1] = x_coord[i-1];
                coord[2] = x_coord[i];
                coord[3] = 2 *  x_coord[i] - x_coord[i-1];
                coord[4] = 2 *  x_coord[i] - x_coord[i-2];
                coord[5] = 2 *  x_coord[i] - x_coord[i-3];
              }
              else
              {
                // point next to right boundary
                if ( (i+2)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-2];
                  val[1] = sol_curr[i-1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+1];
                  coord[0] = x_coord[i-2];
                  coord[1] = x_coord[i-1];
                  coord[2] = x_coord[i];
                  coord[3] = x_coord[i+1];
                  coord[4] = 2 *  x_coord[i+1] - x_coord[i];
                  coord[5] = 2 *  x_coord[i+1] - x_coord[i-1];
                }
                else
                {
                  // point in next layer to right boundary
                  if ( (i+3)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = 2 *  x_coord[i+2] - x_coord[i+1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = sol_curr[i+2];
                    val[5] = sol_curr[i+3];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = x_coord[i+3];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/2.0+val[1] - val[0]/6.0 - val[3]/3.0)/(-coord[2]+coord[1]);
          uhx[1] = (-val[3] +val[2]/2.0+val[1]/3.0 + val[4]/6.0)/(-coord[3]+coord[2]);
          uhx[2] = (-val[5]/3.0+1.5*val[4]-3*val[3] + 11.0*val[2]/6.0)/(-coord[4]+coord[3]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[1] + av[0];
          valy = 3*av[2] - 4*av[1] + av[0];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[3] - 2*av[2] + av[1];
          valy = av[3] - av[1];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[4] - 2*av[3] +av[2];
          valy = av[4] - 4*av[3] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }
          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }                                         // end coeff[1] <=0
        //OutPut(value << " ");
        current_stage_fdm[i] += coeff[1] * (omega[0]*uhx[0] + omega[1]*uhx[1] + omega[2]*uhx[2]);
        // compute the term A_y
        if (coeff[2]  <= 0)
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<=N_x )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+N_x1];
            val[5] = sol_curr[i+2 * N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+3*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[2] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[3] = y_coord[i];
            coord[4] = y_coord[i+N_x1];
            coord[5] = y_coord[i+2 * N_x1];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>=N_x1)&&(i<2*N_x1) )
            {
              val[0] = val[1] = val[2] = sol_curr[i-N_x1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+N_x1];
              val[5] = sol_curr[i+2*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i+N_x1];
              coord[1] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[2] = y_coord[i-N_x1];
              coord[3] = y_coord[i];
              coord[4] = y_coord[i+N_x1];
              coord[5] = y_coord[i+2*N_x1];
            }
            else
            {
              // next layer on bottom boundary
              if ( (i>=2*N_x1)&&(i<3*N_x1) )
              {
                val[0] = val[1] = sol_curr[i-2*N_x1];
                val[2] = sol_curr[i-N_x1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+N_x1];
                val[5] = sol_curr[i+2*N_x1];
                coord[0] = 2 * y_coord[i-2*N_x1]-y_coord[i-N_x1];
                coord[1] = y_coord[i-2*N_x1];
                coord[2] = y_coord[i-N_x1];
                coord[3] = y_coord[i];
                coord[4] = y_coord[i+N_x1];
                coord[5] = y_coord[i+2*N_x1];
              }
              else
              {
                // point on top boundary
                if ( i>= N_y*N_x1)
                {
                  val[0] = sol_curr[i-3*N_x1];
                  val[1] = sol_curr[i-2*N_x1];
                  val[2] = sol_curr[i-N_x1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = y_coord[i-3*N_x1];
                  coord[1] = y_coord[i-2*N_x1];
                  coord[2] = y_coord[i-N_x1];
                  coord[3] = y_coord[i];
                  coord[4] = 2 *  y_coord[i] - y_coord[i-N_x1];
                  coord[5] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                }
                else
                {
                  // point next to top boundary
                  if (( i>= (N_y-1)*N_x1)&&(i<N_y*N_x1))
                  {
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i];
                  }
                  // inner point
                  else
                  {
                    // if (i==14382)
                    //  OutPut("herea "<< coeff[2] << endl);
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+N_x1];
                    val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = y_coord[i+2*N_x1];
                    // if (i==14382)
                    //	OutPut("val " << val[0] -1 << endl);
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/3.0-val[3]/2.0 + val[4] - val[5]/6.0)/(coord[3]-coord[2]);
          uhx[1] = (val[1]/6.0-val[2]+val[3]/2.0 + val[4]/3.0)/(coord[2]-coord[1]);
          uhx[2] = (-val[0]/3.0+1.5*val[1]-3*val[2] + 11.0*val[3]/6.0)/(coord[1]-coord[0]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[3] + av[4];
          valy = 3*av[2] - 4*av[3] + av[4];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[1] - 2*av[2] + av[3];
          valy = av[1] - av[3];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[0] - 2*av[1] +av[2];
          valy = av[0] - 4*av[1] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }
          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }
        else
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<=N_x )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+N_x1];
            val[4] = sol_curr[i+2*N_x1];
            val[5] = sol_curr[i+3*N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[2] = y_coord[i];
            coord[3] = y_coord[i+N_x1];
            coord[4] = y_coord[i+2*N_x1];
            coord[5] = y_coord[i+3*N_x1];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>=N_x1)&&(i<2*N_x1) )
            {
              val[0] = val[1] = sol_curr[i-N_x1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+N_x1];
              val[4] = sol_curr[i+2*N_x1];
              val[5] = sol_curr[i+3*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[1] = y_coord[i-N_x1];
              coord[2] = y_coord[i];
              coord[3] = y_coord[i+N_x1];
              coord[4] = y_coord[i+2*N_x1];
              coord[5] = y_coord[i+3*N_x1];
            }
            else
            {
              // point on top boundary
              if ( i>= N_y*N_x1)
              {
                val[0] = sol_curr[i-2*N_x1];
                val[1] = sol_curr[i-N_x1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = y_coord[i-2*N_x1];
                coord[1] = y_coord[i-N_x1];
                coord[2] = y_coord[i];
                coord[3] = 2 *  y_coord[i] - y_coord[i-N_x1];
                coord[4] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                coord[5] = 2 *  y_coord[i] - y_coord[i-3*N_x1];
              }
              else
              {
                // point next to top boundary
                if (( i>= (N_y-1)*N_x1)&&(i<N_y*N_x1))
                {
                  val[0] = sol_curr[i-2*N_x1];
                  val[1] = sol_curr[i-N_x1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+N_x1];
                  coord[0] = y_coord[i-2*N_x1];
                  coord[1] = y_coord[i-N_x1];
                  coord[2] = y_coord[i];
                  coord[3] = y_coord[i+N_x1];
                  coord[4] = 2 *  y_coord[i+N_x1] - y_coord[i];
                  coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i-N_x1];
                }
                else
                {
                  // point in next layer to top boundary
                  if (( i>= (N_y-2)*N_x1)&&(i<(N_y-1)*N_x1))
                  {
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = 2 *  y_coord[i+2*N_x1] - y_coord[i+N_x1];
                  }
                  // inner point
                  else
                  {
                    // if (i==14382)
                    //  OutPut("here "<< coeff[2] << endl);
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = sol_curr[i+2*N_x1];
                    val[5] = sol_curr[i+3*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = y_coord[i+3*N_x1];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/2.0+val[1] - val[0]/6.0 - val[3]/3.0)/(-coord[2]+coord[1]);
          uhx[1] = (-val[3] +val[2]/2.0+val[1]/3.0 + val[4]/6.0)/(-coord[3]+coord[2]);
          uhx[2] = (-val[5]/3.0+1.5*val[4]-3*val[3] + 11.0*val[2]/6.0)/(-coord[4]+coord[3]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[1] + av[0];
          valy = 3*av[2] - 4*av[1] + av[0];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[3] - 2*av[2] + av[1];
          valy = av[3] - av[1];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[4] - 2*av[3] +av[2];
          valy = av[4] - 4*av[3] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }
          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }                                         // end coeff[2] <=0
        //OutPut(value << endl);
        current_stage_fdm[i] += coeff[2] * (omega[0]*uhx[0] + omega[1]*uhx[1] + omega[2]*uhx[2]);
        break;
    }
    // set Dirichlet boundary conditions
    // on bottom
    if (i<= N_x)
    {
      bound_cond(0, x_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }
    // inflow from the left x = x_min (left)
    if ( i%(N_x1)==0 )
    {
      bound_cond(3, 1.0-y_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }
    // on top
    if (i> (N_x1)*N_y)
    {
      bound_cond(2, 1.0-x_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }
    // on the right hand side
    if ( (i+1)%(N_x1)==0 )
    {
      bound_cond(1, y_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }

  }
  // reset time
  TDatabase::TimeDB->CURRENTTIME -= (TDatabase::TimeDB->RK_c[current_stage] - 1)* tau;
}
#endif

#ifdef __3D__
void ComputeStages_FDM3D(int dim, CoeffFct3D *Coeffs, BoundCondFunct3D *bound_cond,
BoundValueFunct3D *bound_val,
double *sol, double **stages, int current_stage, int N_x, int N_y, int N_z,
int *dof_conversion, double *x_coord, double *y_coord, double *z_coord)

{
  int i, j, N_x1, N_y1,N_z1, N3, N1_[3];;
  double tau, hx_i, hx_i1, hy_i, hy_i1,hz_i, hz_i1, valx, valy, valz, val[6], coord[6], value;
  double *sol_curr, coeff_array[20], *coeff, *current_stage_fdm, *sol_help;
  double d1, d2, d3, d4, uhx[3], omega[3], alpha[3], d[3], beta[3], av[5], c_e = 1e-6;
  double *coordinates[3];
  int *offset_ = NULL, *offset1_;
  BoundCond cond;

  // set time of the stage, in main already new time level
  tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  TDatabase::TimeDB->CURRENTTIME += (TDatabase::TimeDB->RK_c[current_stage] - 1)* tau;
  N_x1 = N_x+1;
  N_y1 = N_y+1;
  N_z1 = N_z+1;
  N3 = (N_x1)*(N_y1)*(N_z1);
  // allocate arrays
  sol_curr = stages[5];
  coeff = coeff_array;

  // compute linear combination of current stages
  // numeration of dof corresponding to fe space
  current_stage_fdm = stages[current_stage];
  memset(current_stage_fdm, 0, N3 * SizeOfDouble);

  // initialize current stage
  memcpy(sol_curr,sol,N3*SizeOfDouble);
  // add previous stages
  for (i=0;i<current_stage;i++)
    Daxpy(N3, tau*TDatabase::TimeDB->RK_A[current_stage][i], stages[i], sol_curr);

  //N1_[0] = N_x1;
  //N1_[1] = N_y1;
  //N1_[2] = N_z1;
  //coordinates[0] = x_coord;
  //coordinates[1] = y_coord;
  //coordinates[2] = z_coord;
  //InitializeConvectiveTermFDM(3, offset_, offset1_, N1_);

  // discretization of pde with FDM
  // loop over all nodes of the FDM grid
  for ( i=0 ; i<N3 ; i++ )
  {
    // coefficients of the equation
    Coeffs(1, &x_coord[i], &y_coord[i], &z_coord[i], NULL, &coeff);
    // compute the coefficients corresponding to the 3d finite-difference-method (9-point)
    if (( i%(N_x1)!=0 )&&( (i+1)%(N_x1)!=0 )&&
      ((i%((N_x1)*(N_y1)))>N_x ) && (((i+N_x1)%((N_x1)*(N_y1)))>N_x) &&
      (i> N_x1 * N_y1-1) && (i< ((N_x1)*(N_y1)*(N_z))))
    {
      // diffusive term wrt x
      val[0] = sol_curr[i-1];
      val[1] = sol_curr[i];
      val[2] = sol_curr[i+1];
      coord[0] = x_coord[i-1];
      coord[1] = x_coord[i];
      coord[2] = x_coord[i+1];
      valx = 2*DividedDifferences(2, val, coord);
      // diffusive term wrt y
      val[0] = sol_curr[i-N_x1];
      val[2] = sol_curr[i+N_x1];
      coord[0] = y_coord[i-N_x1];;
      coord[1] = y_coord[i];
      coord[2] = y_coord[i+N_x1];
      valy = 2*DividedDifferences(2, val, coord);

      // diffusive term wrt z
      val[0] = sol_curr[i-(N_x1)*(N_y1)];
      val[2] = sol_curr[i+(N_x1)*(N_y1)];
      coord[0] = z_coord[i-(N_x1)*(N_y1)];
      coord[1] = z_coord[i];
      coord[2] = z_coord[i+(N_x1)*(N_y1)];
      valz = 2*DividedDifferences(2, val, coord);

      // add also reactive term and right hand side
      current_stage_fdm[i] = coeff[0]*(valx + valy+valz) - coeff[4]*sol_curr[i] + coeff[5];
    }

    //    ConvectiveTermFDM(3, i,
    //		      coeff+1, sol_curr, current_stage_fdm, coordinates,
    //		      offset_, offset1_);
    //    continue;
    switch(TDatabase::ParamDB->DISCTYPE)
    {
      case GALERKIN:
        exit(1);
        break;
        // simple upwind scheme
      case UPWIND:
        // compute the term A_x
        if (coeff[1] >= 0)
        {
          // not on boundary x = x_min -> in the "x-sense" exists a left neighbour
          if ( i%(N_x1)!=0 )
            current_stage_fdm[i] -= coeff[1]*(sol_curr[i]-sol_curr[i-1])/(x_coord[i]-x_coord[i-1]);
        }
        else
        {
          // not on boundary x = x_max -> in the "x-sense" exists a right neighbour
          if ( (i+1)%(N_x1)!=0 )
            current_stage_fdm[i] -= coeff[1]*(sol_curr[i+1]-sol_curr[i])/(x_coord[i+1]-x_coord[i]);
        }

        // compute the term A_y
        if (coeff[2] >= 0)
        {
          // not on boundary y = y_min -> in the "y-sense" exists a left neighbour
          if ((i%((N_x1)*(N_y1)))>N_x )
            current_stage_fdm[i] -= coeff[2]*(sol_curr[i]-sol_curr[i-(N_x1)])/(y_coord[i]-y_coord[i-(N_x1)]);
        }
        else
        {
          // not on boundary y = y_max -> in the "y-sense" exists a right neighbour
          if (((i+N_x1)%((N_x1)*(N_y1)))>N_x)
            current_stage_fdm[i] -= coeff[2]*(sol_curr[i+(N_x1)]-sol_curr[i])/(y_coord[i+(N_x1)]-y_coord[i]);
        }
        if (coeff[3] >= 0)
        {
          // not on boundary z = z_min -> in the "z-sense" exists a left neighbour
          if (i> N_x1 * N_y)
            current_stage_fdm[i] -= coeff[3]*(sol_curr[i]-sol_curr[i-(N_x1)*(N_y1)])/(z_coord[i]-z_coord[i-(N_x1)*(N_y1)]);
        }
        else
        {
          // not on boundary z = z_max -> in the "z-sense" exists a right neighbour
          if (i< ((N_x1)*(N_y1)*(N_z)))
            current_stage_fdm[i] -= coeff[3]*(sol_curr[i+(N_x1)*(N_y1)]-sol_curr[i])/(z_coord[i+(N_x1)*(N_y1)]-z_coord[i]);
        }
        break;
      case ENO_3:
        // compute the term A_x
        if (coeff[1] >  0)
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+1];
            val[5] = sol_curr[i+2];
            coord[0] = 2 * x_coord[i]-x_coord[i+3];
            coord[1] = 2 * x_coord[i]-x_coord[i+2];
            coord[2] = 2 * x_coord[i]-x_coord[i+1];
            coord[3] = x_coord[i];
            coord[4] = x_coord[i+1];
            coord[5] = x_coord[i+2];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )                //gleich zu Zeile 1317
            {
              val[0] = val[1] = val[2] = sol_curr[i-1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+1];
              val[5] = sol_curr[i+2];
              coord[0] = 2 * x_coord[i-1]-x_coord[i+1];
              coord[1] = 2 * x_coord[i-1]-x_coord[i];
              coord[2] = x_coord[i-1];
              coord[3] = x_coord[i];
              coord[4] = x_coord[i+1];
              coord[5] = x_coord[i+2];
            }
            else
            {
              // next layer on left boundary
              if ( (i-2)%(N_x1)==0 )              //gleich zu zeile 1300 (i-2) MUSS IN 2d auch verbessert werden
              {
                val[0] = val[1] = sol_curr[i-2];
                val[2] = sol_curr[i-1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+1];
                val[5] = sol_curr[i+2];
                coord[0] = 2 * x_coord[i-2]-x_coord[i-1];
                coord[1] = x_coord[i-2];
                coord[2] = x_coord[i-1];
                coord[3] = x_coord[i];
                coord[4] = x_coord[i+1];
                coord[5] = x_coord[i+2];
              }
              else
              {
                // point on right boundary
                if ( (i+1)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-3];
                  val[1] = sol_curr[i-2];
                  val[2] = sol_curr[i-1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = x_coord[i-3];
                  coord[1] = x_coord[i-2];
                  coord[2] = x_coord[i-1];
                  coord[3] = x_coord[i];
                  coord[4] = 2 *  x_coord[i] - x_coord[i-1];
                  coord[5] = 2 *  x_coord[i] - x_coord[i-2];
                }
                else
                {
                  // point next to right boundary
                  if ( (i+2)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+1];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = 2 *  x_coord[i+1] - x_coord[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+1];
                    val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = x_coord[i+2];
                    /*sol_help = sol_curr+i-3;
                                val[0] = *sol_help;
                                val[1] = *(sol_help+1);
                                val[2] = *(sol_help+2);
                                val[3] = *(sol_help+3);
                                val[4] = *(sol_help+4);
                                val[5] = *(sol_help+5);
                                sol_help = x_coord+i-3;
                                coord[0] = *sol_help;
                                coord[1] = *(sol_help+1);
                                coord[2] = *(sol_help+2);
                    coord[3] = *(sol_help+3);
                    coord[4] = *(sol_help+4);
                    coord[5] = *(sol_help+5);*/
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+2,coord+2);
          d2 = DividedDifferences(2,val+1,coord+1);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              //value = DividedDifferences(1,val+1,coord+1);
              value =  (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d3 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3a ");
            }
            else
            {
              //value = DividedDifferences(1,val+2,coord+2);
              valx = coord[3]-coord[2];
              value = (val[3]-val[2])/valx;
              value += d1 * valx;
              value += d4 * valx*(coord[3]-coord[4]);
              //OutPut("a2 ");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              //value = DividedDifferences(1,val,coord);
              value = (val[1]-val[0])/(coord[1]-coord[0]);
              value += DividedDifferences(2,val,coord) * (coord[3]-coord[0] + coord[3]-coord[1]);
              valx = coord[3]-coord[0];
              valy = coord[3]-coord[1];
              valz = coord[3]-coord[2];
              //value += d3 * ((coord[3]-coord[1])*(coord[3]-coord[2])
              //  + (coord[3]-coord[0])*(coord[3]-coord[2]) + (coord[3]-coord[0])*(coord[3]-coord[1]));
              value += d3 * (valy * valz + valx * valz + valx * valy);
              //OutPut("a4 ");
            }
            else
            {
              //value = DividedDifferences(1,val+1,coord+1);
              value = (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3b ");
              // OutPut("a3 " << value << "  " << DividedDifferences(1,val+3,coord+3)<<":");
            }
          }
        }
        else
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+1];
            val[4] = sol_curr[i+2];
            val[5] = sol_curr[i+3];
            coord[0] = 2 * x_coord[i]-x_coord[i+2];
            coord[1] = 2 * x_coord[i]-x_coord[i+1];
            coord[2] = x_coord[i];
            coord[3] = x_coord[i+1];
            coord[4] = x_coord[i+2];
            coord[5] = x_coord[i+3];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )
            {
              val[0] = val[1] = sol_curr[i-1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+1];
              val[4] = sol_curr[i+2];
              val[5] = sol_curr[i+3];
              coord[0] = 2 * x_coord[i-1]-x_coord[i];
              coord[1] = x_coord[i-1];
              coord[2] = x_coord[i];
              coord[3] = x_coord[i+1];
              coord[4] = x_coord[i+2];
              coord[5] = x_coord[i+3];
            }
            else
            {
              // point on right boundary
              if ( (i+1)%(N_x1)==0 )
              {
                val[0] = sol_curr[i-2];
                val[1] = sol_curr[i-1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = x_coord[i-2];
                coord[1] = x_coord[i-1];
                coord[2] = x_coord[i];
                coord[3] = 2 *  x_coord[i] - x_coord[i-1];
                coord[4] = 2 *  x_coord[i] - x_coord[i-2];
                coord[5] = 2 *  x_coord[i] - x_coord[i-3];
              }
              else
              {
                // point next to right boundary
                if ( (i+2)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-2];
                  val[1] = sol_curr[i-1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+1];
                  coord[0] = x_coord[i-2];
                  coord[1] = x_coord[i-1];
                  coord[2] = x_coord[i];
                  coord[3] = x_coord[i+1];
                  coord[4] = 2 *  x_coord[i+1] - x_coord[i];
                  coord[5] = 2 *  x_coord[i+1] - x_coord[i-1];
                }
                else
                {
                  // point in next layer to right boundary
                  if ( (i+3)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = 2 *  x_coord[i+2] - x_coord[i+1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = sol_curr[i+2];
                    val[5] = sol_curr[i+3];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = x_coord[i+3];
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+1,coord+1);
          d2 = DividedDifferences(2,val+2,coord+2);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[2]-coord[0] + coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[0])*(coord[2]-coord[1]);
              //OutPut("A3 ");
            }
            else
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d4 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              //OutPut("A2a ");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              valx = coord[2]-coord[1];
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * valx;
              value += d3 * valx*(coord[2]-coord[3]);
              //OutPut("A2b ");
            }
            else
            {
              valx = coord[2]-coord[3];
              value = DividedDifferences(1,val+2,coord+2);
              value += d2 * valx;
              value += d4 * valx*(coord[2]-coord[4]);
              //comparison with backward difference in node i
              //OutPut("A1b ");
              //OutPut("A1 " << value <<"  " <<   DividedDifferences(1,val+1,coord+1) << ":");
            }
          }
        }                                         // end coeff[1] <=0
        //OutPut(value << " ");
        current_stage_fdm[i] -= coeff[1] * value;
        // compute the term A_y
        if (coeff[2]  > 0)
        {
          // prepare vectors for divided differences
          // point on front boundary
          if (( i % ((N_x1)*(N_y1)) ) <=N_x)      //richtig
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+N_x1];
            val[5] = sol_curr[i+2 * N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+3*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[2] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[3] = y_coord[i];
            coord[4] = y_coord[i+N_x1];
            coord[5] = y_coord[i+2 * N_x1];
            //OutPut("front " << i << " " << y_coord[i] << endl);
          }
          else
          {
            //OutPut(" a " << i << " " << N_x1 << " " << N_y1 << " " <<  i % ((N_x1)*(N_y1)) << endl);
            // point next to front boundary
                                                  //richtig
            if ( (i-N_x1) % ((N_x1)*(N_y1))  <=N_x)
            {
              val[0] = val[1] = val[2] = sol_curr[i-N_x1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+N_x1];
              val[5] = sol_curr[i+2*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i+N_x1];
              coord[1] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[2] = y_coord[i-N_x1];
              coord[3] = y_coord[i];
              coord[4] = y_coord[i+N_x1];
              coord[5] = y_coord[i+2*N_x1];
              //OutPut("front +1 " << i << " " << y_coord[i] <<endl);
            }
            else
            {
              // next layer on front boundary
                                                  //GLEICHe Zeile wie Z 1646 (i-2*N_x1)?
                                                  //richtig
              if (  (i-2*N_x1) % ((N_x1)*(N_y1))  <=N_x )
              {
                val[0] = val[1] = sol_curr[i-2*N_x1];
                val[2] = sol_curr[i-N_x1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+N_x1];
                val[5] = sol_curr[i+2*N_x1];
                coord[0] = 2 * y_coord[i-2*N_x1]-y_coord[i-N_x1];
                coord[1] = y_coord[i-2*N_x1];
                coord[2] = y_coord[i-N_x1];
                coord[3] = y_coord[i];
                coord[4] = y_coord[i+N_x1];
                coord[5] = y_coord[i+2*N_x1];
                //OutPut("front +2  " << i << " " << y_coord[i] <<endl);
              }
              else
              {
                // point on back boundary
                                                  //richtig
                if ( (i+N_x1) % ((N_x1)*(N_y1))  <=N_x)
                {
                  val[0] = sol_curr[i-3*N_x1];
                  val[1] = sol_curr[i-2*N_x1];
                  val[2] = sol_curr[i-N_x1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = y_coord[i-3*N_x1];
                  coord[1] = y_coord[i-2*N_x1];
                  coord[2] = y_coord[i-N_x1];
                  coord[3] = y_coord[i];
                  coord[4] = 2 *  y_coord[i] - y_coord[i-N_x1];
                  coord[5] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                  //OutPut("back  " << i << " " << y_coord[i] <<endl);
                }
                else
                {
                  // point next to back boundary
                                                  //richtig
                  if (( (i+2*N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
                  {
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i];
                    //OutPut("back -1  " << i <<  " " << y_coord[i] << endl);
                  }
                  // inner point              richtig
                  else
                  {
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+N_x1];
                    val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = y_coord[i+2*N_x1];
                    //OutPut("inner  " << i << " " << y_coord[i] << endl);
                  }                               // inner point
                }                                 // point next to back boundary
              }                                   // point on back boundary
            }                                     // next layer on front boundary
          }                                       // point next on front boundary
          // end  prepare vectors for divided differences
          // compute approximation of convective term
          /*d1 = DividedDifferences(2,val+2,coord+2);
          d2 = DividedDifferences(2,val+1,coord+1);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
          value = DividedDifferences(1,val+1,coord+1);
          value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
          value += d3 * (coord[3]-coord[1])*(coord[3]-coord[2]);
          //OutPut("a3a ");
          //OutPut("a3 " << value << " ");
          }
          else
          {
          value = DividedDifferences(1,val+2,coord+2);
          value += d1 * (coord[3]-coord[2]);
          //Habe hier nichts verändert von Zeile 1735 bis 1790
          value += d4 * (coord[3]-coord[2])*(coord[3]-coord[4]);
          //OutPut("a2 ");
          }
          }
          else
          {
          //d3 = DividedDifferences(3,val,coord);
          //d4 = DividedDifferences(3,val+1,coord+1);
          DividedDifferences_3_2(val,coord,d);
          d3 = d[0];
          d4 = d[1];
          if   (fabs(d3) < fabs(d4))
          {
          value = DividedDifferences(1,val,coord);
          value += DividedDifferences(2,val,coord) * (coord[3]-coord[0] + coord[3]-coord[1]);
          value += d3 * ((coord[3]-coord[1])*(coord[3]-coord[2])
          + (coord[3]-coord[0])*(coord[3]-coord[2]) + (coord[3]-coord[0])*(coord[3]-coord[1]));
          // OutPut("a4 ");
          //OutPut("a4 " << coeff[2]  << "  " << DividedDifferences(2,val+2,coord+2) << "  "
          // << DividedDifferences(2,val+1,coord+1) << "  " << DividedDifferences(3,val,coord)
          // << "  " << DividedDifferences(3,val+1,coord+1)
          //<< endl);
          }
          else
          {
          value = DividedDifferences(1,val+1,coord+1);
          value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
          value += d4 * (coord[3]-coord[1])*(coord[3]-coord[2]);
          // comparison with backward difference in node i
          // OutPut("a3b ");
          //OutPut("a3 " << value << "  " << DividedDifferences(1,val+3,coord+3)<<":");
          }
          } */

          // compute approximation of convective term
          d1 = DividedDifferences(2,val+2,coord+2);
          d2 = DividedDifferences(2,val+1,coord+1);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              //value = DividedDifferences(1,val+1,coord+1);
              value =  (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d3 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3a ");
            }
            else
            {
              //value = DividedDifferences(1,val+2,coord+2);
              valx = coord[3]-coord[2];
              value = (val[3]-val[2])/valx;
              value += d1 * valx;
              value += d4 * valx*(coord[3]-coord[4]);
              //OutPut("a2 ");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              //value = DividedDifferences(1,val,coord);
              value = (val[1]-val[0])/(coord[1]-coord[0]);
              value += DividedDifferences(2,val,coord) * (coord[3]-coord[0] + coord[3]-coord[1]);
              valx = coord[3]-coord[0];
              valy = coord[3]-coord[1];
              valz = coord[3]-coord[2];
              //value += d3 * ((coord[3]-coord[1])*(coord[3]-coord[2])
              //  + (coord[3]-coord[0])*(coord[3]-coord[2]) + (coord[3]-coord[0])*(coord[3]-coord[1]));
              value += d3 * (valy * valz + valx * valz + valx * valy);
              //OutPut("a4 ");
            }
            else
            {
              //value = DividedDifferences(1,val+1,coord+1);
              value = (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3b ");
              // OutPut("a3 " << value << "  " << DividedDifferences(1,val+3,coord+3)<<":");
            }
          }
        }
        else
        {
          // prepare vectors for divided differences
          // point on front boundary
          if (( i % ((N_x1)*(N_y1)) ) <= N_x )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+N_x1];
            val[4] = sol_curr[i+2*N_x1];
            val[5] = sol_curr[i+3*N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[2] = y_coord[i];
            coord[3] = y_coord[i+N_x1];
            coord[4] = y_coord[i+2*N_x1];
            coord[5] = y_coord[i+3*N_x1];
          }
          else
          {
            // point next to front boundary
            if (((i-N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
            {
              val[0] = val[1] = sol_curr[i-N_x1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+N_x1];
              val[4] = sol_curr[i+2*N_x1];
              val[5] = sol_curr[i+3*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[1] = y_coord[i-N_x1];
              coord[2] = y_coord[i];
              coord[3] = y_coord[i+N_x1];
              coord[4] = y_coord[i+2*N_x1];
              coord[5] = y_coord[i+3*N_x1];
            }
            else
            {
              // point on back boundary
              if (((i+N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
              {
                val[0] = sol_curr[i-2*N_x1];
                val[1] = sol_curr[i-N_x1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = y_coord[i-2*N_x1];
                coord[1] = y_coord[i-N_x1];
                coord[2] = y_coord[i];
                coord[3] = 2 *  y_coord[i] - y_coord[i-N_x1];
                coord[4] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                coord[5] = 2 *  y_coord[i] - y_coord[i-3*N_x1];
              }
              else
              {
                // point next to back boundary
                if (((i+2*N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
                {
                  val[0] = sol_curr[i-2*N_x1];
                  val[1] = sol_curr[i-N_x1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+N_x1];
                  coord[0] = y_coord[i-2*N_x1];
                  coord[1] = y_coord[i-N_x1];
                  coord[2] = y_coord[i];
                  coord[3] = y_coord[i+N_x1];
                  coord[4] = 2 *  y_coord[i+N_x1] - y_coord[i];
                  coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i-N_x1];
                }
                else
                {
                  // point in next layer to back boundary
                  if (((i+3*N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
                  {
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = 2 *  y_coord[i+2*N_x1] - y_coord[i+N_x1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = sol_curr[i+2*N_x1];
                    val[5] = sol_curr[i+3*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = y_coord[i+3*N_x1];
                  }
                }
              }
            }
          }

          // compute approximation of convective term
          d1 = DividedDifferences(2,val+1,coord+1);
          d2 = DividedDifferences(2,val+2,coord+2);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[2]-coord[0] + coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[0])*(coord[2]-coord[1]);
              // OutPut("A3 ");
            }
            else
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d4 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              // OutPut("A2a ");
              //OutPut("A1 " << value << "  " << DividedDifferences(1,val+1,coord+1)<<":");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              // OutPut("A2b ");
            }
            else
            {
              value = DividedDifferences(1,val+2,coord+2);
              value += d2 * (coord[2]-coord[3]);
              value += d4 * (coord[2]-coord[3])*(coord[2]-coord[4]);
              // comparison with backward difference in node i
              // OutPut("A1b ");
              //OutPut("A1 " << value << "  " << DividedDifferences(1,val+1,coord+1)<<":");
            }
          }
        }                                         // end coeff[2] <=0
        //OutPut(value << " " );
        current_stage_fdm[i] -= coeff[2] * value;
        //OutPut(i << " " <<  current_stage_fdm[i] << " " << value << endl);
        // end  prepare vectors for divided///
        // compute approximation of convective term
        // compute the term A_z
        if (coeff[3]  > 0)
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<= (N_x1)*(N_y))
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+(N_x1)*(N_y1)];
            val[5] = sol_curr[i+2 * (N_x1)*(N_y1)];
            coord[0] = 2 * z_coord[i]-z_coord[i+3*(N_x1)*(N_y1)];
            coord[1] = 2 * z_coord[i]-z_coord[i+2*(N_x1)*(N_y1)];
            coord[2] = 2 * z_coord[i]-z_coord[i+(N_x1)*(N_y1)];
            coord[3] = z_coord[i];
            coord[4] = z_coord[i+(N_x1)*(N_y1)];
            coord[5] = z_coord[i+2 *(N_x1)*(N_y1)];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>= (N_x1)*(N_y1)&&(i<2*(N_x1)*(N_y1) )))
            {
              val[0] = val[1] = val[2] = sol_curr[i-(N_x1)*(N_y1)];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+(N_x1)*(N_y1)];
              val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
              coord[0] = 2 * z_coord[i-(N_x1)*(N_y1)]-z_coord[i+(N_x1)*(N_y1)];
              coord[1] = 2 * z_coord[i-(N_x1)*(N_y1)]-z_coord[i];
              coord[2] = z_coord[i-(N_x1)*(N_y1)];
              coord[3] = z_coord[i];
              coord[4] = z_coord[i+(N_x1)*(N_y1)];
              coord[5] = z_coord[i+2*(N_x1)*(N_y1)];
            }
            else
            {
              // next layer on bottom boundary
              if ( (i>=2*(N_x1)*(N_y1)&&(i<3*(N_x1)*(N_y1))))
              {
                val[0] = val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                val[2] = sol_curr[i-(N_x1)*(N_y1)];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+(N_x1)*(N_y1)];
                val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
                coord[0] = 2 * z_coord[i-2*(N_x1)*(N_y1)]-z_coord[i-(N_x1)*(N_y1)];
                coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                coord[2] = z_coord[i-(N_x1)*(N_y1)];
                coord[3] = z_coord[i];
                coord[4] = z_coord[i+(N_x1)*(N_y1)];
                coord[5] = z_coord[i+2*(N_x1)*(N_y1)];
              }
              else
              {
                // point on top boundary
                if ( i>=N_z*N_y1*N_x1)
                {
                  val[0] = sol_curr[i-3*(N_x1)*(N_y1)];
                  val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                  val[2] = sol_curr[i-(N_x1)*(N_y1)];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = z_coord[i-3*(N_x1)*(N_y1)];
                  coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                  coord[2] = z_coord[i-(N_x1)*(N_y1)];
                  coord[3] = z_coord[i];
                  coord[4] = 2 *  z_coord[i] - z_coord[i-(N_x1)*(N_y1)];
                  coord[5] = 2 *  z_coord[i] - z_coord[i-2*(N_x1)*(N_y1)];
                }
                else
                {
                  // point next to top boundary
                  if (( i>= (N_z-1)*N_x1*N_y1)&&(i<N_z*N_x1*N_y1))
                  {
                    val[0] = sol_curr[i-3*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[2] = sol_curr[i-(N_x1)*(N_y1)];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-3*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[2] = z_coord[i-(N_x1)*(N_y1)];
                    coord[3] = z_coord[i];
                    coord[4] = z_coord[i+(N_x1)*(N_y1)];
                    coord[5] = 2 *  z_coord[i+N_x1*(N_y1)] - z_coord[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-3*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[2] = sol_curr[i-(N_x1)*(N_y1)];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+(N_x1)*(N_y1)];
                    val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-3*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[2] = z_coord[i-(N_x1)*(N_y1)];
                    coord[3] = z_coord[i];
                    coord[4] = z_coord[i+(N_x1)*(N_y1)];
                    coord[5] = z_coord[i+2*(N_x1)*(N_y1)];
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+2,coord+2);
          d2 = DividedDifferences(2,val+1,coord+1);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d3 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              //OutPut("a3a ");
              //OutPut("a3 " << value << " ");
            }
            else
            {
              value = DividedDifferences(1,val+2,coord+2);
              value += d1 * (coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[2])*(coord[3]-coord[4]);
              //OutPut("a2 ");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[3]-coord[0] + coord[3]-coord[1]);
              value += d3 * ((coord[3]-coord[1])*(coord[3]-coord[2])
                + (coord[3]-coord[0])*(coord[3]-coord[2]) + (coord[3]-coord[0])*(coord[3]-coord[1]));
              // OutPut("a4 ");
              //OutPut("a4 " << coeff[2]  << "  " << DividedDifferences(2,val+2,coord+2) << "  "
              // << DividedDifferences(2,val+1,coord+1) << "  " << DividedDifferences(3,val,coord)
              // << "  " << DividedDifferences(3,val+1,coord+1)
              //<< endl);
            }
            else
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[1])*(coord[3]-coord[2]);
              // comparison with backward difference in node i
              // OutPut("a3b ");
              //OutPut("a3 " << value << "  " << DividedDifferences(1,val+3,coord+3)<<":");
            }
          }
        }
        else
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<=(N_y)*(N_x1) )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+(N_x1)*(N_y1)];
            val[4] = sol_curr[i+2*(N_x1)*(N_y1)];
            val[5] = sol_curr[i+3*(N_x1)*(N_y1)];
            coord[0] = 2 * z_coord[i]-z_coord[i+2*(N_x1)*(N_y1)];
            coord[1] = 2 * z_coord[i]-z_coord[i+(N_x1)*(N_y1)];
            coord[2] = z_coord[i];
            coord[3] = z_coord[i+(N_x1)*(N_y1)];
            coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
            coord[5] = z_coord[i+3*(N_x1)*(N_y1)];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>=(N_x1)*(N_y1))&&(i<2*(N_x1)*(N_y1)) )
            {
              val[0] = val[1] = sol_curr[i-(N_x1)*(N_y1)];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+(N_x1)*(N_y1)];
              val[4] = sol_curr[i+2*(N_x1)*(N_y1)];
              val[5] = sol_curr[i+3*(N_x1)*(N_y1)];
              coord[0] = 2 * z_coord[i-(N_x1)*(N_y1)]-z_coord[i];
              coord[1] = z_coord[i-(N_x1)*(N_y1)];
              coord[2] = z_coord[i];
              coord[3] = z_coord[i+(N_x1)*(N_y1)];
              coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
              coord[5] = z_coord[i+3*(N_x1)*(N_y1)];
            }
            else
            {
              // point on top boundary
              if ( i>= N_z*N_y1*N_x1)
              {
                val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                val[1] = sol_curr[i-(N_x1)*(N_y1)];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                coord[1] = z_coord[i-(N_x1)*(N_y1)];
                coord[2] = z_coord[i];
                coord[3] = 2 *  z_coord[i] - z_coord[i-(N_x1)*(N_y1)];
                coord[4] = 2 *  z_coord[i] - z_coord[i-2*(N_x1)*(N_y1)];
                coord[5] = 2 *  z_coord[i] - z_coord[i-3*(N_x1)*(N_y1)];
              }
              else
              {
                // point next to top boundary
                if (( i>= (N_z-1)*(N_y1)*N_x1)&&(i<(N_z)*(N_y1)*N_x1))
                {
                  val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                  val[1] = sol_curr[i-(N_x1)*(N_y1)];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+(N_x1)*(N_y1)];
                  coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                  coord[1] = z_coord[i-(N_x1)*(N_y1)];
                  coord[2] = z_coord[i];
                  coord[3] = z_coord[i+(N_x1)*(N_y1)];
                  coord[4] = 2 *  z_coord[i+(N_x1)*(N_y1)] - z_coord[i];
                  coord[5] = 2 *  z_coord[i+(N_x1)*(N_y1)] - z_coord[i-(N_x1)*(N_y1)];
                }
                else
                {
                  // point in next layer to top boundary
                  if (( i>= (N_z-2)*(N_y1)*(N_x1)&&(i<(N_z-1)*(N_y1)*N_x1)))
                  {
                    val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-(N_x1)*(N_y1)];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+(N_x1)*(N_y1)];
                    val[4] = val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-(N_x1)*(N_y1)];
                    coord[2] = z_coord[i];
                    coord[3] = z_coord[i+(N_x1)*(N_y1)];
                    coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
                    coord[5] = 2 *  z_coord[i+2*(N_x1)*(N_y1)] - z_coord[i+(N_x1)*(N_y1)];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-(N_x1)*(N_y1)];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+(N_x1)*(N_y1)];
                    val[4] = sol_curr[i+2*(N_x1)*(N_y1)];
                    val[5] = sol_curr[i+3*(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-(N_x1)*(N_y1)];
                    coord[2] = z_coord[i];
                    coord[3] = z_coord[i+(N_x1)*(N_y1)];
                    coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
                    coord[5] = z_coord[i+3*(N_x1)*(N_y1)];
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences

          // compute approximation of convective term
          d1 = DividedDifferences(2,val+1,coord+1);
          d2 = DividedDifferences(2,val+2,coord+2);
          if (fabs(d1) < fabs(d2))
          {
            //d3 = DividedDifferences(3,val,coord);
            //d4 = DividedDifferences(3,val+1,coord+1);
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[2]-coord[0] + coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[0])*(coord[2]-coord[1]);
              // OutPut("A3 ");
            }
            else
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d4 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              // OutPut("A2a ");
              //OutPut("A1 " << value << "  " << DividedDifferences(1,val+1,coord+1)<<":");
            }
          }
          else
          {
            //d3 = DividedDifferences(3,val+1,coord+1);
            //d4 = DividedDifferences(3,val+2,coord+2);
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[1])*(coord[2]-coord[3]);
              // OutPut("A2b ");
            }
            else
            {
              value = DividedDifferences(1,val+2,coord+2);
              value += d2 * (coord[2]-coord[3]);
              value += d4 * (coord[2]-coord[3])*(coord[2]-coord[4]);
              // comparison with backward difference in node i
              // OutPut("A1b ");
              //OutPut("A1 " << value << "  " << DividedDifferences(1,val+1,coord+1)<<":");
            }
          }
        }                                         // end coeff[3] <=0
        //OutPut(value << endl);
        current_stage_fdm[i] -= coeff[3] * value;
        break;
      case WENO_5:
        d[0] = 0.3;
        d[1] = 0.6;
        d[2] = 0.1;
        // compute the term A_x
        if (coeff[1] >=  0)
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+1];
            val[5] = sol_curr[i+2];
            coord[0] = 2 * x_coord[i]-x_coord[i+3];
            coord[1] = 2 * x_coord[i]-x_coord[i+2];
            coord[2] = 2 * x_coord[i]-x_coord[i+1];
            coord[3] = x_coord[i];
            coord[4] = x_coord[i+1];
            coord[5] = x_coord[i+2];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )                //gleich zu Zeile 1317
            {
              val[0] = val[1] = val[2] = sol_curr[i-1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+1];
              val[5] = sol_curr[i+2];
              coord[0] = 2 * x_coord[i-1]-x_coord[i+1];
              coord[1] = 2 * x_coord[i-1]-x_coord[i];
              coord[2] = x_coord[i-1];
              coord[3] = x_coord[i];
              coord[4] = x_coord[i+1];
              coord[5] = x_coord[i+2];
            }
            else
            {
              // next layer on left boundary
              if ( (i-2)%(N_x1)==0 )              //gleich zu zeile 1300 (i-2) MUSS IN 2d auch verbessert werden
              {
                val[0] = val[1] = sol_curr[i-2];
                val[2] = sol_curr[i-1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+1];
                val[5] = sol_curr[i+2];
                coord[0] = 2 * x_coord[i-2]-x_coord[i-1];
                coord[1] = x_coord[i-2];
                coord[2] = x_coord[i-1];
                coord[3] = x_coord[i];
                coord[4] = x_coord[i+1];
                coord[5] = x_coord[i+2];
              }
              else
              {
                // point on right boundary
                if ( (i+1)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-3];
                  val[1] = sol_curr[i-2];
                  val[2] = sol_curr[i-1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = x_coord[i-3];
                  coord[1] = x_coord[i-2];
                  coord[2] = x_coord[i-1];
                  coord[3] = x_coord[i];
                  coord[4] = 2 *  x_coord[i] - x_coord[i-1];
                  coord[5] = 2 *  x_coord[i] - x_coord[i-2];
                }
                else
                {
                  // point next to right boundary
                  if ( (i+2)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+1];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = 2 *  x_coord[i+1] - x_coord[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-3];
                    val[1] = sol_curr[i-2];
                    val[2] = sol_curr[i-1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+1];
                    val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-3];
                    coord[1] = x_coord[i-2];
                    coord[2] = x_coord[i-1];
                    coord[3] = x_coord[i];
                    coord[4] = x_coord[i+1];
                    coord[5] = x_coord[i+2];
                    /*sol_help = sol_curr+i-3;
                                val[0] = *sol_help;
                                val[1] = *(sol_help+1);
                                val[2] = *(sol_help+2);
                                val[3] = *(sol_help+3);
                                val[4] = *(sol_help+4);
                                val[5] = *(sol_help+5);
                                sol_help = x_coord+i-3;
                                coord[0] = *sol_help;
                                coord[1] = *(sol_help+1);
                                coord[2] = *(sol_help+2);
                    coord[3] = *(sol_help+3);
                    coord[4] = *(sol_help+4);
                    coord[5] = *(sol_help+5);*/
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/3.0-val[3]/2.0 + val[4] - val[5]/6.0)/(coord[3]-coord[2]);
          uhx[1] = (val[1]/6.0-val[2]+val[3]/2.0 + val[4]/3.0)/(coord[2]-coord[1]);
          uhx[2] = (-val[0]/3.0+1.5*val[1]-3*val[2] + 11.0*val[3]/6.0)/(coord[1]-coord[0]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[3] + av[4];
          valy = 3*av[2] - 4*av[3] + av[4];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[1] - 2*av[2] + av[3];
          valy = av[1] - av[3];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[0] - 2*av[1] +av[2];
          valy = av[0] - 4*av[1] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }

          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }
        else
        {
          // prepare vectors for divided differences
          // point on left boundary
          if ( i%(N_x1)==0 )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+1];
            val[4] = sol_curr[i+2];
            val[5] = sol_curr[i+3];
            coord[0] = 2 * x_coord[i]-x_coord[i+2];
            coord[1] = 2 * x_coord[i]-x_coord[i+1];
            coord[2] = x_coord[i];
            coord[3] = x_coord[i+1];
            coord[4] = x_coord[i+2];
            coord[5] = x_coord[i+3];
          }
          else
          {
            // point next to left boundary
            if ( (i-1)%(N_x1)==0 )
            {
              val[0] = val[1] = sol_curr[i-1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+1];
              val[4] = sol_curr[i+2];
              val[5] = sol_curr[i+3];
              coord[0] = 2 * x_coord[i-1]-x_coord[i];
              coord[1] = x_coord[i-1];
              coord[2] = x_coord[i];
              coord[3] = x_coord[i+1];
              coord[4] = x_coord[i+2];
              coord[5] = x_coord[i+3];
            }
            else
            {
              // point on right boundary
              if ( (i+1)%(N_x1)==0 )
              {
                val[0] = sol_curr[i-2];
                val[1] = sol_curr[i-1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = x_coord[i-2];
                coord[1] = x_coord[i-1];
                coord[2] = x_coord[i];
                coord[3] = 2 *  x_coord[i] - x_coord[i-1];
                coord[4] = 2 *  x_coord[i] - x_coord[i-2];
                coord[5] = 2 *  x_coord[i] - x_coord[i-3];
              }
              else
              {
                // point next to right boundary
                if ( (i+2)%(N_x1)==0 )
                {
                  val[0] = sol_curr[i-2];
                  val[1] = sol_curr[i-1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+1];
                  coord[0] = x_coord[i-2];
                  coord[1] = x_coord[i-1];
                  coord[2] = x_coord[i];
                  coord[3] = x_coord[i+1];
                  coord[4] = 2 *  x_coord[i+1] - x_coord[i];
                  coord[5] = 2 *  x_coord[i+1] - x_coord[i-1];
                }
                else
                {
                  // point in next layer to right boundary
                  if ( (i+3)%(N_x1)==0 )
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = val[5] = sol_curr[i+2];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = 2 *  x_coord[i+2] - x_coord[i+1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2];
                    val[1] = sol_curr[i-1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+1];
                    val[4] = sol_curr[i+2];
                    val[5] = sol_curr[i+3];
                    coord[0] = x_coord[i-2];
                    coord[1] = x_coord[i-1];
                    coord[2] = x_coord[i];
                    coord[3] = x_coord[i+1];
                    coord[4] = x_coord[i+2];
                    coord[5] = x_coord[i+3];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/2.0+val[1] - val[0]/6.0 - val[3]/3.0)/(-coord[2]+coord[1]);
          uhx[1] = (-val[3] +val[2]/2.0+val[1]/3.0 + val[4]/6.0)/(-coord[3]+coord[2]);
          uhx[2] = (-val[5]/3.0+1.5*val[4]-3*val[3] + 11.0*val[2]/6.0)/(-coord[4]+coord[3]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[1] + av[0];
          valy = 3*av[2] - 4*av[1] + av[0];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[3] - 2*av[2] + av[1];
          valy = av[3] - av[1];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[4] - 2*av[3] +av[2];
          valy = av[4] - 4*av[3] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }

          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }                                         // end coeff[1] <=0
        //OutPut(value << " ");
        current_stage_fdm[i] -= coeff[1] * (omega[0]*uhx[0] + omega[1]*uhx[1] + omega[2]*uhx[2]);
        // compute the term A_y
        if (coeff[2]  >= 0)
        {
          // prepare vectors for divided differences
          // point on front boundary
          if (( i % ((N_x1)*(N_y1)) ) <=N_x)      //richtig
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+N_x1];
            val[5] = sol_curr[i+2 * N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+3*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[2] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[3] = y_coord[i];
            coord[4] = y_coord[i+N_x1];
            coord[5] = y_coord[i+2 * N_x1];
            //OutPut("front " << i << " " << y_coord[i] << endl);
          }
          else
          {
            //OutPut(" a " << i << " " << N_x1 << " " << N_y1 << " " <<  i % ((N_x1)*(N_y1)) << endl);
            // point next to front boundary
                                                  //richtig
            if ( (i-N_x1) % ((N_x1)*(N_y1))  <=N_x)
            {
              val[0] = val[1] = val[2] = sol_curr[i-N_x1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+N_x1];
              val[5] = sol_curr[i+2*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i+N_x1];
              coord[1] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[2] = y_coord[i-N_x1];
              coord[3] = y_coord[i];
              coord[4] = y_coord[i+N_x1];
              coord[5] = y_coord[i+2*N_x1];
              //OutPut("front +1 " << i << " " << y_coord[i] <<endl);
            }
            else
            {
              // next layer on front boundary
                                                  //GLEICHe Zeile wie Z 1646 (i-2*N_x1)?
                                                  //richtig
              if (  (i-2*N_x1) % ((N_x1)*(N_y1))  <=N_x )
              {
                val[0] = val[1] = sol_curr[i-2*N_x1];
                val[2] = sol_curr[i-N_x1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+N_x1];
                val[5] = sol_curr[i+2*N_x1];
                coord[0] = 2 * y_coord[i-2*N_x1]-y_coord[i-N_x1];
                coord[1] = y_coord[i-2*N_x1];
                coord[2] = y_coord[i-N_x1];
                coord[3] = y_coord[i];
                coord[4] = y_coord[i+N_x1];
                coord[5] = y_coord[i+2*N_x1];
                //OutPut("front +2  " << i << " " << y_coord[i] <<endl);
              }
              else
              {
                // point on back boundary
                                                  //richtig
                if ( (i+N_x1) % ((N_x1)*(N_y1))  <=N_x)
                {
                  val[0] = sol_curr[i-3*N_x1];
                  val[1] = sol_curr[i-2*N_x1];
                  val[2] = sol_curr[i-N_x1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = y_coord[i-3*N_x1];
                  coord[1] = y_coord[i-2*N_x1];
                  coord[2] = y_coord[i-N_x1];
                  coord[3] = y_coord[i];
                  coord[4] = 2 *  y_coord[i] - y_coord[i-N_x1];
                  coord[5] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                  //OutPut("back  " << i << " " << y_coord[i] <<endl);
                }
                else
                {
                  // point next to back boundary
                                                  //richtig
                  if (( (i+2*N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
                  {
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i];
                    //OutPut("back -1  " << i <<  " " << y_coord[i] << endl);
                  }
                  // inner point              richtig
                  else
                  {
                    val[0] = sol_curr[i-3*N_x1];
                    val[1] = sol_curr[i-2*N_x1];
                    val[2] = sol_curr[i-N_x1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+N_x1];
                    val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-3*N_x1];
                    coord[1] = y_coord[i-2*N_x1];
                    coord[2] = y_coord[i-N_x1];
                    coord[3] = y_coord[i];
                    coord[4] = y_coord[i+N_x1];
                    coord[5] = y_coord[i+2*N_x1];
                    //OutPut("inner  " << i << " " << y_coord[i] << endl);
                  }                               // inner point
                }                                 // point next to back boundary
              }                                   // point on back boundary
            }                                     // next layer on front boundary
          }                                       // point next on front boundary
          // compute values for WENO scheme
          uhx[0] = (-val[2]/3.0-val[3]/2.0 + val[4] - val[5]/6.0)/(coord[3]-coord[2]);
          uhx[1] = (val[1]/6.0-val[2]+val[3]/2.0 + val[4]/3.0)/(coord[2]-coord[1]);
          uhx[2] = (-val[0]/3.0+1.5*val[1]-3*val[2] + 11.0*val[3]/6.0)/(coord[1]-coord[0]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[3] + av[4];
          valy = 3*av[2] - 4*av[3] + av[4];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[1] - 2*av[2] + av[3];
          valy = av[1] - av[3];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[0] - 2*av[1] +av[2];
          valy = av[0] - 4*av[1] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }

          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }
        else
        {
          // prepare vectors for divided differences
          // point on front boundary
          if (( i % ((N_x1)*(N_y1)) ) <= N_x )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+N_x1];
            val[4] = sol_curr[i+2*N_x1];
            val[5] = sol_curr[i+3*N_x1];
            coord[0] = 2 * y_coord[i]-y_coord[i+2*N_x1];
            coord[1] = 2 * y_coord[i]-y_coord[i+N_x1];
            coord[2] = y_coord[i];
            coord[3] = y_coord[i+N_x1];
            coord[4] = y_coord[i+2*N_x1];
            coord[5] = y_coord[i+3*N_x1];
          }
          else
          {
            // point next to front boundary
            if (((i-N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
            {
              val[0] = val[1] = sol_curr[i-N_x1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+N_x1];
              val[4] = sol_curr[i+2*N_x1];
              val[5] = sol_curr[i+3*N_x1];
              coord[0] = 2 * y_coord[i-N_x1]-y_coord[i];
              coord[1] = y_coord[i-N_x1];
              coord[2] = y_coord[i];
              coord[3] = y_coord[i+N_x1];
              coord[4] = y_coord[i+2*N_x1];
              coord[5] = y_coord[i+3*N_x1];
            }
            else
            {
              // point on back boundary
              if (((i+N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
              {
                val[0] = sol_curr[i-2*N_x1];
                val[1] = sol_curr[i-N_x1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = y_coord[i-2*N_x1];
                coord[1] = y_coord[i-N_x1];
                coord[2] = y_coord[i];
                coord[3] = 2 *  y_coord[i] - y_coord[i-N_x1];
                coord[4] = 2 *  y_coord[i] - y_coord[i-2*N_x1];
                coord[5] = 2 *  y_coord[i] - y_coord[i-3*N_x1];
              }
              else
              {
                // point next to back boundary
                if (((i+2*N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
                {
                  val[0] = sol_curr[i-2*N_x1];
                  val[1] = sol_curr[i-N_x1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+N_x1];
                  coord[0] = y_coord[i-2*N_x1];
                  coord[1] = y_coord[i-N_x1];
                  coord[2] = y_coord[i];
                  coord[3] = y_coord[i+N_x1];
                  coord[4] = 2 *  y_coord[i+N_x1] - y_coord[i];
                  coord[5] = 2 *  y_coord[i+N_x1] - y_coord[i-N_x1];
                }
                else
                {
                  // point in next layer to back boundary
                  if (((i+3*N_x1) % ((N_x1)*(N_y1)) ) <=N_x)
                  {
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = val[5] = sol_curr[i+2*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = 2 *  y_coord[i+2*N_x1] - y_coord[i+N_x1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2*N_x1];
                    val[1] = sol_curr[i-N_x1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+N_x1];
                    val[4] = sol_curr[i+2*N_x1];
                    val[5] = sol_curr[i+3*N_x1];
                    coord[0] = y_coord[i-2*N_x1];
                    coord[1] = y_coord[i-N_x1];
                    coord[2] = y_coord[i];
                    coord[3] = y_coord[i+N_x1];
                    coord[4] = y_coord[i+2*N_x1];
                    coord[5] = y_coord[i+3*N_x1];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/2.0+val[1] - val[0]/6.0 - val[3]/3.0)/(-coord[2]+coord[1]);
          uhx[1] = (-val[3] +val[2]/2.0+val[1]/3.0 + val[4]/6.0)/(-coord[3]+coord[2]);
          uhx[2] = (-val[5]/3.0+1.5*val[4]-3*val[3] + 11.0*val[2]/6.0)/(-coord[4]+coord[3]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[1] + av[0];
          valy = 3*av[2] - 4*av[1] + av[0];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[3] - 2*av[2] + av[1];
          valy = av[3] - av[1];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[4] - 2*av[3] +av[2];
          valy = av[4] - 4*av[3] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }

          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }                                         // end coeff[2] <=0
        //OutPut(value << " " );
        current_stage_fdm[i] -= coeff[2] * (omega[0]*uhx[0] + omega[1]*uhx[1] + omega[2]*uhx[2]);
        //OutPut(i << " " <<  current_stage_fdm[i] << " " << value << endl);
        // end  prepare vectors for divided///
        // compute approximation of convective term
        // compute the term A_z
        if (coeff[3]  >= 0)
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<= (N_x1)*(N_y))
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[i+(N_x1)*(N_y1)];
            val[5] = sol_curr[i+2 * (N_x1)*(N_y1)];
            coord[0] = 2 * z_coord[i]-z_coord[i+3*(N_x1)*(N_y1)];
            coord[1] = 2 * z_coord[i]-z_coord[i+2*(N_x1)*(N_y1)];
            coord[2] = 2 * z_coord[i]-z_coord[i+(N_x1)*(N_y1)];
            coord[3] = z_coord[i];
            coord[4] = z_coord[i+(N_x1)*(N_y1)];
            coord[5] = z_coord[i+2 *(N_x1)*(N_y1)];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>= (N_x1)*(N_y1)&&(i<2*(N_x1)*(N_y1) )))
            {
              val[0] = val[1] = val[2] = sol_curr[i-(N_x1)*(N_y1)];
              val[3] = sol_curr[i];
              val[4] = sol_curr[i+(N_x1)*(N_y1)];
              val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
              coord[0] = 2 * z_coord[i-(N_x1)*(N_y1)]-z_coord[i+(N_x1)*(N_y1)];
              coord[1] = 2 * z_coord[i-(N_x1)*(N_y1)]-z_coord[i];
              coord[2] = z_coord[i-(N_x1)*(N_y1)];
              coord[3] = z_coord[i];
              coord[4] = z_coord[i+(N_x1)*(N_y1)];
              coord[5] = z_coord[i+2*(N_x1)*(N_y1)];
            }
            else
            {
              // next layer on bottom boundary
              if ( (i>=2*(N_x1)*(N_y1)&&(i<3*(N_x1)*(N_y1))))
              {
                val[0] = val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                val[2] = sol_curr[i-(N_x1)*(N_y1)];
                val[3] = sol_curr[i];
                val[4] = sol_curr[i+(N_x1)*(N_y1)];
                val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
                coord[0] = 2 * z_coord[i-2*(N_x1)*(N_y1)]-z_coord[i-(N_x1)*(N_y1)];
                coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                coord[2] = z_coord[i-(N_x1)*(N_y1)];
                coord[3] = z_coord[i];
                coord[4] = z_coord[i+(N_x1)*(N_y1)];
                coord[5] = z_coord[i+2*(N_x1)*(N_y1)];
              }
              else
              {
                // point on top boundary
                if ( i>=N_z*N_y1*N_x1)
                {
                  val[0] = sol_curr[i-3*(N_x1)*(N_y1)];
                  val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                  val[2] = sol_curr[i-(N_x1)*(N_y1)];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = z_coord[i-3*(N_x1)*(N_y1)];
                  coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                  coord[2] = z_coord[i-(N_x1)*(N_y1)];
                  coord[3] = z_coord[i];
                  coord[4] = 2 *  z_coord[i] - z_coord[i-(N_x1)*(N_y1)];
                  coord[5] = 2 *  z_coord[i] - z_coord[i-2*(N_x1)*(N_y1)];
                }
                else
                {
                  // point next to top boundary
                  if (( i>= (N_z-1)*N_x1*N_y1)&&(i<N_z*N_x1*N_y1))
                  {
                    val[0] = sol_curr[i-3*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[2] = sol_curr[i-(N_x1)*(N_y1)];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[i+(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-3*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[2] = z_coord[i-(N_x1)*(N_y1)];
                    coord[3] = z_coord[i];
                    coord[4] = z_coord[i+(N_x1)*(N_y1)];
                    coord[5] = 2 *  z_coord[i+N_x1*(N_y1)] - z_coord[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-3*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[2] = sol_curr[i-(N_x1)*(N_y1)];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[i+(N_x1)*(N_y1)];
                    val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-3*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[2] = z_coord[i-(N_x1)*(N_y1)];
                    coord[3] = z_coord[i];
                    coord[4] = z_coord[i+(N_x1)*(N_y1)];
                    coord[5] = z_coord[i+2*(N_x1)*(N_y1)];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/3.0-val[3]/2.0 + val[4] - val[5]/6.0)/(coord[3]-coord[2]);
          uhx[1] = (val[1]/6.0-val[2]+val[3]/2.0 + val[4]/3.0)/(coord[2]-coord[1]);
          uhx[2] = (-val[0]/3.0+1.5*val[1]-3*val[2] + 11.0*val[3]/6.0)/(coord[1]-coord[0]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[3] + av[4];
          valy = 3*av[2] - 4*av[3] + av[4];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[1] - 2*av[2] + av[3];
          valy = av[1] - av[3];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[0] - 2*av[1] +av[2];
          valy = av[0] - 4*av[1] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }
          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }
        else
        {
          // prepare vectors for divided differences
          // point on bottom boundary
          if ( i<=(N_y)*(N_x1) )
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[i+(N_x1)*(N_y1)];
            val[4] = sol_curr[i+2*(N_x1)*(N_y1)];
            val[5] = sol_curr[i+3*(N_x1)*(N_y1)];
            coord[0] = 2 * z_coord[i]-z_coord[i+2*(N_x1)*(N_y1)];
            coord[1] = 2 * z_coord[i]-z_coord[i+(N_x1)*(N_y1)];
            coord[2] = z_coord[i];
            coord[3] = z_coord[i+(N_x1)*(N_y1)];
            coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
            coord[5] = z_coord[i+3*(N_x1)*(N_y1)];
          }
          else
          {
            // point next to bottom boundary
            if ( (i>=(N_x1)*(N_y1))&&(i<2*(N_x1)*(N_y1)) )
            {
              val[0] = val[1] = sol_curr[i-(N_x1)*(N_y1)];
              val[2] = sol_curr[i];
              val[3] = sol_curr[i+(N_x1)*(N_y1)];
              val[4] = sol_curr[i+2*(N_x1)*(N_y1)];
              val[5] = sol_curr[i+3*(N_x1)*(N_y1)];
              coord[0] = 2 * z_coord[i-(N_x1)*(N_y1)]-z_coord[i];
              coord[1] = z_coord[i-(N_x1)*(N_y1)];
              coord[2] = z_coord[i];
              coord[3] = z_coord[i+(N_x1)*(N_y1)];
              coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
              coord[5] = z_coord[i+3*(N_x1)*(N_y1)];
            }
            else
            {
              // point on top boundary
              if ( i>= N_z*N_y1*N_x1)
              {
                val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                val[1] = sol_curr[i-(N_x1)*(N_y1)];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                coord[1] = z_coord[i-(N_x1)*(N_y1)];
                coord[2] = z_coord[i];
                coord[3] = 2 *  z_coord[i] - z_coord[i-(N_x1)*(N_y1)];
                coord[4] = 2 *  z_coord[i] - z_coord[i-2*(N_x1)*(N_y1)];
                coord[5] = 2 *  z_coord[i] - z_coord[i-3*(N_x1)*(N_y1)];
              }
              else
              {
                // point next to top boundary
                if (( i>= (N_z-1)*(N_y1)*N_x1)&&(i<(N_z)*(N_y1)*N_x1))
                {
                  val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                  val[1] = sol_curr[i-(N_x1)*(N_y1)];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[i+(N_x1)*(N_y1)];
                  coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                  coord[1] = z_coord[i-(N_x1)*(N_y1)];
                  coord[2] = z_coord[i];
                  coord[3] = z_coord[i+(N_x1)*(N_y1)];
                  coord[4] = 2 *  z_coord[i+(N_x1)*(N_y1)] - z_coord[i];
                  coord[5] = 2 *  z_coord[i+(N_x1)*(N_y1)] - z_coord[i-(N_x1)*(N_y1)];
                }
                else
                {
                  // point in next layer to top boundary
                  if (( i>= (N_z-2)*(N_y1)*(N_x1)&&(i<(N_z-1)*(N_y1)*N_x1)))
                  {
                    val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-(N_x1)*(N_y1)];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+(N_x1)*(N_y1)];
                    val[4] = val[5] = sol_curr[i+2*(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-(N_x1)*(N_y1)];
                    coord[2] = z_coord[i];
                    coord[3] = z_coord[i+(N_x1)*(N_y1)];
                    coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
                    coord[5] = 2 *  z_coord[i+2*(N_x1)*(N_y1)] - z_coord[i+(N_x1)*(N_y1)];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[i-2*(N_x1)*(N_y1)];
                    val[1] = sol_curr[i-(N_x1)*(N_y1)];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[i+(N_x1)*(N_y1)];
                    val[4] = sol_curr[i+2*(N_x1)*(N_y1)];
                    val[5] = sol_curr[i+3*(N_x1)*(N_y1)];
                    coord[0] = z_coord[i-2*(N_x1)*(N_y1)];
                    coord[1] = z_coord[i-(N_x1)*(N_y1)];
                    coord[2] = z_coord[i];
                    coord[3] = z_coord[i+(N_x1)*(N_y1)];
                    coord[4] = z_coord[i+2*(N_x1)*(N_y1)];
                    coord[5] = z_coord[i+3*(N_x1)*(N_y1)];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/2.0+val[1] - val[0]/6.0 - val[3]/3.0)/(-coord[2]+coord[1]);
          uhx[1] = (-val[3] +val[2]/2.0+val[1]/3.0 + val[4]/6.0)/(-coord[3]+coord[2]);
          uhx[2] = (-val[5]/3.0+1.5*val[4]-3*val[3] + 11.0*val[2]/6.0)/(-coord[4]+coord[3]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[1] + av[0];
          valy = 3*av[2] - 4*av[1] + av[0];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[3] - 2*av[2] + av[1];
          valy = av[3] - av[1];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[4] - 2*av[3] +av[2];
          valy = av[4] - 4*av[3] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }

          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }                                         // end coeff[3] <=0
        //OutPut(value << endl);
        current_stage_fdm[i] -= coeff[3] * (omega[0]*uhx[0] + omega[1]*uhx[1] + omega[2]*uhx[2]);
        break;
    }
    // set Dirichlet boundary conditions
    // inflow from the left x = x_min (left)
    if ( i%(N_x1)==0 )
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }

    // on the right hand side
    if ( (i+1)%(N_x1)==0 )
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }

    // on front
    if ((i%((N_x1)*(N_y1)))<N_x1 )
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }

    // on back
    if (((i+N_x1)%((N_x1)*(N_y1)))<N_x1)
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }

    // on bottom
    if (i< N_x1 * N_y1)
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }

    // on top
    if (i>= ((N_x1)*(N_y1)*(N_z)))
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        current_stage_fdm[i] = 0;
        continue;
      }
    }
  }
  // reset time
  TDatabase::TimeDB->CURRENTTIME -= (TDatabase::TimeDB->RK_c[current_stage] - 1)* tau;
}
#endif
/******************************************************************************/
// ComputeSolutionInNextTime_FDM2D
// compute solution in next time by linear combination of stages
/******************************************************************************/

void ComputeSolutionInNextTime_FDM2D(double *sol, double *sol_e, double **stages, int N_stages,
int N_x, int N_y, int N_z)
{
  int i, N2, j;
  double tau, errors[2];

  tau = TDatabase::TimeDB->TIMESTEPLENGTH;
#ifdef __2D__
  N2 = (N_x+1)*(N_y+1);
#endif
#ifdef __3D__
  N2 = (N_x+1)*(N_y+1)*(N_z+1);
#endif

  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL == 1)
  {
    memcpy(sol_e, sol, N2*SizeOfDouble);
    // compute linear combination of stages
    for (i=0;i<N_stages;i++)
    {
      Daxpy(N2, tau*TDatabase::TimeDB->RK_b[i], stages[i], sol);
      Daxpy(N2, tau*TDatabase::TimeDB->RK_e[i], stages[i], sol_e);
    }
  }
  else
    // compute linear combination of stages
    for (i=0;i<N_stages;i++)
  {
    Daxpy(N2, tau*TDatabase::TimeDB->RK_b[i], stages[i], sol);
  }
}


/******************************************************************************/
// ComputeNewStepsize
// computes new stepsize for next timestep
// mlh
/******************************************************************************/

void ComputeNewStepsize(double *sol, double *sol_e, double *sol_old, int N2, double *err, double *steps_old, int *m, bool *acc, int N_x, int *disc)
{
  double *diff, *sc;
  double val, fac, fac_max, fac_min, tol, atol, rtol, eps, k_I, k_P, k_D, k_E, k_R, h_new;
  int k, order, controller;
  diff = new double[N2];
  sc = new double[N2];
  val = 0.0;
  // Ordnung des Verfahrens
  order = TDatabase::TimeDB->RK_ord;
  fac = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_FAC;
  fac_max = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_FAC_MAX;
  fac_min = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_FAC_MIN;
  tol = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_TOL;
  atol = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_ATOL;
  rtol = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_RTOL;
  eps = fac*tol;
  k_I = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_I/(double)order;
  k_P = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_P/(double)order;
  k_E = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_E/(double)order;
  k_R = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_R/(double)order;
  k_D = TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_D/(double)order;
  *acc = false;

  for(k = 0; k < N2; k++)
  {
    sc[k] = atol + std::max(fabs(sol[k]),fabs(sol_e[k]))*rtol;
    diff[k] = fabs(sol[k]-sol_e[k]);
    val += pow(diff[k]/sc[k],2);
  }

  // Fehler als Euklidische Norm der einzelnen Komponenten
  err[2] = sqrt(val/N2);

  controller = TDatabase::TimeDB->TIMESTEPLENGTH_CONTROLLER;
  // Berechnen der neuen Schrittweite entsprechend der Controller
  switch(controller)
  {
    case 0: h_new = pow((eps/err[2]),k_I)*steps_old[1]; break;
    // I controller
    case 1: h_new = pow((eps/err[2]),k_I)*pow((err[1]/err[2]),k_P)*steps_old[1]; break;
    // PI controller
    case 2: h_new = pow((eps/err[2]),k_E)*pow((err[1]/err[2]),k_R)*pow(steps_old[1],2.0)/steps_old[0]; break;
    // PC controller
    case 3: h_new = pow((eps/err[2]),k_I)*pow((err[1]/err[2]),k_P)*pow((pow(err[1],2.0)/(err[2]*err[0])),k_D)*steps_old[1]; break;
    // PID controller
  }

  TDatabase::TimeDB->TIMESTEPLENGTH = std::min(steps_old[1]*fac_max, std::max(steps_old[1]*fac_min,h_new));
  OutPut(TDatabase::TimeDB->TIMESTEPLENGTH << " ::: ");
  // Ueberpruefen der CFL-Bedingung:
  if(TDatabase::TimeDB->TIMESTEPLENGTH > 0.9/(N_x*sqrt(2)))
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.9/(N_x*sqrt(2));
  OutPut(TDatabase::TimeDB->TIMESTEPLENGTH << endl);
 
  // kein Unterschreiten der minimalen Schrittweite:
  if(TDatabase::TimeDB->TIMESTEPLENGTH < 1e-7)
    TDatabase::TimeDB->TIMESTEPLENGTH = 1e-7;

  // Ueberpruefen, ob Schrittweite angenommen wird.
  if (err[2] <= tol)
  {
    *acc = true;
    // Ueberschreiben der alten Schrittweite und des Fehlers
    err[0] = err[1];
    err[1] = err[2];
    steps_old[0] = steps_old[1];
    steps_old[1] = TDatabase::TimeDB->TIMESTEPLENGTH;
  }
  else
  {
    (*m)--;
    for(k = 0; k < N2; k++)
      // Wieder mit alter Naeherung weiterrechnen
      sol[k] = sol_old[k];
    steps_old[1] = TDatabase::TimeDB->TIMESTEPLENGTH;;
    (*disc)++;
  }
  OutPut("steps: " << steps_old[1] << endl);
}


/******************************************************************************/
// SetBoundaryConditions_FDM
// sets the boundary conditions
// only necessary if the boundary conditions are time-dependent
// otherwise, they are set at the initial time
/******************************************************************************/

#ifdef __2D__
void SetBoundaryConditions_FDM2D(BoundCondFunct2D *bound_cond,
BoundValueFunct2D *bound_val,
double *sol, int N_x, int N_y,
int *dof_conversion,
double *x_coord, double *y_coord)
{
  int i, j, N2;
  BoundCond cond;

  N2 = (N_x+1)*(N_y+1);

  // loop over all nodes of the FDM grid
  for ( i=0 ; i<N2 ; i++ )
  {
    // set Dirichlet boundary conditions
    // on bottom
    if (i<= N_x)
    {
      bound_cond(0, x_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(0, x_coord[i], sol[i]);
        continue;
      }
    }
    // inflow from the left x = x_min (left)
    if ( i%(N_x+1)==0 )
    {
      bound_cond(3, 1-y_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(3, 1-y_coord[i], sol[i]);
        continue;
      }
    }
    // on top
    if (i> (N_x+1)*N_y)
    {
      bound_cond(2, 1.0-x_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(2, 1.0-x_coord[i], sol[i]);
        continue;
      }
    }
    // on the right hand side
    if ( (i+1)%(N_x+1)==0 )
    {
      bound_cond(1, y_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(1, y_coord[i], sol[i]);
        continue;
      }
    }
  }
}


void WriteGnuplotFDM(const char *name, double *sol, double *x_coord,
double *y_coord, int N_x, int N_y)
{
  int i, j, index;
  FILE* out = fopen(name,"w");
  if (out==NULL)
  {
    OutPut("cannot write gnuplot fdm data" << endl);
    return;
  }
  OutPut("write gnuplot " << name << " " << N_x+1 << " " << N_y+1 << endl);

  index = 0;
  for (i=0;i<=N_y;i++)
  {
    for (j=0;j<=N_x;j++)
    {
      fprintf(out,"%f %f %f\n",x_coord[index], y_coord[index], sol[index]);
      index++;
    }
    fprintf(out,"\n");
  }
  fclose(out);
}


void WriteVTKFDM(const char *name, double *sol, double *x_coord,
double *y_coord, int N_x, int N_y)
{
  int i, j, k, N2, index;
  FILE* out = fopen(name,"w");
  if (out==NULL)
  {
    OutPut("cannot write vtk fdm data" << endl);
    return;
  }
  OutPut("write vtk " << name << " " << N_x+1 << " " << N_y+1 <<  endl);

  fprintf(out,"# vtk DataFile Version 4.2\n");
  fprintf(out,"file created by MooNMD\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET UNSTRUCTURED_GRID\n\n");
  N2 = (N_x+1)* (N_y+1);
  fprintf(out,"POINTS %d float\n",N2);

  index = 0;
  for (i=0;i<=N_x;i++)
  {
    for (j=0;j<=N_y;j++)
    {
      fprintf(out,"%f %f %f\n",x_coord[index], y_coord[index], 0.0);
      index++;
    }
  }

  fprintf(out,"\n");
  fprintf(out,"CELLS %d %d\n",N_x*N_y, 5*N_x*N_y);

  for (j=0;j<N_y;j++)
  {
    for (i=0;i<N_x;i++)
    {
      // index of left lower vertex
      index = j * (N_x+1) + i;
      fprintf(out,"4 %d %d %d %d\n",index, index+(N_x+1), index+(N_x+2),index+1);
    }
  }

  fprintf(out,"\n");
  fprintf(out,"CELL_TYPES %d\n",N_x*N_y);
  for (j=0;j<N_y;j++)
  {
    for (i=0;i<N_x;i++)
    {
      fprintf(out,"9 ");
    }
  }

  fprintf(out,"\n\n");

  fprintf(out,"POINT_DATA %d\n",N2);
  fprintf(out,"SCALARS solution float\n");
  fprintf(out,"LOOKUP_TABLE default\n");

  index = 0;
  for (i=0;i<=N_x;i++)
  {
    for (j=0;j<=N_y;j++)
    {
      fprintf(out,"%f\n",sol[index]);
      index++;
    }
  }

  fclose(out);
}
#endif
void InitializeConvectiveTermFDM(int dim, int* &offset_, int* &offset1_, int *N1_)
{
  int l;

  if (offset_==NULL)
  {
    offset_ = new int[2*dim];
    offset1_ = offset_+dim;
  }

  offset_[0] = 1;
  offset1_[0] = N1_[0];
  for(l=1;l<dim;l++)
  {
    offset_[l] = N1_[l-1] * offset_[l-1];
    offset1_[l] = N1_[l] * offset1_[l-1] ;
  }
}


void ClearConvectiveTermFDM(int* &offset)
{
  delete[] offset;
  offset = NULL;
}


/******************************************************************************/
// ConvectiveTermFDM
// computes the convective term in node i
// input:    dim        - dimension in which the problem is defined
//           coeff      - convection vector
//           sol_curr   - current solution
//           current_stage_fdm - array for the current stage
//
//
// output:   current_stage_fdm is updated
/******************************************************************************/

void ConvectiveTermFDM(int dim, int i,
double *coeff, double *sol_curr, double *current_stage_fdm,
double **coordinates, int *offset_, int *offset1_)
{
  int j, k, l, offset, offset1, off1;
  int offsetp2, offsetm2, offsetp1, offsetm1, offsetp3, offsetm3;
  double val[6], coord[6], value, valx, valy, valz;
  double d1, d2, d3, d4, uhx[3], omega[3], alpha[3], d[3], beta[3], av[5], c_e = 1e-6;
  double *coordinate;

  switch(TDatabase::ParamDB->DISCTYPE)
  {
    case GALERKIN:
      exit(1);
      break;
      // simple upwind scheme
    case UPWIND:
      // compute the term \partial_k A
      for (k=0;k<dim;k++)
      {
        // the offsets and the coordinates for the k-th coordinate direction
        offset = offset_[k];
        offset1 = offset1_[k];
        coordinate = coordinates[k];
        // first case
        if (coeff[k] > 0)
        {
          off1 = i-offset;
          // not on upwind boundary
          if (!( i%offset1 < offset))
            current_stage_fdm[i] -= coeff[k]*(sol_curr[i]-sol_curr[off1])/(coordinate[i]-coordinate[off1]);
        }
        else
        {
          off1 = i+offset;
          // not on upwind boundary
          if (!((i+offset )%offset1 < offset ))
            current_stage_fdm[i] -= coeff[k]*(sol_curr[off1]-sol_curr[i])/(coordinate[off1]-coordinate[i]);
        }
      }
      break;
    case ENO_3:
      for (k=0;k<dim;k++)
      {
        offset = offset_[k];
        offset1 = offset1_[k];
        coordinate = coordinates[k];
        offsetp1 = i+offset;
        offsetm1 = i-offset;
        offsetp2 = i+2*offset;
        offsetm2 = i-2*offset;
        offsetp3 = i+3*offset;
        offsetm3 = i-3*offset;
        // compute the term \partial_k A
        if (coeff[k] >  0)
        {
          // prepare vectors for divided differences
          // point on left, front, lower ...  boundary
          if ( i%offset1 < offset )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[offsetp1];
            val[5] = sol_curr[offsetp2];
            coord[0] = 2 * coordinate[i]-coordinate[offsetp3];
            coord[1] = 2 * coordinate[i]-coordinate[offsetp2];
            coord[2] = 2 * coordinate[i]-coordinate[offsetp1];
            coord[3] = coordinate[i];
            coord[4] = coordinate[offsetp1];
            coord[5] = coordinate[offsetp2];
          }
          else
          {
            // point next to left, front, lower ...  boundary
            if ( (offsetm1)%offset1 < offset )
            {
              val[0] = val[1] = val[2] = sol_curr[offsetm1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[offsetp1];
              val[5] = sol_curr[offsetp2];
              coord[0] = 2 * coordinate[offsetm1]-coordinate[offsetp1];
              coord[1] = 2 * coordinate[offsetm1]-coordinate[i];
              coord[2] = coordinate[offsetm1];
              coord[3] = coordinate[i];
              coord[4] = coordinate[offsetp1];
              coord[5] = coordinate[offsetp2];
            }
            else
            {
              // next layer on left, front, lower, ...  boundary
              if ( (offsetm2)%offset1 < offset )
              {
                val[0] = val[1] = sol_curr[offsetm2];
                val[2] = sol_curr[offsetm1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[offsetp1];
                val[5] = sol_curr[offsetp2];
                coord[0] = 2 * coordinate[offsetm2]-coordinate[offsetm1];
                coord[1] = coordinate[offsetm2];
                coord[2] = coordinate[offsetm1];
                coord[3] = coordinate[i];
                coord[4] = coordinate[offsetp1];
                coord[5] = coordinate[offsetp2];
              }
              else
              {
                // point on right, back, upper, ...  boundary
                if ( (offsetp1)%offset1 < offset )
                {
                  val[0] = sol_curr[offsetm3];
                  val[1] = sol_curr[offsetm2];
                  val[2] = sol_curr[offsetm1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = coordinate[offsetm3];
                  coord[1] = coordinate[offsetm2];
                  coord[2] = coordinate[offsetm1];
                  coord[3] = coordinate[i];
                  coord[4] = 2 *  coordinate[i] - coordinate[offsetm1];
                  coord[5] = 2 *  coordinate[i] - coordinate[offsetm2];
                }
                else
                {
                  // point next to right, back, upper, ...  boundary
                  if ( (offsetp2)% offset1 < offset )
                  {
                    val[0] = sol_curr[offsetm3];
                    val[1] = sol_curr[offsetm2];
                    val[2] = sol_curr[offsetm1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[offsetp1];
                    coord[0] = coordinate[offsetm3];
                    coord[1] = coordinate[offsetm2];
                    coord[2] = coordinate[offsetm1];
                    coord[3] = coordinate[i];
                    coord[4] = coordinate[offsetp1];
                    coord[5] = 2 *  coordinate[offsetp1] - coordinate[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[offsetm3];
                    val[1] = sol_curr[offsetm2];
                    val[2] = sol_curr[offsetm1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[offsetp1];
                    val[5] = sol_curr[offsetp2];
                    coord[0] = coordinate[offsetm3];
                    coord[1] = coordinate[offsetm2];
                    coord[2] = coordinate[offsetm1];
                    coord[3] = coordinate[i];
                    coord[4] = coordinate[offsetp1];
                    coord[5] = coordinate[offsetp2];
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+2,coord+2);
          d2 = DividedDifferences(2,val+1,coord+1);
          if (fabs(d1) < fabs(d2))
          {
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value =  (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d3 * (coord[3]-coord[1])*(coord[3]-coord[2]);
            }
            else
            {
              valx = coord[3]-coord[2];
              value = (val[3]-val[2])/valx;
              value += d1 * valx;
              value += d4 * valx*(coord[3]-coord[4]);
            }
          }
          else
          {
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = (val[1]-val[0])/(coord[1]-coord[0]);
              value += DividedDifferences(2,val,coord) * (coord[3]-coord[0] + coord[3]-coord[1]);
              valx = coord[3]-coord[0];
              valy = coord[3]-coord[1];
              valz = coord[3]-coord[2];
              value += d3 * (valy * valz + valx * valz + valx * valy);
            }
            else
            {
              value = (val[2]-val[1])/(coord[2]-coord[1]);
              value += d2 * (coord[3]-coord[1] + coord[3]-coord[2]);
              value += d4 * (coord[3]-coord[1])*(coord[3]-coord[2]);
            }
          }
        }
        else                                      // case coeff[k] <= 0
        {
          // prepare vectors for divided differences
          // point on left, front, lower ... boundary
          if ( i%offset1 < offset)
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[offsetp1];
            val[4] = sol_curr[offsetp2];
            val[5] = sol_curr[offsetp3];
            coord[0] = 2 * coordinate[i]-coordinate[offsetp2];
            coord[1] = 2 * coordinate[i]-coordinate[offsetp1];
            coord[2] = coordinate[i];
            coord[3] = coordinate[offsetp1];
            coord[4] = coordinate[offsetp2];
            coord[5] = coordinate[offsetp3];
          }
          else
          {
            // point next to left boundary
            if ( (offsetm1)%offset1 < offset )
            {
              val[0] = val[1] = sol_curr[offsetm1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[offsetp1];
              val[4] = sol_curr[offsetp2];
              val[5] = sol_curr[offsetp3];
              coord[0] = 2 * coordinate[offsetm1]-coordinate[i];
              coord[1] = coordinate[offsetm1];
              coord[2] = coordinate[i];
              coord[3] = coordinate[offsetp1];
              coord[4] = coordinate[offsetp2];
              coord[5] = coordinate[offsetp3];
            }
            else
            {
              // point on right boundary
              if ( (offsetp1)%offset1 < offset )
              {
                val[0] = sol_curr[offsetm2];
                val[1] = sol_curr[offsetm1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = coordinate[offsetm2];
                coord[1] = coordinate[offsetm1];
                coord[2] = coordinate[i];
                coord[3] = 2 *  coordinate[i] - coordinate[offsetm1];
                coord[4] = 2 *  coordinate[i] - coordinate[offsetm2];
                coord[5] = 2 *  coordinate[i] - coordinate[offsetm3];
              }
              else
              {
                // point next to right boundary
                if ( (offsetp2)%offset1 < offset )
                {
                  val[0] = sol_curr[offsetm2];
                  val[1] = sol_curr[offsetm1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[offsetp1];
                  coord[0] = coordinate[offsetm2];
                  coord[1] = coordinate[offsetm1];
                  coord[2] = coordinate[i];
                  coord[3] = coordinate[offsetp1];
                  coord[4] = 2 *  coordinate[offsetp1] - coordinate[i];
                  coord[5] = 2 *  coordinate[offsetp1] - coordinate[offsetm1];
                }
                else
                {
                  // point in next layer to right boundary
                  if ( (offsetp3)%offset1 < offset )
                  {
                    val[0] = sol_curr[offsetm2];
                    val[1] = sol_curr[offsetm1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[offsetp1];
                    val[4] = val[5] = sol_curr[offsetp2];
                    coord[0] = coordinate[offsetm2];
                    coord[1] = coordinate[offsetm1];
                    coord[2] = coordinate[i];
                    coord[3] = coordinate[offsetp1];
                    coord[4] = coordinate[offsetp2];
                    coord[5] = 2 *  coordinate[offsetp2] - coordinate[offsetp1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[offsetm2];
                    val[1] = sol_curr[offsetm1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[offsetp1];
                    val[4] = sol_curr[offsetp2];
                    val[5] = sol_curr[offsetp3];
                    coord[0] = coordinate[offsetm2];
                    coord[1] = coordinate[offsetm1];
                    coord[2] = coordinate[i];
                    coord[3] = coordinate[offsetp1];
                    coord[4] = coordinate[offsetp2];
                    coord[5] = coordinate[offsetp3];
                  }
                }
              }
            }
          }                                       // end  prepare vectors for divided differences
          // compute approximation of convective term
          d1 = DividedDifferences(2,val+1,coord+1);
          d2 = DividedDifferences(2,val+2,coord+2);
          if (fabs(d1) < fabs(d2))
          {
            DividedDifferences_3_2(val,coord,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              value = DividedDifferences(1,val,coord);
              value += DividedDifferences(2,val,coord) * (coord[2]-coord[0] + coord[2]-coord[1]);
              value += d3 * (coord[2]-coord[0])*(coord[2]-coord[1]);
            }
            else
            {
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * (coord[2]-coord[1]);
              value += d4 * (coord[2]-coord[1])*(coord[2]-coord[3]);
            }
          }
          else
          {
            DividedDifferences_3_2(val+1,coord+1,d);
            d3 = d[0];
            d4 = d[1];
            if   (fabs(d3) < fabs(d4))
            {
              valx = coord[2]-coord[1];
              value = DividedDifferences(1,val+1,coord+1);
              value += d1 * valx;
              value += d3 * valx*(coord[2]-coord[3]);
            }
            else
            {
              valx = coord[2]-coord[3];
              value = DividedDifferences(1,val+2,coord+2);
              value += d2 * valx;
              value += d4 * valx*(coord[2]-coord[4]);
            }
          }
        }                                         // end coeff[k] <=0
        current_stage_fdm[i] -= coeff[k] * value;
      }                                           // end k
      break;
    case WENO_5:
      d[0] = 0.3;
      d[1] = 0.6;
      d[2] = 0.1;
      for (k=0;k<dim;k++)
      {
        offset = offset_[k];
        offset1 = offset1_[k];
        coordinate = coordinates[k];
        offsetp1 = i+offset;
        offsetm1 = i-offset;
        offsetp2 = i+2*offset;
        offsetm2 = i-2*offset;
        offsetp3 = i+3*offset;
        offsetm3 = i-3*offset;
        // compute the term \partial_k A
        if (coeff[k] >  0)
        {
          // prepare vectors for divided differences
          // point on left, front, lower ...  boundary
          if ( i%offset1 < offset )
          {
            val[0] = val [1] = val[2] = val[3] = sol_curr[i];
            val[4] = sol_curr[offsetp1];
            val[5] = sol_curr[offsetp2];
            coord[0] = 2 * coordinate[i]-coordinate[offsetp3];
            coord[1] = 2 * coordinate[i]-coordinate[offsetp2];
            coord[2] = 2 * coordinate[i]-coordinate[offsetp1];
            coord[3] = coordinate[i];
            coord[4] = coordinate[offsetp1];
            coord[5] = coordinate[offsetp2];
          }
          else
          {
            // point next to left, front, lower ...  boundary
            if ( (offsetm1)%offset1 < offset )
            {
              val[0] = val[1] = val[2] = sol_curr[offsetm1];
              val[3] = sol_curr[i];
              val[4] = sol_curr[offsetp1];
              val[5] = sol_curr[offsetp2];
              coord[0] = 2 * coordinate[offsetm1]-coordinate[offsetp1];
              coord[1] = 2 * coordinate[offsetm1]-coordinate[i];
              coord[2] = coordinate[offsetm1];
              coord[3] = coordinate[i];
              coord[4] = coordinate[offsetp1];
              coord[5] = coordinate[offsetp2];
            }
            else
            {
              // next layer on left, front, lower, ...  boundary
              if ( (offsetm2)%offset1 < offset )
              {
                val[0] = val[1] = sol_curr[offsetm2];
                val[2] = sol_curr[offsetm1];
                val[3] = sol_curr[i];
                val[4] = sol_curr[offsetp1];
                val[5] = sol_curr[offsetp2];
                coord[0] = 2 * coordinate[offsetm2]-coordinate[offsetm1];
                coord[1] = coordinate[offsetm2];
                coord[2] = coordinate[offsetm1];
                coord[3] = coordinate[i];
                coord[4] = coordinate[offsetp1];
                coord[5] = coordinate[offsetp2];
              }
              else
              {
                // point on right, back, upper, ...  boundary
                if ( (offsetp1)%offset1 < offset )
                {
                  val[0] = sol_curr[offsetm3];
                  val[1] = sol_curr[offsetm2];
                  val[2] = sol_curr[offsetm1];
                  val[3] = val[4] = val[5] = sol_curr[i];
                  coord[0] = coordinate[offsetm3];
                  coord[1] = coordinate[offsetm2];
                  coord[2] = coordinate[offsetm1];
                  coord[3] = coordinate[i];
                  coord[4] = 2 *  coordinate[i] - coordinate[offsetm1];
                  coord[5] = 2 *  coordinate[i] - coordinate[offsetm2];
                }
                else
                {
                  // point next to right, back, upper, ...  boundary
                  if ( (offsetp2)% offset1 < offset )
                  {
                    val[0] = sol_curr[offsetm3];
                    val[1] = sol_curr[offsetm2];
                    val[2] = sol_curr[offsetm1];
                    val[3] = sol_curr[i];
                    val[4] = val[5] = sol_curr[offsetp1];
                    coord[0] = coordinate[offsetm3];
                    coord[1] = coordinate[offsetm2];
                    coord[2] = coordinate[offsetm1];
                    coord[3] = coordinate[i];
                    coord[4] = coordinate[offsetp1];
                    coord[5] = 2 *  coordinate[offsetp1] - coordinate[i];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[offsetm3];
                    val[1] = sol_curr[offsetm2];
                    val[2] = sol_curr[offsetm1];
                    val[3] = sol_curr[i];
                    val[4] = sol_curr[offsetp1];
                    val[5] = sol_curr[offsetp2];
                    coord[0] = coordinate[offsetm3];
                    coord[1] = coordinate[offsetm2];
                    coord[2] = coordinate[offsetm1];
                    coord[3] = coordinate[i];
                    coord[4] = coordinate[offsetp1];
                    coord[5] = coordinate[offsetp2];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/3.0-val[3]/2.0 + val[4] - val[5]/6.0)/(coord[3]-coord[2]);
          uhx[1] = (val[1]/6.0-val[2]+val[3]/2.0 + val[4]/3.0)/(coord[2]-coord[1]);
          uhx[2] = (-val[0]/3.0+1.5*val[1]-3*val[2] + 11.0*val[3]/6.0)/(coord[1]-coord[0]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[3] + av[4];
          valy = 3*av[2] - 4*av[3] + av[4];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[1] - 2*av[2] + av[3];
          valy = av[1] - av[3];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[0] - 2*av[1] +av[2];
          valy = av[0] - 4*av[1] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }
          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }
        else                                      // case coeff[k] <= 0
        {
          // prepare vectors for divided differences
          // point on left, front, lower ... boundary
          if ( i%offset1 < offset)
          {
            val[0] = val [1] = val[2] = sol_curr[i];
            val[3] = sol_curr[offsetp1];
            val[4] = sol_curr[offsetp2];
            val[5] = sol_curr[offsetp3];
            coord[0] = 2 * coordinate[i]-coordinate[offsetp2];
            coord[1] = 2 * coordinate[i]-coordinate[offsetp1];
            coord[2] = coordinate[i];
            coord[3] = coordinate[offsetp1];
            coord[4] = coordinate[offsetp2];
            coord[5] = coordinate[offsetp3];
          }
          else
          {
            // point next to left boundary
            if ( (offsetm1)%offset1 < offset )
            {
              val[0] = val[1] = sol_curr[offsetm1];
              val[2] = sol_curr[i];
              val[3] = sol_curr[offsetp1];
              val[4] = sol_curr[offsetp2];
              val[5] = sol_curr[offsetp3];
              coord[0] = 2 * coordinate[offsetm1]-coordinate[i];
              coord[1] = coordinate[offsetm1];
              coord[2] = coordinate[i];
              coord[3] = coordinate[offsetp1];
              coord[4] = coordinate[offsetp2];
              coord[5] = coordinate[offsetp3];
            }
            else
            {
              // point on right boundary
              if ( (offsetp1)%offset1 < offset )
              {
                val[0] = sol_curr[offsetm2];
                val[1] = sol_curr[offsetm1];
                val[2] = val[3] = val[4] = val[5] = sol_curr[i];
                coord[0] = coordinate[offsetm2];
                coord[1] = coordinate[offsetm1];
                coord[2] = coordinate[i];
                coord[3] = 2 *  coordinate[i] - coordinate[offsetm1];
                coord[4] = 2 *  coordinate[i] - coordinate[offsetm2];
                coord[5] = 2 *  coordinate[i] - coordinate[offsetm3];
              }
              else
              {
                // point next to right boundary
                if ( (offsetp2)%offset1 < offset )
                {
                  val[0] = sol_curr[offsetm2];
                  val[1] = sol_curr[offsetm1];
                  val[2] = sol_curr[i];
                  val[3] = val[4] = val[5] = sol_curr[offsetp1];
                  coord[0] = coordinate[offsetm2];
                  coord[1] = coordinate[offsetm1];
                  coord[2] = coordinate[i];
                  coord[3] = coordinate[offsetp1];
                  coord[4] = 2 *  coordinate[offsetp1] - coordinate[i];
                  coord[5] = 2 *  coordinate[offsetp1] - coordinate[offsetm1];
                }
                else
                {
                  // point in next layer to right boundary
                  if ( (offsetp3)%offset1 < offset )
                  {
                    val[0] = sol_curr[offsetm2];
                    val[1] = sol_curr[offsetm1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[offsetp1];
                    val[4] = val[5] = sol_curr[offsetp2];
                    coord[0] = coordinate[offsetm2];
                    coord[1] = coordinate[offsetm1];
                    coord[2] = coordinate[i];
                    coord[3] = coordinate[offsetp1];
                    coord[4] = coordinate[offsetp2];
                    coord[5] = 2 *  coordinate[offsetp2] - coordinate[offsetp1];
                  }
                  // inner point
                  else
                  {
                    val[0] = sol_curr[offsetm2];
                    val[1] = sol_curr[offsetm1];
                    val[2] = sol_curr[i];
                    val[3] = sol_curr[offsetp1];
                    val[4] = sol_curr[offsetp2];
                    val[5] = sol_curr[offsetp3];
                    coord[0] = coordinate[offsetm2];
                    coord[1] = coordinate[offsetm1];
                    coord[2] = coordinate[i];
                    coord[3] = coordinate[offsetp1];
                    coord[4] = coordinate[offsetp2];
                    coord[5] = coordinate[offsetp3];
                  }
                }
              }
            }
          }
          // compute values for WENO scheme
          uhx[0] = (-val[2]/2.0+val[1] - val[0]/6.0 - val[3]/3.0)/(-coord[2]+coord[1]);
          uhx[1] = (-val[3] +val[2]/2.0+val[1]/3.0 + val[4]/6.0)/(-coord[3]+coord[2]);
          uhx[2] = (-val[5]/3.0+1.5*val[4]-3*val[3] + 11.0*val[2]/6.0)/(-coord[4]+coord[3]);

          // compute smooth indicators
          for (j=0;j<5;j++)
            av[j] =  (val[j+1]-val[j])/(coord[j+1]-coord[j]);
          valx = av[2] - 2* av[1] + av[0];
          valy = 3*av[2] - 4*av[1] + av[0];
          beta[0] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[3] - 2*av[2] + av[1];
          valy = av[3] - av[1];
          beta[1] = 13.0*valx*valx/12.0 + valy*valy/4.0;
          valx = av[4] - 2*av[3] +av[2];
          valy = av[4] - 4*av[3] +3*av[2];
          beta[2] = 13.0*valx*valx/12.0 + valy*valy/4.0;

          // compute alpha
          if (TDatabase::ParamDB->WENO_TYPE == 0)
          {
            alpha[0] = d[0]/((beta[0] + c_e)*(beta[0] + c_e));
            alpha[1] = d[1]/((beta[1] + c_e)*(beta[1] + c_e));
            alpha[2] = d[2]/((beta[2] + c_e)*(beta[2] + c_e));
          }
          else
          {
            alpha[0] = d[0]/(beta[0] + c_e);
            alpha[1] = d[1]/(beta[1] + c_e);
            alpha[2] = d[2]/(beta[2] + c_e);
          }
          // compute weights
          valx = alpha[0] + alpha[1] + alpha[2];
          omega[0] = alpha[0] / valx;
          omega[1] = alpha[1] / valx;
          omega[2] = alpha[2] / valx;
        }                                         // end coeff[k] <=0
        current_stage_fdm[i] -= coeff[k] * (omega[0]*uhx[0] + omega[1]*uhx[1] + omega[2]*uhx[2]);
      }                                           // end k
      break;
    default:
      OutPut("Discretization " << TDatabase::ParamDB->DISCTYPE << " not defined for convective term !!!" << endl);
      exit(1);
  }
}


#ifdef __3D__

void SetBoundaryConditions_FDM3D(BoundCondFunct3D *bound_cond,
BoundValueFunct3D *bound_val,
double *sol, int N_x, int N_y, int N_z,
int *dof_conversion,
double *x_coord, double *y_coord, double *z_coord)

{
  int i, j, N3, N_x1, N_y1, N_z1;
  BoundCond cond;

  // set boundary conditions at current time
  N_x1 = N_x+1;
  N_y1 = N_y+1;
  N_z1 = N_z+1;
  N3 = (N_x1)*(N_y1)*(N_z1);
  // loop over all nodes of the FDM grid
  for ( i=0 ; i<N3 ; i++ )
  {
    // set Dirichlet boundary conditions
    // inflow from the left x = x_min (left)
    if ( i%(N_x1)==0 )
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(x_coord[i], y_coord[i], z_coord[i], sol[i]);
        continue;
      }
    }

    // on the right hand side
    if ( (i+1)%(N_x1)==0 )
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(x_coord[i], y_coord[i], z_coord[i], sol[i]);
        continue;
      }
    }

    // on front
    if ((i%((N_x1)*(N_y1)))<N_x1 )
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(x_coord[i], y_coord[i], z_coord[i], sol[i]);
        continue;
      }
    }

    // on back
    if (((i+N_x1)%((N_x1)*(N_y1)))<N_x1)
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(x_coord[i], y_coord[i], z_coord[i], sol[i]);
        continue;
      }
    }

    // on bottom
    if (i< N_x1 * N_y1)
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(x_coord[i], y_coord[i], z_coord[i], sol[i]);
        continue;
      }
    }

    // on top
    if (i>= ((N_x1)*(N_y1)*(N_z)))
    {
      bound_cond(x_coord[i], y_coord[i], z_coord[i], cond);
      if (cond == DIRICHLET)
      {
        //j =  dof_conversion[i];
        bound_val(x_coord[i], y_coord[i], z_coord[i], sol[i]);
        continue;
      }
    }
  }
}


void WriteVTKFDM(const char *name, double *sol, double *x_coord,
double *y_coord, double *z_coord, int N_x, int N_y, int N_z)
{
  int i, j, k, N2, N3, index;
  FILE* out = fopen(name,"w");
  if (out==NULL)
  {
    OutPut("cannot write vtk fdm data" << endl);
    return;
  }
  OutPut("write vtk " << name << " " << N_x+1 << " " << N_y+1 << " " << N_z+1 <<  endl);

  fprintf(out,"# vtk DataFile Version 4.2\n");
  fprintf(out,"file created by MooNMD\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET UNSTRUCTURED_GRID\n\n");
  N3 = (N_x+1)* (N_y+1)*(N_z+1);
  fprintf(out,"POINTS %d float\n",N3);

  index = 0;
  for (i=0;i<=N_x;i++)
  {
    for (j=0;j<=N_y;j++)
    {
      for (k=0;k<=N_z;k++)
      {
        fprintf(out,"%f %f %f\n",x_coord[index], y_coord[index], z_coord[index]);
        index++;
      }
    }
  }

  fprintf(out,"\n");
  fprintf(out,"CELLS %d %d\n",N_x*N_y*N_z, 9*N_x*N_y*N_z);
  for (k=0;k<N_z;k++)
  {
    for (j=0;j<N_y;j++)
    {
      N2 = (N_x+1)*(N_y+1);
      for (i=0;i<N_x;i++)
      {
        // index of left lower vertex
        index = k*N2 + j * (N_x+1) + i;
        fprintf(out,"8 %d %d %d %d %d %d %d %d\n",index+N2, index+(N_x+1)+N2, index+(N_x+2)+N2,index+1+N2,
          index, index+(N_x+1), index+(N_x+2),index+1);
      }
    }
  }

  fprintf(out,"\n");
  fprintf(out,"CELL_TYPES %d\n",N_x*N_y*N_z);
  for (k=0;k<N_z;k++)
  {
    for (j=0;j<N_y;j++)
    {
      for (i=0;i<N_x;i++)
      {
        fprintf(out,"12 ");
      }
    }
  }
  fprintf(out,"\n\n");

  N3 = (N_x+1)* (N_y+1)*(N_z+1);
  fprintf(out,"POINT_DATA %d\n",N3);
  fprintf(out,"SCALARS solution float\n");
  fprintf(out,"LOOKUP_TABLE default\n");

  index = 0;
  for (i=0;i<=N_x;i++)
  {
    for (j=0;j<=N_y;j++)
    {
      for (k=0;k<=N_z;k++)
      {
        fprintf(out,"%f\n",sol[index]);
        index++;
      }
    }
  }

  fclose(out);
}
#endif
