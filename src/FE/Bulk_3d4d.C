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
   
/****************************************************************************************
 *                                                                                      *
 *                         Bulk_3d4d.C                                                  *
 *                        -------------                                                 *
 *  routines for simulation of precipitation process in 3d cavity                       *
 *                                                                                      *
 ***************************************************************************************/

#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <FEFunction3D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <Solver.h>
#include <FEM_TVD_FCT.h>
#include <Bulk_3d4d.h>
#include <Bulk.h>
#include <TimeDiscRout.h>
#include <RKV_FDM.h>
#include <Database.h>

#include <string.h>
#include <stdlib.h>

/****************************************************************************************
 *                                                                                      *
 *  The declared functions are sorted in the following way :                            *
 *                                                                                      *
 *    1. general routines (e.g. grid generation,...)                                    *
 *    2. FDM routines                                                                   *
 *    3. FEM routines                                                                   *
 *    4. analysing routines                                                             *
 *                                                                                      *
 ***************************************************************************************/


/****************************************************************************************
 *                                                                                      *
 *                         Part I : general routines                                    *
 *                        ---------------------------                                   *
 *                                                                                      *
 ***************************************************************************************/


/****************************************************************************************
 *                                                                                      *
 * generates non-equidistant 4d grid for fem                                            *
 * input:  x_min, x_max, y_min, y_max, z_min, z_max                                     *
 *         are the intervals in x-, y- and z-direction = 3d space                       *
 *         a_min and a_max are the interval limits for the attribute of the reactant    *
 *         N_x, N_y, N_z : number of intervals in x-, y- and z-direction (mesh cells)   *
 *         in analogy : N_a is the number of mesh cells in the 4th dimension            *
 * output: *x_coord,  *y_coord,  *z_coord,  *a_coord, *a_layers_coord                   *
 *                                                                                      *
 * grid type 0: equidistant grid                                                        *
 * grid type 1: non-equidistant grid based on the cosinus function                      *
 * grid type 2: non-equidistant grid based on the hyperbolic tangent function           *
 *                                                                                      *
 ***************************************************************************************/

void grid_generator_4d(TCollection *coll,
		       int &N_x,
                       int &N_y,
                       int &N_z,
		       double &x_min, double &x_max, double &y_min, double &y_max, 
		       double &z_min, double &z_max, 
                       double a_min, double a_max, int N_a,
                       double* &x_coord, double* &y_coord, double* &z_coord, double* &a_coord,
                       double *a_layers_coord)
{
  // check consistency
  if ((a_min >= a_max )|| (N_a <= 0) )
  {
    OutPut(" error : wrong input in grid_generator_4d " << endl);
    exit(-1);
  }

  // integer variables for the loops
  int i, j, k, N_Cells, N_Vertices, found, Nodes, N2, N3, max_dir = 1025;
  double eps = 1e-5, mini, x, y, z, gamma, a1, a2, a3;
  TBaseCell *cell;

  // arrays for the grid -> this allows to implement different distributions of grid points in every direction
  double *x_help, *y_help, *z_help, *a_help;

  x_help = new double[max_dir];  memset(x_help,0,max_dir*SizeOfDouble);
  y_help = new double[max_dir];  memset(y_help,0,max_dir*SizeOfDouble);
  z_help = new double[max_dir];  memset(z_help,0,max_dir*SizeOfDouble);
  a_help = new double[N_a+1];  memset(a_help,0,(N_a+1)*SizeOfDouble);
// double  input[88]= 
//
//{1.250400e-03 ,2.110869e-03 ,2.564006e-03 ,3.166752e-03 ,3.948493e-03 ,4.948318e-03 ,6.217678e-03 ,7.823151e-03 ,9.849857e-03 ,1.127270e-02 ,1.240580e-02 ,1.336284e-02 ,1.419948e-02 ,1.494769e-02 ,1.562767e-02 ,1.683382e-02 ,1.788820e-02 ,1.883107e-02 ,1.968795e-02 ,2.120783e-02 ,2.253643e-02 ,2.372451e-02 ,2.480420e-02 ,2.671929e-02 ,2.839332e-02 ,2.989028e-02 ,3.125067e-02 ,3.366362e-02 ,3.577283e-02 ,3.765893e-02 ,3.937295e-02 ,4.241314e-02 ,4.507062e-02 ,4.744700e-02 ,4.960655e-02 ,5.159286e-02 ,5.343698e-02 ,5.516191e-02 ,5.678522e-02 ,5.832068e-02 ,5.977928e-02 ,6.116998e-02 ,6.250017e-02 ,6.732623e-02 ,7.154477e-02 ,7.531706e-02 ,7.874517e-02 ,8.482564e-02 ,9.014069e-02 ,9.921263e-02 ,1.135701e-01 ,1.250000e-01 ,1.346522e-01 ,1.430893e-01 ,1.506339e-01 ,1.574902e-01 ,1.696511e-01 ,1.802812e-01 ,1.897868e-01 ,1.984251e-01 ,2.063705e-01 ,2.137470e-01 ,2.206468e-01 ,2.271401e-01 ,2.391164e-01 ,2.500000e-01 ,2.693043e-01 ,2.861786e-01 ,3.012678e-01 ,3.149803e-01 ,3.393022e-01 ,3.605624e-01 ,3.795736e-01 ,3.968503e-01 ,4.274940e-01 ,4.542801e-01 ,4.782328e-01 ,5.000000e-01 ,5.386087e-01 ,5.723571e-01 ,6.025356e-01 ,6.299605e-01 ,6.786044e-01 ,7.211248e-01 ,7.591472e-01 ,7.937005e-01 ,9.085603e-01 ,1.000000e+00 };
 double input[88] = {0.001250000000000 ,
0.002110731579413 ,
0.002563914004088 ,
0.003166693847592 ,
0.003948458933296 ,
0.004948299443979 ,
0.006217671606001 ,
0.007823152131089 ,
0.009849848192817 ,
0.011272700427129 ,
0.012405802335839 ,
0.013362834027725 ,
0.014199484324047 ,
0.014947689230070 ,
0.015627666201516 ,
0.016833818754898 ,
0.017888194844297 ,
0.018831072546427 ,
0.019687946143196 ,
0.021207835206421 ,
0.022536431448077 ,
0.023724508066731 ,
0.024804199643075 ,
0.026719286105474 ,
0.028393316830700 ,
0.029890278429282 ,
0.031250666632101 ,
0.033663616515945 ,
0.035772828811153 ,
0.037658931922342 ,
0.039372952752786 ,
0.042413137150944 ,
0.045070619540390 ,
0.047446991851855 ,
0.049606547407142 ,
0.051592857816434 ,
0.053436976294678 ,
0.055161906663233 ,
0.056785220390605 ,
0.058320678197501 ,
0.059779281606086 ,
0.061169980847062 ,
0.062500166625532 ,
0.067326227648337 ,
0.071544767303250 ,
0.075317060474946 ,
0.078745170560444 ,
0.084825640943857 ,
0.090140678210433 ,
0.099212631825087 ,
0.113570087453620 ,
0.125000041585273 ,
0.134652203946564 ,
0.143089312023620 ,
0.150633920105047 ,
0.157490157382677 ,
0.169651123546876 ,
0.180281216202297 ,
0.189786828693602 ,
0.198425147902268 ,
0.206370468208258 ,
0.213747007445228 ,
0.220646784186759 ,
0.227140086575032 ,
0.239116409077369 ,
0.250000010253906 ,
0.269304345055468 ,
0.286178568401423 ,
0.301267789998808 ,
0.314980268830742 ,
0.339302207508500 ,
0.360562397349922 ,
0.379573625756073 ,
0.396850266867541 ,
0.427493989953308 ,
0.454280151067003 ,
0.478232798228372 ,
0.500000002278646 ,
0.538608674401517 ,
0.572357122891389 ,
0.602535567444792 ,
0.629960526177828 ,
0.678604405120686 ,
0.721124785936175 ,
0.759147243604363 ,
0.793700526500832 ,
0.908560296613240 ,
1.000000000000000 };
  // compute coordinates of the 3d domain
  // initialize
  for (i=0;i<max_dir;i++)
      x_help[i] = -47111174;
  for (i=0;i<max_dir;i++)
      y_help[i] = -47111174;
  for (i=0;i<max_dir;i++)
      z_help[i] = -47111174;

  // number of mesh cells
  N_Cells = coll->GetN_Cells();
  for (i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    // number of vertices
    N_Vertices = cell->GetN_Vertices();      
    for ( j=0; j<N_Vertices; j++)
    {
      // get coordinates of the vertex
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      z = cell->GetVertex(j)->GetZ();
      // check if they are already included in the arrays
      found = 0;
      for (k=0;k<max_dir;k++)
      {
	  // coordinate already included
	  if (fabs(x-x_help[k]) < eps)
	  {
	      found = 1;
	      break;
	  }
	  // new entry
	  if (fabs(x_help[k]+47111174) < eps)
	  {
	      x_help[k] = x;
	      found = 1;
	      break;
	  }
      }
      if (!found)
      {
	  OutPut("Error in grid_generation_4d, x coordinate !!!" << endl);
	  exit(4711);
      }
      found = 0;
      for (k=0;k<max_dir;k++)
      {
	  // coordinate already included
	  if (fabs(y-y_help[k]) < eps)
	  {
	      found = 1;
	      break;
	  }
	  // new entry
	  if (fabs(y_help[k]+47111174) < eps)
	  {
	      y_help[k] = y;
	      found = 1;
	      break;
	  }
      }
      if (!found)
      {
	  OutPut("Error in grid_generation_4d, y coordinate !!!" << endl);
	  exit(4711);
      }
      found = 0;
      for (k=0;k<max_dir;k++)
      {
	  // coordinate already included
	  if (fabs(z-z_help[k]) < eps)
	  {
	      found = 1;
	      break;
	  }
	  // new entry
	  if (fabs(z_help[k]+47111174) < eps)
	  {
	      z_help[k] = z;
	      found = 1;
	      break;
	  }
      }
      if (!found)
      {
	  OutPut("Error in grid_generation_4d, z coordinate !!!" << endl);
	  exit(4711);
      }
    }
  }

  // compute actual length of arrays
  for (i=0;i<max_dir;i++)
  {
      if (fabs(x_help[i]+47111174) < eps)
	  break;
  }
  N_x = i-1;
  for (i=0;i<max_dir;i++)
  {
      if (fabs(y_help[i]+47111174) < eps)
	  break;
  }
  N_y = i-1;
  for (i=0;i<max_dir;i++)
  {
      if (fabs(z_help[i]+47111174) < eps)
	  break;
  }
  N_z = i-1;

  // bring arrays in correct order
  for (i=0;i<N_x+1;i++)
  {
      mini = 47111147.0;
      // find smallest entry
      for (j=i;j<N_x+1;j++)
      {
	  if (x_help[j] < mini)
	  {
	      mini = x_help[j];
	      found = j;
	  }
      }
      // change entries
      mini = x_help[i];
      x_help[i] = x_help[found];
      x_help[found] = mini;
  }
  x_min = x_help[0];
  x_max = x_help[N_x];

  for (i=0;i<N_y+1;i++)
  {
      mini = 47111147.0;
      // find smallest entry
      for (j=i;j<N_y+1;j++)
      {
	  if (y_help[j] < mini)
	  {
	      mini = y_help[j];
	      found = j;
	  }
      }
      // change entries
      mini = y_help[i];
      y_help[i] = y_help[found];
      y_help[found] = mini;
  }
  y_min = y_help[0];
  y_max = y_help[N_y];

  for (i=0;i<N_z+1;i++)
  {
      mini = 47111147.0;
      // find smallest entry
      for (j=i;j<N_z+1;j++)
      {
	  if (z_help[j] < mini)
	  {
	      mini = z_help[j];
	      found = j;
	  }
      }
      // change entries
      mini = z_help[i];
      z_help[i] = z_help[found];
      z_help[found] = mini;
  }
  z_min = z_help[0];
  z_max = z_help[N_z];

  // parameters for the loops
  // total number of nodes of the 4d grid
  Nodes = (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1);
  // number of nodes of the 2d plane (spatial coordinates x and y)
  N2 = (N_x+1)*(N_y+1);
  // number of nodes of the 3d cube = "spatial part" of the 4d grid
  N3 = N2*(N_z+1);

  // stretching parameter for grid type 2
  gamma = TDatabase::ParamDB->CHANNEL_GRID_STRETCH;

  //arrays for the extension from 3d to 4d
  double *x_3d_cube, *y_3d_cube, *z_3d_cube;

  x_3d_cube = new double[N3];  memset(x_3d_cube,0,N3*SizeOfDouble);
  y_3d_cube = new double[N3];  memset(y_3d_cube,0,N3*SizeOfDouble);
  z_3d_cube = new double[N3];  memset(z_3d_cube,0,N3*SizeOfDouble);

  switch(TDatabase::ParamDB->GRID_TYPE)
  {
    // equidistant grid
    case 0:
      // a_help
      for ( i=0 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a_min + ( ((a_max-a_min)/(N_a)) * i );
      }
      break;

    // non-equidistant grid, based on the cosinus function in the a-direction
    case 1:
      // fill auxiliary vectors with a distribution
      // a_help
      for ( i=0 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a_max - (a_max-a_min) * cos(i*Pi/(2*N_a));
      }
      break;

    // non-equidistant grid, based on the hyperbolic tangent function in the z-direction
    // stretching parameter is gamma
    case 2:
      // fill auxiliary vectors with a distribution
      // a_help
      for ( i=0 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a_max + (a_max-a_min) * (tanh( gamma*(((double)i/(double)N_a)-1.)))/(tanh(gamma));
      }
      break;

      // equidistant with respect to mass
      // for urea
    case 3:
	// equidistant with power 3
	a1 = a_min*a_min*a_min;
	a2 = a_max*a_max*a_max;
	a3 = 1.0/3.0;
	for ( i=0 ; i<(N_a+1) ; i++ )
	{
	    a_help[i] = a1 + ( ((a2-a1)/(N_a)) * i );
	    // third power
	    a_help[i] = pow(a_help[i],a3);
	}
	break;
    case 4: 
	// test for windtunnel, first coordinate 0, then to r_min, then equi-distant
	a_help[0] = 0;
	a_min = 5e-6/TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
	for ( i=1 ; i<(N_a+1) ; i++ )
	{
	    a_help[i] = a_min + ( ((a_max-a_min)/(N_a-1)) * (i-1) );
	}
	break;
    case 5:
        // transformed grid(from discretized with respect to diam into disc with respect to mass)
      for ( i=0 ; i<(N_a+1); i++ )
	{
	    a_help[i] = input[i];
	}
       break;
    default:
      OutPut(" error : wrong grid_type in grid_generator_4d " << endl);
      exit(-1);
  }
  // fill the "cube"-vectors
  // x_3d_cube
  for ( i=0 ; i<N3 ; i++ )
  {
    x_3d_cube[i] = x_help[(i%(N_x+1))];
  }

  // y_3d_cube
  for ( k=0 ; k<(N_z+1) ; k++ )
  {
    for ( j=0 ; j<(N_y+1) ; j++ )
    {
      for ( i=0 ; i<(N_x+1) ; i++  )
      {
        y_3d_cube[(i+(j*(N_x+1)))+(k*N2)] = y_help[j];
      }
    }
  }
 
  // z_3d_cube
  for ( j=0 ; j<(N_z+1) ; j++ )
  {
    for ( i=0 ; i<N2 ; i++ )
    {
      z_3d_cube[i+(j*N2)] = z_help[j];
    }
  }

  x_coord = new double[Nodes];
  memset(x_coord,0,Nodes*SizeOfDouble);
  y_coord = new double[Nodes];
  memset(y_coord,0,Nodes*SizeOfDouble);
  z_coord = new double[Nodes];
  memset(z_coord,0,Nodes*SizeOfDouble);
  a_coord = new double[Nodes];
  memset(a_coord,0,Nodes*SizeOfDouble);
  
  // fill the output vectors
  for ( i=0 ; i<Nodes ; i++ )
  {
    x_coord[i] = x_3d_cube[i%N3];
    y_coord[i] = y_3d_cube[i%N3];
    z_coord[i] = z_3d_cube[i%N3];
  }
 
  for ( j=0 ; j<(N_a+1) ; j++ )
  {
    for ( i=0 ; i<N3 ; i++ )
    {
      a_coord[i+(j*N3)] = a_help[j];
    }
  }

  // generate the vector with the coordinates of the layers in a-direction
  // the entries are all positive -> in the dim-less form is a_min = d_p_min and a_max = 1.0
  for ( i=0 ; i<(N_a+1) ; i++ )
  {
    a_layers_coord[i] = a_help[i];
  }

  delete z_3d_cube;
  delete y_3d_cube;
  delete x_3d_cube;

  delete a_help;
  delete z_help;
  delete y_help;
  delete x_help;
}

void grid_generator_4d_pipe(int N_x, int N_y, int N_z, int N_a,
                            double* &x_coord, double* &y_coord, 
                            double* &z_coord, 
                            double* &x_coord_pipe, double* &y_coord_pipe, 
                            double* &z_coord_pipe)
{
  
  int i, j, N3 = (N_x+1)*(N_y+1)*(N_z+1), Nodes = N3*(N_a+1);
  
  x_coord_pipe = new double[Nodes];
  memset(x_coord_pipe,0,Nodes*SizeOfDouble);
  y_coord_pipe = new double[Nodes];
  memset(y_coord_pipe,0,Nodes*SizeOfDouble);
  z_coord_pipe = new double[Nodes];
  memset(z_coord_pipe,0,Nodes*SizeOfDouble);
  
  
  for(i=0 ; i<N3 ; i++)
  {
    CoordsTrafo_SquareToCircle(y_coord[i], z_coord[i], y_coord_pipe[i], z_coord_pipe[i]);
    x_coord_pipe[i] = x_coord[i];
    for (j=1 ; j<(N_a+1) ; j++)
    {
      x_coord_pipe[i+j*N3] = x_coord_pipe[i];
      y_coord_pipe[i+j*N3] = y_coord_pipe[i];
      z_coord_pipe[i+j*N3] = z_coord_pipe[i];
    }
  }
}

// ========================================================================
// calculate circle coordinates for given square coordinates
// ========================================================================

void CoordsTrafo_SquareToCircle(double x, double y, double& nx, double& ny)
{
  double phi, r, pi = 3.1415926535897, eps = 1e-8;

  //calculate radius r and angle phi
  r = fmax(fabs(x),fabs(y));
  // case x = y = 0 
  if (fabs(r) < eps)
  {
    nx = 0;  
    ny = 0;
    return;
  }
  else //r>0
  {
    // 1st sector
    if(fabs(y) <= fabs(x) && x > 0)
    { phi = (pi/4)*y/r; }
    
    // 2nd sector
    if(fabs(y) > fabs(x) && y > 0)
    { phi = pi/2 - (pi/4)*x/r;}
      
    // 3rd sector
    if(fabs(y) <= fabs(x) && x < 0)
    { phi = pi - (pi/4)*y/r;}
    
    // 4th sector
    if(fabs(y) > fabs(x) && y < 0)
    { phi = pi*3/2 + (pi/4)*x/r;}

    // calculate x and y coordinate in the circle
    nx = r*cos(phi);
    ny = r*sin(phi);
    return;
  }
}

void grid_generator_5d(TCollection *coll,
int &N_x,
int &N_y,
int &N_z,
double &x_min, double &x_max, double &y_min, double &y_max,
double &z_min, double &z_max,
double a_min, double a_max, double b_min, double b_max, int N_a, int N_b,
double* &x_coord, double* &y_coord, double* &z_coord, double* &a_coord,double* &b_coord,
double *a_layers_coord, double *b_layers_coord)
{
  // check consistency
  if ((a_min >= a_max )|| (N_a <= 0) )
  {
    OutPut(" error : wrong input in grid_generator_5d " << endl);
    exit(-1);
  }
  if ((b_min >= b_max )|| (N_b <= 0) )
  {
    OutPut(" error : wrong input in grid_generator_5d " << endl);
    //exit(-1);
  }
  OutPut("maximal values of internal coordinates " << a_max << " " << b_max  << endl);

  // integer variables for the loops
  int i, j,j1, k, N_Cells, N_Vertices, found, Nodes, N2, N3,N4, max_dir = 1025;
  double eps = 1e-5, mini, x, y, z, gamma, a1, a2, a3, b1, b2, b3;
  TBaseCell *cell;

  // arrays for the grid -> this allows to implement different distributions of grid points in every direction
  double *x_help, *y_help, *z_help, *a_help, *b_help;

  x_help = new double[max_dir];  memset(x_help,0,max_dir*SizeOfDouble);
  y_help = new double[max_dir];  memset(y_help,0,max_dir*SizeOfDouble);
  z_help = new double[max_dir];  memset(z_help,0,max_dir*SizeOfDouble);
  a_help = new double[N_a+1];  memset(a_help,0,(N_a+1)*SizeOfDouble);
  b_help = new double[N_b+1];  memset(b_help,0,(N_b+1)*SizeOfDouble);
  //neu grid
  //neu grid von ARAM
  double input[94] =
  {
    0,
    0.000488281,
    0.000615194,
    0.000775097,
    0.000976561,
    0.00123039,
    0.00155019,
    0.00195312,
    0.00246078,
    0.00310039,
    0.00390625,
    0.00447154,
    0.00492156,
    0.0053016,
    0.00563379,
    0.00593084,
    0.00620079,
    0.00667959,
    0.00709813,
    0.00747238,
    0.0078125,
    0.00841576,
    0.00894308,
    0.00941462,
    0.00984313,
    0.0106032,
    0.0112676,
    0.0118617,
    0.0124016,
    0.0133592,
    0.0141962,
    0.0149448,
    0.015625,
    0.0168315,
    0.0178862,
    0.0188292,
    0.0196863,
    0.0204745,
    0.0212064,
    0.0218909,
    0.0225352,
    0.0231445,
    0.0237233,
    0.0242752,
    0.0248031,
    0.0267184,
    0.0283925,
    0.0298896,
    0.03125,
    0.0336631,
    0.0357723,
    0.0393725,
    0.0450703,
    0.0496062,
    0.0534368,
    0.056785,
    0.0597791,
    0.0625,
    0.0673261,
    0.0715446,
    0.0753169,
    0.0787451,
    0.0818982,
    0.0848255,
    0.0875637,
    0.0901406,
    0.0948934,
    0.0992125,
    0.106873,
    0.11357,
    0.119558,
    0.125,
    0.134652,
    0.143089,
    0.150634,
    0.15749,
    0.169651,
    0.180281,
    0.189787,
    0.198425,
    0.213747,
    0.22714,
    0.239117,
    0.25,
    0.269304,
    0.286179,
    0.301268,
    0.31498,
    0.360562,
    0.39685,
    0.5,
    0.629961,
    0.793701,
    1
  };

  // compute coordinates of the 3d domain
  // initialize
  for (i=0;i<max_dir;i++)
    x_help[i] = -47111174;
  for (i=0;i<max_dir;i++)
    y_help[i] = -47111174;
  for (i=0;i<max_dir;i++)
    z_help[i] = -47111174;

  // number of mesh cells
  N_Cells = coll->GetN_Cells();
  for (i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    // number of vertices
    N_Vertices = cell->GetN_Vertices();
    for ( j=0; j<N_Vertices; j++)
    {
      // get coordinates of the vertex
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      z = cell->GetVertex(j)->GetZ();
      // check if they are already included in the arrays
      found = 0;
      for (k=0;k<max_dir;k++)
      {
        // coordinate already included
        if (fabs(x-x_help[k]) < eps)
        {
          found = 1;
          break;
        }
        // new entry
        if (fabs(x_help[k]+47111174) < eps)
        {
          x_help[k] = x;
          found = 1;
          break;
        }
      }
      if (!found)
      {
        OutPut("Error in grid_generation_5d, x coordinate !!!" << endl);
        exit(4711);
      }
      found = 0;
      for (k=0;k<max_dir;k++)
      {
        // coordinate already included
        if (fabs(y-y_help[k]) < eps)
        {
          found = 1;
          break;
        }
        // new entry
        if (fabs(y_help[k]+47111174) < eps)
        {
          y_help[k] = y;
          found = 1;
          break;
        }
      }
      if (!found)
      {
        OutPut("Error in grid_generation_5d, y coordinate !!!" << endl);
        exit(4711);
      }
      found = 0;
      for (k=0;k<max_dir;k++)
      {
        // coordinate already included
        if (fabs(z-z_help[k]) < eps)
        {
          found = 1;
          break;
        }
        // new entry
        if (fabs(z_help[k]+47111174) < eps)
        {
          z_help[k] = z;
          found = 1;
          break;
        }
      }
      if (!found)
      {
        OutPut("Error in grid_generation_5d, z coordinate !!!" << endl);
        exit(4711);
      }
    }
  }

  // compute actual length of arrays
  for (i=0;i<max_dir;i++)
  {
    if (fabs(x_help[i]+47111174) < eps)
      break;
  }
  N_x = i-1;
  for (i=0;i<max_dir;i++)
  {
    if (fabs(y_help[i]+47111174) < eps)
      break;
  }
  N_y = i-1;
  for (i=0;i<max_dir;i++)
  {
    if (fabs(z_help[i]+47111174) < eps)
      break;
  }
  N_z = i-1;

  // bring arrays in correct order
  for (i=0;i<N_x+1;i++)
  {
    mini = 47111147.0;
    // find smallest entry
    for (j=i;j<N_x+1;j++)
    {
      if (x_help[j] < mini)
      {
        mini = x_help[j];
        found = j;
      }
    }
    // change entries
    mini = x_help[i];
    x_help[i] = x_help[found];
    x_help[found] = mini;
  }
  x_min = x_help[0];
  x_max = x_help[N_x];

  for (i=0;i<N_y+1;i++)
  {
    mini = 47111147.0;
    // find smallest entry
    for (j=i;j<N_y+1;j++)
    {
      if (y_help[j] < mini)
      {
        mini = y_help[j];
        found = j;
      }
    }
    // change entries
    mini = y_help[i];
    y_help[i] = y_help[found];
    y_help[found] = mini;
  }
  y_min = y_help[0];
  y_max = y_help[N_y];

  for (i=0;i<N_z+1;i++)
  {
    mini = 47111147.0;
    // find smallest entry
    for (j=i;j<N_z+1;j++)
    {
      if (z_help[j] < mini)
      {
        mini = z_help[j];
        found = j;
      }
    }
    // change entries
    mini = z_help[i];
    z_help[i] = z_help[found];
    z_help[found] = mini;
  }
  z_min = z_help[0];
  z_max = z_help[N_z];

  // parameters for the loops
  // total number of nodes of the 4d grid
  Nodes = (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1)*(N_b+1);
  // number of nodes of the 2d plane (spatial coordinates x and y)
  N2 = (N_x+1)*(N_y+1);
  // number of nodes of the 3d cube = "spatial part" of the 4d grid
  N3 = N2*(N_z+1);
  N4 = N3*(N_a+1);
  // stretching parameter for grid type 2
  gamma = TDatabase::ParamDB->CHANNEL_GRID_STRETCH;

  //arrays for the extension from 3d to 4d
  double *x_3d_cube, *y_3d_cube, *z_3d_cube;

  x_3d_cube = new double[N3];  memset(x_3d_cube,0,N3*SizeOfDouble);
  y_3d_cube = new double[N3];  memset(y_3d_cube,0,N3*SizeOfDouble);
  z_3d_cube = new double[N3];  memset(z_3d_cube,0,N3*SizeOfDouble);

  switch(TDatabase::ParamDB->GRID_TYPE_1)
  {
    // equidistant grid
    case 0:
      // a_help
      for ( i=0 ; i<(N_a+1) ; i++ )
      {
        a_min=TDatabase::ParamDB->KDP_D_P_0;
        a_help[i] = a_min + ( ((a_max-a_min)/(N_a)) * i );
        //a_help[i] = a_min/a_max + ( ((a_max-a_min)/a_max*(N_a)) * i );
        //OutPut("a_min " << a_min <<endl);

        //a_help[i] = ( (a_max)/(N_a) * i );
      }
      break;

      // non-equidistant grid, based on the cosinus function in the a-direction
    case 1:
      // fill auxiliary vectors with a distribution
      // a_help
      for ( i=0 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a_max - (a_max-a_min) * cos(i*Pi/(2*N_a));
      }
      break;

      // non-equidistant grid, based on the hyperbolic tangent function in the z-direction
      // stretching parameter is gamma
    case 2:
      // fill auxiliary vectors with a distribution
      // a_help
      for ( i=0 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a_max + (a_max-a_min) * (tanh( gamma*(((double)i/(double)N_a)-1.)))/(tanh(gamma));
      }
      break;

      // equidistant with respect to mass
      // for urea
    case 3:
      // equidistant with power 3
      a1 = a_min*a_min*a_min;
      a2 = a_max*a_max*a_max;
      a3 = 1.0/3.0;
      for ( i=0 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a1 + ( ((a2-a1)/(N_a)) * i );
        // third power
        a_help[i] = pow(a_help[i],a3);
      }
      break;
    case 4:
      // test for windtunnel, first coordinate 0, then to r_min, then equi-distant
      a_help[0] = 0;
      a_min = 5e-6/TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
      for ( i=1 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a_min + ( ((a_max-a_min)/(N_a-1)) * (i-1) );
      }
      break;
    case 5:
      // transformed grid(from discretized with respect to diam into disc with respect to mass)
      for ( i=0 ; i<(N_a+1); i++ )
      {
        a_help[i] = input[i];
      }
      break;
      // equidistant grid save first interval
    case 10:
      // a_help
      a_help[0] = a_min;
      a_help[1] = a_min +  ((a_max-a_min)/(N_a))/2.0;

      for ( i=2 ; i<(N_a+1) ; i++ )
      {
        a_help[i] = a_help[1]  + ( ((a_max-a_min)/(N_a)) * (i-1) );
      }
      break;
    default:
      OutPut(" error : wrong grid_type in grid_generator_5d " << endl);
      exit(-1);
  }
  switch(TDatabase::ParamDB->GRID_TYPE_2)
  {
    // equidistant grid
    case 0:
      // b_help
      for ( i=0 ; i<(N_b+1) ; i++ )
      {
        b_min=TDatabase::ParamDB->KDP_D_P_0_2;
        b_help[i] = b_min + ( ((b_max-b_min)/(N_b)) * i );
        //b_help[i] = b_min/b_max + ( ((b_max-b_min)/(b_max*N_b)) * i );
      }
      break;

      // non-equidistant grid, based on the cosinus function in the a-direction
    case 1:
      // fill auxiliary vectors with a distribution
      // a_help
      for ( i=0 ; i<(N_b+1) ; i++ )
      {
        b_help[i] = b_max - (b_max-b_min) * cos(i*Pi/(2*N_b));
      }
      break;

      // non-equidistant grid, based on the hyperbolic tangent function in the z-direction
      // stretching parameter is gamma
    case 2:
      // fill auxiliary vectors with a distribution
      // a_help
      for ( i=0 ; i<(N_b+1) ; i++ )
      {
        b_help[i] = b_max + (b_max-b_min) * (tanh( gamma*(((double)i/(double)N_b)-1.)))/(tanh(gamma));
      }
      break;

      // equidistant with respect to mass
      // for urea
    case 3:
      // equidistant with power 3
      b1 = b_min*b_min*b_min;
      b2 = b_max*b_max*b_max;
      b3 = 1.0/3.0;
      for ( i=0 ; i<(N_b+1) ; i++ )
      {
        b_help[i] = b1 + ( ((b2-b1)/(N_b)) * i );
        // third power
        b_help[i] = pow(b_help[i],b3);
      }
      break;
    case 4:
      // test for windtunnel, first coordinate 0, then to r_min, then equi-distant
      b_help[0] = 0;
      b_min = 5e-6/TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
      for ( i=1 ; i<(N_b+1) ; i++ )
      {
        b_help[i] = b_min + ( ((b_max-b_min)/(N_b-1)) * (i-1) );
      }
      break;
    case 5:
      // transformed grid(from discretized with respect to diam into disc with respect to mass)
      for ( i=0 ; i<(N_b+1); i++ )
      {
        b_help[i] = input[i];
      }
      break;
      // equidistant grid save first interval
    case 10:
      // b_help
      b_help[0] = b_min;
      b_help[1] = b_min +  ((b_max-b_min)/(N_b))/2.0;

      for ( i=2 ; i<(N_b+1) ; i++ )
      {
        b_help[i] = b_help[1]  + ( ((b_max-b_min)/(N_b)) * (i-1) );
      }
      break;
    default:
      OutPut(" error : wrong grid_type in grid_generator_5d " << endl);
      exit(-1);
  }
  // fill the "cube"-vectors
  // x_3d_cube
  for ( i=0 ; i<N3 ; i++ )
  {
    x_3d_cube[i] = x_help[(i%(N_x+1))];
  }

  // y_3d_cube
  for ( k=0 ; k<(N_z+1) ; k++ )
  {
    for ( j=0 ; j<(N_y+1) ; j++ )
    {
      for ( i=0 ; i<(N_x+1) ; i++  )
      {
        y_3d_cube[(i+(j*(N_x+1)))+(k*N2)] = y_help[j];
      }
    }
  }

  // z_3d_cube
  for ( j=0 ; j<(N_z+1) ; j++ )
  {
    for ( i=0 ; i<N2 ; i++ )
    {
      z_3d_cube[i+(j*N2)] = z_help[j];
    }
  }

  x_coord = new double[Nodes];
  memset(x_coord,0,Nodes*SizeOfDouble);
  y_coord = new double[Nodes];
  memset(y_coord,0,Nodes*SizeOfDouble);
  z_coord = new double[Nodes];
  memset(z_coord,0,Nodes*SizeOfDouble);
  a_coord = new double[Nodes];
  memset(a_coord,0,Nodes*SizeOfDouble);
  b_coord = new double[Nodes];
  memset(b_coord,0,Nodes*SizeOfDouble);

  // fill the output vectors
  for ( i=0 ; i<Nodes ; i++ )
  {
    x_coord[i] = x_3d_cube[i%N3];
    y_coord[i] = y_3d_cube[i%N3];
    z_coord[i] = z_3d_cube[i%N3];
  }

  for ( j1=0 ; j1<(N_b+1) ; j1++ )
  {
    for ( j=0 ; j<(N_a+1) ; j++ )
    {
      for ( i=0 ; i<N3 ; i++ )

      {
        a_coord[i+(j*N3)+(j1*N4)]=a_help[j];
        b_coord[i+(j*N3)+(j1*N4)]=b_help[j1];
        //OutPut("a" << a_coord[i+(j*N3)+(j1*N4)] << endl);
        // OutPut("b" << b_coord[i+(j*N3)+(j1*N4)] << endl);

      }

    }
  }

  // generate the vector with the coordinates of the layers in a-direction
  // the entries are all positive -> in the dim-less form is a_min = d_p_min and a_max = 1.0
  for ( i=0 ; i<(N_a+1) ; i++ )
  {
    a_layers_coord[i] = a_help[i];
    //OutPut("a" << a_layers_coord[i] << endl);
  }
  // generate the vector with the coordinates of the layers in b-direction
  for ( j=0 ; j<(N_b+1) ; j++ )
  {
    b_layers_coord[j] = b_help[j];
    //OutPut("b" << a_layers_coord[i] << endl);
  }

  delete z_3d_cube;
  delete y_3d_cube;
  delete x_3d_cube;

  delete b_help;
  delete a_help;
  delete z_help;
  delete y_help;
  delete x_help;
}




/***************************************************************************/
//
// computes boundary conditions for PSD on the inflow boundary wrt the flow
//
// boundary conditions have to be described on the closure of the inflow
//
/***************************************************************************/

int PSD_bound_cound_from_velo_inflow(double x, double y, double z)
{
    int value=0;
    double eps = 1e-6, r2 = TDatabase::ParamDB->P6;
    
    // inflow from left
    if (fabs(x)<eps)
    {
	if ((fabs(y-0.5)<r2+eps)&&(fabs(z-0.5)<r2+eps))
	    value = 3;
    }
   // inflow from right 
    if (fabs(1 - x)<eps)
    {
       if ((fabs(y-0.5)<r2+eps)&&(fabs(z-0.5)<r2+eps))
	    value = 1;
    }
    return value;
}

/****************************************************************************************
*                                                                                       *
*                          Part II.1 : FDM (explicit with upwindig)                     *
*                         -----------------------------------------                     *
*                                                                                       *
****************************************************************************************/


/****************************************************************************************
 *                                                                                      *
 * assemble the matrix corresponding to the 4d finite-difference-method                 *
 *                                                                                      *
 ***************************************************************************************/

void Bulk_FWE_FDM_Upwind_4D(TCollection *coll,
			    TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            double *f_old,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            double x_min, double x_max, double y_min, double y_max,
                            double z_min, double z_max, double a_min, double a_max,
                            double *velo1, double *velo2, double *velo3,
                            double *concent_C_array,
			    int *correspond_3dgrid
)
{
  int i, ii, j, N2, N3, N4, alpha, beta, gamma, no_of_3dcell;
  int maxind;
  int very_first = 0;

  double B_c_C, concent_C_array_val, maxsol, G_c_C_val, val;
  double velocity1_array_val, velocity2_array_val, velocity3_array_val;
  double values[4];
  double *f_new, *derx_val;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;

  // model constants
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double C_2 = TDatabase::ParamDB->BULK_C_2;
  double d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double k_nuc = TDatabase::ParamDB->BULK_k_nuc;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;
  double factor_G;

  // computed model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);

  // the arrays for the velocity field (velo1, velo2 and velo3) are only
  // filled one-time; at the very first computation
  if (fabs (velo1[0] + 4711) < 1e-6)
  {
    very_first++;
    OutPut("very first computation of f" << endl);
  }

  // number of unknowns in 2D (on one layer)
  N2 = (N_x+1)*(N_y+1);
  // number of unknowns in 3D (one cube)
  N3 = (N_x+1)*(N_y+1)*(N_z+1);
  // number of unknowns in 4D (total grid)
  N4 = (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1);

  // array for concentration of species C
  memset(concent_C_array, 0, N3*SizeOfDouble);
  // array for derivatives of f with respect to x,y,z,a
  derx_val = new double[N4];
  memset(derx_val, 0, N4*SizeOfDouble);
  // array for new values of the particle size distribution f
  f_new = new double[N4];
  memset(f_new, 0, N4*SizeOfDouble);

  // discretization of PBE with FDM
  // loop over all nodes of the FDM grid
  for ( i=0 ; i<N4 ; i++ )
  {
    // node is on the first cube (a_coord=0)
    // in this cube, the velocity vector will be filled
    // since this vector is independent of a
    if (i < N3)
    {
	// CHECKED AT 08/10/09
	// find grid cell of 3d grid
	ii = i;
	// top of "first cube" = z_max
	if( ii>= N2*N_z )
	{
	    ii = ii-N2;
	}
	// treat right d.o.f. = x_max seperately
	if ( ((ii+1)%(N_x+1)==0) )
	{
	    ii = ii-1;
	}
	// treat upper d.o.f. =  y_max seperately
	if ( ii%N2 >= (N_y*(N_x+1)) )
	{
	    ii = ii-(N_x+1);
	}
	// level in z-direction 
	gamma = (int)(ii/N2);
	ii -= gamma*N2;
	// level in y-direction 
	alpha = (int)(ii/(N_x+1)+1e-6 );
	// level in x-direction
	beta = ii%(N_x+1);
	no_of_3dcell = correspond_3dgrid[gamma*N_x*N_y+ alpha * N_x + beta];

      // find velocity
      velocity1->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				   x_coord[i],y_coord[i],z_coord[i],values);
      velo1[i] = values[0];
      velocity1_array_val = velo1[i];

      velocity2->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				   x_coord[i],y_coord[i],z_coord[i],values);
      velo2[i] = values[0];
      velocity2_array_val = velo2[i];
      
      velocity3->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				   x_coord[i],y_coord[i],z_coord[i],values);
      velo3[i] = values[0];
      velocity3_array_val = velo3[i];

      // fill the array for the concentrations of C
      // since this concentration does not depend on a
      // c_C
      concent_C->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				   x_coord[i],y_coord[i],z_coord[i],values);
      concent_C_array[i] = values[0];
      concent_C_array_val = values[0];
    }
    // node is not on the first cube
    else
    {
      // the corresponding node on the first cube
      ii = i - N3*((int)(i/N3));
      // compute the value of the velocity for the corresponding (x,y,z) coordinates
      velocity1_array_val = velo1[ii];
      velocity2_array_val = velo2[ii];
      velocity3_array_val = velo3[ii];
      // compute the value of the concentration of C for the corresponding (x,y,z) coordinates
      concent_C_array_val = concent_C_array[ii];
    }

    if ( TDatabase::ParamDB->BULK_GROWTH_RATE==2 )
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(a_coord[i]*d_p_max)));
    }
    else
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty);
    }
    
    // compute the coefficients corresponding to the 4d finite-difference-method (9-point)
    // simple upwind scheme


    // compute the term A_x
    if (velocity1_array_val >= 0)
    {
      // not on boundary x = x_min -> in the "x-sense" exists a left neighbour
      if ( i%(N_x+1)!=0 )
        derx_val[i] += velocity1_array_val*(f_old[i]-f_old[i-1])/(fabs(x_coord[i]-x_coord[i-1]));
    }
    else
    {
      // not on boundary x = x_max -> in the "x-sense" exists a right neighbour
      if ( (i+1)%(N_x+1)!=0 )
        derx_val[i] += velocity1_array_val*(f_old[i+1]-f_old[i])/(fabs(x_coord[i+1]-x_coord[i]));
    }
    
    // compute the term A_y
    if (velocity2_array_val >= 0)
    {
      // not on boundary y = y_min -> in the "y-sense" exists a left neighbour
      if ((i%((N_x+1)*(N_y+1)))>N_x )
       derx_val[i] += velocity2_array_val*(f_old[i]-f_old[i-(N_x+1)])/(fabs(y_coord[i]-y_coord[i-(N_x+1)]));
    }
    else
    {
      // not on boundary y = y_max -> in the "y-sense" exists a right neighbour
      if ((i%((N_x+1)*(N_y+1)))<((N_x+1)*N_y) )
       derx_val[i] += velocity2_array_val*(f_old[i+(N_x+1)]-f_old[i])/(fabs(y_coord[i+(N_x+1)]-y_coord[i]));
    }
    
    // compute the term A_z
    if (velocity3_array_val >= 0)
    {
      // not on boundary z = z_min -> in the "z-sense" exists a left neighbour
      if ( i%N3 >= N2 )
       derx_val[i] += velocity3_array_val*(f_old[i]-f_old[i-N2])/(fabs(z_coord[i]-z_coord[i-N2]));
    }
    else
    {
      // not on boundary z = z_max -> in the "z-sense" exists a right neighbour
      if ( i%N3 < N2*N_z )
       derx_val[i] += velocity3_array_val*(f_old[i+N2]-f_old[i])/(fabs(z_coord[i+N2]-z_coord[i]));
    }
    
    // compute the term C
    if (G_c_C_val >= 0)
    {
      // not on boundary a = a_min -> in the "a-sense" exists a left neighbour
      if ( i >= N3 )
        derx_val[i] += G_c_C_val*(f_old[i]-f_old[i-N3])/(fabs(a_coord[i]-a_coord[i-N3]));
    }
    else
    {
      // not on boundary a = a_max -> in the "a-sense" exists a left neighbour
      if ( i < N3*N_a )
        derx_val[i] += G_c_C_val*(f_old[i+N3]-f_old[i])/(fabs(a_coord[i+N3]-a_coord[i]));
    }
    
    // compute new particle size distribution
    f_new[i] = f_old[i] - deltat * derx_val[i];
    
    // set Dirichlet boundary conditions
    // if convection is positive at the bottom ( = cube with a = a_min )
    if (i<N3)
    {
      if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
      {
        // G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(z_coord[i]*d_p_max)));
        G_c_C_val = k_g*(c_C_infty*concent_C_array_val- c_C_infty_sat*exp(C_2/d_p_0));
      }
      else
      {
        G_c_C_val = k_g*(c_C_infty*concent_C_array_val-c_C_infty_sat);
      }
      // compute G*n, n=(0,0,0,-1);
      if (G_c_C_val*f_infty > 1e-10)
      {
        // compute rate of nucleation
        B_c_C = k_nuc*pow(c_C_infty*(concent_C_array_val - 1),5);
        // truncate negative values
        if (B_c_C < 0)
          B_c_C = 0;
        // compute new particle size distribution
        f_new[i] = B_c_C/(G_c_C_val*f_infty);
      }
    }
    // set Dirichlet boundary conditions
    // if convection is positive at the top
    if (i>=N3*N_a )
    {
      if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
      {
        //G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(z_coord[i]*d_p_max)));
        G_c_C_val = k_g*(c_C_infty*concent_C_array_val- c_C_infty_sat*exp(C_2/d_p_0));
      }
      else
      {
        G_c_C_val = k_g*(c_C_infty*concent_C_array_val-c_C_infty_sat);
      }
      // compute G*n, n=(0,0,0,1);
      if (G_c_C_val*f_infty < 0)
      {
	  f_new[i] = 0.0;
      }
    }

    // set Dirichlet boundary conditions
    // inflow from the left x = x_min (left) or right x = x_max (right)
    if ((( i%(N_x+1)==0 )  ||  ((i+1)%(N_x+1)==0 ))&&(i>N3))
    {
	val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i], z_coord[i]);
	if (val)
	     f_new[i] = 0.0;
    }
  }

  // copy new particle size distribution into array for old one
  Dcopy(N4, f_new, f_old);

  maxsol =  0;
  maxind = -4711;

  // cut undershoots
  for ( i=0 ; i<N4 ; i++ )
  {
    if (f_old[i] > maxsol)
    {
      maxsol = f_old[i];
      maxind = i;
    }
    if (f_old[i]<0)
    {
      f_old[i] = 0;
    }
  }

  // free allocated memory
  delete f_new;
  delete derx_val;
}

void Bulk_RKV_FDM_4D(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3, TFEFunction3D *concent_C,
double *f_old, double **stages,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
double *a_coord,
double *velo1, double *velo2, double *velo3,
double *concent_C_array, int *correspond_3dgrid)
{
  int i, ii, N2, N3, N4, maxind, val, N1_[4], N_, N_stages;
  int very_first = 0, disctype, time_disc;
  int alpha, beta, gamma, no_of_3dcell;
  int *offset_ = NULL, *offset1_;
  double B_c_C, velocity1_array_val, velocity2_array_val, velocity3_array_val;
  double concent_C_array_val, maxsol, G_c_C_val;
  double values[4], t1, t2, coeff[4];
  double *current_stage_fdm, *sol_curr;
  double *coordinates[4];
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  double time = TDatabase::TimeDB->CURRENTTIME;

  // model constants
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double C_2 = TDatabase::ParamDB->BULK_C_2;
  double d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double k_nuc = TDatabase::ParamDB->BULK_k_nuc;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;
  double factor_G;

  // computed model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);

  // very first computation
  if (fabs (velo1[0] + 4711) < 1e-6)
  {
    very_first++;
    OutPut("very first computation of f" << endl);
  }
  // number of unknowns in 2D (on one layer)
  N2 = (N_x+1)*(N_y+1);
  // number of unknowns in 3D (one cube)
  N3 = (N_x+1)*(N_y+1)*(N_z+1);
  // number of unknowns in 4D (total grid)
  N4 = (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1);

  // save parameters
  disctype = TDatabase::ParamDB->DISCTYPE;
  time_disc = TDatabase::TimeDB->TIME_DISC;

  // array for concentration of species C
  memset(concent_C_array, 0, N3*SizeOfDouble);
  // array for solution of current state
  sol_curr = stages[5];

  TDatabase::ParamDB->DISCTYPE = TDatabase::ParamDB->PB_DISC_TYPE;
  TDatabase::TimeDB->TIME_DISC = TDatabase::ParamDB->PB_TIME_DISC;

  N_stages = GetN_SubSteps();

  N1_[0] = N_x+1;
  N1_[1] = N_y+1;
  N1_[2] = N_z+1;
  N1_[3] = N_a+1;
  coordinates[0] = x_coord;
  coordinates[1] = y_coord;
  coordinates[2] = z_coord;
  coordinates[3] = a_coord;

  InitializeConvectiveTermFDM(4, offset_, offset1_, N1_);

  for (N_ = 0; N_ < N_stages; N_++)
  {
    SetTimeDiscParameters(1);
    for (i=0;i<N_stages;i++)
      OutPut(" A("<<N_<<","<<i<<") = " << TDatabase::TimeDB->RK_A[N_][i]);
    OutPut(" : c("<<N_<<") = " << TDatabase::TimeDB->RK_c[N_]);
    OutPut(" : b("<<N_<<") = " << TDatabase::TimeDB->RK_b[N_] << endl);
    current_stage_fdm = stages[N_];
    memset(current_stage_fdm, 0, N4 * SizeOfDouble);
    // initialize current stage
    memcpy(sol_curr,f_old,N4 * SizeOfDouble);
    // add previous stages
    for (i=0;i<N_;i++)
      Daxpy(N4, deltat*TDatabase::TimeDB->RK_A[N_][i], stages[i], sol_curr);
    // compute next stage
    // discretization of PBE with FDM
    // loop over all nodes of the FDM grid
    for ( i=0 ; i<N4 ; i++ )
    {
      // node is on the first cube (a_coord=0)
      // in this cube, the velocity vector will be filled
      // since this vector is independent of a
      if (i < N3)
      {
        // CHECKED AT 08/10/09
        // find grid cell of 3d grid
        ii = i;
        // top of "first cube" = z_max
        if( ii>= N2*N_z )
        {
          ii = ii-N2;
        }
        // treat right d.o.f. = x_max seperately
        if ( ((ii+1)%(N_x+1)==0) )
        {
          ii = ii-1;
        }
        // treat upper d.o.f. =  y_max seperately
        if ( ii%N2 >= (N_y*(N_x+1)) )
        {
          ii = ii-(N_x+1);
        }
        // level in z-direction
        gamma = (int)(ii/N2);
        ii -= gamma*N2;
        // level in y-direction
        alpha = (int)(ii/(N_x+1)+1e-6 );
        // level in x-direction
        beta = ii%(N_x+1);
        no_of_3dcell = correspond_3dgrid[gamma*N_x*N_y+ alpha * N_x + beta];

        // find velocity
        velocity1->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
          x_coord[i],y_coord[i],z_coord[i],values);
        velo1[i] = values[0];
        velocity1_array_val = velo1[i];

        velocity2->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
          x_coord[i],y_coord[i],z_coord[i],values);
        velo2[i] = values[0];
        velocity2_array_val = velo2[i];

        velocity3->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
          x_coord[i],y_coord[i],z_coord[i],values);
        velo3[i] = values[0];
        velocity3_array_val = velo3[i];

        // fill the array for the concentrations of C
        // since this concentration does not depend on a
        // c_C
        concent_C->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
          x_coord[i],y_coord[i],z_coord[i],values);
        concent_C_array[i] = values[0];
        concent_C_array_val = values[0];
      }
      // node is not on the first cube
      else
      {
        // the corresponding node on the first cube
        ii = i - N3*((int)(i/N3));
        // compute the value of the velocity for the corresponding (x,y,z) coordinates
        velocity1_array_val = velo1[ii];
        velocity2_array_val = velo2[ii];
        velocity3_array_val = velo3[ii];
        // compute the value of the concentration of C for the corresponding (x,y,z) coordinates
        concent_C_array_val = concent_C_array[ii];
      }

      if ( TDatabase::ParamDB->BULK_GROWTH_RATE==2 )
      {
        G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(a_coord[i]*d_p_max)));
      }
      else
      {
        G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty);
      }
      coeff[0] = velocity1_array_val;
      coeff[1] = velocity2_array_val;
      coeff[2] = velocity3_array_val;
      coeff[3] = G_c_C_val;

      // compute convection term
      ConvectiveTermFDM(4, i,
        coeff, sol_curr, current_stage_fdm, coordinates,
        offset_, offset1_);

      // set Dirichlet boundary conditions
      // if convection is positive at the bottom ( = cube with a = a_min )
      if (i<N3)
      {
        if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
        {
          // G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(z_coord[i]*d_p_max)));
          G_c_C_val = k_g*(c_C_infty*concent_C_array_val- c_C_infty_sat*exp(C_2/d_p_0));
        }
        else
        {
          G_c_C_val = k_g*(c_C_infty*concent_C_array_val-c_C_infty_sat);
        }
        // compute G*n, n=(0,0,0,-1);
        if (G_c_C_val*f_infty > 1e-10)
        {
          // compute rate of nucleation
          B_c_C = k_nuc*pow(c_C_infty*(concent_C_array_val - 1),5);
          // truncate negative values
          if (B_c_C < 0)
            B_c_C = 0;
          // compute new particle size distribution
          if (N_ == N_stages -1)
            f_old[i] = B_c_C/(G_c_C_val*f_infty);
          current_stage_fdm[i] = 0.0;
        }
      }
      // set Dirichlet boundary conditions
      // if convection is positive at the top
      if (i>=N3*N_a )
      {
        if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
        {
          //G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(z_coord[i]*d_p_max)));
          G_c_C_val = k_g*(c_C_infty*concent_C_array_val- c_C_infty_sat*exp(C_2/d_p_0));
        }
        else
        {
          G_c_C_val = k_g*(c_C_infty*concent_C_array_val-c_C_infty_sat);
        }
        // compute G*n, n=(0,0,0,1);
        if (G_c_C_val*f_infty < 0)
        {
          current_stage_fdm[i] = 0.0;
          if (N_ == N_stages -1)
            f_old[i] = 0.0;
        }
      }

      // set Dirichlet boundary conditions
      // inflow from the left x = x_min (left) or right x = x_max (right)
      if ((( i%(N_x+1)==0 )  ||  ((i+1)%(N_x+1)==0 ))&&(i>N3))
      {
        val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i], z_coord[i]);
        if (val)
        {
          if (N_ == N_stages -1)
            f_old[i] = 0.0;
          current_stage_fdm[i] = 0.0;
        }
      }
    }
  }

  // compute linear combination of stages
  for (i=0;i<N_stages;i++)
  {
    Daxpy(N4, deltat*TDatabase::TimeDB->RK_b[i], stages[i], f_old);
  }

  maxsol =  0;
  maxind = -4711;

  // cut undershoots
  for ( i=0 ; i<N4 ; i++ )
  {
    if (f_old[i] > maxsol)
    {
      maxsol = f_old[i];
      maxind = i;
    }
    if (f_old[i]<0)
    {
      f_old[i] = 0;
    }
  }
  OutPut(time << " maxsol " << maxsol << " maxind " << maxind << endl);

  TDatabase::ParamDB->DISCTYPE = disctype;
  TDatabase::TimeDB->TIME_DISC = time_disc;
  delete[] offset_;
}



/****************************************************************************************
*                                                                                       *
*                          Part II.2 : FDM (implicit Euler with upwindig)               *
*                         -----------------------------------------------               *
*                                                                                       *
****************************************************************************************/


/****************************************************************************************
 *                                                                                      *
 * assemble the matrix corresponding to the 4d finite-difference-method                 *
 *                                                                                      *
 ***************************************************************************************/

void Bulk_BWE_FDM_Upwind_4D(TCollection *coll,
                            TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            double *sol,
                            int *correspond_3dgrid,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            TSquareMatrix3D *mat)
{
  int i, ii, jj, iq, k, N2, N3, N4, diag_index, indices[9];
  int maxind, index, index1, index2;
  int alpha, beta, gamma, no_of_3dcell, N_Entries, range;
  int *col_ptr, *row_ptr;
  int SC_LDS = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;

  double yq, B_c_C, concent_C_array_val, maxsol, G_c_C_val, val, t3;
  double velocity1_array_val, velocity2_array_val, velocity3_array_val;
  double values[4], C_val[4];
  double *velo1, *velo2, *velo3, *concent_C_array;
  double *entries, *rhs;
  double time = TDatabase::TimeDB->CURRENTTIME;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;

  // model constants
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double C_2 = TDatabase::ParamDB->BULK_C_2;
  double d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double k_nuc = TDatabase::ParamDB->BULK_k_nuc;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;
  double factor_G;

  // computed model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);

  // data of the matrix
  entries = mat->GetEntries();
  N_Entries = mat->GetN_Entries();
  col_ptr = mat->GetKCol();
  row_ptr = mat->GetRowPtr();

  // number of unknowns in 2D (on one layer)
  N2 = (N_x+1)*(N_y+1);
  // number of unknowns in 3D (one cube)
  N3 = (N_x+1)*(N_y+1)*(N_z+1);
  // number of unknowns in 4D (total grid)
  N4 = (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1);

  // array for concentration of species C
  concent_C_array = new double[N3];
  memset(concent_C_array, 0, N3*SizeOfDouble);
  // array for the rhs
  rhs = new double[N4];
  memset(rhs, 0, N4*SizeOfDouble);
  // initialization of the entries array
  memset(entries, 0, N_Entries*SizeOfDouble);

  // arrays for the velocity
  velo1 = new double[N3];
  memset(velo1, 0, N3*SizeOfDouble);
  velo2 = new double[N3];
  memset(velo2, 0, N3*SizeOfDouble);
  velo3 = new double[N3];
  memset(velo3, 0, N3*SizeOfDouble);

  // discretization of PBE with FDM
  // loop over all nodes of the FDM grid
  for ( i=0 ; i<N4 ; i++ )
  {
    // node is on the first cube (a_coord=0)
    // in this cube, the velocity vector will be filled
    // since this vector is independent of a
    if (i < N3)
    {
	/* old version 
      // u1
      velocity1->FindGradient(x_coord[i],y_coord[i],z_coord[i],values);
      velo1[i] = values[0];
      velocity1_array_val = velo1[i];

      //u2
      velocity2->FindGradient(x_coord[i],y_coord[i],z_coord[i],values);
      velo2[i] = values[0];
      velocity2_array_val = velo2[i];

      //u3
      velocity3->FindGradient(x_coord[i],y_coord[i],z_coord[i],values);
      velo3[i] = values[0];
      velocity3_array_val = velo3[i];

      // fill the array for the concentrations of C
      // since this concentration does not depend on a
      // c_C
      concent_C->FindGradient(x_coord[i],y_coord[i],z_coord[i],values);
      // nur zum Vergleichen der Ergebnisse -> Programmtest
      // values[0] = 1+i*1.0/N3;
      concent_C_array[i] = values[0];
      concent_C_array_val = values[0];
	*/
	// THIS HAS TO BE CHECKED !!!
	// find grid cell of 3d grid
	ii = i;
	//OutPut("ii " << ii << endl);
	// top of "first cube" = z_max
	if( ii>= N2*N_z )
	{
	    ii = ii-N2;
	}
	//OutPut("a ii " << ii << endl);
	// treat right d.o.f. = x_max seperately
	if ( ((ii+1)%(N_x+1)==0) )
	{
	    ii = ii-1;
	}
	//OutPut("b ii " << ii << endl);
	// treat upper d.o.f. =  y_max seperately
	if ( ii%N2 >= (N_y*(N_x+1)) )
	{
	    ii = ii-(N_x+1);
	}
	//OutPut("c ii " << ii << endl);
	// level in z-direction 
	gamma = (int)(ii/N2);
	ii -= gamma*N2;
	// level in y-direction 
	alpha = (int)(ii/(N_x+1)+1e-6 );
	// level in x-direction
	beta = ii%(N_x+1);
	//OutPut("ii " << ii << " ga " << gamma << " be " << beta << " ga " <<
	//       gamma << " : " <<  gamma*N_x*N_y+ alpha * N_x + beta << endl);
	no_of_3dcell = correspond_3dgrid[gamma*N_x*N_y+ alpha * N_x + beta];

	// find velocity
	velocity1->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				     x_coord[i],y_coord[i],z_coord[i],values);
	velo1[i] = values[0];
	velocity1_array_val = velo1[i];
	
	velocity2->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				     x_coord[i],y_coord[i],z_coord[i],values);
	velo2[i] = values[0];
	velocity2_array_val = velo2[i];
	
	velocity3->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				     x_coord[i],y_coord[i],z_coord[i],values);
	velo3[i] = values[0];
	velocity3_array_val = velo3[i];
	
	// fill the array for the concentrations of C
	// since this concentration does not depend on a
	// c_C
	concent_C->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				     x_coord[i],y_coord[i],z_coord[i],values);
	concent_C_array[i] = values[0];
	concent_C_array_val = values[0];
    }
    // node is not on the first cube
    else
    {
      // the corresponding node on the first cube
      ii = i - N3*((int)(i/N3));
      // compute the value of the velocity for the corresponding (x,y,z) coordinates
      velocity1_array_val = velo1[ii];
      velocity2_array_val = velo2[ii];
      velocity3_array_val = velo3[ii];
      // compute the value of the concentration of C for the corresponding (x,y,z) coordinates
      concent_C_array_val = concent_C_array[ii];
    }

    if ( TDatabase::ParamDB->BULK_GROWTH_RATE==2 )
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(a_coord[i]*d_p_max)));
    }
    else
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty);
    }

    // compute the coefficients corresponding to the 4d finite-difference-method
    // diagonal entry
    // diagonal entry
    for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
    {
	if (col_ptr[k]==i)
	{
	    entries[k] = 1;
	    indices[0] = k;
	    diag_index = k;
	}
	else
	{
	    if (col_ptr[k]==i-1)
	    {
		// left node
		indices[1] = k;
	    }
	    else
	    {
		if (col_ptr[k]==i+1)
		{
		    // right node
		    indices[2] = k;
		}
	    }
	}
    }

    // compute the term A_x
    if (velocity1_array_val >= 0)
    {
      // not on boundary x = x_min -> in the "x-sense" exists a left neighbour
      if ( i%(N_x+1)!=0 )
      {
	  val =   deltat * velocity1_array_val/fabs(x_coord[i]-x_coord[i-1]);
	  index = indices[1];
	  entries[index] -= val;
	  entries[diag_index] += val;
      }
    }
    else
    {
      // not on boundary x = x_max -> in the "x-sense" exists a right neighbour
      if ( (i+1)%(N_x+1)!=0 )
      {
	  val = deltat * velocity1_array_val/fabs(x_coord[i+1]-x_coord[i]);
	  index = indices[2];
	  entries[index] += val;
	  entries[diag_index] -= val;
        }
    }
  
    // compute the term A_y
    index = i%N2;

    if (velocity2_array_val >= 0)
    {
      // not on boundary y = y_min -> in the "y-sense" exists a left neighbour
      if ( index>N_x )
      {
        index1 = i-(N_x+1);
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          if ( col_ptr[k] == index1 )
          {
	      val = deltat * velocity2_array_val/fabs(y_coord[i]-y_coord[index1]);
	      entries[k] -= val;
	      entries[diag_index] += val;
	    break;
          }
        }
      }
    }
    else
    {
      index2 = (N_x+1)*N_y;
      // not on boundary y = y_max -> in the "y-sense" exists a right neighbour
      if ( index<index2 )
      {
        index1 = i+(N_x+1);
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          if (col_ptr[k] == index1)
          {
	      val = deltat * velocity2_array_val/fabs(y_coord[index1]-y_coord[i]);
	      entries[k] += val;
	      entries[diag_index] -= val;
	      break;
	  }
        }
      }
    }

    // compute the term A_z
    index = i%N3;
    if (velocity3_array_val >= 0)
    {
      // not on boundary z = z_min -> in the "z-sense" exists a left neighbour
      if ( index >= N2 )
      {
        index1 = i-N2;
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          if ( col_ptr[k] == index1 )
          {
	      val = deltat * velocity3_array_val/fabs(z_coord[i]-z_coord[index1]);
	      entries[k] -= val;
	      entries[diag_index] += val;
	      break;
	  }
	}
      }
    }
    else
    {
      index2 = N2*N_z;
      // not on boundary z = z_max -> in the "z-sense" exists a right neighbour
      if ( index < index2 )
      {
        index1 = i+N2;
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          if (col_ptr[k] == index1)
          {
	      val = deltat * velocity3_array_val/fabs(z_coord[index1]-z_coord[i]);
	      entries[k] += val;
	      entries[diag_index] -= val;
	      break;
	  }
	}
      }
    }

    // compute the term C
    if (G_c_C_val >= 0)
    {
      // not on boundary a = a_min -> in the "a-sense" exists a left neighbour
      if ( i >= N3 )
      {
        index1 = i-N3;
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          if (col_ptr[k] == index1)
          {
	      val = deltat * G_c_C_val/fabs(a_coord[i]-a_coord[index1]);
	      entries[k] -= val;
	      entries[diag_index] += val;
	      break;
          }
        }
      }
    }
    else
    {
      // not on boundary a = a_max -> in the "a-sense" exists a right neighbour
      if( i<N3*N_a )
      {
        index1 = i+N3;
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
	    if (col_ptr[k] == index1)
	    {
		val = deltat * G_c_C_val/fabs(a_coord[index1]-a_coord[i]);
		entries[k] += val;
		entries[diag_index] -= val;
		break;
	    }
        }
      }
    }

    // rhs
    rhs[i] = sol[i];
  }

  // set Dirichlet boundary conditions
  // if convection is positive at the bottom ( = cube with a = a_min )
  for ( i=0 ; i<N3 ; i++ )
  {
	// CHECKED AT 08/10/09
	// find grid cell of 3d grid
	ii = i;
	// top of "first cube" = z_max
	if( ii>= N2*N_z )
	{
	    ii = ii-N2;
	}
	// treat right d.o.f. = x_max seperately
	if ( ((ii+1)%(N_x+1)==0) )
	{
	    ii = ii-1;
	}
	// treat upper d.o.f. =  y_max seperately
	if ( ii%N2 >= (N_y*(N_x+1)) )
	{
	    ii = ii-(N_x+1);
	}
	// level in z-direction 
	gamma = (int)(ii/N2);
	ii -= gamma*N2;
	// level in y-direction 
	alpha = (int)(ii/(N_x+1)+1e-6 );
	// level in x-direction
	beta = ii%(N_x+1);
	no_of_3dcell = correspond_3dgrid[gamma*N_x*N_y+ alpha * N_x + beta];

      concent_C->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				   x_coord[i],y_coord[i],z_coord[i],C_val);
      if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
      {
        G_c_C_val = k_g*(c_C_infty*C_val[0] - c_C_infty_sat*exp(C_2/d_p_0));
      }
      else
      {
        G_c_C_val = k_g*(c_C_infty*C_val[0] - c_C_infty_sat);
      }
      // compute G*n, n=(0,0,0,-1);
      if (G_c_C_val*f_infty > 1e-10)
      {
	  // compute rate of nucleation
	  B_c_C = k_nuc*pow(c_C_infty*(C_val[0] - 1),5);
	  // truncate negative values
	  if (B_c_C < 0)
	      B_c_C = 0;
	  // set rhs, sol and matrix row 
	  rhs[i] = B_c_C/ (G_c_C_val*f_infty);
	  sol[i] = rhs[i];
	  ii = row_ptr[i];
	  jj = row_ptr[i+1];
	  // off diagonals
	  for ( iq = ii ; iq < jj ; iq++ )
	  {
	      // diagonal entry
	      if(col_ptr[iq]==i)
		  entries[iq] = 1.0;
	      else
		  entries[iq] = 0;
	  }
      }
      else
      {
	  if (G_c_C_val*f_infty >= 0)
	  {
	      rhs[i] = sol[i] = 0;
	      ii = row_ptr[i];
	      jj = row_ptr[i+1];
	      // off diagonals
	      for ( iq = ii ; iq < jj ; iq++ )
	      {
		  // diagonal entry
		  if(col_ptr[iq]==i)
		    entries[iq] = 1.0;
		  else
		      entries[iq] = 0;
	      }
	  }
      }
	  
      // set Dirichlet boundary conditions
      // if convection is positive at the top
      if (G_c_C_val < 0)
      {
	  k = i + N3*N_a;
	  // compute new particle size distribution
	  rhs[k] = 0;
	  sol[k] = rhs[k];
	  // set matrix rows
	  ii = row_ptr[k];
	  jj = row_ptr[k+1];
	  // off diagonals
	  for ( iq = ii ; iq < jj ; iq++ )
	  {
	      // diagonal entry
	      if(col_ptr[iq]==k)
		  entries[iq] = 1.0;
	      else
		  entries[iq] = 0;
	  }
      }
  }

  // inflow from the velocity field
  for (i=N3;i<N4;i++)
  {
      iq =  PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i], z_coord[i]);
      if (iq)
      {
        rhs[i] =  0.0;
        sol[i] = rhs[i];
        ii = row_ptr[i];
        jj = row_ptr[i+1];
        // off diagonals
        for ( iq = ii ; iq < jj ; iq++ )
        {
          // diagonal entry
          if(col_ptr[iq]==i)
            entries[iq] = 1.0;
          else
            entries[iq] = 0;
        }
      }
  }

  if ( sqrt(Ddot(N4,rhs,rhs)) > 0 )
  {
    if ( N4 < SC_LDS )
    {
      DirectSolver(mat, rhs, sol);
      OutPut("SolveLU MEMORY: " << setw(10) << GetMemory() << endl);
    }
    else
    {
	t3 = TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
	TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = 0;
	Solver(mat,rhs,sol);
 	TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = t3;
	// no output of solver data any longer
	TDatabase::ParamDB->SC_VERBOSE_AMG = 1;
    }
  }
  else
  {
    memset(sol,0,N4*SizeOfDouble);
  }

  // cut undershoots
  for ( i=0 ; i<N4 ; i++ )
  {
    if ( sol[i] < 0 )
    {
      sol[i] = 0;
    }
  }

  // free allocated memory
  delete rhs;
  delete velo1;
  delete velo2;
  delete velo3;
  delete concent_C_array;
}

/****************************************************************************************
 *                                                                                      *
 *                         Part III : FEM                                               *
 *                        ----------------                                              *
 *                                                                                      *
 ***************************************************************************************/


/****************************************************************************************
 *                                                                                      *
 * computes the value of the gradient for Q1                                            *
 * input:  *coeff, x, y, z, a                                                           *
 * output: *val                                                                         *
 *                                                                                      *
 ***************************************************************************************/

void Compute_Q1_Value_Gradient_4D(double *coeff, double x, double y, double z, double a, double *val)
{
  // precomputing of some products
  double xy, xz, xa, yz, ya, za;

  xy = x*y;
  xz = x*z;
  xa = x*a;
  yz = y*z;
  ya = y*a;
  za = z*a;

  // value
  val[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z  + coeff[4] *  a
           + coeff[5] * xy + coeff[6] * xz + coeff[7] * xa + coeff[8] * yz + coeff[9] * ya
           + coeff[10] * za + coeff[11] * xy * z + coeff[12] * xy * a + coeff[13] * xz * a
           + coeff[14] * yz * a + coeff[15] * xz * ya;

  // value of the x derivative
  val[1] = coeff[1] + coeff[5] * y + coeff[6] * z + coeff[7] * a + coeff[11] * yz
           + coeff[12] * ya + coeff[13] * za + coeff[15] * yz * a;

  // value of the y derivative
  val[2] = coeff[2] + coeff[5] * x + coeff[8] * z + coeff[9] * a + coeff[11] * xz
           + coeff[12] * xa + coeff[14] * za + coeff[15] * xz * a;

  // value of the z derivative
  val[3] = coeff[3] + coeff[6] * x + coeff[8] * y + coeff[10] * a + coeff[11] * xy
           + coeff[13] * xa + coeff[14] * ya + coeff[15] * xy * a;

  // value of the a derivative
  val[4] = coeff[4] + coeff[7] * x + coeff[9] * y + coeff[10] * z + coeff[12] * xy
           + coeff[13] * xz + coeff[14] * yz + coeff[15] * xy * z;
}


void Compute_Q1_Value_4D(double *coeff, double x, double y, double z, double a, double *val)
{
  // precomputing of some products
  double xy, xz, xa, yz, ya, za;

  xy = x*y;
  xz = x*z;
  xa = x*a;
  yz = y*z;
  ya = y*a;
  za = z*a;

  // value
  val[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z  + coeff[4] *  a
           + coeff[5] * xy + coeff[6] * xz + coeff[7] * xa + coeff[8] * yz + coeff[9] * ya
           + coeff[10] * za + coeff[11] * xy * z + coeff[12] * xy * a + coeff[13] * xz * a
           + coeff[14] * yz * a + coeff[15] * xz * ya;
}


void Compute_Q1_Gradient_4D(double *coeff, double x, double y, double z, double a, double *val)
{
  // precomputing of some products
  double xy, xz, xa, yz, ya, za;

  xy = x*y;
  xz = x*z;
  xa = x*a;
  yz = y*z;
  ya = y*a;
  za = z*a;

  // value
  val[0] = 0;

  // value of the x derivative
  val[1] = coeff[1] + coeff[5] * y + coeff[6] * z + coeff[7] * a + coeff[11] * yz
           + coeff[12] * ya + coeff[13] * za + coeff[15] * yz * a;

  // value of the y derivative
  val[2] = coeff[2] + coeff[5] * x + coeff[8] * z + coeff[9] * a + coeff[11] * xz
           + coeff[12] * xa + coeff[14] * za + coeff[15] * xz * a;

  // value of the z derivative
  val[3] = coeff[3] + coeff[6] * x + coeff[8] * y + coeff[10] * a + coeff[11] * xy
           + coeff[13] * xa + coeff[14] * ya + coeff[15] * xy * a;

  // value of the a derivative
  val[4] = coeff[4] + coeff[7] * x + coeff[9] * y + coeff[10] * z + coeff[12] * xy
           + coeff[13] * xz + coeff[14] * yz + coeff[15] * xy * z;
}


/****************************************************************************************
 *                                                                                      *
 * this routine computes the corresponding 3d grid                                      *
 * input:  N_x, N_y, N_z, x_coord, y_coord, z_coord, coll                               *
 * output: correspond_3dgrid                                                            *
 *                                                                                      *
 ***************************************************************************************/

void generate_correspond_3d_grid(int N_x, int N_y, int N_z,
                                 double *x_coord, double *y_coord, double *z_coord,
                                 TCollection *coll, int *correspond_3dgrid)
{
  int i, j, N_Vertices, N_Cells, alpha, beta, gamma;
  double  x_low, y_low, z_low, sx, sy, sz, eps = 1e-6;

  TBaseCell *cell;

  // number of mesh cells
  N_Cells = coll->GetN_Cells();

  for ( i=0 ; i<N_Cells ; i++ )
  {
    cell = coll->GetCell(i);
    // number of vertices
    N_Vertices = cell->GetN_Vertices();
    // compute left lower corner
    x_low = 1000;
    y_low = 1000;
    z_low = 1000;

    for ( j=0 ; j<N_Vertices ; j++ )
    {
      sx = cell->GetVertex(j)->GetX();
      if ( sx < x_low )
      {
        x_low = sx;
      }
      sy = cell->GetVertex(j)->GetY();
      if ( sy < y_low )
      {
        y_low = sy;
      }
      sz = cell->GetVertex(j)->GetZ();
      if ( sz < z_low )
      {
        z_low = sz;
      }
    }

    // the coordinates of the left lower vertex are computed
    // compare them with the coordinates of the mesh cells -> the "bottom cubes"
    for ( j=0 ; j<(N_x+1)*(N_y+1)*(N_z+1) ; j++ )
    {
      // check if all three coordinates are the same
      if ( ( fabs(x_coord[j]-x_low) < eps ) && ( fabs(y_coord[j]-y_low) < eps )
           && ( fabs(z_coord[j]-z_low) < eps ) )
      {
        gamma = (int)( j/((N_x+1)*(N_y+1))+1e-6 );
        alpha = (int)( (j-gamma*(N_x+1)*(N_y+1))/(N_x+1)+1e-6 );
        beta = (j-gamma*(N_x+1)*(N_y+1))%(N_x+1);
        correspond_3dgrid[alpha*(N_x)+beta+gamma*(N_x)*(N_y)] = i;
      }
    }
  }
}

void generate_correspond_3d_grid_urea_pipe(int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
TCollection *coll, int *correspond_3dgrid)
{
  int i, j, N_Cells;

  TBaseCell *cell;
  // number of mesh cells
  N_Cells = coll->GetN_Cells();

  for ( j=0 ; j<(N_x+1)*(N_y+1)*(N_z+1) ; j++ )
  {
    for ( i=0 ; i<N_Cells ; i++ )
    {
      cell = coll->GetCell(i);
      if (cell->PointInCell(x_coord[j],y_coord[j],z_coord[j]))
        break;
    }
    correspond_3dgrid[j] = i;
  }
}


/****************************************************************************************
 *                                                                                      *
 * this routine computes the entries for the column and the row pointer                 *
 * input:  Nodes, N_x, N_y, N_z, x_max, x_min, y_max, y_min, z_max, z_min, a_max, a_min *
 *         *x_coord, *y_coord, *z_coord, *a_coord                                       *
 * output: *row_ptr, *col_ptr                                                           *
 *                                                                                      *
 ***************************************************************************************/

void filling_row_and_col_ptr(int *N_Entries, int Nodes,
                             int N_x, int N_y, int N_z,
                             double x_max, double x_min, double y_max, double y_min,
                             double z_max, double z_min, double a_max, double a_min,
                             double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                             int *row_ptr, int *col_ptr)
{
  int i, j;
  int x_left, x_right, y_left, y_right, z_left, z_right, a_left, a_right;
  int N_Neigh, loc_index, index;
  int neigh[81];

  // N2 and N3 are defined to save multiplications during the computation of the loc_indices
  int N2 = (N_x+1)*(N_y+1);
  int N3 = N2 * (N_z+1);

  row_ptr[0] = 0;
  *N_Entries = 0;

  // loop over Nodes
  for ( i=0 ; i<Nodes ; i++ )
  {
    // find all neighbors -> number of entries
    // the if queries check the membership of a node with the plains/cubes
    x_left = x_right = y_left = y_right = z_left = z_right = a_left = a_right = 0;

    if ( fabs(x_coord[i]-x_min) < 1e-6 )
    {
      x_left++;
    }
    if ( fabs(x_coord[i]-x_max) < 1e-6 )
    {
      x_right++;
    }
    if ( fabs(y_coord[i]-y_min) < 1e-6 )
    {
      y_left++;
    }
    if ( fabs(y_coord[i]-y_max) < 1e-6 )
    {
      y_right++;
    }
    if ( fabs(z_coord[i]-z_min) < 1e-6 )
    {
      z_left++;
    }
    if ( fabs(z_coord[i]-z_max) < 1e-6 )
    {
      z_right++;
    }
    if ( fabs(a_coord[i]-a_min) < 1e-6 )
    {
      a_left++;
    }
    if ( fabs(a_coord[i]-a_max) < 1e-6 )
    {
      a_right++;
    }

    switch ( x_left + x_right + y_left + y_right + z_left + z_right + a_left + a_right )
    {
      case 0:
        N_Neigh = 81;      // interior node
        break;
      case 1:
        N_Neigh = 54;      // 1 dim off -> one level off  -> 27 neighbours ==> 81 - 27 = 54
        break;
      case 2:
        N_Neigh = 36;      // 2 dim off -> next level off -> 18 neighbours ==> 54 - 27 + 9 = 36
        break;
      case 3:
        N_Neigh = 24;      // 3 dim off -> next level off -> 12 neighbours ==> 36 - 27 + 9 + 6 = 24
        break;
      case 4:
        N_Neigh = 16;      // 4 dim off -> next level off ->  8 neighbours ==> 24 - 27 + 9 + 6 + 4 = 16
        break;
      default:
        OutPut("Error 1 in filling_row_and_col_ptr!"<<endl);
        exit(1);
    }
    // number of entries in the row i
    *N_Entries += N_Neigh;
    index = row_ptr[i];

    row_ptr[i+1] = index + N_Neigh;

    // diagonal entry
    col_ptr[index] = i;
    index++;

    // the array neigh is set to one
    for ( j=0 ; j<81 ; j++ )
    {
      neigh[j] = 1;
    }

    // i is not an inner point -> less neighbors -> corresponding neigh_entries will be set to zero
    if ( N_Neigh!=81 )
    {
      // on the x_left boundary, set local left neighbors to zero
      if (x_left)
      {
        neigh[0] = neigh[3] = neigh[6] = neigh[9] = neigh[12] = neigh[15]
        = neigh[18] = neigh[21] = neigh[24] = neigh[27] = neigh[30]
        = neigh[33] = neigh[36] = neigh[39] = neigh[42] = neigh[45]
        = neigh[48] = neigh[51] = neigh[54] = neigh[57] = neigh[60]
        = neigh[63] = neigh[66] = neigh[69] = neigh[72] = neigh[75]
        = neigh[78] = 0;
      }
      // on the x_right boundary, set local right neighbors to zero
      if (x_right)
      {
        neigh[2] = neigh[5] = neigh[8] = neigh[11] = neigh[14] = neigh[17]
        = neigh[20] = neigh[23] = neigh[26] = neigh[29] = neigh[32]
        = neigh[35] = neigh[38] = neigh[41] = neigh[44] = neigh[47]
        = neigh[50] = neigh[53] = neigh[56] = neigh[59] = neigh[62]
        = neigh[65] = neigh[68] = neigh[71] = neigh[74] = neigh[77]
        = neigh[80] = 0;
      }
      // on the y_left boundary, set local left neighbors to zero
      if (y_left)
      {
        neigh[0] = neigh[1] = neigh[2] = neigh[9] = neigh[10] = neigh[11]
        = neigh[18] = neigh[19] = neigh[20] = neigh[27] = neigh[28]
        = neigh[29] = neigh[36] = neigh[37] = neigh[38] = neigh[45]
        = neigh[46] = neigh[47] = neigh[54] = neigh[55] = neigh[56]
        = neigh[63] = neigh[64] = neigh[65] = neigh[72] = neigh[73]
        = neigh[74] = 0;
      }
      // on the y_right boundary, set local right neighbors to zero
      if (y_right)
      {
        neigh[6] = neigh[7] = neigh[8] = neigh[15] = neigh[16] = neigh[17]
        = neigh[24] = neigh[25] = neigh[26] = neigh[33] = neigh[34]
        = neigh[35] = neigh[42] = neigh[43] = neigh[44] = neigh[51]
        = neigh[52] = neigh[53] = neigh[60] = neigh[61] = neigh[62]
        = neigh[69] = neigh[70] = neigh[71] = neigh[78] = neigh[79]
        = neigh[80] = 0;
      }
      // on the z_left boundary, set local left neighbors to zero
      if (z_left)
      {
        neigh[0] = neigh[1] = neigh[2] = neigh[3] = neigh[4] = neigh[5]
        = neigh[6] = neigh[7] = neigh[8] = neigh[27] = neigh[28]
        = neigh[29] = neigh[30] = neigh[31] = neigh[32] = neigh[33]
        = neigh[34] = neigh[35] = neigh[54] = neigh[55] = neigh[56]
        = neigh[57] = neigh[58] = neigh[59] = neigh[60] = neigh[61]
        = neigh[62] = 0;
      }
      // on the z_right boundary, set local right neighbors to zero
      if (z_right)
      {
        neigh[18] = neigh[19] = neigh[20] = neigh[21] = neigh[22] = neigh[23]
        = neigh[24] = neigh[25] = neigh[26] = neigh[45] = neigh[46]
        = neigh[47] = neigh[48] = neigh[49] = neigh[50] = neigh[51]
        = neigh[52] = neigh[53] = neigh[72] = neigh[73] = neigh[74]
        = neigh[75] = neigh[76] = neigh[77] = neigh[78] = neigh[79]
        = neigh[80] = 0;
      }
      // on the a_left boundary, set local left neighbors to zero
      if (a_left)
      {
        neigh[0] = neigh[1] = neigh[2] = neigh[3] = neigh[4] = neigh[5]
        = neigh[6] = neigh[7] = neigh[8] = neigh[9] = neigh[10]
        = neigh[11] = neigh[12] = neigh[13] = neigh[14] = neigh[15]
        = neigh[16] = neigh[17] = neigh[18] = neigh[19] = neigh[20]
        = neigh[21] = neigh[22] = neigh[23] = neigh[24] = neigh[25]
        = neigh[26] = 0;
      }
      // on the a_right boundary, set local right neighbors to zero
      if (a_right)
      {
        neigh[54] = neigh[55] = neigh[56] = neigh[57] = neigh[58] = neigh[59]
        = neigh[60] = neigh[61] = neigh[62] = neigh[63] = neigh[64]
        = neigh[65] = neigh[66] = neigh[67] = neigh[68] = neigh[69]
        = neigh[70] = neigh[71] = neigh[72] = neigh[73] = neigh[74]
        = neigh[75] = neigh[76] = neigh[77] = neigh[78] = neigh[79]
        = neigh[80] = 0;
      }
    }
    // off-diagonal entries
    for ( j=0 ; j<81 ; j++ )
    {
      // the node itself
      if ( j==40 )
      {
        continue;
      }
      if ( neigh[j] )
      {
        switch(j)
        {
          case 0: loc_index = i-N3-N2-N_x-2;
          break;
          case 1: loc_index = i-N3-N2-N_x-1;
          break;
          case 2: loc_index = i-N3-N2-N_x;
          break;
          case 3: loc_index = i-N3-N2-1;
          break;
          case 4: loc_index = i-N3-N2;
          break;
          case 5: loc_index = i-N3-N2+1;
          break;
          case 6: loc_index = i-N3-N2+N_x;
          break;
          case 7: loc_index = i-N3-N2+N_x+1;
          break;
          case 8: loc_index = i-N3-N2+N_x+2;
          break;
          case 9: loc_index = i-N3-N_x-2;
          break;
          case 10: loc_index = i-N3-N_x-1;
          break;
          case 11: loc_index = i-N3-N_x;
          break;
          case 12: loc_index = i-N3-1;
          break;
          case 13: loc_index = i-N3;
          break;
          case 14: loc_index = i-N3+1;
          break;
          case 15: loc_index = i-N3+N_x;
          break;
          case 16: loc_index = i-N3+N_x+1;
          break;
          case 17: loc_index = i-N3+N_x+2;
          break;
          case 18: loc_index = i-N3+N2-N_x-2;
          break;
          case 19: loc_index = i-N3+N2-N_x-1;
          break;
          case 20: loc_index = i-N3+N2-N_x;
          break;
          case 21: loc_index = i-N3+N2-1;
          break;
          case 22: loc_index = i-N3+N2;
          break;
          case 23: loc_index = i-N3+N2+1;
          break;
          case 24: loc_index = i-N3+N2+N_x;
          break;
          case 25: loc_index = i-N3+N2+N_x+1;
          break;
          case 26: loc_index = i-N3+N2+N_x+2;
          break;
          case 27: loc_index = i-N2-N_x-2;
          break;
          case 28: loc_index = i-N2-N_x-1;
          break;
          case 29: loc_index = i-N2-N_x;
          break;
          case 30: loc_index = i-N2-1;
          break;
          case 31: loc_index = i-N2;
          break;
          case 32: loc_index = i-N2+1;
          break;
          case 33: loc_index = i-N2+N_x;
          break;
          case 34: loc_index = i-N2+N_x+1;
          break;
          case 35: loc_index = i-N2+N_x+2;
          break;
          case 36: loc_index = i-N_x-2;
          break;
          case 37: loc_index = i-N_x-1;
          break;
          case 38: loc_index = i-N_x;
          break;
          case 39: loc_index = i-1;
          break;
          case 41: loc_index = i+1;
          break;
          case 42: loc_index = i+N_x;
          break;
          case 43: loc_index = i+N_x+1;
          break;
          case 44: loc_index = i+N_x+2;
          break;
          case 45: loc_index = i+N2-N_x-2;
          break;
          case 46: loc_index = i+N2-N_x-1;
          break;
          case 47: loc_index = i+N2-N_x;
          break;
          case 48: loc_index = i+N2-1;
          break;
          case 49: loc_index = i+N2;
          break;
          case 50: loc_index = i+N2+1;
          break;
          case 51: loc_index = i+N2+N_x;
          break;
          case 52: loc_index = i+N2+N_x+1;
          break;
          case 53: loc_index = i+N2+N_x+2;
          break;
          case 54: loc_index = i+N3-N2-N_x-2;
          break;
          case 55: loc_index = i+N3-N2-N_x-1;
          break;
          case 56: loc_index = i+N3-N2-N_x;
          break;
          case 57: loc_index = i+N3-N2-1;
          break;
          case 58: loc_index = i+N3-N2;
          break;
          case 59: loc_index = i+N3-N2+1;
          break;
          case 60: loc_index = i+N3-N2+N_x;
          break;
          case 61: loc_index = i+N3-N2+N_x+1;
          break;
          case 62: loc_index = i+N3-N2+N_x+2;
          break;
          case 63: loc_index = i+N3-N_x-2;
          break;
          case 64: loc_index = i+N3-N_x-1;
          break;
          case 65: loc_index = i+N3-N_x;
          break;
          case 66: loc_index = i+N3-1;
          break;
          case 67: loc_index = i+N3;
          break;
          case 68: loc_index = i+N3+1;
          break;
          case 69: loc_index = i+N3+N_x;
          break;
          case 70: loc_index = i+N3+N_x+1;
          break;
          case 71: loc_index = i+N3+N_x+2;
          break;
          case 72: loc_index = i+N3+N2-N_x-2;
          break;
          case 73: loc_index = i+N3+N2-N_x-1;
          break;
          case 74: loc_index = i+N3+N2-N_x;
          break;
          case 75: loc_index = i+N3+N2-1;
          break;
          case 76: loc_index = i+N3+N2;
          break;
          case 77: loc_index = i+N3+N2+1;
          break;
          case 78: loc_index = i+N3+N2+N_x;
          break;
          case 79: loc_index = i+N3+N2+N_x+1;
          break;
          case 80: loc_index = i+N3+N2+N_x+2;
          break;
        }
        col_ptr[index] = loc_index;
        index++;
      }
    }
  }
}


/****************************************************************************************
 *                                                                                      *
 * computation of the size of a mesh cell in convection direction                       *
 * approximation formula by Tezduyar and Park, CMAME 59, 307 - 325, 1986                *
 *                                                                                      *
 ***************************************************************************************/

double Bulk_mesh_size_in_convection_direction(double hK,
                                              double b1, double b2, double b3, double b4,
                                              double *x, double *y, double *z, double *a_4d)
{
  int i;

  double den, val, norm_b;
  double sx, sy, sz, sa;
  double a[256], b[256];

  // only for hexahedra
  sx = sy = sz = sa = 0;

  for ( i=0 ; i<16 ; i++ )
  {
    sx += x[i];
    sy += y[i];
    sz += z[i];
    sa += a_4d[i];
  }
  // bary centre
  sx /= 16;
  sy /= 16;
  sz /= 16;
  sa /= 16;

  // initialize rhs
  memset(b,0,256*SizeOfDouble);

  // set matrices for computation of the coefficients
  // of the bilinear function
  for ( i=0 ; i<16 ; i++ )
  {
      a[16*i]    = 1;
      a[16*i+1]  = x[i];
      a[16*i+2]  = y[i];
      a[16*i+3]  = z[i];
      a[16*i+4]  = a_4d[i];
      a[16*i+5]  = x[i]*y[i];
      a[16*i+6]  = x[i]*z[i];
      a[16*i+7]  = x[i]*a_4d[i];
      a[16*i+8]  = y[i]*z[i];
      a[16*i+9]  = y[i]*a_4d[i];
      a[16*i+10] = z[i]*a_4d[i];
      a[16*i+11] = x[i]*y[i]*z[i];
      a[16*i+12] = x[i]*y[i]*a_4d[i];
      a[16*i+13] = x[i]*z[i]*a_4d[i];
      a[16*i+14] = y[i]*z[i]*a_4d[i];
      a[16*i+15] = x[i]*y[i]*z[i]*a_4d[i];
      b[17*i]    = 1;
  }

  // solve system for the coefficients of the bilinear function
  SolveMultipleSystems(a,b,16,16,16,16);

  // compute numerator
  norm_b = sqrt(b1*b1 + b2*b2 + b3*b3 + b4*b4);

  // compute denominator
  den = 0;
  for ( i=0 ; i<16 ; i++ )
  {
    // value of gradient basis fct. in bary centre
    val  = b1*(b[16*i+1] + b[16*i+5]*sy + b[16*i+6]*sz + b[16*i+7]*sa + b[16*i+11]*sy*sz
               + b[16*i+12]*sy*sa + b[16*i+13]*sz*sa + b[16*i+15]*sy*sz*sa);
    val += b2*(b[16*i+2] + b[16*i+5]*sx + b[16*i+8]*sz + b[16*i+9]*sa + b[16*i+11]*sx*sz
               + b[16*i+12]*sx*sa + b[16*i+14]*sz*sa + b[16*i+15]*sx*sz*sa);
    val += b3*(b[16*i+3] + b[16*i+6]*sx + b[16*i+8]*sy + b[16*i+10]*sa + b[16*i+11]*sx*sy
               + b[16*i+13]*sx*sa + b[16*i+14]*sy*sa + b[16*i+15]*sx*sy*sa);
    val += b4*(b[16*i+4] + b[16*i+7]*sx + b[16*i+9]*sy + b[16*i+10]*sz + b[16*i+12]*sx*sy
               + b[16*i+13]*sx*sz + b[16*i+14]*sy*sz + b[16*i+15]*sx*sy*sz);

    den += fabs(val);
  }

  // return the mesh size in convection direction
  if ( den < 1e-6 )
  {
    return(hK);
  }
  else
  {
    return(2*norm_b/den);
  }
}

/****************************************************************************************
 *                                                                                      *
 *  assembling of matrix with Q_1 finite elements, supg,                                *
 *  solution of arising linear system                                                   *
 *  OutPut: sol (contains the solution)                                                 *
 *                                                                                      *
 ***************************************************************************************/

void Build_4D_FEM_Matrix_Q1(TCollection *coll,
                            TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            double *sol, double *oldsol,
                            double *lump_mass_PSD, double *matrix_D_Entries_PSD,
                            int *correspond_3dgrid,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            TSquareMatrix3D *mat, TSquareMatrix3D *matM)

{
    OutPut("Build_4D_FEM_Matrix_Q1 has been removed !!!"<< endl);
    OutPut("The latest version is in ~MooNMD/src/FE/OLD/Bulk3d4d.081007.C !!!"<< endl);
    exit(4711);
}

/****************************************************************************************
 *                                                                                      *
 *  computes all dof which are Dirichlet dof                                            *
 *  these has to be treated as Neumann dof in FEM--FCT                                  *
 *                                                                                      *
 ***************************************************************************************/
void Compute_Neum_To_Diri_FEM_FCT(int N_x, int N_y, int N_z, int N_a,
				  double *x_coord, double *y_coord, 
				  double *z_coord, double *a_coord,
				  int &N_neum_to_diri, 
				  int* &neum_to_diri,
				  double* &neum_to_diri_x,
				  double* &neum_to_diri_y,
				  double* &neum_to_diri_z,
				  double* &neum_to_diri_a)
{
    int i, Nodes, range, count = 0, val;
    double yq;

    Nodes = (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1);
    
    // inflow from the velocity field
    for (i=(N_x+1)*(N_y+1)*(N_z+1);i<Nodes;i++)
    {
	val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i],z_coord[i]);
	if (val)
	    count++;
    }

    // total number of Dirichlet nodes
    count += (N_x+1)*(N_y+1)*(N_z+1);
    neum_to_diri = new int[count];
    neum_to_diri_x = new double[4*count];
    neum_to_diri_y = neum_to_diri_x + count;
    neum_to_diri_z = neum_to_diri_y + count;    
    neum_to_diri_a = neum_to_diri_z + count;    
    N_neum_to_diri = count;    
    
    // 3D flow domain
    for ( i=0 ; i< (N_x+1)*(N_y+1)*(N_z+1) ; i++ )
    {
	neum_to_diri[i] = i;
	neum_to_diri_x[i] = x_coord[i];
	neum_to_diri_y[i] = y_coord[i];
	neum_to_diri_z[i] = z_coord[i];
	neum_to_diri_a[i] = a_coord[i];
    }
    // velocity inlets
    count = (N_x+1)*(N_y+1)*(N_z+1);
    for (i=(N_x+1)*(N_y+1)*(N_z+1);i<Nodes;i++)
    {
	val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i], z_coord[i]);
	if (val)
	{
	    neum_to_diri[count] = i;
	    neum_to_diri_x[count] = x_coord[i];
	    neum_to_diri_y[count] = y_coord[i];
	    neum_to_diri_z[count] = z_coord[i];
	    neum_to_diri_a[count] = a_coord[i];
	    count++;
	}
    }
}

/****************************************************************************************
 *                                                                                      *
 *  assembling of mass matrix for Q_1 finite elements                                   *
 *  array with indices for assembling is filled too                                     *
 *                                                                                      *
 ***************************************************************************************/
void Build_4D_FEM_FCT_MassMatrix_Q1(TCollection *coll,
                                    int N_x, int N_y, int N_z, int N_a,
                                    double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                                    int* &index_test_ansatz, 
				    TSquareMatrix3D *matM,
                                    double *lump_mass_PSD)
{
  int locdof[16], *col_ptr, *row_ptr;
  int i, j, k, iq, ii, jj, z, z1, N_Entries, Nodes;
  int test_index, ansatz_index, index, index1, z_iq, z_ii, z_local;
  int diag_index, found, N_cells, range_y, range_z;
  int quad_points = 16;

  // N2 and N3 are defined to save multiplications during the computation of loc_dof
  int N2 = (N_x+1)*(N_y+1);
  int N3 = N2 * (N_z+1);

  double a[256], b[256], C_val[4], u_val[4], val_test[5], val_ansatz[5];
  double x_coord_loc[16], y_coord_loc[16], z_coord_loc[16], a_coord_loc[16];
  double *entriesM;
  double x_max, y_max, z_max, a_max, a_min, smag;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, aq, val, weight_det;
  double t1, t2;

  double weight[16]={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  double qx[16]=
  {
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489
  };
  double qy[16]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489
  };
  double qz[16]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489
  };
  double qa[16] =
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489
  };

  t1 = GetTime();

  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();

  entriesM = matM->GetEntries();
  N_Entries = matM->GetN_Entries();
  col_ptr = matM->GetKCol();
  row_ptr = matM->GetRowPtr();

  Nodes = N3 * (N_a+1);
  x_max = x_coord[Nodes-1];
  y_max = y_coord[Nodes-1];
  z_max = z_coord[Nodes-1];
  a_max = a_coord[Nodes-1];
  a_min = a_coord[0];

  // allocate array for storing indices for assembling
  ii = N_cells * N_a * 256;
  index_test_ansatz = new int[ii];
  OutPut("index_test_ansatz " << ii << endl);

  memset(entriesM,0, N_Entries*SizeOfDouble);

  z = 0;
  t2 = GetTime();
  OutPut("time fctmass (1) " << t2-t1 << endl);

  // loop over Nodes
  for ( i=0 ; i<Nodes ; i++ )
  {
    // assign nodes with mesh cell if it is not at x_max, y_max, z_max or a_max
    if ( (fabs(x_coord[i] - x_max)< 1e-6) || (fabs(y_coord[i] - y_max)< 1e-6)
         || (fabs(z_coord[i] - z_max)< 1e-6) || (fabs(a_coord[i] - a_max)< 1e-6) )
    {
      continue;
    }

    // consider the mesh cell with node i on the left lower corner
    // compute dof on these mesh cell
    locdof[0] = i;
    locdof[1] = i+1;
    locdof[2] = locdof[1] + N_x+1;
    locdof[3] = locdof[0] + N_x+1;
    locdof[4] = locdof[0] + N2;
    locdof[5] = locdof[4] + 1;
    locdof[6] = locdof[5] + N_x+1;
    locdof[7] = locdof[4] + N_x+1;
    locdof[8] = locdof[0] + N3;
    locdof[9] = locdof[1] + N3;
    locdof[10] = locdof[2] + N3;
    locdof[11] = locdof[3] + N3;
    locdof[12] = locdof[4] + N3;
    locdof[13] = locdof[5] + N3;
    locdof[14] = locdof[6] + N3;
    locdof[15] = locdof[7] + N3;

    // volume of the 4d "hexahedron"
    area = (x_coord[i+1]-x_coord[i])*(y_coord[i+N_x+1]-y_coord[i])*(z_coord[i+N2]-z_coord[i])*(a_coord[i+N3]-a_coord[i]);
    detJK = area/16.0;

    // compute basis functions
    // first: filling of the *_coord_loc vectors
    for ( j=0 ; j<16 ; j++ )
    {
      index = locdof[j];
      x_coord_loc[j] = x_coord[index];
      y_coord_loc[j] = y_coord[index];
      z_coord_loc[j] = z_coord[index];
      a_coord_loc[j] = a_coord[index];
    }
    // second : set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<16 ; j++ )
    {
      a[16*j]    = 1;
      a[16*j+1]  = x_coord_loc[j];
      a[16*j+2]  = y_coord_loc[j];
      a[16*j+3]  = z_coord_loc[j];
      a[16*j+4]  = a_coord_loc[j];
      a[16*j+5]  = x_coord_loc[j]*y_coord_loc[j];
      a[16*j+6]  = x_coord_loc[j]*z_coord_loc[j];
      a[16*j+7]  = x_coord_loc[j]*a_coord_loc[j];
      a[16*j+8]  = y_coord_loc[j]*z_coord_loc[j];
      a[16*j+9]  = y_coord_loc[j]*a_coord_loc[j];
      a[16*j+10] = z_coord_loc[j]*a_coord_loc[j];
      a[16*j+11] = x_coord_loc[j]*y_coord_loc[j]*z_coord_loc[j];
      a[16*j+12] = x_coord_loc[j]*y_coord_loc[j]*a_coord_loc[j];
      a[16*j+13] = x_coord_loc[j]*z_coord_loc[j]*a_coord_loc[j];
      a[16*j+14] = y_coord_loc[j]*z_coord_loc[j]*a_coord_loc[j];
      a[16*j+15] = x_coord_loc[j]*y_coord_loc[j]*z_coord_loc[j]*a_coord_loc[j];
    }

    // initialize rhs
    memset(b,0,256*SizeOfDouble);
    for ( j=0 ; j<16 ; j++ )
    {
       b[17*j] = 1;
    }

    // solve system for the coefficients of the bilinear function
    // solution is stored in b, row-wise
    SolveMultipleSystemsLapack(a,b,16,16,16,16);

    // assemble matrix entries
    // first index for array of indices
    z_iq = z;
    // loop over the quadrature points
    for ( iq = 0 ; iq < quad_points ; iq++ )
    {
      // quadrature points -> ONLY FOR 4D-PARALLELEPIPED !!!
      index  = locdof[0];
      index1 = locdof[1];
      xq = x_coord[index] + ( x_coord[index1] - x_coord[index])*(1 + qx[iq])/2;
      index1 = locdof[3];
      yq = y_coord[index] + ( y_coord[index1] - y_coord[index])*(1 + qy[iq])/2;
      index1 = locdof[4];
      zq = z_coord[index] + ( z_coord[index1] - z_coord[index])*(1 + qz[iq])/2;
      index1 = locdof[8];
      aq = a_coord[index] + ( a_coord[index1] - a_coord[index])*(1 + qa[iq])/2;
      weight_det = detJK * weight[iq];

      // loop for test function
      // ii -- test function
      z_ii = 0;
      for ( ii=0 ; ii<16 ; ii++ )
      {
        test_index = locdof[ii];

        // values for test function
        Compute_Q1_Value_4D(b+16*ii, xq, yq, zq, aq, val_test);
        val_test[0] *= weight_det;

        // loop for ansatz functions
        // jj -- ansatz function
        for ( jj=0 ; jj<16 ; jj++ )
        {
          ansatz_index = locdof[jj];

          // compute global index
	  // only for first quad point loop over the row_ptr
	  // index depends only on ii and jj but not on iq
 	  //index = test_index;
          // values for ansatz function
	  if (iq==0)
	  {
	      for ( k=row_ptr[test_index] ; k<row_ptr[test_index+1] ; k++ )
	      {
		  if ( col_ptr[k] == ansatz_index )
		  {
		      index = k;
		      index_test_ansatz[z] = index;
		      z++;
		      break;
		  }
	      }
	  }
	  else
	  {
	      // other quad points
	      z_local = z_iq + z_ii;
	      index = index_test_ansatz[z_local];
	      z_ii++;	      
	  }
	  //index = test_index;
          // values for ansatz function
	  if (test_index!=ansatz_index)//||(1))
	  {	  
	      Compute_Q1_Value_Gradient_4D(b+16*jj, xq, yq, zq, aq, val_ansatz);
	      
	      // mass term
	      val = val_ansatz[0]*val_test[0];
	  }
	  else
	  {
	      // mass term
	      val = val_test[0]*val_test[0]/weight_det; 
	  }

	  // add to M(locdof[ii], locdof[jj])
	  entriesM[index] += val;
        }
      }
    }    // end quad points
  }      // end i

  OutPut("z " << z << endl);
  LumpMassMatrixToVector((TSquareMatrix3D*) matM, lump_mass_PSD);
  t2 = GetTime();
  OutPut("time fctmass (2) " << t2-t1 << endl);
}

void Build_4D_FEM_FCT_Matrix_Q1(TCollection *coll,
				TFEFunction3D *velocity1, TFEFunction3D *velocity2, 
				TFEFunction3D *velocity3,
				TFEFunction3D *concent_C,
				double *sol, double *oldsol,
				double *lump_mass_PSD, double *matrix_D_Entries_PSD,
				int *correspond_3dgrid,
				int N_x, int N_y, int N_z, int N_a,
				double *x_coord, double *y_coord, double *z_coord, double *a_coord,
				TSquareMatrix3D *mat,
				TSquareMatrix3D *matM_cons,
				TSquareMatrix3D *matM,
				int *index_test_ansatz,
				int N_neum_to_diri, 
				int *neum_to_diri,
				double *neum_to_diri_x,
				double *neum_to_diri_y,
				double *neum_to_diri_z,
				double *neum_to_diri_a)
{
  int locdof[16], *col_ptr, *row_ptr;
  int i, j, k, iq, ii, jj, z, z1, N_Entries, Nodes;
  int test_index, ansatz_index, index, index1, alpha, beta, gamma, no_of_3dcell;
  int diag_index, found, N_cells, range_y, range_z, topdiri = 0;
  int quad_points = 16, count = 0, count_iq, count_ii, count_local;
  int SC_LDS =  TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  int *test_cells;

  // N2 and N3 are defined to save multiplications during the computation of loc_dof
  int N2 = (N_x+1)*(N_y+1);
  int N3 = N2 * (N_z+1);

  double a[256], b[256], C_val[4], u_val[6], val_test[5], val_ansatz[5], val_sol[5];
  double x_coord_loc[16], y_coord_loc[16], z_coord_loc[16], a_coord_loc[16];
  double sol_loc[16], b_sol[16];
  double *entries, *entriesM, *oldrhs_fem_fct0, *rhs, *RhsArray, *tilde_u, *u1, *u2, *u3, *G;
  double *bdr_val, *entriesM_cons;
  double x_max, y_max, z_max, a_max, a_min, smag;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, aq, val, weight_det;
  double B_c_C, maxsol, norm_b, al, react;
  double time = TDatabase::TimeDB->CURRENTTIME, t1, t2, t3;

  double weight[16]={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  double qx[16]=
  {
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489
  };
  double qy[16]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489
  };
  double qz[16]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489
  };
  double qa[16] =
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489,
     0.5773502691896257645091489,  0.5773502691896257645091489
  };

  t1 = GetTime();

  //model constants
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double C_2 = TDatabase::ParamDB->BULK_C_2;
  double d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double k_nuc = TDatabase::ParamDB->BULK_k_nuc;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;
  double factor_G, G_c_C_val, c_C_infty_sat0;

  Nodes = N3 * (N_a+1);
  
  // no nucleation and no particles present -> solution is zero
  /*if ((Ddot(Nodes,oldsol,oldsol)==0) && (TDatabase::ParamDB->P4<1))
  {
      memcpy(sol, oldsol, Nodes * SizeOfDouble);
      OutPut("set psd to zero"<<endl);
      return;
      }*/

  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();
  // initialize test_cells
  test_cells = new int[N_cells];
  memset(test_cells, 0, N_cells*SizeOfInt);

  // compute model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);

  entries = mat->GetEntries();
  N_Entries = mat->GetN_Entries();
  memset(entries,0, N_Entries*SizeOfDouble);
  col_ptr = mat->GetKCol();
  row_ptr = mat->GetRowPtr();

  entriesM = matM->GetEntries();
  entriesM_cons = matM_cons->GetEntries();

  // copy entries of mass matrix
  memcpy(entriesM, entriesM_cons, N_Entries*SizeOfDouble);

  x_max = x_coord[Nodes-1];
  y_max = y_coord[Nodes-1];
  z_max = z_coord[Nodes-1];
  a_max = a_coord[Nodes-1];
  a_min = a_coord[0];

  // initialization of u1, u2, u3, G in the vector u1
  // N_x*N_y*N_z - number of mesh cells for flow domain
  // need in each mesh cell 16 values (16 quad points) for 4 functions
  u1 = new double[N_x*N_y*N_z*64];
  memset(u1,0,(N_x*N_y*N_z*64)*SizeOfDouble);
  u2 = u1 + N_x*N_y*N_z*16;
  u3 = u2 + N_x*N_y*N_z*16;
  G  = u3 + N_x*N_y*N_z*16;

  // initialization of rhs and the vectors for FEM_FCT_ForConvDiff in the vector rhs
  rhs = new double[Nodes*4];
  memset(rhs,0,Nodes*4*SizeOfDouble);
  tilde_u  = rhs +  Nodes;
  RhsArray = tilde_u + Nodes;
  oldrhs_fem_fct0 = RhsArray + Nodes;

  bdr_val = new double[N_neum_to_diri];
  memset(bdr_val,0,N_neum_to_diri*SizeOfDouble);

  z = z1 = 0;
  c_C_infty_sat0 = c_C_infty_sat/c_C_infty;
  t2 = GetTime();
  OutPut("time bulkfct (1) " << t2-t1 << endl);

  // loop over Nodes
  for ( i=0 ; i<Nodes ; i++ )
  {
    // assign nodes with mesh cell if it is not at x_max, y_max, z_max or a_max
    if ( (fabs(x_coord[i] - x_max)< 1e-6) || (fabs(y_coord[i] - y_max)< 1e-6)
         || (fabs(z_coord[i] - z_max)< 1e-6) || (fabs(a_coord[i] - a_max)< 1e-6) )
    {
      continue;
    }

    // consider the mesh cell with node i on the left lower corner
    // compute dof on these mesh cell
    locdof[0] = i;
    locdof[1] = i+1;
    locdof[2] = locdof[1] + N_x+1;
    locdof[3] = locdof[0] + N_x+1;
    locdof[4] = locdof[0] + N2;
    locdof[5] = locdof[4] + 1;
    locdof[6] = locdof[5] + N_x+1;
    locdof[7] = locdof[4] + N_x+1;
    locdof[8] = locdof[0] + N3;
    locdof[9] = locdof[1] + N3;
    locdof[10] = locdof[2] + N3;
    locdof[11] = locdof[3] + N3;
    locdof[12] = locdof[4] + N3;
    locdof[13] = locdof[5] + N3;
    locdof[14] = locdof[6] + N3;
    locdof[15] = locdof[7] + N3;

    // volume of the 4d "hexahedron"
    area = (x_coord[i+1]-x_coord[i])*(y_coord[i+N_x+1]-y_coord[i])
	*(z_coord[i+N2]-z_coord[i])*(a_coord[i+N3]-a_coord[i]); 
    detJK = area/16.0;

    // compute basis functions
    // first: filling of the *_coord_loc vectors
    for ( j=0 ; j<16 ; j++ )
    {
      index = locdof[j];
      x_coord_loc[j] = x_coord[index];
      y_coord_loc[j] = y_coord[index];
      z_coord_loc[j] = z_coord[index];
      a_coord_loc[j] = a_coord[index];
      //sol_loc[j] = sol[index];
    }
    // second : set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<16 ; j++ )
    {
      a[16*j]    = 1;
      a[16*j+1]  = x_coord_loc[j];
      a[16*j+2]  = y_coord_loc[j];
      a[16*j+3]  = z_coord_loc[j];
      a[16*j+4]  = a_coord_loc[j];
      a[16*j+5]  = x_coord_loc[j]*y_coord_loc[j];
      a[16*j+6]  = x_coord_loc[j]*z_coord_loc[j];
      a[16*j+7]  = x_coord_loc[j]*a_coord_loc[j];
      a[16*j+8]  = y_coord_loc[j]*z_coord_loc[j];
      a[16*j+9]  = y_coord_loc[j]*a_coord_loc[j];
      a[16*j+10] = z_coord_loc[j]*a_coord_loc[j];
      a[16*j+11] = x_coord_loc[j]*y_coord_loc[j]*z_coord_loc[j];
      a[16*j+12] = x_coord_loc[j]*y_coord_loc[j]*a_coord_loc[j];
      a[16*j+13] = x_coord_loc[j]*z_coord_loc[j]*a_coord_loc[j];
      a[16*j+14] = y_coord_loc[j]*z_coord_loc[j]*a_coord_loc[j];
      a[16*j+15] = x_coord_loc[j]*y_coord_loc[j]*z_coord_loc[j]*a_coord_loc[j];
    }
    // initialize rhs
    memset(b,0,256*SizeOfDouble);
    for ( j=0 ; j<16 ; j++ )
    {
      b[17*j] = 1;
    }
    // solve system for the coefficients of the bilinear function
    // solution is stored in b, row-wise
    SolveMultipleSystemsLapack(a,b,16,16,16,16);

    // find corresponding cell in 3d grid
    if ( i < N3 )
    {
	// CHECKED AT 08/10/09
	// find grid cell of 3d grid
	ii = i;
	// top of "first cube" = z_max
	if( ii>= N2*N_z )
	{
	    ii = ii-N2;
	}
	// treat right d.o.f. = x_max seperately
	if ( ((ii+1)%(N_x+1)==0) )
	{
	    ii = ii-1;
	}
	// treat upper d.o.f. =  y_max seperately
	if ( ii%N2 >= (N_y*(N_x+1)) )
	{
	    ii = ii-(N_x+1);
	}
	// level in z-direction 
	gamma = (int)(ii/N2);
	ii -= gamma*N2;
	// level in y-direction 
	alpha = (int)(ii/(N_x+1)+1e-6 );
	// level in x-direction
	beta = ii%(N_x+1);
	no_of_3dcell = correspond_3dgrid[gamma*N_x*N_y+ alpha * N_x + beta];

      if (no_of_3dcell >= N_cells)
      {
	  OutPut("number of cell " << no_of_3dcell << " too large " << N_cells << endl);
	  exit(4711);
      }
      if (test_cells[no_of_3dcell])
      {
	  OutPut("cell " << no_of_3dcell << " treated a second time" << endl);
	  exit(4711);
      }
      test_cells[no_of_3dcell]++;
    }
    // assemble matrix entries
    // first index for array of indices
    count_iq = count;
    // loop over the quadrature points
    for ( iq = 0 ; iq < quad_points ; iq++ )
    {
      // quadrature points -> ONLY FOR 4D-PARALLELEPIPED !!!
      index = locdof[0];
      index1 = locdof[1];
      xq = x_coord[index] + ( x_coord[index1] - x_coord[index])*(1 + qx[iq])/2;
      index1 = locdof[3];
      yq = y_coord[index] + ( y_coord[index1] - y_coord[index])*(1 + qy[iq])/2;
      index1 = locdof[4];
      zq = z_coord[index] + ( z_coord[index1] - z_coord[index])*(1 + qz[iq])/2;
      index1 = locdof[8];
      aq = a_coord[index] + ( a_coord[index1] - a_coord[index])*(1 + qa[iq])/2;
      weight_det = detJK * weight[iq];

      // "bottom cube", fill vectors u1, u2, u3, G
      // pointer z is increased at the end of the loop
      if ( i < N3 )
      {
        // compute parameters in quadrature points
        // u1 = u_val[0], u2 = u_val[1], u3 = u_val[2]
        velocity1->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,xq,yq,zq,u_val);
        velocity2->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,xq,yq,zq,u_val+1);
        velocity3->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,xq,yq,zq,u_val+2);

        // concentration of species C, C = C_val[0]
        concent_C->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,xq,yq,zq,C_val);

        // G(c_C)
        G_c_C_val = factor_G*(C_val[0]-c_C_infty_sat0);

        u1[z1] = u_val[0];
        u2[z1] = u_val[1];
        u3[z1] = u_val[2];
        G[z1] = G_c_C_val;
        z1++;
      }

      count_ii = 0;
      // loop for test function
      // ii -- test function
      for ( ii=0 ; ii<16 ; ii++ )
      {
        // values for test function
        Compute_Q1_Value_4D(b+16*ii, xq, yq, zq, aq, val_test);
	val_test[0] *= weight_det;

        // loop for ansatz functions
        // jj -- ansatz function
        for ( jj=0 ; jj<16 ; jj++ )
        {
          // ansatz_index = locdof[jj];
	    // indices are different
	    // only first quad point
	    if (iq==0)  
	    {
		index = index_test_ansatz[count];
		count++;
	    }
	    else
	    {
	      // other quad points
	      count_local = count_iq + count_ii;
	      index = index_test_ansatz[count_local];
	      count_ii++;
	    }
          // index = test_index;
          // values for ansatz function
          Compute_Q1_Gradient_4D(b+16*jj, xq, yq, zq, aq, val_ansatz);
          entries[index] += (u1[z]*val_ansatz[1]+u2[z]*val_ansatz[2]+u3[z]*val_ansatz[3]
			     +G[z]*val_ansatz[4])*val_test[0];
        }
      }
      z++;
    }   // end quad points

    if ( z == N_x*N_y*N_z*16 )
    {
      z = 0;
    }
  }
  // end i

  t2 = GetTime();
  OutPut("time bulkfct (2) " << t2-t1 << endl);
  // FEM--FCT
  // set Dirichlet boundary conditions on all (possible) inflow boundaries
  // on the "first cube"
  for ( i=0 ; i<N3 ; i++ )
  {
	// CHECKED AT 08/10/09
	// find grid cell of 3d grid
	ii = i;
	// top of "first cube" = z_max
	if( ii>= N2*N_z )
	{
	    ii = ii-N2;
	}
	// treat right d.o.f. = x_max seperately
	if ( ((ii+1)%(N_x+1)==0) )
	{
	    ii = ii-1;
	}
	// treat upper d.o.f. =  y_max seperately
	if ( ii%N2 >= (N_y*(N_x+1)) )
	{
	    ii = ii-(N_x+1);
	}
	// level in z-direction 
	gamma = (int)(ii/N2);
	ii -= gamma*N2;
	// level in y-direction 
	alpha = (int)(ii/(N_x+1)+1e-6 );
	// level in x-direction
	beta = ii%(N_x+1);
	no_of_3dcell = correspond_3dgrid[gamma*N_x*N_y+ alpha * N_x + beta];
      
	concent_C->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
				     x_coord[i],y_coord[i],z_coord[i],C_val);

	G_c_C_val = k_g*(c_C_infty*C_val[0]-c_C_infty_sat);

	if ( G_c_C_val*f_infty > 1e-10 )
	{
	    // compute rate of nucleation
	    B_c_C = k_nuc*pow(c_C_infty*(C_val[0] - 1),5);
	    // truncate negative values
	    if (B_c_C < 0)
	    {
		B_c_C = 0;
	    }
	    // compute new particle size distribution
	    bdr_val[i] = B_c_C/(G_c_C_val*f_infty);
	}
	else
	{
	    // top Dirichlet condition, bdr_val[i] is already set to be zero
	    if (G_c_C_val*f_infty < 0)
	    {
		j = i + N_a * N3;
		neum_to_diri[i] = j;
		topdiri = 1;
	    }
	}
  }

  // this sets the array rhs
  FEM_FCT_ForConvDiff((TSquareMatrix3D*) matM, (TSquareMatrix3D*) mat,
		      Nodes, Nodes,
		      lump_mass_PSD, matrix_D_Entries_PSD,
		      sol, oldsol,
		      rhs, RhsArray, oldrhs_fem_fct0, tilde_u,
		      N_neum_to_diri, neum_to_diri,
		      NULL,NULL,NULL,
		      1, NULL, bdr_val);

  t2 = GetTime();
  OutPut("time bulkfct (3) " << t2-t1 << endl);
  // build matrix for FEM-FCT
  matM->Reset();
  FEM_FCT_SystemMatrix(matM, mat, lump_mass_PSD, Nodes);

  t2 = GetTime();
  OutPut("time bulkfct (4) " << t2-t1 << endl);

  // set Dirichlet boundary conditions in matrix
  for ( i=0 ; i<N_neum_to_diri ; i++ )
  {
      j = neum_to_diri[i];
      rhs[j] = bdr_val[i]; 
      sol[j] = rhs[j];       
      ii = row_ptr[j];
      jj = row_ptr[j+1];
      // off diagonals
      for ( iq = ii ; iq < jj ; iq++ )
      {
	  // diagonal entry
	  if(col_ptr[iq]==j)
	      entriesM[iq] = 1.0;
	  else
	      entriesM[iq] = 0;
      }
  }
  t2 = GetTime();
  OutPut("time bulkfct (5) " << t2-t1 << endl);
  if (sqrt(Ddot(Nodes,rhs,rhs)) > 0)
  {
    if (Nodes<SC_LDS)
    {
      //DirectSolver(matM, rhs, sol);
      OutPut("SolveLU MEMORY: " << setw(10) << GetMemory() << endl);
      //SolveLU(Nodes, N_Entries, row_ptr, col_ptr, entries, rhs, sol);
      OutPut("done SolveLU MEMORY: " << setw(10) << GetMemory() << endl);

    }
    else
    {
	t3 = TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
	TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = 0;
	Solver(matM,rhs,sol);
	TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = t3;
	// no output of solver data any longer
	TDatabase::ParamDB->SC_VERBOSE_AMG = 1;
    }
  }
  else
  {
    memset(sol,0,Nodes*SizeOfDouble);
  }
  
  maxsol = 0;
  // cut undershoots
  for (i=0;i<Nodes;i++)
  {
      if (sol[i] < 0)
	  sol[i] = 0;
      if (sol[i] > maxsol)
	  maxsol = sol[i];
  }
  memcpy(oldsol, sol, Nodes * SizeOfDouble);

//  OutPut(TDatabase::TimeDB->CURRENTTIME << " Solver done " << sqrt(Ddot(Nodes,sol,sol)) << " max " << maxsol << endl);

  // reset indices if necessary
  if (topdiri)
  {
      for (i=0;i<N3;i++)
	  neum_to_diri[i] = i;  
  }

  t2 = GetTime();
  OutPut("time bulkfct (6) " << t2-t1 << endl);

  // deletion of the arrays
  delete bdr_val;
  delete rhs;
  delete u1;
  delete test_cells;
}

/****************************************************************************************
*                                                                                       *
*                          Part IV : analysing routines                                 *
*                         ------------------------------                                *
*                                                                                       *
****************************************************************************************/


/****************************************************************************************
*                                                                                       *
*  feedback to concentration C                                                          *
*                                                                                       *
****************************************************************************************/

void Integral_For_Particle_Increase_Term(TFESpace3D *fespace, TFEFunction3D *fefct,
                                         int N_x, int N_y, int N_z, int N_a,
                                         double *x_coord, double *y_coord, double *z_coord, 
					 double *a_coord,
                                         double *concent_C_array, double *f)
{
  int i, j, N2, N3, i1, i2, m, k, l, N_Cells, local;
  int N_Vertices, index;
  int *GlobalNumbers, *BeginIndex, *DOF;

  double val, x, y, z, value, eps = 1e-6;
  double *integral;
  double c_C_infty_sat, C_2, d_p_max, d_p_0, d_p_min, c_C_infty;

  TCollection *Coll;
  TBaseCell *cell;
  int *indextest;
 
  switch(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY)
  {
     case 1356:
	 //d_p_max = TDatabase::ParamDB->UREA_D_P_MAX;
	 //d_p_0 = TDatabase::ParamDB->UREA_D_P_0;
	 //d_p_min = d_p_0/d_p_max;
               //C_2 =  c_C_infty_sat = c_C_infty = 1;
               break;
     default:
               c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
               C_2 = TDatabase::ParamDB->BULK_C_2;
               d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
               d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
               d_p_min = d_p_0/d_p_max;
               c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;
               break;
  }  

  N2 = (N_x+1)*(N_y+1);
  // number of grid points in 3D domain
  N3 = N2 * (N_z+1);
  // test array to check if all indices in 3D were computed
 
  indextest = new int[N3];
  memset(indextest,0,N3*SizeOfInt);

  // arrays for d.o.f. of fe space (Q1) which contains the integral
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  // values of the integral in the fe space
  integral = fefct->GetValues();

  // initialize
  for ( i=0 ; i<N3 ; i++ )
  {
    integral[i] = -1;
  }

  // all spaces use same Coll
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
	
  //loop over all mesh cells
  for( i=0 ; i<N_Cells ; i++ )
  {
    // get cell i  
    cell = Coll->GetCell(i);
    // get number of edges (or vertices)
    N_Vertices=cell->GetN_Vertices();
    if (N_Vertices!=8)
    {
	OutPut("Integral_For_Particle_Increase_Term only implemented for hexahedra !!!" << endl);
	exit(4711);
    }
    // get pointer to degrees of freedom connected to this mesh cell
    DOF = GlobalNumbers + BeginIndex[i];
    // loop over all vertices
    for ( j=0 ; j<N_Vertices ; j++ )
    {
      // corners of mesh cell
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      z = cell->GetVertex(j)->GetZ();
      // correspondance of local vertex and local dof
      switch(j)
      {
        case 0:
          local = 0;                   
          break;
        case 1:
          local = 1;
          break;
        case 2:
          local = 3;
          break;
        case 3:
          local = 2;
          break;
        case 4:
          local = 4;
          break;
        case 5:
          local = 5;
          break;
        case 6:
          local = 7;
          break;
        case 7:
          local = 6;
          break;
      }
      // l - global index
      l = DOF[local];
      // already done
      if ( integral[l]>=-0.1 )
      {
        continue;
      }
      // corresponding index of array f - comparision of the x, y and z coordinates
      for ( m=0 ; m<(N_x+1) ; m++ )
      {
        if ( fabs(x-x_coord[m])<eps )
        {
          index = m;
          break;
        }
      }
      //OutPut(x_coord[m]);

      for ( m=0 ; m<=N_y*(N_x+1) ; m+=(N_x+1) )
      {
        if ( fabs(y-y_coord[m])<eps )
        {
          index += m;
          break;
        }
      }

     for ( m=0 ; m<=N_z*N2 ; m+=N2 )
      {
        if ( fabs(z-z_coord[m])<eps )
        {
          index += m;
          break;
        }
      }

     // set test array to 1
      indextest[index]++;
      val = 0.0;

      // growth rate does not depend on PSD
      if (TDatabase::ParamDB->BULK_PB_DISC!=2)
      {
        // compute integral in a-direction, composed trapezoidal rule
        for ( k=1 ; k<=N_a ; k++ )
        {
          i1 = index+k*N3;
          i2 = index+(k-1)*N3;
	  val += (a_coord[k]*a_coord[k]*f[i1] + a_coord[k-1]*a_coord[k-1]*f[i2])
                 * (a_coord[k]-a_coord[k-1]) / 2.0;
        }
	//OutPut(" val " << val << endl);
      }
      else
      {
        // compute concentration of C, this array was set in the routine
	// for computing the psd
        value = concent_C_array[index];
        // compute integral in a-direction, composed trapezoidal rule
        for ( k=1 ; k<=N_a ; k++ )
        {
          i1 = index+k*N3;
          i2 = index+(k-1)*N3;
          val += (value-c_C_infty_sat/c_C_infty*exp(C_2/(a_coord[k]*d_p_max)))
                 * (a_coord[k]*a_coord[k]*f[i1] + a_coord[k-1]*a_coord[k-1]*f[i2])
                 * (a_coord[k]-a_coord[k-1]) / 2.0;
        }
      }

      // insert the value into the fe function
      integral[l] = val;
    }
  }

  // check if intergral has been computed for all dof
  for ( i=0 ; i<N3 ; i++ )
  {
    if (!indextest[i])
    {
      OutPut("Index " << i << " not found " << endl);
      exit(4711);
    }
  }

  OutPut(TDatabase::TimeDB->CURRENTTIME << " integral " << sqrt(Ddot(N3,integral,integral)) << endl);

  delete indextest;
}

#if 0 // this is not compiling
void Integral_For_Particle_Increase_Term_2D(TFESpace3D *fespace, TFEFunction3D *fefct,TFEFunction3D *fefct_2,
int N_x, int N_y, int N_z, int N_a,int N_b,
double *x_coord, double *y_coord, double *z_coord,
double *a_coord, double *b_coord,
double *f)
{
  int i,i1,i2,i3,i4, j, N2, N3,N4,N5, m, k,k1,k2, l, N_Cells, local;
  int N_Vertices, index;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int m1,m2,m3,n1,n2,n3;
  double sa,sb,sa2,sb2,val,val_1,val1,val2,val3, val_2,val1_2,val2_2,val3_2,val_t,x, y, z, value, eps = 1e-6,A1,A2;
  double i1a,i1b,i2a,i2b,i3a,i3b,i4a,i4b;
  double j1a,j1b,j2a,j2b,j3a,j3b,k1a,k1b,k2a,k2b,k3a,k3b,fs,fj1,fj2,fj3,fs2,fk1,fk2,fk3;
  double *integral1,*integral2,*integral,*a,*b;
  double c_C_infty_sat, C_2, d_p_max, d_p_0, d_p_min, c_C_infty;
  //double L_max_1 = TDatabase::ParamDB->UREA_D_P_MAX;
  //double L_max_2 = TDatabase::ParamDB->UREA_D_P_MAX;
  double L_max_1 = TDatabase::ParamDB->KDP_D_P_MAX;
  double L_max_2 = TDatabase::ParamDB->KDP_D_P_MAX_2;
  // double G1=1.0;
  //double G2=1.0;
  double val11, val21,val31, val12, val22, val32;
  double val11_2, val21_2,val31_2, val12_2, val22_2, val32_2;
  TCollection *Coll;
  TBaseCell *cell;
  int *indextest;

  //OutPut("diam_MAX" << d_p_max << endl);
  //OutPut("conc_infty" << c_C_infty  << endl);
  N2 = (N_x+1)*(N_y+1);
  // number of grid points in 3D domain
  N3 = N2 * (N_z+1);
  // test array to check if all indices in 3D were computed
  N4 = N3 * (N_a+1);
  N5 = N4 * (N_b+1);
  indextest = new int[N3];
  memset(indextest,0,N3*SizeOfInt);

  //f = new double[N5];
  // memset(f,0,N5*SizeOfDouble);

  //OutPut(TDatabase::TimeDB->CURRENTTIME << " f ist " << sqrt(Ddot(N5,f,f)) << endl);

  a =new double[N5];
  memset(a,0,N5*SizeOfDouble);
  b =new double[N5];
  memset(b,0,N5*SizeOfDouble);
  // arrays for d.o.f. of fe space (Q1) which contains the integral
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  // values of the integral in the fe space
  integral1 = fefct->GetValues();
  integral2 = fefct_2->GetValues();
  integral = fefct->GetValues();

  //OutPut("Integral bis hier !!!" << endl);
  // OutPut(TDatabase::TimeDB->CURRENTTIME << " integral " << sqrt(Ddot(N3,integral,integral)) << endl);
  // OutPut(TDatabase::TimeDB->CURRENTTIME << " integral1 " << sqrt(Ddot(N3,integral1,integral1)) << endl);
  //OutPut(TDatabase::TimeDB->CURRENTTIME << " integral2 " << sqrt(Ddot(N3,integral2,integral2)) << endl);

  // initialize
  for ( i=0 ; i<N3 ; i++ )
  {
    integral[i] = -1;
  }
  for ( i=0 ; i<N3 ; i++ )
  {
    integral1[i] = -1;
  }
  for ( i=0 ; i<N3 ; i++ )
  {
    integral2[i] = -1;
  }
  // all spaces use same Coll
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  //loop over all mesh cells
  for( i=0 ; i<N_Cells ; i++ )
  {
    // get cell i
    cell = Coll->GetCell(i);
    // get number of edges (or vertices)
    N_Vertices=cell->GetN_Vertices();
    if (N_Vertices!=8 && N_Vertices!=4)
    {
      OutPut("Integral_For_Particle_Increase_Term only implemented for hexahedra and tetrahedra!!!" << endl);
      exit(4711);
    }
    // get pointer to degrees of freedom connected to this mesh cell
    DOF = GlobalNumbers + BeginIndex[i];
    // loop over all vertices
    for ( j=0 ; j<N_Vertices ; j++ )
    {
      // corners of mesh cell
      if (N_Vertices == 8)
      {
        cell->GetVertex(j)->GetCoords(x,y,z);
      }
      else
        cell->GetVertex(j)->GetNormalCoords(x,y,z);
      // correspondance of local vertex and local dof
      // hexahedron
      if (N_Vertices == 8)
      {
        switch(j)
        {
          case 0:
            local = 0;
            break;
          case 1:
            local = 1;
            break;
          case 2:
            local = 3;
            break;
          case 3:
            local = 2;
            break;
          case 4:
            local = 4;
            break;
          case 5:
            local = 5;
            break;
          case 6:
            local = 7;
            break;
          case 7:
            local = 6;
            break;
        }
      }
      else
      // tetrahedron
        local = j;
      // l - global index
      l = DOF[local];
      // already done
      if ( integral[l]>=-0.1 )
      {
        continue;
      }
      if ( integral1[l]>=-0.1 )
      {
        continue;
      }
      if ( integral2[l]>=-0.1 )
      {
        continue;
      }
      // corresponding index of array f - comparision of the x, y and z coordinates
      for ( m=0 ; m<(N_x+1) ; m++ )
      {
        if ( fabs(x-x_coord[m])<eps )
        {
          index = m;
          break;
        }
      }
      //OutPut(x_coord[m]);

      for ( m=0 ; m<=N_y*(N_x+1) ; m+=(N_x+1) )
      {
        if ( fabs(y-y_coord[m])<eps )
        {
          index += m;
          break;
        }
      }

      for ( m=0 ; m<=N_z*N2 ; m+=N2 )
      {
        if ( fabs(z-z_coord[m])<eps )
        {
          index += m;
          break;
        }
      }

      // set test array to 1
      indextest[index]++;

      /*   for ( k1=0 ; k1<=N_b ; k1++)
        {
          for ( k2=0 ; k2<=N_a ; k2++ )
          {
              i1 = index+k2*N3+k1*N4;
               a[i1]=a_coord[k2];
         b[i1]=b_coord[k1];

         f[i1]=-7*a[i1]*b[i1]+2*b[i1]*b[i1];
      }
      }*/

      val=0.0;
      val1= 0.0;
      val2= 0.0;
      val_2 = 0.0;
      // growth rate does not depend on PSD

      // compute integral in a-direction, composed trapezoidal rule
      for ( k1=0 ; k1<N_b ; k1++)
      {
        for ( k2=0 ; k2<N_a ; k2++ )
        {
          // lower left node
          i1 = index+k2*N3+k1*N4;
          // lower right node
          i2 = i1+N3;
          // upper left node
          i3 = i1+N4;
          // uppe right node
          i4=  i1+N3+N4;

          // triangle i1 i4 i3
          //                   P1 i1a             P2 i1b                 P3  i4a                 P4 i4b                  P5 i3a                P6 i3b
          //OutPut(" " << a_coord[k2] << " " << b_coord[k1] <<" "<< a_coord[k2+1] << " " << b_coord[k1+1] << " " << a_coord[k2] << " " << b_coord[k1+1] <<endl);
          // coordinates of vertices
          i1a = a_coord[k2];
          i1b = b_coord[k1];
          i2a = a_coord[k2+1];
          i2b = b_coord[k1];
          i3a = a_coord[k2];
          i3b = b_coord[k1+1];
          i4a = a_coord[k2+1];
          i4b = b_coord[k1+1];
          // bary center
          sa = (i1a + i3a + i4a)/3.0;
          sb = (i1b + i3b + i4b)/3.0;
          // if upper triangle below a=b then also lower triangle
          if (sb < sa)
            continue;
          // mid points of edges
          j1a = (i1a + i3a)/2.0;
          j1b = (i1b + i3b)/2.0;
          j2a = (i1a + i4a)/2.0;
          j2b = (i1b + i4b)/2.0;
          j3a = (i3a + i4a)/2.0;
          j3b = (i3b + i4b)/2.0;

          // function values by linear interpolation
          fs = (f[i1]+f[i4]+f[i3])/3.0;

          fj1 = (f[i1]+f[i3])/2.0;
          fj2 = (f[i1]+f[i4])/2.0;
          fj3 = (f[i3]+f[i4])/2.0;

          //A1=abs((P3*P6 - P5*P4)-(P1*P6 - P5*P2)+(P1*P4 - P3*P2))/2
          // A1=fabs((a_coord[k2+1]*b_coord[k1+1]-a_coord[k2]*b_coord[k1+1])-(a_coord[k2]*b_coord[k1+1]-a_coord[k2]*b_coord[k1])+(a_coord[k2]*b_coord[k1+1]-a_coord[k2+1]*b_coord[k1]))/2.0;
          A1=fabs((i4a*i3b-i3a*i4b) -(i1a*i3b-i3a*i1b)+(i1a*i4b-i4a*i1b))/2.0;
          // values of function in edge mid points
          //integral for G1
          val11 =   2*( L_max_2*j1a*j1b -  L_max_1*j1a*j1a) * fj1;
          val21 =  2*( L_max_2*j2a*j2b - L_max_1*j2a*j2a) * fj2;
          val31 =  2*( L_max_2*j3a*j3b - L_max_1*j3a*j3a) * fj3;
          //integral for G1
          val12 = j1a*j1a * fj1;
          val22 = j2a*j2a * fj2;
          val32 = j3a*j3a * fj3;

          val1 += A1 * (val11 + val21 + val31)/3.0;
          val2 += A1 * (val12 + val22 + val32)/3.0;
          // OutPut(" " << i1 << " " << i2 <<" " << i4 <<endl);
          //triangle i1 i2 i4
          //                   P1 i1a                   P2 i1b           P3 i2a                 P4 i2b                  P5 i4a                P6 i4b
          //OutPut(" " << a_coord[k2] << " " << b_coord[k1] << " " << a_coord[k2+1] << " " << b_coord[k1] << " " << a_coord[k2+1] << " " << b_coord[k1+1] <<" "<<endl);

          sa2=(i1a + i2a + i4a)/3.0;
          sb2=(i1b + i2b + i4b)/3.0;
          // if bottom triangle below a=b then also lower triangle
          if (sb2<sa2)
            continue;

          //A2=abs((P3*P6 - P5*P4)-(P1*P6 - P5*P2)+(P1*P4 - P3*P2))/2
          // A2=fabs((a_coord[k2+1]*b_coord[k1+1]-a_coord[k2+1]* b_coord[k1])-(a_coord[k2]*b_coord[k1+1]-a_coord[k2+1]* b_coord[k1])+(a_coord[k2]*b_coord[k1]-a_coord[k2+1]* b_coord[k1]))/2.0;
          A2=fabs((i2a*i4b-i4a*i2b)-(i1a*i4b-i4a*i1b)+(i1a*i2b -i2a*i1b))/2.0;

          // mid points of edges
          k1a = (i1a + i2a)/2.0;
          k1b = (i1b + i2b)/2.0;
          k2a = (i2a + i4a)/2.0;
          k2b = (i2b + i4b)/2.0;
          k3a = (i1a + i4a)/2.0;
          k3b = (i1b + i4b)/2.0;

          // function values by linear interpolation
          fs2 = (f[i1]+f[i2]+f[i4])/3.0;

          fk1 = (f[i1]+f[i2])/2.0;
          fk2 = (f[i2]+f[i4])/2.0;
          fk3 = (f[i1]+f[i4])/2.0;
          //integral for G1
          val11_2 = 2* (L_max_2*k1a*k1b - L_max_1*k1a*k1a) * fk1;
          val21_2 = 2* (L_max_2*k2a*k2b - L_max_1*k2a*k2a) * fk2;
          val31_2 = 2* (L_max_2*k3a*k3b - L_max_1*j3a*k3a) * fk3;
          //integral for G2
          val12_2 =  k1a*k1a * fk1;
          val22_2 =  k2a*k2a * fk2;
          val32_2 =  k3a*k3a * fk3;

          val1 += A2 * (val11_2 + val21_2 + val31_2)/3.0;
          val2 += A2 * (val12_2 + val22_2 + val32_2)/3.0;

          val=val1+val2;
        }

      }
      //OutPut(" val1 " << val1 << endl);
      //  OutPut(" val2 " << val2 << endl);

      // OutPut(" val " << Func(5.,5.) << endl);
      // OutPut(" val1 " << val1 << endl);
      //OutPut(" val2 " << val2 << endl);
      // OutPut(" val " << val1+val2<< endl);
      // OutPut(" A1 " << A1 << endl);
      //  OutPut(" A2 " << A2 << endl);
      // insert the value into the fe function
      integral[l]=val;
      integral1[l] = val1;
      integral2[l] = val2;
      //integral1[l] = 1.;
      //integral2[l] = 2.;
      // if (abs(integral1[l])>1e-12)
      //  {
      //  OutPut(" j = " << j << " integral1[" << l << "]=" << integral1[l] << endl);
      //  }
      //  if (abs(integral2[l])>1e-12)
      //  {
      // OutPut(" j = " << j <<  " integral2[" << l << "]=" << integral2[l] << endl);
      // }
      // OutPut(" integral2 " << integral2[l] << endl);
      // OutPut(" val1 " << integral1[l] << endl);
      // OutPut(" val2 " << integral2[l] << endl);
    }
  }

  for ( i=0 ; i<N3 ; i++ )
  {
    if (!indextest[i])
    {
      OutPut("Index " << i << " not found " << endl);
      exit(4711);
    }
  }
  //OutPut(TDatabase::TimeDB->CURRENTTIME << " integral " << sqrt(Ddot(N3,integral,integral)) << endl);
  //OutPut( setprecision(7) << endl);
  OutPut(TDatabase::TimeDB->CURRENTTIME << " integral1 " << sqrt(Ddot(N3,integral1,integral1)) << endl);
  OutPut(TDatabase::TimeDB->CURRENTTIME << " integral2 " << sqrt(Ddot(N3,integral2,integral2)) << endl);
  //memset(integral1,0.0,N3*SizeOfDouble);
  //memset(integral2,0.0,N3*SizeOfDouble);
  //exit(1);
  delete indextest;
  delete a;
  delete b;
}
#endif // 0, the above method is not compiling

/****************************************************************************************
 *                                                                                       *
 *  feedback to concentration C  2D                                                        *
 *                                                                                       *
 ****************************************************************************************/

#if 0 // not compiling
void Integral_For_Particle_Increase_Term_2D_old(TFESpace3D *fespace, TFEFunction3D *fefct,
int N_x, int N_y, int N_z,int N_a,int N_b,
double *x_coord, double *y_coord, double *z_coord, double *a_coord,double *b_coord, double *f)
{
  int N_P,N_e,N_E,alpha;
  int N_Vertices, index,index1,index2;
  int *GlobalNumbers, *BeginIndex, *DOF;

  double val,val1,val2,val3, x, y, z, eps = 1e-6,sa,sb,A;
  double *integral,*a,*b;
  double c_C_infty_sat, C_2, d_p_max, d_p_0, d_p_min, c_C_infty;
  int i, j, N2, N3, N4, N5, m, k, l, N_Cells, local;
  TCollection *Coll;
  TBaseCell *cell;
  int *indextest;
  const int dim =3;
  int in3,k1,k2;
  double G1=0.4;
  double G2=0.3;
  int  i1=0;
  int  i2=0;
  int ii;
  int iii=0;
  int ii1;
  int ii2;
  int ii3;
  //N_e=N_a*N_b;
  N2 = (N_x+1)*(N_y+1);
  N3 = N2 * (N_z+1);
  N4 = N3 * (N_a+1);
  N5 = N4 * (N_b+1);

  // array for index
  int **in = new int*[N5];
  for (int i = 0; i < N5 ; i++)
    in[i] = new int[dim];

  int **in1 = new int*[N5];
  for (int i = 0; i < N5 ; i++)
    in1[i] = new int[dim];

  int **in2 = new int*[N5];
  for (int i = 0; i < N5 ; i++)
    in2[i] = new int[dim];

  for (int i = 0; i < N5 ; i++)
  {
    for (int j = 0; j < dim ; j++)
    {
      in[i][j]=0;
      in1[i][j]=0;
      in2[i][j]=0;
    }

  }

  indextest = new int[N3];
  memset(indextest,0,N3*SizeOfInt);

  a =new double[N5];
  memset(a,0,N5*SizeOfDouble);
  b =new double[N5];
  memset(b,0,N5*SizeOfDouble);

  // arrays for d.o.f. of fe space (Q1) which contains the integral
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  // values of the integral in the fe space
  integral = fefct->GetValues();

  // initialize
  for ( i=0 ; i<N3 ; i++ )
  {
    integral[i] = -1;
  }

  // all spaces use same Coll
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_P = (N_a+1)*(N_b+1);
  index=0;
  index1=0;
  index2=0;
  //loop over all mesh cells
  for( i=0 ; i<N_Cells ; i++ )
  {
    // get cell i
    cell = Coll->GetCell(i);
    // get number of edges (or vertices)
    N_Vertices=cell->GetN_Vertices();
    if (N_Vertices!=8)
    {
      OutPut("Integral_For_Particle_Increase_Term only implemented for hexahedra !!!" << endl);
      exit(4711);
    }
    // get pointer to degrees of freedom connected to this mesh cell
    DOF = GlobalNumbers + BeginIndex[i];
    // loop over all vertices
    for ( j=0 ; j<N_Vertices ; j++ )
    {

      // corners of mesh cell
      if (N_Vertices == 8)
      {
        cell->GetVertex(j)->GetCoords(x,y,z);
      }
      else
        cell->GetVertex(j)->GetNormalCoords(x,y,z);

      // correspondance of local vertex and local dof
      //hexahedron
      if (N_Vertices == 8)
      {
        switch(j)
        {
          case 0:
            local = 0;
            break;
          case 1:
            local = 1;
            break;
          case 2:
            local = 3;
            break;
          case 3:
            local = 2;
            break;
          case 4:
            local = 4;
            break;
          case 5:
            local = 5;
            break;
          case 6:
            local = 7;
            break;
          case 7:
            local = 6;
            break;
        }
      }
      else
      //tetrahedron
        local = j;

      // l - global index
      l = DOF[local];
      // already done
      if ( integral[l]>=-0.1 )
      {
        continue;
      }
      // corresponding index of array f - comparision of the x, y and z coordinates
      for ( m=0 ; m<(N_x+1) ; m++ )
      {
        if ( fabs(x-x_coord[m])<eps )
        {
          index = m;
          break;
        }
      }
      //OutPut(x_coord[m]);

      for ( m=0 ; m<=N_y*(N_x+1) ; m+=(N_x+1) )
      {
        if ( fabs(y-y_coord[m])<eps )
        {
          index += m;
          break;
        }
      }

      for ( m=0 ; m<=N_z*N2 ; m+=N2 )
      {
        if ( fabs(z-z_coord[m])<eps )
        {
          index += m;
          break;
        }
      }

      // set test array to 1
      indextest[index]++;

      val = 0.0;
      val1 = 0.0;
      val2 = 0.0;
      val3 = 0.0;

      for ( k1=0 ; k1<=N_a ; k1++ )
      {
        for ( k2=0 ; k2<=N_b ; k2++ )
        {
          ii=index+k1*N3+k2*N4;
          //ii=k1+(N_a+1)*k2;
          a[ii]=a_coord[k1];
          b[ii]=b_coord[k2];
          // OutPut("  " <<ii<<"  "<<a[ii]<<"  " <<b[ii]<<"  "<<endl);
          //index=i;
          //alpha=(ii+1)%(N_a+1);

          alpha=(ii+1)%N4;
          //if ((ii+1< N_P-N_a-1)&&(alpha!=0))
          if ((ii+1< N5-N4)&&(alpha!=0))
          {

            in1[index1][0]=ii;
            in1[index1][1]=ii+1;
            in1[index1][2]=ii+N_a+2;
            //OutPut(" test " <<a[5]<<"  " <<b[5]<<"  "<<endl);
            //OutPut("  " <<index1 <<"  " <<in1[index1][0]<<"  " <<  "  " <<in1[index1][1]<<"  " <<   "  " <<in1[index1][2]<<"  " <<   endl);
            // OutPut("  " <<a[ii]<<"  " <<b[i]<<"  "  <<a[i+1]<<"  " <<  "  " <<b[i+1] <<"  " << "  " <<a[i+N_a+2]<< "  " <<b[i+N_a+2]<<"  " <<  endl);
            //OutPut("  " <<ii<<"  "<<a[ii]<<"  " <<b[ii]<<"  "<<a[ii+1]<<"  " <<b[ii+1]<<"  "<<a[ii+N_a+2]<<"  " <<b[ii+N_a+2]<<"  "<<endl);

            in2[index2][0]=ii;
            in2[index2][1]=ii+N_a+2;
            in2[index2][2]=ii+N_a+1;
            //OutPut("  " <<a_coord[k1]<<"  " <<b_coord[k2] <<"  " << endl);

            //OutPut("  " <<a_coord[i+1]<<"  " <<b_coord[i+1] <<"  " << endl);
            //OutPut("  " <<a_coord[index2]<<"  " <<b_coord[index2] <<"  " << endl);
            //OutPut("  " <<a[i]<<"  " <<b[i] <<"  " << endl);
            //OutPut("  " <<a[index2]<<"  " <<b_coord[index2] <<"  " << endl);
            //OutPut("  " <<a_coord[i+N_a+2]<<"  " <<b_coord[i+N_a+2] <<"  " << endl);
            //OutPut("  " <<a_coord[i+N_a+1]<<"  " <<b_coord[i+N_a+1] <<"  " << endl);
            //OutPut("  " <<index2<<"  " <<in2[index2][0]<<"  " <<  "  " <<in2[index2][1]<<"  " <<  "  " <<in2[index][2]<<"  " <<   endl);
            //OutPut("  " <<index2 <<"  " <<in2[index2][0]<<"  " <<  "  " <<in2[index2][1]<<"  " <<   "  " <<in2[index2][2]<<"  " <<   endl);

            index1++;
            index2++;

          }
          //OutPut("  " <<a[i]<<"  " <<b[i]<<"  "   <<  endl);
        }
      }
      //  OutPut("  " <<in1[2][0]<<"  " <<   endl);
      //  OutPut("  " <<in1[2][1]<<"  " <<   endl);
      // OutPut("  " <<in1[2][2]<<"  " <<   endl);

      N_e=index1;
      // OutPut("  " <<N_e<<"  " <<   endl);
      for (i1=0;i1<index1;i1++)
      {

        in[i1][0]=in1[i1][0];
        in[i1][1]=in1[i1][1];
        in[i1][2]=in1[i1][2];
        //OutPut("  " <<i1<<"  " <<in1[i1][0]<<"  " <<  "  " <<in1[i1][1]<<"  " <<  "  " <<in1[i1][2]<<"  " <<   endl);
        // OutPut("  " <<a[in1[i1][0]]<<"  " <<b[in1[i1][0]]<<"  "  <<a[in1[i1][1]]<<"  " <<  "  " <<b[in1[i1][1]] <<"  " << "  " <<a[in1[i1][2]]<< "  " <<a[in1[i1][2]]<<"  " <<  endl);

      }
      for (i2=0;i2<N_e;i2++)
      {
        in[i1+i2][0]=in2[i2][0];
        in[i1+i2][1]=in2[i2][1];
        in[i1+i2][2]=in2[i2][2];
        //OutPut("  " <<i2<<"  " <<in2[i2][0]<<"  " <<  "  " <<in2[i2][1]<<"  " <<  "  " <<in2[i2][2]<<"  " <<   endl);
        // OutPut("  " <<a[in2[i2][0]]<<"  " <<b[in2[i2][0]]<<"  "  <<a[in2[i2][1]]<<"  " <<  "  " <<b[in2[i2][1]] <<"  " << "  " <<a[in2[i2][2]]<< "  " <<a[in2[i2][2]]<<"  " <<  endl);
      }
      N_E=2*N_e;
      for (ii=0;ii<2*N_e;ii++)
      {
        ii1=in[ii][0];
        ii2=in[ii][1];
        ii3=in[ii][2];

        // OutPut("  " <<i<<"  " <<ii1<<"  " <<  "  " <<ii2<<"  " <<  "  " <<ii3<<"  " <<   endl);
        //OutPut("  " <<P(i,1)<<"  " <<P(i,2)<<"  "  <<P(i,3)<<"  " <<  "  " <<P(i,4) <<"  " << "  " <<P(i,5)<< "  " <<P(i,6)<<"  " <<  endl);
        // OutPut("  " <<a[ii1]<<"  " <<b[ii1]<<"  "  <<a[ii2]<<"  " <<  "  " <<b[ii2] <<"  " << "  " <<a[ii3]<< "  " <<b[ii3]<<"  " <<  endl);

        sa =  (a[ii1]+a[ii2]+a[ii3])/3.0;
        sb =  (b[ii1]+b[ii2]+b[ii3])/3.0;
        if (sa>sb)
          continue;
        //  OutPut("  " <<ii<<"  " <<ii1<<"  " <<  "  " <<ii2<<"  " <<  "  " <<ii3<<"  " <<   endl);
        //OutPut("  " <<a[ii1]<<"  " <<b[ii1]<<"  "  <<a[ii2]<<"  " <<  "  " <<b[ii2] <<"  " << "  " <<a[ii3]<< "  " <<b[ii3]<<"  " <<  endl);
        A = fabs(( a[ii2]*b[ii3]-a[ii3]*b[ii2]) - (a[ii1]*b[ii3]-a[ii3]*b[ii1]) + (a[ii1]*b[ii2]-a[ii2]*b[ii1]))/2;
        val1 = G1*a[ii1]*b[ii1]+G1*b[ii1]*b[ii1]-G2*b[ii1]*b[ii1];
        val2 = G1*a[ii2]*b[ii2]+G1*b[ii2]*b[ii2]-G2*b[ii2]*b[ii2];
        val3 = G1*a[ii3]*b[ii3]+G1*b[ii3]*b[ii3]-G2*b[ii3]*b[ii3];
        //OutPut(" val1  " <<val1<<"  " <<   endl);
        //OutPut(" val2  " <<val2<<"  " <<   endl);
        //val=val+A*(val1*f[ii1]+val2*f[ii2]+val3*f[ii3])/3.0;
        //val=val+A*(val1+val2*+val3)/3.0;
        val=val+A*(f[ii1]+f[ii2]*+f[ii3])/3.0;
        // OutPut(" val  " <<val<<"  " <<   endl);

      }
      // insert the value into the fe function
      integral[l] = val;
      // OutPut(" Integr  " <<l <<"  " <<integral[l]<<"  " <<   endl);
    }

  }
  // check if intergral has been computed for all dof
  for ( i=0 ; i<N3 ; i++ )
  {
    if (!indextest[i])
    {
      OutPut("Index " << i << " not found " << endl);
      exit(4711);
    }
  }

  OutPut(TDatabase::TimeDB->CURRENTTIME << " integral " << sqrt(Ddot(N3,integral,integral)) << endl);

  delete indextest;
  delete a;
  delete b;

  for (int j = 0; j < N5 ; j++)
    delete [] in[j] ;
  delete [] in;

  for (int j = 0; j < N5 ; j++)
    delete [] in1[j] ;
  delete [] in1;

  for (int j = 0; j < N5 ; j++)
    delete [] in2[j] ;
  delete [] in2;
}
#endif // 0, the above method is not compiling


/****************************************************************************************
*                                                                                       *
*  Evaluation of population balance at outflow                                          *
*                                                                                       *
****************************************************************************************/

void Evaluate_f_at_outflow(int N_x, int N_y, int N_z, int N_a, 
			   double *x_coord, double *y_coord, 
			   double *a_coord, double *f,
			   double *average_median, int *average_step)
{
  int k, i0, i;
  double *size, *number, median;

  double d_p_min = TDatabase::ParamDB->BULK_D_P_0/TDatabase::ParamDB->BULK_D_P_MAX;
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;

  // centre of outflow is in (0.5,0.5,0)
  // compute index for value of f on a=0
  // this is correct only if N_x and N_y are even
  i0 = (int) ((N_x+1)*(N_y+1)-1)/2;

  if (( fabs(0.5-x_coord[i0])>1e-6 )||( fabs(0.5-y_coord[i0])>1e-6 ))
  {
      OutPut("Wrong index in Evalute_f_at_outflow " << i0 << endl);
    OutPut("N_x and N_y has to be even such that there is a d.o.f. at the center of the outflow " << endl);
    exit(4711);
  }

  size = new  double[N_a+1];
  number = new  double[N_a+1];
  for ( i=0 ; i<(N_a+1) ; i++ )
  {
      size[i] = a_coord[i]*d_p_max;
      number[i] =  f[i0+i*(N_x+1)*(N_y+1)*(N_z+1)]* f_infty;
      OutPut("PSD " << TDatabase::TimeDB->CURRENTTIME << " " <<
	     a_coord[i] << " " <<  f[i0+i*(N_x+1)*(N_y+1)*(N_z+1)] << " " <<
      size[i] << " " << number[i] << endl);
  }

  median =  calculate_dp_50(N_a+1,size,number);
  OutPut(TDatabase::TimeDB->CURRENTTIME << " median of particle size " << 
	 median);

  if (TDatabase::TimeDB->CURRENTTIME>=TDatabase::TimeDB->T1)
  {
      average_median[0] = median/(average_step[0] +1) + 
	  average_step[0] /(average_step[0] +1.0)*average_median[0] ;
      average_step[0] ++;
      OutPut(" average since " << TDatabase::TimeDB->T1 << " sec. : "<< average_median[0] );
      
  }
  delete size;
  delete number;
  OutPut(endl);
}


/*****************************************************************************************
 *                                                                                       *
 *  .vtk file to visualize the psd for a given cut_coordinate of the internal coordinate *
 *                                                                                       *
 ****************************************************************************************/

void write_vtk_file( int N_x, int N_y, int N_z, int N_a, double cut_coord,
                     double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                     double *f_old, const char *name)
{
  int i;
  int N_Cells = N_x*N_y*N_z;
  int N2 = (N_x+1)*(N_y+1);
  int N3 = N2*(N_z+1);
  // routine should work with the first N3-coordinates of x,y,z because a_coord is const. always in N3 blocks
  int N4 = N3*(N_a+1);


  FILE* out = fopen(name,"w");

  fprintf(out,"%s\n","# vtk DataFile Version 4.2" );
  fprintf(out,"%s\n","file created by MooNMD" );
  fprintf(out,"%s\n","ASCII" );
  fprintf(out,"%s\n","DATASET UNSTRUCTURED_GRID" );
  fprintf(out,"\n");
  fprintf(out,"%s","POINTS " );
  fprintf(out,"%i",N3);
  fprintf(out,"%s\n"," float" );

  for ( i=0 ; i<N3 ; i++ )
  {
    fprintf(out,"%f", x_coord[i]);
    fprintf(out,"%s"," ");
    fprintf(out,"%f", y_coord[i]);
    fprintf(out,"%s"," ");
    fprintf(out,"%f\n", z_coord[i]);
  }

  fprintf(out,"\n");
  fprintf(out,"\n");

  fprintf(out,"%s","CELLS ");
  fprintf(out,"%i",N_Cells);
  fprintf(out,"%s"," ");
  fprintf(out,"%i",N_Cells*9);
  fprintf(out,"\n");

  for ( i=0 ; i<N3 ; i++ )
  {
    if( (i<(N3-N2)) && (((i+1)%(N_x+1))!=0) && ((i%N2)<((N_x+1)*N_y)) )
    {
      fprintf(out,"%i",8);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i+1);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i+(N_x+1)+1);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i+(N_x+1));
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i+N2);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i+1+N2);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i+(N_x+1)+1+N2);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",i+(N_x+1)+N2);
      fprintf(out,"%s"," ");
      fprintf(out,"\n");
    }
  }

  fprintf(out,"\n");
  fprintf(out,"\n");

  fprintf(out,"%s","CELL_TYPES ");
  fprintf(out,"%i\n",N_Cells);

  for ( i=0 ; i<N_Cells ; i++ )
  {
    fprintf(out,"%i",12);
    fprintf(out,"%s"," ");
  }

  fprintf(out,"\n");
  fprintf(out,"\n");

  fprintf(out,"%s","POINT_DATA ");
  fprintf(out,"%i\n",N3);
  fprintf(out,"%s\n","SCALARS f float");
  fprintf(out,"%s\n","LOOKUP_TABLE default");

  for ( i=0 ; i<N4 ; i++ )
  {
    if ( a_coord[i]==cut_coord )
    {
      fprintf(out,"%f\n",f_old[i]);
    }
  }

  fprintf(out,"\n");
  fclose(out);
}


/*****************************************************************************************
 *                                                                                       *
 *  save f_old in txt file                                                               *
 *                                                                                       *
 ****************************************************************************************/

void save_f_old_in_txt_file( int Nodes, double *f_old, const char *name)
{
  FILE* out = fopen(name,"w");

  for ( int i=0 ; i<Nodes ; i++ )
  {
    fprintf(out,"%f\n", f_old[i]);
  }

  fclose(out);
}


/*****************************************************************************************
 *                                                                                       *
 *  function for the read in of f from the corresponding txt file                        *
 *                                                                                       *
 ****************************************************************************************/

void read_in_f_old_from_txt_file(int Nodes, double *f_old, const char *name)
{
  // composition of the path and the (changing) file name in: NAME
  char path[] ="./DATA_TMP/";
  int path_length = strlen(path);
  int name_length = strlen(name);

  //char NAME[path_length+name_length+1];
  char NAME[100];
  strcpy(NAME,"./DATA_TMP/");
  strcat(NAME,name);

  // opening of the data file
  FILE* data = fopen(NAME,"r");
  if ( data == NULL )
  {
    OutPut("Error in read_in_f_from_vtk_file: Can not open the data file" <<endl);
    exit(-1);
  }

  // read the data from the data file and filling of the array f
  char line[28];
  int counter = 0;

  while ( fgets(line,27,data) )
  {
    f_old[counter] = atof(line);
    counter++;
  }
}


void Compute_psd_level(TCollection *coll, TFESpace3D *psd_level_space,
		       TFEFunction3D *psd_level_fe_fct,
		       int N_x, int N_y, int N_z, int N_a, 
		       double *x_coord, double *y_coord, 
		       double *z_coord, double *a_coord, double *f)
{
    int i, j, k, l, local, N_Vertices,  N_Cells, index, start, ende;
    int *GlobalNumbers, *BeginIndex, *DOF, layer, found;
    double x, y, z, *psd_level_fct, eps = 1e-6;
    double f_infty = TDatabase::ParamDB->BULK_f_infty;
    TBaseCell *cell;

    layer = TDatabase::ParamDB->OUTPUT_NODE_LAYER_PSD;
    if ((layer<0)||(layer>TDatabase::ParamDB->N_CELL_LAYERS_PSD))
    {
	OutPut("psd node layer does not exist !!!" << endl);
	return;
    }
    // this is the node layer which will be written in the output
    start = layer * (N_x+1) * (N_y+1) * (N_z+1);
    ende =  start +  (N_x+1) * (N_y+1) * (N_z+1);

    GlobalNumbers = psd_level_space->GetGlobalNumbers();
    BeginIndex = psd_level_space->GetBeginIndex();
    psd_level_fct = psd_level_fe_fct->GetValues();

    N_Cells = coll->GetN_Cells();

    //loop over all mesh cells
    for( i=0 ; i<N_Cells ; i++ )
    {
	cell = coll->GetCell(i);
	// get number of edges (or vertices)
	N_Vertices=cell->GetN_Vertices();
	// get pointer to degrees of freedom connected to this mesh cell
	DOF = GlobalNumbers + BeginIndex[i];
	// loop over all vertices
	for ( j=0 ; j<N_Vertices ; j++ )
	{
	    // corners of mesh cell
	    x = cell->GetVertex(j)->GetX();
	    y = cell->GetVertex(j)->GetY();
	    z = cell->GetVertex(j)->GetZ();
	    // correspondance of local vertex and local dof
	    switch(j)
	    {
		case 0:
		    local = 0;                   
		    break;
		case 1:
		    local = 1;
		    break;
		case 2:
		    local = 3;
		    break;
		case 3:
		    local = 2;
		    break;
		case 4:
		    local = 4;
		    break;
		case 5:
		    local = 5;
		    break;
		case 6:
		    local = 7;
		    break;
		case 7:
		    local = 6;
		    break;
	    }
	    // l - global index
	    l = DOF[local];
	    // find the value in the psd_function f
	    for (k=start;k<ende;k++)
	    {
		found = 0;
		// coordinates are identical 
		if (((x_coord[k]-x)<eps)&&((y_coord[k]-y)<eps)&&((z_coord[k]-z)<eps))
		{
		    psd_level_fct[l] = f[k]*f_infty;
		    found = 1;
		    break;
		}
		// just for safety
		if (!found)
		{
		    OutPut("PSD value not found !!! " << endl);
		    exit(4711);
		}
	    }
	} // end j
    } // end i
}
