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
 *                         Bulk_2d3d.C                                                  *
 *                        -------------                                                 *
 *                                                                                      *
 *  routines for simulation of precipitation process in 2d cavity                       *
 *                                                                                      *
 ***************************************************************************************/
#ifdef __2D__

#include <SquareStructure2D.h>
#include <SquareMatrix2D.h>
#include <FEFunction2D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <Solver.h>
#include <FEM_TVD_FCT.h>
#include <Bulk_2d3d.h>
#include <Bulk.h>

#include <Database.h>
#include <TimeDiscRout.h>
#include <RKV_FDM.h>

#include <string.h>
#include <stdlib.h>

//extern "C" void SolveLU(int Nodes, int N_Entries, int *row_ptr, int *col_ptr, double *entries, double *rhs, double *sol);

/****************************************************************************************
 *                                                                                       *
 *  The declared functions are sorted in the following way :                             *
 *                                                                                       *
 *    1. general routines (e.g. grid generation,...)                                     *
 *    2. FDM routines                                                                    *
 *    3. FEM routines                                                                    *
 *    4. analysing routines                                                              *
 *                                                                                       *
 ****************************************************************************************/

/****************************************************************************************
 *                                                                                       *
 *                          Part I : general routines                                    *
 *                         ---------------------------                                   *
 *                                                                                       *
 ****************************************************************************************/

/****************************************************************************************
 *                                                                                       *
 * generates non-equidistant 3d grid for fem                                             *
 * input:  x_min, x_max, y_min, y_max, z_min, z_max                                      *
 *         are the intervals in x-, y- and z-direction                                   *
 *         N_x, N_y, N_z : number of intervals in x-, y- and z-direction                 *
 * output: *x_coord,  *y_coord,  *z_coord,  *z_layers_coord                              *
 *                                                                                       *
 * grid type 0: equidistant grid                                                         *
 * grid type 1: non-equidistant grid based on the cosinus function                       *
 * grid type 2: non-equidistant grid based on the hyperbolic tangent function            *
 *                                                                                       *
 ****************************************************************************************/

void grid_generator_3d(double x_min, double x_max, int N_x,
double y_min, double y_max, int N_y,
double z_min, double z_max, int N_z,
double *x_coord, double *y_coord, double *z_coord,
double *z_layers_coord)
{
  // check consistency
  if ( (x_min >= x_max) || (y_min >= y_max) || (z_min >= z_max) || (N_x <= 0) || (N_y <= 0) || (N_z <= 0) )
  {
    OutPut(" error : wrong input in grid_generator_3dfem " <<endl);
    exit(-1);
  }

  // parameters for the loops
  int Nodes = (N_x+1)*(N_y+1)*(N_z+1);
  int N2 = (N_x+1)*(N_y+1);

  // stretching parameter for grid type 2
  double gamma = TDatabase::ParamDB->CHANNEL_GRID_STRETCH;

  // arrays for the grid -> this allows to implement different distributions of grid points in every direction
  // eventuell mit switch Abfrage, falls mehrere Varianten erforderlich sind
  double *x_help, *y_help, *z_help;

  x_help = new double[N_x+1];  memset(x_help,0,(N_x+1)*SizeOfDouble);
  y_help = new double[N_y+1];  memset(y_help,0,(N_y+1)*SizeOfDouble);
  z_help = new double[N_z+1];  memset(z_help,0,(N_z+1)*SizeOfDouble);

  // im Funktionsprototyp ist grid_type = 0 gesetzt -> siehe Bulk_2d3d.h Zeile 24
  switch(TDatabase::ParamDB->GRID_TYPE)
  {
    // equidistant grid
    case 0:
      // fill auxiliary vectors with a distribution
      // x_help
      for ( int i=0 ; i<(N_x+1) ; i++ )
      {
        x_help[i] = x_min + ( ((x_max-x_min)/(N_x)) * i );
      }

      // y_help
      for ( int i=0 ; i<(N_y+1) ; i++ )
      {
        y_help[i] = y_min + ( ((y_max-y_min)/(N_y)) * i );
      }

      // z_help
      for ( int i=0 ; i<(N_z+1) ; i++ )
      {
        z_help[i] = z_min + ( ((z_max-z_min)/(N_z)) * i );
      }
      break;

      // non-equidistant grid, based on the cosinus function in the z-direction
    case 1:
      // fill auxiliary vectors with a distribution
      // x_help
      for ( int i=0 ; i<(N_x+1) ; i++ )
      {
        x_help[i] = x_min + ( ((x_max-x_min)/(N_x)) * i );
      }

      // y_help
      for ( int i=0 ; i<(N_y+1) ; i++ )
      {
        y_help[i] = y_min + ( ((y_max-y_min)/(N_y)) * i );
      }

      // z_help
      for ( int i=0 ; i<(N_z+1) ; i++ )
      {
        z_help[i] = z_max - (z_max-z_min) * cos(i*Pi/(2*N_z));
      }
      break;

      // non-equidistant grid, based on the hyperbolic tangent function in the z-direction
    case 2:
      // stretching parameter is gamma
      // fill auxiliary vectors with a distribution
      // x_help
      for ( int i=0 ; i<(N_x+1) ; i++ )
      {
        x_help[i] = x_min + ( ((x_max-x_min)/(N_x)) * i );
      }

      // y_help
      for ( int i=0 ; i<(N_y+1) ; i++ )
      {
        y_help[i] = y_min + ( ((y_max-y_min)/(N_y)) * i );
      }

      // z_help
      for ( int i=0 ; i<(N_z+1) ; i++ )
      {
        z_help[i] = z_max + (z_max-z_min) * (tanh( gamma*(((double)i/(double)N_z)-1.)))/(tanh(gamma));
      }
      break;

    default:
      OutPut(" error : wrong grid_type in grid_generator_3d " <<endl);
      exit(-1);
  }

  // fill the output vectors
  // x_coord
  for ( int i=0 ; i<Nodes ; i++ )
  {
    x_coord[i] = x_help[(i%(N_x+1))];
  }

  // y_coord
  for ( int k=0 ; k<(N_z+1) ; k++ )
  {
    for ( int j=0 ; j<(N_y+1) ; j++ )
    {
      for ( int i=0 ; i<(N_x+1) ; i++  )
      {
        y_coord[(i+(j*(N_x+1)))+(k*N2)] = y_help[j];
      }
    }
  }

  // z_coord
  for ( int j=0 ; j<(N_z+1) ; j++ )
  {
    for ( int i=0 ; i<N2 ; i++ )
    {
      z_coord[i+(j*N2)] = z_help[j];
    }
  }

  // generate the vector with the coordinates of the layers in z-direction
  // the entries are all positive -> in the dim-less form is z_min = d_p_min and z_max = 1.0
  for ( int i=0 ; i<(N_z+1) ; i++ )
  {
    z_layers_coord[i] = z_help[i];
  }

  delete x_help;
  delete y_help;
  delete z_help;
}

/***************************************************************************/
//
// computes boundary conditions for PSD on the inflow boundary wrt the flow
//
// boundary conditions have to be described on the closure of the inflow
//
/***************************************************************************/

int PSD_bound_cound_from_velo_inflow(double x, double y)
{
    int value=0, range;
    double eps = 1e-6;
    
    // check left side
    if  (fabs(x)<1e-6)
    {
	range = (int)TDatabase::ParamDB->P7;
	if ((y>range/32.0-eps)&&(y<(range+1)/32.0+eps))
        {
	    value = 3;
	}	
    }
    // check right side
    if  (fabs(1-x)<1e-6)
    {
	range = (int)TDatabase::ParamDB->P8;
	if ((y>range/32.0-eps)&&(y<(range+1)/32.0+eps))
	{
	    value = 1;
	}
    }
    // academic test example
    if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
    {
	if ((fabs(x)<1e-6)||(fabs(1-x)<1e-6)||(fabs(y)<1e-6)||(fabs(1-y)<1e-6))
	  value = 1;
    }

    return value;
}

/****************************************************************************************
 *                                                                                       *
 *                          Part II : FDM (explicit with upwindig)                       *
 *                         ---------------------------------------                       *
 *                                                                                       *
 ****************************************************************************************/

void Bulk_FWE_FDM_Upwind_3D(TCollection *coll,
TFEFunction2D *velocity1, TFEFunction2D *velocity2,
TFEFunction2D *concent_C,
double *f_old,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
double x_min, double x_max, double y_min, double y_max,
double z_min, double z_max,
double *velo1, double *velo2, double *concent_C_array,
			    int *correspond_2dgrid
)
{
  int i, ii, N2, N3, maxind, index, val;
  int very_first = 0;
  int alpha, beta, no_of_2dcell;
  double B_c_C, velocity1_array_val, velocity2_array_val, concent_C_array_val, maxsol, G_c_C_val;
  double values[3], t1, t2;
  double *f_new, *derx_val;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  double time = TDatabase::TimeDB->CURRENTTIME;

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
  double factor_G,c_C;

//  t1 = GetTime();
  // computed model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);
  // academic test example
  double oldtime;
  if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
  {
      //OutPut("f_infty " << f_infty << " G " << factor_G << endl);
      factor_G = 0.02;
      c_C_infty_sat=0;
      oldtime = time - TDatabase::TimeDB->TIMESTEPLENGTH;
   }
  // very first computation
  if (fabs (velo1[0] + 4711) < 1e-6)
  {
    very_first++;
    OutPut("very first computation of f" << endl);
  }

  //number of unknowns in 2D (on one layer)
  N2 = (N_x+1)*(N_y+1);
  // number of unknowns in 3D
  N3 = (N_x+1)*(N_y+1)*(N_z+1);

  // array for concentration of species C
  memset(concent_C_array, 0, N2*SizeOfDouble);
  // array for derivatives of f with respect to x,y,z
  derx_val = new double[N3];
  memset(derx_val, 0, N3*SizeOfDouble);
  // array for new values of the particle size distribution f
  f_new = new double[N3];
  memset(f_new, 0, N3*SizeOfDouble);

    // academic test example
  if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
  {
      if (fabs(time-TDatabase::TimeDB->TIMESTEPLENGTH) < 1e-4)
	  {
	      for (i=0;i<N3;i++)
	      {
		  f_old[i] = 1e3*(2.00-z_coord[i]*z_coord[i]*z_coord[i]-z_coord[i])
		   *(sin(Pi*x_coord[i])*sin(Pi*y_coord[i])+1)
                   *sin(Pi*0/2);
	      }
	  }
  }

  // discretization of PBE with FDM
  // loop over all nodes of the FDM grid
  for (i=0;i<N3;i++)
  {
    // node is on the first layer (z_coord=0)
    // in this layer, the velocity vector will be filled
    // since this vector is independent of z
    if (i < N2)
    {
      ii = i;
      // treat right d.o.f. seperately
      if (((ii+1)%(N_x+1)==0))
        ii = ii-1;
      // treat upper d.o.f. seperately
      if (ii>=N_y*(N_x+1))
        ii = ii-(N_x+1);
      // right corner
      //if (i==(N_y+1)*(N_x+1)-1)
      //  ii = i-(N_x+1)-1;
      alpha = (int)(ii/(N_x+1)+1e-6);
      beta = ii%(N_x+1);
      no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];

      // u1
      velocity1->FindGradientLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
      velo1[i] = values[0];
      velocity1_array_val = velo1[i];

      //u2
      velocity2->FindGradientLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
      velo2[i] = values[0];
      velocity2_array_val = velo2[i];

      // academic test example
      if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
      {
	  velo1[i] = 2*(2*y_coord[i]-1)*(1-(2*x_coord[i]-1)*(2*x_coord[i]-1)) 
	      * (sin(TDatabase::ParamDB->P1*Pi*time)+1)*TDatabase::ParamDB->P9; 
	  velocity1_array_val = velo1[i];
	  velo2[i] =  -2*(2*x_coord[i]-1)*(1-(2*y_coord[i]-1)*(2*y_coord[i]-1)) 
	      *(sin(TDatabase::ParamDB->P1*Pi*time)+1)*TDatabase::ParamDB->P9;
	  velocity2_array_val = velo2[i];
      }
      // fill the array for the concentrations of C
      // since this concentration does not depend on z
      // c_C
      concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
      concent_C_array[i] = values[0];
      concent_C_array_val = values[0];
    }
    //node is not on the first layer
    else
    {
      //the corresponding node on the first layer
      ii = i - N2*((int)(i/N2));
      // compute the value of the velocity for the corresponding (x,y) coordinates
      velocity1_array_val = velo1[ii];
      velocity2_array_val = velo2[ii];
      // compute the value of the concentration of C for the corresponding (x,y) coordinates
      concent_C_array_val = concent_C_array[ii];
    }
    if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(z_coord[i]*d_p_max)));
    }
    else
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty);
    }

    // compute the coefficients corresponding to the 3d finite-difference-method (7-point)
    // simple upwind scheme

    // compute the term A
    if (velocity1_array_val >= 0)
    {
      // not on boundary x = x_min (left)
      if ( i%(N_x+1)!=0 )
        derx_val[i] += velocity1_array_val*(f_old[i]-f_old[i-1])/(fabs(x_coord[i]-x_coord[i-1]));
    }
    else
    {
      // not on boundary x = x_max (right)
      if ( (i+1)%(N_x+1)!=0 )
        derx_val[i] += velocity1_array_val*(f_old[i+1]-f_old[i])/(fabs(x_coord[i+1]-x_coord[i]));
    }

    index = i%N2;
    // compute the term B
    if (velocity2_array_val >= 0)
    {
      // not on boundary y = y_min (front)
      if (index>N_x )
      {
        derx_val[i] += velocity2_array_val*(f_old[i]-f_old[i-(N_x+1)])/(fabs(y_coord[i]-y_coord[i-(N_x+1)]));
      }
    }
    else
    {
      // not on boundary y = y_max (back)
      if (index<(N_x+1)*N_y )
      {
        derx_val[i] += velocity2_array_val*(f_old[i+(N_x+1)]-f_old[i])/(fabs(y_coord[i+(N_x+1)]-y_coord[i]));
      }
    }

    // compute the term C
    if (G_c_C_val >= 0)
    {
      // not on boundary z = z_min (bottom)
      if ( i>((N_x+1)*(N_y+1)-1) )
        derx_val[i] += G_c_C_val*(f_old[i]-f_old[i-N2])/(fabs(z_coord[i]-z_coord[i-N2]));
    }
    else
    {
      // not on boundary z = z_max (top)
      if ( i<((N_x+1)*(N_y+1)*N_z) )
        derx_val[i] += G_c_C_val*(f_old[i+N2]-f_old[i])/(fabs(z_coord[i+N2]-z_coord[i]));
    }

    // academic test example
    if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
    {
	// -f_t
	derx_val[i] -= 1e3*(2.00-z_coord[i]*z_coord[i]*z_coord[i]-z_coord[i])
	    *(sin(Pi*x_coord[i])*sin(Pi*y_coord[i])+1)
	    *cos(Pi*time/2)*Pi/2.0;
	// - u \cdot \nabla f
	derx_val[i] -= velocity1_array_val * 1e3*(2.00-z_coord[i]*z_coord[i]*z_coord[i]-z_coord[i])
	    *(Pi*cos(Pi*x_coord[i])*sin(Pi*y_coord[i]))
	    *sin(Pi*time/2);
	derx_val[i] -= velocity2_array_val * 1e3*(2.00-z_coord[i]*z_coord[i]*z_coord[i]-z_coord[i])
	    *(Pi*sin(Pi*x_coord[i])*cos(Pi*y_coord[i]))
	    *sin(Pi*time/2);
	// -G f_p
	c_C = 0.6*(sin(2*Pi*x_coord[i])*sin(2*Pi*y_coord[i])+1)*sin(Pi*time/2)*time;
	derx_val[i] -= factor_G*c_C*1e3*(-3*z_coord[i]*z_coord[i]-1)
	    *(sin(Pi*x_coord[i])*sin(Pi*y_coord[i])+1)
	    *sin(Pi*time/2);
    }
/*
    if (i>=N2&& i<2*N2)
    {
	OutPut(i << " " << f_old[i] << " " << deltat << " " << derx_val[i] << " "
	       << f_old[i] - deltat * derx_val[i] << endl);
	       }*/
    // compute new particle size distribution
    f_new[i] = f_old[i] - deltat * derx_val[i];

    // set Dirichlet boundary conditions
    // if convection is positive at the bottom
    if (i<N2)
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
      // compute G*n, n=(0,0,-1);
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
      // academic test example
      if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
      {
	  f_new[i] = 1e3*(sin(Pi*x_coord[i])*sin(Pi*y_coord[i])+1) *
	    (2.00 - z_coord[i]*z_coord[i]*z_coord[i]-z_coord[i]) * sin(Pi*time/2);
      }
    }
    // set Dirichlet boundary conditions
    // if convection is positive at the top
    if (i>=N2*N_z)
    {
      if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
      {
        G_c_C_val = k_g*(c_C_infty*concent_C_array_val- c_C_infty_sat*exp(C_2/d_p_0));
      }
      else
      {
        G_c_C_val = k_g*(c_C_infty*concent_C_array_val-c_C_infty_sat);
      }
      // compute G*n, n=(0,0,1);
      if (G_c_C_val*f_infty < 0)
      {
	  if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)>1e-6)
	      f_new[i] = 0.0;
      }
    }

    // set Dirichlet boundary conditions
    // inflow from the left x = x_min (left) or right x = x_max (right)
    if ((( i%(N_x+1)==0 )  ||  ((i+1)%(N_x+1)==0 ))&&(i>N2))
    {
	val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i]);
	if (val)
	     f_new[i] = 0.0;
    }
    // academic test example
    if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
    {
	val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i]);
	if (val)
	    //if ((( i%(N_x+1)==0 )  ||  ((i+1)%(N_x+1)==0 ))
	    // || (i%N2<=N_x) || (i%N2 >=(N_x+1)*N_y))
	    f_new[i] = 1e3*(sin(Pi*x_coord[i])*sin(Pi*y_coord[i])+1) *
		(2.00 - z_coord[i]*z_coord[i]*z_coord[i]-z_coord[i]) * sin(Pi*time/2);
    }
  }

  // copy new particle size distribution into array for old one
  Dcopy(N3, f_new, f_old);

  maxsol =  0;
  maxind = -4711;

  // cut undershoots
  for (i=0;i<N3;i++)
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
  OutPut(time << " maxsol " << maxsol << endl);
  // free allocated memory
  delete f_new;
  delete derx_val;
}

void Bulk_RKV_FDM_3D(TCollection *coll,
         TFEFunction2D *velocity1, TFEFunction2D *velocity2,
         TFEFunction2D *concent_C,
         double *f_old, double **stages, 
         int N_x, int N_y, int N_z,
         double *x_coord, double *y_coord, double *z_coord,
         double *velo1, double *velo2, double *concent_C_array,
         int *correspond_2dgrid)
{
    int i, ii, N2, N3, maxind, val, N1_[3], N_, N_stages;
    int very_first = 0, disctype, time_disc;
  int alpha, beta, no_of_2dcell;
  int *offset_ = NULL, *offset1_;
  double B_c_C, velocity1_array_val, velocity2_array_val, concent_C_array_val, maxsol, G_c_C_val;
  double values[3], t1, t2, coeff[3];
  double *current_stage_fdm, *sol_curr;
  double *coordinates[3];
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  double time = TDatabase::TimeDB->CURRENTTIME;

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
  double factor_G,c_C;

//  t1 = GetTime();
  // computed model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);

  double oldtime;
  // very first computation
  if (fabs (velo1[0] + 4711) < 1e-6)
  {
    very_first++;
    OutPut("very first computation of f" << endl);
  }
  N2 = (N_x+1)*(N_y+1);
  // number of unknowns in 3D
  N3 = (N_x+1)*(N_y+1)*(N_z+1);

  // save parameters
  disctype = TDatabase::ParamDB->DISCTYPE;
  time_disc = TDatabase::TimeDB->TIME_DISC;

  // array for concentration of species C
  memset(concent_C_array, 0, N2*SizeOfDouble);
  // array for solution of current state
  sol_curr = stages[5];
  
  TDatabase::ParamDB->DISCTYPE = TDatabase::ParamDB->PB_DISC_TYPE;
  TDatabase::TimeDB->TIME_DISC = TDatabase::ParamDB->PB_TIME_DISC;

  N_stages = GetN_SubSteps();

  N1_[0] = N_x+1;
  N1_[1] = N_y+1;
  N1_[2] = N_z+1;
  coordinates[0] = x_coord;
  coordinates[1] = y_coord;
  coordinates[2] = z_coord;
  InitializeConvectiveTermFDM(3, offset_, offset1_, N1_);  

  /*  for (i=0;i<N3;i++)
  {
      f_old[i] = sin(i*1.0/N3);
      }*/

 // loop over the stages
  for (N_ = 0; N_ < N_stages; N_++)
  {
      SetTimeDiscParameters(1);
      for (i=0;i<N_stages;i++)
    OutPut(" A("<<N_<<","<<i<<") = " << TDatabase::TimeDB->RK_A[N_][i]);
      OutPut(" : c("<<N_<<") = " << TDatabase::TimeDB->RK_c[N_]);
      OutPut(" : b("<<N_<<") = " << TDatabase::TimeDB->RK_b[N_] << endl);
      current_stage_fdm = stages[N_];
      memset(current_stage_fdm, 0, N3 * SizeOfDouble);
      // initialize current stage
      memcpy(sol_curr,f_old,N3*SizeOfDouble);
      // add previous stages
      for (i=0;i<N_;i++)
    Daxpy(N3, deltat*TDatabase::TimeDB->RK_A[N_][i], stages[i], sol_curr);
      // compute next stage
      // discretization of PBE with FDM
      // loop over all nodes of the FDM grid
      for (i=0;i<N3;i++)
      {
    // node is on the first layer (z_coord=0)
    // in this layer, the velocity vector will be filled
    // since this vector is independent of z
    if ((N_==0) && (i < N2))
    {
        ii = i;
        // treat right d.o.f. seperately
        if (((ii+1)%(N_x+1)==0))
      ii = ii-1;
        // treat upper d.o.f. seperately
        if (ii>=N_y*(N_x+1))
      ii = ii-(N_x+1);
        // right corner
        //if (i==(N_y+1)*(N_x+1)-1)
        //  ii = i-(N_x+1)-1;
        alpha = (int)(ii/(N_x+1)+1e-6);
        beta = ii%(N_x+1);
        no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];
        
        // u1
        velocity1->FindGradientLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
        velo1[i] = values[0];
        velocity1_array_val = velo1[i];
        
        //u2
        velocity2->FindGradientLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
        velo2[i] = values[0];
        velocity2_array_val = velo2[i];
        
        // fill the array for the concentrations of C
        // since this concentration does not depend on z
        // c_C
        concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
        concent_C_array[i] = values[0];
        concent_C_array_val = values[0];
    }
    //node is not on the first layer or not first stage
    else
    {
        //the corresponding node on the first layer
        ii = i - N2*((int)(i/N2));
        // compute the value of the velocity for the corresponding (x,y) coordinates
        velocity1_array_val = velo1[ii];
        velocity2_array_val = velo2[ii];
        // compute the value of the concentration of C for the corresponding (x,y) coordinates
        concent_C_array_val = concent_C_array[ii];
    }
    //velocity1_array_val = cos(i);
    //velocity2_array_val = cos(i)*sin(i);
    //concent_C_array_val = atan(i)+0.2;
    // growth rate
    if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
    {
        G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(z_coord[i]*d_p_max)));
    }
    else
    {
        G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty);
    }
    coeff[0] = velocity1_array_val;
    coeff[1] = velocity2_array_val;
    coeff[2] = G_c_C_val;
    
    // compute convection term 
    ConvectiveTermFDM(3, i, 
          coeff, sol_curr, current_stage_fdm, coordinates,
          offset_, offset1_);

    // set Dirichlet boundary conditions
    // if convection is positive at the bottom
    if (i<N2)
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
        // compute G*n, n=(0,0,-1);
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
    if (i>=N2*N_z)
    {
        if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
        {
      G_c_C_val = k_g*(c_C_infty*concent_C_array_val- c_C_infty_sat*exp(C_2/d_p_0));
        }
        else
        {
      G_c_C_val = k_g*(c_C_infty*concent_C_array_val-c_C_infty_sat);
        }
        // compute G*n, n=(0,0,1);
        if (G_c_C_val*f_infty < 0)
        {
      current_stage_fdm[i] = 0.0;
      if (N_ == N_stages -1)
          f_old[i] = 0.0;
        }
    }
    
    // set Dirichlet boundary conditions
    // inflow from the left x = x_min (left) or right x = x_max (right)
    if ((( i%(N_x+1)==0 )  ||  ((i+1)%(N_x+1)==0 ))&&(i>N2))
    {
        val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i]);
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
      Daxpy(N3, deltat*TDatabase::TimeDB->RK_b[i], stages[i], f_old);
  }

  maxsol =  0;
  maxind = -4711;

  // cut undershoots
  for (i=0;i<N3;i++)
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
   // restore parameters
  TDatabase::ParamDB->DISCTYPE = disctype;
  TDatabase::TimeDB->TIME_DISC = time_disc;
  delete[] offset_;
}


/*****************************************************************************************
 *                                                                                       *
 *                          Part III : FDM (implicit Euler with upwindig)                *
 *                         ----------------------------------------------                *
 *                                                                                       *
 *****************************************************************************************/

//assemble the matrix corresponding to the 3d finite-difference-method (7-point)
void Bulk_BWE_FDM_Upwind_3D(TCollection *coll,
			    TFEFunction2D *velocity1, TFEFunction2D *velocity2,
			    TFEFunction2D *concent_C,
			    double *sol,
			    int *correspond_2dgrid,
			    int N_x, int N_y, int N_z,
			    double *x_coord, double *y_coord, double *z_coord,
			    TSquareMatrix2D *mat)
{
  int i, ii, N2, N3, maxind, index, index1, index2;
  int very_first = 0, jj, iq, k;
  int alpha, beta, no_of_2dcell, N_Entries, range, diag_index;
  int *col_ptr, *row_ptr, indices[9];
  int SC_LDS =  TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  double B_c_C, velocity1_array_val, velocity2_array_val, concent_C_array_val, maxsol, G_c_C_val;
  double values[3], *velo1, * velo2, *concent_C_array, yq, t1, t2, val, z_min, t3;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  double *entries, *rhs;
  double C_val[3];
  double time = TDatabase::TimeDB->CURRENTTIME;

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
  double factor_G;

  // computed model constants
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);

  // data of the matrix
  entries = mat->GetEntries();
  N_Entries = mat->GetN_Entries();
  col_ptr = mat->GetKCol();
  row_ptr = mat->GetRowPtr();

  //number of unknowns in 2D (on one layer)
  N2 = (N_x+1)*(N_y+1);
  // number of unknowns in 3D
  N3 = (N_x+1)*(N_y+1)*(N_z+1);

  // rhs
  rhs = new double[N3];
  // initialization
  memset(rhs,0,N3*SizeOfDouble);
  memset(entries,0, N_Entries*SizeOfDouble);

  // arrays for velocity
  velo1 = new double[N2];
  memset(velo1, 0, N2*SizeOfDouble);
  velo2 = new double[N2];
  memset(velo2, 0, N2*SizeOfDouble);

  // array for concentration of species C
  concent_C_array = new double[N2];
  memset(concent_C_array, 0, N2*SizeOfDouble);

  z_min = z_coord[0];

  // discretization of PBE with FDM
  // loop over all nodes of the FDM grid
  for (i=0;i<N3;i++)
  {
    // node is on the first layer (z_coord=0)
    // in this layer, the velocity vector will be filled
    // since this vector is independent of z
    if (i < N2)
    {
      ii = i;
      // treat right d.o.f. seperately
      if (((ii+1)%(N_x+1)==0))
        ii = ii-1;
      // treat upper d.o.f. seperately
      if (ii>=N_y*(N_x+1))
        ii = ii-(N_x+1);
      // right corner
      //if (i==(N_y+1)*(N_x+1)-1)
      //  ii = i-(N_x+1)-1;
      alpha = (int)(ii/(N_x+1)+1e-6);
      beta = ii%(N_x+1);
      no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];

      // u1
      velocity1->FindGradientLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
      velo1[i] = values[0];
      velocity1_array_val = velo1[i];

      //u2
      velocity2->FindGradientLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
      velo2[i] = values[0];
      velocity2_array_val = velo2[i];

      // fill the array for the concentrations of C
      // since this concentration does not depend on z
      // c_C
      concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],values);
      concent_C_array[i] = values[0];
      concent_C_array_val = values[0];
    }
    //node is not on the first layer
    else
    {
      //the corresponding node on the first layer
      ii = i - N2*((int)(i/N2));
      // compute the value of the velocity for the corresponding (x,y) coordinates
      velocity1_array_val = velo1[ii];
      velocity2_array_val = velo2[ii];
      // compute the value of the concentration of C for the corresponding (x,y) coordinates
      concent_C_array_val = concent_C_array[ii];
    }
    if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty*exp(C_2/(z_coord[i]*d_p_max)));
    }
    else
    {
      G_c_C_val = factor_G*(concent_C_array_val- c_C_infty_sat/c_C_infty);
    }

    // compute the coefficients corresponding to the 3d finite-difference-method (7-point)
    // simple upwind scheme

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
      
    // compute the term A, convection in x-direction
    if (velocity1_array_val >= 0)
    {
      // not on boundary x = x_min (left)
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
      // not on boundary x = x_max (right)
      if ( (i+1)%(N_x+1)!=0 )
      {
	  val = deltat * velocity1_array_val/fabs(x_coord[i+1]-x_coord[i]);
	  index = indices[2];
	  entries[index] += val;
	  entries[diag_index] -= val;
      }
    }

    index = i%((N_x+1)*(N_y+1));
    // compute the term B, convection in y direction
    if (velocity2_array_val >= 0)
    {
      // not on boundary y = y_min (front)
      if (index>N_x )
      {
        index1 = i-(N_x+1);
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
            if (col_ptr[k] == index1)
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
      // not on boundary y = y_max (back)
      index2 = (N_x+1)*N_y;
      if (index<index2 )
      {
        index1 = i+(N_x+1);
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          // fabs: because of possible negative coordinates
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

    // compute the term C
    if (G_c_C_val >= 0)
    {
      // not on boundary z = z_min (bottom)
      if ( i>((N_x+1)*(N_y+1)-1) )
      {
        index1 = i-N2;
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          // fabs: because of possible negative coordinates
          if (col_ptr[k] == index1)
          {
	      val = deltat * G_c_C_val/fabs(z_coord[i]-z_coord[index1]);
	      entries[k] -= val;
	      entries[diag_index] += val;
	    break;
          }
        }
      }
    }
    else
    {
      // not on boundary z = z_max (top)
      if ( i<((N_x+1)*(N_y+1)*N_z) )
      {
        index1 = i+N2;
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          // fabs: because of possible negative coordinates
          if (col_ptr[k] == index1)
          {
	      val = deltat * G_c_C_val/fabs(z_coord[index1]-z_coord[i]);
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

  // set Dirichlet boundary conditions, on all (possible) inflow boundaries
  for ( i=0 ; i<N2 ; i++ )
  {
    ii = i;
    // treat right d.o.f. seperately
    if (((ii+1)%(N_x+1)==0))
      ii = ii-1;
    // treat upper d.o.f. seperately
    if (ii>=N_y*(N_x+1))
      ii = ii-(N_x+1);
    // right corner
    //if (i==(N_y+1)*(N_x+1)-1)
    //  ii = i-(N_x+1)-1;
    alpha = (int)(ii/(N_x+1)+1e-6);
    beta = ii%(N_x+1);
    no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];
    concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],C_val);
    if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
    {
	G_c_C_val = k_g*(C_val[0]-c_C_infty_sat*exp(C_2/z_min));
    }
    else
    {
	G_c_C_val = k_g*(c_C_infty*C_val[0]-c_C_infty_sat);
    }
    // compute G*n, n=(0,0,-1);
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

    // set top since negative G_c_C_val possible
    if (G_c_C_val < 0)
    {
      k = i + N2*N_z;
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
  for (i=N2;i<N3;i++)
  {
      iq =  PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i]);
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

  if (sqrt(Ddot(N3,rhs,rhs)) > 0)
  {
    if (N3<SC_LDS)
    {
      DirectSolver(mat, rhs, sol);
      OutPut("SolveLU MEMORY: " << setw(10) << GetMemory() << endl);
      //SolveLU(Nodes, N_Entries, row_ptr, col_ptr, entries, rhs, sol);
      OutPut("done SolveLU MEMORY: " << setw(10) << GetMemory() << endl);

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
    memset(sol,0,N3*SizeOfDouble);
  }

  // cut undershoots
  for (i=0;i<N3;i++)
  {
    if (sol[i] < 0)
      sol[i] = 0;
  }

  // free allocated memory
  delete rhs;
  delete velo1;
  delete velo2;
  delete concent_C_array;
}

/****************************************************************************************
 *                                                                                       *
 *                          Part III : FEM                                               *
 *                         ----------------                                              *
 *                                                                                       *
 ****************************************************************************************/

/****************************************************************************************
 *                                                                                       *
 * computes the value of the gradient for Q1                                             *
 * input:  *coeff, x, y, z                                                               *
 * output: *val                                                                          *
 *                                                                                       *
 ****************************************************************************************/

void Compute_Q1_Value_Gradient(double *coeff, double x, double y, double z, double *val)
{
  // precomputing of some products
  double xz, x4, x5, z6, y7;
/*
  xy = x*y;
  xz = x*z;
  yz = y*z;

  // value
  test[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z
    + coeff[4] * xy +  coeff[5] * xz +  coeff[6] * yz + coeff[7] * xy * z;
  // value of the x derivative
  test[1] = coeff[1] + coeff[4] * y + coeff[5] * z + coeff[7] * yz;

  // value of the y derivative
  test[2] = coeff[2] + coeff[4] * x + coeff[6] * z + coeff[7] * xz;

  // value of the z derivative
  test[3] = coeff[3] + coeff[5] * x + coeff[6] * y + coeff[7] * xy;
*/
  x4 =  coeff[4]*x;
  x5 =  coeff[5]*x;
  z6 =  coeff[6]*z;
  y7 =  coeff[7]*y;
  xz = x*z;

  // value
  val[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z
    + x4 * y +  x5 * z +  y * z6 + xz * y7;
  // value of the x derivative
  val[1] = coeff[1] + coeff[4] * y + coeff[5] * z + y7 * z;

  // value of the y derivative
  val[2] = coeff[2] + x4 + z6 + coeff[7] * xz;

  // value of the z derivative
  val[3] = coeff[3] + x5 + coeff[6] * y + x*y7;
/*
  if (fabs(test[0]-val[0])>1e-6)
  {
      OutPut("0 " << test[0] << " " << val[0] << endl);
      exit(1);
  }
  if (fabs(test[1]-val[1])>1e-6)
  {
      OutPut("1 " << test[1] << " " << val[1] << endl);
      exit(1);
  }
  if (fabs(test[2]-val[2])>1e-6)
  {
      OutPut("2 " << test[2] << " " << val[2] << endl);
      exit(1);
  }
  if (fabs(test[3]-val[3])>1e-6)
  {
      OutPut("3 " << test[3] << " " << val[3] << endl);
      exit(1);
      }*/
}

void Compute_Q1_Value(double *coeff, double x, double y, double z, double *val)
{
  // precomputing of some products
  double xz, x4, x5, z6, y7;

  x4 =  coeff[4]*x;
  x5 =  coeff[5]*x;
  z6 =  coeff[6]*z;
  y7 =  coeff[7]*y;
  xz = x*z;

  // value
  val[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z
    + x4 * y +  x5 * z +  y * z6 + xz * y7;
}

void Compute_Q1_Gradient(double *coeff, double x, double y, double z, double *val)
{
  // precomputing of some products
  double xz, x4, x5, z6, y7;

  x4 =  coeff[4]*x;
  x5 =  coeff[5]*x;
  z6 =  coeff[6]*z;
  y7 =  coeff[7]*y;
  xz = x*z;

  // value
  val[0] = 0;

  // value of the x derivative
  val[1] = coeff[1] + coeff[4] * y + coeff[5] * z + y7 * z;

  // value of the y derivative
  val[2] = coeff[2] + x4 + z6 + coeff[7] * xz;

  // value of the z derivative
  val[3] = coeff[3] + x5 + coeff[6] * y + x*y7;

}

/****************************************************************************************
*                                                                                       *
* this routine computes the corresponding 2d grid                                       *
* this is the lexicographical numbering of the mesh cells                               *
* input:  N_x, N_y, x_coord, y_coord, coll                                              *
* output: correspond_2dgrid                                                             *
*                                                                                       *
****************************************************************************************/

void generate_correspond_2d_grid(int N_x, int N_y, double *x_coord, double *y_coord,
TCollection *coll, int *correspond_2dgrid)
{
  int i,j,N_Edges, N_Cells, alpha, beta;
  double  x_low, y_low, sx, sy;
  TBaseCell *cell;               

  // number of mesh cells
  N_Cells = coll->GetN_Cells();

  // loop over the mesh cells
  for (i=0 ; i<N_Cells ; i++ )
  {
    cell = coll->GetCell(i);
    // number of edges
    N_Edges = cell->GetN_Edges();
    // compute left lower vertex
    x_low = 1000;
    y_low = 1000;
    for (j=0 ; j<N_Edges ; j++ )
    {
      // get x coordinate
      sx = cell->GetVertex(j)->GetX();
      if ( sx < x_low )
      {
        x_low = sx;
      }
      // get y coordinate
      sy = cell->GetVertex(j)->GetY();
      if ( sy < y_low )
      {
        y_low = sy;
      }
    }
    // the coordinates of the left lower vertex are computed
    // compare them with the coordinates of the mesh cells
    for (j=0 ; j<(N_x+1)*(N_y+1) ; j++)
    {
      // check if both coordinates are the same
      if ( ( fabs(x_coord[j]-x_low) < 1e-8 ) && ( fabs(y_coord[j]-y_low) < 1e-8 ) )
      {
        alpha = (int)(j/(N_x+1)+1e-6);
        beta = j%(N_x+1);
        // index is the cell number
        correspond_2dgrid[alpha*N_x+beta] = i;
      }
    }
  }
}


/****************************************************************************************
 *                                                                                       *
 * this routine computes the entries for the column and the row pointer                  *
 * input:  Nodes, N_x, N_y, x_max, x_min, y_max, y_min, z_max, z_min                     *
 * output: *row_ptr, *col_ptr                                                            *
 *                                                                                       *
 ****************************************************************************************/

void filling_row_and_col_ptr(int *N_Entries, int Nodes, int N_x, int N_y, double x_max, double x_min,
double y_max, double y_min, double z_max, double z_min,
double *x_coord, double *y_coord, double *z_coord,
int *row_ptr, int *col_ptr)
{
  int left, right, front, back, bottom, top;
  int N_Neigh, loc_index, index;
  int neigh[27];

  row_ptr[0] = 0;
  *N_Entries = 0;

  // loop over Nodes
  for ( int i=0 ; i<Nodes ; i++ )
  {
    // find all neighbors -> number of entries
    // the if queries check the membership of a node with the plains
    left = right = front = back = bottom = top = 0;
    if ( fabs(x_coord[i]-x_min) < 1e-10 )
    {
      left++;
    }
    if ( fabs(x_coord[i]-x_max) < 1e-10 )
    {
      right++;
    }
    if ( fabs(y_coord[i]-y_min) < 1e-10 )
    {
      front++;
    }
    if ( fabs(y_coord[i]-y_max) < 1e-10 )
    {
      back++;
    }
    if ( fabs(z_coord[i]-z_min) < 1e-10 )
    {
      bottom++;
    }
    if ( fabs(z_coord[i]-z_max) < 1e-10 )
    {
      top++;
    }

    switch ( left + right + front + back + bottom + top )
    {
      case 0:
        N_Neigh = 27;            // interior node
        break;
      case 1:
        N_Neigh = 18;            // plain
        break;
      case 2:
        N_Neigh = 12;            // edge
        break;
      case 3:
        N_Neigh = 8;             // vertex
        break;
      default:
        OutPut("Error in matrix structure !"<<endl);
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
    for ( int j=0 ; j<27 ; j++ )
    {
      neigh[j] = 1;
    }

    // i is not an inner point -> less neighbors -> corresponding neigh_entries will be set to zero
    if ( N_Neigh!=27 )
    {
      // on the left boundary, set local left neighbors to zero
      if (left)
      {
        neigh[0] = neigh[3] = neigh[6] = neigh[9] = neigh[12] = neigh[15]
          = neigh[18] = neigh[21] = neigh[24] = 0;
      }
      // on the right boundary, set local right neighbors to zero
      if (right)
      {
        neigh[2] = neigh[5] = neigh[8] = neigh[11] = neigh[14] = neigh[17]
          = neigh[20] = neigh[23] = neigh[26] = 0;
      }
      // on the bottom, set local lower neighbors to zero
      if (bottom)
      {
        neigh[0] = neigh[1] = neigh[2] = neigh[3] = neigh[4] = neigh[5]
          = neigh[6] = neigh[7] = neigh[8] = 0;
      }
      // on the top, set local upper neighbors to zero
      if (top)
      {
        neigh[18] = neigh[19] = neigh[20] = neigh[21] = neigh[22] = neigh[23]
          = neigh[24] = neigh[25] = neigh[26] = 0;
      }
      // on the front side, set local front neighbors to zero
      if (front)
      {
        neigh[0] = neigh[1] = neigh[2] = neigh[9] = neigh[10] = neigh[11]
          = neigh[18] = neigh[19] = neigh[20] = 0;
      }
      // on the back side, set local back neighbors to zero
      if (back)
      {
        neigh[6] = neigh[7] = neigh[8] = neigh[15] = neigh[16] = neigh[17]
          = neigh[24] = neigh[25] = neigh[26] = 0;
      }
    }

    // off-diagonal entries
    for ( int j=0 ; j<27 ; j++ )
    {
      // the node itself
      if ( j==13 )
      {
        continue;
      }
      if ( neigh[j] )
      {
        switch(j)
        {
          case 0: loc_index = i-(N_x+1)*(N_y+1)-(N_x+1)-1;
          break;
          case 1: loc_index = i-(N_x+1)*(N_y+1)-(N_x+1);
          break;
          case 2: loc_index = i-(N_x+1)*(N_y+1)-(N_x+1)+1;
          break;
          case 3: loc_index = i-(N_x+1)*(N_y+1)-1;
          break;
          case 4: loc_index = i-(N_x+1)*(N_y+1);
          break;
          case 5: loc_index = i-(N_x+1)*(N_y+1)+1;
          break;
          case 6: loc_index = i-(N_x+1)*(N_y+1)+N_x;;
          break;
          case 7: loc_index = i-(N_x+1)*(N_y+1)+N_x+1;
          break;
          case 8: loc_index = i-(N_x+1)*(N_y+1)+N_x+2;
          break;
          case 9: loc_index = i-N_x-2;
          break;
          case 10: loc_index = i-(N_x+1);
          break;
          case 11: loc_index = i-N_x;
          break;
          case 12: loc_index = i-1;
          break;
          case 14: loc_index = i+1;
          break;
          case 15: loc_index = i+N_x;
          break;
          case 16: loc_index = i+(N_x+1);
          break;
          case 17: loc_index = i+N_x+2;
          break;
          case 18: loc_index = i+(N_x+1)*(N_y+1)-(N_x+1)-1;
          break;
          case 19: loc_index = i+(N_x+1)*(N_y+1)-(N_x+1);
          break;
          case 20: loc_index = i+(N_x+1)*(N_y+1)-(N_x+1)+1;
          break;
          case 21: loc_index = i+(N_x+1)*(N_y+1)-1 ;
          break;
          case 22: loc_index = i+(N_x+1)*(N_y+1);
          break;
          case 23: loc_index = i+(N_x+1)*(N_y+1)+1;
          break;
          case 24: loc_index = i+(N_x+1)*(N_y+1)+N_x;
          break;
          case 25: loc_index = i+(N_x+1)*(N_y+1)+N_x+1;
          break;
          case 26: loc_index = i+(N_x+1)*(N_y+1)+N_x+2;
          break;
        }
        col_ptr[index] = loc_index;
        index++;
      }
    }
  }
}


/******************************************************************************/
//
// computation of the size of a mesh cell in convection direction
// approximation formula by Tezduyar and Park, CMAME 59, 307 - 325, 1986
//
/******************************************************************************/

double Bulk_mesh_size_in_convection_direction(double hK, double b1, double b2,
double b3, double *x, double *y,
double *z)
{
  int i;
  double sx, sy, sz, a[64], b[64], den, val, norm_b;

  // only for hexahedra
  sx = sy = sz = 0;
  for (i=0;i<8;i++)
  {
    sx += x[i];
    sy += y[i];
    sz += z[i];
  }
  //bary centre
  sx /= 8;
  sy /= 8;
  sz /= 8;
  // initialize rhs
  memset(b,0,64*SizeOfDouble);

  // set matrices for computation of the coefficients
  // of the bilinear function
  for (i=0;i<8;i++)
  {
    a[8*i] = 1;
    a[8*i+1] = x[i];
    a[8*i+2] = y[i];
    a[8*i+3] = z[i];
    a[8*i+4] = x[i]*y[i];
    a[8*i+5] = x[i]*z[i];
    a[8*i+6] = y[i]*z[i];
    a[8*i+7] = x[i]*y[i]*z[i];
    b[9*i] = 1;
  }

  // solve system for the coefficients of the bilinear function
  SolveMultipleSystems(a,b,8,8,8,8);

  // compute numerator
  norm_b = sqrt(b1*b1 + b2*b2 + b3*b3);
  // compute denominator
  den = 0;
  for (i=0;i<8;i++)
  {
    // value of gradient basis fct. in bary centre
    val = b1*(b[8*i+1]+ b[8*i+4]*sy + b[8*i+5]*sz + b[8*i+7]*sy*sz);
    val += b2*(b[8*i+2]+ b[8*i+4]*sx + b[8*i+6]*sz + b[8*i+7]*sx*sz);
    val += b3*(b[8*i+3]+ b[8*i+5]*sx + b[8*i+6]*sy + b[8*i+7]*sx*sy);
    den += fabs(val);
  }
  // return the mesh size in convection direction
  //OutPut(b1 << " " << b2 << " " << fabs(den) << " " << 2*norm_b/fabs(den) << " " );
  if (den<1e-10)
    return(hK);
  else
    return(2*norm_b/den);
}


/****************************************************************************************
 *                                                                                       *
 *  assembling of matrix with Q_1 finite elements, supg,                                 *
 *  solution of arising linear system                                                    *
 *  OutPut: sol (contains the solution)                                                  *
 *                                                                                       *
 ****************************************************************************************/

void Build_3D_FEM_Matrix_Q1(TCollection *coll,
                            TFEFunction2D *velocity1, TFEFunction2D *velocity2,
                            TFEFunction2D *concent_C,
                            double *sol, double *oldsol,
                            double *lump_mass_PSD, double *matrix_D_Entries_PSD,
                            int *correspond_2dgrid,
                            int N_x, int N_y, int N_z,
                            double *x_coord, double *y_coord, double *z_coord,
                            TSquareMatrix2D *mat, TSquareMatrix2D *matM)
{

    OutPut("Build_3D_FEM_Matrix_Q1 has been removed !!!"<< endl);
    OutPut("The latest version is in ~MooNMD/src/FE/OLD/Bulk2d3d.081008.C !!!"<< endl);
    exit(4711);
}

/******************************************************************************/
//
//  assembling of mass matrix for Q_1 finite elements          
//  array with indices for assembling is filled too
//
/******************************************************************************/
void Compute_Neum_To_Diri_FEM_FCT(int N_x, int N_y, int N_z,
					 double *x_coord, double *y_coord, 
					 double *z_coord,
					 int &N_neum_to_diri, 
					 int* &neum_to_diri,
					 double* &neum_to_diri_x,
					 double* &neum_to_diri_y,
					 double* &neum_to_diri_z)
{
    int i, Nodes, range, count = 0, val;
    double yq;

    Nodes = (N_x+1)*(N_y+1)*(N_z+1);
    
    // inflow from the velocity field
    for (i=(N_x+1)*(N_y+1);i<Nodes;i++)
    {
	val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i]);
	if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)>1e-6)
	{
	    if (val)
		count++;
	}
	else
	{
	    if ((val)||(fabs(1-z_coord[i])<1e-6))
		count++;
	}
    }

    // total number of Dirichlet nodes
    count += (N_x+1)*(N_y+1);
    neum_to_diri = new int[count];
    neum_to_diri_x = new double[3*count];
    neum_to_diri_y = neum_to_diri_x + count;
    neum_to_diri_z = neum_to_diri_y + count;    
    N_neum_to_diri = count;    
    
    // bottom
    for ( i=0 ; i< (N_x+1)*(N_y+1) ; i++ )
    {
	neum_to_diri[i] = i;
	neum_to_diri_x[i] = x_coord[i];
	neum_to_diri_y[i] = y_coord[i];
	neum_to_diri_z[i] = z_coord[i];
    }
    // velocity inlets
    count = (N_x+1)*(N_y+1);
    for (i=(N_x+1)*(N_y+1);i<Nodes;i++)
    {
	if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)>1e-6)
	{
	    val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i]);
	    if (val)
	    {
		neum_to_diri[count] = i;
		neum_to_diri_x[count] = x_coord[i];
		neum_to_diri_y[count] = y_coord[i];
		neum_to_diri_z[count] = z_coord[i];
		count++;
	    }
	}
	else
	// academic test example
	{
	    /*if ((fabs(x_coord[i])<1e-6)||(fabs(y_coord[i])<1e-6)||
		(fabs(1-x_coord[i])<1e-6)||(fabs(1-y_coord[i])<1e-6))
	    {
		neum_to_diri[count] = i;
		neum_to_diri_x[count] = x_coord[i];
		neum_to_diri_y[count] = y_coord[i];
		neum_to_diri_z[count] = z_coord[i];
		count++;
		}*/
	    val = PSD_bound_cound_from_velo_inflow(x_coord[i], y_coord[i]);
	    if ((val)||(fabs(1-z_coord[i])<1e-6))
	    {
		neum_to_diri[count] = i;
		neum_to_diri_x[count] = x_coord[i];
		neum_to_diri_y[count] = y_coord[i];
		neum_to_diri_z[count] = z_coord[i];
		count++;
	    }
	}
    }
}

/****************************************************************************************
 *                                                                                      *
 *  assembling of mass matrix for Q_1 finite elements                                   *
 *  array with indices for assembling is filled too                                     *
 *                                                                                      *
 ***************************************************************************************/

void Build_3D_FEM_FCT_MassMatrix_Q1(TCollection *coll,
				    int N_x, int N_y, int N_z,
				    double *x_coord, double *y_coord, double *z_coord,
				    int* &index_test_ansatz, 
//				    int* &index_test_ansatz_diag, 
				    TSquareMatrix2D *matM,
				    double *lump_mass_PSD)
{
  int locdof[8], *col_ptr, *row_ptr;
  int i, j, k, iq, ii, jj, z, z1, N_Entries, Nodes, ii8;
  int test_index, ansatz_index, index, z_local, z_iq, z_ii;
  int quad_points = 8, diag_index, found, N_cells, range, index1;

  double a[64], b[64], u_val[4], val_test[4], val_ansatz[4], C_val[3];
  double x_coord_loc[8], y_coord_loc[8], z_coord_loc[8];
  double *entriesM;
  double x_max, y_max, z_max, z_min, smag;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, val, weight_det;
  double t1, t2;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;

  static double weight[8]={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  static double qx[8]=
  {
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489
  };
  static double qy[8]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489
  };
  static double qz[8]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489
  };

  t1 = GetTime();

  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();

  N_Entries = matM->GetN_Entries();
  col_ptr = matM->GetKCol();
  row_ptr = matM->GetRowPtr();
  entriesM = matM->GetEntries();

  Nodes = (N_x+1)*(N_y+1)*(N_z+1);
  x_max = x_coord[Nodes-1];
  y_max = y_coord[Nodes-1];
  z_max = z_coord[Nodes-1];
  z_min = z_coord[0];

  if (fabs(z_min-d_p_min) > 1e-8)
  {
    OutPut("z_min " << z_min << " does not correspond to d_p_min " << d_p_min << endl);
    exit(4711);
  }

  // allocate array for storing indices for assembling
  ii = N_cells * N_z * 64;
  index_test_ansatz = new int[ii];
  OutPut("index_test_ansatz " << ii);
//  index_test_ansatz_diag = new int[Nodes];
//  OutPut(" index_test_ansatz_diag " << Nodes);
  
  memset(entriesM,0, N_Entries*SizeOfDouble);

  z = 0;
  //t2 = GetTime();
  //OutPut("time bulkmass (1) " << t2-t1 << endl);
  // loop over Nodes
  for ( i=0 ; i<Nodes ; i++ )
  {
    // assign nodes with mesh cell if it is not at right, back or top
    if ( (fabs(x_coord[i] - x_max)< 1e-10) || (fabs(y_coord[i] - y_max)< 1e-10)
      || (fabs(z_coord[i] - z_max)< 1e-10) )
    {
      continue;
    }

    // consider the mesh cell with node i on the left lower corner
    // compute dof on these mesh cell
    locdof[0] = i;
    locdof[1] = i+1;
    locdof[2] = locdof[1] + N_x+1;
    locdof[3] = locdof[0] + N_x+1;
    locdof[4] = locdof[0] + (N_x+1)*(N_y+1);
    locdof[5] = locdof[4] + 1;
    locdof[6] = locdof[5] + N_x+1;
    locdof[7] = locdof[4] + N_x+1;

    // volume of the hexahedron
    area = (x_coord[i+1]-x_coord[i])*(y_coord[i+N_x+1]-y_coord[i])*(z_coord[(N_x+1)*(N_y+1)+i]-z_coord[i]);
    detJK = area/8.0;

    // compute basis functions
    // set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<8 ; j++ )
    {
      index = locdof[j];
      x_coord_loc[j] = x_coord[index];
      y_coord_loc[j] = y_coord[index];
      z_coord_loc[j] = z_coord[index];
    }

    // compute basis functions
    // set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<8 ; j++ )
    {
      a[8*j] = 1;
      a[8*j+1] = x_coord_loc[j];
      a[8*j+2] = y_coord_loc[j];
      a[8*j+3] = z_coord_loc[j];
      a[8*j+4] = x_coord_loc[j]*y_coord_loc[j];
      a[8*j+5] = x_coord_loc[j]*z_coord_loc[j];
      a[8*j+6] = y_coord_loc[j]*z_coord_loc[j];
      a[8*j+7] = x_coord_loc[j]*y_coord_loc[j]*z_coord_loc[j];
    }
    // initialize rhs
    memset(b,0,64*SizeOfDouble);
    for ( j=0 ; j<8 ; j++ )
	b[9*j] = 1;
    // solve system for the coefficients of the bilinear function
    // already checked !!!
    // solution is stored in b, row-wise
    SolveMultipleSystemsLapack(a,b,8,8,8,8);

    // assemble matrix entries
    // first index for array of indices
    z_iq = z;
    // loop over the quadrature points
    for (iq = 0;iq < quad_points; iq++)
    {
      // quadrature points -> ONLY FOR PARALLELEPIPED !!!
      index = locdof[0];
      index1 = locdof[1];
      xq = x_coord[index] + ( x_coord[index1] - x_coord[index])*(1+ qx[iq])/2;
      index1 = locdof[3];
      yq = y_coord[index] + ( y_coord[index1] - y_coord[index])*(1+ qy[iq])/2;
      index1 = locdof[5];
      zq = z_coord[index] + ( z_coord[index1] - z_coord[index])*(1+ qz[iq])/2;
      weight_det = detJK * weight[iq];

      // loop for test function
      // ii -- test function
      z_ii = 0;
      for ( ii=0 ; ii<8 ; ii++ )
      {
	test_index = locdof[ii];
        // values for test function
        Compute_Q1_Value(b+8*ii, xq, yq, zq, val_test);
	val_test[0] *= weight_det;
	// loop for ansatz functions
        // jj -- ansatz function
        for ( jj=0 ; jj<8 ; jj++ )
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

	  if (test_index!=ansatz_index)
	  {	  
	      Compute_Q1_Value_Gradient(b+8*jj, xq, yq, zq, val_ansatz);
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
    }                            // end quad points
  }                              // end i

  OutPut(" z " << z << endl);
  LumpMassMatrixToVector((TSquareMatrix2D*) matM, lump_mass_PSD);  
  t2 = GetTime();
  OutPut("time bulkmass (2) " << t2-t1 << endl);
}

void Build_3D_FEM_FCT_Matrix_Q1(TCollection *coll,
				TFEFunction2D *velocity1, TFEFunction2D *velocity2,
				TFEFunction2D *concent_C,
				double *sol, double *oldsol,
				double *lump_mass_PSD, double *matrix_D_Entries_PSD,
				int *correspond_2dgrid,
				int N_x, int N_y, int N_z,
				double *x_coord, double *y_coord, double *z_coord,
				TSquareMatrix2D *mat, TSquareMatrix2D *matM_cons,
				TSquareMatrix2D *matM, 
				int *index_test_ansatz,
//				int *index_test_ansatz_diag,
				int N_neum_to_diri, 
				int *neum_to_diri,
				double *neum_to_diri_x,
				double *neum_to_diri_y,
				double *neum_to_diri_z)
{
  int locdof[8], *col_ptr, *row_ptr;
  int i, j, k, iq, ii, jj, z, z1, N_Entries, Nodes, N2, ii8;
  int test_index, ansatz_index, index, alpha, beta, no_of_2dcell;
  int count = 0, topdiri = 0, count_local, count_iq, count_ii;
  int quad_points = 8, diag_index, found, N_cells, range, index1;
  int SC_LDS =  TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  int *test_cells;

  double a[64], b[64], u_val[4], val_test[4], val_ansatz[4], C_val[3], val_sol[4];
  double x_coord_loc[8], y_coord_loc[8], z_coord_loc[8], sol_loc[8], b_sol[8];
  double *entries, *entriesM, *oldrhs_fem_fct0, *rhs, *RhsArray, *tilde_u, *u1, *u2, *G;
  double *bdr_val, *bdr_val_inter, *entriesM_cons;
  double x_max, y_max, z_max, z_min, smag;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, val, weight_det;
  double B_c_C, maxsol, norm_b, al, react;
  //double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  double time = TDatabase::TimeDB->CURRENTTIME, t1, t2, t3;

  static double weight[8]={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  static double qx[8]=
  {
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489
  };
  static double qy[8]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489
  };
  static double qz[8]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489
  };

  //t1 = GetTime();
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
  double factor_G, G_c_C_val, c_C_infty_sat0, val00, c_C;

  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();
  test_cells = new int[N_cells];
  memset(test_cells, 0, N_cells*SizeOfInt);

  // \tilde G(\tilde c_C)   
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);
  // academic test example
  double oldtime;
  if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
  {
      //OutPut("f_infty " << f_infty << " G " << factor_G << endl);
      factor_G = 0.02;
      c_C_infty_sat=0;
      oldtime = time - TDatabase::TimeDB->TIMESTEPLENGTH;
  }

  entries = mat->GetEntries();
  N_Entries = mat->GetN_Entries();
  memset(entries,0, N_Entries*SizeOfDouble);
 //OutPut(N_Entries << endl);
  col_ptr = mat->GetKCol();
  row_ptr = mat->GetRowPtr();

  entriesM = matM->GetEntries();
  entriesM_cons = matM_cons->GetEntries();

  // copy entries of mass matrix
  memcpy(entriesM, entriesM_cons, N_Entries*SizeOfDouble);  

  N2 = (N_x+1)*(N_y+1);
  Nodes = N2*(N_z+1);
  x_max = x_coord[Nodes-1];
  y_max = y_coord[Nodes-1];
  z_max = z_coord[Nodes-1];
  z_min = z_coord[0];

  // N_x*N_y - number of mesh cells for flow domain
  // need in each mesh cell 8 values (8 quad points) for 3 functions
  u1 = new double[N_x*N_y*24];
  memset(u1,0,(N_x*N_y*24)*SizeOfDouble);
  u2 = u1 + N_x*N_y*8;
  G = u2 + N_x*N_y*8;
  
  rhs = new double[4*Nodes];
  memset(rhs,0,4*Nodes*SizeOfDouble);
  // vectors for FEM_FCT_ForConvDiff
  tilde_u  = rhs +  Nodes;
  RhsArray = tilde_u + Nodes;
  oldrhs_fem_fct0 = RhsArray + Nodes;

  OutPut("N_neum_to_diri " << N_neum_to_diri << endl);
  bdr_val = new double[N_neum_to_diri];
  memset(bdr_val,0,N_neum_to_diri*SizeOfDouble);

  z = z1 = 0;
  c_C_infty_sat0 = c_C_infty_sat/c_C_infty;	    
  //t2 = GetTime();
  //OutPut("time bulkfct (1) " << t2-t1 << endl);
  // loop over Nodes
  for ( i=0 ; i<Nodes ; i++ )
  {
    // assign nodes with mesh cell if it is not at right, back or top
    if ( (fabs(x_coord[i] - x_max)< 1e-10) || (fabs(y_coord[i] - y_max)< 1e-10)
      || (fabs(z_coord[i] - z_max)< 1e-10) )
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

    // volume of the hexahedron
    area = (x_coord[i+1]-x_coord[i])*(y_coord[i+N_x+1]-y_coord[i])*(z_coord[(N_x+1)*(N_y+1)+i]-z_coord[i]);
    detJK = area/8.0;

    // compute basis functions
    // set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<8 ; j++ )
    {
      index = locdof[j];
      x_coord_loc[j] = x_coord[index];
      y_coord_loc[j] = y_coord[index];
      z_coord_loc[j] = z_coord[index];
      //sol_loc[j] = sol[index];
    }
    // compute basis functions
    // set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<8 ; j++ )
    {
      a[8*j] = 1;
      a[8*j+1] = x_coord_loc[j];
      a[8*j+2] = y_coord_loc[j];
      a[8*j+3] = z_coord_loc[j];
      a[8*j+4] = x_coord_loc[j]*y_coord_loc[j];
      a[8*j+5] = x_coord_loc[j]*z_coord_loc[j];
      a[8*j+6] = y_coord_loc[j]*z_coord_loc[j];
      a[8*j+7] = x_coord_loc[j]* a[8*j+6]; //y_coord_loc[j]*z_coord_loc[j];
    }
    // initialize rhs
    memset(b,0,64*SizeOfDouble);
    for ( j=0 ; j<8 ; j++ )
	b[9*j] = 1;
    // solve system for the coefficients of the bilinear function
    // already checked !!!
    // solution is stored in b, row-wise
    SolveMultipleSystemsLapack(a,b,8,8,8,8);

    // find corresponding cell in 2d grid
    if ( i < N2 )
    {
      alpha = (int)(i/(N_x+1)+1e-6);
      beta = i%(N_x+1);
      no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];
      /*  if (no_of_2dcell >= N_cells)
      {
        OutPut("number of cell " << no_of_2dcell << " too large " << N_cells << endl);
        exit(4711);
      }
      if (test_cells[no_of_2dcell])
      {
        OutPut("cell " << no_of_2dcell << " treated a second time" << endl);
        exit(4711);
	}*/
      test_cells[no_of_2dcell]++;
    }

    // assemble matrix entries
    // first index for array of indices
    count_iq = count;
    // loop over the quadrature points
    for (iq = 0;iq < quad_points; iq++)
    {
      // quadrature points -> ONLY FOR PARALLELEPIPED !!!
      index = locdof[0];
      index1 = locdof[1];
      //xq = x_coord[index] + ( x_coord[index1] - x_coord[index])*(1+ qx[iq])/2;
      xq = (x_coord[index] +  x_coord[index1] + qx[iq]* (x_coord[index1] - x_coord[index]))/2;
      index1 = locdof[3];
      //yq = y_coord[index] + ( y_coord[index1] - y_coord[index])*(1+ qy[iq])/2;
      yq = (y_coord[index] +  y_coord[index1] + qy[iq] * (y_coord[index1]  - y_coord[index]))/2;
      index1 = locdof[5];
      //zq = z_coord[index] + ( z_coord[index1] - z_coord[index])*(1+ qz[iq])/2;
      zq = (z_coord[index] +  z_coord[index1] + qz[iq] * (z_coord[index1] - z_coord[index]))/2;
      weight_det = detJK * weight[iq];

      // bottom, fill vectors u1, u2, G
      // pointer z is increased at the end of the loop
      if ( i < N2 )
      {
        // compute parameters in quadrature points
        // u1 = u_val[0], u2 = u_val[1]
        velocity1->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,xq,yq,u_val);
        velocity2->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,xq,yq,u_val+1);

        // concentration of species C, C = C_val[0]
        concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,xq,yq,C_val);
	
        // G(c_C)
        G_c_C_val = factor_G*(C_val[0]-c_C_infty_sat0);

	// academic test example
	if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
	{
	    u_val[0] = 2*(2*y_coord[i]-1)*(1-(2*x_coord[i]-1)*(2*x_coord[i]-1)) 
		* (sin(TDatabase::ParamDB->P1*Pi*time)+1)*TDatabase::ParamDB->P9;
	    u_val[1] = -2*(2*x_coord[i]-1)*(1-(2*y_coord[i]-1)*(2*y_coord[i]-1)) 
		*(sin(TDatabase::ParamDB->P1*Pi*time)+1)*TDatabase::ParamDB->P9;
	}
	
        u1[z1] = u_val[0];
        u2[z1] = u_val[1];
        G[z1] = G_c_C_val;
        z1++;
      }

      count_ii = 0;
      // loop for test function
      // ii -- test function
      for ( ii=0 ; ii<8 ; ii++ )
      {
        // values for test function
        Compute_Q1_Value(b+8*ii, xq, yq, zq, val_test);
	val_test[0] *= weight_det;
	// academic test example
	if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
	{
	    test_index = locdof[ii];
	    // f_t
	    val00 = 1e3*(2.00-zq*zq*zq-zq)
		*(sin(Pi*xq)*sin(Pi*yq)+1)
		*cos(Pi*time/2)*Pi/2.0;
	    // + u \cdot \nabla f
	    val00 += u1[z] * 1e3*(2.00-zq*zq*zq-zq)
		*(Pi*cos(Pi*xq)*sin(Pi*yq))
		*sin(Pi*time/2);
	    val00 += u2[z] * 1e3*(2.00-zq*zq*zq-zq)
		*(Pi*sin(Pi*xq)*cos(Pi*yq))
		*sin(Pi*time/2);
	    // -G f_p
	    c_C = 0.6*(sin(2*Pi*xq)*sin(2*Pi*yq)+1)*sin(Pi*time/2)*time;
	    val00 += factor_G * c_C*1e3*(-3*zq*zq-1)
		*(sin(Pi*xq)*sin(Pi*yq)+1)
		*sin(Pi*time/2);
	    RhsArray[test_index] += val00*val_test[0];

	    // f_t
	    val00 = 1e3*(2.00-zq*zq*zq-zq)
		*(sin(Pi*xq)*sin(Pi*yq)+1)
		*cos(Pi*oldtime/2)*Pi/2.0;
	    // + u \cdot \nabla f
	    val00 += u1[z] * 1e3*(2.00-zq*zq*zq-zq)
		*(Pi*cos(Pi*xq)*sin(Pi*yq))
		*sin(Pi*oldtime/2);
	    val00 += u2[z] * 1e3*(2.00-zq*zq*zq-zq)
		*(Pi*sin(Pi*xq)*cos(Pi*yq))
		*sin(Pi*oldtime/2);
	    // -G f_p
	    c_C = 0.6*(sin(2*Pi*xq)*sin(2*Pi*yq)+1)*sin(Pi*oldtime/2)*oldtime;
	    val00 += factor_G* c_C*1e3*(-3*zq*zq-1)
		*(sin(Pi*xq)*sin(Pi*yq)+1)
		*sin(Pi*oldtime/2);
	    oldrhs_fem_fct0[test_index] += val00*val_test[0];
	}
        // loop for ansatz functions
        // jj -- ansatz function
        for ( jj=0 ; jj<8 ; jj++ )
        {
	    //ansatz_index = locdof[jj];
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
	
	//index = test_index;
	// values for ansatz function
	Compute_Q1_Gradient(b+8*jj, xq, yq, zq, val_ansatz);
	entries[index] += (u1[z]*val_ansatz[1]+u2[z]*val_ansatz[2]+G[z]*val_ansatz[3])*val_test[0];
	
	/* // mass term
	   val = val_ansatz[0]*val_test[0]; 
	  // add to M(locdof[ii], locdof[jj])
	  entriesM[index] += val;
	  */
	  /***********
	  u1[z] = 0.1;
	  u2[z] = -0.05;
	  G[z] = 1;*/
	  // convective term
          //entries[index] += (u1[z]*val_ansatz[1]+u2[z]*val_ansatz[2]+G[z]*val_ansatz[3])*val_test[0];
          //val *= weight_det;
          // add to A(locdof[ii], locdof[jj])
          //entries[index] += val;
        }
      }
      z++;
    }                            // end quad points

    if ( z == N_x*N_y*8 )
    {
      z = 0;
    }
  }                              // end i
 
  //t2 = GetTime();
  //OutPut("time bulkfct (2) " << t2-t1 << endl);

  // FEM--FCT
  // set Dirichlet boundary conditions, on all (possible) inflow boundaries
  for ( i=0 ; i<(N_x+1)*(N_y+1) ; i++ )
  {
      ii = i;
      // treat right d.o.f. seperately
      if ( ((ii+1)%(N_x+1)==0) )
	  ii = ii-1;
      // treat upper d.o.f. seperately
      if ( ii>=N_y*(N_x+1) )
	  ii = ii-(N_x+1);
      // right corner
      //if ( i==N2-1 )
      //  ii = i-(N_x+1)-1;
      
      alpha = (int)(ii/(N_x+1)+1e-6);
      beta = ii%(N_x+1);
      no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];
      concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],C_val);
      //G_c_C_val = k_g*(C_val[0]-c_C_infty_sat*exp(C_2/z_min));
      G_c_C_val = k_g*(c_C_infty*C_val[0]-c_C_infty_sat);
      // compute G*n, n=(0,0,-1);
      if (G_c_C_val*f_infty > 1e-10)
      {
	  // compute rate of nucleation
	  B_c_C = k_nuc*pow(c_C_infty*(C_val[0] - 1),5);
	  // truncate negative values
	  if (B_c_C < 0)
	      B_c_C = 0;
	  // compute new particle size distribution
	  bdr_val[i] = B_c_C/ (G_c_C_val*f_infty);
      }
      else
      {
	  // not academic test example
	  if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)>1e-6)
	  {
	      // top Dirichlet condition, bdr_val[i] is already set to be zero
	      if (G_c_C_val*f_infty < 0)
	      {
		  j = i + N_z * (N_x+1)*(N_y+1);
		  neum_to_diri[i] = j;
		  topdiri = 1;
	      }
	  }
      }
  }  // the other Dirichlet values are already set to zero


  // academic test example
  if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
  {
      for ( i=0 ; i<N_neum_to_diri ; i++ )
      {
	  j = neum_to_diri[i];
	  bdr_val[i] = 1e3*(2.00-z_coord[j]*z_coord[j]*z_coord[j]-z_coord[j])
	      *(sin(Pi*x_coord[j])*sin(Pi*y_coord[j])+1)
	      *sin(Pi*(time/2-TDatabase::TimeDB->TIMESTEPLENGTH/4.0));
      }
  }
  
  // this sets the array rhs
  FEM_FCT_ForConvDiff((TSquareMatrix2D*) matM, (TSquareMatrix2D*) mat,
		      Nodes, Nodes,
		      lump_mass_PSD, matrix_D_Entries_PSD,
		      sol, oldsol,
		      rhs, RhsArray, oldrhs_fem_fct0, tilde_u,
		      N_neum_to_diri, neum_to_diri,
		      NULL,NULL,
		      1, NULL,bdr_val);

  //t2 = GetTime();
  //OutPut("time bulkfct (3) " << t2-t1 << endl);
  // build matrix for FEM-FCT
  matM->Reset();
  FEM_FCT_SystemMatrix(matM, mat, lump_mass_PSD, Nodes);
  OutPut("entries " << Ddot(matM->GetN_Entries(),matM->GetEntries(),matM->GetEntries()) << 
	 " lump " << Ddot(Nodes,lump_mass_PSD,lump_mass_PSD) << endl);

  // academic test example
  if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)<1e-6)
  {
      for ( i=0 ; i<N_neum_to_diri ; i++ )
      {
	  j = neum_to_diri[i];
	  bdr_val[i] = 1e3*(2.00-z_coord[j]*z_coord[j]*z_coord[j]-z_coord[j])
	      *(sin(Pi*x_coord[j])*sin(Pi*y_coord[j])+1)
	      *sin(Pi*time/2);
      }
  }
  //t2 = GetTime();
  //OutPut("time bulkfct (4) " << t2-t1 << endl);
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


  //OutPut(endl);
  //t2 = GetTime();
  //OutPut("time bulkfct (5) " << t2-t1 << endl);
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

  OutPut(TDatabase::TimeDB->CURRENTTIME << " Solver done " << sqrt(Ddot(Nodes,sol,sol)) << " max " << maxsol << endl);

  // reset indices if necessary
  if (topdiri)
  {
      for (i=0;i<N2;i++)
	  neum_to_diri[i] = i;  
  }

  //for (i=0;i<N2;i++)
  //   OutPut(sol[i] << " " <<  bdr_val[i] << endl);
  //t2 = GetTime();
  //OutPut("time bulkfct (6) " << t2-t1 << endl);

  // deletion of the arrays
  delete bdr_val;
  delete rhs;
  delete u1;
  delete test_cells;
}


/****************************************************************************************
 *                                                                                      *
 *  assembling of matricess for Q_1 group finite elements                               *
 *  array with indices for assembling is filled too                                     *
 *                                                                                      *
 ***************************************************************************************/

void Build_3D_FEM_FCT_Matrices_Q1_GroupFEM_Bulk(TCollection *coll,
            int N_x, int N_y, int N_z,
            double *x_coord, double *y_coord, double *z_coord,
            TSquareMatrix2D *matM,TSquareMatrix2D *matU1, 
                  TSquareMatrix2D *matU2, TSquareMatrix2D *matG,
            double *lump_mass_PSD)
{
  int locdof[8], *col_ptr, *row_ptr;
  int i, j, k, iq, ii, jj, z, z1, N_Entries, Nodes, ii8;
  int test_index, ansatz_index, index, z_local, z_iq, z_ii;
  int quad_points = 8, diag_index, found, N_cells, range, index1, index_test_ansatz_loc[64];

  double a[64], b[64], u_val[4], val_test[4], val_ansatz[4], C_val[3];
  double x_coord_loc[8], y_coord_loc[8], z_coord_loc[8];
  double *entriesM, *entriesU1, *entriesU2, *entriesG, *growth;
  double x_max, y_max, z_max, z_min, smag;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, val[4], weight_det;
  double t1, t2;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;

  static double weight[8]={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  static double qx[8]=
  {
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489,  0.5773502691896257645091489
  };
  static double qy[8]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489
  };
  static double qz[8]=
  {
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    -0.5773502691896257645091489, -0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489,
    0.5773502691896257645091489,  0.5773502691896257645091489
  };

  t1 = GetTime();

  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();

  N_Entries = matM->GetN_Entries();
  col_ptr = matM->GetKCol();
  row_ptr = matM->GetRowPtr();
  entriesM = matM->GetEntries();    
  memset(entriesM,0, N_Entries*SizeOfDouble);
 entriesU1 = matU1->GetEntries();
  memset(entriesU1,0, N_Entries*SizeOfDouble);
  entriesU2 = matU2->GetEntries();
  memset(entriesU2,0, N_Entries*SizeOfDouble);
   entriesG = matG->GetEntries();
  memset(entriesG,0, N_Entries*SizeOfDouble);

  Nodes = (N_x+1)*(N_y+1)*(N_z+1);
  x_max = x_coord[Nodes-1];
  y_max = y_coord[Nodes-1];
  z_max = z_coord[Nodes-1];
  z_min = z_coord[0];

  if (fabs(z_min-d_p_min) > 1e-8)
  {
    OutPut("z_min " << z_min << " does not correspond to d_p_min " << d_p_min << endl);
    exit(4711);
  }

  z = 0;
  //t2 = GetTime();
  //OutPut("time bulkmass (1) " << t2-t1 << endl);
  // loop over Nodes
  for ( i=0 ; i<Nodes ; i++ )
  {
    // assign nodes with mesh cell if it is not at right, back or top
    if ( (fabs(x_coord[i] - x_max)< 1e-10) || (fabs(y_coord[i] - y_max)< 1e-10)
      || (fabs(z_coord[i] - z_max)< 1e-10) )
    {
      continue;
    }

    // consider the mesh cell with node i on the left lower corner
    // compute dof on these mesh cell
    locdof[0] = i;
    locdof[1] = i+1;
    locdof[2] = locdof[1] + N_x+1;
    locdof[3] = locdof[0] + N_x+1;
    locdof[4] = locdof[0] + (N_x+1)*(N_y+1);
    locdof[5] = locdof[4] + 1;
    locdof[6] = locdof[5] + N_x+1;
    locdof[7] = locdof[4] + N_x+1;

    // volume of the hexahedron
    area = (x_coord[i+1]-x_coord[i])*(y_coord[i+N_x+1]-y_coord[i])*(z_coord[(N_x+1)*(N_y+1)+i]-z_coord[i]);
    detJK = area/8.0;

    // compute basis functions
    // set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<8 ; j++ )
    {
      index = locdof[j];
      x_coord_loc[j] = x_coord[index];
      y_coord_loc[j] = y_coord[index];
      z_coord_loc[j] = z_coord[index];
    }

    // compute basis functions
    // set matrices for computation of the coefficients of the bilinear function
    for ( j=0 ; j<8 ; j++ )
    {
      a[8*j] = 1;
      a[8*j+1] = x_coord_loc[j];
      a[8*j+2] = y_coord_loc[j];
      a[8*j+3] = z_coord_loc[j];
      a[8*j+4] = x_coord_loc[j]*y_coord_loc[j];
      a[8*j+5] = x_coord_loc[j]*z_coord_loc[j];
      a[8*j+6] = y_coord_loc[j]*z_coord_loc[j];
      a[8*j+7] = x_coord_loc[j]*y_coord_loc[j]*z_coord_loc[j];
    }
    // initialize rhs
    memset(b,0,64*SizeOfDouble);
    for ( j=0 ; j<8 ; j++ )
  b[9*j] = 1;
    // solve system for the coefficients of the bilinear function
    // already checked !!!
    // solution is stored in b, row-wise
    SolveMultipleSystemsLapack(a,b,8,8,8,8);

    // assemble matrix entries
    // first index for array of indices
    z_iq = z;
    // loop over the quadrature points
    for (iq = 0;iq < quad_points; iq++)
    {
      // quadrature points -> ONLY FOR PARALLELEPIPED !!!
      index = locdof[0];
      index1 = locdof[1];
      xq = x_coord[index] + ( x_coord[index1] - x_coord[index])*(1+ qx[iq])/2;
      index1 = locdof[3];
      yq = y_coord[index] + ( y_coord[index1] - y_coord[index])*(1+ qy[iq])/2;
      index1 = locdof[5];
      zq = z_coord[index] + ( z_coord[index1] - z_coord[index])*(1+ qz[iq])/2;
      weight_det = detJK * weight[iq];

      // loop for test function
      // ii -- test function
      z_ii = 0;
      for ( ii=0 ; ii<8 ; ii++ )
      {
  test_index = locdof[ii];
        // values for test function
        Compute_Q1_Value(b+8*ii, xq, yq, zq, val_test);
  val_test[0] *= weight_det;
  // loop for ansatz functions
        // jj -- ansatz function
        for ( jj=0 ; jj<8 ; jj++ )
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
          index_test_ansatz_loc[z_ii] = index;
          z_ii++;
          break;
      }
        }
    }
    else
    {
        // other quad points
           index = index_test_ansatz_loc[z_ii];
            z_ii++;       
    }
         
    Compute_Q1_Value_Gradient(b+8*jj, xq, yq, zq, val_ansatz);
    // mass term
    val[0] = val_ansatz[0]*val_test[0];
    // convective terms
    val[1] = val_ansatz[1]*val_test[0];
    val[2] = val_ansatz[2]*val_test[0];
    val[3] = val_ansatz[3]*val_test[0];

          // add to M(locdof[ii], locdof[jj])
          entriesM[index] += val[0];
          entriesU1[index] += val[1];
          entriesU2[index] += val[2];
          entriesG[index] += val[3];
  }
      }
    }                            // end quad points
  }                              // end i

  OutPut(" z " << z << endl);
  LumpMassMatrixToVector((TSquareMatrix2D*) matM, lump_mass_PSD);  
    
  t2 = GetTime();
  OutPut("time bulkmass (2) " << t2-t1 << endl);
}
/*******************************************************************************/
//
// FEM_FCT_Matrix_Q1_GroupFEM_3D_Bulk
//
// sol            - at beginning: solution of last time
//                  at end: solution of current time
// oldsol         - same as sol
// lump_mass_PSD  - lumped mass matrix, already computed
// matrix_D_Entries_PSD - set in FEM_FCT_ForConvDiff
// mat            - matrix for assembling the convection term
// matM           - mass matrix with incorporation of Dirichlet values
//                  matM_cons is copied to matM in this routine
// matM_cons      - consistent mass matrix, is not changed
//
// sol and oldsol have at the entrance and at leaving the routine the same
//                  values
//
/*******************************************************************************/

void FEM_FCT_Matrix_Q1_GroupFEM_3D_Bulk(TCollection *coll,
TFEFunction2D *velocity1, TFEFunction2D *velocity2,
TFEFunction2D *concent_C,
double *sol, double *oldsol, 
double *lump_mass_PSD, double *matrix_D_Entries_PSD,
int *correspond_2dgrid,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
TSquareMatrix2D *mat,
TSquareMatrix2D *matM_cons,
TSquareMatrix2D *matM,
TSquareMatrix2D *matU1,
TSquareMatrix2D *matU2,
TSquareMatrix2D *matG,
double *psd_coeff,
int N_neum_to_diri,
int *neum_to_diri,
double *neum_to_diri_x,
double *neum_to_diri_y,
double *neum_to_diri_z)
{
  int *col_ptr, *row_ptr;
  int i, j, k, iq, ii, jj, z, N_Entries, Nodes, ij, start, end;
  int test_index, ansatz_index, index, index1, alpha, beta, gamma, no_of_2dcell;
  int diag_index, found, N_cells, range_y, range_z, topdiri = 0;
  int SC_LDS =  TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;

  // N2 and N3 are defined to save multiplications during the computation of loc_dof
  int N2 = (N_x+1)*(N_y+1);
 
  double a[256], b[256], C_val[4],Temp_val[4], u_val[6], val_test[5], val_ansatz[5], val_sol[5];
  double *entries, *entriesM, *rhs, *tilde_u, *u1, *u2,  *growth;
  double *bdr_val, *entriesM_cons, *entriesU1, *entriesU2, *entriesG;
  double *RhsArray, *oldrhs_fem_fct0;
  double x_max, y_max, z_max, a_max, z_min, smag;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, aq, val, weight_det;
  double B_c_C, maxsol, norm_b, al, react, factor_growth;
  double time = TDatabase::TimeDB->CURRENTTIME, t1, t2, t3;

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
  double factor_G, c_C_infty_sat0, G_c_C_val;
  // \tilde G(\tilde c_C)   
  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);
  c_C_infty_sat0 = c_C_infty_sat/c_C_infty;     
  
  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();

  // matrix objects
  entries = mat->GetEntries();
  N_Entries = mat->GetN_Entries();
  memset(entries,0, N_Entries*SizeOfDouble);
  col_ptr = mat->GetKCol();
  row_ptr = mat->GetRowPtr();
  entriesM = matM->GetEntries();
  entriesM_cons = matM_cons->GetEntries();
  entriesU1 = matU1->GetEntries();
  entriesU2 = matU2->GetEntries();
  entriesG = matG->GetEntries();

  // copy entries of mass matrix
  memcpy(entriesM, entriesM_cons, N_Entries*SizeOfDouble);

  Nodes = N2*(N_z+1);
  x_max = x_coord[Nodes-1];
  y_max = y_coord[Nodes-1];
  z_max = z_coord[Nodes-1];
  z_min = z_coord[0];

  // initialization of u1, u2, u3, growth, reaction in the vector u1
  u1 = psd_coeff;
  u2 = u1 + Nodes;
  growth = u2 + Nodes;

  // initialization of rhs and the vectors for FEM_FCT_ForConvDiff in the vector rhs
   rhs = new double[4*Nodes];
  memset(rhs,0,4*Nodes*SizeOfDouble);
  // vectors for FEM_FCT_ForConvDiff
  tilde_u  = rhs +  Nodes;
  RhsArray = tilde_u + Nodes;
  oldrhs_fem_fct0 = RhsArray + Nodes;
 
  bdr_val = new double[N_neum_to_diri];
  memset(bdr_val,0,N_neum_to_diri*SizeOfDouble);

  z = 0;
// loop over the nodes to fill the arrays for the convection
  for ( i=0 ; i<Nodes ; i++ )
  {
    // find corresponding cell in 2d grid
    if ( i < N2 )
    {
      ij = i;
      // treat nodes on boundaries in a special was
      if (fabs(x_coord[i] - x_max)< 1e-6) 
  ij = ij - 1;
      if (fabs(y_coord[i] - y_max)< 1e-6)
  ij = ij - (N_x+1);
     
     alpha = (int)(ij/(N_x+1)+1e-6);
     beta = ij%(N_x+1);
     no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];
     velocity1->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],u_val);
     velocity2->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],u_val+1);
     // store the values
      u1[i] = u_val[0];
      u2[i] = u_val[1];
     // concentration of species C, C = C_val[0]
     concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],C_val); 
     // G(c_C)
      growth[i] = factor_G*(C_val[0]-c_C_infty_sat0);
    }
   else
    {
  // for larger internal coordinates there are the same velocities
  // find corresponding index of samllest internal coordinate
  j = i%N2;
  u1[i] = u1[j];
  u2[i] = u2[j];
  growth[i] = growth[j];
    }
   }
 // compute matrix and rhs
  // working matrix is mat
  // femrhs is array rhs
  // loop over the rows
  for (i=0; i< Nodes; i++)
  {
      start = row_ptr[i];
      end = row_ptr[i+1];
      for (j=start;j<end;j++)
      {
    // matrix
  entries[j] =  entriesU1[j] * u1[i] + entriesU2[j] * u2[i] + entriesG[j] * growth[i];
    }
      // rhs 
      //RhsArray[i] = oldrhs_fem_fct0[i] = 0.0;
   }
  // end i
 
 // FEM--FCT
  // set Dirichlet boundary conditions, on all (possible) inflow boundaries
  for ( i=0 ; i<(N_x+1)*(N_y+1) ; i++ )
  {
      ii = i;
      // treat right d.o.f. seperately
      if ( ((ii+1)%(N_x+1)==0) )
    ii = ii-1;
      // treat upper d.o.f. seperately
      if ( ii>=N_y*(N_x+1) )
    ii = ii-(N_x+1);
      // right corner
      //if ( i==N2-1 )
      //  ii = i-(N_x+1)-1;
      
      alpha = (int)(ii/(N_x+1)+1e-6);
      beta = ii%(N_x+1);
      no_of_2dcell = correspond_2dgrid[alpha*N_x+beta];
      concent_C->FindValueLocal(coll->GetCell(no_of_2dcell),no_of_2dcell,x_coord[i],y_coord[i],C_val);
      //G_c_C_val = k_g*(C_val[0]-c_C_infty_sat*exp(C_2/z_min));
      G_c_C_val = k_g*(c_C_infty*C_val[0]-c_C_infty_sat);
      // compute G*n, n=(0,0,-1);
      if (G_c_C_val*f_infty > 1e-10)
      {
    // compute rate of nucleation
    B_c_C = k_nuc*pow(c_C_infty*(C_val[0] - 1),5);
    // truncate negative values
    if (B_c_C < 0)
        B_c_C = 0;
    // compute new particle size distribution
    bdr_val[i] = B_c_C/ (G_c_C_val*f_infty);
      }
      else
      {
        // top Dirichlet condition, bdr_val[i] is already set to be zero
        if (G_c_C_val*f_infty < 0)
        {
      j = i + N_z * (N_x+1)*(N_y+1);
      neum_to_diri[i] = j;
      topdiri = 1;
        }
      }
  }  // the other Dirichlet values are already set to zero


  // this sets the array rhs
  FEM_FCT_ForConvDiff((TSquareMatrix2D*) matM, (TSquareMatrix2D*) mat,
          Nodes, Nodes,
          lump_mass_PSD, matrix_D_Entries_PSD,
          sol, oldsol,
          rhs, RhsArray, oldrhs_fem_fct0, tilde_u,
          N_neum_to_diri, neum_to_diri,
          NULL,NULL,
          1, NULL,bdr_val);

  //t2 = GetTime();
  //OutPut("time bulkfct (3) " << t2-t1 << endl);
  // build matrix for FEM-FCT
  matM->Reset();
  FEM_FCT_SystemMatrix(matM, mat, lump_mass_PSD, Nodes);
  OutPut("entries " << Ddot(matM->GetN_Entries(),matM->GetEntries(),matM->GetEntries()) << 
   " lump " << Ddot(Nodes,lump_mass_PSD,lump_mass_PSD) << endl);

  //t2 = GetTime();
  //OutPut("time bulkfct (4) " << t2-t1 << endl);
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


  //OutPut(endl);
  //t2 = GetTime();
  //OutPut("time bulkfct (5) " << t2-t1 << endl);
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

  //OutPut(TDatabase::TimeDB->CURRENTTIME << " Solver done " << sqrt(Ddot(Nodes,sol,sol)) << " max " << maxsol << endl);

  // reset indices if necessary
  if (topdiri)
  {
      for (i=0;i<N2;i++)
    neum_to_diri[i] = i;  
  }

  //for (i=0;i<N2;i++)
  //   OutPut(sol[i] << " " <<  bdr_val[i] << endl);
  //t2 = GetTime();
  //OutPut("time bulkfct (6) " << t2-t1 << endl);

  // deletion of the arrays
  delete bdr_val;
  delete rhs;
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

void Integral_For_Particle_Increase_Term(TFESpace2D *fespace, TFEFunction2D *fefct,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_layers_coord,
double *concent_C_array, double *f)
{
  int i, j, N2, i1, i2, m, k, l, N_Cells, local;
  int N_Edges, index;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double val, x, y, value;
  //  double hz, z, val, hx, x, y, value, c_C_infty ;
  double *integral;
  TCollection *Coll;
  TBaseCell *cell;
  int *indextest;

  double c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
  double C_2 = TDatabase::ParamDB->BULK_C_2;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;

  N2 = (N_x+1)*(N_y+1);

  indextest = new int[N2];
  memset(indextest,0,N2*SizeOfInt);

  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  integral = fefct->GetValues();

  for (i=0;i<N2;i++)
  {
    integral[i] = -1;
  }

  // all spaces use same Coll
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  //loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    // get number of edges (or vertices)
    N_Edges=cell->GetN_Edges();
    // get pointer to degrees of freedom connected to this mesh cell
    DOF = GlobalNumbers + BeginIndex[i];
    // loop over all vertices
    for (j=0;j<N_Edges;j++)
    {
      // corners of mesh cell
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      if (cell->GetN_Vertices()!=4)
      {
	  OutPut("Integral_For_Particle_Increase_Term only implemented for quads !!!" << endl);
	  exit(4711);
      }
      // correspondance of local vertex and local dof
      switch(j)
      {
        case 0:
        case 1:
          local = j;
          break;
        case 2:
          local = 3;
          break;
        case 3:
          local = 2;
          break;
      }
      // l - global index
      l = DOF[local];
      // already done
      if (integral[l]>=-0.1)
      {
        continue;
      }
      // corresponding index of array f - comparision of the x and y coordinates
      for (m=0 ; m<(N_x+1) ; m++)
      {
        if ( fabs(x-x_coord[m])<1e-6 )
        {
          index = m;
          break;
        }
      }
      for (m=0 ; m<=N_y*(N_x+1) ; m+=(N_x+1))
      {
        if ( fabs(y-y_coord[m])<1e-6 )
        {
          index += m;
          break;
        }
      }

      indextest[index]++;
      val = 0.0;

      if (TDatabase::ParamDB->BULK_GROWTH_RATE!=2)
      {
        // compute integral in z-direction, composed trapezoidal rule
        for (k=1 ; k<=N_z ; k++)
        {
          i1 = index+k*N2;
          i2 = index+(k-1)*N2;
          val += (z_layers_coord[k]*z_layers_coord[k]*f[i1] + z_layers_coord[k-1]*z_layers_coord[k-1]*f[i2])
            * (z_layers_coord[k]-z_layers_coord[k-1]) / 2;
        }
      }
      else
      {
        // compute concentration of C
        value = concent_C_array[index];
        // compute integral in z-direction, composed trapezoidal rule
        for (k=1 ; k<=N_z ;k++)
        {
          i1 = index+k*N2;
          i2 = index+(k-1)*N2;
          val += (value-c_C_infty_sat/c_C_infty*exp(C_2/(z_layers_coord[k]*d_p_max)))
            * z_layers_coord[k]*z_layers_coord[k]*f[i1];
          val += (value-c_C_infty_sat/c_C_infty*exp(C_2/(z_layers_coord[k-1]*d_p_max)))
            * z_layers_coord[k-1]*z_layers_coord[k-1]*f[i2];
          val *= (z_layers_coord[k]-z_layers_coord[k-1]) / 2;
        }
      }

      /*
            // compute integral in z-direction, composed trapezoidal rule
            // ONLY EQUIDISTANT MESHES !!!
            // following script from 05/08/23
            if (TDatabase::ParamDB->BULK_PB_DISC!=2)
            {
        val = d_p_min*d_p_min*f[index]/2;
        // ONLY equidistant mesh size !!!
        for (k=1;k<Nz-1;k++)
        {
              z = d_p_min + k*hz;
              ii = index+k*N2;
              val += z*z*f[ii];
          }
          // f should be zero at 1, can be deleted
          ii = index+(Nz-1)*N2;
          val += f[ii]/2;
            }
            else
          // G(c_C,d_p)
            {
          // compute concentration of C
          value = concent_C_array[index];
          val = (value-c_C_infty_sat/c_C_infty*exp(C_2/(d_p_min*d_p_max)))*d_p_min*d_p_min*f[index]/2.0;
          // ONLY equidistant mesh size !!!
          for (k=1;k<Nz-1;k++)
          {
              z = d_p_min + k*hz;
              ii = index+k*N2;
              val += (value-c_C_infty_sat/c_C_infty*exp(C_2/(z*d_p_max)))*z*z*f[ii];
          }
          // f should be zero at 1, can be deleted
          ii = index+(Nz-1)*N2;
          val += (value-c_C_infty_sat/c_C_infty*exp(C_2/(d_p_max)))*f[ii]/2;
            }
            val *= hz;
      */

      // insert the value into the fe function
      integral[l] = val;
    }
  }

  for (i=0 ; i<N2 ; i++)
  {
    if (!indextest[i])
    {
      OutPut("Index " << i << " not found " << endl);
      exit(1);
    }
  }

  OutPut(TDatabase::TimeDB->CURRENTTIME << " integral " << sqrt(Ddot(N2,integral,integral)) << endl);

  delete indextest;
}


/****************************************************************************************
 *                                                                                       *
 *  Evaluation of population balance at outflow                                          *
 *                                                                                       *
 ****************************************************************************************/

void Evalute_f_at_outflow(int N_x, int N_y, int N_z, double *x_coord, double *z_layers_coord, 
			  double *f, double *average_median, int *average_step)
{
  int k, i0, i;
  double xout, *size, *number, median;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;

  // compute centre of outflow
  k = (int)TDatabase::ParamDB->P9;
  xout = (k+1)/32.0;
  OutPut("center of outflow " << xout << endl);

  // compute index for value of f on z=0
  i0 = (int) (xout*N_x);

  if (fabs(xout -x_coord[i0])>1e-6)
  {
    OutPut("Wrong index in Evalute_f_at_outflow " << endl);
    exit(1);
  }
  size = new  double[N_z+1];
  number = new  double[N_z+1];
  for (i=0 ; i<(N_z+1) ; i++)
  {
      size[i] = z_layers_coord[i]*d_p_max;
      number[i] =  f[i0+i*(N_x+1)*(N_y+1)]* f_infty;
      // no output for academic test example
      if (fabs(TDatabase::ParamDB->BULK_c_C_infty_sat-1e-4)>1e-6)
      {
	  OutPut("PSD " << TDatabase::TimeDB->CURRENTTIME << " " <<
		 z_layers_coord[i] << " " << f[i0+i*(N_x+1)*(N_y+1)] << " " <<
		 size[i] << " " << number[i] << endl);
      }
  }
  
  median =  calculate_dp_50(N_z+1,size,number);
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
 *  .vtk file to visualize f                                                             *
 *                                                                                       *
 ****************************************************************************************/

void write_vtk_psd(int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
double *f_old, const char *name)
{
  int i;
  int N_Cells = N_x*N_y*N_z;
  int N2 = (N_x+1)*(N_y+1);
  int N3 = N2*(N_z+1);
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;

  FILE* out = fopen(name,"w");
	d_p_max = f_infty = 1;
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
    fprintf(out,"%e\n", z_coord[i]*d_p_max);
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

  for ( i=0 ; i<N3 ; i++ )
  {
    fprintf(out,"%f\n",f_old[i]*f_infty);
  }

  fprintf(out,"\n");
  fclose(out);
}

/****************************************************************************************
 * write population balance                                                             *
 ****************************************************************************************/

void write_psd(int N_x, int N_y, int N_z,
	       double *x_coord, double *y_coord, double *z_layers_coord,
	       double *f_old, const char *name)
{
  int i;
  int N_Cells = N_x*N_y*N_z;
  int N2 = (N_x+1)*(N_y+1);
  int N3 = N2*(N_z+1);

  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double f_infty = TDatabase::ParamDB->BULK_f_infty;

  OutPut("write psd " << name << endl);
  FILE* out = fopen(name,"w");
 
  for ( i=0 ; i<N2 ; i++ )
  {
    //x
    fprintf(out,"%f",x_coord[i]);
    fprintf(out,"%s"," ");
    //y
    fprintf(out,"%f",y_coord[i]);
    fprintf(out,"%s"," ");
/*    //ecken ausschliessen
    if (i==0)
      	{//links unten
         fprintf(out,"%i",-1);
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+(N_x+1));
	}
    if (i==N_x)
      	{
	//rechts unten
         fprintf(out,"%i",-1);
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",i-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+(N_x+1));
	}
    if (i==N_y*(N_x+1))

      	{
         //links oben
         fprintf(out,"%i",i-(N_x+1));
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
	}
    if (i==(N_y+1)*(N_x+1)-1)
      	{ 
          //rechts oben
         fprintf(out,"%i",i-(N_x+1));
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",i-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
	}
    //kanten ausschliessen
    if ( (i%(N_x+1)==0) && (i!=0) && (i!=N_y*(N_x+1)) )
     	{ 
         //kante links
         fprintf(out,"%i",i-(N_x+1));
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+(N_x+1));
        }
        if ((i%(N_x+1)==N_x) && (i!=(N_y+1)*(N_x+1)-1) && (i!=N_y*(N_x+1)))
     	{
          //kante rechts
         fprintf(out,"%i",i-(N_x+1));
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",i-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+(N_x+1));
         }
    if ((i<N_x) && (i!=0))
     	{
        //kante unten
         fprintf(out,"%i",-1);
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",i-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+(N_x+1));
         }
	
    if ((i>N_y*(N_x+1)) && (i!=(N_y+1)*(N_x+1)-1))
     	{
          //kante oben
         fprintf(out,"%i",i-(N_x+1));
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",i-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",-1);
	}
	 else
        //innere knoten
	{
         fprintf(out,"%i",i-(N_x+1));
	 fprintf(out,"%s"," ");
         fprintf(out,"%i",i-1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+1);
         fprintf(out,"%s"," ");
         fprintf(out,"%i",i+(N_x+1));
 	}
*/
    // left neighbour
    if (i%(N_x+1)==0)
        fprintf(out," %i",-1);
    else
	fprintf(out," %i",i-1);
    // lower neighbour
    if (i<N_x+1)
        fprintf(out," %i",-1);
    else
	fprintf(out," %i",i-(N_x+1));
    // right neighbour
    if ((i+1)%(N_x+1)==0)
        fprintf(out," %i",-1);
    else
	fprintf(out," %i",i+1);
    // upper neighbour
    if (i>N_y*(N_x+1)-1)
        fprintf(out," %i",-1);
    else
	fprintf(out," %i",i+(N_x+1));
	
  fprintf(out,"\n");
  }
  
  fprintf(out,"\n");

  for ( i=0 ; i<N_z+1 ; i++ )
  {
   //z_layers
   fprintf(out,"%e\n",z_layers_coord[i]);
   
  }
  fprintf(out,"\n");

  for ( i=0 ; i<N3 ; i++ )
  {
    fprintf(out,"%f\n",f_old[i]*f_infty);
  }

  fclose(out);
}

#endif // #ifdef __2D__
