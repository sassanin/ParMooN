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
 *                         Windtunnel_3d4d.C                                            *
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
#include <Windtunnel_3d4d.h>

#include <Database.h>

#include <string.h>
#include <stdlib.h>
#include <../../Main_Users/Sashi/TNSE_3D/windtunnel_log_normal_parameters.h>


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
 * DROPS_bound_cound_from_velo_inflow                                                   *
 *                                                                                      *
 * computing the boundary conditions for the PSD on the flow inlet (x=0)                *
 *                                                                                      *
 * experimental mean valus + normally distributed random noise scaled with the          *
 * mean variance                                                                        *
 *                                                                                      *
 ***************************************************************************************/

/*double DROPS_bound_cound_from_velo_inflow(int coord_x, int coord_y, int coord_z, int coord_r,
                                          int N_x, int N_y, int N_z, int N_r)
{
  double eps = 1e-6; 
  double value; 
  double stand_deviation;
  double f_infty = TDatabase::ParamDB->WINDTUNNEL_F_INFTY;
 
  value = 0.;

  // for the flow inlet 
  if(coord_x==0)
  {   
      // read mean values of variance and boundary value
      stand_deviation = sqrt(TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[coord_y][coord_z][coord_r][1]);
      value = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[coord_y][coord_z][coord_r][0] 
	  + normal_rand()* stand_deviation;
      // non-dimensionalize value
      //OutPut(coord_y << " " <<  coord_r << " : " << 
      //     TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[coord_y][coord_z][coord_r][0]  <<  endl);
      //     " " << variance << " " << value << endl);
      value = value/f_infty;
      //if (coord_r==48)
      //  OutPut(TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[coord_y][coord_z][coord_r][0] << endl);
      if (value<0) 
	value = 0;
      else 
	{
	  if (value>0)
	    ;//OutPut("PSD bdr " << r << " "  << value << endl);
	}

   }
  return value;
}*/
double DROPS_bound_cound_from_velo_inflow(int coord_x, int coord_y, int coord_z, double r,
                                          int N_x, int N_y, int N_z, int N_r)
{
  double eps = 1e-6; 
  double value; 
  double stand_deviation;
  double f_infty = TDatabase::ParamDB->WINDTUNNEL_F_INFTY;
  double log_norm_par[46][19][6];
  double  mu_dsd, sigma_dsd, scal_dsd, mu_stand_dev, sigma_stand_dev, scal_stand_dev;
  double mean_value, stand_dev;
  double mu, sigma, scal, rand;
  value = 0.;

  // for the flow inlet 
 // if(coord_x==0)
  {   
    lognormal_parameter(log_norm_par);
    sigma_dsd = log_norm_par[coord_y][coord_z][0];
    mu_dsd = log_norm_par[coord_y][coord_z][1]; 
    scal_dsd = log_norm_par[coord_y][coord_z][2];
    sigma_stand_dev = log_norm_par[coord_y][coord_z][3];
    mu_stand_dev = log_norm_par[coord_y][coord_z][4]; 
    scal_stand_dev = log_norm_par[coord_y][coord_z][5];
    rand=normal_rand();
   /* mu = mu_dsd+ rand* mu_stand_dev;
    sigma = sigma_dsd+ rand* sigma_stand_dev;
    scal = scal_dsd+ rand* scal_stand_dev;
    value  = calc_log_normal( r ,  sigma, mu, scal);*/
    mean_value =calc_log_normal( r ,  sigma_dsd, mu_dsd, scal_dsd);
    stand_dev =calc_log_normal( r ,  sigma_stand_dev, mu_stand_dev, scal_stand_dev);
    value = mean_value /*+ normal_rand() * stand_dev*/;
   
      // non-dimensionalize value
       value = value/f_infty; //ich denke f ist schon entdimensionalisiert 
      if (value<0) 
	value = 0;
      else 
	{
	  if (value>0)
	    ;
	}

   }
  return value;

}

double calc_log_normal(double x, double sigma, double mu, double scal)
{
double factor, exp_term;
factor = scal/(x*sigma*sqrt(2.*Pi));
exp_term = exp(-(log(x)-mu)*(log(x)-mu)/(2.*sigma*sigma));
return factor*exp_term;
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

void Windtunnel_FWE_FDM_Upwind_4D(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
double *f_old,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *a_coord,
double x_min, double x_max, double y_min, double y_max,
double z_min, double z_max, double a_min, double a_max,
double *velo1, double *velo2, double *velo3,
int *correspond_3dgrid, double ***diff_velo_air_drops)
{
  int i, ii, j, N2, N3, N4, alpha, beta, gamma, no_of_3dcell;
  int very_first = 0;

  double G_d_val1,G_d_val2, val;
  double velocity1_array_val, velocity2_array_val, velocity3_array_val;
  double values[4];
  double *f_new, *derx_val, velo_drops_diff;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  int delta;
  // model constants
  double envir_cond = TDatabase::ParamDB->WINDTUNNEL_ENVIR_COND;
  double supersat = TDatabase::ParamDB->WINDTUNNEL_SUPERSAT;
  double u_infty = TDatabase::ParamDB->WINDTUNNEL_U_INFTY;
  double l_infty = TDatabase::ParamDB->WINDTUNNEL_L_INFTY;
  double r_infty = TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
  double f_infty = TDatabase::ParamDB->WINDTUNNEL_F_INFTY;
 
  double factor_growth;
  int x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat; 
 
  // no particles in the domain
  if (TDatabase::TimeDB->CURRENTTIME<TDatabase::TimeDB->T1)
    return;

  // computed model constants
  factor_growth = l_infty * envir_cond * supersat/(u_infty * r_infty * r_infty);

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
     
      // correction for different velocities of air and drops
      compute_coordinate(i, &x_coord_mat, &y_coord_mat, &z_coord_mat, &a_coord_mat, 
       N_x, N_y, N_z, N_a);
      velo_drops_diff = calc_velo_u1(diff_velo_air_drops, x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat,  N_x,  N_y,  N_z,  N_a);
     // cout << "velo_drops_diff: " << velo_drops_diff << endl;
     // double velo_air_u1 = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[y_coord_int][z_coord_int][0];
      velo1[i] =values[0]-velo_drops_diff;
      velocity1_array_val = velo1[i];

      velocity2->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
        x_coord[i],y_coord[i],z_coord[i],values);
      velo2[i] = values[0];
      velocity2_array_val = velo2[i];

      velocity3->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
        x_coord[i],y_coord[i],z_coord[i],values);
      velo3[i] = values[0];
      velocity3_array_val = velo3[i];
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
    }

    G_d_val1 = factor_growth/a_coord[i];
    G_d_val2 = factor_growth/(a_coord[i]*a_coord[i]);
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

    // compute the first growth term
    if (G_d_val1 >= 0)
    {
      // not on boundary a = a_min -> in the "a-sense" exists a left neighbour
      if ( i >= N3 )
        derx_val[i] += G_d_val1*(f_old[i]-f_old[i-N3])/(fabs(a_coord[i]-a_coord[i-N3]));
    }
    else
    {
      // not on boundary a = a_max -> in the "a-sense" exists a left neighbour
      if ( i < N3*N_a )
        derx_val[i] += G_d_val1*(f_old[i+N3]-f_old[i])/(fabs(a_coord[i+N3]-a_coord[i]));
    }

    //compute second growth term
    derx_val[i] -= G_d_val2 * f_old[i];

    // compute new particle size distribution
    f_new[i] = f_old[i] - deltat * derx_val[i];

    // set Dirichlet boundary conditions
    // inflow for a = a_min
   

    // set Dirichlet boundary conditions
    // inflow from the left:experimental data
    if ( i < N3 )
    {
	f_new[i] = 0.;
    }
    if ( i%(N_x+1)==0 && a_coord[i]>0. )
      {
      compute_coordinate(i, &x_coord_mat, &y_coord_mat, &z_coord_mat, &a_coord_mat, 
       N_x, N_y, N_z, N_a);
      val = DROPS_bound_cound_from_velo_inflow(x_coord_mat, y_coord_mat, z_coord_mat, a_coord[i], 
					       N_x,  N_y,  N_z, N_a);
      f_new[i] = val;
      }
 
  }



  // copy new particle size distribution into array for old one
  Dcopy(N4, f_new, f_old);

  // cut undershoots
  for ( i=0 ; i<N4 ; i++ )
  {
    if (f_old[i]<0)
    {
      f_old[i] = 0;
    }
  }

  /*  if (i> ((N_y)/2) * N_x && i< (((N_y)/2)+1)*N_x )
    
      OutPut(TDatabase::TimeDB->CURRENTTIME * l_infty/u_infty << " PSD " 
	     << x_coord[i] * l_infty << " "
	     << y_coord[i] * l_infty << " "  
	     << z_coord[i] * l_infty<< " " 
	     << a_coord[i] * r_infty  << " " 
	     << f_old[i] * f_infty  << endl);*/
   

  OutPut("norm PSD after " << Ddot(N4,f_old,f_old)<<endl);
  // free allocated memory
  delete f_new;
  delete derx_val;
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

void Windtunnel_BWE_FDM_Upwind_4D(TCollection *coll,TFEFunction3D *velocity1, 
			TFEFunction3D *velocity2, TFEFunction3D *velocity3,
			double *sol, int *correspond_3dgrid,
			 int N_x, int N_y, int N_z, int N_a, double x_min, double x_max,
			 double y_min, double y_max, double z_min,
                         double z_max, double a_min, double a_max,
			 double *x_coord, double *y_coord, double *z_coord,
			 double *a_coord, double ***diff_velo_air_drops, TSquareMatrix3D *mat)
{
  int i, ii, jj, iq, k, N2, N3, N4, diag_index, indices[9];
  int maxind, index, index1, index2;
  int alpha, beta, gamma, no_of_3dcell, N_Entries, range;
  int *col_ptr, *row_ptr;
  int SC_LDS = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  int x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat; 

  double G_d_val1,G_d_val2;
  double yq, val, t3;
  double velocity1_array_val, velocity2_array_val, velocity3_array_val;
  double values[4], C_val[4];
  double *velo1, *velo2, *velo3;
  double velo_drops_diff;
  double *entries, *rhs;
  double time = TDatabase::TimeDB->CURRENTTIME;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;

  // model constants
  double envir_cond = TDatabase::ParamDB->WINDTUNNEL_ENVIR_COND;
  double supersat = TDatabase::ParamDB->WINDTUNNEL_SUPERSAT;
  double u_infty = TDatabase::ParamDB->WINDTUNNEL_U_INFTY;
  double l_infty = TDatabase::ParamDB->WINDTUNNEL_L_INFTY;
  double r_infty = TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
  double f_infty = TDatabase::ParamDB->WINDTUNNEL_F_INFTY;
  double factor_growth;

  // no particles in the domain
  if (TDatabase::TimeDB->CURRENTTIME<TDatabase::TimeDB->T1)
    return;
  factor_growth = l_infty * envir_cond * supersat/(u_infty * r_infty * r_infty);

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

  rhs = new double[N4];
  // array for the rhs
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
      //OutPut("b ii " << ii << endl);
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
      //OutPut("ii " << ii << " ga " << gamma << " be " << beta << " ga " <<
      //       gamma << " : " <<  gamma*N_x*N_y+ alpha * N_x + beta << endl);
      no_of_3dcell = correspond_3dgrid[gamma*N_x*N_y+ alpha * N_x + beta];

     
      // correction for different velocities of air and drops
      compute_coordinate(i, &x_coord_mat, &y_coord_mat, &z_coord_mat, &a_coord_mat, 
       N_x, N_y, N_z, N_a);
      velo_drops_diff = calc_velo_u1(diff_velo_air_drops, x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat,  N_x,  N_y,  N_z,  N_a); 
      // find velocity
      velocity1->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
        x_coord[i],y_coord[i],z_coord[i],values);
      velo1[i] = values[0]-velo_drops_diff;
      velocity1_array_val = velo1[i];

      velocity2->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
        x_coord[i],y_coord[i],z_coord[i],values);
      velo2[i] = values[0];
      velocity2_array_val = velo2[i];

      velocity3->FindGradientLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
        x_coord[i],y_coord[i],z_coord[i],values);
      velo3[i] = values[0];
      velocity3_array_val = velo3[i];
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
    }

    G_d_val1 = factor_growth/a_coord[i];
    G_d_val2 = factor_growth/(a_coord[i]*a_coord[i]);

    // compute the coefficients corresponding to the 4d finite-difference-method
    // diagonal entry
    // mass term and second growth term
    for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
    {
      if (col_ptr[k]==i)
      {
        entries[k] = 1 - deltat * G_d_val2;
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
        val =  deltat * velocity1_array_val/fabs(x_coord[i]-x_coord[i-1]);
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

    // compute the growth-term
    if (G_d_val1 >= 0)
    {
      // not on boundary a = a_min -> in the "a-sense" exists a left neighbour
      if ( i >= N3 )
      {
        index1 = i-N3;
        for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
        {
          if (col_ptr[k] == index1)
          {
            val = deltat * G_d_val1/fabs(a_coord[i]-a_coord[index1]);
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
            val = deltat * G_d_val1/fabs(a_coord[index1]-a_coord[i]);
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
    // inflow for a = a_min
   

    // set Dirichlet boundary conditions
    // inflow from the left: experimental data
  for (i=0;i<N4;i++)
  {
   //x==0 or first cube 
   if ((i%(N_x+1)==0)||(i<N3))
    { //inflow for a = a_min
    if (i<N3)
       {
        val = 0.;
       }
     if (i%(N_x+1)==0 && a_coord[i]>0.) 
      {//inflow from the left
       compute_coordinate(i, &x_coord_mat, &y_coord_mat, &z_coord_mat, &a_coord_mat, 
                            N_x, N_y, N_z, N_a);
      val = DROPS_bound_cound_from_velo_inflow(x_coord_mat, y_coord_mat, z_coord_mat, a_coord[i], 
					       N_x,  N_y,  N_z, N_a);
      }
    
      
      rhs[i] = val;
      sol[i] = val;
      // off diagonals
      for ( iq = row_ptr[i];  iq < row_ptr[i+1]; iq++ )
      {
        // diagonal entry
        if(col_ptr[iq]==i)
          entries[iq] = 1.0;
        else
          entries[iq] = 0.;
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
  OutPut("norm PSD after " << Ddot(N4,sol,sol)<<endl);

  // free allocated memory
  delete rhs;
  delete velo1;
  delete velo2;
  delete velo3;
}

void FEM_FCT_Matrix_Q1_4D_Windtunnel(TCollection *coll,
				TFEFunction3D *velocity1, TFEFunction3D *velocity2, 
				TFEFunction3D *velocity3, double ***diff_velo_air_drops,
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
  int N4 = N3 * (N_a+1);

  double a[256], b[256], u_val[6], val_test[5], val_ansatz[5], val_sol[5];
  double x_coord_loc[16], y_coord_loc[16], z_coord_loc[16], a_coord_loc[16];
  double sol_loc[16], b_sol[16];
  double *entries, *entriesM, *oldrhs_fem_fct0, *rhs, *RhsArray, *tilde_u, *u1, *u2, *u3;
  double *bdr_val, *entriesM_cons;
  double x_max, y_max, z_max, a_max, a_min;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, aq, val, weight_det;
  int x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat; 
  double maxsol, norm_b, al, G_d_val1, G_d_val2, eps = 1e-6;
  double time = TDatabase::TimeDB->CURRENTTIME, t1, t2, t3;
  double velo_drops_diff;

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
  if (TDatabase::TimeDB->CURRENTTIME<TDatabase::TimeDB->T1)
    return;

  t1 = GetTime();

 //model constants
  double envir_cond = TDatabase::ParamDB->WINDTUNNEL_ENVIR_COND;
  double supersat = TDatabase::ParamDB->WINDTUNNEL_SUPERSAT;
  double u_infty = TDatabase::ParamDB->WINDTUNNEL_U_INFTY;
  double l_infty = TDatabase::ParamDB->WINDTUNNEL_L_INFTY;
  double r_infty = TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
  double f_infty = TDatabase::ParamDB->WINDTUNNEL_F_INFTY;

 // compute model constants
  double factor_growth;
  factor_growth = l_infty * envir_cond * supersat/(u_infty * r_infty * r_infty);
  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();
  // initialize test_cells
  test_cells = new int[N_cells];
  memset(test_cells, 0, N_cells*SizeOfInt);
  // matrix objects
  entries = mat->GetEntries();
  N_Entries = mat->GetN_Entries();
  memset(entries,0, N_Entries*SizeOfDouble);
  col_ptr = mat->GetKCol();
  row_ptr = mat->GetRowPtr();
  entriesM = matM->GetEntries();
  entriesM_cons = matM_cons->GetEntries();

  // copy entries of mass matrix
  memcpy(entriesM, entriesM_cons, N_Entries*SizeOfDouble);

  Nodes = N3 * (N_a+1);
  x_max = x_coord[Nodes-1];
  y_max = y_coord[Nodes-1];
  z_max = z_coord[Nodes-1];
  a_max = a_coord[Nodes-1];
  a_min = a_coord[0];

  // initialization of u1, u2, u3, G in the vector u1
  // N_x*N_y*N_z - number of mesh cells for flow domain
  // need in each mesh cell 16 values (16 quad points) for 3 functions
  u1 = new double[N_x*N_y*N_z*48];
  memset(u1,0,(N_x*N_y*N_z*48)*SizeOfDouble);
  u2 = u1 + N_x*N_y*N_z*16;
  u3 = u2 + N_x*N_y*N_z*16;

  // initialization of rhs and the vectors for FEM_FCT_ForConvDiff in the vector rhs
  rhs = new double[Nodes*4];
  memset(rhs,0,Nodes*4*SizeOfDouble);
  tilde_u  = rhs +  Nodes;
  RhsArray = tilde_u + Nodes;
  oldrhs_fem_fct0 = RhsArray + Nodes;

  bdr_val = new double[N_neum_to_diri];
  memset(bdr_val,0,N_neum_to_diri*SizeOfDouble);

  z = z1 = 0;
 
  t2 = GetTime();
  OutPut("time windtunnelfct (1) " << t2-t1 << endl);

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
      
      // correction for different velocities of air and drops in u1-direction
       compute_coordinate(i, &x_coord_mat, &y_coord_mat, &z_coord_mat, &a_coord_mat, 
               N_x, N_y, N_z, N_a);
       velo_drops_diff = calc_velo_u1(diff_velo_air_drops, x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat,
               N_x, N_y, N_z,  N_a);
        u1[z1] = u_val[0]- velo_drops_diff;
        u2[z1] = u_val[1];
        u3[z1] = u_val[2];
        z1++;
      }
      // convection for internal coordinate
      G_d_val1 = factor_growth/aq;
      // reaction term 
      G_d_val2 = factor_growth/(aq*aq);
      
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
          Compute_Q1_Value_Gradient_4D(b+16*jj, xq, yq, zq, aq, val_ansatz);
          entries[index] += (u1[z]*val_ansatz[1]+u2[z]*val_ansatz[2]+u3[z]*val_ansatz[3] 
			     + G_d_val1 * val_ansatz[4] - G_d_val2 * val_ansatz[0])*val_test[0];
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
  OutPut("time windtunnelfct (2) " << t2-t1 << endl);
  // FEM--FCT
  // set Dirichlet boundary conditions on all (possible) inflow boundaries
  // on the "first cube" (internal coordinate inlet)
  count = 0;
  for ( i=0 ; i< N3 ; i++ )
  {    
      if (fabs(x_coord[i])<eps && fabs(a_coord[i])>eps)
	  continue;
      val  = 0.;
      bdr_val[count] = val;
      neum_to_diri[count] = i;
      topdiri = 1;
      count++;
  }
  // on the flow inlet
  for (i=0;i<Nodes;i++)
  {	
      if ( i%(N_x+1)==0 && a_coord[i]>0.)
      {
	  neum_to_diri[count] = i;
	  compute_coordinate(i, &x_coord_mat, &y_coord_mat, &z_coord_mat, &a_coord_mat, 
			     N_x, N_y, N_z, N_a);
	  bdr_val[count] = DROPS_bound_cound_from_velo_inflow(x_coord_mat, y_coord_mat, z_coord_mat, a_coord[i], 
							      N_x,  N_y,  N_z, N_a);
	    count++; 
      }
  }
   
  // this sets the array rhs
  FEM_FCT_ForConvDiff((TSquareMatrix3D*) matM, (TSquareMatrix3D*) mat,
		      Nodes, Nodes,
		      lump_mass_PSD, matrix_D_Entries_PSD,
		      sol, oldsol,
		      rhs, RhsArray, oldrhs_fem_fct0, tilde_u,
		      N_neum_to_diri, neum_to_diri,
		      NULL, NULL, NULL,
		      1, NULL, bdr_val);

  t2 = GetTime();
  OutPut("time windtunnelfct (3) " << t2-t1 << endl);
  // build matrix for FEM-FCT
  matM->Reset();
  FEM_FCT_SystemMatrix(matM, mat, lump_mass_PSD, Nodes);

  t2 = GetTime();
  OutPut("time windtunnelfct (4) " << t2-t1 << endl);

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
  OutPut("time windtunnelfct (5) " << t2-t1 << endl);
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
	TDatabase::ParamDB->SC_VERBOSE_AMG = 2;
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
  OutPut("norm PSD after " << Ddot(N4,sol,sol)<<endl);

  /*for ( i=0 ; i<N4 ; i++ )
  {
    if ((i)%(N_x+1)==1)
    {
      OutPut(TDatabase::TimeDB->CURRENTTIME * l_infty/u_infty << " PSD " 
	     << x_coord[i] * l_infty << " "
	     << y_coord[i] * l_infty << " "  
	     << z_coord[i] * l_infty<< " " 
	     << a_coord[i] *  r_infty  << " " 
	     << sol[i] * f_infty  << endl);
    }
  }*/
  memcpy(oldsol, sol, Nodes * SizeOfDouble);

  OutPut(TDatabase::TimeDB->CURRENTTIME << " Solver done " << sqrt(Ddot(Nodes,sol,sol)) << " max " << maxsol << endl);

 // reset indices if necessary
  /*if (topdiri)
  {
      for (i=0;i<N3;i++)
	  neum_to_diri[i] = i;  
  }*/

  t2 = GetTime();
  OutPut("time windtunnelfct (6) " << t2-t1 << endl);

  // deletion of the arrays
  delete bdr_val;
  delete rhs;
  delete u1;
  delete test_cells;
}


void Compute_Neum_To_Diri_FEM_FCT_Windtunnel(int N_x, int N_y, int N_z, int N_a,
				  double *x_coord, double *y_coord, 
				  double *z_coord, double *a_coord,
				  int &N_neum_to_diri, 
				  int* &neum_to_diri,
				  double* &neum_to_diri_x,
				  double* &neum_to_diri_y,
				  double* &neum_to_diri_z,
				  double* &neum_to_diri_a)
{
    int i, Nodes, range, count = 0;
    double yq, eps = 1e-6;
    Nodes = (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1);
   
     // total number of Dirichlet nodes
    if (TDatabase::ParamDB->WINDTUNNEL_R_MIN > 0)
	count = (N_x)*(N_y+1)*(N_z+1);//bottom
    else
	count = (N_x+1)*(N_y+1)*(N_z+1);//bottom
    count += (N_z+1) * (N_y+1) * (N_a); //first layer
    
    neum_to_diri = new int[count];
    neum_to_diri_x = new double[4*count];
    neum_to_diri_y = neum_to_diri_x + count;
    neum_to_diri_z = neum_to_diri_y + count;    
    neum_to_diri_a = neum_to_diri_z + count;    
    N_neum_to_diri = count;    

    count = 0;
    // 3D flow domain
    for ( i=0 ; i< (N_x+1)*(N_y+1)*(N_z+1) ; i++ )
    {
        // bottom
	if ((fabs(x_coord[i])<eps)&&(TDatabase::ParamDB->WINDTUNNEL_R_MIN > 0))
	    continue;
	neum_to_diri[count] = i;
	neum_to_diri_x[count] = x_coord[i];
	neum_to_diri_y[count] = y_coord[i];
	neum_to_diri_z[count] = z_coord[i];
	neum_to_diri_a[count] = a_coord[i];
	count++;
    }

    for (i=(N_x+1)*(N_y+1)*(N_z+1);i<Nodes;i++)
    {	
        // inlet
	if ( i%(N_x+1)==0 )
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

void Build_4D_FEM_FCT_MassMatrix_Q1_Windtunnel(TCollection *coll,
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
  double x_max, y_max, z_max, a_max, a_min;
  double area, detJK, hK, hK_conv, tauK, xq, yq, zq, aq, val, weight_det;
  double t1, t2;
  //double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;

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

  // exception for a_min
  /*if (fabs(a_min-d_p_min) > 1e-8)
  {
    OutPut("a_min " << a_min << " does not correspond to d_p_min " << d_p_min << endl);
    exit(4711);
  }*/

  // allocate array for storing indices for assembling
  ii = N_cells * N_a * 256;
  index_test_ansatz = new int[ii];
  OutPut("index_test_ansatz " << ii << endl);

  memset(entriesM,0, N_Entries*SizeOfDouble);

  z = 0;
  t2 = GetTime();
  OutPut("time windtunnelmass (1) " << t2-t1 << endl);

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
    area = (x_coord[i+1]-x_coord[i])*(y_coord[i+N_x+1]-y_coord[i])*
	(z_coord[i+N2]-z_coord[i])*(a_coord[i+N3]-a_coord[i]);
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
	  if (test_index!=ansatz_index)
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

  LumpMassMatrixToVector((TSquareMatrix3D*) matM, lump_mass_PSD);
  t2 = GetTime();
  OutPut("time windtunnelmass (2) " << t2-t1 << endl);
}

double normal_rand()
{//function produces normally distributed random variable X(0,1)
//BOX-Muller-scheme 
 //drand48 produces uniformly distributed random variable on [0,1)
 //must be initialized with  srand48(time(NULL))
  double r, t;
   r=sqrt(-2*log(1-drand48()));
   t=2*Pi*drand48();
   return r*cos(t);  
}

double log_normal_rand(double mu, double sigma, double scal)
 //mu mean value // sigma standard deviation of a normal distributed random variable
{//function produces log_normal distributed random variable X(mu,sigma^2)
 //drand48 produces uniformly distributed random variable on [0,1)
 //must be initialized with  srand48(time(NULL))
  double norm_rand = mu+normal_rand()*sigma;
  return scal* exp(norm_rand);  
}
/*****************************************************************************************
 *                                                                                       *
 *  .vtk file to visualize f for a given cut_coordinate of a                             *
 *                                                                                       *
 ****************************************************************************************/

void write_vtk_file_yzlayer( int N_x, int N_y, int N_z, int N_a, double x_fix_coord,
                     double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                     double *f_old, const char *name)
{
  int i;
  int N_Cells = N_y * N_z * N_a;
  int N_Nodes = (N_y+1)*(N_z+1)*(N_a+1);
  int N1 = N_x+1;
  int N2 = N1*(N_y+1);
  int N3 = N2*(N_z+1);
  // routine should work with the first N3-coordinates of x,y,z because a_coord is const. always in N3 blocks
  int N4 = N3*(N_a+1);
 
  double x_max = x_coord[N4-1];
  double y_max = y_coord[N4-1];
  double z_max = z_coord[N4-1];
  double a_max = a_coord[N4-1];
 int *help;
  help = new int[N4];


  FILE* out = fopen(name,"w");

  fprintf(out,"%s\n","# vtk DataFile Version 4.2" );
  fprintf(out,"%s\n","file created by MooNMD" );
  fprintf(out,"%s\n","ASCII" );
  fprintf(out,"%s\n","DATASET UNSTRUCTURED_GRID" );
  fprintf(out,"\n");
  fprintf(out,"%s","POINTS " );
  fprintf(out,"%i",N_Nodes);
  fprintf(out,"%s\n"," float" );
 int k=0;
  for ( i=0 ; i<N4 ; i++ )
    {
     help[i]=-1;
      if(fabs(x_coord[i]-x_fix_coord)<1e-6)
      {
      fprintf(out,"%f", y_coord[i]);
      fprintf(out,"%s"," ");
      fprintf(out,"%f", z_coord[i]);
      fprintf(out,"%s"," ");
      fprintf(out,"%f\n", a_coord[i]);
      help[i] = k;
      k++;
      }
    }
  fprintf(out,"\n");
  fprintf(out,"\n");

  fprintf(out,"%s","CELLS ");
  fprintf(out,"%i",N_Cells);
  fprintf(out,"%s"," ");
  fprintf(out,"%i",N_Cells*9);
  fprintf(out,"\n");

  for ( i=0 ; i<N4 ; i++ )
  {
   if(fabs(x_coord[i]-x_fix_coord)<1e-6)
   {
      if (fabs(x_coord[i]-x_max)<1e-6 || fabs(y_coord[i]-y_max)<1e-6
         || fabs(z_coord[i]-z_max)<1e-6 || fabs(a_coord[i]-a_max)<1e-6)
           continue;                                               // assign nodes with mesh cell if it is not at x_max, y_max, 
//if (help[i]<0||help[i+N1]<0||help[i+N2]<0||help[i+N1+N2]<0||help[i+N3]<0||help[i+N2+N3]<0||help[i+N1+N2+N3]<0)
    // exit(-1);
      fprintf(out,"%i",8);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i]);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i+N1]);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i+N1+N3]);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i+N3]);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i+N2]);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i+N1+N2]);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i+N1+N2+N3]);
      fprintf(out,"%s"," ");
      fprintf(out,"%i",help[i+N2+N3]);
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
  fprintf(out,"%i\n", N_Nodes);
  fprintf(out,"%s\n","SCALARS f float");
  fprintf(out,"%s\n","LOOKUP_TABLE default");

 for ( i=0 ; i<N4 ; i++ )
    {
  if(fabs(x_coord[i]-x_fix_coord)<1e-6)
   {
    fprintf(out,"%f\n",f_old[i]);
    }
  }

  fprintf(out,"\n");
  fclose(out);
delete help;
}

//function writes .txt file that contains the mean values 
//of the DSD at the outlet of the wind tunnel
void write_data_file_meanvalue(double ***mean_value_outflow,  double *a_layers_coord, int N_x, int N_y ,int N_z, int N_a, const char *name)  
{
 int i,j,r;
  FILE* out = fopen(name,"w");
  for(i=0;i<N_y+1;i++)
      for(j=0;j<N_z+1;j++)
         {
         fprintf(out,"--------***--------\n");
         fprintf(out," 40 %i %i\n\n", i, j );
         for(r=0;r<N_a+1;r++)
            fprintf(out,"%f %e \n",a_layers_coord[r], mean_value_outflow[i][j][r]);
         }
  fclose(out);
}



//func. finds to a given mesh node the cartesian coordinates
void compute_coordinate(int index, int *coord_x, int *coord_y, 
                   int *coord_z, int *coord_a, int N_x, int N_y, 
                     int N_z, int N_a)
 {
 int help;
 *coord_x = (index%(N_x+1));
 help = (int) (index/(N_x+1));
 *coord_y = (help%(N_y+1));
 help = (int) (help/(N_y+1));
 *coord_z = (help%(N_z+1));
 help = (int) (help/(N_z+1));
 *coord_a = (help%(N_a+1));
  }


void midpointflow(double *x_coord, double *y_coord, 
                   double *z_coord, double *a_coord, int N_x, int N_y, 
                     int N_z, int N_a, double *sol )
{
  int inf_coord, out_coord, N3, ii,j,i;
  double u_infty = TDatabase::ParamDB->WINDTUNNEL_U_INFTY;
  double l_infty = TDatabase::ParamDB->WINDTUNNEL_L_INFTY;
  double r_infty = TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
  double f_infty = TDatabase::ParamDB->WINDTUNNEL_F_INFTY;
  
  // number of nodes in 3D
  N3= (N_x+1)*(N_y+1)*(N_z+1);
  // output of f at the outflow x_max
  for ( i=0 ; i<N3 ; i++ )
  {
      if (fabs(x_coord[i]-0.)< 1e-6 && fabs(y_coord[i]-0.005)<1e-6 && fabs(z_coord[i])< 1e-6)
	  inf_coord = i;
      if (fabs(x_coord[i]-0.4)< 1e-6 && fabs(y_coord[i]-0.005)<1e-6 && fabs(z_coord[i])< 1e-6)
	  out_coord = i;
  }
      
  for (i=0;i<=N_a;i++)
  {
      j = inf_coord + N3*i;
      ii = out_coord + N3*i;
      OutPut(TDatabase::TimeDB->CURRENTTIME * l_infty/u_infty << " psd_inf_out " 
	     << a_coord[j] * r_infty << " " 
	     << sol[j] * f_infty  << " " << sol[ii] * f_infty  << endl);
  }
}


void compute_mean_value_outflow(int ***mean_value_outflow_indices, double ***mean_value_outflow, double *sol, double* x_coord,
                      int N_x, int N_y, int N_z, int N_a, int *only_first_time) 
   { double offset =0.0;
    int i,j,r, number_timesteps;
    int  N4, coord_x_mat,coord_y_mat, coord_z_mat, coord_a_mat; 
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    double time = TDatabase::TimeDB->CURRENTTIME; 
    //calculate number of current time step
    number_timesteps= (int)(time/tau);
    
     //initializations, must only be done once 
    //store node numbers of the outflow cube in mean_value_outflow_indices
    if(*only_first_time==1)
      {
      N4= (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1);
       
      for(i=0;i<N4;i++)
       {
        if (fabs(x_coord[i]-0.0)<1e-6)
          {
           compute_coordinate(i, &coord_x_mat, &coord_y_mat, &coord_z_mat, &coord_a_mat, N_x,  N_y, N_z, N_a);
           mean_value_outflow_indices[coord_y_mat][coord_z_mat][coord_a_mat]=i; 
         //cout << "i: " << i << coord_x_mat <<  " x " << coord_y_mat << " y " <<   coord_z_mat  << " z " << coord_a_mat  << endl;
          }
        }
     //mean_value_outflow is initialized to 0
       for (i=0;i<N_y+1;i++)
          for (j=0;j<N_z+1;j++)
             for(r=0;r<N_a+1;r++)
               {
               mean_value_outflow[i][j][r]=0.;
               } 
       *only_first_time=0;
       }

   if((time-offset)>=tau)// start mean value calculation when droplets reach the outlet(here t=0.2)
    {
     number_timesteps =(int) ((time-offset)/tau); //shift the current time step
      //cout << "tau" << tau  << " " <<number_timesteps;
      for (i=0;i<N_y+1;i++)
        for (j=0;j<N_z+1;j++)
            for(r=0;r<N_a+1;r++)
		{ 
                 // calculate mean value
		 mean_value_outflow[i][j][r]  =  (double)(number_timesteps-1)/(double)(number_timesteps)*mean_value_outflow[i][j][r]
                      + 1./(double)(number_timesteps) * sol[mean_value_outflow_indices[i][j][r]];
                }
//cout << sol[mean_value_outflow_indices[0][0][1]]<<endl;
     }
   }
         
double calc_velo_u1(double ***diff_velo_air_drops, int coord_int_x, int coord_int_y, int coord_int_z,
       int  a_coord_int, int N_x, int N_y, int N_z, int N_a)
{
 
  int coord_x, coord_y, coord_z;
  int delta_x, delta_y, delta_z;
  int  dim_x, dim_y, dim_z;
  double dist_x,dist_y, dist_z;
  double value, value1, value2;
 
  //dim_x Number of available files in x-direction,
  // number of measuered layers here dim_x=3
  dim_x = TDatabase::ParamDB->WINDTUNNEL_LAYER_NUMBER_X;
  dim_y = TDatabase::ParamDB->WINDTUNNEL_DIM_Y;
  dim_z = TDatabase::ParamDB->WINDTUNNEL_DIM_Z;
 
  value = 0.;
  
    //values y,z are in a range of 0-0.45 or 0-0.18
    //for this reason they must be converted in values between 0-45 or 0-18
 
    
    // in x-direction the end of the tunnel is assumed to be 0.40m,
    // althougt the calculated lenght is 0.5m 
    // for this reason the 40 is insertet in x-direction
     
    // the last 0.1m the that is measured at 0.4 is taken
    // To force that the value for the last 10 cm in x-direction is the same 
    // for all x_values, the index of the x_value is cut at dim_x-1, so the algorithm 
    // averages between the same value
    
    //difference between two neighbouring grid points
    delta_x = (int) ((40)/(dim_x-1));
    delta_y = (int) ((N_y)/(dim_y-1));
    delta_z = (int) ((N_z)/(dim_z-1));
     
    //calculate corresponding indices
    coord_x  = (int) (coord_int_x / delta_x);        //left neighbour
    coord_y  = (int) (coord_int_y / delta_y);
    coord_z  = (int) (coord_int_z / delta_z);
   
   
  //distance between coordinate and left neighbour
    dist_x = fabs((double) (coord_int_x - coord_x*delta_x)/(double) (delta_x));
    dist_y = fabs((double) (coord_int_y - coord_y*delta_y)/(double)(delta_y)); 
    dist_z = fabs((double) ((coord_int_z - coord_z*delta_z)/(double) delta_z));
    
    //the x index must be smaller then dim_x-1
    if (coord_x >dim_x-1) coord_x=dim_x-1;
  
    

    //interpolatation in y and z direction 
    value1 = (1-dist_y) * (1-dist_z) * diff_velo_air_drops[coord_x][coord_y][coord_z]
      +  dist_y    * (1-dist_z) *  diff_velo_air_drops[coord_x][coord_y+1][coord_z]
      + (1-dist_y) *  dist_z    *  diff_velo_air_drops[coord_x][coord_y][coord_z+1]
      +  dist_y    *  dist_z    *  diff_velo_air_drops[coord_x][coord_y+1][coord_z+1];

   value2 = (1-dist_y) * (1-dist_z) * diff_velo_air_drops[coord_x+1][coord_y][coord_z]
      +  dist_y    * (1-dist_z) *  diff_velo_air_drops[coord_x+1][coord_y+1][coord_z]
      + (1-dist_y) *  dist_z    *  diff_velo_air_drops[coord_x+1][coord_y][coord_z+1]
      +  dist_y    *  dist_z    *  diff_velo_air_drops[coord_x+1][coord_y+1][coord_z+1];
    //interpolation in x direction
    value = (1-dist_x)*value1 + dist_x * value2;
    //non-dimensionalisation
    value /= TDatabase::ParamDB->WINDTUNNEL_U_INFTY;
  
  //OutPut(coord_int_x<< " :: "<< coord_x <<"  "  << dist_x << " " << value1 << " " << value2<< " "<< value << endl);
 //OutPut(coord_int_y<< " :: "<< coord_y <<"  "  << dist_y << " " << " "<< coord_int_z << " :: "<< coord_z <<"  "  << dist_z <<  endl);
  return value;
}


void alloc_cubix( double ****cubix, int dim_x, int dim_y, int dim_z)
{
int i, j;
*cubix = (double ***) malloc (dim_x *sizeof(double **));
for(i=0;i<dim_x;i++)
   { 
   (*cubix)[i]=(double **) malloc (dim_y* sizeof (double *));
   for(j=0; j< dim_y; j++)
       {
        (*cubix)[i][j]= (double *)malloc (dim_z * sizeof(double));
   
       }
   }
return;
}

void disalloc_cubix( double ***cubix, int dim_x, int dim_y, int dim_z)
{
int i, j;
for(i=0;i<dim_x;i++)
    for(j=0;j<dim_y;j++)
        free(cubix[i][j]);
for(i=0;i<dim_x;i++)
    free(cubix[i]);
free(cubix);
return;
}

void alloc_cubix_int(int ****cubix, int dim_x, int dim_y, int dim_z)
{
int i, j;
*cubix = (int ***) malloc (dim_x *sizeof(int **));
for(i=0;i<dim_x;i++)
   { 
   (*cubix)[i]=(int **) malloc (dim_y* sizeof (int *));
   for(j=0; j< dim_y; j++)
       {
        (*cubix)[i][j]= (int *)malloc (dim_z * sizeof(int));
   
       }
   }
return;
}

void disalloc_cubix(int ***cubix, int dim_x, int dim_y, int dim_z)
{
int i, j;
for(i=0;i<dim_x;i++)
    for(j=0;j<dim_y;j++)
        free(cubix[i][j]);
for(i=0;i<dim_x;i++)
    free(cubix[i]);
free(cubix);
return;
}




void alloc_matrix( double ***matrix, int dim_x, int dim_y)
{
int i;
*matrix = (double **) malloc (dim_x *sizeof(double *));
for(i=0;i<dim_x;i++)
      (*matrix)[i]=(double *) malloc (dim_y* sizeof (double));
 return;
}

void disalloc_matrix( double **matrix, int dim_x, int dim_y)
{
int i;
for(i=0;i<dim_x;i++)
    free(matrix[i]);
free(matrix);
return;
}
