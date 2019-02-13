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
 *                         UREA_3d4d.C                                                  *
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
#include <Urea_3d4d.h>
#include <Bulk.h>
#include <BrAgg.h>

#include <Database.h>
#include <Hexahedron.h>

#include <string.h>
#include <stdlib.h>

/****************************************************************************************
 *                                                                                      *
 *  The declared functions are sorted in the following way :                            *
 *                                                                                      *
 *    1. FDM routines                                                                   *
 *    2. FEM routines                                                                   *
 *    3. analysing routines                                                             *
 *                                                                                      *
 ***************************************************************************************/
int PSD_bound_cond_from_velo_inflow_urea(double x, double y, double z)
{
  int value=0;
  double eps = 1e-8;

  // inflow from left
  if ((fabs(x)<eps)&& (y>=1.0/3-eps) && (y<=2.0/3+eps) &&
    (z>=1.0/3-eps) && (z<=2.0/3+eps))
  {
    value =  1;                                  
    //OutPut(x << " "<< y << " " << z << endl);
  }
  // inflow from right
  // if (fabs(1 - x)<eps)
  //{
  //if ((fabs(y-0.5)<r2+eps)&&(fabs(z-0.5)<r2+eps))
  //  value = 1;
  //}
  return value;
}


/****************************************************************************************
 *                                                                                       *
 *                          Part II.1 : FDM (explicit with upwindig)                     *
 *                         -----------------------------------------                     *
 *                                                                                       *
 ****************************************************************************************/
//boundary condition


//AT THE MOMENT NO GROWTH SINCE BY DEFINITION OF THE DENOMINATOR
double growth_rate(double c, double temp)
{
  double k_g =  TDatabase::ParamDB->UREA_k_g;
  double g =  TDatabase::ParamDB->UREA_g;
  double c_infty =  TDatabase::ParamDB->UREA_c_infty;
  double m_mol =  TDatabase::ParamDB->UREA_m_mol;
  double temp_infty =  TDatabase::ParamDB->UREA_temp_infty;
  double val, val1;

  // without nucleation 
  if (TDatabase::ParamDB->UREA_MODEL == 2)
    return(0);

  val = TDatabase::ParamDB->UREA_rho_sat_1+TDatabase::ParamDB->UREA_rho_sat_2*(temp-273.15);

  val1 = m_mol*c_infty*c - val;
  if (val1 <= 0)
    return(0);

  val1 /= val;
  val1 = pow(val1,g);
  val1 *= k_g;
  return(val1);
    
}

double b_nuc(double c, double temp)
{
  double alpha =  TDatabase::ParamDB->UREA_alfa_nuc;
  double beta =  TDatabase::ParamDB->UREA_beta_nuc;
  double c_infty =  TDatabase::ParamDB->UREA_c_infty;
  double temp_infty =  TDatabase::ParamDB->UREA_temp_infty;
  double m_mol =  TDatabase::ParamDB->UREA_m_mol;
  double rho_urea_c_sat, val;

  if ((c<=0)||(TDatabase::ParamDB->UREA_MODEL == 2))
    return(0.0);

  rho_urea_c_sat = TDatabase::ParamDB->UREA_rho_sat_1+TDatabase::ParamDB->UREA_rho_sat_2*(temp-273.15);
  val = m_mol*c_infty*c/rho_urea_c_sat;
  if (val==0)
    return(0.0);
  val=log(val);
  val *= val;
  val = -beta/val;
  val = exp(val);
  val *= alpha;
  // set constant
  val = TDatabase::ParamDB->UREA_alfa_nuc;
  return(val);
}



double InletPSD(double a)
{
  /* double inlet_coord[55] = 
   {0, 9.9166500e-05, 1.1628750e-04, 1.3340900e-04,1.5053050e-04, 1.6765150e-04,   
   1.8477250e-04, 2.0189400e-04, 2.1901550e-04, 2.3613650e-04, 2.5325750e-04,   
   2.7037850e-04, 2.8750000e-04, 3.0462150e-04, 3.2174250e-04, 3.3886350e-04,   
   3.5598450e-04, 3.7310600e-04, 3.9022750e-04, 4.0734850e-04, 4.2446950e-04,  
   4.4159100e-04, 4.5871250e-04, 4.7583350e-04, 4.9295450e-04, 5.1007550e-04,  
   5.2719700e-04, 5.4431850e-04, 5.6143950e-04, 5.7856050e-04, 5.9568150e-04,   
   6.1280300e-04, 6.2992450e-04, 6.4704550e-04, 6.6416650e-04, 6.8128750e-04, 
   6.9840900e-04, 7.1553050e-04, 7.3265150e-04, 7.4977250e-04, 7.6689400e-04,   
   7.8401550e-04, 8.0113650e-04, 8.1825750e-04, 8.3537850e-04, 8.5250000e-04,  
   8.6962150e-04, 8.8674250e-04, 9.0386350e-04, 9.2098450e-04, 9.3810600e-04,
   9.5522750e-04, 9.7234850e-04, 9.8946950e-04, 1.0065910e-03};

   double inlet_f_L_seed[55] =   
    {0, 7.3261885e+12,  5.6038146e+12, 4.7123546e+12, 3.3738790e+12, 2.8909354e+12,
    2.6575694e+12,  2.4939436e+12, 2.7944043e+12, 5.2429412e+12, 8.8255590e+12,
    1.4622235e+13,  2.0629510e+13, 2.6507704e+13, 3.1684903e+13, 3.2411516e+13,
    2.9314408e+13,  2.6100742e+13, 2.2887775e+13, 1.8152796e+13, 1.3153219e+13,
    1.1061698e+13,  8.4815352e+12, 6.3228082e+12, 4.8351689e+12, 3.8342340e+12,
    2.8232362e+12,  2.1730245e+12, 1.7867633e+12, 1.4955162e+12,  9.5780195e+11,
    8.1155267e+11, 5.8733968e+11,  6.7309592e+11, 4.1687593e+11, 3.3071383e+11,
    2.4774618e+11,  1.5439009e+11, 7.4176972e+10,  1.0178349e+11, 5.8212958e+10,
     6.3479029e+10, 3.0547748e+10, 4.7521176e+10,  0.0000000e+00, 0.0000000e+00,
    3.1012405e+10,   0.0000000e+00, 1.5355082e+10, 1.5167527e+10, 1.4362716e+10,
    0.0000000e+00, 0.0000000e+00, 0.0000000e+00,  0.0000000e+00};
*/
  double inlet_coord[100] ={2.500000e-06,1.356050e-05,3.068150e-05,4.780300e-05,6.492450e-05,8.204550e-05,9.916650e-05,1.162875e-04,1.334090e-04,1.505305e-04,1.676515e-04,1.847725e-04,2.018940e-04,2.190155e-04,2.361365e-04,2.532575e-04,2.703785e-04,2.875000e-04,3.046215e-04,3.217425e-04,3.388635e-04,3.559845e-04,3.731060e-04,3.902275e-04,4.073485e-04,4.244695e-04,4.415910e-04,4.587125e-04,4.758335e-04,4.929545e-04,5.100755e-04,5.271970e-04,5.443185e-04,5.614395e-04,5.785605e-04,5.956815e-04,6.128030e-04,6.299245e-04,6.470455e-04,6.641665e-04,6.812875e-04,6.984090e-04,7.155305e-04,7.326515e-04,7.497725e-04,7.668940e-04,7.840155e-04,8.011365e-04,8.182575e-04,8.353785e-04,8.525000e-04,8.696215e-04,8.867425e-04,9.038635e-04,9.209845e-04,9.381060e-04,9.552275e-04,9.723485e-04,9.894695e-04,1.006591e-03,1.023712e-03,1.040833e-03,1.057955e-03,1.075075e-03,1.092197e-03,1.109318e-03,1.126439e-03,1.143561e-03,1.160682e-03,1.177803e-03,1.194925e-03,1.212045e-03,1.229167e-03,1.246287e-03,1.263409e-03,1.280530e-03,1.297651e-03,1.314773e-03,1.331894e-03,1.349016e-03,1.366137e-03,1.383257e-03,1.400379e-03,1.417500e-03,1.434622e-03,1.451742e-03,1.468863e-03,1.485985e-03,1.503106e-03,1.520228e-03,1.537348e-03,1.554470e-03,1.571591e-03,1.588713e-03,1.605834e-03,1.622954e-03,1.640075e-03,1.657197e-03,1.674318e-03,1.691440e-03};

  double inlet_f_L_seed[100] = {0.000000e+00,1.539144e+09,2.082543e+09,1.056886e+09,1.032820e+09,9.174217e+08,7.339395e+08,6.383108e+08,5.630089e+08,4.776040e+08,3.952118e+08,3.244646e+08,2.668657e+08,2.186269e+08,1.774119e+08,1.430176e+08,1.142823e+08,9.072941e+07,7.094795e+07,5.432124e+07,4.189605e+07,3.280686e+07,2.562972e+07,2.013453e+07,1.620424e+07,1.315460e+07,1.040263e+07,7.846039e+06,5.728553e+06,4.495501e+06,3.800058e+06,3.005911e+06,2.229857e+06,1.675695e+06,1.307405e+06,1.003372e+06,7.486895e+05,6.805768e+05,7.435860e+05,7.812756e+05,7.396855e+05,6.082531e+05,3.731120e+05,1.582579e+05,8.976509e+04,1.072953e+05,9.867716e+04,5.302092e+04,3.384872e+04,7.221184e+04,1.380859e+05,1.713877e+05,1.450725e+05,9.890296e+04,6.954686e+04,3.619885e+04,1.297526e+04,1.931888e+04,5.546763e+04,1.066047e+05,1.070156e+05,5.457797e+04,1.540227e+04,3.078500e+03,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};
 int i,n;
 int found;
 double val = 0,V_inj,m; 
 double L_max = TDatabase::ParamDB->UREA_D_P_MAX;
 found=0;
 n=100;

 //V_inj = TDatabase::ParamDB->UREA_INFLOW_SCALE*10.0/(9.0*TDatabase::ParamDB->UREA_u_infty);
 //V_inj = 1/V_inj;
 //V_inj *= 1e6;

 // flow rate in ml/min 
 V_inj = 20 * TDatabase::ParamDB->UREA_INFLOW_SCALE / (3.0 * TDatabase::ParamDB->UREA_u_infty);
 // injection time
 V_inj *= TDatabase::ParamDB->UREA_inflow_time;
 V_inj = 1/V_inj;
 // scaling to units
 V_inj *= 60 * 1e6;

 Dscal(100,V_inj,inlet_f_L_seed);
 for (i=n;i>0;i--)
   { 
      //OutPut(" a vorher ist !!!"<< a *L_max << "  i  "<<   i  << "  inlet  " <<inlet_coord[i]<<endl);
      //OutPut(i << " inlet " << inlet_coord[i] << " a ist " << a << " " << L_max<< endl);
      //if(a*L_max>=inlet_coord[i])
      if((a*L_max<=inlet_coord[i])&&(a*L_max>inlet_coord[i-1]))
      {
        //OutPut(" a ist !!!"<< a *L_max <<" i ist !!!" <<i << " i ist !!!"<<inlet_coord[i]<<endl);
        m=(inlet_f_L_seed[i]-inlet_f_L_seed[i-1])/(inlet_coord[i] - inlet_coord[i-1]);
        val=m *(a*L_max- inlet_coord[i-1]);
        val+=inlet_f_L_seed[i-1];
        found=1;
        break;
      }
     //OutPut("next " << i << " " << n-1 << endl);
    }
    if (found==0)
    {
      val= 0.;
    }
  //OutPut(" a ist "<< a *L_max << " "  << val <<endl);
return(val);
}

double psd_boundary_urea(double x, double y, double z, double a,
double c, double temp)
{
  double L_min = TDatabase::ParamDB->UREA_D_P_0;
  double L_max = TDatabase::ParamDB->UREA_D_P_MAX;
  double f_infty = TDatabase::ParamDB->UREA_f_infty;
  double eps = 1e-8, val, a_min, g;
  double f_in = TDatabase::ParamDB->UREA_f_infty;
  
  // inflow
  if (PSD_bound_cond_from_velo_inflow_urea(x,y,z))
  {  
	if (TDatabase::TimeDB->CURRENTTIME <= TDatabase::ParamDB->UREA_inflow_time)
        {	
      
          // size of the seed crystals
          val = InletPSD(a)/TDatabase::ParamDB->UREA_f_infty;
          //OutPut(a << " " << val << endl);
          return(val);
         }
         else 
          return(0);
  }

  a_min = L_min/L_max ;

  if (fabs(a-a_min)<eps)
  {
    g = growth_rate(c,temp);
    if (g<=0)
    {
      return(0.0);
    }
    val = b_nuc(c,temp);
    val /= (f_infty*g);
    /*if (val > 1)
        {
         return (1.0);
        }*/
    return(val);
  }
  return(0.0);
}

/****************************************************************************************
 *                                                                                      *
 * assemble the matrix corresponding to the 4d finite-difference-method                 *
 *                                                                                      *
 ***************************************************************************************/

void Urea_FWE_FDM_Upwind_4D(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *f_old, double *f_new,
double *rhs_psd,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *a_coord,
double x_min, double x_max, double y_min, double y_max,
double z_min, double z_max, double a_min, double a_max,
double *velo1, double *velo2, double *velo3,
double *concent_C_array, double *Temp_array,
int *correspond_3dgrid)
{
  int i, ii, j, N2, N3, N4, alpha, beta, gamma, no_of_3dcell;
  int maxind;
  int very_first = 0;

  double B_c_C, concent_C_array_val,Temp_array_val, maxsol, G_c_C_val,rho_sat_val, val;
  double velocity1_array_val, velocity2_array_val, velocity3_array_val;
  double values[4];
  double  *derx_val;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;

  // model constants
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
  double u_infty = TDatabase::ParamDB->UREA_u_infty;
  double c_infty = TDatabase::ParamDB->UREA_c_infty;
  double D_j = TDatabase::ParamDB->UREA_D_J;
  double d_p_0 = TDatabase::ParamDB->UREA_D_P_0;
  double d_p_max = TDatabase::ParamDB->UREA_D_P_MAX;
  double k_v = TDatabase::ParamDB->UREA_k_v;
  double k_g = TDatabase::ParamDB->UREA_k_g;
  double g = TDatabase::ParamDB->UREA_g;
  //double k_nuc = TDatabase::ParamDB->UREA_k_nuc;
  //double d_p_min = TDatabase::ParamDB->UREA_D_P_MIN;
  //double c_C_infty = TDatabase::ParamDB->UREA_c_C_infty;
  double f_infty = TDatabase::ParamDB->UREA_f_infty;
  double temp_infty = TDatabase::ParamDB->UREA_temp_infty;
  double rho_sat_1 = TDatabase::ParamDB->UREA_rho_sat_1;
  double rho_sat_2 = TDatabase::ParamDB->UREA_rho_sat_2;
  double rho =  TDatabase::ParamDB->UREA_rho;
  double m_mol =  TDatabase::ParamDB->UREA_m_mol;
  double alfa_nuc =  TDatabase::ParamDB->UREA_alfa_nuc;
  double beta_nuc =  TDatabase::ParamDB->UREA_beta_nuc;
  double eps, factor_G, T_infty, lambda_nuc;
  // double d_p_min;

  // computed model constants
  factor_G = k_g*l_infty/(u_infty*d_p_max);

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
  concent_C_array = new double[N3];
  Temp_array = new double[N3];
  memset(concent_C_array, 0, N3*SizeOfDouble);
  // array for temperature
  memset(Temp_array, 0, N3*SizeOfDouble);
  // array for derivatives of f with respect to x,y,z,a
  derx_val = new double[N4];
  memset(derx_val, 0, N4*SizeOfDouble);
  // array for new values of the particle size distribution f
  //f_new = new double[N4];
  //memset(f_new, 0, N4*SizeOfDouble);

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

      // fill the array for the concentrations of Temp
      // since this concentration does not depend on a
      // ...
      Temp->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
        x_coord[i],y_coord[i],z_coord[i],values);
      Temp_array[i] = values[0];
      Temp_array_val = values[0];
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
      Temp_array_val = Temp_array[ii];
    }

    G_c_C_val=growth_rate(concent_C_array_val,Temp_array_val);
    G_c_C_val *= l_infty /(u_infty * d_p_max);
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
    f_new[i] = f_old[i] - deltat * derx_val[i] + deltat * rhs_psd[i];

    // set Dirichlet boundary conditions
    // if convection is positive at the bottom ( = cube with a = a_min )
    if (i<N3)
    {
      f_new[i]=psd_boundary_urea(x_coord[i], y_coord[i], z_coord[i], a_coord[i],concent_C_array_val,          Temp_array_val);
    }
    // set Dirichlet boundary conditions
    // if convection is positive at the top
    if ((i>=N3*N_a )&&(G_c_C_val <  0))
    {
      f_new[i] = 0.0;
    }

    // set Dirichlet boundary conditions
    // inflow from the left x = x_min (left) or right x = x_max (right)
    if ( i%(N_x+1)==0 )
    {
      if (PSD_bound_cond_from_velo_inflow_urea(x_coord[i], y_coord[i], z_coord[i]))
      {
	  f_new[i]=psd_boundary_urea(x_coord[i], y_coord[i], z_coord[i], a_coord[i],
          concent_C_array_val, Temp_array_val);
      }
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
 OutPut("psd " << Ddot(N4,f_old,f_old) << endl);

  // free allocated memory
 // delete f_new;
  delete derx_val;
  delete concent_C_array;
  delete Temp_array;
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

void Urea_BWE_FDM_Upwind_4D(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *sol,
double *rhs_psd,
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

  double yq, B_c_C, concent_C_array_val, Temp_array_val, rho_sat_val, maxsol, G_c_C_val, val, t3;
  double velocity1_array_val, velocity2_array_val, velocity3_array_val;
  double values[4], C_val[4], Temp_val[4] ;
  double *velo1, *velo2, *velo3, *concent_C_array, *Temp_array;
  double *entries, *rhs;
  double time = TDatabase::TimeDB->CURRENTTIME;
  double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;

  // model constants
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
  double u_infty = TDatabase::ParamDB->UREA_u_infty;
  double c_infty = TDatabase::ParamDB->UREA_c_infty;
  double D_j = TDatabase::ParamDB->UREA_D_J;
  double d_p_0 = TDatabase::ParamDB->UREA_D_P_0;
  double d_p_max = TDatabase::ParamDB->UREA_D_P_MAX;
  double k_v = TDatabase::ParamDB->UREA_k_v;
  double k_g = TDatabase::ParamDB->UREA_k_g;
  double g = TDatabase::ParamDB->UREA_g;
  //double k_nuc = TDatabase::ParamDB->UREA_k_nuc;
  //double d_p_min = TDatabase::ParamDB->UREA_D_P_MIN;
  //double c_C_infty = TDatabase::ParamDB->UREA_c_C_infty;
  double f_infty = TDatabase::ParamDB->UREA_f_infty;
  double temp_infty = TDatabase::ParamDB->UREA_temp_infty;
  double rho_sat_1 = TDatabase::ParamDB->UREA_rho_sat_1;
  double rho_sat_2 = TDatabase::ParamDB->UREA_rho_sat_2;
  double rho =  TDatabase::ParamDB->UREA_rho;
  double m_mol =  TDatabase::ParamDB->UREA_m_mol;
  double alfa_nuc =  TDatabase::ParamDB->UREA_alfa_nuc;
  double beta_nuc =  TDatabase::ParamDB->UREA_beta_nuc;
  double eps, factor_G, T_infty, lambda_nuc;

  // computed model constants
  factor_G = k_g*l_infty/(u_infty*d_p_max);

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
  //array for temperature
  Temp_array = new double[N3];
  memset(Temp_array, 0, N3*SizeOfDouble);
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
      //fill the array for the temperature
      Temp->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
        x_coord[i],y_coord[i],z_coord[i],values);
      Temp_array[i] = values[0];
      Temp_array_val = values[0];
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
      Temp_array_val = Temp_array[ii];
    }

    G_c_C_val=growth_rate(concent_C_array_val,Temp_array_val);
    G_c_C_val *= l_infty /(u_infty * d_p_max);

    // compute the coefficients corresponding to the 4d finite-difference-method
    for ( k=row_ptr[i] ; k<row_ptr[i+1] ; k++ )
    {
      // diagonal entry
      if (col_ptr[k]==i)
      {
        // term from time derivative
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
    rhs[i] = sol[i]+deltat*rhs_psd[i];
    //rhs[i] = sol[i];
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

    Temp->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
      x_coord[i],y_coord[i],z_coord[i],Temp_val);

    G_c_C_val=growth_rate(C_val[0],Temp_val[0]);
    if (G_c_C_val*f_infty > 1e-10)
    {
      rhs[i] = psd_boundary_urea(x_coord[i],y_coord[i],z_coord[i],a_coord[i],C_val[0],Temp_val[0]);
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
      rhs[k] = psd_boundary_urea(x_coord[k],y_coord[k],z_coord[k],a_coord[k],C_val[0],Temp_val[0]);
      sol[k] = rhs[k];
      // set matrix rows
      ii = row_ptr[k];
      jj = row_ptr[k+1];
      // off diagonals
      for ( iq = ii ; iq < jj ; iq++ )
      {
        // diagonal entry
        if(col_ptr[iq]==k)
          entries[iq] = 1;
        else
          entries[iq] = 0;
      }
    }
  }

  // inflow from the velocity field
  for (i=N3;i<N4;i++)
  {
    if ( i%(N_x+1)==0 )
    {
      iq =  PSD_bound_cond_from_velo_inflow_urea(x_coord[i], y_coord[i], z_coord[i]);
      if (iq)
      {
        rhs[i] = psd_boundary_urea(x_coord[i],y_coord[i],z_coord[i],a_coord[i],0,0);
        sol[i] = rhs[i];
        ii = row_ptr[i];
        jj = row_ptr[i+1];
        // off diagonals
        for ( iq = ii ; iq < jj ; iq++ )
        {
          // diagonal entry
          if(col_ptr[iq]==i)
            entries[iq] = 1;
          else
            entries[iq] = 0;
        }
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
OutPut("psd " << Ddot(N4,sol,sol) << endl);
  // free allocated memory
  delete rhs;
  delete velo1;
  delete velo2;
  delete velo3;
  delete concent_C_array;
  delete Temp_array;
}


/****************************************************************************************
 *                                                                                      *
 *                         Part III : FEM                                               *
 *                        ----------------                                              *
 *                                                                                      *
 ***************************************************************************************/

/****************************************************************************************
 *                                                                                      *
 *  computes all dof which are Dirichlet dof                                            *
 *  these has to be treated as Neumann dof in FEM--FCT                                  *
 *                                                                                      *
 ***************************************************************************************/

void Urea_Compute_Neum_To_Diri_FEM_FCT(int N_x, int N_y, int N_z, int N_a,
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
    val = PSD_bound_cond_from_velo_inflow_urea(x_coord[i], y_coord[i],z_coord[i]);
    if (val)
      count++;
  }
  OutPut("Neum_To_Diri_FEM_FCT " << count << endl);

  // total number of Dirichlet nodes
  count += (N_x+1)*(N_y+1)*(N_z+1);
  neum_to_diri = new int[count];
  neum_to_diri_x = new double[4*count];
  neum_to_diri_y = neum_to_diri_x + count;
  neum_to_diri_z = neum_to_diri_y + count;
  neum_to_diri_a = neum_to_diri_z + count;
  N_neum_to_diri = count;

  // inflow wrt to internal coordinate, 3D flow domain
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
    val = PSD_bound_cond_from_velo_inflow_urea(x_coord[i], y_coord[i], z_coord[i]);
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

void Urea_Build_4D_FEM_FCT_Matrix_Q1(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *sol, double *oldsol,double *rhs_psd,
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

  double a[256], b[256], C_val[4],Temp_val[4], u_val[6], val_test[5], val_ansatz[5], val_sol[5];
  double x_coord_loc[16], y_coord_loc[16], z_coord_loc[16], a_coord_loc[16];
  double sol_loc[16], b_sol[16];
  double *entries, *entriesM, *oldrhs_fem_fct0, *rhs, *rhs_psd_old, *RhsArray, *tilde_u, *u1, *u2, *u3, *G;
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
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
  double u_infty = TDatabase::ParamDB->UREA_u_infty;
  double c_infty = TDatabase::ParamDB->UREA_c_infty;
  double D_j = TDatabase::ParamDB->UREA_D_J;
  double d_p_0 = TDatabase::ParamDB->UREA_D_P_0;
  double d_p_max = TDatabase::ParamDB->UREA_D_P_MAX;
  double k_v = TDatabase::ParamDB->UREA_k_v;
  double k_g = TDatabase::ParamDB->UREA_k_g;
  double g = TDatabase::ParamDB->UREA_g;
  double f_infty = TDatabase::ParamDB->UREA_f_infty;
  double temp_infty = TDatabase::ParamDB->UREA_temp_infty;
  double rho_sat_1 = TDatabase::ParamDB->UREA_rho_sat_1;
  double rho_sat_2 = TDatabase::ParamDB->UREA_rho_sat_2;
  double rho =  TDatabase::ParamDB->UREA_rho;
  double m_mol =  TDatabase::ParamDB->UREA_m_mol;
  double alfa_nuc =  TDatabase::ParamDB->UREA_alfa_nuc;
  double beta_nuc =  TDatabase::ParamDB->UREA_beta_nuc;
  double eps, factor_G, T_infty, lambda_nuc,G_c_C_val;
  //model constants

  // compute coefficients of the equation
  N_cells = coll->GetN_Cells();
  // initialize test_cells
  test_cells = new int[N_cells];
  memset(test_cells, 0, N_cells*SizeOfInt);

  // compute model constants
  factor_G = k_g*l_infty/(u_infty*d_p_max);

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

  rhs_psd_old = new double[Nodes];
  memset(rhs_psd_old,0,Nodes*SizeOfDouble);
 // memcpy(rhs,rhs_psd,Nodes*SizeOfDouble);
 // memcpy(RhsArray,rhs_psd_old,Nodes*SizeOfDouble);
  Dcopy(Nodes, rhs_psd, rhs);
  Dcopy(Nodes, rhs_psd_old, RhsArray);

  bdr_val = new double[N_neum_to_diri];
  memset(bdr_val,0,N_neum_to_diri*SizeOfDouble);

  z = z1 = 0;

  t2 = GetTime();
  OutPut("time ureafct (1) " << t2-t1 << endl);

   
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
        Temp->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,xq,yq,zq,Temp_val);
        // G(c_C)
        G_c_C_val = growth_rate(C_val[0],Temp_val[0]);
	G_c_C_val *= l_infty /(u_infty * d_p_max);
	
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
    }                                             // end quad points

    if ( z == N_x*N_y*N_z*16 )
    {
      z = 0;
    }
  }
  // end i

  t2 = GetTime();
  OutPut("time ureafct (2) " << t2-t1 << endl);
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
    Temp->FindValueLocal(coll->GetCell(no_of_3dcell),no_of_3dcell,
      x_coord[i],y_coord[i],z_coord[i],Temp_val);

    G_c_C_val = growth_rate(C_val[0],Temp_val[0]);
    G_c_C_val *= l_infty /(u_infty * d_p_max);
     if ( G_c_C_val*f_infty > 0 )
    {
      // compute rate of nucleation
      bdr_val[i]=psd_boundary_urea(x_coord[i],y_coord[i],z_coord[i],a_coord[i],C_val[0],Temp_val[0]);
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

  // on the the velocity inlet
  for ( i=N3 ; i<N_neum_to_diri ; i++ )
  {
      index = neum_to_diri[i]; 
      // temperature and concentration not of importance
      bdr_val[i]=psd_boundary_urea(x_coord[index],y_coord[index],z_coord[index],
				   a_coord[index],0,0);
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
  OutPut("time ureafct (3) " << t2-t1 << endl);
  // build matrix for FEM-FCT
  matM->Reset();
  FEM_FCT_SystemMatrix(matM, mat, lump_mass_PSD, Nodes);

  t2 = GetTime();
  OutPut("time ureafct (4) " << t2-t1 << endl);
   //memcpy(rhs_psd_old,rhs_psd,Nodes*SizeOfDouble);
 
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
  OutPut("time ureafct (5) " << t2-t1 << endl);
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
  OutPut("time ureafct (6) " << t2-t1 << endl);

  // deletion of the arrays
  delete bdr_val;
  delete rhs;
  delete rhs_psd_old;
  delete u1;
  delete test_cells;
}

double calculate_q_3(int N, double *a_layers_coord, double *sol_psd)
{
int j ;
  double L3, val_left, val_right,t;
 double integral = 0;
 
       for (j=0;j<N-1;j++)
	  {
           val_left = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j]*sol_psd[j];
           val_right = a_layers_coord[j+1]*a_layers_coord[j+1]*a_layers_coord[j+1]*sol_psd[j];
           integral = integral + (val_left + val_right)*(a_layers_coord[j+1] - a_layers_coord[j])*0.5;
          }
       for (j=0;j<N;j++)
	  {
               L3 = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j];
               sol_psd[j] = sol_psd[j] * L3 / integral;          
               
	  }
  return(0);
}
void Evaluate_f_at_outflow1(int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord,double *z_coord,
double *a_coord, double *f)
{
  int k, i0, i;
  double *size, *number, q3,d_p_max,f_infty;

  // centre of outflow is in (0.5,0.5,0)
  // compute index for value of f on a=0
  // this is correct only if N_x and N_y are even
  i0 = (int) (((N_x+1)*(N_y+1)*(N_z+1)-1)/2)+(int)((N_x+1)/2);

  if (( fabs(210-x_coord[i0])>1e-6 )||( fabs(0.5-y_coord[i0])>1e-6 )||( fabs(0.5-z_coord[i0])>1e-6 ))
  {
    OutPut("Wrong index in Evalute_f_at_outflow " << i0 << endl);
    OutPut("x " << x_coord[i0] << " y " << y_coord[i0] << " z " << z_coord[i0]<<   endl);
    OutPut("N_x and N_y has to be even such that there is a d.o.f. at the center of the outflow " << endl);
    exit(4711);
  }
  //OutPut("FOUND " << i0 << " " <<  x_coord[i0] << " " <<  y_coord[i0] << " " <<  z_coord[i0] << endl);
  // exit(1);

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

 q3= calculate_q_3(N_a+1,size,number);
  OutPut(TDatabase::TimeDB->CURRENTTIME << " psd q3 " << q3 << endl);
 delete size;
  delete number;
}


void Calculate_Volume_Distribution(TCollection *Coll,int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
double *a_layers_coord, double *f)
{
  int j,i,k,m,N_V,N_Cells,N1,N2,N3,alpha,beta,gamma,index;
  double S_x,S_y,S_z,x_min,x_max,volume,vol_int;
  double x[8],y[8],z[8],*integral;
  double f_infty = TDatabase::ParamDB->UREA_f_infty;
  TBaseCell *cell;
  TVertex *vertex[8];
  double eps=1e-6;
  x_min=190.0;
  x_max=200.0;

  N1 = N_x+1;
  N2 = (N_x+1)*(N_y+1);
  N3 = N2*(N_z+1);
  // last entry: volume
  integral = new double[N_a+2];
  memset(integral, 0, (N_a+2)*SizeOfDouble);
  vol_int=0.0;
  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  //loop over the mesh cells of the global grid
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    //volume = cell->GetMeasure();
    if (N_V < 8)
    {
      OutPut(" only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    S_x=S_y=S_z=0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex[j] = cell->GetVertex(j);
      vertex[j]->GetCoords(x[j], y[j], z[j]);
      //baricenter of the mesh cell
      S_x+=x[j];
      S_y+=y[j];
      S_z+=z[j];
    }
    // bary center
    S_x/=N_V;
    S_y/=N_V;
    S_z/=N_V;

    if ((S_x >= x_min)&&(S_x <= x_max))
      //integrate for these mesh cells
    {
      // volume of the mesh cell
      //volume = GetMeasure(vertex);
      volume = cell->GetMeasure();
      // update volume of averaging domain
      integral[N_a+1] += volume;
      for (j=0;j<N_V;j++)
      {
	  /* // compute index of psd array
        for (k=0;k<N1;k++)
        {
          if (fabs(x[j]-x_coord[k])<eps)
          {
            beta=k;
            break;
          }
        }
        for (k=0;k<=N2;k+=N1)
        {
          if (fabs(y[j]-y_coord[k])<eps)
          {
            alpha=k/N1;
            break;
          }
        }
        for (k=0;k<=N3;k+=N2)
        {
          if (fabs(z[j]-z_coord[k])<eps)
          {
            gamma=k/N2;
            break;
          }
        }
        // index of psd for smallest particle size
        index=gamma*N_x*N_y+alpha*N_x+beta;
        index=gamma*N2+alpha*N1+beta;
	  */
	  // corresponding index of array f - comparision of the x, y and z coordinates
	  for ( m=0 ; m<N1 ; m++ )
	  {
	      if ( fabs(x[j]-x_coord[m])<eps )
	      {
		  index = m;
		  break;
	      }
	  }
	  //OutPut(x_coord[m]);
	  
	  for ( m=0 ; m<N2 ; m+=N1 )
	  {
	      if ( fabs(y[j]-y_coord[m])<eps )
	      {
		  index += m;
		  break;
	      }
	  }
	  
	  for ( m=0 ; m<N3 ; m+=N2 )
	  {
	      if ( fabs(z[j]-z_coord[m])<eps )
	      {
		  index += m;
		  break;
	      }
	  }
	  
	  //OutPut(x[j] << " " << y[j] << " " << z[j] << " " << gamma << " "  << alpha << " " << beta << 
	  // " " << index << endl);
	  //OutPut(x[j] << " " << y[j] << " " << z[j] << 
	  // " " << index << endl);
	  for (k=0;k<=N_a;k++)
	  {
	      integral[k]+= volume*f[index+k*N3]/N_V;
	  }
	  //OutPut(endl);
      }
    }
    // OutPut("coords " << x[j] << " " << y[j] << " " << z[j] << endl);
  }
  for (k=0;k<=N_a;k++)
  {
    OutPut(TDatabase::TimeDB->CURRENTTIME << " volume " << a_layers_coord[k] 
	   << " " << f_infty*integral[k]/integral[N_a+1]  << endl);
  }

  delete integral;
}


//***********************************************************
//
// computes PSD and q3 in a node
//
//***********************************************************

void Calculate_PSD_on_node(int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
double *a_layers_coord, double *sol_psd, double x, double y, double z)
{
  int N2, N3, i, j, k;
  double eps = 1e-5, eps1 = 1e-12;
  double L3, val_left, val_right;
  double integral = 0;
  // number of nodes in the flow domain
  N2 = (N_x+1)*(N_y+1);
  N3 = N2*(N_z+1);
  
  // look for node
  for (i=0;i<N3;i++)
  {
      // node for d_p_0 found
      if ((fabs(x_coord[i]-x)<eps) && (fabs(y_coord[i]-y)<eps) 
	  && (fabs(z_coord[i]-z)<eps))
      { 
          // calculate third moment of PSD in that node
          for (j=0;j<N_a;j++)
	  {

             val_left = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j]*sol_psd[i+j*N3];
             val_right = a_layers_coord[j+1]*a_layers_coord[j+1]*a_layers_coord[j+1]*sol_psd[i+(j+1)*N3];
             integral = integral + (val_left + val_right)*(a_layers_coord[j+1] - a_layers_coord[j])*0.5;
          }
          // calculate q3
	  for (j=0;j<=N_a;j++)
	  {   

              if (integral < eps1)
              {
                  OutPut("q3 integral " << integral << endl);
                    break;
              }
               else
               {
                   
                  L3 = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j];
      	          OutPut(TDatabase::TimeDB->CURRENTTIME << " PSD " << "("<<x<<","<<y<<","<<z<<") " << a_layers_coord[j] << " " << sol_psd[i+j*N3] << " q3 " << sol_psd[i+j*N3] * L3/integral << endl);
               }          
	  }
	  break;
      }
  }
}

void Calculate_PSD_outflow(TCollection *Coll,int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
double *a_layers_coord, double *sol_psd, int *step,double x_end)
{
  int N2, N3,N4, i,i1, ii,j, k, found,N_V,N_Cells,index, beta, alpha, gamma, val0, val1,m;
  double eps = 1e-5, eps1 = 1e-6;
  double x[8],y[8],z[8],values[4];
  double L3, val_left, val_right;
  double *average_q3,*PSD_array;
  double integral = 0;
    TBaseCell *cell;
    TVertex *vertex[8];
   int *indextest;
    
  // number of nodes in the flow domain
  N2 = (N_x+1)*(N_y+1);
  N3 = N2*(N_z+1);
  N4 = N3*(N_a+1);
  
  PSD_array= new double[N_a+1];
  memset(PSD_array,0,(N_a+1)*SizeOfDouble); 
   
  for (i=0;i<N3;i++)
  { 
      if ((fabs(x_coord[i]-x_end)<eps1))
	
      { 
	//&&(i>=N_y*(N_x+1))&&(i%(N_x+1)==0)&&((i+1)%(N_x+1)==0)
	//if ((i<(N_x+1)))
	//{
	 
	 //OutPut("step " <<step[0] <<"index " <<index << " " << x_coord[index] << " " << y_coord[index] << " " << z_coord[index] << " "<<endl);
	  step[0] ++;
	  //OutPut("step " <<step[0] <<"index " <<i<< " " << x_coord[i] << " " << y_coord[i] << " " << x_coord[i] << " "<<endl);
	//}
	  for (j=0;j<=N_a;j++)
	  {   
	    // OutPut(TDatabase::TimeDB->CURRENTTIME << " step " <<step[0] <<" index " << i+j*N3 << " " << x_coord[i+j*N3] << " " << y_coord[i+j*N3] 
	      //     << " " << z_coord[i+j*N3] << " PSD "<< sol_psd[i+j*N3]<< " " <<endl);
	        
                  PSD_array[j] += sol_psd[i+j*N3];
  	  }
	
  	  integral = 0;
 	  for (j=0;j<N_a;j++)
	  {

             val_left = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j]*sol_psd[i+j*N3];
             val_right = a_layers_coord[j+1]*a_layers_coord[j+1]*a_layers_coord[j+1]*sol_psd[i+(j+1)*N3];
             integral = integral + (val_left + val_right)*(a_layers_coord[j+1] - a_layers_coord[j])*0.5;
          }
          // calculate q3
	  for (j=0;j<=N_a;j++)
	  {   

              /*if (integral < eps1)
              {
                  OutPut(" q3 integral " << integral << endl);
                    break;
              }
               else*/
               {
                   
                  L3 = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j];
 	    OutPut(TDatabase::TimeDB->CURRENTTIME << " step " <<step[0] <<" index " << i+j*N3 << " " << x_coord[i+j*N3] << " " << y_coord[i+j*N3] 
	           << " " << z_coord[i+j*N3] << " "  << a_layers_coord[j] << " PSD "<< sol_psd[i+j*N3] <<" q3 " << sol_psd[i+j*N3]*L3/integral << endl);
     	         // OutPut(TDatabase::TimeDB->CURRENTTIME << " PSD all " << "("<<x_end <<",y,z) " << a_layers_coord[j] << " " << PSD_array[j] << " q3 " <<  PSD_array[j] * L3/integral << endl);
		 //OutPut("step "<< step[0]<<" "<<TDatabase::TimeDB->CURRENTTIME << " PSD " << "("<<x_end<<","<<y_coord[i] <<","<<z_coord[i]<<") " << a_layers_coord[j] << " " << PSD_array[j]<< " q3 " <<  PSD_array[j] * L3/integral<<  endl);
               }          
	  }
     } 
  }

 
  	  integral = 0;
	  for (j=0;j<N_a;j++)
	  {

             val_left = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j]*PSD_array[j];
             val_right = a_layers_coord[j+1]*a_layers_coord[j+1]*a_layers_coord[j+1]*PSD_array[j+1];
             integral = integral + (val_left + val_right)*(a_layers_coord[j+1] - a_layers_coord[j])*0.5;
          }
          // calculate q3
	  for (j=0;j<=N_a;j++)
	  {   

              if (integral < eps1)
              {
                  OutPut("q3 integral " << integral << endl);
                    break;
              }
               else
               {
                   
                  L3 = a_layers_coord[j]*a_layers_coord[j]*a_layers_coord[j];
      	         OutPut(TDatabase::TimeDB->CURRENTTIME << " PSD all " << "("<<x_end <<",y,z) " << a_layers_coord[j] << " " << PSD_array[j] << " q3 " <<  PSD_array[j] * L3/integral << endl);
		 //OutPut("step "<< step[0]<<" "<<TDatabase::TimeDB->CURRENTTIME << " PSD " << "("<<x_end<<","<<y_coord[i] <<","<<z_coord[i]<<") " << a_layers_coord[j] << " " << PSD_array[j]<< " q3 " <<  PSD_array[j] * L3/integral<<  endl);
               }          
	  }

delete PSD_array;
}



void PrepareAgglomerationBreakage(TCollection *Coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *temperature,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
double *a_layers_coord, double *f, double *rhs_new)
{
  int j,i,ii,k,N_V,N_Cells,N3,index, beta, alpha, gamma, val0, val1;
    double eps=1e-6, fac;
    double x[8],y[8],z[8],values[4];
    double params[7];
    double *velo, *grad_velo, *temp, *layers_mass,*shear_rate;
    double f_infty = TDatabase::ParamDB->UREA_f_infty;
    double l_infty = TDatabase::ParamDB->UREA_l_infty;
    double u_infty = TDatabase::ParamDB->UREA_u_infty;
    double L_max = TDatabase::ParamDB->UREA_D_P_MAX;
   
    TBaseCell *cell;
    TVertex *vertex[8];
 
    params[0]= TDatabase::ParamDB->UREA_AGGR_SPATIAL;//space
    params[1]= TDatabase::ParamDB->UREA_AGGR_BROWNIAN; // brown kernel included or not
    params[2]= TDatabase::ParamDB->UREA_AGGR_POL_ORDER;//pol or constant
    params[3]= TDatabase::ParamDB->UREA_AGGR_BROWNIAN_TEMP;//brown kernel depends of temperature or not
    params[4]= TDatabase::ParamDB->UREA_AGGR_BROWNIAN_SCAL;//scal param for brown kernel
    params[5]= TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR_TYPE;//shear depends of velocity or not
    params[6]= TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR;//param for shear kernel

    
    N3 = (N_x+1)*(N_y+1)*(N_z+1);
    velo = new double[3*N3];
    grad_velo = new double[9*N3];
    shear_rate = new double[N3];
    temp = new double[N3];
    
    N_Cells = Coll->GetN_Cells();
    //loop over the mesh cells of the global grid
    for(i=0;i<N_Cells;i++)
    {
	cell = Coll->GetCell(i);
	N_V = cell->GetN_Vertices();
	
	//volume = cell->GetMeasure();
	if (N_V < 8)
	{
	    OutPut(" only for hexahedral mesh cells implemented !!!"<<endl);
	    exit(4711);
	}
	for (j=0;j<N_V;j++)
	{
	    // read coordinates of the mesh cell
	    vertex[j] = cell->GetVertex(j);
	    vertex[j]->GetCoords(x[j], y[j], z[j]);
	}
	
	for (j=0;j<N_V;j++)
	{
	    // compute index of psd array
	    for (k=0;k<=N_x;k++)
	    {
		if (fabs(x[j]-x_coord[k])<eps)
		{
		    beta=k;
		    break;
		}
	    }
	    val0 = (N_x+1)*(N_y+1);
	    for (k=0;k<val0;k+=N_x+1)
	    {
		if (fabs(y[j]-y_coord[k])<eps)
		{
		  alpha=k/(N_x+1);
		    break;
		}
	    }
	    val1 = (N_x+1)*(N_y+1)*(N_z+1);
	    for (k=0;k<val1;k+=val0)
	    {
	      // OutPut(z[j] << " k " << k << " " << z_coord[k] << " " << fabs(z[j]-z_coord[k]) << endl);
		if (fabs(z[j]-z_coord[k])<eps)
		{
		    gamma=k/val0;
		    break;
		}
	    }
	    // index of psd for smallest particle size
	    index=gamma*(N_x+1)*(N_y+1)+alpha*(N_x+1)+beta;
           // if (fabs(x[j] -200.0) < 1e-6)
	   // {
	    //OutPut(index << " " << x[j] << " " << y[j] << " " << z[j] << " "<< " Nx "<<N_x<< " Ny "<<N_y<<" ");
	    //OutPut(gamma << " " << alpha << " " << beta << endl);
	    //}
	    // compute velocity and gradient
	   // shear_rate[index]=0;
	    velocity1->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    velo[3*index] = values[0]*u_infty;
	    grad_velo[9*index] = values[1]*u_infty/l_infty;
	    grad_velo[9*index+1] = values[2]*u_infty/l_infty;
	    grad_velo[9*index+2] = values[3]*u_infty/l_infty;
	    
	   // shear_rate[index] += values[1] * values[1];
            //shear_rate[index] += values[2] * values[2];
            //shear_rate[index] += values[3] * values[3];
	    
	    velocity2->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    velo[3*index+1] = values[0]*u_infty;
	    grad_velo[9*index+3] = values[1]*u_infty/l_infty;
	    grad_velo[9*index+4] = values[2]*u_infty/l_infty;
	    grad_velo[9*index+5] = values[3]*u_infty/l_infty;
	    
	    
	    //shear_rate[index] += values[1] * values[1];
            //shear_rate[index] += values[2] * values[2];
            //shear_rate[index] += values[3] * values[3];  
	    
	    velocity3->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    velo[3*index+2] = values[0]*u_infty;
	    grad_velo[9*index+6] = values[1]*u_infty/l_infty;
	    grad_velo[9*index+7] = values[2]*u_infty/l_infty;
	    grad_velo[9*index+8] = values[3]*u_infty/l_infty;
	    
	    //shear_rate[index] += values[1] * values[1];
            //shear_rate[index] += values[2] * values[2];
            //shear_rate[index] += values[3] * values[3];
     
            //shear_rate[index] = u_infty/l_infty * sqrt(2*shear_rate[index]);
      
          
      
	
	    temperature->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    temp[index] = values[0] * TDatabase::ParamDB->UREA_temp_infty;
            {
              if (PSD_bound_cond_from_velo_inflow_urea(x[j], y[j], z[j]))
              {
                 for (ii=0;ii<=N_a;ii++)
                 {
                    f[index+ii*N3]=psd_boundary_urea(x[j], y[j], z[j], a_layers_coord[ii],
                                        1, values[0]);
                    //OutPut(ii << " aaa "<< x[j]<< " " << y[j] << " " <<  z[j]<< " : " <<  
                    // " a " << a_layers_coord[ii] << " " << f[index+ii*N3] << endl); 
                 }
              }
            }
	    //OutPut(index << " " << temp[index] << " : ");
	}
    }
    layers_mass = new double[N_a+1];
    for (ii=0;ii<=N_a;ii++)
    {
         layers_mass[ii] = pow(a_layers_coord[ii],3.0);
        /* OutPut(layers_mass[ii] );
         if (ii >0)
            OutPut(" " <<  layers_mass[ii] -  layers_mass[ii-1] << endl)
         else
            OutPut(endl);*/
    	}

    Dscal(N3 * (N_a+1), f_infty, f);
    memset(rhs_new,0,N3*(N_a+1)*SizeOfDouble);
    OutPut("call pbe routines "<< Ddot( (N_x+1)*(N_y+1)*(N_z+1)*(N_a+1),rhs_new,rhs_new) <<endl);
     
    OutPut("Grad velo " << Ddot(9*N3,grad_velo,grad_velo) <<" "
	" temp " << Ddot(N3,temp,temp) << " f " << Ddot(N3* (N_a+1), f,f) << endl);
	
	
    OutPut("MEMORY before integral: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl);

   apply_integral_operator(N_x+1, N_y+1, N_z+1, N_a+1,
			 f,
			 velo, grad_velo, temp,
			 rhs_new, layers_mass,
			 TDatabase::ParamDB->UREA_D_P_MAX,
			 TDatabase::ParamDB->UREA_f_infty,
    			 params	
			 );
    OutPut("MEMORY after integral: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl);

   OutPut("call pbe routines done (dim) "<< sqrt(Ddot( N3*(N_a+1),rhs_new,rhs_new)) <<endl);
   fac = l_infty/(u_infty*f_infty);
    Dscal(N3 * (N_a+1), fac, rhs_new);
   OutPut("call pbe routines done "<< sqrt(Ddot( N3*(N_a+1),rhs_new,rhs_new)) <<endl);

    Dscal(N3 * (N_a+1), 1.0/f_infty, f);
 
 // OutPut("pbe routines done"<< endl);
  delete layers_mass;
  delete temp;
  delete grad_velo;
  delete velo;
}

double kernel_urea_old(double x, double x_bar)
{
double kernel;
double c_agg = TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR;
double eps=1.e-8;
double k_B = 1.3806504e-23;
double mu = 1.074e-3;
double T=298.;
double fac;

fac = 2.*k_B*T/(3*mu);

if ((x> eps)&&(x_bar> eps ))
{
kernel= TDatabase::ParamDB->P9 *fac*(x+x_bar)*(1./x+1./x_bar);
}
else
{
return(0.0);
}
kernel+= c_agg * (x+x_bar)*(x+x_bar)*(x+x_bar);
return(kernel);
}

/*double kernel(TCollection *Coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *temperature,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord, double x, double x_bar)
{
double kernel;
double c_agg = TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR;
double eps=1.e-8;
double k_B = 1.3806504e-23;
double mu = 1.074e-3;
double T=293.;
double fac_brown, fac_shear;


    TBaseCell *cell;
    TVertex *vertex[8];

    fac = 2.*k_B/(3*mu);
N3 = (N_x+1)*(N_y+1)*(N_z+1);
    velo = new double[3*N3];
    grad_velo = new double[9*N3];
    temp = new double[N3];
    
    N_Cells = Coll->GetN_Cells();
    //loop over the mesh cells of the global grid
    for(i=0;i<N_Cells;i++)
    {
	cell = Coll->GetCell(i);
	N_V = cell->GetN_Vertices();
	
	//volume = cell->GetMeasure();
	if (N_V < 8)
	{
	    OutPut(" only for hexahedral mesh cells implemented !!!"<<endl);
	    exit(4711);
	}
	for (j=0;j<N_V;j++)
	{
	    // read coordinates of the mesh cell
	    vertex[j] = cell->GetVertex(j);
	    vertex[j]->GetCoords(x[j], y[j], z[j]);
	}
	
	for (j=0;j<N_V;j++)
	{
	    // compute index of psd array
	    for (k=0;k<=N_x;k++)
	    {
		if (fabs(x[j]-x_coord[k])<eps)
		{
		    beta=k;
		    break;
		}
	    }
	    val0 = (N_x+1)*(N_y+1);
	    for (k=0;k<val0;k+=N_x+1)
	    {
		if (fabs(y[j]-y_coord[k])<eps)
		{
		  alpha=k/(N_x+1);
		    break;
		}
	    }
	    val1 = (N_x+1)*(N_y+1)*(N_z+1);
	    for (k=0;k<val1;k+=val0)
	    {
	      // OutPut(z[j] << " k " << k << " " << z_coord[k] << " " << fabs(z[j]-z_coord[k]) << endl);
		if (fabs(z[j]-z_coord[k])<eps)
		{
		    gamma=k/val0;
		    break;
		}
	    }
	    // index of psd for smallest particle size
	    index=gamma*(N_x+1)*(N_y+1)+alpha*(N_x+1)+beta;
	    //OutPut(index << " " << x[j] << " " << y[j] << " " << z[j] << " ");
	    //OutPut(gamma << " " << alpha << " " << beta << endl);
	    // compute velocity and gradient
	    velocity1->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    velo[3*index] = values[0]*u_infty;
	    grad_velo[9*index] = values[1]*u_infty/l_infty;
	    grad_velo[9*index+1] = values[2]*u_infty/l_infty;
	    grad_velo[9*index+2] = values[3]*u_infty/l_infty;
	    velocity2->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    velo[3*index+1] = values[0]*u_infty;
	    grad_velo[9*index+3] = values[1]*u_infty/l_infty;
	    grad_velo[9*index+4] = values[2]*u_infty/l_infty;
	    grad_velo[9*index+5] = values[3]*u_infty/l_infty;
	    velocity3->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    velo[3*index+2] = values[0]*u_infty;
	    grad_velo[9*index+6] = values[1]*u_infty/l_infty;
	    grad_velo[9*index+7] = values[2]*u_infty/l_infty;
	    grad_velo[9*index+8] = values[3]*u_infty/l_infty;
	    temperature->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    temp[index] = values[0] * TDatabase::ParamDB->UREA_temp_infty;
            fac_brown=fac_brown*temp[index];
            fac_shear=sqrt(2)*sqrt(Ddot(9*N3,grad_velo,grad_velo));
            OutPut("shear " << sqrt(2)*sqrt(Ddot(9*N3,grad_velo,grad_velo)) << endl);
	    //OutPut(index << " " << temp[index] << " : ");
                  if ((x> eps)&&(x_bar> eps ))
                      {
                      kernel=fac_brown* (x+x_bar)*(x+x_bar)/(x*x_bar)+fac_shear * (x+x_bar)*(x+x_bar)*(x+x_bar);
                       
                       }
                         else 

                       kernel=0 ;
	}
    }


}
*/
// compute function value, this has to be improved
double psd_value(double *f, double *a_layers_coord, double val, int index, int N_a, int N3)
{
    int i;
    double fval;

    for (i=0;i<N_a;i++)
    {
	if (a_layers_coord[i] >= val)
	{
	    fval = (f[(i-1)*N3+index] + f[i*N3+index])/2.0;
	    break;
	}
    }
    return(fval);
}

double gauss_source_urea_old(int j, int index, double *a_layers_coordinate, 
		    double a_layers_coordinate_outerLoop, double *f, int N3, int Na)
{
  int i;
  double val, val1, val2, val3, val_outer, integral, source, scal, fval;
  double stuetz[3] ={-sqrt(3./5.) , 0., sqrt(3./5.)};
  double gew[3] ={5./9., 8./9., 5./9.};

  // interpolation of function values
  double funcp1=0.5*(f[index+(j+1)*N3]-f[index+j*N3])*(stuetz[0]+1)+f[index+j*N3];
  double funcp2=0.5*(f[index+(j+1)*N3]-f[index+j*N3])*(stuetz[1]+1)+f[index+j*N3];
  double funcp3=0.5*(f[index+(j+1)*N3]-f[index+j*N3])*(stuetz[1]+1)+f[index+j*N3];
  // quadrature points
  double point1 = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2*(stuetz[0]+1)+a_layers_coordinate[j];
  double point2 = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2*(stuetz[1]+1)+a_layers_coordinate[j];
  double point3 = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2*(stuetz[2]+1)+a_layers_coordinate[j];
  // transform of variables
  scal=(a_layers_coordinate[j+1]-a_layers_coordinate[j])/2;

  //calculate source term
  val_outer= a_layers_coordinate_outerLoop*a_layers_coordinate_outerLoop*a_layers_coordinate_outerLoop;

  // first quadrature point
  val=val_outer-(point1*point1*point1);
  val1=pow(val,1.0/3.0);
  for(i=0;i<Na;i++)
    {
      if (a_layers_coordinate[i] > val1)
	{
	  // find PSD for first factor in integral 
	  fval = (f[index+(i+1)*N3] - f[index+i*N3])/(a_layers_coordinate[i+1]-a_layers_coordinate[i])*(val1 - a_layers_coordinate[i]) + f[index+i*N3];
	  break;
	}
    }
  val1 = kernel_urea_old(val1,point1)/(val1*val1)*fval*funcp1;

  // second quadrature point
  val=val_outer-(point2*point2*point2);
  val2=pow(val,1.0/3.0);
  for(i=1;i<=Na;i++)
    {
      if (a_layers_coordinate[i] > val2)
	{
	  // find PSD for first factor in integral 
	  fval = (f[index+(i+1)*N3] - f[index+i*N3])/(a_layers_coordinate[i+1]-a_layers_coordinate[i])*(val2 - a_layers_coordinate[i]) + f[index+i*N3];
      break;
    }
  }
  val2 = kernel_urea_old(val2,point2)/(val2*val2)*fval*funcp2;

  
  // third quadrature point
  val=val_outer-(point3*point3*point3);
  val3=pow(val,1.0/3.0);
  // find PSD for first factor in integral 
  for(i=1;i<=Na;i++)
  {
    if (a_layers_coordinate[i] > val3)
    {
      fval = (f[index+(i+1)*N3] - f[index+i*N3])/(a_layers_coordinate[i+1]-a_layers_coordinate[i])*(val3 - a_layers_coordinate[i]) + f[index+i*N3];
      break;
    }
  }
  val3 = kernel_urea_old(val3,point3)/(val3*val3)*fval*funcp3;

  integral = gew[0]*val1 + gew[1]*val2 + gew[2] *val3;
  integral = 0.5*scal*integral;

  return integral;
}




double gauss_sink_urea_old(int j, int index, double *a_layers_coordinate, 
		  double a_layers_coordinate_outerLoop, double *f, int N3)
{
  double val1, val2,val3, integral, scal;

  double stuetz[3] ={-sqrt(3./5.) , 0., sqrt(3./5.)};
  double gew[3] ={5./9., 8./9., 5./9.};

  // interpolation of function values
  double funcp1=0.5*(f[index+(j+1)*N3]-f[index+j*N3])*(stuetz[0]+1)+f[index+j*N3];
  double funcp2=0.5*(f[index+(j+1)*N3]-f[index+j*N3])*(stuetz[1]+1)+f[index+j*N3];
  double funcp3=0.5*(f[index+(j+1)*N3]-f[index+j*N3])*(stuetz[2]+1)+f[index+j*N3];
  // quadrature points
  double point1 = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2*(stuetz[0]+1)+a_layers_coordinate[j];
  double point2 = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2*(stuetz[1]+1)+a_layers_coordinate[j];
  double point3 = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2*(stuetz[2]+1)+a_layers_coordinate[j];
  // transform of variables
  scal=(a_layers_coordinate[j+1]-a_layers_coordinate[j])/2;

  //calculate sin term
  val1=kernel_urea_old(a_layers_coordinate_outerLoop, point1)*funcp1;
  val2=kernel_urea_old(a_layers_coordinate_outerLoop, point2)*funcp2;
  val3=kernel_urea_old(a_layers_coordinate_outerLoop, point3)*funcp3;

  integral= gew[0]*val1 + gew[1]*val2 + gew[2] *val3;
  integral =scal*integral;

  return integral;
}




void aggregation_conv_urea_old(TCollection *Coll,int N_x, int N_y, int N_z, int N_a, double *x_coord, double *y_coord, double *z_coord,
		      double *f, double *rhs_new, double *a_layers_coord, TFEFunction3D *temperature)
{
    double source, sink,gauss_sink,gauss_source,sum2, val, fval, val_left, val_right, length, scal;
    double L_infty = TDatabase::ParamDB->UREA_D_P_MAX;
    double f_infty = TDatabase::ParamDB->UREA_f_infty;
    double l_infty = TDatabase::ParamDB->UREA_l_infty;
    double u_infty = TDatabase::ParamDB->UREA_u_infty;

    double eps=1e-6;
    int N3 = (N_x+1)*(N_y+1)*(N_z+1);
    int x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat;
    int index, i, j, ii, N_V, k, N_Cells, beta, alpha, gamma, val0, val1;
    double x[8],y[8],z[8],values[4];
    double *velo, *grad_velo, *temp, *layers_mass;
    TBaseCell *cell;
    TVertex *vertex[8];
    temp = new double[N3];
    cout<<"start aggregation" <<endl;
    N_Cells = Coll->GetN_Cells();

    //loop over the mesh cells of the global grid
    for(i=0;i<N_Cells;i++)
    {
	cell = Coll->GetCell(i);
	N_V = cell->GetN_Vertices();
	
	//volume = cell->GetMeasure();
	if (N_V < 8)
	{
	    OutPut(" only for hexahedral mesh cells implemented !!!"<<endl);
	    exit(4711);
	}
	for (j=0;j<N_V;j++)
	{
	    // read coordinates of the mesh cell
	    vertex[j] = cell->GetVertex(j);
	    vertex[j]->GetCoords(x[j], y[j], z[j]);
	}
	
	for (j=0;j<N_V;j++)
	{
	    // compute index of psd array
	    for (k=0;k<=N_x;k++)
	    {
		if (fabs(x[j]-x_coord[k])<eps)
		{
		    beta=k;
		    break;
		}
	    }
	    val0 = (N_x+1)*(N_y+1);
	    for (k=0;k<val0;k+=N_x+1)
	    {
		if (fabs(y[j]-y_coord[k])<eps)
		{
		  alpha=k/(N_x+1);
		    break;
		}
	    }
	    val1 = (N_x+1)*(N_y+1)*(N_z+1);
	    for (k=0;k<val1;k+=val0)
	    {
	      // OutPut(z[j] << " k " << k << " " << z_coord[k] << " " << fabs(z[j]-z_coord[k]) << endl);
		if (fabs(z[j]-z_coord[k])<eps)
		{
		    gamma=k/val0;
		    break;
		}
	    }
	    // index of psd for smallest particle size
	    index=gamma*(N_x+1)*(N_y+1)+alpha*(N_x+1)+beta;
	    //OutPut(index << " " << x[j] << " " << y[j] << " " << z[j] << " ");
	    //OutPut(gamma << " " << alpha << " " << beta << endl);

            // temperature with dimension	   
	    temperature->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    temp[index] = values[0] * TDatabase::ParamDB->UREA_temp_infty;
            // boundary condition
            if (PSD_bound_cond_from_velo_inflow_urea(x[j], y[j], z[j]))
              {
                 for (ii=0;ii<=N_a;ii++)
                 {
                    f[index+ii*N3]=psd_boundary_urea(x[j], y[j], z[j], a_layers_coord[ii],
                                        1, values[0]);
                   // OutPut(ii << " aaa "<< x[j]<< " " << y[j] << " " <<  z[j]<< " : " //<<  
                   //  " a " << a_layers_coord[ii] << " " << f[index+ii*N3] << endl); 
                 }
              }
            
	    //OutPut(index << " " << temp[index] << " : ");
	}
    }
  
   
     // compute dimensional quantities
    Dscal(N_a+1, TDatabase::ParamDB->UREA_D_P_MAX,a_layers_coord);
    Dscal(N3*(N_a+1), f_infty,f);
    scal=l_infty/(u_infty*f_infty);

  // loop over all points of the flow domain
  for(index=0;index<N3;index++)
  {
    // loop over inner coordinates
    for(i=0;i<=N_a;i++)
    {
      // source term: loop over the intervals
      source = 0.;
      for (j=0;j<i;j++)
      {
        source += gauss_source_urea_old(j, index, a_layers_coord, a_layers_coord[i], f, N3, N_a);
      }
      source=a_layers_coord[i]*a_layers_coord[i]*source;
      
      // sink term: loop over the intervals
      sink=0.;
      for(j=0;j<N_a;j++)
      {
        sink += gauss_sink_urea_old(j, index, a_layers_coord, a_layers_coord[i], f, N3);
	 
      }
      
      sink=f[index+i*N3]*sink;
    
            // new rhs 
      rhs_new[i*N3+index] = scal*(source - sink);
      //OutPut(index << " " << i << " " << sink << " " << source << " " << rhs_new[i*N3+index] << endl);
    // if ((isnan(rhs_new[i*N3+index]))||(isinf(rhs_new[i*N3+index])))
     //       exit(1);	 
    }
    //OutPut(" sink: " << sink  <<endl);
    
  }
   
  // loop through all points of the flow domain

   // OutPut(endl);
    // dimensionless quantities
    Dscal(N_a+1, 1.0/TDatabase::ParamDB->UREA_D_P_MAX,a_layers_coord);
    Dscal(N3*(N_a+1), 1.0/f_infty,f);
    
    //Dscal(N3*(N_a+1), scal, rhs_new);
     
    OutPut("rhs from aggregation " << sqrt(Ddot(N3*(N_a+1),rhs_new,rhs_new)) << endl);

    //for (ii=0;ii<=N_a;ii++)
   // {
	//OutPut(ii << " coor " << a_layers_coord[ii] << " result " << l_infty*f_infty*rhs_new[index+ii*N3]/u_infty << endl);
       
   // } 

 //exit(1);
	


cout<<"aggregation done" <<endl;
 delete temp;
  
}

/****************************************************************************************/
//
// kernel for agglomeration
//
/***************************************************************************************/

double kernel_urea(double x, double x_bar, double shear_rate)
{
  double kernel, val;
  double c_agg = TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR;
  double b_aggr = TDatabase::ParamDB->UREA_AGGR_BROWNIAN;
  //double c_agg = TDatabase::ParamDB->WINDTUNNEL_SHEAR_FACTOR;
  //double visc = TDatabase::ParamDB->WINDTUNNEL_dynamic_viscosity;
  //const double kBoltz = 1.3806504e-23;
  //const double temp = 293.15;
  //double scal = 2. * kBoltz * temp /(3*visc);
  double scal;
  double eps=1.e-8;
  double k_B = 1.3806504e-23;
  double mu = 1.074e-3;
  double T=298.;
  scal=b_aggr*2.*k_B*T/(3*mu);
  val = x+x_bar;
  // part coming from Brownschen motion
  if ((x>0) && (x_bar>0))
    kernel = scal*val*(1/x+1/x_bar);
  else
  {
    kernel = 0;
  }
  // add part coming from shear stress
   //kernel = kernel+ c_agg * shear_rate * val * val * val;
  kernel = kernel+ c_agg * val * val * val;

  return(kernel);
}


/****************************************************************************************/
//
// compute source term for aggregation using Gaussian quadrature
//
/***************************************************************************************/

double gauss_source_urea(int j, int index, double *a_layers_coordinate,
double a_layers_coordinate_outerLoop, double *f, int N3, int Na, double shear_rate)
{
    int i, iter, found;
  double val, val1, val2, val3, val_outer, integral, source, scal, fval, tmp;
  //double stuetz[3] ={-sqrt(3./5.)+1 , 1., sqrt(3./5.)+1};
  double stuetz[3] ={0.225403330758517,1.0,1.774596669241483};
  //double gew[3] ={5./9., 8./9., 5./9.};
  double gew[3] ={0.555555555555556,0.888888888888889,0.555555555555556};
  double funcp1, funcp2, funcp3, point1, point2, point3;

  // interpolation of function values
  tmp = 0.5 * (f[index+(j+1)*N3] - f[index+j*N3]);
  funcp1 = tmp * stuetz[0] + f[index+j*N3];
  funcp2 = tmp * stuetz[1] + f[index+j*N3];
  funcp3 = tmp * stuetz[2] + f[index+j*N3];
  // transform of variables
  scal = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2.0;
  // quadrature points
  point1 = scal * stuetz[0] + a_layers_coordinate[j];
  point2 = scal * stuetz[1] + a_layers_coordinate[j];
  point3 = scal * stuetz[2] + a_layers_coordinate[j];
 
  // volume that corresponds to the outer coordinate
  val_outer= a_layers_coordinate_outerLoop*a_layers_coordinate_outerLoop*a_layers_coordinate_outerLoop;

  // slope
  tmp = (f[index+(j+1)*N3] - f[index+j*N3])/(a_layers_coordinate[j+1]-a_layers_coordinate[j]);
  // calculate source term
  // third quadrature point
  val=val_outer-(point3*point3*point3);
  val3=pow(val,1.0/3.0);
  fval = -1;
  found = 0;
  // find PSD for first factor in integral
  for(i=1;i<=Na;i++)
  {
    if (a_layers_coordinate[i] >= val3)
    {
	// lower bound of coordinate
	iter = i;
	found = 1;
      fval = tmp*(val3 - a_layers_coordinate[i]) + f[index+i*N3];
      break;
    }
  }
  if (!found)
  {
      OutPut("val3 not found" << endl);
      exit(1);
  }
  // correct round-off errors
  if (fval < 0)
      fval = 0;
  val3 = kernel_urea(val3,point3, shear_rate)/(val3*val3)*fval*funcp3;

  // second quadrature point
  val=val_outer-(point2*point2*point2);
  val2=pow(val,1.0/3.0);
  fval = -1;
  found = 0;
  for(i=iter;i<=Na;i++)
  {
    if (a_layers_coordinate[i] >= val2)
    {
      // find PSD for first factor in integral
	iter = i;
	found = 1;
      fval = tmp*(val2 - a_layers_coordinate[i]) + f[index+i*N3];
      break;
    }
  }
  // correct round-off errors
  if (fval < 0)
      fval = 0;
  if (!found)
  {
      OutPut("val2 not found" << endl);
      exit(1);
  }

  val2 = kernel_urea(val2,point2, shear_rate)/(val2*val2)*fval*funcp2;

 // first quadrature point
  val=val_outer-(point1*point1*point1);
  val1=pow(val,1.0/3.0);
  fval = -1;
  found = 0;
  for(i=iter;i<=Na;i++)
  {
    if (a_layers_coordinate[i] >= val1)
    {
	found = 1;
      // find PSD for first factor in integral
      fval = tmp*(val1 - a_layers_coordinate[i]) + f[index+i*N3];
      break;
    }
  }
  // correct round-off errors
  if (fval < 0)
      fval = 0;
  if (!found)
  {
      OutPut("val1 not found" << endl);
      exit(1);
  }

  // apply the kernel 
  val1 = kernel_urea(val1,point1, shear_rate)/(val1*val1)*fval*funcp1;
 
  // integral by Gaussian quadrature rule
  integral = gew[0]*val1 + gew[1]*val2 + gew[2] *val3;
  integral = scal*integral;

  return integral;
}


/****************************************************************************************/
//
//  compute sink term for aggregation using Gaussian quadrature
//
/***************************************************************************************/
double gauss_sink_urea(int j, int index, double *a_layers_coordinate,
double a_layers_coordinate_outerLoop, double *f, int N3, double shear_rate)
{
  double val1, val2,val3, integral, scal,tmp;
  //double stuetz[3] ={-sqrt(3./5.)+1 , 1.0, sqrt(3./5.)+1};
  //double gew[3] ={5./9., 8./9., 5./9.};
  double funcp1, funcp2, funcp3, point1, point2, point3;
  double stuetz[3] ={0.225403330758517,1.0,1.774596669241483};
  double gew[3] ={0.555555555555556,0.888888888888889,0.555555555555556};
 
  // interpolation of function values
  tmp = 0.5*(f[index+(j+1)*N3]-f[index+j*N3]);
  funcp1= tmp * stuetz[0]+f[index+j*N3];
  funcp2= tmp * stuetz[1]+f[index+j*N3];
  funcp3= tmp * stuetz[2]+f[index+j*N3];
  // transform of variables
  scal = (a_layers_coordinate[j+1]-a_layers_coordinate[j])/2.0;
  // quadrature points
  point1 = scal * stuetz[0]+a_layers_coordinate[j];
  point2 = scal * stuetz[1]+a_layers_coordinate[j];
  point3 = scal * stuetz[2]+a_layers_coordinate[j];

  //calculate sink term
  val1=kernel_urea(a_layers_coordinate_outerLoop, point1,shear_rate)*funcp1;
  val2=kernel_urea(a_layers_coordinate_outerLoop, point2,shear_rate)*funcp2;
  val3=kernel_urea(a_layers_coordinate_outerLoop, point3,shear_rate)*funcp3;
  // Gaussian quadrature
  integral = gew[0]*val1 + gew[1]*val2 + gew[2] *val3;
  integral = scal*integral;

  return integral;
}


/****************************************************************************************/
//
// compute aggregation integrals with Gaussian quadrature
//
/***************************************************************************************/
void aggregation_conv_urea(TCollection *Coll,
double *x_coord, double *y_coord, double *z_coord,
int N_x, int N_y, int N_z, int N_a,
double *f, double *rhs_new, double *a_layers_coord,TFEFunction3D *temperature)
{
  double source, sink, scal;
  //double r_infty = TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
  //double f_infty = TDatabase::ParamDB->WINDTUNNEL_F_INFTY;
  //double l_infty = TDatabase::ParamDB->WINDTUNNEL_L_INFTY;
  //double u_infty = TDatabase::ParamDB->WINDTUNNEL_U_INFTY;
  double L_infty = TDatabase::ParamDB->UREA_D_P_MAX;
  double f_infty = TDatabase::ParamDB->UREA_f_infty;
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
  double u_infty = TDatabase::ParamDB->UREA_u_infty;
  int N3=(N_x+1)*(N_y+1)*(N_z+1);

  int x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat;
  int index, i, j, ii, k, N_V, N_Cells, beta, alpha, gamma, val0, val1;
  double x[8],y[8],z[8],values[4];
  double eps=1e-6;
  double grad_velo[9], *shear_rate, *temp, fac;
  TBaseCell *cell;
  TVertex *vertex[8];

  OutPut("start aggregation" <<endl);
  // allocate arrays
 // shear_rate = new double[N3];

  N_Cells = Coll->GetN_Cells();
  temp = new double[N3];
  fac = u_infty/l_infty;
  //loop over the mesh cells of the global grid
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();

    //volume = cell->GetMeasure();
    if (N_V < 8)
    {
      OutPut(" only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex[j] = cell->GetVertex(j);
      vertex[j]->GetCoords(x[j], y[j], z[j]);
    }

    for (j=0;j<N_V;j++)
    {
      // compute index of psd array
      for (k=0;k<=N_x;k++)
      {
        if (fabs(x[j]-x_coord[k])<eps)
        {
          beta=k;
          break;
        }
      }
      val0 = (N_x+1)*(N_y+1);
      for (k=0;k<val0;k+=N_x+1)
      {
        if (fabs(y[j]-y_coord[k])<eps)
        {
          alpha=k/(N_x+1);
          break;
        }
      }
      val1 = (N_x+1)*(N_y+1)*(N_z+1);
      for (k=0;k<val1;k+=val0)
      {
        // OutPut(z[j] << " k " << k << " " << z_coord[k] << " " << fabs(z[j]-z_coord[k]) << endl);
        if (fabs(z[j]-z_coord[k])<eps)
        {
          gamma=k/val0;
          break;
        }
      }
      // index of psd for smallest particle size
      index=gamma*(N_x+1)*(N_y+1)+alpha*(N_x+1)+beta;
      // compute gradient of velocity
      // compute shear rate
      
     //  shear_rate[index]=0;
         // first component
     // velocity1->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
     // shear_rate[index] += values[1] * values[1];
     // shear_rate[index] += values[2] * values[2];
     // shear_rate[index] += values[3] * values[3];
         // second component
     // velocity2->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
     // shear_rate[index] += values[1] * values[1];
     // shear_rate[index] += values[2] * values[2];
     // shear_rate[index] += values[3] * values[3];
         // third component
     // velocity3->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
     // shear_rate[index] += values[1] * values[1];
     // shear_rate[index] += values[2] * values[2];
     // shear_rate[index] += values[3] * values[3];
      //shear_rate[index] = fac * sqrt(2*shear_rate[index]);
      
       // temperature with dimension	   
      temperature->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
      temp[index] = values[0] * TDatabase::ParamDB->UREA_temp_infty;
      
      
         // boundary condition
            if (PSD_bound_cond_from_velo_inflow_urea(x[j], y[j], z[j]))
              {
                 for (ii=0;ii<=N_a;ii++)
                 {
                    f[index+ii*N3]=psd_boundary_urea(x[j], y[j], z[j], a_layers_coord[ii],
                                        1, values[0]);
                   // OutPut(ii << " aaa "<< x[j]<< " " << y[j] << " " <<  z[j]<< " : " //<<  
                   //  " a " << a_layers_coord[ii] << " " << f[index+ii*N3] << endl); 
                 }
              }
            
    }
  }
 

  // PSD with dimension is needed
 // Dscal(N3 * (N_a+1), f_infty, f);
 //Dscal(N_a+1, r_infty,  a_layers_coord);
  Dscal(N_a+1, TDatabase::ParamDB->UREA_D_P_MAX,a_layers_coord);
  Dscal(N3*(N_a+1), f_infty,f);
  // factor for non-dimensionless
  scal=l_infty/(u_infty*f_infty);

  // loop over all points of the flow domain
  for(index=0;index<N3;index++)
  {
    // loop over inner coordinates
    for(i=0;i<=N_a;i++)
    {
      // source term: loop over the intervals
      source = 0.;
      for (j=0;j<i;j++)
      {
        source += gauss_source_urea(j, index, a_layers_coord, a_layers_coord[i], f, N3, N_a, 1.);
      }
      source=a_layers_coord[i]*a_layers_coord[i]*source/2.0;

      // sink term: loop over the intervals
      sink=0.;
      for(j=0;j<N_a;j++)
      {
        sink += gauss_sink_urea(j, index, a_layers_coord, a_layers_coord[i], f, N3, 1.);
      }

      sink=f[index+i*N3]*sink;

      // new rhs
      rhs_new[i*N3+index] = scal*(source - sink);
    }
  }

  // dimensionless quantities
//  Dscal(N3 * (N_a+1), 1.0/f_infty, f);
// Dscal(N_a+1, 1.0/r_infty,  a_layers_coord);
    Dscal(N_a+1, 1.0/TDatabase::ParamDB->UREA_D_P_MAX,a_layers_coord);
    Dscal(N3*(N_a+1), 1.0/f_infty,f);
    
   // delete[] shear_rate;
   delete[] temp;
  OutPut("aggregation done " << sqrt(Ddot(N3 * (N_a+1), rhs_new, rhs_new)) << endl);
} 



void call_apply_integral_operator(int nx, int ny, int nz, int na, double* input, double* v, 
                                  double* grad_v, double* temp, double* output, 
                                  double* grid, double* params)
{
  
            apply_integral_operator(nx, ny, nz, na,
                                    input,
                                    v, grad_v, temp,
                                    output, grid,
                                    TDatabase::ParamDB->UREA_D_P_MAX,
                                    TDatabase::ParamDB->UREA_f_infty,
                                    params);
}
  



/*double kernel_urea(double x, double x_bar)
{
double kernel;
if(x>1.250000e-04 && x_bar>1.250000e-04)
{
kernel=0.4*(x+x_bar)*(x+x_bar)*(x+x_bar);
}
else
{
return(0.0);
}
return(kernel);
}

void aggregation_conv_urea(TCollection *Coll,int N_x, int N_y, int N_z, int N_a, double *x_coord, double *y_coord, double *z_coord,
			double *f, double *rhs_new, double *a_layers_coord, TFEFunction3D *temperature)
{
double source, sink;

int N3 = (N_x+1)*(N_y+1)*(N_z+1);
int x_coord_mat, y_coord_mat, z_coord_mat, a_coord_mat;
int iter, i, j, ii, N_V, k, index, N_Cells, beta, alpha, gamma, val0, val1;
    double eps=1e-6, fac;
    double x[8],y[8],z[8],values[4];
    
    double *velo, *grad_velo, *temp, *layers_mass;
    double f_infty = TDatabase::ParamDB->UREA_f_infty;
    double l_infty = TDatabase::ParamDB->UREA_l_infty;
    double u_infty = TDatabase::ParamDB->UREA_u_infty;
    double L_max = TDatabase::ParamDB->UREA_D_P_MAX;
   
    TBaseCell *cell;
    TVertex *vertex[8];
   temp = new double[N3];


    cout<<"start aggregation" <<endl;
    N_Cells = Coll->GetN_Cells();
    //loop over the mesh cells of the global grid
    for(i=0;i<N_Cells;i++)
    {
	cell = Coll->GetCell(i);
	N_V = cell->GetN_Vertices();
	
	//volume = cell->GetMeasure();
	if (N_V < 8)
	{
	    OutPut(" only for hexahedral mesh cells implemented !!!"<<endl);
	    exit(4711);
	}
	for (j=0;j<N_V;j++)
	{
	    // read coordinates of the mesh cell
	    vertex[j] = cell->GetVertex(j);
	    vertex[j]->GetCoords(x[j], y[j], z[j]);
	}
	
	for (j=0;j<N_V;j++)
	{
	    // compute index of psd array
	    for (k=0;k<=N_x;k++)
	    {
		if (fabs(x[j]-x_coord[k])<eps)
		{
		    beta=k;
		    break;
		}
	    }
	    val0 = (N_x+1)*(N_y+1);
	    for (k=0;k<val0;k+=N_x+1)
	    {
		if (fabs(y[j]-y_coord[k])<eps)
		{
		  alpha=k/(N_x+1);
		    break;
		}
	    }
	    val1 = (N_x+1)*(N_y+1)*(N_z+1);
	    for (k=0;k<val1;k+=val0)
	    {
	      // OutPut(z[j] << " k " << k << " " << z_coord[k] << " " << fabs(z[j]-z_coord[k]) << endl);
		if (fabs(z[j]-z_coord[k])<eps)
		{
		    gamma=k/val0;
		    break;
		}
	    }
	    // index of psd for smallest particle size
	    index=gamma*(N_x+1)*(N_y+1)+alpha*(N_x+1)+beta;
	    //OutPut(index << " " << x[j] << " " << y[j] << " " << z[j] << " ");
	    //OutPut(gamma << " " << alpha << " " << beta << endl);
	   
	    temperature->FindGradientLocal(cell,i,x[j],y[j],z[j],values);
	    temp[index] = values[0] * TDatabase::ParamDB->UREA_temp_infty;
            {
              if (PSD_bound_cond_from_velo_inflow_urea(x[j], y[j], z[j]))
              {
                 for (ii=0;ii<=N_a;ii++)
                 {
                    f[index+ii*N3]=psd_boundary_urea(x[j], y[j], z[j], a_layers_coord[ii],
                                        1, values[0]);
                    //OutPut(ii << " aaa "<< x[j]<< " " << y[j] << " " <<  z[j]<< " : " <<  
                    // " a " << a_layers_coord[ii] << " " << f[index+ii*N3] << endl); 
                 }
              }
            }
	    //OutPut(index << " " << temp[index] << " : ");
	}
    }
  



for(iter=0;iter<=N3;iter++)
{
 for(i=0;i<=N_a;i++)
  { 
     source= 0.5 * (a_layers_coord[1]-a_layers_coord[0])
               * (kernel_urea(a_layers_coord[i]-a_layers_coord[0],a_layers_coord[0])
               * f[i*N3+iter]*f[iter]);
     sink= 0.5 * (a_layers_coord[1]-a_layers_coord[0])
               * kernel_urea(a_layers_coord[i],a_layers_coord[0])
               * f[iter]
         + 0.5 * (a_layers_coord[N_a]-a_layers_coord[N_a-1])*kernel_urea(a_layers_coord[i],a_layers_coord[N_a])
               * f[N_a*N3+iter];
     for(j=1;j<i;j++)
     {
       if(j==N_a) continue;
       
       source += 0.5*(a_layers_coord[j+1]-a_layers_coord[j-1])
            * kernel_urea(a_layers_coord[i]-a_layers_coord[j],a_layers_coord[j])
            * f[(i-j)*N3+iter]*f[j*N3+iter];
       if(i==N_a)
         {
          source += 0.5*(a_layers_coord[N_a]-a_layers_coord[N_a-1])
               * (kernel_urea(a_layers_coord[i]-a_layers_coord[N_a],a_layers_coord[N_a])
               * f[(i-N_a)*N3+iter]*f[N_a*N3+iter]);
         }
    
     }
    for(j=1;j<N_a;j++)
    {
     sink += 0.5*(a_layers_coord[j+1]-a_layers_coord[j-1])
          * kernel_urea(a_layers_coord[i],a_layers_coord[j])
          * f[j*N3+iter];
    }
  sink=f[i*N3+iter]*sink;
  rhs_new[i*N3+iter]=0.04*(.5*source+sink);
//cout<<"rhs:" <<rhs_new[i*N3+iter] << endl;
  }

 
}

cout<<"aggregation done" <<endl;
}
*/
