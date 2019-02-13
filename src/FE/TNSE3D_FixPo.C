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
   
// ======================================================================
// @(#)TNSE3D_FixPo.C        1.3 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>
#include <stdlib.h>
#include <TNSE3D_Routines.h>
#include <MainUtilities.h>
// ======================================================================
// compute turbulent viscosity for LES 
// ======================================================================
double TurbulentViscosity3D(double delta, double* gradU, double* u, 
			    double* uConv, double* x, double* y, double* z, double proj_space)
{
  int nu_type = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
  double nu_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  int nu_tensor =  TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
  int found;
  double Re, zplus, r, eps = 1e-6, lambda, val, x0, y0;
  double c, Alpha, b11, b22, b33, b12, b13, b23, Bbeta;
  double nu_power, nu_sigma;
  double delta_x, delta_y, delta_z, hk;
  double mu_max, invariant_2, invariant_3;
  double frobenius_norm_tensor,nu,a11,a12,a13,a22,a23,a33,sigma;
  
  // van Driest damping, Do 09.02.06
  double A =26.0;
  //OutPut(proj_space << " ");
  // no turbulent viscosity
  if ((proj_space==2)||(nu_type==0))
  {
    nu = 0;
    //OutPut("ps "  << proj_space << " ");
    return(nu);
  }

  // compute square of the Frobenius norm of the tensor
  // use deformation tensor
  switch(nu_tensor)
  {
    case 0:
      // compute (grad(u)+grad(u)^T)/2
      
      a11 = gradU[0]+gradU[0];
      a12 = gradU[1]+gradU[3];
      a13 = gradU[2]+gradU[6];
      a22 = gradU[4]+gradU[4];
      a23 = gradU[5]+gradU[7];
      a33 = gradU[8]+gradU[8];
      frobenius_norm_tensor = 2*(a12*a12 + a13*a13 + a23*a23);
      frobenius_norm_tensor += a11*a11 + a22*a22 + a33*a33;
      frobenius_norm_tensor /= 4;
    
      break;

    case 1:
      // use grad u
      frobenius_norm_tensor =  gradU[0]*gradU[0] + gradU[1]*gradU[1] +
      gradU[2]*gradU[2] + gradU[3]*gradU[3]  + gradU[4]*gradU[4] + 
      gradU[5]*gradU[5] + gradU[6]*gradU[6] + gradU[7]*gradU[7] +
      gradU[8]*gradU[8];
  
      break;

    case 2:
    // deformation tensor of small scales
    // compute (grad(u)+grad(u)^T)/2 - G^H
    // works only with VMS methods
    if (uConv == NULL)
    {
       OutPut("TURBULENT_VISCOSITY_TENSOR 2 works only with VMS methods !!!" << endl);
       exit(4711);
    }
    a11 = (gradU[0]+gradU[0])/2.0 - uConv[0];
    a12 = (gradU[1]+gradU[3])/2.0 - uConv[1];
    a13 = (gradU[2]+gradU[6])/2.0 - uConv[2];
    a22 = (gradU[4]+gradU[4])/2.0 - uConv[3];
    a23 = (gradU[5]+gradU[7])/2.0 - uConv[4];
    a33 = (gradU[8]+gradU[8])/2.0 - uConv[5];
    frobenius_norm_tensor = 2*(a12*a12 + a13*a13 + a23*a23);
    frobenius_norm_tensor += a11*a11 + a22*a22 + a33*a33;
    
    break;
    
    default:
      OutPut("TURBULENT_VISCOSITY_TENSOR " << TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR  <<
             " not implemented !!!" << endl);
      exit(4711);
  }

  // compute turbulent viscosity
  switch(nu_type)
  {
  case 1: // Smagorinsky
  case 17: // Smagorinsky
    nu = nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
    break;
  case 2: // p laplacian
    nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
    nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power/2.0);
    break;
  case 3: // Layton SIAM J. Sci. Comput. 1996 
    // nu = nu_0 |ln(h)|^(-2/3*(p-1)) h^sigma \|\nabla u\|^(p-2)
    nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
    nu_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
    nu = nu_constant * pow(delta,nu_sigma) * 
      pow(frobenius_norm_tensor,(nu_power-2)/2.0) 
      /pow(fabs(log(delta)),2*(nu_power-1)/3);
    break;
  case 4: // \|u-g_\delta \ast u\|_2
    nu = nu_constant * delta * sqrt((u[0]-uConv[0])*(u[0]-uConv[0])+
                                    (u[1]-uConv[1])*(u[1]-uConv[1])
                                    +(u[2]-uConv[2])*(u[2]-uConv[2]));
//    OutPut(" nu " << nu << " " << fabs((u[0]-uConv[0])) << " " <<  fabs((u[1]-uConv[1]))
//           << " " << fabs((u[2]-uConv[2])) << endl);
    break;
  case 5: // \|Du^h-G^H\|_2
      // the values of the off diagonal entries of G^H has to be divided by 2
      // since the basis functions are l^H/2 !!!
      a11 = gradU[0] - uConv[0];
      a12 = (gradU[1]+gradU[3]-uConv[1])/2.0;
      a13 = (gradU[2]+gradU[6]-uConv[2])/2.0;
      a22 = gradU[4]-uConv[3];
      a23 = (gradU[5]+gradU[7]-uConv[4])/2.0;
      a33 = gradU[8]-uConv[5];
      nu = nu_constant * delta * delta * sqrt(a11*a11+2*a12*a12+2*a13*a13+a22*a22+2*a23*a23+a33*a33);
      //   OutPut(" nu " << nu);// << " " << fabs((u[0]-uConv[0])) << " " <<  fabs((u[1]-uConv[1]))
//           << " " << fabs((u[2]-uConv[2])) << endl);
    break;
  case 6: // all parameters free
    nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
    nu_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
    nu = nu_constant * pow(delta,nu_sigma) * 
      pow(frobenius_norm_tensor,nu_power/2.0) ;
    break;
      case 100: // van Driest damping for channel 
 	  Re = TDatabase::ParamDB->RE_NR;
	  // walls at z=0 and z=2
	  zplus = Re*(1-fabs(1-z[0]));
	  if(zplus < 5)
	    {
	      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
		(1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
	    }
	  else
	    {
	      nu  =  nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
	    }
	  break;
      case 101: // van Driest damping for cylinder with squared cross--section
	  // left and right wall at the cylinder
	  zplus = 1000;
	  if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
	  {
	      // distance to the wall
	      if (y[0] > 0.7)
		  y0 = y[0] - 0.75;
	      else
		  y0 = 0.65 - y[0];
	      // wall units 
	      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
	  }
	  if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
	  {
	      // distance to the wall
	      if (x[0] < 0.5)
	      {
		  x0 = 0.45 - x[0];
		  // wall units 
		  zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
	      }
	      else
	      {
		  x0 = x[0] - 0.55;
		  // wall units 
		  zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
	      }
	  }
	  
	  if(zplus < 5)
	    {
	      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
		(1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
	      //OutPut("nu " << x[0] << " " << y[0] << " " << x0 << " " << y0 << " " << zplus << endl);
	    }
	  else
	    {
	      nu  =  nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
	    }
	  /*    
	      
	  r = sqrt(x0*x0+y0*y0);
	  // find the distance of the center to the boundary of the cylinder
	  // in the direction of (x[0],y[0])
	  // compute the parameter of the line from the center to (x[0],y[0])
	  // which determines the intersection with the boundary of the cylinder
	  found = 0;
	  if (fabs(x0)>eps)
	  {
	      lambda = 0.05/x0;
	      val = 0.5 + lambda * y0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(lambda>0))
		  found = 1;
	  }
	  if ((fabs(x0)>eps)&&(!found))
	  {
	      lambda = -0.05/x0;
	      val = 0.5 + lambda * y0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(lambda>0))
		  found = 1;
	  }
	  if ((fabs(y0)>eps)&&(!found))
	  {
	      lambda = 0.05/y0;
	      val = 0.5 + lambda * x0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(lambda>0))
		  found = 1;
	  }
	  if ((fabs(y0)>eps)&&(!found))
	  {
	      lambda = -0.05/y0;
	      val = 0.5 + lambda * x0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(!found))
		  found = 1;
	  }
	  if (!found)
	  {
	      OutPut("intersection not found !" << endl);
	      exit(4711);
	  }
	  // compute distance
	  r = (1-lambda)*r;
	  //OutPut("aft " << x[0] << " " << y[0] << " " << r <<endl);
 	  Re = TDatabase::ParamDB->RE_NR;
	  // THIS HAS TO BE CORRECTED
	  zplus = Re*r;
	  if(zplus < 5)
	    {
	      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
		(1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
	    }
	  else
	    {
	      nu  =  nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
	      }*/
	  break;
  case 102:                                     // van Driest damping for channel flow (continuous)
      Re = TDatabase::ParamDB->RE_NR;
      // walls at z=0 and z=2
      zplus = Re*(1-fabs(1-z[0]));
      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
        (1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);

      /*OutPut("Van Driest: " << (1-exp(-zplus/A)) *
      (1-exp(-zplus/A)) << endl);*/
      break;
    case 103:                                     // van Driest damping for channel flow (paper: Rudman, Blackburn'99)
      Re = TDatabase::ParamDB->RE_NR;
      // walls at z=0 and z=2
      zplus = Re*(1-fabs(1-z[0]));
      nu = nu_constant * delta * delta * (1-exp(-(zplus/A)*(zplus/A)*(zplus/A))) *  sqrt(frobenius_norm_tensor);

      /*OutPut("Van Driest: " << (1-exp(-zplus/A)) *
      (1-exp(-zplus/A)) << endl);*/
      break;
    case 104:                                     // van Driest damping for channel flow (paper: Rudman, Blackburn'99) with diff A+
      Re = TDatabase::ParamDB->RE_NR;
      A =17.0;
      // walls at z=0 and z=2
      zplus = Re*(1-fabs(1-z[0]));
      nu = nu_constant * delta * delta * (1-exp(-(zplus/A)*(zplus/A)*(zplus/A))) *  sqrt(frobenius_norm_tensor);

      /*OutPut("Van Driest: " << (1-exp(-zplus/A)) *
      (1-exp(-zplus/A)) << endl);*/
      break;

    case 105:
      // eddy viscosity model: Vreman, Phys. Fluids 16 (10), 3670 -3681, 2004
      // frobenius norm of gradient of velocity
      // use same notations as in paper
      
      Alpha = gradU[0]*gradU[0]+gradU[1]*gradU[1]+gradU[2]*gradU[2]+gradU[3]*gradU[3]
              +gradU[4]*gradU[4]+gradU[5]*gradU[5]
              +gradU[6]*gradU[6]+gradU[7]*gradU[7]+gradU[8]*gradU[8];
      if (fabs(Alpha)<1e-12)
      {
        nu = 0;
        break;
      }

      // compute filter width in coordinate directions
      // hk is just the value if the filter width would be zero (this should not happen)
      
      hk = delta/2.0;
      delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
      delta_x *=delta_x;
      delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
      delta_y *=delta_y;
      delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
      delta_z *=delta_z;

      // compute second invariant of gradient of velocity, scaled with filter widht in coordinate directions
      b11 = delta_x*gradU[0]*gradU[0]+delta_y*gradU[3]*gradU[3]+delta_z*gradU[6]*gradU[6];
      b22 = delta_x*gradU[1]*gradU[1]+delta_y*gradU[4]*gradU[4]+delta_z*gradU[7]*gradU[7];
      b33 = delta_x*gradU[2]*gradU[2]+delta_y*gradU[5]*gradU[5]+delta_z*gradU[8]*gradU[8];
      b12 = delta_x*gradU[0]*gradU[1]+delta_y*gradU[3]*gradU[4]+delta_z*gradU[6]*gradU[7];
      b13 = delta_x*gradU[0]*gradU[2]+delta_y*gradU[3]*gradU[5]+delta_z*gradU[6]*gradU[8];
      b23 = delta_x*gradU[1]*gradU[2]+delta_y*gradU[4]*gradU[5]+delta_z*gradU[7]*gradU[8];
      Bbeta = b11*b22 - b12*b12 +b11*b33 - b13*b13 + b22*b33 - b23*b23;
      // check for round-off errors
      if (Bbeta<0)
      {
        Bbeta = 0;
      }

      // scale, in Vreman (2004) it is recommended to use 2.5 times Smagorinsky constant
      nu = nu_constant*sqrt(Bbeta/Alpha);

      break;
    case 106:                                     // van Driest damping (continuous, classical) for cylinder with squared cross--section
      // left and right wall at the cylinder
      zplus = 1000;
      if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
      {
        // distance to the wall
        if (y[0] > 0.7)
          y0 = y[0] - 0.75;
        else
          y0 = 0.65 - y[0];
        // wall units
        zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
      }
      if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
      {
        // distance to the wall
        if (x[0] < 0.5)
        {
          x0 = 0.45 - x[0];
          // wall units
          zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
        }
        else
        {
          x0 = x[0] - 0.55;
          // wall units
          zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
        }
      }
      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
        (1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
      //OutPut("nu " << x[0] << " " << y[0] << " " << x0 << " " << y0 << " " << zplus << endl);
      break;
    case 107:                                     // van Driest damping (paper: Rudman, Blackburn'99) for cylinder with squared cross--section
      // left and right wall at the cylinder
      zplus = 1000;
      if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
      {
        // distance to the wall
        if (y[0] > 0.7)
          y0 = y[0] - 0.75;
        else
          y0 = 0.65 - y[0];
        // wall units
        zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
      }
      if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
      {
        // distance to the wall
        if (x[0] < 0.5)
        {
          x0 = 0.45 - x[0];
          // wall units
          zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
        }
        else
        {
          x0 = x[0] - 0.55;
          // wall units
          zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
        }
      }
      nu = nu_constant * delta * delta * (1-exp(-(zplus/A)*(zplus/A)*(zplus/A))) *  sqrt(frobenius_norm_tensor);
      //OutPut("nu " << x[0] << " " << y[0] << " " << x0 << " " << y0 << " " << zplus << endl);
      break;
      
    case 108: /** Verstappen model (J Sci Comput'11) */
      
      /* C = 1/mu_max as on page 107 */
      
      // compute filter width in coordinate directions
      // hk is just the value if the filter width would be zero (this should not happen)
      hk = delta/2.0;
      delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
      delta_x *= delta_x;
      delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
      delta_y *= delta_y;
      delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
      delta_z *= delta_z;
      
      mu_max = 4 * ( 1./delta_x + 1./delta_y + 1./delta_z );
      
      switch(nu_tensor)
      {
  case 0:
    a11 = a11/2.;
    a12 = a12/2.;
    a13 = a13/2.;
    a22 = a22/2.;
    a23 = a23/2.;   
    a33 = a33/2.; 
    break;
    
  case 1:
    OutPut("ERROR: Verstappen model needs a symmetric stress tensor!" << endl);
    exit(0);
      }
      
      //invariant_3 = - det(D(u))
      invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 - a23*a23*a11 - a13*a13*a22);
      invariant_2 = 0.5 * frobenius_norm_tensor;
      
      nu = (1.5 * fabs(invariant_3) ) / (mu_max * invariant_2);      
      break;
      
    case 109:  /** Verstappen model (J Sci Comput'11) */
      
      /* C = (h/pi)^2 as on page as on page 97, where Delta = (h_x*h_y*h_z)^(1/3) */
      
      // compute filter width in coordinate directions
      // hk is just the value if the filter width would be zero (this should not happen)
      hk = delta/2.0;
      delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
      delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
      delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
      
      /* TODO: change cell width hk using CELL_MEASURE (more elegant), now too slow! */
      
      hk = delta_x*delta_y*delta_z;
      hk = pow(hk,1.0/3.0);
      
      switch(nu_tensor)
      {
  case 0:
    a11 = a11/2.;
    a12 = a12/2.;
    a13 = a13/2.;
    a22 = a22/2.;
    a23 = a23/2.;   
    a33 = a33/2.; 
    break;
    
  case 1:
    OutPut("ERROR: Verstappen model needs a symmetric stress tensor!" << endl);
    exit(0);
      }
      
      /* invariant_3 = - det(D(u)) */
      invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 - a23*a23*a11 - a13*a13*a22);
      invariant_2 = 0.5 * frobenius_norm_tensor;
      
      nu = ( 1.5 * hk * hk * fabs(invariant_3) )  / ( Pi * Pi * invariant_2 );      
      
      break;
  default:
    OutPut("This type of turbulent viscosity is not implemented !!!" << endl);
    exit(4711);
  break;
  }
  
  return(nu);
}

// ======================================================================
// compute stabilization for div--div term
// ======================================================================
double DivDivStab3D(double u1, double u2, double u3, double hK, double eps) 
{
  int divdiv_type = TDatabase::ParamDB->DIV_DIV_STAB_TYPE;
  double c_1 = TDatabase::ParamDB->DIV_DIV_STAB_C1;
  double c_2 = TDatabase::ParamDB->DIV_DIV_STAB_C2;
  double tau;
    
  switch(divdiv_type)
    {
      // constant
    case 0:
      tau = c_2;
      break;
      // for non inf-sup stable fe, Codina, Gravemeier
    case 1:
      tau = c_2*sqrt(u1*u1+u2*u2+u3*u3)*hK/c_1;
      tau = tau*tau + eps*eps ;
      tau = sqrt(tau);
      break;
    case 2:
      tau = sqrt(u1*u1+u2*u2+u3*u3)*hK/c_1;
      tau = tau*tau + eps*eps ;
      tau = c_2*sqrt(tau);
      break;
    default:
      OutPut("div-div stabilization " << divdiv_type << " not implemented !!!");
      exit(4711);
    }

  return(tau);
}


/******************************************************************************/
//
// computation of SUPG parameter following 
// Bazilevs, Calo, Cottrell, Hughes, Reali, Scovazzi
//
/******************************************************************************/

void SUPG_Param3D(double u1, double u2, double u3, double* coeff, double* params)
{
    double x1, x2, x3, x0, y1, y2, y3, y0, z1, z2, z3, z0, x4, y4, z4;
    double d11, d12, d13, d21, d22, d23, d31, d32, d33, delta, nu;
    double g11, g12, g13, g22, g23, g33;
    double rec_detjk, tau_c, tau_m;
    double eps  = 1e-12;
    double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;    
    double C_I = TDatabase::ParamDB->DELTA0;
    
    nu = coeff[0];             
    rec_detjk = coeff[19];
    rec_detjk = 1/rec_detjk;

    // tetrahedron
    if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
    {
	x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
	y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
	z0 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[0];
	x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
	y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
	z1 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[1];
	x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
	y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
	z2 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[2];
	x3 = TDatabase::ParamDB->INTERNAL_VERTEX_X[3];
	y3 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[3];
	z3 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[3];
	
	d11 = ((y2-y0)*(z3-z0)+(y3-y0)*(z0-z2)) * rec_detjk;  //dxi/dx
	d12 = ((x3-x0)*(z2-z0)+(x2-x0)*(z0-z3)) * rec_detjk;  //dxi/dy
	d13 = ((x2-x0)*(y3-y0)+(x3-x0)*(y0-y2)) * rec_detjk;  //dxi/dz
	
	d21 = ((y3-y0)*(z1-z0)+(y1-y0)*(z0-z3)) * rec_detjk;  //deta/dx
	d22 = ((x1-x0)*(z3-z0)+(x3-x0)*(z0-z1)) * rec_detjk;  //deta/dy
	d23 = ((x3-x0)*(y1-y0)+(x1-x0)*(y0-y3)) * rec_detjk;  //deta/dz
	
	d31 = ((y1-y0)*(z2-z0)+(y2-y0)*(z0-z1)) * rec_detjk;  //dzeta/dx
	d32 = ((x2-x0)*(z1-z0)+(x1-x0)*(z0-z2)) * rec_detjk;  //dzeta/dy
	d33 = ((x1-x0)*(y2-y0)+(x2-x0)*(y0-y1)) * rec_detjk;  //dzeta/dz	
    }
    else
    {
	// hexahedron
	x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
	y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
	z0 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[0];
	x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
	y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
	z1 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[1];
	x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
	y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
	z2 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[2];
	x4 = TDatabase::ParamDB->INTERNAL_VERTEX_X[4];
	y4 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[4];
	z4 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[4];
	
	d11 = ((y1-y0)*(z4-z0)+(y0-y4)*(z1-z0)) * 0.5 * rec_detjk;  //dxi/dx
	d12 = ((x4-x0)*(z1-z0)+(x1-x0)*(z0-z4)) * 0.5 * rec_detjk;  //dxi/dy
	d13 = ((x1-x0)*(y4-y0)+(x0-x4)*(y1-y0)) * 0.5 * rec_detjk;  //dxi/dz
	
	d21 = ((y4-y0)*(z1-z2)+(y1-y2)*(z0-z4)) * 0.5 * rec_detjk;  //deta/dx
	d22 = ((x1-x2)*(z4-z0)+(x0-x4)*(z1-z2)) * 0.5 * rec_detjk;  //deta/dy
	d23 = ((x4-x0)*(y1-y2)+(x1-x2)*(y0-y4)) * 0.5 * rec_detjk;  //deta/dz
	
	d31 = ((y1-y2)*(z1-z0)+(y1-y0)*(z2-z1)) * 0.5 * rec_detjk;  //dzeta/dx
	d32 = ((x1-x0)*(z1-z2)+(x1-x2)*(z0-z1)) * 0.5 * rec_detjk;  //dzeta/dy
	d33 = ((x1-x2)*(y1-y0)+(x1-x0)*(y2-y1)) * 0.5 * rec_detjk;  //dzeta/dz
    }

    g11 = d11*d11 + d21*d21 + d31*d31;
    g12 = d11*d12 + d21*d22 + d31*d32;
    g13 = d11*d13 + d21*d23 + d31*d33;
    g22 = d12*d12 + d22*d22 + d32*d32;
    g23 = d12*d13 + d22*d23 + d32*d33;
    g33 = d13*d13 + d23*d23 + d33*d33;
    
    tau_m = g11*g11 + 2*g12*g12 + 2*g13*g13 + g22*g22 + 2*g23*g23 + g33*g33; // G:G
    //OutPut("det "  << rec_detjk << " " << tau_m << endl);
   
    tau_m *= C_I*nu*nu;
    tau_m +=  4/(time_step*time_step); 
    tau_m += u1 * (g11*u1+g12*u2+g13*u3) + u2*(g12*u1+g22*u2+g23*u3)
	+ u3*(g13*u1+g23*u2+g33*u3);
    if (tau_m < eps)
    {
	params[0] = 0;
	params[1] = 0;
	return;
    }
    tau_m = 1/sqrt(tau_m); // this is the parameter for the momentum equation
    
    tau_c = (d11+d21+d31)*(d11+d21+d31)+(d12+d22+d32)*(d12+d22+d32)
	+(d13+d23+d33)*(d13+d23+d33);
//    OutPut(" taucbef " << tau_c << " "  << d11+d21+d31 << " " << d12+d22+d32 << " " << d13+d23+d33 <<
//	   " " << rec_detjk << ":");
    tau_c *= tau_m;
    if (tau_c < eps)
    {
	tau_c = 0;
    }
    else
	tau_c = 1/tau_c;

    //  OutPut(" tauc " << tau_c << " "  << tau_m);
    params[0] = tau_m;
    params[1] = tau_c;

/*
  delta = (d11*d11+d21*d21+d31*d31)*(d11*d11+d21*d21+d31*d31) + 2*(d11*d12+d21*d22+d31*d32)*(d11*d12+d21*d22+d31*d32) +  // G:G
          (d12*d12+d22*d22+d32*d32)*(d12*d12+d22*d22+d32*d32) + 2*(d11*d13+d21*d23+d31*d33)*(d11*d13+d21*d23+d31*d33) +
          2*(d12*d13+d22*d23+d32*d33)*(d12*d13+d22*d23+d32*d33) + (d13*d13+d23*d23+d33*d33)*(d13*d13+d23*d23+d33*d33);
  delta *= C_I*nu*nu;         
  delta += 4/(time_step*time_step);  
  delta += u1*u1*(d11*d11+d21*d21+d31*d31)+2*u1*u2*(d11*d12+d21*d22+d31*d32)+2*u1*u3*(d11*d13+d21*d23+d31*d33)+
           u2*u2*(d12*d12+d22*d22+d32d32)+2*u2*u3*(d12*d13+d22*d23+d32*d33)+u3*u3*(d13*d13+d23*d23+d33d33);            // uGu
  delta = sqrt(delta);
  
  delta /= (d11+d21+d31)*(d11+d21+d31)+(d12+d22+d32)*(d12+d22+d32)+(d13+d23+d33)*;   // gg
  
  return delta;
*/
}

// ======================================================================
// compute parameter for Leray-alpha model
// ======================================================================
double LerayAlpha_Param3D(double hK) 
{
  double c = TDatabase::ParamDB->DELTA0;
    
  return(c*hK);
}


// ======================================================================
// Type 1, Standard Galerkin
// Type 1, ClassicalLES
// Type 1, GL00Convolution
// ======================================================================
void TimeNSType1Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;
    
    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
       
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;
    
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
        test001*ansatz001);

      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, Smagorinsky
// the nonlinear viscosity is treated implicitly
// ======================================================================
void TimeNSType1Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, Galdi-Layton 98 Model auxiliary problem
// ======================================================================
void TimeNSType1GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM, **AuxMatrix;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double *AuxMatrixRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, mu2, delta;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  MatrixB3 = LocMatrices[5];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
            
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, Standard Galerkin
// Type 2, ClassicalLES
// Type 2, GL00Convolution
// ======================================================================
void TimeNSType2Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  MatrixB3T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 2 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void TimeNSType2Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  MatrixB3T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  MatrixB3T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i
  
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, GL00AuxProblem
// ======================================================================
void TimeNSType2GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM, **AuxMatrix;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double *AuxMatrixRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, mu2, delta;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  MatrixB3 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];
  MatrixB3T = LocMatrices[8];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
            
    } // endfor j
    
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// Type 3, ClassicalLES, (grad u, grad v)
// Type 3, GL00Convolution, (grad u, grad v)`
// ======================================================================
void TimeNSType3Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// Type 3, ClassicalLES, D(u):D(v)
// Type 3, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = Mult*c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;

    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, mu, delta;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i      
}

// ======================================================================
// Type 3, GL00AuxProblem (grad u, grad v)
// ======================================================================
void TimeNSType3GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *AuxMatrixRow;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;  

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;
 
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType3GL00AuxProblemDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *AuxMatrixRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;
  double delta, mu2, viscosity;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
 
      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
            
   } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, ClassicalLES, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T  = LocMatrices[15];
  MatrixB2T  = LocMatrices[16];
  MatrixB3T  = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, ClassicalLES, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T = LocMatrices[15];
  MatrixB2T = LocMatrices[16];
  MatrixB3T = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void TimeNSType4Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T  = LocMatrices[15];
  MatrixB2T  = LocMatrices[16];
  MatrixB3T  = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void TimeNSType4UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T = LocMatrices[15];
  MatrixB2T = LocMatrices[16];
  MatrixB3T = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T  = LocMatrices[15];
  MatrixB2T  = LocMatrices[16];
  MatrixB3T  = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, delta;
  double u1, u2, u3, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T = LocMatrices[15];
  MatrixB2T = LocMatrices[16];
  MatrixB3T = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, GL00AuxProblem, (grad u, grad v)
// ======================================================================
void TimeNSType4GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *AuxMatrixRow;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;  

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;
 
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType4GL00AuxProblemDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *AuxMatrixRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;
  double delta, mu2, viscosity;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
 
      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;            
   } // endfor j 

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 4, VMS_Projection, D(u):D(v)
// ======================================================================

void TimeNSType4VMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22,   **Matrix_G33;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_L;
  double c0, c1, c2, c3, delta;
  double u1, u2, u3, mu, viscosity;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixL   = LocMatrices[12];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  Matrix_tilde_G11  = LocMatrices[19];
  Matrix_tilde_G22  = LocMatrices[20];
  Matrix_tilde_G33  = LocMatrices[21];
  Matrix_G11  = LocMatrices[22];
  Matrix_G22  = LocMatrices[23];
  Matrix_G33  = LocMatrices[24];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p
  Orig5 = OrigValues[5]; // l
  
  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  
  	  if(TDatabase::ParamDB->TURBULENT_MOD_TYPE ==0)
	  { 
	    delta =  CharacteristicFilterWidth(hK);
	    mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
				param[21]);
	  }
	  else if(TDatabase::ParamDB->TURBULENT_MOD_TYPE ==1)
	  {
	    
	    mu=(TDatabase::ParamDB->viscosity_max + TDatabase::ParamDB->viscosity_min)/2.0;
	    
	  }
	  else if (TDatabase::ParamDB->TURBULENT_MOD_TYPE ==2)
	  {
	    mu=TDatabase::ParamDB->viscosity_max;
	  }
	  else if(TDatabase::ParamDB->TURBULENT_MOD_TYPE ==3)
	  {
	    mu=TDatabase::ParamDB->viscosity_min;
	  }
	  
	  mu = mu/2.0;	     
	  viscosity = c0+mu;	
	

      
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];
      val1 = Mult*ansatz000;
      
      val = -val1*test100;
      MatrixRow1[j] += val;
      val = -val1*test010;
      MatrixRow2[j] += val;
      val = -val1*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];
     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig5[j];
        val =  Mult * 2*mu * ansatz000;
        Matrix11Row[j] -= val * test100;
        Matrix22Row[j] -= val * test010;
        Matrix33Row[j] -= val * test001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
     Matrix33Row = Matrix_G33[i];
     test000 = Orig5[i];     
     val =  Mult * test000;

     for(j=0;j<N_U;j++)
     {        
        ansatz100 = Orig0[j];
        ansatz010 = Orig1[j];
        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz100;
        Matrix22Row[j] -= val * ansatz010;
        Matrix33Row[j] -= val * ansatz001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     test000 = Orig5[i];     
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz000 = Orig5[j];
        MatrixRow1[j] += Mult * test000 * ansatz000;
     }
  }
}
// ======================================================================
// Type 4, VMS_Projection with streamline formulation D(u):D(v)
// ======================================================================

void TimeNSType4VMS_ProjectionStreamlineDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22,   **Matrix_G33;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, test_stream, ansatz_stream, conv;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_L;
  double c0, c1, c2, c3, delta;
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixL   = LocMatrices[12];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  Matrix_tilde_G11  = LocMatrices[19];
  Matrix_tilde_G22  = LocMatrices[20];
  Matrix_tilde_G33  = LocMatrices[21];
  Matrix_G11  = LocMatrices[22];
  Matrix_G22  = LocMatrices[23];
  Matrix_G33  = LocMatrices[24];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p
  Orig5 = OrigValues[5]; // l
  
  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           param[21]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;
    test_stream = (u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      ansatz_stream = (u1*ansatz100+u2*ansatz010+u3*ansatz001);
      conv = ansatz_stream * test000;
      ansatz_stream *= mu*test_stream;     // sd projection term 

      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += conv + ansatz_stream;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += conv + ansatz_stream;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += conv + ansatz_stream;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];
      val1 = Mult*ansatz000;
      
      val = -val1*test100;
      MatrixRow1[j] += val;
      val = -val1*test010;
      MatrixRow2[j] += val;
      val = -val1*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];
     val1 =u1*test100+u2*test010+u3*test001;

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig5[j];

        val =  Mult * mu * ansatz000 *val1;
        Matrix11Row[j] -= val * u1;
        Matrix22Row[j] -= val * u2;
        Matrix33Row[j] -= val * u3;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
     Matrix33Row = Matrix_G33[i];
     test000 = Orig5[i];     
     val =  Mult * test000;

     for(j=0;j<N_U;j++)
     {        
        ansatz100 = Orig0[j];
        ansatz010 = Orig1[j];
        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz100;
        Matrix22Row[j] -= val * ansatz010;
        Matrix33Row[j] -= val * ansatz001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     test000 = Orig5[i];     
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz000 = Orig5[j];
        MatrixRow1[j] += Mult * test000 * ansatz000;
     }
  }
}

// ======================================================================
// Type 4, LerayAlpha, D(u):D(v)
// ======================================================================
void TimeNSType4LerayAlphaDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *AuxMatrixRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, alpha;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3
 
  u1 = param[0]; // u1 filtered
  u2 = param[1]; // u2 filtered
  u3 = param[2]; // u3 filtered
  
  alpha = LerayAlpha_Param3D(hK);
 
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
 
      val  = alpha*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;            
   } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// assemble matrix for auxiliary problem
// ======================================================================

void MatrixAuxiliaryProblem(double Mult, double *coeff, 
                            double *param, double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double *AuxMatrixRow, **AuxMatrix;;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U; 
  double delta, mu2, val;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;
 
  // solution does not need to be convolved
  AuxMatrix = LocMatrices[0];
  N_U = N_BaseFuncts[0];
     
  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
    
  // filter width
  delta =  CharacteristicFilterWidth(hK);
  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;
  
  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
       
      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
      
    } // endfor j
  } // endfor i
}

// ======================================================================
// Assembling routine for all nonlinear matrices
// ======================================================================

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// Type 2, Standard Galerkin, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu
 
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// Type 2, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1_2NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  c0 = coeff[0]; // nu

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// Type 1, ClassicalLES, only nonlinear part
// Type 2, ClassicalLES, only nonlinear part
// Type 1, GL00Convolution, only nonlinear part
// Type 2, GL00Convolution, only nonlinear part
// Type 1, GL00AuxProblem, only nonlinear part
// Type 2, GL00AuxProblem, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0, mu, delta;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu
 
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[12],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val, val1;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = Mult*coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = c0*Orig0[i];
    test010 = c0*Orig1[i];
    test001 = c0*Orig2[i];
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 += (test100*ansatz100+test010*ansatz010 +test001*ansatz001);
      Matrix11Row[j] += test100*ansatz100+val1;
      Matrix22Row[j] += test010*ansatz010+val1;
      Matrix33Row[j] += test001*ansatz001+val1;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// Type 4, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3_4NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = Mult*c0*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// Type 4, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3_4NLUpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      Matrix11Row[j] += Mult * val;
      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      Matrix22Row[j] += Mult * val;
      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 3, ClassicalLES, (grad u, grad v), only nonlinear part
// Type 4, ClassicalLES, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0, mu, delta;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[12],&param[12],&param[13],&param[14],
                           -4711);
 
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 3, ClassicalLES, D(u):D(v), only nonlinear diagonal blocks
// Type 4, ClassicalLES, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, viscosity, delta, dummy_param[6] = {0, 0, 0, 0, 0, 0};
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&dummy_param[0],&param[12],&param[13],&param[14],
                           -4711);
  viscosity = Mult*(mu/2.0+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = viscosity*Orig0[i];
    test010 = viscosity*Orig1[i];
    test001 = viscosity*Orig2[i];
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      //t100a100 = test100*ansatz100;
      //t010a010 = test010*ansatz010;
      //t001a001 = test001*ansatz001;

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      /* val1 += (test100*ansatz100+test010*ansatz010+test001*ansatz001);
      Matrix11Row[j] += test100*ansatz100+val1;
      Matrix12Row[j] += test100*ansatz010; 
      Matrix13Row[j] += test100*ansatz001;
      Matrix21Row[j] += test010*ansatz100;
      Matrix22Row[j] += test010*ansatz010+val1;
      Matrix23Row[j] += test010*ansatz001;
      Matrix31Row[j] += test001*ansatz100;
      Matrix32Row[j] += test001*ansatz010;
      Matrix33Row[j] += test001*ansatz001+val1;*/
      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1;
      Matrix12Row[j] += test010*ansatz100; 
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLVMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4, valu1, valu2, valu3;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33; 
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P, N_L;
  double c0, viscosity, delta;
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_tilde_G33  = LocMatrices[11];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old


	  if(TDatabase::ParamDB->TURBULENT_MOD_TYPE ==0)
	  { 
	    delta =  CharacteristicFilterWidth(hK);
	    mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14], param[21]);
	  }
	  else if(TDatabase::ParamDB->TURBULENT_MOD_TYPE ==1)
	  {
	    
	    mu=(TDatabase::ParamDB->viscosity_max + TDatabase::ParamDB->viscosity_min)/2.0;
	    
	  }
	  else if (TDatabase::ParamDB->TURBULENT_MOD_TYPE ==2)
	  {
	    mu=TDatabase::ParamDB->viscosity_max;
	  }
	  else if(TDatabase::ParamDB->TURBULENT_MOD_TYPE ==3)
	  {
	    mu=TDatabase::ParamDB->viscosity_min;
	  }
  
    

	  viscosity = Mult*(mu/2.0+c0);	
      


  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = viscosity*Orig0[i];
    test010 = viscosity*Orig1[i];
    test001 = viscosity*Orig2[i];
    test000 = Mult*Orig3[i];
    valu1 = u1 * test000;
    valu2 = u2 * test000;
    valu3 = u3 * test000;

    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      //val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 = valu1*ansatz100+valu2*ansatz010+valu3*ansatz001;

      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1;
      Matrix12Row[j] += test010*ansatz100; 
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1;
    } // endfor j
  } // endfor i

  val2 = Mult * mu;
  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {       
	 //ansatz000 = Orig4[j];
        val1 = val2 * Orig4[j];
        Matrix11Row[j] -= val1 * test100;
        Matrix22Row[j] -= val1 * test010;
        Matrix33Row[j] -= val1 * test001;
     }
  }   
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), adaptive coarse space
// Type 4, VMS_Projection, D(u):D(v), adaptive coarse space
// ======================================================================
void TimeNSType3_4NL_Adap_VMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22,   **Matrix_G33;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_L;
  double c0, c1, c2, c3, delta, val, val1;
  double u1, u2, u3, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixL   = LocMatrices[9];
  Matrix_tilde_G11  = LocMatrices[10];
  Matrix_tilde_G22  = LocMatrices[11];
  Matrix_tilde_G33  = LocMatrices[12];
  Matrix_G11  = LocMatrices[13];
  Matrix_G22  = LocMatrices[14];
  Matrix_G33  = LocMatrices[15];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l
  
  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           param[21]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
 
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig4[j];
        val =  Mult * 2*mu * ansatz000;
        Matrix11Row[j] -= val * test100;
        Matrix22Row[j] -= val * test010;
        Matrix33Row[j] -= val * test001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
     Matrix33Row = Matrix_G33[i];
     test000 = Orig4[i];     
     val =  Mult * test000;

     for(j=0;j<N_U;j++)
     {        
        ansatz100 = Orig0[j];
        ansatz010 = Orig1[j];
        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz100;
        Matrix22Row[j] -= val * ansatz010;
        Matrix33Row[j] -= val * ansatz001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     test000 = Orig4[i];     
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz000 = Orig4[j];
        MatrixRow1[j] += Mult * test000 * ansatz000;
     }
  }
}

// ======================================================================
// Type 3, VMS_Projection explicit, only Matrix_tilde_G??
// Type 4, VMS_Projection explicit, only Matrix_tilde_G??
// ======================================================================
void TimeNSType3_4VMS_ProjectionExpl3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double val1;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33; 
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz000;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig4;
  int i,j,N_U, N_L;
  double delta, mu;

  Matrix_tilde_G11  = LocMatrices[0];
  Matrix_tilde_G22  = LocMatrices[1];
  Matrix_tilde_G33  = LocMatrices[2];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig4 = OrigValues[3]; // l
 //OutPut(param[12] << " "<< param[13] << " " <<param[14] << " ");
 //OutPut(param[0] << " "<< param[1] << " " <<param[2] << " "<<endl);
  OutPut("CHECK IF PARAM[21] IS SET !!!"<<endl);
  exit(4711);
  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[15],&param[12],&param[13],&param[14],
                           param[21]);
  
  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig4[j];
        val1 = Mult * mu * ansatz000;
        Matrix11Row[j] -= val1 * test100;
        Matrix22Row[j] -= val1 * test010;
        Matrix33Row[j] -= val1 * test001;
     }
  }   
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// streamline projection
// ======================================================================
void TimeNSType3_4NLVMS_ProjectionStreamlineDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33; 
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P, N_L;
  double c0, delta, test_stream, ansatz_stream;
  double u1, u2, u3, mu, val;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_tilde_G33  = LocMatrices[11];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[12],&param[12],&param[13],&param[14],
                           param[21]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test_stream = Mult*(u1*test100+u2*test010+u3*test001);
    test100 *= Mult*c0;
    test010 *= Mult*c0;
    test001 *= Mult*c0;
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz_stream = (u1*ansatz100+u2*ansatz010+u3*ansatz001); // sd direction
      val1 = ansatz_stream*test000;           // convective term
      ansatz_stream *= mu * test_stream;      // sd projection term 

      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1+ansatz_stream;
      Matrix12Row[j] += test010*ansatz100; 
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1+ansatz_stream;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1+ansatz_stream;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];
     val1 =u1*test100+u2*test010+u3*test001;

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig4[j];
        val = Mult * mu * ansatz000 * val1;
        Matrix11Row[j] -= val * u1;
        Matrix22Row[j] -= val * u2;
        Matrix33Row[j] -= val * u3;
     }
  }   
}

// ======================================================================
// Type 3, div-div stabilization
// Type 4, div-div stabilization
// ======================================================================
void TimeNSType3_4NLDivDivDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double tautest100, tautest010, tautest001;
  double c0test100, c0test010, c0test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, tau;
  double u1, u2, u3, mu;
  double val;
  double theta1 = TDatabase::TimeDB->THETA1;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  tau = DivDivStab3D(u1,u2,u3,hK,c0);
  //OutPut(" "<<tau);
  // for time-dependent problems, the div-div term is incorporated
  // into the matrix A which will be multiplied with theta1
  if (theta1 > 0)
      tau /= theta1;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tautest100 = tau*test100;
    tautest010 = tau*test010;
    tautest001 = tau*test001;
    c0test100 = c0*test100;
    c0test010 = c0*test010;
    c0test001 = c0*test001;


    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      // convection term
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001)+tautest100*ansatz100;
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0test010*ansatz100+tautest100*ansatz010;
      Matrix12Row[j] += Mult * val;

      val  = c0test001*ansatz100+tautest100*ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0test100*ansatz010+tautest010*ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001)+tautest010*ansatz010;
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0test001*ansatz010+tautest010*ansatz001;
      Matrix23Row[j] += Mult * val;

      val  = c0test100*ansatz001+tautest001*ansatz100;
      Matrix31Row[j] += Mult * val;

      val  = c0test010*ansatz001+tautest001*ansatz010;
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001)+tautest001*ansatz001;
      val += val1;
      Matrix33Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 14, Extra terms in Hughes D(u):D(v)
//         div-div, SUPG
// ======================================================================
void TimeNSType14VMS_SUPGDD3D(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **MatrixK11, **MatrixK12, **MatrixK13;
  double **MatrixK21, **MatrixK22, **MatrixK23, **MatrixK31;
  double **MatrixK32, **MatrixK33;
  double **MatrixS11, **MatrixS12, **MatrixS13, **MatrixS21;
  double **MatrixS22, **MatrixS23, **MatrixS31, **MatrixS32, **MatrixS33;
  double **MatrixM11, **MatrixM22, **MatrixM33, **MatrixC;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixBT1, **MatrixBT2,  **MatrixBT3;
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6, *Rhs7, val, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixK11Row, *MatrixK12Row, *MatrixK13Row, *MatrixK21Row;
  double *MatrixK22Row, *MatrixK23Row, *MatrixK31Row, *MatrixK32Row;
  double *MatrixK33Row;
  double *MatrixS11Row, *MatrixS12Row, *MatrixS13Row, *MatrixS21Row;
  double *MatrixS22Row, *MatrixS23Row, *MatrixS31Row, *MatrixS32Row;
  double *MatrixS33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row, *MatrixCRow;
  double *MatrixB1Row, *MatrixB2Row, *MatrixB3Row;
  double *MatrixBT1Row, *MatrixBT2Row, *MatrixBT3Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double ansatz200, ansatz020, ansatz002;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  double sh1, sh2, ah, bth1, h1, norm_u, temp1, sh3, m1, m2, ah2, bh1;
  //OutPut("supg");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7, *Orig8, *Orig9, *Orig10;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double supg_params[2];

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  // matrices for vicous and convective term
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8]; 
  // mass matrix
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  // matrices for SUPG
  MatrixK11 = LocMatrices[12];
  MatrixK12 = LocMatrices[13];
  MatrixK13 = LocMatrices[14];
  MatrixK21 = LocMatrices[15];
  MatrixK22 = LocMatrices[16];
  MatrixK23 = LocMatrices[17];
  MatrixK31 = LocMatrices[18];
  MatrixK32 = LocMatrices[19];
  MatrixK33 = LocMatrices[20];
  // matrices for div-div term + extra terms
  MatrixS11 = LocMatrices[21];
  MatrixS12 = LocMatrices[22];
  MatrixS13 = LocMatrices[23];
  MatrixS21 = LocMatrices[24];
  MatrixS22 = LocMatrices[25];
  MatrixS23 = LocMatrices[26];
  MatrixS31 = LocMatrices[27];
  MatrixS32 = LocMatrices[28];
  MatrixS33 = LocMatrices[29]; 

  MatrixC   = LocMatrices[30];

  // matrices for divergence constraint + extra
  MatrixB1  = LocMatrices[31];
  MatrixB2  = LocMatrices[32];
  MatrixB3  = LocMatrices[33];
  // matrices for pressure term in momentum equations +extra
  MatrixBT1 = LocMatrices[34];
  MatrixBT2 = LocMatrices[35];
  MatrixBT3 = LocMatrices[36];

  // right hand sides
  // for velocity space test functions 
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];
  // for pressure space test functions
  Rhs7 = LocRhs[6];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z
  Orig7 = OrigValues[7];         // p
//   Orig8 = OrigValues[8];         // u_xx
//   Orig9 = OrigValues[9];         // u_yy
//   Orig10 = OrigValues[10];       // u_zz

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z



  
  // second order derivatives in the residual will be neglected
  // method is for flows with small viscosity



  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);


///////////////////////////////////////////////////////////////////////////////////////////////////////  
//calculation of modified stab. parameters  
///////////////////////////////////////////////////////////////////////////////////////////////////////
//    norm_u = sqrt(u1*u1+u2*u2+u3*u3);
//    temp1  = (0.5*hK*norm_u)*(0.5*hK*norm_u);
//    if(TDatabase::TimeDB->CURRENTTIME > 1.0)
//      supg_params[0] = TDatabase::ParamDB->DELTA0*hK*hK;
//    else
//      supg_params[0] = 0;
//    //graddiv par. (John/Kindl,2010)
//    supg_params[1] = 0.5*sqrt(c0*c0 + temp1);
supg_params[0] = hK;
supg_params[1] = 0.0;
////////////////////////////////////////////////////////////////////////////////////////////////////////
//end of calculation of modified stab. parameters  
/////////////////////////////////////////////////////////////////////////////////////////////////////////



   
   
   // assembling for velocity test functions
  
  // dummy parameter
  

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row = MatrixM11[i];
    MatrixM22Row = MatrixM22[i];
    MatrixM33Row = MatrixM33[i];
    MatrixK11Row = MatrixK11[i];
    MatrixK12Row = MatrixK12[i];
    MatrixK13Row = MatrixK13[i];
    MatrixK21Row = MatrixK21[i];
    MatrixK22Row = MatrixK22[i];
    MatrixK23Row = MatrixK23[i];
    MatrixK31Row = MatrixK31[i];
    MatrixK32Row = MatrixK32[i];
    MatrixK33Row = MatrixK33[i];
    MatrixS11Row = MatrixS11[i];
    MatrixS12Row = MatrixS12[i];
    MatrixS13Row = MatrixS13[i];
    MatrixS21Row = MatrixS21[i];
    MatrixS22Row = MatrixS22[i];
    MatrixS23Row = MatrixS23[i];
    MatrixS31Row = MatrixS31[i];
    MatrixS32Row = MatrixS32[i];
    MatrixS33Row = MatrixS33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    //for Rhsi, i=4,5,6
    h1 = c1*test100 + c2*test010 + c3*test001;

    //for matrices Sij, Kij, BTi, Rhs
    sh2 = u1*test100+u2*test010+u3*test001;
    
    m1 = Mult*test000;
    m2 = Mult*supg_params[0];
      
    // rhs, this is the part of the term which will be multiplied by theta4*tau
    Rhs1[i] += m1*c1;
    Rhs2[i] += m1*c2;
    Rhs3[i] += m1*c3;
    Rhs4[i] += m2*(u1*h1 + c1*sh2);
    Rhs5[i] += m2*(u2*h1 + c2*sh2); 
    Rhs6[i] += m2*(u3*h1 + c3*sh2);
   

    // test functions for div-div term
    tautest100 = supg_params[1]*test100;
    tautest010 = supg_params[1]*test010;
    tautest001 = supg_params[1]*test001;
    
    // velocity-velocity block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
//       ansatz200 = Orig8[j];
//       ansatz020 = Orig9[j];
//       ansatz002 = Orig10[j];

      //for matrices Sij, Aij, Bi
      sh1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      
//       //for laplace term in S
//       sh3 = ansatz200 + ansatz020 + ansatz002; 
      
      
      
      // matrices Aij
      // this block will be multiplied with theta1*Delta t
      // convection 
      ah = sh1 * test000;
      ah2 = test100*ansatz100+test010*ansatz010+test001*ansatz001;
      // diffusion
      val  = c0*(test100*ansatz100+ah2);
      // add everything
      val += ah;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test010*ansatz010+ah2);
      val += ah;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(ah2+test001*ansatz001);
      val += ah;
      Matrix33Row[j] += Mult * val;

      // mass matrix M
      val = ansatz000*test000;
      MatrixM11Row[j] += Mult * val;
      MatrixM22Row[j] += Mult * val;
      MatrixM33Row[j] += Mult * val;
      
      // div-div term + extras
      // store in matrices S
      MatrixS11Row[j] += Mult * (tautest100*ansatz100+supg_params[0]*sh1*(u1*test100+sh2));
      MatrixS12Row[j] += Mult * (tautest100*ansatz010+supg_params[0]*sh1*u1*test010);
      MatrixS13Row[j] += Mult * (tautest100*ansatz001+supg_params[0]*sh1*u1*test001);
      MatrixS21Row[j] += Mult * (tautest010*ansatz100+supg_params[0]*sh1*u2*test100);
      MatrixS22Row[j] += Mult * (tautest010*ansatz010+supg_params[0]*sh1*(u2*test010+sh2));
      MatrixS23Row[j] += Mult * (tautest010*ansatz001+supg_params[0]*sh1*u2*test001);
      MatrixS31Row[j] += Mult * (tautest001*ansatz100+supg_params[0]*sh1*u3*test100);
      MatrixS32Row[j] += Mult * (tautest001*ansatz010+supg_params[0]*sh1*u3*test010);
      MatrixS33Row[j] += Mult * (tautest001*ansatz001+supg_params[0]*sh1*(u3*test001+sh2));
      
//       // laplace term in S
//       if(TDatabase::ParamDB->DELTA1 == 100){
//  OutPut("LAPLACEEEEEEEEEEEEEEEEEEEE");
//  MatrixS11Row[j] += Mult * (-supg_params[0]*c0*sh3*(sh2+u1*test100));
//  MatrixS12Row[j] += Mult * (-supg_params[0]*c0*sh3*u1*test010);
//  MatrixS13Row[j] += Mult * (-supg_params[0]*c0*sh3*u1*test001);
//  MatrixS21Row[j] += Mult * (-supg_params[0]*c0*sh3*u2*test100);
//  MatrixS22Row[j] += Mult * (-supg_params[0]*c0*sh3*(sh2+u2*test010));
//  MatrixS23Row[j] += Mult * (-supg_params[0]*c0*sh3*u2*test001);
//  MatrixS31Row[j] += Mult * (-supg_params[0]*c0*sh3*u3*test100);
//  MatrixS32Row[j] += Mult * (-supg_params[0]*c0*sh3*u3*test010);
//  MatrixS33Row[j] += Mult * (-supg_params[0]*c0*sh3*(sh2+u3*test001));
//       }
      
      // SUPG terms 
      // store in matrices K
      MatrixK11Row[j] += m2 * ansatz000 * (u1*test100+sh2);
      MatrixK12Row[j] += m2 * u1 * ansatz000 * test010;
      MatrixK13Row[j] += m2 * u1 * ansatz000 * test001;
      MatrixK21Row[j] += m2 * u2 * ansatz000 * test100;
      MatrixK22Row[j] += m2 * ansatz000 * (u2*test010+sh2);
      MatrixK23Row[j] += m2 * u2 * ansatz000 * test001;
      MatrixK31Row[j] += m2 * u3 * ansatz000 * test100;
      MatrixK32Row[j] += m2 * u3 * ansatz000 * test010;
      MatrixK33Row[j] += m2 * ansatz000 * (u3*test001+sh2);

    }                            // endfor j

    // pressure-velocity block, these blocks will be multiplied with Delta t
    MatrixBT1Row = MatrixBT1[i];
    MatrixBT2Row = MatrixBT2[i];
    MatrixBT3Row = MatrixBT3[i];
    for(j=0;j<N_P;j++)
    {
      // pressure ansatz functions
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      // for matrices BTi
      bth1 = ansatz100*test100 + ansatz010*test010 + ansatz001*test001;
     
      // matrices BTi
      MatrixBT1Row[j] += Mult*(-ansatz000*test100+supg_params[0]*(u1*bth1+ansatz100*sh2));
      MatrixBT2Row[j] += Mult*(-ansatz000*test010+supg_params[0]*(u2*bth1+ansatz010*sh2));
      MatrixBT3Row[j] += Mult*(-ansatz000*test001+supg_params[0]*(u3*bth1+ansatz001*sh2));
      
    }        // endfor j
  }                              // endfor i

  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    MatrixB1Row = MatrixB1[i];
    MatrixB2Row = MatrixB2[i];
    MatrixB3Row = MatrixB3[i];
    MatrixCRow = MatrixC[i];

    test100 = Orig4[i];
    test010 = Orig5[i];
    test001 = Orig6[i];
    test000 = Orig7[i];

    Rhs7[i] += Mult*time_step*supg_params[0]*(((1.0/time_step)*u1 + c1)*test100 
            + ((1.0/time_step)*u2 + c2)*test010 + ((1.0/time_step)*u3 + c3)*test001);

  
    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
//       ansatz200 = Orig8[j];
//       ansatz020 = Orig9[j];
//       ansatz002 = Orig10[j];
      
      sh1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      
      //for laplace term in S
      sh3 = ansatz200 + ansatz020 + ansatz002;
      bh1 = time_step*supg_params[0]*((1./time_step)*ansatz000+sh1);
      
      // matrices Bi
      MatrixB1Row[j] += Mult*(ansatz100*test000+bh1*test100);
      MatrixB2Row[j] += Mult*(ansatz010*test000+bh1*test010);
      MatrixB3Row[j] += Mult*(ansatz001*test000+bh1*test001);
      
      
//       //laplace term in B 
//       if(TDatabase::ParamDB->DELTA1 == 100){
//  MatrixB1Row[j] += Mult*(-supg_params[0]*c0*sh3*test100);
//  MatrixB2Row[j] += Mult*(-supg_params[0]*c0*sh3*test010);
//  MatrixB3Row[j] += Mult*(-supg_params[0]*c0*sh3*test001);
//       }
    
    }
    
    
    

    // pressure-pressure block
    for(j=0;j<N_P;j++)
    {
   
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];

      //matrix C
      MatrixCRow[j] +=  m2 * time_step * (ansatz100*test100 + ansatz010*test010 + ansatz001*test001);
    }                   // endfor j
  }     // endfor i
}


// ======================================================================
// Type 4, Extra terms in Hughes D(u):D(v)
//         div-div, SUPG
// ======================================================================
void TimeNSType14NLVMS_SUPGDD3D(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **MatrixK11, **MatrixK12, **MatrixK13;
  double **MatrixK21, **MatrixK22, **MatrixK23, **MatrixK31;
  double **MatrixK32, **MatrixK33;
  double **MatrixS11, **MatrixS12, **MatrixS13, **MatrixS21;
  double **MatrixS22, **MatrixS23, **MatrixS31, **MatrixS32, **MatrixS33;
  double **MatrixC;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixBT1, **MatrixBT2,  **MatrixBT3;
  double *Rhs4, *Rhs5, *Rhs6, val, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixK11Row, *MatrixK12Row, *MatrixK13Row, *MatrixK21Row;
  double *MatrixK22Row, *MatrixK23Row, *MatrixK31Row, *MatrixK32Row;
  double *MatrixK33Row;
  double *MatrixS11Row, *MatrixS12Row, *MatrixS13Row, *MatrixS21Row;
  double *MatrixS22Row, *MatrixS23Row, *MatrixS31Row, *MatrixS32Row;
  double *MatrixS33Row;
  double *MatrixCRow;
  double *MatrixB1Row, *MatrixB2Row, *MatrixB3Row;
  double *MatrixBT1Row, *MatrixBT2Row, *MatrixBT3Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double ansatz200, ansatz020, ansatz002;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  double sh1, sh2, ah, bth1, h1, norm_u, sh3, m2, ah2, bh1;
  //OutPut("supg");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7, *Orig8, *Orig9, *Orig10;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double supg_params[2]   ;

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  // matrices for vicous and convective term
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8]; 
  // matrix for SUPG
  MatrixK11 = LocMatrices[9];
  MatrixK12 = LocMatrices[10];
  MatrixK13 = LocMatrices[11];
  MatrixK21 = LocMatrices[12];
  MatrixK22 = LocMatrices[13];
  MatrixK23 = LocMatrices[14];
  MatrixK31 = LocMatrices[15];
  MatrixK32 = LocMatrices[16];
  MatrixK33 = LocMatrices[17];
  // matrices for div-div term + extra terms
  MatrixS11 = LocMatrices[18];
  MatrixS12 = LocMatrices[19];
  MatrixS13 = LocMatrices[20];
  MatrixS21 = LocMatrices[21];
  MatrixS22 = LocMatrices[22];
  MatrixS23 = LocMatrices[23];
  MatrixS31 = LocMatrices[24];
  MatrixS32 = LocMatrices[25];
  MatrixS33 = LocMatrices[26]; 

  MatrixC   = LocMatrices[27];

  // matrices for divergence constraint + extra
  MatrixB1  = LocMatrices[28];
  MatrixB2  = LocMatrices[29];
  MatrixB3  = LocMatrices[30];
  // matrices for pressure term in momentum equations +extra
  MatrixBT1 = LocMatrices[31];
  MatrixBT2 = LocMatrices[32];
  MatrixBT3 = LocMatrices[33];

  // right hand sides
  Rhs4 = LocRhs[0];
  Rhs5 = LocRhs[1];
  Rhs6 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z
  Orig7 = OrigValues[7];         // p
//   Orig8 = OrigValues[8];         // u_xx
//   Orig9 = OrigValues[9];         // u_yy
//   Orig10 = OrigValues[10];       // u_zz

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z
  
  // second order derivatives in the residual will be neglected
  // method is for flows with small viscosity

  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);

  
///////////////////////////////////////////////////////////////////////////////////////////////////////  
//calculation of modified stab. parameters  
///////////////////////////////////////////////////////////////////////////////////////////////////////
//    norm_u = sqrt(u1*u1+u2*u2+u3*u3);
//    temp1  = (0.5*hK*norm_u)*(0.5*hK*norm_u);
//    if(TDatabase::TimeDB->CURRENTTIME > 1.0)
//      supg_params[0] = TDatabase::ParamDB->DELTA0*hK*hK;
//    else
//      supg_params[0] = 0;
//    //graddiv par. (John/Kindl,2010)  
//    supg_params[1] = 0.5*sqrt(c0*c0 + temp1);
supg_params[0] = hK;
supg_params[1] = 0.0;
////////////////////////////////////////////////////////////////////////////////////////////////////////
//end of calculation of modified stab. parameters  
/////////////////////////////////////////////////////////////////////////////////////////////////////////



  

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixK11Row = MatrixK11[i];
    MatrixK12Row = MatrixK12[i];
    MatrixK13Row = MatrixK13[i];
    MatrixK21Row = MatrixK21[i];
    MatrixK22Row = MatrixK22[i];
    MatrixK23Row = MatrixK23[i];
    MatrixK31Row = MatrixK31[i];
    MatrixK32Row = MatrixK32[i];
    MatrixK33Row = MatrixK33[i];
    MatrixS11Row = MatrixS11[i];
    MatrixS12Row = MatrixS12[i];
    MatrixS13Row = MatrixS13[i];
    MatrixS21Row = MatrixS21[i];
    MatrixS22Row = MatrixS22[i];
    MatrixS23Row = MatrixS23[i];
    MatrixS31Row = MatrixS31[i];
    MatrixS32Row = MatrixS32[i];
    MatrixS33Row = MatrixS33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
//     // SUPG term for the time derivative, with scaling, theta1 already scaled with Delta t
//     ugradu  = (u1*test100+u2*test010+u3*test001)*supg_params[0]*theta1;

    //for Rhsi, i=4,5,6
    h1 = c1*test100 + c2*test010 + c3*test001;

    //for matrices Sij, Kij, BTi, Rhs, Sij
    sh2 = u1*test100+u2*test010+u3*test001;
    
    m2 = Mult*supg_params[0];   
    // rhs, this is part of the term which will be multiplied with theta4
    
    Rhs4[i] += m2*(u1*h1 + c1*sh2);
    Rhs5[i] += m2*(u2*h1 + c2*sh2); 
    Rhs6[i] += m2*(u3*h1 + c3*sh2);
    
    // test functions for div-div term
    tautest100 = supg_params[1]*test100;
    tautest010 = supg_params[1]*test010;
    tautest001 = supg_params[1]*test001;

    // velocity-velocity block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
//       ansatz200 = Orig8[j];
//       ansatz020 = Orig9[j];
//       ansatz002 = Orig10[j];

      //for matrices Sij, Aij, Bi
      sh1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      
//       //for laplace term in S
//       sh3 = ansatz200 + ansatz020 + ansatz002;

      // matrices Aij
      // this block will be multiplied with theta1*Delta t
      // convection 
      ah = sh1 * test000;
      ah2 = test100*ansatz100+test010*ansatz010+test001*ansatz001;
      // diffusion
      val  = c0*(test100*ansatz100+ah2);
      // add everything
      val += ah;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(ah2 + test010*ansatz010);
      val += ah;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(ah2 + test001*ansatz001);
      val += ah;
      Matrix33Row[j] += Mult * val;
      
      // div-div term + extras
      // store in matrices S
      MatrixS11Row[j] += Mult * (tautest100*ansatz100+supg_params[0]*sh1*(u1*test100+sh2));
      MatrixS12Row[j] += Mult * (tautest100*ansatz010+supg_params[0]*sh1*u1*test010);
      MatrixS13Row[j] += Mult * (tautest100*ansatz001+supg_params[0]*sh1*u1*test001);
      MatrixS21Row[j] += Mult * (tautest010*ansatz100+supg_params[0]*sh1*u2*test100);
      MatrixS22Row[j] += Mult * (tautest010*ansatz010+supg_params[0]*sh1*(u2*test010+sh2));
      MatrixS23Row[j] += Mult * (tautest010*ansatz001+supg_params[0]*sh1*u2*test001);
      MatrixS31Row[j] += Mult * (tautest001*ansatz100+supg_params[0]*sh1*u3*test100);
      MatrixS32Row[j] += Mult * (tautest001*ansatz010+supg_params[0]*sh1*u3*test010);
      MatrixS33Row[j] += Mult * (tautest001*ansatz001+supg_params[0]*sh1*(u3*test001+sh2));
      
//       // laplace term in S
//       if(TDatabase::ParamDB->DELTA1 == 100){
//  MatrixS11Row[j] += Mult * (-supg_params[0]*c0*sh3*(sh2+u1*test100));
//  MatrixS12Row[j] += Mult * (-supg_params[0]*c0*sh3*u1*test010);
//  MatrixS13Row[j] += Mult * (-supg_params[0]*c0*sh3*u1*test001);
//  MatrixS21Row[j] += Mult * (-supg_params[0]*c0*sh3*u2*test100);
//  MatrixS22Row[j] += Mult * (-supg_params[0]*c0*sh3*(sh2+u2*test010));
//  MatrixS23Row[j] += Mult * (-supg_params[0]*c0*sh3*u2*test001);
//  MatrixS31Row[j] += Mult * (-supg_params[0]*c0*sh3*u3*test100);
//  MatrixS32Row[j] += Mult * (-supg_params[0]*c0*sh3*u3*test010);
//  MatrixS33Row[j] += Mult * (-supg_params[0]*c0*sh3*(sh2+u3*test001));
//       }
      
      // SUPG terms 
      // store in matrices K
      MatrixK11Row[j] += m2 * ansatz000 * (u1*test100+sh2);
      MatrixK12Row[j] += m2 * u1 * ansatz000 * test010;
      MatrixK13Row[j] += m2 * u1 * ansatz000 * test001;
      MatrixK21Row[j] += m2 * u2 * ansatz000 * test100;
      MatrixK22Row[j] += m2 * ansatz000 * (u2*test010+sh2);
      MatrixK23Row[j] += m2 * u2 * ansatz000 * test001;
      MatrixK31Row[j] += m2 * u3 * ansatz000 * test100;
      MatrixK32Row[j] += m2 * u3 * ansatz000 * test010;
      MatrixK33Row[j] += m2 * ansatz000 * (u3*test001+sh2);

    }                            // endfor j

    // pressure-velocity block, these blocks will be multiplied with Delta t
    MatrixBT1Row = MatrixBT1[i];
    MatrixBT2Row = MatrixBT2[i];
    MatrixBT3Row = MatrixBT3[i];
    for(j=0;j<N_P;j++)
    {
      // pressure ansatz functions
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      // for matrices BTi
      bth1 = ansatz100*test100 + ansatz010*test010 + ansatz001*test001;
     
      // matrices BTi
      MatrixBT1Row[j] += Mult*(-ansatz000*test100+supg_params[0]*(u1*bth1+ansatz100*sh2));
      MatrixBT2Row[j] += Mult*(-ansatz000*test010+supg_params[0]*(u2*bth1+ansatz010*sh2));
      MatrixBT3Row[j] += Mult*(-ansatz000*test001+supg_params[0]*(u3*bth1+ansatz001*sh2));
      
    }        // endfor j
  }                              // endfor i

  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    MatrixB1Row = MatrixB1[i];
    MatrixB2Row = MatrixB2[i];
    MatrixB3Row = MatrixB3[i];
    MatrixCRow = MatrixC[i];

    test100 = Orig4[i];
    test010 = Orig5[i];
    test001 = Orig6[i];
    test000 = Orig7[i];
  
    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
//       ansatz200 = Orig8[j];
//       ansatz020 = Orig9[j];
//       ansatz002 = Orig10[j];
      
      sh1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      
      //for laplace term in S
      sh3 = ansatz200 + ansatz020 + ansatz002;
      bh1 = time_step*supg_params[0]*((1./time_step)*ansatz000+sh1);
      
      // matrices Bi
      MatrixB1Row[j] += Mult*(ansatz100*test000+bh1*test100);
      MatrixB2Row[j] += Mult*(ansatz010*test000+bh1*test010);
      MatrixB3Row[j] += Mult*(ansatz001*test000+bh1*test001);
  /*    
      //laplace term in B 
      if(TDatabase::ParamDB->DELTA1 == 100){
  MatrixB1Row[j] += Mult*(-supg_params[0]*c0*sh3*test100);
  MatrixB2Row[j] += Mult*(-supg_params[0]*c0*sh3*test010);
  MatrixB3Row[j] += Mult*(-supg_params[0]*c0*sh3*test001);
      }
      */
      
    }      

    // pressure-pressure block
    for(j=0;j<N_P;j++)
    {
   
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];

      //matrix C
      MatrixCRow[j] += Mult * time_step * supg_params[0] * (ansatz100*test100 + ansatz010*test010 + ansatz001*test001);
    }                   // endfor j
  }     // endfor i
}




// ======================================================================
// Type 4, Extra terms in Hughes D(u):D(v)
//         div-div, SUPG
// ======================================================================
void TimeNSType4VMS_SUPGDD3D(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **MatrixK;
  double **MatrixS11, **MatrixS12, **MatrixS13, **MatrixS21;
  double **MatrixS22, **MatrixS23, **MatrixS31, **MatrixS32, **MatrixS33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *MatrixKRow;
  double *MatrixS11Row, *MatrixS12Row, *MatrixS13Row, *MatrixS21Row;
  double *MatrixS22Row, *MatrixS23Row, *MatrixS31Row, *MatrixS32Row;
  double *MatrixS33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  //OutPut("supg");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double supg_params[2], ugradu;

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  theta1 *=time_step;
  theta2 *=time_step;
  theta3 *=time_step;
  theta4 *=time_step;

  // matrices for vicous and convective term
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8]; 
  // mass matrix
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  // matrix for SUPG
  MatrixK = LocMatrices[12];
  // matrices for div-div term + 1st extra term
  MatrixS11 = LocMatrices[13];
  MatrixS12 = LocMatrices[14];
  MatrixS13 = LocMatrices[15];
  MatrixS21 = LocMatrices[16];
  MatrixS22 = LocMatrices[17];
  MatrixS23 = LocMatrices[18];
  MatrixS31 = LocMatrices[19];
  MatrixS32 = LocMatrices[20];
  MatrixS33 = LocMatrices[21]; 

  // matrices for divergence constraint
  MatrixB1  = LocMatrices[22];
  MatrixB2  = LocMatrices[23];
  MatrixB3  = LocMatrices[24];
  // matrices for pressure term in momentum equations
  MatrixB1T = LocMatrices[25];
  MatrixB2T = LocMatrices[26];
  MatrixB3T = LocMatrices[27];

  // right hand sides
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z
  
  // second order derivatives in the residual will be neglected
  // method is for flows with small viscosity

  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);
  //OutPut(supg_params[0] << " " << supg_params[1] << " : " << endl);
  //supg_params[0] = 0;
  // assembling for velocity test functions
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    MatrixKRow  = MatrixK[i];
    MatrixS11Row = MatrixS11[i];
    MatrixS12Row = MatrixS12[i];
    MatrixS13Row = MatrixS13[i];
    MatrixS21Row = MatrixS21[i];
    MatrixS22Row = MatrixS22[i];
    MatrixS23Row = MatrixS23[i];
    MatrixS31Row = MatrixS31[i];
    MatrixS32Row = MatrixS32[i];
    MatrixS33Row = MatrixS33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    // SUPG term for the time derivative, with scaling, theta1 already scaled with Delta t
    ugradu  = (u1*test100+u2*test010+u3*test001)*supg_params[0]*theta1;
  	
    // rhs, this is part of the term which will be multiplied with theta4
    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    // test functions for div-div term
    tautest100 = supg_params[1]*test100;
    tautest010 = supg_params[1]*test010;
    tautest001 = supg_params[1]*test001;
    // velocity-velocity block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      // this block will be multiplied with theta1*Delta t
      // convection 
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // diffusion
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      // add everything
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;

      // mass matrix
      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
      
      // div-div term
      // store in matrices S
      MatrixS11Row[j] += Mult * tautest100*ansatz100;
      MatrixS12Row[j] += Mult * tautest100*ansatz010;
      MatrixS13Row[j] += Mult * tautest100*ansatz001;
      MatrixS21Row[j] += Mult * tautest010*ansatz100;
      MatrixS22Row[j] += Mult * tautest010*ansatz010;
      MatrixS23Row[j] += Mult * tautest010*ansatz001;
      MatrixS31Row[j] += Mult * tautest001*ansatz100;
      MatrixS32Row[j] += Mult * tautest001*ansatz010;
      MatrixS33Row[j] += Mult * tautest001*ansatz001;
      
      // SUPG term for the time derivative
      // store in matrices K
      MatrixKRow[j] += Mult * ansatz000 * ugradu;
      
      // 1st extra term for the time derivative
      // store in matrices S
      // tau_m * (u\nabla v + u) 
      val1 = supg_params[0]*theta1*(u1*ansatz100+u2*ansatz010+u3*ansatz001)+ansatz000;
      val = Mult*val1*u1*test100;
      MatrixS11Row[j] += val;
      val = Mult*val1*u1*test010;
      MatrixS12Row[j] += val;
      val = Mult*val1*u1*test001;
      MatrixS13Row[j] += val;
      val = Mult*val1*u2*test100;
      MatrixS21Row[j] += val;
      val = Mult*val1*u2*test010;
      MatrixS22Row[j] += val;
      val = Mult*val1*u2*test001;
      MatrixS23Row[j] += val;
      val = Mult*val1*u3*test100;
      MatrixS31Row[j] += val;
      val = Mult*val1*u3*test010;
      MatrixS32Row[j] += val;
      val = Mult*val1*u3*test001;
      MatrixS33Row[j] += val;
    }                            // endfor j

    // pressure-velocity block, these blocks will be multiplied with Delta t
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      // pressure ansatz functions
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];
     
      // pressure term 
      val  = -ansatz000 * test100;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz000 * test010;
      MatrixRow2[j] += Mult*val;
	  
      val  = -ansatz000 * test001;
      MatrixRow3[j] += Mult*val;
      
    }
  }                              // endfor i

  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test100 = Orig4[i];
    test010 = Orig5[i];
    test001 = Orig6[i];
    test000 = Orig7[i];
	
    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      // divergence constraint
      val = -test000*ansatz100;
      MatrixRow1[j] += Mult*val;

      val = -test000*ansatz010;
      MatrixRow2[j] += Mult*val;
	  
      val = -test000*ansatz001;
      MatrixRow3[j] += Mult*val;	  	  
    }                            // endfor j
  }                              // endfor i
}

void TimeNSType4VMS_SUPGDD3D_old(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **MatrixK, **MatrixC;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, val, val1, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *MatrixKRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixRowC;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  //OutPut("supg");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7, *Orig8, *Orig9, *Orig10;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, c4, c5, c6;
  double u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double u1_xx, u1_yy, u1_zz;
  double u2_xx, u2_yy, u2_zz;
  double u3_xx, u3_yy, u3_zz;
  double supg_params[2], ugradu;

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  theta1 *=time_step;
  theta2 *=time_step;
  theta3 *=time_step;
  theta4 *=time_step;

  // matrices for vicous and convective term
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8]; 
  // mass matrix
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  // matrix for PSPG
  MatrixC = LocMatrices[12];
  // matrix for SUPG
  MatrixK = LocMatrices[13];
  // matrices for divergence constraint
  MatrixB1  = LocMatrices[14];
  MatrixB2  = LocMatrices[15];
  MatrixB3  = LocMatrices[16];
  // matrices for pressure term in momentum equations
  MatrixB1T = LocMatrices[17];
  MatrixB2T = LocMatrices[18];
  MatrixB3T = LocMatrices[19];

  // right hand sides
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3
  c4 = coeff[4];                 // f1_old
  c5 = coeff[5];                 // f2_old
  c6 = coeff[6];                 // f3_old

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z
  
  // second order derivatives in the residual will be neglected
  // method is for flows with small viscosity

  u1_xx = u1_yy = u1_zz = 0;
  u2_xx = u2_yy = u2_zz = 0;
  u3_xx = u3_yy = u3_zz = 0;

  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);
  //OutPut(supg_params[0] << " " << supg_params[1] << " : " << endl);
  //supg_params[0] = 0;
  // assembling for velocity test functions
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    MatrixKRow  = MatrixK[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    // SUPG term, with scaling, theta1 already with scaled with Delta t
    ugradu  = (u1*test100+u2*test010+u3*test001)*supg_params[0]*theta1;
	
    // rhs, this is part of the term which will be multiplied with theta4
    Rhs1[i] += Mult*(test000+ugradu)*c1;
    Rhs2[i] += Mult*(test000+ugradu)*c2;
    Rhs3[i] += Mult*(test000+ugradu)*c3;

    // test functions for div-div term
    tautest100 = supg_params[1]*test100;
    tautest010 = supg_params[1]*test010;
    tautest001 = supg_params[1]*test001;
    // velocity-velocity block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      // this block will be multiplied with theta1*Delta t
      // convection 
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // diffusion + div-div term
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001)+tautest100*ansatz100;
      // add everything
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100)+tautest100*ansatz010;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100)+tautest100*ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010)+tautest010*ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001)+tautest010*ansatz010;
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010)+tautest010*ansatz001;
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001)+tautest001*ansatz100;
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001)+tautest001*ansatz010;
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001)+tautest001*ansatz001;
      val += val1;
      Matrix33Row[j] += Mult * val;

      // mass matrix
      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;

      // SUPG block, without second order derivative 
      // convection 
      val = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      // test with streamline derivative
      val *= ugradu * Mult;
      // term in the diagonal blocks
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
      // SUPG term for the time derivative
      // store in matrices K
      MatrixKRow[j] += Mult * ansatz000 * ugradu;
    }                            // endfor j

    // pressure-velocity block, these blocks will be multiplied with Delta t
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      // pressure ansatz functions
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];
     
      // pressure term 
      val  = -ansatz000 * test100;
      // SUPG term
      val +=  ansatz100 * ugradu;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz000 * test010;
      // SUPG term
      val +=  ansatz010 * ugradu;
      MatrixRow2[j] += Mult*val;
	  
      val  = -ansatz000 * test001;
      // SUPG term
      val +=  ansatz001 * ugradu;
      MatrixRow3[j] += Mult*val;
    }
  }                              // endfor i
  //supg_params[0] = 0;
  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];
    MatrixRowC = MatrixC[i];

    test100 = Orig4[i];
    test010 = Orig5[i];
    test001 = Orig6[i];
    test000 = Orig7[i];
	
    // rhs
    // thetas already multiplied with time_step
    val = u1-theta2*(u1*u1_x+u2*u1_y+u3*u1_z);
    val += theta3*c4 + theta4*c1;
    Rhs4[i] = val*test100;
    val = u2-theta2*(u1*u2_x+u2*u2_y+u3*u2_z);
    val += theta3*c5 + theta4*c2;
    Rhs4[i] += val*test010;
    val = u3-theta2*(u1*u3_x+u2*u3_y+u3*u3_z);
    val += theta3*c6 + theta4*c3;
    Rhs4[i] += val*test001;
    Rhs4[i] *= Mult*time_step*supg_params[0];

    // pressure-pressure block
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
	
      val = supg_params[0] * time_step * time_step *
	  (ansatz100*test100+ansatz010*test010+ansatz001*test001);
      MatrixRowC[j] += Mult*val;
    }

    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      // divergence constraint
      //val1 =  -c0*(ansatz200+ansatz020+ansatz002);
      val1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      val1 *= theta1;
      val1 += u1;
      val1 *= time_step * supg_params[0];

      val = -test000*ansatz100 + val1*test100;
      MatrixRow1[j] += Mult*val;

      val = -test000*ansatz010 + val1*test010;
      MatrixRow2[j] += Mult*val;
	  
      val = -test000*ansatz001 + val1*test001;
      MatrixRow3[j] += Mult*val;	  	  
    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// Type 4, Extra terms in Hughes D(u):D(v)
//         div-div, SUPG
// ======================================================================
void TimeNSType4VMS_Rhs_SUPGDD3D(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
   double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  //OutPut("supg_rhs");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double supg_params[2];

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  theta1 *=time_step;
  theta2 *=time_step;
  theta3 *=time_step;
  theta4 *=time_step;

  // right hand sides
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z

  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);
  //OutPut(supg_params[0] << " " << supg_params[1] << " : " << endl);
  //supg_params[0] = 0;
  // assembling for velocity test functions
  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
	
    // rhs, this is part of the term which will be multiplied with theta4
    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
  }                              // endfor i
}

void TimeNSType4VMS_Rhs_SUPGDD3D_old(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, val, val1;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  //OutPut("supg_rhs");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, c4, c5, c6;
  double u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double supg_params[2], ugradu;

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  theta1 *=time_step;
  theta2 *=time_step;
  theta3 *=time_step;
  theta4 *=time_step;

  // right hand sides
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3
  c4 = coeff[4];                 // f1_old
  c5 = coeff[5];                 // f2_old
  c6 = coeff[6];                 // f3_old

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z
  
  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);
  //OutPut(supg_params[0] << " " << supg_params[1] << " : " << endl);
  //supg_params[0] = 0;
  // assembling for velocity test functions
  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    // SUPG term, with scaling, theta1 already with scaled with Delta t
    ugradu  = (u1*test100+u2*test010+u3*test001)*supg_params[0]*theta1;
	
    // rhs, this is part of the term which will be multiplied with theta4
    Rhs1[i] += Mult*(test000+ugradu)*c1;
    Rhs2[i] += Mult*(test000+ugradu)*c2;
    Rhs3[i] += Mult*(test000+ugradu)*c3;
  }                              // endfor i
  //supg_params[0] = 0;

  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    test100 = Orig4[i];
    test010 = Orig5[i];
    test001 = Orig6[i];
	
    // rhs
    // thetas already multiplied with time_step
    val = u1-theta2*(u1*u1_x+u2*u1_y+u3*u1_z);
    val += theta3*c4 + theta4*c1;
    Rhs4[i] = val*test100;
    val = u2-theta2*(u1*u2_x+u2*u2_y+u3*u2_z);
    val += theta3*c5 + theta4*c2;
    Rhs4[i] += val*test010;
    val = u3-theta2*(u1*u3_x+u2*u3_y+u3*u3_z);
    val += theta3*c6 + theta4*c3;
    Rhs4[i] += val*test001;
    Rhs4[i] *= Mult*time_step*supg_params[0];
  }                              // endfor i
}

// ======================================================================
// Type 4, Extra terms in Hughes D(u):D(v)
//         div-div, SUPG
// ======================================================================
void TimeNSType4NLVMS_SUPGDD3D(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **MatrixK;
  double **MatrixS11, **MatrixS12, **MatrixS13, **MatrixS21;
  double **MatrixS22, **MatrixS23, **MatrixS31, **MatrixS32, **MatrixS33;
  double val, val1, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *MatrixKRow;
  double *MatrixS11Row, *MatrixS12Row, *MatrixS13Row, *MatrixS21Row;
  double *MatrixS22Row, *MatrixS23Row, *MatrixS31Row, *MatrixS32Row;
  double *MatrixS33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  //OutPut("supg");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double supg_params[2], ugradu;

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  theta1 *=time_step;
  theta2 *=time_step;
  theta3 *=time_step;
  theta4 *=time_step;

  // matrices for vicous and convective term
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8]; 
    // matrix for SUPG
  MatrixK = LocMatrices[12];
  // matrices for div-div term + 1st extra term
  MatrixS11 = LocMatrices[13];
  MatrixS12 = LocMatrices[14];
  MatrixS13 = LocMatrices[15];
  MatrixS21 = LocMatrices[16];
  MatrixS22 = LocMatrices[17];
  MatrixS23 = LocMatrices[18];
  MatrixS31 = LocMatrices[19];
  MatrixS32 = LocMatrices[20];
  MatrixS33 = LocMatrices[21]; 

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z
  
  // second order derivatives in the residual will be neglected
  // method is for flows with small viscosity

  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);
  //OutPut(supg_params[0] << " " << supg_params[1] << " : " << endl);
  //supg_params[0] = 0;
  // assembling for velocity test functions
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixKRow  = MatrixK[i];
    MatrixS11Row = MatrixS11[i];
    MatrixS12Row = MatrixS12[i];
    MatrixS13Row = MatrixS13[i];
    MatrixS21Row = MatrixS21[i];
    MatrixS22Row = MatrixS22[i];
    MatrixS23Row = MatrixS23[i];
    MatrixS31Row = MatrixS31[i];
    MatrixS32Row = MatrixS32[i];
    MatrixS33Row = MatrixS33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    // test functions for div-div term
    tautest100 = supg_params[1]*test100;
    tautest010 = supg_params[1]*test010;
    tautest001 = supg_params[1]*test001;

    // SUPG term for the time derivative, with scaling, theta1 already scaled with Delta t
    ugradu  = (u1*test100+u2*test010+u3*test001)*supg_params[0]*theta1;
    // velocity-velocity block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      // this block will be multiplied with theta1*Delta t
      // convection 
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // diffusion
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      // add everything
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;

      
      // div-div term
      // store in matrices S
      MatrixS11Row[j] += Mult * tautest100*ansatz100;
      MatrixS12Row[j] += Mult * tautest100*ansatz010;
      MatrixS13Row[j] += Mult * tautest100*ansatz001;
      MatrixS21Row[j] += Mult * tautest010*ansatz100;
      MatrixS22Row[j] += Mult * tautest010*ansatz010;
      MatrixS23Row[j] += Mult * tautest010*ansatz001;
      MatrixS31Row[j] += Mult * tautest001*ansatz100;
      MatrixS32Row[j] += Mult * tautest001*ansatz010;
      MatrixS33Row[j] += Mult * tautest001*ansatz001;
      
      // store in matrices K     
      MatrixKRow[j] += Mult * ansatz000 * ugradu;
      
      // 1st extra term for the time derivative
      // store in matrices S
      val1 = supg_params[0]*theta1*(u1*ansatz100+u2*ansatz010+u3*ansatz001)+ansatz000;
      val = Mult*val1*u1*test100;
      MatrixS11Row[j] += val;
      val = Mult*val1*u1*test010;
      MatrixS12Row[j] += val;
      val = Mult*val1*u1*test001;
      MatrixS13Row[j] += val;
      val = Mult*val1*u2*test100;
      MatrixS21Row[j] += val;
      val = Mult*val1*u2*test010;
      MatrixS22Row[j] += val;
      val = Mult*val1*u2*test001;
      MatrixS23Row[j] += val;
      val = Mult*val1*u3*test100;
      MatrixS31Row[j] += val;
      val = Mult*val1*u3*test010;
      MatrixS32Row[j] += val;
      val = Mult*val1*u3*test001;
      MatrixS33Row[j] += val;

    }                            // endfor j
  }                              // endfor i
}

void TimeNSType4NLVMS_SUPGDD3D_old(double Mult, double *coeff,
              double *param, double hK,
              double **OrigValues, int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **MatrixK, **MatrixC;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, val, val1, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *MatrixKRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixRowC;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  double test000, test100, test010, test001;
  double tautest001, tautest100, tautest010;
  //OutPut("supgnl");
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7, *Orig8, *Orig9, *Orig10;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, c4, c5, c6;
  double u1, u2, u3, px, py, pz;
  double u1_x, u1_y, u1_z;
  double u2_x, u2_y, u2_z;
  double u3_x, u3_y, u3_z;
  double u1_xx, u1_yy, u1_zz;
  double u2_xx, u2_yy, u2_zz;
  double u3_xx, u3_yy, u3_zz;
  double supg_params[2], ugradu;

  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  theta1 *=time_step;
  theta2 *=time_step;
  theta3 *=time_step;
  theta4 *=time_step;

  // matrices for vicous and convective term
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8]; 
  // mass matrix
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  // matrix for PSPG
  MatrixC = LocMatrices[12];
  // matrix for SUPG
  MatrixK = LocMatrices[13];
  // matrices for divergence constraint
  MatrixB1  = LocMatrices[14];
  MatrixB2  = LocMatrices[15];
  MatrixB3  = LocMatrices[16];
  // matrices for pressure term in momentum equations
  MatrixB1T = LocMatrices[17];
  MatrixB2T = LocMatrices[18];
  MatrixB3T = LocMatrices[19];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // p_z
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3
  c4 = coeff[4];                 // f1_old
  c5 = coeff[5];                 // f2_old
  c6 = coeff[6];                 // f3_old

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old
  u1_x = param[3];             // u1old_x
  u2_x = param[4];             // u2old_x
  u3_x = param[5];             // u3old_x
  u1_y = param[6];             // u1old_y
  u2_y = param[7];             // u2old_y
  u3_y = param[8];             // u3old_y
  u1_z = param[9];             // u1old_z
  u2_z = param[10];            // u2old_z
  u3_z = param[11];            // u3old_z
  
  // second order derivatives in the residual will be neglected
  // method is for flows with small viscosity

  u1_xx = u1_yy = u1_zz = 0;
  u2_xx = u2_yy = u2_zz = 0;
  u3_xx = u3_yy = u3_zz = 0;

  //SUPG parameter   
  // supg_params[0] -> for momentum balance tau_m
  // supg_params[1] -> for continuum equ.   tau_c
  SUPG_Param3D(u1, u2, u3, coeff, supg_params);
  //OutPut(supg_params[0] << " " << supg_params[1] << " : " << endl);
  //supg_params[0] = 0;
  // assembling for velocity test functions
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    MatrixKRow  = MatrixK[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    // SUPG term, with scaling, theta1 already with scaled with Delta t
    ugradu  = (u1*test100+u2*test010+u3*test001)*supg_params[0]*theta1;
	
    // test functions for div-div term
    tautest100 = supg_params[1]*test100;
    tautest010 = supg_params[1]*test010;
    tautest001 = supg_params[1]*test001;
    // velocity-velocity block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      // this block will be multiplied with theta1*Delta t
      // convection 
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // diffusion + div-div term
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001)+tautest100*ansatz100;
      // add everything
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100)+tautest100*ansatz010;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100)+tautest100*ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010)+tautest010*ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001)+tautest010*ansatz010;
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010)+tautest010*ansatz001;
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001)+tautest001*ansatz100;
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001)+tautest001*ansatz010;
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001)+tautest001*ansatz001;
      val += val1;
      Matrix33Row[j] += Mult * val;

      // mass matrix
      //val = Mult*(ansatz000*test000);
      //MatrixM11Row[j] += val;
      //MatrixM22Row[j] += val;
      //MatrixM33Row[j] += val;

      // SUPG block, without second order derivative 
      // convection 
      val = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      // test with streamline derivative
      val *= ugradu * Mult;
      // term in the diagonal blocks
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
      // SUPG term for the time derivative
      // store in matrices K
      MatrixKRow[j] += Mult * ansatz000 * ugradu;
    }                            // endfor j

    // pressure-velocity block, these blocks will be multiplied with Delta t
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      // pressure ansatz functions
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];
     
      // pressure term 
      val  = -ansatz000 * test100;
      // SUPG term
      val +=  ansatz100 * ugradu;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz000 * test010;
      // SUPG term
      val +=  ansatz010 * ugradu;
      MatrixRow2[j] += Mult*val;
	  
      val  = -ansatz000 * test001;
      // SUPG term
      val +=  ansatz001 * ugradu;
      MatrixRow3[j] += Mult*val;
    }
  }                              // endfor i

  //supg_params[0] = 0;
  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];
    MatrixRowC = MatrixC[i];

    test100 = Orig4[i];
    test010 = Orig5[i];
    test001 = Orig6[i];
    test000 = Orig7[i];
	
    // pressure-pressure block
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
	
      val = supg_params[0] * time_step * time_step *
	  (ansatz100*test100+ansatz010*test010+ansatz001*test001);
      MatrixRowC[j] += Mult*val;
    }

    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      // divergence constraint
      //val1 =  -c0*(ansatz200+ansatz020+ansatz002);
      val1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      val1 *= theta1;
      val1 += u1;
      val1 *= time_step * supg_params[0];

      val = -test000*ansatz100 + val1*test100;
      MatrixRow1[j] += Mult*val;

      val = -test000*ansatz010 + val1*test010;
      MatrixRow2[j] += Mult*val;
	  
      val = -test000*ansatz001 + val1*test001;
      MatrixRow3[j] += Mult*val;	  	  
    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// ROSENBROCK
// ======================================================================
void TimeNSType1GalerkinJ3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  
  MatrixA = LocMatrices[0];
     
  N_U = N_BaseFuncts[0];
  //N_P = N_BaseFuncts[1];
  
  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
     
  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

void TimeNSGalerkinC3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, val;
  double test000;
  double *Orig3;
  int i, N_U;
  double c4, c5, c6;
  
  //cout << "TimeNSType1GalerkinC" << endl;
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  
  N_U = N_BaseFuncts[0];
  
  Orig3 = OrigValues[0]; // u
     
  c4 = coeff[4]; // dot f1
  c5 = coeff[5]; // dot f2
  c6 = coeff[6]; // dot f3
  
  for(i=0;i<N_U;i++)
  {
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c4;
    Rhs2[i] += Mult*test000*c5;
    Rhs3[i] += Mult*test000*c6;
  } 
}
void TimeNSType3GalerkinJ3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13;
  double **MatrixA21, **MatrixA22, **MatrixA23;
  double **MatrixA31, **MatrixA32, **MatrixA33;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0;
  double u1, u2, u3, u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, u3_x, u3_y, u3_z;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  
  N_U = N_BaseFuncts[0];
  //N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  //Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu
  
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u2old
  u1_x = param[3]; // u1old
  u2_x = param[4]; // u2old
  u3_x = param[5]; // u2old
  u1_y = param[6]; // u1old
  u2_y = param[7]; // u2old
  u3_y = param[8]; // u2old
  u1_z = param[9]; // u1old
  u2_z = param[10]; // u2old
  u3_z = param[11]; // u2old
  

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += u1_x*ansatz000*test000;
      Matrix11Row[j] += Mult * val;

      val  = u1_y*ansatz000*test000;
      Matrix12Row[j] += Mult * val;
      
      val  = u1_z*ansatz000*test000;
      Matrix13Row[j] += Mult * val;
      

      val  = u2_x*ansatz000*test000;
      Matrix21Row[j] += Mult * val;
      
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += u2_y*ansatz000*test000;
      Matrix22Row[j] += Mult * val;
      
      val  = u2_z*ansatz000*test000;
      Matrix23Row[j] += Mult * val;
      
      
      val  = u3_x*ansatz000*test000;
      Matrix31Row[j] += Mult * val;
      
      val  = u3_y*ansatz000*test000;
      Matrix32Row[j] += Mult * val;
      
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += u3_z*ansatz000*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Assembling routine for right-hand sides only
// ======================================================================

// ======================================================================
// right-hand side for NSE ONLY
// ======================================================================
void TimeNSRHS3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3;
  double test000;
  double *Orig0;
  int i, N_U;
  double c1, c2, c3;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    test000 = Orig0[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " "; 
  } // endfor i
}

// ======================================================================
// right-hand side for NSE ONLY 
// ClassicalLES model
// ======================================================================
void TimeNSRHSClassicalLES3D(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, val;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i, N_U;
  double c1, c2, c3;

  double delta, ngu, val1, mu1;
  double D1u1, D2u1, D3u1, D1u2, D2u2, D3u2, D1u3, D2u3, D3u3;
  double a[3][3];
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  D1u1 = param[3]; // D1u1
  D1u2 = param[4]; // D1u2;
  D1u3 = param[5]; // D1u3;
  D2u1 = param[6]; // D2u1
  D2u2 = param[7]; // D2u2;
  D2u3 = param[8]; // D2u2;
  D3u1 = param[9]; // D3u1
  D3u2 = param[10]; // D3u2;
  D3u3 = param[11]; // D3u2;

  // compute Du Du^T
  a[0][0] = D1u1*D1u1 + D2u1*D2u1 + D3u1*D3u1;
  a[0][1] = D1u1*D1u2 + D2u1*D2u2 + D3u1*D3u2;
  a[0][2] = D1u1*D1u3 + D2u1*D2u3 + D3u1*D3u3;
  a[1][0] = a[0][1];
  a[1][1] = D1u2*D1u2 + D2u2*D2u2 + D3u2*D3u2;
  a[1][2] = D1u2*D1u3 + D2u2*D2u3 + D3u2*D3u3;
  a[2][0] = a[0][2];
  a[2][1] = a[1][2];
  a[2][2] = D1u3*D1u3 + D2u3*D2u3 + D3u3*D3u3;
  

  // filter width
  delta =  CharacteristicFilterWidth(hK);
 
  // (delta^2)/(2 gamma) 
  mu1 = 0.5*delta*delta/gamma;
  val1 =  Mult * mu1;

  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    // LES term
    Rhs1[i] += val1*( test100 * a[0][0] + test010 * a[0][1] 
                      + test001 * a[0][2]);
    Rhs2[i] += val1*( test100 * a[1][0] + test010 * a[1][1] 
                      + test001 * a[1][2]);
    Rhs3[i] += val1*( test100 * a[2][0] + test010 * a[2][1] 
                      + test001 * a[2][2]);
  } // endfor i
}

// ======================================================================
// right-hand side for NSE ONLY
// Galdi-Layton model with convolution and auxiliary problem
// ======================================================================
void TimeNSRHSLESModel3D(double Mult, double *coeff,
                         double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;

  double delta, val1,  mu1;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;
  double gdT11, gdT12, gdT13, gdT22, gdT23, gdT33;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  gdT11 = param[0]; // gDeltaT11 or 11 - component of auxiliary problem
  gdT12 = param[1]; // gDeltaT12 or 12 - and 21 - component
  gdT13 = param[2]; // gDeltaT13 or 13 - and 31 - component
  gdT22 = param[3]; // gDeltaT22 or 22 - component of auxiliary problem
  gdT23 = param[4]; // gDeltaT23 or 23 - component of auxiliary problem
  gdT33 = param[5]; // gDeltaT33 or 33 - component of auxiliary problem

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  mu1 = 0.5*delta*delta/gamma;
  val1 =  Mult* mu1;

  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    // LES term
    Rhs1[i] += val1*( test100* gdT11 + test010* gdT12 + test001 * gdT13);
    Rhs2[i] += val1*( test100* gdT12 + test010* gdT22 + test001 * gdT23);
    Rhs3[i] += val1*( test100* gdT13 + test010* gdT23 + test001 * gdT33);
  } // endfor i
}

// ======================================================================
// right-hand side for auxiliary problem 
// Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSGL00AuxProblemRHS3D(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6, val;
  double mat11, mat12, mat13, mat22, mat23, mat33;
  double test000;
  double *Orig0;
  int i, N_U;
  double D1u1, D2u1, D3u1, D1u2, D2u2, D3u2, D1u3, D2u3, D3u3;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  D1u1 = param[0]; // D1u1
  D1u2 = param[1]; // D1u2;
  D1u3 = param[2]; // D1u3;
  D2u1 = param[3]; // D2u1
  D2u2 = param[4]; // D2u2;
  D2u3 = param[5]; // D2u2;
  D3u1 = param[6]; // D3u1
  D3u2 = param[7]; // D3u2;
  D3u3 = param[8]; // D3u2;
  
  // compute Du Du^T
  mat11 = D1u1*D1u1 + D2u1*D2u1 + D3u1*D3u1;
  mat12 = D1u1*D1u2 + D2u1*D2u2 + D3u1*D3u2;
  mat13 = D1u1*D1u3 + D2u1*D2u3 + D3u1*D3u3;
  mat22 = D1u2*D1u2 + D2u2*D2u2 + D3u2*D3u2;
  mat23 = D1u2*D1u3 + D2u2*D2u3 + D3u2*D3u3;
  mat33 = D1u3*D1u3 + D2u3*D2u3 + D3u3*D3u3;

  for(i=0;i<N_U;i++)
  {
    test000 = Orig0[i];

    val = Mult*test000;

    Rhs1[i] += val* mat11;
    Rhs2[i] += val* mat12;
    Rhs3[i] += val* mat13;
    Rhs4[i] += val* mat22;
    Rhs5[i] += val* mat23;
    Rhs6[i] += val* mat33;
  } // endfor i
}

// ======================================================================
// right-hand side ONLY for auxiliary problem applied to velocity
// ======================================================================
void TimeNSRHSAuxProblemU(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3;
  double test000;
  double *Orig0;
  int i, N_U;
  double c1, c2, c3;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  c1 = param[0]; // u1
  c2 = param[1]; // u2
  c3 = param[2]; // u3

  for(i=0;i<N_U;i++)
  {
    test000 = Orig0[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " "; 
  } // endfor i
}

// ======================================================================
// right-hand side for additional terms in rhs of small scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_SmallRhs3D(double Mult, double *coeff, 
                           double *param, double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  double u1, u2, u3,  ho_u1, ho_u2, ho_u3;
  double D1u1, D2u1, D3u1, D1u2, D2u2, D3u2, D1u3, D2u3, D3u3;
  double ho_D1u1, ho_D2u1, ho_D3u1, ho_D1u2, ho_D2u2, ho_D3u2;
  double ho_D1u3, ho_D2u3, ho_D3u3, p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  // large scales
  u1 = param[0];
  u2 = param[1];
  u3 = param[2];
  D1u1 = param[3]; // D1u1
  D1u2 = param[4]; // D1u2;
  D1u3 = param[5]; // D1u3;
  D2u1 = param[6]; // D2u1
  D2u2 = param[7]; // D2u2;
  D2u3 = param[8]; // D2u2;
  D3u1 = param[9]; // D3u1
  D3u2 = param[10]; // D3u2;
  D3u3 = param[11]; // D3u2;
  // small scales
  ho_u1 = param[12];
  ho_u2 = param[13];
  ho_u3 = param[14];
  ho_D1u1 = param[15]; // D1u1
  ho_D1u2 = param[16]; // D1u2;
  ho_D1u3 = param[17]; // D1u3;
  ho_D2u1 = param[18]; // D2u1
  ho_D2u2 = param[19]; // D2u2;
  ho_D2u3 = param[20]; // D2u2;
  ho_D3u1 = param[21]; // D3u1
  ho_D3u2 = param[22]; // D3u2;
  ho_D3u3 = param[23]; // D3u2;
  // pressure
  p = param[24];

  c0 = coeff[0];
  
  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += 2*c0*(D1u1*test100+(D1u2+D2u1)*test010/4+(D1u3+D3u1)*test001/4);
    Rhs1[i] += (u1*D1u1+u2*D2u1+u3*D3u1) * test000;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1+u3*ho_D3u1) * test000;
    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1+ho_u3*D3u1) * test000;
    Rhs1[i] -= p * test100;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(D2u2*test010+(D1u2+D2u1)*test100/4+(D2u3+D3u2)*test001/4);
    Rhs2[i] += (u1*D1u2+u2*D2u2+u3*D3u2) * test000;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2+u3*ho_D3u2) * test000;
    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2+ho_u3*D3u2) * test000;
    Rhs2[i] -= p * test010;
    Rhs2[i] *= Mult;
    Rhs3[i] += 2*c0*(D3u3*test001+(D1u3+D3u1)*test100/4+(D2u3+D3u2)*test010/4);
    Rhs3[i] += (u1*D1u3+u2*D2u3+u3*D3u3) * test000;
    Rhs3[i] += (u1*ho_D1u3+u2*ho_D2u3+u3*ho_D3u3) * test000;
    Rhs3[i] += (ho_u1*D1u3+ho_u2*D2u3+ho_u3*D3u3) * test000;
    Rhs3[i] -= p * test001;
    Rhs3[i] *= Mult;
  } // endfor i
}

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// used for : GALERKIN
//            COLETTI (non linear steps)
// ========================================================================
void TimeNSParamsVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old
}

// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping
}

// ========================================================================
// parameters: gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsGradVelo3D(double *in, double *out)
{
  //cout << "GRAD" << endl;
  out[0] = in[6]; // D100(u1old)
  out[1] = in[7]; // D100(u2old)
  out[2] = in[8]; // D100(u3old)
  out[3] = in[9]; // D010(u1old)
  out[4] = in[10]; // D010(u2old)
  out[5] = in[11]; // D010(u3old)
  out[6] = in[12]; // D001(u1old)
  out[7] = in[13]; // D001(u2old)
  out[8] = in[14]; // D001(u3old)
  // cout << in[4] << " E " << in[5] << " " << in[6] << " " << in[7] << endl;
}

// ========================================================================
// parameters: u1old, u2old, u3old
// all partial derivatives
// convolution of u1old, u2old, u3old
// ========================================================================
void TimeNSParamsVelo_GradVelo_ConvVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[15]; // convolution of u1old
  out[13] = in[16]; // convolution of u2old
  out[14] = in[17]; // convolution of u3old

  out[15] = in[2]; // z - coordinate for van Driest damping
}

// ========================================================================
// used for assembling of term coming from the LES model
// GL00AuxProb and GL00Conv
// parameters used in TimeNSRHSLESModel3D
// ========================================================================
void TimeNSParamsRHSLES3D(double *in, double *out)
{
  // components of convolved tensor
  out[0] = in[3]; // g_d\ast(D1u1*D1u1+D2u1*D2u1+D3u1*D3u1)
  out[1] = in[4]; // g_d\ast(D1u1*D1u2+D2u1*D2u2+D3u1*D3u2)
  out[2] = in[5]; // g_d\ast(D1u1*D1u3+D2u1*D2u3+D3u1*D3u3)
  out[3] = in[6]; // g_d\ast(D1u2*D1u2+D2u2*D2u2+D3u2*D3u2)
  out[4] = in[7]; // g_d\ast(D1u2*D1u3+D2u2*D2u3+D3u2*D3u3)
  out[5] = in[8]; // g_d\ast(D3u1*D3u1+D3u2*D3u2+D3u3*D3u3)

}

// ========================================================================
// used for : classical LES, first iteration step, viscosity type = 4
// ========================================================================
void TimeNSParamsVelo_GradVeloNuT4_3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  // convolution of the solution
  out[12] = in[15]; // g_\delta \ast u1
  out[13] = in[16]; // g_\delta \ast u2
  out[14] = in[17]; // g_\delta \ast u3
}

// ========================================================================
// used for : GL00Convolution, rhs assembling, viscosity type = 4
// ========================================================================
void TimeNSParamsRHSGL00ConvolutionNuT4_3D(double *in, double *out)
{
  // grad u
  out[0] = in[6]; // D1u1
  out[1] = in[7]; // D1u2
  out[2] = in[8]; // D1u3
  out[3] = in[9]; // D2u1
  out[4] = in[10]; // D2u2
  out[5] = in[11]; // D2u3
  out[6] = in[12]; // D3u1
  out[7] = in[13]; // D3u2
  out[8] = in[14]; // D3u3

  // components of convolved tensor or solution of auxiliary problem
  out[9] = in[15]; // g_d\ast(D1u1*D1u1+D2u1*D2u1+D3u1*D3u1)
  out[10] = in[16]; // g_d\ast(D1u1*D1u2+D2u1*D2u2+D3u1*D3u2)
  out[11] = in[17]; // g_d\ast(D1u1*D1u3+D2u1*D2u3+D3u1*D3u3)
  out[12] = in[18]; // g_d\ast(D1u2*D1u2+D2u2*D2u2+D3u2*D3u2)
  out[13] = in[19]; // g_d\ast(D1u2*D1u3+D2u2*D2u3+D3u2*D3u3)
  out[14] = in[20]; // g_d\ast(D3u1*D3u1+D3u2*D3u2+D3u3*D3u3)

  // velocity 
  out[15] = in[3]; // u1
  out[16] = in[4]; // u2
  out[17] = in[5]; // u3

  // convolved velocity
  out[18] = in[21]; //  g_d\ast u1
  out[19] = in[22]; //  g_d\ast u2
  out[20] = in[23]; //  g_d\ast u3
  
}

// ========================================================================
// parameters: 
// used for : GL00AuxProblem, viscosity type = 4
// ========================================================================
void TimeNSParamsGL00AuxProblemNuT4_3D(double *in, double *out)
{

  OutPut("TimeNSParamsGL00AuxProblemNuT43D not implemented " << endl);
  exit(4711);
  // \nabla u
  out[0] = in[4]; // D1u1
  out[1] = in[5]; // D1u2
  out[2] = in[6]; // D2u1
  out[3] = in[7]; // D2u2

  //  solution of auxiliary problem
  out[4] = in[ 8]; // sol11
  out[5] = in[ 9]; // sol12 = sol21
  out[6] = in[10]; // sol22

  // solution
  out[7] = in[2]; // u1
  out[8] = in[3]; // u2

  // convolution of the solution
  out[9] = in[11]; // g_\delta \ast u1
  out[10] = in[12];// g_\delta \ast u2
}

// ========================================================================
// used for VMS, assembling of rhs for small scale equation 
// ========================================================================
void TimeNSParams_VMS_SmallRhs3D(double *in, double *out)
{
  // large scales
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  // small scales
  out[12] = in[15]; // u1old
  out[13] = in[16]; // u2old
  out[14] = in[17]; // u3old

  out[15] = in[18]; // D1u1
  out[16] = in[19]; // D1u2
  out[17] = in[20]; // D1u3
  out[18] = in[21]; // D2u1
  out[19] = in[22]; // D2u2
  out[20] = in[23]; // D2u3
  out[21] = in[24]; // D3u1
  out[22] = in[25]; // D3u2
  out[23] = in[26]; // D3u3
  // large pressure
  out[24] = in[27]; // p
}
// ========================================================================
// parameters: u1old, u2old, G^H
// used for : projection-based VMS
// ========================================================================
void TimeNSParamsVelo_GradVelo_LargeScale3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping

  // finest grid
  if (TDatabase::ParamDB->INTERNAL_LEVEL == 1)
  {
      out[15] = in[15]; // G_11
      out[16] = in[16]; // G_12
      out[17] = in[17]; // G_13
      out[18] = in[18]; // G_22
      out[19] = in[19]; // G_23
      out[20] = in[20]; // G_33
    
      out[21] = in[21]; // coarse space
  }
  else
  {
      // coarser grids
      out[15] = out[16] = out[17] = out[18] = out[19] = out[20] = out[21] = 0;
  }
}


void TimeNSParamsVelo_GradVelo_VMS3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3
}

//===================used for ALEVMS3D====================

void TimeNSParamsVelo_GradVelo_VMS3D_ALE(double *in, double *out)
{
  out[0] = in[3]-in[15]; // u1old-w1
  out[1] = in[4]-in[16]; // u2old-w2
  out[2] = in[5]-in[17]; // u3old-w3

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3
}

 
