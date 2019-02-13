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
// @(#)ConvDiff.C
//
// Purpose: contains routines which are relevant for ConvDiff2D.C and
//          ConvDiff3D.C
//
// Author: Volker John
//
// History: start of implementation 09.08.2011
//
// =======================================================================

#include <Constants.h>
#include <Database.h>
#include <MainUtilities.h>
#include <LinAlg.h>
#include <ConvDiff.h>
#ifdef __2D__
  #include <FEFunction2D.h>
#endif

#include <stdlib.h>
// #include <malloc.h>
#include <string.h>
#include <sstream>

/******************************************************************************/
// SetParametersCD
// sets parameters of the data base for the main programs CDAdapt2D.C
// and CD_3D.C
/******************************************************************************/

void SetParametersCD(int &nonlinear_method)
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== KLR02_3)
    TDatabase::ParamDB->SOLD_S = 0;
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== LP96)
  {
    OutPut("SOLD_PARAMETER_TYPE == LP96 should be used with higher quadrature rule,"<<endl);
    OutPut("since right hand side is in general not linear !!!"<<endl);
  }

  if ((TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)&&
    (TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION_2_LEVEL))
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
  }
  else
  {
    if (TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
    {
      // check spaces and change if necessary
      switch(TDatabase::ParamDB->ANSATZ_ORDER)
      {
        case 1:
          TDatabase::ParamDB->ANSATZ_ORDER = 100;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 2:
          TDatabase::ParamDB->ANSATZ_ORDER = 201;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 3:
          TDatabase::ParamDB->ANSATZ_ORDER = 302;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 4:
          TDatabase::ParamDB->ANSATZ_ORDER = 403;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 5:
          TDatabase::ParamDB->ANSATZ_ORDER = 504;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        default:
          break;
      }
    }
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
  {
    if(TDatabase::ParamDB->LP_STREAMLINE)
    {
      TDatabase::ParamDB->LP_STREAMLINE = 0;
      TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
      TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
      OutPut("local projection stabilisation in streamline direction ");
      OutPut("is switched off due to stabilisation of full gradient." << endl);
    }
  }

  if(TDatabase::ParamDB->LP_STREAMLINE)
  {
    if(TDatabase::ParamDB->LP_FULL_GRADIENT)
    {
      TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
      TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
      TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;
      OutPut("local projection stabilisation for gradient ");
      OutPut("is switched off due to stabilisation of streamline direction." << endl);
    }
  }

  if (TDatabase::ParamDB->LP_STREAMLINE && TDatabase::ParamDB->LP_CROSSWIND)
    nonlinear_method = 1;

  if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  // has to be changed for all discretizations
  // otherwise access to not allocated array in error computation
  if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
  {
    TDatabase::ParamDB->SDFEM_TYPE = 2;
    OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
      << " since no adjoint problem is solved !!! "<<endl);
  }
  if ((TDatabase::ParamDB->SDFEM_TYPE != 100)&&(TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
    (TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
    TDatabase::ParamDB->SDFEM_TYPE = 100;
    OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
      << " since adjoint problem is solved !!! "<<endl);
  }

  // SUPG, GLS
  // there is only one minor difference in method, treat them otherwise the same
  if (TDatabase::ParamDB->DISCTYPE==GLS)
  {
    TDatabase::ParamDB->DISCTYPE=SUPG;
    TDatabase::ParamDB->INTERNAL_DISC_FLAG = 1;
  }
  if ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE==0))
  {
    // this excludes some not wished side effects
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
  }
  if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&(TDatabase::ParamDB->SOLD_ADJOINT))
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 100;

  if  ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD))
  {
    TDatabase::ParamDB->SDFEM_TYPE = 0;
    TDatabase::ParamDB->DELTA0 =  TDatabase::ParamDB->DELTA1 = 0;
    OutPut("FEM-TVD: switched stabilization off!" << endl);
  }

  if ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE))
    nonlinear_method = 1;
  if ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&((TDatabase::ParamDB->SDFEM_TYPE==2)
    || (TDatabase::ParamDB->SDFEM_TYPE==100)))
    SetPolynomialDegree();

  if ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SDFEM_TYPE==2))
  {
    if (TDatabase::ParamDB->CELL_MEASURE != 4)
    {
      TDatabase::ParamDB->CELL_MEASURE = 4;
      OutPut("CELL_MEASURE changed to " << TDatabase::ParamDB->CELL_MEASURE << endl);
    }
    if (TDatabase::ParamDB->DELTA0 != 1.0)
    {
      TDatabase::ParamDB->DELTA0 = 1.0;
      OutPut("DELTA0 changed to " << TDatabase::ParamDB->DELTA0 << endl);
    }
  }

  if (TDatabase::ParamDB->DISCTYPE==CIP)
  {
    TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 1;
    if (TDatabase::ParamDB->CIP_TYPE < 10)
    {
      TDatabase::ParamDB->DISCTYPE=GALERKIN;
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
    }
    else
    {
      // add SOLD term to CIP discretization
      nonlinear_method = 1;
      TDatabase::ParamDB->CIP_TYPE -=10;
      TDatabase::ParamDB->DISCTYPE=SUPG;
      if (TDatabase::ParamDB->SOLD_TYPE == 0)
        TDatabase::ParamDB->SOLD_TYPE = 2;
      // switch off SUPG part
      TDatabase::ParamDB->SDFEM_TYPE = 0;
      TDatabase::ParamDB->DELTA0 = 0.0;
      TDatabase::ParamDB->DELTA1 = 0.0;
      // 	TDatabase::ParamDB->SOLD_PARAMETER_TYPE = BE05_1;
      //	TDatabase::ParamDB->SOLD_TYPE = 0;
    }
  }

  if (TDatabase::ParamDB->DISCTYPE==DG)
  {
    TDatabase::ParamDB->DISCTYPE=GALERKIN;
    TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 2;
    if ( TDatabase::ParamDB->ANSATZ_ORDER < 10)
      TDatabase::ParamDB->ANSATZ_ORDER = -TDatabase::ParamDB->ANSATZ_ORDER-10;
    else
      // P elements on quads
      TDatabase::ParamDB->ANSATZ_ORDER = -10*TDatabase::ParamDB->ANSATZ_ORDER;
    if (TDatabase::ParamDB->ESTIMATE_ERRORS)
    {
      TDatabase::ParamDB->ESTIMATE_ERRORS = 0;
      OutPut("Error estimation does not work for DG !!!"<< endl);
    }
  }
  if  (!(TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
    TDatabase::ParamDB->SOLD_TYPE = 0;
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE =0;
  }
  if (TDatabase::ParamDB->DISCTYPE==LOCAL_PROJECTION_2_LEVEL)
  {
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    OutPut("LOCAL_PROJECTION_2_LEVEL only with direct solver !!!" << endl);
    nonlinear_method = 1;
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
  }
}



double Mesh_size_in_convection_direction(double hK, double b1, double b2)
{
  if(TDatabase::ParamDB->INTERNAL_HK_CONVECTION < 0)
  {
    // not yet computed for this mesh cell
    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = 
      Mesh_size_in_convection_direction_without_storing(hK, b1, b2);
  }
  // else: already computed for this mesh cell
  return TDatabase::ParamDB->INTERNAL_HK_CONVECTION;
}


double Mesh_size_in_convection_direction_without_storing(double hK, double b1,
                                                         double b2)
{
  int i;
  double x[4], y[4], sx, sy, a[16], b[16], den, val, norm_b;

  // triangles
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3] == -4711)
  {
    for (i=0;i<3;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
    }
    // initialize rhs
    memset(b,0,9*SizeOfDouble);

    // set matrices for computation of the coefficients
    // of the linear function
    for (i=0;i<3;i++)
    {
      a[3*i] = 1;
      a[3*i+1] = x[i];
      a[3*i+2] = y[i];
      b[4*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,3,3,3,3);

    // compute numerator
    norm_b = sqrt(b1*b1 + b2*b2);
    // compute denominator
    den = 0;
    for (i=0;i<3;i++)
    {
      // value of gradient basis fct. in bary centre
      // is a constant
      den += fabs(b1*b[3*i+1]+b2*b[3*i+2]);
    }
    // return the mesh size in convection direction
    if (den<1e-10)
      return(hK);
    else
      return(2*norm_b/den);
  }
  else
  {                                               // quadrilateral
    sx = sy = 0;
    for (i=0;i<4;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      //OutPut(x[i] <<  " " << y[i] << " ");
      sx += x[i];
      sy += y[i];
    }
    //bary centre
    sx /= 4;
    sy /= 4;
    // initialize rhs
    memset(b,0,16*SizeOfDouble);

    // set matrices for computation of the coefficients
    // of the bilinear function
    for (i=0;i<4;i++)
    {
      a[4*i] = 1;
      a[4*i+1] = x[i];
      a[4*i+2] = y[i];
      a[4*i+3] = x[i]*y[i];
      b[5*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,4,4,4,4);

    // compute numerator
    norm_b = sqrt(b1*b1 + b2*b2);
    // compute denominator
    den = 0;
    for (i=0;i<4;i++)
    {
      // value of gradient basis fct. in bary centre
      val = b1*(b[4*i+1] + b[4*i+3] * sy);
      val += b2*(b[4*i+2] + b[4*i+3] * sx);
      den += fabs(val);
    }
    // return the mesh size in convection direction
    //OutPut(b1 << " " << b2 << " " << fabs(den) << " " << 2*norm_b/fabs(den) << " " );
    if (den<1e-10)
    {
      return(hK);
    }
    else
    {
      return(2*norm_b/den);
    }
  }
}

double Mesh_size_in_convection_direction(double hK, double b1, double b2, 
                                         double b3)
{
  if(TDatabase::ParamDB->INTERNAL_HK_CONVECTION < 0)
  {
    // not yet computed for this mesh cell
    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = 
      Mesh_size_in_convection_direction_without_storing(hK, b1, b2, b3);
  }
  // else: already computed for this mesh cell
  return TDatabase::ParamDB->INTERNAL_HK_CONVECTION;
}

double Mesh_size_in_convection_direction_without_storing(double hK, double b1, double b2, 
           double b3)
{
  int i;
  double x[8], y[8], z[8], sx, sy, sz, a[64], b[64], den, val, norm_b;

  // tetrahedra
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
  {
    for (i=0;i<4;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
    }
    // initialize rhs
    memset(b,0,16*SizeOfDouble);

    // set matrices for computation of the coefficients
    // of the bilinear function
    for (i=0;i<4;i++)
    {
      a[4*i] = 1;
      a[4*i+1] = x[i];
      a[4*i+2] = y[i];
      a[4*i+3] = z[i];
      b[5*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    // which is faster ?
    // SolveMultipleSystemsLapack(a,b,4,4,4,4);
    SolveMultipleSystems(a,b,4,4,4,4);
    
    // compute numerator
    norm_b = sqrt(b1*b1 + b2*b2 + b3*b3);
    // compute denominator 
    den = 0;
    for (i=0;i<4;i++)
    {
      // value of gradient basis fct. in bary centre
      // is a constant
      den += fabs(b1*b[4*i+1]+b2*b[4*i+2]+b3*b[4*i+3]);
    } 
    // return the mesh size in convection direction
    if (den<1e-10)
    {
      return(hK);
    }
    else
    {
      return(2*norm_b/den);
    }
  }
  else
  { // hexahedron 
    sx = sy = sz = 0;
    for (i=0;i<8;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
      //OutPut(x[i] <<  " " << y[i] << " ");
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
    for (i=0;i<8;i++)
    {
      b[9*i] = 1;
    }
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
    }

    // solve system for the coefficients of the bilinear function
    // which is faster ??
    //SolveMultipleSystemsLapack(a,b,8,8,8,8);
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
    {
      return(hK);
    }
    else
    {
      return(2*norm_b/den);
    }
  }
}



/******************************************************************************/
//
// definitions of parameters for the SUPG method
// steady-state and time-dependent convection-diffusion-reaction equation
//
/******************************************************************************/

#ifdef __2D__
double Compute_SDFEM_delta(double hK, double eps, double b1, double b2,
double react, double linfb)
#endif
#ifdef __3D__
double Compute_SDFEM_delta(double hK, double eps, double b1, double b2, double b3,
double react, double linfb)
#endif
{
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  double alpha, alpha2, delta, norm_b, nu, reaction, h_K;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  int i;
  
  if (fabs(eps)<1e-20)
    eps = 1e-20;

  // compute cell diameter in convection direction

  if(TDatabase::ParamDB->CELL_MEASURE==4)
    h_K = Mesh_size_in_convection_direction(hK, b1, b2
      #ifdef __3D__
      , b3
      #endif
    );
  else
    h_K = hK;

  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 0:
    case 1:
    case 2:
    case 5:
    case 6:
    case 7:
    case 8:
    case 11:
      norm_b = b1*b1+b2*b2;
      #ifdef __3D__
      norm_b += b3*b3;
      #endif
      norm_b = sqrt(norm_b);
      break;
  }

  // just for safety
  reaction = fabs(react);

  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 0:
      // version from Roos, Stynes, Tobiska 2006
      if (!TDatabase::ParamDB->SHISHKIN_MESH)
      {
        if(eps < h_K*norm_b)
          delta = delta0 * h_K/norm_b;
        else
          delta = delta1 *h_K*h_K/eps ;
      }
      else
      {
        ErrMsg("Shishkin meshes are currently not supported !!!"
               << TDatabase::ParamDB->SHISHKIN_MESH);
        exit(4711);
      }
      /* else                                        // delta for SDFEM only in coarse part of Shishkin mesh
      {
        if(h_K > bound)
          delta = delta0 * h_K/linfb;
        else
          delta = 0;
      }*/
      break;
    case 1:
    case 2:
      // delta based on 1d Green's formula
      // case 2 is the the standard choice, in this choice
      // - h_K is mesh size in convection direction
      // - delta0 = 1
      // are forced
      nu = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      if (norm_b > 0)
      {
        alpha = nu*norm_b*h_K/(2*eps);
        delta = nu*h_K*(1/tanh(alpha) - 1/alpha)/(2*norm_b);
      }
      else
        delta = 0;
      // scaling
      delta *= delta0;
      // for computation of reference solutions for parameter optimization
      //if (TDatabase::ParamDB->INTERNAL_MOMENT == 10)
    	//  delta = 0;
      break;

      // THESE ARE THE PARAMETERS WHERE REACTION BECOMES IMPORTANT
      // can be used for steady-state or time-dependent case
    case 4:
      // Lube/Rapin (2006)
      // equations with reaction
      // under the assumption that the inverse inequality holds
      // with nu_inv = 0 (e.g. linear and bilinear elements)
      // polynomial degree = 1

      delta = -1;

      if (linfb > 0)
        delta = h_K/linfb;

      if (reaction > 0)
      {
        alpha = 1.0/reaction;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      if (eps>0)
      {
        alpha = h_K*h_K/eps;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      delta *= delta0;
      if (delta<0)
        delta = 0;
      break;
    case 5:
    case 8:
      // Lube/Rapin (2006), slightly modified
      // equations with reaction
      // under the assumption that the inverse inequality holds
      // with nu_inv = 0 (e.g. linear and bilinear elements)
      // polynomial degree = 1
      // !!! case 8: reaction without term 1
      delta = -1;
      if (norm_b > 0)
        delta = h_K/(2*norm_b);

      if (reaction > 0)
      {
        alpha = 1.0/reaction;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      if (eps>0)
      {
        alpha = h_K*h_K/eps;
        if ((delta < 0) || (alpha<delta))
          delta = alpha;
      }

      delta *= delta0;
      if (delta<0)
        delta = 0;
      break;

    case 6:
      // Franca/Valentin (2000)
      // equations with reaction
      // under the assumption that the inverse inequality holds
      // with nu_inv = 1/3 (e.g. linear and bilinear elements)
      // polynomial degree = 1
      if (reaction<1e-20)
        reaction = 1e20;

      alpha  = 6*eps/(h_K*h_K*reaction);
      if(alpha <= 1)
        alpha = 1;

      alpha2 = h_K*norm_b/(3*eps);
      if(alpha2 <= 1)
        alpha2 = 1;

      delta=1/(reaction*h_K*h_K*alpha+6*eps*alpha2);
      delta *= delta0*h_K*h_K;
      break;

    case 7:
      // Codina (2000)
      // equations with reaction
      delta = delta0 * h_K * h_K/(4*eps+2 * h_K * norm_b + h_K*h_K * reaction);
      break;

      // THESE ARE PARAMETERS FOR THE TIME-DEPENDENT CASE
    case 9:
      // first estimate in paper John, Novo, SIMUM 2011, formula (3.3)
      delta = delta0*time_step;
      break;
    case 10:
      // second estimate in paper John, Novo, SINUM 2011, asymptotic of formula (5.1)
      // get unscaled diffusion
      eps /= (theta1*time_step);
      if(eps <= h_K)
        delta = delta0 * h_K;
      else
        delta = delta0 * h_K*h_K/eps;
      break;
    case 11:
      // for estimate in paper John, Novo, SINUM 2011, formula (3.10)
      // parameter depending on square root of time step
      norm_b /= time_step;
      if (norm_b >0)
        delta = delta0 * h_K * sqrt(time_step)/norm_b;
      else
        delta = 0;
      break;
    case 100:
      // take value from piecewise constant field
      // only on finest level available
      i = TDatabase::ParamDB->INTERNAL_LOCAL_DOF;
      //OutPut(i << " " << TDatabase::ParamDB->INTERNAL_P1_Array[i] << endl);
      delta =  TDatabase::ParamDB->INTERNAL_P1_Array[i];
      break;
    default :
      OutPut("SDFEM_TYPE "<<TDatabase::ParamDB->SDFEM_TYPE<<
        " not implemented !!!" << endl);
      exit(4711);
  }

  if (TDatabase::ParamDB->SC_VERBOSE>2)
  {
    OutPut("delta " << delta << " " << endl);
  }
  if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 120814)
  	TDatabase::ParamDB->INTERNAL_P1_Array[TDatabase::ParamDB->INTERNAL_MOMENT] = delta;

  return(delta);
}


// ========================================================================
// compute parameters for SOLD schemes
// 2D stationary cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_xx
//        param[4] = u_yy
//        param[5] = ||u^h||_{H^1,K}
//        param[6] = ||R(u^h)||_{L^2,K}
// 2D time-dependent cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_xx
//        param[4] = u_yy
//        param[5] = u_old
//        param[6] = u_old_x
//        param[7] = u_old_y
//        param[8] = u_old_xx
//        param[9] = u_old_yy
// 3D stationary cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_z
//        param[4] = ||u^h||_{H^1,K}
//        param[5] = ||R(u^h)||_{L^2,K}
// 3D time-dependent cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_z
//        param[4] = u_old
//        param[5] = u_old_x
//        param[6] = u_old_y
//        param[7] = u_old_z
//
// The implementation is for the 3D case. The unnecessary flops for 2D
// should be neglibible in the computing time
// ========================================================================

#ifdef __2D__
double Compute_SOLD_sigma(double hK, double eps, double b1,
double b2, double c, double f,
double linfb, double deltaK, double *param,
double residual, int residual_computed,
int time_dependent_problem)
#endif
#ifdef __3D__
double Compute_SOLD_sigma(double hK, double eps, double b1,
double b2, double b3, double c, double f,
double linfb, double deltaK, double *param,
double residual, int residual_computed,
int time_dependent_problem)
#endif
{
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE, i, N;
  double u_x, u_y, u_z, norm_u, norm_res, sigma, res=0.0, norm_b2, value;
  double b1_orth, b2_orth,  b3_orth, norm_der_u2, linfb_orth, z1, z2, z3, linfz, normz;
  double alpha, beta, gamma, lambda, kappa, omega, rho, norm_b_orth;
  double epsilon= 1e-10, hK_project, y, z, u_xx, u_yy;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
#ifdef __2D__
  static double b3 = 0.0;
  u_z = 0.0;
#endif
#ifdef __3D__
  u_z = param[3];
  u_xx = 0.0;
  u_yy = 0.0;
#endif
  // square of norm of convection
  norm_b2 =  b1*b1 + b2*b2 + b3*b3;

  // the residual is already computed
  // if neither residual nor params are provided: initialization is res=0
  if (residual_computed)
    res = residual;

  // parameters are provided
  if (param!=NULL)
  {
    u_x = param[1];
    u_y = param[2];
    // square of the norm of the derivative of the current solution
    norm_der_u2 = u_x*u_x + u_y*u_y + u_z*u_z;

    // compute the residual if this is not already done
    // in 3D, the Laplacian is not provided and the corresponding factors
    // were set to be zero
    if (!residual_computed)
      // only for stationary problems
      res = - eps*(u_xx + u_yy) + b1*u_x + b2*u_y + +b3 *u_z + c*param[0] - f;
  }

  // compute the parameter for the SOLD scheme
  switch (sold_parameter_type)
  {
    case JSW87:                                   // Johnson,Schatz,Wahlbin (1987) (linear)
#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2);
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2, b3);
#endif
      sigma = sqrt(norm_b2)*hK_project*sqrt(hK_project) - eps;
      if (sigma < 0)
        sigma = 0;
      break;

    case HMM86:                                   // Hughes, Mallet, Mizukami (1986)
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      value = b1*u_x + b2*u_y + b3* u_z;          // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      b3_orth = value * u_z/ norm_der_u2;
      // ||b_orth||_\infty
      if (fabs(b2_orth)< fabs (b1_orth))
        linfb_orth = fabs(b1_orth);
      else
        linfb_orth = fabs(b2_orth);
      if (linfb_orth < fabs(b3_orth))
        linfb_orth = fabs(b3_orth);
      // \tau(\b_orth)
#ifdef __2D__
      value = Compute_SDFEM_delta(hK, eps, b1_orth, b2_orth, c, linfb_orth);
#endif
#ifdef __3D__
      value = Compute_SDFEM_delta(hK, eps, b1_orth, b2_orth, b3_orth, c, linfb_orth);
#endif
      // sigma = max ( \tau(b_orth) - \tau(b))
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
      //OutPut("sigma " << sigma << " orth " << value << " deltaK " << deltaK << endl);
      //OutPut("u_x " << u_x << " u_y " << u_y << endl);
      value = b1*u_x + b2*u_y + b3*u_z;           // b \cdot nabla u^h
      if (norm_der_u2>0)
        // sigma = sigma * residual * (b\cdot u^h)/||\nabla u^h||^2
        sigma *= res*value/norm_der_u2;
      else
        sigma = 0;
      break;

    case TP86_1:                                  // Tezduyar, Park (1986), first parameter
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      alpha = b1*u_x + b2*u_y + b3*u_z;           // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;
      b3_orth = alpha * u_z/ norm_der_u2;

      rho =  sqrt(b1_orth*b1_orth + b2_orth*b2_orth + b3_orth*b3_orth);
      value = rho/sqrt(norm_b2);
      value = 2*value*(1-value);
#ifdef __2D__
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth);
#endif
#ifdef __3D__
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth,b3_orth);
#endif
      lambda = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      sigma = lambda * kappa * value/(2*rho)*res*alpha/norm_der_u2;
      break;

    case TP86_2:                                  // Tezduyar, Park (1986), second parameter
      if (TDatabase::ParamDB->SOLD_U0==0)
      {
        OutPut("Paramter SOLD_U0 is zero " << endl);
        exit(4711);
      }
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      alpha = b1*u_x + b2*u_y + b3*u_z;           // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;
      b3_orth = alpha * u_z/ norm_der_u2;

      rho =  sqrt(b1_orth*b1_orth + b2_orth*b2_orth + b3_orth*b3_orth);
      value = rho/sqrt(norm_b2);
      value = 2*value*(1-value);
#ifdef __2D__
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth);
#endif
#ifdef __3D__
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth,b3_orth);
#endif
      lambda = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      sigma = lambda * kappa * kappa* value/(2*rho)*res*alpha/sqrt(norm_der_u2);
      sigma /= TDatabase::ParamDB->SOLD_U0;
      break;

    case GdC88:                                   // Galeao, do Carmo (1988)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      if (time_dependent_problem)
      {
        z1 *= theta1*time_step;
        z2 *= theta1*time_step;
        z3 *= theta1*time_step;
      }
      if (fabs(z2)< fabs (z1))
        linfz = fabs(z1);
      else
        linfz = fabs(z2);
      if (linfz < fabs(z3))
        linfz = fabs(z3);
#ifdef __2D__
      value = Compute_SDFEM_delta(hK, eps, z1, z2, c, linfz);
#endif
#ifdef __3D__
      value = Compute_SDFEM_delta(hK, eps, z1, z2, z3, c, linfz);
#endif
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      break;

    case dCG91:                                   // do Carmo, Galeao (1991)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      normz= sqrt(z1*z1+z2*z2+z3*z3);
      if (normz==0)
      {
        sigma = 0;
        break;
      }
      if (time_dependent_problem)
        normz *= theta1*time_step;
      sigma = deltaK*(sqrt(norm_b2)/normz-1);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      //OutPut(deltaK << " " << sqrt(norm_b2) << " " << normz << " " << sigma << endl);
      break;

    case dCA03:                                   // do Carmo, Alvarez (2003)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      if (time_dependent_problem)
      {
        z1 *= theta1*time_step;
        z2 *= theta1*time_step;
        z3 *= theta1*time_step;
      }
      if (fabs(z2)< fabs (z1))
        linfz = fabs(z1);
      else
        linfz = fabs(z2);
      if (linfz < fabs(z3))
        linfz = fabs(z3);
#ifdef __2D__
      value = Compute_SDFEM_delta(hK, eps, z1, z2, c, linfz);
#endif
#ifdef __3D__
      value = Compute_SDFEM_delta(hK, eps, z1, z2, z3, c, linfz);
#endif
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }

      // sigma of [GdC88] is computed, now compute rho
      normz = sqrt(z1*z1+z2*z2+z3*z3);

      alpha = normz/sqrt(norm_b2);
      if (alpha >= 1)
      {
        rho = 1;
      }

#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2);
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2,b3);
#endif
      beta = pow(hK_project,1-alpha*alpha);
      if (beta > 1)
        beta = 1;

      gamma = (alpha+beta)/2.0;
      if (gamma > beta)
        gamma = beta;

      if (alpha > fabs (res))
        value = alpha;
      else
        value = fabs(res);
      lambda = pow(value,3+alpha/2+alpha*alpha);
      if (0.25+alpha>0.5)
        value = 0.25+alpha;
      else
        value = 0.5;
      lambda /= pow(gamma, value);
      if (lambda >= 1)
      {
        rho = 1;
      }

      value = (1-lambda)/(1+lambda);
      kappa = pow(fabs(2-lambda),value);
      kappa -= 1;

      value = pow(gamma,2-alpha*alpha);
      omega = alpha*alpha * value/deltaK;

      if ((alpha < 1) && (lambda < 1))
      {
        rho = pow(omega*sigma,kappa);
      }
      sigma *= rho;
      sigma *= res*res/norm_der_u2;
      break;

    case AS97:                                    // Almeida, Silva (1997)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      normz= sqrt(z1*z1+z2*z2+z3*z3);
      if (normz==0)
      {
        sigma = 0;
        break;
      }
      if (time_dependent_problem)                 //OutPut("T");
      {
        normz *= theta1*time_step;
      }
      value = b1*u_x + b2*u_y + b3*u_z;

      if (res==0)
        value = 1;
      else
        value /= res;

      if (value < 1)
        value = 1;
      sigma = deltaK*(sqrt(norm_b2)/normz-value);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      break;

    case C93:                                     // Codina (1993)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      value = b1*u_x + b2*u_y + b3*u_z;
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      b3_orth = value * u_z/ norm_der_u2;
      norm_b_orth = sqrt(b1_orth*b1_orth+b2_orth*b2_orth+b3_orth*b3_orth);
      if (norm_b_orth == 0)
      {
        sigma = 0;
        break;
      }
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(hK * norm_b_orth);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * fabs(res)/(2*sqrt(norm_der_u2));
      break;

    case KLR02_1:                                 // Knopp, Lube, Rapin (2002)
      //case KLR02_3:                      // same as KLR02_1 with SOLD_S = 0, new version see below
#ifdef __2D__
      i = 5;
#endif
#ifdef __3D__
      i = 4;
#endif
      value = TDatabase::ParamDB->SOLD_S + param[i];
      if (value == 0)
      {
        sigma = 0;
        break;
      }
      alpha = param[i+1] / value;
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(alpha*hK);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * alpha / 2;
      break;

    case KLR02_3:                                 //pointwise evaluation of residual and norm of u_h
      if (sqrt(norm_der_u2)<epsilon)
      {
        sigma = 0;
        break;
      }
      sigma =  TDatabase::ParamDB->SOLD_CONST * hK *fabs(res)/(2*sqrt(norm_der_u2)) - eps;
      if (sigma<0)
        sigma = 0;
      break;

    case KLR02_4:                                 //pointwise evaluation of residual and norm of u_h with additive constant
      // THIS IS NOT KLR02_4 FROM John/Knobloch 2007 !!!
      if (sqrt(norm_der_u2)<epsilon)
      {
        sigma = 0;
        break;
      }
      sigma =  TDatabase::ParamDB->SOLD_CONST * hK *fabs(res)/(2*(TDatabase::ParamDB->SOLD_S +
        sqrt(norm_der_u2))) - eps;
      if (sigma<0)
        sigma = 0;
      break;

    case KLR02_2:                                 // similar to KLR02_3
    case CS99:
      // norm of convection
      value = sqrt(norm_der_u2);
      if (value == 0)
      {
        sigma = 0;
        break;
      }
      // Q = res/|b|
      alpha = fabs(res) / value;
      // C - (2 eps)/(h Q)
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(alpha*hK);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * alpha / 2;
      break;

    case J90:                                     // Johnson (1990)
      alpha = TDatabase::ParamDB->SOLD_CONST;
      kappa = TDatabase::ParamDB->SOLD_POWER;
      sigma = alpha * pow(hK,kappa)*fabs(res)-eps;
      if (sigma < 0)
        sigma = 0;
      break;

    case BE02_1:                                  // Burman, Ern 2002
      alpha = Pi/6;
      res = res*tanh(res/2.0);
      value = sqrt(norm_b2)*sqrt(norm_der_u2)+fabs(res);
      if (value > 0)
        sigma = deltaK * norm_b2*fabs(res)/value;
      else
      {
        sigma = 0;
        break;
      }
      if (norm_b2>0)
      {
        z1 = (1-b1*b1/norm_b2)*u_x - b1*(b2*u_y+ b3*u_z)/norm_b2;
        z2 = -b2*(b1*u_x + b3*u_z)/norm_b2 + (1-b2*b2/norm_b2)*u_y;
        z3 = -b3*(b1*u_x + b2*u_y)/norm_b2 + (1-b3*b3/norm_b2)*u_z;
      }
      else
      {
        sigma = 0;
        break;
      }
      value = fabs(res)+tan(alpha)*sqrt(norm_b2)*sqrt(z1*z1+z2*z2+z3*z3);
      if (value > 0)
        sigma *= (sqrt(norm_b2)*sqrt(norm_der_u2) + value)/value;
      else
      {
        sigma = 0;
      }
      break;

    case BE02_2:                                  // modified Burman, Ern 2002
      value = sqrt(norm_b2)*sqrt(norm_der_u2)+fabs(res);
      if (value > 0)
        sigma = deltaK * norm_b2*fabs(res)/value;
      else
        sigma = 0;
      break;

    case BE02_3:                                  // Burman, Ern 2002, formula (29)
      alpha = Pi/6;
      res = res*tanh(res/2.0);
      if (norm_b2>0)
      {
        z1 = (1-b1*b1/norm_b2)*u_x - b1*(b2*u_y+ b3*u_z)/norm_b2;
        z2 = -b2*(b1*u_x + b3*u_z)/norm_b2 + (1-b2*b2/norm_b2)*u_y;
        z3 = -b3*(b1*u_x + b2*u_y)/norm_b2 + (1-b3*b3/norm_b2)*u_z;
      }
      else
      {
        sigma = 0;
        break;
      }
      value = sqrt(res*res + tan(alpha)*tan(alpha)*norm_b2 * (z1*z1+z2*z2+z3*z3));
      if (value > 0)
        sigma = deltaK * norm_b2*fabs(res)/value;
      else
        sigma = 0;
      break;

    case Y_Z_BETA:                                // Tezduyar
      y = fabs(TDatabase::ParamDB->SOLD_U0);
      beta = TDatabase::ParamDB->SOLD_POWER;
      if (y==0)
      {
        sigma = 0;
        break;
      }
      if (norm_der_u2==0)
      {
        sigma = 0;
        break;
      }
      z = fabs(res)/y;
      norm_der_u2 /= (y*y);
      norm_der_u2 = pow(norm_der_u2,beta/2.0-1);
#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction_without_storing(hK, u_x, u_y)/2.0;
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction_without_storing(hK, u_x, u_y, u_z)/2.0;
#endif
      sigma = TDatabase::ParamDB->SOLD_CONST * z * norm_der_u2 * pow(hK_project,beta);
      //OutPut(sigma << " ");
      break;

    case JSW87_1:                                 // Johnson,Schatz,Wahlbin (1987) (linear)
#ifdef __2D__
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2);
#endif
#ifdef __3D__
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2, b3);
#endif
      sigma = sqrt(norm_b2)*hK_project*sqrt(hK_project) - eps;
      sigma /= theta1*time_step;
      if (sigma < 0)
        sigma = 0;
      break;

    case BH04:                                    // Burman, Hansbo 2004, edge stabilization
    case BE05_1:                                  // Burman, Hansbo 2004, edge stabilization
    case BE05_2:                                  // Burman, Hansbo 2004, edge stabilization
    case LP96:                                    // Layton, Polman 1996
    case MH_Kno06:                                // improved Mizukami-Hughes, by Knobloch 2006
    case FEM_TVD:                                 // algebraic flux correction
      sigma = 0;
      break;

    case GENERAL_SOLD:
      // TDatabase::ParamDB->SOLD_CONST is the eta from the SOLD2-paper
      if (norm_der_u2>0)
        sigma = TDatabase::ParamDB->SOLD_CONST*hK*fabs(res)/(2*sqrt(norm_der_u2));
      else
        sigma = 0.0;
      break;

    case 100:
      // take value from piecewise constant field
      // only on finest level available
      i = TDatabase::ParamDB->INTERNAL_LOCAL_DOF;
      N = TDatabase::ParamDB->INTERNAL_ARRAY_LENGTH;
      sigma =  TDatabase::ParamDB->INTERNAL_P1_Array[N+i];
      return(sigma);
      break;

      // for lower levels
    case 101:
      return(0.0);
      break;

    default :
      OutPut("SOLD type " << sold_parameter_type << " not available" << endl);
      exit(4711);
  }
#ifdef __2D__
  i = 5;
#endif
#ifdef __3D__
  i = 4;
#endif
  // scaling of sigma accordingly to Knopp, Lube, Rapin (2002)
  if (TDatabase::ParamDB->SOLD_PARAMETER_SCALING)
  {
    value = TDatabase::ParamDB->SOLD_S + param[i];
    if (value == 0)
      sigma = 0;
    else
    {
      alpha = param[i+1] / value;
      sigma *= alpha*alpha;
    }
  }
  else
    sigma *= TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR;

  return (sigma);
}

/*************************************************************************/
// estimate of coercivity constant 
// used eg for residual based estimator of Verf"uhrt 2005
/*************************************************************************/
double EstimateCoercivityConstant(TCollection *Coll,
#ifdef __2D__
                                  CoeffFct2D *Coeffs
#else // 3D
                                  CoeffFct3D *Coeffs
#endif
                                  )
{
  int i, j, N_Cells, N_V;
  double coerc = 4711.0, x, y, xs, ys, *coeff;
  #ifdef __3D__
  double z, zs;
  #endif
  TBaseCell *cell;
  TVertex *vertex;

  coeff = new double[13];
  // index of the reaction term in coeff array
  #ifdef __2D__
  int reac_index = 3;
  #else
  int reac_index = 4;
  #endif
  
    
  N_Cells = Coll->GetN_Cells();                   // number of mesh cells
  // loop over the cells
  for (i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    // loop over the vertices
    xs = ys = 0;
    for (j=0;j<N_V;j++)
    {
      vertex = cell->GetVertex(j);
      #ifdef __2D__
      vertex->GetCoords(x, y);
      Coeffs(1, &x, &y, NULL, &coeff);
      #else
      vertex->GetCoords(x, y, z);
      Coeffs(1, &x, &y, &z, NULL, &coeff);
      #endif
      // assuming that convection is divergence-free
      if (coeff[reac_index] < coerc)
        coerc = coeff[reac_index];
      xs += x;
      ys += y;
      #ifdef __3D__
      zs += z;
      #endif
    }
    xs /= N_V;
    ys /= N_V;
    #ifdef __3D__
    zs /= N_V;
    #endif
    
    #ifdef __2D__
    Coeffs(1, &xs, &ys, NULL, &coeff);
    #else
    Coeffs(1, &xs, &ys, &zs, NULL, &coeff);
    #endif
    // assuming that convection is divergence-free
    if (coeff[reac_index] < coerc)
      coerc = coeff[reac_index];
  }

  delete coeff;
  if(TDatabase::ParamDB->SC_VERBOSE>1)
    OutPut("coercivity constant (assuming div-free convection): " << coerc << endl);
  return(coerc);
}


void SetSoldParameters(int i)
{
  const int max = 16;
  TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 1.0;
  TDatabase::ParamDB->SOLD_TYPE = 1;
  
  if ((i>max-1) && (i<2*max))
  {
    i -= max;
    TDatabase::ParamDB->SOLD_TYPE = 2;
  }

  switch(i)
  {
    case 0:
      TDatabase::ParamDB->SOLD_TYPE = 0;
      break;
    case 1:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
      break;
    case 2:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 1;
      break;
    case 3:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 3;
      break;
    case 4:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 4;
      break;
    case 5:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 5;
      break;
    case 6:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 6;
      break;
    case 7:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.25;
      break;
    case 8:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      break;
    case 9:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.75;
      break;
    case 10:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.25;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 11:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 12:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.75;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 13:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      TDatabase::ParamDB->SOLD_S = 1.0;
      break;
    case 14:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 9;
      TDatabase::ParamDB->SOLD_CONST = 1;
      TDatabase::ParamDB->SOLD_POWER = 1.0;
      break;
    case 15:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 9;
      TDatabase::ParamDB->SOLD_CONST = 1;
      TDatabase::ParamDB->SOLD_POWER = 2.0;
      break;
    default:
      OutPut("Wrong input in SetSoldParameters" << endl);
      exit(4711);
  }
}


#ifdef __2D__
void EdgeStabilization(TFESpace2D *fespace, TFEFunction2D *u, 
                       CoeffFct2D *Coeffs, double *rhs, int time_dependent,
                       double *time_step, TFEFunction2D *old_u)
{
  int i, j, k, ii, N_Cells, *ColInd, *RowPtr, *GlobalNumbers, *BeginIndex;
  int ActiveBound, *DOF, N_Edges, boundedge, locdof;
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double val[3], val_neigh[3], h, norm_t, x[3], y[3], oldval[3];
  double x0, x1, y0, y1, xs, ys, t1, t2, *coeff, jump, fac0, fac1, fac2;
  double phi0_x, phi0_y, phi1_x, phi1_y, phi2_x, phi2_y, n1, n2, maxjump;
  double sx, sy, tmp, meas, area, rho = 2.0;
  TBaseCell *cell, *neigh;
  TCollection *coll;
  FE2D CurrentElement;
  TJoint *joint;
  TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  const int *TmpEdVer;

  coeff = new double[6];

  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  // get start of dirichlet nodes in dof array
  ActiveBound = fespace->GetActiveBound();
  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // assign a numbering to the cells
  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=coll->GetCell(i);
    cell->SetClipBoard(i);
  }                                               // endfor i

  // loop over all cells for computing the edge stabilization
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    h = cell->GetDiameter();
    meas = cell->GetMeasure();
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];

    // local dofs are arranged as follows
    // local dof 0 on vertex 0 opposite to edge 1
    // local dof 1 on vertex 1 opposite to edge 2
    // local dof 2 on vertex 2 opposite to edge 0

    CurrentElement = fespace->GetFE2D(i, cell);
    if (CurrentElement!=C_P1_2D_T_A)
    {
      if (sold_parameter_type!=BE05_2)
      {
        OutPut("Edge stabilization for element " << CurrentElement <<
          " not implemented !!!"<< endl);
        exit(4711);
      }
      if ((CurrentElement!=C_Q1_2D_Q_A)&&(CurrentElement!=C_Q1_2D_Q_M))
      {
        OutPut("Edge stabilization for element " << CurrentElement <<
          " not implemented !!!"<< endl);
        exit(4711);
      }
    }
    // # of edges
    N_Edges = cell->GetN_Edges();

    sx = sy = 0;
    // compute derivatives for basis functions
    for (j=0;j<N_Edges; j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
      sx += x[j];
      sy += y[j];
      //OutPut(x[j] << " " << y[j] << " ");
      u->FindGradientLocal(cell, i, x[j], y[j], val);
      //OutPut("u"<<j << " " << val[0]<<endl);
    }
    sx /= N_Edges;
    sy /= N_Edges;
    //OutPut(endl);
    // compute twice area of triangle
    if (N_Edges==3)
    {
      area = 2*meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }
    else
    {
      // Q1
      OutPut("Implementation not complete"<<endl);
      exit(4711);
      area = meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }

    /* OutPut("0 " << phi0_x << " " << phi0_y << endl);
     OutPut("1 " << phi1_x << " " << phi1_y << endl);
     OutPut("2 " << phi2_x << " " << phi2_y << endl);
    */
    // get refinement descriptor
    refdesc=cell->GetRefDesc();
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);

    // compute gradient of current solution (constant)
    u->FindGradientLocal(cell, i, sx, sy, val);

    // compute maximal normal jump
    if (sold_parameter_type==BH04)
    {
      maxjump = 0;
      for(j=0;j<N_Edges;j++)                      // loop over all edges of cell
      {
        joint=cell->GetJoint(j);
        ver0=cell->GetVertex(TmpEdVer[2*j]);      // get vertices of face j
        ver1=cell->GetVertex(TmpEdVer[2*j+1]);
        x0 = ver0->GetX();                        // coordinates of face j
        y0 = ver0->GetY();
        x1 = ver1->GetX();
        y1 = ver1->GetY();

        // compute tangential
        t1 = x1 - x0;
        t2 = y1 - y0;
        norm_t = sqrt(t1*t1+t2*t2);
        t1 /= norm_t;
        t2 /= norm_t;
        // compute normal
        n1 = -t2;
        n2 = t1;
        // compute solution (including derivative) in midpoint of tangential
        // from point of view of this mesh cell
        xs = (x1+x0)/2;
        ys = (y1+y0)/2;

        // compute solution (including derivative) in midpoint of tangential
        // from point of view of neighbour mesh cell
        // NO ADAPTIVE MESHES ALLOWED
        neigh=joint->GetNeighbour(cell);          // neighbour cell
        if (neigh!=NULL)
        {
          ii =  neigh->GetClipBoard();
          u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        }
        else
        {
          // boundary edge
          // continue;
          val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
        }
        jump = (n1 * val[1] + n2 * val[2]) - (n1 * val_neigh[1] + n2 * val_neigh[2]);
        jump = fabs(jump);
        if (jump > maxjump)
          maxjump = jump;
      }
      //OutPut(" " << maxjump);
    }

    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      ver0=cell->GetVertex(TmpEdVer[2*j]);        // get vertices of face j
      ver1=cell->GetVertex(TmpEdVer[2*j+1]);
      x0 = ver0->GetX();                          // coordinates of face j
      y0 = ver0->GetY();
      x1 = ver1->GetX();
      y1 = ver1->GetY();
      //OutPut("ed " << j << " " << x0 << " " << y0 << " ; " << x1 << " " <<y1
      //      << endl);
      // compute tangential
      t1 = x1 - x0;
      t2 = y1 - y0;
      norm_t = sqrt(t1*t1+t2*t2);
      t1 /= norm_t;
      t2 /= norm_t;
      //OutPut(t1 << " " << t2 << " " << t1*t1+t2*t2 << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of this mesh cell
      xs = (x1+x0)/2;
      ys = (y1+y0)/2;
      //u->FindGradientLocal(cell, i, xs, ys, val);
      //OutPut("grad_i " << val[1] << " " << val[2] << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of neighbour mesh cell
      // NO ADAPTIVE MESHES ALLOWED
      neigh=joint->GetNeighbour(cell);            // neighbour cell
      if (neigh!=NULL)
      {
        ii =  neigh->GetClipBoard();
        //OutPut("ii " << ii << endl);
        u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        boundedge = 0;
      }
      else
      {
        // boundary edge
        val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
        boundedge = 1;
      }
      //OutPut("grad_ii " << val_neigh[1] << " " << val_neigh[2] << endl);

      // compute factor which defines the sign
      fac0 = t1 * val[1] + t2 * val[2];
      //OutPut("fac 0 " << fac0 << " ");
      /*if (fac0 > 1e-8)
          fac0 = 1;
      else
      {
          if (fac0<-1e-8)
        fac0 = -1;
          else
        fac0 = 0;
        }*/
      //OutPut("fac0 " << fac0<<endl);
      tmp = tanh(fac0/1.0)/fac0;
      fac0 = tanh(fac0/1.0);
      //OutPut(fac0 << endl);
      // compute nonlinear factor depending on u (Psi_K(u))
      switch (sold_parameter_type)
      {
        case BH04:
          // compute coefficients in (xs,ys)
          Coeffs(1, &xs, &ys, NULL, &coeff);
          fac1 = coeff[0] * TDatabase::ParamDB->SOLD_CONST;
          if (time_dependent)
            fac1 *= time_step[0];
          fac1 += h * TDatabase::ParamDB->SOLD_S;
          fac1 *= h *maxjump/meas;
          break;
        case BE05_1:
          // compute coefficients in (xs,ys)
          Coeffs(1, &xs, &ys, NULL, &coeff);
          // norm of convection
          fac1 = sqrt(coeff[1]*coeff[1]+coeff[2]*coeff[2]);
          if (time_dependent)
            fac1 *= time_step[0];
          //OutPut("conv " << fac1 << " " << h << endl);
          fac1 *= h*h;
          // reaction term
          fac2 = rho * fabs(coeff[3]);
          if (time_dependent)
          {
            fac2 *= time_step[0];
            fac2 += 1.0;
            //OutPut(fac2 << " ");
          }
          fac2 *=h*h*h;
          fac1 += fac2;
          // jump of gradient
          jump = (val[1]-val_neigh[1]) * (val[1]-val_neigh[1]);
          jump += (val[2]-val_neigh[2]) * (val[2]-val_neigh[2]);
          fac1 *= sqrt(jump);
          //OutPut("jump " << sqrt(jump) << endl);
          fac1 *= TDatabase::ParamDB->SOLD_CONST/meas;
          break;

        case BE05_2:
          // this case does not depend from formerly computed values
          // Simpson rule
          // compute coefficients in (x0,y0)
          Coeffs(1, &x0, &y0, NULL, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((x0<0.25)||(x0>0.75)||(y0<0.25)||(y0>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 = coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4];
          }
          else
          {
            old_u->FindGradientLocal(cell, i, x0, y0, oldval);
            fac1 = val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4];
          }
          Coeffs(1, &xs, &ys, NULL, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((sx<0.25)||(sx>0.75)||(sy<0.25)||(sy>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 += 4*(coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4]);
          }
          else
          {
            old_u->FindGradientLocal(cell, i, xs, ys, oldval);
            fac1 += 4 * (val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4]);
          }
          // compute coefficients in (x1,y1)
          Coeffs(1, &x1, &y1, NULL, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((x1<0.25)||(x1>0.75)||(y1<0.25)||(y1>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 += coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4];
          }
          else
          {
            old_u->FindGradientLocal(cell, i, x1, y1, oldval);
            fac1 += val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4];
          }

          fac1 = TDatabase::ParamDB->SOLD_CONST * fabs(fac1)/6.0;
          break;
        default :
          OutPut("Edge stabilization " << sold_parameter_type
            << " not implemented !!!" <<endl);
          exit(4711);
      }
      // norm_t is the length of the edge
      //OutPut(tmp*fac1/norm_t <<endl);
      fac1 = fac0*meas*fac1*norm_t;
      //OutPut("a_i " << fac1/norm_t << endl );
      // update the rhs
      switch(j)
      {
        // edge zero, active dof are local 0 and 1
        case 0:
          // local dof 0
          locdof = DOF[0];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi0_x + t2*phi0_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 1
          locdof = DOF[1];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi1_x + t2*phi1_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
          // edge one, active dof are local 1 and 2
        case 1:
          // local dof 1
          locdof = DOF[1];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi1_x + t2*phi1_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 2
          locdof = DOF[2];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi2_x + t2*phi2_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
          // edge two, active dof are local 0 and 2
        case 2:
          // local dof 0
          locdof = DOF[0];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi0_x + t2*phi0_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 2
          locdof = DOF[2];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi2_x + t2*phi2_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
      }
    }
  }                                               // loop over cells
  /*
    // loop over all cells for computing the edge stabilization
    for(i=0;i<N_Cells;i++)
    {
      // next cell
      cell = coll->GetCell(i);
      // pointer to global indices of dof connected with this cell
      DOF = GlobalNumbers + BeginIndex[i];

      // local dofs are arranged as follows
      // local dof 0 on vertex 0 opposite to edge 1
  // local dof 1 on vertex 1 opposite to edge 2
  // local dof 2 on vertex 2 opposite to edge 0

  // # of edges
  N_Edges = cell->GetN_Edges();

  // compute derivatives for basis functions
  for (j=0;j<N_Edges; j++)
  {
  x[j] = cell->GetVertex(j)->GetX();
  y[j] = cell->GetVertex(j)->GetY();
  OutPut(x[j] << " " << y[j] << " ori " << rhsori[DOF[j]] << " update " << rhs[DOF[j]] <<
  " diff " <<  rhsori[DOF[j]] - rhs[DOF[j]] << endl);
  }
  }
  */
  delete coeff;

}

#endif // 2D


double ComputeAlpha(double hK) // copy from TCD2D.C
{
  double alpha;
  
  alpha = TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT*
     pow(hK, TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_POWER);
  return(alpha);

  // this is just for the response to the referee and the special example 
  // in [JKL05]
  double b, eps, Pe, t;

  b = sqrt(5.0);
  eps = 1/TDatabase::ParamDB->RE_NR;
  Pe = b*hK/(2*eps);
  t = 1/tanh(Pe) - 1/Pe;
  alpha = t*hK/(2*b);
  return(alpha);
}
