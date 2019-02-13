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
//  TCD2D.C 
// ======================================================================

#include <TCD2D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <FEFunction2D.h>
#include <ConvDiff.h>
#include <math.h>
#include <stdlib.h>

/*
double ComputeAlpha(double hK)
{
  double alpha;
  
  alpha = TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT*
     pow(hK,TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_POWER);
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
*/







// ==

// ======================================================================
// assemble matrices B,C,M from 
// ( A B )
// ( C M ) 
// ======================================================================
/*
void MatricesAssemble_VMM(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixMcoarse, **MatrixB, **MatrixC, val, *MatrixRowMcoarse;
  double *MatrixRowB, *MatrixRowC;
  double ansatz10, ansatz01 ;
  double test10, test01, h;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_, N_coarse;

  MatrixMcoarse = LocMatrices[0];
  MatrixB = LocMatrices[1];
  MatrixC = LocMatrices[2];

  N_ = N_BaseFuncts[0];
  N_coarse = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // phi_x
  Orig1 = OrigValues[1]; // phi_y
  Orig2 = OrigValues[2]; // psi_x
  Orig3 = OrigValues[3]; // psi_y

  h = ComputeAlpha(hK);

  for(i=0;i<N_coarse;i++)
  {
    MatrixRowMcoarse = MatrixMcoarse[i];
    MatrixRowC = MatrixC[i];

    test10 = Orig2[i];
    test01 = Orig3[i];

    for(j=0;j<N_coarse;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      val = ansatz10*test10 + ansatz01*test01;
      MatrixRowMcoarse[j] += Mult * val;
    } // endfor j
    
    for (j=0;j<N_;j++)
    {
       ansatz10 = Orig0[j];
       ansatz01 = Orig1[j];
       val = ansatz10*test10 + ansatz01*test01;
       MatrixRowC[j] -= Mult * val;       
    }
  } // endfor i
  
  for (i=0;i<N_;i++)
  {
     MatrixRowB = MatrixB[i];
     
     test10 = Orig0[i];
     test01 = Orig1[i];
     for(j=0;j<N_coarse;j++)
     {
        ansatz10 = Orig2[j];
        ansatz01 = Orig3[j];
        val = ansatz10*test10+ansatz01*test01;
        MatrixRowB[j] -=  Mult*val*h;
     }
  }
}
*/
// ======================================================================
// assemble matrices B1, B2, C1, C2, M from 
// ( A     B1 B2 )
// ( C1 c2   M ) 
// ======================================================================

/*
void MatricesAssemble_VMM_KL02(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixMcoarse, **MatrixB1, **MatrixC1, val, *MatrixRowMcoarse;
  double **MatrixB2, **MatrixC2;
  double *MatrixRowB1, *MatrixRowC1, *MatrixRowB2, *MatrixRowC2;
  double ansatz10, ansatz01 ;
  double test10, test01;
  double *Orig0, *Orig1, h;
  int i,j, N_, N_coarse;
 
  MatrixMcoarse = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixC1 = LocMatrices[3];
  MatrixC2 = LocMatrices[4];

  N_ = N_BaseFuncts[0];
  N_coarse = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // phi_x
  Orig1 = OrigValues[1]; // phi_y

  h = ComputeAlpha(hK);

  for(i=0;i<N_coarse;i++)
  {
    MatrixRowMcoarse = MatrixMcoarse[i];
    MatrixRowC1 = MatrixC1[i];
    MatrixRowC2 = MatrixC2[i];

    for(j=0;j<N_coarse;j++)
    {
      MatrixRowMcoarse[j] += Mult;
    } // endfor j
    
    for (j=0;j<N_;j++)
    {
       ansatz10 = Orig0[j];
       ansatz01 = Orig1[j];
       MatrixRowC1[j] -= Mult * ansatz10;       
       MatrixRowC2[j] -= Mult * ansatz01;       
    }
  } // endfor i
  
  for (i=0;i<N_;i++)
  {
     MatrixRowB1 = MatrixB1[i];
     MatrixRowB2 = MatrixB2[i];
     
     test10 = Orig0[i];
     test01 = Orig1[i];
     for(j=0;j<N_coarse;j++)
     {
        MatrixRowB1[j] -=  Mult*test10*h;
        MatrixRowB2[j] -=  Mult*test01*h;       
     } 
  }
}
*/

// ======================================================================
// MATRICES FOR REACTION PART OF BULK PRECIPITATION
// ======================================================================

// ======================================================================
// assemble mass matriX 
// ======================================================================
void MatrixMAssemble_Bulk(double Mult, double *coeff, double *param,
                          double hK, 
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00;
  double test00;
  double *Orig0;
  int i,j, N_;
 
  Matrix = LocMatrices[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];
 
    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];

      if (!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	  MatrixRow[j] += Mult * ansatz00*test00;
      else
      {
	  // add to diagonal
          //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      MatrixRow[i] += Mult * ansatz00*test00;
      }
    } // endfor j
  } // endfor i
}

// ======================================================================
// assemble matrix and rhs for upwind discretization
// ======================================================================
void MatricesA_Assemble_Bulk(double Mult, double *coeff, double *param,
			    double hK, 
			    double **OrigValues, int *N_BaseFuncts,
			    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRowA;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, Pe; 
  double k;

  MatrixA = LocMatrices[0];
 
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration

  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
		
      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  ansatz00 = Orig2[j];
	      
	  val =  c0*(test10*ansatz10+test01*ansatz01);
	  if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	      ||(i==j))
	  {
	      val += c3*ansatz00*test00;
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	  }
	  else
	  {
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	      //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      {
		  val += c3*ansatz00*test00;
		  // add to diagonal
		  MatrixRowA[i] += Mult * val;
	      }	     
	  }
      } // endfor j
  } // endfor i
}

void Rhs_Assemble_Bulk(double Mult, double *coeff, double *param,
			    double hK, 
			    double **OrigValues, int *N_BaseFuncts,
			    double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double test00; 
  double *Orig2;
  int i, N_;
  double c4; 

  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];
  Orig2 = OrigValues[0];
  c4 = coeff[4];

  for(i=0;i<N_;i++)
  {
      test00 = Orig2[i];		
      Rhs[i] += Mult*test00*c4;
  }
}

// ======================================================================
// assemble matrices and rhs for SUPGtest00
// matrix A contains the part of the equation which is multiplied with
//        theta_1 Delta t
// matrix K contains the part of the SUPG term which comes from the 
//        temporal discretization of the concentration
// ======================================================================
void MatricesA_Assemble_SUPG_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA,  **MatrixK, val, *MatrixRowA, *MatrixRowK;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_, sold_type=TDatabase::ParamDB->SOLD_TYPE;
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double c0, c1, c2, c3, c4, Pe; 
  double tau, bgradv, bb, k, c00, c11, c22, c33, c44, param_sold[5];
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1, theta2;
  double sigma, conc_old, reaction_old, conc_old_x, conc_old_y;

  MatrixA = LocMatrices[0];
  MatrixK= LocMatrices[1];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration
  c4 = coeff[4];

  c00 = c0;
  c11 = c1;
  c22 = c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0/(time_step * theta1) + c3;
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);
  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      MatrixRowK = MatrixK[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
      
      // supg test term
      bgradv = tau*(c1*test10+c2*test01);
      
      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  ansatz00 = Orig2[j];
	  
	  val =  c0*(test10*ansatz10+test01*ansatz01);
	  val += (c1*ansatz10+c2*ansatz01)* (test00+bgradv);
	  if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	      ||(i==j))
	  {
	      val += c3*ansatz00*(test00+bgradv);
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	  }
	  else
	  {
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	      //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      {
		  val = c3*ansatz00*(test00+bgradv);
		  // add to diagonal
		  MatrixRowA[i] += Mult * val;
	      }	     
	  }
	  // sdfem part which comes from time derivative 
	  MatrixRowK[j] += Mult * ansatz00*bgradv;
      } // endfor j
  } // endfor i
  
  // only SUPG term
  if  (!sold_type)
      return;

  MatrixA = LocMatrices[2];  // SOLD matrix
  theta2 = TDatabase::TimeDB->THETA2;

  // rhs of time-dependent equation for c_A, c_B
  // current velocity is used (as approximation)
  // strong diffusion term is neglected
  // the right hand side is treated in a simplified way
  conc_old = param[6];
  reaction_old = param[7];  
  conc_old_x = param[8];    
  conc_old_y = param[9];
  // rhs of time-dependent equation for c_C
  c44 = conc_old 
      - theta2 * time_step*(c1*conc_old_x + c2 *conc_old_y + reaction_old * conc_old)+ c4;

  //OutPut(conc_old << " " << reaction_old  << " "  << c44 << endl);


  if (sold_parameter_type > BE02_2)
  {
      OutPut("SOLD_PARAMETER_TYPE "<< sold_parameter_type << " not yet implemented!"<<endl);
      exit(4711);
  }
  
  time_step *= theta1;
  c00 = time_step * c0;
  c11 = time_step * c1;
  c22 = time_step * c2;
  c33 = 1.0 + time_step * c3;
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  tau *=  time_step * theta1;

  param_sold[0] = param[10];  // concentration c
  param_sold[1] = param[11];  // c_x
  param_sold[2] = param[12];  // c_y
  param_sold[3] = 0;
  param_sold[4] = 0;
  sigma = Compute_SOLD_sigma(hK, c00, c11, c22, c33, c44, bb, tau, param_sold, 0, 0,1);
  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];

      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  val = 0;
	  if  (sold_type==1)
	  {
	      // isotropic additional diffusion
	      if (bb >0)
		  val = sigma*(test10*ansatz10+test01*ansatz01);
	  }
	  else
	  {
	      // additional diffusion orthogonal to the streamlines
	      if (bb >0)
		  val = sigma * (-c22*ansatz10+c11*ansatz01)*(-c22*test10+c11*test01)/(c11*c11+c22*c22);
	  }
	  MatrixRowA[j] += Mult * val/time_step;
      } // endfor j
  } // endfor i
}
void Rhs_Assemble_SUPG_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, val, *MatrixRowA, *MatrixRowK;
  double *Rhs;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_, sold_type=TDatabase::ParamDB->SOLD_TYPE;
  double c0, c1, c2, c3, c4, Pe; 
  double tau, bgradv, bb, k, c00, c11, c22, c33;
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double sigma, conc_old;
 
  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration
  c4 = coeff[4];

  c00 = c0;
  c11 = c1;
  c22 = c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0/(time_step * theta1) + c3;
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  for(i=0;i<N_;i++) 
  {
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
      
      // supg test term
      bgradv = tau*(c1*test10+c2*test01);
      Rhs[i] +=  Mult*(test00 + bgradv)*c4;      
  }
}

// ======================================================================
// assemble matrices and rhs for Galerkin (FCT-FEM)
// matrix A contains the part of the equation which is multiplied with
//        theta_1 Delta t
// ======================================================================
void MatricesA_Assemble_Galerkin_Bulk(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRowA;
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4;

  MatrixA = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // other concentration
  c4 = coeff[4];

  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      test00 = Orig2[i];
      
      for(j=0;j<N_;j++)
      {
	  ansatz10 = Orig0[j];
	  ansatz01 = Orig1[j];
	  ansatz00 = Orig2[j];
	  
	  val =  c0*(test10*ansatz10+test01*ansatz01);
	  val += (c1*ansatz10+c2*ansatz01)*  test00;
	  if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	      ||(i==j))
	  {
	      val += c3*ansatz00*test00;
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	  }
	  else
	  {
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	      //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	      {
		  val = c3*ansatz00*test00;
		  // add to diagonal
		  MatrixRowA[i] += Mult * val;
	      }	     
	  }
      } // endfor j
  } // endfor i
}





// this routine gives the fe value back
void TimeCDParamsVeloField(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // value of first input fe function
  out[1] = in[3]; // value of second input fe function
  out[2] = in[0]; // x value for axial symmetric  
  // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

// this routine gives the fe value back
void TimeCDParamsVeloField_ALE(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // value of first input fe function
  out[1] = in[3]; // value of second input fe function
  out[2] = in[0]; // x value for axial symmetric  
  out[3] = in[4] + in[5];  // \nabla\cdot w (divergence of the mesh velo, conservatiove ALE)
//   cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

// this routine gives the fe value back
void TimeCDParamsSOLD(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  // cout << in[0] << " " << in[1] << endl;
  out[0] = in[2]; // current solution u
  out[1] = in[3]; // current solution u_x
  out[2] = in[4]; // current solution u_y
  out[3] = in[5]; // current solution u_xx
  out[4] = in[6]; // current solution u_yy
  out[5] = in[7]; // old solution u
  out[6] = in[8]; // old solution u_x
  out[7] = in[9]; // old solution u_y
  out[8] = in[10]; // old solution u_xx
  out[9] = in[11]; // old solution u_yy
  //cout << out[0] << " " << out[1] <<  " " << out[2] << endl; 
}

void TimeCDParamsBulk(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  if (in[2]>0)
     out[0] = in[2]; // value of first input fe function (other concentration)
  else
     out[0] = 0.0;
  out[1] = in[3]; // value of second input fe function (u1)
  out[2] = in[4]; // value of third input fe function (u2)
  // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

void TimeCDParamsBulk_SOLD(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  if (in[2]>0)
     out[0] = in[2]; // value of first input fe function (other concentration)
  else
     out[0] = 0.0;
  out[1] = in[3]; // value of second input fe function (u1)
  out[2] = in[4]; // value of third input fe function (u2)
  // this have to be the same indices as in TimeCDParamsBulk_SOLD_Cc
  out[6] = in[5]; // value of old concentration
  out[7] = in[6]; // value of old reaction (other concentration)
  out[8] = in[7]; // x-deriv of old concentration
  out[9] = in[8]; // y_deriv old concentration
  out[10] = in[9];  // concentration
  out[11] = in[10]; // x-deriv of concentration
  out[12] = in[11]; // y-deriv of concentration
  // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

void TimeCDParamsBulk_Cc(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[4]; // u1
  out[1] = in[5]; // u2
  out[2] = in[2]; // C_a
  out[3] = in[3]; // C_b
  out[4] = in[6]; // C_c_old
  out[5] = in[7]; // integral_value
}

void TimeCDParamsBulk_SOLD_Cc(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[4]; // u1
  out[1] = in[5]; // u2
  out[2] = in[2]; // C_a
  out[3] = in[3]; // C_b
  out[4] = in[6]; // C_c_old
  out[5] = in[7]; // integral_value
  out[6] = in[8]; // value of old concentration
  out[7] = 0; // explicit treatment
  out[8] = in[11]; // x-deriv of old concentration
  out[9] = in[12]; // y_deriv of old concentration
  out[10] = in[4]; // x-deriv of concentration
  out[11] = in[13]; // x-deriv of concentration
  out[12] = in[14]; // y_deriv of concentration
}

void TimeCDParamsBulk_mom(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // c_C
  out[3] = in[5]; // mom_{k-1}
}

void TimeCDParamsBulk_SOLD_mom(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // C_c
  out[3] = in[5]; // C_c_old
}


void TimeCDParamsUrea(double *in, double *out)
{
  // in[0], in[1]  are coordinates
  out[0] = in[2]; // value of second input fe function (u1)
  out[1] = in[3]; // value of third input fe function (u2)
}

void TimeCDParamsUrea_conc(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // concentration
  out[3] = in[5]; // temp
  out[4] = in[6]; // integral conce
}

void TimeCDParamsUrea_temp(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // concentration
  out[3] = in[5]; // temp
  out[4] = in[6]; // integral conce
}

void TimeCDParamsUrea_conc2(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // concentration
  out[3] = in[5]; // temp
  out[4] = in[6]; // integral conce
  out[5] = in[7]; // integral conce
}

void TimeCDParamsUrea_temp2(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
  out[2] = in[4]; // concentration
  out[3] = in[5]; // temp
  out[4] = in[6]; // integral conce
  out[5] = in[7]; // integral conce
}

void TimeCDParamsUrea_conc_mat(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[2]; // u1
  out[1] = in[3]; // u2
}



void JumpTermsForIMEX_P1(TFESpace2D *fespace,
TFEFunction2D *u,
BoundCondFunct2D *BoundaryConditions,
double *sold_param)
{
  int i, j, k, ii, N_Cells, *ColInd, *RowPtr, *GlobalNumbers, *BeginIndex;
  int ActiveBound, *DOF, *DOF_n, N_Edges, boundedge, locdof, found;
  int com00, com01, com10, com11, com20, com21, com_other0, com_other1;
  int loc_vert_n, comp;
  double val[3], val_neigh[3], h, norm_t, x[3], y[3], oldval[3];
  double x_n[3], y_n[3], eps = 1e-6;
  double x0, x1, y0, y1, xs, ys, t1, t2, *coeff, jump, fac0, fac1, fac2;
  double phi0_x, phi0_y, phi1_x, phi1_y, phi2_x, phi2_y, n1, n2, maxjump;
  double phi0_n_x, phi0_n_y, phi1_n_x, phi1_n_y, phi2_n_x, phi2_n_y;
  double phi_n_other_x, phi_n_other_y, p0, p1;
  double sx, sy, tmp, meas, area, rho = 2.0, ansatz, test, area_n, meas_n;
  double integral, norm_grad_u, ave, integral_ave;
  TBaseCell *cell, *neigh;
  TCollection *coll;
  FE2D CurrentElement;
  TJoint *joint;
  TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  BoundCond BdCond;
  TBoundComp *BdComp;
  TBoundEdge *bound_edge;
  TIsoBoundEdge *isobound_edge;
  const int *TmpEdVer;

  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // assign a numbering to the cells
  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=coll->GetCell(i);
    cell->SetClipBoard(i);
  }                                               // endfor i

  // loop over all cells for computing the jump terms
  for(i=0;i<N_Cells;i++)
  {
    integral = integral_ave = 0;
    // next cell
    cell = coll->GetCell(i);
    h = cell->GetDiameter();
    //meas = cell->GetMeasure();
    // pointer to global indices of dof connected with this cell
    //DOF = GlobalNumbers + BeginIndex[i];

    // local dofs are arranged as follows
    // local dof 0 on vertex 0 opposite to edge 1
    // local dof 1 on vertex 1 opposite to edge 2
    // local dof 2 on vertex 2 opposite to edge 0

    CurrentElement = fespace->GetFE2D(i, cell);
    if (CurrentElement!=C_P1_2D_T_A)
    {
      OutPut("JumpTermsForIMEX_P1 for element " << CurrentElement <<
        " not implemented !!!"<< endl);
      exit(4711);
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
    }
    sx /= N_Edges;
    sy /= N_Edges;
    u->FindGradientLocal(cell, i, sx, sy, val);
    norm_grad_u = sqrt(val[1]*val[1]+val[2]*val[2]);

    // get refinement descriptor
    refdesc=cell->GetRefDesc();
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);

    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      ver0=cell->GetVertex(TmpEdVer[2*j]);        // get vertices of face j
      ver1=cell->GetVertex(TmpEdVer[2*j+1]);
      x0 = ver0->GetX();                          // coordinates of face j
      y0 = ver0->GetY();
      x1 = ver1->GetX();
      y1 = ver1->GetY();
      //OutPut(endl << "ed " << j << " " << x0 << " " << y0 << " ; " << x1 << " " <<y1
      //  << endl);
      // compute tangential
      t1 = x1 - x0;
      t2 = y1 - y0;
      norm_t = sqrt(t1*t1+t2*t2);
      t1 /= norm_t;
      t2 /= norm_t;
      // compute normal
      n1 = t2;
      n2 = -t1;
      //OutPut(t1 << " " << t2 << " " << t1*t1+t2*t2 << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of this mesh cell
      xs = (x1+x0)/2;
      ys = (y1+y0)/2;
      // NO ADAPTIVE MESHES ALLOWED
      neigh=joint->GetNeighbour(cell);            // neighbour cell
      if (neigh!=NULL)
      {
        ii =  neigh->GetClipBoard();
        //OutPut("ii " << ii << endl);
        //DOF_n = GlobalNumbers + BeginIndex[ii];
        u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        jump = (val_neigh[1] - val[1]) * n1 + (val_neigh[2] - val[2]) * n2;
        ave = fabs( val[1]*n1 +   val[2]*n2) +  fabs( val_neigh[1]*n1 +   val_neigh[2]*n2);
        ave /= 2;
        integral += jump * jump * norm_t;
        integral_ave += ave * ave * norm_t;
      }
    }
    /*if (norm_grad_u > 0)
      sold_param[i] = sqrt(integral)/norm_grad_u;
    else
    sold_param[i] = 0;*/
    if (integral_ave > 0)
      sold_param[i] = sqrt(integral/integral_ave);
    else
      sold_param[i] = 0;
  }                                               // loop over cells
}


void ParamsFct_HeatLine(double *in, double *out)
{
  out[0] = in[2]; 
  out[1] = in[3];
  out[2] = in[4]; 
  out[3] = in[5];
  out[4] = in[6]; 
  out[5] = in[7]; 
  out[6] = in[8]; 

// cout<< "ParamsFct_HeatLine out[0]  " << out[0]<<endl;
}
