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
// @(#)NSE2D_FixPo.C        1.3 06/27/00
//
// common declaration for all Navier-Stokes problems
// fix point iteration
//
// ======================================================================
#include <Convolution.h>
#include <Database.h>
#include <TNSE2D_Routines.h>
#include <stdlib.h>

// ======================================================================
// compute parameter for RFB stabilization
// Brezzi, Marini, Russo, CMAME 194 (2005) 127 - 148
// ======================================================================

double RFB_Parameter(double hK, double eps, double* b)
{
    double peclet, normb,tau;

    normb = sqrt(b[0]*b[0]+b[1]*b[1]);
    if (normb==0)
	return(0.0);
    peclet = normb*hK/(2*eps);
    if (peclet > 1)
        tau = hK/(3*normb);
    else
	tau = hK*hK/eps;

    tau = hK*hK;
    return(tau);
	    
}

// ======================================================================
// compute parameter for SUPG stabilization
// ======================================================================

double SUPG_Parameter(double hK, double eps, double b1, double b2, double c)
{
  double delta, val, norm_b, M, r=2.0, c_inv_2 = 24.0, beta_0 = 1.0, c_pf = 1.0;
  double beta_1 = 1.0, C0, C_h = 1.0, delta_0 = 0.1;
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  
  switch (TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case 2:
    case 12:
        if (TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE==3)
        c_inv_2 = 48.0;
  else
        c_inv_2 = 24.0;
      break;
      case 3:
      case 13:
        if (TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE==3)
        c_inv_2 = (435+sqrt(26025.0))/4.0;
  else
        c_inv_2 = (244+sqrt(9136.0))/3.0;
         break; 
    default: 
         c_inv_2 = (435+sqrt(26025.0))/4.0;
         break;     
  }

  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 1:
      // standard parameter 
      delta =  hK*hK;
      delta =  delta0*delta;
       //OutPut(delta << endl);
       break;
    case 2:                                      
      // parameter from [BBJL07]
      delta =  hK*hK/(r*r*(eps+c));
      delta =  delta0*delta;
  //OutPut(delta << endl);
       break;
      // paper with Julia Novo
    case 3:
         delta = hK*hK/(3*eps*c_inv_2);
  //OutPut(delta << " ");
  if (c>0)
  {
    val = 1.0/(3.0*c);
    if (val<delta)
      delta = val;
  }
  //OutPut(delta << " ");
  norm_b= fabs(b1);
  if (fabs(b2)>norm_b)
    norm_b = fabs(b2);
  if (norm_b > 0)
  {
  val = beta_0 *hK/(4*norm_b*sqrt(c_inv_2));
  if (val<delta)
          delta = val;
  }
  //OutPut(delta << " ");
  M = 14*eps;
  val = 21 * c * c_pf *  c_pf ;
  if (val > M)
    M = val;
  val = 14 * delta1;
  if (val > M)
    M = val;
  if (c>0)
  {
    val = 21 * norm_b * norm_b/c;
    if (val > M)
      M = val;
  }
  M /= (beta_0*beta_0);
  val = hK * hK /(24 * M * c_inv_2);
  if (val<delta)
          delta = val;

  //OutPut(delta << endl);
     break;
    case 4:
         delta = hK*hK/(2*eps*c_inv_2);
  // OutPut(delta << " ");
  norm_b= fabs(b1);
  if (fabs(b2)>norm_b)
    norm_b = fabs(b2);
  if (norm_b > 0)
  {
          val = beta_1 * hK/(4*norm_b*sqrt(c_inv_2));
    if (val<delta)
            delta = val;
  }
  //OutPut(delta << " ");
  C0 = delta;
  M = 10 * eps * c_inv_2/(C_h*C_h);
  val = 10.0/(C_h*C_h*delta_0);
  if (val > M)
    M = val;
  val = 10 * delta1 * c_inv_2/(C_h*C_h);
  if (val > M)
    M = val;
  val = 10 * C0 * norm_b * norm_b * c_inv_2/(C_h*C_h);
  if (val > M)
    M = val;  
  M /= (beta_1*beta_1);
  val = beta_1 * hK * hK /(16 * M);
  if (val<delta)
          delta = val;

  //OutPut(delta << " " << TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE << endl);
     break;
   case 11:
      // standard parameter 
      if (eps < hK)
        delta = hK;
      else
        delta = hK * hK;
      delta =  delta0*delta;
        //OutPut(delta << endl);
       break;
   case 12:
       delta =  hK;
       delta =  delta0*delta;
         //OutPut(delta << endl);
        break;

    default:
      OutPut("SDFEM_TYPE " << TDatabase::ParamDB->SDFEM_TYPE << " not defined !!!" << endl);
      exit(4711);
  }
  return(delta);
}


// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 1, (reduced) SDFEM or (simplified RFB)
// ======================================================================
void NSType1SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double tau, ugrad;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  //delta = delta0*hK;
  tau = RFB_Parameter(hK,c0,&param[0]); // stabilization parameter

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = test00 + tau * (u1*test10+u2*test01);
    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*ugrad;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  //  cout << "c3";

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      //    val += c3 * Orig2[j] *test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky
// ======================================================================
void NSType1Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 1, VMSProjection, nabla form
// ======================================================================
void NSType1VMSProjection(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11;
  double **MatrixB1, **MatrixB2;
  double **MatrixL, **Matrix_tilde_G11 ;
  double **Matrix_tilde_G22, **Matrix_G11, **Matrix_G22;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2;
  double u1, u2, mu, delta;
  //cout << "VMS" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixL   = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  Matrix_tilde_G11  = LocMatrices[4];
  Matrix_tilde_G22  = LocMatrices[5];
  Matrix_G11  = LocMatrices[6];
  Matrix_G22  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // l

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
  if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM 
     || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
    param[0] = u1;
    param[1] = u2;
  }

  // compute the size of the mesh lenght for the turbulent viscosity
  delta =  CharacteristicFilterWidth(hK);
  // compute the turbulent viscosity
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // assemble right hand side
    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      // assemble matrix A (viscous + convective + Smagorinsky term)
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      // assemble matrices B1 and B2, divergence constraint	
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i

  // matrices for VMS
  for(i=0;i<N_U;i++)
  {
    // projections in the momentum balance
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row  = Matrix_tilde_G22[i];
    test10 = Orig0[i];  // velo_x
    test01 = Orig1[i];  // velo_y

     for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];  // large scales 
      Matrix11Row[j] -= Mult * mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * mu * ansatz00 * test01;
    }
  }

  for(i=0;i<N_L;i++)
  {
    // elliptic projection
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    test00 = Orig4[i]; // large scales

    for(j=0;j<N_U;j++)
    {
	ansatz10 = Orig0[j];  // velo_x
	ansatz01 = Orig1[j];  // velo_y

      Matrix11Row[j] -= Mult * ansatz10 * test00;
      Matrix22Row[j] -= Mult * ansatz01 * test00;
    }
  }

  for(i=0;i<N_L;i++)
  {
    // mass matrix of large scale space
    test00 = Orig4[i];
    MatrixRow1 = MatrixL[i];
    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      MatrixRow1[j] += Mult * ansatz00 * test00;
    }
  }
}

// ======================================================================
// Type 1, Standard Galerkin + divergence term 
// ======================================================================
void NSType1GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;
    }                                             // endfor j
  }                                               // endfor i
 for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                                             // endfor j

  }                                               // endfor i
  //cout << " end. " << endl;
}



// ======================================================================
// Type 2, Standard Galerkin
// ======================================================================
void NSType2Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;
      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;
      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c;
  double u1, u2;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u_xx
  Orig4 = OrigValues[4];         // u_yy
  Orig5 = OrigValues[5];         // p_x
  Orig6 = OrigValues[6];         // p_y
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  
  // for computational comparisons of Oseen problems
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
    c  = coeff[5];
    param[0] = u1;
    param[1] = u2;
  }

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      // standard terms
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SD term
      val += (-c0*(ansatz20+ansatz02)
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // for computational comparisons of Oseen problems
      if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
      {
        ansatz00 = Orig2[j];
        val += c * ansatz00 * test00;
        val += c * ansatz00 * ugrad;
      }

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;

      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig7[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void NSType2Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void NSType2Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 2, Standard Galerkin + divergence term
// ======================================================================
void NSType2GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;
    }                                             // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;
      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                                               // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;
      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                                             // endfor j

  }                                               // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3GalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, ngu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin + div term, (grad u, grad v)
// ======================================================================
void NSType3GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;
// OutPut("divu =" << divu << endl);
      Matrix11Row[j] += Mult * val;


      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;
      Matrix22Row[j] += Mult * val;

    }                                             // endfor j
  }                                               // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                                             // endfor j

  }                                               // endfor i
}


void CSTGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixS11, **MatrixS12, **MatrixS21, **MatrixS22, **MatrixS23, **MatrixS32, **MatrixS33;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_S;
  double c0, c1, c2, c3;
  double u1,u2, u1x, u1y, u2x, u2y;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix32Row, *Matrix33Row;

  MatrixS11 = LocMatrices[0];
  MatrixS12 = LocMatrices[1];
  MatrixS21 = LocMatrices[2];
  MatrixS22 = LocMatrices[3];
  MatrixS23 = LocMatrices[4];
  MatrixS32 = LocMatrices[5];
  MatrixS33 = LocMatrices[6];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // s_x
  Orig1 = OrigValues[1];         // s_y
  Orig2 = OrigValues[2];         // s

  c0 = coeff[0];                 // 1/Wi
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                //u1xold
  u2x = param[3];                //u2xold
  u1y = param[4];                //u1yold
  u2y = param[5];                //u2yold
  
  for(i=0;i<N_S;i++)
  {
    Matrix11Row = MatrixS11[i];
    Matrix12Row = MatrixS12[i];
    Matrix21Row = MatrixS21[i];
    Matrix22Row = MatrixS22[i];
    Matrix23Row = MatrixS23[i];
    Matrix32Row = MatrixS32[i];
    Matrix33Row = MatrixS33[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*(c1);
    Rhs2[i] += Mult*test00*(c2);
    Rhs3[i] += Mult*test00*(c3);
    
    for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = (c0*ansatz00*test00) - (2.0*ansatz00*u1x*test00) + (((u1*ansatz10) + (u2*ansatz01))*test00)  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00;
      Matrix12Row[j] += Mult * val;
      
      val = -ansatz00*u2x*test00;
      Matrix21Row[j] += Mult * val;
      
      val = (c0*ansatz00*test00) + ((((u1*ansatz10) + (u2*ansatz01))*test00)) ;
      Matrix22Row[j] += Mult * val;
       
      val = -ansatz00*u1y*test00;
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix32Row[j] += Mult * val;

      val  = (c0*ansatz00*test00)  - (2.0*ansatz00*u2y*test00) + (((u1*ansatz10) + (u2*ansatz01))*test00);
      Matrix33Row[j] += Mult * val;
       

    }                            // endfor j

  }                              // endfor i

}

void CST_SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixS11, **MatrixS12, **MatrixS21, **MatrixS22, **MatrixS23, **MatrixS32, **MatrixS33;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_S;
  double c0, c1, c2, c3;
  double u1,u2, u1x, u1y, u2x, u2y, ugrad, delta, norm_u;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix32Row, *Matrix33Row;
  
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  
  MatrixS11 = LocMatrices[0];
  MatrixS12 = LocMatrices[1];
  MatrixS21 = LocMatrices[2];
  MatrixS22 = LocMatrices[3];
  MatrixS23 = LocMatrices[4];
  MatrixS32 = LocMatrices[5];
  MatrixS33 = LocMatrices[6];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // s_x
  Orig1 = OrigValues[1];         // s_y
  Orig2 = OrigValues[2];         // s

  c0 = coeff[0];                 // 1/Wi
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                //u1xold
  u2x = param[3];                //u2xold
  u1y = param[4];                //u1yold
  u2y = param[5];                //u2yold
  
  norm_u = sqrt((u1*u1) + (u2*u2));
  delta = delta0*hK/norm_u;

  for(i=0;i<N_S;i++)
  {
    Matrix11Row = MatrixS11[i];
    Matrix12Row = MatrixS12[i];
    Matrix21Row = MatrixS21[i];
    Matrix22Row = MatrixS22[i];
    Matrix23Row = MatrixS23[i];
    Matrix32Row = MatrixS32[i];
    Matrix33Row = MatrixS33[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;
    Rhs3[i] += Mult*(test00+ugrad)*c3;
    
    for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = (c0*ansatz00*(test00+ugrad)) - (2.0*ansatz00*u1x*(test00+ugrad)) + (((u1*ansatz10) + (u2*ansatz01))*(test00+ugrad))  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*(test00+ugrad);
      Matrix12Row[j] += Mult * val;
      
      val = -ansatz00*u2x*(test00+ugrad);
      Matrix21Row[j] += Mult * val;
      
      val = (c0*ansatz00*(test00+ugrad)) + ((((u1*ansatz10) + (u2*ansatz01))*(test00+ugrad))) ;
      Matrix22Row[j] += Mult * val;
       
      val = -ansatz00*u1y*(test00+ugrad);
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*(test00+ugrad);
      Matrix32Row[j] += Mult * val;

      val  = (c0*ansatz00*(test00+ugrad))  - (2.0*ansatz00*u2y*(test00+ugrad)) + (((u1*ansatz10) + (u2*ansatz01))*(test00+ugrad));
      Matrix33Row[j] += Mult * val;
       

    }                            // endfor j

  }                              // endfor i

}

void CST_SUFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixS11, **MatrixS12, **MatrixS21, **MatrixS22, **MatrixS23, **MatrixS32, **MatrixS33;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_S;
  double c0, c1, c2, c3;
  double u1,u2, u1x, u1y, u2x, u2y, ugrad, delta;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix32Row, *Matrix33Row;
  
  MatrixS11 = LocMatrices[0];
  MatrixS12 = LocMatrices[1];
  MatrixS21 = LocMatrices[2];
  MatrixS22 = LocMatrices[3];
  MatrixS23 = LocMatrices[4];
  MatrixS32 = LocMatrices[5];
  MatrixS33 = LocMatrices[6];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // s_x
  Orig1 = OrigValues[1];         // s_y
  Orig2 = OrigValues[2];         // s

  c0 = coeff[0];                 // 1/Wi
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                //u1xold
  u2x = param[3];                //u2xold
  u1y = param[4];                //u1yold
  u2y = param[5];                //u2yold
  

  delta = 0.5*hK/sqrt((u1*u1) + (u2*u2));

  for(i=0;i<N_S;i++)
  {
    Matrix11Row = MatrixS11[i];
    Matrix12Row = MatrixS12[i];
    Matrix21Row = MatrixS21[i];
    Matrix22Row = MatrixS22[i];
    Matrix23Row = MatrixS23[i];
    Matrix32Row = MatrixS32[i];
    Matrix33Row = MatrixS33[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    Rhs3[i] += Mult*test00*c3;
    
    for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
     val  = (c0*ansatz00*test00) - (2.0*ansatz00*u1x*test00) + (((u1*ansatz10) + (u2*ansatz01))*(test00+ugrad))  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00;
      Matrix12Row[j] += Mult * val;
      
      val = -ansatz00*u2x*test00;
      Matrix21Row[j] += Mult * val;
      
      val = (c0*ansatz00*test00) + ((((u1*ansatz10) + (u2*ansatz01))*(test00+ugrad))) ;
      Matrix22Row[j] += Mult * val;
       
      val = -ansatz00*u1y*test00;
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix32Row[j] += Mult * val;

      val  = (c0*ansatz00*test00)  - (2.0*ansatz00*u2y*test00) + (((u1*ansatz10) + (u2*ansatz01))*(test00+ugrad));
      Matrix33Row[j] += Mult * val;
       

    }                            // endfor j

  }                              // endfor i

}



// D(u):D(v)
void NSEType4_Galerkin_CST_SUPG_DEVSS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixC11, **MatrixC12, **MatrixC22, **MatrixC23;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  double **MatrixH11, **MatrixH22, **MatrixH33;
  double **MatrixE11, **MatrixE12, **MatrixE22, **MatrixE23;
  double **MatrixJ11, **MatrixJ21, **MatrixJ22, **MatrixJ32;
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, val;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P, N_S, N_D;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8;
  double u1,u2, u1x, u1y, u2x, u2y;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y, delta, norm_u, ugrad;
  double delta0 = TDatabase::ParamDB->DELTA0;
    
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  
  MatrixG11 = LocMatrices[4];
  MatrixG12 = LocMatrices[5];
  MatrixG21 = LocMatrices[6];
  MatrixG22 = LocMatrices[7];
  MatrixG23 = LocMatrices[8];
  MatrixG32 = LocMatrices[9];
  MatrixG33 = LocMatrices[10];
  
  MatrixH11 = LocMatrices[11];
  MatrixH22 = LocMatrices[12];
  MatrixH33 = LocMatrices[13];
  
  MatrixB1 = LocMatrices[14];
  MatrixB2 = LocMatrices[15];
  
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  
  MatrixC11 = LocMatrices[18];
  MatrixC12 = LocMatrices[19];
  MatrixC22 = LocMatrices[20];
  MatrixC23 = LocMatrices[21];
  
  MatrixD11 = LocMatrices[22];
  MatrixD12 = LocMatrices[23];
  MatrixD21 = LocMatrices[24];
  MatrixD22 = LocMatrices[25];
  MatrixD31 = LocMatrices[26];
  MatrixD32 = LocMatrices[27];
 
  MatrixE11 = LocMatrices[28];
  MatrixE12 = LocMatrices[29];
  MatrixE22 = LocMatrices[30];
  MatrixE23 = LocMatrices[31];
  
  MatrixJ11 = LocMatrices[32];
  MatrixJ21 = LocMatrices[33];
  MatrixJ22 = LocMatrices[34];
  MatrixJ32 = LocMatrices[35];


  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_S = N_BaseFuncts[2];
  N_D = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // tau_x
  Orig5 = OrigValues[5];         // tau_y
  Orig6 = OrigValues[6];         // tau
  Orig7 = OrigValues[7];         // deformation

   c0 = coeff[0]; // D(u):D(v) term 
   c1 = coeff[1]; // f1
   c2 = coeff[2]; // f2
   c3 = coeff[3]; // Tau term in CST 
   c4 = coeff[4]; // f3
   c5 = coeff[5]; // f4 
   c6 = coeff[6]; // f5
   c7 = coeff[7]; // div Tau term 
   c8 = coeff[8]; // div Deformation term 
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold

  
  
    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      { 
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix C
    Matrix11Row = MatrixC11[i];
    Matrix12Row = MatrixC12[i];
    Matrix22Row = MatrixC22[i];
    Matrix23Row = MatrixC23[i];
     for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig6[j];

      val  = c7*test10*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix12Row[j] += Mult * val;

      val  = c7*test10*ansatz00;
      Matrix22Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix23Row[j] += Mult * val;

    }                            // endfor j

    // Matrix E
    Matrix11Row = MatrixE11[i];
    Matrix12Row = MatrixE12[i];
    Matrix22Row = MatrixE22[i];
    Matrix23Row = MatrixE23[i];
     for(j=0;j<N_D;j++)
    {
      ansatz00 = Orig7[j];

      val  = c8*test10*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c8*test01*ansatz00;
      Matrix12Row[j] += Mult * val;

      val  = c8*test10*ansatz00;
      Matrix22Row[j] += Mult * val;

      val  = c8*test01*ansatz00;
      Matrix23Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix BT
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  norm_u = sqrt((u1*u1) + (u2*u2));
  
  if (norm_u < 0.000001)
    delta=0;
  else
    delta = delta0*hK/norm_u;
  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig4[i];
    test01 = Orig5[i];
    test00 = Orig6[i];
    ugrad  = delta * (u1*test10+u2*test01);

    val = (c4 + (u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*(test00+ugrad);
    Rhs3[i] += Mult*val;
    val = (2*c5 + 2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*(test00+ugrad);
    Rhs4[i] += Mult*val;
    val = (c6 + (u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*(test00+ugrad);
    Rhs5[i] += Mult*val;
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*(test00+ugrad);
   
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*(test00+ugrad);
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*(test00+ugrad);
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*(test00+ugrad);
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*(test00+ugrad);
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*(test00+ugrad);
      Matrix32Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];

     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig6[j];

      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1*ansatz10 + u2*ansatz01)*(test00+ugrad)  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*(test00+ugrad);
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*(test00+ugrad);
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1*ansatz10 + u2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*(test00+ugrad) ;
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*(test00+ugrad);
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*(test00+ugrad);
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1*ansatz10 + u2*ansatz01)*(test00+ugrad);
      Matrix33Row[j] += Mult * val;

    }                            // endfor j

   
  }                              // endfor i

  
      for(i=0;i<N_D;i++)
  {

    test00 = Orig7[i];
    
    // Matrix J
    Matrix11Row = MatrixJ11[i];
    Matrix21Row = MatrixJ21[i];
    Matrix22Row = MatrixJ22[i];
    Matrix32Row = MatrixJ32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = -ansatz10*test00;
      Matrix11Row[j] += Mult * val;

      val  = -ansatz01*test00;
      Matrix21Row[j] += Mult * val;

      val  = -ansatz10*test00;
      Matrix22Row[j] += Mult * val;
      
      val  = -ansatz01*test00;
      Matrix32Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix H
    Matrix11Row = MatrixH11[i];
    Matrix22Row = MatrixH22[i];
    Matrix33Row = MatrixH33[i];
    
     for(j=0;j<N_D;j++)
    {
      ansatz00 = Orig7[j];

      val  = ansatz00*test00;
      Matrix11Row[j] += Mult * val;
      
      val = 2.0*ansatz00*test00;
      Matrix22Row[j] += Mult * val;

      val  = ansatz00*test00;
      Matrix33Row[j] += Mult * val;

    }                            // endfor j

   
  }                              // endfor i

  
  // Matrix B
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];
  
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  
}


// D(u):D(v)
void NSEType4_Galerkin_CST_Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixC11, **MatrixC12, **MatrixC22, **MatrixC23;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, val;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6;
  int i,j,N_U, N_P, N_S, N_D;
  double c0, c1, c2, c3, c4, c5, c6, c7;
  double u1,u2, u1x, u1y, u2x, u2y;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  
  MatrixG11 = LocMatrices[4];
  MatrixG12 = LocMatrices[5];
  MatrixG21 = LocMatrices[6];
  MatrixG22 = LocMatrices[7];
  MatrixG23 = LocMatrices[8];
  MatrixG32 = LocMatrices[9];
  MatrixG33 = LocMatrices[10];
  
  MatrixB1 = LocMatrices[11];
  MatrixB2 = LocMatrices[12];
  
  MatrixB1T = LocMatrices[13];
  MatrixB2T = LocMatrices[14];
  
  MatrixC11 = LocMatrices[15];
  MatrixC12 = LocMatrices[16];
  MatrixC22 = LocMatrices[17];
  MatrixC23 = LocMatrices[18];
  
  MatrixD11 = LocMatrices[19];
  MatrixD12 = LocMatrices[20];
  MatrixD21 = LocMatrices[21];
  MatrixD22 = LocMatrices[22];
  MatrixD31 = LocMatrices[23];
  MatrixD32 = LocMatrices[24];
 
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // tau_x
  Orig5 = OrigValues[5];         // tau_y
  Orig6 = OrigValues[6];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c1 = coeff[1]; // f1
   c2 = coeff[2]; // f2
   c3 = coeff[3]; // Tau term in CST 
   c4 = coeff[4]; // f3
   c5 = coeff[5]; // f4 
   c6 = coeff[6]; // f5
   c7 = coeff[7]; // div Tau term 
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
       if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix C
    Matrix11Row = MatrixC11[i];
    Matrix12Row = MatrixC12[i];
    Matrix22Row = MatrixC22[i];
    Matrix23Row = MatrixC23[i];
     for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig6[j];

      val  = c7*test10*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix12Row[j] += Mult * val;

      val  = c7*test10*ansatz00;
      Matrix22Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix23Row[j] += Mult * val;

    }                            // endfor j

    // Matrix BT
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig4[i];
    test01 = Orig5[i];
    test00 = Orig6[i];

    val = (c4 + (u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00;
    Rhs3[i] += Mult*val;
    val = (2*c5 + 2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00;
    Rhs4[i] += Mult*val;
    val = (c6 + (u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00;
    Rhs5[i] += Mult*val;
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00;
   
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00;
      Matrix32Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];

     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig6[j];

      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1*ansatz10 + u2*ansatz01)*test00  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00;
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1*ansatz10 + u2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*test00 ;
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00;
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1*ansatz10 + u2*ansatz01)*test00;
      Matrix33Row[j] += Mult * val;

    }                            // endfor j

   
  }                              // endfor i

      
  // Matrix B
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];
  
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  
}


 
// D(u):D(v) (Axial 3D)
void NSEType4_Galerkin_CST_Galerkin_Axial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixC11, **MatrixC12, **MatrixC22, **MatrixC23;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, val;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6;
  int i,j,N_U, N_P, N_S, N_D;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9;
  double u1,u2, u1x, u1y, u2x, u2y, r, x;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  
  MatrixG11 = LocMatrices[4];
  MatrixG12 = LocMatrices[5];
  MatrixG21 = LocMatrices[6];
  MatrixG22 = LocMatrices[7];
  MatrixG23 = LocMatrices[8];
  MatrixG32 = LocMatrices[9];
  MatrixG33 = LocMatrices[10];
  
  MatrixB1 = LocMatrices[11];
  MatrixB2 = LocMatrices[12];
  
  MatrixB1T = LocMatrices[13];
  MatrixB2T = LocMatrices[14];
  
  MatrixC11 = LocMatrices[15];
  MatrixC12 = LocMatrices[16];
  MatrixC22 = LocMatrices[17];
  MatrixC23 = LocMatrices[18];
  
  MatrixD11 = LocMatrices[19];
  MatrixD12 = LocMatrices[20];
  MatrixD21 = LocMatrices[21];
  MatrixD22 = LocMatrices[22];
  MatrixD31 = LocMatrices[23];
  MatrixD32 = LocMatrices[24];
 
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // tau_x
  Orig5 = OrigValues[5];         // tau_y
  Orig6 = OrigValues[6];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c1 = coeff[1]; // f1
   c2 = coeff[2]; // f2
   c3 = coeff[3]; // Tau term in CST 
   c4 = coeff[4]; // f3
   c5 = coeff[5]; // f4 
   c6 = coeff[6]; // f5
   c7 = coeff[7]; // div Tau term 
   c8 = coeff[8]; // Fluid type, 1-Newtonian, 2-Oldroyd-B, 3-Giesekus
   c9 = coeff[9];// Tau dot Tau in Giesekus
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold
  
  x = param[15];   // x
  
  r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check NSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1*r;
    Rhs2[i] += Mult*test00*c2*r;
    
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*((2*test10*ansatz10+ test01*ansatz01)*r);
      val += c0*2.*ansatz00*test00/r;
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00*r;
      }
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10)*r;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01)*r;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+ 2*test01*ansatz01)*r;
       if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00*r;
      }
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix C
    Matrix11Row = MatrixC11[i];
    Matrix12Row = MatrixC12[i];
    Matrix22Row = MatrixC22[i];
    Matrix23Row = MatrixC23[i];
     for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig6[j];

      if(c8 == 2 || c8==3)
      {
      val  = c7*test10*ansatz00*r;
      Matrix11Row[j] += Mult * val;

      val  = c7*test01*ansatz00*r;
      Matrix12Row[j] += Mult * val;

      val  = c7*test10*ansatz00*r;
      Matrix22Row[j] += Mult * val;

      val  = c7*test01*ansatz00*r;
      Matrix23Row[j] += Mult * val;
      }
      else if(c8 == 1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
      }
      else
      {
	cout<<"Invalid fluid type \n";
	exit(4711);
      }

    }                            // endfor j

    // Matrix BT
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*(ansatz00*test10*r);
      val += -Mult*(ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01*r;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig4[i];
    test01 = Orig5[i];
    test00 = Orig6[i];

    val = c4*test00*r;
    Rhs3[i] += Mult*val;
    val = 2*c5*test00*r;
    Rhs4[i] += Mult*val;
    val = c6 *test00*r;
    Rhs5[i] += Mult*val;
    
      if(c8 == 2 || c8==3)
    {
    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00*r;
    Rhs3[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00*r;
    Rhs4[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00*r;
    Rhs5[i] += Mult*val;
    }
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];


      if(c8==2 || c8==3)
      {
      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00*r;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00*r;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix32Row[j] += Mult * val;
      }
      else if(c8==1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix31Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
    {
	cout<<"Invalid fluid type \n";
	exit(4711);
    }

    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];

     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig6[j];

    if(c8 == 2 || c8==3)
    {
      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1*ansatz10 + u2*ansatz01)*test00*r  ;
      if(c8==3)
     {
      val += c9*tau1*ansatz00*test00*r;
     }
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
     Matrix12Row[j] += Mult * val;
     
      val = -2.0*ansatz00*u2x*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
     Matrix21Row[j] += Mult * val;
     
      val = (c3*2.0*ansatz00 + 2.0*(u1*ansatz10 + u2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*test00*r ;
      if(c8==3)
     {
      val += c9*(tau1 + tau3)*ansatz00*test00*r;
     }
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00*r;
       if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1*ansatz10 + u2*ansatz01)*test00*r;
      if(c8==3)
     {
      val += c9*tau3*ansatz00*test00*r;
     }
      Matrix33Row[j] += Mult * val;
    }
      
      else if(c8 == 1)
      {
	val  = c3*ansatz00*test00*r;
	Matrix11Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val*2.0;
	Matrix33Row[j] += Mult * val;
	val = 0;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
     {
	cout<<"Invalid fluid type \n";
	exit(4711);
     }

    }                            // endfor j

   
  }                              // endfor i

      
  // Matrix B
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];
  
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = -Mult*(test00*ansatz10*r);
      val += -Mult*(ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01*r;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  
}



void CST_GiesekusGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixS11, **MatrixS12, **MatrixS21, **MatrixS22, **MatrixS23, **MatrixS32, **MatrixS33;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_S;
  double c0, c1, c2, c3, c4;
  double u1,u2, u1x, u1y, u2x, u2y, tau1, tau2, tau3;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix32Row, *Matrix33Row;
  
  MatrixS11 = LocMatrices[0];
  MatrixS12 = LocMatrices[1];
  MatrixS21 = LocMatrices[2];
  MatrixS22 = LocMatrices[3];
  MatrixS23 = LocMatrices[4];
  MatrixS32 = LocMatrices[5];
  MatrixS33 = LocMatrices[6];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // s_x
  Orig1 = OrigValues[1];         // s_y
  Orig2 = OrigValues[2];         // s

  c0 = coeff[0];                 // 1/Wi
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3
  c4 = coeff[4];                 // alpha
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];
  u2x = param[3];
  u1y = param[4];
  u2y = param[5];
  tau1 = param[6];
  tau2 = param[7];
  tau3 = param[8];
  
  for(i=0;i<N_S;i++)
  {
    Matrix11Row = MatrixS11[i];
    Matrix12Row = MatrixS12[i];
    Matrix21Row = MatrixS21[i];
    Matrix22Row = MatrixS22[i];
    Matrix23Row = MatrixS23[i];
    Matrix32Row = MatrixS32[i];
    Matrix33Row = MatrixS33[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*(c1);
    Rhs2[i] += Mult*test00*(c2);
    Rhs3[i] += Mult*test00*(c3);
    
    for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = (c0*(1-(2.0*c4))*ansatz00*test00) - (2.0*ansatz00*u1x*test00) + (((u1*ansatz10) + (u2*ansatz01))*test00) + (c4*c0*tau1*ansatz00*test00)  ;
      Matrix11Row[j] += Mult * val;

      val = (-2.0*ansatz00*u1y*test00) + (c4*c0*tau2*ansatz00*test00);
      Matrix12Row[j] += Mult * val;
      
      val = (-ansatz00*u2x*test00) + (c4*c0*0.5*tau2*ansatz00*test00);
      Matrix21Row[j] += Mult * val;
      
      val = (c0*(1-(2.0*c4))*ansatz00*test00) + ((((u1*ansatz10) + (u2*ansatz01))*test00)) + (c4*c0*0.5*(tau1+tau3)*ansatz00*test00) ;
      Matrix22Row[j] += Mult * val;
       
      val = (-ansatz00*u1y*test00)+ (c4*c0*0.5*tau2*ansatz00*test00);
      Matrix23Row[j] += Mult * val;
      
      val = (-2.0*ansatz00*u2x*test00) + (c4*c0*tau2*ansatz00*test00);
      Matrix32Row[j] += Mult * val;

      val  = (c0*(1-(2.0*c4))*ansatz00*test00)  - (2.0*ansatz00*u2y*test00) + (((u1*ansatz10) + (u2*ansatz01))*test00) + (c4*c0*tau3*ansatz00*test00);
      Matrix33Row[j] += Mult * val;
       

    }                            // endfor j

  }                              // endfor i

}

void DFTGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixS11, **MatrixS22, **MatrixS33;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double ansatz00;
  double test00;
  double *Orig0;
  int i,j,N_S;
  double c1, c2, c3;
  double u1,u2, u1x, u1y, u2x, u2y;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  
  MatrixS11 = LocMatrices[0];
  MatrixS22 = LocMatrices[1];
  MatrixS33 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_S = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // d

  u1x = param[2];                //u1xold
  u2x = param[3];                //u2xold
  u1y = param[4];                //u1yold
  u2y = param[5];                //u2yold
  int c=0;
  for(i=0;i<N_S;i++)
  {
    Matrix11Row = MatrixS11[i];
    Matrix22Row = MatrixS22[i];
    Matrix33Row = MatrixS33[i];

    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*(u1x);
    Rhs2[i] += Mult*test00*(u1y+u2x)*0.5;
    Rhs3[i] += Mult*test00*(u2y);
    
    for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig0[j];
      
      val  =  test00*ansatz00 ;
      Matrix11Row[j] += Mult * val;
    
      val =  test00*ansatz00;
      Matrix22Row[j] += Mult * val;
      

      val  = test00*ansatz00;
      Matrix33Row[j] += Mult * val;
       

    }                            // endfor j

  }                              // endfor i

}



// ======================================================================
// rhs for  NSE due to CST
// ======================================================================
void NSCSTRhs(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test10, test01, test00, tau1, tau2, tau3, c1, c2;
  double *Orig0, *Orig1, *Orig2;
  int i,N_U;
  
  double beta = TDatabase::ParamDB->P3, nu = 1.0/TDatabase::ParamDB->WEI_NR;

  double nondim1 = (beta-1.0)*nu;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  
  tau1 = param[6];
  tau2 = param[7];
  tau3 = param[8];
  
  c1 = coeff[1];                
  c2 = coeff[2];                 
  
  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    Rhs1[i] += Mult*(tau1*test10 + tau2*test01)*nondim1 + Mult*c1*test00;
    Rhs2[i] += Mult*(tau2*test10 + tau3*test01)*nondim1 + Mult*c2*test00;
  }                              // endfor i
}

// ======================================================================
// rhs for  NSE due to CST using DEVSS
// ======================================================================
void NSCSTRhs_DEVSS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test10, test01, test00, tau1, tau2, tau3, c1, c2, d1, d2, d3;
  double *Orig0, *Orig1, *Orig2;
  int i,N_U;

  double  beta = TDatabase::ParamDB->P3, nu = 1/TDatabase::ParamDB->WEI_NR;

  double nondim1 = (beta-1.0)*nu;
  double nondim2 = 2*(1.0-beta);
    

  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  tau1 = param[6];
  tau2 = param[7];
  tau3 = param[8];

  d1 = param[9];
  d2 = param[10];
  d3 = param[11];
  
  c1 = coeff[1];                
  c2 = coeff[2];                 
  

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    Rhs1[i] += Mult*(tau1*test10 + tau2*test01)*nondim1 + Mult*(d1*test10 + d2*test01)*nondim2 + Mult*c1*test00;
    Rhs2[i] += Mult*(tau2*test10 + tau3*test01)*nondim1 + Mult*(d2*test10 + d3*test01)*nondim2 + Mult*c2*test00;
  }  
                            // endfor i
                            
  
                       
}

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType4Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1,u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

//       val  = 0;
//       Matrix12Row[j] += Mult * val;
// 
//       val  = 0;
//       Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4GalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad, r=2.0, maxu, tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u_xx
  Orig4 = OrigValues[4];         // u_yy
  Orig5 = OrigValues[5];         // p_x
  Orig6 = OrigValues[6];         // p_y
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
   if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM
      || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
  }

  maxu = fabs(u1);
  if (fabs(u2)>maxu)
      maxu = fabs(u2);

  // parameter from [BBJL07]
  delta =  hK*hK/(r*r*(c0+1));
  delta =  delta0*delta;

  tau =  delta1 *(1+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      // viscous term
      val  = c0*(test10*ansatz10+test01*ansatz01);
      // convective term
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SUPG term
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // term in both diagonal blocks
      Matrix11Row[j] += Mult * (val+tau*test10*ansatz10);
      Matrix22Row[j] += Mult * (val+tau*test01*ansatz01);
      // div-div term
      val = tau * test10*ansatz01; 
      Matrix12Row[j] += Mult * val;

      val = tau * test01*ansatz10; 
      Matrix21Row[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig7[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SDFEMDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u_xx
  Orig4 = OrigValues[4];         // u_yy
  Orig5 = OrigValues[5];         // p_x
  Orig6 = OrigValues[6];         // p_y
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig7[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void NSType4Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
  
  // including ROSSBY term, ONLY FOR STOKES SO FAR
  if (fabs(TDatabase::ParamDB->ROSSBY_NR)>1e-10)
  {
    for(i=0;i<N_U;i++)
    {
      Matrix12Row = MatrixA12[i];
      Matrix21Row = MatrixA21[i];
      test00 = Orig2[i];
      for(j=0;j<N_U;j++)
      {
        ansatz00 = Orig2[j];
        val  = 2*TDatabase::ParamDB->ROSSBY_NR * ansatz00 * test00;
        Matrix12Row[j] -= Mult * val;
        Matrix21Row[j] += Mult * val;
      }   // endfor j
    }
  }
}


// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void NSType4UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  // including ROSSBY term, ONLY FOR STOKES SO FAR
  if (fabs(TDatabase::ParamDB->ROSSBY_NR)>1e-10)
  {
    for(i=0;i<N_U;i++)
    {
      Matrix12Row = MatrixA12[i];
      Matrix21Row = MatrixA21[i];
      test00 = Orig2[i];
      for(j=0;j<N_U;j++)
      {
        ansatz00 = Orig2[j];
        val  = 2*TDatabase::ParamDB->ROSSBY_NR * ansatz00 * test00;
        Matrix12Row[j] -= Mult * val;
        Matrix21Row[j] += Mult * val;
      }   // endfor j
    }
  }
}


// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType4SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, VMSProjection, D(u):D(v) , non-linear term switched off
// ======================================================================
void NSType4VMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixL, **Matrix_tilde_G11 ;
  double **Matrix_tilde_G22, **Matrix_G11, **Matrix_G22;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixL   = LocMatrices[4];
  MatrixB1 = LocMatrices[5];
  MatrixB2 = LocMatrices[6];
  MatrixB1T = LocMatrices[7];
  MatrixB2T = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_G11  = LocMatrices[11];
  Matrix_G22  = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // l

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old


  mu = (1.0-TDatabase::ParamDB->P3)/TDatabase::ParamDB->RE_NR;
  

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
//       val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
//       val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row  = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      Matrix11Row[j] -= Mult * 2*mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * 2*mu * ansatz00 * test01;
    }
  }

  for(i=0;i<N_L;i++)
  {
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    test00 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      Matrix11Row[j] -= Mult * ansatz10 * test00;
      Matrix22Row[j] -= Mult * ansatz01 * test00;
    }
  }

  for(i=0;i<N_L;i++)
  {
    test00 = Orig4[i];
    MatrixRow1 = MatrixL[i];
    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      MatrixRow1[j] += Mult * ansatz00 * test00;
    }
  }
}




// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// Type 1,(reduced) SDFEM or (simplified) RFB
// ======================================================================
void NSType1NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow, *Rhs1, *Rhs2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, c1, c2;
  double u1, u2;
  double tau, ugrad;

  MatrixA = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  tau = RFB_Parameter(hK,c0,&param[0]); // stabilization parameter

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = test00+tau * (u1*test10+u2*test01);
  
    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*ugrad;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1_2NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// ======================================================================
void NSType1_2NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, VMSProjection, nabla form
// Type 2, VMSProjection, nabla form
// ======================================================================
void NSType1_2NLVMSProjection(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double **Matrix_tilde_G11,  **Matrix_tilde_G22;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_L;
  double c0, viscosity, delta;
  double u1, u2, mu;
  //cout << "VMS" << endl;
  MatrixA11 = LocMatrices[0];
  Matrix_tilde_G11  = LocMatrices[1];
  Matrix_tilde_G22  = LocMatrices[2];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // l

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
  if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM
     || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
      u1 = coeff[3];
      u2 = coeff[4];
      // for the vms coefficient, still the params computed from the code
      // are used
      param[0] = u1;
      param[1] = u2;
      /*param[2] = coeff[5];
      param[3] = coeff[6];
      param[4] = coeff[7];
      param[5] = coeff[8];*/
  }

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  viscosity = mu+c0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig3[j];
      Matrix11Row[j] -= Mult * mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * mu * ansatz00 * test01;
    }
  }
}

// ======================================================================
// Type 1, Standard Galerkin + div term, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, du1x, du2y, divu;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u

  c0 = coeff[0];                                  // nu

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;
    }                                             // endfor j
  }                                               // endfor i
}


// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA = LocMatrices[0];
  MatrixB1T = LocMatrices[1];
  MatrixB2T = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u_xx
  Orig4 = OrigValues[4];         // u_yy
  Orig5 = OrigValues[5];         // p_x
  Orig6 = OrigValues[6];         // p_y
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK ;
  else
    delta = delta1*hK*hK ;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
	      + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

//   cout << " u1 " << u1 << " u2 " << u2 <<endl;
  
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3_4NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j,N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3_4NLUpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j,N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3_4NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U, N_P;
  double c0;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
      /* do not see why these values are added, 06/03/20
      val  = mu*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = mu*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;
      */
      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3_4NLSmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;
  double mu, delta;
  //cout << "hier" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin + div term, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, du1x, du2y, divu;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u

  c0 = coeff[0];                                  // nu

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;
      
      Matrix22Row[j] += Mult * val;

    }                                             // endfor j
  }                                               // endfor i
}


// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad, r=2.0, maxu, tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixB1T = LocMatrices[2];
  MatrixB2T = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u_xx
  Orig4 = OrigValues[4];         // u_yy
  Orig5 = OrigValues[5];         // p_x
  Orig6 = OrigValues[6];         // p_y
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
   if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM
      || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
      u1 = coeff[3];
      u2 = coeff[4];
  }

  maxu = fabs(u1);
  if (fabs(u2)>maxu)
      maxu = fabs(u2);

  // parameter from [BBJL07]
  delta =  hK*hK/(r*r*(c0+1));
  delta =  delta0*delta;

  tau =  delta1 *(1+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      // viscous term
      val  = c0*(test10*ansatz10+test01*ansatz01);
      // convective term
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SUPG term
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // term in both diagonal blocks
      Matrix11Row[j] += Mult * (val+tau*test10*ansatz10);
      Matrix22Row[j] += Mult * (val+tau*test01*ansatz01);
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixB1T = LocMatrices[2];
  MatrixB2T = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u_xx
  Orig4 = OrigValues[4];         // u_yy
  Orig5 = OrigValues[5];         // p_x
  Orig6 = OrigValues[6];         // p_y
  Orig7 = OrigValues[7];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}


// ======================================================================
// auxiliary problem
// ======================================================================
void NSAuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **AuxMatrix;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double mu2, delta, mu;
  double u1, u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  AuxMatrix = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  //cout << " " << delta;
  mu2 = 0.25*delta*delta/gamma;
  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*u1;
    Rhs2[i] += Mult*test00*u2;
    //    Rhs1[i] += Mult*test00*coeff[3];
    //    Rhs2[i] += Mult*test00*coeff[3];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// auxiliary problem for differential filter
// ======================================================================
void Filter_Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **AuxMatrix;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double delta, mu;
  double u1, u2;

  AuxMatrix = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  delta *= delta;

  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*u1;
    Rhs2[i] += Mult*test00*u2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = delta*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// rhs for RFB stabilization
// ======================================================================
void NSRFBRhs(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i,N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];
    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
  }                              // endfor i
}

// ======================================================================
// pressure separation
// ======================================================================
void NSPressSep(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double test00;
  double *Orig0;
  int i,j, N_U;
  double c1, c2;
  double p_x, p_y;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  p_x = param[0];                // p_sep to x
  p_y = param[1];                // p_sep to y

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];
    Rhs1[i] += Mult*test00*(c1-p_x);
    Rhs2[i] += Mult*test00*(c2-p_y);
  }                              // endfor i
}


// ======================================================================
// pressure separation
// ======================================================================
void NSPressSepAuxProb(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, *MatrixRow;
  double *Rhs1, val;
  double test10,test01,ansatz10,ansatz01;
  double *Orig0, *Orig1;
  double u1, u2, d1u1, d2u1, d1u2, d2u2;
  int i,j, N_U;
  double c1, c2;
  double p_x, p_y;

  MatrixA = LocMatrices[0];
  Rhs1 = LocRhs[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  if (TDatabase::ParamDB->PRESSURE_SEPARATION==3)
  {
    for(i=0;i<N_U;i++)
    {
      MatrixRow = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      Rhs1[i] += Mult*(test10*c1+test01*c2);
      for(j=0;j<N_U;j++)
      {
        ansatz10 = Orig0[j];
        ansatz01 = Orig1[j];
        MatrixRow[j] += Mult * (ansatz10*test10+ansatz01*test01);
      }
    }
  }
  else
  {
    u1 = param[0];               // u1old
    u2 = param[1];               // u2old
    d1u1 = param[2];
    d1u2 = param[3];
    d2u1 = param[4];
    d2u2 = param[5];
    for(i=0;i<N_U;i++)
    {
      MatrixRow = MatrixA[i];
      test10 = Orig0[i];
      test01 = Orig1[i];
      Rhs1[i] += Mult*(test10*(c1-u1*d1u1-u2*d2u1)
        +test01*(c2-u1*d1u2-u2*d2u2));
      for(j=0;j<N_U;j++)
      {
        ansatz10 = Orig0[j];
        ansatz01 = Orig1[j];
        MatrixRow[j] += Mult * (ansatz10*test10+ansatz01*test01);
      }
    }
  }
}

// ========================================================================
// routines for FJMT07
// ========================================================================

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1GalerkinFJMT07(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, tau;
  double sigma = 1/TDatabase::TimeDB->TIMESTEPLENGTH;

  //OutPut("FJMT");

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = coeff[3];                 // u1old
  u2 = coeff[4];                 // u2old

   // SD parameter
  tau = TDatabase::ParamDB->DELTA0*hK*hK;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*(test00*c1 + tau*c1*(u1*test10+u2*test01));
    Rhs2[i] += Mult*(test00*c2 + tau*c2*(u1*test10+u2*test01));

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += tau*(u1*ansatz10+u2*ansatz01)*(u1*test10+u2*test01);
      val += sigma * ansatz00 * test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkinFJMT07(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2,tau;
  double sigma = 1/TDatabase::TimeDB->TIMESTEPLENGTH;
  //OutPut("07");
  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = coeff[3];                 // u1old
  u2 = coeff[4];                 // u2old

  // SD parameter
  tau = TDatabase::ParamDB->DELTA0*hK*hK;
  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += sigma*ansatz00*test00;
      val += tau*(u1*ansatz10+u2*ansatz01)*(u1*test10+u2*test01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void NSParamsVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
}

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
}

// =======================================================================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, gradient(tau1), gradient(tau2), gradient(tau3)
// ========================================================================================================================
void NSParamsVelo_GradVelo_CST(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
  out[6] = in[8];                // D00(tau1old)
  out[7] = in[9];                // D00(tau2old)
  out[8] = in[10];               // D00(tau3old)
  out[9]  = in[11];              // D10 (tau1old)
  out[10] = in[12];              // D10 (tau2old)
  out[11] = in[13];              // D10 (tau3old)
  out[12] = in[14];              // D01 (tau1old)
  out[13] = in[15];             // D01 (tau2old)
  out[14] = in[16];             // D01 (tau3old)
}


// =======================================================================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, gradient(tau1), gradient(tau2), gradient(tau3)
// ========================================================================================================================
void NSParamsVelo_GradVelo_CST_Axial3D(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
  out[6] = in[8];                // D00(tau1old)
  out[7] = in[9];                // D00(tau2old)
  out[8] = in[10];               // D00(tau3old)
  out[9]  = in[11];              // D10 (tau1old)
  out[10] = in[12];              // D10 (tau2old)
  out[11] = in[13];              // D10 (tau3old)
  out[12] = in[14];              // D01 (tau1old)
  out[13] = in[15];             // D01 (tau2old)
  out[14] = in[16];             // D01 (tau3old)
  out[15] = in[0];              // x value for axial symmetric
}




// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, D1, D2, D3
// ========================================================================
void NSParamsVelo_GradVelo_DEVSS(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
  out[6] = in[8];                // D00(tau1old)
  out[7] = in[9];                // D00(tau2old)
  out[8] = in[10];               // D00(tau3old)
  out[9]  = in[11];              // D00 (D1old)
  out[10] = in[12];              // D00 (D2old)
  out[11] = in[13];              // D00 (D3old)
  out[12] = in[14];              // D00 (pold)

}

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), grad (tau1, tau2, tau3), D1, D2, D3
// ========================================================================
void NSParamsVelo_GradVelo_CST_DEVSS(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
  out[6] = in[8];                // D00(tau1old)
  out[7] = in[9];                // D00(tau2old)
  out[8] = in[10];               // D00(tau3old)
  out[9] = in[11];               // D10(tau1old)
  out[10] = in[12];              // D10(tau2old)
  out[11] = in[13];              // D10(tau3old)
  out[12] = in[14];              // D01(tau1old)
  out[13] = in[15];              // D01(tau2old)
  out[14] = in[16];              // D01(tau3old)
  out[15] = in[17];              // D00 (D1old)
  out[16] = in[18];              // D00 (D2old)
  out[17] = in[19];              // D00 (D3old)

 }




void NSParams_CST(double *in, double *out)
{
  out[0] = in[2];                // tau1old
  out[1] = in[3];                // tau2old
  out[2] = in[4];                // tau3old
  out[3] = in[5];                // tau1old D10
  out[4] = in[6];                // tau2old D10
  out[5] = in[7];                // tau3old D10
  out[6] = in[8];                // tau1old D01
  out[7] = in[9];                // tau2old D01
  out[8] = in[10];               // tau3old D01
}



void NSParamsVeloExact(double *in, double *out)
{
  double x,y, t1, t2, t3, t4, t5, t6, t7, t9, t10, t13, t15, u1, u2;
  double t20, t23;
  double scale = TDatabase::ParamDB->P7;

  x = in[0];
  y = in[1];
  /*  t1 = x*x;
    t2 = 1.0-x;
    t3 = t2*t2;
    t4 = t1*t3;
    t5 = 1.0-y;
    t6 = t5*t5;
    t7 = y*t6;
    t9 = y*y;
    t10 = t9*t5;
    t13 = x*t3;
    t15 = t1*t2;
    u1 = 2.0*t4*t7-2.0*t4*t10;

    t1 = 1.0-x;
    t2 = t1*t1;
    t3 = x*t2;
    t4 = y*y;
    t5 = 1.0-y;
    t6 = t5*t5;
    t7 = t4*t6;
    t9 = x*x;
    t10 = t9*t1;
    t20 = y*t6;
    t23 = t4*t5;
    u2 = -2.0*t3*t7+2.0*t10*t7;
  */
  u1  =2*Pi*sin(Pi*y)*cos(Pi*y)*sin(Pi*x)*sin(Pi*x)*sin(Pi*x);
  u2 = -3*Pi*sin(Pi*x)*sin(Pi*x)*cos(Pi*x)*sin(Pi*y)*sin(Pi*y);
  out[0] = u1-in[2];             // u1old
  out[1] = u2-in[3];             // u2old
  //out[0] = in[2]; // u1old
  // out[1] = in[3]; // u2old
}


// ========================================================================
// parameters: separated pressure
// ========================================================================
void NSParamsPressSep(double *in, double *out)
{
  out[0] = in[2];                // P_sep to x
  out[1] = in[3];                // P_sep to y
}
