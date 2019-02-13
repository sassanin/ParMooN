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
// @(#)NSE2D_EQU_ORD_FixPo.C        1.3 06/27/00
//
// common declaration for all Navier-Stokes problems
// fix point iteration
// for equal order interpolation
//
// ======================================================================

#include <Database.h>
#include <TNSE2D_Routines.h>

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEMEquOrd(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T, **MatrixC;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRowC;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2, px, py;
  double delta, tau, ugrad, r=2.0, maxu;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixC = LocMatrices[4];
  MatrixB1 = LocMatrices[5];
  MatrixB2 = LocMatrices[6];
  MatrixB1T = LocMatrices[7];
  MatrixB2T = LocMatrices[8];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

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
  maxu = sqrt(u1*u1+u2*u2);

  if (TDatabase::ParamDB->NSTYPE > 10)
  {
      // equal order parameter from [BBJL07]
      delta = r*r*r*r*c0/(hK*hK) + r*maxu/hK;
      delta =  delta0/delta;
      
      tau = r*r*c0 + hK*maxu/r;
      tau *= delta1;
  }
  else
  {
      // inf-sub stable parameter from [BBJL07]
      delta =  hK*hK/(r*r*(c0+1));
      delta =  delta0*delta;
      
      tau =  delta1 *(1+c0);
  }
  
  // assembling for velocity test functions
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

    // rhs
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    // velocity-velocity block
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

    // pressure-velocity block
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      // pressure ansatz functions
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];
     
      // pressure term 
      val  = -ansatz00 * test10;
      // SUPG term
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      // SUPG term
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRowC = MatrixC[i];

    test10 = Orig5[i];
    test01 = Orig6[i];
    test00 = Orig7[i];
    // rhs
    Rhs3[i] += Mult*delta*(c1*test10+c2*test01);

    // pressure-pressure block
    for(j=0;j<N_P;j++)
    {
	ansatz10 = Orig5[j];
	ansatz01 = Orig6[j];
	
	val = delta * (ansatz10*test10+ansatz01*test01);
	MatrixRowC[j] += Mult*val;
    }

    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      // divergence constraint
      val = -test00*ansatz10;
      val -=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * delta * test10;
      MatrixRow1[j] -= Mult*val;

      val = -test00*ansatz01;
      val -=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * delta * test01;
      MatrixRow2[j] -= Mult*val;
    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEMEquOrd(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixA12, **MatrixA21;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *Matrix12Row, *Matrix21Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, tau, ugrad, r=2, maxu;

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

  // DIE PARAMETER AUS DER LITERATUR NOCH EINSETZEN
  delta = r*r*r*r*c0/(hK*hK) + r*maxu/hK;
  delta =  delta0/delta;
  
  tau = r*r*c0 + hK*maxu/r;
  tau *= delta1;

  // assembling for velocity test functions
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
    // rhs
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    // velocity-velocity block
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

    // pressure-velocity block
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      // pressure ansatz functions
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      ansatz00 = Orig7[j];
     
      // pressure term 
      val  = -ansatz00 * test10;
      // SUPG term
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test10 = Orig5[i];
    test01 = Orig6[i];
    test00 = Orig7[i];

    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      // divergence constraint
      val = -test00*ansatz10;
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * delta * test10;
      MatrixRow1[j] -= Mult*val;

      val = -test00*ansatz01;
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * delta * test01;
      MatrixRow2[j] -= Mult*val;
    }                            // endfor j
  }
}

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void NSParamsVelo_Press(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // px
  out[3] = in[5];                // py
}
