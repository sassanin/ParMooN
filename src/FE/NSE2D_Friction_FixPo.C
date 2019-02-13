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
#include <Output2D.h>
#include <stdlib.h>

// ======================================================================
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================
void NSType3GalerkinDDFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3, c4;
  double u1, u2, u;

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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // rhs f1
  c2 = coeff[2]; // rhs f2
  c3 = coeff[3]; // friction f1
  c4 = coeff[4]; // friction f2
  //cout << "Friction  " << c3 << endl;

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  //cout << "Friction u2 " << u2 << endl;
  u = sqrt(u1*u1+u2*u2); // |u|

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
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*test00*ansatz00;
      val += c4*u*test00*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*test00*ansatz00;
      val += c4*u*test00*ansatz00;
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i

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
    } // endfor j

  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================
void NSType3GalerkinDDFrictionLocal(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u;
  double Re, eps;

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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction
  eps=0.57; // porosity

  //cout << "Friction  " << c3 << endl;

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  //cout << "Friction u2 " << u2 << endl;
  u = sqrt(u1*u1+u2*u2); // |u|
 
  
  Re=0.18/c0;
  c3 = 6.8*pow(1-eps,1.2)/(eps*eps*eps)*pow(Re,-0.2)*pow(u,0.8)/0.18; // friction
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
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*test00*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*test00*ansatz00;
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i

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
    } // endfor j

  } // endfor i
}


// ======================================================================
// Type 3, with friction, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3UpwindFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u;

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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u = sqrt(u1*u1+u2*u2); // |u|


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
      val += c3*u*ansatz00*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += c3*u*ansatz00*test00;
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i

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
    } // endfor j

  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction (grad u, grad v)
// ======================================================================
void NSType4SDFEMFriction(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction
    
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
   u = sqrt(u1*u1+u2*u2); // |u|
   //cout << "|u|  " << u << endl;
   
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
      ansatz00 = Orig2[j];
      
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
             + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) )*ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
             + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) )*ugrad; // SD term
     Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i

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
    } // endfor j

  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction (grad u, grad v)
// ======================================================================
void NSType4SDFEMFrictionRST(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]/c0; // friction
    
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
   u = sqrt(u1*u1+u2*u2); // |u|
   //cout << "|u|  " << u << endl;
   
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

    ugrad  = delta/c0 * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      ansatz00 = Orig2[j];
      
      val  = test10*ansatz10+test01*ansatz01;
      val += ((u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00)*test00;
      val += (-ansatz20-ansatz02
             + (u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00)*ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = test10*ansatz10+test01*ansatz01;
      val += ((u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00)*test00;
      val += (-ansatz20-ansatz02
             + (u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00)*ugrad; // SD term
     Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i

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
    } // endfor j

  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction, D(u):D(v)
// ======================================================================
void NSType4SDFEMDDFriction(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction

  u1 = param[0]; // u1old
  //cout << "u1old:  " << u1 << endl;
  u2 = param[1]; // u2old
  //cout << "u2old:  " << u2 << endl;
   u = sqrt(u1*u1+u2*u2);

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
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
              + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
              + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
     Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i

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
    } // endfor j

  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction on rhs, D(u):D(v)
// ======================================================================
void NSType4SDFEMDDFrictionRhs(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
   u = sqrt(u1*u1+u2*u2);

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

    Rhs1[i] += Mult*((test00+ugrad)*(c1-c3*u*u1));
    Rhs2[i] += Mult*((test00+ugrad)*(c2-c3*u*u2));
    

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
              + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
              + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
     Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i

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
    } // endfor j

  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDDFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, c3, c4;
  double u1, u2, u;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c3 = coeff[3]; // friction f1
  c4 = coeff[4]; // friction f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u  = sqrt(u1*u1+u2*u2);
  
  //cout << "u1  " << u1 << endl;
  //  cout << "Friction NL  " << c5 << endl;
  //ut << "Friction NL u2  " << u2 << endl;
  

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
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      val += c4*u*ansatz00*test00;
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      val += c4*u*ansatz00*test00;
      Matrix22Row[j] += Mult * val;


    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDDFrictionLocal(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, c3;
  double u1, u2, u;
  double Re, eps;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c3 = coeff[3]; // friction
  eps=0.57; // porosity

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u  = sqrt(u1*u1+u2*u2);
  //cout << "u1  " << u1 << endl;
  //  cout << "Friction NL  " << c5 << endl;
  //ut << "Friction NL u2  " << u2 << endl;
  Re=0.18/c0;
  c3 = 6.8*pow(1-eps,1.2)/(eps*eps*eps)*pow(Re,-0.2)*pow(u,0.8)/0.18; // friction

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
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      Matrix22Row[j] += Mult * val;


    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3_4NLUpwindFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01, ansatz00;
  double test10, test01, test00;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, c3;
  double u1, u2, u;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c3 = coeff[3]; // friction
  
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u = sqrt(u1*u1+u2*u2); // |u|

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
      val += c3*u*test00*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += c3*u*test00*ansatz00;
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEMFriction(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u  = sqrt(u1*u1+u2*u2); // |u|
  //cout << "|u| NL  " << u << endl;

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
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
            + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
            + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction, (grad u, grad v), Lutz
// ======================================================================
void NSType4NLSDFEMFrictionRST(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]/c0; // friction

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u  = sqrt(u1*u1+u2*u2); // |u|
  //cout << "|u| NL  " << u << endl;

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

    ugrad  = delta/c0 * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      
      val  = test10*ansatz10+test01*ansatz01;
      val += ((u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00)*test00;
      val += (-ansatz20-ansatz02
            + (u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = test10*ansatz10+test01*ansatz01;
      val += ((u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00)*test00;
      val += (-ansatz20-ansatz02
            + (u1*ansatz10+u2*ansatz01)/c0+c3*u*ansatz00) * ugrad; // SD term
      Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDDFriction(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
   u = sqrt(u1*u1+u2*u2);

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
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      // val += (-c0*(ansatz20+ansatz02)
      //	      + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      // val += (-c0*(ansatz20+ansatz02)
      //      + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, with friction on rhs, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDDFrictionRhs(double Mult, double *coeff, 
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
  double c0, c1, c2, c3;
  double u1, u2, u;
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy
  Orig5 = OrigValues[5]; // p_x
  Orig6 = OrigValues[6]; // p_y
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // friction

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
   u = sqrt(u1*u1+u2*u2);

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

    Rhs1[i] += Mult*((test00+ugrad)*(c1-c3*u*u1));
    Rhs2[i] += Mult*((test00+ugrad)*(c2-c3*u*u2));

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      
      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
            + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01+c3*u*ansatz00)*test00;
      val += (-c0*(ansatz20+ansatz02)
            + (u1*ansatz10+u2*ansatz01+c3*u*ansatz00) ) * ugrad; // SD term
      Matrix22Row[j] += Mult * val;

    } // endfor j

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
  } // endfor i
}
