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
// @(#)MovingNavierStokes.C        1.2 04/13/00
//
// Assembling routine for all matrices and right-hand sides
// ======================================================================

#include <Database.h>

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void MovingNSType1Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
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
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
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
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void MovingNSType1Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

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

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
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
      
      val = c0*(test10*ansatz10+test01*ansatz01);
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
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
// Type 2, Standard Galerkin
// ======================================================================
void MovingNSType2Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double g1, g2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
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
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    } // endfor j

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
// for Type 2 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void MovingNSType2Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
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
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    } // endfor j

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
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void MovingNSType3Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
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
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void MovingNSType3GalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
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
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void MovingNSType3Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
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
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void MovingNSType3UpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
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
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void MovingNSType4Galerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
    } // endfor j

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
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void MovingNSType4GalerkinDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
    } // endfor j

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
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void MovingNSType4Upwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
    } // endfor j

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
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void MovingNSType4UpwindDD(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
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
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM11Row[j] += Mult * val;

      val = ansatz00*test00;
      MatrixM22Row[j] += Mult * val;
    } // endfor j

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
// Assembling routine for all nonlinear matrices
// ======================================================================

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void MovingNSType1NLGalerkin(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void MovingNSType1NLUpwind(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

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
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, Standard Galerkin, only nonlinear part
// ======================================================================
void MovingNSType2NLGalerkin(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void MovingNSType2NLUpwind(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void MovingNSType3NLGalerkin(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void MovingNSType3NLGalerkinDD(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void MovingNSType3NLUpwind(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void MovingNSType3NLUpwindDD(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void MovingNSType4NLGalerkin(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear part
// ======================================================================
void MovingNSType4NLGalerkinDD(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old - vx
  u2 = param[1]; // u2old - vy

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

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void MovingNSType4NLUpwind(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void MovingNSType4NLUpwindDD(double Mult, double *coeff, 
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
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Assembling routine for right-hand sides only
// ======================================================================

// ======================================================================
// right-hand side ONLY
// ======================================================================
void MovingNSRHS(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
  } // endfor i
}

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void MovingNSParams2(double *in, double *out)
{
  out[0] = in[2]; // u1old
  out[1] = in[3]; // u2old
}

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void GridParams4(double *in, double *out)
{
  out[0] = in[2] - in[4];
  out[1] = in[3] - in[5];
}

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void GridAssemble(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y

  c0 = coeff[0];

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      
      val  = c0*(test10*ansatz10+test01*ansatz01);
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// // kind of boundary condition (for FE space needed)
// void GridBoundCondition(int BdComp, double t, BoundCond &cond)
// {
//   cond = DIRICHLET;
// }
// 
// // value of boundary condition
// void GridBoundValue(int BdComp, double Param, double &value)
// {
//   value = 0;
// }
// 
// void GridCoeffs(int n_points, double *x, double *y,
//         double **parameters, double **coeffs)
// {
//   int i;
//   double *coeff;
//   double r2;
// 
//   for(i=0;i<n_points;i++)
//   {
//     coeff = coeffs[i];
// 
//     // if( (r2 = (x[i]-0.2)*(x[i]-0.2) + (y[i]-0.2)*(y[i]-0.2)) < 0.01)
//     //   coeff[0] = 10*sqrt(r2);
//     // else
//     //   coeff[0] = 10*0.1;
//     coeff[0] = 1;
// 
//     coeff[1] = 0;
//     coeff[2] = 0;
//     coeff[3] = 0;
// 
//     coeff[4] = 0;
//   }
// }
