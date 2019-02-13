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
// TNSE3D_Newton.C
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>
#include <TNSE3D_FixPo.h>
#include <TNSE3D_Routines.h>

#include <stdlib.h>

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// Type 3, ClassicalLES, (grad u, grad v)
// Type 3, GL00Convolution, (grad u, grad v)`
// ======================================================================
void TimeNSType3GalerkinNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, a000t000;
  double u1, u2, u3;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;

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
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

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
      a000t000 = Mult*ansatz000*test000;
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val+ a000t000*du1x;
      Matrix12Row[j] += a000t000*du1y;
      Matrix13Row[j] += a000t000*du1z;
      Matrix21Row[j] += a000t000*du2x;
      Matrix22Row[j] += val+ a000t000*du2y;
      Matrix23Row[j] += a000t000*du2z;
      Matrix31Row[j] += a000t000*du3x;
      Matrix32Row[j] += a000t000*du3y;
      Matrix33Row[j] += val+ a000t000*du3z;

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
void TimeNSType3GalerkinNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType3UpwindNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, a000t000;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;

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

  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

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
      a000t000 = Mult*ansatz000*test000;
      
      val  = Mult*c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      
      Matrix11Row[j] += val+ a000t000*du1x;
      Matrix12Row[j] += a000t000*du1y;
      Matrix13Row[j] += a000t000*du1z;
      Matrix21Row[j] += a000t000*du2x;
      Matrix22Row[j] += val+ a000t000*du2y;
      Matrix23Row[j] += a000t000*du2z;
      Matrix31Row[j] += a000t000*du3x;
      Matrix32Row[j] += a000t000*du3y;
      Matrix33Row[j] += val+ a000t000*du3z;

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
void TimeNSType3UpwindNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType3SmagorinskyNewton3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
                            param[21]);

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
void TimeNSType3SmagorinskyNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType3GL00AuxProblemNewton3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
                            param[21]);

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
void TimeNSType3GL00AuxProblemNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType4GalerkinNewton3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType4GalerkinNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType4UpwindNewton3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType4UpwindNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
void TimeNSType4SmagorinskyNewton3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
                            param[21]);

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
void TimeNSType4SmagorinskyNewtonDD3D(double Mult, double *coeff, 
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
  double c0, c1, c2, c3, delta;
  double u1, u2, u3, mu, viscosity, a000t000;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;


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
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

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
      a000t000 = ansatz000*test000;

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
     
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1+a000t000*du1x;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      val += a000t000*du1y;
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      val += a000t000*du1z;
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      val += a000t000*du2x;
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1+ a000t000*du2y;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      val += a000t000*du2z;
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      val += a000t000*du3x;
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      val += a000t000*du3y;
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1+ a000t000*du3z;
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
void TimeNSType4GL00AuxProblemNewton3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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
                            param[21]);

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
void TimeNSType4GL00AuxProblemNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);


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

void TimeNSType4VMS_ProjectionNewtonDD3D(double Mult, double *coeff, 
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
  double u1, u2, u3, mu, viscosity, a000t000;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;

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
  N_L = N_BaseFuncts[3];

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
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

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
      a000t000 = ansatz000*test000;

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1+a000t000*du1x;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      val += a000t000*du1y;
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      val += a000t000*du1z;
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      val += a000t000*du2x;
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1+ a000t000*du2y;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      val += a000t000*du2z;
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      val += a000t000*du3x;
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      val += a000t000*du3y;
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1+ a000t000*du3z;
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
// Assembling routine for all nonlinear matrices
// ======================================================================

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkinNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3, a000t000;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;
 
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
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

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
      
      a000t000 = Mult*ansatz000*test000;
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;

      Matrix11Row[j] += val+ a000t000*du1x;
      Matrix12Row[j] += a000t000*du1y;
      Matrix13Row[j] += a000t000*du1z;
      Matrix21Row[j] += a000t000*du2x;
      Matrix22Row[j] += val+ a000t000*du2y;
      Matrix23Row[j] += a000t000*du2z;
      Matrix31Row[j] += a000t000*du3x;
      Matrix32Row[j] += a000t000*du3y;
      Matrix33Row[j] += val+ a000t000*du3z;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);

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
void TimeNSType3_4NLUpwindNewton3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0, a000t000;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;
 
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

  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

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
      
      a000t000 = Mult*ansatz000*test000;
      val  = Mult*c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      
      Matrix11Row[j] += val+ a000t000*du1x;
      Matrix12Row[j] += a000t000*du1y;
      Matrix13Row[j] += a000t000*du1z;
      Matrix21Row[j] += a000t000*du2x;
      Matrix22Row[j] += val+ a000t000*du2y;
      Matrix23Row[j] += a000t000*du2z;
      Matrix31Row[j] += a000t000*du3x;
      Matrix32Row[j] += a000t000*du3y;
      Matrix33Row[j] += val+ a000t000*du3z;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// Type 4, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3_4NLUpwindNewtonDD3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);



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
void TimeNSType3_4NLSmagorinskyNewton3D(double Mult, double *coeff, 
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
  OutPut("Discretization not yet implemented !!!"<<endl);
  exit(4711);
 

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
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[15],&param[12],&param[13],&param[14],
                            param[21]);
 
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
void TimeNSType3_4NLSmagorinskyNewtonDD3D(double Mult, double *coeff, 
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
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, viscosity, delta;
  double u1, u2, u3, mu, a000t000;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;

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
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[15],&param[12],&param[13],&param[14],
                            param[21]);
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
      ansatz000 = Orig3[j];
      a000t000 = ansatz000*test000;
      
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
      Matrix11Row[j] += val2+val1+a000t000*du1x;
      Matrix12Row[j] += test010*ansatz100+a000t000*du1y; 
      Matrix13Row[j] += test001*ansatz100+a000t000*du1z;
      Matrix21Row[j] += test100*ansatz010+a000t000*du2x;
      Matrix22Row[j] += val3+val1+a000t000*du2y;
      Matrix23Row[j] += test001*ansatz010+a000t000*du2z;
      Matrix31Row[j] += test100*ansatz001+a000t000*du3x;
      Matrix32Row[j] += test010*ansatz001+a000t000*du3y;
      Matrix33Row[j] += val4+val1+a000t000*du3z;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v)
// Type 4, VMS_Projection, D(u):D(v)
// ======================================================================
void TimeNSType3_4NLVMS_ProjectionNewtonDD3D(double Mult, double *coeff, 
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
  double c0, viscosity, delta;
  double u1, u2, u3, mu, a000t000;
  double du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;


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
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[15],&param[12],&param[13],&param[14],
                            param[21]);
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
      ansatz000 = Orig3[j];
      a000t000 = ansatz000*test000;

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1+a000t000*du1x;
      Matrix12Row[j] += test010*ansatz100+a000t000*du1y; 
      Matrix13Row[j] += test001*ansatz100+a000t000*du1z;
      Matrix21Row[j] += test100*ansatz010+a000t000*du2x;
      Matrix22Row[j] += val3+val1+a000t000*du2y;
      Matrix23Row[j] += test001*ansatz010+a000t000*du2z;
      Matrix31Row[j] += test100*ansatz001+a000t000*du3x;
      Matrix32Row[j] += test010*ansatz001+a000t000*du3y;
      Matrix33Row[j] += val4+val1+a000t000*du3z;
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
        val1 = Mult * mu * ansatz000;
        Matrix11Row[j] -= val1 * test100;
        Matrix22Row[j] -= val1 * test010;
        Matrix33Row[j] -= val1 * test001;
     }
  }   
}
// ======================================================================
// Assembling routine for right-hand sides only
// ======================================================================

// ======================================================================
// right-hand side for NSE ONLY
// ======================================================================
void TimeNSRHSNewton3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6;
  double test000;
  double *Orig0;
  int i, N_U;
  double c1, c2, c3;
  double u1, u2, u3, du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;
 
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

  for(i=0;i<N_U;i++)
  {
    test000 = Mult*Orig0[i];

    Rhs1[i] += test000*c1;
    Rhs2[i] += test000*c2;
    Rhs3[i] += test000*c3;
    Rhs4[i] += test000*(u1*du1x+u2*du1y+u3*du1z);
    Rhs5[i] += test000*(u1*du2x+u2*du2y+u3*du2z);
    Rhs6[i] += test000*(u1*du3x+u2*du3y+u3*du3z);
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " "; 
  } // endfor i
}
void TimeNSRHSNewtonNL3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double *Rhs4, *Rhs5, *Rhs6;
  double test000;
  double *Orig0;
  int i, N_U;
  double c1, c2, c3;
  double u1, u2, u3, du1x, du1y, du1z, du2x, du2y, du2z, du3x, du3y, du3z;

  Rhs4 = LocRhs[0];
  Rhs5 = LocRhs[1];
  Rhs6 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  du1x = param[3]; // du1old/dx
  du2x = param[4]; // du2old/dx
  du3x = param[5]; // du3old/dx
  du1y = param[6]; // du1old/dy
  du2y = param[7]; // du2old/dy
  du3y = param[8]; // du3old/dy
  du1z = param[9]; // du1old/dz
  du2z = param[10]; // du2old/dz
  du3z = param[11]; // du3old/dz

  for(i=0;i<N_U;i++)
  {
    test000 = Mult*Orig0[i];

    Rhs4[i] += test000*(u1*du1x+u2*du1y+u3*du1z);
    Rhs5[i] += test000*(u1*du2x+u2*du2y+u3*du2z);
    Rhs6[i] += test000*(u1*du3x+u2*du3y+u3*du3z);
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " "; 
  } // endfor i
}

// ======================================================================
// right-hand side for NSE ONLY 
// ClassicalLES model
// ======================================================================
void TimeNSRHSClassicalLESNewton3D(double Mult, double *coeff,
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
void TimeNSRHSLESModelNewton3D(double Mult, double *coeff,
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
void TimeNSGL00AuxProblemRHSNewton3D(double Mult, double *coeff,
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
// right-hand side for additional terms in rhs of small scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_SmallRhsNewton3D(double Mult, double *coeff, 
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
