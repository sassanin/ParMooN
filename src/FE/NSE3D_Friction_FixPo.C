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
// @(#)NSE3D_FixPo.C        1.3 06/27/00
//
// common declaration for all Navier-Stokes problems
// fix point iteration
//
// ======================================================================

#include <Database.h>
#include <MooNMD_Io.h>
#include <TNSE3D_Routines.h>
#include <Convolution.h>

// ======================================================================
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================
void NSType3GalerkinDD3DFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, c4, c5;
  double u1, u2, u3, u;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixB1  = LocMatrices[9];
  MatrixB2  = LocMatrices[10];
  MatrixB3  = LocMatrices[11];

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
  c1 = coeff[1]; // rhs f1
  c2 = coeff[2]; // rhs f2
  c3 = coeff[3]; // rhs f3
  c4 = coeff[4]; // friction f1
  c5 = coeff[5]; // friction f1
  

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  u  = sqrt(u1*u1+u2*u2+u3*u3);

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

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz000 = Orig3[j];
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += c4*test000*ansatz000;
      val += c5*u*test000*ansatz000;

      Matrix11Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += c4*test000*ansatz000;
      val += c5*u*test000*ansatz000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += c4*test000*ansatz000;
      val += c5*u*test000*ansatz000;
      Matrix33Row[j] += Mult * val;
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
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================

void NSType3_4NLGalerkinDD3DFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22,  **MatrixA33;
  double val, val1, val2, val3, val4;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0, c4;
  double u1, u2, u3, u;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

//  c0 = 2*coeff[0]; // 2*nu
  c0 = coeff[0]; // 2*nu
  c4 = coeff[4]; // friction

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  u  = sqrt(u1*u1+u2*u2);

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
      ansatz000 = Orig3[j];
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      /*  val  = c0*(test100*ansatz100+0.5*test010*ansatz010+
        0.5*test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(0.5*test100*ansatz100+test010*ansatz010+
        0.5*test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(0.5*test100*ansatz100+0.5*test010*ansatz010+
        test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;*/

      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;

      val  = c0*(2*val2+val3+val4)+val1;
      val += c4*u*ansatz000*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(val2+2*val3+val4)+val1;
      val += c4*u*ansatz000*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(val2+val3+2*val4)+val1;
      val += c4*u*ansatz000*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}
