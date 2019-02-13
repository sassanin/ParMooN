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
//
//
// common declaration for all Navier-Stokes problems in polar coordinates
// fix point iteration
//
// ======================================================================

#include <Convolution.h>
#include <Database.h>
#include <TNSE2D_Routines.h>


// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType4GalerkinAxialSymm3D(double Mult, double *coeff, 
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
  int i,j,N_U, N_P,sign;
  double c0, c1, c2;
  double r, u1,u2,y;

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
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  
  y  = param[0]; // y
  u1 = param[1]; // u1old
  u2 = param[2]; // u2old
  r = fabs(y);

  if (y>=0)
     sign = 1;
  else
     sign = -1;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i]*sign;
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1*r;
    Rhs2[i] += Mult*test00*c2*r;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j]*sign;
      ansatz00 = Orig2[j];
      
      val  = c0*r*(test10*ansatz10+test01*ansatz01);
      val += r*(u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(r*(test10*ansatz10+test01*ansatz01)
		 -(ansatz10+ansatz01)*test00);
                 //+ansatz00*test01+ansatz01*test00+ansatz00*test00/r);
		 // old +ansatz00*test01-ansatz10*test00);
		 
      val += r*(u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10*r;
      MatrixRow1[j] += val;
      
      val = -Mult*ansatz00*(r*test01+test00);
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
      ansatz01 = Orig1[j]*sign;
      ansatz00 = Orig2[j];

      val = -Mult*test00*ansatz10*r;
      MatrixRow1[j] += val;

      val = -Mult*(r*ansatz01+ansatz00)*test00;
      MatrixRow2[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4GalerkinDDAxialSymm3D(double Mult, double *coeff, 
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
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEMAxialSymm3D(double Mult, double *coeff, 
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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
      
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
              + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
              + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
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
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SDFEMDDAxialSymm3D(double Mult, double *coeff, 
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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
              + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
              + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
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
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================

void NSType4UpwindAxialSymm3D(double Mult, double *coeff, 
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
  int i,j,N_U, N_P, sign;
  double c0, c1, c2;
  double r, y;

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
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  y  = param[0]; // y
  r = fabs(y);

  if (y>=0)
     sign = 1;
  else
     sign = -1;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i]*sign;
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1*r;
    Rhs2[i] += Mult*test00*c2*r;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j]*sign;
      ansatz00 = Orig2[j];

      val  = c0*r*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(r*(test10*ansatz10+test01*ansatz01)
		 -(ansatz10+ansatz01)*test00);
		 //+ansatz00*test01+ansatz01*test00+ansatz00*test00/r);
		 // old +ansatz00*test01-ansatz10*test00);
      Matrix22Row[j] += Mult * val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10*r;
      MatrixRow1[j] += val;
      
      val = -Mult*(ansatz00*test01*r+ansatz00*test00);
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
      ansatz01 = Orig1[j]*sign;
      ansatz00 = Orig2[j];

      val = -Mult*test00*ansatz10*r;
      MatrixRow1[j] += val;

      val = -Mult*(r*test00*ansatz01+ansatz00*test00);
      MatrixRow2[j] += val;
    } // endfor j

  } // endfor i
}

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void NSType4UpwindDDAxialSymm3D(double Mult, double *coeff, 
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

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
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4SmagorinskyAxialSymm3D(double Mult, double *coeff, 
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType4SmagorinskyDDAxialSymm3D(double Mult, double *coeff, 
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
void NSType3_4NLGalerkinAxialSymm3D(double Mult, double *coeff, 
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
  int i,j,N_U, sign;
  double c0;
  double r, u1, u2, y;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

  y  = param[0]; // y
  u1 = param[1]; // u1old
  u2 = param[2]; // u2old
  r = fabs(y);

  if (y>=0)
     sign = 1;
  else
     sign = -1;
  
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i]*sign;
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j]*sign;
      ansatz00 = Orig2[j];
      
      val  = c0*r*(test10*ansatz10+test01*ansatz01);
      val += r*(u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

       val  = c0*(r*(test10*ansatz10+test01*ansatz01)
		 -(ansatz10+ansatz01)*test00);
                  //+ansatz00*test01+ansatz01*test00+ansatz00*test00/r);
		  // old +ansatz00*test01-ansatz10*test00);
       val += r*(u1*ansatz10+u2*ansatz01)*test00;
       Matrix22Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDDAxialSymm3D(double Mult, double *coeff, 
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

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
// HERE
void NSType3_4NLUpwindAxialSymm3D(double Mult, double *coeff, 
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
  int i,j,N_U,sign;
  double c0;
  double r, y;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

  y  = param[0]; // y
  r = fabs(y);

  if (y>=0)
     sign = 1;
  else
     sign = -1;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i]*sign;
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j]*sign;
      ansatz00 = Orig2[j];

      val  = c0*r*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(r*(test10*ansatz10+test01*ansatz01)
		 -(ansatz10+ansatz01)*test00);
		 //+ansatz00*test01+ansatz01*test00+ansatz00*test00/r);
		  // old ++ansatz00*test01-ansatz10*test00);
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3_4NLUpwindDDAxialSymm3D(double Mult, double *coeff, 
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y

  c0 = coeff[0]; // nu


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

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3_4NLSmagorinskyAxialSymm3D(double Mult, double *coeff, 
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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

      val  = mu*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = mu*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3_4NLSmagorinskyDDAxialSymm3D(double Mult, double *coeff, 
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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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

    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEMAxialSymm3D(double Mult, double *coeff, 
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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
      
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
            + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
            + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
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
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDDAxialSymm3D(double Mult, double *coeff, 
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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

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
            + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
            + (u1*ansatz10+u2*ansatz01) ) * ugrad; // SD term
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
// DUESE: Type 3 and 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkinAxialSymm3D_Duese(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P,sign;
  double c0, c1, c2;
  double r, u1,u2,y,mu,delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  
  y  = param[0]; // y
  u1 = param[1]; // u1old
  u2 = param[2]; // u2old
  r = fabs(y);

  if (y>=0)
     sign = 1;
  else
     sign = -1;

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[3],&param[1],&param[1]);
  mu = 0;
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i]*sign;
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j]*sign;
      ansatz00 = Orig2[j];
      
      val  = (c0+mu)*r*(test10*ansatz10+test01*ansatz01);
      val += r*(u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(r*(test10*ansatz10+test01*ansatz01)
     		 -(ansatz10+ansatz01)*test00);
      val += r*(u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10*r;
      MatrixRow1[j] += val;
      
      val = -Mult*ansatz00*(r*test01+test00);
      MatrixRow2[j] += val;
    }
  } // endfor i
}

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: y, u1old, u2old
// ========================================================================
void NSParamsVeloAxialSymm3D(double *in, double *out)
{
  out[0] = in[1]; // y
  out[1] = in[2]; // u1old
  out[2] = in[3]; // u2old  
}

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVeloAxialSymm3D(double *in, double *out)
{
  out[0] = in[1]; // y
  out[1] = in[2]; // u1old
  out[2] = in[3]; // u2old
  out[3] = in[4]; // D10(u1old)
  out[4] = in[5]; // D10(u2old)
  out[5] = in[6]; // D01(u1old)
  out[6] = in[7]; // D01(u2old)
}
