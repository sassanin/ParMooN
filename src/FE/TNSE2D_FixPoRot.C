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
// TNSE2D_FixPoRot.C     06/03/20
//
// common declaration for all time dependent Navier-Stokes problems
// rotation form of nonlinear term
// ======================================================================

#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>
#include <TNSE2D_Routines.h>

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// disc type is changed internally to Smagorinsky without additional term
// ======================================================================

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3SmagorinskyRot(double Mult, double *coeff,
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
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

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

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = -u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = -u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00; 
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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
void TimeNSType3SmagorinskyRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, viscosity, delta;

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
  viscosity = c0+mu;

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

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4SmagorinskyRot(double Mult, double *coeff,
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
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

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

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val = -u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val = -u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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
void TimeNSType4SmagorinskyRotDD(double Mult, double *coeff,
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
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

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

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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
// Assembling routine for all nonlinear matrices
// ======================================================================
// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, viscosity, delta;
  double u1, u2, mu;
  //cout << "Sma" << endl;
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
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);

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
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = mu*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = mu*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, viscosity, delta;
  double u1, u2, mu;
  // cout << "Sma" << endl;
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
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);
  mu = mu/2.0;
  viscosity = mu+c0;

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
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// NEU! NEU! NEU! NEU! NEU! NEU! NEU! NEU! NEU! NEU! NEU! NEU! NEU! NEU!
// ======================================================================
// ======================================================================
// Type 3, Standard Galerkin (grad u, grad v)
// disc type is changed internally to Smagorinsky without additional term
// ======================================================================

// ======================================================================
// Type 3, Galerkin (grad u, grad v)
// ======================================================================
void TimeNSType3GalerkinRot(double Mult, double *coeff,
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
  int i,j, N_U, N_P;
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

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old

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
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = -u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = -u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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


// ======================================================================
// Type 3, Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
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

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old

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
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0)*(0.5*test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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


// ======================================================================
// Type 4, Galerkin, (grad u, grad v)
// ======================================================================
void TimeNSType4GalerkinRot(double Mult, double *coeff,
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
  int i,j, N_U, N_P;
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

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old

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
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val = -u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val = -u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = (c0)*(test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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
// Type 4, Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinRotDD(double Mult, double *coeff,
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
  int i,j, N_U, N_P;
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

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // p

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old

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
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
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
// Assembling routine for all nonlinear matrices
// ======================================================================
// ======================================================================
// Type 3, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u

  c0 = coeff[0];                                  // nu

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old

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

      val  = (c0)*(test10*ansatz10+test01*ansatz01);
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      //      val  = mu*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      //      val  = mu*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = (c0)*(test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

    }                                             // endfor j
  }                                               // endfor i
}


// ======================================================================
// Type 3, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u

  c0 = coeff[0];                                  // nu

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old

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

      val  = 2*(c0)*(test10*ansatz10+0.5*test01*ansatz01);
      val += u2*ansatz01*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0)*(test01*ansatz10);
      val -= u2 * ansatz10*test00;
      Matrix12Row[j] += Mult * val;

      val  = (c0)*(test10*ansatz01);
      val -= u1 * ansatz01 * test00;
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0)*(0.5*test10*ansatz10+test01*ansatz01);
      val += u1*ansatz10*test00;
      Matrix22Row[j] += Mult * val;

    }                                             // endfor j
  }                                               // endfor i
}


// ======================================================================
// Type 1, Standard Galerkin + Div Term
// ======================================================================
void TimeNSType1GalerkinDiv(double Mult, double *coeff,
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
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;
 
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

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
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
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


// ======================================================================
// Type 2, Standard Galerkin + div term
// ======================================================================
void TimeNSType2GalerkinDiv(double Mult, double *coeff,
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
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;

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
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
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
// Type 1, Standard Galerkin + div term, only nonlinear part
// Type 2, Standard Galerkin + div term, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkinDiv(double Mult, double *coeff,
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

double ComputePSPGParamter(double hK, double nu, double alpha)
{
  double delta1 = TDatabase::ParamDB->DELTA1, val, val1;
  double c_inv = 1.0;
  
  //val = delta1*hK*hK/(4*c_inv*nu);
  val = TDatabase::ParamDB->DELTA0 * hK * hK/nu;
  if (alpha > 0)
  {
    val1 = 1/(4.*alpha);
    if (val1 < val)
      val = val1;
  }
  return(val);
}

// ======================================================================
// Type 14, PSPG for transient Stokes
// ======================================================================
void TStokes_PSPG_GRADDIV(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22, **MatrixC;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixRowC;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3, delta, u1, u2, mu;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double implicit_version = TDatabase::ParamDB->P9;
  
  if (implicit_version)
    implicit_version = 1;
  else
    implicit_version = 0;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixC = LocMatrices[6];  
  MatrixB1 = LocMatrices[7];
  MatrixB2 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // u_xx
  Orig4 = OrigValues[4];                          // u_yy
  Orig5 = OrigValues[5];                          // p_x
  Orig6 = OrigValues[6];                          // p_y
  Orig7 = OrigValues[7];                          // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 =  TDatabase::ParamDB->OSEEN_ZERO_ORDER_COEFF; 
  delta = ComputePSPGParamter(hK, c0, c3);
  // grad-div parameter
  mu = TDatabase::ParamDB->DIV_DIV_STAB_C1 * pow(hK,TDatabase::ParamDB->DIV_DIV_STAB_C2);
  
  u1 = param[0];                                  // u1old (previous time)
  u2 = param[1];                                  // u2old (previous time)

// velocity test functions
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

    // velocity basis functions
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      val += mu*test10*ansatz10;
      Matrix11Row[j] += Mult * val;

      val  = mu * test10*ansatz01;
      Matrix12Row[j] += Mult * val;

      val  = mu * test01*ansatz10;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      val += mu*test01*ansatz01;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    // pressure basis functions
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig7[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
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
    // contributions that are needed in the explicit and implicit treatment
    // contribution from right-hand side
    val = c1*test10+c2*test01;
    // contribution from time derivative
    val += (u1*test10+u2*test01)/tau;
    val *= delta;
    // this value will be overloaded later
    Rhs3[i] += Mult*val;

    // pressure-pressure block
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig5[j];
      ansatz01 = Orig6[j];
      // contribution from right-hand side
      val = delta * (ansatz10*test10+ansatz01*test01);      
      MatrixRowC[j] += Mult*val;
    }

    // velocity-pressure block
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
      //OutPut(ansatz20 << " "  << ansatz02 << ":");
      // divergence constraint
      val = -test00*ansatz10;
      // stabilization 
      val -=  (ansatz00/tau-c0*(ansatz20+ansatz02) + c3*ansatz00 ) * delta * test10 * (1-implicit_version);
      MatrixRow1[j] -= Mult*val;

      val = -test00*ansatz01;
      val -=  (ansatz00/tau-c0*(ansatz20+ansatz02) + c3*ansatz00 ) * delta * test01 * (1-implicit_version);
      MatrixRow2[j] -= Mult*val;
    }                                             // endfor j
  }                                               // endfor i
}
 // ======================================================================
// Type 14, PSPG for transient Stokes, nonlinear loop 
// ======================================================================
void TStokes_PSPG_GRADDIV_NL(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22, **MatrixC;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3, u1, u2, delta, mu;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double implicit_version = TDatabase::ParamDB->P9;
  
  if (implicit_version)
    implicit_version = 1;
  else
    implicit_version = 0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // u_xx
  Orig4 = OrigValues[4];                          // u_yy
  Orig5 = OrigValues[5];                          // p_x
  Orig6 = OrigValues[6];                          // p_y
  Orig7 = OrigValues[7];                          // p
  Rhs3 = LocRhs[2];

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 =  TDatabase::ParamDB->OSEEN_ZERO_ORDER_COEFF; 
  u1 = param[0];                                  // u1old (previous time)
  u2 = param[1];                                  // u2old (previous time)
  delta = ComputePSPGParamter(hK, c0, c3);
  // grad-div parameter
  mu = TDatabase::ParamDB->DIV_DIV_STAB_C1 * pow(hK,TDatabase::ParamDB->DIV_DIV_STAB_C2);

  // velocity test functions
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // velocity basis functions
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      val += mu*test10*ansatz10;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //val += (u1*ansatz10+u2*ansatz01)*test00;
      val += c3*ansatz00*test00;
      val += mu*test01*ansatz01;
      Matrix22Row[j] += Mult * val;
    }                            // endfor j
  }
   // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    test10 = Orig5[i];
    test01 = Orig6[i];
    test00 = Orig7[i];
    // contributions that change in time step for implicit version
    // contribution from time derivative
    val = -(u1*test10+u2*test01)/tau * implicit_version;
    // contribution from zeroth order term
    val -= c3*(u1*test10+u2*test01) * implicit_version; 
    val *= delta;
    Rhs3[i] += Mult*val;
  }                                               // endfor i
}
// ======================================================================
// Type 14, PSPG for transient Stokes
//  computes all entries that do not change in the time step (if 
//  implicit version is applied)
// ======================================================================
void TStokes_PSPG_GRADDIV_Rhs(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22, **MatrixC;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixRowC;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3, delta, u1, u2;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double implicit_version = TDatabase::ParamDB->P9;
  
  if (implicit_version)
    implicit_version = 1;
  else
    implicit_version = 0;
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // u_xx
  Orig4 = OrigValues[4];                          // u_yy
  Orig5 = OrigValues[5];                          // p_x
  Orig6 = OrigValues[6];                          // p_y
  Orig7 = OrigValues[7];                          // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 =  TDatabase::ParamDB->OSEEN_ZERO_ORDER_COEFF; 
  delta = ComputePSPGParamter(hK, c0, c3);
  
  u1 = param[0];                                  // u1old (previous time)
  u2 = param[1];                                  // u2old (previous time)

// velocity test functions
  for(i=0;i<N_U;i++)
  {
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
  }                              // endfor i

   // assembling for pressure test functions
  for(i=0;i<N_P;i++)
  {
    test10 = Orig5[i];
    test01 = Orig6[i];
    test00 = Orig7[i];
    // contributions that does not change in time step 
    // contribution from right-hand side
    val = c1*test10+c2*test01;
    // contribution from time derivative, old time 
    val += (u1*test10+u2*test01)/tau;
    //val += (coeff[3]*test10+coeff[4]*test01);
    val *= delta;
    Rhs3[i] += Mult*val;
  }                       // endfor i
}
