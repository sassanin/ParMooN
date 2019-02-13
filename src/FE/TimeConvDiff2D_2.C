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
// @(#)TimeConvDiff2D_2.C        1.0 07/25/00
//
// common declaration for time dependent convection diffusion problems
// with two equations
// ======================================================================

#include <Database.h>

void TimeBilinearAssemble_U(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix_U, **Matrix_V, *Rhs, val, *MatrixRow_U, *MatrixRow_V;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8; 

  Matrix_U = LocMatrices[0];
  Matrix_V = LocMatrices[1];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // Laplace
  c1 = coeff[1]; // u_x
  c2 = coeff[2]; // u_y
  c3 = coeff[3]; // u
  c4 = coeff[4]; // Laplace
  c5 = coeff[5]; // u_x
  c6 = coeff[6]; // u_y
  c7 = coeff[7]; // u
  c8 = coeff[8]; // rhs

  for(i=0;i<N_;i++)
  {
    MatrixRow_U = Matrix_U[i];
    MatrixRow_V = Matrix_V[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c8;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val *=Mult;

      MatrixRow_U[j] += val;
      
      val = c4*(test10*ansatz10+test01*ansatz01);
      val += (c5*ansatz10+c6*ansatz01)*test00;
      val += c7*ansatz00*test00;

      val *=Mult;

      MatrixRow_V[j] += val;
    } // endfor j
  } // endfor i
}

void TimeBilinearAssemble_V(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix_U, **Matrix_V, *Rhs, val, *MatrixRow_U, *MatrixRow_V;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8; 

  Matrix_U = LocMatrices[0];
  Matrix_V = LocMatrices[1];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // Laplace
  c1 = coeff[1]; // u_x
  c2 = coeff[2]; // u_y
  c3 = coeff[3]; // u
  c4 = coeff[4]; // Laplace
  c5 = coeff[5]; // u_x
  c6 = coeff[6]; // u_y
  c7 = coeff[7]; // u
  c8 = coeff[8]; // rhs

  for(i=0;i<N_;i++)
  {
    MatrixRow_U = Matrix_U[i];
    MatrixRow_V = Matrix_V[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c8;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val *=Mult;

      MatrixRow_U[j] += val;
      
      val = c4*(test10*ansatz10+test01*ansatz01);
      val += (c5*ansatz10+c6*ansatz01)*test00;
      val += c7*ansatz00*test00;

      val *=Mult;

      MatrixRow_V[j] += val;
    } // endfor j
  } // endfor i
}

void TimeBilinearAssemble_SD_U(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4, c5; 
  double delta, bgradv;
  static double delta0 = TDatabase::ParamDB->DELTA0;

  double bound = TDatabase::ParamDB->P8;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];
  c1 = coeff[1];
  c2 = coeff[2];
  c3 = coeff[3];
  c4 = coeff[4];
  c5 = coeff[5];

/*
  // normal delta
  if(c0 < hK*c5)
  {
    delta = delta0 * hK/c5;
  }
  else
    delta = 0;
*/

/*
  // delta for SDFEM only in coarse part of Shishkin mesh
  if(hK>bound) 
  {
    delta = delta0 * hK;
    // cout << "SDFEM" << endl;
  }
  else
  {
    // cout << "NO" << endl;
    delta = 0;
  }
  // cout << delta << endl;
*/

// /*
  // delta everywhere
  delta = delta0 * hK;
// */

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    Rhs[i] += Mult*(test00+delta*bgradv)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val += delta * (-c0*(ansatz20+ansatz02)
                      +c1*ansatz10+c2*ansatz01
                      +c3*ansatz00) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    } // endfor j
  } // endfor i
}

/*
void TimeBilinearAssembleRB(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, **MatrixA,  *Rhs, val, *MatrixRow, *MatrixRowA;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4,c5; 

  Matrix = LocMatrices[0];
  MatrixA = LocMatrices[1];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // Q
  c5 = coeff[5]; // dQ/du

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      MatrixRowA[j] += Mult*val;

      val += c5*ansatz00*test00;
      val *=Mult;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

void TimeBilinearAssembleRB1(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4,c5; 

  Matrix = LocMatrices[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // Q
  c5 = coeff[5]; // dQ/du

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      // MatrixRowA[j] += Mult*val;

      val += c5*ansatz00*test00;
      val *=Mult;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}
*/

void TimeRhsAssembleRB_U(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c8; 

  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c8 = coeff[8]; // Q

  for(i=0;i<N_;i++)
  {
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c8;

  } // endfor i
}

void TimeRhsAssembleRB_V(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c8; 

  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c8 = coeff[8]; // Q

  for(i=0;i<N_;i++)
  {
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c8;

  } // endfor i
}
