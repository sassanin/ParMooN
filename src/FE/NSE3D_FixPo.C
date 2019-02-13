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
// compute parameter for RFB stabilization
// Brezzi, Marini, Russo, CMAME 194 (2005) 127 - 148
// not clear if 2d situation can be simply used also in 3d
// ======================================================================

double RFB_Parameter3D(double hK, double eps, double* b)
{
    double peclet, normb,tau;

    normb = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
    if (normb==0)
	return(0.0);
    peclet = normb*hK/(2*eps);
    if (peclet > 1)
        //tau = hK/(3*normb);
	tau = hK;
    else
	tau = hK*hK/eps;
    tau = hK;
    return(tau);
	    
}

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  
  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];

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
    MatrixRow = MatrixA[i];
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
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 1, (reduced) SDFEM or (simplified RFB)
// ======================================================================

void NSType1SDFEM3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double tau, ugrad;
  
  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];

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

  tau = RFB_Parameter3D(hK,c0,&param[0]); // stabilization parameter
    
  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    ugrad  = test000+tau * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*ugrad;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad;

      MatrixRow[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];

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
    MatrixRow = MatrixA[i];
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
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
        test001*ansatz001);

      MatrixRow[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 1, Smagorinsky 
// ======================================================================
void NSType1Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];

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
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
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
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 1, VMSProjection, nabla form
// ======================================================================
void NSType1VMSProjection3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22;
  double **Matrix_tilde_G33, **Matrix_G11, **Matrix_G22, **Matrix_G33;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, delta;
  //cout << "VMS" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixL   = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  Matrix_tilde_G11  = LocMatrices[5];
  Matrix_tilde_G22  = LocMatrices[6];
  Matrix_tilde_G33  = LocMatrices[7];
  Matrix_G11  = LocMatrices[8];
  Matrix_G22  = LocMatrices[9];
  Matrix_G33  = LocMatrices[10];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // p
  Orig5 = OrigValues[5];         // l

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);

  for(i=0;i<N_U;i++)
  {
      Matrix11Row = MatrixA11[i];
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

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
 
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
    }                            // endfor j
  }                              // endfor i

  // matrices for VMS
  for(i=0;i<N_U;i++)
  {
      // projections in the momentum balance
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row = Matrix_tilde_G22[i];
    Matrix33Row = Matrix_tilde_G33[i];
    test100 = Orig0[i];  // velo_x
    test010 = Orig1[i];  // velo_y
    test001 = Orig2[i];  // velo_z

     for(j=0;j<N_L;j++)
    {
	ansatz000 = Orig5[j];  // large scales 
      Matrix11Row[j] -= Mult * mu * ansatz000 * test100;
      Matrix22Row[j] -= Mult * mu * ansatz000 * test010;
      Matrix33Row[j] -= Mult * mu * ansatz000 * test001;
    }
  }

  for(i=0;i<N_L;i++)
  {
      // elliptic projection
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    Matrix33Row = Matrix_G33[i];
    test000 = Orig5[i]; // large scales

    for(j=0;j<N_U;j++)
    {
	ansatz100 = Orig0[j];  // velo_x
	ansatz010 = Orig1[j];  // velo_y
	ansatz001 = Orig2[j];  // velo_z

      Matrix11Row[j] -= Mult * ansatz100 * test000;
      Matrix22Row[j] -= Mult * ansatz010 * test000;
      Matrix33Row[j] -= Mult * ansatz001 * test000;
    }
  }

  for(i=0;i<N_L;i++)
  {
      // mass matrix
    test000 = Orig5[i];
    MatrixRow1 = MatrixL[i];
    for(j=0;j<N_L;j++)
    {
      ansatz000 = Orig5[j];
      MatrixRow1[j] += Mult * ansatz000 * test000;
    }
  }
}

// ======================================================================
// Type 2, Standard Galerkin
// ======================================================================
void NSType2Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];
  MatrixB3T = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
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
    MatrixRow = MatrixA[i];
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
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      
      MatrixRow[j] += Mult * val;
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
}

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2SDFEM3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  
  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];
  MatrixB3T = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;
  // delta nach R. Codina
//  delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
  // delta = 1.0/delta;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*(test000+ugrad);
    
    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs2[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term

      MatrixRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig7[i];
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
}

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void  NSType2Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];
  MatrixB3T = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
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
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      
      MatrixRow[j] += Mult * val;
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
}

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void NSType2Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB3 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];
  MatrixB3T = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
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
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
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
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
 
      MatrixRow[j] += Mult * val;
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
}

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3Galerkin3D(double Mult, double *coeff, 
                       double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
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
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

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
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
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
}


// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
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
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];

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
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      /* val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
                 Matrix33Row[j] += Mult * val;*/
    } // endfor j
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
}

// ======================================================================
// Type 3, Upwind (no convection term), (D(u),D(v))
// ======================================================================
void NSType3UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

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
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      Matrix33Row[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double mu, delta;

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
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);

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
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double mu, viscosity, delta;

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
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);
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
      
      val  = 2*viscosity*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
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
}

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType4Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

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
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

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
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
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

}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

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
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

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
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
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
}

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEM3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;

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
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;

  // delta nach R. Codina
//  delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
  // delta = 1.0/delta;

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
    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*(test000+ugrad);
    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term

      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];
    
    test000 = Orig7[i];
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
}

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SDFEMDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;
  double delta1 = TDatabase::ParamDB->DELTA1;

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
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;
  // delta nach R. Codina
//  delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
//  delta = 1.0/delta;

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
    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*(test000+ugrad);

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term

      Matrix11Row[j] += Mult * val;
      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term
      Matrix33Row[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];
    
    test000 = Orig7[i];
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
}

// ======================================================================
// Type 4, Upwind (no convective term), (grad u, grad v)
// ======================================================================
void NSType4Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixB1  = LocMatrices[9];
  MatrixB2  = LocMatrices[10];
  MatrixB3  = LocMatrices[11];
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

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
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);

      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;
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

}

// ======================================================================
// Type 4, Upwind (no convective term), D(u):D(v)
// ======================================================================
void NSType4UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

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
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

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
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      Matrix33Row[j] += Mult * val;
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
}


// ======================================================================
// Type 4, Standard Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, delta;

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
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

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
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);

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
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = (c0+mu) * (test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;
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

}

// ======================================================================
// Type 4, Standard Smagorinsky, D(u):D(v)
// ======================================================================
void NSType4SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity, delta;

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
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

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
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);
  viscosity = c0 + mu/2.0;

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
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = 2*viscosity*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
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
}

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3;
  
  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
    
  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}
// ======================================================================
// Type 1,(reduced) SDFEM or (simplified) RFB
// ======================================================================
void NSType1NLSDFEM3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double *Rhs1, *Rhs2, *Rhs3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double tau, ugrad;
  
  MatrixA = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  tau = RFB_Parameter3D(hK,c0,&param[0]); // stabilization parameter
    
  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    ugrad  = test000+tau * (u1*test100+u2*test010+u3*test001);

    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;
    Rhs3[i] += Mult*ugrad*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad;

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1_2NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  
  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu
    
  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, Smagorinsky , only nonlinear part
// ======================================================================
void NSType1_2NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3, mu, delta;
  
  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
    
  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);
  c0 += mu;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Mult*Orig0[i];
    test010 = Mult*Orig1[i];
    test001 = Mult*Orig2[i];
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, VMSProjection, nabla form
// Type 2, VMSProjection, nabla form
// ======================================================================
void NSType1_2NLVMSProjection3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11;
  double val, val1;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_L;
  double c0, viscosity, delta;
  double u1, u2, u3, mu;
  //cout << "VMS" << endl;
  MatrixA11 = LocMatrices[0];
  Matrix_tilde_G11  = LocMatrices[1];
  Matrix_tilde_G22  = LocMatrices[2];
  Matrix_tilde_G33  = LocMatrices[3];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u_z
  Orig3 = OrigValues[3];         // u
  Orig4 = OrigValues[4];         // l

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);
  viscosity = mu+c0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val  = viscosity*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

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
      Matrix11Row[j] -= Mult * mu * ansatz000 * test100;
      Matrix22Row[j] -= Mult * mu * ansatz000 * test010;
      Matrix33Row[j] -= Mult * mu * ansatz000 * test001;
    }
  }
}

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2NLSDFEM3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad;
  
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;
  
  MatrixA = LocMatrices[0];
  MatrixB1T = LocMatrices[1];
  MatrixB2T = LocMatrices[2];
  MatrixB3T = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;
  // delta nach R. Codina
  //delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
  //delta = 1.0/delta;
  //OutPut(delta << " ");

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    Rhs1[i] += Mult*(test000+ugrad)*c1;
    Rhs2[i] += Mult*(test000+ugrad)*c2;
    Rhs3[i] += Mult*(test000+ugrad)*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*(test000+ugrad);
   
      MatrixRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];

    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, u3;
  
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
      
      val  = c0*(test100*ansatz100+test010*ansatz010
        + test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22,  **MatrixA33;
  double val, val1, val2, val3, val4;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, u3;
 
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

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

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
      Matrix11Row[j] += Mult * val;

      val  = c0*(val2+2*val3+val4)+val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(val2+val3+2*val4)+val1;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind, (grad u, grad v)
// ======================================================================
void NSType3_4NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
 
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
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010
        + test001*ansatz001);
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind , D(u):D(v)
// ======================================================================
void NSType3_4NLUpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22,  **MatrixA33;
  double val, val1, val2, val3;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
 
  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

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
      
      /* val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010+
        0.5*test001*ansatz001);
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010+
        0.5*test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010+
        test001*ansatz001);
        Matrix33Row[j] += Mult * val;*/
      val1 = test100*ansatz100;
      val2 = test010*ansatz010;
      val3 = test001*ansatz001;

      val  = c0*(2*val1+val2+val3);
      Matrix11Row[j] += Mult * val;

      val  = c0*(val1+2*val2+val3);
      Matrix22Row[j] += Mult * val;

      val  = c0*(val1+val2+2*val3);
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3_4NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, u3, mu, viscosity, delta;

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

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);
  viscosity = c0 + mu;

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
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      
      val  = viscosity*(test100*ansatz100+test010*ansatz010+
                        test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;
      /* do not see why these values are added, 06/03/20      
      val  = mu*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = mu*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = mu*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = mu*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = mu*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = mu*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;
      */
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3_4NLSmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, u3, mu, viscosity, delta;

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

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);
  viscosity = c0 + mu/2.0;

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
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      
      val  = 2*viscosity*(test100*ansatz100+0.5*test010*ansatz010+
        0.5*test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;
      
      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+test010*ansatz010+
        0.5*test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+0.5*test010*ansatz010+
        test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEM3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22,  **MatrixA33;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];
  MatrixB3T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;
  // delta nach R. Codina
  //delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
  //delta = 1.0/delta;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*(test000+ugrad);
    
    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term

      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];
  MatrixB3T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;
  // delta nach R. Codina
  //delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
  //delta = 1.0/delta;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*(test000+ugrad);

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term

      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term
      Matrix22Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term
      Matrix33Row[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i

}

// ======================================================================
// Type 4, SDFEM, (grad u, grad v), with div-div stabilization
// ======================================================================
void NSType4NLSDFEM_DivDiv_3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad, tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  MatrixB3T = LocMatrices[11];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;
  // delta nach R. Codina
  //delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
  //delta = 1.0/delta;
  tau = 1;

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
    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*(test000+ugrad);
    
    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term

      Matrix11Row[j] += Mult * (val+tau*test100*ansatz100);
      Matrix22Row[j] += Mult * (val+tau*test010*ansatz010);
      Matrix33Row[j] += Mult * (val+tau*test001*ansatz001);

      val = tau * test100*ansatz010; // div-div term 
      Matrix12Row[j] += Mult * val;

      val = tau * test100*ansatz001; // div-div term 
      Matrix13Row[j] += Mult * val;

      val = tau * test010*ansatz100; // div-div term 
      Matrix21Row[j] += Mult * val;

      val = tau * test010*ansatz001; // div-div term 
      Matrix23Row[j] += Mult * val;

      val = tau * test001*ansatz100; // div-div term 
      Matrix31Row[j] += Mult * val;

      val = tau * test001*ansatz010; // div-div term 
      Matrix32Row[j] += Mult * val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i
}

// ======================================================================
// Type 4, SDFEM, D(u):D(v) with div-div stabilization
// ======================================================================
void NSType4NLSDFEM_DivDiv_DD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, ugrad,tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double power0 = TDatabase::ParamDB->SDFEM_POWER0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  MatrixB3T = LocMatrices[11];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p_x
  Orig5 = OrigValues[5]; // p_y
  Orig6 = OrigValues[6]; // p_z
  Orig7 = OrigValues[7]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  if(c0 < hK)
    delta = delta0*pow(hK,power0);
  else
    delta = delta1*hK*hK;
  // delta nach R. Codina
  //delta = (4* c0/(hK*hK) + 2*sqrt(u1*u1+u2*u2+u3*u3)/hK);
  //delta = 1.0/delta;
  tau = 1;

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
    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val1 = Mult*(test000+ugrad);

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term
      val += tau*test100*ansatz100; // div-div term 
      Matrix11Row[j] += Mult * val;
      
      val = tau * test100*ansatz010; // div-div term 
      Matrix12Row[j] += Mult * val;

      val = tau * test100*ansatz001; // div-div term 
      Matrix13Row[j] += Mult * val;

      val = tau * test010*ansatz100; // div-div term 
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term
      val += tau*test010*ansatz010; // div-div term 
      Matrix22Row[j] += Mult * val;

      val = tau * test010*ansatz001; // div-div term 
      Matrix23Row[j] += Mult * val;

      val = tau * test001*ansatz100; // div-div term 
      Matrix31Row[j] += Mult * val;

      val = tau * test001*ansatz010; // div-div term 
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugrad; // SD term
      val += tau*test001*ansatz001; // div-div term 
      Matrix33Row[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz100 = Orig4[j];
      ansatz010 = Orig5[j];
      ansatz001 = Orig6[j];
      ansatz000 = Orig7[j];

      val  = -ansatz000 * test100;
      val +=  ansatz100 * ugrad;
      MatrixRow1[j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * ugrad;
      MatrixRow2[j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * ugrad;
      MatrixRow3[j] += Mult*val;
    }
  } // endfor i

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
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *AuxMatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double mu2, delta, mu;
  double u1, u2, u3;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  AuxMatrix = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_x
  Orig3 = OrigValues[3]; // u

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],NULL,NULL,NULL,-4711);

  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*u1;
    Rhs2[i] += Mult*test000*u2;
    Rhs3[i] += Mult*test000*u3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
            
    } // endfor j
  } // endfor i
}
// ======================================================================
// rhs for RFB stabilization
// ======================================================================
void NSRFBRhs3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3;
  double test000;
  double *Orig0;
  int i, N_U;
  double c1, c2, c3;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c3 = coeff[3];                 // f3

  for(i=0;i<N_U;i++)
  {
    test000 = Orig0[i];
    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
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
  double *Rhs1, *Rhs2, *Rhs3, val;
  double test000;
  double *Orig0;
  int i,j, N_U;
  double c1, c2, c3;
  double p_x, p_y, p_z;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  p_x = param[0]; // p_sep to x
  p_y = param[1]; // p_sep to y
  p_z = param[2]; // p_sep to z

  for(i=0;i<N_U;i++)
  {
      test000 = Orig0[i];  
      Rhs1[i] += Mult*test000*(c1-p_x);
      Rhs2[i] += Mult*test000*(c2-p_y);
      Rhs3[i] += Mult*test000*(c3-p_z);
  } // endfor i
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
  double test100,test010,test001,ansatz100,ansatz010,ansatz001;
  double *Orig0, *Orig1, *Orig2;
  double u1, u2,u3, d1u1, d2u1,d3u1, d1u2, d2u2, d3u2, d1u3, d2u3, d3u3;
  int i,j, N_U;
  double c1, c2, c3;
  double p_x, p_y, p_z;

  MatrixA = LocMatrices[0];
  Rhs1 = LocRhs[0];
  
  N_U = N_BaseFuncts[0];
 
  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  if (TDatabase::ParamDB->PRESSURE_SEPARATION==3)
  {
     for(i=0;i<N_U;i++)
     {
        MatrixRow = MatrixA[i];
        test100 = Orig0[i];
        test010 = Orig1[i]; 
	test001 = Orig2[i];
        Rhs1[i] += Mult*(test100*c1+test010*c2+test001*c3);
        for(j=0;j<N_U;j++)
        {
           ansatz100 = Orig0[j];
           ansatz010 = Orig1[j];
	   ansatz001 = Orig2[j];
           MatrixRow[j] += Mult * (ansatz100*test100+ansatz010*test010+ansatz001*test001);
        }
     }
  }
  else
  {
     u1 = param[0]; // u1old
     u2 = param[1]; // u2old
     u3 = param[2]; // u3old
     d1u1 = param[3];
     d1u2 = param[4];
     d1u3 = param[5];
     d2u1 = param[6];
     d2u2 = param[7];
     d2u3 = param[8];
     d3u1 = param[9];
     d3u2 = param[10];
     d3u3 = param[11];
     
     for(i=0;i<N_U;i++)
     {
        MatrixRow = MatrixA[i];
        test100 = Orig0[i];
        test010 = Orig1[i];
	test001 = Orig2[i];
        Rhs1[i] += Mult*(test100*(c1-u1*d1u1-u2*d2u1-u3*d3u1)
                         +test010*(c2-u1*d1u2-u2*d2u2-u3*d3u2)+test001*(c3-u1*d1u3-u2*d2u3-u3*d3u3));
        for(j=0;j<N_U;j++)
        {
           ansatz100 = Orig0[j];
           ansatz010 = Orig1[j];
	   ansatz001 = Orig2[j];
           MatrixRow[j] += Mult * (ansatz100*test100+ansatz010*test010+ansatz001*test001);
        }
     }
  }
}


// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void NSParamsVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old  
  out[2] = in[5]; // u3old  
}

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old
  out[3] = in[6]; // D100(u1old)
  out[4] = in[7]; // D100(u2old)
  out[5] = in[8]; // D100(u3old)
  out[6] = in[9]; // D010(u1old)
  out[7] = in[10]; // D010(u2old)
  out[8] = in[11]; // D010(u3old)
  out[9] = in[12]; // D001(u1old)
  out[10] = in[13]; // D001(u2old)
  out[11] = in[14]; // D001(u3old)
}
void NSParamsVeloExact(double *in, double *out)
{
  double sinx, dsinx, ddsinx, dddsinx, siny, dsiny, ddsiny, dddsiny;
  double polyz, dpolyz, ddpolyz, dddpolyz, u1, u2, u3, x, y, z;

  x = in[0];
  y = in[1];
  z = in[2];

/*  sinx = sin(Pi*x)*sin(Pi*x);
  dsinx = 2*Pi*sin(Pi*x)*cos(Pi*x);
  ddsinx = 2*Pi*Pi*(cos(Pi*x)*cos(Pi*x)-sin(Pi*x)*sin(Pi*x));
  dddsinx = -8*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*x);
  siny = sin(2*Pi*y)*sin(2*Pi*y);
  dsiny = 4*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  ddsiny = 8*Pi*Pi*(cos(2*Pi*y)*cos(2*Pi*y) - sin(2*Pi*y)*sin(2*Pi*y));
  dddsiny = -64*Pi*Pi*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  polyz = z*z*(1-z)*(1-z);
  dpolyz = 2*z*(1-z)*(1-2*z);
  ddpolyz = 2-12*z+12*z*z;
  dddpolyz = 24*z-12;

  u1 = sinx*dsiny*polyz+sinx*siny*dpolyz;
  u2 = -dsinx*siny*polyz+sinx*siny*dpolyz;
  u3 =  -dsinx*siny*polyz-sinx*dsiny*polyz;
*/

  u1 = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
  u2 = cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
  u3 = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
     -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;

  out[0] = u1-in[3]; // u1old
  out[1] = u2-in[4]; // u2old  
  out[2] = u3-in[5]; // u2old  
}
// ========================================================================
// parameters: separated pressure
// ========================================================================
void NSParamsPressSep(double *in, double *out)
{
  out[0] = in[3]; // P_sep to x
  out[1] = in[4]; // P_sep to y
  out[2] = in[5]; // P_sep to z
}
