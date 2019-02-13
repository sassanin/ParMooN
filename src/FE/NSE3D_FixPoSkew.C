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
// NSE3D_FixPoSkew.C 
//
// common declaration for all Navier-Stokes problems
// fix point iteration
// skew symmetric from of the convective term
//
// ATTENTION: ONLY FEW OF THE ROUTINES ARE TESTED
//
// ======================================================================

#include <Database.h>
#include <MooNMD_Io.h>
#include <TNSE3D_Routines.h>
#include <Convolution.h>
#include <NSE3D_FixPo.h>


// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1GalerkinSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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

    val2 = 0.5*(u1*test100+u2*test010+u3*test001);
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= val2*ansatz000;

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

void NSType1SDFEMSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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
    ugrad  = tau * (u1*test100+u2*test010+u3*test001);
    val2 =  0.5*(u1*test100+u2*test010+u3*test001);
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
     
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001);
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*val1*test000;
      val -= val2*ansatz000;
      val += val1*ugrad;

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
// Type 1, upwind is not available
// ======================================================================

// ======================================================================
// Type 1, Smagorinsky 
// ======================================================================
void NSType1SmagorinskySkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;

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
// Type 2, Standard Galerkin
// ======================================================================
void NSType2GalerkinSkew3D(double Mult, double *coeff, 
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
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
      
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
void NSType2SDFEMSkew3D(double Mult, double *coeff, 
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
      ansatz000 = Orig3[j];
     
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001);
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*val1*test000;
      val -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
      val += val1*ugrad; // SD term

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
// Type 2, Upwind is not available
// ======================================================================

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void NSType2SmagorinskySkew3D(double Mult, double *coeff, 
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
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
 
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
void NSType3GalerkinSkew3D(double Mult, double *coeff, 
                       double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
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
void NSType3GalerkinDDSkew3D(double Mult, double *coeff, 
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
  double ansatz100, ansatz010, ansatz001, ansatz000;
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
      ansatz000 = Orig3[j];
      
      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;

      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
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
// Type 3, Upwind is not available
// ======================================================================

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3SmagorinskySkew3D(double Mult, double *coeff, 
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
  double ansatz100, ansatz010, ansatz001, ansatz000;
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
      ansatz000 = Orig3[j];

      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
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
void NSType3SmagorinskyDDSkew3D(double Mult, double *coeff, 
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
  double ansatz100, ansatz010, ansatz001, ansatz000;
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
      ansatz000 = Orig3[j];
      
      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
      
      val  = 2*viscosity*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
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
void NSType4GalerkinSkew3D(double Mult, double *coeff, 
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
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
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
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================
void NSType4GalerkinDDSkew3D(double Mult, double *coeff, 
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
      ansatz000 = Orig3[j];

      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 -= 0.5*(u1*test100+u2*test010+u3*test001)*ansatz000;
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
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
void NSType4SDFEMSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2, val3;
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
    val3 = 0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
     
      val2 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      val1 =  0.5*val2*test000;
      val1 -= 0.5*val3*ansatz000;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += val1;
      val += val2*ugrad; // SD term

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
void NSType4SDFEMDDSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2, val3;
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

    val3 = 0.5*(u1*test100+u2*test010+u3*test001);
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val2 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      val1 =  0.5*val2*test000;
      val1 -= 0.5*val3*ansatz000;

      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      val += val2*ugrad; // SD term

      Matrix11Row[j] += Mult * val;
      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      val += val2*ugrad; // SD term
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      val += val2*ugrad; // SD term
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
// Type 4, Upwind not available
// ======================================================================

// ======================================================================
// Type 4, Standard Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4SmagorinskySkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
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

    val2 = 0.5*(u1*test100+u2*test010+u3*test001);
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu) * (test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= val2*ansatz000;

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
void NSType4SmagorinskyDDSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
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

    val2 = 0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 -= val2*ansatz000;

      val  = 2*viscosity*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*viscosity*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
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
void NSType1_2NLGalerkinSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3, val1, val2;
  
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
    val2 = 0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= val2*ansatz000;

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}
// ======================================================================
// Type 1,(reduced) SDFEM or (simplified) RFB
// ======================================================================
void NSType1NLSDFEMSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double *Rhs1, *Rhs2, *Rhs3;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0, c1, c2, c3;
  double u1, u2, u3, val1, val2;
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

    val2 = 0.5*(u1*test100+u2*test010+u3*test001);
    ugrad  = tau * (u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
     
      val1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*val1*test000;
      val -= val2*ansatz000;
      val += val1*ugrad;

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 1, for upwind not available
// ======================================================================

// ======================================================================
// Type 1, Smagorinsky , only nonlinear part
// ======================================================================
void NSType1_2NLSmagorinskySkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3, mu, delta,val2;
  
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
    val2 = 0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= val2*ansatz000;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2NLSDFEMSkew3D(double Mult, double *coeff, 
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
  double c0, c1, c2, c3, val1, val2;
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

    ugrad  = test000+delta * (u1*test100+u2*test010+u3*test001);
    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;
    Rhs3[i] += Mult*ugrad*c3;
    ugrad  = delta * (u1*test100+u2*test010+u3*test001);
    val2 =  0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001); 
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*val1*test000;
      val -= val2*ansatz000;
      val += val1*ugrad;
   
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
void NSType3_4NLGalerkinSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val, val2;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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

    val2 = 0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010
        + test001*ansatz001);

      val += 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val -= val2*ansatz000;

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
void NSType3_4NLGalerkinDDSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22,  **MatrixA33;
  double val, val1, val2, val3, val4, val5, val6;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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

    val5 =  0.5*(u1*test100+u2*test010+u3*test001);
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      
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
      val6 = val1-val5*ansatz000;

      val  = c0*(2*val2+val3+val4)+val6;
      Matrix11Row[j] += Mult * val;

      val  = c0*(val2+2*val3+val4)+val6;
      Matrix22Row[j] += Mult * val;

      val  = c0*(val2+val3+2*val4)+val6;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind  not available
// ======================================================================

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3_4NLSmagorinskySkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val, val1, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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

    val2 = 0.5*(u1*test100+u2*test010+u3*test001);
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      
      val  = viscosity*(test100*ansatz100+test010*ansatz010+
                        test001*ansatz001);
      val += val1-val2*ansatz000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;
      
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
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3_4NLSmagorinskyDDSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val, val1, val2;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
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

    val2 = 0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val1 = 0.5*(u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 -= val2*ansatz000;

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
void NSType4NLSDFEMSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22,  **MatrixA33;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
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

    val2 =  0.5*(u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
     
      val1 = u1*ansatz100+u2*ansatz010+u3*ansatz001;
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += 0.5*val1*test000;
      val -= val2*ansatz000;
      val += val1*ugrad;

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
void NSType4NLSDFEMDDSkew3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, val2;
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

    val2 =  0.5*(u1*test100+u2*test010+u3*test001);
    
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001);     
      val1 = 0.5*val1-val2*ansatz000+val1*ugrad;        

      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
    
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
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
