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
// TNSE2D_FixPo_SSMUM.C        
//
// declarations for SSMUM
//
// author: Volker John 08/05/22
//
// ======================================================================

#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>

#include <stdlib.h>

// ======================================================================
//
// WITH ROTATING FRAME
//
// ======================================================================

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, Coletti, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
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
  double u1, u2, omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;  

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

  x = param[0];                  // x
  y = param[1];                  // y
  u1 = param[2];                 // u1old
  u2 = param[3];                 // u2old

  // compute angular velocity for current point
  r[0] = x - mp_x;
  r[1] = y - mp_y;
  //OutPut("a:" << r[0] << " " << r[1] << endl);
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      omega = 0.0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
  }

  //omega = 0;
  //OutPut(x << " " << y << " "  << omega  << endl);
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

    //if (TDatabase::TimeDB->CURRENTTIME == 0.0)
    {
	Rhs1[i] += Mult*(c1+omega*omega*r[0])*test00;
	Rhs2[i] += Mult*(c2+omega*omega*r[1])*test00;
	//Rhs1[i] += Mult*(c1+0)*test00;
	//Rhs2[i] += Mult*(c2+1)*test00;
    }
    /*else
    {
	Rhs1[i] += Mult*test00*c1+ omega*omega*r[0]*test00;
	Rhs2[i] += Mult*test00*c2+ omega*omega*r[1]*test00;
	}*/

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = 2*omega * ansatz00*test00;
      Matrix12Row[j] += Mult * val;

      val  = -2*omega * ansatz00*test00;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
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
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, Coletti, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow, **AuxMatrix;;
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

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

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

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
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
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, Coletti, (grad u, grad v), only nonlinear part
// Type 4, Coletti, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;  

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  x = param[0];                  // x
  y = param[1];                  // y
  u1 = param[2];                 // u1old
  u2 = param[3];                 // u2old

  // compute angular velocity for current point
  r[0] = x - mp_x;
  r[1] = y - mp_y;
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      omega = 0.0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
  }

  //omega = 0;
  //OutPut(x << " " << y << " "  << omega  << endl);
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
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = 2* omega * ansatz00*test00;
      Matrix12Row[j] += Mult * val;

      val  = -2* omega * ansatz00*test00;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
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

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// right-hand side ONLY, for SSMUM
// ======================================================================
void TimeNSRHS_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2, val;
  double omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;  
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  x = param[0];                  // x
  y = param[1];                  // y
  
  // compute angular velocity for current point
  r[0] = x - mp_x;
  r[1] = y - mp_y;
  //OutPut("r:" << r[0] << " " << r[1] << endl);
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      omega = 0.0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
  }

  //OutPut(c1 << " " << c2 << " "  << omega  << endl);
  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*(c1+omega*omega*r[0])*test00;
    Rhs2[i] += Mult*(c2+omega*omega*r[1])*test00;
    //Rhs1[i] += Mult*test00*c1;
    //Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " ";
  }                              // endfor i
}
// ======================================================================
// right-hand side ONLY, for SSMUM
// ======================================================================
void TimeNS_REL_VELO_SSMUM_WITH_ROTFRAME(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2, val;
  double omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;  
  

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  x = param[0];                  // x
  y = param[1];                  // y
  
  // compute angular velocity for current point
  r[0] = x - mp_x;
  r[1] = y - mp_y;
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      omega = 0.0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
  }

  //OutPut(x << " " << y << " "  << omega  << endl);
  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];
    //OutPut(Mult*test00*omega << " " << Mult << " " << test00 << " " << omega << endl);

    //Rhs1[i] += Mult*test00*omega*r[1];
    //Rhs2[i] += -Mult*test00*omega*r[0];
    Rhs1[i] += Mult*test00*omega;//*r[1];
    Rhs2[i] += -Mult*test00*omega;//*r[0];
    //Rhs1[i] += Mult*test00*c1;
    //Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " : ";
  }                              // endfor i
}

// ======================================================================
//
// ALE
// all matrices have to be assembled newly if the mesh changes
//
// ======================================================================

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, Coletti, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin_SSMUM_ALE(double Mult, double *coeff,
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
  double u1, u2, omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS; 

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

  x = param[0];                  // x
  y = param[1];                  // y
  u1 = param[2];                 // u1old
  u2 = param[3];                 // u2old

  // compute angular velocity for current point
  r[0] = y - mp_y;
  r[1] = mp_x - x;
  //OutPut("a:" << r[0] << " " << r[1] << endl);
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      r[0] = r[1] = 0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
      r[0] *= omega;
      r[1] *= omega;
  }

  u1 -= r[0];
  u2 -= r[1];
  //u1 = u2 = 0;
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
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, Coletti, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow, **AuxMatrix;;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;  ;

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

  x = param[0];                  // x
  y = param[1];                  // y
  u1 = param[2];                 // u1old
  u2 = param[3];                 // u2old

  // compute angular velocity for current point
  r[0] = y - mp_y;
  r[1] = mp_x - x;
  //OutPut("a:" << r[0] << " " << r[1] << endl);
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      omega = 0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
  }

  r[0] *= omega;
  r[1] *= omega;

  u1 -= r[0];
  u2 -= r[1];
  
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

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
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
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, Coletti, (grad u, grad v), only nonlinear part
// Type 4, Coletti, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin_SSMUM_ALE(double Mult, double *coeff,
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
  int i,j, N_U;
  double c0;
  double u1, u2, omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;  

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  x = param[0];                  // x
  y = param[1];                  // y
  u1 = param[2];                 // u1old
  u2 = param[3];                 // u2old
  
  //OutPut(x << " " << y <<  " " << fabs(u1 - 2*Pi*y) << " : " << fabs(u2 + 2*Pi*x) << endl);

  // compute angular velocity for current point
  r[0] = y - mp_y;
  r[1] = mp_x - x;
  //OutPut("a:" << r[0] << " " << r[1] << endl);
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      r[0] = r[1] = 0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
      r[0] *= omega;
      r[1] *= omega;
  }
  u1 -= r[0];
  u2 -= r[1];
  // u1 = u2 = 0;

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
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2, omega, r[2], rad, x, y;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;  

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  x = param[0];                  // x
  y = param[1];                  // y
  u1 = param[2];                 // u1old
  u2 = param[3];                 // u2old

  // compute angular velocity for current point
  r[0] = y - mp_y;
  r[1] = mp_x - x;
  //OutPut("a:" << r[0] << " " << r[1] << endl);
  rad = sqrt(r[0]*r[0] + r[1]*r[1]);
  if (rad>=  outer_radius )
  {
      r[0] = r[1] = 0;
      
  }
  else
  {
      if (rad <= inner_radius)
      {
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      }
      else
      {
	  val = 1 - (rad-inner_radius)/(outer_radius-inner_radius);
	  omega = 2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND * val;
      }
      r[0] *= omega;
      r[1] *= omega;
  }

  u1 -= r[0];
  u2 -= r[1];

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

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// right-hand side ONLY, for NSE
// ======================================================================
void TimeNSRHS_SSMUM_ALE(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " ";
  }                              // endfor i
}
