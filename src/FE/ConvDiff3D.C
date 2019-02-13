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
// %W% %G%
//
// common declaration for all 3D convection diffusion problems
// ======================================================================

#include <Database.h>
#include <LinAlg.h>
#include <ConvDiff.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>

extern double bound;

/******************************************************************************/
//
// computation of the size of a mesh cell in convection direction
// approximation formula by Tezduyar and Park, CMAME 59, 307 - 325, 1986
//
/******************************************************************************/

void Compute_Crosswind_Plane(double *n, double *v1, double *v2)
{
  double nx, ny, nz, len, nn;

  len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (len==0)
  {
    OutPut("no crosswind plane " << endl);
    exit(4711);
  }

  nx = n[0]/len;
  ny = n[1]/len;
  nz = n[2]/len;

  if ( (fabs(nx)>=0.5) || (fabs(ny)>=0.5))
  {
    nn = sqrt(nx*nx+ny*ny);
    v1[0] = ny/nn;
    v1[1] = -nx/nn;
    v1[2] = 0;
    v2[0] = -v1[1]*nz;
    v2[1] = v1[0]*nz;
    v2[2] = v1[1]*nx-v1[0]*ny;
  }
  else
  {
    nn = sqrt(ny*ny+nz*nz);
    v1[0] = 0;
    v1[1] = -nz/nn;
    v1[2] = ny/nn;
    v2[0] = v1[2]*ny-v1[1]*nz;
    v2[1] = - v1[2]*nx;
    v2[2] = v1[1]*nx;
  }
}


void BilinearAssemble(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0];
  c1 = coeff[1];
  c2 = coeff[2];
  c3 = coeff[3];
  c4 = coeff[4];
  c5 = coeff[5];

  
//   cout << " c0 " << c0 <<" c5 " << c5 << endl;
  
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs[i] += Mult*test000*c5;
    
//      if(TDatabase::ParamDB->P14<Mult) TDatabase::ParamDB->P14=Mult;    

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val *= Mult;

      MatrixRow[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// SDFEM without 2nd derivatives
void BilinearAssemble_SD(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, c6;
  double delta, bgradv;
 
  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // b_1
  c2 = coeff[2];                 // b_2
  c3 = coeff[3];                 // b_3
  c4 = coeff[4];                 // c
  c5 = coeff[5];                 // f
  c6 = coeff[6];                 // ||b||

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c4, c6);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    Rhs[i] += Mult*(test000+delta*bgradv)*c5;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val += delta * (c1*ansatz100+c2*ansatz010+c3*ansatz001
        +c4*ansatz000) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                            // endfor j
  }                              // endfor i
}


/*
void BilinearAssemble_UPW1(double Mult, double *coeff, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j, N_;
  double c0;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];

  c0 = coeff[0];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = c0*(test10*ansatz10+test01*ansatz01);

      val *= Mult;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}
*/

void BilinearAssemble_UPW2(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c4, c5;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0];
  c4 = coeff[4];
  c5 = coeff[5];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs[i] += Mult*test000*c5;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += c4*ansatz000*test000;

      val *= Mult;

      MatrixRow[j] += val;
    }                            // endfor j
  }                              // endfor i
}

// ========================================================================
// SOLD schemes which add isotropic diffusion
// Hughes, Mallet, Mizukami (1986)
// Tezduyar, Park (1986)
// Galeao, do Carmo (1988)
// do Carmo, Galeao (1991)
// do Carmo, Alvarez (2003)
// Almeida, Silva (1997)
// Knopp, Lube, Rapin (2002)
// Johnson (1990)
// Johnson (1992)
// ========================================================================

void BilinearAssemble_SOLD(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, c6;
  double delta, bgradv, sigma, norm_b;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // b_1
  c2 = coeff[2];                 // b_2
  c3 = coeff[3];                 // b_3
  c4 = coeff[4];                 // c
  c5 = coeff[5];                 // f
  c6 = coeff[6];                 // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c4, c6);

  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, c6, delta, param, 0, 0, 0);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = (c0+sigma)*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val += delta * (c1*ansatz100+c2*ansatz010+c3*ansatz001
        +c4*ansatz000) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                            // endfor j
  }                              // endfor i
}


// ========================================================================
// SOLD schemes which add diffusion only orthogonal to convection
// Johnson,Schatz,Wahlbin (1987) (linear)
// Knopp, Lube, Rapin (2002)
// Codina (1993)
// ========================================================================

void BilinearAssemble_SOLD_Orthogonal(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, c6;
  double delta, bgradv, sigma, norm_b;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // b_1
  c2 = coeff[2];                 // b_2
  c3 = coeff[3];                 // b_3
  c4 = coeff[4];                 // c
  c5 = coeff[5];                 // f
  c6 = coeff[6];                 // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c4, c6);

  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, c6, delta, param, 0, 0, 0);

  norm_b = c1*c1+c2*c2+c3*c3;

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val += delta * (c1*ansatz100+c2*ansatz010+c3*ansatz001
        +c4*ansatz000) * bgradv;

      if (norm_b > 0)
        val += sigma*( (c2*c2+c3*c3)*ansatz100*test100 + (c1*c1+c3*c3)*ansatz010*test010
          + (c1*c1+c2*c2)*ansatz001*test001
          -c1*c3*(ansatz001*test100+ansatz100*test001)
          -c1*c2*(ansatz010*test100+ansatz100*test010)
          -c2*c3*(ansatz001*test010+ansatz010*test001))/norm_b;

      val *=Mult;

      MatrixRow[j] += val;

    }                            // endfor j
  }                              // endfor i
}


// ========================================================================
// parameters:  H1 norm and  norm of residual
// ========================================================================
void DC_CD_Params(double *in, double *out)
{
  out[0] = in[3];                // H1 norm
  out[1] = in[4];                // norm of residual
}


// ========================================================================
// parameters:  partial derivatives
// ========================================================================
void MBE_Params(double *in, double *out)
{
  out[0] = in[3];                // u
  out[1] = in[4];                // u_x
  out[2] = in[5];                // u_y
  out[3] = in[6];                // u_z
}


// ========================================================================
// parameters:  SC_2
// ========================================================================
void SC_2_Params(double *in, double *out)
{
  out[0] = in[3];                // H1 norm
  out[1] = in[4];                // norm of residual
  out[2] = in[5];                // u_x
  out[3] = in[6];                // u_y
  out[4] = in[7];                // u_z
}


// ========================================================================
// parameters:  SOLD
// ========================================================================
void SOLD_Params(double *in, double *out)
{
  out[0] = in[3];                // u
  out[1] = in[4];                // u_x
  out[2] = in[5];                // u_y
  out[3] = in[6];                // u_z
  out[4] = in[7];                // ||u^h||_{H^1,K}
  out[5] = in[8];                // ||R(u^h)||_{L^2,K}
}
