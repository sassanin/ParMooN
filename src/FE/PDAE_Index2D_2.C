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
// @(#)PorMedia2D.C        1.0 07/25/00
//
// common declaration for time dependent convection diffusion problems
// with two equations
// ======================================================================

#include <Database.h>

void PDAE2D_2_Params2(double *in, double *out)
{
  out[0] = in[2]; // u1old
  out[1] = in[3]; // u2old
  out[2] = in[4]; // u2old
  out[3] = in[5]; // u1old
  out[4] = in[6]; // u2old
  out[5] = in[7]; // u2old
  
}


void TimeMassAssemble_PDAE2(double Mult, double *coeff, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{ double **Matrix11, **Matrix12;
  double **Matrix21, **Matrix22;
  double *MatrixRow11, *MatrixRow12;
  double *MatrixRow21, *MatrixRow22;
  double ansatz00, ansatz10, ansatz01, val;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3;

  Matrix11 = LocMatrices[0];
  Matrix21 = LocMatrices[1];
  Matrix12 = LocMatrices[2];
  Matrix22 = LocMatrices[3];
  
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  c0 = coeff[0];
  c1 = coeff[1];
  c2 = coeff[2];
  c3 = coeff[3];
  
  for(i=0;i<N_;i++)
  {
    MatrixRow11 = Matrix11[i];
    MatrixRow12 = Matrix12[i];
    MatrixRow21 = Matrix21[i];
    MatrixRow22 = Matrix22[i];
    
    test00 = Orig0[i];

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];
      
      MatrixRow11[j] += c0*Mult*ansatz00*test00;
      MatrixRow12[j] += c1*Mult*ansatz00*test00;
      MatrixRow21[j] += c2*Mult*ansatz00*test00;
      MatrixRow22[j] += c3*Mult*ansatz00*test00;
    } // endfor j
  } // endfor i
}

void TimeBilinearAssemble_PDAE2(double Mult, double *coeff, double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix11, **Matrix12;
  double **Matrix21, **Matrix22;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow11, *MatrixRow12;
  double *MatrixRow21, *MatrixRow22;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  
  Matrix11 = LocMatrices[0];
  Matrix21 = LocMatrices[1];
  Matrix12 = LocMatrices[2];
  Matrix22 = LocMatrices[3];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  for(i=0;i<N_;i++)
  {
    MatrixRow11 = Matrix11[i];
    MatrixRow12 = Matrix12[i];
    MatrixRow21 = Matrix21[i];
    MatrixRow22 = Matrix22[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*coeff[24];
    Rhs2[i] += Mult*test00*coeff[25];
    
    /*
    cout << "coeff " << coeff[24] << " : " << coeff[25] <<endl;
    */
    
    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = (coeff[0]*ansatz10 + coeff[1]*ansatz01)*test10;
      val += (coeff[1]*ansatz10 + coeff[2]*ansatz01)*test01;
      val += (coeff[3]*ansatz10 + coeff[4]*ansatz01 + coeff[5]*ansatz00)*test00;
      val *=Mult;

      MatrixRow11[j] += val;
      
      val  = (coeff[6]*ansatz10 + coeff[7]*ansatz01)*test10;
      val += (coeff[7]*ansatz10 + coeff[8]*ansatz01)*test01;
      val += (coeff[9]*ansatz10 + coeff[10]*ansatz01 + coeff[11]*ansatz00)*test00;
      val *=Mult;

      MatrixRow12[j] += val;
      
      val  = (coeff[12]*ansatz10 + coeff[13]*ansatz01)*test10;
      val += (coeff[13]*ansatz10 + coeff[14]*ansatz01)*test01;
      val += (coeff[15]*ansatz10 + coeff[16]*ansatz01 + coeff[17]*ansatz00)*test00;
      val *=Mult;

      MatrixRow21[j] += val;
      
      val  = (coeff[18]*ansatz10 + coeff[19]*ansatz01)*test10;
      val += (coeff[19]*ansatz10 + coeff[20]*ansatz01)*test01;
      val += (coeff[21]*ansatz10 + coeff[22]*ansatz01 + coeff[23]*ansatz00)*test00;
      val *=Mult;

      MatrixRow22[j] += val;
    } // endfor j
  } // endfor i
}

void TimeBilinearAssembleC_PDAE2(double Mult, double *coeff, double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix11, **Matrix12;
  double **Matrix21, **Matrix22;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow11, *MatrixRow12;
  double *MatrixRow21, *MatrixRow22;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  
  Matrix11 = LocMatrices[0];
  Matrix21 = LocMatrices[1];
  Matrix12 = LocMatrices[2];
  Matrix22 = LocMatrices[3];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  for(i=0;i<N_;i++)
  {
    MatrixRow11 = Matrix11[i];
    MatrixRow12 = Matrix12[i];
    MatrixRow21 = Matrix21[i];
    MatrixRow22 = Matrix22[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*coeff[24];
    Rhs2[i] += Mult*test00*coeff[25];
    
    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = (coeff[0]*ansatz10 + coeff[1]*ansatz01)*test10;
      val += (coeff[1]*ansatz10 + coeff[2]*ansatz01)*test01;
      val += (coeff[3]*ansatz10 + coeff[4]*ansatz01 + coeff[5]*ansatz00)*test00;
      val *=Mult;

      MatrixRow11[j] += val;
      
      val  = (coeff[6]*ansatz10 + coeff[7]*ansatz01)*test10;
      val += (coeff[7]*ansatz10 + coeff[8]*ansatz01)*test01;
      val += (coeff[9]*ansatz10 + coeff[10]*ansatz01 + coeff[11]*ansatz00)*test00;
      val *=Mult;

      MatrixRow12[j] += val;
      
      val  = (coeff[12]*ansatz10 + coeff[13]*ansatz01)*test10;
      val += (coeff[13]*ansatz10 + coeff[14]*ansatz01)*test01;
      val += (coeff[15]*ansatz10 + coeff[16]*ansatz01 + coeff[17]*ansatz00)*test00;
      val *=Mult;

      MatrixRow21[j] += val;
      
      val  = (coeff[18]*ansatz10 + coeff[19]*ansatz01)*test10;
      val += (coeff[19]*ansatz10 + coeff[20]*ansatz01)*test01;
      val += (coeff[21]*ansatz10 + coeff[22]*ansatz01 + coeff[23]*ansatz00)*test00;
      val *=Mult;

      MatrixRow22[j] += val;
    } // endfor j
  } // endfor i
}

void TimeBilinearAssembleJ_PDAE2(double Mult, double *coeff, double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix11, **Matrix12;
  double **Matrix21, **Matrix22;
  double val;
  double *MatrixRow11, *MatrixRow12;
  double *MatrixRow21, *MatrixRow22;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  
  Matrix11 = LocMatrices[0];
  Matrix21 = LocMatrices[1];
  Matrix12 = LocMatrices[2];
  Matrix22 = LocMatrices[3];
  
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  for(i=0;i<N_;i++)
  {
    MatrixRow11 = Matrix11[i];
    MatrixRow12 = Matrix12[i];
    MatrixRow21 = Matrix21[i];
    MatrixRow22 = Matrix22[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = (coeff[0]*ansatz10 + coeff[1]*ansatz01)*test10;
      val += (coeff[1]*ansatz10 + coeff[2]*ansatz01)*test01;
      val += (coeff[3]*ansatz10 + coeff[4]*ansatz01 + coeff[5]*ansatz00)*test00;
      val += coeff[24]*ansatz00*test00;
      val *=Mult;

      MatrixRow11[j] += val;
      
      val  = (coeff[6]*ansatz10 + coeff[7]*ansatz01)*test10;
      val += (coeff[7]*ansatz10 + coeff[8]*ansatz01)*test01;
      val += (coeff[9]*ansatz10 + coeff[10]*ansatz01 + coeff[11]*ansatz00)*test00;
      val += coeff[25]*ansatz00*test00;
      val *=Mult;

      MatrixRow12[j] += val;
      
      val  = (coeff[12]*ansatz10 + coeff[13]*ansatz01)*test10;
      val += (coeff[13]*ansatz10 + coeff[14]*ansatz01)*test01;
      val += (coeff[15]*ansatz10 + coeff[16]*ansatz01 + coeff[17]*ansatz00)*test00;
      val += coeff[26]*ansatz00*test00;
      val *=Mult;

      MatrixRow21[j] += val;
      
      val  = (coeff[18]*ansatz10 + coeff[19]*ansatz01)*test10;
      val += (coeff[19]*ansatz10 + coeff[20]*ansatz01)*test01;
      val += (coeff[21]*ansatz10 + coeff[22]*ansatz01 + coeff[23]*ansatz00)*test00;
      val += coeff[27]*ansatz00*test00;
      val *=Mult;

      MatrixRow22[j] += val;
    } // endfor j
  } // endfor i
}

void TimeRhsAssemble_PDAE2(double Mult, double *coeff, double *param,
                    double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c1 = coeff[24]; // Q
  c2 = coeff[25];
  
  for(i=0;i<N_;i++)
  {
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    
  } // endfor i
}
