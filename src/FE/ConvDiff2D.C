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
// @(#)ConvDiff2D.C        1.4 06/27/00
//
// common declaration for all convection diffusion problems
// ======================================================================

#ifdef __2D__

#include <Database.h>
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>

#include <stdlib.h>

extern double bound;



/** @brief the local assembling routines */

/** ========================================================================= */
/** ========================================================================= */
// CD2D: stationary convection diffusion problems

// Galerkin
void BilinearAssembleGalerkin(double Mult, double *coeff, double* param, 
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs)
{
  double val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  int j;
  
  double **Matrix = LocMatrices[0];
  double *Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  double *Orig0 = OrigValues[0];
  double *Orig1 = OrigValues[1];
  double *Orig2 = OrigValues[2];

  // coefficients
  const double c0 = coeff[0]; // eps
  const double c1 = coeff[1]; // b_1
  const double c2 = coeff[2]; // b_2
  const double c3 = coeff[3]; // c
  const double c4 = coeff[4]; // f

  for(int i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig0[i]; // xi derivative
    test01 = Orig1[i]; // eta derivative
    test00 = Orig2[i]; // function
    
    // assemble rhs
    // quad_weigth * test_function * f
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j]; // xi derivative
      ansatz01 = Orig1[j]; // eta derivative
      ansatz00 = Orig2[j]; // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assemble reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;

      // quad weigth
      val *= Mult;

      // update matrix entry
      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

// SDFEM - Streamline Diffusion Finite Element Method, SUPG - Streamline Upwind 
// Petrov Galerkin
void BilinearAssemble_SD(double Mult, double *coeff, double* param, double hK,
                         double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs)
{
  double val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  int j;
  double bgradv;

  double **Matrix = LocMatrices[0];
  double *Rhs = LocRhs[0];

  const int N_ = N_BaseFuncts[0];

  const double *Orig0 = OrigValues[0];
  const double *Orig1 = OrigValues[1];
  const double *Orig2 = OrigValues[2];
  const double *Orig3 = OrigValues[3];
  const double *Orig4 = OrigValues[4];

  const double c0 = coeff[0]; // nu
  const double c1 = coeff[1]; // b_1
  const double c2 = coeff[2]; // b_2
  const double c3 = coeff[3]; // c
  const double c4 = coeff[4]; // f
  const double c5 = coeff[5]; // \|b\|_infty (in the entire cell)

  const double delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);

  for(int i=0;i<N_;i++)
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
      
      // assemble viscous term
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assemble reactive term
      val += c3*ansatz00*test00;
      // assemble stabilization
      val += delta * ( -c0*(ansatz20 + ansatz02)
                       +c1*ansatz10 + c2*ansatz01
                       +c3*ansatz00 ) * bgradv;
      
      // quad weigth
      val *=Mult;
      
      // update matrix entry
      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

// GLS - Galerkin Least Squares
void BilinearAssemble_GLS(double Mult, double *coeff, double* param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, Lu;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    test20 = Orig3[i];
    test02 = Orig4[i];
      
    Lu = (-c0*(test20+test02)+c1*test10+c2*test01 +c3*test00);
    Rhs[i] += Mult*(test00+delta*Lu)*c4;

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
             +c3*ansatz00) * Lu;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}

void BilinearAssemble_Axial3D(double Mult, double *coeff, double *param,
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, x, r;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  // coefficients
  c0 = coeff[0];                                  // eps
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  x  = coeff[20];
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check BilinearAssemble_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
//    cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);
   

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig0[i];                            // xi derivative
    test01 = Orig1[i];                            // eta derivative
    test00 = Orig2[i];                            // function

    // assemble rhs
    // quad_weigth * test_function * f
    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];                        // xi derivative
      ansatz01 = Orig1[j];                        // eta derivative
      ansatz00 = Orig2[j];                        // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assembel reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;

      // quad weigth
      val *= Mult;
      
      //axial 3d
      val *= r;
      // update matrix entry
      MatrixRow[j] += val;
    }                                             // endfor j
  }                                               // endfor i
}



void BilinearAssemble_UPW1(double Mult, double *coeff, double *param, double hK,
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
    }                                             // endfor j
  }                                               // endfor i
}


void BilinearAssemble_UPW2(double Mult, double *coeff, double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c3, c4;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0];
  c3 = coeff[3];
  c4 = coeff[4];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
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
      val += c3*ansatz00*test00;

      val *= Mult;

      MatrixRow[j] += val;
    }                                             // endfor j
  }                                               // endfor i
}



/** ========================================================================= */
/** ========================================================================= */
// TCD2D: time dependent convection diffusion problems

// Galerkin
void MatrixMRhsAssemble(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, *MatrixRow;
  double ansatz00;
  double test00;
  double *Orig0;
  int i,j, N_;
  double c4; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  c4 = coeff[4]; // f

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];

    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];

      MatrixRow[j] += Mult*ansatz00*test00;
    } // endfor j
  } // endfor i
}
void MatrixARhsAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, val, *MatrixRowA, *MatrixRowK;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, Pe, h; 

  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f

  if ((TDatabase::ParamDB->DISCTYPE==5)||(TDatabase::ParamDB->DISCTYPE==6)
      ||(TDatabase::ParamDB->DISCTYPE==7))
  {
    h = ComputeAlpha(hK);
    c0+= h;  
  }
  for(i=0;i<N_;i++)
  {
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

      MatrixRowA[j] += Mult * val;
                
    } // endfor j
  } // endfor i
}


// SUPG
void MatricesAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, *MatrixRowA, *MatrixRowK;
  double **MatrixS, *MatrixRowS, ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double val, val1, val2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, Pe; 
  double c00, c11, c22, c33;
  double tau, bgradv, bb, res, sigma, norm_b;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  MatrixA = LocMatrices[0];
  MatrixK= LocMatrices[1];
  
  if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)  
   MatrixS= LocMatrices[2];
  
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  // coefficients of the problem
  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f
  
  // coefficients for the stabilization parameter
  val = theta1 * time_step;
  c00 = val * c0;
  c11 = val * c1;
  c22 = val * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + val * c3;
  if (TDatabase::ParamDB->SDFEM_TYPE==8)
  {
    c33 = val * c3;
  }
  if (fabs(c11) > fabs(c22))
    bb = fabs(c11);
  else
    bb = fabs(c22);
  // this is tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  if(TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
  {
    // rhs from previous time step 
    c5 = coeff[5];   
    // compute residual
    res = param[0] + theta1*time_step*(-c0*(param[3]+param[4]) + c1*param[1]
          +c2*param[2] + c3*param[0])
          -param[5] +  theta2*time_step*(-c0*(param[8]+param[9]) + c1*param[6]
          +c2*param[7] + c3*param[5])
          -theta3*time_step*c5 - theta4*time_step*c4;
    /*
    c00 =  time_step * theta1 * c0;
    c11 =  time_step * theta1 * c1;
    c22 =  time_step * theta1 * c2;
    c33 = 1.0 + time_step * theta1 *c3;
    */
    c5 = time_step * theta4 * c4;
    // compute the parameter, c5 is just a dummy
    sigma = Compute_SOLD_sigma(hK, c00, c11, c22, c33, c5, bb, tau, param, res, 1,1);
    //OutPut( param[0] << " " << param[5] <<  " " << res << " " <<  sigma << endl);
    val2 = Mult * sigma;
    if (TDatabase::ParamDB->SOLD_TYPE==2)
    {
      norm_b = c1*c1 + c2*c2;
      if (norm_b >1e-10)
        val2 /= norm_b;
      else
        val2 = 0.0;
    }
  }
  // scale appropriately, after it is used for the SOLD scheme
  // do not apply for paper with J. Novo
  // this is \tilde tau
  if((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)
      &&(TDatabase::ParamDB->SDFEM_TYPE!=11))
    tau *= val;

  // loop over the basis functions
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    MatrixRowK = MatrixK[i];
    if(TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
      MatrixRowS = MatrixS[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    // scaling with the stabilization parameter
    // scaling with the time step is done in the main program
    bgradv *= tau;
    // THIS CHANGEs TEST00 !
    test00 += bgradv;
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];
  
      // Galerkin part of the bilinear form
      val1 = c1*ansatz10+c2*ansatz01;
      val1+= c3*ansatz00;

      val = c0*(test10*ansatz10+test01*ansatz01);
      // val1*test00 includes the SUPG part of the convective and reactive term
      val += val1*test00;
      // diffusion part of the SUPG stabilization
      val -= c0*(ansatz20 + ansatz02) * bgradv;
      
      MatrixRowA[j] += Mult * val;
      // time derivative part of the SUPG stabilization
      MatrixRowK[j] += Mult * ansatz00*bgradv;
      
      // isotropic SOLD method
      if((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) 
          && (TDatabase::ParamDB->SOLD_TYPE==1))
      {
        MatrixRowS[j] += val2 * (test10*ansatz10+test01*ansatz01);
      }
      if((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
         && (TDatabase::ParamDB->SOLD_TYPE==2))
      {
        MatrixRowS[j] += val2 * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01);
      }
    } // endfor j
  } // endfor i
}
void MatrixMRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **Matrix, **MatrixS, *Rhs, val, *MatrixRow, *MatrixSRow;
  double ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c00, c11, c22, c33; 
  double tau, bgradv, bb;
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  Matrix = LocMatrices[0];
  MatrixS = LocMatrices[1];
    
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f
  
  c00 = theta1 * time_step * c0;
  c11 = theta1 * time_step * c1;
  c22 = theta1 * time_step * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + theta1 * time_step * c3;
  if (TDatabase::ParamDB->SDFEM_TYPE==8)
  {
    c33 = theta1 * time_step * c3;
  }
  if(fabs(c11) > fabs(c22))
    bb = fabs(c11);
  else
    bb = fabs(c22);
  // this is \tilde tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);
  // scale appropriately
  //OutPut(tau << " ");
  // do not apply for paper with J. Novo
  if((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)
     &&(TDatabase::ParamDB->SDFEM_TYPE!=11))
    tau *= theta1 * time_step;
  //OutPut(theta1 << " " << tau << " : ");
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    MatrixSRow = MatrixS[i];   
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    Rhs[i] += Mult*(test00+tau*bgradv)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig2[j];

      MatrixRow[j] += Mult * ansatz00*test00;
      MatrixSRow[j] += Mult * ansatz00*bgradv;     
    } // endfor j
  } // endfor i
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
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, bgradv, sigma, norm_b;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);

  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, delta, param, 0, 0, 0);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = (c0+sigma)*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val += delta * (-c0*(ansatz20+ansatz02)
        +c1*ansatz10+c2*ansatz01
        +c3*ansatz00) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
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
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, bgradv, sigma, norm_b,sigma0;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);
  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, delta, param, 0, 0, 0);
  norm_b = c1*c1+c2*c2;

  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==CS99)
  {
    //sigma0 = sigma-delta*hK*sqrt(norm_b)/2;          // this is kappa_sl in [CS99]
    sigma0 = sigma-delta*norm_b;                  // this is kappa_sl in [CS99]
    if (sigma0<0)
      sigma0 = 0;
    sigma -= sigma0;                              // effective orthogonal diffusion
    c0 += sigma0;                                 // effective isotropic diffusion
  }

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

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

      if (norm_b >0)
        val += sigma * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01)/norm_b;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}

/******************************************************************************/
// BilinearAssemble_SD_SOLD
// assembles matrix with SD term and SOLD term
/******************************************************************************/
void BilinearAssemble_SD_SOLD(double Mult, double *coeff, double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, u[5];
  double delta, bgradv, sigma, res, norm_grad_u, norm_b2;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);
  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, delta, NULL, 0, 0, 0);

  // current finite element solution of equation
  u[0] = param[0];                                // u
  u[1] = param[1];                                // u_x
  u[2] = param[2];                                // u_y
  u[3] = param[3];                                // u_xx
  u[4] = param[4];                                // u_yy

  if ((!TDatabase::ParamDB->SOLD_ADJOINT))
  {
    norm_b2 = c1*c1+c2*c2;

    if (norm_b2 > 1e-30)
    {
      sigma /= norm_b2;
    }
    else
      sigma = 0;

    norm_grad_u = sqrt(u[1]*u[1]+u[2]*u[2]);
    res = fabs(-c0 * (param[3]+param[4]) + c1 * param[1] + c2 * param[2] + c3 * param[0] - c4);
    if (norm_grad_u > 1e-30)
    {
      sigma *= res * hK /(2*norm_grad_u);
    }
    else
      sigma = 0;
  }
  else
  {
    //sigma *= hK/2.0;
    sigma *= hK*hK;
  }

  //OutPut("sigma " << sigma << endl);
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = (c1*test10+c2*test01)*delta;

    Rhs[i] += Mult*(test00+bgradv)*c4;

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
      // SUPG term
      val += (-c0*(ansatz20+ansatz02)
        +c1*ansatz10+c2*ansatz01
        +c3*ansatz00) * bgradv;
      // SOLD term
      val += sigma * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01);
      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}


/******************************************************************************/
// BilinearAssemble2LevelLPS_Q0
// assembles matrix and rhs for 2-level LPS with projection on Q0
/******************************************************************************/
void BilinearAssemble2LevelLPS_Q0(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4;
  double lpcoeff, lpexponent, stab_coeff, norm_b, bgradv, uh_proj, uh_proj_cw;
  double lpcoeff_crosswind, stab_coeff_cw, val1, ux, uy;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  // coefficients
  c0 = coeff[0];                                  // eps
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f

  lpcoeff = TDatabase::ParamDB->LP_STREAMLINE_COEFF;
  lpexponent = TDatabase::ParamDB->LP_STREAMLINE_EXPONENT;
  lpcoeff_crosswind = TDatabase::ParamDB->LP_CROSSWIND_COEFF;

  // h^2/epsilon
  stab_coeff = 4*hK*hK/c0;
  norm_b = sqrt(c1*c1 + c2*c2);
  if (norm_b > 1e-10)
  {
    // h/|b|
    val = 2*hK/norm_b;
    if (val < stab_coeff)
      stab_coeff = val;
  }
  stab_coeff *= lpcoeff;
  // projection
  uh_proj = param[2];

  if (TDatabase::ParamDB->LP_CROSSWIND)
  {
    ux = param[0];
    uy = param[1];
    uh_proj_cw = param[3];
    // norm squared of convection
    if (norm_b > 1e-20)
    {
      // h/|b|
      stab_coeff_cw = 2*lpcoeff_crosswind*hK/norm_b;
    }
    else
      stab_coeff_cw = 0.0;
  }

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig0[i];                            // xi derivative
    test01 = Orig1[i];                            // eta derivative
    test00 = Orig2[i];                            // function

    bgradv = c1 * test10 + c2 * test01;
    bgradv *= stab_coeff;
    // assemble rhs
    // quad_weigth * test_function * f
    Rhs[i] += Mult*(test00*c4 + uh_proj * bgradv);
    if (TDatabase::ParamDB->LP_CROSSWIND)
    {
      val1 = -c2 * ux + c1 * uy - uh_proj_cw;
      val = stab_coeff_cw * fabs(val1) * val1 * (-c2 * test10 + c1 * test01);
      Rhs[i] -= Mult*val;
      //OutPut(uh_proj * bgradv << " "<< val << " " << uh_proj * bgradv + val << endl);
    }

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];                        // xi derivative
      ansatz01 = Orig1[j];                        // eta derivative
      ansatz00 = Orig2[j];                        // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      // + LPS term
      val += (c1*ansatz10+c2*ansatz01)*(test00+bgradv);
      // assembel reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;
      // quad weigth
      val *= Mult;

      // update matrix entry
      MatrixRow[j] += val;
    }                                             // endfor j
  }                                               // endfor i
}



// Layton, Polman, SIAM Sci. Comput. 17, 1328 - 1346, 1996
void RhsAssemble_LP96(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, bgradv, u;
  double umin = 0, umax = 1, rho;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  u = param[0];

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);
  rho = hK* TDatabase::ParamDB->SOLD_CONST;
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    val = (test00+delta*bgradv)*c4;
    val = val - (MIN(u-umin,0)+MAX(u-umax,0))* test00/rho;
    Rhs[i] += Mult*val;
  }                                               // endfor i
}

/******************************************************************************/
// RhsAssemble_RhsAdjointEnergyEstimate
// assemble rhs for adjoint problem with energy error estimator of Verfuerth
// (2005)
/******************************************************************************/

void RhsAssemble_RhsAdjointEnergyEstimate(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4;
  double test, ansatz, scal, u[5];

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f
  // current finite element solution of equation
  u[0] = param[0]; // u
  u[1] = param[1]; // u_x
  u[2] = param[2]; // u_y
  u[3] = param[3]; // u_xx
  u[4] = param[4]; // u_yy

  // residual eps Delta u - b* nabla u - c u + f
  ansatz = c0 * (u[3] + u[4]) - c1*u[1] - c2*u[2] - c3*u[0] + c4;

  // compute scaling of strong residual in energy error estimator
  scal = hK*hK/c0;
  if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
  {
    // update weight for energy norm estimator
    if (1.0/TDatabase::ParamDB->INTERNAL_COERCIVITY<scal)
      scal = 1.0/TDatabase::ParamDB->INTERNAL_COERCIVITY; 
  }
  // scal ansatz
  ansatz *= 2 * scal;

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    test20 = Orig3[i];
    test02 = Orig4[i];
    // test function
    test = c0 * (test20 + test02) - c1*test10 - c2*test01 - c3*test00;
    // update array for rhs
    val = ansatz * test;
    Rhs[i] += Mult*val;
  }                                               // endfor i
}

/******************************************************************************/
// RhsAssemble_RhsAdjointTV2
// assemble rhs for adjoint problem with total variation squared
/******************************************************************************/

void RhsAssemble_RhsAdjointTV2(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4;
  double test, ansatz, scal, u[5];

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y

  // current finite element solution of equation
  u[1] = param[1];                                // u_x
  u[2] = param[2];                                // u_y

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];

    // update array for rhs
    val = u[1]*test10 + u[2]*test01;
    // scaling 2 from derivative of the square
    Rhs[i] += 2*Mult*val;
  }                                               // endfor i
}


/******************************************************************************/
// RhsAssemble_RhsAdjointTV
// assemble rhs for adjoint problem with total variation
/******************************************************************************/
void RhsAssemble_RhsAdjointTV(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, norm_grad;
  double test, ansatz, scal, u[5], eps=1e-12;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y

  // current finite element solution of equation
  u[1] = param[1];                                // u_x
  u[2] = param[2];                                // u_y

  norm_grad = sqrt(u[1]*u[1]+u[2]*u[2]);

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];

    // update array for rhs
    val = u[1]*test10 + u[2]*test01;
    if (fabs(norm_grad)>1e-30)
      val /= norm_grad;
    else
    {
      val = 0;
    }

    Rhs[i] += Mult*val;
  }                                               // endfor i
}


/******************************************************************************/
// RhsAssemble_RhsAdjointNormBL1_NormBorthL1
// assemble rhs for adjoint problem with total variation plus crosswind derivative
/******************************************************************************/
void RhsAssemble_RhsAdjointNormBL1_NormBorthL1(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, norm_b, b_gradu, sig_b_gradu, borth_gradu, sig_borth_gradu;
  double test, ansatz, scal, u[5], eps=1e-12;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y

  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2

  // current finite element solution of equation
  u[1] = param[1];                                // u_x
  u[2] = param[2];                                // u_y

  norm_b = sqrt(c1*c1+c2*c2);
  b_gradu = c1*u[1] + c2*u[2];
  sig_b_gradu = SIGN(b_gradu);
  borth_gradu = -c2*u[1] + c1*u[2];
  if (norm_b>1e-30)
    sig_borth_gradu = SIGN(borth_gradu)/norm_b;
  else
    sig_borth_gradu = 0;
  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];

    // update array for rhs
    val = (c1*test10 + c2*test01)*sig_b_gradu;
    val += (-c2*test10+c1*test01)* sig_borth_gradu;

    Rhs[i] += Mult*val;
  }                                               // endfor i
}


/******************************************************************************/
// RhsAssemble_RhsAdjointNormResidualL1_NormBorthL1
// assemble rhs for adjoint problem with total variation plus crosswind derivative
/******************************************************************************/
void RhsAssemble_RhsAdjointNormResidualL1_NormBorthL1(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, res, sig_res, borth_gradu, sig_borth_gradu, norm_b;
  double test, ansatz, scal, u[5], eps=1e-12;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // u_xx
  Orig4 = OrigValues[4];                          // u_yy

  c0 = coeff[0];                                  // eps
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  // current finite element solution of equation
  u[0] = param[0];                                // u
  u[1] = param[1];                                // u_x
  u[2] = param[2];                                // u_y
  u[3] = param[3];                                // u_xx
  u[4] = param[4];                                // u_yy

  res = -c0 * (u[3] + u[4]) + c1*u[1] + c2*u[2] + c3 * u[0] - c4;
  sig_res = SIGN(res);
  borth_gradu = -c2*u[1] + c1*u[2];
  norm_b = sqrt(c1*c1+c2*c2);
  if (norm_b>1e-30)
    sig_borth_gradu = SIGN(borth_gradu)/norm_b;
  else
    sig_borth_gradu = 0;
  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    test20 = Orig3[i];
    test02 = Orig4[i];

    // update array for rhs
    val =  (-c0 * (test20 + test02) + c1*test10 + c2*test01 + c3*test00) * sig_res;
    val += (-c2*test10+c1*test01)* sig_borth_gradu;

    Rhs[i] += Mult*val;
  }                                               // endfor i
}


/******************************************************************************/
// RhsAssemble_RhsAdjointEnergyEstimate
// assemble rhs for adjoint problem with energy error estimator of Verfuerth
// (2005)
/******************************************************************************/

void RhsAssemble_RhsAdjointAll(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val, aux_val;
  double alpha1, alpha2, beta1, beta2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  double residual, sig_residual, norm_grad, b_gradu, sig_b_gradu, borth_gradu;
  double norm_b, sig_borth_gradu, norm_borth_gradu_sqrt;
  double normb, p = TDatabase::ParamDB->RESIDUAL_LP_ADJOINT;
  double valres, valdres, valcw, valdcw;
  int i,j, N_;
  double c0, c1, c2, c3, c4, weight, c1_norm, c2_norm;
  double test, ansatz, scal, u[5];

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u_x
  Orig1 = OrigValues[1];                          // u_y
  Orig2 = OrigValues[2];                          // u
  Orig3 = OrigValues[3];                          // u_xx
  Orig4 = OrigValues[4];                          // u_yy

  c0 = coeff[0];                                  // eps
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  // current finite element solution of equation
  u[0] = param[0];                                // u
  u[1] = param[1];                                // u_x
  u[2] = param[2];                                // u_y
  u[3] = param[3];                                // u_xx
  u[4] = param[4];                                // u_yy

  // convection
  norm_b = sqrt(c1*c1+c2*c2);

  // residual eps Delta u - b* nabla u - c u + f
  residual = c0 * (u[3] + u[4]) - c1*u[1] - c2*u[2] - c3*u[0] + c4;
  sig_residual = SIGN(residual);

  // norm of gradient
  norm_grad = sqrt(u[1]*u[1]+u[2]*u[2]);

  // streamline derivative
  b_gradu = c1*u[1] + c2*u[2];
  sig_b_gradu = SIGN(b_gradu);

  // normalize convection
  if (norm_b > 1e-30)
  {
    c1_norm = c1/norm_b;
    c2_norm = c2/norm_b;
  }
  else
    c1_norm = c2_norm = 0;

  // normalized derivative orthogonal to streamline
  borth_gradu = -c2_norm*u[1] + c1_norm*u[2];
  sig_borth_gradu = SIGN(borth_gradu);
  norm_borth_gradu_sqrt = sqrt(fabs(borth_gradu));

  alpha1 = TDatabase::ParamDB->MIN_MAX_EXPONENT_ONE_ADJOINT;
  alpha2 = TDatabase::ParamDB->MIN_MAX_EXPONENT_TWO_ADJOINT;
  beta1 = TDatabase::ParamDB->MIN_MAX_FACTOR_ONE_ADJOINT;
  beta2 = TDatabase::ParamDB->MIN_MAX_FACTOR_TWO_ADJOINT;

  i = TDatabase::ParamDB->INTERNAL_CELL;
  weight = TDatabase::ParamDB->INTERNAL_WEIGHT_Array[i];
  //weight *= TDatabase::ParamDB->INTERNAL_WEIGHT_NODAL_Array[i];
  // weights are computed in MainRoutines2D.C at line 765 and the following

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    test20 = Orig3[i];
    test02 = Orig4[i];
    // test function
    test = c0 * (test20 + test02) - c1*test10 - c2*test01 - c3*test00;
    // update array for rhs
    valdres = 2 * residual * test;
    val = TDatabase::ParamDB->WEIGHT_RESIDUAL_L2_ADJOINT * valdres;
    val += TDatabase::ParamDB->WEIGHT_GRADIENT_L2_ADJOINT * 2 * (u[1]*test10 + u[2]*test01);
    if (fabs(norm_grad)>1e-30)
      val += TDatabase::ParamDB->WEIGHT_GRADIENT_L1_ADJOINT * (u[1]*test10 + u[2]*test01)/norm_grad;
    val += TDatabase::ParamDB->WEIGHT_STREAM_DER_L1_ADJOINT*(c1*test10 + c2*test01)*sig_b_gradu;
    val += TDatabase::ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT*(-c2_norm*test10+c1_norm*test01)* sig_borth_gradu;
    val += TDatabase::ParamDB->WEIGHT_RESIDUAL_L1_ADJOINT *
      (c0 * (test20 + test02) - c1*test10 - c2*test01 - c3*test00) * sig_residual;
    if (p>=1.0)
      val += TDatabase::ParamDB->WEIGHT_RESIDUAL_LP_ADJOINT * p *pow (fabs(residual),p-1) * sig_residual * test;
    if (TDatabase::ParamDB->WEIGHT_RESIDUAL_EXP_LP_ADJOINT > 0)
      val += TDatabase::ParamDB->WEIGHT_RESIDUAL_EXP_LP_ADJOINT * (exp(pow(fabs(residual),p))-1)
        * p *pow (fabs(residual),p-1) * sig_residual * test;
    if (norm_borth_gradu_sqrt<=1)
    {
      valdcw = 5 * borth_gradu * (-c2_norm*test10+c1_norm*test01)
        - 9/2.0 * sig_borth_gradu * borth_gradu * borth_gradu * (-c2_norm*test10+c1_norm*test01);
      val += TDatabase::ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT*valdcw;
    }
    else
    {
      if  (norm_borth_gradu_sqrt > 1e-30)
      {
        valdcw = (-c2_norm*test10+c1_norm*test01)
          * sig_borth_gradu/(2*norm_borth_gradu_sqrt);
        val += TDatabase::ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT*valdcw;
      }
      else
      {
        valdcw = 0;
      }
    }

    if (TDatabase::ParamDB->WEIGHT_RESIDUAL_CW_ADJOINT > 0)
    {
      if (fabs(residual) < TDatabase::ParamDB->P6)
        valres = residual * residual;
      else
      {
        valres = TDatabase::ParamDB->P6 * TDatabase::ParamDB->P6;
        valdres = 0.0;
      }
      if (norm_borth_gradu_sqrt<=1)
      {
        valcw = 2.5 *norm_borth_gradu_sqrt * norm_borth_gradu_sqrt
          - 1.5*norm_borth_gradu_sqrt * norm_borth_gradu_sqrt * norm_borth_gradu_sqrt;
      }
      else
      {
        valcw = sqrt(norm_borth_gradu_sqrt);
      }
      val += TDatabase::ParamDB->WEIGHT_RESIDUAL_CW_ADJOINT*(valcw * valdres + valres * valdcw);
    }

    if(TDatabase::ParamDB->MIN_MAX_ADJOINT == 3 ||
      TDatabase::ParamDB->MIN_MAX_ADJOINT == 1)
    {
      aux_val = TDatabase::ParamDB->MIN_VAL_ADJOINT - u[0];
      if( aux_val > 0 )                           // values below prescribed minimum
      {
        //if(2*beta2*pow(aux_val,1+alpha2) * beta2*(1+alpha2)*pow(aux_val,alpha2) * test00/val>10 ||
        //   2*beta2*pow(aux_val,1+alpha2) * beta2*(1+alpha2)*pow(aux_val,alpha2) * test00/val<0.1)
        //   OutPut("Unterschwinger " << aux_val << " increase/val "
        //     << 2*beta2*pow(aux_val,1+alpha2) *
        //       beta2*(1+alpha2)*pow(aux_val,alpha2) * test00/val << endl);
        val -= 2*beta2*pow(aux_val,1+alpha2) *
          beta2*(1+alpha2)*pow(aux_val,alpha2)*2* aux_val* test00;
      }
      aux_val = u[0] - TDatabase::ParamDB->MAX_VAL_ADJOINT;
      if( aux_val > 0 )                           // values above prescribed maximum
        val += 2*beta2*pow(aux_val,1+alpha2) *
          beta2*(1+alpha2)*pow(aux_val,alpha2) * 2 * aux_val * test00;
    }
    if(TDatabase::ParamDB->MIN_MAX_ADJOINT == 2 ||
      TDatabase::ParamDB->MIN_MAX_ADJOINT == 1)
    {
      aux_val = TDatabase::ParamDB->MIN_VAL_ADJOINT - u[0];
      if( aux_val > 0 )                           // values below prescribed minimum
      {
        //OutPut("Unterschwinger " << aux_val << " increase/val "
        //  << 2*beta1*(exp(alpha1*aux_val*aux_val)-1) *
        //    beat1*2*alpha1*aux_val*exp(alpha1*aux_val*aux_val) * test00
        //  << " val before " << val << endl );
        val -=            2*beta1*(exp(alpha1*aux_val*aux_val)-1) *
          beta1*2*alpha1*aux_val*exp(alpha1*aux_val*aux_val) * test00;
      }
      aux_val = u[0] - TDatabase::ParamDB->MAX_VAL_ADJOINT;
      if( aux_val > 0 )                           // values above prescribed maximum
        val +=            2*beta1*(exp(alpha1*aux_val*aux_val)-1) *
          beta1*2*alpha1*aux_val*exp(alpha1*aux_val*aux_val) * test00;
    }

    //switch(TDatabase::ParamDB->MIN_MAX_ADJOINT)
    //{
    //  case 0:
    //      aux_val = TDatabase::ParamDB->MIN_VAL_ADJOINT - u[0];
    //      if( aux_val > 0 ) // values below prescribed minimum
    //      {
    //        OutPut("Unterschwinger " << aux_val << " increase/val "
    //            << 2*pow(aux_val,1+exponent) *
    //              (1+exponent)*pow(aux_val,exponent) * test00 *
    //              TDatabase::ParamDB->MIN_MAX_FACTOR_ADJOINT/val
    //            << endl);
    //        val += 2*pow(aux_val,1+exponent) *
    //          (1+exponent)*pow(aux_val,exponent) * test00 *
    //          TDatabase::ParamDB->MIN_MAX_FACTOR_ADJOINT;
    //      }
    //      aux_val = u[0] - TDatabase::ParamDB->MAX_VAL_ADJOINT;
    //      if( aux_val > 0 ) // values above prescribed maximum
    //        val += 2*pow(aux_val,1+exponent) *
    //          (1+exponent)*pow(aux_val,exponent) * test00 *
    //          TDatabase::ParamDB->MIN_MAX_FACTOR_ADJOINT;
    //    break;
    //  case 1:
    //      aux_val = TDatabase::ParamDB->MIN_VAL_ADJOINT - u[0];
    //      if( aux_val > 0 ) // values below prescribed minimum
    //      {
    //        //OutPut("Unterschwinger " <<
    //        //  2*(exp(exponent*aux_val*aux_val)-1) * 2*exponent*aux_val*exp(exponent*aux_val*aux_val)*
    //        //  test00 << " val before " << val << endl );
    //        val +=              2*(exp(exponent*aux_val*aux_val)-1) *
    //            2*exponent*aux_val*exp(exponent*aux_val*aux_val) * test00;
    //      }
    //      aux_val = u[0] - TDatabase::ParamDB->MAX_VAL_ADJOINT;
    //      if( aux_val > 0 ) // values above prescribed maximum
    //        val +=              2*(exp(exponent*aux_val*aux_val)-1) *
    //            2*exponent*aux_val*exp(exponent*aux_val*aux_val) * test00;
    //    break;
    //}
    val *= weight;
    Rhs[i] += Mult*val;
  }                                               // endfor i
}


/******************************************************************************/
// RhsAssemble_RhsAdjointL2Error
// assemble rhs for adjoint problem with L2 error to prescribed solution
/******************************************************************************/
void RhsAssemble_RhsAdjointL2Error(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double ansatz, test, u, u_h;
  double *Orig0;
  int i, N_;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  // current finite element solution of equation
  u_h = param[0]; // u

  // prescribed solution
  u = coeff[10];

  // difference in quadrature point
  ansatz = -2 * (u - u_h);

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test = Orig0[i];
    Rhs[i] += Mult * ansatz * test;
  }                                               // endfor i
}

/******************************************************************************/
// RhsAssemble_RhsAdjointH1Error
// assemble rhs for adjoint problem with H1-semi norm error to prescribed solution
/******************************************************************************/
void RhsAssemble_RhsAdjointH1Error(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double test10, test01;
  double ansatz_x, ansatz_y;
  double *Orig0, *Orig1;
  int i,j, N_;
  double scal, u[2];

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y

  // current finite element solution of equation
  u[0] = param[1]; // u_x
  u[1] = param[2]; // u_y

  // error of first partial derivatives
  ansatz_x = coeff[11] - u[0];
  ansatz_y = coeff[12] - u[1]; 
  scal = -2; //savescu
  ansatz_x *= scal;
  ansatz_y *= scal;

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    // test function
    test10 = Orig0[i];
    test01 = Orig1[i];

    // update array for rhs
    val = ansatz_x * test10 + ansatz_y * test01;
    Rhs[i] += Mult*val;
  }                                               // endfor i
}


/******************************************************************************/
//
// IMPROVED MIZUKAMI-HUGHES METHOD (Knobloch, CMAME 2007)
//
// assembling only Galerkin part
//
/******************************************************************************/

void BilinearAssemble_MH_Kno06(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

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
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}


// ========================================================================
// NEEDED PARAM FUNCTIONS START
// ========================================================================
// ========================================================================

// ========================================================================
// parameters:  H1 norm and  norm of residual
// ========================================================================
void DC_CD_Params(double *in, double *out)
{
  out[0] = in[2];                                 // H1 norm
  out[1] = in[3];                                 // norm of residual
}



// ========================================================================
// parameters:  partial derivatives
// ========================================================================
void Params_Sol(double *in, double *out)
{
  out[0] = in[2];                                 // u
}

// ========================================================================
// parameters:  SC_2
// ========================================================================
void SC_2_Params(double *in, double *out)
{
  out[0] = in[2];                                 // H1 norm
  out[1] = in[3];                                 // norm of residual
  out[2] = in[4];                                 // u_x
  out[3] = in[5];                                 // u_y
}

// ========================================================================
// parameters:  SOLD
// ========================================================================
void SOLD_Params(double *in, double *out)
{
  out[0] = in[2];                                 // u
  out[1] = in[3];                                 // u_x
  out[2] = in[4];                                 // u_y
  out[3] = in[5];                                 // u_xx
  out[4] = in[6];                                 // u_yy
  out[5] = in[7];                                 // ||u^h||_{H^1,K}
  out[6] = in[8];                                 // ||R(u^h)||_{L^2,K}
}


// ========================================================================
// parameters: velocity field
// ========================================================================
void Params_Velo(double *in, double *out)
{
  out[0] = in[2];                                 // velo_1
  out[1] = in[3];                                 // velo_2
}

// ========================================================================
// parameters:  SOLD + velocity field
// ========================================================================
void SOLD_Params_And_Velo(double *in, double *out)
{
  out[0] = in[2];                                 // u
  out[1] = in[3];                                 // u_x
  out[2] = in[4];                                 // u_y
  out[3] = in[5];                                 // u_xx
  out[4] = in[6];                                 // u_yy
  out[5] = in[7];                                 // ||u^h||_{H^1,K}
  out[6] = in[8];                                 // ||R(u^h)||_{L^2,K}
  out[7] = in[9];                                 // velo_1
  out[8] = in[10];                                // velo_2
}



// ========================================================================
// parameters: two fe values
// ========================================================================
void Params_Sol4(double *in, double *out)
{
  out[0] = in[2];                                 // u_x
  out[1] = in[3];                                 // u_y
  out[2] = in[4];
  out[3] = in[5];
}

// ========================================================================
// ========================================================================
// NEEDED PARAM FUNCTIONS END
// ========================================================================
// ========================================================================

#endif // #ifdef __2D__
