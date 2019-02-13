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
// @(#)TimeConvDiff2D.C        1.4 04/13/00
//
// common declaration for time dependent convection diffusion problems
// ======================================================================
#include <Database.h>
#include <ConvDiff.h>

void TimeBilinearAssemble(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0];
  c1 = coeff[1];
  c2 = coeff[2];
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
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val *=Mult;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

void TimeMassAssemble(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];
      
      MatrixRow[j] += Mult*ansatz00*test00;
    } // endfor j
  } // endfor i
}

void TimeBilinearAssemble_SD(double Mult, double *coeff, double hK, 
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


void TimeRhsAssembleRB(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4, c5; 

  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c4 = coeff[4]; // Q

  for(i=0;i<N_;i++)
  {
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c4;

  } // endfor i
}


void MatrixMARhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs)
{
  int i,j, N_T;
  
  double c0, c1, c2, c3, c4, x, r;  
  double **MatrixA, **MatrixM, *Rhs;
  double *Orig0, *Orig1, *Orig2;  
  double *MatrixARow, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, val;
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  
  Rhs = LocRhs[0];

  N_T = N_BaseFuncts[0]; 
  
  Orig0 = OrigValues[0]; // T_x
  Orig1 = OrigValues[1]; // T_y
  Orig2 = OrigValues[2]; // T

  c0 = coeff[0]; // nu
//   c1 = coeff[1];  
//   c2 = coeff[2];  
  c3 = coeff[3]; // concentration coeff
  c4 = coeff[4]; // rhs
 
  c1 = param[0]; // u1-w1
  c2 = param[1]; // u2-w2
 
  x  = param[2]; // x
  r  = fabs(x);  
   
  //cout<< "x: "<< x<<" u1: "<< u1<< " u2: "<< u2<<endl;
    
  if(r<1e-12)
   {
   OutPut("check NSE2D Axisymmetric file x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   

  for(i=0;i<N_T;i++)
   {
    MatrixARow = MatrixA[i];
    MatrixMRow  = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // rhs
    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_T;j++)
     {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];    
    
      //stiffness mat
      val  = r*c0*(test10*ansatz10+test01*ansatz01);
      val  += r* ((c1*ansatz10+c2*ansatz01)*test00);
      val += r*c3*ansatz00*test00;
  
      MatrixARow[j] += val*Mult;

      // mass mat
      val = r*Mult*ansatz00*test00;      
      MatrixMRow[j] += val;
      
     } // for(j=0;j<N_T;j++)     
   } //   for(i=0;i<N_T;i++)
   
}// MatrixMARhsAssemble_Axial3D

void MatrixMARhsAssemble(double Mult, double *coeff, double *param,
                         double hK, double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs)
{
  int i,j, N_T;
  
  double c0, c1, c2, c3, c4;  
  double **MatrixA, **MatrixM, *Rhs;
  double *Orig0, *Orig1, *Orig2;  
  double *MatrixARow, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, val;
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  
  Rhs = LocRhs[0];

  N_T = N_BaseFuncts[0]; 
  
  Orig0 = OrigValues[0]; // T_x
  Orig1 = OrigValues[1]; // T_y
  Orig2 = OrigValues[2]; // T

  c0 = coeff[0]; // nu
  c1 = coeff[1] - param[0]; // u1-w1  
  c2 = coeff[2] - param[1]; // u2-w2  
  c4 = coeff[4]; // rhs  
  
  if(TDatabase::ParamDB->P6==1) // con-ALE
   {
    c3 = coeff[3] - param[3] ; // concentration coeff - div w  
    
    if(TDatabase::ParamDB->REACTOR_P29>param[3]) 
       { TDatabase::ParamDB->REACTOR_P29=param[3]; }

     if(TDatabase::ParamDB->REACTOR_P30<param[3]) 
       {  TDatabase::ParamDB->REACTOR_P30=param[3]; }     
      
   }
  else// non-conservative ALE
   {
    c3 = coeff[3]  ; 
    
/*    if(TDatabase::ParamDB->REACTOR_P29>param[3]) 
       { TDatabase::ParamDB->REACTOR_P29=param[3]; }

     if(TDatabase::ParamDB->REACTOR_P30<param[3]) 
       {  TDatabase::ParamDB->REACTOR_P30=param[3]; }   */    
//     cout<< "div w: "<< c3  <<endl;
   }
//   cout<< " u1: "<< u1<< " u2: "<< u2<<endl;
//  cout<< "div w: "<< param[3]  <<endl;

  for(i=0;i<N_T;i++)
   {
    MatrixARow = MatrixA[i];
    MatrixMRow  = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // rhs
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_T;j++)
     {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];    
    
      //stiffness mat
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;
   
      MatrixARow[j] += val*Mult;

      // mass mat
      val = Mult*ansatz00*test00;      
      MatrixMRow[j] += val;
      
     } // for(j=0;j<N_T;j++)     
   } //   for(i=0;i<N_T;i++)
   
}// MatrixMARhsAssemble

void MatrixMARhsALEAssemble_SUPG(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs)
{
  int i,j, N_T;
  
  double c0, c1, c2, c3, c4;  
  double **MatrixA, **MatrixM, **MatrixK, *Rhs;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;  
  double *MatrixARow, *MatrixMRow, *MatrixKRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, val, val1;
  double c00, c11, c22, c33;
  double tau, bgradv, bb, res, sigma, norm_b;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixK = LocMatrices[2];
  
  Rhs = LocRhs[0];

  N_T = N_BaseFuncts[0]; 
  
  Orig0 = OrigValues[0]; // T_x
  Orig1 = OrigValues[1]; // T_y
  Orig2 = OrigValues[2]; // T
  Orig3 = OrigValues[3]; // 
  Orig4 = OrigValues[4]; //

  c0 = coeff[0]; // nu
  c1 = coeff[1] - param[0]; // u1-w1  
  c2 = coeff[2] - param[1]; // u2-w2  
  
  if(TDatabase::ParamDB->P6==1) // con-ALE
   {
    c3 = coeff[3] - param[3]; // concentration coeff - div w  
    
    if(TDatabase::ParamDB->REACTOR_P29>param[3]) 
       { TDatabase::ParamDB->REACTOR_P29=param[3]; }

     if(TDatabase::ParamDB->REACTOR_P30<param[3]) 
       {  TDatabase::ParamDB->REACTOR_P30=param[3]; }     
          
   }
  else// non-conservative ALE
   {
    c3 = coeff[3]; 
   }  
  
  c4 = coeff[4]; // rhs

//  cout<< " u1: "<< u1<< " u2: "<< u2<<endl;
//  cout<< "div w: "<< param[3]  <<endl;

  // coefficients for the stabilization parameter
  val = theta1 * time_step;
  c00 = val * c0;
  c11 = val * c1;
  c22 = val * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + val * c3;
  // reactive coefficient, inclusive term from the temporal derivative
//   c33 = 1.0/(val) + c3;  
  
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
  // this is tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  for(i=0;i<N_T;i++)
   {
    MatrixARow = MatrixA[i];
    MatrixMRow  = MatrixM[i];
    MatrixKRow  = MatrixK[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    // scaling with the stabilization parameter
    // scaling with the time step is done in the main program
    bgradv *= tau;
   
    // rhs
    Rhs[i] += Mult*(test00 + bgradv)*c4;

    for(j=0;j<N_T;j++)
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
      val += val1*(test00 + bgradv);
      
      // diffusion part of the SUPG stabilization
      val -= c0*(ansatz20 + ansatz02) * bgradv; 
        
      MatrixARow[j] += val*Mult;

      // mass mat
      val = Mult*ansatz00*test00;      
      MatrixMRow[j] += val;
           
      // time consistant term
      val = Mult*ansatz00*bgradv;      
      MatrixKRow[j] += val;   
//             cout <<  val << " : " ;
	    
     } // for(j=0;j<N_T;j++)     
   } //   for(i=0;i<N_T;i++)
   
}// MatrixMARhsALEAssemble_SUPG

void MatrixMRhsALEAssemble_SUPG(double Mult, double *coeff, double *param,
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
  c1 = coeff[1] - param[0]; // u1-w1  
  c2 = coeff[2] - param[1]; // u2-w2  
  
  if(TDatabase::ParamDB->P6==1) // con-ALE
   {
    c3 = coeff[3] - param[3]; // concentration coeff - div w  
   }
  else// non-conservative ALE
   {
    c3 = coeff[3]; 
   }  
  
  c4 = coeff[4]; // rhs

    // coefficients for the stabilization parameter
  val = theta1 * time_step;
  c00 = val * c0;
  c11 = val * c1;
  c22 = val * c2;
  // reactive coefficient, inclusive term from the temporal derivative
  c33 = 1.0 + val * c3;
  // reactive coefficient, inclusive term from the temporal derivative
//   c33 = 1.0/(val) + c3;  
     
  
  if (fabs(c11) > fabs(c22))
      bb = fabs(c11);
  else
      bb = fabs(c22);
 // this is \tilde tau
 tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, bb);

  // scale appropriately
  //OutPut(tau << " ");
  // do not apply for paper with J. Novo
//   if((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)
//      &&(TDatabase::ParamDB->SDFEM_TYPE!=11))
//     tau *= theta1 * time_step;
//   OutPut(theta1 << " " << tau << " : ");
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    MatrixSRow = MatrixS[i];   
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    // scaling with the stabilization parameter
    bgradv *= tau;
    
    Rhs[i] += Mult*(test00+bgradv)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig2[j];

      MatrixRow[j] += Mult * ansatz00*test00;
      MatrixSRow[j] +=  Mult*ansatz00*bgradv;
    } // endfor j
  } // endfor i
}

void MatrixARhsAssembleHeatLine(double Mult, double *coeff, double *param,
                                 double hK, 
                                 double **OrigValues, int *N_BaseFuncts,
                                 double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c1, c2, c3, c4, c5, c6, c7, fact; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

      
  c1 = param[0]; // T
  c2 = param[1]; // T_x
  c3 = param[2]; // T_y
  c4 = param[3]; // u1-w1
  c5 = param[4]; // u2-w2
  c6 = param[5]; // u1y-w1y
  c7 = param[6]; // u2x-w2x
     
  fact = c5*c2 + c1*c7 - c4*c3 -c1*c6;

//     for(i=0;i<7;i++)
//    cout<< "fact: "<< param[i] <<endl;
//    cout <<endl;
   
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*fact;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = Mult*(test10*ansatz10+test01*ansatz01);
      MatrixRow[j] += val;
    } // endfor j
  } // endfor i

}


void MatrixARhsAssembleHeatLine_Axial3D(double Mult, double *coeff, double *param, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c1, c2, c3, c4, c5, c6, c7, fact, x, r; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];


   
  c1 = param[0]; // T
  c2 = param[1]; // T_x
  c3 = param[2]; // T_y
  c4 = param[3]; // u1-w1
  c5 = param[4]; // u2-w2
  c6 = param[5]; // u1y-w1y
  c7 = param[6]; // u2x-w2x
  fact = c5*c2 + c1*c7 - c4*c3 -c1*c6; 
  
  x  = param[7]; // x
  r  = fabs(x);  
  Mult *=r; 
  
//   cout<< "fact: "<< fact<<endl;
  
  if(r<1e-12)
   {
   OutPut("check heat Line Axisymmetric file x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

//     Rhs[i] += Mult*test00*fact;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
//       ansatz00 = Orig2[i];
      
      val = Mult*(test10*ansatz10+test01*ansatz01);
//       cout << "val: "<< val<<endl;
      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}
