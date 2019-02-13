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
//  TCD3D.C 
// ======================================================================

#include <Database.h>
#include <ConvDiff.h>
#include <math.h>
#include <stdlib.h>
#include <Enumerations.h>

void MatrixMRhsAssemble(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, *MatrixRow;
  double ansatz000;
  double test000;
  double *Orig0;
  int i,j, N_;
  double c5; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  c5 = coeff[5]; // f

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test000 = Orig0[i];

    Rhs[i] += Mult*test000*c5;

    for(j=0;j<N_;j++)
    {
      ansatz000 = Orig0[j];

      MatrixRow[j] += Mult*ansatz000*test000;
    } // endfor j
  } // endfor i
}

void MatrixMRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, Pe; 
  double tau, bgradv, bb;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // b_3
  c4 = coeff[4]; // c
  c5 = coeff[5]; // f
  
  bb  = sqrt(c1*c1+c2*c2+c3*c3);
  tau  = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c4, bb);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    Rhs[i] += Mult*(test000+tau*bgradv)*c5;

    for(j=0;j<N_;j++)
    {
      ansatz000 = Orig3[j];

      MatrixRow[j] += Mult * ansatz000*test000;
    } // endfor j
  } // endfor i
}

void MatricesAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, *MatrixRowA, *MatrixRowK;
  double **MatrixS, *MatrixRowS;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double val, val1, val2;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, c6, c00, c11, c22, c33, c44, c12, c13, c23; 
  double tau, bgradv, bb, res, sigma, norm_b;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  MatrixA = LocMatrices[0];
  MatrixK= LocMatrices[1];
  MatrixS= LocMatrices[2];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // b_3
  c4 = coeff[4]; // c
  c5 = coeff[5]; // f
  
  val = theta1 * time_step;
  c00 = val * c0;
  c11 = val * c1;
  c22 = val * c2;
  c33 = val * c3;
  // reactive coefficient, inclusive term from the temporal derivative
  c44 = 1.0 + val * c4;
  bb = fabs(c11);
  if (fabs(c22)>bb)
      bb = fabs(c22);
  if (fabs(c33)>bb)
      bb = fabs(c33);

  // this is \tilde tau
  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, c44, bb);

  if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
  {
      // rhs from previous time step 
      c6 = coeff[6];   
      // compute residual
      res = param[0] + val*(c1*param[1] +c2*param[2] + c3*param[3] + c4*param[0])
	  -param[4] +  theta2*time_step*(c1*param[5] +c2*param[6] + c3*param[7] + c4*param[4])
	  -theta3*time_step*c6 - theta4*time_step*c5;
      c6 = time_step * theta4 * c5;
      // compute the parameter, c6 is just a dummy
      sigma = Compute_SOLD_sigma(hK, c00, c11, c22, c33, c44, c6, bb, tau, param, res, 1,1);
      //OutPut( param[0] << " " << param[4] <<  " " << res << " " <<  sigma << endl);
      val2 = Mult * sigma;
      if (TDatabase::ParamDB->SOLD_TYPE==2)
      {
	  c11 = c1 * c1;
	  c22 = c2 * c2;
	  c33 = c3 * c3;
	  c12 = c1 * c2;
	  c13 = c1 * c3;
	  c23 = c2 * c3;
 	  norm_b = c11+c22+c33;
	  if (norm_b >1e-10)
	      val2 /= norm_b;
	  else
	      val2 = 0.0;
     }
  }
  
  // scale appropriately, after it is used for the SOLD scheme
  tau *= val;

  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    MatrixRowK = MatrixK[i];
    if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
	MatrixRowS = MatrixS[i];
	
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;
    bgradv *= tau;
    // CHANGE TEST000 !
    test000 += bgradv;
    Rhs[i] += Mult*test000*c5;


    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val1 = c1*ansatz100+c2*ansatz010+c3*ansatz001;
      val1 +=  c4*ansatz000;

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += val1* test000;

      MatrixRowA[j] += Mult * val;
                
      MatrixRowK[j] += Mult * ansatz000*bgradv;

      // isotropic SOLD method
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==1))
      {
	  MatrixRowS[j] += val2 * (test100*ansatz100+test010*ansatz010+ansatz001*test001);
      }
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==2))
      {
	  MatrixRowS[j] += val2 * ( (c22+c33)*ansatz100*test100 + (c11+c33)*ansatz010*test010
					    + (c11+c22)*ansatz001*test001
					    -c13*(ansatz001*test100+ansatz100*test001)
					    -c12*(ansatz010*test100+ansatz100*test010)
					    -c23*(ansatz001*test010+ansatz010*test001));
      }
    } // endfor j
  } // endfor i
}

void MatrixARhsAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, val, *MatrixRowA, *MatrixRowK;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, Pe; 

  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // b_3
  c4 = coeff[4]; // c
  c5 = coeff[5]; // f
  
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
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
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      MatrixRowA[j] += Mult * val;
                
    } // endfor j
  } // endfor i
}
void MatrixAUpwindRhsAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, val, *MatrixRowA, *MatrixRowK;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c4, c5; 

  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c4 = coeff[4]; // c
  c5 = coeff[5]; // f
  
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
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

      MatrixRowA[j] += Mult * val;
                
    } // endfor j
  } // endfor i
}
void RhsAssemble_SUPG(double Mult, double *coeff, double *param,
                      double hK, 
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c5, Pe; 
  double tau, bgradv, bb;

  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // b_3
  c5 = coeff[5]; // f
  
  bb  = sqrt(c1*c1+c2*c2+c3*c3);
  Pe = hK*bb / (2*c0);
  if(Pe<1e-2)
    tau = Pe/(6*bb);
  else
    tau = (1/tanh(Pe)-1/Pe)/(2*bb);

  tau *= hK;

  for(i=0;i<N_;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    Rhs[i] += Mult*(test000+tau*bgradv)*c5;
  } // endfor i
}

void RhsAssemble(double Mult, double *coeff, double *param,
                 double hK, 
                 double **OrigValues, int *N_BaseFuncts,
                 double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double test000;
  double *Orig0;
  int i,j, N_;
  double c5; 

  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  c5 = coeff[5]; // f
  
  for(i=0;i<N_;i++)
  {
    test000 = Orig0[i];

    Rhs[i] += Mult*test000*c5;
  } // endfor i
}


// ======================================================================
//  definitions for assembling the matrix M, A group fem-fct
// ======================================================================
// int N_Terms_MatrixMAGroupFEMRhs = 4;
// MultiIndex3D Derivatives_MatrixMAGroupFEMRhs[4] = { D100, D010, D001, D000 };
// int SpacesNumbers_MatrixMAGroupFEMRhs[4] = { 0, 0, 0, 0 };
// int N_Matrices_MatrixMAGroupFEMRhs = 2;
// int RowSpace_MatrixMAGroupFEMRhs[2] = { 0, 0};
// int ColumnSpace_MatrixMAGroupFEMRhs[2] = { 0, 0};
// int N_Rhs_MatrixMAGroupFEMRhs = 1;
// int RhsSpace_MatrixMAGroupFEMRhs = { 0 };

void MatrixMAGroupFEMAssemble(double Mult, double *coeff, double *param,
            double hK, 
            double **OrigValues, int *N_BaseFuncts,
            double ***LocMatrices, double **LocRhs)
{
    double **MatrixM, **MatrixA, val, *MatrixRowM, *MatrixRowA, *Rhs;;
    double ansatz000, ansatz100, ansatz010, ansatz001;
    double test000, test100, test010, test001;
    double *Orig0, *Orig1, *Orig2, *Orig3;
    int i,j, N_;
    double c0, c5; 
    
  MatrixM = LocMatrices[0];
  MatrixA = LocMatrices[1];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c5 = coeff[5]; // f
 
  for(i=0;i<N_;i++)
  {
    MatrixRowM = MatrixM[i];
    MatrixRowA = MatrixA[i];
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

      MatrixRowM[j] += Mult*ansatz000*test000;

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      MatrixRowA[j] += Mult * val;
    } // endfor j
  } // endfor i
}

void MatrixGroupFEMAssemble(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
    double **MatrixC1, **MatrixC2, **MatrixC3, **MatrixR, val;
    double *MatrixRowC1, *MatrixRowC2, *MatrixRowC3, *MatrixRowR, *Rhs;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, c5;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;

  MatrixC1 = LocMatrices[0];
  MatrixC2 = LocMatrices[1];
  MatrixC3 = LocMatrices[2];
  MatrixR = LocMatrices[3];
  Rhs = LocRhs[0];
  
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  
  c5 = coeff[5]; // f

  for(i=0;i<N_;i++)
  {
    MatrixRowC1 = MatrixC1[i];
    MatrixRowC2 = MatrixC2[i];
    MatrixRowC3 = MatrixC3[i];
    MatrixRowR = MatrixR[i];
    test000 = Orig3[i];
    Rhs[i] += Mult*test000*c5;
    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      MatrixRowC1[j] += Mult *ansatz100*test000;
      MatrixRowC2[j] += Mult *ansatz010*test000;
      MatrixRowC3[j] += Mult *ansatz001*test000;      
      MatrixRowR[j] += Mult *ansatz000*test000;
    } // endfor j
  } // endfor i
}


// ======================================================================
// MATRICES FOR REACTION PART OF BULK PRECIPITATION
// ======================================================================

// ======================================================================
// assemble mass matrix
// ======================================================================
void MatrixMAssemble_Bulk3D(double Mult, double *coeff, double *param,
                            double hK,
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz000;
  double test000;
  double *Orig0;
  int i,j, N_;

  Matrix = LocMatrices[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test000 = Orig0[i];
 
    for(j=0;j<N_;j++)
    {
      ansatz000 = Orig0[j];

      if (!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	  MatrixRow[j] += Mult * ansatz000*test000;
      else
	  MatrixRow[i] += Mult * ansatz000*test000;
    } // endfor j
  } // endfor i
}

// ======================================================================
// assemble matrix and rhs for upwind discretization
// ======================================================================
void MatricesA_Assemble_Bulk3D(double Mult, double *coeff, double *param,
			       double hK, 
			       double **OrigValues, int *N_BaseFuncts,
			       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRowA;
  double ansatz000, ansatz100, ansatz010, ansatz001 ;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c4, Pe;
  double k;

  MatrixA = LocMatrices[0];
 
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // u_3
  c4 = coeff[4]; // other concentration or 0

  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test100 = Orig0[i];
      test010 = Orig1[i];
      test001 = Orig2[i];
      test000 = Orig3[i];

      for(j=0;j<N_;j++)
      {
	  ansatz100 = Orig0[j];
	  ansatz010 = Orig1[j];
	  ansatz001 = Orig2[j];
	  ansatz000 = Orig3[j];

	  val =  c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
	  if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
	      ||(i==j))
	  {
	      val += c4*ansatz000*test000;
	      // stiffness matrix
	      MatrixRowA[j] += Mult * val;
	  }
	  else
	  {
	     // stiffness matrix
	     MatrixRowA[j] += Mult * val;
	     val += c4*ansatz000*test000;
	     // add to diagonal
	     MatrixRowA[i] += Mult * val;
	  }
      } // endfor j
  } // endfor i
}

void Rhs_Assemble_Bulk3D(double Mult, double *coeff, double *param,
			 double hK,
			 double **OrigValues, int *N_BaseFuncts,
			 double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double test000;
  double *Orig3;
  int i, N_;
  double c5;
  Rhs = LocRhs[0];
 
  N_ = N_BaseFuncts[0];
  Orig3 = OrigValues[0];
  c5 = coeff[5];  

  for(i=0;i<N_;i++)
  {
      test000 = Orig3[i];
      Rhs[i] += Mult*test000*c5;
  }
}

// ======================================================================
// assemble matrices and rhs for SUPG
// matrix A contains the part of the equation which is multiplied with
//        theta_1 Delta t
// matrix K contains the part of the SUPG term which comes from the 
//        temporal discretization of the concentration
// ======================================================================
void MatricesA_Assemble_SUPG_Bulk3D(double Mult, double *coeff, double *param,
				    double hK,
				    double **OrigValues, int *N_BaseFuncts,
				    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA,  **MatrixK, val, *MatrixRowA, *MatrixRowK;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_, sold_type=TDatabase::ParamDB->SOLD_TYPE;
  double c0, c1, c2, c3, c4, c5, Pe;
  double tau, bgradv, bb, k, c00, c11, c22, c33, c44;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double sigma, conc_old;

  MatrixA = LocMatrices[0];
  MatrixK= LocMatrices[1];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // u_3
  c4 = coeff[4]; // other concentration or 0
  c5 = coeff[5];

  c00 = c0;
  c11 = c1;
  c22 = c2;
  c33 = c3;
  // reactive coefficient, inclusive term from the temporal derivative
  c44 = 1.0/(time_step * theta1) + c4;
  if ( fabs(c11) > fabs(c22) )
  {
    if ( fabs(c11) > fabs(c33) )
    {
      bb = fabs(c11);
    }
    else
    {
      bb = fabs(c33);
    }
  }
  else
  {
    if ( fabs(c22) > fabs(c33) )
    {
      bb = fabs(c22);
    }
    else
    {
      bb = fabs(c33);
    }
  }

  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, c44, bb);
/*
  if (!sold_type)
      sigma = 0;
  else
      sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, bb, tau, param);
*/
 
  if (sold_type < 2)
  {
      for(i=0;i<N_;i++)
      {
	  MatrixRowA = MatrixA[i];
	  MatrixRowK = MatrixK[i];
	  test100 = Orig0[i];
	  test010 = Orig1[i];
	  test001 = Orig2[i];
	  test000 = Orig3[i];

	  // supg test term
	  bgradv = tau*(c1*test100+c2*test010+c3*test001);
	  
	  for(j=0;j<N_;j++)
	  {
	      ansatz100 = Orig0[j];
	      ansatz010 = Orig1[j];
	      ansatz001 = Orig2[j];
	      ansatz000 = Orig3[j];

	      val =  c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
	      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001) * (test000+bgradv);
	      if ((!TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
		  ||(i==j))
	      {
		  val += c4*ansatz000*(test000+bgradv);
		  // stiffness matrix
		  MatrixRowA[j] += Mult * val;
	      }
	      else
	      {
		  // stiffness matrix
		  MatrixRowA[j] += Mult * val;
		  val = c4*ansatz000*(test000+bgradv);
		  // add to diagonal
		  MatrixRowA[i] += Mult * val;
	      }
	      // sdfem part which comes from time derivative 
	      MatrixRowK[j] += Mult * ansatz000*bgradv;
	  } // endfor j
      } // endfor i
  }
/*  else
  {
      // NOT YET CORRECTED !!!
      OutPut("NOT YET CORRECTED"<<endl);
      exit(4711);
      for(i=0;i<N_;i++)
      {
	  MatrixRowA = MatrixA[i];
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
	      
	      val += tau * (c1*ansatz100+c2*ansatz010+c3*ansatz001
			    +c4*ansatz000) * bgradv;
	      if (bb > 0)
		  val += sigma * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01)       /(bb*bb);
	      
	      MatrixRowA[j] += Mult * val;


	  } // endfor j
      } // endfor i
      }*/
}
void Rhs_Assemble_SUPG_Bulk3D(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, val, *MatrixRowA, *MatrixRowK;
  double *Rhs;
  double ansatz000, ansatz100, ansatz010, ansatz001 ;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_, sold_type=TDatabase::ParamDB->SOLD_TYPE;
  double c0, c1, c2, c3, c4, c5, Pe;
  double tau, bgradv, bb, k, c00, c11, c22, c33, c44;
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double sigma, conc_old;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // u_3
  c4 = coeff[4]; // other concentration
  c5 = coeff[5];

  c00 = c0;
  c11 = c1;
  c22 = c2;
  c33 = c3;
  // reactive coefficient, inclusive term from the temporal derivative
  c44 = 1.0/(time_step * theta1) + c4;
  if ( fabs(c11) > fabs(c22) )
  {
    if ( fabs(c11) > fabs(c33) )
    {
      bb = fabs(c11);
    }
    else
    {
      bb = fabs(c33);
    }
  }
  else
  {
    if ( fabs(c22) > fabs(c33) )
    {
      bb = fabs(c22);
    }
    else
    {
      bb = fabs(c33);
    }
  }

  tau = Compute_SDFEM_delta(hK, c00, c11, c22, c33, c44, bb);

/*  if (!sold_type)
      sigma = 0;
  else
      sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, bb, tau, param);
*/ 
  if (sold_type < 2)
  {
      for(i=0;i<N_;i++)
      {
	  test100 = Orig0[i];
	  test010 = Orig1[i];
	  test001 = Orig2[i];
          test000 = Orig3[i];
	  
	  // supg test term
	  bgradv = tau*(c1*test100+c2*test010+c3*test001);
	  Rhs[i] +=  Mult*(test000 + bgradv)*c5;
	  
      }
  }
  else
      exit(4711);
}

// ======================================================================
// assemble matrices and rhs for Galerkin (FCT-FEM)
// matrix A contains the part of the equation which is multiplied with
//        theta_1 Delta t
// ======================================================================

void MatricesA_Assemble_Galerkin_Bulk3D(double Mult, double *coeff, double *param,
				  double hK, 
				  double **OrigValues, int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, val, *MatrixRowA;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, Pe;

  MatrixA = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // u_1
  c2 = coeff[2]; // u_2
  c3 = coeff[3]; // u_3
  c4 = coeff[4]; // other concentration or 0
  c5 = coeff[5];

  for(i=0;i<N_;i++)
  {
      MatrixRowA = MatrixA[i];
      test100 = Orig0[i];
      test010 = Orig1[i];
      test001 = Orig2[i];
      test000 = Orig3[i];
      
      for(j=0;j<N_;j++)
      {
	  ansatz100 = Orig0[j];
	  ansatz010 = Orig1[j];
	  ansatz001 = Orig2[j];
	  ansatz000 = Orig3[j];
	  
	  val =  c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
	  val += (c1*ansatz100+c2*ansatz010+c3*ansatz001) * test000;
	  val += c4*ansatz000*test000;
	  // stiffness matrix
	  MatrixRowA[j] += Mult * val;
      } // endfor j
  } // endfor i
}

// this routine gives the fe value back
void TimeCDParamsVeloField(double *in, double *out)
{
  // in[0], in[1] and in[2] are coordinates 
  out[0] = in[3]; // value of first input fe function
  out[1] = in[4]; // value of second input fe function
  out[2] = in[5]; // value of third input fe function
  // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

void TimeCDParamsBulk(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates
  if (in[3]>0)
     out[0] = in[3]; // value of first input fe function (other concentration)
  else
     out[0] = 0.0;
  out[1] = in[4]; // value of second input fe function (u1)
  out[2] = in[5]; // value of third input fe function (u2)
  out[3] = in[6]; // value of third input fe function (u3)
 // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

void TimeCDParamsBulk_Cc(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[5]; // u1
  out[1] = in[6]; // u2
  out[2] = in[7]; // u3
  out[3] = in[3]; // C_a
  out[4] = in[4]; // C_b
  out[5] = in[8]; // C_c_old
  out[6] = in[9]; // integral_value
}
// this routine gives the fe value back
void TimeCDParamsSOLD(double *in, double *out)
{
  // in[0], in[1] are coordinates 
  // cout << in[0] << " " << in[1] << endl;
  out[0] = in[3]; // current solution u
  out[1] = in[4]; // current solution u_x
  out[2] = in[5]; // current solution u_y
  out[3] = in[6]; // current solution u_z
  out[4] = in[7]; // old solution u
  out[5] = in[8]; // old solution u_x
  out[6] = in[9]; // old solution u_y
  out[7] = in[10]; // old solution u_z
}

// this routine gives the fe value back
void TimeCDParamsSolution(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[3]; // current solution u
}

void TimeCDParamsUrea(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates
  out[0] = in[3]; // value of second input fe function (u1)
  out[1] = in[4]; // value of third input fe function (u2)
  out[2] = in[5]; // value of third input fe function (u3)
 // cout << out[0] << " " << out[1] <<  " " << out[2] << endl;
}

void TimeCDParamsUrea_conc(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[3]; // u1
  out[1] = in[4]; // u2
  out[2] = in[5]; // u3
  out[3] = in[6]; // concentration
  out[4] = in[7]; // temp
  out[5] = in[8]; // integral conce
}

void TimeCDParamsUrea_temp(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[3]; // u1
  out[1] = in[4]; // u2
  out[2] = in[5]; // u3
  out[3] = in[6]; // concentration
  out[4] = in[7]; // temp
  out[5] = in[8]; // integral conce 
}

void TimeCDParamsUrea_conc_mat(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[3]; // u1
  out[1] = in[4]; // u2
  out[2] = in[5]; // u3
}

void TimeCDParamsUrea_conc2(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[3]; // u1
  out[1] = in[4]; // u2
  out[2] = in[5]; // u3
  out[3] = in[6]; // concentration
  out[4] = in[7]; // temp
  out[5] = in[8]; // integral conce
  out[6] = in[9]; // integral conce
}
void TimeCDParamsUrea_temp2(double *in, double *out)
{
  // in[0], in[1], in[2] are coordinates 
  out[0] = in[3]; // u1
  out[1] = in[4]; // u2
  out[2] = in[5]; // u3
  out[3] = in[6]; // concentration
  out[4] = in[7]; // temp
  out[5] = in[8]; // integral conce
  out[6] = in[9]; // integral conce 
}
