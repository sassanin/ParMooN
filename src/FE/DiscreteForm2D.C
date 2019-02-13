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
   
// =======================================================================
// @(#)DiscreteForm2D.C        1.6 10/18/99
// 
// Class:       TDiscreteForm2D
// Purpose:     assemble a couple of matrices and right-hand side at once
//
// Author:      Gunar Matthies (02.10.98)
//
// History:     start of implementation 02.10.98 (Gunar Matthies)
//        :     added moving domains and 2phaseflows 20.09.09 (Sashikumaar Ganesan)
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Database.h>
#include <FEDatabase2D.h>
#include <DiscreteForm2D.h>
#include <string.h>
#include <stdlib.h>
#include <MovingNavierStokes.h>

#ifdef __2D__
#include <NSE2D_Param.h>
#include <NSE2D_EquOrd_FixPo.h>
#include <NSE2D_FixPo.h>
#include <NSE2D_Friction_FixPo.h>
#include <NSE2D_FixPoRot.h>
#include <NSE2D_FixPoSkew.h>
#include <NSE2D_Newton.h>
#include <NSE2D_AxialSymm3D_FixPo.h>
#include <TNSE2D_FixPo.h>
#include <TNSE2D_FixPoRot.h>
#include <TNSE2D_FixPo_SSMUM.h>
#include <TNSE2D_Routines.h>
#include <TCD2D.h>

#include <MainUtilities.h>
#include <ConvDiff.h>
#include <ConvDiff2D.h>
#include <TimeConvDiff2D.h>

#endif


/** constructor with vector initialization */
TDiscreteForm2D::TDiscreteForm2D(char *name, char *description,
        int n_terms, MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFct2D *assemble, CoeffFct2D *coeffs,
        ManipulateFct2D *manipulate)
{
  int i, j, max;
  MultiIndex2D alpha;

  Name = strdup(name);
  Description = strdup(description);

  N_Terms = n_terms;
  Derivatives = derivatives;
  FESpaceNumber = fespacenumber;

  N_Matrices = n_matrices;
  N_Rhs = n_rhs;
  RowSpace = rowspace;
  ColumnSpace = columnspace;
  RhsSpace = rhsspace;

  Coeffs = coeffs;

  Assemble = assemble;
  AssembleParam = NULL;

  Manipulate = manipulate;

  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  // find number of spaces
  max = -1;
  for(i=0;i<N_Terms;i++)
  {
    j = FESpaceNumber[i];
    if(j > max) max = j;
  }

  N_Spaces = max+1;

  Needs2ndDerivatives = new bool[N_Spaces];
  for(i=0;i<N_Spaces;i++)
    Needs2ndDerivatives[i] = FALSE;

  for(i=0;i<N_Terms;i++)
  {
    alpha = Derivatives[i];
    j = FESpaceNumber[i];
    if(alpha == D20 || alpha == D11 || alpha == D02)
      Needs2ndDerivatives[j] = TRUE;
  }

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
  #endif 
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    cout << "---------------------" << endl;
    cout << "number of spaces: " << N_Spaces << endl;
    for(i=0;i<N_Spaces;i++)
      cout << i << " " << Needs2ndDerivatives[i] << endl;
    cout << "---------------------" << endl;
  }

}

/** constructor with assembling using parameters */
TDiscreteForm2D::TDiscreteForm2D(char *name, char *description,
        int n_terms, MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFctParam2D *assembleparam, CoeffFct2D *coeffs,
        ManipulateFct2D *manipulate)
{
  int i, j, max;
  MultiIndex2D alpha;

  Name = strdup(name);
  Description = strdup(description);

  N_Terms = n_terms;
  Derivatives = derivatives;
  FESpaceNumber = fespacenumber;

  N_Matrices = n_matrices;
  N_Rhs = n_rhs;
  RowSpace = rowspace;
  ColumnSpace = columnspace;
  RhsSpace = rhsspace;

  Coeffs = coeffs;

  Assemble = NULL;
  AssembleParam = assembleparam;

  Manipulate = manipulate;

  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  // find number of spaces
  max = -1;
  for(i=0;i<N_Terms;i++)
  {
    j = FESpaceNumber[i];
    if(j > max) max = j;
  }

  N_Spaces = max+1;

  Needs2ndDerivatives = new bool[N_Spaces];
  for(i=0;i<N_Spaces;i++)
    Needs2ndDerivatives[i] = FALSE;

  for(i=0;i<N_Terms;i++)
  {
    alpha = Derivatives[i];
    j = FESpaceNumber[i];
    if(alpha == D20 || alpha == D11 || alpha == D02)
      Needs2ndDerivatives[j] = TRUE;
  }

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
  #endif 
  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
    cout << "---------------------" << endl;
    cout << "number of spaces: " << N_Spaces << endl;
    for(i=0;i<N_Spaces;i++)
      cout << i << " " << Needs2ndDerivatives[i] << endl;
    cout << "---------------------" << endl;
  }
}

TDiscreteForm2D::~TDiscreteForm2D()
{
  delete AllOrigValues;
  delete OrigValues;
  delete Needs2ndDerivatives;
  delete Name;
  delete Description;
}

void TDiscreteForm2D::GetLocalForms(int N_Points, double *weights, 
                                    double *AbsDetjk, double hK, 
                                    double *X, double *Y,
                                    int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                                    double **Parameters, double **AuxArray,
                                    TBaseCell *Cell, int n_matrices, int n_rhs,
                                    double ***LocMatrix, double **LocRhs,
                                    double factor)
{
  int i,j,k,l, N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  double Mult, *Coeff, *Param;

  // cout << "in TDiscreteForm2D::GetLocalForms" << endl;

  for(i=0;i<n_matrices;i++)
  {
    CurrentMatrix = LocMatrix[i];
    N_Rows = N_BaseFuncts[RowSpace[i]];
    N_Columns = N_BaseFuncts[ColumnSpace[i]];
    for(j=0;j<N_Rows;j++)
    {
      MatrixRow = CurrentMatrix[j];
      memset(MatrixRow, 0, SizeOfDouble*N_Columns);
    } // endfor j
  } // endfor i

  for(i=0;i<n_rhs;i++)
  {
    N_Rows = N_BaseFuncts[RhsSpace[i]];
    memset(LocRhs[i], 0, SizeOfDouble*N_Rows);
  }

// *****************************************************
// for 2Phase flow problems (Sashikumaar Ganesan)
  AuxArray[0][0] = Cell->GetRegionID();
// *****************************************************

  if(Coeffs)
    Coeffs(N_Points, X, Y, Parameters, AuxArray);

  if(Manipulate)
    Manipulate(N_Points, AuxArray, Parameters, Cell);

  for(i=0;i<N_Terms;i++)
  {
    AllOrigValues[i] = 
      TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]], 
                                        Derivatives[i]);
  }

  for(i=0;i<N_Points;i++)
  {
    Mult = weights[i] * AbsDetjk[i] * factor;
    Coeff = AuxArray[i];
    Coeff[19] = AbsDetjk[i];
    
    if(TDatabase::ParamDB->Axial3DAxis==1)
    {
     Coeff[20] = Y[i];  // r in axial3D (X: symmetric) problems (Sashikumaar Ganesan)      
    }
    else
    {
    Coeff[20] = X[i];  // r in axial3D (Y: symmetric) problems (Sashikumaar Ganesan)      
    }

    Param = Parameters[i];

    for(j=0;j<N_Terms;j++)
      OrigValues[j] = AllOrigValues[j][i];

    if(Assemble)
      Assemble(Mult, Coeff, hK, OrigValues, N_BaseFuncts, 
               LocMatrix, LocRhs);

    if(AssembleParam)
      AssembleParam(Mult, Coeff, Param, hK, OrigValues, N_BaseFuncts, 
                    LocMatrix, LocRhs);
  } // endfor i
}

void TDiscreteForm2D::GetLocalForms(int N_Points, double *weights, 
                        double *AbsDetjk, double hK, 
                        double *X, double *Y,
                        int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                        TBaseCell *Cell,
                        double ***LocMatrix, double **LocRhs)
{
  double Mult;
  double *Coefficients[N_Points];
  double *aux = new double [N_Points*20]; // do not change below 20
  for(int j=0;j<N_Points;j++)
    Coefficients[j] = aux + j*20;
  
  if(Coeffs)
    Coeffs(N_Points, X, Y, NULL, Coefficients);

  if(Manipulate)
    Manipulate(N_Points, Coefficients, NULL, Cell);
  for(int j=0;j<N_Terms;j++)
  {
    AllOrigValues[j] = 
      TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[j]], 
                                        Derivatives[j]);
  }
  
  for(int i=0;i<N_Points;i++)
  {
 
    Mult = weights[i]*AbsDetjk[i];
    Coefficients[i][19] = AbsDetjk[i];
    
    for(int j=0;j<N_Terms;j++) {
       OrigValues[j] = AllOrigValues[j][i];
    }

 
    if(Assemble)
      Assemble(Mult, Coefficients[i], hK, OrigValues, N_BaseFuncts, 
               LocMatrix, LocRhs);
    if(AssembleParam)
      AssembleParam(Mult, Coefficients[i], NULL, hK, OrigValues, N_BaseFuncts, 
                    LocMatrix, LocRhs);
  } // endfor i
  delete [] aux;
}

/******************************************************************************/
//
// Routines for initializing discrete forms
//
/******************************************************************************/

// since this routine will also compiled for 3d:
#ifdef __2D__



// part for standard Galerkin
int N_Terms = 3;
MultiIndex2D Derivatives[3] = { D10, D01, D00 };
int SpacesNumbers[3] = { 0, 0, 0 };

// part for SDFEM
int N_Terms_SD = 5;
MultiIndex2D Derivatives_SD[5] = { D10, D01, D00, D20, D02 };
int SpacesNumbers_SD[5] = { 0, 0, 0, 0, 0 };

// part for UPWIND with lumping of reaction term and rhs
int N_Terms_UPW1 = 2;
MultiIndex2D Derivatives_UPW1[2] = { D10, D01 };
int SpacesNumbers_UPW1[2] = { 0, 0 };

// part for UPWIND without lumping of reaction term and rhs
int N_Terms_UPW2 = 3;
MultiIndex2D Derivatives_UPW2[3] = { D10, D01, D00 };
int SpacesNumbers_UPW2[3] = { 0, 0, 0 };

// part for rhs
int N_Terms_rhs = 1;
MultiIndex2D Derivatives_rhs[1] = { D00 };
int SpacesNumbers_rhs[1] = { 0 };

// part for all
int CD_N_Matrices = 1;
int CD_RowSpace[1] = { 0 };
int CD_ColumnSpace[1] = { 0 };
int CD_N_Rhs = 1;
int CD_RhsSpace[1] = { 0 };

MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
MultiIndex2D ZeroDerivative[1] = { D00};

// part for all
int N_Matrices = 1;
int RowSpace[1] = { 0 };
int ColumnSpace[1] = { 0 };
int N_Rhs = 1;
int RhsSpace[1] = { 0 };


// parameters for DC/CD shock capturing scheme
void DC_CD_Params(double *in, double *out);



// parameters: velocity field
void Params_Velo(double *in, double *out);


int DC_CD_N_FESpaces = 1;
int DC_CD_N_Fct = 2;
int DC_CD_N_ParamFct = 1;
int DC_CD_N_FEValues = 2;
int DC_CD_N_Params = 2;
int DC_CD_FEFctIndex[2] = { 0, 1 };
MultiIndex2D DC_CD_FEMultiIndex[2] = { D00, D00 };
ParamFct *DC_CD_Fct[1] = { DC_CD_Params };
int DC_CD_BeginParam[1] = { 0 };

// parameters for SC_2 shock capturing scheme
void SC_2_Params(double *in, double *out);

int SC_2_N_FESpaces = 2;
int SC_2_N_Fct = 3;
int SC_2_N_ParamFct = 1;
int SC_2_N_FEValues = 4;
int SC_2_N_Params = 4;
int SC_2_FEFctIndex[4] = { 0, 1, 2, 2 };
MultiIndex2D SC_2_FEMultiIndex[4] = { D00, D00, D10, D01 };
ParamFct *SC_2_Fct[1] = { SC_2_Params };
int SC_2_BeginParam[1] = { 0 };

// parameters for SOLD schemes
void SOLD_Params(double *in, double *out);

// parameters:  SOLD + velocity field
void SOLD_Params_And_Velo(double *in, double *out);


int SOLD_N_FESpaces = 2;
int SOLD_N_Fct = 3;
int SOLD_N_ParamFct = 1;
int SOLD_N_FEValues = 7;
int SOLD_N_Params = 7;
int SOLD_FEFctIndex[7] = { 0, 0, 0, 0, 0, 1, 2 };
MultiIndex2D SOLD_FEMultiIndex[7] = { D00, D10, D01, D20, D02, D00, D00 };
ParamFct *SOLD_Fct[1] = { SOLD_Params };
int SOLD_BeginParam[1] = { 0 };


// ======================================================================
// definitions for assembling the matrix A and rhs
// ======================================================================
int N_Terms_MatrixARhs = 3;
MultiIndex2D Derivatives_MatrixARhs[3] = { D10, D01, D00};
int SpacesNumbers_MatrixARhs[3] = { 0, 0, 0  };
int N_Matrices_MatrixARhs = 1;
int RowSpace_MatrixARhs[1] = { 0 };
int ColumnSpace_MatrixARhs[1] = { 0 };
int N_Rhs_MatrixARhs = 1;
int RhsSpace_MatrixARhs[1] = { 0 };



void MatrixARhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, val, *MatrixRowA, *MatrixRowK;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, Pe, h, x, r; 

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
  
  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
//   cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);  */   
  
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      MatrixRowA[j] += r*Mult * val;
                
    } // endfor j
  } // endfor i
}


// ======================================================================
// definitions for assembling the matrices A, K and rhs
// ======================================================================
int N_Terms_MatricesAKRhs_SUPG = 5;
MultiIndex2D Derivatives_MatricesAKRhs_SUPG[5] = { D10, D01, D00, D20, D02 };
int SpacesNumbers_MatricesAKRhs_SUPG[5] = { 0, 0, 0, 0, 0  };
int N_Matrices_MatricesAKRhs_SUPG = 2;
int RowSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
int ColumnSpace_MatricesAKRhs_SUPG[2] = { 0, 0 };
int N_Rhs_MatricesAKRhs_SUPG = 1;
int RhsSpace_MatricesAKRhs_SUPG[1] = { 0 };

int N_Matrices_MatricesAKRhs_SOLD = 3;
int RowSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };
int ColumnSpace_MatricesAKRhs_SOLD[3] = { 0, 0, 0 };

void MatricesAKRhsAssemble_SUPG_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixK, *Rhs, *MatrixRowA, *MatrixRowK;
  double **MatrixS, *MatrixRowS, ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, r, x;
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

  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_SUPG_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
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

  if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
  {
      // rhs from previous time step 
      c5 = coeff[5];   
      // compute residual
      res = param[0] + theta1*time_step*(-c0*(param[3]+param[4]) + c1*param[1]
           +c2*param[2] + c3*param[0])
    -param[5] +  theta2*time_step*(-c0*(param[8]+param[9]) + c1*param[6]
           +c2*param[7] + c3*param[5])
    -theta3*time_step*c5 - theta4*time_step*c4;
      /*c00 =  time_step * theta1 * c0;
      c11 =  time_step * theta1 * c1;
      c22 =  time_step * theta1 * c2;
      c33 = 1.0 + time_step * theta1 *c3;*/
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
  if ((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)&&(TDatabase::ParamDB->SDFEM_TYPE!=11))
    tau *= val;

  // loop over the basis functions
  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    MatrixRowK = MatrixK[i];
    if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
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
    Rhs[i] += r*Mult*test00*c4;

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
      
      MatrixRowA[j] += r*Mult * val;
      // time derivative part of the SUPG stabilization          
      MatrixRowK[j] += r*Mult * ansatz00*bgradv;
      
      // isotropic SOLD method
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==1))
      {
    MatrixRowS[j] += val2 * (test10*ansatz10+test01*ansatz01);
      }
      if ((TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE) && (TDatabase::ParamDB->SOLD_TYPE==2))
      {

    MatrixRowS[j] += val2 * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01);
      }
     
    } // endfor j
  } // endfor i
}



int N_Terms_MatrixMRhs_SUPG = 3;
MultiIndex2D Derivatives_MatrixMRhs_SUPG[3] = { D10, D01, D00 };
int SpacesNumbers_MatrixMRhs_SUPG[3] = { 0, 0, 0 };
int N_Matrices_MatrixMRhs_SUPG = 1;
int RowSpace_MatrixMRhs_SUPG[1] = { 0 };
int ColumnSpace_MatrixMRhs_SUPG[1] = { 0 };
int N_Rhs_MatrixMRhs_SUPG = 1;
int RhsSpace_MatrixMRhs_SUPG[1] = { 0 };

void MatrixMRhsAssemble_SUPG_Axial3D(double Mult, double *coeff, double *param,
                            double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **Matrix, **MatrixS, *Rhs, val, *MatrixRow, *MatrixSRow;
  double ansatz00, x, r;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c00, c11, c22, c33; 
  double tau, bgradv, bb;
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  Matrix = LocMatrices[0];
  MatrixS= LocMatrices[1];
  
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
  
  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_SUPG_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
   cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);  
  
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
  if ((TDatabase::ParamDB->SDFEM_TYPE!=9)&&(TDatabase::ParamDB->SDFEM_TYPE!=10)&&(TDatabase::ParamDB->SDFEM_TYPE!=11))
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

    Rhs[i] += r*Mult*(test00+tau*bgradv)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig2[j];

      MatrixRow[j] += r*Mult * ansatz00*test00;
      MatrixSRow[j] += r*Mult * ansatz00*bgradv;    
      
    } // endfor j
  } // endfor i
}


int N_Terms_MatrixMRhs = 1;
MultiIndex2D Derivatives_MatrixMRhs[1] = { D00 };
int SpacesNumbers_MatrixMRhs[1] = { 0 };
int N_Matrices_MatrixMRhs = 1;
int RowSpace_MatrixMRhs[1] = { 0 };
int ColumnSpace_MatrixMRhs[1] = { 0 };
int N_Rhs_MatrixMRhs = 1;
int RhsSpace_MatrixMRhs[1] = { 0 };



void MatrixMRhsAssemble_Axial3D(double Mult, double *coeff, double *param,
                           double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, *MatrixRow;
  double ansatz00, x, r;
  double test00;
  double *Orig0;
  int i,j, N_;
  double c4; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  c4 = coeff[4]; // f
  x  = coeff[20]; // see DiscreteForm2D.C
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check MatrixMRhsAssemble_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
//    cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];

    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];

      MatrixRow[j] += r*Mult*ansatz00*test00;
    } // endfor j
  } // endfor i
}




//////////////////////////////////////////////////////////////////////////////
// ======================================================================
// ======================================================================
// ======================================================================
// ======================================================================




void InitializeDiscreteFormsScalar(TDiscreteForm2D *&DiscreteFormMatrixMRhs, 
                                   TDiscreteForm2D *&DiscreteFormMatrixARhs,
                                   TDiscreteForm2D *&DiscreteFormMatrixMRhs_SUPG,
                                   TDiscreteForm2D *&DiscreteFormMatrixARhs_SUPG,
                                   CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SDFEMString[] = "SDFEM";

  if(TDatabase::ParamDB->Axial3D)
   {
    // discrete form for assembling mass matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixMRhs = new TDiscreteForm2D(GalerkinString, allString, N_Terms_MatrixMRhs, Derivatives_MatrixMRhs,
      SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
      RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
      MatrixMRhsAssemble_Axial3D, LinCoeffs, NULL);

    // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixARhs = new TDiscreteForm2D
      (GalerkinString, allString, N_Terms_MatrixARhs, Derivatives_MatrixARhs,
      SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
      RowSpace_MatrixARhs, ColumnSpace_MatrixARhs, RhsSpace_MatrixARhs,
      MatrixARhsAssemble_Axial3D, LinCoeffs, NULL);

    // discrete form for assembling mass matrix and rhs (Galerkin FEM)
//     DiscreteFormMatrixMRhs_SUPG = NULL;
    DiscreteFormMatrixMRhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatrixMRhs_SUPG, Derivatives_MatrixMRhs_SUPG,
      SpacesNumbers_MatrixMRhs_SUPG, N_Matrices_MatrixMRhs_SUPG, N_Rhs_MatrixMRhs_SUPG,
      RowSpace_MatrixMRhs_SUPG, ColumnSpace_MatrixMRhs_SUPG, RhsSpace_MatrixMRhs_SUPG,
      MatrixMRhsAssemble_SUPG_Axial3D, LinCoeffs, NULL);
      
      
    // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixARhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
      SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SUPG, N_Rhs_MatricesAKRhs_SUPG,
      RowSpace_MatricesAKRhs_SUPG, ColumnSpace_MatricesAKRhs_SUPG, RhsSpace_MatricesAKRhs_SUPG,
      MatricesAKRhsAssemble_SUPG_Axial3D, LinCoeffs, NULL);  
 
//     OutPut("InitializeDiscreteFormsScalar !!! Axial3D SUPG not implemented yet !!!!!!!!" << endl);
//     exit(0);
   }
  else
   {
    // discrete form for assembling mass matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixMRhs = new TDiscreteForm2D(GalerkinString, allString, N_Terms_MatrixMRhs, Derivatives_MatrixMRhs,
      SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
      RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
      MatrixMRhsAssemble, LinCoeffs, NULL);

    // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixARhs = new TDiscreteForm2D
      (GalerkinString, allString, N_Terms_MatrixARhs, Derivatives_MatrixARhs,
      SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
      RowSpace_MatrixARhs, ColumnSpace_MatrixARhs, RhsSpace_MatrixARhs,
      MatrixARhsAssemble, LinCoeffs, NULL);

    // discrete form for assembling mass matrix and rhs ( )
    DiscreteFormMatrixMRhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatrixMRhs_SUPG, Derivatives_MatrixMRhs_SUPG,
      SpacesNumbers_MatrixMRhs_SUPG, N_Matrices_MatrixMRhs_SUPG, N_Rhs_MatrixMRhs_SUPG,
      RowSpace_MatrixMRhs_SUPG, ColumnSpace_MatrixMRhs_SUPG, RhsSpace_MatrixMRhs_SUPG,
      MatrixMRhsAssemble_SUPG, LinCoeffs, NULL);

    // discrete form for assembling stiffness matrix and rhs ( )
    DiscreteFormMatrixARhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
      SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SUPG, N_Rhs_MatricesAKRhs_SUPG,
      RowSpace_MatricesAKRhs_SUPG, ColumnSpace_MatricesAKRhs_SUPG, RhsSpace_MatricesAKRhs_SUPG,
      MatricesAKRhsAssemble_SUPG, LinCoeffs, NULL);  
   }
  
}


void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormSDFEM,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLSDFEM,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormPressSep,
  TDiscreteForm2D *&DiscreteFormPressSepAuxProb,
  TDiscreteForm2D *&DiscreteFormNSRFBRhs,
  TDiscreteForm2D *&DiscreteFormNSCSTRhs,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";


  DiscreteFormVMSProjection = NULL;
  DiscreteFormNLVMSProjection = NULL;

  DiscreteFormPressSep = new TDiscreteForm2D(GalerkinString, allString,
					     NSPressSepN_Terms, NSPressSepDerivatives, 
					     NSPressSepSpaceNumbers,
					     NSPressSepN_Matrices, NSPressSepN_Rhs, 
					     NSPressSepRowSpace, NSPressSepColumnSpace,
					     NSPressSepRhsSpace, NSPressSep, LinCoeffs, NULL);

  DiscreteFormPressSepAuxProb = new TDiscreteForm2D(GalerkinString, allString,
                                                    NSPressSepAuxProbN_Terms, NSPressSepAuxProbDerivatives, 
                                                    NSPressSepAuxProbSpaceNumbers,
                                                    NSPressSepAuxProbN_Matrices, NSPressSepAuxProbN_Rhs,
                                                    NSPressSepAuxProbRowSpace, NSPressSepAuxProbColumnSpace,
                                                    NSPressSepAuxProbRhsSpace, NSPressSepAuxProb, LinCoeffs, NULL);

  DiscreteFormNSRFBRhs = new TDiscreteForm2D(GalerkinString, allString,
					   NSRFBRhsN_Terms, NSRFBRhsDerivatives, 
					   NSRFBRhsSpaceNumbers,
					   NSRFBRhsN_Matrices, NSRFBRhsN_Rhs, 
					   NSRFBRhsRowSpace, NSRFBRhsColumnSpace,
					   NSRFBRhsRhsSpace, NSRFBRhs, LinCoeffs, NULL);
  if (TDatabase::ParamDB->TENSOR_TYPE == 1)
  {
   DiscreteFormNSCSTRhs = new TDiscreteForm2D(GalerkinString, allString,
					   NSCSTRhsN_Terms, NSCSTRhsDerivatives, 
					   NSCSTRhsSpaceNumbers,
					   NSCSTRhsN_Matrices, NSCSTRhsN_Rhs, 
					   NSCSTRhsRowSpace, NSCSTRhsColumnSpace,
					   NSCSTRhsRhsSpace, NSCSTRhs, LinCoeffs, NULL);
  }
  else if (TDatabase::ParamDB->TENSOR_TYPE == 2)
  {
      DiscreteFormNSCSTRhs = new TDiscreteForm2D(GalerkinString, allString,
					   NSCSTRhsN_Terms, NSCSTRhsDerivatives, 
					   NSCSTRhsSpaceNumbers,
					   NSCSTRhsN_Matrices, NSCSTRhsN_Rhs, 
					   NSCSTRhsRowSpace, NSCSTRhsColumnSpace,
					   NSCSTRhsRhsSpace, NSCSTRhs_DEVSS, LinCoeffs, NULL);
  }
  else 
  {
   DiscreteFormNSCSTRhs = NULL; 
  }
  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
  {
    switch(NSTYPE)
    {
      case 1:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Galerkin, LinCoeffs, NULL);
                //NSType1RhsSpace, NSType1GalerkinFJMT07, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SDFEM, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Smagorinsky, LinCoeffs, NULL);

         DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  NSType1VMSProjectionN_Terms, NSType1VMSProjectionDerivatives, 
                  NSType1VMSProjectionSpaceNumbers,
                  NSType1VMSProjectionN_Matrices, NSType1N_Rhs, 
                  NSType1VMSProjectionRowSpace, NSType1VMSProjectionColumnSpace,
                  NSType1RhsSpace, NSType1VMSProjection, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);
                //NSType1NLRhsSpace, NSType1_2NLGalerkinFJMT07, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLSDFEMN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLSDFEMRhsSpace, NSType1NLSDFEM, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);

         DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  NSType1_2NLVMSProjectionN_Terms, NSType1_2NLVMSProjectionDerivatives, 
                  NSType1_2NLVMSProjectionSpaceNumbers,
                  NSType1_2NLVMSProjectionN_Matrices, NSType1NLN_Rhs,
                  NSType1_2NLVMSProjectionRowSpace, NSType1_2NLVMSProjectionColumnSpace,
                  NSType1NLRhsSpace, NSType1_2NLVMSProjection, LinCoeffs, NULL);
        }
	else
	    // skew-symmetric case
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SDFEMSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLSDFEMN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLSDFEMRhsSpace, NSType1NLSDFEMSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 2:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Galerkin, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEM, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEM, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);
        }
        else 
        { // skew symmetric case
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEMSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEMSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Galerkin, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);
          }
	if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
	{ // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
	if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
	{ // rotation form
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyRot, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRot, LinCoeffs, NULL);
	}
	}
	else
        {
          // (D(u):D(v))
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
          { // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
          { // rotation form
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRotDD, LinCoeffs, NULL);
	  }
	}
        break;

      case 4:
	case 14:
	  
	 
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        { 
	    if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Galerkin, LinCoeffs, NULL);
	    if (TDatabase::ParamDB->NSTYPE!=14)
		DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
							NSType4SDN_Terms, NSType4SDDerivatives, 
							NSType4SDSpaceNumbers,
							NSType4SDN_Matrices, NSType4SDN_Rhs, 
							NSType4SDRowSpace, NSType4SDColumnSpace,
							NSType4SDRhsSpace, NSType4SDFEM, LinCoeffs, NULL);
	    else
		DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
							NSType4EquOrdN_Terms, NSType4EquOrdDerivatives, 
							NSType4EquOrdSpaceNumbers,
							NSType4EquOrdN_Matrices, NSType4EquOrdN_Rhs, 
							NSType4EquOrdRowSpace, NSType4EquOrdColumnSpace,
							NSType4EquOrdRhsSpace, NSType4SDFEMEquOrd, LinCoeffs, NULL);
	      
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
	    if (TDatabase::ParamDB->NSTYPE!=14)
		DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
							  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
							  NSType4NLSDSpaceNumbers,
							  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
							  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
							  NSType4NLSDRhsSpace, NSType4NLSDFEM, LinCoeffs, NULL);
	    else
		DiscreteFormNLSDFEM = DiscreteFormSDFEM;
/*
	    new TDiscreteForm2D(SDFEMString, nonlinearString,
							  NSType4NLEquOrdN_Terms, NSType4NLEquOrdDerivatives, 
							  NSType4NLEquOrdSpaceNumbers,
							  NSType4NLEquOrdN_Matrices, NSType4NLEquOrdN_Rhs, 
							  NSType4NLEquOrdRowSpace, NSType4NLEquOrdColumnSpace,
							  NSType4NLEquOrdRhsSpace, NSType4NLSDFEMEquOrd, LinCoeffs, NULL);
*/
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);

	    // not yet implemented 
	    if (TDatabase::ParamDB->NSTYPE==14)
	    {
		DiscreteFormGalerkin = NULL;
		DiscreteFormUpwind = NULL;
		DiscreteFormSmagorinsky = NULL;
		DiscreteFormVMSProjection = NULL;
		DiscreteFormNLGalerkin = NULL;
		DiscreteFormNLUpwind = NULL;
		DiscreteFormNLSmagorinsky = NULL;
		DiscreteFormNLVMSProjection = NULL;		
	    }

          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMRot, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyRot, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMRot, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRot, LinCoeffs, NULL);
          }
        }
        else
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDD, LinCoeffs, NULL);

            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDD, LinCoeffs, NULL);
	  if (TDatabase::ParamDB->TENSOR_TYPE ==1)
	  {
	  	  DiscreteFormVMSProjection = new TDiscreteForm2D(UpwindString, allString,
                  NSType4VMSProjection_CST_N_Terms, NSType4VMSProjection_CST_Derivatives, 
                  NSType4VMSProjection_CST_SpaceNumbers,
                  NSType4VMSProjection_CST_N_Matrices, NSType4N_Rhs, 
                  NSType4VMSProjection_CST_RowSpace, NSType4VMSProjection_CST_ColumnSpace,
                  NSType4RhsSpace, NSType4VMSProjectionDD, LinCoeffs, NULL);
	    
	  }
	  else
	  {
	  DiscreteFormVMSProjection = new TDiscreteForm2D(UpwindString, allString,
                  NSType4VMSProjectionN_Terms, NSType4VMSProjectionDerivatives, 
                  NSType4VMSProjectionSpaceNumbers,
                  NSType4VMSProjectionN_Matrices, NSType4N_Rhs, 
                  NSType4VMSProjectionRowSpace, NSType4VMSProjectionColumnSpace,
                  NSType4RhsSpace, NSType4VMSProjectionDD, LinCoeffs, NULL);
	  }
	    
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
          {
            
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkewDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkewDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
          {
            
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMRotDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMRotDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRotDD, LinCoeffs, NULL);
          }
        }
        break;
    } // endswitch
  } 
  else //Newton iteration
  {
    switch(NSTYPE)
    {
      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
        // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDDNewton, LinCoeffs, NULL);
  

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyNewton, LinCoeffs, NULL);
 
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinNewton, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDNewton, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDNewton, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  }
}
void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLGalerkinDuese,
  CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";
  char UpwindString[] = "Upwind";

  DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                             NSType4N_Terms, NSType4Derivatives, 
                                             NSType4SpaceNumbers,
                                             NSType4N_Matrices, NSType4N_Rhs, 
                                             NSType4RowSpace, NSType4ColumnSpace,
                                             NSType4RhsSpace, NSType4GalerkinAxialSymm3D, LinCoeffs, NULL);

  DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                                               NSType4NLN_Terms, NSType4NLDerivatives, 
                                               NSType4NLSpaceNumbers,
                                               NSType4NLN_Matrices, NSType4NLN_Rhs, 
                                               NSType4NLRowSpace, NSType4NLColumnSpace,
                                               NSType4NLRhsSpace, NSType3_4NLGalerkinAxialSymm3D, LinCoeffs, NULL);

  DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
					   NSType4N_Terms, NSType4Derivatives, 
					   NSType4SpaceNumbers,
					   NSType4N_Matrices, NSType4N_Rhs, 
					   NSType4RowSpace, NSType4ColumnSpace,
					   NSType4RhsSpace, NSType4UpwindAxialSymm3D, LinCoeffs, NULL);
  
  DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
					     NSType4NLN_Terms, NSType4NLDerivatives, 
					     NSType4NLSpaceNumbers,
					     NSType4NLN_Matrices, NSType4NLN_Rhs, 
					     NSType4NLRowSpace, NSType4NLColumnSpace,
					     NSType4NLRhsSpace, NSType3_4NLUpwindAxialSymm3D, LinCoeffs, NULL);

  DiscreteFormNLGalerkinDuese = new TDiscreteForm2D(GalerkinString, nonlinearString,
                                               NSType4N_Terms,NSType4Derivatives,
					       NSType4SpaceNumbers,
                                               NSType4NLSDNewtonN_Matrices, NSType4NLN_Rhs, 
					       NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                                               NSType4NLRhsSpace, NSType3_4NLGalerkinAxialSymm3D_Duese, LinCoeffs, NULL);
}

// friction
void InitializeDiscreteFormsFriction(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormSDFEM,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLSDFEM,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormGalerkinFriction, 
  TDiscreteForm2D *&DiscreteFormNLGalerkinFriction,
  TDiscreteForm2D *&DiscreteFormGalerkinFrictionLocal, 
  TDiscreteForm2D *&DiscreteFormNLGalerkinFrictionLocal,
  TDiscreteForm2D *&DiscreteFormSDFEMFriction,
  TDiscreteForm2D *&DiscreteFormNLSDFEMFriction,
  TDiscreteForm2D *&DiscreteFormSDFEMDDFriction,
  TDiscreteForm2D *&DiscreteFormNLSDFEMDDFriction,
  TDiscreteForm2D *&DiscreteFormSDFEMDDFrictionRhs,
  TDiscreteForm2D *&DiscreteFormNLSDFEMDDFrictionRhs,
  TDiscreteForm2D *&DiscreteFormSDFEMFrictionRST,
  TDiscreteForm2D *&DiscreteFormNLSDFEMFrictionRST,
  TDiscreteForm2D *&DiscreteFormUpwindFriction,
  TDiscreteForm2D *&DiscreteFormNLUpwindFriction,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";


  DiscreteFormUpwindFriction = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives,
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindFriction, LinCoeffs, NULL);

  DiscreteFormNLUpwindFriction = new TDiscreteForm2D(UpwindString, 
                  nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindFriction, LinCoeffs, NULL);
  
  DiscreteFormGalerkinFrictionLocal 
           = new TDiscreteForm2D(GalerkinString, 
                  allString,
		  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDFrictionLocal, LinCoeffs, NULL);
  
  DiscreteFormNLGalerkinFrictionLocal 
           = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDDFrictionLocal, 
                  LinCoeffs, NULL);


  DiscreteFormGalerkinFriction 
           = new TDiscreteForm2D(GalerkinString, 
                  allString,
		  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDFriction, LinCoeffs, NULL);
  
  DiscreteFormNLGalerkinFriction 
           = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDDFriction, 
                  LinCoeffs, NULL);

  DiscreteFormSDFEMFriction = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMFriction, LinCoeffs, NULL);

  DiscreteFormNLSDFEMFriction = new TDiscreteForm2D(SDFEMString, 
		  nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMFriction, LinCoeffs, 
                  NULL);
  DiscreteFormSDFEMDDFriction = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDFriction, LinCoeffs, NULL);

  DiscreteFormNLSDFEMDDFriction = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDFriction, LinCoeffs, NULL);
    
  DiscreteFormSDFEMDDFrictionRhs = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDFrictionRhs, LinCoeffs, NULL);

  DiscreteFormNLSDFEMDDFrictionRhs = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDFrictionRhs, LinCoeffs, NULL);


  DiscreteFormSDFEMFrictionRST = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMFrictionRST, LinCoeffs, NULL);

  DiscreteFormNLSDFEMFrictionRST = new TDiscreteForm2D(SDFEMString, 
		  nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMFrictionRST, LinCoeffs, 
                  NULL);
  
  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
  {
    switch(NSTYPE)
    {
      case 1:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Galerkin, LinCoeffs, NULL);

	  DiscreteFormUpwind= new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);

	  DiscreteFormNLUpwindFriction = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);
        }
	// skew--symmetric
        else
	{
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 2:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Galerkin, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEM, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEM, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);
        }
        else 
        { // skew symmetric case
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEMSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEMSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Galerkin, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);
          }
          else
          { // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
        }
        else
        {
          // (D(u):D(v))
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          else
          { // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Galerkin, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEM, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEM, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);
          }
          else
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
        }
        else
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          else
          {
            
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkewDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkewDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
        }
        break;
    } // endswitch
  } 
  else //Newton iteration
  {
    switch(NSTYPE)
    {
      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
        // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDDNewton, LinCoeffs, NULL);
  

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyNewton, LinCoeffs, NULL);
 
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinNewton, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDNewton, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDNewton, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  }
}

void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormColetti,
  TDiscreteForm2D *&DiscreteFormGL00Convolution, 
  TDiscreteForm2D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind, 
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLColetti, 
  TDiscreteForm2D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm2D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHSColetti,
  TDiscreteForm2D *&DiscreteFormRHSLESModel,
  TDiscreteForm2D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm2D *&DiscreteFormRHSSmagorinskyExpl,
  TDiscreteForm2D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm2D *&DiscreteFormRHSAuxProblemU,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";
  char rhsString[] = "rhs";
  char ColettiString[] = "Coletti";
  char GaldiLaytonString[] = "Galdi-Layton";
  char auxprobString[] = "aux prob";
  char Layton96String[] = "Layton96";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS, LinCoeffs, NULL);

  DiscreteFormRHSAuxProblemU = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHSAuxProblemU, LinCoeffs, NULL);

  DiscreteFormMatrixAuxProblemU = new TDiscreteForm2D(auxprobString, allString,
            MatrixAuxiliaryProblemN_Terms, MatrixAuxiliaryProblemDerivatives, 
            MatrixAuxiliaryProblemSpaceNumbers,
            MatrixAuxiliaryProblemN_Matrices, MatrixAuxiliaryProblemN_Rhs, 
            MatrixAuxiliaryProblemRowSpace, MatrixAuxiliaryProblemColumnSpace,
            MatrixAuxiliaryProblemRhsSpace, MatrixAuxiliaryProblem, LinCoeffs, NULL);

  DiscreteFormRHSColetti = new TDiscreteForm2D(ColettiString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSColetti, LinCoeffs, NULL);


   DiscreteFormRHSLESModel =  new TDiscreteForm2D(GaldiLaytonString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSLESModel, LinCoeffs, NULL);

  DiscreteFormGL00AuxProblemRHS = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSGL00AuxProblemRHSN_Terms, TimeNSGL00AuxProblemRHSDerivatives, 
            TimeNSGL00AuxProblemRHSSpaceNumbers,
            TimeNSGL00AuxProblemRHSN_Matrices, TimeNSGL00AuxProblemRHSN_Rhs, 
            TimeNSGL00AuxProblemRHSRowSpace, TimeNSGL00AuxProblemRHSColumnSpace,
            TimeNSGL00AuxProblemRHSRhsSpace, TimeNSGL00AuxProblemRHS,
            LinCoeffs, NULL);

  // this is temporarily changed to the discrete form for the defect correction 
	    if ( TDatabase::ParamDB->DEFECT_CORRECTION_TYPE==0)
		DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
								    TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
								    TimeNSRHSColSpaceNumbers,
								    TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
								    TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
								    TimeNSRHSColRhsSpace, TimeNSRHSSmagorinskyExplicit,
								    LinCoeffs, NULL);
	    if ( TDatabase::ParamDB->DEFECT_CORRECTION_TYPE==1)
		DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
								    TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
								    TimeNSRHSColSpaceNumbers,
								    TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
								    TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
								    TimeNSRHSColRhsSpace, TimeNSRHSDefectCorrectionU2,
								    LinCoeffs, NULL);
	    
	    if ( TDatabase::ParamDB->DEFECT_CORRECTION_TYPE==2)
		DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
								    TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
								    TimeNSRHSColSpaceNumbers,
								    TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
								    TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
								    TimeNSRHSColRhsSpace, TimeNSRHSDefectCorrectionU2_1,
								    LinCoeffs, NULL);
		

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
  
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);
     
       DiscreteFormSmagorinsky = new TDiscreteForm2D(Layton96String, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, 
                TimeNSType1SpaceNumbers,
                TimeNSType1GL00AuxProblemN_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1GL00AuxProblemRowSpace, TimeNSType1GL00AuxProblemColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GL00AuxProblem, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        // same as Smagorinsky
        DiscreteFormNLColetti =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky, LinCoeffs, NULL);


        // same as Smagorinsky
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2GL00AuxProblemN_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2GL00AuxProblemRowSpace, TimeNSType2GL00AuxProblemColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2GL00AuxProblem, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
       break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

	  if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
	      DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
		  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
 
          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
          
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

	  if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
	      DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyRot, LinCoeffs, NULL);

        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
 
          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblemDD, LinCoeffs, NULL);

         DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3VMSProjectionN_Terms, TimeNSType3VMSProjectionDerivatives, 
                  TimeNSType3VMSProjectionSpaceNumbers,
                  TimeNSType3VMSProjectionN_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3VMSProjectionRowSpace, TimeNSType3VMSProjectionColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
 
         DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        
       }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky, LinCoeffs, NULL);


          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
        else
        { 
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;
        
          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblemDD, LinCoeffs, NULL);

          DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4VMSProjectionN_Terms, TimeNSType4VMSProjectionDerivatives, 
                  TimeNSType4VMSProjectionSpaceNumbers,
                  TimeNSType4VMSProjectionN_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4VMSProjectionRowSpace, TimeNSType4VMSProjectionColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
   
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

          DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        }
        break;
    } // endswitch

  else
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}
void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormColetti,
  TDiscreteForm2D *&DiscreteFormGL00Convolution, 
  TDiscreteForm2D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind, 
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLColetti, 
  TDiscreteForm2D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm2D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHSColetti,
  TDiscreteForm2D *&DiscreteFormRHSLESModel,
  TDiscreteForm2D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm2D *&DiscreteFormRHSSmagorinskyExpl,
  TDiscreteForm2D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm2D *&DiscreteFormRHSAuxProblemU,
  TDiscreteForm2D *&DiscreteFormC,
  TDiscreteForm2D *&DiscreteFormJ,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";
  char rhsString[] = "rhs";
  char ColettiString[] = "Coletti";
  char GaldiLaytonString[] = "Galdi-Layton";
  char auxprobString[] = "aux prob";
  char Layton96String[] = "Layton96";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS, LinCoeffs, NULL);

  DiscreteFormRHSAuxProblemU = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHSAuxProblemU, LinCoeffs, NULL);

  DiscreteFormMatrixAuxProblemU = new TDiscreteForm2D(auxprobString, allString,
            MatrixAuxiliaryProblemN_Terms, MatrixAuxiliaryProblemDerivatives, 
            MatrixAuxiliaryProblemSpaceNumbers,
            MatrixAuxiliaryProblemN_Matrices, MatrixAuxiliaryProblemN_Rhs, 
            MatrixAuxiliaryProblemRowSpace, MatrixAuxiliaryProblemColumnSpace,
            MatrixAuxiliaryProblemRhsSpace, MatrixAuxiliaryProblem, LinCoeffs, NULL);

  DiscreteFormRHSColetti = new TDiscreteForm2D(ColettiString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSColetti, LinCoeffs, NULL);


   DiscreteFormRHSLESModel =  new TDiscreteForm2D(GaldiLaytonString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSLESModel, LinCoeffs, NULL);

  DiscreteFormGL00AuxProblemRHS = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSGL00AuxProblemRHSN_Terms, TimeNSGL00AuxProblemRHSDerivatives, 
            TimeNSGL00AuxProblemRHSSpaceNumbers,
            TimeNSGL00AuxProblemRHSN_Matrices, TimeNSGL00AuxProblemRHSN_Rhs, 
            TimeNSGL00AuxProblemRHSRowSpace, TimeNSGL00AuxProblemRHSColumnSpace,
            TimeNSGL00AuxProblemRHSRhsSpace, TimeNSGL00AuxProblemRHS,
            LinCoeffs, NULL);

  DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
            TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSSmagorinskyExplicit,
            LinCoeffs, NULL);

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        cout << "Super!!!" << endl;
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);
     
       DiscreteFormSmagorinsky = new TDiscreteForm2D(Layton96String, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, 
                TimeNSType1SpaceNumbers,
                TimeNSType1GL00AuxProblemN_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1GL00AuxProblemRowSpace, TimeNSType1GL00AuxProblemColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GL00AuxProblem, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        // same as Smagorinsky
        DiscreteFormNLColetti =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        
        // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                0, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                0, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GalerkinJ, LinCoeffs, NULL);
        
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky, LinCoeffs, NULL);


        // same as Smagorinsky
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2GL00AuxProblemN_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2GL00AuxProblemRowSpace, TimeNSType2GL00AuxProblemColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2GL00AuxProblem, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        
        // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                0, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                0, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType1GalerkinJ, LinCoeffs, NULL);
        
       break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
          
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
 
          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblemDD, LinCoeffs, NULL);

         DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3VMSProjectionN_Terms, TimeNSType3VMSProjectionDerivatives, 
                  TimeNSType3VMSProjectionSpaceNumbers,
                  TimeNSType3VMSProjectionN_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3VMSProjectionRowSpace, TimeNSType3VMSProjectionColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
 
         DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        
       }
       
         // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                0, TimeNSType3N_Rhs, 
                TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                TimeNSType3RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                0, TimeNSType3N_Rhs, 
                TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                TimeNSType3RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GalerkinJ, LinCoeffs, NULL);
       
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky, LinCoeffs, NULL);


          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
        else
        {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;
        
          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblemDD, LinCoeffs, NULL);

          DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4VMSProjectionN_Terms, TimeNSType4VMSProjectionDerivatives, 
                  TimeNSType4VMSProjectionSpaceNumbers,
                  TimeNSType4VMSProjectionN_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4VMSProjectionRowSpace, TimeNSType4VMSProjectionColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
   
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

          DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        }
        
        // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                0, TimeNSType4N_Rhs, 
                TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                TimeNSType4RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                0, TimeNSType4N_Rhs, 
                TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                TimeNSType4RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                TimeNSType4NLN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                TimeNSType4RhsSpace, TimeNSType3GalerkinJ, LinCoeffs, NULL); 
        
        break;
    } // endswitch
  else
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}
void InitializeDiscreteFormsPaper2(  
  TDiscreteForm2D *&DiscreteFormRHSGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHSPaper2,
  CoeffFct2D *LinCoeffs, CoeffFct2D *Coeffs)
{
  char GaldiLaytonString[] = "GaldiLayton";
  char auxprobString[] = "aux prob";

  DiscreteFormRHSGL00AuxProblem = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSGL00AuxProblemPaper2, LinCoeffs, NULL);

  DiscreteFormGL00AuxProblemRHSPaper2 = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSGL00AuxProblemRHSN_Terms, TimeNSGL00AuxProblemRHSDerivatives, 
            TimeNSGL00AuxProblemRHSSpaceNumbers,
            TimeNSGL00AuxProblemRHSN_Matrices, TimeNSGL00AuxProblemRHSN_Rhs, 
            TimeNSGL00AuxProblemRHSRowSpace, TimeNSGL00AuxProblemRHSColumnSpace,
            TimeNSGL00AuxProblemRHSRhsSpace,  TimeNSGL00AuxProblemRHSPaper2 ,
            Coeffs, NULL);
}
void InitializeDiscreteFormsVMS(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind, 
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteForm_ho_RHS,
  TDiscreteForm2D *&DiscreteForm_ls_RHS,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";
  char rhsString[] = "rhs";
  char ColettiString[] = "Coletti";
  char GaldiLaytonString[] = "Galdi-Layton";
  char auxprobString[] = "aux prob";
  char Layton96String[] = "Layton96";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS, LinCoeffs, NULL);

  DiscreteForm_ho_RHS = new TDiscreteForm2D(
    GalerkinString, rhsString,
    TimeNS_ho_RHSN_Terms, TimeNS_ho_RHSDerivatives, TimeNS_ho_RHSSpaceNumbers,
    TimeNS_ho_RHSN_Matrices, TimeNS_ho_RHSN_Rhs, 
    TimeNS_ho_RHSRowSpace, TimeNS_ho_RHSColumnSpace,
    TimeNS_ho_RHSRhsSpace, TimeNS_VMS_SmallRhs2D, LinCoeffs, NULL);

  // THIS SHOULD BE DELETED IN THE MAIN PROGRAM
  DiscreteForm_ls_RHS = NULL;

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);
     
       DiscreteFormSmagorinsky = new TDiscreteForm2D(Layton96String, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
        }
        else
        {
	    //switch(TDatabase::ParamDB->VMS_TYPE)
	    //{
	    //case 0:
                 // (D(u):D(v))
                 DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                                            TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                                                            TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                                                            TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                                                            TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
                 
                 DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                                                          TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                                                          TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                                                          TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                                                          TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                                                               TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                                                               TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                                                               TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                                                               TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                                                              TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                              TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                                                              TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                                                              TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                                                            TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                            TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                                                            TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                                                            TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                                                                 TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                                 TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                                                                 TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                                                                 TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
		 /*     break;
             case 1:
                 // (D(u):D(v))
                 DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                         TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                         TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                         TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                         TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
                 
                 DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                         TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                         TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                         TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                         TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                         TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                         TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                         TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                         TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                         TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                         TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                         TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                         TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin_VMS_1_DD, LinCoeffs, NULL);
                 
                 DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                                                            TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                            TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                                                            TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                                                            TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                                                                 TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                                 TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                                                                 TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                                                                 TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
                 break;
		 }*/
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  else
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}

void InitializeDiscreteForms_SSMUM(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHS1,
  CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";
  char rhsString[] = "rhs";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS_SSMUM_ALE, LinCoeffs, NULL);

  DiscreteFormRHS1 = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNS_REL_VELO_SSMUM_WITH_ROTFRAME, LinCoeffs, NULL);

  DiscreteFormGalerkin = 
      new TDiscreteForm2D(GalerkinString, allString,
			  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
			  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
			  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
			  TimeNSType4RhsSpace, TimeNSType4Galerkin_SSMUM_ALE, LinCoeffs, NULL);
  
  DiscreteFormNLGalerkin = 
      new TDiscreteForm2D(GalerkinString, nonlinearString,
			  TimeNSType4_SSMUM_NLN_Terms, TimeNSType4_SSMUM_NLDerivatives, TimeNSType4_SSMUM_NLSpaceNumbers,
			  TimeNSType4_SSMUM_NLN_Matrices, TimeNSType4_SSMUM_NLN_Rhs, 
			  TimeNSType4_SSMUM_NLRowSpace, TimeNSType4_SSMUM_NLColumnSpace,
			  TimeNSType4_SSMUM_NLRhsSpace, TimeNSType3_4NLGalerkin_SSMUM_ALE, LinCoeffs, NULL);
}

/******************************************************************************/
// InitializeDiscreteFormsCDAdapt2D()
// initializes parameters of the data base for the main program CDAdap_2D.C
/******************************************************************************/
void InitializeDiscreteFormsCDAdapt2D(TDiscreteForm2D **DiscreteForms,
CoeffFct2D *BilinearCoeffs)
{
  char CdString[] = "Conv-Diff";
  char GalString[] = "Galerkin";
  char SDFEMString[] = "SDFEM";
  char UpwString[] = "Upwind";
  char SoldString[] = "SOLD";
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormGalerkinMatrix;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormSOLD;
  TDiscreteForm2D *DiscreteFormSOLD_Orthogonal;
  TDiscreteForm2D *DiscreteFormRhsLP96;
  TDiscreteForm2D *DiscreteFormMH_Kno06;
  TDiscreteForm2D *DiscreteFormSD_SOLD;
  TDiscreteForm2D *DiscreteFormRhsAdjointEnergyErrorEstimate;
  TDiscreteForm2D *DiscreteFormRhsAdjointTV;
  TDiscreteForm2D *DiscreteFormRhsAdjointTV2;
  TDiscreteForm2D *DiscreteFormRhsAdjointNormBL1_NormBorthL1;
  TDiscreteForm2D *DiscreteFormRhsAdjointNormResidualL1_NormBorthL1;
  TDiscreteForm2D *DiscreteFormRhsAdjointAll;
  TDiscreteForm2D *DiscreteFormRhsAdjointL2Error;
  TDiscreteForm2D *DiscreteFormRhsAdjointH1Error;
  TDiscreteForm2D *DiscreteForm2LevelLPS_Q0;
  int i;

  ManipulateFct2D *manipulate;
  if (TDatabase::ParamDB->SDFEM_NORM_B==0)
    manipulate = linfb;
  else
    manipulate = ave_l2b_quad_points;

  DiscreteFormGalerkin = new TDiscreteForm2D
    (CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    (AssembleFctParam2D*)BilinearAssembleGalerkin, BilinearCoeffs, NULL);

  DiscreteFormSDFEM = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble_SD, BilinearCoeffs,  manipulate);

  DiscreteFormUpwind = new TDiscreteForm2D
    (CdString, UpwString , N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    (AssembleFctParam2D*)BilinearAssemble_UPW2, BilinearCoeffs, NULL);

  DiscreteFormSOLD = new TDiscreteForm2D
    (CdString, SoldString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_SOLD, BilinearCoeffs,  manipulate);

  DiscreteFormSOLD_Orthogonal = new TDiscreteForm2D
    (CdString, SoldString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_SOLD_Orthogonal, BilinearCoeffs,  manipulate);

  DiscreteFormRhsLP96 = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms, Derivatives, SpacesNumbers,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_LP96, BilinearCoeffs,  manipulate);

  DiscreteFormMH_Kno06 = new TDiscreteForm2D
    (CdString, SoldString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_MH_Kno06, BilinearCoeffs,  manipulate);

  DiscreteFormSD_SOLD = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble_SD_SOLD, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointEnergyErrorEstimate = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointEnergyEstimate, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointTV = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointTV, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointTV2 = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointTV2, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointNormBL1_NormBorthL1 = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointNormBL1_NormBorthL1, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointNormResidualL1_NormBorthL1 = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointNormResidualL1_NormBorthL1, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointAll = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointAll, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointL2Error = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_rhs,
    Derivatives_rhs, SpacesNumbers_rhs,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointL2Error, BilinearCoeffs,  manipulate);

  DiscreteFormRhsAdjointH1Error = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointH1Error, BilinearCoeffs,  manipulate);

  DiscreteFormGalerkinMatrix = new TDiscreteForm2D
    (CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, 0, RowSpace, ColumnSpace, NULL,
    (AssembleFctParam2D*)BilinearAssembleGalerkin, BilinearCoeffs, NULL);

  DiscreteForm2LevelLPS_Q0 = new TDiscreteForm2D
    (CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble2LevelLPS_Q0, BilinearCoeffs, NULL);

  // initializing
  for (i=0;i<=17;i++)
    DiscreteForms[i] = NULL;

  DiscreteForms[0] = DiscreteFormGalerkin;
  DiscreteForms[1] = DiscreteFormSDFEM;
  DiscreteForms[2] = DiscreteFormUpwind;
  DiscreteForms[3] = DiscreteFormSOLD;
  DiscreteForms[4] = DiscreteFormSOLD_Orthogonal;
  DiscreteForms[5] = DiscreteFormRhsLP96;
  DiscreteForms[6] = DiscreteFormMH_Kno06;
  DiscreteForms[7] = DiscreteFormSD_SOLD;
  DiscreteForms[8] = DiscreteFormRhsAdjointEnergyErrorEstimate;
  DiscreteForms[9] = DiscreteFormRhsAdjointTV;
  DiscreteForms[10] = DiscreteFormRhsAdjointTV2;
  DiscreteForms[11] = DiscreteFormRhsAdjointNormBL1_NormBorthL1;
  DiscreteForms[12] = DiscreteFormRhsAdjointNormResidualL1_NormBorthL1;
  DiscreteForms[13] = DiscreteFormRhsAdjointAll;
  DiscreteForms[14] = DiscreteFormRhsAdjointL2Error;
  DiscreteForms[15] = DiscreteFormRhsAdjointH1Error;
  DiscreteForms[16] = DiscreteFormGalerkinMatrix;
  DiscreteForms[17] = DiscreteForm2LevelLPS_Q0;
}

void  InitializeDiscreteForms_Moving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormNLGalerkin,
                                     TDiscreteForm2D *&DiscreteFormGrid, CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs)
{

  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";
  char GridString[] = "Grid";
  
  DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
                         GridN_Terms, GridDerivatives, GridSpaceNumbers,
                         GridN_Matrices, GridN_Rhs,
                         GridRowSpace, GridColumnSpace,
                         GridRhsSpace, GridAssemble4,
                         GridCoeffs, NULL);


 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                             TimeNSType4N_Matrices, TimeNSType4N_Rhs,
                             TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                             TimeNSType4RhsSpace, TimeNSType4GalerkinDD_Axial3D, LinCoeffs, NULL);

   DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                               TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                               TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs,
                               TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                               TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD_Axial3D, LinCoeffs, NULL);
  } 
 else
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                             TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                             TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                             TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);

   DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                               TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                               TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                               TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                               TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);

 }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int N_Terms_MatrixMARhs = 3;
MultiIndex2D Derivatives_MatrixMARhs[3] = { D10, D01, D00};
int SpacesNumbers_MatrixMARhs[3] = { 0, 0, 0  };
int N_Matrices_MatrixMARhs = 2;
int RowSpace_MatrixMARhs[2] = { 0, 0 };
int ColumnSpace_MatrixMARhs[2] = { 0, 0};
int N_Rhs_MatrixMARhs = 1;
int RhsSpace_MatrixMARhs[1] = { 0 };


int N_Terms_MatrixMARhs_SUPG = 5;
MultiIndex2D Derivatives_MatrixMARhs_SUPG[5] = { D10, D01, D00, D20, D02 };
int SpacesNumbers_MatrixMARhs_SUPG[5] = { 0, 0, 0, 0, 0  };
int N_Matrices_MatrixMARhs_SUPG = 3;
int RowSpace_MatrixMARhs_SUPG[3] = { 0, 0, 0 };
int ColumnSpace_MatrixMARhs_SUPG[3] = { 0, 0, 0};
int N_Rhs_MatrixMARhs_SUPG = 1;
int RhsSpace_MatrixMARhs_SUPG[1] = { 0 };



 


void InitializeDiscreteForms_ScalarMoving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormGrid,
                                    TDiscreteForm2D *&DiscreteFormMatrixMRhs_SUPG, TDiscreteForm2D *&DiscreteFormSUPG, CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SUPGString[] = "SUPG";
  char GridString[] = "Grid";
  
    DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
                         GridN_Terms, GridDerivatives, GridSpaceNumbers,
                         GridN_Matrices, GridN_Rhs,
                         GridRowSpace, GridColumnSpace,
                         GridRhsSpace, GridAssemble4,
                         GridCoeffs, NULL);
  
  
  
 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble_Axial3D, LinCoeffs, NULL);  
   
   
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;  
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented InitializeDiscreteForms !!!!!!"<<endl;     
    DiscreteFormSUPG = NULL; 
    DiscreteFormMatrixMRhs_SUPG  = NULL;
  } 
 else
 {   
  DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble, LinCoeffs, NULL);
  
    // discrete form for assembling mass matrix and rhs ( )
    DiscreteFormMatrixMRhs_SUPG = new TDiscreteForm2D(SUPGString, allString, N_Terms_MatrixMRhs_SUPG, Derivatives_MatrixMRhs_SUPG,
                                         SpacesNumbers_MatrixMRhs_SUPG, N_Matrices_MatrixMRhs_SUPG, N_Rhs_MatrixMRhs_SUPG,
                                         RowSpace_MatrixMRhs_SUPG, ColumnSpace_MatrixMRhs_SUPG, RhsSpace_MatrixMRhs_SUPG,
                                         MatrixMRhsALEAssemble_SUPG, LinCoeffs, NULL);
      
      
    DiscreteFormSUPG = new TDiscreteForm2D(SUPGString, allString,
                             N_Terms_MatrixMARhs_SUPG, Derivatives_MatrixMARhs_SUPG, SpacesNumbers_MatrixMARhs_SUPG,
                             N_Matrices_MatrixMARhs_SUPG, N_Rhs_MatrixMARhs_SUPG,
                             RowSpace_MatrixMARhs_SUPG, ColumnSpace_MatrixMARhs_SUPG,
                             RhsSpace_MatrixMARhs_SUPG, MatrixMARhsALEAssemble_SUPG, LinCoeffs, NULL);
 }
  
  
} // InitializeDiscreteForms
   
   
void InitializeDiscreteForms_HeatLine(TDiscreteForm2D *&DiscreteFormHeatLine, CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";

  
 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormHeatLine = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs, Derivatives_MatrixARhs, SpacesNumbers_MatrixARhs,
                             N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
                             RowSpace_MatrixARhs, ColumnSpace_MatrixARhs,
                             RhsSpace_MatrixARhs, MatrixARhsAssembleHeatLine_Axial3D, LinCoeffs, NULL);
  } 
 else
 {   
  DiscreteFormHeatLine = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs, Derivatives_MatrixARhs, SpacesNumbers_MatrixARhs,
                             N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
                             RowSpace_MatrixARhs, ColumnSpace_MatrixARhs,
                             RhsSpace_MatrixARhs, MatrixARhsAssembleHeatLine, LinCoeffs, NULL);
 }
  
  
} // InitializeDiscreteForms   
   
 int N_Terms_MatrixARhs_CST = 3;
 MultiIndex2D Derivatives_MatrixARhs_CST[3] = { D10, D01, D00};
 int SpacesNumbers_MatrixARhs_CST[3] = { 2, 2, 2  };
 int N_Matrices_MatrixARhs_CST = 7;
 int RowSpace_MatrixARhs_CST[7] = { 2,2,2,2,2,2,2 };
 int ColumnSpace_MatrixARhs_CST[7] = { 2,2,2,2,2,2,2 };
 int N_Rhs_MatrixARhs_CST = 3;
 int RhsSpace_MatrixARhs_CST[3] = { 2,2,2 };
   
   
void InitializeDiscreteForms_CST(TDiscreteForm2D *&DiscreteFormGalerkin, 
				 TDiscreteForm2D *&DiscreteFormSDFEM, 
				 CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SDString[] = "SUPG";

   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs_CST, Derivatives_MatrixARhs_CST, 
			     SpacesNumbers_MatrixARhs_CST,N_Matrices_MatrixARhs_CST,
                              N_Rhs_MatrixARhs_CST,
                             RowSpace_MatrixARhs_CST, ColumnSpace_MatrixARhs_CST,
                             RhsSpace_MatrixARhs_CST, CSTGalerkin, LinCoeffs, NULL);
   
      DiscreteFormSDFEM = new TDiscreteForm2D(SDString, allString,
                             N_Terms_MatrixARhs_CST, Derivatives_MatrixARhs_CST, 
			     SpacesNumbers_MatrixARhs_CST,N_Matrices_MatrixARhs_CST,
                              N_Rhs_MatrixARhs_CST,
                             RowSpace_MatrixARhs_CST, ColumnSpace_MatrixARhs_CST,
                             RhsSpace_MatrixARhs_CST, CST_SDFEM, LinCoeffs, NULL);
  
}   
  
  void InitializeDiscreteForms_CST_Giesekus(TDiscreteForm2D *&DiscreteFormGalerkin, CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";

   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs_CST, Derivatives_MatrixARhs_CST, 
			     SpacesNumbers_MatrixARhs_CST,N_Matrices_MatrixARhs_CST,
                              N_Rhs_MatrixARhs_CST,
                             RowSpace_MatrixARhs_CST, ColumnSpace_MatrixARhs_CST,
                             RhsSpace_MatrixARhs_CST, CST_GiesekusGalerkin, LinCoeffs, NULL);
  
} 

 int N_Terms_MatrixARhs_DFT = 1;
 MultiIndex2D Derivatives_MatrixARhs_DFT[1] = {D00};
 int SpacesNumbers_MatrixARhs_DFT[1] = { 3 };
 int N_Matrices_MatrixARhs_DFT = 3;
 int RowSpace_MatrixARhs_DFT[3] = { 3,3,3};
 int ColumnSpace_MatrixARhs_DFT[3] = { 3,3,3};
 int N_Rhs_MatrixARhs_DFT = 3;
 int RhsSpace_MatrixARhs_DFT[3] = { 3,3,3 };


void InitializeDiscreteForms_DFT(TDiscreteForm2D *&DiscreteFormGalerkin, CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";


   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs_DFT, Derivatives_MatrixARhs_DFT, 
			     SpacesNumbers_MatrixARhs_DFT,N_Matrices_MatrixARhs_DFT,
                              N_Rhs_MatrixARhs_DFT,
                             RowSpace_MatrixARhs_DFT, ColumnSpace_MatrixARhs_DFT,
                             RhsSpace_MatrixARhs_DFT, DFTGalerkin, LinCoeffs, NULL);
   
}   

 int N_Terms_MatrixARhs_NSECST_DEVSS = 8;
 MultiIndex2D Derivatives_MatrixARhs_NSECST_DEVSS[8] = { D10, D01, D00, D00, D10, D01, D00, D00};
 int SpacesNumbers_MatrixARhs_NSECST_DEVSS[8] = { 0, 0, 0, 1, 2, 2, 2, 3  };
 int N_Matrices_MatrixARhs_NSECST_DEVSS = 36;
//  int RowSpace_MatrixARhs_NSECST_DEVSS[36] =    {0,0,0,0,  1,1,  0,0,  2,2,2,2,2,2,2,  0,0,0,0,  2,2,2,2,2,2,  3,3,3,  0,0,0,0,  3,3,3,3};
//  int ColumnSpace_MatrixARhs_NSECST_DEVSS[36] = {0,0,0,0,  0,0,  1,1,  2,2,2,2,2,2,2,  2,2,2,2,  0,0,0,0,0,0,  3,3,3,  3,3,3,3,  0,0,0,0};
 int RowSpace_MatrixARhs_NSECST_DEVSS[36] =    {0,0,0,0, 2,2,2,2,2,2,2, 3,3,3,  1,1,  0,0,  0,0,0,0,  2,2,2,2,2,2,  0,0,0,0,  3,3,3,3};
 int ColumnSpace_MatrixARhs_NSECST_DEVSS[36] = {0,0,0,0, 2,2,2,2,2,2,2, 3,3,3,  0,0,  1,1,  2,2,2,2,  0,0,0,0,0,0,  3,3,3,3,  0,0,0,0};
 int N_Rhs_MatrixARhs_NSECST_DEVSS = 5;
 int RhsSpace_MatrixARhs_NSECST_DEVSS[5] = { 0,0,2,2,2 };


  int N_Terms_MatrixARhs_NSECST_LPS = 7;
 MultiIndex2D Derivatives_MatrixARhs_NSECST_LPS[7] = { D10, D01, D00, D00, D10, D01, D00};
 int SpacesNumbers_MatrixARhs_NSECST_LPS[7] = { 0, 0, 0, 1, 2, 2, 2  };
 int N_Matrices_MatrixARhs_NSECST_LPS = 25;
 int RowSpace_MatrixARhs_NSECST_LPS[25] =    {0,0,0,0, 2,2,2,2,2,2,2, 1,1,  0,0,  0,0,0,0,  2,2,2,2,2,2 };
 int ColumnSpace_MatrixARhs_NSECST_LPS[25] = {0,0,0,0, 2,2,2,2,2,2,2, 0,0,  1,1,  2,2,2,2,  0,0,0,0,0,0 };
 int N_Rhs_MatrixARhs_NSECST_LPS = 5;
 int RhsSpace_MatrixARhs_NSECST_LPS[5] = { 0,0,2,2,2 };
 
 
void InitializeDiscreteForms_NSECST(TDiscreteForm2D *&DiscreteFormGalerkinSUPG, TDiscreteForm2D *&DiscreteFormLPS, CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  
     DiscreteFormGalerkinSUPG = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs_NSECST_DEVSS, Derivatives_MatrixARhs_NSECST_DEVSS, 
			     SpacesNumbers_MatrixARhs_NSECST_DEVSS,N_Matrices_MatrixARhs_NSECST_DEVSS,
                              N_Rhs_MatrixARhs_NSECST_DEVSS,
                             RowSpace_MatrixARhs_NSECST_DEVSS, ColumnSpace_MatrixARhs_NSECST_DEVSS,
                             RhsSpace_MatrixARhs_NSECST_DEVSS, NSEType4_Galerkin_CST_SUPG_DEVSS, LinCoeffs, NULL);
     
     if(TDatabase::ParamDB->Axial3D==0)
     {
          DiscreteFormLPS = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs_NSECST_LPS, Derivatives_MatrixARhs_NSECST_LPS, 
			     SpacesNumbers_MatrixARhs_NSECST_LPS,N_Matrices_MatrixARhs_NSECST_LPS,
                              N_Rhs_MatrixARhs_NSECST_LPS,
                             RowSpace_MatrixARhs_NSECST_LPS, ColumnSpace_MatrixARhs_NSECST_LPS,
                             RhsSpace_MatrixARhs_NSECST_LPS, NSEType4_Galerkin_CST_Galerkin, LinCoeffs, NULL);
     }
     else if(TDatabase::ParamDB->Axial3D==1)
     {
          DiscreteFormLPS = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixARhs_NSECST_LPS, Derivatives_MatrixARhs_NSECST_LPS, 
			     SpacesNumbers_MatrixARhs_NSECST_LPS,N_Matrices_MatrixARhs_NSECST_LPS,
                              N_Rhs_MatrixARhs_NSECST_LPS,
                             RowSpace_MatrixARhs_NSECST_LPS, ColumnSpace_MatrixARhs_NSECST_LPS,
                             RhsSpace_MatrixARhs_NSECST_LPS, NSEType4_Galerkin_CST_Galerkin_Axial3D, LinCoeffs, NULL);
     }
}

 int Time_N_Terms_MatrixARhs_NSECST_DEVSS = 8;
 MultiIndex2D Time_Derivatives_MatrixARhs_NSECST_DEVSS[8] = { D10, D01, D00, D00, D10, D01, D00, D00};
 int Time_SpacesNumbers_MatrixARhs_NSECST_DEVSS[8] = { 0, 0, 0, 1, 2, 2, 2, 3  };
 int Time_N_Matrices_MatrixARhs_NSECST_DEVSS = 41;
 int Time_RowSpace_MatrixARhs_NSECST_DEVSS[41] =    {0,0,0,0, 0,0, 2,2,2,2,2,2,2, 2,2,2, 3,3,3,  1,1,  0,0,  0,0,0,0,  2,2,2,2,2,2,  0,0,0,0,  3,3,3,3};
 int Time_ColumnSpace_MatrixARhs_NSECST_DEVSS[41] = {0,0,0,0, 0,0, 2,2,2,2,2,2,2, 2,2,2, 3,3,3,  0,0,  1,1,  2,2,2,2,  0,0,0,0,0,0,  3,3,3,3,  0,0,0,0};
 int Time_N_Rhs_MatrixARhs_NSECST_DEVSS = 8;
 int Time_RhsSpace_MatrixARhs_NSECST_DEVSS[8] = { 0,0,2,2,2, 2,2,2 };

  int Time_N_Terms_MatrixARhs_NSECST_LPS = 7;
 MultiIndex2D Time_Derivatives_MatrixARhs_NSECST_LPS[7] = { D10, D01, D00, D00, D10, D01, D00};
 int Time_SpacesNumbers_MatrixARhs_NSECST_LPS[7] = { 0, 0, 0, 1, 2, 2, 2  };
 int Time_N_Matrices_MatrixARhs_NSECST_LPS = 30;
 int Time_RowSpace_MatrixARhs_NSECST_LPS[30] =    {0,0,0,0,  0,0,  2,2,2,2,2,2,2,  2,2,2,  1,1,  0,0,  0,0,0,0,  2,2,2,2,2,2 };
 int Time_ColumnSpace_MatrixARhs_NSECST_LPS[30] = {0,0,0,0,  0,0,  2,2,2,2,2,2,2,  2,2,2,  0,0,  1,1,  2,2,2,2,  0,0,0,0,0,0 };
 int Time_N_Rhs_MatrixARhs_NSECST_LPS = 8;
 int Time_RhsSpace_MatrixARhs_NSECST_LPS[8] = { 0,0,2,2,2, 2,2,2 };

  int Time_N_Terms_Matrix_NL_NSECST = 6;
 MultiIndex2D Time_Derivatives_Matrix_NL_NSECST[6] = { D10, D01, D00, D10, D01, D00};
 int Time_SpacesNumbers_Matrix_NL_NSECST[6] = { 0, 0, 0, 2, 2, 2};
 int Time_N_Matrices_Matrix_NL_NSECST = 15;
 int Time_RowSpace_Matrix_NL_NSECST[15] =    {0,0, 2,2,2,2,2,2,2, 2,2,2,2,2,2};
 int Time_ColumnSpace_Matrix_NL_NSECST[15] = {0,0, 2,2,2,2,2,2,2, 0,0,0,0,0,0};
 int Time_N_Rhs_MatrixA_NL_NSECST = 3;
 int Time_RhsSpace_Matrix_NL_NSECST[3] = {2,2,2};
 
   int Time_N_Terms_Matrix_Rhs_NSECST = 4;
 MultiIndex2D Time_Derivatives_Matrix_Rhs_NSECST[4] = { D00, D10, D01, D00};
 int Time_SpacesNumbers_Matrix_Rhs_NSECST[4] = { 0, 2, 2, 2};
 int Time_N_Matrices_Matrix_Rhs_NSECST = 0;
 int *Time_RowSpace_Matrix_Rhs_NSECST =  NULL;
 int *Time_ColumnSpace_Matrix_Rhs_NSECST = NULL;
 int Time_N_Rhs_MatrixA_Rhs_NSECST = 5;
 int Time_RhsSpace_Matrix_Rhs_NSECST[5] = {0,0, 2,2,2};

void InitializeDiscreteForms_TNSECST(TDiscreteForm2D *&DiscreteFormGalerkinSUPG, TDiscreteForm2D *&DiscreteFormLPS, 
				     TDiscreteForm2D *&DiscreteFormNLGalerkinSUPG, TDiscreteForm2D *&DiscreteFormNLLPS, 
				     TDiscreteForm2D *&DiscreteFormRHSGalerkinSUPG,  TDiscreteForm2D *&DiscreteFormRHSLPS,
				     CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char nonlinearString[] = "nonlinear";
  char allString[] = "all";
  char rhsString[] = "rhs";
  
     
     DiscreteFormGalerkinSUPG = new TDiscreteForm2D(GalerkinString, allString,
                             Time_N_Terms_MatrixARhs_NSECST_DEVSS, Time_Derivatives_MatrixARhs_NSECST_DEVSS, 
			     Time_SpacesNumbers_MatrixARhs_NSECST_DEVSS, Time_N_Matrices_MatrixARhs_NSECST_DEVSS,
                              Time_N_Rhs_MatrixARhs_NSECST_DEVSS,
                             Time_RowSpace_MatrixARhs_NSECST_DEVSS, Time_ColumnSpace_MatrixARhs_NSECST_DEVSS,
                             Time_RhsSpace_MatrixARhs_NSECST_DEVSS, Time_NSEType4_Galerkin_CST_SUPG_DEVSS, LinCoeffs, NULL);
      
     DiscreteFormLPS = new TDiscreteForm2D(GalerkinString, allString,
                             Time_N_Terms_MatrixARhs_NSECST_LPS, Time_Derivatives_MatrixARhs_NSECST_LPS, 
			     Time_SpacesNumbers_MatrixARhs_NSECST_LPS, Time_N_Matrices_MatrixARhs_NSECST_LPS,
                              Time_N_Rhs_MatrixARhs_NSECST_LPS,
                             Time_RowSpace_MatrixARhs_NSECST_LPS, Time_ColumnSpace_MatrixARhs_NSECST_LPS,
                             Time_RhsSpace_MatrixARhs_NSECST_LPS, Time_NSEType4_Galerkin_CST_Galerkin, LinCoeffs, NULL);
     
      DiscreteFormNLGalerkinSUPG = new TDiscreteForm2D(GalerkinString, nonlinearString,
                             Time_N_Terms_Matrix_NL_NSECST, Time_Derivatives_Matrix_NL_NSECST, 
			     Time_SpacesNumbers_Matrix_NL_NSECST, Time_N_Matrices_Matrix_NL_NSECST,
                              Time_N_Rhs_MatrixA_NL_NSECST,
                             Time_RowSpace_Matrix_NL_NSECST, Time_ColumnSpace_Matrix_NL_NSECST,
                             Time_RhsSpace_Matrix_NL_NSECST, Time_NSEType4_Galerkin_CST_SUPG_DEVSS_NLTerms, LinCoeffs, NULL);
      
      DiscreteFormNLLPS = new TDiscreteForm2D(GalerkinString, nonlinearString,
                             Time_N_Terms_Matrix_NL_NSECST, Time_Derivatives_Matrix_NL_NSECST, 
			     Time_SpacesNumbers_Matrix_NL_NSECST, Time_N_Matrices_Matrix_NL_NSECST,
                              Time_N_Rhs_MatrixA_NL_NSECST,
                             Time_RowSpace_Matrix_NL_NSECST, Time_ColumnSpace_Matrix_NL_NSECST,
                             Time_RhsSpace_Matrix_NL_NSECST, Time_NSEType4_Galerkin_CST_Galerkin_NLTerms, LinCoeffs, NULL);
      
      DiscreteFormRHSGalerkinSUPG = new TDiscreteForm2D(GalerkinString, rhsString,
                             Time_N_Terms_Matrix_Rhs_NSECST, Time_Derivatives_Matrix_Rhs_NSECST, 
			     Time_SpacesNumbers_Matrix_Rhs_NSECST, Time_N_Matrices_Matrix_Rhs_NSECST,
                              Time_N_Rhs_MatrixA_Rhs_NSECST,
                             Time_RowSpace_Matrix_Rhs_NSECST, Time_ColumnSpace_Matrix_Rhs_NSECST,
                             Time_RhsSpace_Matrix_Rhs_NSECST, Time_NSEType4_Galerkin_CST_SUPG_RhsOnly, LinCoeffs, NULL);
     
     DiscreteFormRHSLPS = new TDiscreteForm2D(GalerkinString, rhsString,
                             Time_N_Terms_Matrix_Rhs_NSECST, Time_Derivatives_Matrix_Rhs_NSECST, 
			     Time_SpacesNumbers_Matrix_Rhs_NSECST, Time_N_Matrices_Matrix_Rhs_NSECST,
                              Time_N_Rhs_MatrixA_Rhs_NSECST,
                             Time_RowSpace_Matrix_Rhs_NSECST, Time_ColumnSpace_Matrix_Rhs_NSECST,
                             Time_RhsSpace_Matrix_Rhs_NSECST, Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly, LinCoeffs, NULL);
     
}

void InitializeDiscreteForms_2PhaseAxial3D_TNSECST(TDiscreteForm2D *&DiscreteFormLPS, 
				    TDiscreteForm2D *&DiscreteFormNLLPS, 
				     TDiscreteForm2D *&DiscreteFormRHSLPS,
				     TDiscreteForm2D *&DiscreteFormGrid,
				     CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char nonlinearString[] = "nonlinear";
  char allString[] = "all";
  char rhsString[] = "rhs";
    char GridString[] = "Grid";
  
     DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
                            GridN_Terms, GridDerivatives, GridSpaceNumbers,
                            GridN_Matrices, GridN_Rhs,
                            GridRowSpace, GridColumnSpace,
                            GridRhsSpace, GridAssemble4,
                            GridCoeffs, NULL);

     DiscreteFormLPS = new TDiscreteForm2D(GalerkinString, allString,
                             Time_N_Terms_MatrixARhs_NSECST_LPS, Time_Derivatives_MatrixARhs_NSECST_LPS, 
			     Time_SpacesNumbers_MatrixARhs_NSECST_LPS, Time_N_Matrices_MatrixARhs_NSECST_LPS,
                              Time_N_Rhs_MatrixARhs_NSECST_LPS,
                             Time_RowSpace_MatrixARhs_NSECST_LPS, Time_ColumnSpace_MatrixARhs_NSECST_LPS,
                             Time_RhsSpace_MatrixARhs_NSECST_LPS, Time_NSEType4_Galerkin_CST_Galerkin_2PhaseAxial3D , LinCoeffs, NULL);

      DiscreteFormNLLPS = new TDiscreteForm2D(GalerkinString, nonlinearString,
                             Time_N_Terms_Matrix_NL_NSECST, Time_Derivatives_Matrix_NL_NSECST, 
			     Time_SpacesNumbers_Matrix_NL_NSECST, Time_N_Matrices_Matrix_NL_NSECST,
                              Time_N_Rhs_MatrixA_NL_NSECST,
                             Time_RowSpace_Matrix_NL_NSECST, Time_ColumnSpace_Matrix_NL_NSECST,
                             Time_RhsSpace_Matrix_NL_NSECST, Time_NSEType4_Galerkin_CST_Galerkin_NLTerms_2PhaseAxial3D, LinCoeffs, NULL);
      
     DiscreteFormRHSLPS = new TDiscreteForm2D(GalerkinString, rhsString,
                             Time_N_Terms_Matrix_Rhs_NSECST, Time_Derivatives_Matrix_Rhs_NSECST, 
			     Time_SpacesNumbers_Matrix_Rhs_NSECST, Time_N_Matrices_Matrix_Rhs_NSECST,
                              Time_N_Rhs_MatrixA_Rhs_NSECST,
                             Time_RowSpace_Matrix_Rhs_NSECST, Time_ColumnSpace_Matrix_Rhs_NSECST,
                             Time_RhsSpace_Matrix_Rhs_NSECST, Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly_2PhaseAxial3D, LinCoeffs, NULL);
     
}


void InitializeDiscreteForms_ImpDropAxial3D_TNSECST(TDiscreteForm2D *&DiscreteFormLPS, 
				    TDiscreteForm2D *&DiscreteFormNLLPS, 
				     TDiscreteForm2D *&DiscreteFormRHSLPS,
				     TDiscreteForm2D *&DiscreteFormGrid,
				     CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char nonlinearString[] = "nonlinear";
  char allString[] = "all";
  char rhsString[] = "rhs";
    char GridString[] = "Grid";
  
     DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
                            GridN_Terms, GridDerivatives, GridSpaceNumbers,
                            GridN_Matrices, GridN_Rhs,
                            GridRowSpace, GridColumnSpace,
                            GridRhsSpace, GridAssemble4,
                            GridCoeffs, NULL);

     DiscreteFormLPS = new TDiscreteForm2D(GalerkinString, allString,
                             Time_N_Terms_MatrixARhs_NSECST_LPS, Time_Derivatives_MatrixARhs_NSECST_LPS, 
			     Time_SpacesNumbers_MatrixARhs_NSECST_LPS, Time_N_Matrices_MatrixARhs_NSECST_LPS,
                              Time_N_Rhs_MatrixARhs_NSECST_LPS,
                             Time_RowSpace_MatrixARhs_NSECST_LPS, Time_ColumnSpace_MatrixARhs_NSECST_LPS,
                             Time_RhsSpace_MatrixARhs_NSECST_LPS, Time_NSEType4_Galerkin_CST_Galerkin_ImpDropAxial3D, LinCoeffs, NULL);

      DiscreteFormNLLPS = new TDiscreteForm2D(GalerkinString, nonlinearString,
                             Time_N_Terms_Matrix_NL_NSECST, Time_Derivatives_Matrix_NL_NSECST, 
			     Time_SpacesNumbers_Matrix_NL_NSECST, Time_N_Matrices_Matrix_NL_NSECST,
                              Time_N_Rhs_MatrixA_NL_NSECST,
                             Time_RowSpace_Matrix_NL_NSECST, Time_ColumnSpace_Matrix_NL_NSECST,
                             Time_RhsSpace_Matrix_NL_NSECST, Time_NSEType4_Galerkin_CST_Galerkin_NLTerms_ImpDropAxial3D, LinCoeffs, NULL);
      
     DiscreteFormRHSLPS = new TDiscreteForm2D(GalerkinString, rhsString,
                             Time_N_Terms_Matrix_Rhs_NSECST, Time_Derivatives_Matrix_Rhs_NSECST, 
			     Time_SpacesNumbers_Matrix_Rhs_NSECST, Time_N_Matrices_Matrix_Rhs_NSECST,
                              Time_N_Rhs_MatrixA_Rhs_NSECST,
                             Time_RowSpace_Matrix_Rhs_NSECST, Time_ColumnSpace_Matrix_Rhs_NSECST,
                             Time_RhsSpace_Matrix_Rhs_NSECST, Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly_ImpDropAxial3D, LinCoeffs, NULL);
     
}




void InitializeDiscreteForms_Moving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormSUPG, CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SUPGString[] = "SUPG";
  
 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble_Axial3D, LinCoeffs, NULL);  
   
   
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;  
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented InitializeDiscreteForms !!!!!!"<<endl;     
    DiscreteFormSUPG = NULL; 
  } 
 else
 {   
  DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble, LinCoeffs, NULL);
  
  
  DiscreteFormSUPG = new TDiscreteForm2D(SUPGString, allString,
                             N_Terms_MatrixMARhs_SUPG, Derivatives_MatrixMARhs_SUPG, SpacesNumbers_MatrixMARhs_SUPG,
                             N_Matrices_MatrixMARhs_SUPG, N_Rhs_MatrixMARhs_SUPG,
                             RowSpace_MatrixMARhs_SUPG, ColumnSpace_MatrixMARhs_SUPG,
                             RhsSpace_MatrixMARhs_SUPG, MatrixMARhsALEAssemble_SUPG, LinCoeffs, NULL);
 }
  
  
} // InitializeDiscreteForms
   
   
void InitializeDiscreteForms_Stationary( TDiscreteForm2D *&DiscreteFormUpwind,  TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormSDFEM,
                                TDiscreteForm2D *&DiscreteFormGLS, CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SUPGString[] = "SUPG";
  
    ManipulateFct2D *manipulate;
  if (TDatabase::ParamDB->SDFEM_NORM_B==0)
    manipulate = linfb;
  else
    manipulate = ave_l2b_quad_points;
  
 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                  N_Terms, Derivatives, SpacesNumbers,
                                  CD_N_Matrices, CD_N_Rhs, CD_RowSpace,
                                  CD_ColumnSpace, CD_RhsSpace, 
                                  (AssembleFctParam2D*)BilinearAssemble_Axial3D,
                                  LinCoeffs, NULL); 
   
  
   cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;
   cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;   
   cout<< "InitializeDiscreteForms_Stationary"<<endl;  
   DiscreteFormUpwind = NULL;
   DiscreteFormSDFEM = NULL;
   DiscreteFormGLS = NULL;
  } 
 else
 {
   
     DiscreteFormUpwind = new TDiscreteForm2D(GalerkinString, allString,
                                  N_Terms, Derivatives, SpacesNumbers,
                                  CD_N_Matrices, CD_N_Rhs, CD_RowSpace,
                                  CD_ColumnSpace, CD_RhsSpace,
                                  (AssembleFctParam2D*)BilinearAssemble_UPW2,
                                  LinCoeffs, NULL);  
   
     DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                  N_Terms, Derivatives, SpacesNumbers,
                                  CD_N_Matrices, CD_N_Rhs, CD_RowSpace,
                                  CD_ColumnSpace, CD_RhsSpace, 
                                  (AssembleFctParam2D*)BilinearAssembleGalerkin,
                                  LinCoeffs, NULL);
     
     DiscreteFormSDFEM = new TDiscreteForm2D(SUPGString, SUPGString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD, 
                                             CD_N_Matrices, CD_N_Rhs, CD_RowSpace, CD_ColumnSpace, CD_RhsSpace,
                                             BilinearAssemble_SD, LinCoeffs,  manipulate);
     
     DiscreteFormGLS = new TDiscreteForm2D(SUPGString, SUPGString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD, 
                                             CD_N_Matrices, CD_N_Rhs, CD_RowSpace, CD_ColumnSpace, CD_RhsSpace,
                                             BilinearAssemble_GLS, LinCoeffs,  manipulate);     
     
     
 }  
} // InitializeDiscreteForms
     

void  InitializeDiscreteForms(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormNLGalerkin,
                              CoeffFct2D *LinCoeffs)
{

  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";




//   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
//                                              NSType4N_Terms, NSType4Derivatives, 
//                                              NSType4SpaceNumbers,
//                                              NSType4N_Matrices, NSType4N_Rhs, 
//                                              NSType4RowSpace, NSType4ColumnSpace,
//                                              NSType4RhsSpace, NSType4GalerkinAxialSymm3D, LinCoeffs, NULL);
// 
//   DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
//                                                NSType4NLN_Terms, NSType4NLDerivatives, 
//                                                NSType4NLSpaceNumbers,
//                                                NSType4NLN_Matrices, NSType4NLN_Rhs, 
//                                                NSType4NLRowSpace, NSType4NLColumnSpace,
//                                                NSType4NLRhsSpace, NSType3_4NLGalerkinAxialSymm3D, LinCoeffs, NULL);
// 
//   DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
// 					   NSType4N_Terms, NSType4Derivatives, 
// 					   NSType4SpaceNumbers,
// 					   NSType4N_Matrices, NSType4N_Rhs, 
// 					   NSType4RowSpace, NSType4ColumnSpace,
// 					   NSType4RhsSpace, NSType4UpwindAxialSymm3D, LinCoeffs, NULL);
//   
//   DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
// 					     NSType4NLN_Terms, NSType4NLDerivatives, 
// 					     NSType4NLSpaceNumbers,
// 					     NSType4NLN_Matrices, NSType4NLN_Rhs, 
// 					     NSType4NLRowSpace, NSType4NLColumnSpace,
// 					     NSType4NLRhsSpace, NSType3_4NLUpwindAxialSymm3D, LinCoeffs, NULL);
}


void InitializeDiscreteFormGrid(TDiscreteForm2D *&DiscreteFormGrid, CoeffFct2D *GridCoeffs)
{
  char GridString[] = "Grid";
  char allString[] = "all";
  
  DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
            GridN_Terms, GridDerivatives, GridSpaceNumbers,
            GridN_Matrices, GridN_Rhs,
            GridRowSpace, GridColumnSpace,
            GridRhsSpace, GridAssemble4,
             GridCoeffs, NULL);    
    
    
}
    
void InitializeDiscreteForms_2PhaseAxial3D(
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormGrid,
  CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs, int NSTYPE)
{

  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char nonlinearString[] = "nonlinear";
  char rhsString[] = "rhs";
  char GridString[] = "Grid";

DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs,
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, MovingNSRHS, LinCoeffs, NULL);

 DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
            GridN_Terms, GridDerivatives, GridSpaceNumbers,
            GridN_Matrices, GridN_Rhs,
            GridRowSpace, GridColumnSpace,
            GridRhsSpace, GridAssemble4,
             GridCoeffs, NULL);


  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    case 2:
    case 3:
    cout<< " DiscreteFormGalerkin_Psep not implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
    cout<< " DiscreteFormUpwind_Psep not implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;


    cout<< " DiscreteForm noting is implemented for 2Phase NSTYPE 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
    exit(-1);
  
    break;

    case 4:
      if(TDatabase::ParamDB->LAPLACETYPE == 0)
      {

         cout<< " DiscreteForm noting is implemented for 2Phase NSTYPE 4 grad:grad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
        exit(-1);
  
      }
      else
      {

       DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                   TimeNSType4N_Matrices, TimeNSType4N_Rhs,
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD_2PhaseAxial3D, LinCoeffs, NULL);

       DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                   TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs,
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD_2PhaseAxial3D, LinCoeffs, NULL);

      }
    break;
  } // endswitch
}



#endif
