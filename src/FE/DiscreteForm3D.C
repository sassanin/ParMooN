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
// %W% %G%
// 
// Class:       TDiscreteForm3D
// Purpose:     assemble a couple of matrices and right-hand side at once
//
// Author:      Gunar Matthies (02.10.98)
//
// History:     start of implementation 02.10.98 (Gunar Matthies)
//
// =======================================================================

#include <Database.h>
#include <FEDatabase3D.h>
#include <DiscreteForm3D.h>
#include <NSE3D_Param.h>
#include <NSE3D_FixPo.h>
#include <NSE3D_FixPoSkew.h>
#include <NSE3D_Friction_FixPo.h>
#include <NSE3D_Newton.h>
#include <TNSE3D_FixPo.h>
#include <TNSE3D_Newton.h>
//#include <TNSE3D_Routines.h>
#include <string.h>
#include <stdlib.h>
#include <ConvDiff3D.h>
#include <TCD3D.h>
#include <ALE_3D.h>

/** constructor with vector initialization */
TDiscreteForm3D::TDiscreteForm3D(char *name, char *description, int n_terms,
        MultiIndex3D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFct3D *assemble, CoeffFct3D *coeffs,
        ManipulateFct3D *manipulate)
{
  int i, j, max;
  MultiIndex3D alpha;

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
    Needs2ndDerivatives[i] = false;

  for(i=0;i<N_Terms;i++)
    {
      alpha = Derivatives[i];
      j = FESpaceNumber[i];
      if(alpha == D200 || alpha == D110 || alpha == D101 || 
         alpha == D020 || alpha == D011 || alpha == D002)
        Needs2ndDerivatives[j] = true;
    }

#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==TDatabase::ParamDB->Par_P0)
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
TDiscreteForm3D::TDiscreteForm3D(char *name, char *description, int n_terms,
        MultiIndex3D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFctParam3D *assembleparam, CoeffFct3D *coeffs,
        ManipulateFct3D *manipulate)
{
  int i, j, max;
  MultiIndex3D alpha;

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
    Needs2ndDerivatives[i] = false;

  for(i=0;i<N_Terms;i++)
  {
    alpha = Derivatives[i];
    j = FESpaceNumber[i];
    if(alpha == D200 || alpha == D110 || alpha == D101 || 
       alpha == D020 || alpha == D011 || alpha == D002) 
      Needs2ndDerivatives[j] = true;
  }


#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==TDatabase::ParamDB->Par_P0)
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

TDiscreteForm3D::~TDiscreteForm3D()
{
  delete AllOrigValues;
  delete OrigValues;
  delete Needs2ndDerivatives;
  delete Name;
  delete Description;
}

void TDiscreteForm3D::GetLocalForms(int N_Points, double *weights, 
                                    double *AbsDetjk, double hK, 
                                    double *X, double *Y, double *Z,
                                    int *N_BaseFuncts, BaseFunct3D *BaseFuncts, 
                                    double **Parameters, double **AuxArray,
                                    TBaseCell *Cell,
                                    int n_matrices, int n_rhs,
                                    double ***LocMatrix, double **LocRhs)
{
  int i,j,k,l, N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  double Mult, *Coeff, *Param;
 
//   cout << "in TDiscreteForm::GetLocalForms " << n_matrices<< endl;

  for(i=0;i<n_matrices;i++)
  {
      
    N_Rows =-50;
    N_Columns =-50;
      
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
//   AuxArray[0][0] = Cell->GetPhase_ID();
  AuxArray[0][0] = Cell->GetRegionID();
  AuxArray[0][1] = hK;
  
//   AuxArray[0][2] = Cell->GetVertex(0)-> GetX(); 
//   AuxArray[0][3] = Cell->GetVertex(1)-> GetX();  
//   AuxArray[0][4] = Cell->GetVertex(2)-> GetX();
//   AuxArray[0][5] = Cell->GetVertex(3)-> GetX();   
// *****************************************************

  if(Coeffs)
    Coeffs(N_Points, X, Y, Z, Parameters, AuxArray);

  if(Manipulate)
    Manipulate(N_Points, AuxArray, Parameters, Cell);

  for(i=0;i<N_Terms;i++)
  {
    AllOrigValues[i] = 
      TFEDatabase3D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]], 
                                          Derivatives[i]);
    if(!(AllOrigValues[i]))
    {
      ErrMsg("second derivatives not yet supported. Exiting");
      exit(1);
    }
  }  
  
  for(i=0;i<N_Points;i++)
  {
    Mult = weights[i]*AbsDetjk[i];
       
    Coeff = AuxArray[i];
    Coeff[19] = AbsDetjk[i];
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

/******************************************************************************/
//
// Routines for initializing discrete forms
//
/******************************************************************************/

void InitializeDiscreteForms(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormSDFEM,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormSmagorinsky,
  TDiscreteForm3D *&DiscreteFormVMSProjection,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLSDFEM,
  TDiscreteForm3D *&DiscreteFormNLUpwind,
  TDiscreteForm3D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm3D *&DiscreteFormNLVMSProjection,
  TDiscreteForm3D *&DiscreteFormNLSDFEM_DivDiv,
  TDiscreteForm3D *&DiscreteFormPressSep,
  TDiscreteForm3D *&DiscreteFormPressSepAuxProb,
  TDiscreteForm3D *&DiscreteFormNSRFBRhs,
  CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char Galerkin[] = "Galerkin";
  char all[] = "all";
  char Upwind[] = "Upwind";
  char Smagorinsky[] = "Smagorinsky";
  char nonlinear[] = "nonlinear";
  char Sdfem[] = "SDFEM";

  DiscreteFormNLSDFEM_DivDiv = NULL;
  DiscreteFormVMSProjection = NULL;
  DiscreteFormNLVMSProjection = NULL;

  DiscreteFormPressSep = new TDiscreteForm3D(Galerkin, all,
					     NSPressSepN_Terms, NSPressSepDerivatives, 
					     NSPressSepSpaceNumbers,
					     NSPressSepN_Matrices, NSPressSepN_Rhs, 
					     NSPressSepRowSpace, NSPressSepColumnSpace,
					     NSPressSepRhsSpace, NSPressSep, LinCoeffs, NULL);
  DiscreteFormPressSepAuxProb = new TDiscreteForm3D(Galerkin, all,
					     NSPressSepAuxProbN_Terms, NSPressSepAuxProbDerivatives, 
					     NSPressSepAuxProbSpaceNumbers,
					     NSPressSepAuxProbN_Matrices, NSPressSepAuxProbN_Rhs, 
					     NSPressSepAuxProbRowSpace, NSPressSepAuxProbColumnSpace,
					     NSPressSepAuxProbRhsSpace, NSPressSepAuxProb, LinCoeffs, NULL);
 
  DiscreteFormNSRFBRhs = new TDiscreteForm3D(Galerkin, all,
					     NSRFBRhsN_Terms, NSRFBRhsDerivatives, 
					     NSRFBRhsSpaceNumbers,
					     NSRFBRhsN_Matrices, NSRFBRhsN_Rhs, 
					     NSRFBRhsRowSpace, NSRFBRhsColumnSpace,
					     NSRFBRhsRhsSpace, NSRFBRhs3D, LinCoeffs, NULL);
 //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
  {
    switch(NSTYPE)
    {
	case 1:
	    if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
	    {
		DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
							   NSType1N_Terms, NSType1Derivatives, 
							   NSType1SpaceNumbers,
							   NSType1N_Matrices, NSType1N_Rhs, 
							   NSType1RowSpace, NSType1ColumnSpace,
							   NSType1RhsSpace, NSType1Galerkin3D, LinCoeffs, NULL);
		
		DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
							NSType1N_Terms, NSType1Derivatives, 
							NSType1SpaceNumbers,
							NSType1N_Matrices, NSType1N_Rhs, 
							NSType1RowSpace, NSType1ColumnSpace,
							NSType1RhsSpace, NSType1SDFEM3D, LinCoeffs, NULL);
		
		DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
							 NSType1N_Terms, NSType1Derivatives, 
							 NSType1SpaceNumbers,
							 NSType1N_Matrices, NSType1N_Rhs, 
							 NSType1RowSpace, NSType1ColumnSpace,
							 NSType1RhsSpace, NSType1Upwind3D, LinCoeffs, NULL);
		
		DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
							      NSType1N_Terms, NSType1Derivatives, 
							      NSType1SpaceNumbers,
							      NSType1N_Matrices, NSType1N_Rhs, 
							      NSType1RowSpace, NSType1ColumnSpace,
							      NSType1RhsSpace, NSType1Smagorinsky3D, LinCoeffs, NULL);
		
		DiscreteFormVMSProjection = new TDiscreteForm3D(Smagorinsky, all, 
								NSType1VMSProjectionN_Terms, NSType1VMSProjectionDerivatives, 
								NSType1VMSProjectionSpaceNumbers,
								NSType1VMSProjectionN_Matrices, NSType1N_Rhs, 
								NSType1VMSProjectionRowSpace, NSType1VMSProjectionColumnSpace,
								NSType1RhsSpace, NSType1VMSProjection3D, LinCoeffs, NULL);

		DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
							     NSType1NLN_Terms, NSType1NLDerivatives, 
							     NSType1NLSpaceNumbers,
							     NSType1NLN_Matrices, NSType1NLN_Rhs, 
							     NSType1NLRowSpace, NSType1NLColumnSpace,
							     NSType1NLRhsSpace, NSType1_2NLGalerkin3D, LinCoeffs, NULL);
		
		DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
							  NSType1NLN_Terms, NSType1NLDerivatives, 
							  NSType1NLSpaceNumbers,
							  NSType1NLN_Matrices, NSType1NLSDFEMN_Rhs, 
							  NSType1NLRowSpace, NSType1NLColumnSpace,
							  NSType1NLSDFEMRhsSpace, NSType1NLSDFEM3D, LinCoeffs, NULL);
		
		DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
							   NSType1NLN_Terms, NSType1NLDerivatives, 
							   NSType1NLSpaceNumbers,
							   NSType1NLN_Matrices, NSType1NLN_Rhs, 
							   NSType1NLRowSpace, NSType1NLColumnSpace,
							   NSType1NLRhsSpace, NSType1_2NLUpwind3D, LinCoeffs, NULL);
		
		DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
								NSType1NLN_Terms, NSType1NLDerivatives, 
								NSType1NLSpaceNumbers,
								NSType1NLN_Matrices, NSType1NLN_Rhs, 
								NSType1NLRowSpace, NSType1NLColumnSpace,
								NSType1NLRhsSpace, NSType1_2NLSmagorinsky3D, LinCoeffs, NULL);
		DiscreteFormNLVMSProjection = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType1_2NLVMSProjectionN_Terms, NSType1_2NLVMSProjectionDerivatives, 
                  NSType1_2NLVMSProjectionSpaceNumbers,
                  NSType1_2NLVMSProjectionN_Matrices, NSType1NLN_Rhs,
                  NSType1_2NLVMSProjectionRowSpace, NSType1_2NLVMSProjectionColumnSpace,
                  NSType1NLRhsSpace, NSType1_2NLVMSProjection3D, LinCoeffs, NULL);

	    }
	    else // skew--symmetric
	    {
		DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
							   NSType1N_Terms, NSType1Derivatives, 
							   NSType1SpaceNumbers,
							   NSType1N_Matrices, NSType1N_Rhs, 
							   NSType1RowSpace, NSType1ColumnSpace,
							   NSType1RhsSpace, NSType1GalerkinSkew3D, LinCoeffs, NULL);
		
		DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
							NSType1N_Terms, NSType1Derivatives, 
							NSType1SpaceNumbers,
							NSType1N_Matrices, NSType1N_Rhs, 
							NSType1RowSpace, NSType1ColumnSpace,
							NSType1RhsSpace, NSType1SDFEMSkew3D, LinCoeffs, NULL);
		
		DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
							 NSType1N_Terms, NSType1Derivatives, 
							 NSType1SpaceNumbers,
							 NSType1N_Matrices, NSType1N_Rhs, 
							 NSType1RowSpace, NSType1ColumnSpace,
							 NSType1RhsSpace, NSType1Upwind3D, LinCoeffs, NULL);
		
		DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
							      NSType1N_Terms, NSType1Derivatives, 
							      NSType1SpaceNumbers,
							      NSType1N_Matrices, NSType1N_Rhs, 
							      NSType1RowSpace, NSType1ColumnSpace,
							      NSType1RhsSpace, NSType1SmagorinskySkew3D, LinCoeffs, NULL);
		
		DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
							     NSType1NLN_Terms, NSType1NLDerivatives, 
							     NSType1NLSpaceNumbers,
							     NSType1NLN_Matrices, NSType1NLN_Rhs, 
							     NSType1NLRowSpace, NSType1NLColumnSpace,
							     NSType1NLRhsSpace, NSType1_2NLGalerkinSkew3D, LinCoeffs, NULL);
		
		DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
							  NSType1NLN_Terms, NSType1NLDerivatives, 
							  NSType1NLSpaceNumbers,
							  NSType1NLN_Matrices, NSType1NLSDFEMN_Rhs, 
							  NSType1NLRowSpace, NSType1NLColumnSpace,
							  NSType1NLSDFEMRhsSpace, NSType1NLSDFEMSkew3D, LinCoeffs, NULL);
		
		DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
							   NSType1NLN_Terms, NSType1NLDerivatives, 
							   NSType1NLSpaceNumbers,
							   NSType1NLN_Matrices, NSType1NLN_Rhs, 
							   NSType1NLRowSpace, NSType1NLColumnSpace,
							   NSType1NLRhsSpace, NSType1_2NLUpwind3D, LinCoeffs, NULL);
		
		DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
								NSType1NLN_Terms, NSType1NLDerivatives, 
								NSType1NLSpaceNumbers,
								NSType1NLN_Matrices, NSType1NLN_Rhs, 
								NSType1NLRowSpace, NSType1NLColumnSpace,
								NSType1NLRhsSpace, NSType1_2NLSmagorinskySkew3D, LinCoeffs, NULL);
	    }

        break;
      case 2:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
        DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Galerkin3D, LinCoeffs, NULL);

        DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEM3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind3D, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Smagorinsky3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkin3D, LinCoeffs, NULL);
      
        DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEM3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind3D, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinsky3D, LinCoeffs, NULL);
	}
	else // skew-symmetric
         {
        DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2GalerkinSkew3D, LinCoeffs, NULL);

        DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEMSkew3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind3D, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2SmagorinskySkew3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkinSkew3D, LinCoeffs, NULL);
      
        DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEMSkew3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind3D, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinskySkew3D, LinCoeffs, NULL);
	 }
       break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
        // (grad, grad)
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Galerkin3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Smagorinsky3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky3D, LinCoeffs, NULL);
	  }
	  else // skew-symmetric
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkew3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew3D, LinCoeffs, NULL);
	  }
        }
        else  
        {
          // (D(u):D(v))
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);
	  }
	  else // skew-symmetric
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDSkew3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDDSkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDDSkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDDSkew3D, LinCoeffs, NULL);
	  }

        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Galerkin3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEM3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Smagorinsky3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEM3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky3D, LinCoeffs, NULL);

	  DiscreteFormNLSDFEM_DivDiv = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSD_DivDiv_N_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSD_DivDiv_RowSpace, NSType4NLSD_DivDiv_ColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEM_DivDiv_3D, LinCoeffs, NULL);
	  }
	  else // skew-symmetric
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkew3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkew3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew3D, LinCoeffs, NULL);
	  //
	  }

        }
        else
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDD3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDD3D, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);

          DiscreteFormNLSDFEM_DivDiv = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSD_DivDiv_N_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSD_DivDiv_RowSpace, NSType4NLSD_DivDiv_ColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEM_DivDiv_DD3D, LinCoeffs, NULL);
	  }
	  else // skew-symmetric
          {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDDSkew3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDSkew3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDDSkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinDDSkew3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDSkew3D, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDDSkew3D, LinCoeffs, NULL);
	  }
	      
        }
        break;
    } // endswitch
  } 
  else //Newton iteration
  {
    switch(NSTYPE)
    {
      case 1:
      case 2:
         OutPut("Wrong NSTYPE " << NSTYPE << " for Newton's method !!!");
         OutPut("Use NSTYPE 3 or 4 !!!");
         exit(4711);
      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
        // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinNewton3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindNewton3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyNewton3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLGalerkinNewton3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLUpwindNewton3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyNewton3D, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLGalerkinNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLUpwindNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyNewtonDD3D, LinCoeffs, NULL);
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinNewton3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMNewton3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindNewton3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyNewton3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLGalerkinNewton3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMNewton3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLUpwindNewton3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyNewton3D, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMNewtonDD3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLGalerkinNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMNewtonDD3D, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType3_4NLNewtonN_Matrices, NSType3_4NLNewtonN_Rhs, 
                  NSType3_4NLNewtonRowSpace, NSType3_4NLNewtonColumnSpace,
                  NSType3_4NLNewtonRhsSpace, NSType3_4NLUpwindNewtonDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyNewtonDD3D, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  }
}

void InitializeDiscreteFormsFriction(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormSDFEM,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormSmagorinsky,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLSDFEM,
  TDiscreteForm3D *&DiscreteFormNLUpwind,
  TDiscreteForm3D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm3D *&DiscreteFormGalerkinFriction,
  TDiscreteForm3D *&DiscreteFormNLGalerkinFriction,
  CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char Galerkin[] = "Galerkin";
  char all[] = "all";
  char Upwind[] = "Upwind";
  char Smagorinsky[] = "Smagorinsky";
  char nonlinear[] = "nonlinear";
  char Sdfem[] = "SDFEM";

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
  {
    DiscreteFormGalerkinFriction = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDD3DFriction, 
		  LinCoeffs, NULL);

     DiscreteFormNLGalerkinFriction = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDD3DFriction, LinCoeffs, NULL);

    switch(NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Galerkin3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind3D, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Smagorinsky3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkin3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind3D, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinsky3D, LinCoeffs, NULL);
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Galerkin3D, LinCoeffs, NULL);

        DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEM3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind3D, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Smagorinsky3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkin3D, LinCoeffs, NULL);
      
        DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEM3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind3D, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinsky3D, LinCoeffs, NULL);
        break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
        // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Galerkin3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Smagorinsky3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky3D, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Galerkin3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEM3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Smagorinsky3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEM3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky3D, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm3D(Sdfem, all,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDD3D, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Upwind, all, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm3D(Sdfem, nonlinear,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDD3D, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Upwind, nonlinear, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  } 
  else //Newton iteration
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}

void InitializeDiscreteForms(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormUpwindNC,
  TDiscreteForm3D *&DiscreteFormSmagorinsky,
  TDiscreteForm3D *&DiscreteFormClassicalLES,
  TDiscreteForm3D *&DiscreteFormGL00Convolution, 
  TDiscreteForm3D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm3D *&DiscreteFormVMS_Projection, 
  TDiscreteForm3D *&DiscreteFormVMS_SUPG,
  TDiscreteForm3D *&DiscreteFormLerayAlpha,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLUpwind, 
  TDiscreteForm3D *&DiscreteFormNLUpwindNC, 
  TDiscreteForm3D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm3D *&DiscreteFormNLClassicalLES, 
  TDiscreteForm3D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm3D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm3D *&DiscreteFormNLVMS_Projection, 
  TDiscreteForm3D *&DiscreteFormNLVMS_ProjectionExpl, 
  TDiscreteForm3D *&DiscreteFormNLVMSRFBExplRhs,
  TDiscreteForm3D *&DiscreteFormNLVMS_SUPG, 
  TDiscreteForm3D *&DiscreteFormNLLerayAlpha,
  TDiscreteForm3D *&DiscreteFormRHS,
  TDiscreteForm3D *&DiscreteFormRHSClassicalLES,
  TDiscreteForm3D *&DiscreteFormRHSLES,
  TDiscreteForm3D *&DiscreteFormRHSSUPG,
  TDiscreteForm3D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm3D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm3D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm3D *&DiscreteFormRHSAuxProblemU,
  TDiscreteForm3D *&DiscreteFormRHSNewton,
  TDiscreteForm3D *&DiscreteFormRHSNewtonNL,			
  TDiscreteForm3D *&DiscreteFormC, 
  TDiscreteForm3D *&DiscreteFormJ,
  CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char rhs[] = "rhs";
  char ClassicalLES[] = "ClassicalLES";
  char GaldiLayton[] = "Galdi-Layton";
  char auxprob[] = "aux prob";
  char Upwind[] = "Upwind";
  char Smagorinsky[] = "Smagorinsky";
  char nonlinear[] = "nonlinear";
  char all[] = "all";
  char Layton96[] = "Layton96";
  char Leray[] = "Leray";

  DiscreteFormGalerkin = NULL;
  DiscreteFormUpwind = NULL;
  DiscreteFormUpwindNC = NULL;
  DiscreteFormSmagorinsky = NULL;
  DiscreteFormClassicalLES = NULL;
  DiscreteFormGL00Convolution = NULL; 
  DiscreteFormGL00AuxProblem = NULL; 
  DiscreteFormVMS_Projection = NULL; 
  DiscreteFormVMS_SUPG = NULL;
  DiscreteFormLerayAlpha = NULL;
  DiscreteFormNLGalerkin = NULL;
  DiscreteFormNLUpwind = NULL; 
  DiscreteFormNLUpwindNC = NULL; 
  DiscreteFormNLSmagorinsky = NULL;
  DiscreteFormNLClassicalLES = NULL; 
  DiscreteFormNLGL00Convolution = NULL;
  DiscreteFormNLGL00AuxProblem = NULL;
  DiscreteFormNLVMS_Projection = NULL; 
  DiscreteFormNLVMS_ProjectionExpl = NULL; 
  DiscreteFormNLVMS_SUPG = NULL; 
  DiscreteFormNLVMSRFBExplRhs = NULL;
  DiscreteFormNLLerayAlpha = NULL;
  DiscreteFormRHS = NULL;
  DiscreteFormRHSClassicalLES = NULL;
  DiscreteFormRHSSUPG = NULL;
  DiscreteFormRHSLES = NULL;
  DiscreteFormMatrixGL00AuxProblem = NULL;
  DiscreteFormGL00AuxProblemRHS = NULL;
  DiscreteFormMatrixAuxProblemU = NULL;
  DiscreteFormRHSAuxProblemU = NULL;
  DiscreteFormRHSNewton = NULL;
  DiscreteFormRHSNewtonNL = NULL;			
  DiscreteFormC = NULL; 
  DiscreteFormJ = NULL;

  DiscreteFormRHS = new TDiscreteForm3D(GalerkinString, rhs,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS3D, LinCoeffs, NULL);

  DiscreteFormRHSAuxProblemU = new TDiscreteForm3D(GalerkinString, rhs,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHSAuxProblemU, LinCoeffs, NULL);

  DiscreteFormMatrixAuxProblemU = new TDiscreteForm3D(auxprob, all,
            MatrixAuxiliaryProblemN_Terms, MatrixAuxiliaryProblemDerivatives, 
            MatrixAuxiliaryProblemSpaceNumbers,
            MatrixAuxiliaryProblemN_Matrices, MatrixAuxiliaryProblemN_Rhs, 
            MatrixAuxiliaryProblemRowSpace, MatrixAuxiliaryProblemColumnSpace,
            MatrixAuxiliaryProblemRhsSpace, MatrixAuxiliaryProblem, LinCoeffs, NULL);

  DiscreteFormRHSClassicalLES = new TDiscreteForm3D(ClassicalLES, rhs,
            TimeNSRHSLESN_Terms, TimeNSRHSLESDerivatives, TimeNSRHSLESSpaceNumbers,
            TimeNSRHSLESN_Matrices, TimeNSRHSLESN_Rhs, 
            TimeNSRHSLESRowSpace, TimeNSRHSLESColumnSpace,
            TimeNSRHSLESRhsSpace, TimeNSRHSClassicalLES3D, LinCoeffs, NULL);

  DiscreteFormRHSLES = new TDiscreteForm3D(GaldiLayton, rhs,
            TimeNSRHSLESN_Terms, TimeNSRHSLESDerivatives, TimeNSRHSLESSpaceNumbers,
            TimeNSRHSLESN_Matrices, TimeNSRHSLESN_Rhs, 
            TimeNSRHSLESRowSpace, TimeNSRHSLESColumnSpace,
            TimeNSRHSLESRhsSpace, TimeNSRHSLESModel3D, LinCoeffs, NULL);

  DiscreteFormRHSSUPG = new TDiscreteForm3D(GalerkinString, rhs,
            TimeNSVMS_Rhs_SUPGN_Terms, TimeNSVMS_Rhs_SUPGDerivatives, TimeNSVMS_Rhs_SUPGSpaceNumbers,
            TimeNSVMS_Rhs_SUPGN_Matrices, TimeNSVMS_Rhs_SUPGN_Rhs, 
            TimeNSVMS_Rhs_SUPGRowSpace, TimeNSVMS_Rhs_SUPGColumnSpace,
            TimeNSVMS_Rhs_SUPGRhsSpace, TimeNSType4VMS_Rhs_SUPGDD3D, LinCoeffs, NULL);

  DiscreteFormGL00AuxProblemRHS = new TDiscreteForm3D(GaldiLayton, auxprob,
            TimeNSGL00AuxProblemRHSN_Terms, TimeNSGL00AuxProblemRHSDerivatives, 
            TimeNSGL00AuxProblemRHSSpaceNumbers,
            TimeNSGL00AuxProblemRHSN_Matrices, TimeNSGL00AuxProblemRHSN_Rhs, 
            TimeNSGL00AuxProblemRHSRowSpace, TimeNSGL00AuxProblemRHSColumnSpace,
            TimeNSGL00AuxProblemRHSRhsSpace, TimeNSGL00AuxProblemRHS3D,
            LinCoeffs, NULL);

 DiscreteFormRHSNewton = new TDiscreteForm3D(GalerkinString, rhs,
            TimeNSRHSNewtonN_Terms, TimeNSRHSNewtonDerivatives, TimeNSRHSNewtonSpaceNumbers,
            TimeNSRHSNewtonN_Matrices, TimeNSRHSNewtonN_Rhs, 
            TimeNSRHSNewtonRowSpace, TimeNSRHSNewtonColumnSpace,
            TimeNSRHSNewtonRhsSpace, TimeNSRHSNewton3D, LinCoeffs, NULL);

  DiscreteFormRHSNewtonNL = new TDiscreteForm3D(GalerkinString, rhs,
	    TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
	    TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHSNewtonNL3D, LinCoeffs, NULL);

  DiscreteFormC = new TDiscreteForm3D(GalerkinString, rhs,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSGalerkinC3D, LinCoeffs, NULL);

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind3D, LinCoeffs, NULL);
	DiscreteFormUpwindNC = DiscreteFormUpwind;

        DiscreteFormSmagorinsky = new TDiscreteForm3D(Layton96, all,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky3D, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormClassicalLES = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLGalerkin3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLUpwind3D, LinCoeffs, NULL);
	DiscreteFormNLUpwindNC = DiscreteFormNLUpwind;
        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLSmagorinsky3D, LinCoeffs, NULL);

        DiscreteFormGL00AuxProblem =  new TDiscreteForm3D(GaldiLayton, all,
                TimeNSType1N_Terms, TimeNSType1Derivatives, 
                TimeNSType1SpaceNumbers,
                TimeNSType1GL00AuxProblemN_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1GL00AuxProblemRowSpace, TimeNSType1GL00AuxProblemColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GL00AuxProblem3D, LinCoeffs, NULL);

        // same as Smagorinsky
        DiscreteFormNLClassicalLES =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;	

        // Discrete Forms for Rosenbrock Methods
	DiscreteFormJ = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1GalerkinJ3D, LinCoeffs, NULL);
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType2N_Terms, TimeNSType3Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                TimeNSType2N_Terms, TimeNSType3Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind3D, LinCoeffs, NULL);
	DiscreteFormUpwindNC = DiscreteFormUpwind;

        DiscreteFormSmagorinsky = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType2N_Terms, TimeNSType3Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky3D, LinCoeffs, NULL);


        // same as Smagorinsky
        DiscreteFormClassicalLES = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        DiscreteFormGL00AuxProblem =  new TDiscreteForm3D(GaldiLayton, all,
                TimeNSType2N_Terms, TimeNSType3Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2GL00AuxProblemN_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2GL00AuxProblemRowSpace, TimeNSType2GL00AuxProblemColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2GL00AuxProblem3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLGalerkin3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLUpwind3D, LinCoeffs, NULL);
	DiscreteFormNLUpwindNC = DiscreteFormNLUpwind;

        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLSmagorinsky3D, LinCoeffs, NULL);

        // same as Smagorinsky
        DiscreteFormNLClassicalLES = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
 
        // Discrete Forms for Rosenbrock Methods
        DiscreteFormJ = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType1GalerkinJ3D, LinCoeffs, NULL);
	break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind3D, LinCoeffs, NULL);
	  DiscreteFormUpwindNC = DiscreteFormUpwind;

          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky3D, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormClassicalLES = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;
          DiscreteFormNLVMSRFBExplRhs = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm3D(GaldiLayton, all,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblem3D, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D, LinCoeffs, NULL);
	  DiscreteFormNLUpwindNC = DiscreteFormNLUpwind;

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinsky3D, LinCoeffs, NULL);
          
          DiscreteFormNLClassicalLES = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

          // Discrete Forms for Rosenbrock Methods
          DiscreteFormJ = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinJ3D, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3UpwindDD3D, LinCoeffs, NULL);

	  // use nabla-nabla for nonconformings
          DiscreteFormUpwindNC = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind3D, LinCoeffs, NULL);
	  	  
         DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD3D, LinCoeffs, NULL);
 
          // same as Smagorinsky
          DiscreteFormClassicalLES = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm3D(GaldiLayton, all,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblemDD3D, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwindDD3D, LinCoeffs, NULL);

          DiscreteFormNLUpwindNC = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D, LinCoeffs, NULL);
	    
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLSmagorinskyN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLSmagorinskyRowSpace, TimeNSType3_4NLSmagorinskyColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);
          
          // same as Smagorinsky
          DiscreteFormNLClassicalLES = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

          // Discrete Forms for Rosenbrock Methods
          DiscreteFormJ = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinJ3D, LinCoeffs, NULL);
       }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind3D, LinCoeffs, NULL);

	  DiscreteFormUpwindNC = DiscreteFormUpwind;

         DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky3D, LinCoeffs, NULL);


          // same as Smagorinsky
          DiscreteFormClassicalLES = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm3D(GaldiLayton, all,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblem3D, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D, LinCoeffs, NULL);
	  DiscreteFormNLUpwindNC = DiscreteFormNLUpwind;

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinsky3D, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormNLClassicalLES = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD3D, LinCoeffs, NULL);
	
	  // use nabla-nabla for nonconformings
	  DiscreteFormUpwindNC = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind3D, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD3D, LinCoeffs, NULL);

          DiscreteFormVMS_SUPG = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4VMS_SUPGN_Terms, TimeNSType4VMS_SUPGDerivatives, 
                  TimeNSType4VMS_SUPGSpaceNumbers,
                  TimeNSType4VMS_SUPGN_Matrices, TimeNSType4VMS_SUPGN_Rhs, 
                  TimeNSType4VMS_SUPGRowSpace, TimeNSType4VMS_SUPGColumnSpace,
                  TimeNSType4VMS_SUPGRhsSpace, TimeNSType4VMS_SUPGDD3D, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormClassicalLES = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;
        
          DiscreteFormGL00AuxProblem =  new TDiscreteForm3D(GaldiLayton, all,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblemDD3D, LinCoeffs, NULL);

          DiscreteFormVMS_Projection = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4VMS_ProjectionN_Terms, TimeNSType4VMS_ProjectionDerivatives, 
                  TimeNSType4VMS_ProjectionSpaceNumbers,
                  TimeNSType4VMS_ProjectionN_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4VMS_ProjectionRowSpace, TimeNSType4VMS_ProjectionColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4VMS_ProjectionDD3D, LinCoeffs, NULL);

	  if (TDatabase::ParamDB->DISCTYPE==VMS_PROJECTION_SD)
	      DiscreteFormVMS_Projection = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4VMS_ProjectionN_Terms, TimeNSType4VMS_ProjectionDerivatives, 
                  TimeNSType4VMS_ProjectionSpaceNumbers,
                  TimeNSType4VMS_ProjectionN_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4VMS_ProjectionRowSpace, TimeNSType4VMS_ProjectionColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4VMS_ProjectionStreamlineDD3D, LinCoeffs, NULL);
        
          DiscreteFormLerayAlpha=  new TDiscreteForm3D(Leray, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, 
                  TimeNSType4SpaceNumbers,
                  TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4LerayAlphaDD3D, LinCoeffs, NULL);

	      
          DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwindDD3D, LinCoeffs, NULL);

          DiscreteFormNLUpwindNC = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLSmagorinskyN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLSmagorinskyRowSpace, TimeNSType3_4NLSmagorinskyColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);
   
           DiscreteFormNLVMS_SUPG = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4NLVMS_SUPGN_Terms, TimeNSType4NLVMS_SUPGDerivatives, 
                  TimeNSType4NLVMS_SUPGSpaceNumbers,
                  TimeNSType4NLVMS_SUPGN_Matrices, TimeNSType4NLVMS_SUPGN_Rhs, 
                  TimeNSType4NLVMS_SUPGRowSpace, TimeNSType4NLVMS_SUPGColumnSpace,
                  TimeNSType4NLVMS_SUPGRhsSpace, TimeNSType4NLVMS_SUPGDD3D, LinCoeffs, NULL);

         // same as Smagorinsky
          DiscreteFormNLClassicalLES = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
          if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE==17)
          {
           DiscreteFormNLVMS_Projection = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NL_Adap_VMS_ProjectionN_Terms, TimeNSType3_4NL_Adap_VMS_ProjectionDerivatives, 
                  TimeNSType3_4NL_Adap_VMS_ProjectionSpaceNumbers,
                  TimeNSType3_4NL_Adap_VMS_ProjectionN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NL_Adap_VMS_ProjectionRowSpace, TimeNSType3_4NL_Adap_VMS_ProjectionColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NL_Adap_VMS_ProjectionDD3D, LinCoeffs, NULL);
           }
           else
           {
             DiscreteFormNLVMS_Projection = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLVMS_ProjectionN_Terms, TimeNSType3_4NLVMS_ProjectionDerivatives, 
                  TimeNSType3_4NLVMS_ProjectionSpaceNumbers,
                  TimeNSType3_4NLVMS_ProjectionN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLVMS_ProjectionRowSpace, TimeNSType3_4NLVMS_ProjectionColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLVMS_ProjectionDD3D, LinCoeffs, NULL);
	   }
	  if (TDatabase::ParamDB->DISCTYPE==VMS_PROJECTION_SD)
	      DiscreteFormNLVMS_Projection = new TDiscreteForm3D(Smagorinsky, nonlinear,
		  TimeNSType3_4NLVMS_ProjectionN_Terms, TimeNSType3_4NLVMS_ProjectionDerivatives, 
                  TimeNSType3_4NLVMS_ProjectionSpaceNumbers,
                  TimeNSType3_4NLVMS_ProjectionN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLVMS_ProjectionRowSpace, TimeNSType3_4NLVMS_ProjectionColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLVMS_ProjectionStreamlineDD3D, LinCoeffs, NULL);

          DiscreteFormNLVMS_ProjectionExpl = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLVMS_ProjectionExplN_Terms, TimeNSType3_4NLVMS_ProjectionExplDerivatives, 
                  TimeNSType3_4NLVMS_ProjectionExplSpaceNumbers,
                  TimeNSType3_4NLVMS_ProjectionExplN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLVMS_ProjectionExplRowSpace, TimeNSType3_4NLVMS_ProjectionExplColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4VMS_ProjectionExpl3D, LinCoeffs, NULL);

          DiscreteFormNLVMSRFBExplRhs = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLSmagorinskyN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLSmagorinskyRowSpace, TimeNSType3_4NLSmagorinskyColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLDivDivDD3D, LinCoeffs, NULL);

          DiscreteFormNLLerayAlpha = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
         
       // Discrete Forms for Rosenbrock Methods
        DiscreteFormJ = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                TimeNSType4RhsSpace, TimeNSType3GalerkinJ3D, LinCoeffs, NULL);
	  
        break;
      case 14:
	      DiscreteFormVMS_SUPG = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType14VMS_SUPGN_Terms, TimeNSType14VMS_SUPGDerivatives, 
                  TimeNSType14VMS_SUPGSpaceNumbers,
                  TimeNSType14VMS_SUPGN_Matrices, TimeNSType14VMS_SUPGN_Rhs, 
                  TimeNSType14VMS_SUPGRowSpace, TimeNSType14VMS_SUPGColumnSpace,
                  TimeNSType14VMS_SUPGRhsSpace, TimeNSType14VMS_SUPGDD3D, LinCoeffs, NULL);
        DiscreteFormNLVMS_SUPG = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType14NLVMS_SUPGN_Terms, TimeNSType14NLVMS_SUPGDerivatives, 
                  TimeNSType14NLVMS_SUPGSpaceNumbers,
                  TimeNSType14NLVMS_SUPGN_Matrices, TimeNSType14NLVMS_SUPGN_Rhs, 
                  TimeNSType14NLVMS_SUPGRowSpace, TimeNSType14NLVMS_SUPGColumnSpace,
                  TimeNSType14NLVMS_SUPGRhsSpace, TimeNSType14NLVMS_SUPGDD3D, LinCoeffs, NULL);
	  break;
    } // endswitch
  else // Newton's method
  {

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 3:
	if(TDatabase::ParamDB->LAPLACETYPE == 0)
	  {
	    // (grad, grad)
	    DiscreteFormGalerkin = 
	      new TDiscreteForm3D(GalerkinString, all,
				  TimeNSType3NewtonN_Terms, TimeNSType3NewtonDerivatives, TimeNSType3NewtonSpaceNumbers,
				  TimeNSType3NewtonN_Matrices, TimeNSType3NewtonN_Rhs, 
				  TimeNSType3NewtonRowSpace, TimeNSType3NewtonColumnSpace,
				  TimeNSType3NewtonRhsSpace, TimeNSType3GalerkinNewton3D, LinCoeffs, NULL);
	    
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
				  TimeNSType3NewtonN_Terms, TimeNSType3NewtonDerivatives, TimeNSType3NewtonSpaceNumbers,
				  TimeNSType3NewtonN_Matrices, TimeNSType3NewtonN_Rhs, 
				  TimeNSType3NewtonRowSpace, TimeNSType3NewtonColumnSpace,
				  TimeNSType3NewtonRhsSpace, TimeNSType3UpwindNewton3D, LinCoeffs, NULL);
          
          DiscreteFormUpwindNC = DiscreteFormUpwind;
          DiscreteFormNLGalerkin = 
	    new TDiscreteForm3D(GalerkinString, nonlinear,
				TimeNSType3_4NLNewtonN_Terms, TimeNSType3_4NLNewtonDerivatives, TimeNSType3_4NLNewtonSpaceNumbers,
				TimeNSType3_4NLNewtonN_Matrices, TimeNSType3_4NLNewtonN_Rhs, 
				TimeNSType3_4NLNewtonRowSpace, TimeNSType3_4NLNewtonColumnSpace,
				TimeNSType3_4NLNewtonRhsSpace, TimeNSType3_4NLGalerkinNewton3D, LinCoeffs, NULL);

          DiscreteFormNLUpwind = 
	    new TDiscreteForm3D(Upwind, nonlinear,
				TimeNSType3_4NLNewtonN_Terms, TimeNSType3_4NLNewtonDerivatives, TimeNSType3_4NLNewtonSpaceNumbers,
				TimeNSType3_4NLNewtonN_Matrices, TimeNSType3_4NLNewtonN_Rhs, 
				TimeNSType3_4NLNewtonRowSpace, TimeNSType3_4NLNewtonColumnSpace,
				TimeNSType3_4NLNewtonRhsSpace, TimeNSType3_4NLUpwindNewton3D, LinCoeffs, NULL);
          DiscreteFormNLUpwindNC = DiscreteFormNLUpwind;
	  }
	else
	  {
	    OutPut("Newton method in this case not implemented yet !!!!!!!!" << endl);
	    exit(4711);
	  }
	break;
      case 4:
	DiscreteFormSmagorinsky = 
	  new TDiscreteForm3D(Smagorinsky, all,
			      TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
			      TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
			      TimeNSType4RowSpace, TimeNSType4ColumnSpace,
			      TimeNSType4RhsSpace, TimeNSType4SmagorinskyNewtonDD3D, LinCoeffs, NULL);
	
	DiscreteFormVMS_Projection = 
	  new TDiscreteForm3D(Smagorinsky, all,
			      TimeNSType4VMS_ProjectionN_Terms, TimeNSType4VMS_ProjectionDerivatives, 
			      TimeNSType4VMS_ProjectionSpaceNumbers,
			      TimeNSType4VMS_ProjectionN_Matrices, TimeNSType4N_Rhs, 
			      TimeNSType4VMS_ProjectionRowSpace, TimeNSType4VMS_ProjectionColumnSpace,
			      TimeNSType4RhsSpace, TimeNSType4VMS_ProjectionNewtonDD3D, LinCoeffs, NULL);
	DiscreteFormNLSmagorinsky = 
	  new TDiscreteForm3D(Smagorinsky, nonlinear,
			      TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
			      TimeNSType3_4NLSmagorinskyN_Matrices, TimeNSType3_4NLN_Rhs, 
			      TimeNSType3_4NLSmagorinskyRowSpace, TimeNSType3_4NLSmagorinskyColumnSpace,
			      TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinskyNewtonDD3D, LinCoeffs, NULL);
	
	DiscreteFormNLVMS_Projection = 
	  new TDiscreteForm3D(Smagorinsky, nonlinear,
			      TimeNSType3_4NLVMS_ProjectionN_Terms, TimeNSType3_4NLVMS_ProjectionDerivatives, 
			      TimeNSType3_4NLVMS_ProjectionSpaceNumbers,
			      TimeNSType3_4NLVMS_ProjectionN_Matrices, TimeNSType3_4NLN_Rhs, 
			      TimeNSType3_4NLVMS_ProjectionRowSpace, TimeNSType3_4NLVMS_ProjectionColumnSpace,
			      TimeNSType3_4NLRhsSpace, TimeNSType3_4NLVMS_ProjectionNewtonDD3D, LinCoeffs, NULL);   
	break;
    default: 
      OutPut("Newton method in this case not implemented yet !!!!!!!!" << endl);
      exit(4711);
      
    } 
  }
}

// set all needed discrete forms
void InitializeDiscreteFormsForFreeSurface(  
  TDiscreteForm3D* &DiscreteFormGalerkin,
  TDiscreteForm3D* &DiscreteFormUpwind,
  TDiscreteForm3D* &DiscreteFormNLGalerkin,
  TDiscreteForm3D* &DiscreteFormNLUpwind, 
  CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char Galerkin[] = "Galerkin";
  char rhs[] = "rhs";
  char Upwind[] = "Upwind";
  char nonlinear[] = "nonlinear";
  char all[] = "all";

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                TimeNSType1N_Terms, TimeNSType1Derivatives,
                TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                TimeNSType1N_Terms, TimeNSType1Derivatives,
                TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives,
                TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLGalerkin3D,
                LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives,
                TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLUpwind3D,
                LinCoeffs, NULL);
      break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                TimeNSType2N_Terms, TimeNSType3Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                TimeNSType2N_Terms, TimeNSType3Derivatives,
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives,
                TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLGalerkin3D,
                LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives,
                TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLUpwind3D,
                LinCoeffs, NULL);
      break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives,
                  TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives,
                  TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkin3D,
                  LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D,
                  LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives,
                  TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinDD3D,
                  LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives,
                  TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3UpwindDD3D,
                  LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D,
                  LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwindDD3D,
                  LinCoeffs, NULL);
        }
      break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives,
                  TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin3D,
                  LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives,
                  TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind3D,
                  LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkin3D,
                  LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D,
                  LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(Galerkin, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives,
                  TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD3D,
                  LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives,
                  TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD3D,
                  LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm3D(Galerkin, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D,
                  LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives,
                  TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwindDD3D,
                  LinCoeffs, NULL);
        }
      break;
    } // endswitch
  else
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}

void InitializeDiscreteFormsVMS(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormSmagorinsky,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLUpwind, 
  TDiscreteForm3D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm3D *&DiscreteFormRHS,
  TDiscreteForm3D *&DiscreteForm_ho_RHS,
  CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char rhs[] = "rhs";
  char Upwind[] = "Upwind";
  char Smagorinsky[] = "Smagorinsky";
  char nonlinear[] = "nonlinear";
  char all[] = "all";
  char Layton96[] = "Layton96";
  
  DiscreteFormRHS = new TDiscreteForm3D(
    GalerkinString, rhs,
    TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
    TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
    TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
    TimeNSRHSRhsSpace, TimeNSRHS3D, LinCoeffs, NULL);
  
  DiscreteForm_ho_RHS = new TDiscreteForm3D(
    GalerkinString, rhs,
    TimeNS_ho_RHSN_Terms, TimeNS_ho_RHSDerivatives, TimeNS_ho_RHSSpaceNumbers,
    TimeNS_ho_RHSN_Matrices, TimeNS_ho_RHSN_Rhs, 
    TimeNS_ho_RHSRowSpace, TimeNS_ho_RHSColumnSpace,
    TimeNS_ho_RHSRhsSpace, TimeNS_VMS_SmallRhs3D, LinCoeffs, NULL);
  
   //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind3D, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm3D(Layton96, all,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLGalerkin3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLUpwind3D, LinCoeffs, NULL);
     
        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLSmagorinsky3D, LinCoeffs, NULL);
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType2N_Terms, TimeNSType3Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin3D, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                TimeNSType2N_Terms, TimeNSType3Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind3D, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm3D(GalerkinString, all,
                TimeNSType2N_Terms, TimeNSType3Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky3D, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLGalerkin3D, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLUpwind3D, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm3D(GalerkinString, nonlinear,
                TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
                TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, 
                TimeNSType1_2NLRowSpace, TimeNSType1_2NLColumnSpace,
                TimeNSType1_2NLRhsSpace, TimeNSType1_2NLSmagorinsky3D, LinCoeffs, NULL);

       break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky3D, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinsky3D, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3UpwindDD3D, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD3D, LinCoeffs, NULL);
 
         DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLSmagorinskyN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLSmagorinskyRowSpace, TimeNSType3_4NLSmagorinskyColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);
          
       }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind3D, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky3D, LinCoeffs, NULL);

         DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkin3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwind3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinsky3D, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm3D(Upwind, all, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD3D, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm3D(Smagorinsky, all,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD3D, LinCoeffs, NULL);
        
          DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm3D(Upwind, nonlinear, 
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLRowSpace, TimeNSType3_4NLColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLUpwindDD3D, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm3D(Smagorinsky, nonlinear,
                  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
                  TimeNSType3_4NLSmagorinskyN_Matrices, TimeNSType3_4NLN_Rhs, 
                  TimeNSType3_4NLSmagorinskyRowSpace, TimeNSType3_4NLSmagorinskyColumnSpace,
                  TimeNSType3_4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD3D, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  else
  {
     OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
 
  }
}

void InitializeDiscreteForms
     ( TDiscreteForm3D *&DiscreteFormGalerkin,
       TDiscreteForm3D *&DiscreteFormNLGalerkin,
       CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "galerkin";
  char all[]      = "all";
  
  DiscreteFormGalerkin = NULL;
  DiscreteFormNLGalerkin = NULL;
  
  switch (NSTYPE)
  {
    case 1:
      DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                NSType1N_Terms, NSType1Derivatives, NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Galerkin3D, LinCoeffs, NULL);
		
      DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, all,
                NSType1NLN_Terms, NSType1NLDerivatives, NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkin3D, LinCoeffs, NULL);
      
      break;
      
    case 2:
      DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,
                NSType2N_Terms, NSType2Derivatives, NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Galerkin3D, LinCoeffs, NULL);
		
      DiscreteFormNLGalerkin = new TDiscreteForm3D(GalerkinString, all,
                NSType2NLN_Terms, NSType2NLDerivatives, NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkin3D, LinCoeffs, NULL);
      
      break;
      
    case 3:
      
      break;
    
    case 4:
      
      break;
      
    default:
      OutPut("Unknown NSTYPE " << NSTYPE << endl);
  }
}

void InitializeDiscreteForms
     ( TDiscreteForm3D *&DiscreteFormGalerkin,
       TDiscreteForm3D *&DiscreteFormNLGalerkin,
       TDiscreteForm3D *&DiscreteFormRHS,
       CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "galerkin";
  char Description[] = "desc";
  
  DiscreteFormGalerkin = NULL;
  DiscreteFormNLGalerkin = NULL;
  DiscreteFormRHS = NULL;
  
  switch (NSTYPE)
  {      
    case 2:
      DiscreteFormGalerkin = new TDiscreteForm3D (GalerkinString, Description,
		  TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
		  TimeNSType2N_Matrices, TimeNSType2N_Rhs, TimeNSType2RowSpace,
		  TimeNSType2ColumnSpace, TimeNSType2RhsSpace, TimeNSType2Galerkin3D,
		  LinCoeffs, NULL);
		  
      DiscreteFormNLGalerkin = new TDiscreteForm3D (GalerkinString, Description,
	TimeNSType1_2NLN_Terms, TimeNSType1_2NLDerivatives, TimeNSType1_2NLSpaceNumbers,
	TimeNSType1_2NLN_Matrices, TimeNSType1_2NLN_Rhs, TimeNSType1_2NLRowSpace,
	TimeNSType1_2NLColumnSpace, TimeNSType1_2NLRhsSpace, TimeNSType1_2NLGalerkin3D,
	LinCoeffs, NULL);
	
      DiscreteFormRHS = new TDiscreteForm3D (GalerkinString, Description,
		  TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
		  TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, TimeNSRHSRowSpace,
		  TimeNSRHSColumnSpace, TimeNSRHSRhsSpace, TimeNSRHS3D,
		  LinCoeffs, NULL);
      
      break;
      
    case 3:
      break;
      
    case 4:
      if ( TDatabase::ParamDB->LAPLACETYPE == 0) // (grad u, grad v)
      {
	DiscreteFormGalerkin = new TDiscreteForm3D (GalerkinString, Description,
		  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
		  TimeNSType4N_Matrices, TimeNSType4N_Rhs, TimeNSType4RowSpace,
		  TimeNSType4ColumnSpace, TimeNSType4RhsSpace, TimeNSType4Galerkin3D,
		  LinCoeffs, NULL);
		  
	DiscreteFormNLGalerkin = new TDiscreteForm3D (GalerkinString, Description,
	  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
	  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, TimeNSType3_4NLRowSpace,
	  TimeNSType3_4NLColumnSpace, TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkin3D,
	  LinCoeffs, NULL);
	  
	DiscreteFormRHS = new TDiscreteForm3D (GalerkinString, Description,
		  TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
		  TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, TimeNSRHSRowSpace,
		  TimeNSRHSColumnSpace, TimeNSRHSRhsSpace, TimeNSRHS3D,
		  LinCoeffs, NULL);
      }
      else // ( D(u), D(v) )
      {
	DiscreteFormGalerkin = new TDiscreteForm3D (GalerkinString, Description,
		  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
		  TimeNSType4N_Matrices, TimeNSType4N_Rhs, TimeNSType4RowSpace,
		  TimeNSType4ColumnSpace, TimeNSType4RhsSpace, TimeNSType4GalerkinDD3D,
		  LinCoeffs, NULL);
		  
	DiscreteFormNLGalerkin = new TDiscreteForm3D (GalerkinString, Description,
	  TimeNSType3_4NLN_Terms, TimeNSType3_4NLDerivatives, TimeNSType3_4NLSpaceNumbers,
	  TimeNSType3_4NLN_Matrices, TimeNSType3_4NLN_Rhs, TimeNSType3_4NLRowSpace,
	  TimeNSType3_4NLColumnSpace, TimeNSType3_4NLRhsSpace, TimeNSType3_4NLGalerkinDD3D,
	  LinCoeffs, NULL);
	  
	DiscreteFormRHS = new TDiscreteForm3D (GalerkinString, Description,
		  TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
		  TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, TimeNSRHSRowSpace,
		  TimeNSRHSColumnSpace, TimeNSRHSRhsSpace, TimeNSRHS3D,
		  LinCoeffs, NULL);
      }
      break;
    
    default:
      OutPut("NSTYPE " << NSTYPE << " not yet supported" << endl);
  }
}

void InitializeDiscreteForms(TDiscreteForm3D *&DiscreteForm, CoeffFct3D *LinCoeff)
{
  
  char CdString[] = "Conv-Diff";
  char GalString[] = "Galerkin";
  
  
  DiscreteForm =  new TDiscreteForm3D(CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
                                      N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
                                      BilinearAssemble, LinCoeff, NULL);
  
 
}

void InitializeDiscreteFormsScalar(TDiscreteForm3D *&DiscreteFormMRhs_Galerkin, TDiscreteForm3D *&DiscreteFormARhs_Galerkin, 
                                   TDiscreteForm3D *&DiscreteFormRhs, CoeffFct3D *LinCoeff)
{
  
  char MMString[] = "Mass matrix";
  
   DiscreteFormMRhs_Galerkin = new TDiscreteForm3D(MMString, MMString, N_Terms_MatrixMRhs,
                                 Derivatives_MatrixMRhs, SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
                                 RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
                                 MatrixMRhsAssemble, LinCoeff, NULL);

   DiscreteFormARhs_Galerkin = new TDiscreteForm3D(MMString, MMString, N_Terms_MatrixARhs,
                                 Derivatives_MatrixARhs, SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs,
                                 N_Rhs_MatrixARhs, RowSpace_MatrixARhs, ColumnSpace_MatrixARhs,
                                 RhsSpace_MatrixARhs, MatrixARhsAssemble, LinCoeff, NULL);

   DiscreteFormRhs = new TDiscreteForm3D(MMString, MMString, N_Terms_Rhs, Derivatives_Rhs,
                          SpacesNumbers_Rhs, N_Matrices_Rhs, N_Rhs_Rhs,
                          RowSpace_Rhs, ColumnSpace_Rhs, RhsSpace_Rhs,
                          RhsAssemble, LinCoeff, NULL);
}

void InitializeDiscreteFormGrid(TDiscreteForm3D *&DiscreteFormGrid,
				CoeffFct3D *GridCoeffs)
{
  char GridString_3D[] = "Grid";
  char allString_3D[] = "all";
  
  DiscreteFormGrid = new TDiscreteForm3D(GridString_3D, allString_3D,
            GridN_Terms_3D, GridDerivatives_3D, GridSpaceNumbers_3D,
            GridN_Matrices_3D, GridN_Rhs_3D,
            GridRowSpace_3D, GridColumnSpace_3D,
            GridRhsSpace_3D, GridAssemble4,
             GridCoeffs, NULL);    
//  cout << " InitializeDiscreteFormGrid not yet implemented" <<endl;
//  exit(0);
    
}

/*
void InitializeDiscreteFormsOS_ST(
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormStiffRhsOS_ST,
  CoeffFct3D *LinCoeffs, int NSTYPE)
{
  char AllString[] = "mass, stiffness, rhs";
  char StiffRhsString[] = "stiffness, rhs";
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      DiscreteFormGalerkin = new TDiscreteForm3D(
        AllString, AllString, TimeOSType1N_Terms,
        TimeOSType1Derivatives, TimeOSType1SpaceNumbers,
        TimeOSType1N_Matrices, TimeOSType1N_Rhs, 
        TimeOSType1RowSpace, TimeOSType1ColumnSpace,
        TimeOSType1RhsSpace, TimeOSType1Galerkin3D,
        LinCoeffs, NULL);
      DiscreteFormStiffRhsOS_ST = new TDiscreteForm3D(
        StiffRhsString, StiffRhsString, 
        TimeOSType1_2N_TermsStiffRhs, TimeOSType1_2DerivativesStiffRhs,
        TimeOSType1_2SpaceNumbersStiffRhs, TimeOSType1_2N_MatricesStiffRhs,
        TimeOSType1_2N_RhsStiffRhs, TimeOSType1_2RowSpaceStiffRhs,
        TimeOSType1_2ColumnSpaceStiffRhs, TimeOSType1_2RhsSpaceStiffRhs,
        TimeOSType1_2StiffRhsAssemble3D, LinCoeffs, NULL);
      break;
    case 2:
      DiscreteFormGalerkin = new TDiscreteForm3D(
        AllString, AllString, TimeOSType2N_Terms,
        TimeOSType2Derivatives, TimeOSType2SpaceNumbers,
        TimeOSType2N_Matrices, TimeOSType2N_Rhs, 
        TimeOSType2RowSpace, TimeOSType2ColumnSpace,
        TimeOSType2RhsSpace, TimeOSType2Galerkin3D,
        LinCoeffs, NULL);
       
       DiscreteFormStiffRhsOS_ST = new TDiscreteForm3D(
        StiffRhsString, StiffRhsString, 
        TimeOSType1_2N_TermsStiffRhs, TimeOSType1_2DerivativesStiffRhs,
        TimeOSType1_2SpaceNumbersStiffRhs, TimeOSType1_2N_MatricesStiffRhs,
        TimeOSType1_2N_RhsStiffRhs, TimeOSType1_2RowSpaceStiffRhs,
        TimeOSType1_2ColumnSpaceStiffRhs, TimeOSType1_2RhsSpaceStiffRhs,
        TimeOSType1_2StiffRhsAssemble3D, LinCoeffs, NULL);
       
       break;
    case 4:
      OutPut("not implemented ");
      exit(-1);
    break;
  }
}
*/
