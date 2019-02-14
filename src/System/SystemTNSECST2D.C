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
   
 
/** ************************************************************************ 
* @brief     source file for TSystemTNSECST2D
* @author    Jagannath Venkatesan, 
* @date      19.08.16
* @History    
*            10.09.16 : DEVSS working (Sequential only)
*            12.09.16 : LPS working (Sequential only)
*            15.09.16 - Implemented SMPI for DEVSS and LPS
*            Slip with friction implemented
* 
 ************************************************************************  */
#ifdef __2D__
#include <Database.h>
#include <SystemNSECST2D.h>
#include <SystemTNSECST2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
#include <MainUtilities.h>

#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif

TSystemTNSECST2D::TSystemTNSECST2D(int N_FESpaces, TFESpace2D **FE_Spaces, int N_FEFunctions, 
				 TFEFunction2D **FeFunctions, int disctype_NSECST,  int solver) : 
				 TSystemNSECST2D(N_FESpaces, FE_Spaces, N_FEFunctions, FeFunctions, disctype_NSECST,  solver)
{

   B = new double[N_TotalDof];
   defect = new double[N_TotalDof];
   NLRhs = new double[3*N_S];
   gamma =0.; 
   
   // allocate the mass matrices in addition
   SqmatrixMU11 = new TSquareMatrix2D(sqstructureA);
   SqmatrixMU12 = new TSquareMatrix2D(sqstructureA);
   SqmatrixMU21 = new TSquareMatrix2D(sqstructureA);
   SqmatrixMU22 = new TSquareMatrix2D(sqstructureA);
   
   SqmatrixMS11 = new TSquareMatrix2D(sqstructureG);  
   SqmatrixMS12 = new TSquareMatrix2D(sqstructureG); 
   SqmatrixMS21 = new TSquareMatrix2D(sqstructureG); 
   SqmatrixMS22 = new TSquareMatrix2D(sqstructureG); 
   SqmatrixMS23 = new TSquareMatrix2D(sqstructureG);
   SqmatrixMS32 = new TSquareMatrix2D(sqstructureG); 
   SqmatrixMS33 = new TSquareMatrix2D(sqstructureG);

   SystMatAssembled  = FALSE;
   
   olderror_l_2_l_2u = 0;
   olderror_l_2_l_2s = 0;
   olderror_l_2_h_1u = 0;
   olderror_l_2_h_1s = 0;
   olderror_l_2_h_1p = 0;
   olderror_l_2_h_1p = 0;

   
} // TSystemTNSECST2D::TSystemTNSECST2D

TSystemTNSECST2D::~TSystemTNSECST2D()
{
  delete SqmatrixMU11; delete SqmatrixMU12;  delete SqmatrixMU21; delete SqmatrixMU22; 
  delete SqmatrixMS11; delete SqmatrixMS12;  delete SqmatrixMS21; delete SqmatrixMS22; 
  delete SqmatrixMS23; delete SqmatrixMS32;  delete SqmatrixMS33; 
}
  
 
void TSystemTNSECST2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D **BoundCond, BoundValueFunct2D **BoundValue,  
			   TAuxParam2D *TNSECSTaux, TAuxParam2D *TNSEauxerror, TAuxParam2D *TCSTauxerror)
{ int i;
  TDiscreteForm2D *DiscreteFormGalerkinSUPG, *DiscreteFormLPS, *DiscreteFormNLGalerkinSUPG, *DiscreteFormNLLPS, *DiscreteFormRHSGalerkinSUPG, *DiscreteFormRHSLPS;
  
  for (i=0;i<8;i++) 
  { 
  // save the boundary condition
  BoundaryConditions[i] = BoundCond[i];
  // save the boundary values  
  BoundaryValues[i] = BoundValue[i];
  }
  
  // save the bilinear coefficients
  LinCoeffs[0] = lincoeffs;
  
  // aux for assembling and error calculations
  aux = TNSECSTaux;
  auxerror_NSE = TNSEauxerror;
  auxerror_CST = TCSTauxerror;
  
  // set the Discreteforms
   InitializeDiscreteForms_TNSECST(DiscreteFormGalerkinSUPG, DiscreteFormLPS, 
				   DiscreteFormNLGalerkinSUPG, DiscreteFormNLLPS,
				   DiscreteFormRHSGalerkinSUPG, DiscreteFormRHSLPS,
				   LinCoeffs[0]);
  
     if (Disctype_NSECST == 1)
     {  
       DiscreteFormARhs = DiscreteFormGalerkinSUPG;
       DiscreteFormNL   = DiscreteFormNLGalerkinSUPG;
       DiscreteFormRhs  = DiscreteFormRHSGalerkinSUPG;
     }
    else if (Disctype_NSECST == 2)
    {
      DiscreteFormARhs = DiscreteFormLPS;
      DiscreteFormNL   = DiscreteFormNLLPS;
      DiscreteFormRhs  = DiscreteFormRHSLPS;
    }
    else 
    {
      cout<<"Invalid Disctype_NSECST !!!\n";
      exit(0);
    }
  
     // LPS assembling terms
         if (Disctype_NSECST == 2)
      {
       clock_t Time_Assemble = clock();
       
	double delta0 = TDatabase::ParamDB->DELTA0;
	double delta1 = TDatabase::ParamDB->DELTA1 * (1.0 - TDatabase::ParamDB->P3);
	double delta2 = TDatabase::ParamDB->P4;

        AddDeformationTensorTerm(SqmatrixA11_LPS, SqmatrixA12_LPS, SqmatrixA21_LPS, SqmatrixA22_LPS, delta1, 1.0, 1);
        AddTauTerm(SqmatrixG11_LPS,SqmatrixG12_LPS,SqmatrixG21_LPS,SqmatrixG22_LPS,SqmatrixG23_LPS,SqmatrixG32_LPS,SqmatrixG33_LPS,delta0,delta2,1.0,1); 
  
	 cout << "LPS Assemble time : "<<double( clock() - Time_Assemble ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
       }
       
       // Matrix initialization due to SMPI (Otherwise not required)
   int N_SquareMatrices, N_RectMatrices, N_Rhs;

   SQMATRICES[0] = SqmatrixMU11;
   SQMATRICES[1] = SqmatrixMU12;
   SQMATRICES[2] = SqmatrixMU21;
   SQMATRICES[3] = SqmatrixMU22;
   
   SQMATRICES[4] = SqmatrixMS11;
   SQMATRICES[5] = SqmatrixMS12;
   SQMATRICES[6] = SqmatrixMS21;
   SQMATRICES[7] = SqmatrixMS22;
   SQMATRICES[8] = SqmatrixMS23;
   SQMATRICES[9] = SqmatrixMS32;
   SQMATRICES[10] = SqmatrixMS33;
  
   N_SquareMatrices = 11;
   
   SQMATRICES[0]->Reset();
   SQMATRICES[1]->Reset();
   SQMATRICES[2]->Reset();
   SQMATRICES[3]->Reset();
   SQMATRICES[4]->Reset();
   SQMATRICES[5]->Reset();
   SQMATRICES[6]->Reset();
   SQMATRICES[7]->Reset();
   SQMATRICES[8]->Reset();
   SQMATRICES[9]->Reset();
   SQMATRICES[10]->Reset();
   
   if (Disctype_NSECST == 1)
   {
   SQMATRICES[11] = SqmatrixH11;
   SQMATRICES[12] = SqmatrixH22;
   SQMATRICES[13] = SqmatrixH33;
   N_SquareMatrices += 3;
   
   SQMATRICES[11]->Reset();
   SQMATRICES[12]->Reset();
   SQMATRICES[13]->Reset();
   }
   
   
   MATRICES[0] = MatrixB1;
   MATRICES[1] = MatrixB2;
   
   MATRICES[2] = MatrixB1T;
   MATRICES[3] = MatrixB2T;
   
   MATRICES[4] = MatrixC11;
   MATRICES[5] = MatrixC12;
   MATRICES[6] = MatrixC22;
   MATRICES[7] = MatrixC23;
   
   MATRICES[8] = MatrixD11;
   MATRICES[9] = MatrixD12;
   MATRICES[10] = MatrixD21;
   MATRICES[11] = MatrixD22;
   MATRICES[12] = MatrixD31;
   MATRICES[13] = MatrixD32;
   
   N_RectMatrices = 14;
   
   MATRICES[0]->Reset();
   MATRICES[1]->Reset();
   MATRICES[2]->Reset();
   MATRICES[3]->Reset();
   MATRICES[4]->Reset();
   MATRICES[5]->Reset();
   MATRICES[6]->Reset();
   MATRICES[7]->Reset();
   MATRICES[8]->Reset();
   MATRICES[9]->Reset();
   MATRICES[10]->Reset();
   MATRICES[11]->Reset();
   MATRICES[12]->Reset();
   MATRICES[13]->Reset();
   
    if (Disctype_NSECST == 1)
   {
   MATRICES[14] = MatrixE11;
   MATRICES[15] = MatrixE12;
   MATRICES[16] = MatrixE22;
   MATRICES[17] = MatrixE23;
   
   MATRICES[18] = MatrixJ11;
   MATRICES[19] = MatrixJ21;
   MATRICES[20] = MatrixJ22;
   MATRICES[21] = MatrixJ32;
   N_RectMatrices += 8;
   
   MATRICES[14]->Reset();
   MATRICES[15]->Reset();
   MATRICES[16]->Reset();
   MATRICES[17]->Reset();
   MATRICES[18]->Reset();
   MATRICES[19]->Reset();
   MATRICES[20]->Reset();
   MATRICES[21]->Reset();
   }
  
   #ifdef _SMPI
    if(SOLVER == DIRECT)
     {
      
      if(Disctype_NSECST == 2)
      {
      P_DS = new TSeqParDirectSolver(N_U,N_P,N_S,0,SQMATRICES,MATRICES);
      }
      else if (Disctype_NSECST == 1)
      {
	P_DS = new TSeqParDirectSolver(N_U,N_P,N_S,N_D,SQMATRICES,MATRICES);
      }
     }
   #endif
       
       

  
} // void TSystemTNSECST2D::Init
  
  
 void TSystemTNSECST2D::Assemble(double *sol, double *rhs)
 {
     double *RHSs[8];
    TFESpace2D *fesprhs[8];
    int N_SquareMatrices, N_RectMatrices, N_Rhs;
   
   SQMATRICES[0] = SqmatrixA11;
   SQMATRICES[1] = SqmatrixA12;
   SQMATRICES[2] = SqmatrixA21;
   SQMATRICES[3] = SqmatrixA22;
   SQMATRICES[4] = SqmatrixMU11;
   SQMATRICES[5] = SqmatrixMU22;
   
   SQMATRICES[6] = SqmatrixG11;
   SQMATRICES[7] = SqmatrixG12;
   SQMATRICES[8] = SqmatrixG21;
   SQMATRICES[9] = SqmatrixG22;
   SQMATRICES[10] = SqmatrixG23;
   SQMATRICES[11] = SqmatrixG32;
   SQMATRICES[12] = SqmatrixG33;
   SQMATRICES[13] = SqmatrixMS11;
   SQMATRICES[14] = SqmatrixMS22;
   SQMATRICES[15] = SqmatrixMS33;
   
   N_SquareMatrices = 16;
   
   SQMATRICES[0]->Reset();
   SQMATRICES[1]->Reset();
   SQMATRICES[2]->Reset();
   SQMATRICES[3]->Reset();
   SQMATRICES[4]->Reset();
   SQMATRICES[5]->Reset();
   SQMATRICES[6]->Reset();
   SQMATRICES[7]->Reset();
   SQMATRICES[8]->Reset();
   SQMATRICES[9]->Reset();
   SQMATRICES[10]->Reset();
   SQMATRICES[11]->Reset();
   SQMATRICES[12]->Reset();
   SQMATRICES[13]->Reset();
   SQMATRICES[14]->Reset();
   SQMATRICES[15]->Reset();
     
    if (Disctype_NSECST == 1)
   {
   SQMATRICES[16] = SqmatrixH11;
   SQMATRICES[17] = SqmatrixH22;
   SQMATRICES[18] = SqmatrixH33;
   N_SquareMatrices += 3;
   
   SQMATRICES[16]->Reset();
   SQMATRICES[17]->Reset();
   SQMATRICES[18]->Reset();
   }
   
   MATRICES[0] = MatrixB1;
   MATRICES[1] = MatrixB2;
   
   MATRICES[2] = MatrixB1T;
   MATRICES[3] = MatrixB2T;
   
   MATRICES[4] = MatrixC11;
   MATRICES[5] = MatrixC12;
   MATRICES[6] = MatrixC22;
   MATRICES[7] = MatrixC23;
   
   MATRICES[8] = MatrixD11;
   MATRICES[9] = MatrixD12;
   MATRICES[10] = MatrixD21;
   MATRICES[11] = MatrixD22;
   MATRICES[12] = MatrixD31;
   MATRICES[13] = MatrixD32;
   
   N_RectMatrices = 14;
   
   MATRICES[0]->Reset();
   MATRICES[1]->Reset();
   MATRICES[2]->Reset();
   MATRICES[3]->Reset();
   MATRICES[4]->Reset();
   MATRICES[5]->Reset();
   MATRICES[6]->Reset();
   MATRICES[7]->Reset();
   MATRICES[8]->Reset();
   MATRICES[9]->Reset();
   MATRICES[10]->Reset();
   MATRICES[11]->Reset();
   MATRICES[12]->Reset();
   MATRICES[13]->Reset();
   
     if (Disctype_NSECST == 1)
   {
   MATRICES[14] = MatrixE11;
   MATRICES[15] = MatrixE12;
   MATRICES[16] = MatrixE22;
   MATRICES[17] = MatrixE23;
   
   MATRICES[18] = MatrixJ11;
   MATRICES[19] = MatrixJ21;
   MATRICES[20] = MatrixJ22;
   MATRICES[21] = MatrixJ32;
   N_RectMatrices += 8;
   
   MATRICES[14]->Reset();
   MATRICES[15]->Reset();
   MATRICES[16]->Reset();
   MATRICES[17]->Reset();
   MATRICES[18]->Reset();
   MATRICES[19]->Reset();
   MATRICES[20]->Reset();
   MATRICES[21]->Reset();
   }
   
   N_Rhs = 8;
   RHSs[0] = rhs;
   RHSs[1] = rhs +   N_U;
   RHSs[2] = rhs + 2*N_U;
   RHSs[3] = rhs + 2*N_U +   N_S;
   RHSs[4] = rhs + 2*N_U + 2*N_S;
   RHSs[5] = NLRhs;
   RHSs[6] = NLRhs +   N_S;
   RHSs[7] = NLRhs + 2*N_S;
   
   memset(rhs, 0, N_TotalDof*SizeOfDouble);
   memset(NLRhs, 0, 3*N_S*SizeOfDouble);

     
      fesprhs[0] = FeSpace[0];
      fesprhs[1] = FeSpace[0];
      
      fesprhs[2] = FeSpace[2];
      fesprhs[3] = FeSpace[2];
      fesprhs[4] = FeSpace[2];
      
      fesprhs[5] = FeSpace[2];
      fesprhs[6] = FeSpace[2];
      fesprhs[7] = FeSpace[2];

    // assemble
      Assemble2D(N_FeSpace, FeSpace,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormARhs,
        BoundaryConditions,
        BoundaryValues,
        aux);

      // Add the LPS terms
      if (Disctype_NSECST == 2)
     {
        MatAdd(SQMATRICES[0], SqmatrixA11_LPS,1.);
	MatAdd(SQMATRICES[1], SqmatrixA12_LPS,1.);
	MatAdd(SQMATRICES[2], SqmatrixA21_LPS,1.);
	MatAdd(SQMATRICES[3], SqmatrixA22_LPS,1.);
	
	MatAdd(SQMATRICES[6],  SqmatrixG11_LPS,1.);
	MatAdd(SQMATRICES[7],  SqmatrixG12_LPS,1.);
	MatAdd(SQMATRICES[8],  SqmatrixG21_LPS,1.);
	MatAdd(SQMATRICES[9],  SqmatrixG22_LPS,1.);
	MatAdd(SQMATRICES[10], SqmatrixG23_LPS,1.);
	MatAdd(SQMATRICES[11], SqmatrixG32_LPS,1.);
	MatAdd(SQMATRICES[12], SqmatrixG33_LPS,1.);
     }
     
          int N_RecMat_Slip;
        // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
        SQMATRICES[2] = SqmatrixA12;
        SQMATRICES[3] = SqmatrixA21;
	SQMATRICES[4] = SqmatrixMU11;
        SQMATRICES[5] = SqmatrixMU22;
        SQMATRICES[6] = SqmatrixMU12;
        SQMATRICES[7] = SqmatrixMU21;

        MATRICES[0] = MatrixB1T;
        MATRICES[1] = MatrixB2T;
	N_RecMat_Slip=2;
	
	if (Disctype_NSECST == 2)
	 {
	   MATRICES[2] = MatrixC11;
           MATRICES[3] = MatrixC12;
           MATRICES[4] = MatrixC22;
           MATRICES[5] = MatrixC23;
	   N_RecMat_Slip+=4;
	 }
	 else if(Disctype_NSECST == 1)
	 {
	  cout<<"Not yet implemented \n";
	  exit(4711);
	 }
	 else
	 {
	  cout<<"Invalid Disctype_NSECST \n";
	  exit(4711);
	 }

        Assemble2DSlipBC(1, FeSpace,
                         8, SQMATRICES,
                         N_RecMat_Slip, MATRICES,
                         2, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         aux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=     
 
   
 
   
     // set rhs for Dirichlet nodes
     memcpy(sol +                         N_Active_U, rhs +                         N_Active_U, N_DirichletDof_U*SizeOfDouble);
     memcpy(sol +   N_U +                 N_Active_U, rhs +   N_U +                 N_Active_U, N_DirichletDof_U*SizeOfDouble); 
     memcpy(sol + 2*N_U +                 N_Active_S, rhs + 2*N_U +                 N_Active_S, N_DirichletDof_S*SizeOfDouble);
     memcpy(sol + 2*N_U +   N_S +         N_Active_S, rhs + 2*N_U +   N_S +         N_Active_S, N_DirichletDof_S*SizeOfDouble);
     memcpy(sol + 2*N_U + 2*N_S +         N_Active_S, rhs + 2*N_U + 2*N_S +         N_Active_S, N_DirichletDof_S*SizeOfDouble);
   
     if (Disctype_NSECST == 1)
     {
     memcpy(sol + 2*N_U + 3*N_S +         N_Active_D, rhs + 2*N_U + 3*N_S +         N_Active_D, N_DirichletDof_D*SizeOfDouble);
     memcpy(sol + 2*N_U + 3*N_S +   N_D + N_Active_D, rhs + 2*N_U + 3*N_S +   N_D + N_Active_D, N_DirichletDof_D*SizeOfDouble);
     memcpy(sol + 2*N_U + 3*N_S + 2*N_D + N_Active_D, rhs + 2*N_U + 3*N_S + 2*N_D + N_Active_D, N_DirichletDof_D*SizeOfDouble);
     }
   

 } // void TSystemTNSECST2D::Assemble
  
  void TSystemTNSECST2D::AssembleRhs(double *sol, double *rhs)
 {
  int N_SquareMatrices, N_RectMatrices, N_Rhs;
  
  double *RHSs[5];
  
  TFESpace2D *fesprhs[5];
  
  
      N_Rhs = 5;
      N_SquareMatrices = 0;
      N_RectMatrices = 0;
  
      RHSs[0] = rhs;
      RHSs[1] = rhs +   N_U;
      RHSs[2] = rhs + 2*N_U;
      RHSs[3] = rhs + 2*N_U +   N_S;
      RHSs[4] = rhs + 2*N_U + 2*N_S;

      memset(rhs, 0, N_TotalDof*SizeOfDouble);
     
      fesprhs[0] = FeSpace[0];
      fesprhs[1] = FeSpace[0];
      
      fesprhs[2] = FeSpace[2];
      fesprhs[3] = FeSpace[2];
      fesprhs[4] = FeSpace[2];
  
    // assemble
      Assemble2D(N_FeSpace, FeSpace,
        N_SquareMatrices, NULL,
        N_RectMatrices, NULL,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormRhs,
        BoundaryConditions,
        BoundaryValues,
        aux);    
  }
  
  
 void TSystemTNSECST2D:: AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol)
 {
  
   double tau, val = TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT, temp;
  
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  memset(B, 0, N_TotalDof*SizeOfDouble);

  // old rhs multiplied with current subtime step and theta3 on B
  temp = tau*TDatabase::TimeDB->THETA3; 
  Daxpy(N_Active_U, temp, oldrhs, B);
  Daxpy(N_Active_U, temp, oldrhs+N_U, B+N_U); 
  Daxpy(N_Active_S, temp, oldrhs+2*N_U, B+2*N_U);
  Daxpy(N_Active_S, temp, oldrhs+2*N_U+N_S, B+2*N_U+N_S);
  Daxpy(N_Active_S, temp, oldrhs+2*N_U+2*N_S, B+2*N_U+2*N_S);

  // add rhs from current sub time step to rhs array B
  temp = tau*TDatabase::TimeDB->THETA4;
  Daxpy(N_Active_U, temp, rhs, B);
  Daxpy(N_Active_U, temp, rhs+N_U, B+N_U);
  Daxpy(N_Active_S, temp, rhs+2*N_U, B+2*N_U);
  Daxpy(N_Active_S, temp, rhs+2*N_U+N_S, B+2*N_U+N_S);
  Daxpy(N_Active_S, temp, rhs+2*N_U+2*N_S, B+2*N_U+2*N_S);
  
  // scale the BT, C, E, B, matrices, not in nonlinear step 
  if (scale != 1.0)
  {
         Dscal(MatrixB1T->GetN_Entries(), scale, MatrixB1T->GetEntries());
         Dscal(MatrixB2T->GetN_Entries(), scale, MatrixB2T->GetEntries());
	 
	 temp = scale*TDatabase::TimeDB->THETA1;
	 Dscal(MatrixC11->GetN_Entries(), temp, MatrixC11->GetEntries());
	 Dscal(MatrixC12->GetN_Entries(), temp, MatrixC12->GetEntries());
	 Dscal(MatrixC22->GetN_Entries(), temp, MatrixC22->GetEntries());
	 Dscal(MatrixC23->GetN_Entries(), temp, MatrixC23->GetEntries());	 
       
	 if (Disctype_NSECST == 1)
        { 
         Dscal(MatrixE11->GetN_Entries(), temp, MatrixE11->GetEntries());
	 Dscal(MatrixE12->GetN_Entries(), temp, MatrixE12->GetEntries());
	 Dscal(MatrixE22->GetN_Entries(), temp, MatrixE22->GetEntries());
	 Dscal(MatrixE23->GetN_Entries(), temp, MatrixE23->GetEntries());
        }
      // scale divergence constraint
      if(val>0) 
       {
        Dscal(MatrixB1->GetN_Entries(), val*scale, MatrixB1->GetEntries());
        Dscal(MatrixB2->GetN_Entries(), val*scale, MatrixB2->GetEntries());
       }    
  
  }
 
         temp = tau*TDatabase::TimeDB->THETA1;
	 Dscal(MatrixD11->GetN_Entries(), temp, MatrixD11->GetEntries());
	 Dscal(MatrixD12->GetN_Entries(), temp, MatrixD12->GetEntries());
	 Dscal(MatrixD21->GetN_Entries(), temp, MatrixD21->GetEntries());
	 Dscal(MatrixD22->GetN_Entries(), temp, MatrixD22->GetEntries());
	 Dscal(MatrixD31->GetN_Entries(), temp, MatrixD31->GetEntries());
	 Dscal(MatrixD32->GetN_Entries(), temp, MatrixD32->GetEntries());

       temp = - tau*TDatabase::TimeDB->THETA2;
       MatAdd(SqmatrixMU11, SqmatrixA11, temp);
       MatAdd(SqmatrixMU12, SqmatrixA12, temp);
       MatAdd(SqmatrixMU21, SqmatrixA21, temp);
       MatAdd(SqmatrixMU22, SqmatrixA22, temp);   

       gamma = - tau*TDatabase::TimeDB->THETA2;
       
       memset(defect, 0, N_TotalDof*SizeOfDouble);

       MatVectActive(SqmatrixMU11, sol, defect);
       Daxpy(N_Active_U, 1, defect, B);
       MatVectActive(SqmatrixMU12, sol+N_U, defect);
       Daxpy(N_Active_U, 1, defect, B);
       MatVectActive(SqmatrixMU21, sol, defect+N_U);
       Daxpy(N_Active_U, 1, defect+N_U, B+N_U);
       MatVectActive(SqmatrixMU22, sol+N_U, defect+N_U);
       Daxpy(N_Active_U, 1, defect+N_U, B+N_U);
       
       // Assume we dont use forward euler
       temp = - TDatabase::TimeDB->THETA2/TDatabase::TimeDB->THETA1;
       MatVect1(MatrixC11, sol+2*N_U, defect);
       Daxpy(N_Active_U, temp, defect, B);
       MatVect1(MatrixC12, sol+2*N_U+N_S, defect);
       Daxpy(N_Active_U, temp, defect, B);
       MatVect1(MatrixC22, sol+2*N_U+N_S, defect+N_U);
       Daxpy(N_Active_U, temp, defect+N_U, B+N_U);
       MatVect1(MatrixC23, sol+2*N_U+2*N_S, defect+N_U);
       Daxpy(N_Active_U, temp, defect+N_U, B+N_U);
      
       if (Disctype_NSECST == 1)
      {// Assume we dont use forward euler
       MatVect1(MatrixE11, sol+2*N_U+3*N_S, defect);
       Daxpy(N_Active_U, temp, defect, B);
       MatVect1(MatrixE12, sol+2*N_U+3*N_S+N_D, defect);
       Daxpy(N_Active_U, temp, defect, B);
       MatVect1(MatrixE22, sol+2*N_U+3*N_S+N_D, defect+N_U);
       Daxpy(N_Active_U, temp, defect+N_U, B+N_U);
       MatVect1(MatrixE23, sol+2*N_U+3*N_S+2*N_D, defect+N_U);
       Daxpy(N_Active_U, temp, defect+N_U, B+N_U);
     }
     
       temp = - tau*TDatabase::TimeDB->THETA2;
       MatAdd(SqmatrixMS11, SqmatrixG11, temp);
       MatAdd(SqmatrixMS12, SqmatrixG12, temp);
       MatAdd(SqmatrixMS21, SqmatrixG21, temp);
       MatAdd(SqmatrixMS22, SqmatrixG22, temp);
       MatAdd(SqmatrixMS23, SqmatrixG23, temp);
       MatAdd(SqmatrixMS32, SqmatrixG32, temp);
       MatAdd(SqmatrixMS33, SqmatrixG33, temp);
            
       MatVectActive(SqmatrixMS11, sol+2*N_U, defect+2*N_U);
       Daxpy(N_Active_S, 1, defect+2*N_U, B+2*N_U);
       MatVectActive(SqmatrixMS12, sol+2*N_U+N_S, defect+2*N_U);
       Daxpy(N_Active_S, 1, defect+2*N_U, B+2*N_U);
       MatVectActive(SqmatrixMS21, sol+2*N_U, defect+2*N_U+N_S);
       Daxpy(N_Active_S, 1, defect+2*N_U+N_S, B+2*N_U+N_S);
       MatVectActive(SqmatrixMS22, sol+2*N_U+N_S, defect+2*N_U+N_S);
       Daxpy(N_Active_S, 1, defect+2*N_U+N_S, B+2*N_U+N_S);
       MatVectActive(SqmatrixMS23, sol+2*N_U+2*N_S, defect+2*N_U+N_S);
       Daxpy(N_Active_S, 1, defect+2*N_U+N_S, B+2*N_U+N_S);
       MatVectActive(SqmatrixMS32, sol+2*N_U+N_S, defect+2*N_U+2*N_S);
       Daxpy(N_Active_S, 1, defect+2*N_U+2*N_S, B+2*N_U+2*N_S);
       MatVectActive(SqmatrixMS33, sol+2*N_U+2*N_S, defect+2*N_U+2*N_S);
       Daxpy(N_Active_S, 1, defect+2*N_U+2*N_S, B+2*N_U+2*N_S);
         
       temp = tau*TDatabase::TimeDB->THETA1;
       Daxpy(N_Active_S, temp, NLRhs,       B+2*N_U      );
       Daxpy(N_Active_S, temp, NLRhs+  N_S, B+2*N_U+  N_S);
       Daxpy(N_Active_S, temp, NLRhs+2*N_S, B+2*N_U+2*N_S);
       
       //assembling system matrix
      temp =  -gamma + tau*TDatabase::TimeDB->THETA1;
      
       MatAdd(SqmatrixMU11, SqmatrixA11, temp);
       MatAdd(SqmatrixMU12, SqmatrixA12, temp);
       MatAdd(SqmatrixMU21, SqmatrixA21, temp);
       MatAdd(SqmatrixMU22, SqmatrixA22, temp);  
       
       MatAdd(SqmatrixMS11, SqmatrixG11, temp);
       MatAdd(SqmatrixMS12, SqmatrixG12, temp);
       MatAdd(SqmatrixMS21, SqmatrixG21, temp);
       MatAdd(SqmatrixMS22, SqmatrixG22, temp);
       MatAdd(SqmatrixMS23, SqmatrixG23, temp);
       MatAdd(SqmatrixMS32, SqmatrixG32, temp);
       MatAdd(SqmatrixMS33, SqmatrixG33, temp);
          
       gamma = tau*TDatabase::TimeDB->THETA1;

     // set rhs for Dirichlet nodes
   memcpy(B +                         N_Active_U, rhs +                         N_Active_U, N_DirichletDof_U*SizeOfDouble);
   memcpy(B +   N_U +                 N_Active_U, rhs +   N_U +                 N_Active_U, N_DirichletDof_U*SizeOfDouble); 
   memcpy(B + 2*N_U +                 N_Active_S, rhs + 2*N_U +                 N_Active_S, N_DirichletDof_S*SizeOfDouble);
   memcpy(B + 2*N_U +   N_S +         N_Active_S, rhs + 2*N_U +   N_S +         N_Active_S, N_DirichletDof_S*SizeOfDouble);
   memcpy(B + 2*N_U + 2*N_S +         N_Active_S, rhs + 2*N_U + 2*N_S +         N_Active_S, N_DirichletDof_S*SizeOfDouble);

      if (Disctype_NSECST == 1)
     {
     memcpy(B + 2*N_U + 3*N_S +         N_Active_D, rhs + 2*N_U + 3*N_S +         N_Active_D, N_DirichletDof_D*SizeOfDouble);
     memcpy(B + 2*N_U + 3*N_S +   N_D + N_Active_D, rhs + 2*N_U + 3*N_S +   N_D + N_Active_D, N_DirichletDof_D*SizeOfDouble);
     memcpy(B + 2*N_U + 3*N_S + 2*N_D + N_Active_D, rhs + 2*N_U + 3*N_S + 2*N_D + N_Active_D, N_DirichletDof_D*SizeOfDouble);
     }
   
   SystMatAssembled  = TRUE;

 }  // void TSystemTNSECST2D:: AssembleSystMat
  
  void TSystemTNSECST2D::AssembleNonLinear(double *sol, double *rhs)
 {
   int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
    double *RHSs[3];
    TFESpace2D *fesprhs[3];
      
      SQMATRICES[0] = SqmatrixA11;
      SQMATRICES[1] = SqmatrixA22;
      SQMATRICES[2] = SqmatrixG11;
      SQMATRICES[3] = SqmatrixG12;
      SQMATRICES[4] = SqmatrixG21;
      SQMATRICES[5] = SqmatrixG22;
      SQMATRICES[6] = SqmatrixG23;
      SQMATRICES[7] = SqmatrixG32;
      SQMATRICES[8] = SqmatrixG33;
      
      N_SquareMatrices = 9;
      
      SQMATRICES[0]->Reset();
      SQMATRICES[1]->Reset();
      SQMATRICES[2]->Reset();
      SQMATRICES[3]->Reset();
      SQMATRICES[4]->Reset();
      SQMATRICES[5]->Reset();
      SQMATRICES[6]->Reset();
      SQMATRICES[7]->Reset();
      SQMATRICES[8]->Reset();
      
         
      MATRICES[0] = MatrixD11;
      MATRICES[1] = MatrixD12;
      MATRICES[2] = MatrixD21;
      MATRICES[3] = MatrixD22;
      MATRICES[4] = MatrixD31;
      MATRICES[5] = MatrixD32;
   
      N_RectMatrices = 6;
   
      MATRICES[0]->Reset();
      MATRICES[1]->Reset();
      MATRICES[2]->Reset();
      MATRICES[3]->Reset();
      MATRICES[4]->Reset();
      MATRICES[5]->Reset();
      
      memset(NLRhs, 0, 3*N_S*SizeOfDouble);
      
      N_Rhs = 3;
      RHSs[0] = NLRhs;
      RHSs[1] = NLRhs +   N_S;
      RHSs[2] = NLRhs + 2*N_S;
      fesprhs[0] = FeSpace[2];
      fesprhs[1] = FeSpace[2];
      fesprhs[2] = FeSpace[2];
      
         // assemble the nonlinear part of matrices
      Assemble2D(N_FeSpace, FeSpace,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormNL,
        BoundaryConditions+5,
        BoundaryValues+5,
        aux);  
      
      // Add the LPS terms
      if (Disctype_NSECST == 2)
     {
        MatAdd(SQMATRICES[0], SqmatrixA11_LPS,1.);
	MatAdd(SQMATRICES[1], SqmatrixA22_LPS,1.);
	MatAdd(SQMATRICES[2], SqmatrixG11_LPS,1.);
	MatAdd(SQMATRICES[3], SqmatrixG12_LPS,1.);
	MatAdd(SQMATRICES[4], SqmatrixG21_LPS,1.);
	MatAdd(SQMATRICES[5], SqmatrixG22_LPS,1.);
	MatAdd(SQMATRICES[6], SqmatrixG23_LPS,1.);
	MatAdd(SQMATRICES[7], SqmatrixG32_LPS,1.);
	MatAdd(SQMATRICES[8], SqmatrixG33_LPS,1.);
     }
     

        // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;

        Assemble2DSlipBC(1, FeSpace,
                         2, SQMATRICES,
                         0, MATRICES,
                         0, NULL, NULL,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         aux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=     
 
     
     
      
     // set rhs for Dirichlet nodes
     memcpy(sol +                         N_Active_U, rhs +                         N_Active_U, N_DirichletDof_U*SizeOfDouble);
     memcpy(sol +   N_U +                 N_Active_U, rhs +   N_U +                 N_Active_U, N_DirichletDof_U*SizeOfDouble); 
     memcpy(sol + 2*N_U +                 N_Active_S, rhs + 2*N_U +                 N_Active_S, N_DirichletDof_S*SizeOfDouble);
     memcpy(sol + 2*N_U +   N_S +         N_Active_S, rhs + 2*N_U +   N_S +         N_Active_S, N_DirichletDof_S*SizeOfDouble);
     memcpy(sol + 2*N_U + 2*N_S +         N_Active_S, rhs + 2*N_U + 2*N_S +         N_Active_S, N_DirichletDof_S*SizeOfDouble);
   
     if (Disctype_NSECST == 1)
     {
     memcpy(sol + 2*N_U + 3*N_S +         N_Active_D, rhs + 2*N_U + 3*N_S +         N_Active_D, N_DirichletDof_D*SizeOfDouble);
     memcpy(sol + 2*N_U + 3*N_S +   N_D + N_Active_D, rhs + 2*N_U + 3*N_S +   N_D + N_Active_D, N_DirichletDof_D*SizeOfDouble);
     memcpy(sol + 2*N_U + 3*N_S + 2*N_D + N_Active_D, rhs + 2*N_U + 3*N_S + 2*N_D + N_Active_D, N_DirichletDof_D*SizeOfDouble);
     }
 }
  
  /* assemble only LHS, and non-linear term in rhs */
  void TSystemTNSECST2D::AssembleSystMatNonLinear()
 {
   double temp;
    if(SystMatAssembled)
    {
     cout << "Restore System mat before calling AssembleSystMatNonLinear" <<endl;
    }
    
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    
    temp = tau*TDatabase::TimeDB->THETA1;
    
    Dscal(MatrixD11->GetN_Entries(), temp, MatrixD11->GetEntries());
    Dscal(MatrixD12->GetN_Entries(), temp, MatrixD12->GetEntries());
    Dscal(MatrixD21->GetN_Entries(), temp, MatrixD21->GetEntries());
    Dscal(MatrixD22->GetN_Entries(), temp, MatrixD22->GetEntries());
    Dscal(MatrixD31->GetN_Entries(), temp, MatrixD31->GetEntries());
    Dscal(MatrixD32->GetN_Entries(), temp, MatrixD32->GetEntries());

    MatAdd(SqmatrixMU11, SqmatrixA11, temp);
    MatAdd(SqmatrixMU12, SqmatrixA12, temp);
    MatAdd(SqmatrixMU21, SqmatrixA21, temp);
    MatAdd(SqmatrixMU22, SqmatrixA22, temp);  
       
    MatAdd(SqmatrixMS11, SqmatrixG11, temp);
    MatAdd(SqmatrixMS12, SqmatrixG12, temp);
    MatAdd(SqmatrixMS21, SqmatrixG21, temp);
    MatAdd(SqmatrixMS22, SqmatrixG22, temp);
    MatAdd(SqmatrixMS23, SqmatrixG23, temp);
    MatAdd(SqmatrixMS32, SqmatrixG32, temp);
    MatAdd(SqmatrixMS33, SqmatrixG33, temp);
    
    Daxpy(N_Active_S, temp, NLRhs,       B+2*N_U      );
    Daxpy(N_Active_S, temp, NLRhs+  N_S, B+2*N_U+  N_S);
    Daxpy(N_Active_S, temp, NLRhs+2*N_S, B+2*N_U+2*N_S);
    
    gamma = tau*TDatabase::TimeDB->THETA1;
    
    SystMatAssembled  = TRUE;
 }
  
  void TSystemTNSECST2D::RestoreMassMat()
{

  if(SystMatAssembled)
   {
       //assembling system matrix
       MatAdd(SqmatrixMU11, SqmatrixA11, -gamma);
       MatAdd(SqmatrixMU12, SqmatrixA12, -gamma);
       MatAdd(SqmatrixMU21, SqmatrixA21, -gamma);
       MatAdd(SqmatrixMU22, SqmatrixA22, -gamma); 
       
       MatAdd(SqmatrixMS11, SqmatrixG11, -gamma);
       MatAdd(SqmatrixMS12, SqmatrixG12, -gamma);
       MatAdd(SqmatrixMS21, SqmatrixG21, -gamma);
       MatAdd(SqmatrixMS22, SqmatrixG22, -gamma);
       MatAdd(SqmatrixMS23, SqmatrixG23, -gamma);
       MatAdd(SqmatrixMS32, SqmatrixG32, -gamma);
       MatAdd(SqmatrixMS33, SqmatrixG33, -gamma);
       
       Daxpy(N_Active_S, -gamma, NLRhs,       B+2*N_U      );
       Daxpy(N_Active_S, -gamma, NLRhs+  N_S, B+2*N_U+  N_S);
       Daxpy(N_Active_S, -gamma, NLRhs+2*N_S, B+2*N_U+2*N_S);
       
       gamma = 0.;     
    
    SystMatAssembled  = FALSE;  
   }
  else
  {
    cout << "System is not assembled to restore " <<endl;
  }
   
  
} //void TSystemTNSECST2D::RestoreMassMat()
  
void TSystemTNSECST2D::Solve(double *sol)
{  
  if(!SystMatAssembled)
  {
    cout << "System Matrix is not assembled to solve " <<endl;
    exit(0);
  }
  
    switch(SOLVER)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
	if (Disctype_NSECST == 1)
        {
   #ifdef _SEQ
   #ifdef _SMPI
      P_DS->Solve(sol,B,true);
    #else        
      
	  	DirectSolver(SqmatrixMU11, SqmatrixMU12, SqmatrixMU21, SqmatrixMU22, 
		     SqmatrixMS11, SqmatrixMS12, SqmatrixMS21, SqmatrixMS22, SqmatrixMS23, SqmatrixMS32, SqmatrixMS33,
                     SqmatrixH11, SqmatrixH22, SqmatrixH33,
	             MatrixB1,  MatrixB2, MatrixB1T, MatrixB2T, 
	             MatrixC11, MatrixC12, MatrixC22, MatrixC23,
	             MatrixD11, MatrixD12, MatrixD21, MatrixD22, MatrixD31, MatrixD32,
	             MatrixE11, MatrixE12, MatrixE22, MatrixE23,
	             MatrixJ11, MatrixJ21, MatrixJ22, MatrixJ32,
	             B, sol); 
		  
   #endif
   #endif
        }
       else if(Disctype_NSECST == 2)
         {
   #ifdef _SEQ
   #ifdef _SMPI
      P_DS->Solve(sol,B,true);
    #else  
	     DirectSolver(SqmatrixMU11, SqmatrixMU12, SqmatrixMU21, SqmatrixMU22, 
		     SqmatrixMS11, SqmatrixMS12, SqmatrixMS21, SqmatrixMS22, SqmatrixMS23, SqmatrixMS32, SqmatrixMS33,                  
	             MatrixB1,  MatrixB2, MatrixB1T, MatrixB2T, 
	             MatrixC11, MatrixC12, MatrixC22, MatrixC23,
	             MatrixD11, MatrixD12, MatrixD21, MatrixD22, MatrixD31, MatrixD32,
	             B, sol); 
	     
  #endif
  #endif
	 }
	     break;
 
      default:

            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    

} // void TSystemTNSECST2D::Solve
  
  void TSystemTNSECST2D::GetTNSECSTResidual(double *sol, double *res)
{
  
    if(!SystMatAssembled)
   {
    cout << "System Matrix is not assembled to calculate residual " <<endl;
    exit(0);
   }
   double *defect = new double[N_TotalDof];
   
   memset(defect, 0, N_TotalDof*SizeOfDouble);
   memset(res, 0, N_TotalDof*SizeOfDouble);
   
   // Block A
   MatVect(SqmatrixMU11, sol, defect);
   Daxpy(N_U, 1, defect, res);
   MatVectActive(SqmatrixMU12, sol+N_U, defect);
   Daxpy(N_Active_U, 1, defect, res);
   MatVectActive(SqmatrixMU21, sol, defect+N_U);
   Daxpy(N_Active_U, 1, defect+N_U, res+N_U);
   MatVect(SqmatrixMU22, sol+N_U, defect+N_U);
   Daxpy(N_U, 1, defect+N_U, res+N_U);
    
   // Block G
   MatVect(SqmatrixMS11, sol+2*N_U, defect+2*N_U);
   Daxpy(N_S, 1, defect+2*N_U, res+2*N_U);
   MatVectActive(SqmatrixMS12, sol+2*N_U+N_S, defect+2*N_U);
   Daxpy(N_Active_S, 1, defect+2*N_U, res+2*N_U);
   MatVectActive(SqmatrixMS21, sol+2*N_U, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVect(SqmatrixMS22, sol+2*N_U+N_S, defect+2*N_U+N_S);
   Daxpy(N_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVectActive(SqmatrixMS23, sol+2*N_U+2*N_S, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVectActive(SqmatrixMS32, sol+2*N_U+N_S, defect+2*N_U+2*N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+2*N_S, res+2*N_U+2*N_S);
   MatVect(SqmatrixMS33, sol+2*N_U+2*N_S, defect+2*N_U+2*N_S);
   Daxpy(N_S, 1, defect+2*N_U+2*N_S, res+2*N_U+2*N_S);
   
   // Block C
   MatVect1(MatrixC11, sol+2*N_U, defect);
   Daxpy(N_Active_U, 1, defect, res);
   MatVect1(MatrixC12, sol+2*N_U+N_S, defect);
   Daxpy(N_Active_U, 1, defect, res);
   MatVect1(MatrixC22, sol+2*N_U+N_S, defect+N_U);
   Daxpy(N_Active_U, 1, defect+N_U, res+N_U);
   MatVect1(MatrixC23, sol+2*N_U+2*N_S, defect+N_U);
   Daxpy(N_Active_U, 1, defect+N_U, res+N_U);

   
   // Block D
   MatVect1(MatrixD11, sol, defect+2*N_U);
   Daxpy(N_Active_S, 1, defect+2*N_U, res+2*N_U);
   MatVect1(MatrixD12, sol+N_U, defect+2*N_U);
   Daxpy(N_Active_S, 1, defect+2*N_U, res+2*N_U);
   MatVect1(MatrixD21, sol, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVect1(MatrixD22, sol+N_U, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVect1(MatrixD31, sol, defect+2*N_U+2*N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+2*N_S, res+2*N_U+2*N_S);
   MatVect1(MatrixD32, sol+N_U, defect+2*N_U+2*N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+2*N_S, res+2*N_U+2*N_S);   
      
        if (Disctype_NSECST == 1)
     {
       // Block H
       MatVect(SqmatrixH11, sol+2*N_U+3*N_S, defect+2*N_U+3*N_S);
       Daxpy(N_D, 1, defect+2*N_U+3*N_S, res+2*N_U+3*N_S);
       MatVect(SqmatrixH22, sol+2*N_U+3*N_S+N_D, defect+2*N_U+3*N_S+N_D);
       Daxpy(N_D, 1, defect+2*N_U+3*N_S+N_D, res+2*N_U+3*N_S+N_D);
       MatVect(SqmatrixH33, sol+2*N_U+3*N_S+2*N_D, defect+2*N_U+3*N_S+2*N_D);
       Daxpy(N_D, 1, defect+2*N_U+3*N_S+2*N_D, res+2*N_U+3*N_S+2*N_D);
       
       // Block E
       MatVect1(MatrixE11, sol+2*N_U+3*N_S, defect);
       Daxpy(N_Active_U, 1, defect, res);
       MatVect1(MatrixE12, sol+2*N_U+3*N_S+N_D, defect);
       Daxpy(N_Active_U, 1, defect, res);
       MatVect1(MatrixE22, sol+2*N_U+3*N_S+N_D, defect+N_U);
       Daxpy(N_Active_U, 1, defect+N_U, res+N_U);
       MatVect1(MatrixE23, sol+2*N_U+3*N_S+2*N_D, defect+N_U);
       Daxpy(N_Active_U, 1, defect+N_U, res+N_U);
       
       // Block J
       MatVect1(MatrixJ11, sol, defect+2*N_U+3*N_S);
       Daxpy(N_Active_D, 1, defect+2*N_U+3*N_S, res+2*N_U+3*N_S);
       MatVect1(MatrixJ21, sol, defect+2*N_U+3*N_S+N_D);
       Daxpy(N_Active_D, 1, defect+2*N_U+3*N_S+N_D, res+2*N_U+3*N_S+N_D);
       MatVect1(MatrixJ22, sol+N_U, defect+2*N_U+3*N_S+N_D);
       Daxpy(N_Active_D, 1, defect+2*N_U+3*N_S+N_D, res+2*N_U+3*N_S+N_D);
       MatVect1(MatrixJ32, sol+N_U, defect+2*N_U+3*N_S+2*N_D);
       Daxpy(N_Active_D, 1, defect+2*N_U+3*N_S+2*N_D, res+2*N_U+3*N_S+2*N_D);
       
      // Block BT
      MatVect1(MatrixB1T, sol+2*N_U+3*N_S+3*N_D, defect);
      Daxpy(N_Active_U, 1, defect, res); 
      MatVect1(MatrixB2T, sol+2*N_U+3*N_S+3*N_D, defect+N_U);
      Daxpy(N_Active_U, 1, defect+N_U, res+N_U); 
       
      // Block B
      MatVect1(MatrixB1, sol, defect+2*N_U+3*N_S+3*N_D);
      Daxpy(N_P, 1, defect+2*N_U+3*N_S+3*N_D, res+2*N_U+3*N_S+3*N_D);
      MatVect1(MatrixB2, sol+N_U, defect+2*N_U+3*N_S+3*N_D);
      Daxpy(N_P, 1, defect+2*N_U+3*N_S+3*N_D, res+2*N_U+3*N_S+3*N_D); 
     }
     else if (Disctype_NSECST == 2)
     {
      // Block BT
      MatVect1(MatrixB1T, sol+2*N_U+3*N_S, defect);
      Daxpy(N_Active_U, 1, defect, res); 
      MatVect1(MatrixB2T, sol+2*N_U+3*N_S, defect+N_U);
      Daxpy(N_Active_U, 1, defect+N_U, res+N_U); 
      
      // Block B
      MatVect1(MatrixB1, sol, defect+2*N_U+3*N_S);
      Daxpy(N_P, 1, defect+2*N_U+3*N_S, res+2*N_U+3*N_S); 
      MatVect1(MatrixB2, sol+N_U, defect+2*N_U+3*N_S);
      Daxpy(N_P, 1, defect+2*N_U+3*N_S, res+2*N_U+3*N_S);       
     }

     Daxpy(N_TotalDof, -1., B, res);
     

}
  
void TSystemTNSECST2D::MeasureTNSECSTErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3,
				    double *AllErrors)
{
  double errors[4], u_error[4], s_error[6];
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
     // errors in first velocity component
     FeFct[0]->GetErrors(ExactU1, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, auxerror_NSE, 1, FeSpace, errors);
      u_error[0] = errors[0];
      u_error[1] = errors[1];
      
      
     // errors in second velocity component
     FeFct[1]->GetErrors(ExactU2, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, auxerror_NSE, 1, FeSpace, errors);
     u_error[2] = errors[0];
     u_error[3] = errors[1]; 
     
      
     AllErrors[0] = sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]);
     AllErrors[1] = sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]);
     
      
      // errors in pressure
     FeFct[2]->GetErrors(ExactP, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, auxerror_NSE, 1, FeSpace+1, errors);     
     AllErrors[2] = errors[0];
     AllErrors[3] = errors[1];

     
     // errors in first stress component
     FeFct[3]->GetErrors(ExactS1, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, auxerror_CST, 1, FeSpace+2, errors);
      s_error[0] = errors[0];
      s_error[1] = errors[1];
      
     // errors in second stress component
     FeFct[4]->GetErrors(ExactS2, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, auxerror_CST, 1, FeSpace+2, errors);
     s_error[2] = errors[0];
     s_error[3] = errors[1];      
      
      // errors in third stress component
     FeFct[5]->GetErrors(ExactS3, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, auxerror_CST, 1, FeSpace+2, errors);     
     s_error[4] = errors[0];
     s_error[5] = errors[1]; 
  
     
     AllErrors[4] = sqrt(s_error[0]*s_error[0]+(2*s_error[2]*s_error[2])+s_error[4]*s_error[4]);
     AllErrors[5] = sqrt(s_error[1]*s_error[1]+(2*s_error[3]*s_error[3])+s_error[5]*s_error[5]);
  
     
      // error in L^infty(0,t,L^2) for velocity
      if(AllErrors[0] > AllErrors[7])
       {
        AllErrors[7]  = AllErrors[0];
        AllErrors[6]  =  TDatabase::TimeDB->CURRENTTIME;
      }
   
      // error in L^2(0,t,L^2)  for velocity  
      AllErrors[8] += (u_error[0]*u_error[0] + u_error[2]*u_error[2] +olderror_l_2_l_2u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_l_2u = u_error[0]*u_error[0] + u_error[2]*u_error[2];
     
      
      
      // error in L^infty(0,t,L^2) for stress
      if(AllErrors[4] > AllErrors[10])
       {
        AllErrors[10]  = AllErrors[4];
        AllErrors[9]  =  TDatabase::TimeDB->CURRENTTIME;
      }
      
       // error in L^2(0,t,L^2)  for  stress  
      AllErrors[11] += ( s_error[0]*s_error[0]+(2*s_error[2]*s_error[2])+s_error[4]*s_error[4] + olderror_l_2_l_2s)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_l_2s = s_error[0]*s_error[0]+(2*s_error[2]*s_error[2])+s_error[4]*s_error[4] ;
     
      
      // error in L^infty(0,t,H^1) for velocity
      if(AllErrors[1] > AllErrors[13])
       {
        AllErrors[13]  = AllErrors[1];
        AllErrors[12]  =  TDatabase::TimeDB->CURRENTTIME;
      }
   
      // error in L^2(0,t,H^1)    
      AllErrors[14] += (u_error[1]*u_error[1] + u_error[3]*u_error[3] +olderror_l_2_h_1u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_h_1u = u_error[1]*u_error[1] + u_error[3]*u_error[3];
      
        // error in L^infty(0,t,H^1) for stress
      if(AllErrors[5] > AllErrors[10])
       {
        AllErrors[16]  = AllErrors[5];
        AllErrors[15]  =  TDatabase::TimeDB->CURRENTTIME;
      }
      
       // error in L^2(0,t,H^1)  for  stress  
      AllErrors[17] += ( s_error[1]*s_error[1]+(2*s_error[3]*s_error[3])+s_error[5]*s_error[5] + olderror_l_2_h_1s)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_h_1s = s_error[1]*s_error[1]+(2*s_error[3]*s_error[3])+s_error[5]*s_error[5] ;
      
      
      // error in L^infty(0,t,L^2) for pressure
      if(AllErrors[2] > AllErrors[19])
       {
        AllErrors[19]  = AllErrors[2];
        AllErrors[18]  =  TDatabase::TimeDB->CURRENTTIME;
      }
      
      // error in L^2(0,t,L^2)  for pressure  
      AllErrors[20] += (AllErrors[2]*AllErrors[2] +olderror_l_2_l_2p)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_l_2p = AllErrors[2]*AllErrors[2];
      
      // error in L^infty(0,t,H^1) for pressure
      if(AllErrors[3] > AllErrors[22])
       {
        AllErrors[22]  = AllErrors[3];
        AllErrors[21]  =  TDatabase::TimeDB->CURRENTTIME;
      }
      
      // error in L^2(0,t,H^1)  for pressure  
      AllErrors[23] += (AllErrors[3]*AllErrors[3] +olderror_l_2_h_1p)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_h_1p = AllErrors[3]*AllErrors[3];
      
      
     
}  // void TSystemTNSECST2D::MeasureTNSECSTErrors
  
  
  
#endif   // #ifdef __2D__
