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
* @brief     source file for TSystemNSECST2D
* @author    Jagannath Venkatesan, 
* @date      23.02.16
* @History   31.08.16 : Included GetResidual
*            13.09.16 : Checked - Both Sequential and SMPI working
 ************************************************************************  */
#ifdef __2D__
#include <Database.h>
#include <SystemNSECST2D.h>
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

TSystemNSECST2D::TSystemNSECST2D(int N_FESpaces, TFESpace2D **FE_Spaces, int N_FEFunctions, 
				 TFEFunction2D **FeFunctions, int disctype_NSECST, 
				 int solver)
{
  int i;
  N_FeSpace = N_FESpaces;
  N_FeFunction = N_FEFunctions;

  for (i=0;i<N_FeSpace;i++)
  FeSpace[i] = FE_Spaces[i];
  
  for (i=0;i<N_FeFunction;i++)
  FeFct[i] = FeFunctions[i];
  
  //set the discretization type (1-DEVSS/SUPG, 2-LPS)
  Disctype_NSECST = disctype_NSECST;
  
  //set the solver type
  SOLVER = solver;
 
  // velocity and pressure
  N_U = FeSpace[0]->GetN_DegreesOfFreedom();
  N_P = FeSpace[1]->GetN_DegreesOfFreedom(); 
  N_Active_U =  FeSpace[0]->GetActiveBound();
  N_DirichletDof_U = N_U - N_Active_U;  
    
  // conformation stress
  N_S = FeSpace[2]->GetN_DegreesOfFreedom();
  N_Active_S =  FeSpace[2]->GetActiveBound();
  N_DirichletDof_S = N_S - N_Active_S;
   
  N_TotalDof = 2*N_U + 3*N_S + N_P;

  // Deformation tensor
  if (Disctype_NSECST == 1)
  {
  N_D = FeSpace[3]->GetN_DegreesOfFreedom();
  N_Active_D =  FeSpace[3]->GetActiveBound();
  N_DirichletDof_D = N_D - N_Active_D;
  N_TotalDof += 3*N_D;
  }
  
  // build matrices
  // first build matrix structure
  sqstructureA = new TSquareStructure2D(FeSpace[0]);
  sqstructureA->Sort();  // sort column numbers: numbers are in increasing order
  
  structureB = new TStructure2D(FeSpace[1], FeSpace[0]);
  structureBT = new TStructure2D(FeSpace[0], FeSpace[1]);
  structureB->Sort();
  structureBT->Sort();
  
// As of now considering NSEType 4 only
  SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
  SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
  SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
  SqmatrixA22 = new TSquareMatrix2D(sqstructureA);
  
  MatrixB1 = new TMatrix2D(structureB);
  MatrixB2 = new TMatrix2D(structureB);
  MatrixB1T = new TMatrix2D(structureBT);
  MatrixB2T = new TMatrix2D(structureBT);

  sqstructureG = new TSquareStructure2D(FeSpace[2]);
  sqstructureG->Sort();  
  SqmatrixG11 = new TSquareMatrix2D(sqstructureG);  
  SqmatrixG12 = new TSquareMatrix2D(sqstructureG); 
  SqmatrixG21 = new TSquareMatrix2D(sqstructureG); 
  SqmatrixG22 = new TSquareMatrix2D(sqstructureG); 
  SqmatrixG23 = new TSquareMatrix2D(sqstructureG);
  SqmatrixG32 = new TSquareMatrix2D(sqstructureG); 
  SqmatrixG33 = new TSquareMatrix2D(sqstructureG);
  
  structureC = new TStructure2D(FeSpace[0], FeSpace[2]);
  structureC->Sort();
  MatrixC11 = new TMatrix2D(structureC);
  MatrixC12 = new TMatrix2D(structureC);
  MatrixC22 = new TMatrix2D(structureC);
  MatrixC23 = new TMatrix2D(structureC);
  
  structureD = new TStructure2D(FeSpace[2], FeSpace[0]);
  structureD->Sort();
  MatrixD11 = new TMatrix2D(structureD);
  MatrixD12 = new TMatrix2D(structureD);
  MatrixD21 = new TMatrix2D(structureD);
  MatrixD22 = new TMatrix2D(structureD);
  MatrixD31 = new TMatrix2D(structureD);
  MatrixD32 = new TMatrix2D(structureD);
  
    if (Disctype_NSECST == 1)
    {
      sqstructureH = new TSquareStructure2D(FeSpace[3]);
      sqstructureH->Sort();
      SqmatrixH11 = new TSquareMatrix2D(sqstructureH);  
      SqmatrixH22 = new TSquareMatrix2D(sqstructureH); 
      SqmatrixH33 = new TSquareMatrix2D(sqstructureH);
  
      structureE = new TStructure2D(FeSpace[0], FeSpace[3]);
      structureE->Sort();
      MatrixE11 = new TMatrix2D(structureE);
      MatrixE12 = new TMatrix2D(structureE);
      MatrixE22 = new TMatrix2D(structureE);
      MatrixE23 = new TMatrix2D(structureE);
  
      structureJ = new TStructure2D(FeSpace[3], FeSpace[0]);
      structureJ->Sort();
      MatrixJ11 = new TMatrix2D(structureJ);
      MatrixJ21 = new TMatrix2D(structureJ);
      MatrixJ22 = new TMatrix2D(structureJ);
      MatrixJ32 = new TMatrix2D(structureJ);
    }
  
           if (Disctype_NSECST == 2)
      {   // Matrices for LPS terms
          // As of now considering NSEType 4 only
       SqmatrixA11_LPS = new TSquareMatrix2D(sqstructureA);
       SqmatrixA12_LPS = new TSquareMatrix2D(sqstructureA);
       SqmatrixA21_LPS = new TSquareMatrix2D(sqstructureA);
       SqmatrixA22_LPS = new TSquareMatrix2D(sqstructureA);
       
       SqmatrixG11_LPS = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG12_LPS = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG21_LPS = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG22_LPS = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG23_LPS = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG32_LPS = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG33_LPS = new TSquareMatrix2D(sqstructureG);  
       
       SqmatrixG11_LPS_Streamline = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG22_LPS_Streamline = new TSquareMatrix2D(sqstructureG); 
       SqmatrixG33_LPS_Streamline = new TSquareMatrix2D(sqstructureG); 
      }
  
    
 // matrices for methods
   sqmatrices = (TSquareMatrix **)SQMATRICES;
   matrices = (TMatrix **)MATRICES;

   aux = NULL;
   auxerror_NSE = NULL;
   auxerror_CST = NULL;
   
} // TSystemNSECST2D::TSystemNSECST2D

TSystemNSECST2D::~TSystemNSECST2D()
{
   delete sqstructureA; delete sqstructureG; 
   delete structureB; delete structureBT; delete structureC;
    delete structureD; 
   delete SqmatrixA11; delete SqmatrixA12;
   delete SqmatrixA21; delete SqmatrixA22;
   delete SqmatrixG11; delete SqmatrixG12; delete SqmatrixG21;
   delete SqmatrixG22; delete SqmatrixG23; delete SqmatrixG32;
   delete SqmatrixG33;  delete MatrixB1; delete MatrixB2;
   delete MatrixB1T; delete MatrixB2T; delete MatrixC11;
   delete MatrixC12; delete MatrixC22; delete MatrixC23;
   delete MatrixD11; delete MatrixD12;
   delete MatrixD21; delete MatrixD22; delete MatrixD31; delete MatrixD32; 
      if (Disctype_NSECST == 1)
    {  
    delete sqstructureH;
    delete structureJ;
    delete structureE;
   delete MatrixJ11; delete MatrixJ21;
   delete MatrixJ22; delete MatrixJ32; 
   delete SqmatrixH11; delete SqmatrixH22; delete SqmatrixH33; 
   delete MatrixE11; delete MatrixE12; delete MatrixE22;
   delete MatrixE23;
    }
}
  
  
void TSystemNSECST2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D **BoundCond, BoundValueFunct2D **BoundValue,  
			   TAuxParam2D *NSECSTaux, TAuxParam2D *NSEauxerror, TAuxParam2D *CSTauxerror)
{
  int i;
  TDiscreteForm2D *DiscreteFormGalerkinSUPG, *DiscreteFormLPS;
  
  for (i=0;i<N_FeFunction-1;i++)
  { 
  // save the boundary condition
  BoundaryConditions[i] = BoundCond[i];
  // save the boundary values  
  BoundaryValues[i] = BoundValue[i];
  }
  
  // save the bilinear coefficients
  LinCoeffs[0] = lincoeffs;
  
  // aux for assembling and error calculation
  aux = NSECSTaux;
  auxerror_NSE = NSEauxerror;
  auxerror_CST = CSTauxerror;
  
  // set the Discreteforms
   InitializeDiscreteForms_NSECST(DiscreteFormGalerkinSUPG, DiscreteFormLPS, LinCoeffs[0]);
   
     if (Disctype_NSECST == 1)
      DiscreteFormARhs = DiscreteFormGalerkinSUPG;
    else if (Disctype_NSECST == 2)
      DiscreteFormARhs = DiscreteFormLPS;
    else 
    {
      cout<<"Invalid Disctype_NSECST !!!\n";
      exit(0);
    }

   // Matrix initialization due to SMPI (Otherwise not required)
   int N_SquareMatrices, N_RectMatrices, N_Rhs;

   SQMATRICES[0] = SqmatrixA11;
   SQMATRICES[1] = SqmatrixA12;
   SQMATRICES[2] = SqmatrixA21;
   SQMATRICES[3] = SqmatrixA22;
   
   SQMATRICES[4] = SqmatrixG11;
   SQMATRICES[5] = SqmatrixG12;
   SQMATRICES[6] = SqmatrixG21;
   SQMATRICES[7] = SqmatrixG22;
   SQMATRICES[8] = SqmatrixG23;
   SQMATRICES[9] = SqmatrixG32;
   SQMATRICES[10] = SqmatrixG33;
  
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
   
} 


void TSystemNSECST2D::Assemble(double *sol, double *rhs)
{

    double *RHSs[9];
    TFESpace2D *fesprhs[9];
   
    
    int N_SquareMatrices, N_RectMatrices, N_Rhs;

   SQMATRICES[0] = SqmatrixA11;
   SQMATRICES[1] = SqmatrixA12;
   SQMATRICES[2] = SqmatrixA21;
   SQMATRICES[3] = SqmatrixA22;
   
   SQMATRICES[4] = SqmatrixG11;
   SQMATRICES[5] = SqmatrixG12;
   SQMATRICES[6] = SqmatrixG21;
   SQMATRICES[7] = SqmatrixG22;
   SQMATRICES[8] = SqmatrixG23;
   SQMATRICES[9] = SqmatrixG32;
   SQMATRICES[10] = SqmatrixG33;
  
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
    
   N_Rhs = 5;
   RHSs[0] = rhs;
   RHSs[1] = rhs +   N_U;
   RHSs[2] = rhs + 2*N_U;
   RHSs[3] = rhs + 2*N_U +   N_S;
   RHSs[4] = rhs + 2*N_U + 2*N_S;
   RHSs[5] = rhs + 2*N_U + 3*N_S;
   
   if (Disctype_NSECST == 1)
   {
   RHSs[6] = rhs + 2*N_U + 3*N_S +   N_D;
   RHSs[7] = rhs + 2*N_U + 3*N_S + 2*N_D;
   RHSs[8] = rhs + 2*N_U + 3*N_S + 3*N_D;
   }
   
   memset(rhs, 0, N_TotalDof*SizeOfDouble);

     
      fesprhs[0] = FeSpace[0];
      fesprhs[1] = FeSpace[0];
      fesprhs[2] = FeSpace[2];
      fesprhs[3] = FeSpace[2];
      fesprhs[4] = FeSpace[2];
      
      if (Disctype_NSECST == 1)
   {
      fesprhs[5] = FeSpace[3];
      fesprhs[6] = FeSpace[3];
      fesprhs[7] = FeSpace[3];
      fesprhs[8] = FeSpace[1];
   }
   else if(Disctype_NSECST == 2)
   {
     fesprhs[5] = FeSpace[1];
   }
  
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
	
	MatAdd(SQMATRICES[4], SqmatrixG11_LPS,1.);
	MatAdd(SQMATRICES[5], SqmatrixG12_LPS,1.);
	MatAdd(SQMATRICES[6], SqmatrixG21_LPS,1.);
	MatAdd(SQMATRICES[7], SqmatrixG22_LPS,1.);
	MatAdd(SQMATRICES[8], SqmatrixG23_LPS,1.);
	MatAdd(SQMATRICES[9], SqmatrixG32_LPS,1.);
	MatAdd(SQMATRICES[10],SqmatrixG33_LPS,1.);
     }
           
      int N_RecMat_Slip;
        // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
        SQMATRICES[2] = SqmatrixA12;
        SQMATRICES[3] = SqmatrixA21;

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
	
	 if (Disctype_NSECST == 1)
	 {
	   MATRICES[2] = MatrixC11;
           MATRICES[3] = MatrixC12;
           MATRICES[4] = MatrixC22;
           MATRICES[5] = MatrixC23;
	   MATRICES[6] = MatrixE11;
           MATRICES[7] = MatrixE12;
           MATRICES[8] = MatrixE22;
           MATRICES[9] = MatrixE23;
	   N_RecMat_Slip+=8;
	 }
	

        Assemble2DSlipBC(1, FeSpace,
                         4, SQMATRICES,
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
     
     
//      double *Entries = SqmatrixA11->GetEntries(), sum=0;
//      int *RowPtr = SqmatrixA11->GetRowPtr();
//      for(int i=0;i<N_U;i++)
//   {
//     for(int j=RowPtr[i];j<RowPtr[i+1];j++)
//       sum += Entries[j]*Entries[j];
//   }
//       cout<< sum<<"\n";
//      
    


} // void TSystemNSECST2D::Assemble(

  

void TSystemNSECST2D::Solve(double *sol, double *rhs)
{
  
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
      P_DS->Solve(sol,rhs,true);
    #else     
	DirectSolver(SqmatrixA11, SqmatrixA12, SqmatrixA21, SqmatrixA22, 
		     SqmatrixG11, SqmatrixG12, SqmatrixG21, SqmatrixG22, SqmatrixG23, SqmatrixG32, SqmatrixG33,
                     SqmatrixH11, SqmatrixH22, SqmatrixH33,
	             MatrixB1,  MatrixB2, MatrixB1T, MatrixB2T, 
	             MatrixC11, MatrixC12, MatrixC22, MatrixC23,
	             MatrixD11, MatrixD12, MatrixD21, MatrixD22, MatrixD31, MatrixD32,
	             MatrixE11, MatrixE12, MatrixE22, MatrixE23,
	             MatrixJ11, MatrixJ21, MatrixJ22, MatrixJ32,
	             rhs, sol); 	
#endif
#endif
}
else if(Disctype_NSECST == 2)
{
   #ifdef _SEQ
   #ifdef _SMPI
      P_DS->Solve(sol,rhs,true);
    #else  
 	DirectSolver(SqmatrixA11, SqmatrixA12, SqmatrixA21, SqmatrixA22, 
		     SqmatrixG11, SqmatrixG12, SqmatrixG21, SqmatrixG22, SqmatrixG23, SqmatrixG32, SqmatrixG33,                  
	             MatrixB1,  MatrixB2, MatrixB1T, MatrixB2T, 
	             MatrixC11, MatrixC12, MatrixC22, MatrixC23,
	             MatrixD11, MatrixD12, MatrixD21, MatrixD22, MatrixD31, MatrixD32,
	             rhs, sol); 
 #endif
#endif
}
	break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
}    // void TSystemNSECST2D::Solve


  void TSystemNSECST2D::GetResidual(double *sol, double *rhs, double *res)
{

   double *defect = new double[N_TotalDof];
   
   memset(defect, 0, N_TotalDof*SizeOfDouble);
   memset(res, 0, N_TotalDof*SizeOfDouble);
   
   // Block A
   MatVect(SqmatrixA11, sol, defect);
   Daxpy(N_U, 1., defect, res);
   MatVectActive(SqmatrixA12, sol+N_U, defect);
   Daxpy(N_Active_U, 1., defect, res);
   MatVectActive(SqmatrixA21, sol, defect+N_U);
   Daxpy(N_Active_U, 1., defect+N_U, res+N_U);
   MatVect(SqmatrixA22, sol+N_U, defect+N_U);
   Daxpy(N_U, 1., defect+N_U, res+N_U);

    
   // Block G
   MatVect(SqmatrixG11, sol+2*N_U, defect+2*N_U);
   Daxpy(N_S, 1, defect+2*N_U, res+2*N_U);
   MatVectActive(SqmatrixG12, sol+2*N_U+N_S, defect+2*N_U);
   Daxpy(N_Active_S, 1, defect+2*N_U, res+2*N_U);
   MatVectActive(SqmatrixG21, sol+2*N_U, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVect(SqmatrixG22, sol+2*N_U+N_S, defect+2*N_U+N_S);
   Daxpy(N_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVectActive(SqmatrixG23, sol+2*N_U+2*N_S, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVectActive(SqmatrixG32, sol+2*N_U+N_S, defect+2*N_U+2*N_S);
   Daxpy(N_Active_S, 1, defect+2*N_U+2*N_S, res+2*N_U+2*N_S);
   MatVect(SqmatrixG33, sol+2*N_U+2*N_S, defect+2*N_U+2*N_S);
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
   Daxpy(N_Active_S, 1., defect+2*N_U, res+2*N_U);
   MatVect1(MatrixD12, sol+N_U, defect+2*N_U);
   Daxpy(N_Active_S, 1., defect+2*N_U, res+2*N_U);
   MatVect1(MatrixD21, sol, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1., defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVect1(MatrixD22, sol+N_U, defect+2*N_U+N_S);
   Daxpy(N_Active_S, 1., defect+2*N_U+N_S, res+2*N_U+N_S);
   MatVect1(MatrixD31, sol, defect+2*N_U+2*N_S);
   Daxpy(N_Active_S, 1., defect+2*N_U+2*N_S, res+2*N_U+2*N_S);
   MatVect1(MatrixD32, sol+N_U, defect+2*N_U+2*N_S);
   Daxpy(N_Active_S, 1., defect+2*N_U+2*N_S, res+2*N_U+2*N_S);   
      
        if (Disctype_NSECST == 1)
     {
       // Block H
       MatVect(SqmatrixH11, sol+2*N_U+3*N_S, defect+2*N_U+3*N_S);
       Daxpy(N_D, 1., defect+2*N_U+3*N_S, res+2*N_U+3*N_S);
       MatVect(SqmatrixH22, sol+2*N_U+3*N_S+N_D, defect+2*N_U+3*N_S+N_D);
       Daxpy(N_D, 1., defect+2*N_U+3*N_S+N_D, res+2*N_U+3*N_S+N_D);
       MatVect(SqmatrixH33, sol+2*N_U+3*N_S+2*N_D, defect+2*N_U+3*N_S+2*N_D);
       Daxpy(N_D, 1., defect+2*N_U+3*N_S+2*N_D, res+2*N_U+3*N_S+2*N_D);
       
       // Block E
       MatVect1(MatrixE11, sol+2*N_U+3*N_S, defect);
       Daxpy(N_Active_U, 1., defect, res);
       MatVect1(MatrixE12, sol+2*N_U+3*N_S+N_D, defect);
       Daxpy(N_Active_U, 1., defect, res);
       MatVect1(MatrixE22, sol+2*N_U+3*N_S+N_D, defect+N_U);
       Daxpy(N_Active_U, 1., defect+N_U, res+N_U);
       MatVect1(MatrixE23, sol+2*N_U+3*N_S+2*N_D, defect+N_U);
       Daxpy(N_Active_U, 1., defect+N_U, res+N_U);
       
       // Block J
       MatVect1(MatrixJ11, sol, defect+2*N_U+3*N_S);
       Daxpy(N_Active_D, 1., defect+2*N_U+3*N_S, res+2*N_U+3*N_S);
       MatVect1(MatrixJ21, sol, defect+2*N_U+3*N_S+N_D);
       Daxpy(N_Active_D, 1., defect+2*N_U+3*N_S+N_D, res+2*N_U+3*N_S+N_D);
       MatVect1(MatrixJ22, sol+N_U, defect+2*N_U+3*N_S+N_D);
       Daxpy(N_Active_D, 1., defect+2*N_U+3*N_S+N_D, res+2*N_U+3*N_S+N_D);
       MatVect1(MatrixJ32, sol+N_U, defect+2*N_U+3*N_S+2*N_D);
       Daxpy(N_Active_D, 1., defect+2*N_U+3*N_S+2*N_D, res+2*N_U+3*N_S+2*N_D);
       
      // Block BT
      MatVect1(MatrixB1T, sol+2*N_U+3*N_S+3*N_D, defect);
      Daxpy(N_Active_U, 1., defect, res); 
      MatVect1(MatrixB2T, sol+2*N_U+3*N_S+3*N_D, defect+N_U);
      Daxpy(N_Active_U, 1., defect+N_U, res+N_U); 
      
      // Block B
      MatVect1(MatrixB1, sol, defect+2*N_U+3*N_S+3*N_D);
      Daxpy(N_P, 1., defect+2*N_U+3*N_S+3*N_D, res+2*N_U+3*N_S+3*N_D);
      MatVect1(MatrixB2, sol+N_U, defect+2*N_U+3*N_S+3*N_D);
      Daxpy(N_P, 1., defect+2*N_U+3*N_S+3*N_D, res+2*N_U+3*N_S+3*N_D); 
      }
     else if (Disctype_NSECST == 2)
     {
      // Block BT
      MatVect1(MatrixB1T, sol+2*N_U+3*N_S, defect);
      Daxpy(N_Active_U, 1., defect, res); 
      MatVect1(MatrixB2T, sol+2*N_U+3*N_S, defect+N_U);
      Daxpy(N_Active_U, 1., defect+N_U, res+N_U); 
      
      // Block B
      MatVect1(MatrixB1, sol, defect+2*N_U+3*N_S);
      Daxpy(N_P, 1., defect+2*N_U+3*N_S, res+2*N_U+3*N_S); 
      MatVect1(MatrixB2, sol+N_U, defect+2*N_U+3*N_S);
      Daxpy(N_P, 1., defect+2*N_U+3*N_S, res+2*N_U+3*N_S);       
     }

     Daxpy(N_TotalDof, -1., rhs, res);

}  // void TSystemNSECST2D::GetResidual

void TSystemNSECST2D::MeasureErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3,
				    double *u_error, double *p_error, double *s_error)
{
  double errors[4];
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
      
      // errors in pressure
     FeFct[2]->GetErrors(ExactP, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, auxerror_NSE, 1, FeSpace+1, errors);     
     p_error[0] = errors[0];
     p_error[1] = errors[1]; 

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

     
} // void TSystemNSECST2D::MeasureErrors
    

#endif // #ifdef __2D__
