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
* @date      21.03.17
* @History    
* 
 ************************************************************************  */

#ifdef __2D__

#include <Database.h>
#include <SystemNSECST2D.h>
#include <SystemTNSECST2D.h>
#include <SystemTNSECST2D_ALE.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <FreeSurface2D.h>
#include <IsoInterfaceJoint.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
#include <MainUtilities.h>

#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif


TSystemTNSECST2D_ALE::TSystemTNSECST2D_ALE(int N_FESpaces, TFESpace2D **FE_Spaces, int N_FEFunctions, 
				 TFEFunction2D **FeFunctions, int disctype_NSECST,  int solver,
				 TFESpace2D *gridFESpace, TFEVectFunct2D *MeshVelocity ) : 
				 TSystemTNSECST2D(N_FESpaces, FE_Spaces, N_FEFunctions, FeFunctions, disctype_NSECST,  solver)
{  
  char WString[] = "w";  
  TFESpace2D *fesp[1]; 
  
  GridFESpace = gridFESpace;
  N_GridDOFs = gridFESpace->GetN_DegreesOfFreedom();
  
  N_GridActive = gridFESpace->GetActiveBound();
  FeSpace[3] = GridFESpace;
  // grid 
  SquareStructureGrid= new TSquareStructure2D(GridFESpace); 
  SquareStructureGrid->Sort();
   
  // for mesh
  SqmatrixGrid11 = new TSquareMatrix2D(SquareStructureGrid); // Grid11
  SqmatrixGrid12 = new TSquareMatrix2D(SquareStructureGrid); // Grid12
  SqmatrixGrid21 = new TSquareMatrix2D(SquareStructureGrid); // Grid21
  SqmatrixGrid22 = new TSquareMatrix2D(SquareStructureGrid); // Grid22
  
  if(TDatabase::ParamDB->TWO_PHASE_FLOW==1 || TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
  {
   SqmatrixF11 = new TSquareMatrix2D(sqstructureA);
   SqmatrixF22 = new TSquareMatrix2D(sqstructureA);
  }

  SQMATRICES_GRID[0] = SqmatrixGrid11;
  SQMATRICES_GRID[1] = SqmatrixGrid12;
  SQMATRICES_GRID[2] = SqmatrixGrid21;
  SQMATRICES_GRID[3] = SqmatrixGrid22;

  Entries[0] = SqmatrixGrid11->GetEntries();
  Entries[1] = SqmatrixGrid12->GetEntries();
  Entries[2] = SqmatrixGrid21->GetEntries();
  Entries[3] = SqmatrixGrid22->GetEntries();

  GridKCol = SquareStructureGrid->GetKCol();
  GridRowPtr = SquareStructureGrid->GetRowPtr();
     
  fesp[0] = GridFESpace;
  Meshaux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  

  MeshVectFct = MeshVelocity;
  
  MeshVeloFct[0] = MeshVelocity->GetComponent(0);
  MeshVeloFct[1] = MeshVelocity->GetComponent(1);
  MeshVelo =  MeshVelocity->GetValues();
  
  gridpos = new double[2*N_GridDOFs];
  gridpos_old = new double[2*N_GridDOFs];   
  gridpos_ref = new double[2*N_GridDOFs];
  gridpos_aux = new double[2*N_GridDOFs];
  griddisp = new double[2*N_GridDOFs];   
  
 
   memset(gridpos, 0, 2*N_GridDOFs*SizeOfDouble);
   GridPos = new TFEVectFunct2D(GridFESpace, WString, WString, gridpos, N_GridDOFs, 2);  
   GridPos->GridToData();
   
   GridRhs = new double[2*N_GridDOFs];
   
   RefGridPos = new TFEVectFunct2D(GridFESpace, WString, WString, gridpos_ref, N_GridDOFs, 2);
   AuxGridPos = new TFEVectFunct2D(GridFESpace, WString, WString, gridpos_aux, N_GridDOFs, 2); 
 
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
   memcpy(gridpos_ref, gridpos, 2*N_GridDOFs*SizeOfDouble);
   memcpy(gridpos_aux, gridpos, 2*N_GridDOFs*SizeOfDouble);
   
   Aux_ALE = NULL;
   SolveLinearElastic = TRUE;
   
    tmp_Gd = new double[2*N_GridDOFs];
    tmp_Gsol = new double[2*N_GridDOFs];

}

TSystemTNSECST2D_ALE::~TSystemTNSECST2D_ALE()
{
  delete RefGridPos; delete GridPos;  delete AuxGridPos;
  delete [] GridRhs;
  delete SqmatrixGrid11; delete SqmatrixGrid12; delete SqmatrixGrid21; delete SqmatrixGrid22;
    if(TDatabase::ParamDB->TWO_PHASE_FLOW==1 || TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
  {
  delete SqmatrixF11; delete SqmatrixF22;
  }
  delete SquareStructureGrid;
  delete Meshaux;
  delete [] gridpos; delete [] gridpos_old; delete [] gridpos_ref; delete [] griddisp;delete [] gridpos_aux;
  delete [] tmp_Gd; delete [] tmp_Gsol;
}

void TSystemTNSECST2D_ALE::Init(CoeffFct2D *LinCoeffs_Monolithic, CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D **BoundCond, BoundValueFunct2D **BoundValue,  
			   TAuxParam2D *TNSECSTaux, TAuxParam2D *TNSEauxerror, TAuxParam2D *TCSTauxerror,
			      BoundCondFunct2D *GridBoundCond, BoundValueFunct2D *gridBoundValue)
{ int i;
  TDiscreteForm2D *DiscreteFormGrid, *DiscreteFormGalerkinSUPG, *DiscreteFormLPS, *DiscreteFormNLGalerkinSUPG, *DiscreteFormNLLPS, *DiscreteFormRHSGalerkinSUPG, *DiscreteFormRHSLPS;
  
  for (i=0;i<8;i++) 
  { 
  // save the boundary condition
  BoundaryConditions[i] = BoundCond[i];
  // save the boundary values  
  BoundaryValues[i] = BoundValue[i];
  }
  
  GridBoundaryConditions[0] = GridBoundCond;
  GridBoundValue[0] = gridBoundValue;
  
  // save the bilinear coefficients
  LinCoeffs[0] = LinCoeffs_Monolithic;
  
  // aux for assembling and error calculations
  aux = TNSECSTaux;
  auxerror_NSE = TNSEauxerror;
  auxerror_CST = TCSTauxerror;
  
  if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
{
      // set the Discreteforms
   InitializeDiscreteForms_2PhaseAxial3D_TNSECST(DiscreteFormLPS, DiscreteFormNLLPS,
				   DiscreteFormRHSLPS, 
				   DiscreteFormGrid, 
				   LinCoeffs[0], GridBilinearCoeffs);
//     cout<<"Correct loop entered !!!!!!!!!!!!! \n";
 
}
else if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
{
    InitializeDiscreteForms_ImpDropAxial3D_TNSECST(DiscreteFormLPS, DiscreteFormNLLPS,
				 DiscreteFormRHSLPS, DiscreteFormGrid, 
				 LinCoeffs[0], GridBilinearCoeffs);
        cout<<"Correct loop entered for initialization of discrete forms of imp drop !!!!!!!!!!!!! \n";
	
}
else
{
    OutPut("Implemented only for two-phase and free surface flows \n");
    OutPut("Change TWO_PHASE_FLOW or FREE_SURFACE_FLOW  to 1 in dat file !!!!! " << endl);
    exit(4711);
  
}

  
     if (Disctype_NSECST == 1)
     {  
      OutPut("Implemented only using LPS \n");
      OutPut("Change Disctype_NSECST to 2  !!!!! " << endl);
      exit(4711);
     }
    else if (Disctype_NSECST == 2)
    {
      DiscreteFormARhs = DiscreteFormLPS;
      DiscreteFormNL   = DiscreteFormNLLPS;
      DiscreteFormRhs  = DiscreteFormRHSLPS;
      DiscreteFormMesh = DiscreteFormGrid;
    }
    else 
    {
      cout<<"Invalid Disctype_NSECST !!!\n";
      exit(0);
    }
  
     // LPS assembling terms
        if (Disctype_NSECST == 2)
      {
	SqmatrixA11_LPS->Reset();
	SqmatrixA12_LPS->Reset();
	SqmatrixA21_LPS->Reset();
	SqmatrixA22_LPS->Reset();
	
	SqmatrixG11_LPS->Reset();
	SqmatrixG12_LPS->Reset();
	SqmatrixG21_LPS->Reset();
	SqmatrixG22_LPS->Reset();
	SqmatrixG23_LPS->Reset();
	SqmatrixG32_LPS->Reset();
	SqmatrixG33_LPS->Reset();
	
	SqmatrixG11_LPS_Streamline->Reset();
	SqmatrixG22_LPS_Streamline->Reset();
	SqmatrixG33_LPS_Streamline->Reset();
	
        clock_t Time_Assemble = clock();
       	
	if(TDatabase::ParamDB->Axial3D=1)
	{
	AddDeformationTensorTerm_2PhaseOrImpDropFlow_3DAxial(SqmatrixA11_LPS, SqmatrixA12_LPS, SqmatrixA21_LPS, SqmatrixA22_LPS, 1.0, 1);
	
	if (TDatabase::ParamDB->LP_FULL_GRADIENT == 1 && TDatabase::ParamDB->LP_DIVERGENCE == 1)
        AddTauTerm_2PhaseOrImpDropFlow_3DAxial(SqmatrixG11_LPS,SqmatrixG12_LPS,SqmatrixG21_LPS,SqmatrixG22_LPS,SqmatrixG23_LPS,SqmatrixG32_LPS,SqmatrixG33_LPS,TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF,TDatabase::ParamDB->LP_DIVERGENCE_COEFF,1.0,1); 
	else if(TDatabase::ParamDB->LP_STREAMLINE == 1 && TDatabase::ParamDB->LP_DIVERGENCE == 1)
	{
	  AddStreamlineTerm_2PhaseOrImpDropFlow_3DAxial(SqmatrixG11_LPS_Streamline, SqmatrixG22_LPS_Streamline, SqmatrixG33_LPS_Streamline, FeFct[0], FeFct[1], TDatabase::ParamDB->LP_STREAMLINE_COEFF, 1.0, 1);
	  AddTauTerm_2PhaseOrImpDropFlow_3DAxial(SqmatrixG11_LPS,SqmatrixG12_LPS,SqmatrixG21_LPS,SqmatrixG22_LPS,SqmatrixG23_LPS,SqmatrixG32_LPS,SqmatrixG33_LPS,0.0,TDatabase::ParamDB->LP_DIVERGENCE_COEFF,1.0,1); 
	}
	else
	{
	  cout<<"Invalid LPS type !!!!!!! \n"; exit(0);
	}
	}
	else
	{
	AddDeformationTensorTerm_2PhaseOrImpDropFlow(SqmatrixA11_LPS, SqmatrixA12_LPS, SqmatrixA21_LPS, SqmatrixA22_LPS, 1.0, 1);
        AddTauTerm_2PhaseOrImpDropFlow(SqmatrixG11_LPS,SqmatrixG12_LPS,SqmatrixG21_LPS,SqmatrixG22_LPS,SqmatrixG23_LPS,SqmatrixG32_LPS,SqmatrixG33_LPS,TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF,TDatabase::ParamDB->LP_DIVERGENCE,1.0,1); 
	}
	 cout << "LPS Assemble time : "<<double( clock() - Time_Assemble ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
       }  
       
       
  // Matrix initialization due to SMPI (Otherwise not required)
   int N_SquareMatrices, N_RectMatrices;

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
   

  
   #ifdef _SMPI
    if(SOLVER == DIRECT)
     {
      if(Disctype_NSECST == 2)
      {
      P_DS = new TSeqParDirectSolver(N_U,N_P,N_S,0,SQMATRICES,MATRICES);
      }

     }
   #endif
       
  
} // void TSystemTNSECST2D_ALE::Init


void TSystemTNSECST2D_ALE::AssembleMeshMat()
{
  
 TFESpace2D *fesp[1];

   fesp[0] = GridFESpace; 
   
   SQMATRICES_GRID[0]->Reset();
   SQMATRICES_GRID[1]->Reset();
   SQMATRICES_GRID[2]->Reset();
   SQMATRICES_GRID[3]->Reset();    

     Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormMesh,
             GridBoundaryConditions,
             GridBoundValue,
             Meshaux);

  // for Dirichlet rows in off-diagonal matrices
  memset(Entries[1] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);     
  
} //AssembleMeshMat
  

 
 void TSystemTNSECST2D_ALE::Assemble(double *sol, double *rhs, double dt)
 {
     double *RHSs[8], Params[10];
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
   
   SqmatrixMU12->Reset();
   SqmatrixMU21->Reset();
   SqmatrixMS12->Reset();
   SqmatrixMS21->Reset();
   SqmatrixMS23->Reset();
   SqmatrixMS32->Reset();
     
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
      Assemble2D(N_FeSpace+1, FeSpace,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormARhs,
        BoundaryConditions,
        BoundaryValues,
        aux);
      
        if(TDatabase::ParamDB->TWO_PHASE_FLOW == 1)
      {
        SqmatrixF11->Reset();
        SqmatrixF22->Reset();
        Interface_2PhaseSurfAxial3D(SqmatrixF11, SqmatrixF22, rhs, rhs + N_U, BoundaryConditions[0], dt);
//        Adding interface entries to A11 and A22
      MatAdd(SqmatrixA11, SqmatrixF11, 1.);
      MatAdd(SqmatrixA22, SqmatrixF22, 1.);
      }   
       else if(TDatabase::ParamDB->FREE_SURFACE_FLOW == 1)
      {
        SqmatrixF11->Reset();
        SqmatrixF22->Reset();
	FreeSurf_axial3D_new(SqmatrixF11, SqmatrixF22,  rhs, rhs + N_U, BoundaryConditions[0], dt,
                          FeFct[0]->GetValues(), NULL, Params);
//        Adding free surface entries to A11 and A22
      MatAdd(SqmatrixA11, SqmatrixF11, 1.);
      MatAdd(SqmatrixA22, SqmatrixF22, 1.);
      }  

      // Add the LPS terms
      if (Disctype_NSECST == 2)
     {
        MatAdd(SqmatrixA11, SqmatrixA11_LPS,1.);
	MatAdd(SqmatrixA12, SqmatrixA12_LPS,1.);
	MatAdd(SqmatrixA21, SqmatrixA21_LPS,1.);
	MatAdd(SqmatrixA22, SqmatrixA22_LPS,1.);
	
	MatAdd(SqmatrixG11,  SqmatrixG11_LPS,1.);
	MatAdd(SqmatrixG12,  SqmatrixG12_LPS,1.);
	MatAdd(SqmatrixG21,  SqmatrixG21_LPS,1.);
	MatAdd(SqmatrixG22,  SqmatrixG22_LPS,1.);
	MatAdd(SqmatrixG23, SqmatrixG23_LPS,1.);
	MatAdd(SqmatrixG32, SqmatrixG32_LPS,1.);
	MatAdd(SqmatrixG33, SqmatrixG33_LPS,1.);
	
	if(TDatabase::ParamDB->LP_STREAMLINE == 1)
	{
	  SqmatrixG11_LPS_Streamline->Reset();
	  SqmatrixG22_LPS_Streamline->Reset();
	  SqmatrixG33_LPS_Streamline->Reset();
	  
	  AddStreamlineTerm_2PhaseOrImpDropFlow_3DAxial(SqmatrixG11_LPS_Streamline, SqmatrixG22_LPS_Streamline, SqmatrixG33_LPS_Streamline, FeFct[0], FeFct[1], TDatabase::ParamDB->LP_STREAMLINE_COEFF, 1.0, 1); 
	  
	  MatAdd(SqmatrixG11,  SqmatrixG11_LPS_Streamline,1.);
	  MatAdd(SqmatrixG22,  SqmatrixG22_LPS_Streamline,1.);
	  MatAdd(SqmatrixG33,  SqmatrixG33_LPS_Streamline,1.);
	}
     }
    
         int N_RecMat_Slip;
        // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
	 double *RHSs_Slip[2];
         TFESpace2D *fesprhs_slip[2];
	
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

      N_Rhs = 2;
      RHSs_Slip[0] = rhs;
      RHSs_Slip[1] = rhs +   N_U;
     
      fesprhs_slip[0] = FeSpace[0];
      fesprhs_slip[1] = FeSpace[0];
            
        Assemble2DSlipBC(1, FeSpace,
                         8, SQMATRICES,
                         N_RecMat_Slip, MATRICES,
                         N_Rhs, RHSs_Slip, fesprhs_slip,
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

 } // void TSystemTNSECST2D_ALE::Assemble
  
  
 
 void TSystemTNSECST2D_ALE:: AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol)
 {
  
   double tau, val = TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT, temp;
  
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  memset(B, 0, N_TotalDof*SizeOfDouble);

    if(TDatabase::ParamDB->TWO_PHASE_FLOW == 1 || TDatabase::ParamDB->FREE_SURFACE_FLOW == 1)
   { 
     // since rhs depends on moving grid so its better to use current geo
     Daxpy(N_Active_U, tau, rhs, B);
     Daxpy(N_Active_U, tau, rhs+N_U, B+N_U);
     Daxpy(N_Active_S, tau, rhs+2*N_U, B+2*N_U);
     Daxpy(N_Active_S, tau, rhs+2*N_U+N_S, B+2*N_U+N_S);
     Daxpy(N_Active_S, tau, rhs+2*N_U+2*N_S, B+2*N_U+2*N_S);
   }
    else
  {
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
  }
  
  // scale the BT, C, B, matrices, not in nonlinear step 
  if (scale != 1.0)
  {
         Dscal(MatrixB1T->GetN_Entries(), scale, MatrixB1T->GetEntries());
         Dscal(MatrixB2T->GetN_Entries(), scale, MatrixB2T->GetEntries()); 
       
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
       

       temp = - TDatabase::TimeDB->THETA2;
       MatVect1(MatrixC11, sol+2*N_U, defect);
       Daxpy(N_Active_U, temp, defect, B);
       MatVect1(MatrixC12, sol+2*N_U+N_S, defect);
       Daxpy(N_Active_U, temp, defect, B);
       MatVect1(MatrixC22, sol+2*N_U+N_S, defect+N_U);
       Daxpy(N_Active_U, temp, defect+N_U, B+N_U);
       MatVect1(MatrixC23, sol+2*N_U+2*N_S, defect+N_U);
       Daxpy(N_Active_U, temp, defect+N_U, B+N_U);
       
         temp = tau*TDatabase::TimeDB->THETA1;
	 Dscal(MatrixC11->GetN_Entries(), temp, MatrixC11->GetEntries());
	 Dscal(MatrixC12->GetN_Entries(), temp, MatrixC12->GetEntries());
	 Dscal(MatrixC22->GetN_Entries(), temp, MatrixC22->GetEntries());
	 Dscal(MatrixC23->GetN_Entries(), temp, MatrixC23->GetEntries());
      
     
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

   
   SystMatAssembled  = TRUE;

 }  // void TSystemTNSECST2D_ALE:: AssembleSystMat
  
  
    
  void TSystemTNSECST2D_ALE::GetTNSECSTResidual(double *sol, double *res)
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
//    cout<<"MU11 : "<<Ddot(N_U, defect, defect)<<"\n";
   Daxpy(N_U, 1, defect, res);
//    cout<<"MU11 : "<<Ddot(N_U, res, res)<<"\n";
   MatVectActive(SqmatrixMU12, sol+N_U, defect);
//    cout<<"MU12 : "<<Ddot(N_U, defect, defect)<<"\n";
   Daxpy(N_Active_U, 1, defect, res);
//    cout<<"MU12 : "<<Ddot(N_U, res, res)<<"\n";
   MatVectActive(SqmatrixMU21, sol, defect+N_U);
//    cout<<"MU21 : "<<Ddot(N_U, defect+N_U, defect+N_U)<<"\n";
   Daxpy(N_Active_U, 1, defect+N_U, res+N_U);
//    cout<<"MU21 : "<<Ddot(N_U, res+N_U, res+N_U)<<"\n";
   MatVect(SqmatrixMU22, sol+N_U, defect+N_U);
//    cout<<"MU22 : "<<Ddot(N_U, defect+N_U, defect+N_U)<<"\n";
   Daxpy(N_U, 1, defect+N_U, res+N_U);
//    cout<<"MU22 : "<<Ddot(N_U, res+N_U, res+N_U)<<"\n";
   
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

      // Block BT
      MatVect1(MatrixB1T, sol+2*N_U+3*N_S, defect);
//       cout<<"B1T : "<<Ddot(N_U, defect, defect)<<"\n";
      Daxpy(N_Active_U, 1, defect, res); 
//       cout<<"B1T : "<<Ddot(N_U, res, res)<<"\n";
      MatVect1(MatrixB2T, sol+2*N_U+3*N_S, defect+N_U);
//       cout<<"B2T : "<<Ddot(N_U, defect+N_U, defect+N_U)<<"\n";
      Daxpy(N_Active_U, 1, defect+N_U, res+N_U);
//       cout<<"B2T : "<<Ddot(N_U, res+N_U, res+N_U)<<"\n";
      
      // Block B
      MatVect1(MatrixB1, sol, defect+2*N_U+3*N_S);
      Daxpy(N_P, 1, defect+2*N_U+3*N_S, res+2*N_U+3*N_S); 
      MatVect1(MatrixB2, sol+N_U, defect+2*N_U+3*N_S);
      Daxpy(N_P, 1, defect+2*N_U+3*N_S, res+2*N_U+3*N_S);       

//       cout<<"Before rhs : "<<Ddot(2*N_U, res, res)<<"\n";
//              cout<<"B1 : "<<Ddot(N_U, B, B)<<"\n";
//        cout<<"B2 : "<<Ddot(N_U, B+N_U, B+N_U)<<"\n";
     Daxpy(N_TotalDof, -1., B, res);
//      cout<<"After rhs : "<<Ddot(2*N_U, res, res)<<"\n";

}
  
  void TSystemTNSECST2D_ALE::Solve(double *sol)
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
	     
		     
      break;
 
      default:

            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    

} // void TSystemTNSECST2D_ALE::Solve
 
   void TSystemTNSECST2D_ALE::RestoreMassMat()
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
   
  
} //void TSystemTNSECST2D_ALE::RestoreMassMat()
  
  void TSystemTNSECST2D_ALE::AssembleNonLinear(double *sol, double *rhs)
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
      Assemble2D(N_FeSpace+1, FeSpace,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormNL,
        BoundaryConditions+5,
        BoundaryValues+5,
        aux);  
      
       if(TDatabase::ParamDB->TWO_PHASE_FLOW == 1 || TDatabase::ParamDB->FREE_SURFACE_FLOW == 1)
      {
       // Adding interface entries to A11 and A22
      MatAdd(SqmatrixA11, SqmatrixF11, 1.);
      MatAdd(SqmatrixA22, SqmatrixF22, 1.);
      }
      
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
	
          if(TDatabase::ParamDB->LP_STREAMLINE == 1)
	{
	  SqmatrixG11_LPS_Streamline->Reset();
	  SqmatrixG22_LPS_Streamline->Reset();
	  SqmatrixG33_LPS_Streamline->Reset();
	  
	  AddStreamlineTerm_2PhaseOrImpDropFlow_3DAxial(SqmatrixG11_LPS_Streamline, SqmatrixG22_LPS_Streamline, SqmatrixG33_LPS_Streamline, FeFct[0], FeFct[1], TDatabase::ParamDB->LP_STREAMLINE_COEFF, 1.0, 1); 
	  
	  MatAdd(SQMATRICES[2],  SqmatrixG11_LPS_Streamline,1.);
	  MatAdd(SQMATRICES[5],  SqmatrixG22_LPS_Streamline,1.);
	  MatAdd(SQMATRICES[8],  SqmatrixG33_LPS_Streamline,1.);
	}
     }
      
      int N_RecMat_Slip;
        // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
	RHSs[0] = rhs;
	RHSs[1] = rhs + N_U;
	
	fesprhs[0] = FeSpace[0];
	fesprhs[1] = FeSpace[0];	
	
        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
	
        Assemble2DSlipBC(1, FeSpace,
                         2, SQMATRICES,
                         0, NULL,
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
   
 }
 
  /* assemble only LHS, and non-linear term in rhs */
  void TSystemTNSECST2D_ALE::AssembleSystMatNonLinear()
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
 
 
// for free-surface flows (impinging drop)
void TSystemTNSECST2D_ALE::GetMeshVelo(double tau, TVertex ***MovBoundVert, int *N_MovVert, 
				    TBaseCell **Free_Cells, int **IsoCellEdgeNos, bool &reparam, 
				    TFEVectFunct2D **FEVectFuncts)
{
 
    GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, GridRhs,
                        GridKCol, GridRowPtr,
                        GridPos, AuxGridPos,
                        FEVectFuncts[0], tau,
                        MeshVectFct, MovBoundVert, N_MovVert,
                        Free_Cells, IsoCellEdgeNos, reparam, RefGridPos, FEVectFuncts[1]);
  
}

// for free-surface flows (impinging drop)
void TSystemTNSECST2D_ALE::MoveMesh(double tau, TVertex ***MovBoundVert, int *N_MovVert, 
				 TBaseCell **Free_Cells, int **IsoCellEdgeNos, bool &reparam, 
				 int N_ReParam, TFEVectFunct2D **FEVectFuncts)
{
       MoveGrid_imping(Entries, tmp_Gsol, tmp_Gd, GridRhs,
                  GridKCol, GridRowPtr,
                  GridPos, FEVectFuncts[0], tau,
                  AuxGridPos, 
                  MovBoundVert, N_MovVert,
                  Free_Cells, IsoCellEdgeNos, reparam,
                  N_ReParam, FEVectFuncts[1]);
}
 
 
 
 

#endif   // #ifdef __2D__
