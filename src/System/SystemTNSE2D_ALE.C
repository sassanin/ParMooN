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
* @brief     source file for TSystemTNSE2D_ALE
* @author    Sashikumaar Ganesan
* @date      23.02.15

 ************************************************************************  */
#ifdef __2D__

#include <Database.h>
#include <SystemTNSE2D_ALE.h>
// #include <SquareStructure2D.h>
// #include <SquareMatrix2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <FreeSurface2D.h>
#include <MainUtilities.h>
#include <Upwind.h>
#include <IsoInterfaceJoint.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include<Constants.h>
#include <VMS.h>

                                   
TSystemTNSE2D_ALE::TSystemTNSE2D_ALE(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                                           TFEFunction2D *p, double *sol, double *rhs,  int disctype, int nsetype, int solver,
#ifdef __PRIVATE__                                            
                                     TFESpace2D *Projection_space, TFESpace2D *Stress_FeSpace,   TFESpace2D *Deformation_FeSpace,
#endif   
                                           TFESpace2D *gridFESpace, TFEVectFunct2D *MeshVelocity, bool conservativeale)
                                           : TSystemTNSE2D(velocity_fespace, presssure_fespace, Velocity, p, sol, rhs, disctype, nsetype, solver
#ifdef __PRIVATE__  
                                        ,Projection_space, Stress_FeSpace, Deformation_FeSpace
#endif   
                                     )
{
  char WString[] = "w";  
  TFESpace2D *fesp[1]; 
  
//   /** old M mass matrix */
//   OldSqmatrixM11 = new TSquareMatrix2D(sqstructureA);
//   OldSqmatrixM12 = new TSquareMatrix2D(sqstructureA);
//   OldSqmatrixM21 = new TSquareMatrix2D(sqstructureA); 
//   OldSqmatrixM22 = new TSquareMatrix2D(sqstructureA);
  
  GridFESpace = gridFESpace;
  N_GridDOFs = gridFESpace->GetN_DegreesOfFreedom();
  N_GridActive = gridFESpace->GetActiveBound();
  FeSpaces[2] = GridFESpace;
  // grid 
  SquareStructureG= new TSquareStructure2D(GridFESpace); 
  SquareStructureG->Sort();
   
  // for mesh
  SqmatrixG11 = new TSquareMatrix2D(SquareStructureG); // G11
  SqmatrixG12 = new TSquareMatrix2D(SquareStructureG); // G12
  SqmatrixG21 = new TSquareMatrix2D(SquareStructureG); // G21
  SqmatrixG22 = new TSquareMatrix2D(SquareStructureG); // G22
  
  if(TDatabase::ParamDB->TWO_PHASE_FLOW==1 || TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
  {
   SqmatrixF11 = new TSquareMatrix2D(sqstructureA);
   SqmatrixF22 = new TSquareMatrix2D(sqstructureA);
  }
  
  SQMATRICES_GRID[0] = SqmatrixG11;
  SQMATRICES_GRID[1] = SqmatrixG12;
  SQMATRICES_GRID[2] = SqmatrixG21;
  SQMATRICES_GRID[3] = SqmatrixG22;
  
  Entries[0] = SqmatrixG11->GetEntries();
  Entries[1] = SqmatrixG12->GetEntries();
  Entries[2] = SqmatrixG21->GetEntries();
  Entries[3] = SqmatrixG22->GetEntries();

  GridKCol = SquareStructureG->GetKCol();
  GridRowPtr = SquareStructureG->GetRowPtr();
     
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
   CONSERVATIVEALE = conservativeale;
   
    tmp_Gd = new double[2*N_GridDOFs];
    tmp_Gsol = new double[2*N_GridDOFs];
   
}

TSystemTNSE2D_ALE::~TSystemTNSE2D_ALE()
{
  delete RefGridPos; delete GridPos; delete AuxGridPos;
  delete [] GridRhs;
  delete SqmatrixG11; delete SqmatrixG12; delete SqmatrixG21; delete SqmatrixG22;
  
  if(TDatabase::ParamDB->TWO_PHASE_FLOW==1 || TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
  {
  delete SqmatrixF11; delete SqmatrixF22;
  }
  
  delete SquareStructureG;
  delete Meshaux;
  delete [] gridpos; delete [] gridpos_old; delete [] gridpos_ref; delete [] griddisp; delete [] gridpos_aux;
  delete [] tmp_Gd; delete [] tmp_Gsol;
}


void TSystemTNSE2D_ALE::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
                                CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D *GridBoundCond, BoundValueFunct2D *gridBoundValue,
                                TAuxParam2D *aux, TAuxParam2D *nseaux_error)
{
//  TFESpace2D *fesp[2];  
//  TFEFunction2D *fefct[4]; 
  
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormColetti;
  TDiscreteForm2D *DiscreteFormGL00Convolution;
  TDiscreteForm2D *DiscreteFormGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMSProjection;

  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLColetti;
  TDiscreteForm2D *DiscreteFormNLGL00Convolution;
  TDiscreteForm2D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;
  TDiscreteForm2D *DiscreteFormNLVMSProjection;

  TDiscreteForm2D *DiscreteFormRHS;
  TDiscreteForm2D *DiscreteFormRHSColetti;
  TDiscreteForm2D *DiscreteFormRHSSmagorinskyExpl;
  TDiscreteForm2D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm2D *DiscreteFormRHSLESModel;
  TDiscreteForm2D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm2D *DiscreteFormMatrixAuxProblemU;
  
  // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  
  GridBoundaryConditions[0] = GridBoundCond;
  
  // save the boundary values  
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
  GridBoundValue[0] = gridBoundValue; 
  
  // save the nse bilinear coefficient   
  LinCoeffs[0] = lincoeffs;

//   default, i.e., velocity for nonlinear term
  NSEaux = aux;
// aux for calculating the error
  NSEaux_error = nseaux_error;

if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
{
InitializeDiscreteForms_2PhaseAxial3D(DiscreteFormGalerkin,
                                  DiscreteFormNLGalerkin,
                                  DiscreteFormRHS,
                                  DiscreteFormGrid, LinCoeffs[0], GridBilinearCoeffs,
                                  TDatabase::ParamDB->NSTYPE);
}
else if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
{
    InitializeDiscreteForms_Moving(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
                                 DiscreteFormGrid, LinCoeffs[0], GridBilinearCoeffs);
}
else
{
   // set the Discreteform for grid
  InitializeDiscreteFormGrid(DiscreteFormGrid, GridBilinearCoeffs);

  // set the Discreteforms
  InitializeDiscreteForms(DiscreteFormGalerkin, DiscreteFormUpwind,
              DiscreteFormSmagorinsky,DiscreteFormColetti,
              DiscreteFormGL00Convolution,DiscreteFormGL00AuxProblem,
              DiscreteFormVMSProjection,
              DiscreteFormNLGalerkin,
              DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky,
              DiscreteFormNLColetti,DiscreteFormNLGL00Convolution,
              DiscreteFormNLGL00AuxProblem,
              DiscreteFormNLVMSProjection,
              DiscreteFormRHS,
              DiscreteFormRHSColetti,
              DiscreteFormRHSLESModel,
              DiscreteFormMatrixGL00AuxProblem,
              DiscreteFormGL00AuxProblemRHS,
              DiscreteFormRHSSmagorinskyExpl,
              DiscreteFormMatrixAuxProblemU,
              DiscreteFormRHSAuxProblemU,
              LinCoeffs[0], NSEType);
}
  
  
    // find discrete form
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
	    if(TDatabase::ParamDB->FREE_SURFACE_FLOW==0)
            DiscreteFormRhs = DiscreteFormRHS;
          break;

          case UPWIND:
          if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
	  { 
	    Error("Not implemented for DISCTYPE UPWIND" << endl);
            exit(-1);
	  }
            DiscreteFormARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;    
            DiscreteFormRhs = DiscreteFormRHS;
            break;

//           case SMAGORINSKY:
//             DiscreteFormARhs = DiscreteFormSmagorinsky;
//             DiscreteFormNL = DiscreteFormNLSmagorinsky;              
//             break;

#ifdef __PRIVATE__ 
          case VMS_PROJECTION:
    
            DiscreteFormARhs = DiscreteFormVMSProjection;
            DiscreteFormNL = DiscreteFormNLVMSProjection;
            DiscreteFormRhs = DiscreteFormRHS;    
            break;
#endif
	    
          default:
            Error("Unknown DISCTYPE" << endl);
            exit(-1);
        } 
     
     // set the discrete form for the Stokes equation
      if (TDatabase::ParamDB->PROBLEM_TYPE == 3)
       {
        DiscreteFormARhs = DiscreteFormUpwind;     
        DiscreteFormNL = NULL;
       }
 
} // TSystemTNSE2D_ALE::Init

 
/* Assemble M, A and rhs */ 
void TSystemTNSE2D_ALE::Assemble(double *sol, double *rhs, double dt)
{
  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  
  double *RHSs[3], Params[10];

  TFESpace2D *fesprhs[3];
  
  // all matrices are assembled, so 
  SystMatAssembled = FALSE;

  N_Rhs = 2;
  N_FESpaces = 3;   

     // initialize matrices
     switch(NSEType)
      {
        case 1:
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixM11;
         MATRICES[0] = MatrixB1;
         MATRICES[1] = MatrixB2;

         SQMATRICES[0]->Reset();
         SQMATRICES[1]->Reset();
         MATRICES[0]->Reset();
         MATRICES[1]->Reset();

         N_SquareMatrices = 2;
         N_RectMatrices = 2;

         break;

        case 2:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[0] = SqmatrixM11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 2;
          N_RectMatrices = 4;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          SQMATRICES[4] = SqmatrixM11;
          SQMATRICES[5] = SqmatrixM22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          // assembling all matrices, so reset all
          SqmatrixM12->Reset();
          SqmatrixM21->Reset();

          N_SquareMatrices = 6;
          N_RectMatrices = 2;
	  
#ifdef __PRIVATE__ 
        if(Disctype == VMS_PROJECTION)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] =  MatricesL;
          SQMATRICES[6]->Reset();

          N_RectMatrices = 6;
          MATRICES[2] = Matrices_tilde_G11;
          MATRICES[3] = Matrices_tilde_G22;
          MATRICES[4] = Matrices_G11;
          MATRICES[5] = Matrices_G22;
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();

          N_FESpaces = 4;
        }
  
        break;
#endif
        case 4:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          SQMATRICES[4] = SqmatrixM11;
          SQMATRICES[5] = SqmatrixM22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          // assembling all matrices, so reset all
          SqmatrixM12->Reset();
          SqmatrixM21->Reset();
  
          N_SquareMatrices = 6;
          N_RectMatrices = 4;
#ifdef __PRIVATE__ 
        if(Disctype == VMS_PROJECTION)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] =  MatricesL;
          SQMATRICES[6]->Reset();

          N_RectMatrices = 8;
          MATRICES[4] = Matrices_tilde_G11;
          MATRICES[5] = Matrices_tilde_G22;
          MATRICES[6] = Matrices_G11;
          MATRICES[7] = Matrices_G22;
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();
          MATRICES[6]->Reset();
          MATRICES[7]->Reset();

          N_FESpaces = 4;
        }
#endif        
       break;
      } //  switch(NSEType)

      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
     
      fesprhs[0] = FeSpaces[0];
      fesprhs[1] = FeSpaces[0];
      fesprhs[2] = FeSpaces[1];

      // assemble
      Assemble2D(N_FESpaces, FeSpaces,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormARhs,
        BoundaryConditions,
        BoundaryValues,
        NSEaux);

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
        
      if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == 3) )
       {
        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[3], FeFct[0], FeFct[1]);
           break;
         }                        // endswitch
       }                          // endif     
          
          
          
      // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {

        if(NSEType <4)
         {
          OutPut("For slip with friction bc NSTYPE 4 is ");
          OutPut("necessary !!!!! " << endl);
          exit(4711);
         }
 
        N_FESpaces = 1;
        N_SquareMatrices = 8;
        N_RectMatrices = 2;
        N_Rhs = 2;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
        SQMATRICES[2] = SqmatrixA12;
        SQMATRICES[3] = SqmatrixA21;
        SQMATRICES[4] = SqmatrixM11;
        SQMATRICES[5] = SqmatrixM22;
        SQMATRICES[6] = SqmatrixM12;
        SQMATRICES[7] = SqmatrixM21;

        MATRICES[0] = MatrixB1T;
        MATRICES[1] = MatrixB2T;


        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, MATRICES,
                         N_Rhs, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);


      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=      
    

     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 

#ifdef __PRIVATE__ 
      // update matrices
      if (Disctype == VMS_PROJECTION)
        {
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixA12;
         SQMATRICES[2] = SqmatrixA21;
         SQMATRICES[3] = SqmatrixA22;
         SQMATRICES[6] =  MatricesL;
         MATRICES[2] = Matrices_tilde_G11;
         MATRICES[3] = Matrices_tilde_G22;
         MATRICES[4] = Matrices_G11;
         MATRICES[5] = Matrices_G22;

         VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
        }
#endif
} // TSystemTNSE2D_ALE::Assemble(T


void TSystemTNSE2D_ALE::GetMeshVeloAndMove(int *N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
                                             double * Iso_refX,  double Currtime, double tau)
{
 int i, N_GridBDDOFs;
 
  GridPos->GridToData();   
  memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);  
 
  // modyfy the boundary 
  RefGridPos->DataToGrid();  
  ModifyBoudary(N_MovVert, MovBoundVert, Free_Joint, Iso_refX, Currtime);    
  
  // data with updated BD values
  GridPos->GridToData();  

  N_GridBDDOFs = N_GridDOFs - N_GridActive;
  
  memset(GridRhs, 0, 2*N_GridDOFs*SizeOfDouble);     
  memcpy(GridRhs+N_GridActive, gridpos+N_GridActive, N_GridBDDOFs*SizeOfDouble);     //rhs1  
  memcpy(GridRhs+(N_GridDOFs+N_GridActive), gridpos+(N_GridDOFs+N_GridActive), N_GridBDDOFs*SizeOfDouble);   //rhs2       
    
  Daxpy(N_GridBDDOFs, -1., gridpos_old+N_GridActive, GridRhs+N_GridActive);
  Daxpy(N_GridBDDOFs, -1., gridpos_old+(N_GridDOFs+N_GridActive), GridRhs+(N_GridDOFs+N_GridActive));    
     
  memcpy(MeshVelo, GridRhs, 2*N_GridDOFs*SizeOfDouble);   
    
  SolveGridEquation(Entries, MeshVelo, GridRhs, GridKCol, GridRowPtr, N_GridDOFs);
         
  memcpy(gridpos, gridpos_old, 2*N_GridDOFs*SizeOfDouble);
  Daxpy(2*N_GridDOFs, 1., MeshVelo, gridpos);

  //move the mesh
  GridPos->DataToGrid();   

  // mesh velocity
  Dscal(2*N_GridDOFs, 1./tau, MeshVelo);
//    cout<< "MeshVelo " <<Ddot((2*N_GridDOFs), MeshVelo, MeshVelo)<<endl; 
} // TSystemTNSE2D_ALE::GetMeshVeloAndMove


void TSystemTNSE2D_ALE::GetMeshVeloAndMove(double Currtime, double tau)
{
 int i;
  
  GridPos->GridToData();   
  memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);  

   //  move velo in current time  
   for(i=0;i<N_GridDOFs;i++)
     ModifyCoord(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], Currtime);   

   //compute mesh velocity
   memcpy(MeshVelo, gridpos, 2*N_GridDOFs*SizeOfDouble);     
   Daxpy(2*N_GridDOFs, -1., gridpos_old, MeshVelo);        
   Dscal(2*N_GridDOFs, -1./tau, MeshVelo); // - sign du*/ //e to -w\cdot\nabla C in the equation   
   
   //move the mesh
   GridPos->DataToGrid();   
//    memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble); 
} //TSystemTNSE2D_ALE::GetMeshVeloAndMove(
 
  
// for free-surface flows (impinging drop)
void TSystemTNSE2D_ALE::GetMeshVelo(double tau, TVertex ***MovBoundVert, int *N_MovVert, 
				    TBaseCell **Free_Cells, int **IsoCellEdgeNos, bool &reparam)
{
 
    GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, GridRhs,
                        GridKCol, GridRowPtr,
                        GridPos, AuxGridPos,
                        VelocityFct, tau,
                        MeshVectFct, MovBoundVert, N_MovVert,
                        Free_Cells, IsoCellEdgeNos, reparam, RefGridPos, NULL);

 
}


// for free-surface flows (impinging drop)
void TSystemTNSE2D_ALE::MoveMesh(double tau, TVertex ***MovBoundVert, int *N_MovVert, 
				 TBaseCell **Free_Cells, int **IsoCellEdgeNos, bool &reparam, int N_ReParam)
{
       MoveGrid_imping(Entries, tmp_Gsol, tmp_Gd, GridRhs,
                  GridKCol, GridRowPtr,
                  GridPos, VelocityFct, tau,
                  AuxGridPos, 
                  MovBoundVert, N_MovVert,
                  Free_Cells, IsoCellEdgeNos, reparam,
                  N_ReParam, NULL);
}

  
void TSystemTNSE2D_ALE::AssembleMeshMat()
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
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValue,
             Meshaux);
     
  // for Dirichlet rows in off-diagonal matrices
  memset(Entries[1] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);     
  
} //AssembleMeshMat
  

void TSystemTNSE2D_ALE::AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol)
{
 double tau, val = TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT;
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  memset(B, 0, (2*N_U+N_P)*SizeOfDouble);    
     
        if(TDatabase::ParamDB->TWO_PHASE_FLOW == 1 || TDatabase::ParamDB->FREE_SURFACE_FLOW == 1)
      { 
	// since rhs depends on moving grid so its better to use current geo
	Daxpy(N_Active, tau, rhs, B);
        Daxpy(N_Active, tau, rhs+N_U, B+N_U);    
      }
      else
      {
          // old rhs multiplied with current subtime step and theta3 on B
         Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs, B);
         Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs+N_U, B+N_U);   
    
        // add rhs from current sub time step to rhs array B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs+N_U, B+N_U);   
      }
  // scale the pressure matrices, not in nonlinear step 
  if (scale != 1.0)
  {
   switch(NSEType)
    {
     case 1:
     case 3:
        Dscal(MatrixB1->GetN_Entries(), scale, MatrixB1->GetEntries());
        Dscal(MatrixB2->GetN_Entries(), scale, MatrixB2->GetEntries());
     break;

    case 2:
    case 4:     
         Dscal(MatrixB1T->GetN_Entries(), scale, MatrixB1T->GetEntries());
         Dscal(MatrixB2T->GetN_Entries(), scale, MatrixB2T->GetEntries());

      // scale divergence constraint
      if(val>0) 
       {
        Dscal(MatrixB1->GetN_Entries(), val*scale, MatrixB1->GetEntries());
        Dscal(MatrixB2->GetN_Entries(), val*scale, MatrixB2->GetEntries());
       }
      break;
    } // switch(NSETyp
  } //  if (scale != 1.0)
    
 
   // Also currently : M := M + gamma A
   // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A 
   // defect = M * sol
   // B:= B + defect 
   memset(defect, 0, (2*N_U+N_P)*SizeOfDouble);  
   switch(NSEType)
    {
     case 1:
     case 2:
//        MatAdd(SqmatrixM11, SqmatrixA11, -tau*TDatabase::TimeDB->THETA2);          
       gamma = - tau*TDatabase::TimeDB->THETA2;
   
       MatVectActive(SqmatrixM11, sol, defect);
       MatVectActive(SqmatrixM11, sol+N_U, defect+N_U);
       Daxpy(N_Active, 1, defect, B);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);
 
       // assembling of system matrix       
//        MatAdd(SqmatrixM11, SqmatrixA11, -gamma + tau*TDatabase::TimeDB->THETA1);   
       gamma = tau*TDatabase::TimeDB->THETA1;
     break;

     case 3:
     case 4:
       MatAdd(SqmatrixM11, SqmatrixA11, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM12, SqmatrixA12, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM21, SqmatrixA21, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM22, SqmatrixA22, - tau*TDatabase::TimeDB->THETA2);       
       gamma = - tau*TDatabase::TimeDB->THETA2;

       MatVectActive(SqmatrixM11, sol, defect);
       Daxpy(N_Active, 1, defect, B);
       MatVectActive(SqmatrixM12, sol+N_U, defect);
       Daxpy(N_Active, 1, defect, B);
       MatVectActive(SqmatrixM21, sol, defect+N_U);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);
       MatVectActive(SqmatrixM22, sol+N_U, defect+N_U);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);

       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM12, SqmatrixA12, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM21, SqmatrixA21, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM22, SqmatrixA22, -gamma + tau*TDatabase::TimeDB->THETA1);       
       gamma = tau*TDatabase::TimeDB->THETA1;     
     break;
    } 
  
   // set rhs for Dirichlet nodes
   memcpy(B+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
   memcpy(B+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
//                 cout<< "B " <<Ddot((2*N_U+N_P), B, B)<<endl; 
   SystMatAssembled  = TRUE;
//    exit(0);
   
} // AssembleSystMat

/* assemble only LHS, not rhs */
void TSystemTNSE2D_ALE::AssembleSystMatNonLinear()
{
 double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  if(SystMatAssembled)
   {
    cout << "Restore System mat before calling AssembleSystMat" <<endl;
   }

   switch(NSEType)
    {
     case 1:
     case 2: 
       // assembling of system matrix       
       MatAdd(SqmatrixM11, SqmatrixA11,   tau*TDatabase::TimeDB->THETA1);   
       gamma = tau*TDatabase::TimeDB->THETA1;
     break;

     case 3:
     case 4:       
       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM12, SqmatrixA12,  tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM21, SqmatrixA21,  tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM22, SqmatrixA22,  tau*TDatabase::TimeDB->THETA1);       
       gamma = tau*TDatabase::TimeDB->THETA1;     
     break;
    } 
  
   SystMatAssembled  = TRUE;
} // AssembleSystMatNonLinear


void TSystemTNSE2D_ALE::RestoreMassMat()
{

//   cout << "RestoreMassMat  gamma " << gamma << endl;
  if(SystMatAssembled)
   {
    // restore the mass matrix
    switch(NSEType)
     {
      case 1:
      case 2:
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma);          
       gamma = 0.;
      break;

     case 3:
     case 4:
       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma);
       MatAdd(SqmatrixM12, SqmatrixA12, -gamma);
       MatAdd(SqmatrixM21, SqmatrixA21, -gamma);
       MatAdd(SqmatrixM22, SqmatrixA22, -gamma);       
       gamma = 0.;     
     break;
    } 
    
    SystMatAssembled  = FALSE;  
   }
  else
  {
    cout << "System is not assembled to restore " <<endl;
  }
   
//   cout << "RestoreMassMat" << endl;
//   exit(0);
  
}

void TSystemTNSE2D_ALE::AssembleANonLinear(double *sol, double *rhs)
{
 int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces, last_sq;
TFESpace2D *fesprhs[2];
     N_RectMatrices = 0;          
     N_Rhs = 0;
     N_FESpaces = 3;
 
     // set the nonliner matrices
      switch(TDatabase::ParamDB->NSTYPE)
       {
        case 1:
        case 2:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[0]->Reset();

          N_SquareMatrices = 1;
        break;

        case 3:
        case 4:
          if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
           {
            SQMATRICES[0] = SqmatrixA11;
            SQMATRICES[1] = SqmatrixA22;
            SQMATRICES[0]->Reset();
            SQMATRICES[1]->Reset();

            N_SquareMatrices = 2;
            last_sq = 1;
#ifdef __PRIVATE__ 
            if (Disctype == VMS_PROJECTION)
              {
               SQMATRICES[0] = SqmatrixA11;
               SQMATRICES[1] = SqmatrixA12;
               SQMATRICES[2] = SqmatrixA21;
               SQMATRICES[3] = SqmatrixA22;
               SQMATRICES[0]->Reset();
               SQMATRICES[1]->Reset();
               SQMATRICES[2]->Reset();
               SQMATRICES[3]->Reset();

               N_SquareMatrices = 4;
               last_sq = 3;

               N_RectMatrices = 2;
               MATRICES[0] = Matrices_tilde_G11;
               MATRICES[1] = Matrices_tilde_G22;
               MATRICES[0]->Reset();
               MATRICES[1]->Reset();
       
               N_FESpaces = 4;
              }  
#endif              
           }
          else
           {
            // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }

         break;
        } // switch(TDatabase::ParamDB->NSTYPE)
         

      // assemble the nonlinear part of NSE
      Assemble2D(N_FESpaces, FeSpaces,
                 N_SquareMatrices, SQMATRICES,
                 N_RectMatrices, MATRICES,
                 N_Rhs, NULL, NULL,
                 DiscreteFormNL,
                 BoundaryConditions,
                 BoundaryValues,
                 NSEaux); 

        if(TDatabase::ParamDB->TWO_PHASE_FLOW == 1 || TDatabase::ParamDB->FREE_SURFACE_FLOW == 1)
      {
       // Adding interface entries to A11 and A22
      MatAdd(SqmatrixA11, SqmatrixF11, 1.);
      MatAdd(SqmatrixA22, SqmatrixF22, 1.);
      }
      

       // apply upwind disc
      if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == 3) )
       {
        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[last_sq], FeFct[0], FeFct[1]);
          break;
         }                        // endswitch
       }                          // endif     
       
       
      // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      { 
        N_FESpaces = 1;
        N_SquareMatrices = 2;
        N_RectMatrices = 0;
        N_Rhs = 2;

	fesprhs[0] = FeSpaces[0];
	fesprhs[1] = FeSpaces[0];
	
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;

        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, NULL,
                         N_Rhs, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=         

     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble);       
#ifdef __PRIVATE__ 
      // update matrices
      if (Disctype == VMS_PROJECTION)
        {
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixA12;
         SQMATRICES[2] = SqmatrixA21;
         SQMATRICES[3] = SqmatrixA22;
         SQMATRICES[6] =  MatricesL;
         MATRICES[2] = Matrices_tilde_G11;
         MATRICES[3] = Matrices_tilde_G22;
         MATRICES[4] = Matrices_G11;
         MATRICES[5] = Matrices_G22;

         VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
        }
#endif
} //TSystemTNSE2D_ALE::AssembleNonLinear(

  

 
void TSystemTNSE2D_ALE::Solve(double *sol)
{  
  if(!SystMatAssembled)
  {
    cout << "System Matrix is not assembled to solve " <<endl;
    exit(0);
  }
  
    switch(Solver)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
        switch(NSEType)
         {
          case 1:
            DirectSolver(SqmatrixM11, MatrixB1,  MatrixB2, B, sol);
          break;

          case 2:
             DirectSolver(SqmatrixA11, MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, B, sol);
          break;

          case 3:
           cout << "Direct solver not yet implemented for NSTYPE 3 " <<endl;
          break;

          case 4:
             DirectSolver(SqmatrixM11, SqmatrixM12, SqmatrixM21, SqmatrixM22, 
                          MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, B, sol);      
          break;
      } //  switch(NSEType) 

      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    

     
//      cout << "NSEType " << NSEType << endl;
//      exit(0);
     
}



void TSystemTNSE2D_ALE::GetTNSEResidual(double *sol, double *res)
{
  
  if(!SystMatAssembled)
   {
    cout << "System Matrix is not assembled to calculate residual " <<endl;
    exit(0);
   }
     
   
     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixM11;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;  
         break;

        case 2:
          SQMATRICES[0] = SqmatrixM11;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixM11;
          SQMATRICES[1] = SqmatrixM12;
          SQMATRICES[2] = SqmatrixM21;
          SQMATRICES[3] = SqmatrixM22;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;  
  
        break;

        case 4:
          SQMATRICES[0] = SqmatrixM11;
          SQMATRICES[1] = SqmatrixM12;
          SQMATRICES[2] = SqmatrixM21;
          SQMATRICES[3] = SqmatrixM22;

          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;
	  
 
        break;
      } //  switch(NSEType) 
      
//    Defect(sqmatrices, matrices, sol, B, res); 
   

   memset(defect, 0, 2*N_U+N_P*SizeOfDouble);
   memset(res, 0, 2*N_U+N_P*SizeOfDouble);
   
   // Block A
   MatVect(SqmatrixM11, sol, defect);
//    cout<<"MU11 : "<<Ddot(N_U, defect, defect)<<"\n";
   Daxpy(N_U, 1, defect, res);
//     cout<<"MU11 : "<<Ddot(N_U, res, res)<<"\n";
   MatVectActive(SqmatrixM12, sol+N_U, defect);
//    cout<<"MU12 : "<<Ddot(N_U, defect, defect)<<"\n";
   Daxpy(N_Active, 1, defect, res);
//    cout<<"MU12 : "<<Ddot(N_U, res, res)<<"\n";
   MatVectActive(SqmatrixM21, sol, defect+N_U);
//    cout<<"MU21 : "<<Ddot(N_U, defect+N_U, defect+N_U)<<"\n";
   Daxpy(N_Active, 1, defect+N_U, res+N_U);
//    cout<<"MU21 : "<<Ddot(N_U, res+N_U, res+N_U)<<"\n";
   MatVect(SqmatrixM22, sol+N_U, defect+N_U);
//    cout<<"MU22 : "<<Ddot(N_U, defect+N_U, defect+N_U)<<"\n";
   Daxpy(N_U, 1, defect+N_U, res+N_U);
//    cout<<"MU22 : "<<Ddot(N_U, res+N_U, res+N_U)<<"\n";

      // Block BT
      MatVect1(MatrixB1T, sol+2*N_U, defect);
//       cout<<"B1T : "<<Ddot(N_U, defect, defect)<<"\n";
      Daxpy(N_Active, 1, defect, res);
//       cout<<"B1T : "<<Ddot(N_U, res, res)<<"\n";
      MatVect1(MatrixB2T, sol+2*N_U, defect+N_U);
//       cout<<"B2T : "<<Ddot(N_U, defect+N_U, defect+N_U)<<"\n";
      Daxpy(N_Active, 1, defect+N_U, res+N_U); 
//       cout<<"B2T : "<<Ddot(N_U, res+N_U, res+N_U)<<"\n";
      
      // Block B
      MatVect1(MatrixB1, sol, defect+2*N_U);
      Daxpy(N_P, 1, defect+2*N_U, res+2*N_U); 
      MatVect1(MatrixB2, sol+N_U, defect+2*N_U);
      Daxpy(N_P, 1, defect+2*N_U, res+2*N_U); 

//        cout<<"Before rhs : "<<Ddot(2*N_U, res, res)<<"\n";
//        cout<<"B1 : "<<Ddot(N_U, B, B)<<"\n";
//        cout<<"B2 : "<<Ddot(N_U, B+N_U, B+N_U)<<"\n";
      Daxpy(2*N_U+N_P, -1., B, res);
//    cout<<"After rhs : "<<Ddot(2*N_U, res, res)<<"\n";
   
      
      

} // TSystemTNSE2D_ALE::GetResidual


// void TSystemTNSE2D_ALE::MeasureTNSEErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
//                                     double *AllErrors)
// {
//   double errors[4],  u_error[4];    
//     
//      // errors in first velocity component
//      FeFct[0]->GetErrors(ExactU1, 3, TimeNSAllDerivatives, 2,
//                          L2H1Errors,
//                          NULL, NSEaux_error, 1, FeSpaces, errors);
//       u_error[0] = errors[0];
//       u_error[1] = errors[1];
//       
//      // errors in second velocity component
//      FeFct[1]->GetErrors(ExactU2, 3, TimeNSAllDerivatives, 2,
//                          L2H1Errors,
//                          NULL, NSEaux_error, 1, FeSpaces, errors);
//      u_error[2] = errors[0];
//      u_error[3] = errors[1];      
//       
//       // errors in pressure
//      FeFct[2]->GetErrors(ExactP, 3, TimeNSAllDerivatives, 2,
//                          L2H1Errors,
//                          NULL, NSEaux_error, 1, FeSpaces+1, errors);     
// 
//      
//      // calculate all errors
//      AllErrors[0] = sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]);
//      AllErrors[1] = sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]);
//      AllErrors[2] = errors[0];
//      AllErrors[3] = errors[1];    
//      
//       // error in L^infty(0,t,L^2)
//       if(AllErrors[0] > AllErrors[5])
//        {
//         AllErrors[5]  = AllErrors[0];
//         AllErrors[4]  =  TDatabase::TimeDB->CURRENTTIME;
//       }
//       
//       // error in L^2(0,t,L^2)    
//       AllErrors[6] += (u_error[0]*u_error[0] + u_error[2]*u_error[2] +olderror_l_2_l_2u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
//       olderror_l_2_l_2u = u_error[0]*u_error[0] + u_error[2]*u_error[2];
//      
// }
 
#endif //#ifdef __2D__
