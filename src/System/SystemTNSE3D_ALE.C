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
* @brief     source file for TSystemTNSE3D_ALE
* @author    Sashikumaar Ganesan
* @date       6.4.2016
* @History    MoveBDwithVelo added 23.11.2017
 ************************************************************************  */
#include <Database.h>
#include <SystemNSE3D.h>
#include <SystemTNSE3D.h>
#include <SystemTNSE3D_ALE.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <FEVectFunct3D.h>
#include <AuxParam3D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <TNSE3D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Upwind3D.h>
#include <AssembleMat3D.h>
#include <ALE_3D.h>
#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif
#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

#include <stdlib.h>
#include <string.h>

void merge(double *a, int low, int mid, int high,int * a2,int max)
{
    double * b;
	b = new double[max];
    int *b2;
	b2 = new int[max];
    int i = low, j = mid + 1, k = 0;
  
    while (i <= mid && j <= high) {
        if (a[i] <= a[j])
	{
            b[k] = a[i];
	    b2[k] = a2[i];
	    k++;
	    i++;
	}
        else
	{
            b[k] = a[j];
	    b2[k] = a2[j];
	    k++;
	    j++;
	}
    }
    while (i <= mid)
    {
        b[k] = a[i];
	b2[k] = a2[i];
	k++;
	i++;
    }
    while (j <= high)
    {
	b[k] = a[j];
	b2[k] = a2[j];
	k++;
	j++;
    }
    k--;
    while (k >= 0) 
    {
        a[low + k] = b[k];
	a2[low + k] = b2[k];
        k--;
    }

delete []b;
delete []b2;
}
  
void mergesort(double *a, int low, int high, int * a2,int max)
{
    if (low < high) {
        int m = (high + low)/2;
        mergesort(a, low, m,a2,max);
        mergesort(a, m + 1, high,a2,max);
        merge(a, low, m, high,a2,max);
    }
}



TSystemTNSE3D_ALE::TSystemTNSE3D_ALE(int N_levels, TFESpace3D **velocity_fespace, TFESpace3D **presssure_fespace, TFEVectFunct3D **velocity, 
                                    TFEFunction3D **pressure, double **sol, double **rhs,  int disctype, int nsetype, int solver,
                                   TFESpace3D **Projection_space, TFESpace3D ** gridFESpace, TFEVectFunct3D **meshVelocity, bool conservativeale)
                                   :TSystemNSE3D(N_levels,velocity_fespace, presssure_fespace, velocity, pressure, sol, rhs, disctype,nsetype, solver)
{
 int i;  
  
  Projection_Spaces = Projection_space;
  B = new double[3*N_U+N_P];
  defect = new double[3*N_U+N_P];

  gamma =0.;  
  
  SqmatrixM11 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixM12 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixM13 = new TSquareMatrix3D*[N_levels]; 	
  SqmatrixM21 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixM22 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixM23 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixM31 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixM32 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixM33 = new TSquareMatrix3D*[N_levels]; 
  
  Mesh_restrict_aux = new double*[N_levels];

  for(i=Start_Level;i<N_Levels;i++)
   {
    // allocate the mass matrices in addition
     cout<<"NSEType::"<<NSEType<<endl;
    switch(NSEType)
     {
      case 1:
      case 2:
        SqmatrixM11[i] = new TSquareMatrix3D(sqstructureA[i]);
      break;

      case 3:
      case 4:
        SqmatrixM11[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixM12[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixM13[i] = new TSquareMatrix3D(sqstructureA[i]);	
        SqmatrixM21[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixM22[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixM23[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixM31[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixM32[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixM33[i] = new TSquareMatrix3D(sqstructureA[i]);
      break;
      
      default:
            OutPut("Unknown NSETYPE, it must be 1 to 4" << endl);
            exit(4711);;      
     }    
        
   if (Disctype == VMS_PROJECTION)
   { 
     if(NSEType==1 || NSEType==2)
     {
      Error("NSETYPE should be 3 or 4 for VMS_PROJECTION !!!" << endl);
      exit(-1);          
     }
      sqstructureL = new TSquareStructure3D(Projection_Spaces[i]);     
      sqstructureL->Sort();
      structure_tilde_G = new TStructure3D(U_Space[i], Projection_Spaces[i]); 
      structure_G = new TStructure3D(Projection_Spaces[i], U_Space[i]);
   
      MatricesL = new TSquareMatrix3D*[N_levels];
      Matrices_tilde_G11 = new TMatrix3D*[N_levels];
      Matrices_tilde_G22 = new TMatrix3D*[N_levels];
      Matrices_tilde_G33 = new TMatrix3D*[N_levels];
      Matrices_G11 = new TMatrix3D*[N_levels];
      Matrices_G22 = new TMatrix3D*[N_levels];
      Matrices_G33 = new TMatrix3D*[N_levels];
      

      sqmatrixL = new TSquareMatrix3D(sqstructureL);
      MatricesL[i] = sqmatrixL;
      LumpMassMatrixToDiag(MatricesL[i]);
      matrix_tilde_G11 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G11[i] = matrix_tilde_G11;
      matrix_tilde_G22 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G22[i] = matrix_tilde_G22;
      matrix_tilde_G33 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G33[i] = matrix_tilde_G33;
      matrix_G11 = new TMatrix3D(structure_G);
      Matrices_G11[i] = matrix_G11;
      matrix_G22 = new TMatrix3D(structure_G);
      Matrices_G22[i] = matrix_G22;
      matrix_G33 = new TMatrix3D(structure_G);
      Matrices_G33[i] = matrix_G33;
//       OutPut("dof projection : "     << endl);
//       exit(0);
    } //  if (Disctype == VMS_PROJECTION)
   } //  for(i=Start_Level;i<N_Levels;i++)
   
  NSE_Rhsaux = NULL;
  SystMatAssembled  = FALSE;
  olderror_l_2_l_2u = 0.;
 
 //=======================used for ALE=============================== 
  char WString[] = "w";  
  TFESpace3D *fesp_grid[1]; 
  MeshVelo = new double* [N_levels];

  //----------------for finest level--------------------------
  
  GridFESpace = gridFESpace;   
  
  // grid 
  SquareStructureG= new TSquareStructure3D(GridFESpace[N_Levels-1]); 
  SquareStructureG->Sort();
  
  // for mesh
  SqmatrixG11 = new TSquareMatrix3D(SquareStructureG); // G11
  SqmatrixG12 = new TSquareMatrix3D(SquareStructureG); // G12
  SqmatrixG13 = new TSquareMatrix3D(SquareStructureG); // G13
  SqmatrixG21 = new TSquareMatrix3D(SquareStructureG); // G21
  SqmatrixG22 = new TSquareMatrix3D(SquareStructureG); // G22
  SqmatrixG23 = new TSquareMatrix3D(SquareStructureG); // G23
  SqmatrixG31 = new TSquareMatrix3D(SquareStructureG); // G31
  SqmatrixG32 = new TSquareMatrix3D(SquareStructureG); // G32
  SqmatrixG33 = new TSquareMatrix3D(SquareStructureG); // G33
  
  SQMATRICES_GRID[0] = SqmatrixG11;
  SQMATRICES_GRID[1] = SqmatrixG12;
  SQMATRICES_GRID[2] = SqmatrixG13;
  SQMATRICES_GRID[3] = SqmatrixG21;
  SQMATRICES_GRID[4] = SqmatrixG22;
  SQMATRICES_GRID[5] = SqmatrixG23;
  SQMATRICES_GRID[6] = SqmatrixG31;
  SQMATRICES_GRID[7] = SqmatrixG32;
  SQMATRICES_GRID[8] = SqmatrixG33;
   
  
  Entries[0] = SqmatrixG11->GetEntries();
  Entries[1] = SqmatrixG12->GetEntries();
  Entries[2] = SqmatrixG13->GetEntries();
  Entries[3] = SqmatrixG21->GetEntries();
  Entries[4] = SqmatrixG22->GetEntries();
  Entries[5] = SqmatrixG23->GetEntries();
  Entries[6] = SqmatrixG31->GetEntries();
  Entries[7] = SqmatrixG32->GetEntries();
  Entries[8] = SqmatrixG33->GetEntries();

  GridKCol = SquareStructureG->GetKCol();
  GridRowPtr = SquareStructureG->GetRowPtr();
   
  Aux_ALE = NULL;
  SolveLinearElastic = TRUE;
  CONSERVATIVEALE = conservativeale;
  
 //---------------------for finest level---------------------------------- 
  MeshVelocity= meshVelocity;
  
  for(i=0;i<N_Levels;i++)  
   {
    N_GridDOFs = GridFESpace[i]->GetN_DegreesOfFreedom();
    N_GridActive = GridFESpace[i]->GetActiveBound();
    
    Mesh_restrict_aux[i] = new double[N_GridDOFs];
    
    MeshVelo[i] =  MeshVelocity[i]->GetValues();  
   }  
//     MeshVeloFct[0] = MeshVelocity[i]->GetComponent(0);
//     MeshVeloFct[1] = MeshVelocity[i]->GetComponent(1);
//     MeshVeloFct[2] = MeshVelocity[i]->GetComponent(2);

   fesp_grid[0] = GridFESpace[N_Levels-1];
   
    GridRhs = new double[3*N_GridDOFs];
    gridpos = new double[3*N_GridDOFs];
    gridpos_old = new double[3*N_GridDOFs];   
    gridpos_ref = new double[3*N_GridDOFs];   

    memset(gridpos, 0, 3*N_GridDOFs*SizeOfDouble);
    GridPos = new TFEVectFunct3D(fesp_grid[0], WString, WString, gridpos, N_GridDOFs, 3);  
    GridPos->GridToData();
   
    RefGridPos = new TFEVectFunct3D(fesp_grid[0], WString, WString, gridpos_ref, N_GridDOFs, 3);     
 
    memcpy(gridpos_old, gridpos, 3*N_GridDOFs*SizeOfDouble);
    memcpy(gridpos_ref, gridpos, 3*N_GridDOFs*SizeOfDouble); 
 

   // aux for Mesh Mat assemble
   Meshaux = new TAuxParam3D(1, 0, 0, 0, fesp_grid, NULL, NULL, NULL, NULL, 0, NULL); 
} // SystemTNSE3D_ALE::TSystemTNSE3D_ALE(

TSystemTNSE3D_ALE::~TSystemTNSE3D_ALE()
{
  
  
}  

void TSystemTNSE3D_ALE::Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
                             BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue, BoundCondFunct3D *GridBoundCondition,
                             BoundValueFunct3D *GridBoundValue,CoeffFct3D *GridBilinearCoeffs)
{ 
  int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  int N_U_Current;
  int velocity_space_code, pressure_space_code;
  int mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
    
  double alpha[2];  
  
  TDiscreteForm3D *DiscreteFormGalerkin, *DiscreteFormSDFEM, *DiscreteFormUpwind, *DiscreteFormSmagorinsky;
  TDiscreteForm3D *DiscreteFormVMSProjection, *DiscreteFormNLGalerkin, *DiscreteFormNLSDFEM, *DiscreteFormNLUpwind;
  TDiscreteForm3D *DiscreteFormNLSmagorinsky, *DiscreteFormNLVMSProjection,*DiscreteFormPressSep, *DiscreteFormAuxProbPressSep;
  TDiscreteForm3D *DiscreteFormNSRFBRhs, *DiscreteFormNLSDFEM_DivDiv;
  TDiscreteForm3D *DiscreteFormUpwindNC, *DiscreteFormClassicalLES, *DiscreteFormGL00Convolution, *DiscreteFormGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormVMS_Projection, *DiscreteFormVMS_SUPG, *DiscreteFormLerayAlpha, *DiscreteFormNLUpwindNC, *DiscreteFormNLClassicalLES;
  TDiscreteForm3D *DiscreteFormNLGL00Convolution, *DiscreteFormNLGL00AuxProblem, *DiscreteFormNLVMS_Projection,  *DiscreteFormNLVMS_ProjectionExpl;
  TDiscreteForm3D *DiscreteFormNLVMSRFBExplRhs, *DiscreteFormNLVMS_SUPG, *DiscreteFormNLLerayAlpha, *DiscreteFormRHS, *DiscreteFormRHSClassicalLES, *DiscreteFormRHSLES;
  TDiscreteForm3D *DiscreteFormRHSSUPG, *DiscreteFormMatrixGL00AuxProblem, *DiscreteFormGL00AuxProblemRHS, *DiscreteFormMatrixAuxProblemU;
  TDiscreteForm3D *DiscreteFormRHSAuxProblemU, *DiscreteFormRHSNewton, *DiscreteFormRHSNewtonNL, *DiscreteFormC, *DiscreteFormJ;
  
  
  /**  save the boundary condition */
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  
  BoundaryConditions[2] = BoundCond;
  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundaryConditions[1] = GridBoundCondition;
  GridBoundaryConditions[2] = GridBoundCondition;
  
  /**  save the boundary values   */
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
  BoundaryValues[2] = U3BoundValue;
  GridBoundValues[0] = GridBoundValue;
  GridBoundValues[1] = GridBoundValue;
  GridBoundValues[2] = GridBoundValue;
  
  /** save the nse bilinear coefficient   */
  LinCoeffs[0] = lincoeffs;
  
    // set the Discreteform for grid
  InitializeDiscreteFormGrid(DiscreteFormGrid, GridBilinearCoeffs);
  
//   /** set the Discreteforms */
  InitializeDiscreteForms(DiscreteFormGalerkin, DiscreteFormUpwind, DiscreteFormUpwindNC, DiscreteFormSmagorinsky,
                          DiscreteFormClassicalLES, DiscreteFormGL00Convolution, DiscreteFormGL00AuxProblem,
                          DiscreteFormVMS_Projection, DiscreteFormVMS_SUPG,DiscreteFormLerayAlpha, DiscreteFormNLGalerkin,
                          DiscreteFormNLUpwind, DiscreteFormNLUpwindNC, DiscreteFormNLSmagorinsky,
                          DiscreteFormNLClassicalLES, DiscreteFormNLGL00Convolution, DiscreteFormNLGL00AuxProblem,
                          DiscreteFormNLVMS_Projection, DiscreteFormNLVMS_ProjectionExpl, DiscreteFormNLVMSRFBExplRhs,
                          DiscreteFormNLVMS_SUPG, DiscreteFormNLLerayAlpha, DiscreteFormRHS, DiscreteFormRHSClassicalLES, DiscreteFormRHSLES,
                          DiscreteFormRHSSUPG, DiscreteFormMatrixGL00AuxProblem, DiscreteFormGL00AuxProblemRHS,
                          DiscreteFormMatrixAuxProblemU, DiscreteFormRHSAuxProblemU, DiscreteFormRHSNewton,
                          DiscreteFormRHSNewtonNL, DiscreteFormC, DiscreteFormJ, LinCoeffs[0], NSEType);
 
    /** find discrete form */
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
            DiscreteFormRhs = DiscreteFormRHS;    
          break;
 

          case UPWIND:
            DiscreteFormARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;   
            DiscreteFormRhs = DiscreteFormRHS;	     
            break;

          case SMAGORINSKY:
            DiscreteFormARhs = DiscreteFormSmagorinsky;
            DiscreteFormNL = DiscreteFormNLSmagorinsky;  
            DiscreteFormRhs = DiscreteFormRHS; 
            break;
 
          case VMS_PROJECTION:
            DiscreteFormARhs = DiscreteFormVMS_Projection;
            DiscreteFormNL = DiscreteFormNLVMS_Projection;
            DiscreteFormRhs = DiscreteFormRHS; 
            break;
    
          default:
            Error("Unknown DISCTYPE" << endl);
            exit(-1);
        } 
     
     // set the discrete form for the Stokes equation
      if (TDatabase::ParamDB->PROBLEM_TYPE==STOKES)
       {
        DiscreteFormARhs = DiscreteFormUpwind;     
        DiscreteFormNL = NULL;
        DiscreteFormRhs = DiscreteFormRHS; 
       }

#ifdef _OMPONLY
     if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
      {
       cout<<"NOT YET IMPLEMENTED !!!"<<endl;
       exit(0);
      }
#endif  
 
    // initilize the assemble    
    for(i=Start_Level;i<N_Levels;i++)
    {    
     // initialize the aux
     fesp_aux[0] =  U_Space[i];
     fesp_aux[1] =  P_Space[i];
     fesp_aux[2] =  GridFESpace[i];

     fefct_aux[i*6 ] = Velocity[i]->GetComponent(0);
     fefct_aux[i*6+1] = Velocity[i]->GetComponent(1);
     fefct_aux[i*6+2] = Velocity[i]->GetComponent(2);

     fefct_aux[i*6+3] = MeshVelocity[i]->GetComponent(0);
     fefct_aux[i*6+4] = MeshVelocity[i]->GetComponent(1);
     fefct_aux[i*6+5] = MeshVelocity[i]->GetComponent(2);

      switch(Disctype)
       {
         // turbulent viscosity must be computed
        case SMAGORINSKY:
        case CLASSICAL_LES:
        case GL00_CONVOLUTION:
        case GL00_AUX_PROBLEM:
        case SDFEM:

         cout<< " aux not yet implemented  for this NSType!!! " << endl;
         exit(0);
 
        break; 

        case VMS_PROJECTION:
        case VMS_PROJECTION_EXPL:

         NSEaux[i] =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo_VMS3D_ALE, TimeNSN_FctVelo_GradVelo_VMS3D_ALE,
                                   TimeNSN_ParamFctVelo_GradVelo_VMS3D_ALE, TimeNSN_FEValuesVelo_GradVelo_VMS3D_ALE,
                                   fesp_aux, fefct_aux+(i*6), TimeNSFctVelo_GradVelo_VMS3D_ALE,
                                   TimeNSFEFctIndexVelo_GradVelo_VMS3D_ALE, TimeNSFEMultiIndexVelo_GradVelo_VMS3D_ALE,
                                   TimeNSN_ParamsVelo_GradVelo_VMS3D_ALE, TimeNSBeginParamVelo_GradVelo_VMS3D_ALE);
        break; 

        case UPWIND:
          NSEaux[i] = new TAuxParam3D(1, 0, 0, 0, fesp_aux, NULL, NULL, NULL, NULL, 0, NULL);
          break; 
   
        default:
            NSEaux[i] =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo, TimeNSN_ParamFctVelo, TimeNSN_FEValuesVelo,
                                      fesp_aux, fefct_aux+(i*6), TimeNSFctVelo, TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
                                      TimeNSN_ParamsVelo, TimeNSBeginParamVelo);
        } 

      
      // initialize matrices
      N_Rhs = 3;
      N_FESpaces = 2;   
     
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();   
      
      fesp[0] =  U_Space[i];
      fesp[1] =  P_Space[i];
      
      fefct[0] = Velocity[i]->GetComponent(0);
      fefct[1] = Velocity[i]->GetComponent(1);
      fefct[2] = Velocity[i]->GetComponent(2);  
      
        //-----------used for ALE-------    
     if (i==N_Levels - 1)   
      {   
       fesp[2] = GridFESpace[i];
     
       fefct[3] = MeshVelocity[i]->GetComponent(0);
       fefct[4] = MeshVelocity[i]->GetComponent(1);
       fefct[5] = MeshVelocity[i]->GetComponent(2);    
      }  
  //-------------------------------------  
     
      fesprhs[0] =  U_Space[i];
      fesprhs[1] =  U_Space[i];
      fesprhs[2] =  U_Space[i];

      RHSs[0] = RhsArray[i];
      RHSs[1] = RhsArray[i] + N_U_Current;
      RHSs[2] = RhsArray[i] + 2*N_U_Current;
      RHSs[3] = RhsArray[i] + 3*N_U_Current;    
      
      
     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixA11[i];
          SQMATRICES[1] = SqmatrixM11[i];
          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];
 
          N_SquareMatrices = 2;
          N_RectMatrices = 3;
        break;

        case 2:
          SQMATRICES[0] = SqmatrixA11[i];
          SQMATRICES[1] = SqmatrixM11[i];
          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];
          MATRICES[3] = MatrixB1T[i];
          MATRICES[4] = MatrixB2T[i];
          MATRICES[5] = MatrixB3T[i];
  
          N_SquareMatrices = 2;
          N_RectMatrices = 6;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11[i];
          SQMATRICES[1] = SqmatrixA12[i];
          SQMATRICES[2] = SqmatrixA13[i];  
          SQMATRICES[3] = SqmatrixA21[i];
          SQMATRICES[4] = SqmatrixA22[i];
          SQMATRICES[5] = SqmatrixA23[i]; 
          SQMATRICES[6] = SqmatrixA31[i];
          SQMATRICES[7] = SqmatrixA32[i];
          SQMATRICES[8] = SqmatrixA33[i];  
          SQMATRICES[9] = SqmatrixM11[i];
          SQMATRICES[10] = SqmatrixM22[i];
          SQMATRICES[11] = SqmatrixM33[i];

          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];

          N_SquareMatrices = 12;
          N_RectMatrices = 3;
  
          Error("NSETYPE  3 is not yet implemented for VMS_PROJECTION !!!" << endl);
          exit(-1);          

        break;

        case 4:
          SQMATRICES[0] = SqmatrixA11[i];
          SQMATRICES[1] = SqmatrixA12[i];
          SQMATRICES[2] = SqmatrixA13[i];	  
          SQMATRICES[3] = SqmatrixA21[i];
          SQMATRICES[4] = SqmatrixA22[i];
          SQMATRICES[5] = SqmatrixA23[i]; 
          SQMATRICES[6] = SqmatrixA31[i];
          SQMATRICES[7] = SqmatrixA32[i];
          SQMATRICES[8] = SqmatrixA33[i];  
          SQMATRICES[9] = SqmatrixM11[i];
          SQMATRICES[10] = SqmatrixM22[i];
          SQMATRICES[11] = SqmatrixM33[i];
 
          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];
          MATRICES[3] = MatrixB1T[i];
          MATRICES[4] = MatrixB2T[i];
          MATRICES[5] = MatrixB3T[i];

          N_SquareMatrices = 12;
          N_RectMatrices = 6;
    
        if(Disctype == VMS_PROJECTION)
         {
          N_SquareMatrices = 13;
          SQMATRICES[12] =  MatricesL[i];
          SQMATRICES[12]->Reset();

          N_RectMatrices = 12;
          MATRICES[6] = Matrices_tilde_G11[i];
          MATRICES[7] = Matrices_tilde_G22[i];
          MATRICES[8] = Matrices_tilde_G33[i];
          MATRICES[9] = Matrices_G11[i];
          MATRICES[10] = Matrices_G22[i];
          MATRICES[11] = Matrices_G33[i];
          MATRICES[6]->Reset();
          MATRICES[7]->Reset();
          MATRICES[8]->Reset();
          MATRICES[9]->Reset();
          MATRICES[10]->Reset();
          MATRICES[11]->Reset();

          
          N_FESpaces = 4;
          fesp[3] = Projection_Spaces[i];
         }
        break;
      } //  switch(NSEType)
       
     // array of assemble objects
     AMatRhsAssemble[i] = new TAssembleMat3D(N_FESpaces, fesp, N_SquareMatrices, SQMATRICES, N_RectMatrices,
                                             MATRICES, N_Rhs, RHSs, fesprhs, DiscreteFormARhs, BoundaryConditions,
                                             BoundaryValues, NSEaux[i]);                            
     AMatRhsAssemble[i]->Init();    

      // set Rhs Only Assemble
     if(i == N_Levels-1)
      {
       //aux for error calculation
       NSEaux_error[i] =  new TAuxParam3D(TimeNSN_FESpacesVelo, 
		 TimeNSN_FctVelo, 
		 TimeNSN_ParamFctVelo,
		 TimeNSN_FEValuesVelo,
                                          fesp,
		 fefct, 
		 TimeNSFctVelo,
		 TimeNSFEFctIndexVelo,
		 TimeNSFEMultiIndexVelo,
		 TimeNSN_ParamsVelo,
		 TimeNSBeginParamVelo);       
       
       
       // it is suffucuent to assemble only at finnest level, as the coarse level rhs' are not used in multigrid solver
       N_FESpaces = 1;
       N_SquareMatrices = 0;
       N_RectMatrices = 0;
      
       RhsOnlyAssemble  = new TAssembleMat3D(N_FESpaces,
		    fesp, 
		    N_SquareMatrices,
		    NULL, 
		    N_RectMatrices,
		    NULL,
                                             N_Rhs, 
		    RHSs, 
		    fesprhs, 
		    DiscreteFormRhs, 
		    BoundaryConditions,
		    BoundaryValues,
		    NSEaux[i]);
       RhsOnlyAssemble->Init();   
      }
     
#ifdef _MPI
   if(i == N_Levels-1) {
     
	       SQMATRICES[0] = SqmatrixM11[i];
               SQMATRICES[1] = SqmatrixM12[i];
               SQMATRICES[2] = SqmatrixM13[i];	  
               SQMATRICES[3] = SqmatrixM21[i];
               SQMATRICES[4] = SqmatrixM22[i];
               SQMATRICES[5] = SqmatrixM23[i]; 
               SQMATRICES[6] = SqmatrixM31[i];
               SQMATRICES[7] = SqmatrixM32[i];
               SQMATRICES[8] = SqmatrixM33[i];
     
     
    if(SOLVER == DIRECT)
     {
      DS = new TParDirectSolver(ParComm_U[N_Levels-1],ParComm_P[N_Levels-1],SQMATRICES,MATRICES);
     }
   }
#endif
#ifdef _SMPI
   if(i == N_Levels-1) 
   {
     
	       SQMATRICES[0] = SqmatrixM11[i];
               SQMATRICES[1] = SqmatrixM12[i];
               SQMATRICES[2] = SqmatrixM13[i];	  
               SQMATRICES[3] = SqmatrixM21[i];
               SQMATRICES[4] = SqmatrixM22[i];
               SQMATRICES[5] = SqmatrixM23[i]; 
               SQMATRICES[6] = SqmatrixM31[i];
               SQMATRICES[7] = SqmatrixM32[i];
               SQMATRICES[8] = SqmatrixM33[i];
     
     
    if(SOLVER == DIRECT)
     {
      P_DS = new TSeqParDirectSolver(3,N_U,N_P,SQMATRICES,MATRICES);
     }
   }
#endif
 
     
 //    set the nonliner matrices
     N_RectMatrices = 0;          
     N_Rhs = 0;
     N_FESpaces = 1;
     
     switch(NSEType)
       {
        case 1:
        case 2:
          SQMATRICES[0] = SqmatrixA11[i];
          N_SquareMatrices = 1;
        break;

        case 3:
        case 4:
          if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
           {
            SQMATRICES[0] = SqmatrixA11[i];
            SQMATRICES[1] = SqmatrixA22[i];
            SQMATRICES[2] = SqmatrixA33[i];
            N_SquareMatrices = 3;

            if (Disctype == VMS_PROJECTION)
              {
               SQMATRICES[0] = SqmatrixA11[i];
               SQMATRICES[1] = SqmatrixA12[i];
               SQMATRICES[2] = SqmatrixA13[i];	  
               SQMATRICES[3] = SqmatrixA21[i];
               SQMATRICES[4] = SqmatrixA22[i];
               SQMATRICES[5] = SqmatrixA23[i]; 
               SQMATRICES[6] = SqmatrixA31[i];
               SQMATRICES[7] = SqmatrixA32[i];
               SQMATRICES[8] = SqmatrixA33[i];

               SQMATRICES[0]->Reset();
               SQMATRICES[1]->Reset();
               SQMATRICES[2]->Reset();
               SQMATRICES[3]->Reset();
               SQMATRICES[4]->Reset();
               SQMATRICES[5]->Reset();
               SQMATRICES[6]->Reset();
               SQMATRICES[7]->Reset();
               SQMATRICES[8]->Reset();

               N_SquareMatrices = 9;

               N_RectMatrices = 3;
               MATRICES[0] = Matrices_tilde_G11[i];
               MATRICES[1] = Matrices_tilde_G22[i];
               MATRICES[2] = Matrices_tilde_G33[i];
               MATRICES[0]->Reset();
               MATRICES[1]->Reset();
               MATRICES[2]->Reset();
       
                N_FESpaces = 4;
                fesp[3] = Projection_Spaces[i];
              }  
           }
          else
           {
            // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }

         break;
        } // switch(NSEType)
                
     AMatAssembleNonLinear[i] = new TAssembleMat3D(N_FESpaces, fesp, 
 N_SquareMatrices, SQMATRICES,
 N_RectMatrices, MATRICES,
                                                   N_Rhs, NULL, NULL,
 DiscreteFormNL,
 BoundaryConditions, 
 BoundaryValues,
 NSEaux[i]);
     AMatAssembleNonLinear[i]->Init();      
     
     //===============================================================================================================
// //      // set the slip with friction assemble matrices
// //        if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
// //         {
// //          if(NSEType <4)
// //           {
// //            OutPut("For slip with friction bc NSTYPE 4 is necessary !!!!! " << endl);
// //            exit(4711);
// //           }
//           
// //           // prepare everything for the assembling of slip with friction bc
// //           // on all levels
// //           N_FESpaces = 1;
// //           N_SquareMatrices = 9;
// //           N_RectMatrices = 0;
// //           N_Rhs = 3; 
// //       
// //           SQMATRICES[0] = SqmatrixA11[i];
// //           SQMATRICES[1] = SqmatrixA22[i];
// //           SQMATRICES[2] = SqmatrixA33[i];
// //           SQMATRICES[3] = SqmatrixA12[i];
// //           SQMATRICES[4] = SqmatrixA13[i];
// //           SQMATRICES[5] = SqmatrixA21[i];
// //           SQMATRICES[6] = SqmatrixA23[i];
// //           SQMATRICES[7] = SqmatrixA31[i];
// //           SQMATRICES[8] = SqmatrixA32[i];
// // 
// //           Assemble3DSlipBC(N_FESpaces, fesp,
// //             N_SquareMatrices, SQMATRICES,
// //             N_RectMatrices, NULL,
// //             N_Rhs, RHSs, fesprhs,
// //             NULL,
// //             BoundaryConditions,
// //             BoundaryValues,
// //             NSEaux);     
// //      
//      
// // 	} // if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
     
     //===============================================================================================================     
     // initialize solver
       if(SOLVER==GMG)
        {       
         //setup the multigrid solver
         alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
         alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;  
         velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
         pressure_space_code = TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE;

         if(mg_type==1)
          {
           if(i==0)
            {
             alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
             alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
            }     
           else if(i==N_Levels-1)
            {
             alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
             alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE;  
            }
          
           if(i<N_Levels-1)
            {
             velocity_space_code = -1;
             pressure_space_code = 0; 
            }          
          }  
  
         switch(NSEType)
          {
           case 1:
              MGLevel = new TNSE_MGLevel1(i, SqmatrixM11[i],  MatrixB1[i], MatrixB2[i], MatrixB3[i],
                                             structureBT[i], RhsArray[i], SolArray[i], N_aux, alpha,
                                             velocity_space_code, pressure_space_code,
                                             NULL, NULL);
              MG->AddLevel(MGLevel);
          break;
          case 2:
              MGLevel = new TNSE_MGLevel2(i, SqmatrixM11[i], MatrixB1[i], MatrixB2[i], MatrixB3[i],
                                             MatrixB1T[i], MatrixB2T[i], MatrixB3T[i], RhsArray[i], SolArray[i], N_aux, alpha,
                                             velocity_space_code, pressure_space_code,
                                             NULL, NULL);
              MG->AddLevel(MGLevel);
          break;
          case 3:
              MGLevel = new TNSE_MGLevel3(i, SqmatrixM11[i], SqmatrixM12[i], SqmatrixM13[i],
                                             SqmatrixM21[i], SqmatrixM22[i], SqmatrixM23[i],
                                             SqmatrixM31[i], SqmatrixM32[i], SqmatrixM33[i],
                                             MatrixB1[i], MatrixB2[i], MatrixB3[i],
                                             structureBT[i], RhsArray[i], SolArray[i], N_aux, alpha,
                                             velocity_space_code, pressure_space_code,
                                             NULL, NULL);
              MG->AddLevel(MGLevel);
          break;
          case 4:
              MGLevel = new TNSE_MGLevel4(i, SqmatrixM11[i], SqmatrixM12[i], SqmatrixM13[i],
                                             SqmatrixM21[i], SqmatrixM22[i], SqmatrixM23[i],
                                             SqmatrixM31[i], SqmatrixM32[i], SqmatrixM33[i],
                                             MatrixB1[i], MatrixB2[i], MatrixB3[i],
                                             MatrixB1T[i], MatrixB2T[i], MatrixB3T[i], RhsArray[i], SolArray[i], N_aux, alpha,
                                             velocity_space_code, pressure_space_code,
                                             NULL, NULL
		#ifdef _MPI
	  , ParComm_U[i], ParComm_P[i]
#endif		
	      );
              MG->AddLevel(MGLevel);
          break;	
        } //  switch(NSEType)
       }  // if(SOLVER==GMG)     
     } // for(i=Start_Level;i<N_Levels;i++)      
         
         
   // Mesh mat assemble
   SQMATRICES_GRID[0]->Reset();
   SQMATRICES_GRID[1]->Reset();
   SQMATRICES_GRID[2]->Reset();
   SQMATRICES_GRID[3]->Reset();
   SQMATRICES_GRID[4]->Reset();
   SQMATRICES_GRID[5]->Reset();
   SQMATRICES_GRID[6]->Reset();
   SQMATRICES_GRID[7]->Reset();
   SQMATRICES_GRID[8]->Reset();
      
   fesp[0] = GridFESpace[N_Levels-1];    
   
   MeshMatAssemble = new TAssembleMat3D(1, fesp, 9, SQMATRICES_GRID, 0, NULL, 0, NULL, NULL, DiscreteFormGrid, 
                                        GridBoundaryConditions, GridBoundValues, Meshaux);                 
   MeshMatAssemble->Init();
   
//  cout << " TSystemTNSE3D_ALE::Init done ! " << endl;
//          exit(0);
	 
} // TSystemNSE3D::Init

/** add info of moving boundaries move with velocity*/
void TSystemTNSE3D_ALE::AddMoveBoundFunction(MoveBound_3D *movebdwithvelo, int n_movVert, TVertex **movboundVert,
                                int n_movfaces, int *movfaces, TBaseCell ** movcells)
{ 
  MoveBDwithVelo = movebdwithvelo;
  N_MovVert     = n_movVert;
  N_Movfaces    = n_movfaces;
  Movfaces      = movfaces;
  MovBoundVert  = movboundVert;
  MovCells      = movcells;
  
  SolveLinearElastic = TRUE;
  BdMoveWithVelo = TRUE;
 }  // 


/** add info of moving boundaries */
void TSystemTNSE3D_ALE::AddBoundModifyFunction(ModifyBoundCoords_3D *modifyboudary, int n_movVert, TVertex **movboundVert,
                                int n_movfaces, int *movfaces, TBaseCell ** movcells)
{ 
  ModifyBoudary = modifyboudary;
  N_MovVert     = n_movVert;
  N_Movfaces    = n_movfaces;
  Movfaces      = movfaces;
  MovBoundVert  = movboundVert;
  MovCells      = movcells;
  
  SolveLinearElastic = TRUE;
  BdMoveWithVelo  = FALSE;
 }  // TSystemTNSE3D_ALE::AddBoundModifyFunction(

 /** assemble only rhs */
void TSystemTNSE3D_ALE::AssembleRhs()
 { 
   // it is sufficient to assemble only the finnest level RHS, as lo RHS are not necessary for multigrid solver
   
   // initialize rhs
   RhsOnlyAssemble->Reset();
    
   // cout << " RhsOnlyAssemble " << endl;
   // asssemble rhs
   RhsOnlyAssemble->Assemble3D();
   
   // set rhs for Dirichlet nodes
//    memcpy(SolArray[N_Levels-1]+N_Active, RhsArray[N_Levels-1]+N_Active, N_DirichletDof*SizeOfDouble);
//    memcpy(SolArray[N_Levels-1]+N_U+N_Active, RhsArray[N_Levels-1]+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
//    memcpy(SolArray[N_Levels-1]+2*N_U+N_Active, RhsArray[N_Levels-1]+2*N_U+N_Active, N_DirichletDof*SizeOfDouble);     
 }
   
   
/** assemble entire matrices and rhs */
void TSystemTNSE3D_ALE::Assemble()
 { 
  int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  int N_U_Current, N_P_Current, N_Active_Current, N_DirichletDof;
    
  double alpha[2];

   // before assemble NSE mat, the mesh velocity has to be calculated
   this->AssembleMeshMat();

   //Get the mesh velocity
   this->GetMeshVelo(FALSE, 0);
 
   for(i=Start_Level;i<N_Levels;i++)
    {     
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();
      N_Active_Current  = U_Space[i]->GetActiveBound();     
      N_DirichletDof = N_U_Current - N_Active_Current;
      N_P_Current = P_Space[i]->GetN_DegreesOfFreedom();      
      
      // initialize matrices
      AMatRhsAssemble[i]->Reset();
      
      //these matrices are not part of AMatRhsAssemble, so reset here
      if(NSEType==3 || NSEType==4)
       {
        SqmatrixM12[i]->Reset();
        SqmatrixM13[i]->Reset();  
        SqmatrixM21[i]->Reset();
        SqmatrixM23[i]->Reset(); 
        SqmatrixM31[i]->Reset();
        SqmatrixM32[i]->Reset();
       }

      /** assemble */
      AMatRhsAssemble[i]->Assemble3D();

      fefct[0] = Velocity[i]->GetComponent(0);
      fefct[1] = Velocity[i]->GetComponent(1);
      fefct[2] = Velocity[i]->GetComponent(2);
      
       if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE==3) )
        {
         switch(NSEType)
          {
           case 1:
           case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes3D(SqmatrixA11[i], fefct[0], fefct[1], fefct[2]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with three matrices
            UpwindForNavierStokes3D(SqmatrixA11[i], fefct[0], fefct[1], fefct[2]);
            UpwindForNavierStokes3D(SqmatrixA22[i], fefct[0], fefct[1], fefct[2]);
            UpwindForNavierStokes3D(SqmatrixA33[i], fefct[0], fefct[1], fefct[2]);
            cout << "UPWINDING DONE : level " << endl;
           break;
         }                        // endswitch
        }                          // endif     
       

        if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
         if(NSEType <4)
          {
           OutPut("For slip with friction bc NSTYPE 4 is ");
           OutPut("necessary !!!!! " << endl);
           exit(4711);
          }
          
//          AMatRhsAssemble[i]->AssembleNavierSlip(); 
         // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 9;
          N_RectMatrices = 0;
          N_Rhs = 3; 

          fesp[0] =  U_Space[i];
	  
          fesprhs[0] =  U_Space[i];
          fesprhs[1] =  U_Space[i];
          fesprhs[2] =  U_Space[i];
      
          SQMATRICES[0] = SqmatrixA11[i];
          SQMATRICES[1] = SqmatrixA22[i];
          SQMATRICES[2] = SqmatrixA33[i];
          SQMATRICES[3] = SqmatrixA12[i];
          SQMATRICES[4] = SqmatrixA13[i];
          SQMATRICES[5] = SqmatrixA21[i];
          SQMATRICES[6] = SqmatrixA23[i];
          SQMATRICES[7] = SqmatrixA31[i];
          SQMATRICES[8] = SqmatrixA32[i];

          RHSs[0] = RhsArray[i];
          RHSs[1] = RhsArray[i] + N_U_Current;
          RHSs[2] = RhsArray[i] + 2*N_U_Current;
      
          Assemble3DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, NULL,
                           N_Rhs, RHSs, fesprhs,
                           NULL,
                           BoundaryConditions,
                           BoundaryValues,
                           NSEaux[i]);

        } //  if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRIC

// //       delete NSEaux; 

      // set rhs for Dirichlet nodes
      memcpy(SolArray[i]+N_Active_Current, RhsArray[i]+N_Active_Current, N_DirichletDof*SizeOfDouble);
      memcpy(SolArray[i]+N_U_Current+N_Active_Current, RhsArray[i]+N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble); 
      memcpy(SolArray[i]+2*N_U_Current+N_Active_Current, RhsArray[i]+2*N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble);     

#ifdef _MPI
ParComm_U[i]->CommUpdate(SolArray[i]);
ParComm_P[i]->CommUpdate(SolArray[i] + 3*N_U);
#endif
 
    if (Disctype == VMS_PROJECTION)
    {     
     if (i==N_Levels-1)
      {    
       LumpMassMatrixToDiag(MatricesL[i]);
  
       SQMATRICES[0] = SqmatrixA11[i];
       SQMATRICES[1] = SqmatrixA12[i];
       SQMATRICES[2] = SqmatrixA13[i];
       SQMATRICES[3] = SqmatrixA21[i];
       SQMATRICES[4] = SqmatrixA22[i];
       SQMATRICES[5] = SqmatrixA23[i];
       SQMATRICES[6] = SqmatrixA31[i];
       SQMATRICES[7] = SqmatrixA32[i];
       SQMATRICES[8] = SqmatrixA33[i];
       SQMATRICES[9] =  MatricesL[i];
       MATRICES[0] = Matrices_tilde_G11[i];
       MATRICES[1] = Matrices_tilde_G22[i];
       MATRICES[2] = Matrices_tilde_G33[i];
       MATRICES[3] = Matrices_G11[i];
       MATRICES[4] = Matrices_G22[i];
       MATRICES[5] = Matrices_G33[i];

       VMS_ProjectionUpdateMatrices(U_Space[i]->GetN_DegreesOfFreedom(), U_Space[i]->GetActiveBound(), 
                                    Projection_Spaces[i]->GetN_DegreesOfFreedom(), SQMATRICES, MATRICES);
       //  OutPut("update done"<<endl);
       }      
   
      }
     } // for(i=Start_Level;i<N_Levels;i++)  
   
//     OutPut("necessary !!!!! " << endl);
//     exit(0);
 }
 
/** Get Mesh Velo */ 
void TSystemTNSE3D_ALE::GetMeshVelo(bool MoveMesh, int solver_flag)
{
 int i, j, N_GridBDDOFs, k=0;
 double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
 double  Currtime = TDatabase::TimeDB->CURRENTTIME;
     
    N_GridDOFs = GridFESpace[N_Levels-1]->GetN_DegreesOfFreedom();
    N_GridActive = GridFESpace[N_Levels-1]->GetActiveBound();

    GridPos->GridToData();   
    memcpy(gridpos_old, gridpos, 3*N_GridDOFs*SizeOfDouble);  
    
    if(BdMoveWithVelo)
     {
      // move the BDpos with Velocity
      MoveBDwithVelo(Velocity[N_Levels-1], GridPos, gridpos_old, N_MovVert, MovBoundVert, N_Movfaces, Movfaces, MovCells, tau);
     }
    else
     {
      // if the ALE reference domain is initial domain, then uncomment the below line
      //RefGridPos->DataToGrid();

      // modify the boundary vertices
      ModifyBoudary(N_MovVert, MovBoundVert, N_Movfaces, Movfaces,MovCells, Currtime);    
      // data with updated BD values
      GridPos->GridToData();
     }

    N_GridBDDOFs = N_GridDOFs - N_GridActive;
    memset(GridRhs, 0, 3*N_GridDOFs*SizeOfDouble);    

    memcpy(GridRhs+N_GridActive, gridpos+N_GridActive, N_GridBDDOFs*SizeOfDouble);     //rhs1  
    memcpy(GridRhs+(N_GridActive+N_GridDOFs), gridpos+(N_GridActive+N_GridDOFs), N_GridBDDOFs*SizeOfDouble); //rhs2 
    memcpy(GridRhs+(N_GridActive+2*N_GridDOFs),gridpos+(N_GridActive+2*N_GridDOFs), N_GridBDDOFs*SizeOfDouble);   //rhs3  

    Daxpy(N_GridBDDOFs, -1., gridpos_old+N_GridActive, GridRhs+N_GridActive);
    Daxpy(N_GridBDDOFs, -1., gridpos_old+(N_GridDOFs+N_GridActive), GridRhs+(N_GridDOFs+N_GridActive));
    Daxpy(N_GridBDDOFs, -1., gridpos_old+(2*N_GridDOFs+N_GridActive), GridRhs+(2*N_GridDOFs+N_GridActive));

    //call direct solver

    DirectSolver(SQMATRICES_GRID, 3, 3, MeshVelo[N_Levels-1], GridRhs, MeshMatSolver_Values, MeshMatSolver_KCol,
                 MeshMatSolver_Row, MeshMatSolver_Symbolic, MeshMatSolver_Numeric,  solver_flag);

    memcpy(gridpos, gridpos_old, 3*N_GridDOFs*SizeOfDouble);
    //move the mesh back to old position
    if(!MoveMesh)
     { 
      GridPos->DataToGrid(); 
     }
    
    Daxpy(3*N_GridDOFs, 1., MeshVelo[N_Levels-1], gridpos);
 
    // mesh velocity
    Dscal(3*N_GridDOFs, 1./tau, MeshVelo[N_Levels-1]);

    //populate mesh velo to lower levels
    if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
     {
      for(i=N_Levels-1;i>0;i++)
       RestrictFunction(GridFESpace[i-1],GridFESpace[i],3,MeshVelo[i-1],MeshVelo[i],Mesh_restrict_aux[i-1]);
     }
 
    //move the mesh
    if(MoveMesh)
     { 
      GridPos->DataToGrid(); 
      // free the memory, as the matrices will be assembled, rb_flag = 4 ==> only free up memory
      DirectSolver(SQMATRICES_GRID, 3, 3, MeshVelo[N_Levels-1], GridRhs, MeshMatSolver_Values, MeshMatSolver_KCol,
                   MeshMatSolver_Row, MeshMatSolver_Symbolic, MeshMatSolver_Numeric,  4);   
     }     

   // cout<< "MoveMesh  " <<Ddot((3*N_GridDOFs), MeshVelo[N_Levels-1], MeshVelo[N_Levels-1])<<endl;
} // TSystemTNSE3D_ALE::GetMeshVelo

void TSystemTNSE3D_ALE::GetMeshVeloAndMove(double Currtime, double tau)
{
 int i,k;
 
//   for( k=0;k<N_Levels;k++){
//      N_GridDOFs = GridFESpace[k]->GetN_DegreesOfFreedom();
//     GridPos->GridToData();   
//     memcpy(gridpos_old[k], gridpos[k], 3*N_GridDOFs*SizeOfDouble);  
// 
//     //  move velo in current time  
//     for(i=0;i<N_GridDOFs;i++)
//       ModifyCoord(gridpos_ref[i][k], gridpos_ref[i+N_GridDOFs][k],
// 		  gridpos_ref[i+2*N_GridDOFs][k], gridpos[i][k],
// 		  gridpos[i+N_GridDOFs][k],  gridpos[i+2*N_GridDOFs][k],Currtime);   
// 
//     //compute mesh velocity
//     memcpy(MeshVelo[k], gridpos[k], 3*N_GridDOFs*SizeOfDouble);     
//     Daxpy(3*N_GridDOFs, -1., gridpos_old[k], MeshVelo[k]);        
//     Dscal(3*N_GridDOFs, -1./tau, MeshVelo[k]); // - sign du*/ //e to -w\cdot\nabla C in the equation   
//     
//     //move the mesh
//     GridPos->DataToGrid(); 
//   }
//    memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble); 
cout<< "TSystemTNSE3D_ALE::GetMeshVeloAndMove "  <<endl; 
 exit(0);
 
} //TSystemTNSE2D_ALE::GetMeshVelo
 
void TSystemTNSE3D_ALE::AssembleMeshMat()
{      
  //assemble the mesh matrices
  MeshMatAssemble->Reset(); 

  TDatabase::ParamDB->ASSEMBLEMESHMAT=TRUE;
  MeshMatAssemble->Assemble3D();
  TDatabase::ParamDB->ASSEMBLEMESHMAT=FALSE;

  // set Dirichlet rows in off-diagonal mesh matrices
  memset(Entries[1] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble); 
  memset(Entries[3] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[5] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[6] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[7] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  
  if (TDatabase::ParamDB->MESH_SLIP_WITH_FRICTION >= 1)
   {
	// set SLIP_WITH_FRICTION rows for no-penetration
//     MeshMatAssemble->SetNoPenetration();

      fesp[0] =  GridFESpace[N_Levels-1];

      fesprhs[0] =  GridFESpace[N_Levels-1];
      fesprhs[1] =  GridFESpace[N_Levels-1];
      fesprhs[2] =  GridFESpace[N_Levels-1];

      RHSs[0] = GridRhs;
      RHSs[1] = GridRhs + N_GridDOFs;
      RHSs[2] = GridRhs + 2*N_GridDOFs;

      SQMATRICES_GRID[0] = SqmatrixG11;
      SQMATRICES_GRID[1] = SqmatrixG12;
      SQMATRICES_GRID[2] = SqmatrixG13;
      SQMATRICES_GRID[3] = SqmatrixG21;
      SQMATRICES_GRID[4] = SqmatrixG22;
      SQMATRICES_GRID[5] = SqmatrixG23;
      SQMATRICES_GRID[6] = SqmatrixG31;
      SQMATRICES_GRID[7] = SqmatrixG32;
      SQMATRICES_GRID[8] = SqmatrixG33;

      Assemble3DSlipBC(1, fesp,
                       9, SQMATRICES_GRID,
                       0, NULL,
                       3, RHSs, fesprhs,
                       NULL,
					   GridBoundaryConditions,
                       GridBoundValues, Meshaux);


   } //  if (TDatabase::ParamDB->MESH_SLIP_WITH_FRICTION >= 1)

//   cout << " AssembleMeshMat done " <<endl;
//   exit(0);
  
} //AssembleMeshMat
 

 
void TSystemTNSE3D_ALE::AssembleNonLinear()
{
 int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
 int N_U_Current, N_Active_Current, N_DirichletDof;

  if (TDatabase::ParamDB->PROBLEM_TYPE==STOKES)
   {
    cout<<"AssembleNonLinear cannot be called for STOKES problem !!!"<<endl;
    exit(0);
   }     
  
   if(SystMatAssembled)
   {
    cout << "Restore System mat before calling AssembleSystMat" <<endl;
    exit(0);
   }
   
   
//    if(TDatabase::TimeDB->CURRENTTIME>4.0)
    {  
     // before assemble the mesh velocity has to be calculated     
     this->GetMeshVelo(FALSE, 1);
    }
 
 
   for(i=Start_Level;i<N_Levels;i++)
    {    
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();
      N_Active_Current  = U_Space[i]->GetActiveBound();     
      N_DirichletDof = N_U_Current - N_Active_Current;
      
      // reset the nonliner matrices
      AMatAssembleNonLinear[i]->Reset();
      
      // assemble the nonlinear matrix */      
      AMatAssembleNonLinear[i]->Assemble3D();      
      
       // apply upwind disc
      if( (Disctype==UPWIND) && (TDatabase::ParamDB->PROBLEM_TYPE!=STOKES) )
       {
        fefct[0] = Velocity[i]->GetComponent(0);
        fefct[1] = Velocity[i]->GetComponent(1);
        fefct[2] = Velocity[i]->GetComponent(2);

        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes3D(SqmatrixA11[i], fefct[0], fefct[1], fefct[2]);	    
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes3D(SqmatrixA11[i], fefct[0], fefct[1], fefct[2]);
            UpwindForNavierStokes3D(SqmatrixA22[i], fefct[0], fefct[1], fefct[2]);
            UpwindForNavierStokes3D(SqmatrixA33[i], fefct[0], fefct[1], fefct[2]);    
          break;
         }                        // endswitch
       }                          // endif     
        
       // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      { 
       // 	AMatRhsAssemble[i]->AssembleNavierSlip();
            
        fesp[0] =  U_Space[i];
        N_FESpaces = 1;
        N_SquareMatrices = 3;
        N_RectMatrices = 0;
        N_Rhs = 3;

        SQMATRICES[0] = SqmatrixA11[i];
        SQMATRICES[1] = SqmatrixA22[i];
        SQMATRICES[2] = SqmatrixA33[i];

        RHSs[0] = RhsArray[i];
        RHSs[1] = RhsArray[i] + N_U_Current;
        RHSs[2] = RhsArray[i] + 2*N_U_Current;

        fesprhs[0] =  U_Space[i];
        fesprhs[1] =  U_Space[i];
        fesprhs[2] =  U_Space[i];
  
        Assemble3DSlipBC(N_FESpaces, fesp,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, NULL,
                         N_Rhs, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux[i]);
       }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=         
          
    if (Disctype == VMS_PROJECTION)
    {     
     if (i==N_Levels-1)
      {    
//        LumpMassMatrixToDiag(MatricesL[i]);
      
       SQMATRICES[0] = SqmatrixA11[i];
       SQMATRICES[1] = SqmatrixA12[i];
       SQMATRICES[2] = SqmatrixA13[i];
       SQMATRICES[3] = SqmatrixA21[i];
       SQMATRICES[4] = SqmatrixA22[i];
       SQMATRICES[5] = SqmatrixA23[i];
       SQMATRICES[6] = SqmatrixA31[i];
       SQMATRICES[7] = SqmatrixA32[i];
       SQMATRICES[8] = SqmatrixA33[i];
       SQMATRICES[9] =  MatricesL[i];
       MATRICES[0] = Matrices_tilde_G11[i];
       MATRICES[1] = Matrices_tilde_G22[i];
       MATRICES[2] = Matrices_tilde_G33[i];
       MATRICES[3] = Matrices_G11[i];
       MATRICES[4] = Matrices_G22[i];
       MATRICES[5] = Matrices_G33[i];
  
       VMS_ProjectionUpdateMatrices(U_Space[i]->GetN_DegreesOfFreedom(), U_Space[i]->GetActiveBound(), 
                                    Projection_Spaces[i]->GetN_DegreesOfFreedom(), SQMATRICES, MATRICES);
       }      
     }               
   } // for(i=Start_Level;i<N_Levels;i+
      
} //TSystemTNSE3D::AssembleNonLinear(
 
 
 
/** prepare the system matrix */
void TSystemTNSE3D_ALE::AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol)
{
 int i;
 double tau, val = TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT;
 
   if(SystMatAssembled)
   {
    cout << "Restore System mat before calling AssembleSystMat" <<endl;
    exit(0);
   }
 
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  // reset working rhs
  memset(B, 0, N_TotalDOF*SizeOfDouble);    
     
  // old rhs multiplied with current subtime step and theta3 on B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs+N_U, B+N_U);   
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs+2*N_U, B+2*N_U);  
  
  // add rhs from current sub time step to rhs array B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs+N_U, B+N_U);   
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs+2*N_U, B+2*N_U);     
  
  // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A 
  // defect = M * sol
  // B:= B + defect 
  memset(defect, 0, N_TotalDOF*SizeOfDouble);  

  //working rhs needed only at finnest level
  i = N_Levels-1;

  switch(NSEType)
    {
     case 1:
     case 2:
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], -tau*TDatabase::TimeDB->THETA2);          

         MatVectActive(SqmatrixM11[i], sol, defect);
         MatVectActive(SqmatrixM11[i], sol+N_U, defect+N_U);
         MatVectActive(SqmatrixM11[i], sol+2*N_U, defect+2*N_U); 
         Daxpy(N_Active, 1, defect, B);
         Daxpy(N_Active, 1, defect+N_U, B+N_U);
         Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);
     break;
     
     case 3:
     case 4:
        MatAdd(SqmatrixM11[i], SqmatrixA11[i], - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqmatrixM12[i], SqmatrixA12[i], - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqmatrixM13[i], SqmatrixA13[i], - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqmatrixM21[i], SqmatrixA21[i], - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqmatrixM22[i], SqmatrixA22[i], - tau*TDatabase::TimeDB->THETA2); 
        MatAdd(SqmatrixM23[i], SqmatrixA23[i], - tau*TDatabase::TimeDB->THETA2); 
        MatAdd(SqmatrixM31[i], SqmatrixA31[i], - tau*TDatabase::TimeDB->THETA2);
        MatAdd(SqmatrixM32[i], SqmatrixA32[i], - tau*TDatabase::TimeDB->THETA2); 
        MatAdd(SqmatrixM33[i], SqmatrixA33[i], - tau*TDatabase::TimeDB->THETA2); 
  
        MatVectActive(SqmatrixM11[i], sol, defect);
        Daxpy(N_Active, 1, defect, B);
        MatVectActive(SqmatrixM12[i], sol+N_U, defect);
        Daxpy(N_Active, 1, defect, B);
        MatVectActive(SqmatrixM13[i], sol+2*N_U, defect);
        Daxpy(N_Active, 1, defect, B);
        MatVectActive(SqmatrixM21[i], sol, defect+N_U);
        Daxpy(N_Active, 1, defect+N_U, B+N_U);
        MatVectActive(SqmatrixM22[i], sol+N_U, defect+N_U);
        Daxpy(N_Active, 1, defect+N_U, B+N_U);
        MatVectActive(SqmatrixM23[i], sol+2*N_U, defect+N_U);
        Daxpy(N_Active, 1, defect+N_U, B+N_U);  
        MatVectActive(SqmatrixM31[i], sol, defect+2*N_U);
        Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);
        MatVectActive(SqmatrixM32[i], sol+N_U, defect+2*N_U);
        Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);
        MatVectActive(SqmatrixM33[i], sol+2*N_U, defect+2*N_U);
        Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);   
     break;     
    } // switch(NSETyp
  
   // set rhs for Dirichlet Values
   memcpy(B+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
   memcpy(B+(N_U+N_Active), rhs+(N_U+N_Active), N_DirichletDof*SizeOfDouble); 
   memcpy(B+(2*N_U+N_Active), rhs+(2*N_U+N_Active), N_DirichletDof*SizeOfDouble);     

  //     OutPut("B Dot: "<< Ddot(3*N_U, B, B) << endl); 
   
   // generate the system matrix 
   for(i=Start_Level;i<N_Levels;i++)  
    {  
      if(i==(N_Levels-1)) 
       { gamma = - tau*TDatabase::TimeDB->THETA2; }
      else
       { gamma = 0.; }
      
      switch(NSEType)
      {
       case 1:
       case 2:
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], -gamma+tau*TDatabase::TimeDB->THETA1);    
       break;

       case 3:
       case 4:    
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], -gamma+tau*TDatabase::TimeDB->THETA1);
         MatAdd(SqmatrixM12[i], SqmatrixA12[i], -gamma+tau*TDatabase::TimeDB->THETA1);
         MatAdd(SqmatrixM13[i], SqmatrixA13[i], -gamma+tau*TDatabase::TimeDB->THETA1);
         MatAdd(SqmatrixM21[i], SqmatrixA21[i], -gamma+tau*TDatabase::TimeDB->THETA1);
         MatAdd(SqmatrixM22[i], SqmatrixA22[i], -gamma+tau*TDatabase::TimeDB->THETA1); 
         MatAdd(SqmatrixM23[i], SqmatrixA23[i], -gamma+tau*TDatabase::TimeDB->THETA1); 
         MatAdd(SqmatrixM31[i], SqmatrixA31[i], -gamma+tau*TDatabase::TimeDB->THETA1);
         MatAdd(SqmatrixM32[i], SqmatrixA32[i], -gamma+tau*TDatabase::TimeDB->THETA1); 
         MatAdd(SqmatrixM33[i], SqmatrixA33[i], -gamma+tau*TDatabase::TimeDB->THETA1);         
       break;
      } // switch(NSETyp    
      
      
      switch(NSEType)
      {
       case 1:
       case 3:
         if(scale != 1.0)
          {
           Dscal(MatrixB1[i]->GetN_Entries(), scale, MatrixB1[i]->GetEntries());
           Dscal(MatrixB2[i]->GetN_Entries(), scale, MatrixB2[i]->GetEntries());
           Dscal(MatrixB3[i]->GetN_Entries(), scale, MatrixB3[i]->GetEntries());
          }
       break;
       
       case 2:
       case 4:   
         if(scale != 1.0)
          { 
           Dscal(MatrixB1T[i]->GetN_Entries(), scale, MatrixB1T[i]->GetEntries());
           Dscal(MatrixB2T[i]->GetN_Entries(), scale, MatrixB2T[i]->GetEntries());
           Dscal(MatrixB3T[i]->GetN_Entries(), scale, MatrixB3T[i]->GetEntries());

           // scale divergence constraint
           if(val>0) 
            {
             Dscal(MatrixB1[i]->GetN_Entries(), val*scale, MatrixB1[i]->GetEntries());
             Dscal(MatrixB2[i]->GetN_Entries(), val*scale, MatrixB2[i]->GetEntries());
             Dscal(MatrixB3[i]->GetN_Entries(), val*scale, MatrixB3[i]->GetEntries()); 
            }
          } //  if(scale != 1.0) 
       break; 
       
      } // switch(NSETyp    
    } //  for(i=Start_Level;i<N_Levels-1;i++)  
  
   gamma = tau*TDatabase::TimeDB->THETA1;      
   SystMatAssembled  = TRUE;   
      
      
} // TSystemTNSE3D::AssembleSystMat(do

void TSystemTNSE3D_ALE::GetResidual(double *sol, double &impuls_residual, double &residual)
{
#ifdef _MPI   
  int i,j,rank, *master = ParComm_U[N_Levels-1]->GetMaster();  ;
  double residual_scalar = 0.0, sum =0.;
 
   MPI_Comm_rank(Comm, &rank);
#endif
  
     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixM11[N_Levels-1];
          MATRICES[0] = MatrixB1[N_Levels-1];
          MATRICES[1] = MatrixB2[N_Levels-1];
          MATRICES[2] = MatrixB3[N_Levels-1];
        break;

        case 2:
          SQMATRICES[0] = SqmatrixM11[N_Levels-1];
          MATRICES[0] = MatrixB1[N_Levels-1];
          MATRICES[1] = MatrixB2[N_Levels-1];
          MATRICES[2] = MatrixB3[N_Levels-1];
          MATRICES[3] = MatrixB1T[N_Levels-1];
          MATRICES[4] = MatrixB2T[N_Levels-1];
          MATRICES[5] = MatrixB3T[N_Levels-1];
        break;

        case 3:
          SQMATRICES[0] = SqmatrixM11[N_Levels-1];
          SQMATRICES[1] = SqmatrixM12[N_Levels-1];
          SQMATRICES[2] = SqmatrixM13[N_Levels-1];	  
          SQMATRICES[3] = SqmatrixM21[N_Levels-1];
          SQMATRICES[4] = SqmatrixM22[N_Levels-1];
          SQMATRICES[5] = SqmatrixM23[N_Levels-1]; 
          SQMATRICES[6] = SqmatrixM31[N_Levels-1];
          SQMATRICES[7] = SqmatrixM32[N_Levels-1];
          SQMATRICES[8] = SqmatrixM33[N_Levels-1];  

          MATRICES[0] = MatrixB1[N_Levels-1];
          MATRICES[1] = MatrixB2[N_Levels-1];
          MATRICES[2] = MatrixB3[N_Levels-1];
        break;

        case 4:
          SQMATRICES[0] = SqmatrixM11[N_Levels-1];
          SQMATRICES[1] = SqmatrixM12[N_Levels-1];
          SQMATRICES[2] = SqmatrixM13[N_Levels-1];	  
          SQMATRICES[3] = SqmatrixM21[N_Levels-1];
          SQMATRICES[4] = SqmatrixM22[N_Levels-1];
          SQMATRICES[5] = SqmatrixM23[N_Levels-1]; 
          SQMATRICES[6] = SqmatrixM31[N_Levels-1];
          SQMATRICES[7] = SqmatrixM32[N_Levels-1];
          SQMATRICES[8] = SqmatrixM33[N_Levels-1];  
          MATRICES[0] = MatrixB1[N_Levels-1];
          MATRICES[1] = MatrixB2[N_Levels-1];
          MATRICES[2] = MatrixB3[N_Levels-1];
          MATRICES[3] = MatrixB1T[N_Levels-1];
          MATRICES[4] = MatrixB2T[N_Levels-1];
          MATRICES[5] = MatrixB3T[N_Levels-1];
        break;
      } //  switch(NSEType)  

#ifdef _MPI      
    ParComm_U[N_Levels-1]->CommUpdate(sol);
    ParComm_P[N_Levels-1]->CommUpdate(sol+3*N_U);
#endif

    memset(defect, 0, N_TotalDOF*SizeOfDouble);  
    
    Defect(sqmatrices, matrices, sol, B, defect); 
  
     //correction due to L^2_O Pressure space 
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector3D(defect+3*N_U, N_P,  TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
     
#ifdef _MPI   
   for(i=0;i<N_U;i++)
   {
     if(master[i]!=rank)    continue;
      
      residual_scalar += defect[i      ]*defect[i      ];
      residual_scalar += defect[i+  N_U]*defect[i+  N_U];
      residual_scalar += defect[i+2*N_U]*defect[i+2*N_U];

    }

    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
    impuls_residual = (sum);

    master = ParComm_P[N_Levels-1]->GetMaster();
    for(i=0;i<N_P;i++)
    {
      if(master[i]!=rank)    continue;
      
      residual_scalar += defect[i+3*N_U]*defect[i+3*N_U];

    }
    
    sum = 0;
    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
    residual = (sum);

#else
    impuls_residual  =  Ddot(3*N_U, defect, defect);
    residual         =  Ddot(N_TotalDOF, defect, defect); 

#endif
} // TSystemTNSE3D::GetResidual


void TSystemTNSE3D_ALE::Solve(double *sol)
{
 int i,j,rank,*master, N_LinIter=0; 
 double summ = 0., residual,residual_scalar = 0., sum =0.;
 
    switch(SOLVER)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:

          if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
           {
            memcpy(Itmethod_sol, sol, N_TotalDOF*SizeOfDouble);
            memcpy(Itmethod_rhs, B, N_TotalDOF*SizeOfDouble);
           }
          // solve the linear system
          N_LinIter += Itmethod->Iterate(sqmatrices, matrices, Itmethod_sol, Itmethod_rhs);

          if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
           {
            memcpy(sol, Itmethod_sol, N_TotalDOF*SizeOfDouble);
            memcpy(B, Itmethod_rhs, N_TotalDOF*SizeOfDouble);
           }

      break;

      case DIRECT:
 
        switch(NSEType)
         {
          case 1:
           cout << "Solver not included for NSTYPE 1 in this version" <<endl;
            cout << "try NSTYPE 2 or 4 " <<endl;   
	    exit(0);
          break;

          case 2:    
#ifdef _MPI
             DS->Solve(sol, B, true);    
#else	        
             DirectSolver(SqmatrixM11[N_Levels-1], MatrixB1T[N_Levels-1], MatrixB2T[N_Levels-1], MatrixB3T[N_Levels-1],
                          MatrixB1[N_Levels-1], MatrixB2[N_Levels-1], MatrixB3[N_Levels-1], B, sol);

#endif    
          break;

          case 3:
           cout << "Solver not included for NSTYPE 3 in this version" <<endl;
            cout << "try NSTYPE 2 or 4 " <<endl;   
            exit(0);
          break;

          case 4:
#ifdef _MPI
        	DS->Solve(sol, B, true);
#endif
	
#ifdef _OMPONLY
	       if(TDatabase::ParamDB->DSType == 1)
	         DS->Solve(sol, B, true);
	       else{
	         OutPut("Select Proper Solver" << endl);
	         exit(0);
	       }
#endif

#ifdef _SEQ
 #ifdef _SMPI

	  P_DS->Solve(sol,B,true);

#else
             DirectSolver(SqmatrixM11[N_Levels-1], SqmatrixM12[N_Levels-1], SqmatrixM13[N_Levels-1], 
                          SqmatrixM21[N_Levels-1], SqmatrixM22[N_Levels-1], SqmatrixM23[N_Levels-1],  
                          SqmatrixM31[N_Levels-1], SqmatrixM32[N_Levels-1], SqmatrixM33[N_Levels-1],  
                          MatrixB1T[N_Levels-1], MatrixB2T[N_Levels-1], MatrixB3T[N_Levels-1],
                          MatrixB1[N_Levels-1], MatrixB2[N_Levels-1], MatrixB3[N_Levels-1], B, sol,0);

#endif

#endif
          break;
         } //  switch(NSEType) 

      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
} // TSystemTNSE3D::Solve(do


void TSystemTNSE3D_ALE::RestoreMassMat()
{
  int i;
  
   if(SystMatAssembled)
   {
    // restore the mass matrix
    for(i=Start_Level;i<N_Levels;i++)  
     {  
      switch(NSEType)
      {
       case 1:
       case 2:
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], -gamma);    
       break;

       case 3:
       case 4:    
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], -gamma);
         MatAdd(SqmatrixM22[i], SqmatrixA22[i], -gamma); 
         MatAdd(SqmatrixM33[i], SqmatrixA33[i], -gamma);   

        //these mass mattrices are zero
         SqmatrixM12[i]->Reset();
         SqmatrixM13[i]->Reset();  
         SqmatrixM21[i]->Reset();
         SqmatrixM23[i]->Reset(); 
         SqmatrixM31[i]->Reset();
         SqmatrixM32[i]->Reset();
        break;
       } // switch(NSETyp     
      } //for(i=Start_Level;i<N_Levels-1;i++)  
    
    gamma = 0.;      
    SystMatAssembled  = FALSE;  
   }
  else
   {
    cout << "System is not assembled to restore !! " <<endl;
   } 
  
} // TSystemTNSE3D::RestoreMassMatNonLinear()

void TSystemTNSE3D_ALE::RestoreMassMatNonLinear()
{
  int i;
  
   if(SystMatAssembled)
   {
    // restore the mass matrix
    for(i=Start_Level;i<N_Levels;i++)  
     {  
      switch(NSEType)
      {
       case 1:
       case 2:
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], -gamma);    
       break;

       case 3:
       case 4:    
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], -gamma);
         MatAdd(SqmatrixM22[i], SqmatrixA22[i], -gamma); 
         MatAdd(SqmatrixM33[i], SqmatrixA33[i], -gamma);  
         
         if (Disctype == VMS_PROJECTION)
          {     
           MatAdd(SqmatrixM12[i], SqmatrixA12[i], -gamma);
           MatAdd(SqmatrixM13[i], SqmatrixA13[i], -gamma); 
           
           MatAdd(SqmatrixM21[i], SqmatrixA21[i], -gamma);         
           MatAdd(SqmatrixM23[i], SqmatrixA23[i], -gamma);
           
           MatAdd(SqmatrixM31[i], SqmatrixA31[i], -gamma); 
           MatAdd(SqmatrixM32[i], SqmatrixA32[i], -gamma);  
          }
        break;
       } // switch(NSETyp     
      } //for(i=Start_Level;i<N_Levels-1;i++)  
    
    gamma = 0.;      
    SystMatAssembled  = FALSE;  
   }
  else
   {
    cout << "System is not assembled to restore !! " <<endl;
   } 
  
} // TSystemTNSE3D::RestoreMassMatNonLinear()

void TSystemTNSE3D_ALE::AssembleSystMatNonLinear()
{
  int i;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  
   if(SystMatAssembled)
   {
    cout << "Restore System mat before calling AssembleSystMat" <<endl;
    exit(0);
   }
   
   gamma = tau*TDatabase::TimeDB->THETA1; 
  
    // generate the system matrix 
    for(i=Start_Level;i<N_Levels;i++)  
     {  
      switch(NSEType)
      {
       case 1:
       case 2:
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], gamma);    
       break;

       case 3:
       case 4:    
         MatAdd(SqmatrixM11[i], SqmatrixA11[i], gamma);
         MatAdd(SqmatrixM22[i], SqmatrixA22[i], gamma); 
         MatAdd(SqmatrixM33[i], SqmatrixA33[i], gamma);   
         
         if (Disctype == VMS_PROJECTION)
          {     
           MatAdd(SqmatrixM12[i], SqmatrixA12[i], gamma);
           MatAdd(SqmatrixM13[i], SqmatrixA13[i], gamma); 
           
           MatAdd(SqmatrixM21[i], SqmatrixA21[i], gamma);         
           MatAdd(SqmatrixM23[i], SqmatrixA23[i], gamma);
           
           MatAdd(SqmatrixM31[i], SqmatrixA31[i], gamma); 
           MatAdd(SqmatrixM32[i], SqmatrixA32[i], gamma);  
          }
          
        break;
       } // switch(NSETyp     
      }  //  for(i=Start_Level;i<N_Levels-1;i++)  
    
    SystMatAssembled  = TRUE;  
} // TSystemTNSE3D::AssembleSystMatNonLinear()


void TSystemTNSE3D_ALE::MeasureTNSEErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2,  DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP, double *AllErrors)
{
  int i;
  double errors[9];
  
   MultiIndex3D TimeNSAllDerivatives[4] = {D000, D100, D010, D001}; 
   
   
      //error at finnest level
      i = N_Levels -1;
   
      fesp[0] =  U_Space[i];
      fesp[1] =  P_Space[i];
      
      fefct[0] = Velocity[i]->GetComponent(0);
      fefct[1] = Velocity[i]->GetComponent(1);
      fefct[2] = Velocity[i]->GetComponent(2);    
      
      if(NSEaux_error[i]==NULL)
       {
        NSEaux_error[i] =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo, TimeNSN_ParamFctVelo, TimeNSN_FEValuesVelo,
                                        fesp, fefct, TimeNSFctVelo, TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo, TimeNSN_ParamsVelo, TimeNSBeginParamVelo); 
       }

      // L2: error[0], H1-semi: error[1]
      // errors in first velocity component
      fefct[0]->GetErrors(ExactU1, 4, TimeNSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[i], 1, fesp, errors);    
      
      // errors in second velocity component
      fefct[1]->GetErrors(ExactU2, 4, TimeNSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[i], 1, fesp, errors+2);  
                              
      // errors in third velocity component
      fefct[2]->GetErrors(ExactU3, 4, TimeNSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[i], 1, fesp, errors+4);  
      
      // errors in pressure
      Pressure[i]->GetErrors(ExactP, 4, TimeNSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[i], 1,  P_Space+i, errors+6);       
     
      // "L2(u),  H1-semi(u), "L2(p),  H1-semi(p)
       AllErrors[0] = sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4]);
       AllErrors[1] = sqrt(errors[1]*errors[1]+errors[3]*errors[3]+errors[5]*errors[5]); 
       AllErrors[2] = errors[6];
       AllErrors[3] = errors[7];  
     
       // L^infty(0,t,L^2(u)) & time
       if(AllErrors[0] > AllErrors[5])
       {
        AllErrors[4]  =  TDatabase::TimeDB->CURRENTTIME;	 
        AllErrors[5]  = AllErrors[0];
       }       
              
      // error in L^2(0,t,L^2)
      AllErrors[6] += ( (errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4]) + olderror_l_2_l_2u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror_l_2_l_2u = errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4];    
       
} //TSystemTNSE3D::MeasureTNSEErrors(Do



void printall(TSquareMatrix * A11,
//    TSquareMatrix * A12, TSquareMatrix * A13,TSquareMatrix * A21,
   TSquareMatrix * A22,
//    TSquareMatrix * A23, TSquareMatrix * A31, TSquareMatrix * A32,
   TSquareMatrix * A33
//          TSquareMatrix * L,
//                           TMatrix* B1T, TMatrix* B2T, TMatrix* B3T,
//                           TMatrix* B1, TMatrix* B2, TMatrix* B3
                    )
      {
double * entries;
std::ofstream myfile;
myfile.open("entries/A11.txt");
entries = A11->GetEntries();
for(int ii =0; ii< A11->GetN_Entries() ; ii++ )
 myfile << " " << entries[ii] << endl;
myfile.close();

// myfile.open("entries/A12.txt");
// entries = A12->GetEntries();
// for(int ii =0; ii< A12->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/A13.txt");
// entries = A13->GetEntries();
// for(int ii =0; ii< A13->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/A21.txt");
// entries = A21->GetEntries();
// for(int ii =0; ii< A21->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
myfile.open("entries/A22.txt");
entries = A22->GetEntries();
for(int ii =0; ii< A22->GetN_Entries() ; ii++ )
 myfile << " " << entries[ii] << endl;
myfile.close();
// myfile.open("entries/A23.txt");
// entries = A23->GetEntries();
// for(int ii =0; ii< A23->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/A31.txt");
// entries = A31->GetEntries();
// for(int ii =0; ii< A31->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/A32.txt");
// entries = A32->GetEntries();
// for(int ii =0; ii< A32->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
myfile.open("entries/A33.txt");
entries = A33->GetEntries();
for(int ii =0; ii< A33->GetN_Entries() ; ii++ )
 myfile << " " << entries[ii] << endl;
myfile.close();
// myfile.open("L.txt");
// entries = L->GetEntries();
// for(int ii =0; ii< L->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
// myfile.open("entries/tilde_G11.txt");
// entries = B1T->GetEntries();
// for(int ii =0; ii< B1T->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
// 
// myfile.open("entries/tilde_G22.txt");
// entries = B2T->GetEntries();
// for(int ii =0; ii< B2T->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/tilde_G33.txt");
// entries = B3T->GetEntries();
// for(int ii =0; ii< B3T->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/G11.txt");
// entries = B1->GetEntries();
// for(int ii =0; ii< B1->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/G22.txt");
// entries = B2->GetEntries();
// for(int ii =0; ii< B2->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
//
// myfile.open("entries/G33.txt");
// entries = B3->GetEntries();
// for(int ii =0; ii< B3->GetN_Entries() ; ii++ )
//  myfile << " " << entries[ii] << endl;
// myfile.close();
 }

void TSystemTNSE3D::printall_matrix()
 {
  cout << "started printing" << endl;
     
   printall( SqmatrixA11[N_Levels-1], SqmatrixA22[N_Levels-1], SqmatrixA33[N_Levels-1]);
      
    cout << "finished printing" << endl;
     
} 

void TSystemTNSE3D_ALE::CheckAllMat()
{
        char key[14];
      sprintf(key,"SqmatrixA11");
      DOF_stats((TMatrix3D*)SqmatrixA11[N_Levels-1],'U',0,key);
      sprintf(key,"SqmatrixA12");
      DOF_stats((TMatrix3D*)SqmatrixA12[N_Levels-1],'U',1,key);
      sprintf(key,"SqmatrixA13");
      DOF_stats((TMatrix3D*)SqmatrixA13[N_Levels-1],'U',2,key);
      sprintf(key,"SqmatrixA21");
      DOF_stats((TMatrix3D*)SqmatrixA21[N_Levels-1],'U',0,key);
      sprintf(key,"SqmatrixA22");
      DOF_stats((TMatrix3D*)SqmatrixA22[N_Levels-1],'U',1,key);
      sprintf(key,"SqmatrixA23");
      DOF_stats((TMatrix3D*)SqmatrixA23[N_Levels-1],'U',2,key);
      sprintf(key,"SqmatrixA31");
      DOF_stats((TMatrix3D*)SqmatrixA31[N_Levels-1],'U',0,key);
      sprintf(key,"SqmatrixA32");
      DOF_stats((TMatrix3D*)SqmatrixA32[N_Levels-1],'U',1,key);
      sprintf(key,"SqmatrixA33");
      DOF_stats((TMatrix3D*)SqmatrixA33[N_Levels-1],'U',2,key);
      sprintf(key,"MatrixB1");
      DOF_stats(MatrixB1[N_Levels-1],'P',0,key);
      sprintf(key,"MatrixB2");
      DOF_stats(MatrixB2[N_Levels-1],'P',1,key);
      sprintf(key,"MatrixB3");
      DOF_stats(MatrixB3[N_Levels-1],'P',2,key);
      sprintf(key,"MatrixB1T");
      DOF_stats(MatrixB1T[N_Levels-1],'U',3,key);
      sprintf(key,"MatrixB2T");
      DOF_stats(MatrixB2T[N_Levels-1],'U',3,key);
      sprintf(key,"MatrixB3T");
      DOF_stats(MatrixB3T[N_Levels-1],'U',3,key);
}
void TSystemTNSE3D_ALE::DOF_stats(TMatrix3D* MatB, char M , int k,char * name)
{

int n_row;
double cou;
int * row_ptr;
int * Kcol;
int *Vect,*tVect;
double* nnz;
int nnzs=0;
int rank,begin,end;
double* range;
char str[40];
double *col_nnz;
double d_x,d_y,d_z;
std::ofstream myfile;
int* GlNr;

int nrs;


double * solu = SolArray[N_Levels-1];
TFESpace3D* fspace;

if(M == 'U')
  fspace = U_Space[N_Levels-1];
else
  fspace = P_Space[N_Levels-1];

//nrs = fspace->GetN_DegreesOfFreedom();

//nrs = Pfespace->GetN_DegreesOfFreedom();

#ifdef _MPI
int *Master;

if(M=='U')
 Master = ParComm_U[N_Levels-1]->GetMaster();
 else
  Master = ParComm_P[N_Levels-1]->GetMaster();
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 //GlNr = ParComm_U[N_Levels-1]->Get_Local2Global();
#endif

 #ifdef _MPI
 sprintf(str,"All_Matrices/%s_%d.txt",name,rank);
#else
 sprintf(str,"All_Matrices/%s.txt",name);
#endif

myfile.open(str);
                                   

// //TMatrix3D* MatB =  MatrixB1[N_Levels-1];
// TSquareMatrix* MatB = SqmatrixA12[N_Levels-1];
// //TMatrix3D* MatB =  MatrixB1T[N_Levels-1];

int * rowp =  MatB->GetRowPtr();
int *colp = MatB->GetKCol();
double * entp = MatB->GetEntries();

int r_p =0;
double *row_ent;
int *d_id;
double summer = 0.0;
double tsummer = 0.0;
int max=0;

nrs = fspace->GetN_DegreesOfFreedom();



for(int ii=0; ii< nrs; ii++)
{
#ifdef _MPI
        //if(MasterP[ii]==rank)
         if(Master[ii]==rank)
#endif
        {
             // Pfespace->GetDOFPosition(ii,d_x,d_y,d_z);
              fspace->GetDOFPosition(ii,d_x,d_y,d_z);

              r_p = rowp[ii + 1] - rowp[ii];

              if(max<r_p)
                max =r_p;


        }
}
row_ent = new double[max];
d_id = new int[max];
//cout << "Max here at " << str << " is :: " << max << endl;
for(int ii=0; ii< nrs; ii++)
{
#ifdef _MPI
        //if(MasterP[ii]==rank)
         if(Master[ii]==rank)
#endif
        {
             // Pfespace->GetDOFPosition(ii,d_x,d_y,d_z);
              fspace->GetDOFPosition(ii,d_x,d_y,d_z);

// #ifdef _MPI
              //myfile << " " << solu[ii]  << " " << GlNr[ii] << " " << d_x << " "<< d_y << "   "<< d_z << endl;

              myfile   << " " << d_x << " "<< d_y << "   "<< d_z;
              r_p = 0;
              for(int jj=rowp[ii]; jj < rowp[ii + 1]; jj++)
              {
                row_ent[r_p] = entp[jj] ;
                d_id[r_p] = colp[jj];
                 summer += ( (solu[colp[jj]+ k*N_U]* entp[jj]))*( (solu[colp[jj]+ k*N_U]* entp[jj]));

                r_p++;
              }


              mergesort(row_ent,0,r_p-1,d_id,max);

              for(int jj=0; jj < r_p; jj++)
                myfile << " " <<  row_ent[jj] << " s:" << solu[d_id[jj] + k*N_U];

              myfile << "" << endl;           

        }


}

#ifdef _MPI
//printf("our parallel output :: %lf rank ::%d \n",summer, rank);
//cout << "our parallel output :: " << summer <<  endl;
   MPI_Allreduce(&summer,&tsummer,1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if(rank ==0)
     cout << "our output :: " << tsummer << endl;
#else
     tsummer = summer;
     cout << "our output :: " << tsummer << endl;
#endif

     delete []d_id;
     delete []row_ent;

myfile.close();

}


void TSystemTNSE3D_ALE::RHS_stats()
{

char str[40];
int rank;
int U_nrs, P_nrs;
TFESpace3D *Uspace,*Pspace;
double d_x,d_y,d_z;
std::ofstream myfile;

#ifdef _MPI
  int *MasterU,*MasterP;

  MasterU = ParComm_U[N_Levels-1]->GetMaster();
  MasterP = ParComm_P[N_Levels-1]->GetMaster();

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#endif

double * solu = SolArray[N_Levels-1];

#ifdef _MPI
 sprintf(str,"RHS/rhs_%d.txt",rank);
#else
 sprintf(str,"RHS/rhs.txt");
#endif

myfile.open(str);



  Uspace = U_Space[N_Levels-1];
  Pspace = P_Space[N_Levels-1];

  U_nrs = Uspace->GetN_DegreesOfFreedom();
  P_nrs = Pspace->GetN_DegreesOfFreedom();

        for(int ii=0; ii< U_nrs; ii++)
          {
          #ifdef _MPI
                  //if(MasterP[ii]==rank)
                  if(MasterU[ii]==rank)
          #endif
                  {
                      // Pfespace->GetDOFPosition(ii,d_x,d_y,d_z);
                        Uspace->GetDOFPosition(ii,d_x,d_y,d_z);
                        myfile << "" << d_x<< " " << d_y << " " << d_z << " " << solu[ii] << " " << solu[ii + N_U]
                        << " " << solu[ii + 2*N_U] << endl;
                  }
          }

        for(int ii=0; ii< P_nrs; ii++)
          {
          #ifdef _MPI
                  //if(MasterP[ii]==rank)
    if(MasterP[ii]==rank)
          #endif
                  {
                      // Pfespace->GetDOFPosition(ii,d_x,d_y,d_z);
                        Pspace->GetDOFPosition(ii,d_x,d_y,d_z);
                        myfile << "" << d_x<< " " << d_y << " " << d_z << " " << solu[ii + 3*N_U]  << endl;
                  }
          }

myfile.close();

}
