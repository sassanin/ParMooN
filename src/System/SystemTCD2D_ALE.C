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
* @brief     source file for TSystemTCD2D_ALE
* @author    Sashikumaar Ganesan
* @date      19.02.15
* @History 
 ************************************************************************  */
#ifdef __2D__

#include <Database.h>
#include <SystemTCD2D_ALE.h>
#include <SystemTCD2D.h>
#include <FEVectFunct2D.h>
// #include <TimeConvDiff2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
#include <LinAlg.h>
#include <FreeSurface2D.h>
#include <CDSystemTimeDG_1.h>

TSystemTCD2D_ALE::TSystemTCD2D_ALE(TFESpace2D *fespace, int disctype, int solver, TFESpace2D *gridFESpace, TFEVectFunct2D *MeshVelocity, 
                                                       bool conservativeale) : TSystemTCD2D(fespace,  disctype, solver)
{
  char WString[] = "w";  
  TFESpace2D *fesp[1];
    
  
  /** old M mass matrix */
  sqmatrixM_old = new TSquareMatrix2D(sqstructure);  
  
  GridFESpace = gridFESpace;
  N_GridDOFs = gridFESpace->GetN_DegreesOfFreedom();
  N_GridActive = gridFESpace->GetActiveBound();
  
  // grid 
  SquareStructureG= new TSquareStructure2D(GridFESpace); 
  SquareStructureG->Sort();
   
  // for mesh
  SqmatrixG11 = new TSquareMatrix2D(SquareStructureG); // G11
  SqmatrixG12 = new TSquareMatrix2D(SquareStructureG); // G12
  SqmatrixG21 = new TSquareMatrix2D(SquareStructureG); // G21
  SqmatrixG22 = new TSquareMatrix2D(SquareStructureG); // G22
  
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
    
  MeshVeloFct[0] = MeshVelocity->GetComponent(0);
  MeshVeloFct[1] = MeshVelocity->GetComponent(1);
  MeshVelo =  MeshVelocity->GetValues();
  
  gridpos = new double[2*N_GridDOFs];
  gridpos_old = new double[2*N_GridDOFs];   
  gridpos_ref = new double[2*N_GridDOFs];   
  griddisp = new double[2*N_GridDOFs];   
  
  
   memset(gridpos, 0, 2*N_GridDOFs*SizeOfDouble);
   GridPos = new TFEVectFunct2D(GridFESpace, WString, WString, gridpos, N_GridDOFs, 2);  
   GridPos->GridToData();
   
   GridRhs = new double[2*N_GridDOFs];
   
   RefGridPos = new TFEVectFunct2D(GridFESpace, WString, WString, gridpos_ref, N_GridDOFs, 2);     
 
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
   memcpy(gridpos_ref, gridpos, 2*N_GridDOFs*SizeOfDouble); 
   
   Aux_ALE = NULL;
   SolveLinearElastic = TRUE;
   CONSERVATIVEALE = conservativeale;  
   NeedInterMassMat = TRUE;
     
   if(!CONSERVATIVEALE && TDatabase::TimeDB->TIME_DISC==1)
     NeedInterMassMat = FALSE;  
 
  /**matrices at time quad points for dG time steppings */
  if(TDatabase::TimeDB->DG_TimeDisc)
   {    
     
     TimeDG->SetALEForm(CONSERVATIVEALE);
      
     switch(TDatabase::TimeDB->DG_Order)
      {    
       case 1: // dG(1)
     
         /** additional rhs for Qp1 */ 
	 rhs_Qp1 = new double[N_DOF];

         /** M mass matrix */
         sqmatrixM_Qp1 = new TSquareMatrix2D(sqstructure);  
         N_Matrices++;   
         /** A stifness matrix */
         sqmatrixA_Qp1 = new TSquareMatrix2D(sqstructure);  
         N_Matrices++; 
 
         // add additional QuadPt matrics to dG
         TimeDG->AddQp1Matrices(sqmatrixM_Qp1, sqmatrixA_Qp1, rhs_Qp1); 
       break;      
 
       default:
            OutPut("time dG order" << TDatabase::TimeDB->DG_Order <<" not yet implemented " << endl);
            exit(4711);;
     }     
   }   
  
    
} // TSystemTCD2D_ALE


void TSystemTCD2D_ALE::Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue, 
                                      CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D *GridBoundCond, BoundValueFunct2D *gridBoundValue,
                                      TAuxParam2D *aux      
                                      )
{
  
  Aux_ALE = aux;
     
  BoundaryConditions[0] =  BoundCond;
  BoundaryValues[0] = BoundValue;
  
  GridBoundValue[0] = gridBoundValue;
  GridBoundaryConditions[0] = GridBoundCond;
  
  TDiscreteForm2D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm2D *DiscreteFormARhs_Galerkin; 
  TDiscreteForm2D *DiscreteFormMRhs_SUPG_NULL;
  TDiscreteForm2D *DiscreteFormMRhs_SUPG;
  TDiscreteForm2D *DiscreteFormARhs_SUPG;
  
  TDiscreteForm2D *DiscreteFormMARhs_Galerkin;
  TDiscreteForm2D *DiscreteFormMatrixMARhs_SUPG;
     
  InitializeDiscreteForms_ScalarMoving(DiscreteFormMARhs_Galerkin, DiscreteFormGrid, DiscreteFormMRhs_SUPG, DiscreteFormMatrixMARhs_SUPG,
                                        BilinearCoeffs, GridBilinearCoeffs);
   
  
  InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormMRhs_SUPG_NULL,
                                  DiscreteFormARhs_SUPG, BilinearCoeffs);
  
    switch(Disctype)
     {
      case GALERKIN:
      case LOCAL_PROJECTION:
//            DiscreteFormARhs = DiscreteFormARhs_Galerkin;
           DiscreteFormMRhs = DiscreteFormMRhs_Galerkin;
           DiscreteFormMARhs = DiscreteFormMARhs_Galerkin;
      break;
      
      case SUPG:
//            DiscreteFormARhs = DiscreteFormARhs_SUPG;
           DiscreteFormMRhs = DiscreteFormMRhs_SUPG;
           DiscreteFormMARhs = DiscreteFormMatrixMARhs_SUPG;          
      break;
      
      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }  
     
    /** dG time steppings */
    if(TDatabase::TimeDB->DG_TimeDisc)
     {
      TimeDG->AddBilinear(BilinearCoeffs);
     }//  if(TDatabase::TimeDB->DG_TimeDisc)
   
} // Init


void TSystemTCD2D_ALE::AddBoundModifyFunction(ModifyBoundCoords *modifyboudary, int *n_MovVert, TVertex **movBoundVert, 
                                              TIsoBoundEdge **free_Joint, double *iso_refX)
 { 
   ModifyBoudary = modifyboudary; 
   
   SolveLinearElastic = TRUE; 
   
   N_MovVert = n_MovVert;
   MovBoundVert = movBoundVert;
   Free_Joint = free_Joint;
   Iso_refX = iso_refX;  
} 


void TSystemTCD2D_ALE::StoreMmat()
 { 
   sqmatrixM_old->Reset(); 
   MatAdd(sqmatrixM_old, sqmatrixM, 1.); 
 
 } // TSystemTCD2D_ALE::StoreMmat()

 
 /** move the bound and mesh to the given time */
void TSystemTCD2D_ALE::MoveMesh(double Currtime)
{
 int i, N_GridBDDOFs;
 
 
 
   if(SolveLinearElastic)
    {
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
     
      memcpy(griddisp, GridRhs, 2*N_GridDOFs*SizeOfDouble);   
      SolveGridEquation(Entries, griddisp, GridRhs, GridKCol, GridRowPtr, N_GridDOFs);
           
      memcpy(gridpos, gridpos_old, 2*N_GridDOFs*SizeOfDouble);
      Daxpy(2*N_GridDOFs, 1., griddisp, gridpos);
    }
   else 
    {
     // domain at the given time  
     for(i=0;i<N_GridDOFs;i++)
       ModifyCoord(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], Currtime);   
    }   
           
   // move the mesh
   GridPos->DataToGrid();    
  
} // MoveMesh




void ModifyCoord_Exact(double x, double y, double &X, double &Y, double t, double &W1, double &W2)
{
 double d;
 
  d = 2. - cos(20.*Pi*t);
 
 X = x*d;
 Y = y*d;
   
 // mesh velocity
 W1 = (  20.*Pi*x*sin(20.*Pi*t))/d;
 W2 = (  20.*Pi*y*sin(20.*Pi*t))/d;
 
//  double tn, tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
 
//  tn =  t - tau ;
  // mesh velocity
//  W1 = -(x/tau)*(cos(20.*Pi*t) + cos(20.*Pi*tn));
//  W2 = -(y/tau)*(cos(20.*Pi*t) + cos(20.*Pi*tn));
//  cout << "test mesh velo " <<endl;
}



void TSystemTCD2D_ALE::GetMeshVelo(double Currtime, double tau, bool MoveMesh)
{
 int i, N_GridBDDOFs;
 
   GridPos->GridToData();   
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);  
      
   if(SolveLinearElastic)
    { 
      // assemble the mesh matrix, that is, the previous time-step domain is the reference domain
      this->AssembleMeshMat();

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
     
      //   memcpy(MeshVelo, GridRhs, 2*N_GridDOFs*SizeOfDouble);   
      memset(MeshVelo, 0, 2*N_GridDOFs*SizeOfDouble); 
      
//       OutPut("MeshVelo " << Ddot(2*N_GridDOFs, GridRhs, GridRhs ) << endl);  
      
      SolveGridEquation(Entries, MeshVelo, GridRhs, GridKCol, GridRowPtr, N_GridDOFs);        
  
      memcpy(gridpos, gridpos_old, 2*N_GridDOFs*SizeOfDouble);
      if(MoveMesh)      
       { Daxpy(2*N_GridDOFs, 1., MeshVelo, gridpos); } // add the displacement
       
      Dscal(2*N_GridDOFs, 1./tau, MeshVelo);   
      // move the mesh back to original pos, if we calculate only the mesh velocity otherwise
      GridPos->DataToGrid();   
    }
   else 
    { 
//        for(i=0;i<N_GridDOFs;i++)
//           ModifyCoord_Exact(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], 
//                             Currtime, MeshVelo[i], MeshVelo[i+N_GridDOFs]); 
 
      //  move velo in current time  
      for(i=0;i<N_GridDOFs;i++)
        ModifyCoord(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], Currtime);   

      //compute mesh velocity
      memset(MeshVelo, 0, 2*N_GridDOFs*SizeOfDouble); 
      memcpy(MeshVelo, gridpos, 2*N_GridDOFs*SizeOfDouble);     
      Daxpy(2*N_GridDOFs, -1., gridpos_old, MeshVelo);        
      Dscal(2*N_GridDOFs, 1./tau, MeshVelo);   
 
       // move the mesh back to original pos, as we calculate only the mesh velocity
      if(MoveMesh)
      { GridPos->DataToGrid(); }      
    }
   

   
} // TSystemTCD2D_ALE::GetMeshVelo
  
  
void TSystemTCD2D_ALE::AssembleMeshMat()
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
  
}

void TSystemTCD2D_ALE::AssembleMRhs(double *sol, double *rhs)
{
  int N_DOF, N_Active, N_SquareMatrices;
  double *RHSs[1];

  TFESpace2D *fesp[2], *ferhs[1];
  TFEFunction2D  *fefct[4];
  TAuxParam2D *aux;
  
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    fesp[1] = GridFESpace;    
    ferhs[0] = FeSpace;
    
    // initialize matrices
//     TSquareMatrix2D *SQMATRICES[2] = { sqmatrixM, NULL };
    SQMATRICES[0] = sqmatrixM;
    SQMATRICES[0]->Reset();

    N_SquareMatrices =1;
    
//     aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); 
    
    if(Disctype == SDFEM)
    {
     N_SquareMatrices = 2;  
     SQMATRICES[1] = sqmatrixS;
     SQMATRICES[1]->Reset();  
    }

    fefct[0] = MeshVeloFct[0];
    fefct[1] = MeshVeloFct[1];
    fefct[2] = MeshVeloFct[0]; //for calculating divergence of w
    fefct[3] = MeshVeloFct[1]; //for calculating divergence of w   

    // assemble
    Assemble2D(2, fesp,
               N_SquareMatrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormMRhs,
               BoundaryConditions,
               BoundaryValues,
               Aux_ALE);
 
    
   // copy Dirichlet values from rhs into sol
   memcpy(sol+N_Active, rhs+N_Active, (N_DOF - N_Active)*SizeOfDouble);  

   if(Disctype==SDFEM && TDatabase::ParamDB->REACTOR_P28==1)
    {     
     MatAdd(sqmatrixM, sqmatrixS, 1.);
    }
     

     //  memset(defect, 0, N_DOF*SizeOfDouble);   
//          MatVectActive(SQMATRICES[0], sol, defect);
//         OutPut("MatM-Sol " << Ddot(N_DOF, defect, defect ) << endl);  
//        OutPut("MatM-Sol " << Ddot(N_DOF, sol, sol ) << endl); 
  
//     cout << " AssembleMRhs " << endl;  
     
} // TSystemMatScalar2D::AssembleMRhs 

void TSystemTCD2D_ALE::AssembleMARhs(double *sol, double *rhs)
{
  int N_DOF, N_Active, N_SquareMatrices;
  double *RHSs[1];

  TFESpace2D *fesp[2], *ferhs[1];
  TFEFunction2D  *fefct[4];
    
    // assemble the mass mat and rhs 
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    fesp[1] = GridFESpace;    
    ferhs[0] = FeSpace;
    
    // initialize matrices
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[1] = sqmatrixM;    
    SQMATRICES[0]->Reset();
    SQMATRICES[1]->Reset();
    
    N_SquareMatrices =2;
    
   if(Disctype == SDFEM)
    {
     N_SquareMatrices = 3;
     SQMATRICES[2] = sqmatrixS;
     SQMATRICES[2]->Reset();  
    }      

    fefct[0] = MeshVeloFct[0];
    fefct[1] = MeshVeloFct[1];
    fefct[2] = MeshVeloFct[0]; //for calculating divergence of w
    fefct[3] = MeshVeloFct[1]; //for calculating divergence of w    

    // assemble
    Assemble2D(2, fesp,
               N_SquareMatrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormMARhs,
               BoundaryConditions,
               BoundaryValues,
               Aux_ALE);
     
      // copy Dirichlet values from rhs into sol
      memcpy(sol+N_Active, rhs+N_Active, (N_DOF - N_Active)*SizeOfDouble);  
      
      
    if(Disctype==SUPG && TDatabase::ParamDB->REACTOR_P28==1)
     {     
      MatAdd(sqmatrixM, sqmatrixS, 1.);
     }
     
//         memset(defect, 0, N_DOF*SizeOfDouble);   
//         MatVectActive(SQMATRICES[1], sol, defect);
//         OutPut("MatM-Sol " << Ddot(N_DOF, defect, defect ) << endl);  
//        OutPut("MatM-Sol " << Ddot(N_DOF, sol, sol ) << endl); 
  
//     cout << " AssembleMARhs " << endl;
//       exit(0);  
  
} // TSystemMatScalar2D::AssembleMARhs 



void TSystemTCD2D_ALE::AssembleALEMat(double *sol, double *rhs, double tau)
{
  double CurrTime;
  
  // assemble M & A Mat in "CurrTime" mesh
  TDatabase::ParamDB->REACTOR_P29=1.e10; // get min \div w
  TDatabase::ParamDB->REACTOR_P30=-1.e10;   // get max \div w 
     
   //w^{n+1}, must be computed before moving the mesh in this time step
   this->GetMeshVelo(TDatabase::TimeDB->CURRENTTIME, tau, FALSE);   
       
  /** dG time steppings */
  if(TDatabase::TimeDB->DG_TimeDisc)
   {
      
    this->StoreMmat();  
     
     switch(TDatabase::TimeDB->DG_Order)
      {    
       case 1: // dG(1)
     
          CurrTime = TDatabase::TimeDB->CURRENTTIME - 2*tau/3.;
	  
// 	  cout << TDatabase::TimeDB->CURRENTTIME << " CurrTime " << CurrTime <<endl;
	  
          //w^{n+1/3}, must be computed before moving the mesh in this time step
          // and move the mesh
//           this->GetMeshVelo(CurrTime, tau/3., FALSE); 
          this->MoveMesh(CurrTime);  
          // assemble M_n+1/3 & A_n+1/3 at current domain
          memset(rhs_Qp1, 0, N_DOF*SizeOfDouble);  
          this->AssembleMARhs(sol, rhs_Qp1);  
          
          // store the current matrices
          sqmatrixM_Qp1->Reset(); 
          MatAdd(sqmatrixM_Qp1, sqmatrixM, 1.); 
          sqmatrixA_Qp1->Reset(); 
          MatAdd(sqmatrixA_Qp1, sqmatrixA, 1.); 

          CurrTime = TDatabase::TimeDB->CURRENTTIME;  
// 	  	  cout << tau/3 << " CurrTime " << CurrTime <<endl;
          //w^{n+1}, must be computed before moving the mesh in this time step
          // and move the mesh
//           this->GetMeshVelo(TDatabase::TimeDB->CURRENTTIME, 2*tau/3., TRUE); 
          // move the mesh  
          this->MoveMesh(CurrTime);    
          // assemble M_n+1 & A_n+1 at current domain
          memset(rhs, 0, N_DOF*SizeOfDouble);  
          this->AssembleMARhs(sol, rhs); 

       break;      
 
       default:
            OutPut("time dG order" << TDatabase::TimeDB->DG_Order <<" not yet implemented " << endl);
            exit(4711);;
        }  //  switch(TDatabase::TimeDB->DG_Order)  

   }
  else // BE, CN
   {   
    if(NeedInterMassMat)
      { CurrTime = TDatabase::TimeDB->CURRENTTIME - tau/2.; }  
     else
      { CurrTime = TDatabase::TimeDB->CURRENTTIME; }    

  
     // move the mesh  
     this->MoveMesh(CurrTime); 
//     cout << "test" <<endl;
//      this->GetMeshVelo(CurrTime = TDatabase::TimeDB->CURRENTTIME, tau  , FALSE);   
     // assemble M_n+1/2 & A_n+1/2 at current domain
     this->AssembleMARhs(sol, rhs);       
    
     if(NeedInterMassMat)
      { 
       CurrTime = TDatabase::TimeDB->CURRENTTIME;
   
       // move the mesh to CurrTime
       this->MoveMesh(CurrTime);
       
       // assemble M_n+1  Mat in "CurrTime" mesh, only M mat so mesh velo is not needed
       this->AssembleMRhs(sol, rhs);
      }
   } // else // BE, CN
  
}



void TSystemTCD2D_ALE::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol)
{
  int N_DOF, N_Active, N_SquareMatrices;
  double tau, fact, CurrTime, Oldtime;
    
  N_DOF = FeSpace->GetN_DegreesOfFreedom();
  memset(defect, 0, N_DOF*SizeOfDouble);  
  
   /** dG time steppings */
  if(TDatabase::TimeDB->DG_TimeDisc)
   {
     // M_n*u_n
    MatVectActive(sqmatrixM_old, oldsol, defect); 
//      OutPut("B  " << Ddot(N_DOF, defect, defect ) << endl);  
    // all sys mat entries will be reset to zero
    TimeDG->ResetSysMat(); 
    
    //pass Mu_old and current rhs
    TimeDG->AssembleALESysMat_Qp1(defect, rhs);
//      cout<<"AssembleALESysMat_Qp1    " <<   endl ;
//     exit(0);
   }
  else
   {    
    N_Active =  FeSpace->GetActiveBound();   
    
    tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  
    memset(B, 0, N_DOF*SizeOfDouble);      
    
    // old rhs multiplied with current subtime step and theta3 on B
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3,  oldrhs, B);    

    // add rhs from current sub time step to rhs array B
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4,  rhs, B);    


    if(CONSERVATIVEALE)        
     {
      MatAdd(sqmatrixM_old, sqmatrixA, - tau*TDatabase::TimeDB->THETA2);
      gamma =0.;  
      
      //defect = M * oldsol
      MatVectActive(sqmatrixM_old, oldsol, defect);
          
      //store the mass matrix for next time step (only for linear problem)
      this->StoreMmat();    
     }
    else
     {
//       CurrTime = TDatabase::TimeDB->CURRENTTIME  ;
//       Oldtime = TDatabase::TimeDB->CURRENTTIME - TDatabase::TimeDB->CURRENTTIMESTEPLENGTH ;
//       fact = ( 2. - cos(20.*Pi*Oldtime))/( 2. - cos(20.*Pi*CurrTime));  // expanding squaredomain
// //       fact *=fact;   
// //       fact = 1.;
//       // cout << "CurrTime " << CurrTime << " " << "Oldtime  " << Oldtime << " "   << "fact " << fact<<endl;
//       Dscal(sqmatrixM->GetN_Entries(), fact, sqmatrixM->GetEntries());    
         
      MatAdd(sqmatrixM, sqmatrixA, - tau*TDatabase::TimeDB->THETA2);
      gamma = -tau*TDatabase::TimeDB->THETA2;  
     
      MatVectActive(sqmatrixM, oldsol, defect);       

/*      MatAdd(sqmatrixM, sqmatrixA, -gamma);
      gamma = 0;
      Dscal(sqmatrixM->GetN_Entries(), 1./fact, sqmatrixM->GetEntries());    */     
     }

     // B:= B + defec    
     Daxpy(N_Active, 1, defect, B);
     
     // set Dirichlet values
     memcpy(B+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
     memcpy(sol+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
     
     //assemble the system matrix
     MatAdd(sqmatrixM, sqmatrixA, -gamma + tau*TDatabase::TimeDB->THETA1);

     memset(defect, 0, N_DOF*SizeOfDouble);  
     MatVectActive(sqmatrixM, oldsol, defect);
    } // else  if(TDatabase::TimeDB->DG_TimeDisc)
        
 
// 	       exit(0);
	       
} // AssembleSystMat
 
void TSystemTCD2D_ALE::Solve(double *sol, double *rhs)
{
  
   if(TDatabase::TimeDB->DG_TimeDisc)
     {     
       // all sys mat entries will be reset to zero
       TimeDG->SoveTimedG(sol); 
     }
    else
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
        DirectSolver(sqmatrixM, B, sol);
       break;      
 
       default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
      }    
     } //   else if(TDatabase::TimeDB->DG_TimeDisc)
     
//            OutPut("Sol " << Ddot(N_DOF, sol, sol) << endl);   
     
}


// double TSystemTCD2D_ALE::GetResidual(double *sol)
// {
//   int N_DOF = FeSpace->GetN_DegreesOfFreedom(); 
//   double residual_scalar;
//   
// //   if(SystMatAssembled)
//    {
//     memset(defect, 0, N_DOF*SizeOfDouble);         
//     ScalarDefect(sqmatrixM, sol, B, defect, residual_scalar);   
//    }
// //   else
//    {
// //     OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
// //     exit(4711);;   
//    }
//    
//    return residual_scalar;
//       
// }
// 
// 

#endif // #ifdef __2D__
