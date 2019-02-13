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
* @brief     source file for TSystemNSE3D
* @author    Sashikumaar Ganesan, 
* @date      27.01.15
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemNSE3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <FEVectFunct3D.h>
#include <AuxParam3D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <NSE3D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Upwind3D.h>
#include <AssembleMat3D.h>

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
// #include <sstream>
// #include <MooNMD_Io.h>

TSystemNSE3D::TSystemNSE3D(int N_levels, TFESpace3D **velocity_fespace, TFESpace3D **presssure_fespace, TFEVectFunct3D **velocity, 
                     TFEFunction3D **pressure, double **sol, double **rhs, int disctype, int nsetype, int solver)
{
  int i, zerostart;
  int profiling = TDatabase::ParamDB->timeprofiling;
 
  /** need it for solver */
  sqmatrices = (TSquareMatrix **)SQMATRICES;
  matrices = (TMatrix **)MATRICES;
  
   if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE) || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
    { N_aux= 4; }
   else
    { N_aux= 2; }
    
  //set number of multigrid levels
  N_Levels = N_levels;
   
  //set the discretization type
  Disctype = disctype;
  
  // NSE type
  NSEType = nsetype;
  
#ifdef _MPI
  Comm = TDatabase::ParamDB->Comm;
  if(!(NSEType==2 || NSEType==4))
  {
    printf("parallel implemented for NSTYPE 2 and 4 only\n");
    MPI_Finalize();
    exit(0);
  }
#endif  
  
  //set the solver type
  SOLVER = solver;
  
  Velocity = velocity;
  Pressure = pressure;
 
  SolArray = sol;
  RhsArray = rhs;
  
  U_Space = velocity_fespace;
  P_Space = presssure_fespace;
  
  N_U = velocity_fespace[N_levels-1]->GetN_DegreesOfFreedom();
  N_P = presssure_fespace[N_levels-1]->GetN_DegreesOfFreedom();
  N_TotalDOF = 3*N_U + N_P;
  
  N_Active =  velocity_fespace[N_levels-1]->GetActiveBound();
  N_DirichletDof = N_U - N_Active;  
 
  sqstructureA = new TSquareStructure3D *[N_levels];
  structureB = new TStructure3D *[N_levels];
  structureBT = new TStructure3D *[N_levels];
  
  SqmatrixA11 = new TSquareMatrix3D*[N_levels];
  SqmatrixA12 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA13 = new TSquareMatrix3D*[N_levels];  
  SqmatrixA21 = new TSquareMatrix3D*[N_levels];
  SqmatrixA22 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA23 = new TSquareMatrix3D*[N_levels];    
  SqmatrixA31 = new TSquareMatrix3D*[N_levels];
  SqmatrixA32 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA33 = new TSquareMatrix3D*[N_levels];    
  
  MatrixB1 = new TMatrix3D*[N_levels];
  MatrixB2 = new TMatrix3D*[N_levels];
  MatrixB3 = new TMatrix3D*[N_levels];  
  
  MatrixB1T = new TMatrix3D*[N_levels];
  MatrixB2T = new TMatrix3D*[N_levels];
  MatrixB3T = new TMatrix3D*[N_levels];    
 
  AMatRhsAssemble = new TAssembleMat3D*[N_levels];  
  AMatAssembleNonLinear = new TAssembleMat3D*[N_levels];  
  
  if(TDatabase::ParamDB->INTERFACE_FLOW)
   {
    if(NSEType!=4)
    {
     OutPut(" FREE SURFACE or INTERFACE_FLOW needs NSEType 4"<<endl);
     exit(0);
    }    
    SqmatrixF11 = new TSquareMatrix3D*[N_levels]; 
    SqmatrixF22 = new TSquareMatrix3D*[N_levels]; 
    SqmatrixF33 = new TSquareMatrix3D*[N_levels]; 
    
    N_FreeSurfFaces = new int[N_levels];
    FreeSurfCellNumbers = new int*[N_levels];
    FreeSurfJointNumbers = new int*[N_levels];
   }
  
  
  if(SOLVER==AMG_SOLVE || SOLVER==DIRECT)
   {
    Start_Level=N_Levels-1;
   }
  else 
   {
    Start_Level=0;
   }
   
  // build matrices   
  for(i=Start_Level;i<N_Levels;i++)
   {
    if(SOLVER==GMG)
     OutPut("MULTIGRID LEVEL : " << i<<endl;)
    
     // first build matrix structure
     sqstructureA[i] = new TSquareStructure3D(U_Space[i]);
     sqstructureA[i]->Sort();  // sort column numbers: numbers are in increasing order
  
     structureB[i] = new TStructure3D(P_Space[i], U_Space[i]);
     structureBT[i] = new TStructure3D(U_Space[i], P_Space[i]);   

    switch(NSEType)
     {
      case 1:
        MatrixB1[i] = new TMatrix3D(structureB[i]);
        MatrixB2[i] = new TMatrix3D(structureB[i]);
        MatrixB3[i] = new TMatrix3D(structureB[i]);

        SqmatrixA11[i] = new TSquareMatrix3D(sqstructureA[i]);

        MatVect = MatVect_NSE1;
        Defect = Defect_NSE1;
      break;

      case 2:
        MatrixB1[i] = new TMatrix3D(structureB[i]);
        MatrixB2[i] = new TMatrix3D(structureB[i]);
        MatrixB3[i] = new TMatrix3D(structureB[i]);
        MatrixB1T[i] = new TMatrix3D(structureBT[i]);
        MatrixB2T[i] = new TMatrix3D(structureBT[i]);
        MatrixB3T[i] = new TMatrix3D(structureBT[i]);

        SqmatrixA11[i] = new TSquareMatrix3D(sqstructureA[i]);

        MatVect = MatVect_NSE2;
        Defect = Defect_NSE2;
      break;

      case 3:
        MatrixB1[i] = new TMatrix3D(structureB[i]);
        MatrixB2[i] = new TMatrix3D(structureB[i]);
        MatrixB3[i] = new TMatrix3D(structureB[i]);
       
        SqmatrixA11[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA12[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA13[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA21[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA22[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA23[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA31[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA32[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA33[i] = new TSquareMatrix3D(sqstructureA[i]);

        MatVect = MatVect_NSE3;
        Defect = Defect_NSE3;
      break;

      case 4:
        MatrixB1[i] = new TMatrix3D(structureB[i]);
        MatrixB2[i] = new TMatrix3D(structureB[i]);
        MatrixB3[i] = new TMatrix3D(structureB[i]);
        MatrixB1T[i] = new TMatrix3D(structureBT[i]);
        MatrixB2T[i] = new TMatrix3D(structureBT[i]);
        MatrixB3T[i] = new TMatrix3D(structureBT[i]);

        SqmatrixA11[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA12[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA13[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA21[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA22[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA23[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA31[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA32[i] = new TSquareMatrix3D(sqstructureA[i]);
        SqmatrixA33[i] = new TSquareMatrix3D(sqstructureA[i]);

        MatVect = MatVect_NSE4;
        Defect = Defect_NSE4;

        if(TDatabase::ParamDB->INTERFACE_FLOW)
         {
          SqmatrixF11[i] = new TSquareMatrix3D(sqstructureA[i]);
          SqmatrixF22[i] = new TSquareMatrix3D(sqstructureA[i]);
          SqmatrixF33[i] = new TSquareMatrix3D(sqstructureA[i]);
         }
      break;
      
      default:
            OutPut("Unknown NSETYPE " << NSEType <<"  it must be 1 to 4" << endl);
            exit(4711);;      
      
     }  
     OutPut(endl;)
    } //  for(i=Start_Level;i<N_Levels;i++)
    
#ifdef _MPI
  i = N_Levels-1; 
  switch(NSEType)
      {
        case 1:
          printf("parallel implemented for NSTYPE 2 only\n");
	  MPI_Finalize();
	  exit(0);
        break;

//         case 2:
// 	  SQMATRICES[0] = SqmatrixA11[i];
//   
//           MATRICES[0] = MatrixB1[i];
//           MATRICES[1] = MatrixB2[i];
//           MATRICES[2] = MatrixB3[i];
//           MATRICES[3] = MatrixB1T[i];
//           MATRICES[4] = MatrixB2T[i];
//           MATRICES[5] = MatrixB3T[i];
          
//         break;

        case 3:
          printf("parallel implemented for NSTYPE 2 only\n");
	  MPI_Finalize();
	  exit(0);
        break;

//         case 4:
//             SQMATRICES[0] = SqmatrixA11[i];
// 	    SQMATRICES[1] = SqmatrixA12[i];
//             SQMATRICES[2] = SqmatrixA13[i];
// 	    SQMATRICES[3] = SqmatrixA21[i];
// 	    SQMATRICES[4] = SqmatrixA22[i];
// 	    SQMATRICES[5] = SqmatrixA23[i];
// 	    SQMATRICES[6] = SqmatrixA31[i];
// 	    SQMATRICES[7] = SqmatrixA32[i];
// 	    SQMATRICES[8] = SqmatrixA33[i];
//           MATRICES[0] = MatrixB1[i];
//           MATRICES[1] = MatrixB2[i];
//           MATRICES[2] = MatrixB3[i];
//           MATRICES[3] = MatrixB1T[i];
//           MATRICES[4] = MatrixB2T[i];
//           MATRICES[5] = MatrixB3T[i];

//           break;
      } //  switch(NSEType)

  double t1,t2,tdiff;
  if(profiling)  t1 = MPI_Wtime();

  ParMapper_U = new TParFEMapper3D*[N_levels]; 
  ParMapper_P = new TParFEMapper3D*[N_levels];

  ParComm_U = new TParFECommunicator3D*[N_levels];
  ParComm_P = new TParFECommunicator3D*[N_levels];
  
  for(i=Start_Level;i<N_levels;i++)
  {
    ParMapper_U[i] = new TParFEMapper3D(3, U_Space[i], sqstructureA[i]->GetRowPtr(),sqstructureA[i]->GetKCol());
    ParMapper_P[i] = new TParFEMapper3D(1, P_Space[i], structureBT[i]->GetRowPtr(),  structureBT[i]->GetKCol());
      
    ParComm_U[i] = new TParFECommunicator3D(ParMapper_U[i]);
    ParComm_P[i] = new TParFECommunicator3D(ParMapper_P[i]);
  }



  if(profiling)
  {
     t2 = MPI_Wtime();
     tdiff = t2-t1;
     int out_rank=TDatabase::ParamDB->Par_P0;
     int rank;
     MPI_Comm_rank(Comm, &rank);
     MPI_Reduce(&tdiff, &t1, 1, MPI_DOUBLE, MPI_MIN, out_rank, Comm);
     if(rank == out_rank)
     {
      printf("Time taken for ParFEMapper3D for all levels is %e\n", t1);
     }
  }
#endif
   //initialize multigrid solver
   if(SOLVER==GMG)
   {    
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
    MG = new TNSE_MultiGrid(1, 2, Parameters);
    
     switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
      {
       case 11:
          zerostart = 0;
       break;
       case 16:
           zerostart = 1;
        break;
      default:
         zerostart = 1;
       }    
    
    // build preconditioner
    switch (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
     {
      case 5:
       prec = new TMultiGridIte(MatVect, Defect, NULL, 0, N_TotalDOF, MG, zerostart);
       Itmethod_sol = new double[N_TotalDOF];
       Itmethod_rhs = new double[N_TotalDOF];
      break;
      default:
         OutPut("Unknown preconditioner !!!" << endl);
         exit(4711);
     }
     
    // build solver
    switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
     {
        case 11:
            Itmethod = new TFixedPointIte(MatVect, Defect, prec, 0, N_TotalDOF, 0);
         break;
         case 16:
            Itmethod = new TFgmresIte(MatVect, Defect, prec, 0, N_TotalDOF, 0
            #ifdef _MPI 
,ParComm_U[N_levels-1], ParComm_P[N_levels-1]
#endif
);
         break;
         default:
            OutPut("Unknown solver !!!" << endl);
            exit(4711);
      }         
     
        
     } //  if(SOLVER==GMG)  
   

#ifdef _MPI
  int rank;
  TCollection *coll = U_Space[N_levels-1]->GetCollection();
  int owncells = coll->GetN_OwnCells();
  int problemSize=0;
  
  MPI_Comm_rank(Comm, &rank);
    
  MPI_Reduce(&owncells, &problemSize, 1, MPI_INT, MPI_SUM, 0, Comm);
  if(rank==0)
    OutPut( "total own cells over all sub domains : " << problemSize << endl);
  
  problemSize=0;
  int n_master = ParComm_U[N_levels-1]->GetN_Master();
  MPI_Reduce(&n_master, &problemSize, 1, MPI_INT, MPI_SUM, 0, Comm);
  if(rank==0)
    OutPut( "total own dofs over all sub domains : " << problemSize*3 << endl);
  
    coll = P_Space[N_levels-1]->GetCollection();
  owncells = coll->GetN_OwnCells();
  problemSize=0;
  MPI_Reduce(&owncells, &problemSize, 1, MPI_INT, MPI_SUM, 0, Comm);
  if(rank==0)
    OutPut( "pressure :: total own cells over all sub domains : " << problemSize << endl);
  
  problemSize=0;
  n_master = ParComm_P[N_levels-1]->GetN_Master();
  MPI_Reduce(&n_master, &problemSize, 1, MPI_INT, MPI_SUM, 0, Comm);
  if(rank==0)
    OutPut( "pressure :: total own dofs over all sub domains : " << problemSize << endl);
#endif
    

 
    fefct_aux = new  TFEFunction3D*[N_levels*6];  
    NSEaux = new TAuxParam3D*[N_levels];
    NSEaux_error = new TAuxParam3D*[N_levels];    

   for(i=Start_Level;i<N_Levels;i++)
   {
    NSEaux_error[i] = NULL;
    NSEaux[i] =  NULL;         
   }
   
}

// TSystemNSE3D::~TSystemNSE3D()
// {
//     delete NSEaux; 
//        
//     if(NSEaux_error)
//       delete NSEaux_error;
// 
// }


void TSystemNSE3D::Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
                           BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue)
{ 

  int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  int N_U_Current;
  int velocity_space_code, pressure_space_code;
  int mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
    
  double alpha[2];  
  
  TDiscreteForm3D *DiscreteFormGalerkin, *DiscreteFormSDFEM, *DiscreteFormUpwind, *DiscreteFormSmagorinsky;
  TDiscreteForm3D *DiscreteFormVMSProjection, *DiscreteFormNLGalerkin, *DiscreteFormNLSDFEM, *DiscreteFormNLUpwind;
  TDiscreteForm3D *DiscreteFormNLSmagorinsky, *DiscreteFormNLVMSProjection, *DiscreteFormPressSep, *DiscreteFormAuxProbPressSep;
  TDiscreteForm3D *DiscreteFormNSRFBRhs, *DiscreteFormNLSDFEM_DivDiv;
    
  /**  save the boundary condition */
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  
  BoundaryConditions[2] = BoundCond;  
  
  /**  save the boundary values   */
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
  BoundaryValues[2] = U3BoundValue;
  
  /** save the nse bilinear coefficient   */
  LinCoeffs[0] = lincoeffs;
  
  /** set the Discreteforms */
  InitializeDiscreteForms(DiscreteFormGalerkin, DiscreteFormSDFEM,
                          DiscreteFormUpwind, DiscreteFormSmagorinsky,
                          DiscreteFormVMSProjection,
                          DiscreteFormNLGalerkin, DiscreteFormNLSDFEM,
                          DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky, 
                          DiscreteFormNLVMSProjection,
                          DiscreteFormNLSDFEM_DivDiv,
                          DiscreteFormPressSep,
                          DiscreteFormAuxProbPressSep,
                          DiscreteFormNSRFBRhs,
                          LinCoeffs[0], NSEType);
  
  
    /** find discrete form */
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
          break;

          case SDFEM:
            DiscreteFormARhs = DiscreteFormSDFEM;
            DiscreteFormNL = DiscreteFormNLSDFEM; 
          break;

          case UPWIND:
            DiscreteFormARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;    
            break;

          case SMAGORINSKY:
            DiscreteFormARhs = DiscreteFormSmagorinsky;
            DiscreteFormNL = DiscreteFormNLSmagorinsky;              
            break;

          case VMS_PROJECTION:
	      DiscreteFormARhs = DiscreteFormVMSProjection;
	      if (NSEType != 1)
	       {
                OutPut("VMS only for NSTYPE 1 implemented !!!"<<endl);
		exit(4711);
	       }

            DiscreteFormNL = DiscreteFormNLVMSProjection;
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
     // initialize matrices
     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixA11[i];
          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];

          N_SquareMatrices = 1;
          N_RectMatrices = 3;
        break;

        case 2:
          SQMATRICES[0] = SqmatrixA11[i];
          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];
          MATRICES[3] = MatrixB1T[i];
          MATRICES[4] = MatrixB2T[i];
          MATRICES[5] = MatrixB3T[i];
  
          N_SquareMatrices = 1;
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

          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];

          N_SquareMatrices = 9;
          N_RectMatrices = 3;
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
          MATRICES[0] = MatrixB1[i];
          MATRICES[1] = MatrixB2[i];
          MATRICES[2] = MatrixB3[i];
          MATRICES[3] = MatrixB1T[i];
          MATRICES[4] = MatrixB2T[i];
          MATRICES[5] = MatrixB3T[i];

          N_SquareMatrices = 9;
          N_RectMatrices = 6;

          break;
      } //  switch(NSEType)
      
      N_Rhs = 3;
      N_FESpaces = 2;   
     
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();   
      
      fesp_aux[0] =  U_Space[i];
      fesp_aux[1] =  P_Space[i];
      
      fefct_aux[i*6 ] = Velocity[i]->GetComponent(0);
      fefct_aux[i*6+1] = Velocity[i]->GetComponent(1);
      fefct_aux[i*6+2] = Velocity[i]->GetComponent(2);

      
      NSEaux[i] =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo, NSN_FEValuesVelo,
                                fesp_aux, fefct_aux+(i*6), NSFctVelo, NSFEFctIndexVelo, NSFEMultiIndexVelo,
                                 NSN_ParamsVelo, NSBeginParamVelo);

      fesp[0] =  U_Space[i];
      fesp[1] =  P_Space[i];
      
      fefct[0] = Velocity[i]->GetComponent(0);
      fefct[1] = Velocity[i]->GetComponent(1);
      fefct[2] = Velocity[i]->GetComponent(2);      
     
      fesprhs[0] =  U_Space[i];
      fesprhs[1] =  U_Space[i];
      fesprhs[2] =  U_Space[i];

      RHSs[0] = RhsArray[i];
      RHSs[1] = RhsArray[i] + N_U_Current;
      RHSs[2] = RhsArray[i] + 2*N_U_Current;
      RHSs[3] = RhsArray[i] + 3*N_U_Current;
      
     // array of assemble objects
     AMatRhsAssemble[i] = new TAssembleMat3D(N_FESpaces, fesp, N_SquareMatrices, SQMATRICES, N_RectMatrices, MATRICES,
                              N_Rhs, RHSs, fesprhs, DiscreteFormARhs, BoundaryConditions, BoundaryValues, NSEaux[i]);
     AMatRhsAssemble[i]->Init();    
     
     if(TDatabase::ParamDB->INTERFACE_FLOW)
      {
       this->FindFreeSurfJoints(i, 0);
      }
        
#ifdef _MPI
   if(i == N_Levels-1) {
    if(SOLVER == DIRECT)
     {
      DS = new TParDirectSolver(ParComm_U[N_Levels-1],ParComm_P[N_Levels-1],SQMATRICES,MATRICES);
     }
   }
#endif

 //    ===============================================================================================================
 //    set the nonliner matrices
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
           }
          else
           {
            // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }

         break;
        } // switch(NSEType)
            
      N_RectMatrices = 0;          
      N_Rhs = 0;
      N_FESpaces = 1;
     
     AMatAssembleNonLinear[i] = new TAssembleMat3D(N_FESpaces, fesp, N_SquareMatrices, SQMATRICES, N_RectMatrices, NULL,
                              N_Rhs, NULL, NULL, DiscreteFormNL, BoundaryConditions, BoundaryValues, NSEaux[i]);
     AMatAssembleNonLinear[i]->Init();       
     
     //===============================================================================================================
//      // set the slip with friction assemble matrices
//        if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
//         {
//          if(NSEType <4)
//           {
//            OutPut("For slip with friction bc NSTYPE 4 is necessary !!!!! " << endl);
//            exit(4711);
//           }
          
//           // prepare everything for the assembling of slip with friction bc
//           // on all levels
//           N_FESpaces = 1;
//           N_SquareMatrices = 9;
//           N_RectMatrices = 0;
//           N_Rhs = 3; 
//       
//           SQMATRICES[0] = SqmatrixA11[i];
//           SQMATRICES[1] = SqmatrixA22[i];
//           SQMATRICES[2] = SqmatrixA33[i];
//           SQMATRICES[3] = SqmatrixA12[i];
//           SQMATRICES[4] = SqmatrixA13[i];
//           SQMATRICES[5] = SqmatrixA21[i];
//           SQMATRICES[6] = SqmatrixA23[i];
//           SQMATRICES[7] = SqmatrixA31[i];
//           SQMATRICES[8] = SqmatrixA32[i];
// 
//           Assemble3DSlipBC(N_FESpaces, fesp,
//             N_SquareMatrices, SQMATRICES,
//             N_RectMatrices, NULL,
//             N_Rhs, RHSs, fesprhs,
//             NULL,
//             BoundaryConditions,
//             BoundaryValues,
//             NSEaux);     
//      
     
// 	} // if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
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
              MGLevel = new TNSE_MGLevel1(i, SqmatrixA11[i],  MatrixB1[i], MatrixB2[i], MatrixB3[i],
                                             structureBT[i], RhsArray[i], SolArray[i], N_aux, alpha,
                                             velocity_space_code, pressure_space_code,
                                             NULL, NULL);
              MG->AddLevel(MGLevel);
          break;
          case 2:
              MGLevel = new TNSE_MGLevel2(i, SqmatrixA11[i], MatrixB1[i], MatrixB2[i], MatrixB3[i],
                                             MatrixB1T[i], MatrixB2T[i], MatrixB3T[i], RhsArray[i], SolArray[i], N_aux, alpha,
                                             velocity_space_code, pressure_space_code,
                                             NULL, NULL);
              MG->AddLevel(MGLevel);
          break;
          case 3:
              MGLevel = new TNSE_MGLevel3(i, SqmatrixA11[i], SqmatrixA12[i], SqmatrixA13[i],
                                             SqmatrixA21[i], SqmatrixA22[i], SqmatrixA23[i],
                                             SqmatrixA31[i], SqmatrixA32[i], SqmatrixA33[i],
                                             MatrixB1[i], MatrixB2[i], MatrixB3[i],
                                             structureBT[i], RhsArray[i], SolArray[i], N_aux, alpha,
                                             velocity_space_code, pressure_space_code,
                                             NULL, NULL);
              MG->AddLevel(MGLevel);
          break;
          case 4:
              MGLevel = new TNSE_MGLevel4(i, SqmatrixA11[i], SqmatrixA12[i], SqmatrixA13[i],
                                             SqmatrixA21[i], SqmatrixA22[i], SqmatrixA23[i],
                                             SqmatrixA31[i], SqmatrixA32[i], SqmatrixA33[i],
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
              
//  cout << " TSystemNSE3D::Init done ! " << endl; 
              
} // TSystemNSE3D::Init


void TSystemNSE3D::UpdateUpwind(int i)
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
            // do upwinding with three matrices
            UpwindForNavierStokes3D(SqmatrixA11[i], fefct[0], fefct[1], fefct[2]);
            UpwindForNavierStokes3D(SqmatrixA22[i], fefct[0], fefct[1], fefct[2]);
            UpwindForNavierStokes3D(SqmatrixA33[i], fefct[0], fefct[1], fefct[2]);
            cout << "UPWINDING DONE : level " << endl;
           break;
         }                        // endswitch      
      
      
}


void TSystemNSE3D::Assemble()
{
  int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  int N_U_Current, N_P_Current, N_Active_Current, N_DirichletDof;
    
  double alpha[2];

  
   for(i=Start_Level;i<N_Levels;i++)
    {     
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();
      N_Active_Current  = U_Space[i]->GetActiveBound();     
      N_DirichletDof = N_U_Current - N_Active_Current;
      N_P_Current = P_Space[i]->GetN_DegreesOfFreedom();      
      
      // initialize matrices
      AMatRhsAssemble[i]->Reset();
 
      /** assemble */
      AMatRhsAssemble[i]->Assemble3D();
 
      /** free surface/interface integration */
//       if(TDatabase::ParamDB->INTERFACE_FLOW)
//        {      
//         FreeSurfInt(U_Space[i]->GetCollection(), N_FreeSurfFaces[i], FreeSurfCellNumbers[i], FreeSurfJointNumbers[i],
    
	 
/*		      if(N_FreeJoints)
   {
    FreeSurfCellNumbers[level] = new int[N_FreeJoints];
    FreeSurfJointNumbers[level] = new int[N_FreeJoints];*/  

//   void FreeSurfInt(TCollection *Coll, int N_BoundFaces,
//                  int *CellNumbers, int *JointNumbers,
//                  TFEFunction3D *potential, double dt,
//                  TSquareMatrix3D *Aii,
//                  double *rhs1, double *rhs2, double *rhs3)
  
//        }

       /** upwind */
      if( (Disctype==UPWIND) && (TDatabase::ParamDB->PROBLEM_TYPE!=STOKES) )
        {
         this->UpdateUpwind(i);
        }    
       

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
    
      // set rhs for Dirichlet nodes
      memcpy(SolArray[i]+N_Active_Current, RhsArray[i]+N_Active_Current, N_DirichletDof*SizeOfDouble);
      memcpy(SolArray[i]+N_U_Current+N_Active_Current, RhsArray[i]+N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble); 
      memcpy(SolArray[i]+2*N_U_Current+N_Active_Current, RhsArray[i]+2*N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble);     
        
     } // for(i=Start_Level;i<N_Levels;i++)
      

//     cout << "Test Assemble " << endl; 
} // TSystemNSE3D::Assemble(T

void TSystemNSE3D::AssembleNonLinear(double **sol, double **rhs)
{
 int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
 int N_U_Current, N_Active_Current, N_DirichletDof;
  
   for(i=Start_Level;i<N_Levels;i++)
    {    
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();
      N_Active_Current  = U_Space[i]->GetActiveBound();     
      N_DirichletDof = N_U_Current - N_Active_Current;
      
      // reset the nonliner matrices
      AMatAssembleNonLinear[i]->Reset();
            
      // assemble the nonlinear matrix */      
      AMatAssembleNonLinear[i]->Assemble3D();      
      
       /** upwind */
      if( (Disctype==UPWIND) && (TDatabase::ParamDB->PROBLEM_TYPE!=STOKES) )
        {
         this->UpdateUpwind(i);
        }     
        
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
     
     /** no change in rhs, so no need to update */
      // set rhs for Dirichlet nodes 
      memcpy(SolArray[i]+N_Active_Current, RhsArray[i]+N_Active_Current, N_DirichletDof*SizeOfDouble);
      memcpy(SolArray[i]+N_U_Current+N_Active_Current, RhsArray[i]+N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble); 
      memcpy(SolArray[i]+2*N_U_Current+N_Active_Current, RhsArray[i]+2*N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble);     
        
      } //
      
} //TSystemNSE3D::AssembleNonLinear(


void TSystemNSE3D::GetResidual(double *sol, double *rhs, double *res, double &impuls_residual, double &residual)
{

     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixA11[N_Levels-1];
          MATRICES[0] = MatrixB1[N_Levels-1];
          MATRICES[1] = MatrixB2[N_Levels-1];
          MATRICES[2] = MatrixB3[N_Levels-1];
        break;

        case 2:
          SQMATRICES[0] = SqmatrixA11[N_Levels-1];
          MATRICES[0] = MatrixB1[N_Levels-1];
          MATRICES[1] = MatrixB2[N_Levels-1];
          MATRICES[2] = MatrixB3[N_Levels-1];
          MATRICES[3] = MatrixB1T[N_Levels-1];
          MATRICES[4] = MatrixB2T[N_Levels-1];
          MATRICES[5] = MatrixB3T[N_Levels-1];
	  
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11[N_Levels-1];
          SQMATRICES[1] = SqmatrixA12[N_Levels-1];
          SQMATRICES[2] = SqmatrixA13[N_Levels-1];	  
          SQMATRICES[3] = SqmatrixA21[N_Levels-1];
          SQMATRICES[4] = SqmatrixA22[N_Levels-1];
          SQMATRICES[5] = SqmatrixA23[N_Levels-1]; 
          SQMATRICES[6] = SqmatrixA31[N_Levels-1];
          SQMATRICES[7] = SqmatrixA32[N_Levels-1];
          SQMATRICES[8] = SqmatrixA33[N_Levels-1];  

          MATRICES[0] = MatrixB1[N_Levels-1];
          MATRICES[1] = MatrixB2[N_Levels-1];
          MATRICES[2] = MatrixB3[N_Levels-1];
        break;

        case 4:
          SQMATRICES[0] = SqmatrixA11[N_Levels-1];
          SQMATRICES[1] = SqmatrixA12[N_Levels-1];
          SQMATRICES[2] = SqmatrixA13[N_Levels-1];	  
          SQMATRICES[3] = SqmatrixA21[N_Levels-1];
          SQMATRICES[4] = SqmatrixA22[N_Levels-1];
          SQMATRICES[5] = SqmatrixA23[N_Levels-1]; 
          SQMATRICES[6] = SqmatrixA31[N_Levels-1];
          SQMATRICES[7] = SqmatrixA32[N_Levels-1];
          SQMATRICES[8] = SqmatrixA33[N_Levels-1];  
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
       

    
    Defect(sqmatrices, matrices, sol, rhs, res); 
  

     
#ifdef _MPI
   double residual_scalar = 0.0;
   double sum =0.;
   int i,j,rank;
   MPI_Comm_rank(Comm, &rank);
   int *master = ParComm_U[N_Levels-1]->GetMaster();
   
   for(i=0;i<N_U;i++)
   {
     if(master[i]!=rank)    continue;
      
      residual_scalar += res[i      ]*res[i      ];
      residual_scalar += res[i+  N_U]*res[i+  N_U];
      residual_scalar += res[i+2*N_U]*res[i+2*N_U];

    }

    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
    impuls_residual = (sum);

    master = ParComm_P[N_Levels-1]->GetMaster();
    for(i=0;i<N_P;i++)
    {
      if(master[i]!=rank)    continue;
      
      residual_scalar += res[i+3*N_U]*res[i+3*N_U];

    }
    
    sum = 0;
    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
    residual = (sum);

#else
    impuls_residual  =  Ddot(3*N_U, res, res);
    residual         =  Ddot(N_TotalDOF, res, res); 

#endif
   
} // TSystemNSE3D::GetResidual

void TSystemNSE3D::Solve(double *sol, double *rhs)
{
  int N_LinIter=0;
  double summ = 0;
  double residual,residual_scalar = 0.0;
  double sum =0.;
  int i,j,rank;
  int *master;
   
    switch(SOLVER)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:

          if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
           {
            memcpy(Itmethod_sol, sol, N_TotalDOF*SizeOfDouble);
            memcpy(Itmethod_rhs, rhs, N_TotalDOF*SizeOfDouble);
           }
          // solve the linear system
          N_LinIter += Itmethod->Iterate(sqmatrices, matrices, Itmethod_sol, Itmethod_rhs);
	  
          if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
           {
            memcpy(sol, Itmethod_sol, N_TotalDOF*SizeOfDouble);
            memcpy(rhs, Itmethod_rhs, N_TotalDOF*SizeOfDouble);
           }
          MG->RestrictToAllGrids();

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
	    DS->Solve(sol, rhs, true);
	    
#else	    

	    
	    
	    DirectSolver(SqmatrixA11[N_Levels-1], MatrixB1T[N_Levels-1], MatrixB2T[N_Levels-1], MatrixB3T[N_Levels-1],
                          MatrixB1[N_Levels-1], MatrixB2[N_Levels-1], MatrixB3[N_Levels-1], rhs, sol);

#endif

	    
          break;

          case 3:
           cout << "Solver not included for NSTYPE 3 in this version" <<endl;
            cout << "try NSTYPE 2 or 4 " <<endl;   
            exit(0);
          break;

          case 4:
#ifdef _MPI
        	DS->Solve(sol, rhs, true);
#endif
	
#ifdef _OMPONLY
	       if(TDatabase::ParamDB->DSType == 1)
	         DS->Solve(sol, rhs, true);
	       else{
	         OutPut("Select Proper Solver" << endl);
	         exit(0);
	       }
#endif

#ifdef _SEQ

         
     
             DirectSolver(SqmatrixA11[N_Levels-1], SqmatrixA12[N_Levels-1], SqmatrixA13[N_Levels-1], 
                          SqmatrixA21[N_Levels-1], SqmatrixA22[N_Levels-1], SqmatrixA23[N_Levels-1],  
                          SqmatrixA31[N_Levels-1], SqmatrixA32[N_Levels-1], SqmatrixA33[N_Levels-1],  
                          MatrixB1T[N_Levels-1], MatrixB2T[N_Levels-1], MatrixB3T[N_Levels-1],
                          MatrixB1[N_Levels-1], MatrixB2[N_Levels-1], MatrixB3[N_Levels-1], rhs, sol,0);

#endif
          break;
         } //  switch(NSEType) 

      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
}

void TSystemNSE3D::MeasureErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2, DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP,
                                  double *u_error, double *p_error)
{
  double errors[4];
  fesp[0] =  U_Space[N_Levels-1];

    fefct[0] = Velocity[N_Levels-1]->GetComponent(0);
    fefct[1] = Velocity[N_Levels-1]->GetComponent(1);
    fefct[2] = Velocity[N_Levels-1]->GetComponent(2);     

  
     if(NSEaux_error[N_Levels-1]==NULL)
      NSEaux_error[N_Levels-1] =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo, NSN_FEValuesVelo,
                                    fesp, fefct, NSFctVelo, NSFEFctIndexVelo, NSFEMultiIndexVelo,
                                    NSN_ParamsVelo, NSBeginParamVelo);

     // errors in first velocity component
      fefct[0]->GetErrors(ExactU1, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1, U_Space+(N_Levels-1), errors);
   
      u_error[0] = errors[0];
      u_error[1] = errors[1];
      
     // errors in second velocity component
      fefct[1]->GetErrors(ExactU2, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1, U_Space+(N_Levels-1), errors);
      u_error[2] = errors[0];
      u_error[3] = errors[1];     
     
      // errors in third velocity component
      fefct[2]->GetErrors(ExactU3, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1, U_Space+(N_Levels-1), errors);
      u_error[4] = errors[0];
      u_error[5] = errors[1];       
     
     
      // errors in pressure
      Pressure[N_Levels-1]->GetErrors(ExactP, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1,  P_Space+(N_Levels-1), errors);     
      p_error[0] = errors[0];
      p_error[1] = errors[1]; 
}
    

void TSystemNSE3D::FindFreeSurfJoints(int level, int Region_ID)
{
  int i,j, N_Joints;
  int N_Cells, N_FreeJoints;
  
  TBaseCell *cell, *neigh, *neigh0, *neigh1;
  TJoint *joint;
  TCollection *Coll;
  TIsoJointEqN *isojoint;

  Coll = U_Space[level]->GetCollection();  
  N_Cells = Coll->GetN_Cells();
  N_FreeJoints=0;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->GetRegionID() == Region_ID)
    {
      N_Joints = cell->GetN_Joints();
      for(j=0;j<N_Joints;j++)
       {
        joint = cell->GetJoint(j);
        if(joint->GetType()==InterfaceJoint3D || joint->GetType()==IsoInterfaceJoint3D || joint->GetType()==IsoBoundFace)
         {  
          N_FreeJoints++;
         }
       } // endfor j
    } // end ID == 1
  } // endfor i

  N_FreeSurfFaces[level] = N_FreeJoints;
  
  if(N_FreeJoints)
   {
    FreeSurfCellNumbers[level] = new int[N_FreeJoints];
    FreeSurfJointNumbers[level] = new int[N_FreeJoints];    
   }
   
  N_FreeJoints = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->GetRegionID() == Region_ID)
    {
      N_Joints = cell->GetN_Joints();
      for(j=0;j<N_Joints;j++)
       {
        joint = cell->GetJoint(j);
        if(joint->GetType()==InterfaceJoint3D || joint->GetType()==IsoInterfaceJoint3D || joint->GetType()==IsoBoundFace)
         { 
          FreeSurfCellNumbers[level][N_FreeJoints] = i;
          FreeSurfJointNumbers[level][N_FreeJoints] = j;
          N_FreeJoints++;
         }
       } // endfor j
    } // end ID == 1
  } // endfor i
  
  
//   cout << " FindFreeSurfJoints done " << N_FreeJoints << endl;
//   exit(0);  
}
