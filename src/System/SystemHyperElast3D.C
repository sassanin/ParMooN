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
#    along with ParMooN. If not, see <http://www.gnu.org/licenses/>.
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
#include <SystemHyperElast3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <FEVectFunct3D.h>
#include <AuxParam3D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
// #include <NSE3D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Upwind3D.h>
#include <AssembleMat3D.h>
#include <LinAlg.h>

// #include <NSE_MultiGrid.h>
// #include <NSE_MGLevel1.h>
// #include <NSE_MGLevel2.h>
// #include <NSE_MGLevel3.h>
// #include <NSE_MGLevel4.h>

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

void Galerkin3D(double Mult, double *coeff, double *param, double hK, double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
void HyperParamsVelo(double *in, double *out);
double Piola_Kir(double *param, double test100, double test010, double test001, double test000, double ansatz100, double ansatz010, double ansatz001, double ansatz000, int row, int col);

/**  parameters for Auxpaam */
 int Hyper_FEFctIndex[9] = { 0, 1, 2, 0, 1, 2, 0, 1, 2};
 int Hyper_BeginParam[1] = { 0 };  
 MultiIndex3D Hyper_FEMultiIndex[9] = { D100, D100, D100, D010, D010, D010, D001, D001, D001 };  
 ParamFct *Hyper_ParamFct[1] = { HyperParamsVelo };


/** parameters for Discrete form */
  MultiIndex3D Derivatives[4] = { D100, D010, D001, D000};
  int SpaceNumbers[4] = { 0, 0, 0, 0};
  int RowSpace[9]    = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int ColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int RhsSpace[3] = { 0, 0, 0 };

TSystemHyperElast3D::TSystemHyperElast3D(int N_levels, TFESpace3D **disp_fespace, TFEVectFunct3D **displacement, double **sol, double **rhs, int disctype, int solver)
{
  int i, zerostart;
  int profiling = TDatabase::ParamDB->timeprofiling;
 
  /** need it for solver */
  sqmatrices = (TSquareMatrix **)SQMATRICES;
      
  //set number of multigrid levels
  N_Levels = N_levels;
   
//   cout << " N_levels " << N_levels <<endl;
//   exit(0);
//   
  //set the discretization type
  Disctype = disctype;
  
#ifdef _MPI
  Comm = TDatabase::ParamDB->Comm;
#endif  
  
  //set the solver type
  SOLVER = solver;
  
  Displacement = displacement;
 
  SolArray = sol;
  RhsArray = rhs;
  
  U_Space = disp_fespace;
  
  N_U = disp_fespace[N_levels-1]->GetN_DegreesOfFreedom();
  N_TotalDOF = 3*N_U;
       
  N_Active =  disp_fespace[N_levels-1]->GetActiveBound();
  N_DirichletDof = N_U - N_Active;  
 
  sqstructureA = new TSquareStructure3D *[N_levels];
   
  SqmatrixA11 = new TSquareMatrix3D*[N_levels];
  SqmatrixA12 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA13 = new TSquareMatrix3D*[N_levels];  
  SqmatrixA21 = new TSquareMatrix3D*[N_levels];
  SqmatrixA22 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA23 = new TSquareMatrix3D*[N_levels];    
  SqmatrixA31 = new TSquareMatrix3D*[N_levels];
  SqmatrixA32 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA33 = new TSquareMatrix3D*[N_levels];    
  
  AMatRhsAssemble = new TAssembleMat3D*[N_levels];  
  AMatAssembleNonLinear = new TAssembleMat3D*[N_levels];  
  
  if(TDatabase::ParamDB->INTERFACE_FLOW)
   {  
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

     SqmatrixA11[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA12[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA13[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA21[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA22[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA23[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA31[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA32[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA33[i] = new TSquareMatrix3D(sqstructureA[i]);

//         MatVect = MatVect_NSE4;
//         Defect = Defect_NSE4;

        if(TDatabase::ParamDB->INTERFACE_FLOW)
         {
          SqmatrixF11[i] = new TSquareMatrix3D(sqstructureA[i]);
          SqmatrixF22[i] = new TSquareMatrix3D(sqstructureA[i]);
          SqmatrixF33[i] = new TSquareMatrix3D(sqstructureA[i]);
         }
    } //  for(i=Start_Level;i<N_Levels;i++)
    
#ifdef _MPI
  i = N_Levels-1; 
 

  double t1,t2,tdiff;
  if(profiling)  t1 = MPI_Wtime();

  ParMapper_U = new TParFEMapper3D*[N_levels]; 
  ParComm_U = new TParFECommunicator3D*[N_levels];
  
  for(i=Start_Level;i<N_levels;i++)
  {
    ParMapper_U[i] = new TParFEMapper3D(3, U_Space[i], sqstructureA[i]->GetRowPtr(),sqstructureA[i]->GetKCol());

    ParComm_U[i] = new TParFECommunicator3D(ParMapper_U[i]);
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
//    if(SOLVER==GMG)
//    {    
//     Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
//     Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
//     MG = new TNSE_MultiGrid(1, 2, Parameters);
//     
//      switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
//       {
//        case 11:
//           zerostart = 0;
//        break;
//        case 16:
//            zerostart = 1;
//         break;
//       default:
//          zerostart = 1;
//        }    
//     
//     // build preconditioner
//     switch (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
//      {
//       case 5:
//        prec = new TMultiGridIte(MatVect, Defect, NULL, 0, N_TotalDOF, MG, zerostart);
//        Itmethod_sol = new double[N_TotalDOF];
//        Itmethod_rhs = new double[N_TotalDOF];
//       break;
//       default:
//          OutPut("Unknown preconditioner !!!" << endl);
//          exit(4711);
//      }
//      
//     // build solver
//     switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
//      {
//         case 11:
//             Itmethod = new TFixedPointIte(MatVect, Defect, prec, 0, N_TotalDOF, 0);
//          break;
//          case 16:
//             Itmethod = new TFgmresIte(MatVect, Defect, prec, 0, N_TotalDOF, 0
//             #ifdef _MPI 
// ,ParComm_U[N_levels-1], ParComm_P[N_levels-1]
// #endif
// );
//          break;
//          default:
//             OutPut("Unknown solver !!!" << endl);
//             exit(4711);
//       }         
//      
//         
//      } //  if(SOLVER==GMG)  
   

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
    
  
    fefct_aux = new  TFEFunction3D*[N_levels*3]; 
    Hyperaux = new TAuxParam3D*[N_levels];
    Hyperaux_error = new TAuxParam3D*[N_levels];    

   for(i=Start_Level;i<N_Levels;i++)
   {
    Hyperaux_error[i] = NULL;
    Hyperaux[i] =  NULL;         
   }
   
} // TSystemHyperElast3D

// TSystemHyperElast3D::~TSystemHyperElast3D()
// {
//     delete [] Hyperaux; 
//        
//     if(Hyperaux_error)
//       delete [] Hyperaux_error;
// 
// }


void TSystemHyperElast3D::Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
                           BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue)
{ 
  int i, N_SquareMatrices, N_Rhs, N_FESpaces;
  int N_U_Current;
  int N_Terms = 4;
  int N_Matrices = 9;
   
  char GalerkinString[] = "Galerkin";
  char rhs[] = "rhs";
  char all[] = "all";
  
//   int mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
//     
//   double alpha[2];  
  
  TDiscreteForm3D *DiscreteFormGalerkin;
 
    
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
  
  N_Rhs = 3;
  DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,  N_Terms, Derivatives, SpaceNumbers,
                                             N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace, Galerkin3D, LinCoeffs[0], NULL); 
  
    /** find discrete form */
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
          break;
 
          default:
            Error("Unknown DISCTYPE" << endl);
            exit(-1);
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
     SQMATRICES[0] = SqmatrixA11[i];
     SQMATRICES[1] = SqmatrixA12[i];
     SQMATRICES[2] = SqmatrixA13[i];	  
     SQMATRICES[3] = SqmatrixA21[i];
     SQMATRICES[4] = SqmatrixA22[i];
     SQMATRICES[5] = SqmatrixA23[i]; 
     SQMATRICES[6] = SqmatrixA31[i];
     SQMATRICES[7] = SqmatrixA32[i];
     SQMATRICES[8] = SqmatrixA33[i];  
 
     N_SquareMatrices = 9;

      N_Rhs = 3;
      N_FESpaces = 1;   
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();   
      
      this->InitHyperAuxParm(i);
           
      fesp[0] =  U_Space[i];
 
      fefct[0] = Displacement[i]->GetComponent(0);
      fefct[1] = Displacement[i]->GetComponent(1);
      fefct[2] = Displacement[i]->GetComponent(2);      
     
      fesprhs[0] =  U_Space[i];
      fesprhs[1] =  U_Space[i];
      fesprhs[2] =  U_Space[i];

      RHSs[0] = RhsArray[i];
      RHSs[1] = RhsArray[i] + N_U_Current;
      RHSs[2] = RhsArray[i] + 2*N_U_Current;
//    cout << " Hyper_FEMultiIndex   3D: " << Hyper_FEMultiIndex[0] << endl;
  
     // array of assemble objects
     AMatRhsAssemble[i] = new TAssembleMat3D(N_FESpaces, fesp, N_SquareMatrices, SQMATRICES, 0, NULL,
                                             N_Rhs, RHSs, fesprhs, DiscreteFormARhs, BoundaryConditions, BoundaryValues, Hyperaux[i]);
     AMatRhsAssemble[i]->Init();    
        
//      if(TDatabase::ParamDB->INTERFACE_FLOW)
//       {
//        this->FindFreeSurfJoints(i, 0);
//       }
        
#ifdef _MPI
   if(i == N_Levels-1) {
    if(SOLVER == DIRECT)
     {
      DS = new TParDirectSolver(ParComm_U[N_Levels-1], ParComm_P[N_Levels-1], SQMATRICES, NULL);
     }
   }
#endif
    
     // initialize solver
       if(SOLVER==GMG)
        {       
//          //setup the multigrid solver
//          alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
//          alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;  
//          velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
//          pressure_space_code = TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE;
// 
//          if(mg_type==1)
//           {
//            if(i==0)
//             {
//              alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
//              alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
//             }     
//            else if(i==N_Levels-1)
//             {
//              alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
//              alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE;  
//             }
//           
//            if(i<N_Levels-1)
//             {
//              velocity_space_code = -1;
//              pressure_space_code = 0; 
//             }          
//           }  
//   
//          switch(NSEType)
//           {
//            case 1:
//               MGLevel = new TNSE_MGLevel1(i, SqmatrixA11[i],  MatrixB1[i], MatrixB2[i], MatrixB3[i],
//                                              structureBT[i], RhsArray[i], SolArray[i], N_aux, alpha,
//                                              velocity_space_code, pressure_space_code,
//                                              NULL, NULL);
//               MG->AddLevel(MGLevel);
//           break;
//           case 2:
//               MGLevel = new TNSE_MGLevel2(i, SqmatrixA11[i], MatrixB1[i], MatrixB2[i], MatrixB3[i],
//                                              MatrixB1T[i], MatrixB2T[i], MatrixB3T[i], RhsArray[i], SolArray[i], N_aux, alpha,
//                                              velocity_space_code, pressure_space_code,
//                                              NULL, NULL);
//               MG->AddLevel(MGLevel);
//           break;
//           case 3:
//               MGLevel = new TNSE_MGLevel3(i, SqmatrixA11[i], SqmatrixA12[i], SqmatrixA13[i],
//                                              SqmatrixA21[i], SqmatrixA22[i], SqmatrixA23[i],
//                                              SqmatrixA31[i], SqmatrixA32[i], SqmatrixA33[i],
//                                              MatrixB1[i], MatrixB2[i], MatrixB3[i],
//                                              structureBT[i], RhsArray[i], SolArray[i], N_aux, alpha,
//                                              velocity_space_code, pressure_space_code,
//                                              NULL, NULL);
//               MG->AddLevel(MGLevel);
//           break;
//           case 4:
//               MGLevel = new TNSE_MGLevel4(i, SqmatrixA11[i], SqmatrixA12[i], SqmatrixA13[i],
//                                              SqmatrixA21[i], SqmatrixA22[i], SqmatrixA23[i],
//                                              SqmatrixA31[i], SqmatrixA32[i], SqmatrixA33[i],
//                                              MatrixB1[i], MatrixB2[i], MatrixB3[i],
//                                              MatrixB1T[i], MatrixB2T[i], MatrixB3T[i], RhsArray[i], SolArray[i], N_aux, alpha,
//                                              velocity_space_code, pressure_space_code,
//                                              NULL, NULL
// #ifdef _MPI
// 	  , ParComm_U[i], ParComm_P[i]
// #endif
// 	  );
//               MG->AddLevel(MGLevel);
//           break;	
//         } //  switch(NSEType)
       }  // if(SOLVER==GMG)     
     } // for(i=Start_Level;i<N_Levels;i++)      
//               
//  cout << " TSystemHyperElast3D::Init done ! " << endl; 
              
} // TSystemHyperElast3D::Init


// void TSystemNSE3D::UpdateUpwind(int i)
// {
//        
//     fefct[0] = Velocity[i]->GetComponent(0);
//     fefct[1] = Velocity[i]->GetComponent(1);
//     fefct[2] = Velocity[i]->GetComponent(2);
//     
//          switch(NSEType)
//           {
//            case 1:
//            case 2:
//             // do upwinding with one matrix
//             UpwindForNavierStokes3D(SqmatrixA11[i], fefct[0], fefct[1], fefct[2]);
//             cout << "UPWINDING DONE : level " << endl;
//             break;
// 
//           case 3:
//           case 4:
//             // do upwinding with three matrices
//             UpwindForNavierStokes3D(SqmatrixA11[i], fefct[0], fefct[1], fefct[2]);
//             UpwindForNavierStokes3D(SqmatrixA22[i], fefct[0], fefct[1], fefct[2]);
//             UpwindForNavierStokes3D(SqmatrixA33[i], fefct[0], fefct[1], fefct[2]);
//             cout << "UPWINDING DONE : level " << endl;
//            break;
//          }                        // endswitch      
//       
//       
// }
// 
// 
void TSystemHyperElast3D::Assemble()
{
  int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  int N_U_Current,  N_Active_Current, N_DirichletDof;
    
  double alpha[2];

//   
   for(i=Start_Level;i<N_Levels;i++)
    {     
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();
      N_Active_Current  = U_Space[i]->GetActiveBound();     
      N_DirichletDof = N_U_Current - N_Active_Current;
      
      // initialize matrices
      AMatRhsAssemble[i]->Reset();

      /** assemble */
      AMatRhsAssemble[i]->Assemble3D();
 
//       /** free surface/interface integration */
// //       if(TDatabase::ParamDB->INTERFACE_FLOW)
// //        {      
// //         FreeSurfInt(U_Space[i]->GetCollection(), N_FreeSurfFaces[i], FreeSurfCellNumbers[i], FreeSurfJointNumbers[i],
//     
// 	 
// /*		      if(N_FreeJoints)
//    {
//     FreeSurfCellNumbers[level] = new int[N_FreeJoints];
//     FreeSurfJointNumbers[level] = new int[N_FreeJoints];*/  
// 
// //   void FreeSurfInt(TCollection *Coll, int N_BoundFaces,
// //                  int *CellNumbers, int *JointNumbers,
// //                  TFEFunction3D *potential, double dt,
// //                  TSquareMatrix3D *Aii,
// //                  double *rhs1, double *rhs2, double *rhs3)
//   
// //        }
// 
//        /** upwind */
//       if( (Disctype==UPWIND) && (TDatabase::ParamDB->PROBLEM_TYPE!=STOKES) )
//         {
//          this->UpdateUpwind(i);
//         }    
//        
// 
//        if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
//         {
//          if(NSEType <4)
//           {
//            OutPut("For slip with friction bc NSTYPE 4 is ");
//            OutPut("necessary !!!!! " << endl);
//            exit(4711);
//           }
//           
// //          AMatRhsAssemble[i]->AssembleNavierSlip(); 
//           
//           // prepare everything for the assembling of slip with friction bc
//           // on all levels
//           N_FESpaces = 1;
//           N_SquareMatrices = 9;
//           N_RectMatrices = 0;
//           N_Rhs = 3; 
// 
//           fesp[0] =  U_Space[i];
// 	  
//           fesprhs[0] =  U_Space[i];
//           fesprhs[1] =  U_Space[i];
//           fesprhs[2] =  U_Space[i];
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
//           RHSs[0] = RhsArray[i];
//           RHSs[1] = RhsArray[i] + N_U_Current;
//           RHSs[2] = RhsArray[i] + 2*N_U_Current;
//       
//           Assemble3DSlipBC(N_FESpaces, fesp,
//                            N_SquareMatrices, SQMATRICES,
//                            N_RectMatrices, NULL,
//                            N_Rhs, RHSs, fesprhs,
//                            NULL,
//                            BoundaryConditions,
//                            BoundaryValues,
//                            NSEaux[i]);
//         } //  if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRIC  
//     
//       // set rhs for Dirichlet nodes
//       memcpy(SolArray[i]+N_Active_Current, RhsArray[i]+N_Active_Current, N_DirichletDof*SizeOfDouble);
//       memcpy(SolArray[i]+N_U_Current+N_Active_Current, RhsArray[i]+N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble); 
//       memcpy(SolArray[i]+2*N_U_Current+N_Active_Current, RhsArray[i]+2*N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble);     
        
     } // for(i=Start_Level;i<N_Levels;i++)
//     cout << "Test Assemble " << endl; 
//     exit(0);
    
} // TSystemHyperElast3D::Assemble(T

void TSystemHyperElast3D::InitHyperAuxParm(int i)
{
  int Hyper_N_FESpace = 1;
  int Hyper_N_FEFunction = 3;
  int Hyper_N_ParamFct = 1;
  int Hyper_N_FEValues = 9;
  int Hyper_N_ParamValues = 99;                                      
                                            
  fesp_aux[0] =  U_Space[i];
     
  fefct_aux[i*3 ] = Displacement[i]->GetComponent(0);
  fefct_aux[(i*3)+1] = Displacement[i]->GetComponent(1);
  fefct_aux[(i*3)+2] = Displacement[i]->GetComponent(2);

  Hyperaux[i] =  new TAuxParam3D(Hyper_N_FESpace, Hyper_N_FEFunction, Hyper_N_ParamFct, Hyper_N_FEValues,
                                fesp_aux, fefct_aux+(i*3), Hyper_ParamFct, Hyper_FEFctIndex, Hyper_FEMultiIndex,
                                Hyper_N_ParamValues, Hyper_BeginParam);
}



void Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21, **MatrixA22, **MatrixA23; 
  double **MatrixA31, **MatrixA32, **MatrixA33;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRowA11, *MatrixRowA12, *MatrixRowA13, *MatrixRowA21, *MatrixRowA22, *MatrixRowA23;
  double *MatrixRowA31, *MatrixRowA32, *MatrixRowA33;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];    
  MatrixA31 = LocMatrices[6];  
  MatrixA32 = LocMatrices[7];  
  MatrixA33 = LocMatrices[8];


  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // phi_x
  Orig1 = OrigValues[1]; // phi_y
  Orig2 = OrigValues[2]; // phi_z
  Orig3 = OrigValues[3]; // phi


  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

/*  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  u1x = param[3]; // u1old
  u2x = param[4]; // u2old
  u3x = param[5]; // u3old 
  u1y = param[6]; // u1old
  u2y = param[7]; // u2old
  u3y = param[8]; // u3old  
  u1z = param[9]; // u1old
  u2z = param[10]; // u2old
  u3z = param[11]; // u3old */ 
 
  double Sder[9];
  double *g = new double[9];
  double h[9];
  double F[9];
  double v1[3]; 
  for(i=0;i<N_U;i++)
  {
    MatrixRowA11 = MatrixA11[i];
    MatrixRowA12 = MatrixA12[i];
    MatrixRowA13 = MatrixA13[i];
    MatrixRowA21 = MatrixA21[i];    
    MatrixRowA22 = MatrixA22[i];    
    MatrixRowA23 = MatrixA23[i];    
    MatrixRowA31 = MatrixA31[i];    
    MatrixRowA32 = MatrixA32[i];
    MatrixRowA33 = MatrixA33[i];    
    
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    
    val1 = Mult*test000;
    
    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;
    
    
    for(j=0;j<N_U;j++)
    { 
      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      for(int i =0; i<9; i++)
         F[i] = param[90+i]; 
      int d;
      //===============Blocks A11, A21, and A31 ================//
      //Construction of derivative of (F^t)F matrix ; H;
      v1[0] = param[90]; v1[1] = param[93]; v1[2] = param[96];
      
      for(int k =0; k < 3; k++){
       int d = 3*k;
       h[d] = ansatz100*v1[k];
       h[d+1] = ansatz010*v1[k];
       h[d+2] = ansatz001*v1[k];
       }
      for( int k =0; k < 3; k++){
       h[k] += ansatz100*v1[k];
       h[k+3] += ansatz010*v1[k];
       h[k+6] += ansatz001*v1[k];
       }
      //=====================================================//
      for(int l =0; l < 9; l++){
        Sder[l] = 0;
        d = 9*(l+1);
        for(int k =0 ; k< 9; k++){ 
        Sder[l] += h[k]*param[d+k];
         }
       }
      MatrixMult(F, Sder,g, 'n', 'n') ;
      g[0] += (ansatz100*param[0] + ansatz010*param[1] + ansatz001*param[2]);
      g[3] += (ansatz100*param[3] + ansatz010*param[4] + ansatz001*param[5]);
      g[6] += (ansatz100*param[6] + ansatz010*param[7] + ansatz001*param[8]);
      MatrixRowA11[j] += Mult * (test100*g[0] + test010*g[3] + test001*g[6]);
      MatrixRowA21[j] += Mult * (test100*g[1] + test010*g[4] + test001*g[7]);
      MatrixRowA31[j] += Mult * (test100*g[2] + test010*g[5] + test001*g[8]);

      //===============Blocks A12, A22, and A32 ================//      
      //Construction of derivative of (F^t)F matrix ; H;
      v1[0] = param[91]; v1[1] = param[94]; v1[2] = param[97];
      for(int k =0; k < 3; k++){
       int d = 3*k;
       h[d] = ansatz100*v1[k];
       h[d+1] = ansatz010*v1[k];
       h[d+2] = ansatz001*v1[k];
       }
      for( int k =0; k < 3; k++){
       h[k] += ansatz100*v1[k];
       h[k+3] += ansatz010*v1[k];
       h[k+6] += ansatz001*v1[k];
       }
      //=====================================================//
      for(int l =0; l < 9; l++){
        Sder[l] = 0;
        d = 9*(l+1);
        for(int k =0 ; k< 9; k++){ 
        Sder[l] += h[k]*param[d+k];
         }
       }
       MatrixMult(F, Sder,g, 'n', 'n') ;
       g[1] += (ansatz100*param[0] + ansatz010*param[1] + ansatz001*param[2]);
       g[4] += (ansatz100*param[3] + ansatz010*param[4] + ansatz001*param[5]);
       g[7] += (ansatz100*param[6] + ansatz010*param[7] + ansatz001*param[8]);
       MatrixRowA12[j] += Mult * (test100*g[0] + test010*g[3] + test001*g[6]);
       MatrixRowA22[j] += Mult * (test100*g[1] + test010*g[4] + test001*g[7]);
       MatrixRowA32[j] += Mult * (test100*g[2] + test010*g[5] + test001*g[8]);

      //===============Blocks A13, A23, and A33 ================//

      v1[0] = param[92]; v1[1] = param[95]; v1[2] = param[98];
      for(int k =0; k < 3; k++){
       int d = 3*k;
       h[d] = ansatz100*v1[k];
       h[d+1] = ansatz010*v1[k];
       h[d+2] = ansatz001*v1[k];
       }
      for( int k =0; k < 3; k++){
       h[k] += ansatz100*v1[k];
       h[k+3] += ansatz010*v1[k];
       h[k+6] += ansatz001*v1[k];
       }
      //=====================================================//
      for(int l =0; l < 9; l++){
        Sder[l] = 0;
        d = 9*(l+1);
        for(int k =0 ; k< 9; k++){ 
        Sder[l] += h[k]*param[d+k];
         }
       }
       MatrixMult(F, Sder,g, 'n', 'n') ;
       g[2] += (ansatz100*param[0] + ansatz010*param[1] + ansatz001*param[2]);
       g[5] += (ansatz100*param[3] + ansatz010*param[4] + ansatz001*param[5]);
       g[8] += (ansatz100*param[6] + ansatz010*param[7] + ansatz001*param[8]);
       MatrixRowA13[j] += Mult * (test100*g[0] + test010*g[3] + test001*g[6]);
       MatrixRowA23[j] += Mult * (test100*g[1] + test010*g[4] + test001*g[7]);
       MatrixRowA33[j] += Mult * (test100*g[2] + test010*g[5] + test001*g[8]);
     
      
    } // endfor j
  } // endfor i

//     for (i = 0; i < 99; i++)
//   cout << "param " << i << ": " << param[i] <<endl;
   
// cout << " Galerkin3D " <<endl;
// exit(0);
 
 
}


void HyperParamsVelo(double *in, double *out)
{  
  /** based on the elastic model, S_ij, 0< i,j < 9 vales will be returned */

  int i, j, k, l, J, K, L, disp;
  const double c1=1.0;
  const double c2=1.0;
  const double D1=1.0;
 
  double I1, I3, I3pow, k1, k2, I1_bar, *OUT;
  double  F[9], C[9];
  double cij, cik, cjk, ckl, cil, cjl; 
  
  // grad U
  for(int i =0; i<9; i++)
    F[i] = in[3+i];
  
  // grad U + I
  F[0] += 1.0; 
  F[4] += 1.0;
  F[8] += 1.0;
  
  // C = F^T F
  MatrixMult(F, F, C, 't', 'n');
  
  // I1 = tr(C)
  I1 = C[0] +C[4] + C[8];
  
  // I3=det(C), Blas will retun LU decomposed value in C
  I3 = MatrixDeterminant(C, 3);
  I3pow = pow(I3, -1.0/3.0);
  
  // compute C^{-1}
  MatrixInverse(C, 3, 1);
  
  I1_bar = I1 * I3pow;
  
  for(i=0; i<9; i++)
   out[i] = 0.;
       
  k1 = c1 + c2*I1_bar - 6.*c2;
  k2 = 2 * (I3-1.) * I3;
  for(i=0; i<9; i++)
   {
    out[i] = k1*(I3pow - (I1_bar*C[i])/3) + (k2* C[i])/D1;;
   }
  
  
  OUT = out+9;
  k1 = I1_bar*(c1 + 2.*c2*I1_bar- 6.*c2)/6.0;
  for (i = 0; i < 3; i++)
   {
    disp = 3*i;
    for(j = 0; j < 3; j++)
     {
      J = 3*j;
      disp +=J;     
      cij = C[J+i];
      k2 = 2*c2*(I3pow - I1_bar*cij/3.0);
      for(k =0; k< 3; k++)
       {
        K = 3*k;
        disp +=K;         
        cik= C[K +i];
        cjk= C[K +j];
        for(l=0; l<3 ;l++)
         {
          L = 3*l;
          ckl= C[L +k];
          cjl= C[L +j];
          cil= C[L +i];
        
         OUT[disp + l] = k2*(I3pow - I1_bar * ckl/3.0)  +
                  (-I3pow * ckl)/3.0 - (I3pow - I1_bar * ckl/3.0)*cij - k1*(cik*cjl + cil*cjk) +
                  (2./D1) * ( I3*(I3-1)*ckl*cij + I3*I3*ckl*cij + I3*(I3-1)*(cik*cjl + cil*cjk)/2.0);
         }
       }
    }
  }

    for(i =0; i<9; i++)
      out[90+i] = F[i];
    
//   out[9] = in[0]; // x - coordinate  
//   out[10] = in[1]; // y - coordinate  
//   out[11] = in[2]; // z - coordinate 
//   for (i = 0; i < 90; i++)
//     cout << "out " << i << ": " << out[i] <<endl;
//     
//  cout << " stress cal " <<endl;
//  exit(0);
 
}

double Piola_Kir(double *param, double test100, double test010, double test001, double test000, double ansatz100, double ansatz010, double ansatz001, double ansatz000, int row, int col)
{
 const int c1=1.0;
 const int c2=1.0;
 const int D1=1.0;
	
 double I[9] = {1, 0, 0, 0,1,0, 0,0,1};
 double *F = new double[9];
 double * F_der = new double[9]; 
 double *C = new double[9];
 double *C_inv = new double[9];
 double I_first, I_third, I_inv, I_1_bar, cij,S_ij =0, val;
 double *temp_1 = new double[9];
 double *temp_2 = new double[9];
 double *S= new double[9];
 double *S_der = new double[9];
 double ckl, cik, cjl, cil, cjk;
 double entry;

 for(int i =0; i<9; i++)
  F[i] = param[3+i] + I[i] ;
 
 MatrixMult(F, F, C, 't', 'n');
 I_first = C[0] +C[4] + C[8];
 
 memcpy (temp_1, C, sizeof(temp_1)); 	
//  I_third = MatrixDeterminant(temp_1, 3);
 I_inv = pow(I_third, -1.0/3.0);
 memcpy (C_inv, C, sizeof(temp_1));
//  MatrixInverse(C_inv, 3, 1);
 I_1_bar = I_first * I_inv;
	
 for (int i = 0; i < 3; i++)
 {
  for(int j = 0; j < 3; j++)
  {
   cij = temp_1[3*j + i];
   val = (c1 + c2*I_1_bar - 6*c2)*(I_inv - (I_1_bar*cij)/3) + (2 * (I_third-1) * I_third * cij)/D1;
   S[3*j + i] = val;
  }
 }
 
 if(row==0)
 {
  for(int j =0 ; j <3; j++)
   {
    switch(j)
    {			
     case 0 : val = ansatz100; 
     break;
     case 1 : val = ansatz010; 
     break;				
     case 2 : val = ansatz001; 
     break;			 
    }
   F_der[3*j + row] = val;
   }
 }
 else if(row==1)
 {
  for(int j =0 ; j <3; j++)
  {
   switch(j)
   {
    case 0 : val = ansatz100; 
    break;
    case 1 : val = ansatz010; 
    break;
    case 2 : val = ansatz001; 
    break;
   }
   F_der[3*j + row] = val;
  }
 }	
else if(row==2)
 {
  for(int j =0 ; j <3; j++)
  {
   switch(j)
   {
   case 0 : val = ansatz100; 
   break;
   case 1 : val = ansatz010; 
   break;
   case 2 : val = ansatz001; 
   break;
   }
   F_der[3*j + row] = val;
  }
 }
 MatrixMult(F_der, F, temp_1, 't', 'n');
 MatrixMult(F, F_der, temp_2, 't', 'n');
 double *H = new double[9];
 for(int i =0; i<=8; i++) H[i] = temp_1[i] + temp_2[i] ;
 
 for (int i = 0; i < 3; i++)
 {
  for(int j = 0; j < 3; j++)
  {
   cij = C_inv[3*j +i];
   S_ij =0;
   for(int k =0; k< 3; k++)
   {
    for(int l=0; l<3 ;l++)
    {
     ckl= C_inv[3*l +k];
     cik= C_inv[3*k +i];
     cjl= C_inv[3*l +j];
     cil= C_inv[3*l +i];
     cjk= C_inv[3*k +j];
    double a = 2*c2*(I_inv - I_1_bar*cij/3.0)*(I_inv - I_1_bar * ckl/3.0);
    double b = (-1.0/3.0) * ( I_inv * ckl - (I_inv - (-1.0/3.0)*I_1_bar * ckl)*cij - (1.0/2.0)*(cik*cjl + cil*cjk)) * (c1 + 2*c2*I_1_bar- 6*c2);
    double c = (2/D1)*(I_third*(I_third-1)*ckl*cij + pow(I_third, 2.0)*ckl*cij +(I_third/2.0)*(I_third-1)*(cik*cjl + cil*cjk));
    val = a + b + c;
    S_ij += val * H[3*j +i];
    }
   }
   S_der[3*j +i] = S_ij;
  }
 }
 MatrixMult(F_der, S, temp_1, 'n', 'n');
 MatrixMult(F, S_der, temp_2, 'n', 'n');
 for(int i =0; i<=8; i++) H[i] = temp_1[i] + temp_2[i] ;
 switch(col)
 {
  case 0: entry = H[3*0 +0]*test100 + H[3*1 +0]*test010 + H[3*2 +0]*test001;				
  break;
  case 1: entry = H[3*0 +1]*test100 + H[3*1 +1]*test010 + H[3*2 +1]*test001;
  break;
  case 2: entry = H[3*0 +2]*test100 + H[3*1 +2]*test010 + H[3*2 +2]*test001;	
  break;		
 }
 delete [] H;
 delete [] F;   delete [] F_der;
 delete [] C;  delete [] C_inv; 
 delete [] S;   delete [] S_der;
 return entry;
}

//End of Piola_Kir function

