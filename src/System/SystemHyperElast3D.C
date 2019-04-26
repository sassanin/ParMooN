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
 int Hyper_FEFctIndex[12] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
 int Hyper_BeginParam[1] = { 0 };  
 MultiIndex3D Hyper_FEMultiIndex[12] = { D000, D000, D000, D100, D100, D100,
                                         D010, D010, D010, D001, D001, D001 };  
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
  
  /** set the Discreteforms */
//   this->InitHyperDiscreteForms(DiscreteFormGalerkin);
  
       

  
  DiscreteFormGalerkin = NULL; 
  int N_Terms = 4;
  N_Rhs = 3;
  int N_Matrices = 9;
 
  char GalerkinString[] = "Galerkin";
  char rhs[] = "rhs";
  char all[] = "all";    
  
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
       
                 
//       bool *SecondDer;
//       SecondDer = DiscreteFormARhs->GetNeeds2ndDerivatives();
//       cout << " SecondDer " << SecondDer[0] << endl; 
//       exit(0);   
       
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
 
//           bool *SecondDer;
//       SecondDer = DiscreteFormARhs->GetNeeds2ndDerivatives();
//       cout << " SecondDer " << SecondDer[0] << endl; 
      
      
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



// void TSystemNSE3D::AssembleNonLinear(double **sol, double **rhs)
// {
//  int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
//  int N_U_Current, N_Active_Current, N_DirichletDof;
//   
//    for(i=Start_Level;i<N_Levels;i++)
//     {    
//       N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();
//       N_Active_Current  = U_Space[i]->GetActiveBound();     
//       N_DirichletDof = N_U_Current - N_Active_Current;
//       
//       // reset the nonliner matrices
//       AMatAssembleNonLinear[i]->Reset();
//             
//       // assemble the nonlinear matrix */      
//       AMatAssembleNonLinear[i]->Assemble3D();      
//       
//        /** upwind */
//       if( (Disctype==UPWIND) && (TDatabase::ParamDB->PROBLEM_TYPE!=STOKES) )
//         {
//          this->UpdateUpwind(i);
//         }     
//         
//       // slip with boundary condition
//       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
//       { 
// // 	AMatRhsAssemble[i]->AssembleNavierSlip(); 
//            
//         fesp[0] =  U_Space[i];
//         N_FESpaces = 1;
//         N_SquareMatrices = 3;
//         N_RectMatrices = 0;
//         N_Rhs = 3;
// 
//         SQMATRICES[0] = SqmatrixA11[i];
//         SQMATRICES[1] = SqmatrixA22[i];
//         SQMATRICES[2] = SqmatrixA33[i];
// 
//         RHSs[0] = RhsArray[i];
//         RHSs[1] = RhsArray[i] + N_U_Current;
//         RHSs[2] = RhsArray[i] + 2*N_U_Current;
// 
//         fesprhs[0] =  U_Space[i];
//         fesprhs[1] =  U_Space[i];
//         fesprhs[2] =  U_Space[i];
//   
//         Assemble3DSlipBC(N_FESpaces, fesp,
//                          N_SquareMatrices, SQMATRICES,
//                          N_RectMatrices, NULL,
//                          N_Rhs, RHSs, fesprhs,
//                          NULL,
//                          BoundaryConditions,
//                          BoundaryValues,
//                          NSEaux[i]);
// 
//        }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=         
//      
//      /** no change in rhs, so no need to update */
//       // set rhs for Dirichlet nodes 
//       memcpy(SolArray[i]+N_Active_Current, RhsArray[i]+N_Active_Current, N_DirichletDof*SizeOfDouble);
//       memcpy(SolArray[i]+N_U_Current+N_Active_Current, RhsArray[i]+N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble); 
//       memcpy(SolArray[i]+2*N_U_Current+N_Active_Current, RhsArray[i]+2*N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble);     
//         
//       } //
//       
// } //TSystemNSE3D::AssembleNonLinear(
// 
// 
// void TSystemNSE3D::GetResidual(double *sol, double *rhs, double *res, double &impuls_residual, double &residual)
// {
// 
//      switch(NSEType)
//       {
//         case 1:
//           SQMATRICES[0] = SqmatrixA11[N_Levels-1];
//           MATRICES[0] = MatrixB1[N_Levels-1];
//           MATRICES[1] = MatrixB2[N_Levels-1];
//           MATRICES[2] = MatrixB3[N_Levels-1];
//         break;
// 
//         case 2:
//           SQMATRICES[0] = SqmatrixA11[N_Levels-1];
//           MATRICES[0] = MatrixB1[N_Levels-1];
//           MATRICES[1] = MatrixB2[N_Levels-1];
//           MATRICES[2] = MatrixB3[N_Levels-1];
//           MATRICES[3] = MatrixB1T[N_Levels-1];
//           MATRICES[4] = MatrixB2T[N_Levels-1];
//           MATRICES[5] = MatrixB3T[N_Levels-1];
// 	  
//         break;
// 
//         case 3:
//           SQMATRICES[0] = SqmatrixA11[N_Levels-1];
//           SQMATRICES[1] = SqmatrixA12[N_Levels-1];
//           SQMATRICES[2] = SqmatrixA13[N_Levels-1];	  
//           SQMATRICES[3] = SqmatrixA21[N_Levels-1];
//           SQMATRICES[4] = SqmatrixA22[N_Levels-1];
//           SQMATRICES[5] = SqmatrixA23[N_Levels-1]; 
//           SQMATRICES[6] = SqmatrixA31[N_Levels-1];
//           SQMATRICES[7] = SqmatrixA32[N_Levels-1];
//           SQMATRICES[8] = SqmatrixA33[N_Levels-1];  
// 
//           MATRICES[0] = MatrixB1[N_Levels-1];
//           MATRICES[1] = MatrixB2[N_Levels-1];
//           MATRICES[2] = MatrixB3[N_Levels-1];
//         break;
// 
//         case 4:
//           SQMATRICES[0] = SqmatrixA11[N_Levels-1];
//           SQMATRICES[1] = SqmatrixA12[N_Levels-1];
//           SQMATRICES[2] = SqmatrixA13[N_Levels-1];	  
//           SQMATRICES[3] = SqmatrixA21[N_Levels-1];
//           SQMATRICES[4] = SqmatrixA22[N_Levels-1];
//           SQMATRICES[5] = SqmatrixA23[N_Levels-1]; 
//           SQMATRICES[6] = SqmatrixA31[N_Levels-1];
//           SQMATRICES[7] = SqmatrixA32[N_Levels-1];
//           SQMATRICES[8] = SqmatrixA33[N_Levels-1];  
//           MATRICES[0] = MatrixB1[N_Levels-1];
//           MATRICES[1] = MatrixB2[N_Levels-1];
//           MATRICES[2] = MatrixB3[N_Levels-1];
//           MATRICES[3] = MatrixB1T[N_Levels-1];
//           MATRICES[4] = MatrixB2T[N_Levels-1];
//           MATRICES[5] = MatrixB3T[N_Levels-1];
//         break;
//       } //  switch(NSEType)  
// 
// #ifdef _MPI      
//     ParComm_U[N_Levels-1]->CommUpdate(sol);
//     ParComm_P[N_Levels-1]->CommUpdate(sol+3*N_U);
// #endif
//        
// 
//     
//     Defect(sqmatrices, matrices, sol, rhs, res); 
//   
// 
//      
// #ifdef _MPI
//    double residual_scalar = 0.0;
//    double sum =0.;
//    int i,j,rank;
//    MPI_Comm_rank(Comm, &rank);
//    int *master = ParComm_U[N_Levels-1]->GetMaster();
//    
//    for(i=0;i<N_U;i++)
//    {
//      if(master[i]!=rank)    continue;
//       BdComp
//       residual_scalar += res[i      ]*res[i      ];
//       residual_scalar += res[i+  N_U]*res[i+  N_U];
//       residual_scalar += res[i+2*N_U]*res[i+2*N_U];
// 
//     }
// 
//     MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
//     impuls_residual = (sum);
// 
//     master = ParComm_P[N_Levels-1]->GetMaster();
//     for(i=0;i<N_P;i++)
//     {
//       if(master[i]!=rank)    continue;
//       
//       residual_scalar += res[i+3*N_U]*res[i+3*N_U];
// 
//     }
//     
//     sum = 0;
//     MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
//     residual = (sum);
// 
// #else
//     impuls_residual  =  Ddot(3*N_U, res, res);
//     residual         =  Ddot(N_TotalDOF, res, res); 
// 
// #endif
//    
// } // TSystemNSE3D::GetResidual
// 
// void TSystemNSE3D::Solve(double *sol, double *rhs)
// {
//   int N_LinIter=0;
//   double summ = 0;
//   double residual,residual_scalar = 0.0;
//   double sum =0.;
//   int i,j,rank;
//   int *master;
//    
//     switch(SOLVER)
//      {
//       case AMG_SOLVE:
//         cout << "AMG_SOLVE not yet implemented " <<endl;
//       break;
// 
//       case GMG:
// 
//           if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
//            {
//             memcpy(Itmethod_sol, sol, N_TotalDOF*SizeOfDouble);
//             memcpy(Itmethod_rhs, rhs, N_TotalDOF*SizeOfDouble);
//            }
//           // solve the linear system
//           N_LinIter += Itmethod->Iterate(sqmatrices, matrices, Itmethod_sol, Itmethod_rhs);
// 	  
//           if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
//            {
//             memcpy(sol, Itmethod_sol, N_TotalDOF*SizeOfDouble);
//             memcpy(rhs, Itmethod_rhs, N_TotalDOF*SizeOfDouble);
//            }
//           MG->RestrictToAllGrids();
// 
//       break;
// 
//       
//       case DIRECT:
//  
//         switch(NSEType)
//          {
//           case 1:
//            cout << "Solver not included for NSTYPE 1 in this version" <<endl;
//             cout << "try NSTYPE 2 or 4 " <<endl;   
// 	    exit(0);
//           break;
// 
//           case 2:
// 
// 	    
// #ifdef _MPI
// 	    DS->Solve(sol, rhs, true);
// 	    
// #else	    
// 
// 	    
// 	    
// 	    DirectSolver(SqmatrixA11[N_Levels-1], MatrixB1T[N_Levels-1], MatrixB2T[N_Levels-1], MatrixB3T[N_Levels-1],
//                           MatrixB1[N_Levels-1], MatrixB2[N_Levels-1], MatrixB3[N_Levels-1], rhs, sol);
// 
// #endif
// 
// 	    
//           break;
// 
//           case 3:
//            cout << "Solver not included for NSTYPE 3 in this version" <<endl;
//             cout << "try NSTYPE 2 or 4 " <<endl;   
//             exit(0);
//           break;
// 
//           case 4:
// #ifdef _MPI
//         	DS->Solve(sol, rhs, true);
// #endif
// 	
// #ifdef _OMPONLY
// 	       if(TDatabase::ParamDB->DSType == 1)
// 	         DS->Solve(sol, rhs, true);
// 	       else{
// 	         OutPut("Select Proper Solver" << endl);
// 	         exit(0);
// 	       }
// #endif
// 
// #ifdef _SEQ
// 
//          
//      
//              DirectSolver(SqmatrixA11[N_Levels-1], SqmatrixA12[N_Levels-1], SqmatrixA13[N_Levels-1], 
//                           SqmatrixA21[N_Levels-1], SqmatrixA22[N_Levels-1], SqmatrixA23[N_Levels-1],  
//                           SqmatrixA31[N_Levels-1], SqmatrixA32[N_Levels-1], SqmatrixA33[N_Levels-1],  
//                           MatrixB1T[N_Levels-1], MatrixB2T[N_Levels-1], MatrixB3T[N_Levels-1],
//                           MatrixB1[N_Levels-1], MatrixB2[N_Levels-1], MatrixB3[N_Levels-1], rhs, sol,0);
// 
// #endif
//           break;
//          } //  switch(NSEType) 
// 
//       break;      
//  
//       default:
//             OutPut("Unknown Solver" << endl);
//             exit(4711);;
//      }    
//   
// }
// 
// void TSystemNSE3D::MeasureErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2, DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP,
//                                   double *u_error, double *p_error)
// {
//   double errors[4];
//   fesp[0] =  U_Space[N_Levels-1];
// 
//     fefct[0] = Velocity[N_Levels-1]->GetComponent(0);
//     fefct[1] = Velocity[N_Levels-1]->GetComponent(1);
//     fefct[2] = Velocity[N_Levels-1]->GetComponent(2);     
// 
//   
//      if(NSEaux_error[N_Levels-1]==NULL)
//       NSEaux_error[N_Levels-1] =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo, NSN_FEValuesVelo,
//                                     fesp, fefct, NSFctVelo, NSFEFctIndexVelo, NSFEMultiIndexVelo,
//                                     NSN_ParamsVelo, NSBeginParamVelo);
// 
//      // errors in first velocity component
//       fefct[0]->GetErrors(ExactU1, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1, U_Space+(N_Levels-1), errors);
//    
//       u_error[0] = errors[0];
//       u_error[1] = errors[1];
//       
//      // errors in second velocity component
//       fefct[1]->GetErrors(ExactU2, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1, U_Space+(N_Levels-1), errors);
//       u_error[2] = errors[0];
//       u_error[3] = errors[1];     
//      
//       // errors in third velocity component
//       fefct[2]->GetErrors(ExactU3, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1, U_Space+(N_Levels-1), errors);
//       u_error[4] = errors[0];
//       u_error[5] = errors[1];       
//      
//      
//       // errors in pressure
//       Pressure[N_Levels-1]->GetErrors(ExactP, 4, NSAllDerivatives, 2, L2H1Errors, NULL, NSEaux_error[N_Levels-1], 1,  P_Space+(N_Levels-1), errors);     
//       p_error[0] = errors[0];
//       p_error[1] = errors[1]; 
// }
//     
// 
// void TSystemNSE3D::FindFreeSurfJoints(int level, int Region_ID)
// {
//   int i,j, N_Joints;
//   int N_Cells, N_FreeJoints;
//   
//   TBaseCell *cell, *neigh, *neigh0, *neigh1;
//   TJoint *joint;
//   TCollection *Coll;
//   TIsoJointEqN *isojoint;
// 
//   Coll = U_Space[level]->GetCollection();  
//   N_Cells = Coll->GetN_Cells();
//   N_FreeJoints=0;
// 
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = Coll->GetCell(i);
//     if(cell->GetRegionID() == Region_ID)
//     {
//       N_Joints = cell->GetN_Joints();
//       for(j=0;j<N_Joints;j++)
//        {
//         joint = cell->GetJoint(j);
//         if(joint->GetType()==InterfaceJoint3D || joint->GetType()==IsoInterfaceJoint3D || joint->GetType()==IsoBoundFace)
//          {  
//           N_FreeJoints++;
//          }
//        } // endfor j
//     } // end ID == 1
//   } // endfor i
// 
//   N_FreeSurfFaces[level] = N_FreeJoints;
//   
//   if(N_FreeJoints)
//    {
//     FreeSurfCellNumbers[level] = new int[N_FreeJoints];
//     FreeSurfJointNumbers[level] = new int[N_FreeJoints];    
//    }
//    
//   N_FreeJoints = 0;
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = Coll->GetCell(i);
//     if(cell->GetRegionID() == Region_ID)
//     {
//       N_Joints = cell->GetN_Joints();
//       for(j=0;j<N_Joints;j++)
//        {
//         joint = cell->GetJoint(j);
//         if(joint->GetType()==InterfaceJoint3D || joint->GetType()==IsoInterfaceJoint3D || joint->GetType()==IsoBoundFace)
//          { 
//           FreeSurfCellNumbers[level][N_FreeJoints] = i;
//           FreeSurfJointNumbers[level][N_FreeJoints] = j;
//           N_FreeJoints++;
//          }
//        } // endfor j
//     } // end ID == 1
//   } // endfor i
//   
//   
// //   cout << " FindFreeSurfJoints done " << N_FreeJoints << endl;
// //   exit(0);  
// }



void TSystemHyperElast3D::InitHyperDiscreteForms(TDiscreteForm3D *DiscreteFormGalerkin)
{
/*       
  char GalerkinString[] = "Galerkin";
  char rhs[] = "rhs";
  char all[] = "all";
  
  DiscreteFormGalerkin = NULL; 
   
  int N_Terms = 4;
  MultiIndex3D Derivatives[4] = { D100, D010, D001, D000};
  int SpaceNumbers[4] = { 0, 0, 0, 0};
  int N_Matrices = 9;
  int RowSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int ColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int N_Rhs = 3;
  int RhsSpace[3] = { 0, 0, 0 };
  
  
  DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,  N_Terms, Derivatives, SpaceNumbers,
                                             N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace, Galerkin3D, LinCoeffs[0], NULL);         
  
  
      bool *SecondDer;
      SecondDer = DiscreteFormGalerkin->GetNeeds2ndDerivatives();
    
                 cout << "DiscreteFormGalerkin SecondDer " << SecondDer[0] << endl; */
          
                 
}


void TSystemHyperElast3D::InitHyperAuxParm(int i)
{
  int Hyper_N_FESpace = 1;
  int Hyper_N_FEFunction = 3;
  int Hyper_N_ParamFct = 1;
  int Hyper_N_FEValues = 12;
  int Hyper_N_ParamValues = 15;                                      
                                            
  fesp_aux[0] =  U_Space[i];
     
  fefct_aux[i*3 ] = Displacement[i]->GetComponent(0);
  fefct_aux[(i*3)+1] = Displacement[i]->GetComponent(1);
  fefct_aux[(i*3)+2] = Displacement[i]->GetComponent(2);

  Hyperaux[i] =  new TAuxParam3D(Hyper_N_FESpace, Hyper_N_FEFunction, Hyper_N_ParamFct, Hyper_N_FEValues,
                                fesp_aux, fefct_aux+(i*3), Hyper_ParamFct, Hyper_FEFctIndex, Hyper_FEMultiIndex,
                                Hyper_N_ParamValues, Hyper_BeginParam);
    
  
  
//         cout << " Hyper_FEMultiIndex   3D: " << Hyper_FEMultiIndex[0] << endl; 
  
  
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
  double u1, u2, u3, u1x, u2x, u3x;
  double u1y, u2y, u3y, u1z, u2z, u3z;

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

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u


  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
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
  u3z = param[11]; // u3old  

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
            
      MatrixRowA11[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 0, 0);
      MatrixRowA12[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 1, 0);
      MatrixRowA13[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 2, 0);     
      MatrixRowA21[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 0, 1);          
      MatrixRowA22[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 1, 1);          
      MatrixRowA23[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 2, 1);          
      MatrixRowA31[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 0, 2);          
      MatrixRowA32[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 1, 2);          
      MatrixRowA33[j] += Mult * Piola_Kir(param, test100, test010, test001, test000, ansatz100,ansatz010, ansatz001, ansatz000, 2, 2);     
      
    } // endfor j
  } // endfor i

}

void HyperParamsVelo(double *in, double *out)
{
    
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[0]; // x - coordinate  
  out[13] = in[1]; // y - coordinate  
  out[14] = in[2]; // z - coordinate 
    
}


//Evaluation of Matrix Entries 

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
 double I_first, I_third, I_inv, I_1_bar, c_ij,S_ij =0, val;
 double *temp_1 = new double[9];
 double *temp_2 = new double[9];
 double *S= new double[9];
 double *S_der = new double[9];
 double c_kl, c_ik, c_jl, c_il, c_jk;
 double entry;

 for(int i =0; i<=8; i++)
  F[i] = param[3+i] + I[i] ;
 
 MatrixMult(F, F, C, 't', 'n');
 I_first = C[0] +C[4] + C[8];
 
 memcpy (temp_1, C, sizeof(temp_1)); 	
 I_third = MatrixDeterminant(temp_1);
 I_inv = pow(I_third, -1.0/3.0);
 memcpy (C_inv, C, sizeof(temp_1));
 MatrixInverse(C_inv);
 I_1_bar = I_first * I_inv;
	
 for (int i = 0; i < 3; i++)
 {
  for(int j = 0; j < 3; j++)
  {
   c_ij = temp_1[3*j + i];
   val = (c1 + c2*I_1_bar - 6*c2)*(I_inv - (I_1_bar*c_ij)/3) + (2 * (I_third-1) * I_third * c_ij)/D1;
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
   c_ij = C_inv[3*j +i];
   S_ij =0;
   for(int k =0; k< 3; k++)
   {
    for(int l=0; l<3 ;l++)
    {
     c_kl= C_inv[3*l +k];
     c_ik= C_inv[3*k +i];
     c_jl= C_inv[3*l +j];
     c_il= C_inv[3*l +i];
     c_jk= C_inv[3*k +j];
    double a = 2*c2*(I_inv - I_1_bar*c_ij/3.0)*(I_inv - I_1_bar * c_kl/3.0);
    double b = (-1.0/3.0) * ( I_inv * c_kl - (I_inv - (-1.0/3.0)*I_1_bar * c_kl)*c_ij - (1.0/2.0)*(c_ik*c_jl + c_il*c_jk)) * (c1 + 2*c2*I_1_bar- 6*c2);
    double c = (2/D1)*(I_third*(I_third-1)*c_kl*c_ij + pow(I_third, 2.0)*c_kl*c_ij +(I_third/2.0)*(I_third-1)*(c_ik*c_jl + c_il*c_jk));
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

