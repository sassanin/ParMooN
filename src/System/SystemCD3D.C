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
* @brief     source file for TSystemCD3D
* @author    Sashikumaar Ganesan
* @date      23.01.15
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemCD3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <AuxParam3D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <AssembleMat3D.h>

#include <stdlib.h>
#include <string.h>

#include <Solver.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

//#ifdef _MPI
TSystemCD3D::TSystemCD3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver)
{
  int i;
  int profiling = TDatabase::ParamDB->timeprofiling;
  /** need it for solver */
  sqmatrices = (TSquareMatrix **)SQMATRICES;
  
  //set number of multigrid levels
  N_Levels = N_levels;
  
  //store the FEspace
  FeSpaces = fespaces;
  
  SolArray = sol;
  RhsArray = rhs;
  
#ifdef _MPI
  Comm = TDatabase::ParamDB->Comm;
#endif
  //set the discretization type
  Disctype = disctype;
  
  //set the solver type
  SOLVER = solver;

  N_DOF = FeSpaces[N_Levels-1]->GetN_DegreesOfFreedom();  
  
  // build matrices
  // first build matrix structure
  sqstructure = new TSquareStructure3D*[N_Levels];
  sqmatrixA = new TSquareMatrix3D*[N_Levels];
  
  //instance of the Assemble class
  AMatRhsAssemble = new TAssembleMat3D *[N_Levels];
  
  if(SOLVER==AMG_SOLVE || SOLVER==DIRECT)
   {
    Start_Level=N_Levels-1;
   }
  else 
   {
    Start_Level=0;
   }
  
  for(i=Start_Level;i<N_Levels;i++)
  {
//    if(SOLVER==GMG)
//     OutPut("MULTIGRID LEVEL : " << i<<endl;)
    
   sqstructure[i] = new TSquareStructure3D(FeSpaces[i]);
   
   if(SOLVER==DIRECT || SOLVER==GMG)
   { 
     sqstructure[i]->Sort(); 
     
  } // sort column numbers: numbers are in increasing order
   else if(SOLVER==AMG_SOLVE)
   { sqstructure[i]->SortDiagFirst(); }
   
   /** A is the stiffness/system mat for stationary problem   */
   sqmatrixA[i] = new TSquareMatrix3D(sqstructure[i]);  
   N_Matrices = 1;
      
   OutPut(endl);   
  }// for(;i<N_Levels;i++)
  
  
#ifdef _MPI
  double t1,t2,tdiff;
  
  ParMapper = new TParFEMapper3D*[N_levels]; 
  ParComm   = new TParFECommunicator3D*[N_levels];
  
  if(profiling)  t1 = MPI_Wtime();
  for(i=Start_Level;i<N_levels;i++)
   {   
        ParMapper[i] = new TParFEMapper3D(1, FeSpaces[i], sqstructure[i]->GetRowPtr(), sqstructure[i]->GetKCol());
        ParComm[i]   = new TParFECommunicator3D(ParMapper[i]);
   }// for(i=0;i<N_levels;i++)

   int out_rank=TDatabase::ParamDB->Par_P0;
   int rank;
   MPI_Comm_rank(Comm, &rank);
   
   if(profiling)
   {
     t2 = MPI_Wtime();
     tdiff = t2-t1;
     MPI_Reduce(&tdiff, &t1, 1, MPI_DOUBLE, MPI_MIN, out_rank, Comm);
     if(rank == out_rank)
     {
      printf("Time taken for FeSpace SubDomain dof mapping for all levels is %e\n", t1);
     }
   }
   
  TCollection *coll = fespaces[N_levels-1]->GetCollection();
  int owncells = coll->GetN_OwnCells();
  int problemSize=0;
  MPI_Reduce(&owncells, &problemSize, 1, MPI_INT, MPI_SUM, out_rank, Comm);
  if(rank==0)
    OutPut( "total own cells over all sub domains : " << problemSize << endl);
  
  problemSize=0;
  int n_master = ParComm[N_levels-1]->GetN_Master();
  MPI_Reduce(&n_master, &problemSize, 1, MPI_INT, MPI_SUM, out_rank, Comm);
  if(rank==0)
    OutPut( "total own dofs over all sub domains : " << problemSize << endl);
  
#endif
 
   //initialize multigrid solver
   if(SOLVER==GMG)
   {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    MG = new TMultiGrid3D(1, 2, Parameters);
    
    // determine number of auxiliary arrays
    if ( (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
         || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR) )
     {  N_aux= 4; }
        else
     {  N_aux= 2; }   
   
    // build preconditioner
    switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
     {
      case 1:
            prec = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL, 0, N_DOF, 1
#ifdef _MPI   
                                  ,ParComm[N_Levels-1]
#endif    
                                  );
	    break;	    
      case 5:
            prec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL, 0, N_DOF, MG, 1);
            Itmethod_sol = new double[N_DOF];
            Itmethod_rhs = new double[N_DOF];    
            break;
      default:
            OutPut("Unknown preconditioner !!!" << endl);
            exit(4711);
     }     

       switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
        {
          case 11:
            Itmethod = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, prec, 0, N_DOF, 1
#ifdef _MPI   
                               , ParComm[N_Levels-1]
#endif
	);
          break;
          case 16:
            Itmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, prec, 0, N_DOF, 1
#ifdef _MPI   
                               , ParComm[N_Levels-1]
#endif
	);
          break;
          default:
            OutPut("Unknown solver !!!" << endl);
            exit(4711);
        }     
     
   }// if(solver==GMG)  
//   exit(0);
} //TSystemCD3D::TSystemCD3D
//#endif

TSystemCD3D::~TSystemCD3D()
{
  int i;
  
  for(i=Start_Level;i<N_Levels;i++)
   {
    delete sqstructure[i];
    delete sqmatrixA[i];   
   }
   
    delete [] sqstructure;
    delete [] sqmatrixA;
  
  if (SOLVER==GMG && TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
   {
    delete [] Itmethod_sol;
    delete [] Itmethod_rhs;
   }
}
  
  
void TSystemCD3D::Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                              TAuxParam3D *aux)
{
#ifdef _MPI
     if(SOLVER == DIRECT)
    {
     SQMATRICES[0] = sqmatrixA[N_Levels-1];
     DS = new TParDirectSolver(ParComm[N_Levels-1],NULL,SQMATRICES,NULL);
    }
#endif

#ifdef _OMPONLY
   if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
   {
     DS = new TParDirectSolver(sqmatrixA[N_Levels-1]);
   }
#endif 

 int i;
 
  BoundaryConditions[0] = BoundCond;
  BoundaryValues[0] = BoundValue;
  
//   TDiscreteForm3D *DiscreteFormUpwind;  
  TDiscreteForm3D *DiscreteFormGalerkin;
//   TDiscreteForm3D *DiscreteFormSDFEM;
//   TDiscreteForm3D *DiscreteFormGLS;  

   InitializeDiscreteForms(DiscreteFormGalerkin, BilinearCoeffs);  

    switch(Disctype)
     {
      case GALERKIN:
//       case LOCAL_PROJECTION:
           DiscreteFormARhs = DiscreteFormGalerkin;
      break;

//       case SUPG:
//            DiscreteFormARhs = DiscreteFormSDFEM;
//       break;
// 
//       case UPWIND:
//            DiscreteFormARhs = DiscreteFormUpwind;
//       break;      
//       
//       case GLS:
//            DiscreteFormARhs = DiscreteFormGLS;
//       break;
      
      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }  
     
   // initilize the assemble class  
   if(aux==NULL)
    { aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }   
   
    for(i=Start_Level;i<N_Levels;i++)
    { 
     fesp[0] = FeSpaces[i];
     ferhs[0] = FeSpaces[i];  
     
     RHSs[0] = RhsArray[i];  
     SQMATRICES[0] = sqmatrixA[i];
     
     // array of assemble objects
     AMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs, 
                              DiscreteFormARhs, BoundaryConditions, BoundaryValues, aux);
     AMatRhsAssemble[i]->Init();
  
     
     //setup the multigrid solver
     if(SOLVER==GMG)
      {
#ifdef _MPI  
       MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], ParComm[i], ParMapper[i], N_aux, NULL);
#else
       MGLevel = new TMGLevel3D(i, SQMATRICES[0], RHSs[0], SolArray[i], N_aux, NULL);
#endif
       MG->AddLevel(MGLevel);      
      }       
    } // for(i=Star
    
} // TSystemCD3D::Init


void TSystemCD3D::Assemble()
{
  int i, N_DOF_low, N_Active;

   for(i=Start_Level;i<N_Levels;i++)
    {    
     N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
     N_Active =  FeSpaces[i]->GetActiveBound();
    
     // initialize matrices and rhs
     AMatRhsAssemble[i]->Reset(); 

     // assemble
     AMatRhsAssemble[i]->Assemble3D();
  
     // set rhs for Dirichlet nodes
     memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);   
    } //  for(i=Start_Level;i<N_Levels;i++)    

//have to shift this in pardirectsolver    
#ifdef _OMPONLY     
    if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
      DS->AssembleMatrix(sqmatrixA[N_Levels-1]);
#endif
    
} // void TSystemCD3D::Assemble(T


void TSystemCD3D::Solve(double *sol, double *rhs)
{ 
    switch(SOLVER)
     {
      case AMG_SOLVE:
        Solver(sqmatrixA[N_Levels-1], rhs, sol);
      break;

      case GMG:
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          memcpy(Itmethod_sol, sol, N_DOF*SizeOfDouble);
          memcpy(Itmethod_rhs, rhs, N_DOF*SizeOfDouble);
         }
        else
         {
          Itmethod_sol = sol;
          Itmethod_rhs = rhs;
         }
         
         // solve linear system
        Itmethod->Iterate(sqmatrices, NULL, Itmethod_sol, Itmethod_rhs);

        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          memcpy(sol, Itmethod_sol, N_DOF*SizeOfDouble);
          memcpy(rhs, Itmethod_rhs, N_DOF*SizeOfDouble);
         }
      break;

      case DIRECT:
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
	DirectSolver((TSquareMatrix*)sqmatrixA[N_Levels-1], rhs, sol);
#endif
      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
}
