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
* @brief     source file for TSystemTCD3D
* @author    Sashikumaar Ganesan
* @date      24.01.15
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemTCD3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <AuxParam3D.h>
#include <MultiGridScaIte.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <Solver.h>
#include <AssembleMat3D.h>

#include <stdlib.h>
#include <string.h>

TSystemTCD3D::TSystemTCD3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver):TSystemCD3D(N_levels, fespaces, sol, rhs, disctype, solver)
{
  int i;
  
  /** M mass and system matrix */
  sqmatrixM = new TSquareMatrix3D*[N_Levels];
  N_Matrices++;
  
  MMatRhsAssemble = new TAssembleMat3D*[N_Levels];
  
  for(i=Start_Level;i<N_Levels;i++)
   {
    sqmatrixM[i] = new TSquareMatrix3D(sqstructure[i]);       
   }
  
  /** working rhs, used in AssembleSystMat() */
  B = new double[N_DOF];
  defect = new double[N_DOF];
  
  gamma =0.;
  
//   /** time-consistent part of the SUPG matrix */
//   if(Disctype==SDFEM || Disctype==SUPG)
//    {
//     sqmatrixS = new TSquareMatrix3D(sqstructure); 
//     N_Matrices++;
//     
//     sqmatrixK = new TSquareMatrix3D(sqstructure); 
//     N_Matrices++; 
//    }
   
  SystMatAssembled  = FALSE;
    
  
  
  Values = NULL; KCol = NULL; Row = NULL; Symbolic=NULL, Numeric=NULL;  
  UpdateStiffnessMat=FALSE;
  FirstTime   = TRUE;
   
} // constructor


TSystemTCD3D::~TSystemTCD3D()
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
  
  delete [] B;
  delete [] defect;
  
//   delete sqmatrixM;
//   if(Disctype==SDFEM || Disctype==SUPG)
//    {
//      delete sqmatrixS;
//      delete sqmatrixK;    
//    }
}


void TSystemTCD3D::Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue,
                              TAuxParam3D *aux )
{
#ifdef _MPI
   if(SOLVER == DIRECT)
   {
     SQMATRICES[0] = sqmatrixM[N_Levels-1];
     TDS = new TParDirectSolver(ParComm[N_Levels-1],NULL,SQMATRICES,NULL);
   }
#endif

#ifdef _OMPONLY
   if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
   {
     DS = new TParDirectSolver(sqmatrixM[N_Levels-1]);
   }
#endif 

 int i; 
  
  BoundaryConditions[0] = BoundCond;
  BoundaryValues[0] = BoundValue;
  
  TDiscreteForm3D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm3D *DiscreteFormARhs_Galerkin; 
//   TDiscreteForm3D *DiscreteFormMRhs_SUPG;
//   TDiscreteForm3D *DiscreteFormARhs_SUPG;

  
  InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormRhs, BilinearCoeffs);
  
    switch(Disctype)
     {
      case GALERKIN:
//       case LOCAL_PROJECTION:
           DiscreteFormARhs = DiscreteFormARhs_Galerkin;
           DiscreteFormMRhs = DiscreteFormMRhs_Galerkin;
      break;
      
//       case SUPG:
//            DiscreteFormARhs = DiscreteFormARhs_SUPG;
//            DiscreteFormMRhs = DiscreteFormMRhs_SUPG;
//       break;
      
      default:
            OutPut("Unknown or not yet implemented DISCTYPE" << endl);
            exit(4711);;
     }  
   
   // initialize the assemble 
   if(aux==NULL)
    { aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }
    
    for(i=Start_Level;i<N_Levels;i++)
    { 
     fesp[0] = FeSpaces[i];
     ferhs[0] = FeSpaces[i];  
     
     RHSs[0] = RhsArray[i];
  
     // A Matrix     
     SQMATRICES[0] = sqmatrixA[i];
     AMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs, 
                              DiscreteFormARhs, BoundaryConditions, BoundaryValues, aux);
     AMatRhsAssemble[i]->Init();
  
     // M matrix     
     SQMATRICES[0] = sqmatrixM[i];
     MMatRhsAssemble[i] = new TAssembleMat3D(1, fesp, 1, SQMATRICES, 0, NULL, 1, RHSs, ferhs, 
                              DiscreteFormMRhs, BoundaryConditions, BoundaryValues, aux);
     MMatRhsAssemble[i]->Init();   
      
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
} // Init


void TSystemTCD3D::AssembleMRhs()
{
  //this is set to true for direct solver factorization
  factorize = true;
  
  int i, N_DOF_low, N_Active;

   for(i=Start_Level;i<N_Levels;i++)
    {      
     N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
     N_Active =  FeSpaces[i]->GetActiveBound();
  
     /** initialize matrices and rhs */
     MMatRhsAssemble[i]->Reset(); 
     
     // assemble
     MMatRhsAssemble[i]->Assemble3D();
     
     /** free the Mass mat array, no need in time loop */
     MMatRhsAssemble[i]->DeAllocate();
     
     /** set rhs for Dirichlet nodes */
     memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);         
    } //  for(i=Start_Level;i<N_Levels;i++)

} // TSystemMatScalar3D::AssembleMRhs 


void TSystemTCD3D::AssembleARhs()
{
  //this is set to true for direct solver factorization
  factorize = true;
  
  int i, N_DOF_low, N_Active;

   for(i=Start_Level;i<N_Levels;i++)
    {    
     N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
     N_Active =  FeSpaces[i]->GetActiveBound();

     /** reset the matrix and rhs */
     AMatRhsAssemble[i]->Reset(); 
    
     // assemble
     AMatRhsAssemble[i]->Assemble3D();     

     /** set rhs for Dirichlet nodes */
     memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);           
    }//   for(i=Start_Level;i<N_Le  
    
} // TSystemMatScalar3D::AssembleARhs 

void TSystemTCD3D::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol
#ifdef _MPI
                                             , double **Rhs_array
#endif
                                             )
{
    int i, N_Active;
    double tau;
    
    if(SystMatAssembled)
     {
      OutPut("System is has to be restored before calling AssembleSystMat! " <<endl);
      exit(0);
     }
    
    SQMATRICES[0] = sqmatrixM[N_Levels-1];

    N_Active =  FeSpaces[N_Levels-1]->GetActiveBound();     
    tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;  
    
    memset(B, 0, N_DOF*SizeOfDouble); 

    /** old rhs multiplied with current subtime step and theta3 on B */
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3,  oldrhs, B);    

    /** add rhs from current sub time step to rhs array B */
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4,  rhs, B);    

    /** M = M + (- tau*THETA2)A */
     MatAdd(sqmatrixM[N_Levels-1], sqmatrixA[N_Levels-1], - tau*TDatabase::TimeDB->THETA2);
     gamma = -tau*TDatabase::TimeDB->THETA2;  // set current factor of steady state matrix

     /** defect = M * oldsol */
     memset(defect, 0, N_DOF*SizeOfDouble);  
     MatVectActive(sqmatrixM[N_Levels-1], oldsol, defect); 
    //cout << "defect " << Ddot(N_Active, sol, sol)<< endl;      
    
     /** B:= B + defec  */
     Daxpy(N_Active, 1, defect, B);

     /** set Dirichlet values */
     memcpy(B+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);  
     memcpy(sol+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
          
     /** assemble the system matrix */
     for(i=Start_Level;i<N_Levels;i++)   
      {
       if(i==N_Levels-1)
         { MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma + tau*TDatabase::TimeDB->THETA1);}
        else
         { MatAdd(sqmatrixM[i], sqmatrixA[i], tau*TDatabase::TimeDB->THETA1);} 
         
#ifdef _MPI  
       SQMATRICES[0] = sqmatrixM[i];  
#endif
       }
     gamma = tau*TDatabase::TimeDB->THETA1;
     
//have to shift this in pardirectsolver     
#ifdef _OMPONLY     
    if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
      DS->AssembleMatrix(sqmatrixM[N_Levels-1]);
#endif
     
     SystMatAssembled  = TRUE;

} // AssembleSystMat

void TSystemTCD3D::RestoreMassMat()
{
 int i;

  if(SystMatAssembled)
   {
     // restore the mass matrix
     for(i=Start_Level;i<N_Levels;i++)  
      MatAdd(sqmatrixM[i], sqmatrixA[i], -gamma);
     
     gamma = 0.;
     SystMatAssembled  = FALSE;
   }
  else
  {
    cout << "System is not assembled to restore " <<endl;
    exit(0);
  }

}

void TSystemTCD3D::Solve(double *sol)
{  
    switch(SOLVER)
     {
      case AMG_SOLVE:
         Solver(sqmatrixM[N_Levels-1], B, sol);
      break;

      case GMG:
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          memcpy(Itmethod_sol, sol, N_DOF*SizeOfDouble);
          memcpy(Itmethod_rhs, B, N_DOF*SizeOfDouble);
         }
        else
         {
          Itmethod_sol = sol;
          Itmethod_rhs = B;
         }
      
         /** solve linear system */
        Itmethod->Iterate(sqmatrices, NULL, Itmethod_sol, Itmethod_rhs);
#ifdef _MPI
    if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==6)
         ParComm[N_Levels-1]->CommUpdateH2(Itmethod_sol);
#endif
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          memcpy(sol, Itmethod_sol, N_DOF*SizeOfDouble);
         }
      break;

      case DIRECT:
#ifdef _MPI
	TDS->Solve(sol, B, factorize);
// 	exit(0);
#endif

#ifdef _OMPONLY
	if(TDatabase::ParamDB->DSType == 1)
	  DS->Solve(sol, B, factorize);
	else{
	  OutPut("Select Proper Solver" << endl);
	  exit(0);
	}
#endif

#ifdef _SEQ


//   Values = NULL; KCol = NULL; Row = NULL; Symbolic=NULL, Numeric=NULL;  
//   UpdateStiffnessMat=FALSE;
//   FirstTime   = TRUE;
    if(FirstTime)
      { /** allocation and LU-decomposition forward/backward */
       DirectSolver((TSquareMatrix*)sqmatrixM[N_Levels-1], B, sol, Values, KCol, Row, Symbolic, Numeric, 0);
       FirstTime = FALSE;
      }
    else if(UpdateStiffnessMat)
      { 
       OutPut("Not yet implemented" << endl);
       exit(0);
      }
    else
     { /** only forward/backward*/
      DirectSolver((TSquareMatrix*)sqmatrixM[N_Levels-1], B, sol, Values, KCol, Row, Symbolic, Numeric, 1);   
     }
          
//         DirectSolver((TSquareMatrix*)sqmatrixM[N_Levels-1], B, sol);
#endif
	//this is set to false for direct solver factorization
        factorize = false;
      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }
     
}

double TSystemTCD3D::GetResidual(double *sol)
{
  double residual_scalar=0.0;
  
  if(SystMatAssembled)
   {
    memset(defect, 0, N_DOF*SizeOfDouble);         
    ScalarDefect(sqmatrixM[N_Levels-1], sol, B, defect, residual_scalar);
    
#ifdef _MPI 
    residual_scalar = 0.0;
    double sum =0.;
    int i,rank;
    MPI_Comm_rank(Comm, &rank);
    int *master = ParComm[N_Levels-1]->GetMaster();
    for(i=0;i<N_DOF;i++)
    {
      if(master[i]!=rank)    continue;
      residual_scalar += defect[i]*defect[i];
    }
   MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
   residual_scalar = sqrt(sum);
#endif
   }
  else
   {
    OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
    exit(4711);;   
   }
   return residual_scalar;    
}
