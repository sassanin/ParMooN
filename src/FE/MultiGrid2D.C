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
   
// =======================================================================
// @(#)MultiGrid2D.C        1.13 06/27/00
//
// Class:       TMultiGrid2D
// Purpose:     store all data for a multi grid method
//
// Author:      Gunar Matthies 02.11.1998
//
// History:     02.11.1998 start of implementation
//
// =======================================================================

#ifdef __2D__

#include <MultiGrid2D.h>
#include <FEDatabase2D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <LinAlg.h>

#include <stdlib.h>
#include <string.h>

#ifdef _MPI  
#include <ParFECommunicator2D.h>
#include <FEFunction2D.h>
#endif


/** constructor */
TMultiGrid2D::TMultiGrid2D(int n_problems, int n_parameters, 
                           double *parameters)
{
  N_Levels = 0;
  
  N_Problems = n_problems;

  N_Parameters = n_parameters;

  Parameters = parameters;
}

/** add new level as finest */
void TMultiGrid2D::AddLevel(TMGLevel2D *MGLevel)
{ 
  MultiGridLevels[N_Levels] = MGLevel;

  FESpaces[N_Levels] = MGLevel->GetFESpace();

  N_Levels++;
}

/** add new level as finest */
void TMultiGrid2D::ReplaceLevel(int i,TMGLevel2D *MGLevel)
{
  TMGLevel2D *ret;

  ret = MultiGridLevels[i];
  MultiGridLevels[i] = MGLevel;

//  MultiGridLevels[N_Levels] = MGLevel;

  FESpaces[i] = MGLevel->GetFESpace();

  if (i>=N_Levels)
    N_Levels = i+1;
}

/** restrict solution from finest grid to all coarser grids */
void TMultiGrid2D::RestrictToAllGrids()
{
  int lev, j;
  TMGLevel2D *CurrentLevel, *CoarserLevel;
  double *X1, *R1;
  double *X2, *R2, *Aux;

  for(lev=N_Levels-1;lev>0;lev--)
  {
    CurrentLevel = MultiGridLevels[lev];
    X1 = CurrentLevel->GetSolution();

    CoarserLevel = MultiGridLevels[lev-1];
    X2 = CoarserLevel->GetSolution();
    Aux = CoarserLevel->GetAuxVector(0);

    RestrictFunction(FESpaces[lev-1], FESpaces[lev], X2, 
                     X1, Aux);
  } // endfor lev
} // RestrictToAllGrids



// #ifdef _MPI  
// /***************************************************************************************
// /** Copy Values (own+Hallo) to Own FeFunct sol  
// ************************************************************************************** */
// void CopyValuesToOwnFeFunct(TFEFunction2D *C, TFEFunction2D *OwnC, double *value)
// {
// int i,j,N_owncells,N_cells,N_LocDof,*OwnBeginIndex,*BeginIndex,*OwnGN,*GN,*DOF,*OwnDOF;
// double *Solution,*OwnSolution;	
// TCollection *coll,*owncoll;
// TBaseCell *cell,*owncell;
// TFESpace2D *ScalarSpace,*OwnScalarSpace;
// 
// ScalarSpace = C->GetFESpace2D();
// OwnScalarSpace = OwnC->GetFESpace2D();
// coll = ScalarSpace->GetCollection();
// owncoll = OwnScalarSpace->GetCollection();
// 
// BeginIndex = ScalarSpace->GetBeginIndex();		
// GN = ScalarSpace->GetGlobalNumbers();         
// 
// OwnBeginIndex = OwnScalarSpace->GetBeginIndex();		
// OwnGN = OwnScalarSpace->GetGlobalNumbers();         
// 
// N_owncells = owncoll->GetN_OwnCells();
// N_cells = coll->GetN_Cells();
// 
// 
//   for(i=0;i<N_cells;i++)
//    {
//    cell = coll->GetCell(i);
//     cell->SetLocalCellNo(i);
//    }
// 
//   Solution = C->GetValues();
//   OwnSolution = OwnC->GetValues();
//   N_LocDof = BeginIndex[1]-BeginIndex[0];
//   //printf("NLOC %d\n",OwnScalarSpace->GetN_DegreesOfFreedom());
// 
//   for(i=0;i<N_owncells;i++)
//   {
//     owncell = owncoll->GetCell(i);
//     DOF = GN + BeginIndex[owncell->GetLocalCellNo()];
//     OwnDOF = OwnGN + OwnBeginIndex[i];
//     for(j=0;j<N_LocDof;j++)
//      OwnSolution[OwnDOF[j]] = value[DOF[j]];		
//     } 
// } // CopyValuesToOwnFeFunct 
// 
// 
// /** Copy Own FeFunct sol To Values (own+Hallo) */
// void CopyOwnFeFunctToValues(TFEFunction2D *C,TFEFunction2D *OwnC,double *value, double *aux, int flag)
// {
// int i,j,N_owncells,N_cells,N_LocDof,*OwnBeginIndex,*BeginIndex,*OwnGN,*GN,*DOF,*OwnDOF, N_FineDOFs;
// double *Solution,*OwnSolution;	
// TCollection *coll,*owncoll;
// TBaseCell *cell,*owncell;
// TFESpace2D *ScalarSpace,*OwnScalarSpace;
// 
// ScalarSpace = C->GetFESpace2D();
// OwnScalarSpace = OwnC->GetFESpace2D();
// coll = ScalarSpace->GetCollection();
// owncoll = OwnScalarSpace->GetCollection();
// 
// BeginIndex = ScalarSpace->GetBeginIndex();		
// GN = ScalarSpace->GetGlobalNumbers();         
// 
// OwnBeginIndex = OwnScalarSpace->GetBeginIndex();		
// OwnGN = OwnScalarSpace->GetGlobalNumbers();         
// 
// N_owncells = owncoll->GetN_OwnCells();
// N_cells = coll->GetN_Cells();
// N_FineDOFs = ScalarSpace->GetN_DegreesOfFreedom();
//  
// for(i=0;i<N_cells;i++)
// {
//   cell = coll->GetCell(i);
//   cell->SetLocalCellNo(i);
// }
// 
// Solution = C->GetValues();
// OwnSolution = OwnC->GetValues();
// N_LocDof = BeginIndex[1]-BeginIndex[0];
// //printf("NLOC %d\n",OwnScalarSpace->GetN_DegreesOfFreedom());
// 
// if(flag!=0)
//   memset(aux, 0, SizeOfDouble*N_FineDOFs);
// 
// for(i=0;i<N_owncells;i++)
// {
// 	owncell = owncoll->GetCell(i);
// 	DOF = GN + BeginIndex[owncell->GetLocalCellNo()];
// 	OwnDOF = OwnGN + OwnBeginIndex[i];
// 	for(j=0;j<N_LocDof;j++)
// 	 {
// 	   if(flag==0)
// 	   {
// 	    value[DOF[j]] = OwnSolution[OwnDOF[j]];
// 
// //             printf("OwnDof %d DOF %d GlobalCellNo %d\n", OwnDOF[j], DOF[j], owncell->GetGlobalCellNo());
// 	   }
//            else
// 	   {
// 	     if( fabs(aux[DOF[j]])<1.e-10)
// 	     {
//               value[DOF[j]] += TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR*OwnSolution[OwnDOF[j]];     
// 	      aux[DOF[j]] += 1.;
// 	     }
// 	   }
// 	 }	
// }
// 
// 
// 
// }// 
//  
// /** defect restriction from level+1 to level */
// void CellsPerDOF(TFESpace2D *FineSpace, double *aux)
// {
//   int i,j;
//   TBaseCell *cell;
//   TCollection *CoarseColl, *FineColl;
//   FE2D CoarseId, FineId;
//   TFE2D *CoarseElement, *FineElement;
//   BaseFunct2D FineBF;
//   int N_FineCells;
//   int N_FineDOFs;
//   int  *FineBeginIndex;
//   int *FineGlobalNumbers;
//   int *FineDOF;
//   int N_Fine;
//   int *DOF;
// 
//   
//   FineColl = FineSpace->GetCollection();
//   N_FineCells = FineColl->GetN_Cells();
//   FineBeginIndex = FineSpace->GetBeginIndex();
//   FineGlobalNumbers = FineSpace->GetGlobalNumbers();
//   N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();
// 
//    memset(aux, 0, S  for(i=0; i<N_DOF; i++)
//     if(MasterOfDof[i] == rank)
//       res += d[i]*d[i];
//     
//   MPI_Allreduce(&res, &res_global, 1, MPI_DOUBLE, MPI_SUM, ParComm->GetComm());
//   res = sqrt(res_global); izeOfDouble*N_FineDOFs);
// 
//   // set fine grid clipboard to -1
//   for(i=0;i<N_FineCells;i++)
//   {
//     cell = FineColl->GetCell(i);
//   
//     DOF = FineGlobalNumbers+FineBeginIndex[i];
//     FineId = FineSpace->GetFE2D(i, cell);
//     FineElement = TFEDatabase2D::GetFE2D(FineId);
//     FineBF = FineElement->GetBaseFunct2D_ID();
//     N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();
//     for(j=0;j<N_Fine;j++)
//       aux[DOF[j]] += 1;
//   }
// }
// 
// #endif  



/** one cycle on level i */
void TMultiGrid2D::Cycle(int i, double &res)
{
  double s;
  
  TMGLevel2D *CurrentLevel, *CoarserLevel;
  double *CurrentSol, *CoarserSol, *CoarserRhs;
  double *CurrentRhs, *CurrentDefect, *CurrentAux;
  double *CurrentAux2, *OldSol, *OldDefect, res_global;
  double oldres,reduction, alpha;
  int j, N_DOF, maxit, it, slc,gam;

  CurrentLevel = MultiGridLevels[i];
  CurrentDefect = CurrentLevel->GetAuxVector(0);
  CurrentAux = CurrentLevel->GetAuxVector(1);
  slc =0;
  if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)||
      (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR))
    slc = 1;
  if (slc)
  {
    OldSol = CurrentLevel->GetAuxVector(2);
    OldDefect = CurrentLevel->GetAuxVector(3);
  }
  CurrentAux2 = CurrentDefect;
  CurrentSol = CurrentLevel->GetSolution();
  CurrentRhs = CurrentLevel->GetRhs();
  N_DOF = CurrentLevel->GetN_DOF(); 
  
#ifdef _MPI  
  TParFECommunicator2D *ParComm, *CoarseParComm; 
  int rank;
 
//   ParComm = CurrentLevel->GetParComm();  

//   MPI_Comm_rank(ParComm->GetComm(), &rank); 
#endif  
      
  
  if(i==0)
  {
    // coarse grid
    res = 1;
    maxit =  TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
    it = 0;
           
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
           
    if (TDatabase::ParamDB->SC_VERBOSE>=2      
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif        
       )
    {
//        OutPut("rank: " <<rank <<" level 0 res before smoothing: " << res << endl);
            printf("level 0 res before smoothing: %e\n",  res);    
//          printf("Multigrid rank %d \n", rank);
    }
    
    
    reduction = TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR*res;
    while ((res>reduction)&&(it<maxit))
    {
      switch(TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR)
      {
 
        case 1: // Jacobi
          CurrentLevel->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
          break;
        case 2: // SOR
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
          break;
        case 3: // SSOR
          CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
          break;
        case 4: // ILU
          CurrentLevel->ILU(CurrentSol, CurrentRhs, CurrentDefect,
            N_Parameters, Parameters);
          break;
        case 17: // solution with Gaussian elimination
          CurrentLevel->SolveExact(CurrentSol, CurrentRhs);
          break;
        default:
           OutPut("Coarse smoother not implemented !! Use coarse smoother 3" << endl);
           CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                              N_Parameters, Parameters);
      } // endswitch SC_COARSE_SMOOTHER
      
#ifdef _MPI  
      // communicate the values (sol & rhs) to the slave DOFs from master DOF
//       ParComm->CommUpdate(CurrentSol);
#endif   
      
      CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
     
      it++;
      if (TDatabase::ParamDB->SC_VERBOSE>=2      
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif        
       )
      {
//          OutPut("No. Itr " << it << "  level 0 res after smoothing: " << res << endl);
        printf("No. Itr %d level 0 res after smoothing: %e\n", it, res);    
      }   
     
    }
/*    
printf("Multigrid rank %d \n", rank);
MPI_Finalize();  
exit(0);    */
    
  }
  else
  {
    slc =0;

    if (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
      slc = 1;
    else if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR)
             &&(i==N_Levels-1))
      slc = 1;

    CoarserLevel = MultiGridLevels[i-1];
    CoarserSol = CoarserLevel->GetSolution();
    CoarserRhs = CoarserLevel->GetRhs();
#ifdef _MPI  
//   CoarseParComm = CoarserLevel->GetParComm();      
#endif       
    
    // smoothing
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);

    if (TDatabase::ParamDB->SC_VERBOSE>=2      
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif        
       )
    {
       OutPut("level " << i << " ");
       OutPut("res before presmoothing : " << oldres << endl);
    }
    if (slc)
    {
      memcpy(OldSol, CurrentSol, N_DOF*SizeOfDouble);
      memcpy(OldDefect, CurrentDefect, N_DOF*SizeOfDouble);
    }

    switch(TDatabase::ParamDB->SC_SMOOTHER_SCALAR)
    {
// #ifdef _MPI  
//       ParComm->CommUpdate(CurrentSol, CurrentRhs);
// #endif
      
      case 1: // Jacobi
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
//          ParComm->CommUpdate(CurrentSol);
#endif
         }
  
        break;
      case 2: // SOR
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
         {
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
// #ifdef _MPI  
//          ParComm->CommUpdate(CurrentSol);
// #endif
         }

        break;
      case 3: // SSOR
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
          CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
        break;
      case 4: // ILU
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
        {
          CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);
 
          CurrentLevel->ILU(CurrentSol, CurrentRhs, CurrentDefect,
                N_Parameters, Parameters);
        }
        break;
      default:
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
    } // endswitch SC_SMOOTHER

    // calculate defect
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);

    
    if (TDatabase::ParamDB->SC_VERBOSE>=2      
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif        
       )
      {
        OutPut("level " << i << " ");
        OutPut("res after presmoothing: " << oldres << endl);
      }
      
/* int ii;
for(ii=0;ii<N_DOF;ii++)
	CurrentDefect[ii] = 2.0;   */   
 
#ifdef _MPI  
//         ParComm->CommUpdate(CurrentDefect);
//         CellsPerDOF(CurrentLevel->GetFESpace(),CurrentAux);
//         CopyValuesToOwnFeFunct(CurrentLevel->GetFEFunction(),  CurrentLevel->GetOwnFEFunction(), CurrentAux);
//         memcpy(CurrentAux,CurrentLevel->GetOwnFEFunction()->GetValues(),(CurrentLevel->GetOwnFESpace()->GetN_DegreesOfFreedom())*sizeof(double));
// 
//         CopyValuesToOwnFeFunct(CurrentLevel->GetFEFunction(),  CurrentLevel->GetOwnFEFunction(), CurrentDefect);

        /** first restrict defect to coarse sol, then copy coarse own sol to all (hallo) coarse rhs **/ 
/*         DefectRestriction(CoarserLevel->GetOwnFESpace(), CurrentLevel->GetOwnFESpace(),
                           CoarserLevel->GetOwnSolution(), CurrentLevel->GetOwnSolution(),
                           CurrentAux);

         memset(CoarserRhs, 0, SizeOfDouble*CoarserLevel->GetN_DOF());
         CopyOwnFeFunctToValues(CoarserLevel->GetFEFunction(),  CoarserLevel->GetOwnFEFunction(), CoarserRhs, CurrentAux, 0); 
         CoarseParComm->CommUpdateReduce(CoarserSol,CoarserRhs);  */     
#else      
    // restrict defect
    DefectRestriction(FESpaces[i-1], FESpaces[i], CoarserRhs, CurrentDefect, CurrentAux);
    
#endif    
    //OutPut("a"<<endl);
    CoarserLevel->CorrectDefect(CoarserRhs);
    CoarserLevel->Reset(CoarserSol);
    // coarse grid correction, apply mg recursively
     
    
    
//      CoarserLevel->Defect(CoarserSol, CoarserRhs, CoarserLevel->GetAuxVector(0), oldres);
//             OutPut("res Before Coarse Solve : " << oldres << endl);
//      
     
     
    for(j=0;j<mg_recursions[i];j++)
      Cycle(i-1, res);
    if (TDatabase::ParamDB->SC_MG_CYCLE_SCALAR<1) mg_recursions[i] = 1;              /* F--cycle */


    // prolongate correction
#ifdef _MPI  
//       // copy sol from own+Hallo fespace to own fespace
//       CopyValuesToOwnFeFunct(CoarserLevel->GetFEFunction(),  CoarserLevel->GetOwnFEFunction(),  CoarserSol);
      
      /** first restrict to coarse sol, then copy coarse sol to all (hallo) coarse rhs **/
//       Prolongate(CoarserLevel->GetOwnFESpace(), CurrentLevel->GetOwnFESpace(),
//                        CoarserLevel->GetOwnSolution(), CurrentLevel->GetOwnSolution(),
//                        CurrentLevel->GetAuxVector(1));
// 
//       CopyOwnFeFunctToValues(CurrentLevel->GetFEFunction(),  CurrentLevel->GetOwnFEFunction(), CurrentSol, CurrentLevel->GetAuxVector(1), 1);    
//       
//       ParComm->CommUpdate(CurrentSol);
#else  
    Prolongate(FESpaces[i-1], FESpaces[i], CoarserSol, CurrentAux2, CurrentAux);
    
    CurrentLevel->CorrectNodes(CurrentAux2);
    CurrentLevel->Update(CurrentSol, CurrentAux2);    
#endif     

//   OutPut("MGlelvel Coarse norm " <<  sqrt(Ddot(FESpaces[i-1]->GetN_DegreesOfFreedom(), 
//                                      CoarserSol,CoarserSol))  << endl);
//      
//    OutPut("MGlelvel Fine norm " <<  sqrt(Ddot(FESpaces[i]->GetN_DegreesOfFreedom() ,CurrentAux2,CurrentAux2))  << endl);

    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);

    if (TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif        
       )
      {
        OutPut("level " << i << " ");
//         OutPut("res before postsmoothing: " << oldres << endl);
        printf("level %d, res before postsmoothing: %e \n",  i, oldres);	
      }

    // smoothing
    switch(TDatabase::ParamDB->SC_SMOOTHER_SCALAR)
    {
      case 1: // Jacobi
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
         {
          CurrentLevel->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
//           ParComm->CommUpdate(CurrentSol);
#endif     
	 }
        break;
      case 2: // SOR
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
// #ifdef _MPI  
//          ParComm->CommUpdate(CurrentSol);
// #endif
         }
  
        break;
      case 3: // SSOR
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
          CurrentLevel->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
        break;
      case 4: // ILU
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
        {
          CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);

          CurrentLevel->ILU(CurrentSol, CurrentRhs, CurrentDefect,
                N_Parameters, Parameters);
        }
        break;
      default:
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          CurrentLevel->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
#ifdef _MPI  
//          ParComm->CommUpdate(CurrentSol);
#endif
         }
	  
    } // endswitch SC_SMOOTHER

   
    
    // calculate defect

    if (slc)
    {
      alpha = CurrentLevel->StepLengthControl(CurrentSol, OldSol,
                          OldDefect,                  
                          N_Parameters,Parameters);       
      
      for (j=0;j<N_DOF;j++)
        CurrentSol[j] = OldSol[j] + alpha *( CurrentSol[j]-OldSol[j]);
    }

    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
     
    if (TDatabase::ParamDB->SC_VERBOSE>=2      
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif        
       )
      {
//         OutPut("level " << i << " ");
//         OutPut("res after postsmoothing: " << res << endl);
        printf("level %d, res after postsmoothing: %e \n",  i, res);	
      }
      
/*//          if(rank==TDatabase::ParamDB->Par_P0   )    
// printf("mg_recursions %d, Multigrid rank %d \n",  i, rank);
MPI_Finalize();  
exit(0); */     
      
  }
}

void TMultiGrid2D::SetDirichletNodes(int i)
{
  int HangingNodeBound, N_Dirichlet;
  TMGLevel2D *CurrentLevel;
  double *X, *R;

  if(i>=N_Levels) return;

  CurrentLevel = MultiGridLevels[i];

  X = CurrentLevel->GetSolution();
  R = CurrentLevel->GetRhs();

  HangingNodeBound = CurrentLevel->GetHangingNodeBound();
  N_Dirichlet = CurrentLevel->GetN_Dirichlet();

  memcpy(X+HangingNodeBound, R+HangingNodeBound,
         SizeOfDouble*N_Dirichlet);
}

/** set recursion for multigrid */ 
void TMultiGrid2D::SetRecursion(int levels)
{
  int gam = TDatabase::ParamDB->SC_MG_CYCLE_SCALAR,k;

  // coarsest grid 
  mg_recursions[1] = 1;
  if (gam>0)
    for (k=2;k<=levels;k++)        
      mg_recursions[k] = gam;
  else                /* F -- cycle */
    for (k=2;k<=levels;k++)        
      mg_recursions[k] = 2;
}

#endif // #ifdef __2D__
