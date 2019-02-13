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
// %W% %G%
//
// Class:       TMultiGrid3D
// Purpose:     store all data for a multi grid method in 3d
//
// Author:      Gunar Matthies 26.06.2000
//
// History:     26.06.2000 start of implementation
//
// =======================================================================

#include <MultiGrid3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <LinAlg.h>

#include <stdlib.h>
#include <string.h>
#ifdef _MPI  
#include <ParFECommunicator3D.h>
#include <FEFunction3D.h>
#endif

#define PreSmooth -1
#define CoarseSmooth 0
#define PostSmooth 1

double tSmoother = 0.0;
/** constructor */
TMultiGrid3D::TMultiGrid3D(int n_problems, int n_parameters, 
                       double *parameters)
{
  N_Levels = 0;
  
  N_Problems = n_problems;

  N_Parameters = n_parameters;

  Parameters = parameters;
}

/** add new level as finest */
void TMultiGrid3D::AddLevel(TMGLevel3D *MGLevel)
{ 
  MultiGridLevels[N_Levels] = MGLevel;

  FESpaces[N_Levels] = MGLevel->GetFESpace();

  N_Levels++;
}
/** add new level as finest */
void TMultiGrid3D::ReplaceLevel(int i,TMGLevel3D *MGLevel)
{ 
  TMGLevel3D *ret;

  ret = MultiGridLevels[i];
  MultiGridLevels[i] = MGLevel;

//  MultiGridLevels[N_Levels] = MGLevel;

  FESpaces[i] = MGLevel->GetFESpace();

  if (i>=N_Levels)
    N_Levels = i+1;
}

/** restrict solution from finest grid to all coarser grids */
void TMultiGrid3D::RestrictToAllGrids()
{
  int lev, j;
  TMGLevel3D *CurrentLevel, *CoarserLevel;
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

void TMultiGrid3D::Smooth(int smoother_type, TMGLevel3D *Level, 
#ifdef _MPI
			  TParFECommunicator3D *ParComm, 
#endif
			  double &oldres)
{
  int i,j,k;
  double *CurrentSol, *CurrentRhs, *CurrentDefect, *CurrentAux;
#ifdef _HYBRID
  bool firstTime, LastTime;
#endif
  
  CurrentSol    = Level->GetSolution();
  CurrentRhs    = Level->GetRhs();
  CurrentDefect = Level->GetAuxVector(0);
  CurrentAux    = Level->GetAuxVector(1);
  
  if(smoother_type == PreSmooth)			//presmoothing
  {
    switch(TDatabase::ParamDB->SC_SMOOTHER_SCALAR)
    {
      case 1: // Jacobi
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          Level->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol);
#endif
	}
        break;
      case 2: // SOR
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol);
#endif
	}
        break;
      case 3: // SSOR
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          Level->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
         ParComm->CommUpdate(CurrentSol);
#endif
	}
        break;
      case 4: // ILU
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
        {
          Level->Defect(CurrentSol, CurrentRhs, CurrentDefect, 
                oldres);
          Level->ILU(CurrentSol, CurrentRhs, CurrentDefect,
                N_Parameters, Parameters);
        }
        break;
#ifdef _MPI
  #ifdef _HYBRID
	case 5: //SOR_Reorder
	printf("Not Working\n");
	MPI_Finalize();
	exit(0);
        break;
  #else	
	case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #endif
#endif

#ifdef _MPI
  #ifdef _HYBRID
	case 6: //SOR Reorder and color
	  Level->SOR_Re_Color(CurrentSol, CurrentRhs, CurrentAux, N_Parameters, Parameters, PreSmooth);
        break;
  #else	
	case 6: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          Level->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #endif
#ifdef _MPI
	case 7: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          Level->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
#endif
#endif	
      default:
        for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
          Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
    } // endswitch SC_SMOOTHER
  }
  else if(smoother_type == CoarseSmooth)
  {
#ifdef _MPI  
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
#endif  
    double res;
    int it = 0;
    int maxit =  TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
    
    Level->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
     if(TDatabase::ParamDB->SC_VERBOSE>=2 
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif  
     )
     {
      OutPut("residual before on coarse "<<res << endl);
     }
     
    double reduction = TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR*res;
  
    while ((res>reduction)&&(it<maxit))
    {  
      switch(TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR)
            {
                 case 1: // Jacobi
                         Level->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
                                                       N_Parameters, Parameters);
#ifdef _MPI  
                         ParComm->CommUpdate(CurrentSol);
#endif
                         break;
		       
                 case 2: // SOR
                         Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
                                                     N_Parameters, Parameters);
#ifdef _MPI  
                         ParComm->CommUpdate(CurrentSol);
#endif
                         break;
		       
                 case 3: // SSOR
                         Level->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                                                    N_Parameters, Parameters);
#ifdef _MPI  
                         ParComm->CommUpdate(CurrentSol);
#endif
                         break;
        
	         case 4: // ILU
                         Level->ILU(CurrentSol, CurrentRhs, CurrentDefect,
                                                        N_Parameters, Parameters);
                         break;
		       
                 case 17: // solution with Gaussian elimination
                          Level->SolveExact(CurrentSol, CurrentRhs);
                          break;
#ifdef _MPI
  #ifdef _HYBRID
	         case 5: //SOR_Reorder
	                 printf("Not Working\n");
	                 MPI_Finalize();
	                 exit(0);
                         break;
  #else
                 case 5: //SOR_Reorder
	                 for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	                 {
                               Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
                                                            N_Parameters, Parameters);
	                 }
                         break;
  #endif
#endif

#ifdef _MPI
  #ifdef _HYBRID
	  case 6: //SOR_Reorder
	    Level->SOR_Re_Color(CurrentSol, CurrentRhs, CurrentAux, N_Parameters, Parameters, CoarseSmooth);
          break;
  #else	
	  case 6: //SOR_Reorder
            Level->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                  N_Parameters, Parameters);
          break;
  #endif
		
	   case 7: //SOR_Reorder
            Level->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                  N_Parameters, Parameters);
          break;
#endif	  
               default:
                       OutPut("Coarse smoother not implemented !! Use coarse smoother 3" << endl);
                       Level->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                                                    N_Parameters, Parameters);
#ifdef _MPI  
                       ParComm->CommUpdate(CurrentSol);
#endif
          } // endswitch SC_COARSE_SMOOTHER_SCALAR
          
       Level->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
       it++;
       if(TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif
        )
         OutPut("itr no. :: "<<it-1<<"        res on coarse: " << res << endl);   
    }//endwhile
    oldres = res;
  }
  else if(smoother_type == PostSmooth)
  {
    switch(TDatabase::ParamDB->SC_SMOOTHER_SCALAR)
    {
      case 1: // Jacobi
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          Level->Jacobi(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol);
#endif 
	}
        break;
      case 2: // SOR
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol);
#endif 
	}
        break;
      case 3: // SSOR
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          Level->SSOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol);
#endif 
	}
        break;
      case 4: // ILU
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
        {
          Level->Defect(CurrentSol, CurrentRhs, CurrentDefect, 
                oldres);
          Level->ILU(CurrentSol, CurrentRhs, CurrentDefect,
                N_Parameters, Parameters);
        }
        break;
#ifdef _MPI
  #ifdef _HYBRID
	case 5: //SOR_Reorder
	printf("Not Working\n");
	MPI_Finalize();
	exit(0);
        break;
  #else	
	case 5: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;j++)
	{
          Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #endif
#endif

#ifdef _MPI
  #ifdef _HYBRID
	case 6: //SOR Reorder and color

	  Level->SOR_Re_Color(CurrentSol, CurrentRhs, CurrentAux, N_Parameters, Parameters, PostSmooth);

        break;
  #else	
	case 6: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          Level->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
  #endif	

	case 7: //SOR_Reorder
	for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          Level->SOR_Re(CurrentSol, CurrentRhs, CurrentAux,
                N_Parameters, Parameters);
	}
        break;
#endif	
      default:
        for(j=0;j<TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;j++)
	{
          Level->SOR(CurrentSol, CurrentRhs, CurrentAux,
            N_Parameters, Parameters);
#ifdef _MPI  
          // communicate the values (sol & rhs) to the slave DOFs from master DOF
          ParComm->CommUpdate(CurrentSol);
#endif 
	}
    } // endswitch SC_SMOOTHER_SCALAR
  }
  else{
    printf("Wrong smoother_type \n");
#ifdef _MPI
    MPI_Finalize();
#endif
    exit(0);
    }
}

/** one cycle on level i */
void TMultiGrid3D::Cycle(int i, double &res)
{
  double s,t1,t2;
  
  TMGLevel3D *CurrentLevel, *CoarserLevel;
  double *CurrentSol, *CoarserSol, *CoarserRhs;
  double *CurrentRhs, *CurrentDefect, *CurrentAux;
  double *CurrentAux2, *OldSol, *OldDefect;
  double oldres,reduction, alpha;
  int j, N_DOF, maxit, it, slc,gam;
  double initres, normsol, firstres;

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
  TParFECommunicator3D *ParComm, *CoarseParComm; 
  int rank;
 
   ParComm = CurrentLevel->GetParComm();  

   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
#endif  
 // OutPut("Norm of B rhs in cycle " <<  sqrt(Ddot(N_DOF,CurrentRhs,CurrentRhs)) <<"i is:"<<i <<endl); 
  
  if(i==0)
  {
    // coarse grid
//     cout << "coarse grid" << endl;
//     res = 1;
//     maxit =  TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
//     it = 0;
//     CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
    
//    if(TDatabase::ParamDB->SC_VERBOSE>=2 
// #ifdef _MPI  
//         && rank==TDatabase::ParamDB->Par_P0
// #endif  
//      )
//      {
//       OutPut("residual before on coarse "<<res << endl);
//      }
     
//     reduction = TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR*res;
//     while ((res>reduction)&&(it<maxit))
//     {
#ifdef _MPI
     t1 = MPI_Wtime();
     Smooth(CoarseSmooth, CurrentLevel, ParComm, res); 
     t2 = MPI_Wtime();
#else
     t1 = GetTime();
     Smooth(CoarseSmooth, CurrentLevel, res);
     t2 = GetTime();
#endif
     tSmoother += t2-t1 ;
     
//       CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, res);
//       it++;
//       if(TDatabase::ParamDB->SC_VERBOSE>=2
// #ifdef _MPI  
//         && rank==TDatabase::ParamDB->Par_P0
// #endif
//         )
//          OutPut("itr no. :: "<<it-1<<"        res on coarse: " << res << endl);
//     }//end while
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
  CoarseParComm = CoarserLevel->GetParComm();   
#endif  
    // smoothing
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);  
    firstres = initres = oldres;  
    normsol = sqrt(Ddot(N_DOF, CurrentSol, CurrentSol));
  
    if(TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif   
      )
      {
       OutPut("level " << i << " ");
       OutPut("res before presmoothing: " << oldres << endl);
      }

    if (slc)
    {
      memcpy(OldSol, CurrentSol, N_DOF*SizeOfDouble);
      memcpy(OldDefect, CurrentDefect, N_DOF*SizeOfDouble);
    }

#ifdef _MPI
     t1 = MPI_Wtime();
     Smooth(PreSmooth, CurrentLevel, ParComm, oldres); 
     t2 = MPI_Wtime();
#else
     t1 = GetTime();
     Smooth(PreSmooth, CurrentLevel, oldres);
     t2 = GetTime();
#endif
     tSmoother += t2-t1 ;
    // calculate defect
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);
    
    if (TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif 
       )
      {
         OutPut("normsol: " << normsol << " oldres: " << oldres << endl);
        OutPut("level " << i << " ");
        OutPut("res after presmoothing: " << oldres << endl);
         OutPut("Smoothing (" << i << "): " << oldres/normsol << endl);
      }
    // restrict defect
//     exit(0);
#ifdef _MPI  
        ParComm->CommUpdate(CurrentDefect);
  
	memcpy(CurrentLevel->GetTemp_arr(),CurrentDefect,CurrentLevel->GetN_DOF());
    
	DefectRestriction(CoarserLevel->GetFESpace(), CurrentLevel->GetFESpace(),
                           CoarserRhs, CurrentDefect,
                           CurrentAux);
// 	OutPut("restriction "<<endl);
// 	exit(0);

	  CoarseParComm->CommUpdateReduce(CoarserRhs);	 

// 	OutPut("2.restriction "<<endl);
	//exit(0);
#else
        DefectRestriction(FESpaces[i-1], FESpaces[i],CoarserRhs, CurrentDefect, CurrentAux);
#endif

    CoarserLevel->CorrectDefect(CoarserRhs);  //non-active part set to 0
    CoarserLevel->Reset(CoarserSol);          //all set to 0
     
    // coarse grid correction
    // coarse grid correction, apply mg recursively*/
    for(j=0;j<mg_recursions[i];j++)
       Cycle(i-1, res);
    if (TDatabase::ParamDB->SC_MG_CYCLE_SCALAR<1) mg_recursions[i] = 1;              // F--cycle 

    // prolongate correction
#ifdef _MPI  
      Prolongate(CoarserLevel->GetFESpace(), CurrentLevel->GetFESpace(),
                       CoarserSol, CurrentLevel->GetTemp_arr(),
                       CurrentLevel->GetAuxVector(1));


	ParComm->CommUpdateReduce(CurrentLevel->GetTemp_arr());

      
      CurrentLevel->Update(CurrentSol, CurrentLevel->GetTemp_arr());
      
      ParComm->CommUpdate(CurrentSol);
#else 
    Prolongate(FESpaces[i-1], FESpaces[i], 
                   CoarserSol, CurrentAux2, CurrentAux);

    CurrentLevel->CorrectNodes(CurrentAux2);

    CurrentLevel->Update(CurrentSol, CurrentAux2);
#endif  
    
    CurrentLevel->Defect(CurrentSol, CurrentRhs, CurrentDefect, oldres);
  
    initres = oldres;
    normsol = sqrt(Ddot(N_DOF, CurrentSol, CurrentSol));
    
    if (TDatabase::ParamDB->SC_VERBOSE>=2
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif  
       )
      {
        OutPut("level " << i << " ");
        OutPut("res before postsmoothing: " << oldres << endl);
      }
//        exit(0);
    // smoothing
#ifdef _MPI
     t1 = MPI_Wtime();
     Smooth(PostSmooth, CurrentLevel, ParComm, oldres); 
     t2 = MPI_Wtime();
#else
     t1 = GetTime();
     Smooth(PostSmooth, CurrentLevel, oldres);
     t2 = GetTime();
#endif
     tSmoother += t2-t1 ;
     
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
        OutPut("level " << i << " ");
        OutPut("res after postsmoothing: " << res);
        OutPut(" rate: " << res/firstres << endl);
        // OutPut("Smoothing2 (" << i << "): " << initres/normsol << endl);
      }
  }
}

void TMultiGrid3D::SetDirichletNodes(int i)
{
  int HangingNodeBound, N_Dirichlet;
  TMGLevel3D *CurrentLevel;
  double *X, *R;

  if(i>=N_Levels) return;

  CurrentLevel = MultiGridLevels[i];

  X = CurrentLevel->GetSolution();
  R = CurrentLevel->GetRhs();

  HangingNodeBound = CurrentLevel->GetHangingNodeBound();
  N_Dirichlet = CurrentLevel->GetN_Dirichlet();

  memcpy(X+HangingNodeBound, R+HangingNodeBound, SizeOfDouble*N_Dirichlet);
}

/** set recursion for MultiGrid3D */ 
void TMultiGrid3D::SetRecursion(int levels)
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
