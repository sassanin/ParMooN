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
// @(#)FixedPointIte.C        1.24 06/27/00
//
// Class:       TFixedPointIte
// Purpose:     iteration methods
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#ifdef _MPI
# include "mpi.h"
#endif

#include <ItMethod.h>
#include <MultiGridScaIte.h>
#include <MultiGrid2D.h>
#include <FixedPointIte.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>
#ifdef _MPI 
#include <ParFECommunicator3D.h>
#endif

#include <stdlib.h>
#include <stdio.h>

extern double tSor,tCyc,tP,tR,tD,tS,tSmoother;
double tt=0.0;
/** constructor with initialization */
TFixedPointIte::TFixedPointIte(MatVecProc *MatVec, 
                               DefectProc *Defect, 
                               TItMethod *Prec,
                               int n_aux, int n_dof,
                               int scalar                         
                               ): TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
                                    
{
  int i;
  double *aux;

  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;

  if (scalar)
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  }
  else
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
  }

  prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT;
  div_factor = TDatabase::ParamDB->SC_DIV_FACTOR;
  minit = TDatabase::ParamDB->SC_MINIT;
  restart = TDatabase::ParamDB->SC_GMRES_RESTART;

  defect = new double[N_DOF];
  
  N_Aux = n_aux;
  AuxArray = new double* [n_aux]; 
  aux = new double[n_aux*N_DOF];
  for(i=0;i<n_aux;i++)
    AuxArray[i] = aux+i*N_DOF;  
  
#ifdef _MPI   
     ParComm = NULL;
#endif     
}


#ifdef _MPI   
TFixedPointIte::TFixedPointIte(MatVecProc *MatVec, 
                               DefectProc *Defect, 
                               TItMethod *Prec,
                               int n_aux, int n_dof,
                               int scalar, 	
	#ifdef  __3D__
			       TParFECommunicator3D *parcomm
	#else		       
				TParFECommunicator2D *parcomm
	#endif 
                                ): TItMethod(MatVec, Defect, Prec, n_aux, n_dof)


                              
{
  int i;
  double *aux;

  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;

  if (scalar)
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  }
  else
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
  }

  prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT;
  div_factor = TDatabase::ParamDB->SC_DIV_FACTOR;
  minit = TDatabase::ParamDB->SC_MINIT;
  restart = TDatabase::ParamDB->SC_GMRES_RESTART;

  defect = new double[N_DOF];
  
  N_Aux = n_aux;
  AuxArray = new double* [n_aux]; 
  aux = new double[n_aux*N_DOF];
  for(i=0;i<n_aux;i++)
    AuxArray[i] = aux+i*N_DOF;  

   ParComm = parcomm;
     
}

#endif   



TFixedPointIte::~TFixedPointIte()
{
  delete defect;
  if (N_Aux>0)
  {
    delete AuxArray;
    delete AuxArray[0];
  }
}

int TFixedPointIte::Iterate (TSquareMatrix **sqmat,
                             TMatrix **mat, double *sol, 
                             double *rhs)
{
  int i=0, maxite= maxit;
  double res, res0, reslast, t1, t2;

  
#ifdef _MPI  
  int ii, rank, *MasterOfDof;  
  double  resglobal;
  
   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
   MasterOfDof = ParComm->GetMaster();
   t1 = MPI_Wtime();
#else
   t1 = GetTime();
#endif    

  matvecdefect(sqmat, mat, sol, rhs, defect);
   

#ifdef _MPI
   res=0.;
   for(ii=0; ii<N_DOF; ii++)
    if(MasterOfDof[ii] == rank){
      res += defect[ii]*defect[ii];
    }
    
   MPI_Allreduce(&res, &resglobal, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);

   res=res0=reslast = sqrt(resglobal); 
#else  
  
  res=res0=reslast=sqrt(Ddot(N_DOF,defect,defect));     
#endif  

   if (TDatabase::ParamDB->SC_VERBOSE>0     
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif      
       )   
     printf("Fixed Point Iteration 0: %e\n",  res);    
  
  // residual small enough before iteration
  if (res<=res_norm_min)
    maxite=minit;

  for (i=1; i<=maxite; i++)
  { 
#ifdef _MPI 
//   ParComm->CommUpdate(defect, defect);
    if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR != 5)
      ParComm->CommUpdate(defect);
#endif  

    // solve defect equation with the preconditioner
    // rhs and result are on last component 
// #ifdef _MPI
//     if(rank==0)
// #endif
//      printf("Fixed Point Iteration %d: %e\n", i, res); 
    
    prec->Iterate(sqmat, mat, defect, defect); 
   
    // add to sol and damp
    Daxpy(N_DOF,1.0,defect,sol);
    // compute new defect
    matvecdefect(sqmat, mat, sol, rhs, defect);
    
#ifdef _MPI
   res=0.;
   for(ii=0; ii<N_DOF; ii++)
    if(MasterOfDof[ii] == rank){
      res += defect[ii]*defect[ii];
    }
    
   MPI_Allreduce(&res, &resglobal, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);

   res = sqrt(resglobal); 
#else  
   res = sqrt(Ddot(N_DOF,defect,defect));     
#endif 
 
   if (TDatabase::ParamDB->SC_VERBOSE>0     
#ifdef _MPI  
        && rank==TDatabase::ParamDB->Par_P0
#endif        
        ) 
       printf("Fixed Point Iteration %d res %e rate: %e\n", i, res, res/reslast);    
//       OutPut("Fixed Point Iteration " << i << " " << res << " " << 
//             res/reslast << endl); 

    reslast=res;
    //  if (i<sc->maxit)
    //    {
    //      residuals[residual_cnt%AMG_CONV_RATE_BACK]=res;
    //      residual_cnt++;
    //    }  
    if (i<minit) continue;
    // relative convergence criterion fulfilled
    if (res<=res0*red_factor) break;
    // absolute convergence criterion fulfilled
    if (res<=res_norm_min) break;
    // divergence of iteration  
    if (res>div_factor*res0)
    {
      OutPut("Fixed Point iteration diverges !!!\n" << endl);
      exit(4711);
    }
  }
#ifdef _MPI
  t2 = MPI_Wtime();
#else
  t2 = GetTime();  
#endif 
  tt += t2-t1 ;
  // iteration stopped because of reaching maximal iteration number
 /* if (i>=maxite && res>res_norm_min && res>res0*red_factor )
    {
      OutPut("Fixed Point Iteration not converged !!!" << endl);
      //    iteration_cnt=i;
      //    end_residual=res;                
    }
  
  if (i>maxite) 
    i=maxite;*/
 if(1)
#ifdef _MPI
  if(rank==0)
#endif
  {
  OutPut("----------------------------------------------------------------------------------------"<<endl);  
  OutPut(" FP ITE: " << setw(4) << i <<endl);
  OutPut(" t (time taken for this FP) : " << setw(6) << t2-t1 <<endl);
  OutPut(" t (total time taken for FP) : " << setw(6) << tt <<endl);
  OutPut(" t/cyc (time taken per FP cycle) : " << setw(6) << (t2-t1)/i <<endl<<endl);

  OutPut(" tSmoother: " << setw(6) << tSmoother <<endl);
  OutPut(" tCyc (time taken for all MG cycles in FP) : " << setw(6) << tCyc <<endl<<endl);

  OutPut(" res : "  << setw(8) << res <<endl);
  OutPut(" rate: " << pow(res/res0,1.0/i) << endl<<endl);

  OutPut(" tSD (time taken for scalar defect computation): " << setw(6) << tS <<endl);
  OutPut(" tDef (time taken for defect computation): " << setw(6) << tD <<endl <<endl);
  OutPut(" tP (Prolongation): " << setw(6) << tP <<endl);
  OutPut(" tR (Restriction) : " << setw(6) << tR<<endl);
  OutPut("----------------------------------------------------------------------------------------"<<endl);
  }
  // iteration_cnt=i;
  // end_residual=res;                
  return(i);
}                        
