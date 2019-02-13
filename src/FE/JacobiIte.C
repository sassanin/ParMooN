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
// @(#)JacobiIte.C        1.24 06/27/00
//
// Class:       TJacobiIte
// Purpose:     iteration methods
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#include <ItMethod.h>
#include <JacobiIte.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/** constructor with initialization */
TJacobiIte::TJacobiIte(MatVecProc *MatVec, 
                       DefectProc *Defect, 
                       TItMethod *Prec,
                       int n_aux, int n_dof,
                       int scalar
#ifdef _MPI   
                               ,TParFECommunicator3D *ParComm
#endif      
)
  : TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
  
{
  int i;
  double *aux;
  
  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;
  oldSol=new double[N_DOF];
#ifdef _MPI
  this->ParComm=ParComm;
#endif
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
  
  N_Aux = n_aux;
  if (n_aux>0)
  {
    AuxArray = new double* [n_aux]; 
    aux = new double[n_aux*N_DOF];
    for(i=0;i<n_aux;i++)
      AuxArray[i] = aux+i*N_DOF;
  }
}

// destructor
TJacobiIte::~TJacobiIte()
{
  if (N_Aux>0)
  {
    delete AuxArray;
    delete AuxArray[0];
  }
}

// Jacobi iteration
int TJacobiIte::Iterate (TSquareMatrix **sqmat,
                         TMatrix **mat, double *sol, 
                         double *rhs)
{
  int i, j, k, l, *ARowPtr, *AKCol;
  double *AEntries;


  ARowPtr = sqmat[0]->GetRowPtr();
  AKCol = sqmat[0]->GetKCol();
  AEntries = sqmat[0]->GetEntries();
  double om = TDatabase::ParamDB->SC_SOR_OMEGA;

  for (i=0; i<N_DOF;i++)
  {
    j = ARowPtr[i];
    l = ARowPtr[i+1];
    for (k=j;k<l;k++)
      if (AKCol[k]==i)          // diagonal entry
      {
        sol[i] = om*rhs[i]/AEntries[k];
        break;
      }
  }
  return(i);
}


void TJacobiIte::Iterate_p(TSquareMatrix **sqmat, TMatrix **mat, double *sol, double *rhs
#ifdef _MPI   
                               ,TParFECommunicator3D *ParComm
#endif
                                                              )
{
  int i, j, k, l,Diagonal, *ARowPtr, *AKCol,iter=0,rank;
  double *AEntries,*d,sum,res=1.0,error=1.e-3;
  double tstrt=0,tend=0,tcomp=0,tcomm=0,comp,comm;
  
  ARowPtr = sqmat[0]->GetRowPtr();
  AKCol = sqmat[0]->GetKCol();
  AEntries = sqmat[0]->GetEntries();
  d = new double[N_DOF];
  memcpy(oldSol, sol, N_DOF*SizeOfDouble);
 
#ifdef _HYBRID
    int thid;
    omp_set_num_threads(2);
   
    #pragma omp parallel default(shared)
    {
      printf("No.of threads is %d--------\n",omp_get_num_threads());
      thid=omp_get_thread_num();
      printf("hello thid %d---\n",thid);
    }
#endif


#ifdef _HYBRID
 #pragma omp parallel default(shared) private(thid,i,j,k,l,sum,Diagonal)
 {
#endif
  while(res>error)
  {
    iter++;res=0.0;
    

#ifdef _HYBRID
   if(thid==0)  
#endif
    {
#ifdef _MPI
      tstrt = MPI_Wtime();
#endif
     Defect(sqmat,sol,rhs,d,res);
#ifdef _MPI   
     MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
     
     ParComm->CommUpdate(sol);
     tend = MPI_Wtime();
     tcomm+=(tend-tstrt);
   //ParComm->CommUpdateAlltoAllv(sol,rhs);
   // if(rank==0) printf("residue after %d iteration is %lf--------\n",iter,res);
//#else
   //  printf("residue after %d iteration is %lf--------\n",iter,res);
#endif
    }
   
    
#ifdef _MPI
    tstrt = MPI_Wtime();
#endif
    
#ifdef _HYBRID
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0; i<N_DOF;i++)
    {
      j = ARowPtr[i];
      l = ARowPtr[i+1];
      sum=rhs[i];
    
      for (k=j;k<l;k++)
      {
        sum -= AEntries[k]*sol[AKCol[k]];
        if (AKCol[k]==i)          // diagonal entry
           Diagonal=k;
      }
      sol[i]+=sum/AEntries[Diagonal];
     }
#ifdef _MPI
     tend = MPI_Wtime();
     tcomp+=(tend-tstrt);
#endif



 }
#ifdef _HYBRID
  }
#endif  

#ifdef _MPI
   tcomp=tcomp/iter; tcomm=tcomm/iter;
   MPI_Allreduce (&tcomp,&comp,1,MPI_DOUBLE,MPI_MAX,TDatabase::ParamDB->Comm);
   MPI_Allreduce (&tcomm,&comm,1,MPI_DOUBLE,MPI_MAX,TDatabase::ParamDB->Comm);
   if(rank==0) printf("no.of iterations %d comp %lf comm %lf per iter\n",iter,comp,comm);
#else
     printf("no.of iterations----%d\n",iter);
#endif
 
}


// calculate defect d=f-Ax
void TJacobiIte::Defect(TSquareMatrix **A, double *sol, double *f, double *d, double &res)
{
  ScalarDefect(A[0], sol, f, d, res);
#ifdef _MPI  
  int i, rank, *MasterOfDof, dof;
  double res_global;
  
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
  MasterOfDof = ParComm->GetMaster();
  
  for(i=0; i<N_DOF; i++)
    if(MasterOfDof[i] == rank)
      res += d[i]*d[i];
    
  MPI_Allreduce(&res, &res_global, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  res = sqrt(res_global); 
 
#endif  
} // end Defect
