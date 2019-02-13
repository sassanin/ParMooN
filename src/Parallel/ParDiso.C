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
// @(#)ParDiso.C
//
// Class:      TParDiso
// Purpose:    Solve equation system by ParDiso routines
//
// Author:     Sashikumaar Ganesan & Abdus Shamim (01.05.15)
//
// History:    Start of implementation 01.05.15 (Sashikumaar Ganesan & Abdus Shamim) 
//
// =======================================================================

#ifdef _OMPONLY
#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Database.h>

#include <ParDiso.h>

#define iparam(x) iparam[x-1]

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);

static void pardiso_error(int error , int *iparam)
{
  if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
	else
	   printf("error = %d \n",error);
    }
    else
        printf("[PARDISO]: License check was successful ... \n");

    exit(0);
    
}

/** constructor */
TParDiso::TParDiso(int neqns, int nnz, int* rowptr, int* kcol)
{
  int i, j, k;
  double t1;
  
  t1 = GetTime();
  
  int numThreads = TDatabase::ParamDB->OMPNUMTHREADS;
  omp_set_num_threads(numThreads);

  ierror = 0;

  nrhs          = 1;
  Nmax_system   = 1;
  matrix_number = 1;
  matrix_type   = 11;			//real and unsymmetric
  N_Eqns        = neqns;
  N_Nz          = nnz;
  msglvl        = 0;
  //csr format
  RowPtr        = new int[N_Eqns+1];
  KCol          = new int[N_Nz];
  
  for(i=0;i<N_Eqns+1;i++)
    RowPtr[i] = rowptr[i] + 1;
  
  for(i=0;i<N_Nz;i++)
    KCol[i] = kcol[i] + 1;
  
  //0 - sparse direct solver
  //1 - multi recursive iterative solver
  Solver        = 0;
  /* Dummy pointer for user permutation. */
  perm_user     = 0;
    
  memset(dparam , 0, sizeof(dparam));
  memset(pt     , 0, sizeof(pt));

  pardisoinit (pt,  &matrix_type, &Solver, iparam, dparam, &ierror);  

  if (ierror)
      pardiso_error(ierror , iparam);

    iparam(1) = 1;		             //0-Set all entries to their default values except IPARM(3)
    iparam(2) = 2;			     //use metis
    iparam(3) = numThreads;                  /* maximal number of OpenMP tasks */
    iparam(4) = 0;
    iparam(5) = 0;
    iparam(6) = 0;			     /* 1-sol is written to rhs**/                                            
    iparam(8) = 0;			     /* no iterative refinement */
    iparam(10) = 13;
    iparam(11) = 1;
    iparam(12) = 0;
    iparam(13) = 1;
    iparam(18) = -1;
    iparam(19) = 0;
    iparam(24) = 1;
    iparam(25) = 1;
    iparam(26) = 0;
    iparam(28) = 1;
    iparam(29) = 0;
    iparam(30) = 0;
    iparam(31) = 0;
    iparam(32) = 0;
    iparam(33) = 0;
    iparam(34) = 0;
    iparam(36) = 0;
    iparam(37) = 0;
    iparam(38) = 0;
    iparam(51) = 0;			     /* 0-Use OpenMP-threaded solver 
                                                1-Use Mixed OpenMP-MPI solver*/
    iparam(52) = 1;			     /* 1 - For OpenMP-threaded solver
                                                p - Use p compute nodes*/
                                                
                                                
//     dparam(4) = 10;                          /* Maximum Number of Grid Levels */
    
    printf("iparam(3)::%d\tNumber of OpenMP tasks per host\n", iparam(3));
    printf("iparam(52)::%d\tNumber of hosts\n", iparam(52));

    OutPut("Time taken for initialization & Analysis: "<< GetTime()- t1 << endl);
}


void TParDiso::FactorizeAndSolve(double *Mat, double *rhs, double *sol,bool Factorize)
{
    if(0)//verification
    {
    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
	
	pardiso_chkmatrix(&matrix_type, &N_Eqns, Mat, RowPtr, KCol, &ierror);
	if (ierror != 0) {
	    printf("\nERROR in consistency of matrix: %d", ierror);
	    exit(1);
	}
	
    /* -------------------------------------------------------------------- */
    /* ..  pardiso_chkvec(...)                                              */
    /*     Checks the given vectors for infinite and NaN values             */
    /*     Input parameters (see PARDISO user manual for a description):    */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */

	pardiso_chkvec(&N_Eqns, &nrhs, rhs, &ierror);
	if (ierror != 0) {
	    printf("\nERROR  in right hand side: %d", ierror);
	    exit(1);
	}
	
    /* -------------------------------------------------------------------- */
    /* ..  pardiso_chkvec(...)                                              */
    /*     Checks the given vectors for infinite and NaN values             */
    /*     Input parameters (see PARDISO user manual for a description):    */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */

	pardiso_chkvec(&N_Eqns, &nrhs, sol, &ierror);
	if (ierror != 0) {
	    printf("\nERROR  in sol vector: %d", ierror);
	    exit(1);
	}    
	
    /* -------------------------------------------------------------------- */
    /* .. pardiso_printstats(...)                                           */
    /*    prints information on the matrix to STDOUT.                       */
    /*    Use this functionality only for debugging purposes                */
    /* -------------------------------------------------------------------- */

	pardiso_printstats(&matrix_type, &N_Eqns, Mat, RowPtr, KCol, &nrhs, rhs, &ierror);
	if (ierror != 0) {
	    printf("\nERROR right hand side: %d", ierror);
	    exit(1);
	}    
	printf("\n``````````````````````````````Matrix check done: \n");
	
	
    }//verification
    
    double tf,temp;
    double t1 = GetTime();
    int i,j,k;
    
    if(Factorize)
    {
      phase = 11;
      pardiso(pt , &Nmax_system , &matrix_number , &matrix_type , &phase  , &N_Eqns,
	      Mat, RowPtr       , KCol           , &perm_user   , &nrhs   , iparam ,
	      &msglvl           , &ddum            , &ddum          , &ierror , dparam);
    
      if (ierror)
        pardiso_error(ierror , iparam);
            
      phase = 22;
      pardiso(pt , &Nmax_system , &matrix_number , &matrix_type , &phase  , &N_Eqns,
	      Mat, RowPtr       , KCol           , &perm_user   , &nrhs   , iparam ,
	      &msglvl           , &ddum            , &ddum          , &ierror , dparam);
    
      if (ierror)
        pardiso_error(ierror , iparam);
      
      tf = GetTime()- t1;
      OutPut("Time taken for factorization: "<< tf << endl);
      t1 = GetTime();
      
    }
   
    /*
    * Start forward/backward substitution of the PARDISO solver.
    */
    phase = 33;
    pardiso(pt , &Nmax_system , &matrix_number , &matrix_type , &phase  , &N_Eqns,
	       Mat, RowPtr       , KCol           , &perm_user   , &nrhs   , iparam ,
	       &msglvl           , rhs            , sol          , &ierror , dparam);
    
    if (ierror)
      pardiso_error(ierror , iparam);
    
    OutPut("Time taken for solving: "<< GetTime()- t1 << endl);
    OutPut("Time taken for factorization & solving: "<< tf + GetTime()- t1 << endl);
    
//     printf("\nSolve completed ... ");
//     printf("\nThe solution of the system is: ");
//     for (i = 0; i < N_Eqns; i++) {
//         printf("\n x [%d] = % lf", i, sol[i] );
//     }
    
    rhsptr = rhs;
    solptr = sol;
}

void TParDiso::Solve(double *Mat, double *rhs, double *sol)
{
  FactorizeAndSolve(Mat, rhs, sol,false);
}

void TParDiso::Clean(double *Mat)
{
  delete [] RowPtr;
  delete [] KCol;
  
    /* Free PARDISO internal memory*/
    phase = 0;
    pardiso(pt , &Nmax_system , &matrix_number, &matrix_type , &phase, &N_Eqns,
	     Mat, RowPtr       , KCol          , &perm_user   , &nrhs , iparam ,
	     &msglvl , rhsptr     , solptr           , &ierror      , dparam);
    
    if (ierror)
      pardiso_error(ierror , iparam);

}

#endif
