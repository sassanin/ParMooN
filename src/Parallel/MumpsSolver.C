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
// @(#) MumpsSolver.C
//
// Purpose:    Solve a system by MPI based MUMPS parallel direct solvers
//
// Author:     Sashikumaar Ganesan (30.09.09)
//
// History:    Start of implementation 30.09.09 (Sashikumaar Ganesan)
//             Major changes 21.10.2010 (Sashi)
// =======================================================================

#if defined(_MPI) || defined(_SMPI)
#  include "mpi.h"

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <Database.h>
 #include <MumpsSolver.h>
#ifdef __2D__
 #include <ParFECommunicator2D.h>
#endif
#ifdef __3D__
 #include <ParFECommunicator3D.h>
#endif
 
//  #include <ParVector.h>
//  #include <ParVectorNSE3D.h>

extern "C"
 {
  #include "dmumps_c.h"
 }

#  define JOB_INIT -1
#  define JOB_ANALYSIS 1
#  define JOB_FACTORIZ 2
#  define JOB_SOLVE 3
#  define JOB_FACTORIZANDSOLVE 5
#  define JOB_END -2
#  define USE_COMM_WORLD -987654

#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]
#define INFOG(I) infog[(I)-1]
#define INFO(I) info[(I)-1]
#define RINFOG(I) rinfog[(I)-1]
#define RINFO(I) rinfo[(I)-1]

/** constructor */
TMumpsSolver::TMumpsSolver(int N_Eqns, int M_dist_Nz, int *M_dist_Irn, int *M_dist_Jcn, int N_Rhs)
{
   int i, j, k, rank;
   double t1;

#ifdef _SMPI
  
  Comm = MPI_COMM_WORLD;
  MPI_Comm_rank(Comm, &rank);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
//    Comm = TDatabase::ParamDB->Comm;
  MPI_Bcast(&N_Eqns,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&M_dist_Nz,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&N_Rhs,1,MPI_INT,0,MPI_COMM_WORLD);

#endif
   
#ifdef _MPI
   Comm = TDatabase::ParamDB->Comm;
   MPI_Comm_rank(Comm, &rank);
#endif
   
   t1 = MPI_Wtime();

   // calls initmumps for initializing mumps solver
   id.job = JOB_INIT;
   
   // host will take part in computations
   id.par = 1;  
   
   id.comm_fortran = MPI_Comm_c2f(Comm);
//    id.comm_fortran = 1;
   
   // unsymmetric matrix
   id.sym = 0;       

   //   id.ooc_tmpdir=/home/sashikum/MumpsOutOfCore
   //  id.comm_fortran=USE_COMM_WORLD;

  //init the MUMPS solver
   dmumps_c(&id);

   id.ICNTL(1) = 1; 
   id.ICNTL(2) = -1; 
   id.ICNTL(3) = -1; 
   id.ICNTL(4) = 1; 

   // martix will be given in assembled form
   id.ICNTL(5) = 0; 
   
   id.ICNTL(6) = 1; 

   // 3 is good or 5 is slightly better than 3
   // 0 AMD, 2 AMF, 3 SCOTCH, 4 PORD, 5  METIS, 6 QAMD, 7  Automatic choice (default)
   id.ICNTL(7) = 5;
   
//   id.ICNTL(8)=4; 
//   id.ICNTL(10)=5;  // maximum number of allowed iterative refinement steps
//   id.ICNTL(13)=1; 
   id.ICNTL(14)=50; //the percentage increase in the estimated working space
   
   // 3 - structure and values on slave process
#ifdef _MPI
   id.ICNTL(18) = 3;
#endif
   
#ifdef _SMPI
    id.ICNTL(18) = 0;
#endif
   
   //id.ICNTL(22)=1; // out of core
  //  id.ICNTL(23) = 100;
   // parallel analysis: 0 automatic, 1 sequential, 2 parallel
   id.ICNTL(28)=0; 
   
   
//   id.ICNTL(29)=1; // 1 - PT-SCOTCH, 2 ParMetis parallel ordering tool
//   id.CNTL(1)=0.00001; //threshold for numerical pivoting a larger value of CNTL(1) increases 
                      // fill-in but leads to a more accurate factorization
//   id.CNTL(2)=1.e-6;  // stoping criterion for iterative refinement steps

   //set the structure into the MUMPS solver
#ifdef _SMPI
   if(rank==0)  
   {
	id.nz   = M_dist_Nz;
	id.irn  = M_dist_Irn;
	id.jcn  = M_dist_Jcn;

	if(rank == 0)
	  id.n  = N_Eqns;

	if(N_Rhs>1)
	{
	  id.nrhs = N_Rhs;
	  id.lrhs = N_Eqns;
	}
	
	if(rank==0)
	OutPut("MUMPS Analysis : ");
      
	// call the MUMPS for analysis
   }
#endif
 
#ifdef _MPI 
   id.nz_loc   = M_dist_Nz;
   id.irn_loc  = M_dist_Irn;
   id.jcn_loc  = M_dist_Jcn;

   if(rank == 0)
     id.n  = N_Eqns;

   if(N_Rhs>1)
   {
     id.nrhs = N_Rhs;
     id.lrhs = N_Eqns;
   }
   
  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE != 24)
   OutPut("MUMPS Analysis : ");
#endif
  // call the MUMPS for analysis
  id.job=JOB_ANALYSIS;
  dmumps_c(&id);

// printf("Mumps Rank %d Nz %d N_Eqns %d  \n", rank, M_dist_Nz, N_Eqns);

 #ifdef _SMPI
   if(rank==0)
#endif
     {	
	    OutPut("Time taken : "<< MPI_Wtime()- t1 << endl);
     } 

#ifdef _MPI
if(rank==0&& TDatabase::ParamDB->SC_VERBOSE != 24)
   OutPut("Time taken : "<< MPI_Wtime()- t1 << endl);
#endif
// MPI_Finalize();
// exit(0);
}


void TMumpsSolver::FactorizeAndSolve(double *Mat_loc, double *rhs)
{               	  	        
 
  int rank;
  MPI_Comm_rank(Comm, &rank);
  double t1 = MPI_Wtime();
  // for centralized input to the host
  if(rank == 0)
    id.rhs = rhs;
 
  // define the local system matrix
  
#ifdef _SMPI
  if(rank==0)
  {
  id.a  = Mat_loc;  
  }  
#endif 

#ifdef _MPI
  id.a_loc  = Mat_loc;
#endif
  

//   if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE != 24)
//    OutPut("MUMPS Factorization : ");
  
  // use iterative approach to residual 
  id.ICNTL(10)=1;
  // residual value: 2 - machine precision, else proved residual 
  id.ICNTL(11)= 2;
   
  id.job=JOB_FACTORIZ;
  dmumps_c(&id);
  
  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE != 24)
   OutPut("Time taken for Factorization: "<< MPI_Wtime()- t1 << endl);
  
  FactorFlag = TRUE;   

#ifdef _SMPI
  if(rank==0)
#endif
  {
    if(id.INFOG(1)<0)
    {
      printf("MUMPS Factorize failed INFOG(1) %d \n",  id.INFOG(1));
      MPI_Finalize();
      exit(0); 
    }
  }

  id.job=JOB_SOLVE;
  dmumps_c(&id);

#ifdef _SMPI
  if(rank==0)
#endif
  { 
    if(id.INFOG(1)<0)
    {
      printf("MUMPS Solve failed INFOG(1) %d \n",  id.INFOG(1));
      MPI_Finalize();
      exit(0); 
    }
  }
}

void TMumpsSolver::Solve(double *Mat_loc, double *rhs)
{
  int rank;
  MPI_Comm_rank(Comm, &rank);

  // for centralized input to the host
  if(rank == 0)
    id.rhs = rhs;

  // define the local system matrix
#ifdef _MPI
  id.a_loc  = Mat_loc;
#endif
#ifdef _SMPI
  id.a  = Mat_loc;
#endif


  id.job=JOB_SOLVE;
  dmumps_c(&id);
  
  if(id.INFOG(1)<0)
  {
    printf("MUMPS Solve failed INFOG(1) %d \n",  id.INFOG(1));
    MPI_Finalize();
    exit(0); 
  }

//   if(id.INFO(1)==8) // more than ICNTL(10) iterative refinement is required
//   {
//     if(rank==0)
//     printf("MUMPS: Number of iterative refinement steps: %d \n", id.INFOG(15));
// 
// //     id.job=JOB_FACTORIZ;
// //     dmumps_c(&id);
// // 
// //     id.job=JOB_SOLVE;
// //     dmumps_c(&id);
// 
// //     if(rank==0)
// //      printf("MUMPS: Number of iterative refinement steps (after factorize): %d \n", id.INFOG(15));
//    }
}

void TMumpsSolver::Clean()
{
  id.job=JOB_END;
  dmumps_c(&id);
  printf("MUMPS: Exiting\n");
}




#endif
