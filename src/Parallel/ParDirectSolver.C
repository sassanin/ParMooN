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
// @(#)ParDirectSolver.h
//
// Class:      TParDirectSolver
// Purpose:    Class for interfacing ParMooN with MumpsSolver
//
// Author:     Sashikumaar Ganesan & Abdus Shamim (27.04.15)
//
// History:    Start of implementation 27.04.15 (Sashikumaar Ganesan & Abdus Shamim)
//	       Adding INTERNAL_PROJECT_PRESSURE (Raviteja Meesala)
// ======================================================================= 
#ifdef _MPI
#  include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <LinAlg.h>

#include <ParDirectSolver.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>
#include <MumpsSolver.h>

TParDirectSolver::TParDirectSolver(TParFECommunicator3D *parcomm,TParFECommunicator3D *parcomm_p,TSquareMatrix3D **mat,TMatrix3D **matB)
{
  int rank,size;
  NSEType = TDatabase::ParamDB->NSTYPE;  
  
  ParComm   = parcomm;
  NDof      = ParComm->GetNDof();
  Comm      = TDatabase::ParamDB->Comm;
  N_rhs     = 1;
  Mat       = mat[0];
  DSType    = TDatabase::ParamDB->DSType;
  
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);  
  
  if(matB != NULL)
   { SystemType = NSE; }
  else
   { SystemType = SCALAR; }   
  
  if(SystemType==NSE)
  {
    ParComm_P   = parcomm_p;
    NDof_P      = ParComm_P->GetNDof();
    NDof_U      = ParComm->GetNDof();
    NDof        = 3*NDof_U + NDof_P;   
    
   switch(NSEType)
     {
      case 1:
             if(rank == 0)
              printf("ParDirectSolver NSEType %d not yet implemented\n",NSEType);
              MPI_Finalize();
              exit(0);
      break;

      case 2:
             MatB  = matB[0];   
             MatBT = matB[3];
            
             MatA11 = mat[0];    
	     MatB1  = matB[0];   MatB2  = matB[1];   MatB3  = matB[2];      
	     MatBT1 = matB[3];   MatBT2 = matB[4];   MatBT3 = matB[5]; 
	
      break;
    
      case 3:
             if(rank == 0)
              printf("ParDirectSolver NSEType %d not yet implemented\n",NSEType);
              MPI_Finalize();
              exit(0);
      break;
    
      case 4:
             MatB  = matB[0];   
             MatBT = matB[3];
    
             MatA11 = mat[0];    MatA12 = mat[1];    MatA13 = mat[2];	MatBT1 = matB[3];
             MatA21 = mat[3];    MatA22 = mat[4];    MatA23 = mat[5];	MatBT2 = matB[4];
             MatA31 = mat[6];    MatA32 = mat[7];    MatA33 = mat[8];	MatBT3 = matB[5];
             MatB1  = matB[0];   MatB2  = matB[1];   MatB3  = matB[2];      //[0]
       break;
      
      default:
            OutPut("Unknown NSETYPE " << NSEType <<"  it must be 1 to 4" << endl);
            exit(4711);;      
      
     }  
   }
  else //scalar problem
   {
    MatB = NULL;
   }
  
  if(DSType == 1)
   {
    if(rank == 0)
      printf("ParDirectSolver Type ---(Select _OMPONLY in makefile)---> ParDiso\n");
      MPI_Finalize();
      exit(0);
   }
  else if(DSType == 2) // Mumps
   {
    if(rank == 0 && TDatabase::ParamDB->SC_VERBOSE != 24)
      printf("ParDirectSolver Type ------> Mumps\n");
    
    if(SystemType==NSE)
     { 
      switch(NSEType)
       {
        case 1:
             if(rank == 0)
              printf("ParDirectSolver NSEType %d not yet implemented\n",NSEType);
              MPI_Finalize();
              exit(0);
        break;

        case 2:
              InitMumps_NSE2();
        break;
    
        case 3:
            if(rank == 0)
             printf("ParDirectSolver NSEType %d not yet implemented\n",NSEType);
             MPI_Finalize();
             exit(0);
        break;
    
        case 4:
              InitMumps_NSE4();
        break; 
       } //  switch(NSEType)     
     }
    else
     {
      InitMumps_Scalar();
     }
    
    //initilize the Mumps solver
    Mumps = new TMumpsSolver(N_Eqns, N_Nz, I_rn, J_cn, N_rhs);
   }
  else
   {
    printf("Select ParDirectSolver Type\n");
    MPI_Finalize();
    exit(0);
   }
  
}

TParDirectSolver::~TParDirectSolver()
{
  if(DSType == 2)
    Mumps->Clean();

  delete [] MatLoc;
  delete [] OwnRhs;
  delete [] I_rn;
  delete [] J_cn;
  delete [] GlobalRhs;
//   delete Mumps;
}

void TParDirectSolver::InitMumps_Scalar()
{
  int i,j,k,l,m,t;
  int *Master_P,*local2global_P;
  int *Master         = ParComm->GetMaster();
  int *local2global   = ParComm->Get_Local2Global();

  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
    RowPtr = Mat->GetRowPtr();
    KCol   = Mat->GetKCol();
    
    N_Nz = 0;
    for(i=0;i<NDof;i++)
      if(Master[i] == rank)
	N_Nz += RowPtr[i+1] - RowPtr[i];
    
    N_Master = ParComm->GetN_Master();
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
    k = 0;
    for(i=0;i<NDof;i++)
    {
      if(Master[i] == rank)
      {
	for(j=RowPtr[i];j<RowPtr[i+1];j++)
	{
	  I_rn[k] = local2global[i] + 1;              //fortran format
	  J_cn[k] = local2global[KCol[j]] + 1;        //fortran format
	  k++;
	}
      }
    }

    OwnRhs = new double[N_Master];
    
    MPI_Allreduce(&N_Master, &N_Eqns, 1, MPI_INT, MPI_SUM, Comm);
    
    GlobalRhsSize = N_Eqns;
    
    if(rank == 0)
      GlobalRhs = new double[GlobalRhsSize];
}

void TParDirectSolver::InitMumps_NSE2()
{
             
  int i,j,k,l,m,t,N_Active, row_index, index_disp2U, index_disp3U;
  int *Master_P,*local2global_P;
  int *Master         = ParComm->GetMaster();
  int *local2global   = ParComm->Get_Local2Global();

  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  
    Master_P       = ParComm_P->GetMaster();
    local2global_P = ParComm_P->Get_Local2Global();
  
    N_Active = Mat->GetActiveBound();
     
    //compute N_Nz in sqmatrices i.e. A blocks
    RowPtr = Mat->GetRowPtr();
    KCol   = Mat->GetKCol();
    
    RowPtr_B = MatB->GetRowPtr();
    KCol_B   = MatB->GetKCol();

    RowPtr_BT = MatBT->GetRowPtr();
    KCol_BT   = MatBT->GetKCol();
    
    N_Nz_U = 0;
    N_Nz_P = 0;
    
   for(i=0;i<NDof_U;i++)
    {
      if(Master[i] == rank)
      {
	N_Nz_U += RowPtr[i+1] - RowPtr[i];		//for each A type blocks
	N_Nz_P += RowPtr_BT[i+1] - RowPtr_BT[i];	//for each BT type blocks
      }
    }
    
    for(i=0;i<NDof_P;i++)
    {
      if(Master_P[i] == rank)
      {
	N_Nz_P += RowPtr_B[i+1] - RowPtr_B[i];		//for each B type blocks
      }
    }
    
    N_Nz = 3*N_Nz_P + 3*N_Nz_U; 	
      
    N_Master_U  = ParComm->GetN_Master();
    N_Master_P  = ParComm_P->GetN_Master();
    N_Master    = N_Master_P + 3*N_Master_U;
    
    MPI_Allreduce(&N_Master_U, &Global_N_DOF_U, 1, MPI_INT, MPI_SUM, Comm);
    
    OwnRhs = new double[N_Master];
    
    MPI_Allreduce(&N_Master, &N_Eqns, 1, MPI_INT, MPI_SUM, Comm);
    
    GlobalRhsSize = N_Eqns;
    
    if(rank == 0)
      GlobalRhs = new double[GlobalRhsSize];
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
    
    index_disp3U = 3*Global_N_DOF_U;
    index_disp2U = 2*Global_N_DOF_U;
 /** U1 component */   
    k = 0;
    for(i=0;i<NDof_U;i++)
    {
      if(Master[i] == rank)
       {
        row_index = local2global[i] + 1; //fortran format

	 for(j=RowPtr[i];j<RowPtr[i+1];j++)
	  {
	   I_rn[k] =  row_index;     
	   J_cn[k] =  local2global[KCol[j]] + 1;        //fortran format
	   k++;
	  }//for(j=
	 
	
	 //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	   {
	    I_rn[k] = row_index;                
	    J_cn[k] = index_disp3U + local2global_P[KCol_BT[j]] + 1;       //fortran format
	    k++;
	   }//for(j=
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
    
  /** U2 component */      
   for(i=0;i<NDof_U;i++)
    {
     if(Master[i] == rank)
      {
       row_index = Global_N_DOF_U + local2global[i] + 1;              //fortran format

       for(j=RowPtr[i];j<RowPtr[i+1];j++)
        {
         I_rn[k] = row_index;
         J_cn[k] = Global_N_DOF_U + local2global[KCol[j]] + 1;        //fortran format
         k++;
       }//for(j=
	 
	
         //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] = row_index;  
	    J_cn[k] = index_disp3U + local2global_P[KCol_BT[j]] + 1;       //fortran format
	    k++;
	  }//for(j=
	  
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
   
       
  /** U3 component */    
    for(i=0;i<NDof_U;i++)
    {
      if(Master[i] == rank)
      {
       row_index = index_disp2U + local2global[i] + 1;              //fortran format
       
	for(j=RowPtr[i];j<RowPtr[i+1];j++)
	 {
	  I_rn[k] = row_index;
          J_cn[k] = index_disp2U + local2global[KCol[j]] + 1;        //fortran format
           k++;
	  }//for(j=
	 
	
	  //Mat BT(m)
        if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] = row_index;
	    J_cn[k] = index_disp3U + local2global_P[KCol_BT[j]] + 1;       //fortran format
	    k++;
	  }//for(j=  
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
     
    /** pressure */
    for(i=0;i<NDof_P;i++)
    {
      
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      {  
        printf("ParDirectSolver->InitMumps_NSE2(): NSEType 2  INTERNAL_PROJECT_PRESSURE Not yet implemented !!!!\n");
        MPI_Finalize();
        exit(0); 
      }
      
      if(Master_P[i] == rank)
      {
       row_index = index_disp3U + local2global_P[i] + 1;                //fortran format
	//Mat B(l)
        for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;
	    J_cn[k] = local2global[KCol_B[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
 
         for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;
	    J_cn[k] = Global_N_DOF_U + local2global[KCol_B[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
	  
	 for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;
	    J_cn[k] = index_disp2U + local2global[KCol_B[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
    /** since the BT entries of Dirichlet U DOF row is not added */
    N_Nz = k ; 
   
}


void TParDirectSolver::InitMumps_NSE4()
{
  int i,j,k,l,m,t,N_Active, row_index, index_disp2U, index_disp3U;
  int *Master_P,*local2global_P;
  int *Master         = ParComm->GetMaster();
  int *local2global   = ParComm->Get_Local2Global();

  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

    Master_P       = ParComm_P->GetMaster();
    local2global_P = ParComm_P->Get_Local2Global();
    
    N_Active = Mat->GetActiveBound();
    
    //compute N_Nz in sqmatrices i.e. A blocks
    RowPtr = Mat->GetRowPtr();
    KCol   = Mat->GetKCol();
    
    RowPtr_B = MatB->GetRowPtr();
    KCol_B   = MatB->GetKCol();

    RowPtr_BT = MatBT->GetRowPtr();
    KCol_BT   = MatBT->GetKCol();
    
    N_Nz_U = 0;
    N_Nz_P = 0;
    
    for(i=0;i<NDof_U;i++)
    {
      if(Master[i] == rank)
      {
	N_Nz_U += RowPtr[i+1] - RowPtr[i];		//for each A type blocks
	N_Nz_P += RowPtr_BT[i+1] - RowPtr_BT[i];	//for each BT type blocks
      }
    }

    for(i=0;i<NDof_P;i++)
    {
      if(Master_P[i] == rank)
      {
	N_Nz_P += RowPtr_B[i+1] - RowPtr_B[i];		//for each B type blocks
      }
    }

    N_Nz = 3*N_Nz_P + 9*N_Nz_U; 	//total non zeros in system matrix ( 3*(B + BT) + 9*A )
      
    N_Master_U  = ParComm->GetN_Master();
    N_Master_P  = ParComm_P->GetN_Master();
    N_Master    = N_Master_P + 3*N_Master_U;
    
    OwnRhs = new double[N_Master];
    
    MPI_Allreduce(&N_Master_U, &Global_N_DOF_U, 1, MPI_INT, MPI_SUM, Comm);
    
    MPI_Allreduce(&N_Master, &N_Eqns, 1, MPI_INT, MPI_SUM, Comm);
    
    GlobalRhsSize = N_Eqns;
    
    if(rank == 0)
      GlobalRhs = new double[GlobalRhsSize];
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
    index_disp3U = 3*Global_N_DOF_U;
    index_disp2U = 2*Global_N_DOF_U;
    
    k = 0;
    for(i=0;i<NDof_U;i++)
    {
      
      if(Master[i] == rank)
      {
	row_index = local2global[i] + 1; //fortran format
	
	
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  local2global[KCol[j]] + 1;        //fortran format
	      k++;
	      
	      if(i<N_Active)
	    {	
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  Global_N_DOF_U + local2global[KCol[j]] + 1;        //fortran format
	      k++;
	    }
	    if(i<N_Active)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  index_disp2U + local2global[KCol[j]] + 1;        //fortran format
	      k++;     
	    } 
	    }//for(j=
	
	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = index_disp3U + local2global_P[KCol_BT[j]] + 1;       //fortran format
	    k++;
	  }//for(j=
	  
      
      }//if(Master[i] == rank)
    }//for(i=0;i<
    
        for(i=0;i<NDof_U;i++)
    {
      
      if(Master[i] == rank)
      {
	row_index = Global_N_DOF_U + local2global[i] + 1; //fortran format
	
	
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	          if(i<N_Active)
	    {	
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  local2global[KCol[j]] + 1;        //fortran format
	      k++;
	    }
	  
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  Global_N_DOF_U + local2global[KCol[j]] + 1;        //fortran format
	      k++;
	  
	    if(i<N_Active)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  index_disp2U + local2global[KCol[j]] + 1;        //fortran format
	      k++;     
	    } 
	    }//for(j=
	
	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = index_disp3U + local2global_P[KCol_BT[j]] + 1;       //fortran format
	    k++;
	  }//for(j=
	  
      
      }//if(Master[i] == rank)
    }//for(i=0;i<
    
            for(i=0;i<NDof_U;i++)
    {
      
      if(Master[i] == rank)
      {
	row_index = index_disp2U + local2global[i] + 1; //fortran format
	
	
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	          if(i<N_Active)
	    {	
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  local2global[KCol[j]] + 1;        //fortran format
	      k++;
	    }
	    
	   if(i<N_Active)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  Global_N_DOF_U + local2global[KCol[j]] + 1;        //fortran format
	      k++;
	    }
	   
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  index_disp2U + local2global[KCol[j]] + 1;        //fortran format
	      k++;     
	   
	    }//for(j=
	
	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = index_disp3U + local2global_P[KCol_BT[j]] + 1;       //fortran format
	    k++;
	  }//for(j=
	  
      
      }//if(Master[i] == rank)
    }
    
    
    for(i=0;i<NDof_P;i++)
    {
//       if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
//       {  
//         printf("ParDirectSolver->InitMumps_NSE2(): NSEType 4  INTERNAL_PROJECT_PRESSURE Not yet implemented !!!!\n");
//         MPI_Finalize();
//         exit(0); 
//       }
      
      if(Master_P[i] == rank)
      {
	row_index = index_disp3U + local2global_P[i] + 1;                //fortran format
	
	//Mat B1
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = local2global[KCol_B[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
	  
	  	//Mat B2
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = Global_N_DOF_U + local2global[KCol_B[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
	  
	  	//Mat B3
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = index_disp2U + local2global[KCol_B[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
	  
	  

      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
     /** since the BT entries of Dirichlet U DOF row is not added */
    N_Nz = k ; 

  
  
}




void TParDirectSolver::AssembleMatrix()
{ 
  int i,j,k,l,m,t;
  int *Master = ParComm->GetMaster();;
  double *EntriesA;
 
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  if(MatB == NULL)
  {
    EntriesA = Mat->GetEntries();
    k = 0;
    for(i=0;i<NDof;i++)
    {
      if(Master[i] == rank)
      {
        for(j=RowPtr[i];j<RowPtr[i+1];j++)
        {
          MatLoc[k] = EntriesA[j];
	  k++;
        }
      }
    }
  }
  else
  {
    switch(NSEType)
       {
        case 1:
             if(rank == 0)
               printf("ParDirectSolver NSEType %d not yet implemented\n",NSEType);
             MPI_Finalize();
             exit(0);
        break;

        case 2:
              AssembleMatrix_NSE2();
        break;
    
        case 3:
            if(rank == 0)
              printf("ParDirectSolver NSEType %d not yet implemented\n",NSEType);
            MPI_Finalize();
            exit(0);
        break;
    
        case 4:
	      AssembleMatrix_NSE4();
        break; 
       } //  switch(NSEType)   
  }
 
}

void TParDirectSolver::AssembleMatrix_NSE2()
{
  int rank,size, N_Active;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  int i,j,k,l,m,t;
  int *Master   = ParComm->GetMaster();
  int *Master_P = ParComm_P->GetMaster();
  
  double *EntriesA = MatA11->GetEntries();
  double *EntriesB = MatBT1->GetEntries();
  
  N_Active = Mat->GetActiveBound();
  

  
  k = 0;
  for(i=0;i<NDof_U;i++)
  {
    if(Master[i] == rank)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
	MatLoc[k] = EntriesA[j];
	
	k++;
      }//for(j=
      
      //Mat B1T(m)
      if(i<N_Active)
      for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
      {
	MatLoc[k] = EntriesB[j];

	k++;
      }//for(j=
    }//if(Master_P[i] == rank)
  }//for(i=0;i<
  
  EntriesB = MatBT2->GetEntries();
  for(i=0;i<NDof_U;i++)
  {
    if(Master[i] == rank)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
	MatLoc[k] = EntriesA[j];

	k++;
      }//for(j=

      //Mat B2T(m)
      if(i<N_Active)
      for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
      {
	MatLoc[k] = EntriesB[j];

	k++;
      }//for(j=
    }//if(Master_P[i] == rank)
  }//for(i=0;i<

  EntriesB = MatBT3->GetEntries();
  for(i=0;i<NDof_U;i++)
  {
    if(Master[i] == rank)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
       MatLoc[k] = EntriesA[j];

	k++;
      }//for(j=

      //Mat BT(m)
      if(i<N_Active)
      for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
      {
	MatLoc[k] = EntriesB[j];

	k++;
      }//for(j=
    }//if(Master_P[i] == rank)  
  }//for(i=0;i<
  
  EntriesB = MatB1->GetEntries();
  double *EntriesB2 = MatB2->GetEntries();
  double *EntriesB3 = MatB3->GetEntries();
  
  for(i=0;i<NDof_P;i++)
  {
    if(Master_P[i] == rank)
    {
      //Mat B(l)
      for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
      {
	MatLoc[k] = EntriesB[j];

	k++;
      }//for(j=
      
      for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
      {
	MatLoc[k] = EntriesB2[j];

	k++;
      }//for(j=

      for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
      {
	MatLoc[k] = EntriesB3[j];

	k++;
      }//for(j=
    }//if(Master_P[i] == rank)  
  }//for(i=0;i<


  
}

void TParDirectSolver::AssembleMatrix_NSE4()
{
  int rank,size, N_Active;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  int i,j,k,l,m,t;
  int *Master   = ParComm->GetMaster();
  int *Master_P = ParComm_P->GetMaster();
  
  double* EntriesA[9];
  double* EntriesB[6];
  
  EntriesA[0] = MatA11->GetEntries();
  EntriesA[1] = MatA12->GetEntries();
  EntriesA[2] = MatA13->GetEntries();
  
  EntriesA[3] = MatA21->GetEntries();
  EntriesA[4] = MatA22->GetEntries();
  EntriesA[5] = MatA23->GetEntries();
  
  EntriesA[6] = MatA31->GetEntries();
  EntriesA[7] = MatA32->GetEntries();
  EntriesA[8] = MatA33->GetEntries();
  
  EntriesB[0] = MatBT1->GetEntries();
  EntriesB[1] = MatBT2->GetEntries();
  EntriesB[2] = MatBT3->GetEntries();
  
  EntriesB[3] = MatB1->GetEntries();
  EntriesB[4] = MatB2->GetEntries();
  EntriesB[5] = MatB3->GetEntries();
  
  N_Active = Mat->GetActiveBound();

     k = 0;
    for(i=0;i<NDof_U;i++)
    {
      
      if(Master[i] == rank)
      {

	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      MatLoc[k] = EntriesA[0][j]; 
	      k++;
	      
	      if(i<N_Active)
	      {
	       MatLoc[k] = EntriesA[1][j];
	      k++;
	      }
	      
	       if(i<N_Active)
	       {
		 MatLoc[k] = EntriesA[2][j];
	         k++;
	       }

	    }


	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[0][j];
	    k++;
	  }//for(j=
	  
     
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
    
        for(i=0;i<NDof_U;i++)
    {
      
      if(Master[i] == rank)
      {

	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      if(i<N_Active)
	      {
	       MatLoc[k] = EntriesA[3][j];
	       k++;
	      }
	        MatLoc[k] = EntriesA[4][j];
	        k++;
	      
	       if(i<N_Active)
	       {
		 MatLoc[k] = EntriesA[5][j];
	          k++; 
	      }

	      
	    }

	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[1][j];     
	    k++;
	  }//for(j=
	  
     
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
    
        for(i=0;i<NDof_U;i++)
    {
      
      if(Master[i] == rank)
      {

	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      if(i<N_Active)
	      {
		MatLoc[k] = EntriesA[6][j];
	      k++;
		
	      }
	      
	      if(i<N_Active)
	      {
		MatLoc[k] = EntriesA[7][j];
	      k++;
	      }
	      
	      MatLoc[k] = EntriesA[8][j];
	      k++;  
	    }

	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[2][j]; 
	    k++;
	  }//for(j=
	  
     
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
    
    
    
    for(i=0;i<NDof_P;i++)
    {
      if(Master_P[i] == rank)
      {
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[3][j];
	    k++;
	  }
	  
	  // if(i<N_Active)
	   for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[4][j];
	    k++;
	  }
	  
	//   if(i<N_Active)
	   for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[5][j];
	    k++;
	  }
	
      }//if(Master_P[i] == rank)
    }//for(i=0;i<

}

void TParDirectSolver::GetRhs(double *Rhs)
{
  int i,j,k,t;
  int *Master = ParComm->GetMaster();
  
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  if(MatB == NULL)
  {
    t = 0;
    for(i=0;i<NDof;i++)
    {
      if(Master[i] == rank)
      {
	OwnRhs[t++] = Rhs[i];
      }
    }
    
    ParComm->GatherToRoot(GlobalRhs, GlobalRhsSize, OwnRhs, N_Master, 0);
  }
  else
  {
    double *temp,*temp2;
    int *Master_P = ParComm_P->GetMaster();
    t = 0;
    for(i=0;i<NDof_U;i++)
    {
      if(Master[i] == rank)
      {
	OwnRhs[t               ] = Rhs[i           ];
	OwnRhs[t +   N_Master_U] = Rhs[i +   NDof_U];
	OwnRhs[t + 2*N_Master_U] = Rhs[i + 2*NDof_U];
	t++;
      }
    }
    
    //send Ux
    ParComm->GatherToRoot(GlobalRhs, GlobalRhsSize, OwnRhs, N_Master_U, 0);
    
    //send Uy
    temp  = GlobalRhs+Global_N_DOF_U;
    temp2 = OwnRhs+N_Master_U;
    ParComm->GatherToRoot(temp, GlobalRhsSize, temp2, N_Master_U, 0);
    
    //send Uz
    temp  = GlobalRhs+2*Global_N_DOF_U;
    temp2 = OwnRhs+2*N_Master_U;
    ParComm->GatherToRoot(temp, GlobalRhsSize, temp2, N_Master_U, 0);
    
    //Send P
    k = 0;
    for(i=0;i<NDof_P;i++)
    {
      if(Master_P[i] == rank)
      {
	OwnRhs[k + 3*N_Master_U] = Rhs[i + 3*NDof_U];
	k++;
      }
    }
    
    temp  = GlobalRhs+3*Global_N_DOF_U;
    temp2 = OwnRhs+3*N_Master_U;
    ParComm->GatherToRoot(temp	, GlobalRhsSize, temp2, N_Master_P, 0);

  }
}

void TParDirectSolver::UpdateSol(double *Sol)
{
  int i,j,k,t;
  int *Master = ParComm->GetMaster();
  
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  if(MatB == NULL)
  {
    ParComm->ScatterFromRoot(GlobalRhs, GlobalRhsSize, OwnRhs, N_Master, 0);
  
    t = 0;
    for(i=0;i<NDof;i++)
    {
      if(Master[i] == rank)
      {
	Sol[i] = OwnRhs[t++];
      }
    }
    ParComm->CommUpdate(Sol);
  }
  else
  {
//     if(rank == 0)
//     {
//       double sum=0;
//       for(i=0;i<GlobalRhsSize;i++)
//       {
// 	sum += GlobalRhs[i]*GlobalRhs[i];
// // 	GlobalRhs[i] = i;
//       }
// //       cout<<"norm of sol :: "<<sqrt(sum)<<endl;
// //       cout<<"GlobalRhsSize :: "<<GlobalRhsSize<<endl;
//     }
// //     MPI_Finalize();
// //     exit(0);
    
    
    int *Master_P = ParComm_P->GetMaster();
    double *temp,*temp2;
    
    //Recv Ux
    ParComm->ScatterFromRoot(GlobalRhs, GlobalRhsSize, OwnRhs, N_Master_U, 0);
    
    //Recv Uy
    temp  = GlobalRhs+Global_N_DOF_U;
    temp2 = OwnRhs+N_Master_U;
    ParComm->ScatterFromRoot(temp, GlobalRhsSize, temp2, N_Master_U, 0);
    
    //Recv Uz
    temp  = GlobalRhs+2*Global_N_DOF_U;
    temp2 = OwnRhs+2*N_Master_U;
    ParComm->ScatterFromRoot(temp, GlobalRhsSize, temp2, N_Master_U, 0);
  
    t = 0;
    for(i=0;i<NDof_U;i++)
    {
      if(Master[i] == rank)
      {
	Sol[i           ] = OwnRhs[t               ];
	Sol[i +   NDof_U] = OwnRhs[t +   N_Master_U];
	Sol[i + 2*NDof_U] = OwnRhs[t + 2*N_Master_U];
	t++;
      }
    }
    ParComm->CommUpdate(Sol);
    
    temp  = GlobalRhs + 3*Global_N_DOF_U;
    temp2 = OwnRhs + 3*N_Master_U;
    ParComm_P->ScatterFromRoot(temp, GlobalRhsSize, temp2, N_Master_P, 0);
    
    t=0;
    for(i=0;i<NDof_P;i++)
    {
      if(Master_P[i] == rank)
      {
	Sol[i + 3*NDof_U] = OwnRhs[t + 3*N_Master_U];
	t++;
      }
    }
    ParComm_P->CommUpdate(Sol+3*NDof_U);
    
  }
}

void TParDirectSolver::Solve(double *Sol, double *Rhs, bool Factorize)
{
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  
  double t = MPI_Wtime();
  int i,j;
 
  if(DSType == 2)
  {
       
           
    if(NSEType == 4){
      if(rank == 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE == 0 )
	Rhs[3*NDof_U] = 0.0;
    }
    GetRhs(Rhs);
    
    if(Factorize)
    {
      AssembleMatrix();
 
      Mumps->FactorizeAndSolve(MatLoc,GlobalRhs);
 
    }
    else
    {
      Mumps->Solve(MatLoc,GlobalRhs);
    }
    
    UpdateSol(Sol);
  }
  else
  {
    printf("Select ParDirectSolver Type\n");
    MPI_Finalize();
    exit(0);
  }
  
//   printf("time taken for solving::%lf\n",MPI_Wtime()-t);
}

#else

//--------------------------------------------------------------------------------------------------------//
//						OMP ONLY
//--------------------------------------------------------------------------------------------------------//


#ifdef _OMPONLY

#include "omp.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <ParDirectSolver.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParDiso.h>

TParDirectSolver::TParDirectSolver(TSquareMatrix3D *mat)
{
  N_rhs   = 1;
  Mat     = mat;
  DSType  = TDatabase::ParamDB->DSType;

  if(DSType == 1)
  {
      printf("ParDirectSolver Type ------> ParDiso\n");
    
    InitPardiso();
    ParDiso = new TParDiso(N_Eqns, N_Nz, RowPtr, KCol);
  }
  else
  {
    printf("Select ParDirectSolver Type as 1 for ParDiso\n");
    exit(0);
  }
  
}

TParDirectSolver::~TParDirectSolver()
{
  ParDiso->Clean(Mat_Values);
  delete [] Mat_Values; Mat_Values = NULL;
}

void TParDirectSolver::InitPardiso()
{
  int i,j,k,t;
   
  RowPtr = Mat->GetRowPtr();
  KCol   = Mat->GetKCol();

  N_Nz   = Mat->GetN_Entries();
  
  N_Eqns = Mat->GetN_Rows();
  
  Mat_Values = new double[N_Nz];
}

void TParDirectSolver::AssembleMatrix(TSquareMatrix3D *matrix)
{
  int i;
  double *EntriesA = matrix->GetEntries();

  for(i=0;i<N_Nz;i++)
    Mat_Values[i] = EntriesA[i];
}

void TParDirectSolver::Solve(double *Sol, double *Rhs, bool Factorize)
{
  int i,j;
  double t = GetTime();
  
  if(DSType == 1)
  {
    ParDiso->FactorizeAndSolve(Mat_Values, Rhs, Sol, Factorize);
  }
  else
  {
    printf("Select ParDirectSolver Type\n");
    exit(0);
  }
  
  printf("time taken for solving::%lf\n",GetTime()-t);
 
}

#endif

#endif















   
