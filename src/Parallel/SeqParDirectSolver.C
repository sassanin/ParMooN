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
// @(#)SeqParDirectSolver.h
//
// Class:      TSeqParDirectSolver.h
// Purpose:    Class for interfacing ParMooN with MumpsSolver
//
// Author:     Raviteja
//
// History:    Start of implementation 05.05.16 
// History:    Implement SMPI for TNSE2D - Ankit
//
// ======================================================================= 
#ifdef _SMPI
#  include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <LinAlg.h>

#include <SeqParDirectSolver.h>
#include <Database.h>

#ifdef __2D__
#include <SquareMatrix2D.h>
#else 
#include <SquareMatrix3D.h>
#endif

#include <MumpsSolver.h>
#define FORTRAN 1

#ifdef __3D__
TSeqParDirectSolver::TSeqParDirectSolver(int dim,int N_U,int N_P,TSquareMatrix3D **mat,TMatrix3D **matB)
{
  int rank,size;
  NSEType = TDatabase::ParamDB->NSTYPE;  
  
//   ParComm   = parcomm;
//   NDof      = ParComm->GetNDof();

  Comm      = TDatabase::ParamDB->Comm;
  N_rhs     = 1;
  Mat       = mat[0];
  DSType    = TDatabase::ParamDB->DSType;
  
  NDof = N_U;
  
/*  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);*/  
  
  if(matB != NULL)
   { 
     SystemType = NSE;      
  }
  else
   { SystemType = SCALAR; }   
  
  if(SystemType==NSE)
  {
 
    NDof_P      = N_P;
    NDof_U      = N_U;
    NDof        = dim*NDof_U + NDof_P;   
    
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
  
  if(DSType == 4)
   {
    
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
    printf("Select ParDirectSolver Type\n, DSType %d", DSType);
    MPI_Finalize();
    exit(0);
   }
  
}

void TSeqParDirectSolver::InitMumps_Scalar()
{
  
}

void TSeqParDirectSolver::InitMumps_NSE2()
{
 
}

void TSeqParDirectSolver::InitMumps_NSE4()
{
  int i,j,k,l,m,t,N_Active, row_index,index_dispU, index_disp2U, index_disp3U;

  int rank=0,size;
  int pr_p;
  
    
    N_Active = Mat->GetActiveBound();
    
    //compute N_Nz in sqmatrices i.e. A blocks
    RowPtr = Mat->GetRowPtr();
    KCol   = Mat->GetKCol();
    
    RowPtr_B = MatB->GetRowPtr();
    KCol_B   = MatB->GetKCol();

    RowPtr_BT = MatBT->GetRowPtr();
    KCol_BT   = MatBT->GetKCol();
    
    N_Nz_U = Mat->GetN_Entries();
    N_Nz_P = MatB->GetN_Entries();
    N_Nz = 6*N_Nz_P + 9*N_Nz_U; 	//total non zeros in system matrix ( 3*(B + BT) + 9*A )

    N_Eqns = NDof;
    GlobalRhsSize = NDof;
    

      GlobalRhs = new double[GlobalRhsSize];
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
    index_dispU = NDof_U;
    index_disp3U = 3*NDof_U;
    index_disp2U = 2*NDof_U;
    
    k = 0;
    for(i=0;i<NDof_U;i++)
    {      
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      I_rn[k] =  i + FORTRAN;              //fortran format
	      J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
	      k++;
	      
		  if(i<N_Active)
		  {	
		    I_rn[k] =  i + FORTRAN;             //fortran format
		    J_cn[k] =  index_dispU + KCol[j] + FORTRAN;           //fortran format
		    k++;
		  }
		  if(i<N_Active)
		  {
		    I_rn[k] =  i + FORTRAN;             //fortran format
		    J_cn[k] =  index_disp2U + KCol[j] + FORTRAN;           //fortran format
		    k++;     
		  } 
	    }//for(j=
	
	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  i + FORTRAN;             //fortran format
	    J_cn[k] = index_disp3U + KCol_BT[j] + FORTRAN;      //fortran format
	    k++;
	  }//for(j=
	  
      

    }//for(i=0;i<
    
        for(i=0;i<NDof_U;i++)
    {
      
	row_index = index_dispU + i + FORTRAN; //fortran format
	
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
		if(i<N_Active)
		{	
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
		  k++;
		}
	    
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  index_dispU + KCol[j] + FORTRAN;        //fortran format
	      k++;
	  
	    if(i<N_Active)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  index_disp2U + KCol[j] + FORTRAN;        //fortran format
	      k++;     
	    } 
	    }//for(j=
	
	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = index_disp3U + KCol_BT[j] + FORTRAN;       //fortran format
	    k++;
	  }//for(j=

    }//for(i=0;i<
    
            for(i=0;i<NDof_U;i++)
    {      
	row_index = index_disp2U + i + FORTRAN; //fortran format
	
	
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	          if(i<N_Active)
	    {	
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
	      k++;
	    }
	    
	   if(i<N_Active)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  index_dispU + KCol[j] + FORTRAN;        //fortran format
	      k++;
	    }
	   
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  index_disp2U + KCol[j] + FORTRAN;        //fortran format
	      k++;     
	   
	    }//for(j=
	
	  //Mat BT(m)
	 if(i<N_Active)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = index_disp3U + KCol_BT[j] + FORTRAN;       //fortran format
	    k++;
	  }//for(j=
    }
    
    
    for(i=0;i<NDof_P;i++)
    {
  
	row_index = index_disp3U + i + FORTRAN;                //fortran format
	pr_p = k;
	
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  	//Mat B2
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = index_dispU + KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  	//Mat B3
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = index_disp2U + KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
      if(i== 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	{
	  if(k>pr_p)
	  {
	    I_rn[k-1] = row_index;                //for project pressure
	    J_cn[k-1] = row_index;        //for project pressure
	  }
	  else
	  {
	    cout << "fehler in ParDirectSolver" << endl;
	  }
	}

    }//for(i=0;i<
    
     /** since the BT entries of Dirichlet U DOF row is not added */
    N_Nz = k ; 

  
  
}

void TSeqParDirectSolver::Solve(double *Sol, double *Rhs, bool Factorize)
{
  
  double t = MPI_Wtime();
  int i,j;
  int rank=0; 
   //printf("are going in here %d\n",DSType);
  if(DSType == 4)
  {
    GetRhs(Rhs);
    if(NSEType == 4)
    {
      if( rank==0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE == 0 )
	GlobalRhs[3*NDof_U] = 0;
    }
    
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
  
}


void TSeqParDirectSolver::AssembleMatrix()
{ 
  int i,j,k,l,m,t;
  double *EntriesA;
  int rank=0;
  
  if(MatB == NULL)
  {
    EntriesA = Mat->GetEntries();
    k = 0;
    for(i=0;i<NDof;i++)
    {

        for(j=RowPtr[i];j<RowPtr[i+1];j++)
        {
          MatLoc[k] = EntriesA[j];
	  k++;
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

void TSeqParDirectSolver::AssembleMatrix_NSE2()
{
  
}

void TSeqParDirectSolver::AssembleMatrix_NSE4()
{
  int N_Active;

  int i,j,k,l,m,t;

  int pr_p;
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
    }//for(i=0;i<
    
    
        for(i=0;i<NDof_U;i++)
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
    }//for(i=0;i<
    
    
        for(i=0;i<NDof_U;i++)
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
    }//for(i=0;i<
    
    
    
    
    for(i=0;i<NDof_P;i++)
    {
      	  pr_p = k;

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
	  
		
	  if(i == 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	  {
	    if(k>pr_p)
	    {
	      for(int ij = pr_p; ij<k ;ij++ )
		MatLoc[ij] = 0.0;
	      MatLoc[k-1] = 1.0;
	    }
	    else
	      cout << "fehler in ParDirectSolver" << endl;
	  }
    }//for(i=0;i<

}




#endif

#ifdef __2D__
TSeqParDirectSolver::TSeqParDirectSolver(int N_U,int N_P, int N_S, int N_D, TSquareMatrix2D **sqmat,TMatrix2D **recmat)
{
    int rank,size;

    Comm      = TDatabase::ParamDB->Comm;
    N_rhs     = 1;

  
    DSType    = TDatabase::ParamDB->DSType;

    NDof_P      = N_P;
    NDof_U      = N_U;
    NDof_S      = N_S;
    NDof_D      = N_D;
    NDof        = 2*NDof_U + NDof_P + 3*NDof_S + 3*NDof_D;   
    
    MatA     = sqmat[0];
    MatB    = recmat[0];
    MatBT   = recmat[2];
   
    SqmatrixA11 = sqmat[0];
    SqmatrixA12 = sqmat[1];
    SqmatrixA21 = sqmat[2];
    SqmatrixA22 = sqmat[3];
    
    if(NDof_S!=0)
    {
        MatC    = recmat[4];
        MatD    = recmat[8];
//         if(NDof_D!=0)
            MatG     = sqmat[4];
    }
    
    if(NDof_S!=0)
    {
        SqmatrixG11 = sqmat[4];
        SqmatrixG12 = sqmat[5];
        SqmatrixG21 = sqmat[6];
        SqmatrixG22 = sqmat[7];
        SqmatrixG23 = sqmat[8];
        SqmatrixG32 = sqmat[9];
        SqmatrixG33 = sqmat[10];
    }
    
    if(NDof_D!=0)
    {
        MatH = sqmat[11];
        SqmatrixH11 = sqmat[11];
        SqmatrixH22 = sqmat[12];
        SqmatrixH33 = sqmat[13];
    }
    
    MatrixB1    = recmat[0];
    MatrixB2    = recmat[1];
    MatrixB1T   = recmat[2];
    MatrixB2T   = recmat[3];
    
    if(NDof_S!=0)
    {
        MatrixC11   = recmat[4];   
        MatrixC12   = recmat[5];
        MatrixC22   = recmat[6];
        MatrixC23   = recmat[7];
        
        MatrixD11   = recmat[8];
        MatrixD12   = recmat[9];
        MatrixD21   = recmat[10];
        MatrixD22   = recmat[11];
        MatrixD31   = recmat[12];
        MatrixD32   = recmat[13];
    }
    
    if(NDof_D!=0)
    {
        MatE    = recmat[14];
        MatJ    = recmat[18];
        MatrixE11  = recmat[14];
        MatrixE12  = recmat[15];
        MatrixE22  = recmat[16];
        MatrixE23  = recmat[17];
        
        MatrixJ11  = recmat[18];
        MatrixJ21  = recmat[19];
        MatrixJ22  = recmat[20];
        MatrixJ32  = recmat[21];
    }
	     
  
    if(DSType == 4)
    {
        if(NDof_D==0 && NDof_S !=0) 
        {
            InitMumps_NSEType4_CST_2D();
        }
        else if(NDof_D!=0)
        {
            InitMumps_NSEType4_CST_DEVSS_2D();
        }
        else if(NDof_D==0 && NDof_S ==0) // Only NSE
        {
            InitMumps_NSEType4();
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

void TSeqParDirectSolver::InitMumps_NSEType4()
{
  
    int i,j,k,l,m,t,N_Active_U,  row_index;

    int rank=0,size;
    int pr_p;  
    
    N_Active_U = MatA->GetActiveBound();

    RowPtr = MatA->GetRowPtr();
    KCol   = MatA->GetKCol();
    
    RowPtr_B = MatB->GetRowPtr();
    KCol_B   = MatB->GetKCol();

    RowPtr_BT = MatBT->GetRowPtr();
    KCol_BT   = MatBT->GetKCol();
    
    N_Nz_U = MatA->GetN_Entries();
    N_Nz_P = MatB->GetN_Entries();

    N_Nz = 4*N_Nz_U + 4*N_Nz_P ; //total non zeros in system matrix 

    N_Eqns = NDof;
    GlobalRhsSize = NDof;
    
    GlobalRhs = new double[GlobalRhsSize];
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
    k = 0;
    for(i=0;i<NDof_U;i++)
    {      // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      I_rn[k] =  i + FORTRAN;              //fortran format
	      J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
	      k++;
	      
		  if(i<N_Active_U)
		  {	
		    I_rn[k] =  i + FORTRAN;             //fortran format
		    J_cn[k] =  NDof_U + KCol[j] + FORTRAN;           //fortran format
		    k++;
		  }
	    }//for(j=
	
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  i + FORTRAN;             //fortran format
	    J_cn[k] = 2*NDof_U + KCol_BT[j] + FORTRAN;      //fortran format
	    k++;
	  }//for(j=


    }//for(i=0;i<
    
    for(i=0;i<NDof_U;i++)
    {
      
	row_index = NDof_U + i + FORTRAN; //fortran format
	     // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
		if(i<N_Active_U)
		{	
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
		  k++;
		}
	    
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	 
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = 2*NDof_U +  KCol_BT[j] + FORTRAN;       //fortran format
	    k++;
	  }//for(j=

    }//for(i=0;i<
    
  
    for(i=0;i<NDof_P;i++)
    {
  
	row_index = 2*NDof_U  + i + FORTRAN;                //fortran format
	pr_p = k;
	
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  	//Mat B2
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = NDof_U + KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  
      if(i== 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	{
	  if(k>pr_p)
	  {
	    I_rn[k-1] = row_index;                //for project pressure
	    J_cn[k-1] = row_index;        //for project pressure
	  }
	  else
	  {
	    cout << "fehler in ParDirectSolver" << endl;
	  }
	}

    }//for(i=0;i<
  
       /** since the BT entries of Dirichlet U DOF row is not added */
    N_Nz = k ; 
  
}


void TSeqParDirectSolver::InitMumps_NSEType4_CST_2D()
{
  
    int i,j,k,l,m,t,N_Active_U, N_Active_S, row_index;

  int rank=0,size;
  int pr_p;
  
    
    N_Active_U = MatA->GetActiveBound();
    N_Active_S = MatG->GetActiveBound();
    

    RowPtr = MatA->GetRowPtr();
    KCol   = MatA->GetKCol();
    
    RowPtr_G = MatG->GetRowPtr();
    KCol_G   = MatG->GetKCol();
    
    RowPtr_B = MatB->GetRowPtr();
    KCol_B   = MatB->GetKCol();

    RowPtr_BT = MatBT->GetRowPtr();
    KCol_BT   = MatBT->GetKCol();
    
    RowPtr_C = MatC->GetRowPtr();
    KCol_C   = MatC->GetKCol();
    
    RowPtr_D = MatD->GetRowPtr();
    KCol_D   = MatD->GetKCol();
    
    N_Nz_U = MatA->GetN_Entries();
    N_Nz_G = MatG->GetN_Entries();
    N_Nz_P = MatB->GetN_Entries();
    N_Nz_C = MatC->GetN_Entries();
    N_Nz_D = MatD->GetN_Entries();
    
    N_Nz = 4*N_Nz_U + 7*N_Nz_G + 4*N_Nz_P + 4*N_Nz_C + 6*N_Nz_D; //total non zeros in system matrix 

    N_Eqns = NDof;
    GlobalRhsSize = NDof;
    

      GlobalRhs = new double[GlobalRhsSize];
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
  
    
     k = 0;
    for(i=0;i<NDof_U;i++)
    {      // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      I_rn[k] =  i + FORTRAN;              //fortran format
	      J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
	      k++;
	      
		  if(i<N_Active_U)
		  {	
		    I_rn[k] =  i + FORTRAN;             //fortran format
		    J_cn[k] =  NDof_U + KCol[j] + FORTRAN;           //fortran format
		    k++;
		  }
	    }//for(j=
	
	// Mat C
	 if(i<N_Active_U)
	 {
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	      I_rn[k] =  i + FORTRAN;              //fortran format
	      J_cn[k] =  2*NDof_U + KCol_C[j] + FORTRAN;        //fortran format
	      k++;
	      
              I_rn[k] =  i + FORTRAN;             //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_C[j] + FORTRAN;           //fortran format
              k++;

	    }//for(j=
	   
	 }	
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  i + FORTRAN;             //fortran format
	    J_cn[k] = 2*NDof_U + 3*NDof_S + KCol_BT[j] + FORTRAN;      //fortran format
	    k++;
	  }//for(j=


    }//for(i=0;i<
    
  for(i=0;i<NDof_U;i++)
    {
      
	row_index = NDof_U + i + FORTRAN; //fortran format
	     // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
		if(i<N_Active_U)
		{	
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
		  k++;
		}
	    
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	
	// Mat C
	 if(i<N_Active_U)
	 {
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	      I_rn[k] =  row_index;             //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_C[j] + FORTRAN;        //fortran format
	      k++;
	      
              I_rn[k] =  row_index;          //fortran format
	      J_cn[k] =  2*NDof_U + 2*NDof_S + KCol_C[j] + FORTRAN;           //fortran format
              k++;

	    }//for(j=
	   
	 }	
	 
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = 2*NDof_U + 3*NDof_S + KCol_BT[j] + FORTRAN;       //fortran format
	    k++;
	  }//for(j=

    }//for(i=0;i<
    
    for(i=0;i<NDof_S;i++)
    {
      
	row_index = 2*NDof_U + i + FORTRAN; //fortran format
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_D[j] + FORTRAN;        //fortran format
		  k++;
 
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol_D[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  2*NDof_U + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      if(i<N_Active_S)
	      {
              I_rn[k] =  row_index ;          //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_G[j] + FORTRAN;           //fortran format
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
    
      for(i=0;i<NDof_S;i++)
    {
      
	row_index = 2*NDof_U + NDof_S + i + FORTRAN; //fortran format
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_D[j] + FORTRAN;        //fortran format
		  k++;
 
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol_D[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  2*NDof_U + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      }
	      
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      
	      if(i<N_Active_S)
	      {
              I_rn[k] =  row_index;           //fortran format
	      J_cn[k] =  2*NDof_U + 2*NDof_S + KCol_G[j] + FORTRAN;           //fortran format
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
  
        for(i=0;i<NDof_S;i++)
    {
      
	row_index = 2*NDof_U + 2*NDof_S + i + FORTRAN; //fortran format
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_D[j] + FORTRAN;        //fortran format
		  k++;
 
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol_D[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	      I_rn[k] =  row_index ;             //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      }
	      
	      I_rn[k] =  row_index ;            //fortran format
	      J_cn[k] =  2*NDof_U + 2*NDof_S + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
      
	    }//for(j=

    }//for(i=0;i<
  
  for(i=0;i<NDof_P;i++)
    {
  
	row_index = 2*NDof_U + 3*NDof_S + i + FORTRAN;                //fortran format
	pr_p = k;
	
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  	//Mat B2
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = NDof_U + KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  
      if(i== 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	{
	  if(k>pr_p)
	  {
	    I_rn[k-1] = row_index;                //for project pressure
	    J_cn[k-1] = row_index;        //for project pressure
	  }
	  else
	  {
	    cout << "fehler in ParDirectSolver" << endl;
	  }
	}

    }//for(i=0;i<
  
       /** since the BT entries of Dirichlet U DOF row is not added */
    N_Nz = k ; 
  
}

void TSeqParDirectSolver::InitMumps_NSEType4_CST_DEVSS_2D()
{
        int i,j,k,l,m,t,N_Active_U, N_Active_S, N_Active_D, row_index;

  int rank=0,size;
  int pr_p;
  
    
    N_Active_U = MatA->GetActiveBound();
    N_Active_S = MatG->GetActiveBound();
    N_Active_D = MatH->GetActiveBound();

    RowPtr = MatA->GetRowPtr();
    KCol   = MatA->GetKCol();
    
    RowPtr_G = MatG->GetRowPtr();
    KCol_G   = MatG->GetKCol();
    
    RowPtr_H = MatH->GetRowPtr();
    KCol_H   = MatH->GetKCol();
    
    RowPtr_B = MatB->GetRowPtr();
    KCol_B   = MatB->GetKCol();

    RowPtr_BT = MatBT->GetRowPtr();
    KCol_BT   = MatBT->GetKCol();
    
    RowPtr_C = MatC->GetRowPtr();
    KCol_C   = MatC->GetKCol();
    
    RowPtr_D = MatD->GetRowPtr();
    KCol_D   = MatD->GetKCol();
    
    RowPtr_E = MatE->GetRowPtr();
    KCol_E   = MatE->GetKCol();
    
    RowPtr_J = MatJ->GetRowPtr();
    KCol_J   = MatJ->GetKCol();
    
    N_Nz_U = MatA->GetN_Entries();
    N_Nz_G = MatG->GetN_Entries();
    N_Nz_H = MatH->GetN_Entries();
    N_Nz_P = MatB->GetN_Entries();
    N_Nz_C = MatC->GetN_Entries();
    N_Nz_D = MatD->GetN_Entries();
    N_Nz_E = MatE->GetN_Entries();
    N_Nz_J = MatJ->GetN_Entries();
    
    N_Nz = 4*N_Nz_U + 7*N_Nz_G + 4*N_Nz_P + 4*N_Nz_C + 6*N_Nz_D + 3*N_Nz_H + 4*N_Nz_E + 4*N_Nz_J; //total non zeros in system matrix 

    N_Eqns = NDof;
    GlobalRhsSize = NDof;
    

      GlobalRhs = new double[GlobalRhsSize];
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
         k = 0;
    for(i=0;i<NDof_U;i++)
    {      // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      I_rn[k] =  i + FORTRAN;              //fortran format
	      J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
	      k++;
	      
		  if(i<N_Active_U)
		  {	
		    I_rn[k] =  i + FORTRAN;             //fortran format
		    J_cn[k] =  NDof_U + KCol[j] + FORTRAN;           //fortran format
		    k++;
		  }
	    }//for(j=
	
	
	 if(i<N_Active_U)
	 {// Mat C
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	      I_rn[k] =  i + FORTRAN;              //fortran format
	      J_cn[k] =  2*NDof_U + KCol_C[j] + FORTRAN;        //fortran format
	      k++;
	      
              I_rn[k] =  i + FORTRAN;             //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_C[j] + FORTRAN;           //fortran format
              k++;

	    }//for(j=
	    
	    // Mat E
	    for(j=RowPtr_E[i];j<RowPtr_E[i+1];j++)
	    {
	      I_rn[k] =  i + FORTRAN;              //fortran format
	      J_cn[k] =  2*NDof_U + 3*NDof_S + KCol_E[j] + FORTRAN;        //fortran format
	      k++;
	      
              I_rn[k] =  i + FORTRAN;             //fortran format
	      J_cn[k] =  2*NDof_U + 3*NDof_S + NDof_D + KCol_E[j] + FORTRAN;           //fortran format
              k++;

	    }//for(j=
	    
	   
	 }	
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  i + FORTRAN;             //fortran format
	    J_cn[k] = 2*NDof_U + 3*NDof_S + 3*NDof_D + KCol_BT[j] + FORTRAN;      //fortran format
	    k++;
	  }//for(j=

    }//for(i=0;i<
    
    for(i=0;i<NDof_U;i++)
    {
      
	row_index = NDof_U + i + FORTRAN; //fortran format
	     // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
		if(i<N_Active_U)
		{	
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol[j] + FORTRAN;        //fortran format
		  k++;
		}
	    
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	
	
	 if(i<N_Active_U)
	 {// Mat C
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	      I_rn[k] =  row_index;             //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_C[j] + FORTRAN;        //fortran format
	      k++;
	      
              I_rn[k] =  row_index;          //fortran format
	      J_cn[k] =  2*NDof_U + 2*NDof_S + KCol_C[j] + FORTRAN;           //fortran format
              k++;

	    }//for(j=
	   
	  // Mat E
	    for(j=RowPtr_E[i];j<RowPtr_E[i+1];j++)
	    {
	      I_rn[k] =  row_index ;            //fortran format
	      J_cn[k] =  2*NDof_U + 3*NDof_S + NDof_D + KCol_E[j] + FORTRAN;        //fortran format
	      k++;
	      
              I_rn[k] = row_index  ;          //fortran format
	      J_cn[k] =  2*NDof_U + 3*NDof_S + 2*NDof_D + KCol_E[j] + FORTRAN;           //fortran format
              k++;

	    }//for(j=
	   
	 }	
	 
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    I_rn[k] =  row_index;                //fortran format
	    J_cn[k] = 2*NDof_U + 3*NDof_S + 3*NDof_D + KCol_BT[j] + FORTRAN;       //fortran format
	    k++;
	  }//for(j=

    }//for(i=0;i<
    
    
    for(i=0;i<NDof_S;i++)
    {
      
	row_index = 2*NDof_U + i + FORTRAN; //fortran format
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_D[j] + FORTRAN;        //fortran format
		  k++;
 
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol_D[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  2*NDof_U + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      if(i<N_Active_S)
	      {
              I_rn[k] =  row_index ;          //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_G[j] + FORTRAN;           //fortran format
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
    
     for(i=0;i<NDof_S;i++)
    {
      
	row_index = 2*NDof_U + NDof_S + i + FORTRAN; //fortran format
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_D[j] + FORTRAN;        //fortran format
		  k++;
 
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol_D[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  2*NDof_U + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      }
	      
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      
	      if(i<N_Active_S)
	      {
              I_rn[k] =  row_index;           //fortran format
	      J_cn[k] =  2*NDof_U + 2*NDof_S + KCol_G[j] + FORTRAN;           //fortran format
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
    
    
    for(i=0;i<NDof_S;i++)
    {
      
	row_index = 2*NDof_U + 2*NDof_S + i + FORTRAN; //fortran format
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_D[j] + FORTRAN;        //fortran format
		  k++;
 
	      I_rn[k] =  row_index;              //fortran format
	      J_cn[k] =  NDof_U + KCol_D[j] + FORTRAN;        //fortran format
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	      I_rn[k] =  row_index ;             //fortran format
	      J_cn[k] =  2*NDof_U + NDof_S + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
	      }
	      
	      I_rn[k] =  row_index ;            //fortran format
	      J_cn[k] =  2*NDof_U + 2*NDof_S + KCol_G[j] + FORTRAN;        //fortran format
	      k++;
      
	    }//for(j=

    }//for(i=0;i<
    
    
      for(i=0;i<NDof_D;i++)
    {
      
	row_index = 2*NDof_U + 3*NDof_S + i + FORTRAN; //fortran format
	     // Mat J
	     if(i<N_Active_D)
	     {
	    for(j=RowPtr_J[i];j<RowPtr_J[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_J[j] + FORTRAN;        //fortran format
		  k++;
		  
	    }//for(j=
	   } // if(i<N_Active_D)
	   
	// Mat H
	   for(j=RowPtr_H[i];j<RowPtr_H[i+1];j++)
	    { 
	      I_rn[k] =  row_index ;            //fortran format
	      J_cn[k] =  2*NDof_U + 3*NDof_S + KCol_H[j] + FORTRAN;        //fortran format
	      k++;
      
	    }//for(j=

    }//for(i=0;i<
    
    for(i=0;i<NDof_D;i++)
    {
      
	row_index = 2*NDof_U + 3*NDof_S + NDof_D  + i + FORTRAN; //fortran format
	     // Mat J
	     if(i<N_Active_D)
	     {
	    for(j=RowPtr_J[i];j<RowPtr_J[i+1];j++)
	    {
	          I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  KCol_J[j] + FORTRAN;        //fortran format
		  k++;
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  NDof_U + KCol_J[j] + FORTRAN;        //fortran format
		  k++;
		  
	    }//for(j=
	   } // if(i<N_Active_D)
	   
	// Mat H
	   for(j=RowPtr_H[i];j<RowPtr_H[i+1];j++)
	    { 
	      I_rn[k] =  row_index ;            //fortran format
	      J_cn[k] =  2*NDof_U + 3*NDof_S + NDof_D + KCol_H[j] + FORTRAN;        //fortran format
	      k++;
      
	    }//for(j=

    }//for(i=0;i<
    
          for(i=0;i<NDof_D;i++)
    {
      
	row_index = 2*NDof_U + 3*NDof_S + 2*NDof_D + i + FORTRAN; //fortran format
	     // Mat J
	     if(i<N_Active_D)
	     {
	    for(j=RowPtr_J[i];j<RowPtr_J[i+1];j++)
	    {
			
		  I_rn[k] =  row_index;              //fortran format
		  J_cn[k] =  NDof_U + KCol_J[j] + FORTRAN;        //fortran format
		  k++;
		  
	    }//for(j=
	   } // if(i<N_Active_D)
	   
	// Mat H
	   for(j=RowPtr_H[i];j<RowPtr_H[i+1];j++)
	    { 
	      I_rn[k] =  row_index ;            //fortran format
	      J_cn[k] =  2*NDof_U + 3*NDof_S + 2*NDof_D + KCol_H[j] + FORTRAN;        //fortran format
	      k++;
      
	    }//for(j=

    }//for(i=0;i<
  
  
  for(i=0;i<NDof_P;i++)
    {
  
	row_index = 2*NDof_U + 3*NDof_S + 3*NDof_D + i + FORTRAN;                //fortran format
	pr_p = k;
	
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  	//Mat B2
	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    I_rn[k] = row_index;                //fortran format
	    J_cn[k] = NDof_U + KCol_B[j] + FORTRAN;        //fortran format
	    k++;
	  }//for(j=
	  
	  
      if(i== 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	{
	  if(k>pr_p)
	  {
	    I_rn[k-1] = row_index;                //for project pressure
	    J_cn[k-1] = row_index;        //for project pressure
	  }
	  else
	  {
	    cout << "fehler in ParDirectSolver" << endl;
	  }
	}

    }//for(i=0;i<
  
     /** since the BT entries of Dirichlet U DOF row is not added */
    N_Nz = k ; 
  
}


void TSeqParDirectSolver::AssembleMatrix()
{ 

    if(NDof_D != 0 && NDof_S != 0)
    {    
        AssembleMatrix_NSEType4_CST_DEVSS2D();
    }
    else if(NDof_D == 0 && NDof_S != 0)
    {
        AssembleMatrix_NSEType4_CST_2D();
    }
    else if(NDof_D == 0 && NDof_S == 0)
    {
        AssembleMatrix_NSEType4D();
    }
}

void TSeqParDirectSolver::AssembleMatrix_NSEType4D()
{
  int N_Active_U;

  int i,j,k,l,m,t;

  int pr_p;
  double* EntriesA[4];
  double* EntriesB[4];

  
  EntriesA[0] = SqmatrixA11->GetEntries();
  EntriesA[1] = SqmatrixA12->GetEntries();
  EntriesA[2] = SqmatrixA21->GetEntries();
  EntriesA[3] = SqmatrixA22->GetEntries();
  
  EntriesB[0] = MatrixB1T->GetEntries();
  EntriesB[1] = MatrixB2T->GetEntries();
  EntriesB[2] = MatrixB1->GetEntries();
  EntriesB[3] = MatrixB2->GetEntries();

  N_Active_U = MatA->GetActiveBound();

  
     k = 0;
    for(i=0;i<NDof_U;i++)
    {       // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      MatLoc[k] = EntriesA[0][j]; 
	      k++;
	      
	      if(i<N_Active_U)
	      {
	       MatLoc[k] = EntriesA[1][j];
	      k++;
	      }      
	    }

	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[0][j];
	    k++;
	  }//for(j=
    }//for(i=0;i<
    

        for(i=0;i<NDof_U;i++)
    {       // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      if(i<N_Active_U)
	      {
	       MatLoc[k] = EntriesA[2][j];
	       k++;
	      }
	        MatLoc[k] = EntriesA[3][j];
	        k++;      
	    }
	    
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[1][j];     
	    k++;
	  }//for(j=
    }//for(i=0;i<
    
   
    
    for(i=0;i<NDof_P;i++)
    {
      	  pr_p = k;

	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[2][j];
	    k++;
	  }

	   for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[3][j];
	    k++;
	  } 
		
	  if(i == 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	  {
	    if(k>pr_p)
	    {
	      for(int ij = pr_p; ij<k ;ij++ )
		MatLoc[ij] = 0.0;
	      MatLoc[k-1] = 1.0;
	    }
	    else
	      cout << "fehler in ParDirectSolver" << endl;
	  }
    }//for(i=0;i<

}


void TSeqParDirectSolver::AssembleMatrix_NSEType4_CST_2D()
{
  int N_Active_U, N_Active_S;

  int i,j,k,l,m,t;

  int pr_p;
  double* EntriesA[4];
  double* EntriesB[4];
  double* EntriesG[7];
  double* EntriesC[4];
  double* EntriesD[6];
  
  EntriesA[0] = SqmatrixA11->GetEntries();
  EntriesA[1] = SqmatrixA12->GetEntries();
  EntriesA[2] = SqmatrixA21->GetEntries();
  EntriesA[3] = SqmatrixA22->GetEntries();
  
  EntriesB[0] = MatrixB1T->GetEntries();
  EntriesB[1] = MatrixB2T->GetEntries();
  EntriesB[2] = MatrixB1->GetEntries();
  EntriesB[3] = MatrixB2->GetEntries();

  EntriesC[0] = MatrixC11->GetEntries();
  EntriesC[1] = MatrixC12->GetEntries();
  EntriesC[2] = MatrixC22->GetEntries();
  EntriesC[3] = MatrixC23->GetEntries();
  
  EntriesD[0] = MatrixD11->GetEntries();
  EntriesD[1] = MatrixD12->GetEntries();
  EntriesD[2] = MatrixD21->GetEntries();
  EntriesD[3] = MatrixD22->GetEntries();
  EntriesD[4] = MatrixD31->GetEntries();
  EntriesD[5] = MatrixD32->GetEntries();
  
  EntriesG[0] = SqmatrixG11->GetEntries();
  EntriesG[1] = SqmatrixG12->GetEntries();
  EntriesG[2] = SqmatrixG21->GetEntries();
  EntriesG[3] = SqmatrixG22->GetEntries();
  EntriesG[4] = SqmatrixG23->GetEntries();
  EntriesG[5] = SqmatrixG32->GetEntries();
  EntriesG[6] = SqmatrixG33->GetEntries();
  
  N_Active_U = MatA->GetActiveBound();
  N_Active_S = MatG->GetActiveBound();
  
     k = 0;
    for(i=0;i<NDof_U;i++)
    {       // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      MatLoc[k] = EntriesA[0][j]; 
	      k++;
	      
	      if(i<N_Active_U)
	      {
	       MatLoc[k] = EntriesA[1][j];
	      k++;
	      }      
	    }
          // Mat C
	 if(i<N_Active_U)
	 {
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	      MatLoc[k] = EntriesC[0][j]; 
	      k++;
	      MatLoc[k] = EntriesC[1][j]; 
              k++;

	    }//for(j=
	 }	

	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[0][j];
	    k++;
	  }//for(j=
    }//for(i=0;i<
    

        for(i=0;i<NDof_U;i++)
    {       // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      if(i<N_Active_U)
	      {
	       MatLoc[k] = EntriesA[2][j];
	       k++;
	      }
	        MatLoc[k] = EntriesA[3][j];
	        k++;      
	    }
          // Mat C
	 if(i<N_Active_U)
	 {
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	   
	      MatLoc[k] = EntriesC[2][j]; 
	      k++;
	      MatLoc[k] = EntriesC[3][j]; 
              k++;

	    }//for(j=
	   
	 }	
	    
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[1][j];     
	    k++;
	  }//for(j=
    }//for(i=0;i<
    
     for(i=0;i<NDof_S;i++)
    {
      
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
               MatLoc[k] = EntriesD[0][j]; 
	       k++;
               MatLoc[k] = EntriesD[1][j]; 
	       k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
              MatLoc[k] = EntriesG[0][j]; 
	      k++;
	      if(i<N_Active_S)
	      {
              MatLoc[k] = EntriesG[1][j]; 
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
    
    
       for(i=0;i<NDof_S;i++)
    {
      
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
                  MatLoc[k] = EntriesD[2][j]; 
		  k++;
                  MatLoc[k] = EntriesD[3][j]; 
	          k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	       MatLoc[k] = EntriesG[2][j]; 
	      k++;
	      }
	      
	       MatLoc[k] = EntriesG[3][j]; 
	      k++;
	      
	      if(i<N_Active_S)
	      {
               MatLoc[k] = EntriesG[4][j]; 
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
    
     for(i=0;i<NDof_S;i++)
    {
      
	
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		 MatLoc[k] = EntriesD[4][j]; 
		  k++;
 
	      MatLoc[k] = EntriesD[5][j]; 
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	      MatLoc[k] = EntriesG[5][j]; 
	      k++;
	      }
	      
	     MatLoc[k] = EntriesG[6][j]; 
	      k++;
      
	    }//for(j=

    }//for(i=0;i<
    
    
    for(i=0;i<NDof_P;i++)
    {
      	  pr_p = k;

	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[2][j];
	    k++;
	  }

	   for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[3][j];
	    k++;
	  } 
		
	  if(i == 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	  {
	    if(k>pr_p)
	    {
	      for(int ij = pr_p; ij<k ;ij++ )
		MatLoc[ij] = 0.0;
	      MatLoc[k-1] = 1.0;
	    }
	    else
	      cout << "fehler in ParDirectSolver" << endl;
	  }
    }//for(i=0;i<

}


void TSeqParDirectSolver::AssembleMatrix_NSEType4_CST_DEVSS2D()
{
  int N_Active_U, N_Active_S, N_Active_D;

  int i,j,k,l,m,t;

  int pr_p;
  double* EntriesA[4];
  double* EntriesB[4];
  double* EntriesG[7];
  double* EntriesC[4];
  double* EntriesD[6];
  double* EntriesE[4];
  double* EntriesJ[4];
  double* EntriesH[3];
  
  EntriesA[0] = SqmatrixA11->GetEntries();
  EntriesA[1] = SqmatrixA12->GetEntries();
  EntriesA[2] = SqmatrixA21->GetEntries();
  EntriesA[3] = SqmatrixA22->GetEntries();
  
  EntriesB[0] = MatrixB1T->GetEntries();
  EntriesB[1] = MatrixB2T->GetEntries();
  EntriesB[2] = MatrixB1->GetEntries();
  EntriesB[3] = MatrixB2->GetEntries();

  EntriesC[0] = MatrixC11->GetEntries();
  EntriesC[1] = MatrixC12->GetEntries();
  EntriesC[2] = MatrixC22->GetEntries();
  EntriesC[3] = MatrixC23->GetEntries();
  
  EntriesD[0] = MatrixD11->GetEntries();
  EntriesD[1] = MatrixD12->GetEntries();
  EntriesD[2] = MatrixD21->GetEntries();
  EntriesD[3] = MatrixD22->GetEntries();
  EntriesD[4] = MatrixD31->GetEntries();
  EntriesD[5] = MatrixD32->GetEntries();
  
  EntriesG[0] = SqmatrixG11->GetEntries();
  EntriesG[1] = SqmatrixG12->GetEntries();
  EntriesG[2] = SqmatrixG21->GetEntries();
  EntriesG[3] = SqmatrixG22->GetEntries();
  EntriesG[4] = SqmatrixG23->GetEntries();
  EntriesG[5] = SqmatrixG32->GetEntries();
  EntriesG[6] = SqmatrixG33->GetEntries();
  
  EntriesH[0] = SqmatrixH11->GetEntries();
  EntriesH[1] = SqmatrixH22->GetEntries();
  EntriesH[2] = SqmatrixH33->GetEntries();
  
  EntriesE[0] = MatrixE11->GetEntries();
  EntriesE[1] = MatrixE12->GetEntries();
  EntriesE[2] = MatrixE22->GetEntries();
  EntriesE[3] = MatrixE23->GetEntries();
  
  EntriesJ[0] = MatrixJ11->GetEntries();
  EntriesJ[1] = MatrixJ21->GetEntries();
  EntriesJ[2] = MatrixJ22->GetEntries();
  EntriesJ[3] = MatrixJ32->GetEntries();
  
  N_Active_U = MatA->GetActiveBound();
  N_Active_S = MatG->GetActiveBound();
  N_Active_D = MatH->GetActiveBound();
  
   k = 0;
    for(i=0;i<NDof_U;i++)
    {       // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      MatLoc[k] = EntriesA[0][j]; 
	      k++;
	      
	      if(i<N_Active_U)
	      {
	       MatLoc[k] = EntriesA[1][j];
	      k++;
	      }      
	    }
          
	 if(i<N_Active_U)
	 {// Mat C
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	      MatLoc[k] = EntriesC[0][j]; 
	      k++;
	      MatLoc[k] = EntriesC[1][j]; 
              k++;

	    }//for(j=
	    
	    // Mat E
	      for(j=RowPtr_E[i];j<RowPtr_E[i+1];j++)
	    {
	      MatLoc[k] = EntriesE[0][j]; 
	      k++;
	      MatLoc[k] = EntriesE[1][j]; 
              k++;

	    }//for(j=
	        
	 }	
	 
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[0][j];
	    k++;
	  }//for(j=
    }//for(i=0;i<
  
   for(i=0;i<NDof_U;i++)
    {       // Mat A
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      if(i<N_Active_U)
	      {
	       MatLoc[k] = EntriesA[2][j];
	       k++;
	      }
	        MatLoc[k] = EntriesA[3][j];
	        k++;      
	    }
          
	 if(i<N_Active_U)
	 { // Mat C
	   for(j=RowPtr_C[i];j<RowPtr_C[i+1];j++)
	    {
	   
	      MatLoc[k] = EntriesC[2][j]; 
	      k++;
	      MatLoc[k] = EntriesC[3][j]; 
              k++;

	    }//for(j=
	   
	   // Mat E
	      for(j=RowPtr_E[i];j<RowPtr_E[i+1];j++)
	    {
	      MatLoc[k] = EntriesE[2][j]; 
	      k++;
	      MatLoc[k] = EntriesE[3][j]; 
              k++;

	    }//for(j=
	   
	 }	
	    
	  //Mat BT(m)
	 if(i<N_Active_U)
	  for(j=RowPtr_BT[i];j<RowPtr_BT[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[1][j];     
	    k++;
	  }//for(j=
    }//for(i=0;i<
  
   for(i=0;i<NDof_S;i++)
    {
      
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
               MatLoc[k] = EntriesD[0][j]; 
	       k++;
               MatLoc[k] = EntriesD[1][j]; 
	       k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
              MatLoc[k] = EntriesG[0][j]; 
	      k++;
	      if(i<N_Active_S)
	      {
              MatLoc[k] = EntriesG[1][j]; 
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
    
    
       for(i=0;i<NDof_S;i++)
    {
      
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
                  MatLoc[k] = EntriesD[2][j]; 
		  k++;
                  MatLoc[k] = EntriesD[3][j]; 
	          k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	       MatLoc[k] = EntriesG[2][j]; 
	      k++;
	      }
	      
	       MatLoc[k] = EntriesG[3][j]; 
	      k++;
	      
	      if(i<N_Active_S)
	      {
               MatLoc[k] = EntriesG[4][j]; 
              k++;
	      }
	      
	    }//for(j=

    }//for(i=0;i<
    
     for(i=0;i<NDof_S;i++)
    {
      
	
	     // Mat D
	     if(i<N_Active_S)
	     {
	    for(j=RowPtr_D[i];j<RowPtr_D[i+1];j++)
	    {
			
		 MatLoc[k] = EntriesD[4][j]; 
		  k++;
 
	      MatLoc[k] = EntriesD[5][j]; 
	      k++;
	    }//for(j=
	   } // if(i<N_Active_S)
	   
	// Mat G
	   for(j=RowPtr_G[i];j<RowPtr_G[i+1];j++)
	    {
	      if(i<N_Active_S)
	      {
	      MatLoc[k] = EntriesG[5][j]; 
	      k++;
	      }
	      
	     MatLoc[k] = EntriesG[6][j]; 
	      k++;
      
	    }//for(j=

    }//for(i=0;i<
    
    
    for(i=0;i<NDof_D;i++)
    {
      
	     // Mat J
	     if(i<N_Active_D)
	     {
	    for(j=RowPtr_J[i];j<RowPtr_J[i+1];j++)
	    {
			
		   MatLoc[k] = EntriesJ[0][j]; 
		  k++;
		  
	    }//for(j=
	   } // if(i<N_Active_D)
	   
	// Mat H
	   for(j=RowPtr_H[i];j<RowPtr_H[i+1];j++)
	    { 
	       MatLoc[k] = EntriesH[0][j]; 
		  k++;
      
	    }//for(j=

    }//for(i=0;i<
    
     for(i=0;i<NDof_D;i++)
    {
     
	     // Mat J
	     if(i<N_Active_D)
	     {
	    for(j=RowPtr_J[i];j<RowPtr_J[i+1];j++)
	    {
			
		   MatLoc[k] = EntriesJ[1][j]; 
		  k++;
		  
		  MatLoc[k] = EntriesJ[2][j]; 
		  k++;
		  
	    }//for(j=
	   } // if(i<N_Active_D)
	   
	// Mat H
	   for(j=RowPtr_H[i];j<RowPtr_H[i+1];j++)
	    { 
	       MatLoc[k] = EntriesH[1][j]; 
		  k++;
      
	    }//for(j=

    }//for(i=0;i<
    
         for(i=0;i<NDof_D;i++)
    {
     
	     // Mat J
	     if(i<N_Active_D)
	     {
	    for(j=RowPtr_J[i];j<RowPtr_J[i+1];j++)
	    {
			
		   MatLoc[k] = EntriesJ[3][j]; 
		  k++;
		  
	    }//for(j=
	   } // if(i<N_Active_D)
	   
	// Mat H
	   for(j=RowPtr_H[i];j<RowPtr_H[i+1];j++)
	    { 
	       MatLoc[k] = EntriesH[2][j]; 
		  k++;
      
	    }//for(j=

    }//for(i=0;i<
    
    
        for(i=0;i<NDof_P;i++)
    {
      	  pr_p = k;

	  for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[2][j];
	    k++;
	  }

	   for(j=RowPtr_B[i];j<RowPtr_B[i+1];j++)
	  {
	    MatLoc[k] = EntriesB[3][j];
	    k++;
	  } 
		
	  if(i == 0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
	  {
	    if(k>pr_p)
	    {
	      for(int ij = pr_p; ij<k ;ij++ )
		MatLoc[ij] = 0.0;
	      MatLoc[k-1] = 1.0;
	    }
	    else
	      cout << "fehler in ParDirectSolver" << endl;
	  }
    }//for(i=0;i<
  
  
}


void TSeqParDirectSolver::Solve(double *Sol, double *Rhs, bool Factorize)
{
  
  double t = MPI_Wtime();
  int i,j;
  int rank=0; 
   //printf("are going in here %d\n",DSType);
  if(DSType == 4)
  {
    GetRhs(Rhs);
           
       
      if( rank==0 && TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE != 0 )
      { if(NDof_D==0 && NDof_S!=0)
	GlobalRhs[2*NDof_U+3*NDof_S] = 1.0;
      else if(NDof_D!=0 && NDof_S!=0)
	GlobalRhs[2*NDof_U+3*NDof_S+3*NDof_D] = 1.0;
      else if(NDof_D==0 && NDof_S==0)
	GlobalRhs[2*NDof_U] = 1.0;
      }
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
  
}




#endif





TSeqParDirectSolver::~TSeqParDirectSolver()
{
  if(DSType == 4)
    Mumps->Clean();

  delete [] MatLoc;
  delete [] OwnRhs;
  delete [] I_rn;
  delete [] J_cn;
  delete [] GlobalRhs;
//   delete Mumps;
}


void TSeqParDirectSolver::GetRhs(double *Rhs)
{

  for(int i=0;i<NDof;i++)
    GlobalRhs[i] = Rhs[i];
}

void TSeqParDirectSolver::UpdateSol(double *Sol)
{
 for(int i=0;i<NDof;i++)
    Sol[i] = GlobalRhs[i];
}



#endif









   
