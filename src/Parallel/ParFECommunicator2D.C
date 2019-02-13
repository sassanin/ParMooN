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
// @(#)ParFECommunicator2D.C
//
// Class:      TParFECommunicator2D
// Purpose:    Class containing all info needed for communication between subdomains
//
// Author:     Sashikumaar Ganesan (01.10.09)
//
// History:    Start of implementation 01.10.09 (Sashikumaar Ganesan)
//
// ======================================================================= 
#ifdef _MPI

#  include "mpi.h"

#include <ParFECommunicator2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <MeshPartition.h>
#include <Database.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <Joint.h>
#include <SubDomainJoint.h>

TParFECommunicator2D::TParFECommunicator2D(MPI_Comm comm, TFESpace2D *FESpace)
{
 int N_U, rank;
 const int root = 0;

 Comm    = comm;
 fespace = FESpace;
 MaxSubDomainPerDof = fespace->GetMaxSubDomainPerDof();
 if(MaxSubDomainPerDof<0)
  {
   printf("Error: SetMaxSubDomainPerDof in FeSpace before calling ParFECommunicator2D \n");
   MPI_Finalize();
   exit(0);
  }

 N_DofRankIndex = NULL;
 DofRankIndex = NULL;
 GlobalDofOFLocalDof = NULL;
 OwnDofIndex = NULL;
 GlobalDofOFLocalDofAllRank = NULL;
 N_LocalDofAllRank = NULL;
 N_OwnDof = -1;
 N_ActiveOwnDof = -1;
 N_DependentCells = 0;
 DependentCellIndex = NULL;
 ActiveOwnDofIndex = NULL;
 OwnDofNeibPtr = NULL;
 OwnDofNeibPtrList = NULL;
 NeedAtNeib = NULL;
 N_DofNeibs = 0;
 DofNeibIDs = NULL;
 Max_CommunicationSteps=0;
 ReceiveID_Index = NULL;
 SendID_Index = NULL;


  N_U=0;
  MPI_Comm_rank(Comm, &rank);
  if(TDatabase::ParamDB->Par_P1) // root take part in computation
   {
    N_U = fespace->GetN_DegreesOfFreedom();
   }
  else if(rank!=root)
   {
    N_U = fespace->GetN_DegreesOfFreedom();
   }
 MPI_Reduce(&N_U, &MaxN_LocalDofAllRank, 1, MPI_INT, MPI_MAX, 0, Comm);

 ConstructDofRankIndex();

//  SetFENeibCommunicationSteps();

 MakeDofMappingFromRoot();
//  MakeDofMappingFromNeib();
} 

int TParFECommunicator2D::MakeDofMappingFromNeib()
{
 int i, j, k, l, n, rank, size, N, P;
 int N_U, GlobalDofBegin, *N_OwnDofAll, *displs, *recvcount;
 int N_LocDof, N_JointDOF, *DOF, *GlobalNumbers, *BeginIndex;
 int *N_DepCellAllNeibs, N_CellNeib, *CellNeibList, MaxN_DepCellAllNeibs;
 int *DepCellGlobalNoAllNeibs, *DepCellPtrAllNeibs, *DepLocDofPosAllNeibs, *N_DepLocValuesAllNeibs;
 int MaxN_DepLocValuesAllNeibs, *DepLocGlobalDofAllNeibs;
 int *temp1, *temp2;

 bool UPDATE;

 FE2D FeId;
 TFEDesc2D *FeDesc;
 TCollection *Coll;
 TBaseCell *cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = fespace->GetCollection();
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  N_U = fespace->GetN_DegreesOfFreedom();

  GlobalDofOFLocalDof = new int[N_U];
  N_OwnDofAll = new int[size];
  displs = new int[size];
  recvcount = new int[size];

  for(i=0; i<N_U; i++)
   GlobalDofOFLocalDof[i] = -1;

  for(i=0; i<size; i++)
   {
    displs[i]=i;
    recvcount[i]=1;
   }

  MPI_Allgatherv(&N_OwnDof, 1, MPI_INT,  N_OwnDofAll, recvcount, displs, MPI_INT,  Comm);

  GlobalDofBegin = 0;
  for(i=0; i<rank; i++)
   GlobalDofBegin +=N_OwnDofAll[i];

  // the global dof  will be assigened only for own dof, other dofs will get from neibs
  for(i=0; i<N_OwnDof; i++)
   {
    j=OwnDofIndex[i];
    GlobalDofOFLocalDof[j] = GlobalDofBegin + j;
   }

  N_DepCellAllNeibs = new int[N_DofNeibs];
  temp1 = new int[N_DofNeibs];
  temp2 = new int[N_DofNeibs];
  memset(N_DepCellAllNeibs, 0, N_DofNeibs*SizeOfInt);

  N_DepLocValuesAllNeibs =  new int[N_DofNeibs];
  memset(N_DepLocValuesAllNeibs, 0, N_DofNeibs*SizeOfInt);

  // find which neibs need info from this rank and what are the values
  for(i=0; i<N_DependentCells; i++)
   {
    N = DependentCellIndex[i];
    cell = Coll->GetCell(N);
    FeId = fespace->GetFE2D(N, cell);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_LocDof = FeDesc->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[N];
    N_CellNeib = cell->GetN_NeibProcesses();
    CellNeibList = cell->GetNeibProcessesIds();

    for(j=0; j<N_CellNeib; j++)
     {
      k =  IndexOfNeibRank[ CellNeibList[j] ];
      N_DepCellAllNeibs[k]++;

      for(l=0; l<N_LocDof; l++)
       {
        P = DOF[l];
        if( DofRankIndex[P*MaxSubDomainPerDof] == rank)
         {
          // this dof is own and need to be send to all of its neibs
          // may be accounted more than once, but this is the case (for each cell) with DepLocDofPosAllNeibs
          N_DepLocValuesAllNeibs[k] ++;
         }
       }
     }
  } // for(i=0; i<N_DependentCells; i++)

  MaxN_DepCellAllNeibs = -1;
  for(i=0; i<N_DofNeibs; i++)
   if(MaxN_DepCellAllNeibs<N_DepCellAllNeibs[i]) 
     MaxN_DepCellAllNeibs = N_DepCellAllNeibs[i];

  MaxN_DepLocValuesAllNeibs = -1;
  for(i=0; i<N_DofNeibs; i++)
   if(MaxN_DepLocValuesAllNeibs<N_DepLocValuesAllNeibs[i]) 
     MaxN_DepLocValuesAllNeibs = N_DepLocValuesAllNeibs[i];


  DepCellGlobalNoAllNeibs = new int[N_DofNeibs*MaxN_DepCellAllNeibs];
  DepCellPtrAllNeibs = new int[N_DofNeibs*(MaxN_DepCellAllNeibs+1)];
  DepLocDofPosAllNeibs = new int[N_DofNeibs*MaxN_DepLocValuesAllNeibs];
  DepLocGlobalDofAllNeibs = new int[N_DofNeibs*MaxN_DepLocValuesAllNeibs];

  for(i=0; i<N_DofNeibs; i++)
   DepCellPtrAllNeibs[i*MaxN_DepCellAllNeibs] = 0;

  memset(N_DepCellAllNeibs, 0, N_DofNeibs*SizeOfInt);
  memset(N_DepLocValuesAllNeibs, 0, N_DofNeibs*SizeOfInt);

  for(i=0; i<N_DependentCells; i++)
   {
    N = DependentCellIndex[i];
    cell = Coll->GetCell(N);
    FeId = fespace->GetFE2D(N, cell);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_LocDof = FeDesc->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[N];
    N_CellNeib = cell->GetN_NeibProcesses();
    CellNeibList = cell->GetNeibProcessesIds();

    for(j=0; j<N_CellNeib; j++)
     {
      k =  IndexOfNeibRank[ CellNeibList[j] ];
      DepCellGlobalNoAllNeibs[k*MaxN_DepCellAllNeibs + N_DepCellAllNeibs[k] ] = cell->GetGlobalCellNo();

      for(l=0; l<N_LocDof; l++)
       {
        P = DOF[l];

        if( DofRankIndex[P*MaxSubDomainPerDof] == rank)
         {

          if(GlobalDofOFLocalDof[P]==-1)
           {
            printf(" Error! Own rank Global Dof not assigned properly %d !!\n", GlobalDofOFLocalDof[P]);
            MPI_Abort(Comm, 0);
           }

          n = N_DepLocValuesAllNeibs[k];
          DepLocDofPosAllNeibs[k*MaxN_DepLocValuesAllNeibs + n] = l;
          DepLocGlobalDofAllNeibs[k*MaxN_DepLocValuesAllNeibs + n] = GlobalDofOFLocalDof[P];
          N_DepLocValuesAllNeibs[k]++;
         } //  if( DofRankIndex
        } // for(l=0; l

      DepCellPtrAllNeibs[k*MaxN_DepCellAllNeibs + N_DepCellAllNeibs[k] + 1 ]  =  n;
      N_DepCellAllNeibs[k] ++;
     } //  for(j=0; j<N_CellNeib; j++)
  } // for(i=0; i<N_DependentCells; i++)

 int *N_HaloCellAllNeibs, *HaloCellGlobalNoAllNeibs, *N_HaloLocValuesAllNeib;
 int MaxN_HaloCellAllNeibs, **sendarray, **recvarray, MaxN_HaloLocValuesAllNeibs;
 int *sendarraydisp, *recvarraydisp, *HaloCellPtrAllNeibs, *HaloLocDofPosAllNeibs, *HaloLocGlobalDofAllNeibs;
 int **sendarraylen, **recvarraylen;

  N_HaloCellAllNeibs = new int[N_DofNeibs];
  memset(N_HaloCellAllNeibs, 0, N_DofNeibs*SizeOfInt);

  N_HaloLocValuesAllNeib =  new int[N_DofNeibs];
  memset(N_HaloLocValuesAllNeib, 0, N_DofNeibs*SizeOfInt);

  // communicate to the neibs accoring to the predefined schedule
  // sendbuf[k], value to "k" the neib rank
  // recievebuf[k], recieved value from the neib rank "k"
  sendarray = new int*[4];
  recvarray = new int*[4];
  sendarraydisp = new int[4];
  recvarraydisp = new int[4];
  sendarraylen = new int*[4];
  recvarraylen = new int*[4];

  sendarray[0] = N_DepCellAllNeibs;
  sendarray[1] = N_DepLocValuesAllNeibs;

  recvarray[0] = N_HaloCellAllNeibs;
  recvarray[1] = N_HaloLocValuesAllNeib;

  this->MooNMD_FECommunicateNeib(sendarray, recvarray, 2);

  MaxN_HaloCellAllNeibs = -1;
  for(i=0; i<N_DofNeibs; i++)
   if(MaxN_HaloCellAllNeibs<N_HaloCellAllNeibs[i]) 
     MaxN_HaloCellAllNeibs = N_HaloCellAllNeibs[i];

  MaxN_HaloLocValuesAllNeibs = -1;
  for(i=0; i<N_DofNeibs; i++)
   if(MaxN_HaloLocValuesAllNeibs<N_HaloLocValuesAllNeib[i]) 
     MaxN_HaloLocValuesAllNeibs = N_HaloLocValuesAllNeib[i];

//   int outrank=3;
//   if(rank==outrank)
//    {
//     for(i=0; i<N_DofNeibs; i++)
//      printf(" rank %d DofNeibIDs %d N_DepCellAllNeibs %d\n", rank,  DofNeibIDs[i],  N_DepCellAllNeibs[i]);
// 
//     for(i=0; i<N_DofNeibs; i++)
//      printf(" rank %d  N_HaloLocValuesAllNeib %d\n", rank,   N_HaloCellAllNeibs[i]);
// 
//     printf(" rank %d MaxN_DepCellAllNeibs %d MaxN_HaloCellAllNeibs %d\n", rank, MaxN_DepCellAllNeibs, MaxN_HaloCellAllNeibs);
//    }

//    if(IndexOfNeibRank[outrank]!=-1)
//     printf(" rank %d neib rank %d N_HaloLocValuesAllNeib %d\n", rank,  DofNeibIDs[IndexOfNeibRank[outrank]],  N_HaloCellAllNeibs[IndexOfNeibRank[outrank]]);
// 

//     printf(" rank %d MaxN_DepCellAllNeibs %d MaxN_HaloCellAllNeibs %d\n", rank, MaxN_DepCellAllNeibs, MaxN_HaloCellAllNeibs);
//   MPI_Finalize();
//   exit(0); 

  HaloCellGlobalNoAllNeibs = new int[N_DofNeibs*MaxN_HaloCellAllNeibs];
  HaloCellPtrAllNeibs = new int[N_DofNeibs*(MaxN_HaloCellAllNeibs+1)];
  HaloLocDofPosAllNeibs = new int[N_DofNeibs*MaxN_HaloLocValuesAllNeibs];
  HaloLocGlobalDofAllNeibs = new int[N_DofNeibs*MaxN_HaloLocValuesAllNeibs];

  // pack-->send and recv-->unpack all data between neibs
  sendarraydisp[0] = MaxN_DepCellAllNeibs;
  sendarraydisp[1] = MaxN_DepCellAllNeibs+1;
  sendarraydisp[2] = MaxN_DepLocValuesAllNeibs;
  sendarraydisp[3] = MaxN_DepLocValuesAllNeibs;

  recvarraydisp[0] = MaxN_HaloCellAllNeibs;
  recvarraydisp[1] = MaxN_HaloCellAllNeibs+1;
  recvarraydisp[2] = MaxN_HaloLocValuesAllNeibs;
  recvarraydisp[3] = MaxN_HaloLocValuesAllNeibs;

  for(i=0; i<N_DofNeibs; i++)
   temp1[i] = N_DepCellAllNeibs[i]+1;

  sendarraylen[0] = N_DepCellAllNeibs;
  sendarraylen[1] = temp1;
  sendarraylen[2] = N_DepLocValuesAllNeibs;
  sendarraylen[3] = N_DepLocValuesAllNeibs;

  for(i=0; i<N_DofNeibs; i++)
   temp2[i] = N_HaloCellAllNeibs[i]+1;

  recvarraylen[0] = N_HaloCellAllNeibs;
  recvarraylen[1] = temp2;
  recvarraylen[2] = N_HaloLocValuesAllNeib;
  recvarraylen[3] = N_HaloLocValuesAllNeib;

  sendarray[0] = DepCellGlobalNoAllNeibs;
  sendarray[1] = DepCellPtrAllNeibs;
  sendarray[2] = DepLocDofPosAllNeibs;
  sendarray[3] = DepLocGlobalDofAllNeibs;

  recvarray[0] = HaloCellGlobalNoAllNeibs;
  recvarray[1] = HaloCellPtrAllNeibs;
  recvarray[2] = HaloLocDofPosAllNeibs;
  recvarray[3] = HaloLocGlobalDofAllNeibs;

  this->MooNMD_FECommunicateNeib(sendarray,  sendarraydisp, sendarraylen, recvarray, recvarraydisp, recvarraylen, 4);

//   for(i=0; i<N_DofNeibs; i++)
//    HaloCellPtrAllNeibs[i*MaxN_HaloCellAllNeibs] = 0;


    printf(" rank %d N_HaloLocValuesAllNeib %d\n", rank, MaxN_HaloLocValuesAllNeibs);
//     printf(" rank %d GlobalDofBegin %d\n", rank, GlobalDofBegin);


  delete [] N_OwnDofAll;
  delete [] displs;
  delete [] recvcount;
  delete [] N_DepCellAllNeibs;
  delete [] DepCellGlobalNoAllNeibs;
  delete [] DepLocDofPosAllNeibs;
  delete [] N_DepLocValuesAllNeibs;
  delete [] DepLocGlobalDofAllNeibs;
  delete [] sendarray;
  delete [] recvarray;

  delete [] N_HaloCellAllNeibs;
  delete [] N_HaloLocValuesAllNeib;
  delete [] HaloCellGlobalNoAllNeibs;
  delete [] HaloCellPtrAllNeibs;
  delete [] HaloLocDofPosAllNeibs;
  delete [] HaloLocGlobalDofAllNeibs;


  MPI_Finalize();
  exit(0); 


  return 0;
} 


int TParFECommunicator2D::MakeDofMappingFromRoot()
{
 int i, j, k, l, M, rank, size, N_Cells;
 int N_U, N_LocDof;
 int *IntArray = new int[2], *SubDomainParentCellNo, *SubDomainDofGlobalNumbers;
 int *SubDomainGlobalNumbers, *SubDomainBeginIndex;
 int *GlobalNumbers, *BeginIndex, *DOF, *SubDomainDOF, P, N;
 int Type = TDatabase::ParamDB->Par_P4, Disp, ID, begin;
 const int root = 0;

 bool UPDATE;

 double x, y;

 TCollection *Coll;
 TBaseCell *cell;

 MPI_Status status;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();



//   N_U=0;
// 
//   if(TDatabase::ParamDB->Par_P1) // root take part in computation
//    {
//     N_U = fespace->GetN_DegreesOfFreedom();
//    }
//   else if(rank!=root)
//    {
//     N_U = fespace->GetN_DegreesOfFreedom();
//    }
// 
//   MPI_Reduce(&N_U, &MaxN_LocalDofAllRank, 1, MPI_INT, MPI_MAX, root, Comm);

//   if(rank==root)
   N_U = fespace->GetN_DegreesOfFreedom();

  N_LocalDofAllRank = new int[size];
  N_LocalDofAllRank[rank] =  N_U;


  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();
  N_LocDof = BeginIndex[1] - BeginIndex[0];  // assume that  all cells in fespace have same FE

   if(rank==root)
    {
     int N_SubDomainCells, N_SubDomainDOF;

     N_LocalDofAllRank[0] = N_U;
 
     if(TDatabase::ParamDB->Par_P1) // root take part in computation
      {
       GlobalDofOFLocalDofAllRank = new int[size*MaxN_LocalDofAllRank];
       for(i=0;i<MaxN_LocalDofAllRank;i++)
         GlobalDofOFLocalDofAllRank[i] = -1;

        cout << "not yet implemented !!!!!!!!!!!!!" << endl;
        cout << "Check ParCommunicator !!!!!!!!!!!!!" << endl;
        MPI_Finalize();
        exit(0);
      }
     else
      {
       GlobalDofOFLocalDofAllRank = new int[(size-1)*MaxN_LocalDofAllRank];
      }

     for(j=1;j<size;j++)
      {
       MPI_Recv(IntArray, 2, MPI_INT, j, 100, Comm, &status);

       N_SubDomainCells = IntArray[0]; // including halo cells if any
       N_SubDomainDOF = IntArray[1];   // including halo cells DOF if any   
       N_LocalDofAllRank[j] = N_SubDomainDOF;

       if(j>1)
        {
         delete [] SubDomainParentCellNo;
         delete [] SubDomainGlobalNumbers;
         delete [] SubDomainBeginIndex;
        }

       /** get the subdomain cells' global cell numbers */
       SubDomainParentCellNo = new int[N_SubDomainCells];
       MPI_Recv(SubDomainParentCellNo, N_SubDomainCells, MPI_INT, j, 200, Comm, &status);

       SubDomainGlobalNumbers = new int[N_SubDomainCells*N_LocDof];
       MPI_Recv(SubDomainGlobalNumbers, N_SubDomainCells*N_LocDof, MPI_INT, j, 300, Comm, &status);

       SubDomainBeginIndex = new int[N_SubDomainCells];
       MPI_Recv(SubDomainBeginIndex, N_SubDomainCells, MPI_INT, j, 400, MPI_COMM_WORLD, &status);

//       if(j==2)
//        {
//         for(i=0; i<N_SubDomainCells; i++)
//          printf("Localcell %d:  SubDomainParentCellNo  : %d \n",i,  SubDomainParentCellNo[i]);
//        }

      /** mapping begin*/
      for(k=0; k<N_SubDomainCells; k++)
       {
        M = SubDomainParentCellNo[k];
        DOF = GlobalNumbers + BeginIndex[M];
        SubDomainDOF = SubDomainGlobalNumbers + SubDomainBeginIndex[k];
//        if(j==1)
//         printf("Localcell %d:  SubDomainParentCellNo  : %d \n",k,  M);
        for(l=0; l<N_LocDof; l++)
         {
          N=DOF[l];
          P = SubDomainDOF[l];

//           if(j==2)
//            printf("Localcell %d: LocDof%d,  GlobalDofOF : %d \n",k, P,  N);

          if(TDatabase::ParamDB->Par_P1) // root take part in computation
            {
             GlobalDofOFLocalDofAllRank[j*MaxN_LocalDofAllRank +  P] = N;
            }
          else
           {
            GlobalDofOFLocalDofAllRank[(j-1)*MaxN_LocalDofAllRank +  P] = N;
//             if(j==1)
//               {
//               printf("Localcell %d: LocDof%d,  GlobalDofOFLocalDof  : %d \n",k, P,  GlobalDofOFLocalDofAllRank[(j-1)*MaxN_LocalDofAllRank +  P]);
//              }
           }
         } // for(l=0; l<N_LocDof; l++)
        } // for(i=0; i<N_SubDomainCells; i++)

//       if(j==2)
//        {
//         for(i=0; i<N_SubDomainDOF; i++) //  if(TDatabase::ParamDB->Par_P1) // root take part in computation
//          printf("LocalDof %d:  GlobalDofOFLocalDof  : %d \n",i,  GlobalDofOFLocalDofAllRank[(j-1)*MaxN_LocalDofAllRank+ i]);
//        }

      if(j>1)
       MPI_Wait(&request004, MPI_STATUS_IGNORE);

      if(TDatabase::ParamDB->Par_P1) // root take part in computation
       {
        MPI_Isend(GlobalDofOFLocalDofAllRank+j*MaxN_LocalDofAllRank, N_SubDomainDOF, MPI_INT, j, 500, Comm, &request004);
       }
      else
       {
        MPI_Isend(GlobalDofOFLocalDofAllRank+(j-1)*MaxN_LocalDofAllRank, N_SubDomainDOF, MPI_INT, j, 500, Comm, &request004);
       }

      } // for(j=1;j<size;j++)

     delete [] SubDomainParentCellNo;
     delete [] SubDomainGlobalNumbers;
     delete [] SubDomainBeginIndex;
    }
   else
    {
     IntArray[0] = N_Cells;
     IntArray[1] = N_U;
     MPI_Send(IntArray, 2, MPI_INT, 0, 100, Comm);

     SubDomainParentCellNo = Coll->GetGlobalIndex();
     MPI_Send(SubDomainParentCellNo, N_Cells, MPI_INT, 0, 200, Comm);

     MPI_Send(GlobalNumbers, N_Cells*N_LocDof, MPI_INT, 0, 300, Comm);
     MPI_Send(BeginIndex, N_Cells, MPI_INT, 0, 400, Comm);

     GlobalDofOFLocalDof = new int[N_U];
     MPI_Irecv(GlobalDofOFLocalDof, N_U, MPI_INT, 0, 500, Comm, &request004);


//       WaitForMakeDofMapping();
//       if(rank==2)
//        {
//         for(i=0; i<N_U; i++)
//          printf("In rank LocalDof %d:  GlobalDofOFLocalDof  : %d \n",i,  GlobalDofOFLocalDof[i]);
//        }

    }


//     printf("Main Parcommunicator2D.C \n" );
//     MPI_Finalize();
//     exit(0); 

  delete [] IntArray;
  return 0;
} 


int TParFECommunicator2D::ConstructDofRankIndex()
{
 int i, j, k, N_Cells, ID, M, N, N_LocDof, N_U, Disp, rank;
 int *DOF, *GlobalNumbers, *BeginIndex;
 int Type = TDatabase::ParamDB->Par_P4, N_Active;
 int N_Joints, N_JointDOF, N_InnerDOF;

 FE2D FeId;
 TFEDesc2D *FeDesc;
 TJoint *Joint;
 TSubDomainJoint *SubDJoint;

 bool UPDATE;

  MPI_Comm_rank(Comm, &rank);

  TCollection *Coll;
  TBaseCell *cell;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = fespace->GetN_DegreesOfFreedom();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();
  N_LocDof = BeginIndex[1] - BeginIndex[0];  // assume that  all cells in fespace have same FE

  N_DofRankIndex = new int[N_U];
  DofRankIndex = new int[MaxSubDomainPerDof*N_U];
  memset(N_DofRankIndex, 0, N_U*SizeOfInt);

  Disp = N_U*MaxSubDomainPerDof;
  for(i=0; i<Disp; i++)
   DofRankIndex[i] = -1;



     // find how many SubDomain contain each dof and their corresponding IDs
     for(i=0; i<N_Cells; i++)
      {
       cell = Coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       DOF = GlobalNumbers + BeginIndex[i];

       for(j=0; j<N_LocDof; j++)
        { 
         N = DOF[j];
         M = N_DofRankIndex[N];
         UPDATE = TRUE;

         for(k=0; k<M; k++)
          if(ID==DofRankIndex[N*MaxSubDomainPerDof + k ])
           {
            UPDATE = FALSE;
            break;
           } // for(k=0; k<M; k++)

         if(UPDATE)
          {
           DofRankIndex[N*MaxSubDomainPerDof + M ] = ID;
           N_DofRankIndex[N]++;
          }
         }  // for(j=0; j<N_DOF; j++)
        } // for(i=0; i<N_Cells; i++)

//      // find how many SubDomain contain each dof and their corresponding IDs
//      for(i=0; i<N_Cells; i++)
//       {
//        cell = Coll->GetCell(i);
//        N_Joints = cell->GetN_Joints();
//        FeId = fespace->GetFE2D(i, cell);
//        FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
//        N_LocDof = FeDesc->GetN_DOF();
//        N_JointDOF = FeDesc->GetN_JointDOF();
//        N_InnerDOF = FeDesc->GetN_InnerDOF();
// 
//        DOF = GlobalNumbers + BeginIndex[i];
//        ID = cell->GetSubDomainNo();
// 
//        if(N_LocDof== N_InnerDOF)
//         {
//          //disc space 
//          for(j=0; j<N_LocDof; j++)
//           {
//            N = DOF[j];
//            M = N_DofRankIndex[N];
//            DofRankIndex[N*MaxSubDomainPerDof + M ] = ID;
//            N_DofRankIndex[N]++;
//           }
// 
//         }
//        else if(N_LocDof == N_Joints*N_JointDOF)
//         {
//          //non-conf space
//          for(j=0; j<N_Joints; j++)
//           {
//            Joint = cell->GetJoint(j);
//            if(Joint->GetType() == SubDomainJoint )
//             {
//              SubDJoint = (TSubDomainJoint *)Joint;
//             }
//    
//           }
//  
//         }
//        else
//         {
//          //cont space 
//         }
// 
//       } // for(i=0; i<N_Cells; i++)
// 

//     for(i=0; i<N_U; i++)
//      {
//       N = N_DofRankIndex[i];
//        for(j=0; j<N; j++)
//          printf("dof %d: rank index : %d \n",i, DofRankIndex[i*MaxSubDomainPerDof + j]);
//       fespace->GetDOFPosition(i, x, y);
//          printf("dof %d: x : %8.2f, y : %8.2f \n",i,x, y);
//       printf(" \n");
//      }

//    printf("ConstructDofRankIndex \n" );

//    MPI_Finalize();
//    exit(0);
    // set the SubDomain interface dof's rank index by various methods
    for(i=0; i<N_U; i++)
     {
      N = N_DofRankIndex[i];

      if(N<1)
       {
        cout<<" Error in  GetDofRankIndex !!" <<endl;
        MPI_Abort(Comm, 0);
       } 
      else if(N>1)
       {
       // more than one subdomain contain this dof
         switch(Type)
         {
          case 0: 
            // Min(rank ID) will take this dof
            for(k=1; k<N; k++)
             if(DofRankIndex[i*MaxSubDomainPerDof] > DofRankIndex[i*MaxSubDomainPerDof + k ])
              {
               M = DofRankIndex[i*MaxSubDomainPerDof + k ];
               DofRankIndex[i*MaxSubDomainPerDof + k ] = DofRankIndex[i*MaxSubDomainPerDof];
               DofRankIndex[i*MaxSubDomainPerDof] = M;
              } 

          break;
          case 1: 
            // Max(rank ID) will take this dof
            for(k=1; k<N; k++)
             if(DofRankIndex[i*MaxSubDomainPerDof] < DofRankIndex[i*MaxSubDomainPerDof + k ])
              {
               M = DofRankIndex[i*MaxSubDomainPerDof + k ];
               DofRankIndex[i*MaxSubDomainPerDof + k ] = DofRankIndex[i*MaxSubDomainPerDof];
               DofRankIndex[i*MaxSubDomainPerDof] = M;
              } 
          break;

          default: 
            cout << "Wrong type in GetDofRankIndex" << endl;
            MPI_Abort(Comm, 0);
          break;
        } // switch

       } // else if(N>1)

      } // for(i=0; i<N_U; i++)


//     for(i=0; i<N_U; i++)
//      printf("dof %d: rank index : %d \n",i, DofRankIndex[i*MaxSubDomainPerDof]);

  /** needed for SetOwnDof() */
  OwnDofIndex = new int[N_U];

  N_OwnDof = 0;
  N_ActiveOwnDof = 0;
  N_Active = fespace->GetActiveBound();

  ActiveOwnDofIndex = new int[N_Active];
  for(i=0; i<N_U; i++)
   {
    if(DofRankIndex[i*MaxSubDomainPerDof] ==rank)
     {
      OwnDofIndex[N_OwnDof] = i;
      N_OwnDof ++;  

      if(i<N_Active)
       {
        ActiveOwnDofIndex[N_ActiveOwnDof] = i;
        N_ActiveOwnDof++;
       }
     }
   }

 return 0;
}

int TParFECommunicator2D::WaitForMakeDofMapping()
{
//   MPI_Status status;

//   MPI_Wait(&request001, MPI_STATUS_IGNORE);
//   MPI_Wait(&request002, MPI_STATUS_IGNORE);
//   MPI_Wait(&request003, MPI_STATUS_IGNORE);
     MPI_Wait(&request004, MPI_STATUS_IGNORE);
//   MPI_Wait(&request005, MPI_STATUS_IGNORE);
//   MPI_Wait(&request006, MPI_STATUS_IGNORE);

  return 0;
}

// set the own dof excluding the halo dof (in slave process)
void TParFECommunicator2D::SetOwnDof()
{
 int i, N_U, rank, size, N_Active;


    printf("No need to call SetOwnDof()! already did in ConstructDofRankIndex() \n" );
    MPI_Finalize();
    exit(0); 

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  N_U = fespace->GetN_DegreesOfFreedom();

  N_OwnDof = 0;
  N_ActiveOwnDof = 0;
  N_Active = fespace->GetActiveBound();

  ActiveOwnDofIndex = new int[N_Active];
  for(i=0; i<N_U; i++)
   {
    if(DofRankIndex[i*MaxSubDomainPerDof] ==rank)
     {
      OwnDofIndex[N_OwnDof] = i;
      N_OwnDof ++;  

      if(i<N_Active)
       {
        ActiveOwnDofIndex[N_ActiveOwnDof] = i;
        N_ActiveOwnDof++;
       }
     }
   }
//     printf(" rank %d, value %d\n",rank, N_ActiveOwnDof);

}

int TParFECommunicator2D::GetN_OwnDof()
 {
  if(N_OwnDof<0)
   {
    cout << "Call SetOwnDof() before calling GetN_OwnDof()" << endl;
    MPI_Abort(Comm, 0);
   }
  else
   {
    return N_OwnDof;
   }
 }

int* TParFECommunicator2D::GetOwnDofIndex()
{
  if(N_OwnDof<0)
   {
    cout << "Call SetOwnDof() before calling GetN_OwnDof()" << endl;
    MPI_Abort(Comm, 0);
   }
  else
   {
    return OwnDofIndex;
   }
}

int TParFECommunicator2D::GetActiveN_OwnDof()
 {
  if(N_ActiveOwnDof<0)
   {
    cout << "Call SetOwnDof() before calling GetN_OwnDof()" << endl;
    MPI_Abort(Comm, 0);
   }
  else
   {
    return N_ActiveOwnDof;
   }
 }

int* TParFECommunicator2D::GetActiveOwnDofIndex()
{
  if(N_ActiveOwnDof<0)
   {
    cout << "Call SetOwnDof() before calling GetN_OwnDof()" << endl;
    MPI_Abort(Comm, 0);
   }
  else
   {
    return ActiveOwnDofIndex;
   }
}

// set steps and order for communication between processes
// set the scheduling algorightm, which rank communicates to which rank and when and how 
void TParFECommunicator2D::SetFENeibCommunicationSteps()
{
 int i, j, k, l, m, N, P, Q, N_JointDOF, NeibID, N_Neibs, rank, size, N_OwnCells;
 int N_LocDof, *N_Neibs_All, *NeibList_All, N_U;
 int *NeibList, *neibsHalloCellIndex, *JointDOF, *GlobalNumbers, *BeginIndex;
 int *DOF, *NeibRanks, total, *dofNeibIDs, N_DofNeibs_MaxAll;
 int *pos, *N_DofNeibs_Array, *DofNeibIDs_Array, N_Entries;


 bool UPDATE;

 FE2D FeId;
 TFEDesc2D *FeDesc;
 TCollection *Coll;
 TBaseCell *cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = fespace->GetCollection();
  N_OwnCells = Coll->GetN_OwnCells();
  neibsHalloCellIndex = new int [N_OwnCells];
  N_DependentCells = 0;
  for(i=0; i<N_OwnCells; i++)
   {
    cell = Coll->GetCell(i);
    if(cell->IsDependentCell())
     {
      neibsHalloCellIndex[N_DependentCells] = i;
      N_DependentCells++;
    }
   } //  for(i=0; i<N_Cells; i++)

  DependentCellIndex = new int[N_DependentCells];
  memcpy(DependentCellIndex, neibsHalloCellIndex, N_DependentCells*SizeOfInt);
  delete [] neibsHalloCellIndex;


  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  N_U = fespace->GetN_DegreesOfFreedom();

  N_DependentCellNeibs = new int [size];
//   DependentCellNeibIDs = new int [size*N_DependentCells];
  NeibList = new int [size];
  N_Neibs_All = new int[N_U];
  NeibList_All =  new int[N_U*MaxSubDomainPerDof];
  memset(N_Neibs_All, 0, N_U*SizeOfInt);
  memset(N_DependentCellNeibs, 0, size*SizeOfInt);



  // find which neibs need info from this rank
  for(i=0; i<N_DependentCells; i++)
   {
    //find number of neibs associated with this cell which are required info from this cell DOFs
    //check neibs of each DOF of this cell
    //enough to check DOFs associated on vertices of this cell, i.e, P1 part only
    //in case of disc-fespace no neibs
    N = DependentCellIndex[i];
    cell = Coll->GetCell(N);
    k = cell->GetN_Joints();
    FeId = fespace->GetFE2D(N, cell);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_LocDof = FeDesc->GetN_DOF();
    N_JointDOF = FeDesc->GetN_JointDOF();
    DOF = GlobalNumbers+BeginIndex[N];

    N_Neibs = 0;
    if(N_JointDOF) // no nebs for disc fespaces
    {
     for(j=0;j<k;j++)
      {
       JointDOF = FeDesc->GetJointDOF(j);
       P = DOF[JointDOF[0]];  // only first DOF of each joint

       Q = N_DofRankIndex[P];
       NeibRanks = DofRankIndex + P*MaxSubDomainPerDof;

       for(l=0;l<Q;l++)
        {
         NeibID = NeibRanks[l];
         if(NeibID==rank)
          continue;

         UPDATE = TRUE;
         for(m=0;m<N_Neibs;m++)
          {
           if(NeibID==NeibList[m])
            {
             UPDATE = FALSE;
             break;
            }
          }

         if(UPDATE)
          {
           NeibList[N_Neibs] = NeibID;
           N_Neibs++;
          }
        }
      } //  for(j=0; j<k ; j++

    // set info of how many neibs need this cell values and their IDs
    cell->SetN_NeibProcesses(N_Neibs);
    cell->SetNeibProcessesIds(NeibList);

     // in this cell, now collect info for only own DOFs, which values are needed in neib ranks
     for(j=0;j<N_LocDof;j++)
      {
       P = DOF[j];
       if( DofRankIndex[P*MaxSubDomainPerDof] == rank)
        {
         // this DOF is own Dof, and its value has to be send to its neib ranks, which are having this cell as a Halo cell
         for(l=0;l<N_Neibs;l++)
          {
           UPDATE = TRUE;
           Q = N_Neibs_All[P];

           for(m=0;m<Q;m++)
            {
             if(NeibList_All[P*MaxSubDomainPerDof + m]==NeibList[l])
              {
               UPDATE = FALSE;
               break;
              }
            }

           if(UPDATE)
            {
             NeibList_All[P*MaxSubDomainPerDof + Q] = NeibList[l];
             N_Neibs_All[P]++;
            }     
           } //  for(l=0;l<N_Neibs;l++)
         } // if( DofRankIndex[P*MaxSubD
        } // for(j=0;j<N_LocDof;j

     } // if(N_JointDOF)
   } // for(i=0; i<N_DependentCells; i++)

  total = 0;
  N_DofNeibs=0;
  dofNeibIDs= new int[size];
  for(i=0;i<N_U;i++)
   {
    N = N_Neibs_All[i];
    total += N;

    // find all neibs of this rank for setting communication steps
    for(j=0;j<N;j++)
     {
      NeibID = NeibList_All[i*MaxSubDomainPerDof + j];
      UPDATE = TRUE;
      for(m=0;m<N_DofNeibs;m++)
       {
        if(dofNeibIDs[m]==NeibID)
         {
          UPDATE=FALSE;
          break;
         }
       }

      if(UPDATE)
       {
        dofNeibIDs[N_DofNeibs]=NeibID;
        N_DofNeibs++;
       }
     }
   }

  DofNeibIDs = new int[N_DofNeibs];
  memcpy(DofNeibIDs, dofNeibIDs, N_DofNeibs*SizeOfInt);
  delete [] dofNeibIDs;

  IndexOfNeibRank = new int[size];

  for(i=0; i<size; i++)
   IndexOfNeibRank[i] = -1;

  for(i=0; i<N_DofNeibs; i++)
   IndexOfNeibRank[DofNeibIDs[i]] = i;

//   printf("Number of exact DOF Neib-processes of process %d is %d\n",rank, N_DofNeibs);

  MPI_Allreduce(&N_DofNeibs, &N_DofNeibs_MaxAll, 1, MPI_INT, MPI_MAX, Comm); 

  if(rank==0)
   {
    pos = new int[size];
    for(i=0;i<size;i++)
     pos[i] = i*N_DofNeibs_MaxAll;

    N_DofNeibs_Array = new int[size];
    DofNeibIDs_Array= new int[size*N_DofNeibs_MaxAll];
   }


   MPI_Gather(&N_DofNeibs, 1, MPI_INT, N_DofNeibs_Array, 1, MPI_INT, 0, Comm);
   MPI_Gatherv(DofNeibIDs, N_DofNeibs, MPI_INT, DofNeibIDs_Array, N_DofNeibs_Array, pos, 
               MPI_INT, 0, Comm);


   if(rank==0)
    {
     int step, *Entries, *ColInd, *RowPtr, *Collection, MaxPossibleSteps;

     N_Entries = 0;

     for(i=0;i<size;i++)
      N_Entries +=N_DofNeibs_Array[i];

     // adjacency matrix memory allocation
     Entries = new int[N_Entries];
     ColInd     = new int[N_Entries];
     RowPtr = new int[size+1];

     RowPtr[0] = 0;
     m = 0;
     for(i=0;i<size;i++)
      {
       N =  N_DofNeibs_Array[i];
       RowPtr[i+1] = RowPtr[i] + N;

       for(j=0;j<N;j++)
        {
         ColInd[m] = DofNeibIDs_Array[i*N_DofNeibs_MaxAll + j];
         // assmeble the matrix
         Entries[m] =-1;
         m++;
        }
      }

    /* print adjacency matrix
    for(i=0;i<size;i++)
     for(j=RowPtr[i];j<RowPtr[i+1];j++)
       printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);
    */

     // find number of steps to communicate between neibs (edge coloring alogrithm)
     Collection = new int[size];
     Max_CommunicationSteps = -1;
     MaxPossibleSteps = N_DofNeibs_MaxAll+1;

     for(step=0;step<MaxPossibleSteps;step++)
      {
       m=0;
       for(k=0;k<size;k++)
        {
         UPDATE = TRUE;
         for(l=0;l<m;l++)
           if(Collection[l] == k)
            {
             UPDATE = FALSE;
             break;
            }

         if(UPDATE)
          {
           Collection[m] = k;
            m++;
           }
          else
           {
            continue;
           }
 
         for(l=RowPtr[k];l<RowPtr[k+1];l++)
          { 
           // only upper tirangular martix is enough (due to symmetry)
           if( ColInd[l]<k) continue;

           UPDATE = TRUE;
           for(P=0;P<m;P++)
            if(Collection[P] == ColInd[l])
             {
              // ColInd[l] is already in the collection
              UPDATE = FALSE;
              break;
             } 
       
            if(UPDATE)
             {
              if(Entries[l]<0)
               {
                Entries[l] = step;
                if(Max_CommunicationSteps<step+1) Max_CommunicationSteps=step+1;
                Collection[m] = ColInd[l];
                m++;
                break;
               }
             }
          } //  for(l=RowPtr[k];
        } // for(k=0;
       } //  for(step=0;step<

     printf("Total number of steps needed to communicate is %d\n",Max_CommunicationSteps);

    /* print adjacency matrix
    for(i=0;i<size;i++)
     for(j=RowPtr[i];j<RowPtr[i+1];j++)
       printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);
    */
    N_CommunicationProcesses = new int[Max_CommunicationSteps];
    SendID_Index= new int[Max_CommunicationSteps*size];
    ReceiveID_Index = new int[Max_CommunicationSteps*size];

    for(j=0;j<Max_CommunicationSteps;j++)
     N_CommunicationProcesses[j] =0;

    for(i=0;i<Max_CommunicationSteps*size;i++)
     {
      SendID_Index[i] = -1;
      ReceiveID_Index[i] = -1;
     }

    for(i=0;i<size;i++)
     for(j=RowPtr[i];j<RowPtr[i+1];j++)
       if(i<ColInd[j])
        {
         if( Entries[j]<0)
          { 
           printf("Error in finding DOF communication steps %d\n",Entries[j]);
           MPI_Abort(MPI_COMM_WORLD,  0);
          }
         k =  Entries[j];
         SendID_Index[k*Max_CommunicationSteps + N_CommunicationProcesses[k]] = i;
         ReceiveID_Index[k*Max_CommunicationSteps + N_CommunicationProcesses[k]] =  ColInd[j];
         N_CommunicationProcesses[k]++;
        }

     for(i=0;i<Max_CommunicationSteps;i++)
      {
       printf("Step %d:\n",i);
       N= N_CommunicationProcesses[i];
       for(j=0;j<N;j++)
        printf("Send from process %d to %d -------- and onto\n",
                   SendID_Index[i*Max_CommunicationSteps +j], ReceiveID_Index[i*Max_CommunicationSteps +j]);
       printf("\n");
      }

    delete [] Collection;
    } // if(rank==0)

   MPI_Bcast(&Max_CommunicationSteps, 1, MPI_INT, 0, Comm);
   MPI_Barrier(Comm);

   if(rank!=0)
    {
     N_CommunicationProcesses = new int[Max_CommunicationSteps];
     SendID_Index= new int[Max_CommunicationSteps*size];
     ReceiveID_Index = new int[Max_CommunicationSteps*size];
    }

   MPI_Bcast(N_CommunicationProcesses, Max_CommunicationSteps, MPI_INT, 0, Comm);
   MPI_Bcast(SendID_Index, Max_CommunicationSteps*size, MPI_INT, 0, Comm);
   MPI_Bcast(ReceiveID_Index, Max_CommunicationSteps*size, MPI_INT,  0, Comm);

   delete [] NeibList;
   delete [] N_Neibs_All;
   delete [] NeibList_All;

//            printf("Parallel SetIntraCommunicator ranks  %d, N_DofNeibs %d\n", rank, N_DofNeibs);
//            MPI_Finalize();
//            exit(0);
}

// sendbuf[k], value to "k" the neib rank
// recievebuf[k], recieved value from the neib rank "k"
//send and recieve one value between the neibs
void TParFECommunicator2D::MooNMD_FECommunicateNeib(int **sendbuf, int **recvbuf, int N_array)
{
 int i, j, k, l, N, rank;
 int sendID, recvID, *val_send, *val_recv;

  val_send = new int[N_array];
  val_recv = new int[N_array];

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);

  for(i=0; i<Max_CommunicationSteps; i++)
   {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*Max_CommunicationSteps +j];
      recvID = ReceiveID_Index[i*Max_CommunicationSteps +j];

      if(rank==sendID)
       {
        for(l=0;l<N_array; l++)
          val_send[l] = sendbuf[l][ IndexOfNeibRank[recvID] ];

        MPI_Send(val_send, N_array, MPI_INT,  recvID, 100, Comm);
       }
      else if(rank==recvID)
       {
        MPI_Recv(val_recv, N_array, MPI_INT,  sendID, 100, Comm, &status);

        for(l=0;l<N_array; l++)
         recvbuf[l][IndexOfNeibRank[sendID]] = val_recv[l];
       }

      if(rank==recvID)
       {
        for(l=0;l<N_array; l++)
        val_send[l] = sendbuf[l][ IndexOfNeibRank[sendID] ];

        MPI_Send(val_send, N_array, MPI_INT,  sendID, 200, Comm);
       }
      else if(rank==sendID)
       {
        MPI_Recv(val_recv, N_array, MPI_INT,  recvID, 200, Comm, &status);

        for(l=0;l<N_array; l++)
         recvbuf[l][IndexOfNeibRank[recvID]] = val_recv[l];
       }
     } //for(j=0;j<N
   } //for(i=0; i<Ma

  delete [] val_send;
  delete [] val_recv;
}

 void TParFECommunicator2D::MooNMD_FECommunicateNeib(int **sendbuf, int *senddisp, int **sendlen, int **recvbuf, int *recvdisp, int **recvlen, int N_array)
 {
 int i, j, k, l, m, N, M, rank, totsendlen=0, totrecvlen=0, pos;
 int sendID, recvID, *buffer, *Sendbuf, *Recvbuf, *tmp;
 int len, max=0;

 for(l=0;l<N_array; l++)
  {
   totrecvlen +=senddisp[l];
   totsendlen +=recvdisp[l];
  }
 
  if(totsendlen<totrecvlen)
    max = totrecvlen;
  else
    max = totsendlen;


  buffer = new int[max];
  tmp = new int[max];

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);

  for(i=0; i<Max_CommunicationSteps; i++)
   {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*Max_CommunicationSteps +j];
      recvID = ReceiveID_Index[i*Max_CommunicationSteps +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
        pos = 0;
//         len = 0;
        for(l=0;l<N_array; l++)
         {
          Sendbuf = sendbuf[l];
          MPI_Pack(Sendbuf+M*senddisp[l], sendlen[l][M], MPI_INT, buffer, max*SizeOfInt,  &pos, Comm);
//           len += sendlen[l][M];
         }

//         printf(" rank %d  pos   %d len %d\n", rank, pos, len*SizeOfInt);
        MPI_Send(buffer, pos, MPI_PACKED,  recvID, 100, Comm);
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        len = 0;
        for(l=0;l<N_array; l++)
         len += recvlen[l][M];

//          printf(" rank %d  recvd from %d len %d\n", rank, sendID, len);
        MPI_Recv(buffer, len*SizeOfInt, MPI_PACKED,  sendID, 100, Comm, &status);

        pos = 0;
        for(l=0;l<N_array; l++)
         {
          MPI_Unpack(buffer, len*SizeOfInt, &pos, tmp, recvlen[l][M], MPI_INT, Comm);

          Recvbuf = recvbuf[l];
          for(m=0;m<recvlen[l][M]; m++)
           {
            Recvbuf[M*recvdisp[l] + m] = tmp[m];
//              printf(" rank %d  recvd from %d tmp[m] %d\n", rank, sendID, tmp[m]);
           }
         }
       }

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        pos = 0;

        for(l=0;l<N_array; l++)
         {
          Sendbuf = sendbuf[l];
          MPI_Pack(Sendbuf+M*senddisp[l], sendlen[l][M], MPI_INT, buffer, max*SizeOfInt,  &pos, Comm);
//           len += sendlen[l][M];
         }

//         printf(" rank %d  pos   %d len %d\n", rank, pos, len*SizeOfInt);
        MPI_Send(buffer, pos, MPI_PACKED,  sendID, 100, Comm);
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
        len = 0;
        for(l=0;l<N_array; l++)
         len += recvlen[l][M];

//          printf(" rank %d  recvd from %d len %d\n", rank, sendID, len);
        MPI_Recv(buffer, len*SizeOfInt, MPI_PACKED,  recvID, 100, Comm, &status);

        pos = 0;
        for(l=0;l<N_array; l++)
         {
          MPI_Unpack(buffer, len*SizeOfInt, &pos, tmp, recvlen[l][M], MPI_INT, Comm);

          Recvbuf = recvbuf[l];
          for(m=0;m<recvlen[l][M]; m++)
           {
            Recvbuf[M*recvdisp[l] + m] = tmp[m];
//              printf(" rank %d  recvd from %d tmp[m] %d\n", rank, sendID, tmp[m]);
           }
         }
       }
     } //for(j=0;j<N
   } //for(i=0; i<Ma

  delete [] buffer;
  delete [] tmp;
//     printf(" rank %d test sendlen %d recvlen %d\n", rank, totsendlen, totrecvlen);
//   MPI_Finalize();
//   exit(0); 
}


TParFECommunicator2D::~TParFECommunicator2D()
{
 if(N_DofRankIndex != NULL)
    delete [] N_DofRankIndex;
 if(DofRankIndex != NULL)
    delete [] DofRankIndex;
 if(GlobalDofOFLocalDof  != NULL)
    delete [] GlobalDofOFLocalDof;
 if(GlobalDofOFLocalDofAllRank  != NULL)
    delete [] GlobalDofOFLocalDofAllRank;
 if(N_LocalDofAllRank  != NULL)
    delete [] N_LocalDofAllRank;
}

#endif
