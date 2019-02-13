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
// @(#)ParFEMapper3D.C
//
// Class:      TParFEMapper3D
// Purpose:    Class containing all info needed for communication between subdomains
//
// Author:     Sashikumaar Ganesan (24.04.15)
//
// History:    Start of implementation 24.04.15 (Sashikumaar Ganesan)
//
// ======================================================================= 
#ifdef _MPI

#include "mpi.h"
#include <ParFEMapper3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <SubDomainJoint.h>
#include <Edge.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#define GLOBAL_NO 0
#define DOF_NO 1

#define HALOCELL 0
#define NONHALO  1

extern double timeC;

TParFEMapper3D::TParFEMapper3D(int N_dim, TFESpace3D *fespace, int *rowptr, int *kcol)
{
  N_Dim       = N_dim; 
  Comm        = TDatabase::ParamDB->Comm;
  FESpace     = fespace;
  RowPtr      = rowptr;
  KCol        = kcol;
 
  N_Dof = FESpace->GetN_DegreesOfFreedom();
  
  MaxSubDomainPerDof = fespace->GetMaxSubDomainPerDof();
  if(MaxSubDomainPerDof<0)
  {
    printf("Error: SetMaxSubDomainPerDof in FeSpace before calling ParFECommunicator3D \n");
    MPI_Finalize();
    exit(0);
  }
  
  if(TDatabase::ParamDB->MapperType != 2)
  {
    ConstructDofMap_Master_Halo();
#ifdef _HYBRID
    Color(N_CInt,ptrCInt,'i');
    Color(N_CMaster,ptrCMaster,'m');
    Color(N_CDept1,ptrCDept1,'D');
    Color(N_CDept2,ptrCDept2,'d');
#endif
  }
  else
  {
    ConstructDofMap();
  }
  
//   if(TDatabase::ParamDB->SOLVER_TYPE == DIRECT)
    Assign_GlobalDofNo();

}

static int GetLocalIndex(int N, int *array, int val)
{
 int m=0;
 while(array[m] != val)
    m++;
 return m;
}

int TParFEMapper3D::find_min(int *arr, int N, char *temp_arr){
  int i,j;
  int min = 1000000; 
  j = 0;
  for(i=0;i<N;i++){
    if(temp_arr[i] == 'x')   continue;
    if(arr[i]<min){
      min = arr[i];
      j = i;
    }
  }
  if(min == 1000000){
    printf("\n................this shudnt happen...................\n");
    MPI_Finalize();
    exit(0);
  }
  return j;
}

void TParFEMapper3D::ConstructDofMap_Master_Halo()
{ 
  int aa, bb, i, ii, j, jj, k, kk, l, ll, m, mm;
  int M, N, N_Vert, ID;
  int *DOF, *GlobalNumbers, *BeginIndex, *JointDof, *LocalIndex;
  int N_Cells, N_OwnCells, N_Dof, N_Active, N_LocDof;
  int temp,temp_globalno,temp_dofno;
 
  double start_time, end_time, temp_time;
  bool UPDATE;

  TCollection *Coll;
  TBaseCell *cell;

  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  start_time = MPI_Wtime();
  temp_time  = start_time;
  
  Coll       = FESpace->GetCollection();
  N_Cells    = Coll->GetN_Cells();
  N_Dof      = FESpace->GetN_DegreesOfFreedom();
  N_Active   = FESpace->GetN_ActiveDegrees();
  N_OwnCells = Coll->GetN_OwnCells();

  for(i=0;i<N_OwnCells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }
  for(i=N_OwnCells;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(!cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }

  /** *************************************************************************************/
  /** ARRAY CONTAINING GLOBAL DOF numbers (GLOBAL means LOCALLY GLOBAL for the subdomain) */ 
  /** BeginIndex: START ID IN GlobalNumbers ARRAY Cellwise*********************************/
  BeginIndex    = FESpace->GetBeginIndex();		
  GlobalNumbers = FESpace->GetGlobalNumbers();          

  /** Array containing global number of all Local cells */
  LocalIndex = new int[N_Cells]; 
  
  int **MappingData;
  MappingData = new int*[2];
  for(i=0;i<2;i++)
  {
    MappingData[i] = new int[N_Dof];
    memset(MappingData[i], 0, N_Dof*SizeOfInt);
  }
  
  /** *************************************************************************************/
  /** Master   :: Array containing the rank to which the Dof belongs                      */ 
  /** DofMarker:: Array to help finding Master, Slave, Dependent, Halo                    */
  /** DofMarker:: ALL the Dofs are marked 'i'                                             */
  /** *************************************************************************************/
  Master = new int[N_Dof];
  for(i=0;i<N_Dof;i++) 
    Master[i] = rank;
    
  DofMarker = new char[N_Dof];
  memset(DofMarker, 'i', N_Dof*sizeof(char));
  
  /** *************************************************************************************/
  /** Local Index Array is filled with GlobalNumbers of each cell                         */ 
  /** DofMarker::Dofs of Dependent cells are marked 'd'                                   */
  /** *************************************************************************************/
  for(i=0; i<N_Cells; i++)
  {
     cell          = Coll->GetCell(i);
     LocalIndex[i] = cell->GetGlobalCellNo();
     
     if(cell->IsDependentCell() && !cell->IsHaloCell())
     {
      DOF      = GlobalNumbers   + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N            = DOF[j];
	DofMarker[N] = 'd';
      }  // for(j=0; j<N_DOF; j++)   
     }
  } //  for(i=0; i<N_Cells
    
  /** *************************************************************************************/
  /** MappingData for dofs to be updated across domains are stored                        */ 
  /** DofMarker::Dofs of Halo cells not a part of dependent cells are marked 'h'          */
  /** *************************************************************************************/  
  for(i=N_OwnCells;i<N_Cells;i++)
  {
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();
     DOF = GlobalNumbers + BeginIndex[i];
     N_LocDof = BeginIndex[i+1] - BeginIndex[i];

     for(j=0; j<N_LocDof; j++)
     {
       N = DOF[j];
       if(DofMarker[N] == 'i') 
       {
	 DofMarker[N] = 'h';
	 Master[N]    = ID;
      }
      
      if(ID < Master[N]) Master[N] = ID;
      MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
      MappingData[DOF_NO][N]       = j; 
    } 
  }
  
  //------------------------------------------//
  /** Master DOF verification by other ranks **/
  //------------------------------------------//
 
  int *N_DOFtobeverified_otherRank;		//Array containing how many DOF's need to be verified from other processors(ranks) 
  int Total_DOFtobeverified_otherRank=0;		//Array containing how many DOF's need to be verified by this rank for other processors(ranks)
  int *N_DOFtobeverified_thisRank;		//Total no of DOF's that need to be verified from other ranks
  int Total_DOFtobeverified_thisRank=0;		//Total no of DOF verified by this rank
  int **sendbuf;
  int **recvbuf;
  int *verrecvbuf;
  int *versendbuf;
  int *temp_arr;
 
  sendbuf  = new int*[2];
  recvbuf  = new int*[2];
  sdispl   = new int[size];
  rdispl   = new int[size];
  temp_arr = new int[size];
  
  N_DOFtobeverified_otherRank = new int[size];
  N_DOFtobeverified_thisRank  = new int[size];
  memset(N_DOFtobeverified_otherRank,0,size*SizeOfInt);
  memset(N_DOFtobeverified_thisRank, 0,size*SizeOfInt); 
  
   //only the dofs in HALO cells (not dependent cells) needs to be verified
  for(i=0;i<N_Dof;i++){
    if(DofMarker[i] == 'h'){
      N_DOFtobeverified_otherRank[Master[i]]++;
      Total_DOFtobeverified_otherRank++;
    }
  }
  
  MPI_Alltoall(N_DOFtobeverified_otherRank,1, MPI_INT,N_DOFtobeverified_thisRank ,1, MPI_INT, Comm);
  for(i=0;i<size;i++)
    Total_DOFtobeverified_thisRank +=  N_DOFtobeverified_thisRank[i];
  
   for(i=0;i<2;i++){
     sendbuf[i] = new int[Total_DOFtobeverified_otherRank];
     recvbuf[i] = new int[Total_DOFtobeverified_thisRank];
   }

  //Send dof info to other processors and verify for master
  sdispl[0] = 0;
  for(i=1;i<size;i++){
    sdispl[i] = N_DOFtobeverified_otherRank[i-1] + sdispl[i-1];
  }
  
  rdispl[0] = 0;
  for(i=1;i<size;i++)
    rdispl[i] = N_DOFtobeverified_thisRank[i-1] + rdispl[i-1];
  
  memcpy (temp_arr,sdispl, size*SizeOfInt );
  
  for(i=0;i<N_Dof;i++){
     if(DofMarker[i] == 'h'){
       //if(rank==0) printf("%d ",temp_arr[Master[i]]);
       sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
       sendbuf[DOF_NO][temp_arr[Master[i]]]    = MappingData[DOF_NO][i];
       temp_arr[Master[i]]++;
     }
  }
  
  MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified_otherRank, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFtobeverified_thisRank, rdispl, MPI_INT, Comm);
  MPI_Alltoallv(sendbuf[DOF_NO],   N_DOFtobeverified_otherRank, sdispl, MPI_INT,recvbuf[DOF_NO],    N_DOFtobeverified_thisRank, rdispl, MPI_INT, Comm);
  
  verrecvbuf  	= new int[Total_DOFtobeverified_thisRank];
  versendbuf  	= new int[Total_DOFtobeverified_otherRank];

  for(i=0;i<Total_DOFtobeverified_thisRank;i++){
    temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
    verrecvbuf[i] = Master[GlobalNumbers[temp * N_LocDof + recvbuf[DOF_NO][i]]];
  }
  
  MPI_Alltoallv(verrecvbuf, N_DOFtobeverified_thisRank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified_otherRank, sdispl, MPI_INT, Comm);

  for(i=0;i<Total_DOFtobeverified_otherRank;i++){
    temp_globalno = sendbuf[GLOBAL_NO][i];
    temp_dofno    = sendbuf[DOF_NO][i];
    temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
    temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
    if(DofMarker[temp] != 'h')
      printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
    Master[temp] = versendbuf[i];
  }  

  //-----------------------------------------------------//
 /** Master DOF verification by other ranks completed  **/
 //-----------------------------------------------------//
  
  delete [] verrecvbuf;                       verrecvbuf = NULL;
  delete [] versendbuf;                       versendbuf = NULL; 
  delete [] N_DOFtobeverified_otherRank;      N_DOFtobeverified_otherRank = NULL;
  delete [] N_DOFtobeverified_thisRank;       N_DOFtobeverified_thisRank  = NULL;
  
  for(i=0;i<2;i++){
    delete [] sendbuf[i];                     sendbuf[i] = NULL;
    delete [] recvbuf[i];                     recvbuf[i] = NULL;
  }
  delete [] sendbuf;                          sendbuf = NULL;
  delete [] recvbuf;                          recvbuf = NULL;  
  
  end_time = MPI_Wtime();
  if(rank == 0)
    printf("total time taken for master verification = %lf\n",end_time-temp_time);
  temp_time = end_time;
//#################################################################  Master verification  ##################################################################################//   

//################################################################# Redistribution of interface dofs #######################################################################//
//*/
if(TDatabase::ParamDB->Par_P4){      
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//  
  int total_interface_dofs        = 0;                    //these are the dofs which lie on the interface of sub domains (total over all sub domains)
  int N_own_interface_dofs        = 0;                    //these are the own interface dofs
  int N_interface_dofs            = 0;                    //these are the interface dofs including shared ones
  int T_interface_dofs            = 0;                    //these are the total of N_interface_dofs over all sub domains
  int total_own_dofs              = 0;                    //these are the dofs which belongs to my rank (total implies total over all sub domains)
  int start = 0,           end    = 0;                    //temporary variables used for numbering the dofs
  int max_n_ranks_interface_dofs  = 0;                    //maximum number of ranks sharing any interface dof
  
  int *GlobalDofNo                = new int[N_Dof];       //this array is used to assign unique global dof no. over all sub domains 
      memset(GlobalDofNo, -1, N_Dof*sizeof(int));         //all is set to -1 for a default value
  int *N_Dof_Slave                = new int[size];        //array of number of slave interface dofs to be verified by other ranks
      memset(N_Dof_Slave, 0, size*sizeof(int));
  int *N_Dof_Master               = new int[size];        //array of interface dofs to be verified by this rank for other ranks 
  int **MasterPos                 = new int*[2];          //a double array for storing the position info for master interface dofs
  int **SlavePos                  = new int*[2];          //a double array for storing the position info for slave interface dofs 
  int *masterInfo, *slaveInfo;                            //these are used to update the global dof no for interface master-slave dofs
  int *GlobalDofNo_interface;                             //this array contains only the globaldof no. of interface dofs
  int *all_own_dofs_info          = new int[size];        //array of N_own dofs by each rank
  int *all_interface_dofs_info    = new int[size];        //array of N_own interface dofs by each rank
  int *all_T_interface_dofs_info  = new int[size];        //array of interface dofs(including shared ones) by each rank
  int *all_GlobalDofNo;                                   //array of global interface dof no from all ranks
  int *N_ranks_per_interface_dofs;                        //array of count of number of ranks sharing an interface dof
  int *N_allocated_masters        = new int[size];        //array containing number of masters(interface dofs) allocated to each processors
  char *tempc                     = new char[size];       //temporary array
  char **Master_Table             = new char*[size];      //a Table to mark the shared interface dofs

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//   
  
  //mark all (interface) dofs as 'z'
  for(i=N_OwnCells;i<N_Cells;i++){
    cell = Coll->GetCell(i);
    ID = cell->GetSubDomainNo();
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocDof = BeginIndex[i+1] - BeginIndex[i];
    
    //z is both interface and dependent dof, while d is only dependent dof
    for(j=0; j<N_LocDof; j++){
      N = DOF[j];
      if(DofMarker[N] == 'd'){
	DofMarker[N] = 'z';
      }
    }
  }
  
  //count total number of OWN interface_dofs and own dofs
  N_OwnDof = 0;
  for(i=0;i<N_Dof;i++){
    if(DofMarker[i] == 'z')
      if(Master[i] == rank)
        N_own_interface_dofs++;
    
    if(Master[i] == rank)
      N_OwnDof++;
  }
  
  //Send the count info to all ranks
  MPI_Allgather(&N_own_interface_dofs, 1, MPI_INT, all_interface_dofs_info, 1, MPI_INT, Comm);
  MPI_Allgather(&N_OwnDof,             1, MPI_INT, all_own_dofs_info,       1, MPI_INT, Comm);
    
  for(aa=0;aa<size;aa++){
    total_interface_dofs += all_interface_dofs_info[aa];
    total_own_dofs       += all_own_dofs_info[aa];
  }
  
  //start- numbering position for interface dofs
  aa=0;
  while(aa<rank){
    start += all_interface_dofs_info[aa];
    aa++;
  } 
  //printf("###########   rank = %d       start=%d\n",rank,start);
  
  Total_DOFtobeverified_otherRank = 0;                                               //this counts the total number of interface dofs(slave) that needs to be verified from other ranks(master) 
  Total_DOFtobeverified_thisRank  = 0;                                               //this counts the total number of interface dofs(masters) that needs to be verified for other ranks(slaves)
  
  //number the own interface dofs
  //count the total dofs to be verified by this and other rank
  for(i=0;i<N_Dof;i++){
    if(DofMarker[i]=='z'){
      if(Master[i]==rank){
	GlobalDofNo[i] = start;
	start++;
      }
      else{
	N_Dof_Slave[Master[i]]++;
	Total_DOFtobeverified_otherRank++;
      }
    }
  }
 
 MPI_Alltoall(N_Dof_Slave, 1, MPI_INT, N_Dof_Master, 1, MPI_INT, Comm);
 for(i=0;i<size;i++)
    Total_DOFtobeverified_thisRank +=  N_Dof_Master[i];
 
  for(i=0;i<2;i++){
   MasterPos[i] = new int[Total_DOFtobeverified_thisRank];
   SlavePos[i]  = new int[Total_DOFtobeverified_otherRank];
  }
  
  masterInfo  = new int[Total_DOFtobeverified_thisRank];
  slaveInfo   = new int[Total_DOFtobeverified_otherRank];
  
  rdispl[0] = 0;
  sdispl[0] = 0;
  for(i=1;i<size;i++){
    rdispl[i] = rdispl[i-1] + N_Dof_Master[i-1];
    sdispl[i] = sdispl[i-1] + N_Dof_Slave[i-1];                                 //slaves will be sent to gather info about their global dof no.
  }
 
  memcpy (temp_arr, sdispl, size*SizeOfInt );
 
  //store the position info of the slave dofs
  for(i=0;i<N_Dof;i++){
    if(DofMarker[i]=='z'){
      if(Master[i]!=rank){
       SlavePos[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
       SlavePos[DOF_NO][temp_arr[Master[i]]]    = MappingData[DOF_NO][i];
       temp_arr[Master[i]]++;
      }
    }
  }
  
//   for(aa=0;aa<size;aa++){
//     if(rank==aa){
//       printf("rank = %d\n",rank);	
//       for(i=0;i<size;i++){
// 	printf("N_Dof_Master[%d] = %d\t N_Dof_Slave[%d]=%d\n",i,N_Dof_Master[i],i,N_Dof_Slave[i]);
//       }
//     }
//     MPI_Barrier(Comm);
//   }
//   exit(0);

  //send the slave info to masters
  MPI_Alltoallv(SlavePos[GLOBAL_NO], N_Dof_Slave, sdispl, MPI_INT, MasterPos[GLOBAL_NO], N_Dof_Master, rdispl, MPI_INT, Comm);
  MPI_Alltoallv(SlavePos[DOF_NO],    N_Dof_Slave, sdispl, MPI_INT, MasterPos[DOF_NO],    N_Dof_Master, rdispl, MPI_INT, Comm);
 
  start = 0;
  for(aa=0;aa<size;aa++){
    end = start + N_Dof_Master[aa];
    for(i=start;i<end;i++){
      //get the cell with global cell no = MasterPos[GLOBAL_NO][i]
      temp = GetLocalIndex(N_Cells,LocalIndex,MasterPos[GLOBAL_NO][i]);
      
      //the master of this dof, i should be rank
      if(Master[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]] != rank){
	printf("....................Wrong Master.....................\n");
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should be a category z dof
      if(DofMarker[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]] != 'z'){
	printf("....................Wrong Dof.........................\n");
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should have a GlobalDofNo
      if(GlobalDofNo[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]] == -1){
	printf("....................Wrong GlobalDofNo.........................\n");
	MPI_Finalize();
	exit(0);
      }
      
      masterInfo[i] = GlobalDofNo[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]];
    }
    start = end;
  }

  MPI_Alltoallv(masterInfo, N_Dof_Master, rdispl, MPI_INT, slaveInfo, N_Dof_Slave, sdispl, MPI_INT, Comm);

  start = 0;
  for(aa=0;aa<size;aa++){
    end = start + N_Dof_Slave[aa];
    for(i=start;i<end;i++){
      temp_globalno = SlavePos[GLOBAL_NO][i];
      temp_dofno    = SlavePos[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp * N_LocDof + temp_dofno];
      
      //the master of this dof, i should be aa
      if(Master[temp] != aa){
	printf("...................2.Wrong Master.....master[%d]=%d....rank=%d............\n",temp,Master[temp],rank);
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should be a category z dof
      if(DofMarker[temp] != 'z'){
	printf("...................2.Wrong Dof.........................\n");
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should not have a GlobalDofNo
      if(GlobalDofNo[temp] != -1){
	printf("...................2.Wrong GlobalDofNo.........................\n");
	MPI_Finalize();
	exit(0);
      }
      //assign the global dof no to the slave obtained from master
      GlobalDofNo[temp] = slaveInfo[i];
    }
    start = end;
  }
  
  //now we mark the own dofs apart from the interface ones
  aa=0;
  start = total_interface_dofs;
  while(aa<rank){
    start += all_own_dofs_info[aa];
    aa++;
  }

  for(i=0;i<N_Dof;i++){
    if(DofMarker[i] != 'z'){
      if(Master[i]==rank){
	GlobalDofNo[i] = start;
	start++;
      }
    }
    
    if(Master[i]==rank){
      if(GlobalDofNo[i] == -1){
	printf("......................This shudnt happen.................\n");
	MPI_Finalize();
	exit(0);
      }
    }
    
    if(DofMarker[i]=='z'){
      if(GlobalDofNo[i]>=total_interface_dofs){
	printf("......................This shudnt happen.................\n");
	MPI_Finalize();
	exit(0);
      }
      
      N_interface_dofs++;
    }
  }
  
  GlobalDofNo_interface = new int[N_interface_dofs];
  aa = 0;
  for(i=0;i<N_Dof;i++){
     if(DofMarker[i]=='z'){
      GlobalDofNo_interface[aa] = GlobalDofNo[i];
      aa++;
     }
  }
  
//   MPI_Allreduce(&N_interface_dofs, &T_interface_dofs, 1, MPI_INT, MPI_SUM, Comm);
  
  delete [] N_Dof_Master;
  delete [] N_Dof_Slave;
  
  for(i=0;i<2;i++){
    delete [] MasterPos[i];                    MasterPos[i] = NULL;
    delete [] SlavePos[i];                     SlavePos[i]  = NULL;
  }
  delete [] MasterPos;                         MasterPos = NULL;
  delete [] SlavePos;                          SlavePos  = NULL;  
  
  delete [] masterInfo;                        masterInfo = NULL;
  delete [] slaveInfo;                         slaveInfo  = NULL;

  MPI_Allgather(&N_interface_dofs, 1, MPI_INT, all_T_interface_dofs_info, 1, MPI_INT, Comm);
  
  T_interface_dofs = 0;
  for(aa=0;aa<size;aa++)
    T_interface_dofs += all_T_interface_dofs_info[aa];
  
  all_GlobalDofNo = new int[T_interface_dofs];
  
//   for(aa=0;aa<size;aa++){
//     if(rank==aa){
//       printf("rank = %d      #total interface(own) dofs = %d       #interface dofs(own) = %d\n",rank,total_interface_dofs,N_own_interface_dofs);
//       printf("              #total interface dofs      = %d      #interface dofs      = %d\n",T_interface_dofs,N_interface_dofs);
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//   }
  
  rdispl[0] = 0;
  for(aa=1;aa<size;aa++)
    rdispl[aa] = rdispl[aa-1] + all_T_interface_dofs_info[aa-1];

  MPI_Allgatherv(GlobalDofNo_interface, N_interface_dofs, MPI_INT, all_GlobalDofNo, all_T_interface_dofs_info, rdispl, MPI_INT, Comm);

  //table is created only over the total own interface dofs over the sub domains
  for(i=0;i<size;i++){
    Master_Table[i] = new char[total_interface_dofs];
    memset(Master_Table[i], 'x', total_interface_dofs*sizeof(char));
  }
  
  start = 0;end = 0;
  for(i=0;i<size;i++){
    end += all_T_interface_dofs_info[i];
    for(aa=start;aa<end;aa++){
      Master_Table[i][all_GlobalDofNo[aa]] = 'y';
    }
    start = end; 
  }
  
  N_ranks_per_interface_dofs = new int[total_interface_dofs];
  memset(N_ranks_per_interface_dofs, 0, total_interface_dofs*sizeof(int));
  
  for(j=0;j<total_interface_dofs;j++){
    //check how many ranks share the interface dof
    for(i=0;i<size;i++){
      if(Master_Table[i][j] == 'y'){
	N_ranks_per_interface_dofs[j]++;
      }
    }
    //compute the max munber of rank sharing any interface dof
    if(max_n_ranks_interface_dofs<N_ranks_per_interface_dofs[j])
      max_n_ranks_interface_dofs = N_ranks_per_interface_dofs[j];
  }
  
//   if(rank==1)
//   for(i=0;i<size;i++){
//     for(aa=0;aa<20;aa++){
//       printf("%c\t",Master_Table[i][aa]);
//     }
//     printf("\n");
//   }
  
  for(j=0;j<total_interface_dofs;j++){
  //a interface dof should have at least one neighbour from other rank
    if(N_ranks_per_interface_dofs[j] == 1){
	printf("......................This shudnt happen.....j=%d............\n",j);
	MPI_Finalize();
	exit(0);
    }
  }
 
  memset(N_allocated_masters, 0, size*sizeof(int));
  for(i=0; i<N_Dof; i++){
    if(DofMarker[i]=='z'){
      N_allocated_masters[Master[i]]++;
    }
  }
  
  
//   for(aa=0;aa<size;aa++){
//      if(rank==aa){
//        printf("\nrank=%d\n",rank);
//        for(i=0;i<size;i++)
// 	 printf("Before Redistribution::             N_allocated_masters[%d] = %d\n",i,N_allocated_masters[i]);
//      }
//      MPI_Barrier(Comm);
//      //break;
//   }
  
  memset(N_allocated_masters, 0, size*sizeof(int));
  int min = 0;
  //a interface dof should have at least one neighbour from other rank
  for(i=2; i<=max_n_ranks_interface_dofs;i++){
    for(j=0;j<total_interface_dofs;j++){
      if(N_ranks_per_interface_dofs[j] == i){
	
	for(aa=0;aa<size;aa++){
	  tempc[aa] = Master_Table[aa][j];
	}
	
	min = find_min(N_allocated_masters,size,tempc);
	N_allocated_masters[min]++;
	
	for(aa=0;aa<N_Dof;aa++){
	  if(GlobalDofNo[aa] == j){
	    Master[aa] = min;
	  }
	}
	N_ranks_per_interface_dofs[j] = -1;
      }
    }
  }
  
//   for(aa=0;aa<size;aa++){
//      if(rank==aa){
//        printf("\nrank=%d\n",rank);
//        for(i=0;i<size;i++)
// 	 printf("After Redistribution::             N_allocated_masters[%d] = %d\n",i,N_allocated_masters[i]);
//      }
//      MPI_Barrier(Comm);
//   }
  
  memset(N_allocated_masters, 0, size*sizeof(int));
  for(i=0; i<N_Dof; i++){
    if(DofMarker[i]=='z'){
      N_allocated_masters[Master[i]]++;
      
      if(GlobalDofNo[i] == -1){
	printf("......................This shudnt happen.................\n");
	MPI_Finalize();
	exit(0);
      }
    }
  }
  
//   for(aa=0;aa<size;aa++){
//      if(rank==aa){
//        printf("\nrank=%d\n",rank);
//        for(i=0;i<size;i++)
// 	 printf("N_allocated_masters[%d] = %d\n",i,N_allocated_masters[i]);
//      }
//      MPI_Barrier(Comm);
//      //break;
//   }
  
  delete [] all_GlobalDofNo;                  all_GlobalDofNo            = NULL;
  delete [] GlobalDofNo;                      GlobalDofNo                = NULL;
  delete [] GlobalDofNo_interface;            GlobalDofNo_interface      = NULL;
  for(i=0;i<size;i++) 
    delete [] Master_Table[i];                Master_Table[i]            = NULL;
  delete [] N_ranks_per_interface_dofs;       N_ranks_per_interface_dofs = NULL;     
  delete [] N_allocated_masters;              N_allocated_masters        = NULL;
  delete [] temp_arr;                         temp_arr                   = NULL;
  
  end_time = MPI_Wtime();
  if(rank == 0)
    printf("Total Time Taken for Redistribution of master dofs = %lf\n",end_time-temp_time);
  temp_time = end_time;

    } //if(TDatabase::ParamDB->Par_P4)  
  //*/

//################################################################# Redistribution of interface dofs #######################################################################//


  
//#################################################################  Marking The Dofs ######################################################################################//     
  //------------------------------------------------------------------//
  /** Gather information about other DOFs ::                          */
  /**    ## dofs marked as 'i'             -->independent             */
  /**    ## dofs marked as 's'             -->slave                   */
  /**    ## dofs marked as 'D'             -->dependent_type1         */
  /**    ## dofs marked as 'd'             -->dependent_type2         */
  /**    ## dofs marked as 'H'             -->halo_type1              */
  /**    ## dofs marked as 'h'             -->halo_type2              */
  //------------------------------------------------------------------// 
  
  memset(DofMarker, 'i', N_Dof*sizeof(char));		//all dofs marked as independent
  //all dofs in halo cell marked as halo_type2(unused)
  for(i=N_OwnCells;i<N_Cells;i++){
    cell     = Coll->GetCell(i);
    DOF      = GlobalNumbers + BeginIndex[i];
    N_LocDof = BeginIndex[i+1] - BeginIndex[i];
   
    for(j=0;j<N_LocDof;j++){
      N            = DOF[j];
      DofMarker[N] = 'h';
    }
  }
  
  for(i=0;i<N_OwnCells;i++){
    cell = Coll->GetCell(i);
    //now mark dofs in dependent cells
    if(cell->IsDependentCell()){
      DOF      = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];
     
      for(j=0;j<N_LocDof;j++){
        N = DOF[j];
        if(DofMarker[N] != 'h')	//if dof is not marked halo then it is dependent dof
	  DofMarker[N] = 'd';
      }
    }
  }
  
  for(i=0;i<N_OwnCells;i++){
    cell = Coll->GetCell(i);
    //now mark dofs in dependent cells
    if(cell->IsDependentCell()){
      DOF      = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];
     
      for(j=0;j<N_LocDof;j++){
        N = DOF[j];
        if(DofMarker[N] == 'h'){
	  if(Master[N] != rank)
	    DofMarker[N] = 's';	//if dof is marked halo and master of it is some other proc then it is a slave
	 else
	    DofMarker[N] = 'm';	//else it is a master
        }
      }
    }
  }
  
  //mark dependent type1 & type2
  bool flag = false;
  for(i=0;i<N_OwnCells;i++){
    cell = Coll->GetCell(i);
    //now mark dofs in dependent cells
    if(cell->IsDependentCell()){
      DOF      = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];
     
      for(j=0;j<N_LocDof;j++){
        N = DOF[j];
        if(Master[N] != rank){
	  flag = true;
	  break;
        }
      }
     //dependent dofs connected to slave dofs are marked as type1(D) else type2(d)
     //type1 dofs must be smoothed first as they are the halo dofs (type1 i.e. useful) to other procs
     if(flag == true){
       for(j=0;j<N_LocDof;j++){
	 N = DOF[j];
	 if(DofMarker[N] == 'd')
	   DofMarker[N] = 'D';
         }
         flag = false;
      }
    }
  }

#ifdef _HYBRID
if(TDatabase::ParamDB->Par_P5 == 1)
{  
   N_InterfaceM = 0;
   N_Halo1 = 0;
   N_Dept2 = 0;
   N_Int = 0;
  
  for(i=0;i<N_Dof;i++)
  {
    if(DofMarker[i] == 'H')
      N_Halo1++;
    else if(DofMarker[i] == 'd')
      N_Dept2++;
    else if(DofMarker[i] == 'i')
      N_Int++;
    else if(DofMarker[i] == 'm')
      N_InterfaceM++;
    else if(DofMarker[i] == 's')
      N_InterfaceS++;
    else if(DofMarker[i] == 'D')
      N_Dept1++;
  }
  
  double reqd_ratio = (N_InterfaceM+N_InterfaceS)/(N_Dept1+N_Halo1);;
  double ratio;
  int ctr = 0;
  
  N_Dept2 = 0;
  N_Int = 0;
  for(i=0;i<N_Dof;i++)
  {
    if(DofMarker[i] == 'd')
      N_Dept2++;
    else if(DofMarker[i] == 'i')
      N_Int++;
  }
  
  if(N_Dept2!=0)
    ratio = N_Int/N_Dept2;
  else
    ratio=0;
  
  while(ratio<reqd_ratio && ctr<4)
  {
    ctr++;
      flag = false;
      for(i=0;i<N_OwnCells;i++){
	cell = Coll->GetCell(i);
	//now mark dofs in independent cells
	if( !(cell->IsDependentCell()) ){
	  DOF      = GlobalNumbers + BeginIndex[i];
	  N_LocDof = BeginIndex[i+1] - BeginIndex[i];
	
	  for(j=0;j<N_LocDof;j++){
	    N = DOF[j];
	    if(DofMarker[N] == 'd' || DofMarker[N] == 'D'){
	      flag = true;
	      break;
	    }
	  }
	//dependent dofs connected to slave dofs are marked as type1(D) else type2(d)
	//here we increase the ratio of dependent_type2 dofs to independent dofs
	if(flag == true){
	  for(j=0;j<N_LocDof;j++){
	    N = DOF[j];
	    if(DofMarker[N] == 'i')
	      DofMarker[N] = 't';
	    }
	    flag = false;
	  }
	}
      }
      
      for(i=0;i<N_Dof;i++)
      {
	if(DofMarker[i] == 't')
	  DofMarker[i] = 'd';
      }
    
      N_Dept2 = 0;
      N_Int = 0;
      for(i=0;i<N_Dof;i++)
      {
	if(DofMarker[i] == 'd')
	  N_Dept2++;
	else if(DofMarker[i] == 'i')
	  N_Int++;
      }
      if(N_Dept2!=0)
	ratio = N_Int/N_Dept2;
      else
	ratio=0;
//       ratio = N_Int/N_Dept2;
  }
}
#endif
  
   //mark halo type1(H)-->useful & type2(h)
  flag = false;
  for(i=N_OwnCells;i<N_Cells;i++){
    DOF      = GlobalNumbers + BeginIndex[i];
    N_LocDof = BeginIndex[i+1] - BeginIndex[i];
   
    for(j=0;j<N_LocDof;j++){
      N = DOF[j];
      if(Master[N]==rank){
        flag = true;
        break;
      }
    }
   //halo dofs connected to master dofs are marked as type1(H) else type2(h)
   //type2 halo dofs are not required in smoothing operations
    if(flag == true){
      for(j=0;j<N_LocDof;j++){
        N = DOF[j];
        if(DofMarker[N] == 'h')
	  DofMarker[N] = 'H';  
      }
      flag = false;
    }
  }
  
//   for(i=0;i<N_Dof;i++)
//   {
//     if(DofMarker[i]!='m')	continue;
//     
//     for(j=RowPtr[i];j<RowPtr[i+1];j++)
//     {
//       if(DofMarker[KCol[j]] == 'h')
//       {
// 	DofMarker[KCol[j]] = 'H';
//       }
//     }
//   }

  end_time = MPI_Wtime();
  if(rank == 0)
    printf("Total Time Taken for marking the dofs = %lf\n",end_time-temp_time);
  temp_time = end_time;
//#################################################################  Marking The Dofs ######################################################################################// 
  
//#############################################################  Mapper for Master Dof are set ##############################################################################//  
  int **SlaveBuf,**MasterBuf;
  int *temp_arrH1,*temp_arrH2;
 
  N_InterfaceM = 0;      N_InterfaceS  = 0;
  N_Slave      = 0;     
  N_OwnDof     = 0; 
  N_Master     = 0;
  N_Int        = 0; 
  N_Dept       = 0;      N_Dept1 = 0;    N_Dept2 = 0;    //N_Dept3 = 0;
  N_Halo       = 0;      N_Halo1 = 0;    N_Halo2 = 0;    
 
  N_DofSend    = new int[size];
  N_DofSendMS  = new int[size]; 
  N_DofSendH2  = new int[size];
  N_DofSendH1  = new int[size];
 
  N_DofRecv    = new int[size];
  N_DofRecvMS  = new int[size]; 
  N_DofRecvH2  = new int[size];
  N_DofRecvH1  = new int[size];
  
  memset(N_DofRecvMS  ,0,size*SizeOfInt);
  memset(N_DofRecvH1,0,size*SizeOfInt);
  memset(N_DofRecvH2,0,size*SizeOfInt);
  
  temp_arr   = new int[size];
  temp_arrH1 = new int[size];
  temp_arrH2 = new int[size];
  
  SlaveBuf  = new int*[2];
  MasterBuf = new int*[2]; 

  for(N=0;N<N_Dof;N++){
   
    if(Master[N] != rank){
      N_Slave++;
      N_DofRecv[Master[N]]++;
     
      if(DofMarker[N] == 's'){
        N_InterfaceS++;
        N_DofRecvMS[Master[N]]++;
      }
      else if(DofMarker[N] == 'H'){
        N_Halo1++;       
        N_DofRecvH1[Master[N]]++;
      }
      else{
        N_Halo2++;
        N_DofRecvH2[Master[N]]++;
      }
    }
    else{
      N_OwnDof++;
      N_Master++;
      if(DofMarker[N] == 'm')
        N_InterfaceM++;
      else if(DofMarker[N] == 'D')
        N_Dept1++;
      else if(DofMarker[N] == 'd')
        N_Dept2++;
    // else if(DofMarker[N] == 'x')
       //N_Dept3++;
      else
        N_Int++;
    }
  }
 
  N_Halo = N_Halo1 + N_Halo2;
  N_Dept = N_Dept1 + N_Dept2;// + N_Dept3;
  
  //DEBUG
  if(size<9)
  for(aa=0;aa<size;aa++){
     if(rank==aa){
          printf("\nRank::%d\n",rank);
	  printf("N_Dof               = %d\t N_OwnDof           = %d\t N_Active           = %d\n", N_Dof, N_OwnDof,N_Active);
	  printf("N_Master(Total)     = %d\t N_Slave(Total)     = %d\n", N_Master, N_Slave);
	  printf("N_Master(Interface) = %d\t N_Slave(Interface) = %d\n", N_InterfaceM, N_InterfaceS);
	  printf("N_Halo              = %d\t N_Halo1            = %d\t   N_Halo2  = %d\n",N_Halo,N_Halo1,N_Halo2);
	  printf("N_Depndt            = %d\t N_Dept1            = %d\t   N_Dept2  = %d\n",N_Dept,N_Dept1,N_Dept2);
	  printf("N_Indpdt            = %d\n",N_Int); 
     }
     MPI_Barrier(MPI_COMM_WORLD);
  }//verified
  
  rdispl[0] = 0; sdispl[0] = 0; 
  sdisplMS  = new int[size];   sdisplMS[0] = 0;  
  rdisplMS  = new int[size];   rdisplMS[0] = 0;
  sdisplH1  = new int[size];   sdisplH1[0] = 0;  
  rdisplH1  = new int[size];   rdisplH1[0] = 0;
  sdisplH2  = new int[size];   sdisplH2[0] = 0;
  rdisplH2  = new int[size];   rdisplH2[0] = 0;
  
  OwnDofs = new int[N_OwnDof];
  
  N_SendDofMS = 0; 
  N_SendDofH1 = 0; 
  N_SendDofH2 = 0;  

  MPI_Alltoall(N_DofRecv  , 1, MPI_INT,N_DofSend   , 1, MPI_INT, Comm);
  MPI_Alltoall(N_DofRecvMS, 1, MPI_INT,N_DofSendMS , 1, MPI_INT, Comm);
  MPI_Alltoall(N_DofRecvH1, 1, MPI_INT,N_DofSendH1 , 1, MPI_INT, Comm);
  MPI_Alltoall(N_DofRecvH2, 1, MPI_INT,N_DofSendH2 , 1, MPI_INT, Comm);
  
  for(i=1;i<size;i++){
    rdispl[i]   = rdispl[i-1]   + N_DofRecv[i-1];
    rdisplMS[i] = rdisplMS[i-1] + N_DofRecvMS[i-1];
    rdisplH1[i] = rdisplH1[i-1] + N_DofRecvH1[i-1];
    rdisplH2[i] = rdisplH2[i-1] + N_DofRecvH2[i-1];
   
    sdispl[i]   = sdispl[i-1]   + N_DofSend[i-1];
    sdisplMS[i] = sdisplMS[i-1] + N_DofSendMS[i-1];
    sdisplH1[i] = sdisplH1[i-1] + N_DofSendH1[i-1];
    sdisplH2[i] = sdisplH2[i-1] + N_DofSendH2[i-1];
  }
 
  for(i=0;i<size;i++){
    N_SendDofMS += N_DofSendMS[i];
    N_SendDofH1 += N_DofSendH1[i];
    N_SendDofH2 += N_DofSendH2[i];
  }
  N_SendDof = N_SendDofMS + N_SendDofH1 + N_SendDofH2;

  memcpy (temp_arr,   rdisplMS, size*SizeOfInt );
  memcpy (temp_arrH1, rdisplH1, size*SizeOfInt );
  memcpy (temp_arrH2, rdisplH2, size*SizeOfInt );
 
  DofSend    = new int[N_SendDofMS+N_SendDofH1+N_SendDofH2];
  DofSendMS  = DofSend;
  DofSendH1  = DofSend + N_SendDofMS;
  DofSendH2  = DofSend + N_SendDofMS + N_SendDofH1;
  
  DofRecv    = new int[N_InterfaceS+N_Halo1+N_Halo2];
  DofRecvMS  = DofRecv;
  DofRecvH1  = DofRecv + N_InterfaceS;
  DofRecvH2  = DofRecv + N_InterfaceS + N_Halo1;
  
  for(i=0;i<2;i++){
    if(N_InterfaceS>0)
      SlaveBuf[i]  = new int[N_InterfaceS];
    if(N_SendDofMS>0)
      MasterBuf[i] = new int[N_SendDofMS];
    //memset (SlaveBuf[i] , 0, size*SizeOfInt);
    //memset (MasterBuf[i], 0, size*SizeOfInt);
  } 
  
  Reorder = new int[N_Dof];
  NewGN = new int[N_Dof];
  m = 0;
 
  int Mstr  = 0;                  Reorder_M  = Reorder;
  int Indpt = N_InterfaceM;       Reorder_I  = Reorder + Indpt;
  int Dept1 = Indpt + N_Int;      Reorder_D1 = Reorder + Dept1;
  int Dept2 = Dept1 + N_Dept1;    Reorder_D2 = Reorder + Dept2;
  //int Dept3 = Dept2 + N_Dept2;    Reorder_D3 = Reorder + Dept3;
  int Slv   = Dept2 + N_Dept2;
  int Hl1   = Slv   + N_InterfaceS;
  int Hl2   = Hl1   + N_Halo1;
  int ts = 0, th1 = 0, th2 = 0, ti = 0, tm = 0, td1 = 0, td2 = 0;// td3 = 0;
  
  for(i=0;i<N_Dof;i++){
   //slave dofs
   if(Master[i] != rank){
     if(DofMarker[i] == 's'){
       SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
       SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
   
//        if(SlaveBuf[0][temp_arr[Master[i]]]<0 || Master[i]>=size || Master[i]<0 || temp_arr[Master[i]]>N_Dof){
// 	printf("\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"); 
// 	printf("\n.................This shudnt happen!!!(rank=%d,cell=%d,dof=%d)\n",rank,MappingData[GLOBAL_NO][i],i);
// 	exit(0);
//        }
   
       DofRecvMS[temp_arr[Master[i]]] = i;
       temp_arr[Master[i]]++;
       
       {
         Reorder[Slv] = i;
	 NewGN[i]     = ts++;//Slv;
         Slv++;
       }
     }
     else if(DofMarker[i] == 'H'){
       {
	 Reorder[Hl1] = i;
	 NewGN[i]     = th1++;//Hl1; 
	 Hl1++;
       }
     }
     else{
       {
	 Reorder[Hl2] = i;
	 NewGN[i]     = th2++;//Hl2;
	 Hl2++;
       }
     }
   }
   else{
     OwnDofs[m++]=i;
       {
       if(DofMarker[i]=='m'){
	 Reorder[Mstr] = i;
	 NewGN[i]      = tm++;//Mstr;
	 Mstr++;
       }
       else if(DofMarker[i]=='D'){
	 Reorder[Dept1] = i;
	 NewGN[i]       = td1++;//Dept1;
	 Dept1++;
       }
       else if(DofMarker[i]=='d'){
	 Reorder[Dept2] = i;
	 NewGN[i]       = td2++;//Dept2;
	 Dept2++;
       }
       //else if(DofMarker[i]=='x'){
	 //Reorder[Dept3] = i;
	 //NewGN[i]       = td3++;//Dept3;
	 //Dept3++;
       //}
       else{
	 Reorder[Indpt] = i;
	 NewGN[i]       = ti++;//Indpt;
	 Indpt++;
       }
     }
   }
 }
 
  MPI_Alltoallv(SlaveBuf[GLOBAL_NO], N_DofRecvMS, rdisplMS, MPI_INT, MasterBuf[GLOBAL_NO], N_DofSendMS, sdisplMS, MPI_INT, Comm);	
  MPI_Alltoallv(SlaveBuf[DOF_NO],    N_DofRecvMS, rdisplMS, MPI_INT, MasterBuf[DOF_NO],    N_DofSendMS, sdisplMS, MPI_INT, Comm);
 
  for(i=0;i<N_SendDofMS;i++){
    temp_globalno = MasterBuf[GLOBAL_NO][i];
    temp_dofno    = MasterBuf[DOF_NO][i];
    temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
    temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
    DofSendMS[i] = temp; 
  }

  for(i=0;i<2;i++){
    if(N_InterfaceS>0){
      delete [] SlaveBuf[i];
      SlaveBuf[i] = NULL;
    }
    if(N_SendDofMS>0){
      delete [] MasterBuf[i];
      MasterBuf[i] = NULL;
    }
  }
  delete [] SlaveBuf;    SlaveBuf  = NULL;
  delete [] MasterBuf;   MasterBuf = NULL;
  delete [] temp_arr;    temp_arr  = NULL;  
  
 //*********************************************************************************************************************************//
 //                                                  Mapper for Halo Dof are set                                              //
 //*********************************************************************************************************************************//
  int **SlaveBufH1,**MasterBufH1;
  int **SlaveBufH2,**MasterBufH2;
 
  SlaveBufH1  = new int*[2];
  MasterBufH1 = new int*[2];
 
  SlaveBufH2  = new int*[2];
  MasterBufH2 = new int*[2];
 
  for(i=0;i<2;i++){
    if(N_Halo1>0)
      SlaveBufH1[i]  = new int[N_Halo1];
    if(N_SendDofH1>0)
      MasterBufH1[i] = new int[N_SendDofH1];
     //memset (SlaveBufH1[i] , 0, size*SizeOfInt);
     //memset (MasterBufH1[i], 0, size*SizeOfInt);
    if(N_Halo2>0)
      SlaveBufH2[i]  = new int[N_Halo2];
    if(N_SendDofH2>0)
      MasterBufH2[i] = new int[N_SendDofH2];
      //memset (SlaveBufH2[i] , 0, size*SizeOfInt);
      //memset (MasterBufH2[i], 0, size*SizeOfInt);
  }
  
  for(i=0;i<N_Dof;i++){
    if(DofMarker[i] == 'H'){
      if(Master[i] == rank){
        printf("\nThis should be a HALO_1 dof\n");
        exit(0);
      }
     
      SlaveBufH1[GLOBAL_NO][temp_arrH1[Master[i]]] = MappingData[GLOBAL_NO][i];
      SlaveBufH1[DOF_NO]   [temp_arrH1[Master[i]]] = MappingData[DOF_NO]   [i];
	  
      DofRecvH1[temp_arrH1[Master[i]]] = i;
      temp_arrH1[Master[i]]++;
    }
    else if(DofMarker[i] == 'h'){
      if(Master[i] == rank){
        printf("\nThis should be a HALO_2 dof\n");
        exit(0);
      }
     
      SlaveBufH2[GLOBAL_NO][temp_arrH2[Master[i]]] = MappingData[GLOBAL_NO][i];
      SlaveBufH2[DOF_NO]   [temp_arrH2[Master[i]]] = MappingData[DOF_NO]   [i];
	  
      DofRecvH2[temp_arrH2[Master[i]]] = i;
      temp_arrH2[Master[i]]++;
    }
  }
  
  MPI_Alltoallv(SlaveBufH1[GLOBAL_NO], N_DofRecvH1, rdisplH1, MPI_INT, MasterBufH1[GLOBAL_NO], N_DofSendH1, sdisplH1, MPI_INT, Comm);	
  MPI_Alltoallv(SlaveBufH1[DOF_NO],    N_DofRecvH1, rdisplH1, MPI_INT, MasterBufH1[DOF_NO],    N_DofSendH1, sdisplH1, MPI_INT, Comm);
 
  for(i=0;i<N_SendDofH1;i++){
    temp_globalno = MasterBufH1[GLOBAL_NO][i];
    temp_dofno    = MasterBufH1[DOF_NO][i];
    temp          = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
    temp          = GlobalNumbers[temp*N_LocDof + temp_dofno];
    DofSendH1[i]  = temp; 
  }
  
  for(i=0;i<2;i++){
    if(N_Halo1>0){
      delete [] SlaveBufH1[i];
      SlaveBufH1[i] = NULL;
    }
    if(N_SendDofH1>0){
      delete [] MasterBufH1[i];
      MasterBufH1[i] = NULL;
    }
  }
  delete [] SlaveBufH1;    SlaveBufH1  = NULL;
  delete [] MasterBufH1;   MasterBufH1 = NULL;
  delete [] temp_arrH1;    temp_arrH1  = NULL;
  
  MPI_Alltoallv(SlaveBufH2[GLOBAL_NO], N_DofRecvH2, rdisplH2, MPI_INT, MasterBufH2[GLOBAL_NO], N_DofSendH2, sdisplH2, MPI_INT, Comm);	
  MPI_Alltoallv(SlaveBufH2[DOF_NO],    N_DofRecvH2, rdisplH2, MPI_INT, MasterBufH2[DOF_NO],    N_DofSendH2, sdisplH2, MPI_INT, Comm);

  for(i=0;i<N_SendDofH2;i++){
    temp_globalno = MasterBufH2[GLOBAL_NO][i];
    temp_dofno    = MasterBufH2[DOF_NO][i];
    temp          = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
    temp          = GlobalNumbers[temp*N_LocDof + temp_dofno];
    DofSendH2[i]  = temp; 
  }
 
  for(i=0;i<2;i++){
    if(N_Halo2>0){
      delete [] SlaveBufH2[i];
      SlaveBufH2[i] = NULL;
    }
    if(N_SendDofH2>0){
       delete [] MasterBufH2[i];
       MasterBufH2[i] = NULL;
    }
  }
  delete [] SlaveBufH2;    SlaveBufH2  = NULL;
  delete [] MasterBufH2;   MasterBufH2 = NULL;
  delete [] temp_arrH2;    temp_arrH2  = NULL;
  
  if(N_Dim>1)
  {
     for(i=0;i<size;i++)
     {
        sdispl[i]   *= N_Dim;
	sdisplMS[i] *= N_Dim;
	sdisplH1[i] *= N_Dim;
	sdisplH2[i] *= N_Dim;
	
	rdispl[i]     *= N_Dim;
	rdisplMS[i]   *= N_Dim;
	rdisplH1[i]   *= N_Dim;
	rdisplH2[i]   *= N_Dim;
	
	N_DofSend[i]   *= N_Dim;
        N_DofSendMS[i] *= N_Dim;
        N_DofSendH1[i] *= N_Dim;
        N_DofSendH2[i] *= N_Dim;
     
        N_DofRecv[i]    *= N_Dim;
        N_DofRecvMS[i]  *= N_Dim;
        N_DofRecvH1[i]  *= N_Dim;
        N_DofRecvH2[i]  *= N_Dim;
     }
  }
	
  
  if(N_SendDof>0) Send_Info   = new double[N_SendDof*N_Dim];
  Send_InfoMS = Send_Info;
  Send_InfoH1 = Send_Info + N_SendDofMS*N_Dim;
  Send_InfoH2 = Send_Info + N_SendDofMS*N_Dim + N_SendDofH1*N_Dim;
 
  if(N_Slave>0) Recv_Info   = new double[N_Slave*N_Dim];
  Recv_InfoMS = Recv_Info;
  Recv_InfoH1 = Recv_Info + N_InterfaceS*N_Dim;
  Recv_InfoH2 = Recv_Info + N_InterfaceS*N_Dim + N_Halo1*N_Dim;
 
  end_time = MPI_Wtime();
  if(rank == 0)
    printf("Total Time Taken for mapping the dofs = %lf\n",end_time-temp_time);
  
  //if(TDatabase::ParamDB->SC_VERBOSE>2)
    if(rank==TDatabase::ParamDB->Par_P0)
      printf("\n################       Mapping for slave-master dofs and halo_1, halo_2 dofs done !!!    ################\n");
  
  if(rank == 0)
      printf("total time taken by the ConstructDofMap_light() = %lf\n",MPI_Wtime()-start_time);  

 //only for checking purpose   
 if(0){
	  int Total_Own = 0;   
	  for(i=0;i<N_Dof;i++)
	  {
	    if(Master[i]!=rank){
	      if(!(DofMarker[i]=='s' || DofMarker[i]=='h' || DofMarker[i]=='H')){
		printf("\n................1.This shudnt happen...................\n");
		MPI_Finalize();
		exit(0);
	      }
	    }
	    else if(Master[i]==rank){
	      Total_Own++;
	      if(!(DofMarker[i]=='m' || DofMarker[i]=='d' || DofMarker[i]=='D' || DofMarker[i]=='i' || DofMarker[i]=='x')){
		printf("\n................2.This shudnt happen...........DofMarker[%d]=%c........\n",i,DofMarker[i]);
		MPI_Finalize();
		exit(0);
	      }
	    }
	    else{
		printf("\n................3.This shudnt happen...................\n");
		MPI_Finalize();
		exit(0);
	    } 
	  }
	  
 }
//     MPI_Barrier(MPI_COMM_WORLD);
    
//  MPI_Finalize();   
//  exit(0);
 
}

void TParFEMapper3D::ConstructDofMap()
{
   
 int rank, size, i,jj,ii, j, k, l, m, m1, P, N_Cells, N_U,N_Active, N_LocDof, ID;
 int M, N, N_Vert, N_Joints, N_JointDOF, Neib_ID;
 int *DOF, *GlobalNumbers, *BeginIndex, *JointDof, Disp, N_EdgeDOF, *EdgeDof,*LocalIndex,*Verify;
 int N_CrossEdgeNeibs, *CrossEdgeNeibsRank, N_Edges, N_VertInCell,N_OwnCells;
 int N_VertCrossNeibs, *VertCrossNeibs, VertDof;
 int test, N_Dof;

 double x,y,z;

 bool UPDATE;

 TCollection *Coll;
 TBaseCell *cell, *Neib_cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = FESpace->GetN_DegreesOfFreedom();
  N_Active = FESpace->GetN_ActiveDegrees();
  N_OwnCells = Coll->GetN_OwnCells();
  
 // printf("N_U is %d----Rank %d\n",N_U,rank);
  
  
  for(i=0;i<N_OwnCells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }
  for(i=N_OwnCells;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(!cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }

  /** *************************************************************************************/
  /** ARRAY CONTAINING GLOBAL DOF numbers (GLOBAL means LOCALLY GLOBAL for the subdomain) */ 
  /** BeginIndex: START ID IN GlobalNumbers ARRAY Cellwise*********************************/
  BeginIndex = FESpace->GetBeginIndex();		
  GlobalNumbers = FESpace->GetGlobalNumbers();          
   
  /** Array containing number of ranks(process ID) associated (surrounding) with each DOF (including halo DOF) */ 
//   N_DofRankIndex = new int[N_U];
 
//   memset(N_DofRankIndex, 1, N_U*SizeOfInt);
  /** Array containing global number of all Local cells */
  LocalIndex = new int[N_Cells]; 
  
  int **MappingData;
  MappingData = new int*[2];
  for(i=0;i<2;i++)
  {
    MappingData[i] = new int[N_U];
    memset(MappingData[i], 0, N_U*SizeOfInt);
  }
  
  Master = new int[N_U];
  for(i=0;i<N_U;i++) Master[i]=rank;
    
  Verify = new int[N_U];
  memset(Verify, 0, N_U*SizeOfInt);
  
 /** START ---> [ DofRankIndex ------- N_DofRankIndex -------- N_DependentCells ] **/ 
 for(i=0; i<N_Cells; i++)
  {
     cell = Coll->GetCell(i);
     LocalIndex[i] = cell->GetGlobalCellNo();
     
     if(cell->IsDependentCell() && !cell->IsHaloCell())
     {
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	Verify[N] = -1;
      }  // for(j=0; j<N_DOF; j++)
      
     }
  } //  for(i=0; i<N_Cells
    

 for(i=N_OwnCells;i<N_Cells;i++)
 {
      cell = Coll->GetCell(i);
      ID = cell->GetSubDomainNo();
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	//if(N>=N_Active) continue;
	if(Verify[N] == 0) 
	{
	  Verify[N] = 1;
	  Master[N] = ID;
	}
	if(ID < Master[N]) Master[N] = ID;
	MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
	MappingData[DOF_NO][N]       = j; 
      }
      
  }
 
   
    int temp,temp_globalno,temp_dofno;
  
  
   /** END ---> [ Master ------------ Verify] **/ 
  
    int *N_DOFtobeverified;	  /** Array containing how many DOF's need to be verified from other processors(ranks) */
    int *N_DOFverifiedbythisrank; /** Array containing how many DOF's need to be verified by this rank for other processors(ranks) */
    int VerifiedDOF = 0;	  /** Total no of DOF's that need to be verified from other ranks */
    int Verifiedbythisrank =0;    /** Total no of DOF verified by this rank */
    int **sendbuf,**recvbuf,*verrecvbuf,*versendbuf;
     
    /** PROTECTED VARIABLES  * int *sdispl,*rdispl;*/
    int *temp_arr;
	
    sendbuf 	= new int*[2];
    recvbuf 	= new int*[2];
    sdispl  	= new int[size];
    rdispl  	= new int[size];
    temp_arr	= new int[size];
	
    N_DOFtobeverified = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
    N_DOFverifiedbythisrank = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
	
    for(i=0;i<N_U;i++)
    {
      if(Verify[i] == 1)
      {
	N_DOFtobeverified[Master[i]]++;
	VerifiedDOF++;
      }
    }
		
    MPI_Alltoall(N_DOFtobeverified,1, MPI_INT,N_DOFverifiedbythisrank ,1, MPI_INT, Comm);
					
    for(i=0;i<size;i++)
      Verifiedbythisrank +=  N_DOFverifiedbythisrank[i];
		
    
    for(i=0;i<2;i++)
    {
      sendbuf[i] = new int[VerifiedDOF];
      recvbuf[i] = new int[Verifiedbythisrank];
    }
		
    //printf("Rank F %d----------%d\n",rank,Verifiedbythisrank);
    sdispl[0] = 0;
    for(i=1;i<size;i++)
    {
      sdispl[i] = N_DOFtobeverified[i-1] + sdispl[i-1];
//       if(rank == 0)
      //printf(" %d \n",sdispl[i]);
    }
					
    rdispl[0] = 0;
    for(i=1;i<size;i++)
      rdispl[i] = N_DOFverifiedbythisrank[i-1] + rdispl[i-1];
					
					
    memcpy (temp_arr,sdispl, size*SizeOfInt );
					
    for(i=0;i<N_U;i++)
     {
      if(Verify[i] == 1)
      {
	//printf("%d ",temp_arr[Master[i]]);
	sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	sendbuf[DOF_NO][temp_arr[Master[i]]] = MappingData[DOF_NO][i];
	temp_arr[Master[i]]++;
      }
    }

    MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
    MPI_Alltoallv(sendbuf[DOF_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[DOF_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
	      
    verrecvbuf  	= new int[Verifiedbythisrank];
    versendbuf  	= new int[VerifiedDOF];

    for(i=0;i<Verifiedbythisrank;i++)
    {
      temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
      verrecvbuf[i] = Master[GlobalNumbers[temp * N_LocDof + recvbuf[DOF_NO][i]]];
    }
	
    MPI_Alltoallv(verrecvbuf, N_DOFverifiedbythisrank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified, sdispl, MPI_INT, Comm);
					
    for(i=0;i<VerifiedDOF;i++)
    {
      temp_globalno = sendbuf[GLOBAL_NO][i];
      temp_dofno    = sendbuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      if(Verify[temp] != 1)
	printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
      Master[temp] = versendbuf[i];
    }
					
/**===================================== END : MASTER FOR ALL DEGREES OF FREEDOM HAS BEEN DECIDED ========================*/

    delete [] verrecvbuf;
    delete [] versendbuf;
    delete [] N_DOFtobeverified;
    delete [] N_DOFverifiedbythisrank;
    delete [] temp_arr;

		
    for(i=0;i<2;i++)
    {
      delete [] sendbuf[i];
      delete [] recvbuf[i];
    }

    delete [] sendbuf;
    delete [] recvbuf;
		
				
/**==================================START : NOW FIRST GATHERING INFORMATION REGARDING WHICH 
 *                                    DOF's have to be received from OTHER PROCESSORS==================*/
					
    int **SlaveBuf,**MasterBuf,SizeDofSend = 0;
    /** PROTECTED VARIABLES * int *N_DofSend,*N_DofRecv,*DofSend,*DofRecv;*/
    N_Slave = 0;
    N_OwnDof = 0;
					
    N_DofSend = new int[size];
    N_DofRecv = new int[size];
    temp_arr  = new int[size];
    SlaveBuf  = new int*[2];
    MasterBuf = new int*[2];
	
	//for(i=N_Active;i<N_U;i++) Master[i]=rank;
	
    memset(N_DofRecv,0,size*SizeOfInt);
			
    for(i=0;i<N_U;i++)
    {
     if(Master[i] != rank)
     {  
      N_Slave++;
      N_DofRecv[Master[i]]++;
     } 
     else 
       N_OwnDof++;
    }
  //  printf("Rank %d      OwnDofs is %d---------\n",rank,N_OwnDof);

    OwnDofs = new int[N_OwnDof];
    rdispl[0] = 0; 

    MPI_Alltoall(N_DofRecv,1, MPI_INT,N_DofSend ,1, MPI_INT, Comm);

    for(i=1;i<size;i++)
      rdispl[i] = rdispl[i-1] + N_DofRecv[i-1];
	
    sdispl[0] = 0;
    for(i=0;i<size;i++)
    {
      SizeDofSend += N_DofSend[i];
      sdispl[i] = sdispl[i-1] + N_DofSend[i-1];
    }


    DofSend  = new int[SizeDofSend];
    DofRecv  = new int[N_Slave];
    N_SendDof = SizeDofSend;
	
    for(i=0;i<2;i++)
    {
      SlaveBuf[i]  = new int[N_Slave];
      MasterBuf[i] = new int[SizeDofSend];
    }
	
    memcpy (temp_arr,rdispl, size*SizeOfInt );
	
    m=0;
    for(i=0;i<N_U;i++)
    {
      if(Master[i] != rank)
      {
	SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
	DofRecv[temp_arr[Master[i]]] = i;
	temp_arr[Master[i]]++;
      }
      else 
      {
	OwnDofs[m]=i;
	m++;
      }
    }
			
    MPI_Alltoallv(SlaveBuf[GLOBAL_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[GLOBAL_NO], N_DofSend, sdispl, MPI_INT, Comm);
    MPI_Alltoallv(SlaveBuf[DOF_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[DOF_NO], N_DofSend, sdispl, MPI_INT, Comm);
		
    for(i=0;i<SizeDofSend;i++)
    {
      temp_globalno = MasterBuf[GLOBAL_NO][i];
      temp_dofno    = MasterBuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      DofSend[i] = temp;  
    }
  
  //  if(TDatabase::ParamDB->SC_VERBOSE>1)   
	//printf(" Rank %d ------ NUMBER OF DOF's to be sent = %d -------- NUMBER OF DOF's to be recv = %d\n",rank,N_SendDof,N_Slave); 
    for(i=0;i<2;i++)
    {
      delete [] SlaveBuf[i];
      delete [] MasterBuf[i];
    }

    delete [] SlaveBuf;
    delete [] MasterBuf;
    
   if(N_Dim>1)
   {
     for(i=0;i<size;i++)
     {
        sdispl[i]     *= N_Dim;
	rdispl[i]     *= N_Dim;
	
	N_DofSend[i]  *= N_Dim;
        N_DofRecv[i]  *= N_Dim; 
     }
   }
	
  if(N_SendDof>0)
   Send_Info = new double[N_SendDof*N_Dim];
  if(N_Slave>0)
   Recv_Info = new double[N_Slave*N_Dim];
    
    if(TDatabase::ParamDB->SC_VERBOSE>2)
      if(rank==TDatabase::ParamDB->Par_P0)
	printf("ConstructDofMap done !!!\n");
      
       //DEBUG
 int aa;
 if(size<9)
 for(aa=0;aa<size;aa++){
     if(rank==aa){
          printf("\nRank::%d\n",rank);
	  printf("N_OwnDof     = %d\t Not_Own = %d\n",N_OwnDof,N_Slave);
    }
     MPI_Barrier(MPI_COMM_WORLD);
  }//verified
  
//    printf("N_SendDofMS = %d\n",N_SendDofMS);
}


void TParFEMapper3D::Assign_GlobalDofNo()
{ 
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  int *Send_data = new int[N_SendDof];
  int *Recv_data = new int[N_Slave];
  
  int *Send_dataMS, *Send_dataH1, *Send_dataH2;
  int *Recv_dataMS, *Recv_dataH1, *Recv_dataH2;
  
  int i,j,k,start;
  double t1 = MPI_Wtime();
  int *all_own_dofs_info = new int[size];

  Local2Global = new int[N_Dof];
  for(i=0;i<N_Dof;i++)
    Local2Global[i] = -1;
  
  MPI_Allgather(&N_OwnDof, 1, MPI_INT, all_own_dofs_info, 1, MPI_INT, Comm);
  
  k     = 0;
  start = 0;
  while(k<rank){
    start += all_own_dofs_info[k];
    k++;
  }
  
  for(i=0;i<N_Dof;i++)
  {
    if(Master[i] != rank) continue;
    
    Local2Global[i] = start++;
  }
  
//   cout<<"..................rank:: "<<rank<<"  start::"<<start<<"    N_OwnDof::"<<N_OwnDof<<endl;
     if(N_Dim>1)
  {
     for(i=0;i<size;i++)
     {
        sdispl[i]   /= N_Dim;
	sdisplMS[i] /= N_Dim;
	sdisplH1[i] /= N_Dim;
	sdisplH2[i] /= N_Dim;
	
	rdispl[i]     /= N_Dim;
	rdisplMS[i]   /= N_Dim;
	rdisplH1[i]   /= N_Dim;
	rdisplH2[i]   /= N_Dim;
	
	N_DofSend[i]   /= N_Dim;
        N_DofSendMS[i] /= N_Dim;
        N_DofSendH1[i] /= N_Dim;
        N_DofSendH2[i] /= N_Dim;
     
        N_DofRecv[i]    /= N_Dim;
        N_DofRecvMS[i]  /= N_Dim;
        N_DofRecvH1[i]  /= N_Dim;
        N_DofRecvH2[i]  /= N_Dim;
     }
  }
   
  if(TDatabase::ParamDB->MapperType != 2)
  {
    Send_dataMS = Send_data;
    Send_dataH1 = Send_data + N_SendDofMS;
    Send_dataH2 = Send_data + N_SendDofMS +N_SendDofH1;
  
    Recv_dataMS = Recv_data;
    Recv_dataH1 = Recv_data + N_InterfaceS;
    Recv_dataH2 = Recv_data + N_InterfaceS +N_Halo1;
  
//---------------------------------------------------------------------  
    for(i=0;i<N_SendDofMS;i++)
    {
      Send_dataMS[i] = Local2Global[DofSendMS[i]];
    }

    MPI_Alltoallv(Send_dataMS,N_DofSendMS,sdisplMS,MPI_INT,Recv_dataMS,N_DofRecvMS,rdisplMS,MPI_INT,Comm);
  
    for(i=0;i<N_InterfaceS;i++)  
    {
      Local2Global[DofRecvMS[i]] = Recv_dataMS[i];
    }
//----------------------------------------------------------------------    
    for(i=0;i<N_SendDofH1;i++)
    {
      Send_dataH1[i] = Local2Global[DofSendH1[i]];
    }

    MPI_Alltoallv(Send_dataH1,N_DofSendH1,sdisplH1,MPI_INT,Recv_dataH1,N_DofRecvH1,rdisplH1,MPI_INT,Comm);
  
    for(i=0;i<N_Halo1;i++)  
    {
      Local2Global[DofRecvH1[i]] = Recv_dataH1[i];
    }
//------------------------------------------------------------------------
    for(i=0;i<N_SendDofH2;i++)
    {
      Send_dataH2[i] = Local2Global[DofSendH2[i]];
    }
  
    MPI_Alltoallv(Send_dataH2,N_DofSendH2,sdisplH2,MPI_INT,Recv_dataH2,N_DofRecvH2,rdisplH2,MPI_INT,Comm);
  
    for(i=0;i<N_Halo2;i++)  
    {
      Local2Global[DofRecvH2[i]] = Recv_dataH2[i];
    }
//-------------------------------------------------------------------------    
  }
  else
  {
    int *Send_data = new int[N_SendDof];
    int *Recv_data = new int[N_Slave];
    
    for(i=0;i<N_SendDof;i++)
    {
	Send_data[i] = Local2Global[DofSend[i]];
    }
  
    MPI_Alltoallv(Send_data,N_DofSend,sdispl,MPI_INT,Recv_data,N_DofRecv,rdispl,MPI_INT,Comm);
  
    for(i=0;i<N_Slave;i++)  
    {
        Local2Global[DofRecv[i]] = Recv_data[i];
    }
      
  }
  
  if(N_Dim>1)
  {
     for(i=0;i<size;i++)
     {
        sdispl[i]   *= N_Dim;
	sdisplMS[i] *= N_Dim;
	sdisplH1[i] *= N_Dim;
	sdisplH2[i] *= N_Dim;
	
	rdispl[i]     *= N_Dim;
	rdisplMS[i]   *= N_Dim;
	rdisplH1[i]   *= N_Dim;
	rdisplH2[i]   *= N_Dim;
	
	N_DofSend[i]   *= N_Dim;
        N_DofSendMS[i] *= N_Dim;
        N_DofSendH1[i] *= N_Dim;
        N_DofSendH2[i] *= N_Dim;
     
        N_DofRecv[i]    *= N_Dim;
        N_DofRecvMS[i]  *= N_Dim;
        N_DofRecvH1[i]  *= N_Dim;
        N_DofRecvH2[i]  *= N_Dim;
     }
  }
  
  
  delete [] Send_data;
  delete [] Recv_data;
  
  //verification
  if(1)
  for(i=0;i<N_Dof;i++)
  {
    if(Local2Global[i] == -1)
    {
      printf("..........all global dof number not assigned..............\n");
      MPI_Finalize();
      exit(0);
    }
  }
  
  if(rank == 0)
    printf("Time taken for Assign_GlobalDofNo :: %lf\n",MPI_Wtime()-t1);
  
}


#ifdef _HYBRID
void TParFEMapper3D::Color(int &numColors, int *&ptrColors, char type)
{
  int i,j,k;
  
  int *myReorder;
  int ndof;   

  //find the type of dof that needs to be colored
  if(type == 'm'){
    myReorder = Reorder_M;	ndof = N_InterfaceM; }
  else if(type == 'i'){
    myReorder = Reorder_I;	ndof = N_Int;        }
  else if(type == 'D'){
    myReorder = Reorder_D1;	ndof = N_Dept1;      }
  else if(type == 'd'){
    myReorder = Reorder_D2;	ndof = N_Dept2;      }
  //else if(type == 'x'){
    //myReorder = Reorder_D3;	ndof = N_Dept3;      }
  else{
    printf("wrong type!!!\n"); exit(0);              }
  
  //this stores the total number of colrs used for coloring
  numColors = 0;
  if(!ndof)	return;
  
  int max, temp, t;
  //this stores the color number for dofs
  int *allocatedColor = new int[ndof];
  //initialize all colors with default
  for(i=0;i<ndof;i++)
    allocatedColor[i] = -1;
  
  //now we start coloring the dofs
  for(i=0;i<ndof;i++)
  {
    temp = -1;
    k = myReorder[i];
    for(j=RowPtr[k];j<RowPtr[k+1];j++)
    {
      if(KCol[j] >= k || DofMarker[KCol[j]] != type)	continue;
      //if(DofMarker[KCol[j]] != type)	continue;
      
      t = NewGN[KCol[j]];		    //locate the pos of the dof in myreorder to identify its color in allocatedColor
      if(temp < allocatedColor[t])
	temp = allocatedColor[t];
    }
    
    allocatedColor[i] = temp+1;
    
    if(numColors < allocatedColor[i])
      numColors = allocatedColor[i];
  }
  //colors were numbered from 0, hece total will be 1 more
  numColors++;
  
  ptrColors = new int[numColors+1];
  temp = 0; k = 0;
  //arrange the dofs, such that same color dofs are together
  for(i=0;i<numColors;i++)
  {
    ptrColors[i] = k;
    for(j=0;j<ndof;j++)
    {
      if(allocatedColor[j] == i)
      {
	temp = myReorder[j];
	myReorder[j] = myReorder[k];
	myReorder[k] = temp;
	k++;
	
	allocatedColor[j] = -1;
      }
    }
  }
  
  ptrColors[numColors] = k;
  printf("numcolors (type = %c):: %d\t total dof = %d\n",type,numColors,ndof);
//   exit(0);
}
#endif

#endif
