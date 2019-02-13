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
// @(#)ParFECommunicator2D.h
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

#include <FESpace2D.h>

#ifndef __PARFECOMMUNICATOR2D__
#define __PARFECOMMUNICATOR2D__

class TParFECommunicator2D 
{
  protected:

    /** MPI_Comm for which the fespace communications are needed */
     MPI_Comm Comm;

    /** fespace for which the communications are needed */
    TFESpace2D *fespace;

    /**number of own dof (excluding Halo dof) */
    int N_OwnDof;

    /**number of Active own dof (excluding Dirichlet and Halo dof) */
    int N_ActiveOwnDof;

    /**number of dependent Cells*/
    int N_DependentCells;

    /**dependent Cells indices*/
    int *DependentCellIndex;

    /**dependent Cells indices*/
    int *N_DependentCellNeibs;

    /**dependent Cells indices*/
    int *DependentCellNeibIDs;

    /**maximum number(possible) of subdomains among all dofs */
    int MaxSubDomainPerDof;

    /**array containing list of all rank(process ID) associated with each dof (including Halo dof) */
    int *N_DofRankIndex;

    /**array rank(process) ID of each dof (including Halo dof) */
    int *DofRankIndex;

    /**array containing the corresponding global dof of local dof */
    int *GlobalDofOFLocalDof;

    /**own dof of size  N_OwnDof, excluding the Halo dofs  */
    int *OwnDofIndex;

    /**Active own dof of size  N_OwnDof, excluding Dirichlet and Halo dofs  */
    int *ActiveOwnDofIndex;

    /**list of start and end of OwnDofNeibPtrList of each dof  */
    int *OwnDofNeibPtr;

    /**list of neib rank IDs of each own dof, which value has to be send to the neibs */
    int *OwnDofNeibPtrList;

    /**is this dof value need for any neib ranks*/
    bool *NeedAtNeib;

    /**total number of DOF neib ranks for this rank */
    int N_DofNeibs;

    /**DOF neib ranks' IDs of this rank */
    int *DofNeibIDs;

    /**index in the neib list (DofNeibIDs) for the given neib rank */
    int *IndexOfNeibRank;

    /** Max number of Communication Steps for all ranks*/
    int Max_CommunicationSteps;

    /** Number of Communication Processes in each step of  communication*/
    int *N_CommunicationProcesses;

    /** Receive rank ID in each Communication Process*/
    int *ReceiveID_Index;

    /** send rank ID in each Communication Process*/
    int *SendID_Index;

    /**request tags for non-blocking communication */
    MPI_Request request001, request002, request003, request004, request005, request006;
    
    MPI_Request requestGatherVectorAtRoot;


    /** ************************************************************ */
    /** values which are available only at root */
    /** ************************************************************ */

    /**number of local dof (including Halo) among all ranks, except root */
    int *N_LocalDofAllRank;

    /**maximum number of local dof (including Halo) among all ranks, except root */
    int MaxN_LocalDofAllRank;

    /**array containing the correponding global dof of local dof for all ranks*/
    int *GlobalDofOFLocalDofAllRank;

  private:
    /** collect all communication info for the given fespace */
    /** assumed that the root contains the Global fespace */
    int MakeDofMappingFromRoot();
    int MakeDofMappingFromNeib();
    int ConstructDofRankIndex();

  public:
    /** constructor with a fespace */
    TParFECommunicator2D(MPI_Comm comm, TFESpace2D *FESpace);

    int WaitForMakeDofMapping();

    void SetOwnDof();

    /** set the scheduling algorightm, which rank communicates to which rank and when and how */
    void SetFENeibCommunicationSteps();

    /** send and recieve N_array of values between the neibs */
    void MooNMD_FECommunicateNeib(int **sendbuf, int **recvbuf, int N_array);

    /** send and recieve N_array of arrays with given displacements between the neibs */
    void MooNMD_FECommunicateNeib(int **sendbuf, int *senddisp, int **sendlen, int **recvbuf, int *recvdisp, int **recvlen, int N_array);

    MPI_Comm GetComm()
     { return Comm;}

    int GetN_OwnDof();       

    int *GetOwnDofIndex();

    int GetActiveN_OwnDof();       

    int *GetActiveOwnDofIndex();

    int GetMaxSubDomainPerDof()
        { return MaxSubDomainPerDof; }

    int *GetDofRankIndex()
        { return DofRankIndex; }

    int *GetGlobalDofOFLocalDof()
        { return GlobalDofOFLocalDof; }

    TFESpace2D *Getfespace()
      {return fespace; }

    int GetMaxN_LocalDofAllRank()
      {return MaxN_LocalDofAllRank; }


    int *GetN_LocalDofAllRank()
      {return N_LocalDofAllRank; }


    int *GetGlobalDofOFLocalDofAllRank()
      {return GlobalDofOFLocalDofAllRank; }

    /** destructor */
    ~TParFECommunicator2D();

};
#endif
#endif
