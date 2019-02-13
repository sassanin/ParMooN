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
// @(#)ParFECommunicator3D.h
//
// Class:      TParFECommunicator3D
// Purpose:    Class containing communication routines
//
// Author:     Sashikumaar Ganesan (24.04.15)
//
// History:    Start of implementation 24.04.15 (Sashikumaar Ganesan)
//
// ======================================================================= 

#ifdef _MPI

#ifndef __PARFECOMMUNICATOR3D__
#define __PARFECOMMUNICATOR3D__

#include "mpi.h"

#include <ParFEMapper3D.h>
#include <FESpace3D.h>
#include <SquareStructure.h>
#include <SquareStructure3D.h>

class TParFECommunicator3D
{
  protected:
    TParFEMapper3D *Mapper;
    
    int N_Dim;
    
    int N_Dof;
    
    int N_SendDof, N_SendDofMS, N_SendDofH1, N_SendDofH2;
    
    double *Send_Info, *Send_InfoMS, *Send_InfoH1, *Send_InfoH2;
    
    double *Recv_Info, *Recv_InfoMS, *Recv_InfoH1, *Recv_InfoH2;
    
    int *N_DofSend, *N_DofSendMS, *N_DofSendH1, *N_DofSendH2; 
    
    int *N_DofRecv, *N_DofRecvMS, *N_DofRecvH1, *N_DofRecvH2;
    
    int *sdispl, *sdisplMS, *sdisplH1, *sdisplH2;
    
    int *rdispl, *rdisplMS, *rdisplH1, *rdisplH2; 
    
    int *DofSend, *DofSendMS, *DofSendH1, *DofSendH2;
    
    int *DofRecv, *DofRecvMS, *DofRecvH1, *DofRecvH2;
    
    int N_Slave, N_InterfaceS, N_Halo1, N_Halo2;
    
    MPI_Comm Comm;
    
    int *Master;
    
    
  public:
    TParFECommunicator3D(TParFEMapper3D *mapper);
    
    TParFECommunicator3D();
    
    ~TParFECommunicator3D();
    
    void CommUpdateMS(double *sol);
    
    void CommUpdateH1(double *sol);
    
    void CommUpdateH2(double *sol);
    
    void CommUpdate_M_H1(double *sol);
    
    void CommUpdate(double *sol);
    
    void CommUpdateReduce(double *rhs);
    
    void CommUpdate(double *sol, double *rhs);
    
    void GatherToRoot(double *&GlobalArray, int &GolabalSize, double *LocalArray, int LocalSize, int root);
    
    void ScatterFromRoot(double *GlobalArray, int &GlobalSize, double *LocalArray, int LocalSize, int root);
    
    void CommAverageUpdate(double *sol);
    void CommUpdateMS_Average(double *buffer, int *count);
    void CommUpdateH1_Average(double *buffer, int *count);
    void CommUpdateH2_Average(double *buffer, int *count);
    
    
    int *GetMaster()
    {return Master;}
    
    int GetNDof()
    {return N_Dof;}
    
    int* Get_Local2Global()
    { return Mapper->Get_Local2Global();}
    
    int GetN_Master()
    { return Mapper->GetN_Master();}
    
    char* Get_DofMarker()
    {return Mapper->Get_DofMarker();}
    
};

#endif
#endif
