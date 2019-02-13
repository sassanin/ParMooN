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
// @(#)ParFEMapper3D.h
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

#ifndef __TParFEMapper3D__
#define __TParFEMapper3D__

#include "mpi.h"

#include <FESpace3D.h>
#include <SquareStructure.h>
#include <SquareStructure3D.h>

class TParFEMapper3D
{
  protected:
    /** MPI_Comm for which the fespace communications are needed */
    MPI_Comm Comm;
     
    /** Number of dimensions 2d/3d **/
    int N_Dim;
    
    int N_Dof;
    
    /** row ptr and col ptr of A aatrix**/
    int *RowPtr, *KCol;
    
    /** fespace for which the communications are needed */
    TFESpace3D *FESpace;
    
    /**maximum number(possible) of subdomains among all dofs */
    int MaxSubDomainPerDof;
    
    /** array containing the master of the dofs **/ 
    int *Master, *OwnDofs;
    
    /**array containing the dof type**/
    char *DofMarker;
    
    /**array which carries info across procs**/
    double *Send_Info,*Recv_Info; 
    
    /**array to store globalDofNo of all its dofs**/
    int *Local2Global;
    
    /** info of number of dofs to be send/recv from other procs**/
    int *sdispl, *rdispl;
    
    int N_InterfaceM, N_InterfaceS, N_Halo1, N_Halo2, N_Dept1, N_Dept2;
    
    int N_Slave, N_Halo, N_Master, N_Dept, N_Int, N_OwnDof;
    
    int *N_DofSend, *N_DofRecv, *N_DofSendMS, *N_DofSendH1, *N_DofSendH2, *N_DofRecvMS, *N_DofRecvH1, *N_DofRecvH2;
     
    int *sdisplMS, *sdisplH1, *sdisplH2, *rdisplMS, *rdisplH1, *rdisplH2;
     
    int N_SendDof, N_SendDofMS, N_SendDofH1, N_SendDofH2;
     
    int *DofSend, *DofSendMS, *DofSendH1, *DofSendH2, *DofRecv, *DofRecvMS, *DofRecvH1, *DofRecvH2;
     
    double *Send_InfoMS, *Send_InfoH1, *Send_InfoH2, *Recv_InfoMS, *Recv_InfoH1, *Recv_InfoH2; 
    
    int *NewGN, *Reorder, *Reorder_M, *Reorder_I, *Reorder_D1, *Reorder_D2, *Reorder_D3;
     
    int N_CMaster, N_CDept1, N_CDept2, N_CInt, *ptrCMaster, *ptrCDept1, *ptrCDept2, *ptrCInt;
     
  public:
    TParFEMapper3D(int N_dim, TFESpace3D *fespace, int *rowptr, int *kcol);
    
    ~TParFEMapper3D();
    
    void ConstructDofMap_Master_Halo();		//MapperType 1
    
    void ConstructDofMap();			//MapperType 2
    
    int find_min(int *arr, int N, char *temp_arr);
    
    void GetCommInfo(int &n_Dim, int &n_Dof,
		     int &n_SendDof, int &n_SendDofMS, int &n_SendDofH1, int &n_SendDofH2,
		     double *&send_Info, double *&send_InfoMS, double *&send_InfoH1, double *&send_InfoH2,
		     double *&recv_Info, double *&recv_InfoMS, double *&recv_InfoH1, double *&recv_InfoH2,
		     int *&n_DofSend, int *&n_DofSendMS, int *&n_DofSendH1, int *&n_DofSendH2,
		     int *&n_DofRecv, int *&n_DofRecvMS, int *&n_DofRecvH1, int *&n_DofRecvH2,
		     int *&Sdispl, int *&SdisplMS, int *&SdisplH1, int *&SdisplH2,
		     int *&Rdispl, int *&RdisplMS, int *&RdisplH1, int *&RdisplH2,
		     int *&dofSend, int *&dofSendMS, int *&dofSendH1, int *&dofSendH2,
		     int *&dofRecv, int *&dofRecvMS, int *&dofRecvH1, int *&dofRecvH2,
		     int &n_Slave, int &n_InterfaceS, int &n_Halo1, int &n_Halo2)
    {
      n_Dim        = N_Dim; 
      n_Dof        = N_Dof;
      
      n_SendDof    = N_SendDof;
      n_SendDofMS  = N_SendDofMS;
      n_SendDofH1  = N_SendDofH1;
      n_SendDofH2  = N_SendDofH2;
      
      send_Info    = Send_Info;
      send_InfoMS  = Send_InfoMS;
      send_InfoH1  = Send_InfoH1;
      send_InfoH2  = Send_InfoH2;
      
      recv_Info    = Recv_Info;
      recv_InfoMS  = Recv_InfoMS;
      recv_InfoH1  = Recv_InfoH1;
      recv_InfoH2  = Recv_InfoH2;
      
      n_DofSend    = N_DofSend;
      n_DofSendMS  = N_DofSendMS;
      n_DofSendH1  = N_DofSendH1;
      n_DofSendH2  = N_DofSendH2;
    
      n_DofRecv    = N_DofRecv;
      n_DofRecvMS  = N_DofRecvMS;
      n_DofRecvH1  = N_DofRecvH1;
      n_DofRecvH2  = N_DofRecvH2;
    
      Sdispl       = sdispl;
      SdisplMS     = sdisplMS;
      SdisplH1     = sdisplH1;
      SdisplH2     = sdisplH2;
      
      Rdispl       = rdispl;
      RdisplMS     = rdisplMS;
      RdisplH1     = rdisplH1;
      RdisplH2     = rdisplH2;
    
      dofSend      = DofSend;
      dofSendMS    = DofSendMS;
      dofSendH1    = DofSendH1;
      dofSendH2    = DofSendH2;
      
      dofRecv      = DofRecv;
      dofRecvMS    = DofRecvMS;
      dofRecvH1    = DofRecvH1;
      dofRecvH2    = DofRecvH2;
    
      n_Slave      = N_Slave;
      n_InterfaceS = N_InterfaceS;
      n_Halo1      = N_Halo1;
      n_Halo2      = N_Halo2;
    }
    
    int GetN_Master()
    {return N_Master;}
    
    int *GetMaster()
    {return Master;}
    
    char *Get_DofMarker()
    {  return DofMarker; }
    
    int* GetReorder_M()
    {return Reorder_M;}
    
    int* GetReorder_I()
    {return Reorder_I;}
    
    int* GetReorder_D1()
    {return Reorder_D1;}
    
    int* GetReorder_D2()
    {return Reorder_D2;}
    
    int* GetReorder_D3()
    {return Reorder_D3;}
    
    int GetN_InterfaceM()
    {return N_InterfaceM;}
    
    int GetN_Int_light()
    {return N_Int;}
    
    int GetN_Dept1()
    {return N_Dept1;}
    
    int GetN_Dept2()
    {return N_Dept2;}
    
    int* Get_Local2Global()
    { return Local2Global;}

#ifdef _HYBRID
    void Color(int &numColors, int *&ptrColors, char type);
    
    int GetN_CMaster()
    {return N_CMaster;}
    int* GetptrCMaster()
    {return ptrCMaster;}
    
    int GetN_CDept1()
    {return N_CDept1;}
    int* GetptrCDept1()
    {return ptrCDept1;}
    
    int GetN_CDept2()
    {return N_CDept2;}
    int* GetptrCDept2()
    {return ptrCDept2;}
    
    int GetN_CInt()
    {return N_CInt;}
    int* GetptrCInt()
    {return ptrCInt;}
    
#endif
  
    void Assign_GlobalDofNo();
};

#endif
#endif
