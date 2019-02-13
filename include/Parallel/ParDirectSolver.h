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
//
// =======================================================================
#ifdef _MPI
#include "mpi.h"
#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>

#include <MumpsSolver.h>

#ifndef __DIRECTSOLVER_MUMPS__
#define __DIRECTSOLVER_MUMPS__

class TParDirectSolver
{
protected:
  TParFECommunicator3D *ParComm;
  
  TParFECommunicator3D *ParComm_P;
  
  MPI_Comm Comm;
  
  TMumpsSolver *Mumps;
 
  //problem type: NSE 0, SCALAR 100
  int SystemType;
  
  //Select ParDirectSolver type
  int DSType;
  
  int NSEType;
  
  //Mumps Solver Parameters
  int N_Master, NDof;
  int N_Master_U, N_Master_P, NDof_U, NDof_P;
  
  //row ptr and col ptr of system matrix
  int *RowPtr, *KCol, *RowPtr_global;
  int *RowPtr_B, *KCol_B, *RowPtr_BT, *KCol_BT;
  
  // number of non zeroes in system matrix
  int N_Nz, N_Nz_global;
  int N_Nz_U, N_Nz_P;
  
  // number of global degrees of freedom over all ranks
  int N_Eqns;
  
  // number of rhs, OwnRhs(only master dofs)
  int N_rhs;
  double *OwnRhs;
  int Global_N_DOF_U;
  double *GlobalRhs,*GlobalSol;
  int GlobalRhsSize;
  
  //(I,J) for MUMPS input
  int *I_rn, *J_cn;
  
  //SqMatrix
  TSquareMatrix3D *Mat;
  TMatrix3D *MatB, *MatBT;
  
  //for NSTYPE 4
  TSquareMatrix3D *MatA11,*MatA12,*MatA13,
                  *MatA21,*MatA22,*MatA23,
		  *MatA31,*MatA32,*MatA33;
		  
  TMatrix3D       *MatB1, *MatB2, *MatB3,
                  *MatBT1,*MatBT2,*MatBT3;
  
  //MatLoc is the system matrix entry values(only master rows)
  double *MatLoc;
  
  double *MatGlobal;
  
  int *all_Nnz;
  
  int offset;
  
public:
 
  TParDirectSolver(TParFECommunicator3D *parcomm,TParFECommunicator3D *parcomm_p,TSquareMatrix3D **mat,TMatrix3D **matB);
  
  ~TParDirectSolver();
  
  void AssembleMatrix();
  
  void AssembleMatrix_NSE2();
  
  void AssembleMatrix_NSE4();
  
  void InitMumps_Scalar();
  
  void InitMumps_NSE2();  
  
  void InitMumps_NSE4();  
  
  void InitPardiso();
  
  void GetRhs(double *Rhs);
  
  void UpdateSol(double *Sol);
  
  void Solve(double *Sol, double *Rhs, bool Factorize);
  
  void InitPardiso_test();
};
 
#endif
#else

//------------------------------------------------------------------------------------------------------//
#ifdef _OMPONLY
#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <ParDiso.h>

#ifndef __DIRECTSOLVER_MUMPS__
#define __DIRECTSOLVER_MUMPS__

class TParDirectSolver
{
protected:
  TParDiso *ParDiso;
  
  //Select ParDirectSolver type
  int DSType;
  
  //Pardiso Solver Parameters
  int NDof;
  
  //row ptr and col ptr of system matrix
  int *RowPtr, *KCol;
  
  // number of non zeroes in system matrix
  int N_Nz;
  
  // number of global degrees of freedom over all ranks
  int N_Eqns;
  
  // number of rhs, OwnRhs(only master dofs)
  int N_rhs;
  
  //SqMatrix
  TSquareMatrix3D *Mat;
  
  //MatLoc is the system matrix entry values(only master rows)
  double *Mat_Values;
  
public:
  TParDirectSolver(TSquareMatrix3D *mat);
  
  ~TParDirectSolver();
  
  void AssembleMatrix(TSquareMatrix3D *matrix);

  void InitPardiso();

  void Solve(double *Sol, double *Rhs, bool Factorize);
  
};
 
#endif
#endif
#endif
