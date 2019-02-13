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
#ifdef _SMPI
#include "mpi.h"

#ifdef __2D__
#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#else  
#include <SquareStructure2D.h>
#include <SquareMatrix2D.h>
#endif

#include <MumpsSolver.h>

#ifndef __SEQDIRECTSOLVER_MUMPS__
#define __SEQDIRECTSOLVER_MUMPS__

class TSeqParDirectSolver
{
protected:
  
  MPI_Comm Comm;
  
  TMumpsSolver *Mumps;
 
  //problem type: NSE 0, SCALAR 100
  int SystemType;
  
  //Select ParDirectSolver type
  int DSType;
  
  int NSEType;
  
  //Mumps Solver Parameters
  int N_Master, NDof;
  int N_Master_U, N_Master_P, NDof_U, NDof_P, NDof_S, NDof_D;
  
  //row ptr and col ptr of system matrix
  int *RowPtr, *KCol, *RowPtr_global;
  int *RowPtr_B, *KCol_B, *RowPtr_BT, *KCol_BT;
 
  
  // number of non zeroes in system matrix
  int N_Nz, N_Nz_global;
  int N_Nz_U, N_Nz_P, N_Nz_G, N_Nz_C, N_Nz_D, N_Nz_H, N_Nz_E, N_Nz_J;
  
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
#ifdef __3D__  
  //SqMatrix
  TSquareMatrix3D *Mat;
  TMatrix3D *MatB, *MatBT;
  
  
  //for NSTYPE 4 - 3D
  TSquareMatrix3D *MatA11,*MatA12,*MatA13,
                  *MatA21,*MatA22,*MatA23,
		  *MatA31,*MatA32,*MatA33;
		  
  TMatrix3D       *MatB1, *MatB2, *MatB3,
                  *MatBT1,*MatBT2,*MatBT3;
#endif
		  
#ifdef __2D__
  //for NSECST - 2D 
		  
int *RowPtr_G, *KCol_G, *RowPtr_H, *KCol_H, *RowPtr_C, *KCol_C, *RowPtr_D, *KCol_D;
int *RowPtr_E, *KCol_E, *RowPtr_J, *KCol_J;

  TSquareMatrix2D *MatA, *MatG, *MatH;
  TMatrix2D *MatB, *MatBT, *MatC, *MatE, *MatD, *MatJ;  
  
    TSquareMatrix2D *SqmatrixA11, *SqmatrixA12, 
                    *SqmatrixA21, *SqmatrixA22,
                    *SqmatrixG11, *SqmatrixG12, 
		    *SqmatrixG21, *SqmatrixG22, *SqmatrixG23, 
		    *SqmatrixG32, *SqmatrixG33,
                    *SqmatrixH11, *SqmatrixH22, *SqmatrixH33;  
		    
    TMatrix2D *MatrixB1, *MatrixB2, *MatrixB1T, *MatrixB2T,
              *MatrixC11, *MatrixC12, *MatrixC22, *MatrixC23,
              *MatrixE11, *MatrixE12, *MatrixE22, *MatrixE23,
              *MatrixD11, *MatrixD12, *MatrixD21, *MatrixD22, *MatrixD31, *MatrixD32,
              *MatrixJ11, *MatrixJ21, *MatrixJ22, *MatrixJ32;  
#endif
  //MatLoc is the system matrix entry values(only master rows)
  double *MatLoc;
  
  double *MatGlobal;
  
  int *all_Nnz;
  
  int offset;
  
public:
 
  #ifdef __3D__  
  TSeqParDirectSolver(int dim,int N_U,int N_P,TSquareMatrix3D **mat,TMatrix3D **matB);
  
  void AssembleMatrix();
  
  void AssembleMatrix_NSE2();
  
  void AssembleMatrix_NSE4();
  
  void InitMumps_Scalar();
  
  void InitMumps_NSE2();  
  
  void InitMumps_NSE4();
  
  void Solve(double *Sol, double *Rhs, bool Factorize);
  #endif
  
  
  #ifdef __2D__
  // for NSE-CST-2D
  TSeqParDirectSolver(int N_U,int N_P, int N_S, int N_D, TSquareMatrix2D **sqmat,TMatrix2D **recmat);
  
  void InitMumps_NSEType4();
  
   void InitMumps_NSEType4_CST_2D();
  
  void InitMumps_NSEType4_CST_DEVSS_2D();
  
   void AssembleMatrix();
   
   void AssembleMatrix_NSEType4D();
    
    void AssembleMatrix_NSEType4_CST_2D();
      
    void AssembleMatrix_NSEType4_CST_DEVSS2D();
  
  void Solve(double *Sol, double *Rhs, bool Factorize);
  
  #endif
  
  ~TSeqParDirectSolver();
  
  void GetRhs(double *Rhs);
  
  void UpdateSol(double *Sol);
  


};
 
#endif
#endif
