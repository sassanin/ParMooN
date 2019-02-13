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
// @(#)ParDiso.h
//
// Class:      TParDiso
// Purpose:    Solve equation system by ParDiso routines
//
// Author:     Sashikumaar Ganesan & Abdus Shamim (01.05.15)
//
// History:    Start of implementation 01.05.15 (Sashikumaar Ganesan & Abdus Shamim) 
//
// =======================================================================

#ifdef _OMPONLY

#ifndef __PARDISO__
#define __PARDISO__

#include <Database.h>
// #include <SquareMatrix.h>
// #include <Matrix.h>
// #include <SquareMatrix3D.h>
// #include <Matrix3D.h>
// #include <ParFECommunicator3D.h>

class TParDiso
{
  protected:
    
  void *pt[64];
  int iparam[64];
  double dparam[64];
  int phase, nrhs, Nmax_system, matrix_number, matrix_type, N_Eqns, N_Nz;
  //csr format
  int *RowPtr, *KCol;
  int Solver;
  int perm_user, msglvl, ierror;
  
  double idum, ddum;
  
  double *rhsptr,*solptr;
  
  public:
   /** constructor */
   TParDiso(int neqns, int nnz, int* rowptr, int* kcol);
   
   void FactorizeAndSolve(double *Mat, double *rhs, double *sol,bool Factorize);

   void Solve(double *Mat, double *rhs, double *sol);
   
   void Clean(double *Mat);
  
    /** destructor */
    ~TParDiso();
};
#endif


#endif
