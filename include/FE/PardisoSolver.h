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
// Class:      TPardisoSolver
// Purpose:    Solve equation system by OpenMP based direct
//	       solver PARDISO	       
//
// Author:     
//
// History:    
//
// =======================================================================

#ifndef __PARDIRECTSOLVER__
#define __PARDIRECTSOLVER__

#include <Database.h>
#include <MooNMD_Io.h>
#include <DirectSolver.h>

#ifdef __2D__
  #include <SquareMatrix2D.h>
  #include <Matrix2D.h>
#else
  #include <SquareMatrix3D.h>
  #include <Matrix3D.h>
#endif

#ifdef _OPENMP
  #include <omp.h>
#endif

class TPardisoSolver : public TDirectSolver
{
  protected:
    /** data address pointer */
    void **pt;

    /** controll parameters */
    int *iparam;
    
    /** controll parameters */
    double *dparam;

    /** Column indicies */
    int *KCol;
    int KCol_size;

    /** Row pointer */
    int *RowPtr;
    int RowPtr_size;

    /** Entries */
    double *Entries;
    int Entries_size;

    /** Number of equations */
    int N_Eq;
  
    /** Number of non-zeros */
    int N_Entries;
    

    /** number of processors */
    int num_prc;

    /** matrix type */
    int mtype;

    /** pressureprojection */
    bool pp;
    int rhs_index;
    
    /** Benchmark */
    std::ofstream *dat;
    double time_a;
	double time_f;
	double time_s;
    int runs_a;
	int runs_f;
	int runs_s;

  protected:
    void ShiftIndicies();
    void PrintIparam();
//     void PrintDparam();
    void ErrorMsg(int);
    void Sort();
    void FreeMemory();
    void AllocMemory();

    int SetMsgLvl()
    { return TDatabase::ParamDB->SC_VERBOSE > 0 ? 1 : 0;};

#ifdef __2D__
    // nstype 1
    void FillRowPtr(int N_U, int N_P, int N_Active,
		    int *RowPtrA, int *RowPtrB, int *ColPtrB);

    void GetTransposedArrays(TMatrix *B, int *&KRowB,
			     int *&ColPtrB, int *&MapB);
#endif

  public:
    TPardisoSolver();

    ~TPardisoSolver();

#ifdef __2D__
    /** scalar */
    void SetMatrix(TSquareMatrix2D *Matrix);
    
    /** NSE_TYPE 1 */
    void SetMatrix(TSquareMatrix2D *A, TMatrix *B1, TMatrix *B2);

    void SetMatrixPar(TSquareMatrix2D *A, TMatrix *B1, TMatrix *B2);

   /** NSE_TYPE 2 */
   void SetMatrix(TSquareMatrix2D *sqmatrixA,
		  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
		  TMatrix2D *matrixB1,  TMatrix2D *matrixB2);

   /** NSE_TYPE 4 */
   void SetMatrix(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
		  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
		  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
		  TMatrix2D *matrixB1,  TMatrix2D *matrixB2);
#endif

#ifdef __3D__
    /** NSTYPE 2 */
    void SetMatrix(TSquareMatrix3D *sqmatrixA,
		   TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		   TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3);
		
    /** NSTYPE 4 */
    void SetMatrix(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
		      TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
		      TSquareMatrix3D *sqmatrixA22, TSquareMatrix3D *sqmatrixA23,
		      TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
		      TSquareMatrix3D *sqmatrixA33,
		      TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		      TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3);
#endif

    void Analyse();
    void Factorize();
    void Solve(double *sol, double *rhs);
    void FactorizeSolve(double *sol, double *rhs);

    /** Benchmark */

    void BenchReset();

};

#endif
