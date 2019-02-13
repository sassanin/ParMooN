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
// @(#)DirectSolver.h
// 
// Purpose:     solve equation system by direct solver
//
// Author:      Gunar Matthies (06.09.05)
//
// History:     start of implementation 06.09.05 (Gunar Matthies)
//
// =======================================================================

#ifndef __DIRECTSOLVER__
#define __DIRECTSOLVER__

#include <SquareMatrix2D.h>
#include <Matrix2D.h>

#ifdef __3D__
#include <SquareMatrix3D.h>
#include <Matrix3D.h>
#endif

class TDirectSolver
{
  public:
    TDirectSolver() {};
    
    virtual ~TDirectSolver() {};
    
#ifdef __3D__
    /** NSTYPE 2 */
    virtual void SetMatrix(TSquareMatrix3D *sqmatrixA,
		   TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		   TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3) = 0;
		   
    /** NSTYPE 4 */
    virtual void SetMatrix(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
		      TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
		      TSquareMatrix3D *sqmatrixA22, TSquareMatrix3D *sqmatrixA23,
		      TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
		      TSquareMatrix3D *sqmatrixA33,
		      TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		      TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3) = 0;
#endif
    
    virtual void Analyse() = 0;
    virtual void Factorize() = 0;
    virtual void Solve(double *sol, double *rhs) = 0;
    virtual void FactorizeSolve(double *sol, double *rhs) = 0;
    
};

/** solve equation system */

void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol);
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag);

void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp);
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag);

void DirectSolverLong(TSquareMatrix *matrix, double *rhs, double *sol);

void DirectSolver(TSquareMatrix2D *sqmatrixA,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  double *rhs1, double *rhs2, double *sol1, double *sol2, int rb_flag=3);

/*void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag=3);*/

void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag);

//****************************************************************************/
// for NSTYPE == 1
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

void DirectSolver(TSquareMatrix2D *sqmatrixA, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag);

//****************************************************************************/
// for NSTYPE == 2
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  TMatrix2D *matrixC,
                  double *rhs, double *sol);

//****************************************************************************/
// for NSTYPE == 4
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

//****************************************************************************/
// for NSTYPE == 14
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
		  TSquareMatrix2D *sqmatrixC,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

//****************************************************************************/
// for Darcy
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA, TSquareMatrix2D *sqmatrixC,
                  TMatrix2D *matrixBT, TMatrix2D *matrixB,
                  double *rhs, double *sol);

//****************************************************************************/
// for Conformation stress tensor
//****************************************************************************/

void DirectSolver(TSquareMatrix2D *sqmatrixS11, TSquareMatrix2D *sqmatrixS12, 
                  TSquareMatrix2D *sqmatrixS21, TSquareMatrix2D *sqmatrixS22, TSquareMatrix2D *sqmatrixS23, 
		   TSquareMatrix2D *sqmatrixS32, TSquareMatrix2D *sqmatrixS33,
                  double *rhs, double *sol);

//****************************************************************************/
// for Deformation tensor in DEVSS
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixS11,  TSquareMatrix2D *sqmatrixS22, 
		   TSquareMatrix2D *sqmatrixS33, double *rhs, double *sol);


//****************************************************************************/
// Monolithic NSE(Type4)_CST_2D with DEVSS
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12, TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
		  TSquareMatrix2D *SqmatrixG11, TSquareMatrix2D *SqmatrixG12, TSquareMatrix2D *SqmatrixG21 , TSquareMatrix2D *SqmatrixG22, TSquareMatrix2D *SqmatrixG23 , TSquareMatrix2D *SqmatrixG32 , TSquareMatrix2D *SqmatrixG33,
                  TSquareMatrix2D *SqmatrixH11, TSquareMatrix2D *SqmatrixH22, TSquareMatrix2D *SqmatrixH33,
		  TMatrix2D *matrixB1,  TMatrix2D *matrixB2, TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *MatrixC11, TMatrix2D *MatrixC12, TMatrix2D *MatrixC22, TMatrix2D *MatrixC23,
		  TMatrix2D *MatrixD11, TMatrix2D *MatrixD12, TMatrix2D *MatrixD21, TMatrix2D *MatrixD22, TMatrix2D *MatrixD31, TMatrix2D *MatrixD32,
	          TMatrix2D *MatrixE11, TMatrix2D *MatrixE12, TMatrix2D *MatrixE22, TMatrix2D *MatrixE23,
	          TMatrix2D *MatrixJ11, TMatrix2D *MatrixJ21, TMatrix2D *MatrixJ22, TMatrix2D *MatrixJ32,
                  double *rhs, double *sol);


//****************************************************************************/
// Monolithic NSE(Type4)_CST_2D with LPS
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12, TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
		  TSquareMatrix2D *SqmatrixG11, TSquareMatrix2D *SqmatrixG12, TSquareMatrix2D *SqmatrixG21 , TSquareMatrix2D *SqmatrixG22, TSquareMatrix2D *SqmatrixG23 , TSquareMatrix2D *SqmatrixG32 , TSquareMatrix2D *SqmatrixG33,                
		  TMatrix2D *matrixB1,  TMatrix2D *matrixB2, TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *MatrixC11, TMatrix2D *MatrixC12, TMatrix2D *MatrixC22, TMatrix2D *MatrixC23,
		  TMatrix2D *MatrixD11, TMatrix2D *MatrixD12, TMatrix2D *MatrixD21, TMatrix2D *MatrixD22, TMatrix2D *MatrixD31, TMatrix2D *MatrixD32,
                  double *rhs, double *sol);


#ifdef __3D__
//****************************************************************************/
// for NSTYPE == 2
//****************************************************************************/
void DirectSolver(TSquareMatrix3D *sqmatrixA,
                  TMatrix3D *matrixB1T, TMatrix3D *matrixB2T,
                  TMatrix3D *matrixB3T,
                  TMatrix3D *matrixB1,  TMatrix3D *matrixB2,
                  TMatrix3D *matrixB3,
                  double *rhs, double *sol);
//****************************************************************************/
// for NSTYPE == 4
//****************************************************************************/
void DirectSolver(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
		  TSquareMatrix3D *sqmatrixA13,
                  TSquareMatrix3D *sqmatrixA21, TSquareMatrix3D *sqmatrixA22,
		  TSquareMatrix3D *sqmatrixA23,
                  TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
		  TSquareMatrix3D *sqmatrixA33,
                  TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T, 
                  TMatrix3D *matrixB1,  TMatrix3D *matrixB2, TMatrix3D *matrixB3,
                  double *rhs, double *sol, int flag);


void DirectSolver(TSquareMatrix3D **sqmatrices, int n_row, int n_column,
		  double *sol, double *rhs);

void DirectSolver(TSquareMatrix3D **sqmatrices, int n_row, int n_column,
                  double *sol, double *rhs, double *&Entries,
                   int *&KCol, int *&RowPtr, void *&Symbolic, void *&Numeric, int rb_flag);

#endif
#endif
