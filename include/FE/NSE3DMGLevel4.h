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
// @(#)NSE3DMGLevel4.h        1.6 07/03/00
//
// Class:       TNSE3DMGLevel4
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system of
//              type 4 (all Aij, B1, B2, B3, B1T, B2T, B3T)
//
// Author:      Volker John 25.07.2000
//
// History:     25.07.2000 start of implementation
//
// =======================================================================

#ifndef __NSE3DMGLEVEL4__
#define __NSE3DMGLEVEL4__

#include <NSE3DMGLevel.h>

class TNSE3DMGLevel4 : public TNSE3DMGLevel
{
  protected:
    /** matrix A11 */
    TSquareMatrix3D *A11;

    /** matrix A12 */
    TSquareMatrix3D *A12;

    /** matrix A12 */
    TSquareMatrix3D *A13;

    /** matrix A21 */
    TSquareMatrix3D *A21;

    /** matrix A22 */
    TSquareMatrix3D *A22;

    /** matrix A23 */
    TSquareMatrix3D *A23;

    /** matrix A31 */
    TSquareMatrix3D *A31;

    /** matrix A32 */
    TSquareMatrix3D *A32;

    /** matrix A33 */
    TSquareMatrix3D *A33;

    /** structure of matrix A */
    TSquareStructure3D *StructureA;

    /** row pointer for matrix A */
    int *ARowPtr;

    /** column number vector for matrix A */
    int *AKCol;

    /** matrix entries of matrix A */
    double *A11Entries;

    /** matrix entries of matrix A */
    double *A12Entries;

    /** matrix entries of matrix A */
    double *A13Entries;

    /** matrix entries of matrix A */
    double *A21Entries;

    /** matrix entries of matrix A */
    double *A22Entries;

    /** matrix entries of matrix A */
    double *A23Entries;

    /** matrix entries of matrix A */
    double *A31Entries;

    /** matrix entries of matrix A */
    double *A32Entries;

    /** matrix entries of matrix A */
    double *A33Entries;

    /** matrix B1 */
    TMatrix3D *B1;

    /** matrix B2 */
    TMatrix3D *B2;

    /** matrix B3 */
    TMatrix3D *B3;

    /** structure of matrix B */
    TStructure3D *StructureB;

    /** row pointer for matrix B */
    int *BRowPtr;

    /** column number vector for matrix B */
    int *BKCol;

    /** matrix entries of matrix B1 */
    double *B1Entries;

    /** matrix entries of matrix B2 */
    double *B2Entries;

    /** matrix entries of matrix B2 */
    double *B3Entries;

    /** matrix B1T */
    TMatrix3D *B1T;

    /** matrix B2T */
    TMatrix3D *B2T;

    /** matrix B3T */
    TMatrix3D *B3T;

    /** structure of matrix BT */
    TStructure3D *StructureBT;

    /** row pointer for matrix BT */
    int *BTRowPtr;

    /** column number vector for matrix BT */
    int *BTKCol;

    /** matrix entries of matrix B1T */
    double *B1TEntries;

    /** matrix entries of matrix B2T */
    double *B2TEntries;
 
    /** matrix entries of matrix B3T */
    double *B3TEntries;

  public:
    /** constructor */
    TNSE3DMGLevel4(int level, 
                   TSquareMatrix3D *A11, TSquareMatrix3D *A12, 
                   TSquareMatrix3D *A13, 
                   TSquareMatrix3D *A21, TSquareMatrix3D *A22, 
                   TSquareMatrix3D *A23, 
                   TSquareMatrix3D *A31, TSquareMatrix3D *A32, 
                   TSquareMatrix3D *A33, 
                   TMatrix3D *B1, TMatrix3D *B2, TMatrix3D *B3,
                   TMatrix3D *B1T, TMatrix3D *B2T, TMatrix3D *B3T,
                   double *f1, double *f2, double *f3, double *g,
                   double *u1, double *u2, double *u3, double *p,
                   int n_aux, double al, int VelocitySpace, 
                   int PressureSpace, TCollection *coll);

    /** destructor */
    ~TNSE3DMGLevel4();

    /** calculate defect */
    virtual void Defect(double *u1, double *u2, double *u3, double *p,
                        double *f1, double *f2, double *f3, double *g,
                        double *d1, double *d2, double *d3, double *d4,
                        double &res);

    /** correct Dirichlet and hanging nodes */
    virtual void CorrectNodes(double *u1, double *u2, double *u3, double *p);

    /** Vanka smoother */
    virtual void CellVanka(double *u1, double *u2,  double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *def1, double *def2, double *def3, double *def4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

    /** Vanka smoother */
    virtual void NodalVanka(double *u1, double *u2, double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *def1, double *def2, double *def3, double *def4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

    /** solve exact on this level */
    virtual void SolveExact(double *u1, double *u2, double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters);

    /** Braess Sarazin smoother */
    virtual void BraessSarazin(double *u1, double *u2, double *u3, double *p,
        double *rhs1, double *rhs2, double *rhs3, double *rhs4,
        double *aux, double *Counters,
        int N_Parameters, double *Parameters,int N_Levels);

    /** step length control for Vanka */
    virtual double StepLengthControl(double *u1, double *u2, double *u3, double *p,
           double *u1old, double *u2old, double *u3old, double *pold,
           double *def1, double *def2, double *def3, double *def4,
           int N_Parameters, double *Parameter);

    /** print all matrices and both right hand sides */
    virtual void PrintAll();

};

#endif
