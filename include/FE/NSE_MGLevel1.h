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
// @(#)NSE_MGLevel1.h        1.6 07/03/00
//
// Class:       TNSE_MGLevel1
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system of
//              type 1 (A, B1, B2, B3)
//
// Author:      Volker John 25.07.2000
//
// History:      25.07.2000 start of implementation
//
// =======================================================================

#ifndef __NSE_MGLEVEL1__
#define __NSE_MGLEVEL1__

#include <NSE_MGLevel.h>

class TNSE_MGLevel1 : public TNSE_MGLevel
{
  protected:
#ifdef __2D__
    /** matrix A */
    TSquareMatrix2D *A;

    /** structure of matrix A */
    TSquareStructure2D *StructureA;

    /** matrix B1 */
    TMatrix2D *B1;

    /** matrix B2 */
    TMatrix2D *B2;

    /** structure of matrix B */
    TStructure2D *StructureB;

    /** structure of matrix BT */
    TStructure2D *StructureBT;
#endif  
#ifdef __3D__
    /** matrix A */
    TSquareMatrix3D *A;

    /** structure of matrix A */
    TSquareStructure3D *StructureA;

    /** matrix B1 */
    TMatrix3D *B1;

    /** matrix B2 */
    TMatrix3D *B2;

    /** matrix B3 */
    TMatrix3D *B3;

    /** structure of matrix B */
    TStructure3D *StructureB;

     /** matrix entries of matrix B3 */
    double *B3Entries;

    /** structure of matrix BT */
    TStructure3D *StructureBT;
#endif  

   /** row pointer for matrix A */
    int *ARowPtr;

    /** column number vector for matrix A */
    int *AKCol;

    /** matrix entries of matrix A */
    double *AEntries;

    /** row pointer for matrix B */
    int *BRowPtr;

    /** column number vector for matrix B */
    int *BKCol;

    /** matrix entries of matrix B1 */
    double *B1Entries;

    /** matrix entries of matrix B2 */
    double *B2Entries;

    /** row pointer for matrix BT */
    int *BTRowPtr;

    /** column number vector for matrix BT */
    int *BTKCol;
    

  public:
    /** constructor */
#ifdef __2D__
    TNSE_MGLevel1(int level, TSquareMatrix2D *A, 
                  TMatrix2D *B1, TMatrix2D *B2,
                  TStructure2D *structureBT,
                  double *f1, double *u1,
                  int n_aux, double *al,  int VelocitySpace, 
                  int PressureSpace, TCollection *coll,
                  int *dw);
#endif  
#ifdef __3D__
    TNSE_MGLevel1(int level, TSquareMatrix3D *A, 
                  TMatrix3D *B1, TMatrix3D *B2, TMatrix3D *B3, 
                  TStructure3D *structureBT,
                  double *f1, double *u1,
                  int n_aux, double *al,  int VelocitySpace, 
                  int PressureSpace, TCollection *coll,
                  int *dw);
#endif  

    /** destructor */
    ~TNSE_MGLevel1();

    /** calculate defect */
    virtual void Defect(double *u1, double *f1, double *d1, double &res);

    /** correct Dirichlet and hanging nodes */
    virtual void CorrectNodes(double *u1);

    /** Vanka smoother */
    virtual void  CellVanka(double *u1, double *rhs1, double *aux,
                            int N_Parameters, double *Parameters, 
                            int smoother, int N_Levels);

     /** Vanka smoother */
    virtual void NodalVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

   /** solve exact on this level */
    virtual void SolveExact(double *u1, double *rhs1);

    /** solve exact on this level */
    virtual void SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag);

   /** Braess--Sarazin smoother */
    virtual void BraessSarazin(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters,int N_Levels);

    /** step length control for Vanka */
    virtual double  StepLengthControl(double *u1, double *u1old, double *def1, 
           int N_Parameters, double *Parameter);

    /** print all matrices and both right hand sides */
    virtual void PrintAll();
    #ifdef _MPI
    virtual void UpdateHaloRhs(double*a, double*b)
    {
      
    }; 
    #endif
    
};

#endif
