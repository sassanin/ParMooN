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
// @(#)MGLevel2D.h        1.3 08/24/99
//
// Class:       TMGLevel2D
// Purpose:     store all data for one level in a multi grid method
//
// Author:      Gunar Matthies 02.11.1998
//
// History:     02.11.1998 start of implementation
//
// =======================================================================

#ifndef __MGLEVEL2D__
#define __MGLEVEL2D__

#include <SquareMatrix2D.h>

class TMGLevel2D
{
  protected:
    /** level number */
    int Level;

    /** FE space */
    TFESpace2D *FESpace;

    /** permutation vector */
    int *Permutation;

    /** number of active nodes */
    int N_Active;

    /** upper bound for hanging node number */
    int HangingNodeBound;

    /** number of Dirichlet nodes */
    int N_Dirichlet;

    /** number of all degrees of freedom */
    int N_DOF;

    /** used matrix */
    TSquareMatrix2D *A;

    /** structure of used matrix */
    TSquareStructure2D *MatrixStructure;

    /** row pointer for matrix */
    int *RowPtr;

    /** column number vector */
    int *KCol;

    /** matrix entries */
    double *Entries;

    /** array with right-hand sides */
    double *Rhs;

    /** array with approximate solution */
    double *X;

    /** number of auxiliary vectors */
    int N_Aux;

    /** array of auxiliary vectors */
    double **Aux;

    /** array for additional data, e.g. ILU decomposition */
    double *Additional;

    /** generate ILU decomposition */
    void ILUDecomposition();

#ifdef _MPI
     /** number of all degrees of freedom in own cells*/
     int OwnN_DOF;   
    
     TParFECommunicator2D *ParComm;

     /** FEFunction in own+Hallo fe space */
     TFEFunction2D *C;     
     
     /** Own FE space */
     TFESpace2D *OwnScalarSpace;

     /** FEFunction in own fe space */
     TFEFunction2D *OwnC;
 
     /** Own solution */
     double *OwnSolArray;
#endif    
    
    
  public:
    /** constructor */
    TMGLevel2D(int level, TSquareMatrix2D *A,
               double *rhs, double *sol, int n_aux,
               int *permutation);

#ifdef _MPI   
     TMGLevel2D(int level, TSquareMatrix2D *a, double *rhs, double *sol, 
                TFEFunction2D *c, TParFECommunicator2D *parComm,TFESpace2D *ownScalarSpace, int n_aux,
                int *permutation);  
#endif    
    
    /** destructor */
    ~TMGLevel2D();

    /** return i-th auxiliary vector */
    double *GetAuxVector(int i);

    /** return FunctionVectors */
    double *GetSolution()
    { return X; }

    /** return Rhs */
    double *GetRhs()
    { return Rhs; }

    /** return AuxVectors */
    double **GetAuxVectors()
    { return Aux; }

    /** return number of degrees of freedom */
    int GetN_DOF()
    { return N_DOF; }

    /** get HangingNodeBound */
    int GetHangingNodeBound()
    { return HangingNodeBound; }

    /** get number of Dirichlet nodes */
    int GetN_Dirichlet()
    { return N_Dirichlet; }

    /** calculate defect */
    void Defect(double *sol, double *f, double *d, double &res);

    /** update solution */
    void Update(double *sol, double *upd);

    /** correct Dirichlet and hanging nodes */
    void CorrectNodes(double *vect);

    /** correct defect */
    void CorrectDefect(double *vect);

    /** reset vector to zero */
    void Reset(double *vect);

    /** return FE space */
    TFESpace2D *GetFESpace()
    { return FESpace; }

    /** smoother */
    void ILU(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void SOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void SSOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void Jacobi(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void Block2x2(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** solve exact on this level */
    void SolveExact(double *u1, double *rhs1);

    /** step length control */
    double StepLengthControl(double *u, 
                         double *uold, 
                         double *def,
                         int N_Parameters, 
                         double *Parameters);
    
#ifdef _MPI       
    TFESpace2D *GetOwnFESpace()
       { return OwnScalarSpace; }

    TParFECommunicator2D *GetParComm()
      { return ParComm; }

    double *GetOwnSolution()
      { return OwnSolArray; }  

     TFEFunction2D *GetFEFunction()
      { return C; }

     TFEFunction2D *GetOwnFEFunction()
      { return OwnC; }

    
#endif

};

#endif
