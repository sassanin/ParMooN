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
   
/** ************************************************************************ 
*
* @class     TSystemNSECST2D
* @brief     stores the information of a monolithic 2D NSE-CST system matrix 
* @author    Jagannath Venkatesan, 
* @date      23.02.16
* @History   
 ************************************************************************  */


#ifndef __SYSTEMNSECST2D__
#define __SYSTEMNSECST2D__

 #ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif
    

#include <SquareMatrix2D.h>

/**class for 2D NSE-CST system matrix */
class TSystemNSECST2D
{
  protected:

    /** DOFs of velocity, pressure, stress and deformation tensor space */    
    int N_U, N_P, N_S, N_D, N_TotalDof;
    int N_Active_U, N_Active_S, N_Active_D;
    int N_DirichletDof_U, N_DirichletDof_S, N_DirichletDof_D;
    int N_FeSpace, N_FeFunction;
    
      /** fespace */
    TFESpace2D *FeSpace[4];
      
        /** Fe functions */
    TFEFunction2D *FeFct[9];    
  
    /** Discretization type */
    int Disctype_NSECST;
    
    /** Solver type */
    int SOLVER;
               
    /** Bilinear coefficient   */
    CoeffFct2D *LinCoeffs[1];    

    /** sqstructure of the system matrices */
    TSquareStructure2D *sqstructureA, *sqstructureG, *sqstructureH;

    /** structure of the system matrices */
    TStructure2D *structureB, *structureBT, *structureC, *structureE,  *structureD,  *structureJ;
    
    
    TSquareMatrix2D *SqmatrixA11, *SqmatrixA12, *SqmatrixA21, *SqmatrixA22;
    TSquareMatrix2D *SqmatrixG11, *SqmatrixG12, *SqmatrixG21, *SqmatrixG22, *SqmatrixG23, *SqmatrixG32, *SqmatrixG33;
    TSquareMatrix2D *SqmatrixH11, *SqmatrixH22, *SqmatrixH33;
    
    TSquareMatrix2D *SQMATRICES[19];
    TSquareMatrix **sqmatrices;
    
    // LPS terms  
    TSquareMatrix2D *SqmatrixA11_LPS, *SqmatrixA12_LPS, *SqmatrixA21_LPS, *SqmatrixA22_LPS;
    TSquareMatrix2D *SqmatrixG11_LPS, *SqmatrixG22_LPS, *SqmatrixG33_LPS;
    TSquareMatrix2D *SqmatrixG12_LPS, *SqmatrixG21_LPS, *SqmatrixG23_LPS, *SqmatrixG32_LPS;
    TSquareMatrix2D *SqmatrixG11_LPS_Streamline, *SqmatrixG22_LPS_Streamline, *SqmatrixG33_LPS_Streamline;
   
    TMatrix2D *MatrixB1, *MatrixB2, *MatrixB1T, *MatrixB2T;
    TMatrix2D *MatrixC11, *MatrixC12, *MatrixC22, *MatrixC23;
    TMatrix2D *MatrixE11, *MatrixE12, *MatrixE22, *MatrixE23;
    TMatrix2D *MatrixD11, *MatrixD12, *MatrixD21, *MatrixD22, *MatrixD31, *MatrixD32;
    TMatrix2D *MatrixJ11, *MatrixJ21, *MatrixJ22, *MatrixJ32;
    TMatrix2D *MATRICES[22];
    TMatrix **matrices;
    
    /** aux is used to pass additional fe functions that is nedded for assembling */
    TAuxParam2D *aux, *auxerror_NSE, *auxerror_CST;
     
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[8];
  
     /** Boundary values */   
    BoundValueFunct2D *BoundaryValues[8];
        
    /** Discrete form for the equation */
    TDiscreteForm2D *DiscreteFormARhs, *DiscreteFormNL, *DiscreteFormRhs;

    
#ifdef _SMPI
   TSeqParDirectSolver *P_DS;
#endif
 

  public:
    /** constructor */
     TSystemNSECST2D(int N_FESpaces, TFESpace2D **FE_Spaces, int N_FEFunctions, TFEFunction2D **FeFunctions, int disctype_NSECST, int solver);

    /** destructor */
    ~TSystemNSECST2D();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D **BoundCond, BoundValueFunct2D **BoundValue,  TAuxParam2D *NSECSTaux, TAuxParam2D *NSEauxerror, TAuxParam2D *CSTauxerror);
 
    /** assemble the system matrix */
    void Assemble(double *sol, double *rhs);

    
    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
    
    // Residual calculation
      void GetResidual(double *sol, double *rhs, double *res);
    
    /** measure the errors */
    void MeasureErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP, DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3, double *u_error, double *p_error, double *s_error);

   
    
};

#endif
