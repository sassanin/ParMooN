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
* @class     TSystemCST2D_Giesekus
* @brief     stores the information of a 2D CST system matrix of Giesekus type
* @author    Jagannath Venkatesan, 
* @date      03.09.15
* @History  
 ************************************************************************  */


#ifndef __SYSTEMCST2DGiesekus__
#define __SYSTEMCST2DGiesekus__

#include <SquareMatrix2D.h>

/**class for 2D  CST system matrix */
class TSystemCST2D_Giesekus
{
  protected:

    /** DOFs of stress space */    
    int N_S, N_Active, N_DirichletDof;
    
      /** fespace */
    TFESpace2D *FeSpace;
    
      /** Velocity fespace */
    TFESpace2D *FeSpaces_All[2];
    
        /** Fe functions */
    TFEFunction2D *FeFct[3];    
    
      /** Fe functions of Velocity */
    TFEFunction2D *FeFct_All[5]; 
    
    /** Discretization type */
    int Disctype;
    
    /** Solver type */
    int SOLVER;
    
       /** number of matrices in the system matrix*/
    int N_Matrices;
           
    /** Bilinear coefficient   */
    CoeffFct2D *LinCoeffs[1];    

    /** sqstructureA of the system matrix */
    TSquareStructure2D *sqstructure;
  
    /** S is the stiffness/system mat */
    TSquareMatrix2D *SqmatrixS11, *SqmatrixS12, *SqmatrixS21, *SqmatrixS22, *SqmatrixS23, *SqmatrixS32, *SqmatrixS33, *SQMATRICES[8];
    TSquareMatrix **sqmatrices;
  
    TAuxParam2D *CSTNSEaux, *CSTaux_error;
     
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[3];
  
     /** Boundary values */   
    BoundValueFunct2D *BoundaryValues[3];
        
    /** Discrete form for the equation */
    TDiscreteForm2D *DiscreteFormARhs;

    
  public:
    /** constructor */
     TSystemCST2D_Giesekus(TFESpace2D *stress_fespace, TFEVectFunct2D *Stress, int disctype, int solver);

    /** destructor */
    ~TSystemCST2D_Giesekus();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *S1BoundValue, BoundValueFunct2D *S2BoundValue, BoundValueFunct2D *S3BoundValue, TAuxParam2D *aux, TAuxParam2D *auxerror);
 
    /** assemble the system matrix */
    void Assemble(double *sol, double *rhs);
    
    /** get the resudual of the system */
    void GetResidual(double *sol, double *rhs, double *res);

    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
    
    /** measure the error */
    void MeasureErrors(DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3, double *s_error);

};

#endif
