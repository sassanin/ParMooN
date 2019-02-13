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
* @class     TSystemNSE2D
* @brief     stores the information of a 2D NSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History   Added methods (Sashi, 23.08.14)
 ************************************************************************  */


#ifndef __SYSTEMNSE2D__
#define __SYSTEMNSE2D__

#include <SquareMatrix2D.h>
#include <AssembleMat2D.h>

/**class for 2D  NSE system matrix */
class TSystemNSE2D
{
  protected:

    /** DOFs of velocity and pressure spaces */    
    int N_U, N_P, N_Active, N_DirichletDof;
    
    /** sol, rhs arrays */
    double *Sol, *Rhs, *RHSs[3];
           
    /** velocity fespace */
    TFESpace2D *FeSpaces[5];
    
    TFEVectFunct2D *VelocityFct;
    
    /** Fe functions of NSE */
    TFEFunction2D *FeFct[5];    
    
    /** Discretization type */
    int Disctype;
       
    /** NSE type */
    int NSEType; 
           
    /** Bilinear coefficient   */
    CoeffFct2D *LinCoeffs[1];    
    
    /** NSEaux is used to pass additional fe functions (eg. mesh velocity) that is nedded for assembling */
    TAuxParam2D *NSEaux, *NSEaux_error;
    
    /** method for resudual calculation */
    DefectProc *Defect;   
    
    /** Solver type */
    int Solver;
       
    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructureA of the system matrix */
    TSquareStructure2D *sqstructureA;

    /** structure of the system matrix */
    TStructure2D *structureB, *structureBT;
    
    /** A is the stiffness/system mat for NSE velocity component   */
    TSquareMatrix2D *SqmatrixA11, *SqmatrixA12, *SqmatrixA21, *SqmatrixA22, *SQMATRICES[9];
    TSquareMatrix **sqmatrices;
  
    /** B is the  system mat for NSE pressure component   */
    TMatrix2D *MatrixB1, *MatrixB2, *MatrixB1T, *MatrixB2T, *MATRICES[8];
    TMatrix **matrices;
    
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[2];
  
     /** Boundary values */   
    BoundValueFunct2D *BoundaryValues[2];
        
    /** Discrete form for the equation */
    TDiscreteForm2D *DiscreteFormARhs, *DiscreteFormNL, *DiscreteFormRhs;

    /** Assembling */  
    TAssembleMat2D *AMatRhsAssemble, *AMatAssembleNonLinear;  
   
   
        /** these private routines are not available for public */
#ifdef __PRIVATE__  
    /** sqstructureG of the  vms projection matrix */
    TSquareStructure2D *sqstructureL;

    /** structure of the vms projection  matrix */
    TStructure2D *structure_G, *structure_tilde_G;      
    
    /** L -  mat for VMS   */
    TSquareMatrix2D *MatricesL;
    
    /** G -  mat for VMS   */
    TMatrix2D *Matrices_tilde_G11, *Matrices_tilde_G22, *Matrices_G11, *Matrices_G22;  
#endif 
 
  private:
   void UpdateUpwind();
    
   void UpdateLPS();
  public:
    /** constructor */
     TSystemNSE2D(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                  TFEFunction2D *p, double *sol, double *rhs, int disctype, int nsetype, int solver
#ifdef __PRIVATE__  
                  ,TFESpace2D *Projection_space, TFESpace2D *Stress_FeSpace, TFESpace2D *Deformation_FeSpace
#endif         
                );

    /** destrcutor */
    ~TSystemNSE2D();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
              TAuxParam2D *aux, TAuxParam2D *auxerror) ;
    
    
    /** assemble the system matrix */
    void Assemble(double *sol, double *rhs);
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(double *sol, double *rhs);
    
    /** assemble only the rhs when coupled with CST */
    void AssembleRhsOnly(double *sol, double *rhs);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double *rhs, double *res);
    
    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
  
    /** measure the error in the NSE */
    void MeasureErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    double *u_error, double *p_error);
};

#endif
