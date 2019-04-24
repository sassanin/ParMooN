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
* @class     TSystemHyperElast3D
* @brief     stores the information of a hyperelastic system matrix 
* @author    Sashikumaar Ganesan, 
* @date      02.03.19
* @History    
 ************************************************************************  */


#ifndef __SYSTEMTHYPERELAST3D3D__
#define __SYSTEMTHYPERELAST3D3D__

#include <SquareMatrix3D.h>
#include <FEVectFunct3D.h>
// #include <NSE_MultiGrid.h>
#include <ItMethod.h>
#include <AssembleMat3D.h>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#include <ParDirectSolver.h>
#endif

#ifdef _OMPONLY
#include <ParDirectSolver.h>
#endif

/**class for 3D  hyperelastic system matrix */
class TSystemHyperElast3D
{
  protected:

#ifdef _MPI
    TParFEMapper3D **ParMapper_U, **ParMapper_P;
    TParFECommunicator3D **ParComm_U, **ParComm_P;
    MPI_Comm Comm;
    TParDirectSolver *DS;
#endif
    
#ifdef _OMPONLY
    TParDirectSolver *DS;
#endif
    
    /** Number of multigrid levels */
    int N_Levels;
    
    /** starting level for the solver, e.g. Start_Level = N_Levels-1 for direct solver */
    int Start_Level;  
    
    /** DOFs of displacement spaces */    
    int N_TotalDOF, N_U, N_Active, N_DirichletDof;
    
    /** Number of free surface/interface faces */
    int *N_FreeSurfFaces;
    
    /** Cell numbers and Joint numbers of free surface/interface faces */    
    int **FreeSurfCellNumbers, **FreeSurfJointNumbers;
    
    /** disp fespace */
    TFESpace3D **U_Space, *FeSpaces[5];
    
    /** disp FE function */
    TFEVectFunct3D **Displacement;
    
    double **SolArray, **RhsArray, *RHSs[3];

    TFESpace3D *fesp[1], *fesprhs[3], *fesp_aux[1];
    TFEFunction3D *fefct[7], **fefct_aux;
  
    // Assembling */
    TAssembleMat3D **AMatRhsAssemble, **AMatAssembleNonLinear;
    
    /** Discretization type */
    int Disctype;
                 
    /** Bilinear coefficient   */
    CoeffFct3D *LinCoeffs[1];    
    
    /** aux is used to pass additional fe functions (eg. mesh velocit, nonlinear functions) that is nedded for assembling */
    TAuxParam3D **Hyperaux, **Hyperaux_error;
    
    /** method for matrix vector mult */
    MatVecProc *MatVect;  
    
    /** method for resudual calculation */
    DefectProc *Defect;   
    
    /** Solver type */
    int SOLVER;
       
    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructureA of the system matrix */
    TSquareStructure3D **sqstructureA;
    
    /** A is the stiffness/system mat for NSE velocity component   */
    TSquareMatrix3D **SqmatrixA11, **SqmatrixA12, **SqmatrixA13;
    TSquareMatrix3D **SqmatrixA21, **SqmatrixA22, **SqmatrixA23;
    TSquareMatrix3D **SqmatrixA31, **SqmatrixA32, **SqmatrixA33, *SQMATRICES[15];
    TSquareMatrix **sqmatrices;
    TSquareMatrix3D **SqmatrixF11, **SqmatrixF22, **SqmatrixF33;
        
    /** Boundary conditon */
    BoundCondFunct3D *BoundaryConditions[3];
  
     /** Boundary values */   
    BoundValueFunct3D *BoundaryValues[3];
        
    /** Discrete form for the equation */
    TDiscreteForm3D *DiscreteFormARhs, *DiscreteFormNL;

//     /** variables for multigrid */
//     int N_aux;
//     double Parameters[2], *Itmethod_sol, *Itmethod_rhs;
//     TNSE_MultiGrid *MG;
//     TNSE_MGLevel *MGLevel;
//     TItMethod *Itmethod, *prec;
 
    
  private:
//    void UpdateUpwind(int i);
    
//    void UpdateLPS();
      void InitHyperDiscreteForms(TDiscreteForm3D *DiscreteFormGalerkin);
      void InitHyperAuxParm(int i);
  public:
    /** constructor */
     TSystemHyperElast3D(int N_levels, TFESpace3D **disp_fespace, TFEVectFunct3D **displacement, 
                     double **sol, double **rhs, int disctype, int solver);

    /** destrcutor */
//     ~TSystemMatNSE3D ();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
              BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue);
    
    /** assemble the system matrix */
    void Assemble();
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(double **sol, double **rhs);
    
//     /** get the resudual of the NSE system */
//     void GetResidual(double *sol, double *rhs, double *res, double &impuls_residual, double &residual);
//     
//     /** solve the system matrix */
//     void  Solve(double *sol, double *rhs);
//   
//     /** measure the error in the NSE */
//     void MeasureErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2, DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP,
//                        double *u_error, double *p_error);
//     
//     /** find all joints which form the free surface/interface */   
//     void FindFreeSurfJoints(int level, int Region_ID);
};

#endif
