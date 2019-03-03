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
* @brief     source file for TSystemTNSE3D_ALE
* @author    Sashikumaar Ganesan
* @date       6.4.2016
* @History    MoveBDwithVelo added 23.11.2017
 ************************************************************************  */

#ifndef __SYSTEMTNSE3D_ALE__
#define __SYSTEMTNSE3D_ALE__

#include <SystemNSE3D.h>
#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif
/**class for 3D  TNSE system matrix */
class TSystemTNSE3D_ALE : public TSystemNSE3D
{
  protected:
    
    /** vms projection fespace */
    TFESpace3D **Projection_Spaces;
    
    /** M - mass/system mat for TNSE velocity component   */
    TSquareMatrix3D **SqmatrixM11, **SqmatrixM12, **SqmatrixM13, 
                    **SqmatrixM21, **SqmatrixM22, **SqmatrixM23,
                    **SqmatrixM31, **SqmatrixM32, **SqmatrixM33;
  
    /** sqstructureG of the  vms projection matrix */
    TSquareStructure3D *sqstructureL;

    /** structure of the vms projection  matrix */
    TStructure3D *structure_G, *structure_tilde_G; 
    
    /** G -  mat for VMS   */    
    TMatrix3D *matrix_tilde_G11, 
              *matrix_tilde_G22,
	      *matrix_tilde_G33,
              *matrix_G11,
	      *matrix_G22,
	      *matrix_G33;
	      
    TMatrix3D **Matrices_tilde_G11,
              **Matrices_tilde_G22,
	      **Matrices_tilde_G33, 
	      **Matrices_G11,
	      **Matrices_G22,
	      **Matrices_G33;  

    /** L -  mat for VMS   */
    TSquareMatrix3D *sqmatrixL, **MatricesL;
     
    // Assembling rhs*/
    TAssembleMat3D *RhsOnlyAssemble;
    
#ifdef _SMPI
   TSeqParDirectSolver* P_DS;
#endif
    
     /** working rhs, used in AssembleSystMat() */
    double *B;   
   
    /** to store defect */
    double *defect;   
    
    /** factor that multplied with Mat A in working rhs */
    double gamma;      

    /** Discrete form of the M and rhs matrics */
    TDiscreteForm3D *DiscreteFormRhs; 
    
    /** NSE_Rhsaux is used to for assembling rhs only*/
    TAuxParam3D *NSE_Rhsaux;
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
    
    /** needed for error calculation in time */
    double olderror_l_2_l_2u;
    
      
    //-------------------------used for ALE-------------------------------
     /** No. of Grid DOFs */
    int N_GridDOFs, N_GridActive, *GridKCol, *GridRowPtr;
    
    TFEVectFunct3D **MeshVelocity;
    
    /** Grid Posistions */
    double **MeshVelo, *gridpos, *gridpos_old, *gridpos_ref,  *GridRhs, *Entries[9],**Mesh_restrict_aux;   
    
    /** grid fespace */
    TFESpace3D **GridFESpace;    
    
     
    /** mesh matrices assembling */
    TAssembleMat3D *MeshMatAssemble;    
    
    /** grid pos vector */
    TFEVectFunct3D *GridPos, *RefGridPos;
     
//     /** Fe functions of NSE */
//     TFEFunction3D *MeshVeloFct[3];    
    
    /** Discrete form for moving mesh */
    TDiscreteForm3D *DiscreteFormMARhs, *DiscreteFormGrid;  
           
    /** marices for the moving grid **/
    TSquareMatrix3D *SqmatrixG11, *SqmatrixG12, *SqmatrixG13,
                    *SqmatrixG21, *SqmatrixG22, *SqmatrixG23,
		    *SqmatrixG31, *SqmatrixG32, *SqmatrixG33,
		    *SQMATRICES_GRID[9];
 
    /** structure for the moving grid **/
    TSquareStructure3D *SquareStructureG;
    
    /** aux for mesh */
    TAuxParam3D *Aux_ALE, *Meshaux;

    /** Grid bounadry conditions **/ 
    BoundCondFunct3D *GridBoundaryConditions[3];
    BoundValueFunct3D *GridBoundValues[3];
      
    /** method for Modify Mesh Coords */
    ModifyMeshCoords_3D *ModifyCoord;
    
     /** method for Modify Mesh Coords */
    ModifyBoundCoords_3D *ModifyBoudary;   
    MoveBound_3D  *MoveBDwithVelo;
    
    
    /** */
    bool SolveLinearElastic, CONSERVATIVEALE, BdMoveWithVelo;
    
    /** info of moving boundaries */
    int N_MovVert, N_Movfaces, *Movfaces;
    TVertex **MovBoundVert;
    TBaseCell ** MovCells;
    

    /**mesh mat LU decomp values */
    void *MeshMatSolver_Symbolic, *MeshMatSolver_Numeric;
    double *MeshMatSolver_Values;
    int *MeshMatSolver_KCol, *MeshMatSolver_Row;
    //------------------------------------------------------------------------------
    
  public:
    /** constructor */
     TSystemTNSE3D_ALE(int N_levels, TFESpace3D **velocity_fespace, TFESpace3D **presssure_fespace, TFEVectFunct3D **velocity, 
                       TFEFunction3D **pressure, double **sol, double **rhs,  int disctype, int nsetype, int solver,
                       TFESpace3D **Projection_space, 
                       TFESpace3D ** gridFESpace, TFEVectFunct3D **meshVelocity, bool conservativeale);
     
    /** destrcutor */
    ~TSystemTNSE3D_ALE();

    /** methods */
    
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
              BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue, BoundCondFunct3D *GridBoundCondition,
              BoundValueFunct3D *GridBoundValue,CoeffFct3D *GridBilinearCoeffs);
    
    void AddMeshModifyFunction(ModifyMeshCoords_3D *modifyCoord)
       {  ModifyCoord = modifyCoord; SolveLinearElastic = FALSE; }
  
    /** add info of moving boundaries */
    void AddBoundModifyFunction(ModifyBoundCoords_3D *modifyboudary, int n_movVert, TVertex **movboundVert,
                                int n_movfaces, int *movfaces, TBaseCell ** movcells);

    void AddMoveBoundFunction(MoveBound_3D *movebdwithvelo, int n_movVert, TVertex **movboundVert,
                                int n_movfaces, int *movfaces, TBaseCell ** movcells);
    
    /** assemble the M, A and rhs */
    void Assemble();
        
    /** assemble only the rhs of NSE system */
    void AssembleRhs(); 
    
    /** Get Mesh Velo */
    void GetMeshVelo(bool MoveMesh, int solver_flag);
    
    /** Get Mesh Velo */
    void GetMeshVeloAndMove(double Currtime, double tau);
    
    /** assemble the Mesh mat and rhs */ 
    void AssembleMeshMat();
    
    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol);

    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMatNonLinear();
    
     /** restoring the mass matrix */
    void RestoreMassMat();
    
    /** restoring the mass matrix */
    void RestoreMassMatNonLinear();
      
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear();

    /** solve the system matrix */
    void  Solve(double *sol);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double &impuls_residual, double &residual);
 
    /** measure the error in the NSE */
    void MeasureTNSEErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2, DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP, double *AllErrors);

    /** print the matrices in a file */
    void printall_matrix();
    
    void GetMesh();
    
    TFEVectFunct3D *GetGridPosFeVect()
    { return GridPos; }
    
    void All_levels_check();
    
    double BlockMatVect(TSquareMatrix *A);
    
    double BlockMatVect(TMatrix *A);

    void CheckAllMat();

    void DOF_stats(TMatrix3D* MatB, char M , int k,char * name);

    void RHS_stats();
    
    void SetSystMatAssembled(bool val)
     {SystMatAssembled=val;}
    
};

#endif
