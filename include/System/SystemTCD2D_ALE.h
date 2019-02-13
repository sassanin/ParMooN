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
* @class     TSystemTCD2D_ALE
* @brief     stores the information of a timedependent part of a 2D scalar system matrix 
* @author    Sashikumaar Ganesan
* @date      08.08.14
* @History 
 ************************************************************************  */


#ifndef __SYSTEMCD2D_ALE__
#define __SYSTEMCD2D_ALE__

#include <SquareMatrix2D.h>
#include <SystemTCD2D.h>

/**class for 2D scalar system matrix */
class TSystemTCD2D_ALE : public TSystemTCD2D
{
  protected:
    /** No. of Grid DOFs */
    int N_GridDOFs, N_GridActive, *GridKCol, *GridRowPtr;
    
    /** additional rhs for dG time disc */       
    double *rhs_Qp1; // dG(1)  
    
    /** Grid Posistions */
    double *MeshVelo, *gridpos, *gridpos_old, *gridpos_ref, *griddisp, *GridRhs, *Entries[4];   
    
    /** grid fespace */
    TFESpace2D *GridFESpace;    
    
    /** Fgrid BC */    
    BoundValueFunct2D *GridBoundValue[1];
     
    /** grid pos vector */
    TFEVectFunct2D *GridPos, *RefGridPos;
     
    /** Fe functions of NSE */
    TFEFunction2D *MeshVeloFct[2];    
    
    /** Discrete form for moving mesh */
    TDiscreteForm2D *DiscreteFormMARhs, *DiscreteFormGrid;  
    
    /** Old M mass matrix */
    TSquareMatrix2D *sqmatrixM_old;  
    
    /** matrices at time quad points for dG time disc */
    TSquareMatrix2D *sqmatrixM_Qp1, *sqmatrixA_Qp1;   // dG(1)   
    
    /** marices for the moving grid **/
    TSquareMatrix2D *SqmatrixG11, *SqmatrixG12, *SqmatrixG21, *SqmatrixG22, *SQMATRICES_GRID[4];
    
    /** structure for the moving grid **/
    TSquareStructure2D *SquareStructureG;
    
    /** aux for mesh */
    TAuxParam2D *Aux_ALE, *Meshaux;

    /** Grid bounadry conditions **/ 
    BoundCondFunct2D *GridBoundaryConditions[1];
    BoundValueFunct2D *GridBoundValues[1];
      
    /** method for Modify Mesh Coords */
    ModifyMeshCoords *ModifyCoord;
    
     /** method for Modify Bound Coords */
    ModifyBoundCoords *ModifyBoudary;   
    
    /** data for Modify Bound Coords  */
    int *N_MovVert;
    TVertex **MovBoundVert;
    TIsoBoundEdge **Free_Joint;
    double *Iso_refX;
    
    /** for mesh displacement */
    bool SolveLinearElastic, CONSERVATIVEALE, NeedInterMassMat;
    
  public:
    /** constructor */
    TSystemTCD2D_ALE(TFESpace2D *fespace, int disctype, int solver, TFESpace2D *gridFESpace, TFEVectFunct2D *MeshVelocity, 
                                                       bool conservativeale);

//     /** destrcutor */
//     ~TSystemMatTimeScalar2D();

    /** methods */
    void Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue, 
              CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D *GridBoundCond, BoundValueFunct2D *gridBoundValue,
              TAuxParam2D *aux);

    void AddMeshModifyFunction(ModifyMeshCoords *modifyCoord)
    {  ModifyCoord = modifyCoord; SolveLinearElastic = FALSE; }
    
    void AddBoundModifyFunction(ModifyBoundCoords *modifyboudary, int *n_MovVert, TVertex **movBoundVert, 
				TIsoBoundEdge **free_Joint, double * iso_refX);  
    
    
//     /** return the stiffness matric */
//     TSquareMatrix2D *GetAMatrix()
//     { return sqmatrixA; }

    /** store M mat for next time step */
    void StoreMmat();
    
    /** move mesh  */
    void MoveMesh(double Currtime);
    
    /** Get Mesh Velo */
    void GetMeshVelo(double Currtime, double tau, bool MoveMesh);
    
    /** assemble the Mesh mat and rhs */ 
    void AssembleMeshMat();
    
    /** assemble the Mass mat and rhs */
    void AssembleMRhs(double *sol, double *rhs);   
    
    /** assemble the Mass, stiff mat and rhs */
    void AssembleMARhs(double *sol, double *rhs); 

    /** assemble the M, A and rhs for different time disc schemes */
    void AssembleALEMat(double *sol, double *rhs, double tau);
    
    /** M = M + (tau*TDatabase::TimeDB->THETA1)*A */ 
    /** B = (tau*TDatabase::TimeDB->THETA1)*rhs +(tau*TDatabase::TimeDB->THETA2)*oldrhs + [ M - (tau*TDatabase::TimeDB->THETA2)A]*oldsol */  
    void AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol);
 
     /** solve the system matrix */
    void Solve(double *sol, double *rhs);  
    
//     /** return the residual of the system for the given sol*/
//     double GetResidual(double *sol);
    
};

#endif
