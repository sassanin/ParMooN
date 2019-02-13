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
* @class     TSystemTNSECST2D_ALE
* @brief     stores the information of a monolithic 2D T-NSE-CST-ALE system matrix 
* @author    Jagannath Venkatesan, 
* @date      21.03.17
* @History   
 ************************************************************************  */

#ifndef __SYSTEMTNSECST2D_ALE__
#define __SYSTEMTNSECST2D_ALE__

#include <SystemTNSECST2D.h>

/**class for 2D  TNSE system matrix */
class TSystemTNSECST2D_ALE : public TSystemTNSECST2D
{
  protected:
    /** No. of Grid DOFs */
    int N_GridDOFs, N_GridActive, *GridKCol, *GridRowPtr;
    
    /** Grid Posistions */
    double *MeshVelo, *gridpos, *gridpos_old, *gridpos_ref, *gridpos_aux, *griddisp, *GridRhs, *Entries[4], *tmp_Gd, *tmp_Gsol;     
    
    /** grid fespace */
    TFESpace2D *GridFESpace;    
    
    /** Fgrid BC */    
    BoundValueFunct2D *GridBoundValue[1];
     
    /** grid pos vector */
    TFEVectFunct2D *GridPos, *RefGridPos, *AuxGridPos, *MeshVectFct;
     
    /** Fe functions of NSE */
    TFEFunction2D *MeshVeloFct[2];    
    
    /** Discrete form for moving mesh */
    TDiscreteForm2D *DiscreteFormMARhs, *DiscreteFormMesh;  
   
       
    /** matices for the moving grid **/
    TSquareMatrix2D *SqmatrixGrid11, *SqmatrixGrid12, *SqmatrixGrid21, *SqmatrixGrid22, *SQMATRICES_GRID[4];
 
    /** matrices for adding freesurf/interface entries to A11 and A22 **/
    TSquareMatrix2D *SqmatrixF11, *SqmatrixF22;
    
    
    /** structure for the moving grid **/
    TSquareStructure2D *SquareStructureGrid;
    
    /** aux for mesh */
    TAuxParam2D *Aux_ALE, *Meshaux;

    /** Grid bounadry conditions **/ 
    BoundCondFunct2D *GridBoundaryConditions[1];
    BoundValueFunct2D *GridBoundValues[1];
      
    /** method for Modify Mesh Coords */
    ModifyMeshCoords *ModifyCoord;
    
     /** method for Modify Mesh Coords */
    ModifyBoundCoords *ModifyBoudary;   
    
    bool SolveLinearElastic;
    
    
  public:
    
    /** constructor */
    TSystemTNSECST2D_ALE(int N_FESpaces, TFESpace2D **FE_Spaces, int N_FEFunctions, TFEFunction2D **FeFunctions, int disctype_NSECST, int solver, TFESpace2D *gridFESpace, TFEVectFunct2D *MeshVelocity);

    /** destrcutor */
    ~TSystemTNSECST2D_ALE();
    
     /** methods */
    /** Initilize the discrete forms and the matrices */ 
    void Init(CoeffFct2D *LinCoeffs_Monolithic, CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D **BoundCond, BoundValueFunct2D **BoundValue,  
			   TAuxParam2D *TNSECSTaux, TAuxParam2D *TNSEauxerror, TAuxParam2D *TCSTauxerror, BoundCondFunct2D *GridBoundCond, BoundValueFunct2D *gridBoundValue);
       
    void AssembleMeshMat();

     void Assemble(double *sol, double *rhs, double dt);
     
     void AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol);
     
     void GetTNSECSTResidual(double *sol, double *res);
     
     void Solve(double *sol);
     
     void RestoreMassMat();
     
     void AssembleNonLinear(double *sol, double *rhs);
     
     void AssembleSystMatNonLinear();
     
     void GetMeshVelo(double tau, TVertex ***MovBoundVert, int *N_MovVert, TBaseCell **Free_Cells, int **IsoCellEdgeNos, bool &reparam, TFEVectFunct2D **VelocityFct);

     void MoveMesh(double tau, TVertex ***MovBoundVert, int *N_MovVert, TBaseCell **Free_Cells, int **IsoCellEdgeNos, bool &reparam, int N_ReParam, TFEVectFunct2D **VelocityFct);
  
};

#endif
