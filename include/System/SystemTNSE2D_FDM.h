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
* @class     TSystemTNSE2D_FDM
* @brief     stores the information of a 2D TNSE system matrix for Fictitious domain method
* @author    Sashikumaar Ganesan
* @date      06.07.16
* @History    
 ************************************************************************  */

#ifndef __SYSTEMTNSE2D_FDM__
#define __SYSTEMTNSE2D_FDM__

#include <SystemTNSE2D.h>

/**class for 2D TNSE system matrix */
class TSystemTNSE2D_FDM : public TSystemTNSE2D
{
  protected:
    /** No. of Grid DOFs */
    int N_SolidCells, N_SolidDOF;
      
    /** Level Set fespace */
    TFESpace2D *LevelSetFESpace;    
 
    /** Level Set functions of NSE */
    TFEFunction2D *LevelSetFunction;    
    
    /** Solid domain cells */
    TBaseCell **Solid_Cells;
    TCollection *Solid_Coll;
    
    /** Lagrange fespace */
    TFESpace2D *SolidVelocity_FeSpace, *LagrangeFESpace;     
    
    /** structure for the moving grid **/
    TStructure2D *structureCT, *structureC;
    
    /** lagrangian matrices */
    TMatrix2D *MatrixC1, *MatrixC1T, *MatrixC2, *MatrixC2T;

    /** other row vectors */ 
    double *E1, *E2, *F1, *F2, D1, D2;
  
//     /** */
//     bool SolveLinearElastic, CONSERVATIVEFDM, NeedInterMassMat;
    
  public:
    
    /** constructor */
    TSystemTNSE2D_FDM(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                                     TFEFunction2D *p, double *sol, double *rhs,  int disctype, int nsetype, int solver,                                    
                                     TFESpace2D *levelSetFESpace,   TFEFunction2D *levelSetFunction);
// 
//     /** destrcutor */
//     ~TSystemTNSE2D_FDM();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
              TAuxParam2D *aux, TAuxParam2D *nseaux_error, TFESpace2D *solidVelocity_FeSpace, TFESpace2D *lagrangeFESpace);
    
    /** partition the mesh for solid domain based on Level set function */
    TCollection *GetSolidColl();
    
//     void AddMeshModifyFunction(ModifyMeshCoords *modifyCoord)
//     {  ModifyCoord = modifyCoord; SolveLinearElastic = FALSE; }
//     
//     void AddBoundModifyFunction(ModifyBoundCoords *modifyboudary)
//     {  ModifyBoudary = modifyboudary; SolveLinearElastic = TRUE; }  
 
    /** assemble the M, A and rhs */
    void Assemble(double *sol, double *rhs);
  
// //     /** store M mat for next time step */
// //     void StoreMmat();
//     
// //     /** move mesh  */   
// //     void MoveMesh(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
// //                                              double * Iso_refX, double Currtime);
// //     
// //     /** move mesh  */
// //     void MoveMesh(double Currtime);
//     
//     /** Get Mesh Velo */
//     void GetMeshVeloAndMove(int *N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
//                      double * Iso_refX,  double Currtime, double tau);
//     
//     /** Get Mesh Velo */
//     void GetMeshVeloAndMove(double Currtime, double tau);
//     
//     /** assemble the Mesh mat and rhs */ 
//     void AssembleMeshMat();
//     
//     /** scale B matices and assemble rhs based on the \theta scheme  */
//     void AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol);
// 
//     /** scale B matices and assemble rhs based on the \theta scheme  */
//     void AssembleSystMatNonLinear();
//      
//     /** restoring the mass matrix */
//     void RestoreMassMat();
//         
//     /** assemble the nonlinear part of the NSE system */
//     void AssembleANonLinear(double *sol, double *rhs);
// 
//     /** solve the system matrix */
//     void  Solve(double *sol);
//     
//     /** get the resudual of the NSE system */
//     void GetTNSEResidual(double *sol, double *res);   
//  
// //     /** measure the error in the NSE */
// //     void MeasureTNSEErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP, double *AllErrors);
// //     
// //     void GridMovingScheme(TFESpace2D *GridFESpace);
};

#endif
