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
   
#ifdef __2D__
  #include <FESpace2D.h>
  #include <TriaAffin.h>
  #include <QuadAffin.h>
  #include <NodalFunctional2D.h>
  #include <BoundEdge.h>
#endif
#ifdef __3D__
  #include <FESpace3D.h>
  #include <FEVectFunct3D.h>
#endif

#include <BaseCell.h>
#include <JointEqN.h>

#include <string.h>

double GetTime();
int GetMemory();

void SwapDoubleArray(double *doublearray, int length);
void SwapIntArray(int *intarray, int length);


#ifdef __2D__
// ====================================================================
// calculate the streamfunction from the (u1, u2) velocity field
// ====================================================================
void StreamFunction(TFESpace2D *velo, double *u1, double *u2,
                    TFESpace2D *stream, double *psi);

void ComputeVorticityDivergence(TFESpace2D *velo, TFEFunction2D *u1, TFEFunction2D *u2,
                      TFESpace2D *vorticity, double *vort, double *div);

// determine L2 and H1 error
void L2H1Errors(int N_Points, double *X, double *Y, double *AbsDetjk, 
                double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);

// determine L2-error, divergence error and H1 error, 2D
void L2DivH1Errors(int N_Points, double *X, double *Y, double *AbsDetjk, 
                   double *Weights, double hK, 
                   double **Der, double **Exact,
                   double **coeffs, double *LocError);

// determine L1 error, 2D
void L1Error(int N_Points, double *X, double *Y, double *AbsDetjk, 
	     double *Weights, double hK, 
	     double **Der, double **Exact,
	     double **coeffs, double *LocError);

// determine L2, H1 and SDFEM error
void SDFEMErrors(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError);

// determine L2, H1 and SDFEM error, in (0,P6)^2
void SDFEMErrorsSmooth(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError);

// determine L2, H1 and SDFEM error for smooth region in the
// example JohnMaubachTobiska1997 (x-0.5)^2+(y-0.5)^2 > r^2 
void SDFEMErrorsSmooth_JohnMaubachTobiska1997
(int N_Points, double *X, double *Y, double *AbsDetjk, 
 double *Weights, double hK, double **Der, double **Exact,
 double **coeffs, double *LocError);

// determine errors to interpolant
// paper with Julia Novo
 void SDFEMErrorsInterpolant(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError);
     
 // determine L2, H1 and SDFEM error for NSE
void SPGErrorsOseen(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError);

// determine L2, H1 and pressure part of SUPG error for Oseen
void SPGErrorsOseenPressure(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError);

// determine deformation tensor error
void DeformationTensorError(int N_Points, double *X, double *Y, double *AbsDetjk, 
                double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);
// compute the subgrid dissipation 
void SubGridDissipation(int N_Points, double *X, double *Y, double *AbsDetjk, 
                double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);

// determine H1 norm
void H1Norm(int N_Points, double *X, double *Y, double *AbsDetjk, 
            double *Weights, double hK, 
            double **Der, double **Exact,
            double **coeffs, double *LocError);

// compute the error in the divergence
void DivergenceError(int N_Points, double *X, double *Y,
		     double *AbsDetjk, double *Weights, double hK, 
		     double **Der, double **Exact,
		     double **coeffs, double *LocError);

// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *X, double *Y, double *AbsDetjk, 
           double *Weights, double hK, 
           double **Der, double **Exact,
           double **coeffs, double *LocError);

void DivergenceErrorGradDivOseen(int N_Points, double *X, double *Y,
         double *AbsDetjk, double *Weights, double hK, 
         double **Der, double **Exact,
         double **coeffs, double *LocError);
         
void Parameters_Gradient_Residual(int N_Points, double *X, double *Y, double *AbsDetjk,
           double *Weights, double hK,
           double **Der, double **Exact,
           double **coeffs, double *LocError);
#endif

double graddiv_parameterOseen(double hK, double nu, double b1, double b2);

#ifdef __3D__
void ComputeVorticityDivergence(TFESpace3D *velo, TFEFunction3D *u1, 
                                TFEFunction3D *u2, TFEFunction3D *u3,
                                TFESpace3D *vorticity_space, 
                                TFEFunction3D *vort1, 
                                TFEFunction3D *vort2, TFEFunction3D *vort3,
                                double *div);

// determine L2 and H1 error
void L2H1Errors(int N_Points, double *X, double *Y, double *Z,
                double *AbsDetjk, 
                double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);
void L2H1ErrorsSmooth(int N_Points, double *X, double *Y, double *Z,
                double *AbsDetjk, 
                double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError);

// compute L2 error, L2 error of divergence, and H1 error for vector valued
// basis functions (Raviart-Thomas or Brezzi-Douglas-Marini)
void L2DivH1Errors(int N_Points, double *X, double *Y, double *Z, 
                   double *AbsDetjk, double *Weights, double hK, double **Der,
                   double **Exact, double **coeffs, double *LocError);

// determine L1 error
void L1Error(int N_Points, double *X, double *Y,  double *Z, double *AbsDetjk, 
	     double *Weights, double hK, 
	     double **Der, double **Exact,
	     double **coeffs, double *LocError);

// determine deformation tensor error
void DeformationTensorError(int N_Points, double *X, double *Y, double *Z, 
                            double *AbsDetjk, 
                            double *Weights, double hK, 
                            double **Der, double **Exact,
                            double **coeffs, double *LocError);
// compute the subgrid dissipation 
void SubGridDissipation(int N_Points, double *X, double *Y, double *Z, 
                        double *AbsDetjk, 
                        double *Weights, double hK, 
                        double **Der, double **Exact,
                        double **coeffs, double *LocError);


// compute the error in the divergence
void DivergenceError(int N_Points, double *X, double *Y, double *Z,
		     double *AbsDetjk, double *Weights, double hK, 
		     double **Der, double **Exact,
		     double **coeffs, double *LocError);

// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *X, double *Y, double *Z,
                      double *AbsDetjk, 
                      double *Weights, double hK, 
                      double **Der, double **Exact,
                      double **coeffs, double *LocError);

void Q_criterion(TCollection *Coll,
                 TFEFunction3D *velocity1, TFEFunction3D *velocity2,
                 TFEFunction3D *velocity3, double *Qcrit);
#endif

// ====================================================================
// auxiliary routines
// ====================================================================
void linfb(int N_Points, double **Coeffs, double ** Params,
           TBaseCell *cell);
           
void ave_l2b_quad_points(int N_Points, double **Coeffs, double **Params,
           TBaseCell *cell);

void CFPM2D_linfb(int N_Points, double **Coeffs, double ** Params,
           TBaseCell *cell);

void LInfU(int N_Points, double **Coeffs, double **Params, 
           TBaseCell *cell);


#ifdef __2D__

// ====================================================================
// get velocity and pressure space from user defined parameters  
// ====================================================================
int GetVelocityAndPressureSpace(TCollection *coll,
                                BoundCondFunct2D *BoundCondition,
                                TCollection *mortarcoll,
                                TFESpace2D* &velocity_space,
                                TFESpace2D* &pressure_space,
                                int* pressure_space_code,
                                int velo_order, int pressure_order);


#endif

// ====================================================================
// read in solution from a Grape file 
// ====================================================================
int ReadGrapeFile(char *name, int N_FEFct, int N_FEVectFct,
                   TFEFunction2D **FEFct, TFEVectFunct2D **FEVectFct);

#ifdef __3D__

// ====================================================================
// get velocity and pressure space from user defined parameters  
// ====================================================================
int GetVelocityAndPressureSpace3D(TCollection *coll,
                                BoundCondFunct3D *BoundCondition,
                                TFESpace3D* &velocity_space,
                                TFESpace3D* &pressure_space,
                                int* pressure_space_code,
                                int velo_order, int pressure_order);

// ====================================================================
// get velocity and pressure space for MG_TYPE 2 
// ====================================================================
int GetVelocityAndPressureSpaceLow3D(TCollection *coll,
                                     BoundCondFunct3D *BoundCondition,
                                     TFESpace3D* &velocity_space,
                                     int *velocity_space_code,
                                     TFESpace3D* &pressure_space,
                                     int *pressure_space_code,
                                     int velo_order, int pressure_order);
// ====================================================================
// read in solution from a Grape file 
// ====================================================================
int ReadGrapeFile3D(char *name, int N_FEFct, int N_FEVectFct,
                    TFEFunction3D **FEFct, TFEVectFunct3D **FEVectFct);

#endif
void ExactNull(double x, double y, double z, double *values);
void ExactNull(double x, double y, double *values);
int ComputeNewTimeStep(double err); 
void BoundConditionNoBoundCondition(int BdComp, double t, BoundCond &cond);
void BoundConditionNoBoundCondition(double x, double y, double z, BoundCond &cond);
void BoundaryValueNoBoundaryValue(int BdComp, double Param, double &value);
void BoundaryValueHomogenous(int BdComp, double Param, double &value);
void BoundaryValueHomogenous(double x, double y, double z, double &value);
void BoundConditionVMM(int BdComp, double t, BoundCond &cond);
void BoundConditionNSE(int BdComp, double t, BoundCond &cond);
void BoundaryConditionPressSep(int i, double t, BoundCond &cond);
void BoundaryValuePressSep(int BdComp, double Param, double &value);
void BoundaryConditionPressSep3D(double x, double y, double z, BoundCond &cond);
void BoundaryValuePressSep3D(double x, double y, double z, double &value);
void BoundaryConditionNewton(double x, double y, double z, BoundCond &cond);
void BoundaryValueNewton(double x, double y, double z, double &value);
void BoundCondition_FEM_FCT(int i, double t, BoundCond &cond);
void BoundValue_FEM_FCT(int BdComp, double Param, double &value);
void BoundCondition_FEM_FCT(double x, double y, double z, BoundCond &cond);
void BoundValue_FEM_FCT(double x, double y, double z, double &value);
void BoundConditionAuxProblem(int i, double t, BoundCond &cond);
void BoundValueAuxProblem(int BdComp, double Param, double &value);
void ho_BoundCondition(int i, double t, BoundCond &cond);
void ho_BoundValue(int BdComp, double Param, double &value);

void SetPolynomialDegree();
void CheckMaximumPrinciple(TSquareMatrix *A, double *sol, int N_Active,
			   double *errors);
void SaveData(char *name, int N_Array, double **sol, int *N_Unknowns);
void ReadData(char *name, int N_Array, double **sol, int *N_Unknowns);

void SaveData(std::string basename, double *sol, int nDOF);
void ReadData(std::string filename, double *sol, int nDOF);
/******************************************************************************/
//
// sets the nodes with global dof neum_to_diri to Dirichlet nodes 
// by setting the matrix and the rhs
//
/******************************************************************************/
 
#ifdef __2D__
void SetDirichletNodesFromNeumannNodes(TSquareMatrix2D **SQMATRICES, 
				       double *rhs, double *sol,
				       int N_neum_to_diri,
				       int *neum_to_diri,
				       int *neum_to_diri_bdry,
				       double *neum_to_diri_param,
				       BoundValueFunct2D *BoundaryValue);
#endif    
#ifdef __3D__
void SetDirichletNodesFromNeumannNodes(TSquareMatrix3D **SQMATRICES, 
				       double *rhs, double *sol,
				       int N_neum_to_diri,
				       int *neum_to_diri,
				       double *neum_to_diri_x,
				       double *neum_to_diri_y,
				       double *neum_to_diri_z,
				       BoundValueFunct3D *BoundaryValue);
#endif    
