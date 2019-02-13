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
// @(#)NS2DErrorEstimator.h        1.3 11/15/99
// 
// Class:       TNS2DErrorEstimator
//
// Purpose:     define error estimator
//
// Author:      Volker John
//
// History:     29.01.98 start implementation
//
// =======================================================================
#ifdef __2D__

#ifndef __NS_ERROR_ESTIMATOR__
#define __NS_ERROR_ESTIMATOR__

#include <FEFunction2D.h>
#include <FEVectFunct2D.h>

#define ns_gradient_indicator                      0
#define ns_residual_estimator_h1                   1
#define ns_residual_estimator_l2                   2
#define ns_residual_estimator_energy_quasi_robust  3
#define ns_gradient_recovery                       4
#define ns_implicit_estimator_neumann              5

class TNS2DErrorEstimator
{
  protected:
    /** cell collection which the spaces are based on */
    TCollection *Collection_U;
    TCollection *Collection_P;
    
    /** ansatz space */
    TFESpace2D *FESpace2D_U;
    TFESpace2D *FESpace2D_P;

    /** FE function */
    TFEVectFunct2D *U;
    TFEFunction2D *P;

    /** error estimator routine */
    int FELocalEstimator;

    /** do or not do error control */
    int ErrorControl;
 
    /** Stokes or Navier-Stokes equation */
    int NavierStokes;

public:
    /** initialize error estimator */
    TNS2DErrorEstimator(int fe_local_estimator, 
                      TFEVectFunct2D *u, 
                      TFEFunction2D *p, 
                      int error_control,
                      int navierstokes);

    /** return Collection */
    TCollection *GetCollection_U()
      { return Collection_U; };
    TCollection *GetCollection_P()
      { return Collection_P; };

    /** return */
    TFEVectFunct2D *GetU()
      { return U;};

    TFEFunction2D *GetP()
      { return P;};

    /** return element error estimator */
    int GetFELocalEstimator()
      { return FELocalEstimator;};

    /** return error control */
    int GetErrorControl()
      { return ErrorControl;}

    /** return problem type */
    int GetNavierStokes()
      { return NavierStokes;}

    /** the estimator */
    void GetErrorEstimate(int N_Derivatives,
                          MultiIndex2D *NeededDerivatives,
                          int N_DerivativesP,
                          MultiIndex2D *NeededDerivativesP,
                          CoeffFct2D *Coeff, 
                          BoundCondFunct2D **BoundaryConds,
                          BoundValueFunct2D **BoundaryValues,
                          TAuxParam2D *Aux,
                          int n_fespaces,
                          TFESpace2D **fespaces,
                          double *eta_K, 
                          double *eta_max,
                          double *estimated_global_error);
      
    /** elementwise estimator */
    void  EstimateCellError(TFESpace2D **fespaces,
                            TBaseCell *cell,
                            int N_Points, 
                            double *X, 
                            double *Y, 
                            double *AbsDetjk, 
                            double *weights,
                            double **Derivatives, 
                            double **AuxArray, 
                            BoundCondFunct2D **BoundaryConds,
                            BoundValueFunct2D **BoundaryValues,
                            int N_Points1D, 
                            double *zeta,
                            double *X1D[4], 
                            double *Y1D[4], 
                            double *weights1D,
                            double *xyval_ref1D[4],
                            double *xderiv_ref1D[4],
                            double *yderiv_ref1D[4],
                            int *GlobalNumbers, 
                            int *BeginIndex,
                            int *DOF,
                            double *Values,                                           
                            int *GlobalNumbersP, 
                            int *BeginIndexP,
                            int *DOFP,
                            double *ValuesP,
                            double *local_error);

};

#endif

#endif // #ifdef __2D__
