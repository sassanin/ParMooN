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
// LocalProjection.h
//
// Purpose:   routines for local projection stabilization
//
// Author:    Gunar Matthies  2007/03/06
//
// =======================================================================

#ifndef __LOCAL_PROJECTION__
#define __LOCAL_PROJECTION__

#include <LinAlg.h>

#ifdef __2D__
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *b, double *r);

void Defect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x,
                  double *b, double *r);

void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *y);

void MatVect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void UltraLocalProjection(void* A, bool ForPressure, CoeffFct2D *Coeff);
void UltraLocalProjection(void* A, bool ForPressure);
void AddDeformationTensorTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff);
void UltraLocalProjectionSD(void* A, bool ForPressure);

double UltraLocalError(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorDivergence(TFEFunction2D *uh1, TFEFunction2D *uh2,
                       DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2,
                       double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorStreamline(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorStreamlinePWConst(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff);

void AddStreamlineTerm(TSquareMatrix2D* A, TFEFunction2D *uh1,
                       TFEFunction2D *uh2,
                       double lpcoeff, double lpexponent, int OrderDiff);

void AddStreamlineTermPWConst(TSquareMatrix2D* A, TFEFunction2D *uh1,
                              TFEFunction2D *uh2,
                              double lpcoeff, double lpexponent, int OrderDiff);

void AddDivergenceTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff);

// for conformation stress tensor equation (div of stress tensor)
void AddDivergenceTerm(TSquareMatrix2D **SQMATRICES,
                       double lpcoeff, double lpexponent, int OrderDiff);

// for conformation stress tensor equation (grad of stress tensor)
void AddGradTauTerm(TSquareMatrix2D *G11,TSquareMatrix2D *G22,
                       TSquareMatrix2D *G33, double lpcoeff, double lpexponent, int OrderDiff);

// for conformation stress tensor equation (grad and div of stress tensor)
void AddTauTerm(TSquareMatrix2D *G11, TSquareMatrix2D *G12, TSquareMatrix2D *G21,
		TSquareMatrix2D *G22, TSquareMatrix2D *G23, TSquareMatrix2D *G32,
                       TSquareMatrix2D *G33, double lpcoeff_grad, double lpcoeff_div, double lpexponent, int OrderDiff);

// for conformation stress tensor equation in 2 phase flows (grad and div of stress tensor)
void AddTauTerm_2PhaseOrImpDropFlow(TSquareMatrix2D *G11, TSquareMatrix2D *G12, TSquareMatrix2D *G21,
		TSquareMatrix2D *G22, TSquareMatrix2D *G23, TSquareMatrix2D *G32,
                       TSquareMatrix2D *G33, double lpcoeff_grad, double lpcoeff_div, double lpexponent, int OrderDiff);

// for conformation stress tensor equation in 3D-Axis-symmetric form (grad and div of stress tensor)
void AddTauTerm_2PhaseOrImpDropFlow_3DAxial(TSquareMatrix2D *G11, TSquareMatrix2D *G12, TSquareMatrix2D *G21,
		TSquareMatrix2D *G22, TSquareMatrix2D *G23, TSquareMatrix2D *G32, TSquareMatrix2D *G33, 
		       double lpcoeff_grad, double lpcoeff_div, double lpexponent, int OrderDiff);

void AddDeformationTensorTerm_2PhaseOrImpDropFlow(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpexponent, int OrderDiff);

void AddDeformationTensorTerm_2PhaseOrImpDropFlow_3DAxial(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpexponent, int OrderDiff);


void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    TMatrix *C,
                    double *x, double *y);

void MatVect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   TMatrix *C,
                   double *x, double *b, double *r);

void Defect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

void AddStreamlineTerm(TSquareMatrix2D *G11,TSquareMatrix2D *G22,
                       TSquareMatrix2D *G33, TFEFunction2D *u1, TFEFunction2D *u2,
                       double lpcoeff, double lpexponent, int OrderDiff);

// for conformation stress tensor equation
void AddStreamlineTerm_2PhaseOrImpDropFlow_3DAxial(TSquareMatrix2D *G11,TSquareMatrix2D *G22,
                       TSquareMatrix2D *G33, TFEFunction2D *u1, TFEFunction2D *u2,
                       double lpcoeff, double lpexponent, int OrderDiff);

void AddStretchingTerm(TSquareMatrix2D **SQMATRICES, TFEFunction2D *u1, TFEFunction2D *u2,
                       double lpcoeff, double lpexponent, int OrderDiff);


double UltraLocalErrorSmooth(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                             double lpcoeff, double lpexponent, int OrderDiff);
bool TestCell(TBaseCell *cell);

void UltraLocalProjectionFunction(void* A, bool ForPressure);

void UltraLocalProjectionStreamlinePLaplacian(TSquareMatrix2D* A, 
                                              TFEFunction2D *uh, 
                                              CoeffFct2D *Coeff);

void LocalProjectionCoarseGridQ0(TFEFunction2D *uh, TFEFunction2D *uh_proj,
                                 CoeffFct2D *Coeff, int convection_flag);  
         
void LocalProjectionCrossWindCoarseGridQ0(TDomain *Domain, int mg_level,
                                          TFEFunction2D *uh,
                                          TFEFunction2D *uh_proj,
                                          CoeffFct2D *Coeff,
                                          double *rhs, int convection_flag); 


void AdaptivePostProcess(TFEFunction2D *FeFunction, double *PostSol, bool DirichletBC);


void AddALEStreamlineLPS(TSquareMatrix2D* A, int N_FeFunct, TFEFunction2D **FeFunct,
                         double lpcoeff, double lpexponent, int OrderDiff);
#else // __3D__ 

void AddStreamlineTerm(TSquareMatrix3D* A, TFEFunction3D *uh1,
                       TFEFunction3D *uh2, TFEFunction3D *uh3,
                       double lpcoeff, double lpexponent, int OrderDiff); 

void UltraLocalProjection(TSquareMatrix3D* A, 
                          double lpcoeff, double lpexponent, int OrderDiff);

FE3D GetElement3D(TBaseCell *cell, int CoarseOrder);

void UltraLocalProjection3D(void* A, bool ForPressure);

double UltraLocalError3D(TFEFunction3D *uh, DoubleFunct3D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff);


#endif // __3D__

#endif // __LOCAL_PROJECTION__
