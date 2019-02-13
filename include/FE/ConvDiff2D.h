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
   
// ======================================================================
// @(#)ConvDiff2D.h        1.13 10/19/99
//
// common declaration for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF2D__
#define __CONVDIFF2D__

#ifdef __2D__
/** the local assembling routines. Each of them corresponds to one 
 * LocalAssembling2D_type */

/** ========================================================================= */
/** ========================================================================= */
// CD2D: stationary convection diffusion problems

// Galerkin
void BilinearAssembleGalerkin(double Mult, double *coeff, double* param, 
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);
// SDFEM - Streamline Diffusion Finite Element Method, SUPG - Streamline Upwind 
// Petrov Galerkin
void BilinearAssemble_SD(double Mult, double *coeff, double* param, double hK,
                         double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs);
// GLS - Galerkin Least Squares
void BilinearAssemble_GLS(double Mult, double *coeff, double* param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs);

void BilinearAssemble_Axial3D(double Mult, double *coeff, double* param,
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);


/** ========================================================================= */
/** ========================================================================= */
// TCD2D: time dependent convection diffusion problems

// Galerkin
void MatrixARhsAssemble(double Mult, double *coeff, double *param, double hK, 
                        double **OrigValues, int *N_BaseFuncts, 
                        double ***LocMatrices, double **LocRhs);
void MatrixMRhsAssemble(double Mult, double *coeff, double *param, double hK, 
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// SUPG
void MatricesAKRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                                double hK, double **OrigValues,
                                int *N_BaseFuncts,double ***LocMatrices,
                                double **LocRhs);
void MatrixMRhsAssemble_SUPG(double Mult, double *coeff, double *param,
                             double hK, double **OrigValues,int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);


/** ========================================================================= */
/** ========================================================================= */
/** ========================================================================= */
/** ========================================================================= */


void BilinearAssemble_UPW1(double Mult, double *coeff, double* param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

void BilinearAssemble_UPW2(double Mult, double *coeff, double* param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD(double Mult, double *coeff, double *param,
                           double hK, double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD_Orthogonal(double Mult, double *coeff,
                                      double *param, double hK,
                                      double **OrigValues, int *N_BaseFuncts,
                                      double ***LocMatrices, double **LocRhs);

void RhsAssemble_LP96(double Mult, double *coeff, double *param, double hK,
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);

void BilinearAssemble_MH_Kno06(double Mult, double *coeff, double *param,
                               double hK,
                               double **OrigValues, int *N_BaseFuncts,
                               double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SD_SOLD(double Mult, double *coeff, double *param, 
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);
            
void BilinearAssemble2LevelLPS_Q0(double Mult, double *coeff, double *param,
                                  double hK, double **OrigValues,
                                  int *N_BaseFuncts, double ***LocMatrices,
                                  double **LocRhs);

void RhsAssemble_RhsAdjointEnergyEstimate(double Mult, double *coeff, 
                                          double *param, double hK,
                                          double **OrigValues, 
                                          int *N_BaseFuncts,
                                          double ***LocMatrices,
                                          double **LocRhs);
void RhsAssemble_RhsAdjointTV(double Mult, double *coeff, double *param,
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

void RhsAssemble_RhsAdjointTV2(double Mult, double *coeff, double *param,
                               double hK, double **OrigValues,
                               int *N_BaseFuncts, double ***LocMatrices,
                               double **LocRhs);

void RhsAssemble_RhsAdjointNormBL1_NormBorthL1(double Mult, double *coeff,
                                               double *param, double hK,
                                               double **OrigValues,
                                               int *N_BaseFuncts,
                                               double ***LocMatrices,
                                               double **LocRhs);

void RhsAssemble_RhsAdjointNormResidualL1_NormBorthL1(double Mult,
                                                      double *coeff,
                                                      double *param, double hK,
                                                      double **OrigValues,
                                                      int *N_BaseFuncts,
                                                      double ***LocMatrices,
                                                      double **LocRhs);
void RhsAssemble_RhsAdjointAll(double Mult, double *coeff, double *param,
                               double hK, double **OrigValues,
                               int *N_BaseFuncts, double ***LocMatrices,
                               double **LocRhs);


void RhsAssemble_RhsAdjointL2Error(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);

void RhsAssemble_RhsAdjointH1Error(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues, 
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);


// ========================================================================
// parameter functions
// ========================================================================



#endif // __CONVDIFF2D__

#endif // #ifdef __2D__
