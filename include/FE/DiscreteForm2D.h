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
// @(#)DiscreteForm2D.h        1.6 10/18/99
// 
// Class:       TDiscreteForm2D
// Purpose:     assemble on one cell a couple of bilinear and linear forms
//
// Author:      Gunar Matthies (02.10.98)
//
// History:     start of implementation 02.10.98 (Gunar Matthies)
//              add new data 14.04.99 (Gunar Matthies)
//
// =======================================================================

#ifndef __DISCRETEFORM2D__
#define __DISCRETEFORM2D__

#include <Enumerations.h>
#include <Constants.h>

/** a function from a finite element space */
class TDiscreteForm2D
{
  protected:
    /** name */
    char *Name;

    /** some more describing words */
    char *Description;

    /** number of terms */
    int N_Terms;

    /** number of involved spaces */
    int N_Spaces;

    /** Are second derivatives from space i needed */
    bool *Needs2ndDerivatives;

    /** multiindices for derivatives of ansatz functions */
    MultiIndex2D *Derivatives;

    /** number of FESpace2D which is used for a derivative */
    int *FESpaceNumber;

    /** number of matrices */
    int N_Matrices;

    /** number of right-hand sides */
    int N_Rhs;

    /** which FE space corresponds to each row */
    int *RowSpace;

    /** which FE space corresponds to each column */
    int *ColumnSpace;

    /** which FE space corresponds to each right-hand side */
    int *RhsSpace;

    /** function for calculating the coefficients */
    CoeffFct2D *Coeffs;

    /** function doing the real assembling */
    AssembleFct2D *Assemble;

    /** function doing the real assembling using parameters from 
        argument list */
    AssembleFctParam2D *AssembleParam;

    /** function for manipulating the coefficients */
    ManipulateFct2D *Manipulate;

    /** memory for storing the original value arrays */
    double ***AllOrigValues;

    /** memory for storing the original value arrays at one point */
    double **OrigValues;

  public:
    /** constructor */
    TDiscreteForm2D(char *name, char *description, int n_terms,
        MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFct2D *assemble, CoeffFct2D *coeffs,
        ManipulateFct2D *manipulate);

    /** constructor with assembling using parameters */
    TDiscreteForm2D(char *name, char *description, int n_terms,
        MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFctParam2D *assembleparam, CoeffFct2D *coeffs,
        ManipulateFct2D *manipulate);

    /** destructor */
    ~TDiscreteForm2D();

    /** return name */
    char *GetName()
    { return Name; }

    /** return description */
    char *GetDescription()
    { return Description; }

    /** return local stiffness matrix */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                       double hK, double *X, double *Y,
                       int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                       double **Parameters, double **AuxArray,
                       TBaseCell *Cell, int N_Matrices, int N_Rhs,
                       double ***LocMatrix, double **LocRhs,
                       double factor = 1.);
    
    /** assemble local matrices and right hand sides 
     * 
     * This is a simplified version of the above GetLocalForms(...).
     */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                        double hK, double *X, double *Y,
                        int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                        TBaseCell *Cell,double ***LocMatrix, double **LocRhs);

    /** return array Needs2ndDerivatives */
    bool *GetNeeds2ndDerivatives() const
    { return Needs2ndDerivatives; };

    /** function for calculating the coefficients */
    CoeffFct2D *GetCoeffFct() const
    { return Coeffs; }
    
    /** return the index of the row space of the i-th matrix */
    int rowSpaceOfMat(int i) const
    { return RowSpace[i]; }
    
    /** return the index of the column space of the i-th matrix */
    int colSpaceOfMat(int i) const
    { return ColumnSpace[i]; }
};

/******************************************************************************/
//
// Routines for initializing discrete forms
//
/******************************************************************************/

#ifdef __2D__
void InitializeDiscreteFormsScalar(TDiscreteForm2D *&DiscreteFormMatrixMRhs, TDiscreteForm2D *&DiscreteFormMatrixARhs, TDiscreteForm2D *&DiscreteFormMatrixMRhs_SUPG, TDiscreteForm2D *&DiscreteFormMatrixARhs_SUPG, CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormSDFEM,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLSDFEM,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormPressSep,
  TDiscreteForm2D *&DiscreteFormPressSepAuxProb,
  TDiscreteForm2D *&DiscreteFormNSRFBRhs,
  TDiscreteForm2D *&DiscreteFormNSCSTRhs,
  CoeffFct2D *LinCoeffs, int NSTYPE);

void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLGalerkinDuese,
  CoeffFct2D *LinCoeffs);

void InitializeDiscreteFormsFriction(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormSDFEM,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLSDFEM,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormGalerkinFriction, 
  TDiscreteForm2D *&DiscreteFormNLGalerkinFriction,
  TDiscreteForm2D *&DiscreteFormGalerkinFrictionLocal, 
  TDiscreteForm2D *&DiscreteFormNLGalerkinFrictionLocal,
  TDiscreteForm2D *&DiscreteFormSDFEMFriction,
  TDiscreteForm2D *&DiscreteFormNLSDFEMFriction,
  TDiscreteForm2D *&DiscreteFormSDFEMDDFriction,
  TDiscreteForm2D *&DiscreteFormNLSDFEMDDFriction,
  TDiscreteForm2D *&DiscreteFormSDFEMDDFrictionRhs,
  TDiscreteForm2D *&DiscreteFormNLSDFEMDDFrictionRhs,
  TDiscreteForm2D *&DiscreteFormSDFEMDDFrictionRST,
  TDiscreteForm2D *&DiscreteFormNLSDFEMDDFrictionRST,
  TDiscreteForm2D *&DiscreteFormUpwindFriction,
  TDiscreteForm2D *&DiscreteFormNLUpwindFriction,
  CoeffFct2D *LinCoeffs, int NSTYPE);

void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormColetti,
  TDiscreteForm2D *&DiscreteFormGL00Convolution, 
  TDiscreteForm2D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLColetti, 
  TDiscreteForm2D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm2D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHSColetti,
  TDiscreteForm2D *&DiscreteFormRHSLESModel,
  TDiscreteForm2D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm2D *&DiscreteFormRHSSmagorinskyExpl,
  TDiscreteForm2D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm2D *&DiscreteFormRHSAuxProblemU,
  CoeffFct2D *LinCoeffs, int NSTYPE);
  
void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormColetti,
  TDiscreteForm2D *&DiscreteFormGL00Convolution, 
  TDiscreteForm2D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLColetti, 
  TDiscreteForm2D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm2D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHSColetti,
  TDiscreteForm2D *&DiscreteFormRHSLESModel,
  TDiscreteForm2D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm2D *&DiscreteFormRHSSmagorinskyExpl,
  TDiscreteForm2D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm2D *&DiscreteFormRHSAuxProblemU,
  TDiscreteForm2D *&DiscreteFormC,
  TDiscreteForm2D *&DiscreteFormJ,
  CoeffFct2D *LinCoeffs, int NSTYPE);

void InitializeDiscreteFormsPaper2(  
  TDiscreteForm2D *&DiscreteFormRHSGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHSPaper2,
   CoeffFct2D *LinCoeffs, CoeffFct2D *Coeffs);

void InitializeDiscreteFormsVMS(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind, 
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteForm_ho_RHS,
  TDiscreteForm2D *&DiscreteForm_ls_RHS,
  CoeffFct2D *LinCoeffs, int NSTYPE);

void InitializeDiscreteForms_SSMUM(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHS1,
  CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms_NSECST(TDiscreteForm2D *&DiscreteFormGalerkinSUPG, TDiscreteForm2D *&DiscreteFormLPS, CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms_TNSECST(TDiscreteForm2D *&DiscreteFormGalerkinSUPG, TDiscreteForm2D *&DiscreteFormLPS, TDiscreteForm2D *&DiscreteFormNLGalerkinSUPG, TDiscreteForm2D *&DiscreteFormNLLPS, TDiscreteForm2D *&DiscreteFormRHSGalerkinSUPG, TDiscreteForm2D *&DiscreteFormRHSLPS, CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms_CST(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormSDFEM, CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms_CST_Giesekus(TDiscreteForm2D *&DiscreteFormGalerkin, CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms_DFT(TDiscreteForm2D *&DiscreteFormGalerkin, CoeffFct2D *LinCoeffs);

void InitializeDiscreteFormsCDAdapt2D(TDiscreteForm2D **DiscreteForms,
              CoeffFct2D *BilinearCoeffs);

void  InitializeDiscreteForms_Moving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormNLGalerkin,
                                  TDiscreteForm2D *&DiscreteFormGrid, CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs);

void InitializeDiscreteForms_ScalarMoving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormGrid,
                                    TDiscreteForm2D *&DiscreteFormMatrixMRhs_SUPG, TDiscreteForm2D *&DiscreteFormSUPG, CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs);

void InitializeDiscreteForms_Moving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormSUPG, CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms_HeatLine(TDiscreteForm2D *&DiscreteFormHeatLine, CoeffFct2D *LinCoeffs);

void InitializeDiscreteForms_Stationary(TDiscreteForm2D *&DiscreteFormUpwind, TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormSDFEM,
                                        TDiscreteForm2D *&DiscreteFormGLS, CoeffFct2D *LinCoeffs);

void InitializeDiscreteFormGrid(TDiscreteForm2D *&DiscreteFormGrid, CoeffFct2D *GridCoeffs);

void InitializeDiscreteForms_2PhaseAxial3D(
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormGrid,
  CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs, int NSTYPE);

void InitializeDiscreteForms_2PhaseAxial3D_TNSECST(TDiscreteForm2D *&DiscreteFormLPS, 
				    TDiscreteForm2D *&DiscreteFormNLLPS, 
				     TDiscreteForm2D *&DiscreteFormRHSLPS,
				     TDiscreteForm2D *&DiscreteFormGrid,
				     CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs);

void InitializeDiscreteForms_ImpDropAxial3D_TNSECST(TDiscreteForm2D *&DiscreteFormLPS, 
				    TDiscreteForm2D *&DiscreteFormNLLPS, 
				     TDiscreteForm2D *&DiscreteFormRHSLPS,
				     TDiscreteForm2D *&DiscreteFormGrid,
				     CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs);



#endif

#endif
