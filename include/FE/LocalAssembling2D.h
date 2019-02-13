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
   
/** =======================================================================
 * @(#)LocalAssembling2D.h        10.03.2015
 * 
 * @Class:    LocalAssembling2D
 * Purpose:   Assemble on one cell a couple of bilinear and linear forms. That
 *            means the loop over all quadrature points is done within this 
 *            class.
 * 
 *            Furthermore this class includes the computation of function values
 *            at quadrature points, where the function is given by some finite
 *            element function. 
 * 
 * @note This class replaces former TDiscreteForm2D and TAuxParam2D
 * 
 * @Author:      Ulrich Wilbrandt (10.03.2015)
 * 
 * History:     start of implementation 10.03.2015 (Ulrich Wilbrandt)
 * 
 * =======================================================================
 */
#ifndef __LOCAL_ASSEMBLING_2D__
#define __LOCAL_ASSEMBLING_2D__

#include <Enumerations.h>
#include <Constants.h>
#include <string>
#include <vector>

enum LocalAssembling2D_type { CD2D_Galerkin,
                              CD2D_SUPG,
                              CD2D_GLS,
                              CD2D_Axiax3D_Galerkin,
                              TCD2D_Mass_Rhs_Galerkin, // mass matrix and rhs
                              TCD2D_Stiff_Rhs_Galerkin,// stiffness matrix + rhs
                              TCD2D_Mass_Rhs_SUPG, // mass matrix and rhs
                              TCD2D_Stiff_Rhs_SUPG, // stiffness matrix + rhs
                              NSE2D_Galerkin,
                              NSE2D_Galerkin_Nonlinear,
                              Darcy2D_Galerkin
};

/** a function from a finite element space */
class LocalAssembling2D
{
  protected:
    /** name */
    std::string name;

    /** @brief number of terms */
    int N_Terms;

    /** @brief number of involved spaces (typically one or two) */
    int N_Spaces;

    /** @brief for each space we store a bool indicatin if second derivatives 
     *         are needed */
    bool *Needs2ndDerivatives;

    /** @brief multiindices for derivatives of ansatz and test functions 
     * 
     * This is an array of size N_Terms.
     */
    std::vector<MultiIndex2D> Derivatives;

    /** @brief for each term, there is one FESpace2D asociated with that term */
    std::vector<int> FESpaceNumber;

    /** @brief which FE space corresponds to each row */
    std::vector<int> RowSpace;

    /** @brief which FE space corresponds to each column */
    std::vector<int> ColumnSpace;

    /** @brief which FE space corresponds to each right-hand side */
    std::vector<int> RhsSpace;

    /** function for calculating the coefficients */
    CoeffFct2D *Coeffs;

    /** @brief function doing the real assembling using parameters from 
     *         argument list */
    AssembleFctParam2D *AssembleParam;

    /** function for manipulating the coefficients */
    ManipulateFct2D *Manipulate;

    /** memory for storing the original value arrays */
    double ***AllOrigValues;

    /** memory for storing the original value arrays at one point */
    double **OrigValues;
    
    int N_Matrices;
    int N_Rhs;
    
    
    /** number of stored parameter functions (ParamFct) */
    int N_ParamFct;
    
    /** array of stored parameter function */
    std::vector<ParamFct*> ParameterFct;
    
    /** index of first parameter produced by parameter function i */
    std::vector<int> BeginParameter;
    
    // number of parameters
    int N_Parameters;
    
    /** array of stored FEFunction2D */
    TFEFunction2D **FEFunctions2D;
    
    /** number of FE values */
    int N_FEValues;
    
    /** index of FEFunction2D used for FE value i */
    std::vector<int> FEValue_FctIndex;
    
    /** which multiindex is used for FE value i */
    std::vector<MultiIndex2D> FEValue_MultiIndex;

    /** Depending on the NSTYPE and the NSE_NONLINEAR_FORM all parameters are 
     * set within this function. This function is called from the constructor 
     * in case of Navier-Stokes problems. It only exists in order to not make 
     * the constructor huge. 
     * 
     * Basically this function implements three nested switches.
     */
    void set_parameters_for_nse(LocalAssembling2D_type type);
  public:
    /** constructor */
    LocalAssembling2D(LocalAssembling2D_type type, TFEFunction2D **fefunctions2d,
                      CoeffFct2D *coeffs);

    /** destructor */
    ~LocalAssembling2D();
    
    

    /** return local stiffness matrix */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                       double *X, double *Y,
                       int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                       double **Parameters, double **AuxArray,
                       TBaseCell *Cell, int N_Matrices, int N_Rhs,
                       double ***LocMatrix, double **LocRhs,
                       double factor = 1.);
    
    /** return all parameters at all quadrature points */
    void GetParameters(int n_points, TCollection *Coll,
                       TBaseCell *cell, int cellnum,
                       double *x, double *y,
                       double **Parameters);

    /** return all parameters at boundary points */
    void GetParameters(int N_Points, TCollection *Coll,
                       TBaseCell *cell, int cellnum,
                       double *s, int joint,
                       double **Parameters);
    
    

    /** return name */
    const std::string& get_name() const
    { return name; }

    /** return array Needs2ndDerivatives */
    bool *GetNeeds2ndDerivatives() const
    { return Needs2ndDerivatives; }

    /** function for calculating the coefficients */
    CoeffFct2D *GetCoeffFct() const
    { return Coeffs; }
    
    /** return the index of the row space of the i-th matrix */
    int rowSpaceOfMat(int i) const
    { return RowSpace[i]; }
    
    /** return the index of the column space of the i-th matrix */
    int colSpaceOfMat(int i) const
    { return ColumnSpace[i]; }
    
    int GetN_ParamFct() const
    { return N_ParamFct; }
    
    int GetN_Parameters() const
    { return N_Parameters; }
};


#endif // __LOCAL_ASSEMBLING_2D__
