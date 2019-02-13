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
// %W% %G%
// 
// Class:       TFEFunction3D
// Purpose:     a function from a finite element space in 3D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================

#ifndef __FEFUNCTION3D__
#define __FEFUNCTION3D__

#include <AllClasses.h>
#include <FESpace3D.h>
#include <AuxParam3D.h>
#include <Constants.h>
#include <vector>

/** a function from a finite element space */
class TFEFunction3D
{
  protected:
    /** name of the function */
    char *Name;

    /** some more words describing the function */
    char *Description;

    /** space to which this function belongs to */
    TFESpace3D *FESpace3D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /** constructor with vector initialization */
    TFEFunction3D(TFESpace3D *fespace3D, char *name, char *description,
                  double *values, int length);

    /** destructor */
    ~TFEFunction3D();

    /** return name */
    char *GetName()
    { return Name; }

    /** return description */
    char *GetDescription()
    { return Description; }

    /** return fe space */
    TFESpace3D *GetFESpace3D()
    { return FESpace3D; }

    /** return length */
    int GetLength()
    { return Length; }

    /** return vector of data */
    double *GetValues()
    { return Values; }

    /** calculate errors to given function 
     * NOTE: errors must be of length N_Errors+1 !!!! */
    void GetErrors(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod3D *ErrorMeth, 
                   CoeffFct3D *Coeff, TAuxParam3D *Aux,
                   int n_fespaces, TFESpace3D **fespaces,
                   double *errors);
    
    void GetErrorsForVectorValuedFunction(DoubleFunct3D ** const Exact,
                                          ErrorMethod3D * const ErrMeth,
                                          double * const errors);
    
    void GetErrorsAdapt(DoubleFunct3D *Exact, int N_Derivatives,
			MultiIndex3D *NeededDerivatives,
			int N_Errors, ErrorMethod3D *ErrorMeth, 
			CoeffFct3D *Coeff, 
			TAuxParam3D *Aux,
			int n_fespaces, TFESpace3D **fespaces,
			double *errors);
    
    /** calculate errors to given function taylored to OPTPDE */
    void GetErrorsOPTPDE(DoubleFunct3D *Exact, int N_Derivatives,
		   MultiIndex3D *NeededDerivatives,
		   int N_Errors, ErrorMethod3D *ErrorMeth, 
		   CoeffFct3D *Coeff, TAuxParam3D *Aux,
		   int n_fespaces, TFESpace3D **fespaces,
		   double radius, double upper, double lower, double *errors);
    

    /** calculate errors to given function */
    void GetMeshCellParams(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod3D *ErrorMeth, 
                   CoeffFct3D *Coeff, TAuxParam3D *Aux,
                   int n_fespaces, TFESpace3D **fespaces,
                   double *errors, double *cell_parameters);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradient(double x, double y, double z, double *values);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradientLocal(TBaseCell *cell, int cell_no, 
                           double x, double y, double z, 
                           double *values);

    /** determine the value of function
        the given point */
    void FindValueLocal(TBaseCell *cell, int cell_no, 
                           double x, double y, double z, 
                           double *values);

    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct3D *Exact);

   /** calculate the super-convergence interpolation of an exact function */
    void InterpolateSuper(DoubleFunct3D *Exact);
    
    /**  interpolation of an exact function with give FeFunction values */
    void Interpolate(int N_Coord, double *Coords, int N_AuxFeFcts, 
                     TFEFunction3D **AuxFeFcts, DoubleFunctND *Exact);
    
    /** @brief interpolate a vector valued function
     *
     * @param[in] Exact must be of length 3 (= space dimension)
     * @warning EvalAll must be correctly implemented for the used finite element
     */
    void Interpolate_vector_valued_function(std::vector<DoubleFunct3D*> Exact);
    
    /**
     * @brief project this functions into the space L20 (having zero mean value)
     *
     * After a call to this function the mean value (integral of this function
     * devided by the measure of its domain) has the value a. This is for
     * example needed for the pressure in a Stokes problem with Dirichlet
     * conditions on all boundaries.
     *
     * @param[in] a set mean value of this FEFunction3D to a
     */
    void project_into_L20(double a = 0.0);

    /**
     * @brief find the integral of this function and the measure of its domain
     *
     * @param[out] integral double value for the integral of this TFEFunction2D
     * @param[out] measure double value for the measure of its domain
     */
    void compute_integral_and_measure(double& integral, double& measure) const;

  /** @brief Set Dirichlet values according to boundary conditions */
    void SetDirichletBC(BoundCondFunct3D *BoundaryCondition,
                                   BoundValueFunct3D *BoudaryValue);
    
    /** @brief find the largest and smallest element in the vector of this FE 
     *         function */
   void MinMax(double & min, double & max);

   /** @brief print the largest and smallest element in the vector of this FE 
    *         function */
   void PrintMinMax();
};

#endif
