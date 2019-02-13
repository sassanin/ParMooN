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
// @(#)FEFunction2D.h        1.3 04/13/00
// 
// Class:       TFEFunction2D
// Purpose:     a function from a finite element space in 2D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================

#ifndef __FEFUNCTION2D__
#define __FEFUNCTION2D__

#include <AllClasses.h>
#include <FESpace2D.h>
#include <AuxParam2D.h>
#include <Constants.h>

/** a function from a finite element space */
class TFEFunction2D
{
  protected:
    /** name of the function */
    char *Name;

    /** some more words describing the function */
    char *Description;

    /** space to which this function belongs to */
    TFESpace2D *FESpace2D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /** constructor with vector initialization */
    TFEFunction2D(TFESpace2D *fespace2D, char *name, char *description,
                  double *values, int length);

    /** destructor */
    ~TFEFunction2D();

    /** return name */
    char *GetName()
    { return Name; }

    /** return description */
    char *GetDescription()
    { return Description; }

    /** return fe space */
    TFESpace2D *GetFESpace2D()
    { return FESpace2D; }

    /** return length */
    int GetLength()
    { return Length; }

    /** return vector of data */
    double *GetValues()
    { return Values; }

    /** calculate errors to given function */
    void GetErrors(DoubleFunct2D *Exact, int N_Derivatives,
                   MultiIndex2D *NeededDerivatives,
                   int N_Errors, ErrorMethod2D *ErrorMeth, 
                   CoeffFct2D *Coeff, TAuxParam2D *Aux,
                   int n_fespaces, TFESpace2D **fespaces,
                   double *errors);

    /** @brief use this for vector valued basis functions (Raviart-Thomas (RT)
     *         or Brezzi-Douglas-Marini (BDM) elements) */
    void  GetErrorsForVectorValuedFunction(
                  DoubleFunct2D * const * const Exact, 
                  ErrorMethod2D * const ErrorMeth, 
                  double * const errors);
    
    void GetErrorsAdapt(DoubleFunct2D *Exact, int N_Derivatives,
		   MultiIndex2D *NeededDerivatives,
		   int N_Errors, ErrorMethod2D *ErrorMeth, 
		   CoeffFct2D *Coeff, TAuxParam2D *Aux,
		   int n_fespaces, TFESpace2D **fespaces,
		   double *errors);
    
    void GetErrorsOPTPDE(DoubleFunct2D *Exact, int N_Derivatives,
		   MultiIndex2D *NeededDerivatives,
		   int N_Errors, ErrorMethod2D *ErrorMeth, 
		   CoeffFct2D *Coeff, TAuxParam2D *Aux,
		   int n_fespaces, TFESpace2D **fespaces,
		   int& kink, double upper, double lower, double *errors);
    
    void GetErrorsAdaptOPTPDE(DoubleFunct2D *Exact, int N_Derivatives,
			MultiIndex2D *NeededDerivatives,
			int N_Errors, ErrorMethod2D *ErrorMeth, 
			CoeffFct2D *Coeff, TAuxParam2D *Aux,
			int n_fespaces, TFESpace2D **fespaces,
			double radius, double upper, double lower,double *errors);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradient(double x, double y, double *values);

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradientLocal(TBaseCell *cell, int cell_no, double x, double y, double *values);

    /** determine the value of function at
        the given point */
    void FindValueLocal(TBaseCell *cell, int cell_no, double x, double y, double *values);

    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct2D *Exact);
    /** interpolate the old mesh fe function values to the new fe function
    * 
    * Note that this is rather slow, because no further information is 
    * required. The function 'OldFeFunction' could even live on a larger domain.
    */
    void Interpolate(TFEFunction2D *F);
    
    /**
     * @brief project this functions into the space L20 (having zero mean value)
     * 
     * After a call to this function the mean value (integral of this function
     * devided by the measure of its domain) has the value a. This is for 
     * example needed for the pressure in a Stokes problem with Dirichlet 
     * conditions on all boundaries.
     * 
     * @param a set mean value of this FEFunction2D to a
     */
    void project_into_L20(double a = 0.0);
    
    /**
     * @brief find the integral of this function and the measure of its domain
     * 
     * @param integral double value for the integral of this TFEFunction2D
     * @param measure double value for the measure of its domain 
     */
    void compute_integral_and_measure(double& integral, double& measure) const;
    
    /**
     * @brief compute the mean value of this TFEFunction2D
     * 
     * this functions uses 'compute_integral_and_measure'. Then the mean is the
     * integral divided by the measure.
     */
    double compute_mean() const;    

    /** calculate parameters which are connected to a mesh cell */
    void GetMeshCellParams(DoubleFunct2D *Exact, int N_Derivatives,
                           MultiIndex2D *NeededDerivatives,
                           int N_Errors, ErrorMethod2D *ErrorMeth, 
                           CoeffFct2D *Coeff, 
                           TAuxParam2D *Aux,
                           int n_fespaces, TFESpace2D **fespaces,
                           double *errors, double *parameters);

    /** calculate the super-convergence interpolation of an exact function */
    void InterpolateSuper(DoubleFunct2D *Exact);
    
    /** set Dirichlet values according to boundary conditions */
    void SetDirichletBC(BoundCondFunct2D *BoundaryCondition,
                        BoundValueFunct2D *BoundaryValue);
    
   /** write the solution into a data file - written by Sashi **/
   void WriteSol();

   /** Read the solution from a given data file - written by Sashi **/
   void ReadSol(char *BaseName);
   

   /** sol will be correct to conserve the Old_Mass (remessing, temp, surfact, psd, etc) - added by sashi */   
   void CorrectMass(double OldMass);

   /** Retun the mass, domain volume and mean values of the function - added by sashi */
   void GetMassAndMean(double *OutVal);
   
   
   /** multiply function with a scalar alpha. Only non-Dirichlet dofs are 
       multiplied! */
   TFEFunction2D& operator*=(double alpha);
   
   /** add one TFEFunction2D to another one. Both have to be defined on the 
       same space. Only non-Dirichlet dofs are added!  */
   TFEFunction2D & operator+=(const TFEFunction2D & rhs);
   
   /** copy one TFEFunction2D to another one. Both have to be defined on the 
       same space */
   TFEFunction2D & operator=(const TFEFunction2D & rhs);
   /** find the largest and smallest element in the vector of this FE function
   */
   void MinMax(double & min, double & max);

   /** print the largest and smallest element in the vector of this FE 
       function 
   */
   void PrintMinMax();

};
#endif
