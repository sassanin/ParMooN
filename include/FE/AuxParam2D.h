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
// @(#)AuxParam2D.h        1.1 10/30/98
// 
// Class:       TAuxParam2D
// Purpose:     store parameter functions and FE functions
//
// Author:      Gunar Matthies (06.08.98)
//
// History:     start of implementation 06.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __AUXPARAM2D__
#define __AUXPARAM2D__

#include <Constants.h>
#include <FESpace2D.h>
#include <FEFunction2D.h>
#include <string>

/** store parameter functions and FE functions */
class TAuxParam2D
{
  public:
// =======================================================================
//  numbers of stored objects
// =======================================================================
    /** number of stored FESpace2D */
    int N_FESpace2D;

    /** number of stored FEFunction2D */
    int N_FEFunction2D;

    /** number of stored parameter function (ParamFct) */
    int N_ParamFct;

// =======================================================================
//  array of pointers to stored objects
// =======================================================================
    /** array of stored FESpace2D */
    TFESpace2D **FESpaces2D;

    /** array of stored FEFunction2D */
    TFEFunction2D **FEFunctions2D;

    /** array of stored parameter function */
    ParamFct **ParameterFct;

// =======================================================================
//  information of FE values used by parameter functions
// =======================================================================
    /** number of FE values */
    int N_FEValues;

    /** index of FEFunction2D used for FE value i */
    int *FEValue_FctIndex;

    /** which multiindex is used for FE value i */
    MultiIndex2D *FEValue_MultiIndex;

// =======================================================================
//  information of parameter functions
// =======================================================================
    /** number of all parameters */
    int N_Parameters;

    /** index of first parameter produced by parameter function i */
    int *BeginParameter;

// =======================================================================
//  information of parameter functions
// =======================================================================
    /** storage for temporary FE values */
    double *Temp;

    double **Values;
    double ***OrigValues;
    int **Index;
    int *N_BaseFunct;

  public:
    /** constructor */
    TAuxParam2D(int n_fespace2d, int n_fefunction2d, int n_paramfct,
                int n_fevalues,
                TFESpace2D **fespaces2d, TFEFunction2D **fefunctions2d,
                ParamFct **parameterfct,
                int *fevalue_fctindex, MultiIndex2D *fevalue_multiindex,
                int n_parameters, int *beginparameter);

    /** @brief standard constructor
     * 
     * If you don't need values of a finite element function in your assembling,
     * choose this constructor. This is equivalent to calling 
     * TAuxParam2D(0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
     */
    TAuxParam2D();
    /** @brief constructor used if finite element function values are needed 
     *         during assembling
     * 
     * Depending on the given parameter 'name', this object will be initialized
     * properly. Currently supported values for 'name' are:
     *  - "velocity", this is used for Navier-Stokes problems
     * 
     * You can achieve the same behavior using the first constructor above, but 
     * this constructor is easier.
     */
    TAuxParam2D(std::string name, TFEFunction2D **fefunctions2d);

    
    /** destructor */
    ~TAuxParam2D();

    /** return all parameters at all quadrature points */
    void GetParameters(int n_points, TCollection *Coll,
                       TBaseCell *cell, int cellnum,
                       double *xi, double *eta,
                       double *x, double *y,
                       double **Parameters);

    /** return all parameters at boundary points */
    void GetParameters(int N_Points, TCollection *Coll,
                       TBaseCell *cell, int cellnum,
                       double *s, int joint,
                       double **Parameters);

    int GetN_Parameters() const
    { return N_Parameters; }
    int GetN_ParamFct() const
    { return N_ParamFct; }

};

// standard function to use for Navier-Stokes
void Velocity_Fct(double *inputList, double *outputValues);


#endif
