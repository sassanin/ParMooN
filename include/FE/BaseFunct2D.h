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
// @(#)BaseFunct2D.h        1.5 02/08/00
//
// Class:      TBaseFunct2D
//
// Purpose:    represents the set of base functions for a finite element
//             in two dimensions
//
// Author:     Gunar Matthies
//
// History:    08.07.98 restart implementation
// //     :    10.06.2010 methods for vector basis functions (Sashikumaar Ganesan)
// =======================================================================

#ifndef __BASEFUNCT2D__
#define __BASEFUNCT2D__

#include <QuadFormula1D.h>
#include <QuadFormula2D.h>
#include <Constants.h>
#include <GridCell.h>

#include <Enumerations.h>

/** set of all base function on the reference element for a finite 
    element in two dimensions */
class TBaseFunct2D
{
  protected:
    /** number of base functions = dimension of local space */
    int Dimension;

    /** Id for this set of base functions */
    BaseFunct2D BaseFunct;

    /** array for all functions and derivatives */
    DoubleFunct2D *Functions[N_MultiIndices2D];

    /** status of changability of entries */
    bool changable;

    /** reference element used for this set of base functions */
    BF2DRefElements RefElement;

    /** polynomial degree */
    int PolynomialDegree;

    /** accuracy */
    int Accuracy;

    /** number of basis functions per joint where the sign
        has to be changed if needed */
    int N_BF2Change;

    /** indices of basis functions with changeable sign,
        sorted by joints */
    int **BF2Change;    
 
    /** Dimension of the vector basis function */
    int BaseVectDim;
    
    /** Space dependent basis function, (LPS with Exp bubble) */
    bool SpaceDeptBasis;    
    
  public:
    /** constructor, fill in all information */
    TBaseFunct2D(int dimension, BaseFunct2D basefunct,
                 BF2DRefElements refelement,
                 DoubleFunct2D* functions, 
                 DoubleFunct2D* derivativesxi,
                 DoubleFunct2D* derivativeseta,
                 DoubleFunct2D* derivativesxixi,
                 DoubleFunct2D* derivativesxieta,
                 DoubleFunct2D* derivativesetaeta,
                 int polynomialdegree,
                 int accuracy,
                 int n_bf2change,
                 int **bf2change
                );

    /** constructor, fill in all information  with scalar basis function dimension*/
    TBaseFunct2D(int dimension, BaseFunct2D basefunct,
                 BF2DRefElements refelement,
                 DoubleFunct2D* functions, 
                 DoubleFunct2D* derivativesxi,
                 DoubleFunct2D* derivativeseta,
                 DoubleFunct2D* derivativesxixi,
                 DoubleFunct2D* derivativesxieta,
                 DoubleFunct2D* derivativesetaeta,
                 int polynomialdegree,
                 int accuracy,
                 int n_bf2change,
                 int **bf2change,
                 int baseVectDim
                );


    TBaseFunct2D(int dimension, BaseFunct2D basefunct,
                 BF2DRefElements refelement,
                 DoubleFunct2D* functions, 
                 DoubleFunct2D* derivativesxi,
                 DoubleFunct2D* derivativeseta,
                 DoubleFunct2D* derivativesxixi,
                 DoubleFunct2D* derivativesxieta,
                 DoubleFunct2D* derivativesetaeta,
                 int polynomialdegree,
                 int accuracy,
                 int n_bf2change,
                 int **bf2change,
                 bool spaceDeptBasis
                );



    /** constructor without filling data structure */
    TBaseFunct2D(int dimension);

    /** return the dimension of local space */
    int GetDimension() const
    { return Dimension; }

    /** return BaseFunct_ID */
    BaseFunct2D GetID() const
    { return BaseFunct; }

    /** return the values for derivative MultiIndex at (xi,eta) */
    void GetDerivatives(MultiIndex2D MultiIndex, double xi,
                        double eta, double *values)
    { Functions[MultiIndex](xi, eta, values); };

    /** return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex2D MultiIndex, 
                        TQuadFormula2D *formula, double **values);

    /** return values on joint i */
    void GetValues(int N_Points, double *zeta, int i, double **Values);

    /** return derivatives on joint i */
    void GetDerivatives(MultiIndex2D MultiIndex, int N_Points,
                        double *zeta, int i, double **Values);

   /** return values of derivative index on joint */
   void GetValues(int N_Points, double *zeta, int i, 
                  MultiIndex2D index, double **Values);

    /** set status to unchangable */
    void SetUnchangable()
      { changable = false; };

    /** set function for derivative MultiIndex */
    void SetFunction(MultiIndex2D MultiIndex, DoubleFunct2D* function);

    /** make date on reference element */
    void MakeRefElementData(QuadFormula1D QuadFormula);

    /** make date on reference element */
    void MakeRefElementData(QuadFormula2D QuadFormula);

    /** generate reference element */
    TGridCell *GenerateRefElement();

    /** return reference element */
    BF2DRefElements GetRefElement() const
      { return RefElement; };

    /** return polynomial degree */
    int GetPolynomialDegree() const
      { return PolynomialDegree; };

    /** return accuracy */
    int GetAccuracy() const
      { return Accuracy; };

    /** return number of changeable basis functions per joint */
    int GetN_BF2Change() const
      { return N_BF2Change; }

    /** return array with basis function indices */
    int **GetBF2Change() const
      { return BF2Change; }

    /** change basis functions on cell if needed */
    void ChangeBF(TCollection *Coll, TBaseCell *Cell, double *Values);

    /** change basis functions on cell in all points if needed */
    void ChangeBF(TCollection *Coll, TBaseCell *Cell, int N_Points, double **Values);
    
    /** return the dimension of the vector basis function */
    int GetBaseVectDim() const
    { return BaseVectDim; }
    

};

#endif
