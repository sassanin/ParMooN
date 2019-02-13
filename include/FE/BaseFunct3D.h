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
// Class:      TBaseFunct3D
//
// Purpose:    represents the set of base functions for a finite element
//             in three dimensions
//
// Author:     Gunar Matthies
//
// History:    start of implementation 19.11.99
// 
// =======================================================================

#ifndef __BASEFUNCT3D__
#define __BASEFUNCT3D__

#include <QuadFormula2D.h>
#include <QuadFormula3D.h>
#include <Constants.h>
#include <GridCell.h>

#include <Enumerations.h>

/** set of all base function on the reference element for a finite 
    element in two dimensions */
class TBaseFunct3D
{
  protected:
    /** number of base functions = dimension of local space */
    int Dimension;

    /** Id for this set of base functions */
    BaseFunct3D BaseFunct;

    /** array for all functions and derivatives */
    DoubleFunct3D *Functions[N_MultiIndices3D];

    /** status of changability of entries */
    bool changable;

    /** reference element used for this set of base functions */
    BF3DRefElements RefElement;

    /** polynomial degree */
    int PolynomialDegree;

    /** accuracy */
    int Accuracy;

    /** number of basis functions per joint where the sign
        has to be changed if needed */
    int N_BF2Change;

    /** indices of basis functions with changeable sign,
        sorted by joints */
    int ***BF2Change;
    
    /** Dimension of the vector basis function */
    int BaseVectDim;

  public:
    /** constructor, fill in all information */
    TBaseFunct3D(int dimension, BaseFunct3D basefunct,
                 BF3DRefElements refelement,
                 DoubleFunct3D* functions, 
                 DoubleFunct3D* derivativesXi,
                 DoubleFunct3D* derivativesEta,
                 DoubleFunct3D* derivativesZeta,
                 DoubleFunct3D* derivativesXiXi,
                 DoubleFunct3D* derivativesXiEta,
                 DoubleFunct3D* derivativesXiZeta,
                 DoubleFunct3D* derivativesEtaEta,
                 DoubleFunct3D* derivativesEtaZeta,
                 DoubleFunct3D* derivativesZetaZeta,
                 int polynomialdegree,
                 int accuracy,
                 int n_bf2change,
                 int ***bf2change,
                 int baseVectDim = 1
                );

    /** constructor without filling data structure */
    TBaseFunct3D(int dimension);

    /** return the dimension of local space */
    int GetDimension() const
    { return Dimension; }
    
    /** return the values for derivative MultiIndex at (xi,eta) */
    void GetDerivatives(MultiIndex3D MultiIndex, double xi,
                        double eta, double zeta, double *values)
    { Functions[MultiIndex](xi, eta, zeta, values); };

    /** return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex3D MultiIndex, 
                        TQuadFormula3D *formula, double **values);

    /** return values on joint i */
    void GetValues(int N_Points, double *t, double *s, 
                   int i, double **Values);

   /** return values of derivative index on joint */
   void GetValues(int N_Points, double *t, double *s, int i, 
                  MultiIndex3D index, double **Values);

    /** set status to unchangable */
    void SetUnchangable()
      { changable = false; };

    /** set function for derivative MultiIndex */
    void SetFunction(MultiIndex3D MultiIndex, DoubleFunct3D* function);

    /** make date on reference element */
    void MakeRefElementData(QuadFormula2D QuadFormula);

    /** make date on reference element */
    void MakeRefElementData(QuadFormula3D QuadFormula);

    /** generate reference element */
    TGridCell *GenerateRefElement();

    /** return reference element */
    BF3DRefElements GetRefElement() const
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
    int ***GetBF2Change() const
      { return BF2Change; }

    /** change basis functions on cell if needed */
    void ChangeBF(TCollection *Coll, TBaseCell *Cell, double *Values);

    /** change basis functions on cell in all points if needed */
    void ChangeBF(TCollection *Coll, TBaseCell *Cell, int N_Points, double **Values);

    /** return BaseFunct_ID */
    BaseFunct3D GetID() const
    { return BaseFunct; }
    
    /** return the dimension of the vector basis function */
    int GetBaseVectDim() const
    { return BaseVectDim; }
    
};

#endif
