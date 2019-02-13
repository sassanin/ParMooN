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
// @(#)TriaIsoparametric.h        1.5 04/13/00
//
// Class:      TTriaIsoparametric
//
// Purpose:    isoparametric reference transformations for triangle
//
// Author:     Gunar Matthies
//
// History:    29.04.98 start implementation
// 
// =======================================================================

#ifndef __TRIAISOPARAMETRIC__
#define __TRIAISOPARAMETRIC__

#include <Enumerations.h>
#include <RefTrans2D.h>

/** reference transformations for triangle */
class TTriaIsoparametric : public TRefTrans2D
{
  protected:
    /** x coordinate */
    double x[3];

    /** y coordinate */
    double y[3];

    /** x parameters for reference transformation */
    double xc0, xc1, xc2;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2;

    /** number of additional points */
    int N_AuxPoints;

    /** distance in x direction between real auxiliary point and 
        its position after a affine mapping */
    double XDistance[MaxN_BaseFunctions2D];

    /** distance in y direction between real auxiliary point and 
        its position after a affine mapping */
    double YDistance[MaxN_BaseFunctions2D];

    /** order of approximation */
    int ApproximationOrder;

    /** values of corresponding base function at quadpoints */
    double FctValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** xi-derivatives of corresponding base function at quadpoints */
    double XiDerValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** eta-derivatives of corresponding base function at quadpoints */
    double EtaDerValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** base function type for each order of approximation */
    static BaseFunct2D BaseFunctFromOrder[]; 

    /** base function type for each order of approximation */
    static FEDesc2D FEDescFromOrder[]; 

    /** auxiliary array */
    double DoubleAux[MaxN_BaseFunctions2D];

    /** auxiliary array */
    int IntAux[MaxN_BaseFunctions2D];

    /** used quadrature rule */
    QuadFormula2D QuadFormula;
  
    /** for data from quadrature formula */
    double *XI, *ETA, *W;

    /** number of quadrature points */
    int N_QuadPoints;

  public:
    /** constuctor */
    TTriaIsoparametric();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double &x, double &y);

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *eta, double *xi, 
                        double *x, double *y, double *absdetjk);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double &eta, double &xi);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct2D BaseFunct,
                       int N_Points, double *xi, double *eta,
                       int N_Functs, QuadFormula2D formula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct2D *BaseFunct,
                       int N_Points, double *xi, double *eta,
                       QuadFormula2D formula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, int N_BaseFunct,
                double *uref, double *uxiref, double *uetaref,
                double *uorig, double *uxorig, double *uyorig);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                double *uref, double *uxiref, double *uetaref,
                double *uorig, double *uxorig, double *uyorig);

    /** set element to cell */
    void SetCell(TBaseCell * cell);

    /** set order of approximation */
    void SetApproximationOrder(int order)
    {
      if(order <= 0)
        ApproximationOrder = 1;
       else
        ApproximationOrder = order;
    }

    /** set used quadrature formula */
    void SetQuadFormula(QuadFormula2D formula)
    { QuadFormula = formula; }

    /** return outer normal vector */
    void GetOuterNormal(int j, double zeta,
                                double &n1, double &n2);

    /** return tangent */
    void GetTangent(int j, double zeta,
                                double &t1, double &t2);

    /** return volume of cell according to reference transformation */
    double GetVolume();

    /** return boundary vertices */
    void GetOrigBoundFromRef( int joint, int N_Points, double *zeta, double *X, double *Y);

};

#endif
