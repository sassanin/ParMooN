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
// @(#)QuadAffin.h        1.5 04/13/00
//
// Class:      TQuadAffin
//
// Purpose:    affin reference transformations for parallelogram
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#ifndef __QUADAFFIN__
#define __QUADAFFIN__

#include <Enumerations.h>
#include <RefTrans2D.h>

/** reference transformations for triangle */
class TQuadAffin : public TRefTrans2D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3;

    /** y coordinate */
    double y0, y1, y2, y3;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    TQuadAffin();

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
                       int N_Functs, QuadFormula2D QuadFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct2D *BaseFunct,
                       int N_Points, double *xi, double *eta,
                       QuadFormula2D QuadFormula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, int N_BaseFunct,
                       double *uref, double *uxiref, double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1);
    
    void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                       double *uref, double *uxiref, double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1);

    /** set element to cell */
    void SetCell(TBaseCell * cell);

    /** return outer normal vector */
    void GetOuterNormal(int j, double zeta,
                                double &n1, double &n2);

    /** return tangent */
    void GetTangent(int j, double zeta,
                                double &t1, double &t2);

    /** return volume of cell according to reference transformation */
    double GetVolume();
    
    
    /** transfer a  set of boundary points from reference to original element */
    void GetOrigBoundFromRef(int joint, int N_Points, double *zeta, double *X, double *Y);   
    
    /** @brief Piola transformation for vector valued basis functions */
    void PiolaMapOrigFromRef(int N_Functs, double *refD00, double *origD00);
    /** @brief Piola transformation for the derivatives of vector valued basis 
     *         functions */
    void PiolaMapOrigFromRef(int N_Functs, double *refD10, double *refD01, 
                             double *origD10, double *origD01);
};

#endif // __QUADAFFIN__
