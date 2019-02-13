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
// @(#)TetraAffin.h        1.3 08/18/99
//
// Class:      TTetraAffin
//
// Purpose:    reference transformations for tetrahedron
//
// Author:     Daniel Quoos
//
// History:    01.02.00 start implementation
// 
// =======================================================================

#ifndef __TETRAAFFIN__
#define __TETRAAFFIN__

#include <Enumerations.h>
#include <RefTrans3D.h>

/** reference transformations for tetrahedron */
class TTetraAffin : public TRefTrans3D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3;

    /** y coordinate */
    double y0, y1, y2, y3;

    /** z coordinate */
     double z0, z1, z2, z3;   

    /** x parameters for reference transformation */
    double xc0, xc1, xc2, xc3;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2, yc3;

    /** z parameters for reference transformation */
    double zc0, zc1, zc2, zc3;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    TTetraAffin();

    /** transfer a point from reference face to original element face */
    void GetOrigBoundFromRef(int Joint, double xi, double eta,
			     double &X, double &Y, double &Z);
    
    /** transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double zeta, double &x, double &y, double &z);

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *eta, double *xi, double *zeta, 
                        double *x, double *y, double *z, double *absdetjk);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z, double &eta, double &xi, double &zeta);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct3D BaseFunct,
                       int N_Points, double *xi, double *eta, double *zeta,
                       int N_Functs, QuadFormula3D QuadFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct3D *BaseFunct,
                       int N_Points, double *xi, double *eta, double *zeta,
                       QuadFormula3D QuadFormula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, double zeta, int N_BaseFunct,
                       double *uref, double *uxiref, double *uetaref, 
                       double *zetaref,
                       double *uorig, double *uxorig, double *uyorig, 
                       double *uzorig, int _BaseVectDim = 1);

    /** calculate functions and derivatives from reference element
        to original element on joint, parameters on joint are p1, p2 */
    void GetOrigValuesJoint(int JointNr, double p1, double p2, int N_BaseFunct,
          double *uref, double *uxiref, double *uetaref, double *uzetaref,
          double *uorig, double *uxorig, double *uyorig, double *uzorig);
	    
    // for compatibility 
    void GetOrigValues(int JointNr, double p1, double p2, int N_BaseFunct,
          double *uref, double *uxiref, double *uetaref, double *uzetaref,
          double *uorig, double *uxorig, double *uyorig, double *uzorig);
		
    /** set element to cell */
    void SetCell(TBaseCell * cell);

    /** return outer normal unit vector */
    void GetOuterNormal(int j, double s, double t,
                        double &n1, double &n2, double &n3);

    /** return two tangent vectors */
    void GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23);
    
    /** @brief Piola transformation for vector basis */
    void PiolaMapOrigFromRef(int N_Functs, double *refD000, double *origD000);
};

#endif
