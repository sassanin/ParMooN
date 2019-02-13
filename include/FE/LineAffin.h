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
// @(#)LineAffin.h
//
// Class:      LineAffin
//
// Purpose:    reference transformations for Line
//
// Author:     Sashikumaar Ganesan
//
// History:    17.05.2007 start implementation
// 
// =======================================================================

#ifndef __LINEAFFIN__
#define __LINEAFFIN__

#include <Enumerations.h>
#include <RefTrans1D.h>

/** reference transformations for line */
class TLineAffin : public TRefTrans1D
{
  protected:
    /** x coordinate */
    double x0, x1;

    /** x coordinate (useful in surface finite element) */
    double y0, y1;

    /** x parameters for reference transformation */
    double xc0, xc1;

    /** x parameters for reference transformation */
    double yc0, yc1;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    TLineAffin();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double xi, double &X
#ifdef __2D__
, double &Y
#endif
                                );

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *xi, double *X, double *Y, double *absdetjk);
    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, double *xi, double *X, double *absdetjk);

    void GetOrigValues(int N_Sets, BaseFunct1D *BaseFuncts, int N_Points, double *zeta,
                       QuadFormula1D QuadFormula, bool *Needs2ndDer);

//     /** set element to cell */
    void SetCell(TBaseCell * cell);
    
    double Getrec_detjk()
    { return rec_detjk; }
};

#endif
 
