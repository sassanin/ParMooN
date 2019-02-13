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
// @(#)RefTrans2D.h        1.4 04/13/00
//
// Class:      TRefTrans2D
//
// Purpose:    reference transformations for 2D geometric objects
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#ifndef __REFTRANS2D__
#define __REFTRANS2D__

#include <Constants.h>
#include <Enumerations.h>
#include <BaseCell.h>

/** reference transformations for 2D geometric objects */
class TRefTrans2D
{
  protected:
    TBaseCell *Cell;

  public:
    /** constuctor */
    TRefTrans2D() {};

    /** transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double &x, double &y);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double &eta, double &xi);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(TBaseCell *cell);

    /** set original element to cell */
    virtual void SetCell(TBaseCell *cell)
    {  Cell = cell; }

    static RefTrans2D FindRefTrans2D
        (int N_LocalUsedElements, FE2D *LocalUsedElements);

    /** return outer normal vector */
    virtual void GetOuterNormal(int j, double zeta,
                                double &n1, double &n2) = 0;

    /** return tangent */
    virtual void GetTangent(int j, double zeta,
                                double &t1, double &t2) = 0;

    /** return volume of cell according to reference transformation */
    double GetVolume();
    
    // piola map, needed for vector values basis functions such as 
    // Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
    virtual void PiolaMapOrigFromRef(int N_Functs, double *refD00, double *origD00 ) 
    { cout << " Piola Map not defined for this element " << endl; };
    
    // piola map for derivatives
    virtual void PiolaMapOrigFromRef(int N_Functs, double *refD10, 
                                     double *refD01, double *origD10, 
                                     double *origD01)
    { cout << " Piola Map not defined for this element " << endl; };
};

#endif
