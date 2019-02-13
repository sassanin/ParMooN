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
// @(#)RefTrans1D.h        
//
// Class:      TRefTrans1D
//
// Purpose:    reference transformations for 1D geometric objects
//
// Author:     Sashikumaar Ganesan
//
// History:    17.05.2007 start implementation
// 
// =======================================================================

#ifndef __REFTRANS1D__
#define __REFTRANS1D__

#include <Constants.h>
#include <Enumerations.h>
#include <BaseCell.h>

/** reference transformations for 1D geometric objects */
class TRefTrans1D
{
  protected:
    TBaseCell *Cell;

  public:
    /** constuctor */
    TRefTrans1D() {};

    /** transfer form reference element to original element */
    void GetOrigFromRef(double xi, double &x);


    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double &eta);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(TBaseCell *cell);

    /** set original element to cell */
    virtual void SetCell(TBaseCell *cell)
    {  Cell = cell; }

    static RefTrans1D FindRefTrans1D
        (int N_LocalUsedElements, FE1D *LocalUsedElements);


    /** return volume of cell according to reference transformation */
    double GetVolume();
};

#endif
