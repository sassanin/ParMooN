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
// @(#)It_Mortar.h        1.2 02/08/99
// 
// Class:       TIt_Mortar
// Purpose:     iterator which gets all cells on a given mortar face
//
// Author:      Volker Behns  20.03.98
//
// =======================================================================

#ifndef __IT_MORTAR__
#define __IT_MORTAR__

#include <It_Search.h>

/** iterator which gets all cells on a given mortar face */
class TIt_Mortar : public TIt_Search
{
  protected:
    /** local number of mortar face on each level */
    int LocMortarFace[MAX_ItLevel];

    /** indicator for mortar side */
    bool MortarSide;

  public:
    // Constructors

    // Methods
    /** initialize the iterator */
    virtual int Init(int level);

    /** return the maximum level */
    virtual int GetMaxLevel();

    /** return the next cell */
    virtual TBaseCell *Next(int &info);

    /** return the previous cell */
    virtual TBaseCell *Prev()
    { return NULL; }

    /** return number of moratr faces */
    int GetN_MortarFace()
    { return Domain->GetN_MortarFace(); }

    /** return coordinates of 'start'-point of mortar edge */
    void GetPoint(double &X, double &Y)
    {
      TVertex *Vert = ActiveCell->GetVertex(LocMortarFace[ActiveLevel]);
      X = Vert->GetX();
      Y = Vert->GetY();
    }
};

#endif
