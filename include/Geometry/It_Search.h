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
// @(#)It_Search.h        1.1 10/30/98
// 
// Class:       TIt_Search
// Purpose:     class of iterators which search the cell tree by 
//              graph-theoretic methods
//
// Author:      Volker Behns  01.09.97
//
// =======================================================================

#ifndef __IT_SEARCH__
#define __IT_SEARCH__

#include <Iterator.h>

/** class of iterators which search the cell tree by
    graph-theoretic methods */
class TIt_Search : public TIterator
{
  protected:
    /** status vector */
    struct {int N_Children, CurrentChild;} Status[MAX_ItLevel];

    /** active level in status vector */
    int ActiveLevel;
    /** active root cell */
    int ActiveRootCell;

    /** active cell */
    TBaseCell *ActiveCell;

  public:
    // Constructors

    // Methods
    /** Initialize on level */
    virtual int Init(int level);

    /** return the maximum level */
    virtual int GetMaxLevel();
};

#endif
