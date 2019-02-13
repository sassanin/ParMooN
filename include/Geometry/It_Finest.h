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
// @(#)It_Finest.h        1.1 10/30/98
// 
// Class:       TIt_Finest
// Purpose:     iterator to produce a series of cells which lie
//              on top of cell tree
//
// Author:      Volker Behns  05.08.97
//
// =======================================================================

#ifndef __IT_Finest__
#define __IT_Finest__

#include <It_LE.h>

/** iterator to produce a series of cells which lie on
    top of cell tree */
class TIt_Finest : public TIt_LE
{
  protected:

  public:
    // Constructors

    // Methods
    /** initialize the iterator */
    virtual int Init(int level);

    /** return the maximum level */
    virtual int GetMaxLevel();
};

#endif
