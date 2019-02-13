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
// @(#)RefNoRef.h        1.1 10/30/98
//
// Class:       TRefNoRef
// Purpose:     no refinement - only path to shape descriptor
//
// Author:      Volker Behns  21.07.97
//
// =======================================================================

#ifndef __REFNOREF__
#define __REFNOREF__

#include <RefDesc.h>

/** no refinement - only path to shape descriptor */
class TRefNoRef : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a nonrefined object */
    TRefNoRef(TShapeDesc *shape);

    // Methods
    /** return TRUE because object is to refine */
    virtual int IsToRefine()
    { return FALSE; }

};

#endif
