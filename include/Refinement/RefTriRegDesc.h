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
// @(#)RefTriRegDesc.h        1.1 10/30/98
//
// Class:       TRefTriRegDesc
// Purpose:     refinement descriptor for regular refinement of a triangle
//
// Author:      Volker Behns  17.07.97
//
// =======================================================================

#ifndef __REFTRIREGDESC__
#define __REFTRIREGDESC__

#include <RefDesc.h>

#define TRIRRN_E         3
#define TRIRRMAXN_VpC    3
#define TRIRRMAXN_CpV    3
#define TRIRRMAXN_EpC    3
#define TRIRRMAXN_CpE    2
#define TRIRRMAXN_EpV    4
#define TRIRRMAXN_iVpE   1
#define TRIRRMAXN_nVpoE  3
#define TRIRRMAXN_nEpoE  2

/** refinement descriptor for regular refinement of a triangle */
class TRefTriRegDesc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a triangle */
    TRefTriRegDesc(TShapeDesc *shape);

    // Methods
};

#endif
