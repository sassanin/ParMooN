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
// @(#)RefQuadToTri1Desc.h        1.1 10/30/98
//
// Class:       TRefQuadToTri1Desc
// Purpose:     refinement descriptor for refinement of a quadrangle
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//
// =======================================================================

#ifndef __REFQUADTOTRI1DESC__
#define __REFQUADTOTRI1DESC__

#include <RefDesc.h>

#ifndef __REFQUADTOTRI0DESC__
  #define QUADTTN_E         3
  #define QUADTTMAXN_VpC    3
  #define QUADTTMAXN_CpV    2
  #define QUADTTMAXN_EpC    3
  #define QUADTTMAXN_CpE    2
  #define QUADTTMAXN_EpV    3
  #define QUADTTMAXN_nVpoE  2
  #define QUADTTMAXN_nEpoE  1
#endif

/** refinement descriptor for refinement of a quadrangle into
    two triangles */
class TRefQuadToTri1Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor */
    TRefQuadToTri1Desc(TShapeDesc *shape);

    // Methods
};

#endif
