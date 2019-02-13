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
// @(#)RefTriBis1Desc.h        1.1 10/30/98
//
// Class:       TRefTriBis1Desc
// Purpose:     refinement descriptor for bisection of a triangle
//              bisection of edge 1
//
// Author:      Volker Behns  18.07.97
//
// =======================================================================

#ifndef __REFTRIBIS1DESC__
#define __REFTRIBIS1DESC__

#include <RefDesc.h>

#define TRIBI1MAXN_VpC    3
#define TRIBI1MAXN_CpV    2
#define TRIBI1MAXN_EpC    3
#define TRIBI1MAXN_CpE    2
#define TRIBI1MAXN_EpV    3
#define TRIBI1MAXN_iVpE   1
#define TRIBI1MAXN_nVpoE  3
#define TRIBI1MAXN_nEpoE  2
#define TRIBI1N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis1Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 1 */
    TRefTriBis1Desc(TShapeDesc *shape);

    // Methods
};

#endif
