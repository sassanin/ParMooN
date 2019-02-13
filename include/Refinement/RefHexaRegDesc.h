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
// @(#)RefHexaRegDesc.h        1.3 10/18/99
//
// Class:       TRefHexaRegDesc
// Purpose:     refinement descriptor for regular refinement of a
//              hexahedron
//
// Author:      Volker Behns  18.07.97
//
// =======================================================================

#ifndef __REFHEXAREGDESC__
#define __REFHEXAREGDESC__

#include <RefDesc.h>

#define REFHEXAREGMAXN_EpV      6
#define REFHEXAREGMAXN_FpV      12
#define REFHEXAREGMAXN_VpF      4
#define REFHEXAREGMAXN_VpC      8
#define REFHEXAREGMAXN_CpV      8
#define REFHEXAREGMAXN_EpC      12
#define REFHEXAREGMAXN_CpE      4
#define REFHEXAREGMAXN_FpC      6
#define REFHEXAREGMAXN_CpF      2
#define REFHEXAREGMAXN_EpF      4
#define REFHEXAREGMAXN_FpE      4
#define REFHEXAREGMAXN_iVpE     1
#define REFHEXAREGMAXN_oVpoF    4
#define REFHEXAREGMAXN_nVpoE    3
#define REFHEXAREGMAXN_nEpoE    2
#define REFHEXAREGMAXN_nFpoF    4
#define REFHEXAREGMAXN_nVpoF    9
#define REFHEXAREGMAXN_nEpoF    12
#define REFHEXAREGMAXN_iVpF     1
#define REFHEXAREGMAXN_iEpF     4
#define HEXAN_V                 8
#define HEXAN_E                 12
#define HEXAN_F                 6



/** refinement descriptor for regular refinement of a hexahedron */
class TRefHexaRegDesc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a hexahedron */
    TRefHexaRegDesc(TShapeDesc *shape);

    // Methods
};

#endif
