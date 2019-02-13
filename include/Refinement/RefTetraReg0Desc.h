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
// @(#)RefTetraReg0Desc.h        1.2 10/18/99
//
// Class:       TRefTetraRegDesc
// Purpose:     refinement descriptor for regular refinement of a
//              tetrahedron (variant 0)
//
// Authors:      Volker Behns  18.07.97
//               Matthias Ebelibg 06.09.99
//
// =======================================================================

#ifndef __REFTETRAREG0DESC__
#define __REFTETRAREG0DESC__

#include <RefDesc.h>

#define REFTETRAREGMAXN_EpV      7
#define REFTETRAREGMAXN_FpV     12
#define REFTETRAREGMAXN_VpF      3
#define REFTETRAREGMAXN_VpC      4
#define REFTETRAREGMAXN_CpV      6
#define REFTETRAREGMAXN_EpC      6
#define REFTETRAREGMAXN_CpE      4
#define REFTETRAREGMAXN_FpC      4
#define REFTETRAREGMAXN_CpF      2
#define REFTETRAREGMAXN_EpF      3
#define REFTETRAREGMAXN_FpE      4
#define REFTETRAREGMAXN_iVpE     1
#define REFTETRAREGMAXN_nVpoE    3
#define REFTETRAREGMAXN_nEpoE    2
#define REFTETRAREGMAXN_nFpoF    4
#define REFTETRAREGMAXN_nVpoF    6
#define REFTETRAREGMAXN_oVpoF    3
#define REFTETRAREGMAXN_nEpoF    9
#define REFTETRAREGMAXN_iEpF     3
#define TETRAN_V                 4
#define TETRAN_E                 6
#define TETRAN_F                 4

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraReg0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraReg0Desc(TShapeDesc *shape);

    // Methods
};

#endif
