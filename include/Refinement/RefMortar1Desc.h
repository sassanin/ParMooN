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
// @(#)RefMortar1Desc.h        1.1 10/30/98
//
// Class:       TRefMortar1Desc
// Purpose:     refinement descriptor for a mortar cell (i.e. generate a
//              Nx1 grid with base edge 1)
//
// Author:      Volker Behns  30.12.97
//
// =======================================================================

#ifndef __REFMORTAR1DESC__
#define __REFMORTAR1DESC__

#include <RefDesc.h>

#ifndef __REFMORTAR0DESC__
  #define MORTARRMAXN_VpC    4
  #define MORTARRMAXN_CpV    2
  #define MORTARRMAXN_EpC    4
  #define MORTARRMAXN_CpE    2
  #define MORTARRMAXN_EpV    3
#endif

/** refinement descriptor for a mortar cell */
class TRefMortar1Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for mortar refinement of a quadrangle */
    TRefMortar1Desc(TShapeDesc *shape, int Mortar_Ni, int N);

    // Methods
};

#endif
