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
// @(#)RefQuad2Conf0Desc.h        1.1 08/26/99
//
// Class:       TRefQuad2Conf0Desc
// Purpose:     refinement descriptor for conforming closure of a
//              quadrangle with 2 hanging node
//
// Author:      Matthias Ebeling  20.08.99
//
// =======================================================================

#ifndef __REFQUAD2CONF0DESC__
#define __REFQUAD2CONF0DESC__

#include <RefDesc.h>

#define QUADConfN_E          4
#define QUAD2ConfMAXN_VpC    4
#define QUAD2ConfMAXN_CpV    3
#define QUAD2ConfMAXN_EpC    4
#define QUAD2ConfMAXN_CpE    2
#define QUAD2ConfMAXN_EpV    3
#define QUADN_V              4
#define QUAD2ConfMAXN_iVpE   1
#define QUAD2ConfMAXN_nVpoE  3
#define QUAD2ConfMAXN_nEpoE  2

/** refinement descriptor for conforming closure of a quadrangle */
class TRefQuad2Conf0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for conforming closure of a quadrangle */
    TRefQuad2Conf0Desc(TShapeDesc *shape);

    // Methods
};

#endif
