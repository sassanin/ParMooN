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
// @(#)InnerEdge.h  
// 
// Class:       TInnerEdge
// Purpose:     class for inner edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#ifndef __INNEREDGE__
#define __INNEREDGE__

#include <Edge.h>

/** an edge in a 3D grid */
class TInnerEdge : public TEdge
{
  protected:

  public:

  /** constructor with neighbours */
   TInnerEdge(int n_Neibs, TBaseCell **neighbs);

};

#endif
