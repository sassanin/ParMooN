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
// @(#)MortarJoint.C        1.1 10/30/98
// 
// Class:       TMortarJoint
// Purpose:     indicates a mortar joint
//
// Author:      Volker Behns  19.03.98
//
// =======================================================================

#ifndef __MORTAR__
#define __MORTAR__
#endif

#include <MortarJoint.h>

// Constructors
TMortarJoint::TMortarJoint()
{
  ID = MortarJoint;
}

// Methods
int TMortarJoint::CheckMatchingRef(TBaseCell *Me, int J_i,
                    struct StoreGeom &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#ifdef __MORTAR__
int TMortarJoint::CheckMatchingRef(TBaseCell *Me, int J_i,
                    StoreGeomMortar &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#endif
