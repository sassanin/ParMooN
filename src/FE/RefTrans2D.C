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
// @(#)RefTrans2D.C        1.1 10/30/98
//
// Class:      TRefTrans2D
//
// Purpose:    reference transformations for 2D geometric objects
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#include <RefTrans2D.h>

RefTrans2D TRefTrans2D::FindRefTrans2D
        (int N_LocalUsedElements, FE2D *LocalUsedElements)
{
  RefTrans2D rf;

  if( ((int)(LocalUsedElements[0])) < 5)
    rf = TriaAffin;
  else
    if( ((int)(LocalUsedElements[0])) < 10)
      rf = QuadAffin;
    else
      rf = QuadBilinear;

  return rf;
}
