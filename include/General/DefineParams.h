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
// @(#)DefineParams.h        1.3 10/18/99
// 
// Purpose:     defines for mortar
//
// Author:      Volker Behns 13.09.99
//
// =======================================================================

#ifdef __MORTAR__
  // additional link term for SDFEM
  //#define __ADD_LINK_SDFEM__
  // additional link term for upwind
  //#define __ADD_LINK_UPW__

  // continuity in cross points
  #define __CONNECT_CROSSPOINTS__
#endif // __MORTAR__

// correct new midpoints (regular refinement of quadrangles)
#define __CORRECT_MP__

// -------------- end of changable parameters ------------------

#ifdef __ADD_LINK_SDFEM__
#ifndef __ADD_LINK__
// extension of matrix structure
#define __ADD_LINK__
#endif
#endif

#ifdef __ADD_LINK_UPW__
#ifndef __ADD_LINK__
// extension of matrix structure
#define __ADD_LINK__
#endif
#endif
