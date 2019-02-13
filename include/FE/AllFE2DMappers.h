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
   
// ============================================================================
// FE mapper for local regular triangulation, same pattern on considered edge
// ============================================================================
#include <C0_2_C0_2.h>
#include <C1_2_C1_2.h>
#include <C2_2_C2_2.h>
#include <C3_2_C3_2.h>
#include <C4_2_C4_2.h>
#include <C5_2_C5_2.h>
#include <C6_2_C6_2.h>
#include <C7_2_C7_2.h>
#include <C8_2_C8_2.h>
#include <C9_2_C9_2.h>

#include <N1_2_N1_2.h>

#include <N2_2_N2_2.h>
#include <N3_2_N3_2.h>
#include <N4_2_N4_2.h>
#include <N5_2_N5_2.h>

// ============================================================================
// FE mapper for local one regular triangulation, same pattern on considered edge
// ============================================================================
#include <C0_2_C0_2_1Reg.h>
#include <C1_2_C1_2_1Reg.h>
#include <C2_2_C2_2_1Reg.h>
#include <C3_2_C3_2_1Reg.h>
#include <C4_2_C4_2_1Reg.h>
#include <C5_2_C5_2_1Reg.h>

#include <N1_2_N1_2_1Reg.h>
#include <N2_2_N2_2_1Reg.h>
#include <N3_2_N3_2_1Reg.h>
#include <N4_2_N4_2_1Reg.h>
#include <N5_2_N5_2_1Reg.h>



#include <N2_2_N1_2_1Reg.h>
#include <N3_2_N2_2_1Reg.h>
#include <N4_2_N3_2_1Reg.h>
#include <N5_2_N4_2_1Reg.h>
