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
#include <P1_P1.h>
#include <P2_P2.h>
#include <P3_P3.h>

#include <Q1_Q1.h>
#include <Q2_Q2.h>
#include <Q3_Q3.h>
#include <Q4_Q4.h>

#include <N1_N1.h>
#include <NP1_NP1.h>
#include <NP2_NP2.h>

#include <NP1_NP1_1Reg.h>
#include <NP2_NP2_1Reg.h>

#include <D_D.h>

#include <Q1_Q1_1Reg.h>
#include <Q2_Q2_1Reg.h>
#include <Q3_Q3_1Reg.h>

#include <P1_P1_1Reg.h>
#include <P2_P2_1Reg.h>
#include <P3_P3_1Reg.h>

#include <N2_N2.h>
#include <N3_N3.h>
#include <N4_N4.h>

#include <P2B_P2B.h>
