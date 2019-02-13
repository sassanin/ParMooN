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
   
#include <BaseFunct2D.h>

#include <BF_C_T_P00_2D.h>
#include <BF_C_T_P0_2D.h>
#include <BF_C_T_P1_2D.h>
#include <BF_C_T_P2_2D.h>
#include <BF_C_T_P3_2D.h>
#include <BF_C_T_P4_2D.h>
#include <BF_C_T_P5_2D.h>
#include <BF_C_T_P6_2D.h>

#include <BF_C_Q_Q00_2D.h>
#include <BF_C_Q_Q0_2D.h>
#include <BF_C_Q_Q1_2D.h>
#include <BF_C_Q_Q2_2D.h>
#include <BF_C_Q_Q3_2D.h>
#include <BF_C_Q_Q4_2D.h>
#include <BF_C_Q_Q5_2D.h>
#include <BF_C_Q_Q6_2D.h>
#include <BF_C_Q_Q7_2D.h>
#include <BF_C_Q_Q8_2D.h>
#include <BF_C_Q_Q9_2D.h>

#include <BF_N_T_P1_2D.h>
#include <BF_N_Q_Q1_2D.h>

#include <BF_D_Q_P1_2D.h>
#include <BF_D_Q_P2_2D.h>
#include <BF_D_Q_P3_2D.h>

#include <BF_C_T_B2_2D.h>
#include <BF_C_T_B3_2D.h>
#include <BF_C_T_B4_2D.h>

#include <BF_D_T_P1_2D.h>
#include <BF_D_T_P2_2D.h>

#include <BF_N_Q_Q2_2D.h>
#include <BF_N_Q_Q3_2D.h>
#include <BF_N_Q_Q4_2D.h>
#include <BF_N_Q_Q5_2D.h>

#include <BF_D_Q_P4_2D.h>
#include <BF_D_Q_P5_2D.h>
#include <BF_D_Q_P6_2D.h>
#include <BF_D_Q_P7_2D.h>

#include <BF_N_T_P1MOD_2D.h>

#include <BF_C_T_P1MINI_2D.h>

#include <BF_N_T_P2_2D.h>
#include <BF_N_T_P3_2D.h>
#include <BF_N_T_P4_2D.h>
#include <BF_N_T_P5_2D.h>

#include <BF_D_T_P3_2D.h>
#include <BF_D_T_P4_2D.h>

#include <BF_B_Q_IB2_2D.h>

#include <BF_D_Q_Q1_2D.h>
#include <BF_D_Q_Q2_2D.h>
#include <BF_D_Q_Q3_2D.h>
#include <BF_D_Q_Q4_2D.h>

#include <BF_D_Q_D2_2D.h>

#include <BF_D_T_SV1_2D.h>
#include <BF_C_T_SV2_2D.h>

#include <BF_C_T_RR2_2D.h>

//==== localProjection BaseFunctions
#include <BF_C_Q_UL1_2D.h>
#include <BF_C_Q_UL2_2D.h>
#include <BF_C_Q_UL3_2D.h>
#include <BF_C_Q_UL4_2D.h>
#include <BF_C_Q_UL5_2D.h>

#include <BF_C_T_UL1_2D.h>
#include <BF_C_T_UL2_2D.h>
#include <BF_C_T_UL3_2D.h>
#include <BF_C_T_UL4_2D.h>
#include <BF_C_T_UL5_2D.h>

#include <BF_C_Q_UL2S_2D.h>
#include <BF_C_Q_UL3S_2D.h>
#include <BF_C_Q_UL4S_2D.h>
#include <BF_C_Q_UL5S_2D.h>
#include <BF_C_Q_UL6S_2D.h>
#include <BF_C_Q_UL7S_2D.h>
#include <BF_C_Q_UL8S_2D.h>
#include <BF_C_Q_UL9S_2D.h>

#include <BF_C_Q_UL2SE_2D.h>
#include <BF_C_Q_UL3SE_2D.h>
#include <BF_C_Q_UL4SE_2D.h>
#include <BF_C_Q_UL5SE_2D.h>
#include <BF_C_Q_UL6SE_2D.h>
#include <BF_C_Q_UL7SE_2D.h>
#include <BF_C_Q_UL8SE_2D.h>
#include <BF_C_Q_UL9SE_2D.h>

#include <BF_C_Q_M2_2D.h>
#include <BF_C_Q_M3_2D.h>
#include <BF_C_Q_M4_2D.h>
#include <BF_C_Q_M5_2D.h>
#include <BF_C_Q_M6_2D.h>
#include <BF_C_Q_M7_2D.h>
#include <BF_C_Q_M8_2D.h>
#include <BF_C_Q_M9_2D.h>

#include <BF_C_Q_EL1_2D.h>

// Raviart-Thomas (RT) basis functions, vector valued basis functions
#include <BF_N_Q_RT0_2D.h>
#include <BF_N_Q_RT1_2D.h>
#include <BF_N_Q_RT2_2D.h>
#include <BF_N_Q_RT3_2D.h>
#include <BF_N_T_RT0_2D.h>
#include <BF_N_T_RT1_2D.h>
#include <BF_N_T_RT2_2D.h>
#include <BF_N_T_RT3_2D.h>
// Brezzi-Douglas-Marini (BDM) basis functions, vector valued basis functions
#include <BF_N_Q_BDM1_2D.h>
#include <BF_N_Q_BDM2_2D.h>
#include <BF_N_Q_BDM3_2D.h>
#include <BF_N_T_BDM1_2D.h>
#include <BF_N_T_BDM2_2D.h>
#include <BF_N_T_BDM3_2D.h>
