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
   
#include <NodalFunctional3D.h>

#include <NF_C_T_P00_3D.h>
#include <NF_C_T_P0_3D.h>
#include <NF_C_T_P1_3D.h>
#include <NF_C_T_P2_3D.h>
#include <NF_C_T_P3_3D.h>

#include <NF_C_T_B2_3D.h>

#include <NF_C_H_Q00_3D.h>
#include <NF_C_H_Q0_3D.h>
#include <NF_C_H_Q1_3D.h>
#include <NF_C_H_Q2_3D.h>
#include <NF_C_H_Q3_3D.h>
#include <NF_C_H_Q4_3D.h>

#include <NF_N_H_Q1_3D.h>
#include <NF_N_T_P1_3D.h>

#include <NF_D_T_P1_3D.h>
#include <NF_D_T_P2_3D.h>
#include <NF_D_T_P3_3D.h>

#include <NF_D_H_P1_3D.h>
#include <NF_D_H_P2_3D.h>
#include <NF_D_H_P3_3D.h>

#include <NF_D_H_Q1_3D.h>
#include <NF_D_H_Q2_3D.h>

#include <NF_B_H_IB2_3D.h>

// Superconvergence
#include <NF_S_H_Q2_3D.h>

#include <NF_N_T_P2_3D.h>

#include <NF_N_H_Q2_3D.h>
#include <NF_N_H_Q3_3D.h>
#include <NF_N_H_Q4_3D.h>

#include <NF_C_H_UL1_3D.h>
#include <NF_C_H_UL2_3D.h>
#include <NF_C_H_UL3_3D.h>

// Raviart-Thomas (RT) nodal functional, vector valued basis functions
#include <NF_N_T_RT0_3D.h>
#include <NF_N_T_RT1_3D.h>
#include <NF_N_T_RT2_3D.h>
#include <NF_N_T_RT3_3D.h>
#include <NF_N_H_RT0_3D.h>
#include <NF_N_H_RT1_3D.h>
#include <NF_N_H_RT2_3D.h>

// Brezzi-Douglas-Marini (BDM) nodal functional, vector valued basis functions
#include <NF_N_T_BDDF1_3D.h>
#include <NF_N_T_BDDF2_3D.h>
#include <NF_N_T_BDDF3_3D.h>
#include <NF_N_H_BDDF1_3D.h>
#include <NF_N_H_BDDF2_3D.h>
#include <NF_N_H_BDDF3_3D.h>
