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
   
// hanging node for conforming P3/Q3 elements

static int C_P3_2D_N_Nodes_0 = 4;
static double C_P3_2D_Coupling_0[] = {  0.0625, -0.3125,  0.9375,  0.3125 };

static int C_P3_2D_N_Nodes_1 = 4;
static double C_P3_2D_Coupling_1[] = { -0.0625,  0.5625,  0.5625, -0.0625 };

static int C_P3_2D_N_Nodes_2 = 4;
static double C_P3_2D_Coupling_2[] = {  0.3125,  0.9375, -0.3125,  0.0625 };

THNDesc *HN_C_P3_2D_0_Obj = new THNDesc(C_P3_2D_N_Nodes_0, C_P3_2D_Coupling_0);
THNDesc *HN_C_P3_2D_1_Obj = new THNDesc(C_P3_2D_N_Nodes_1, C_P3_2D_Coupling_1);
THNDesc *HN_C_P3_2D_2_Obj = new THNDesc(C_P3_2D_N_Nodes_2, C_P3_2D_Coupling_2);
