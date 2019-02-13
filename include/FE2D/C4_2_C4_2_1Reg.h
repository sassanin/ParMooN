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
   
// P4 and Q4 like behaviour on considered edge
static char C4_2_C4_2_1Reg_Name[] = "C4_2_C4_2_1Reg";
static char C4_2_C4_2_1Reg_Desc[] = "conforming P4 or Q4 element, one regular grid";
static int C4_2_C4_2_1Reg_N0 = 5;
static int C4_2_C4_2_1Reg_N1 = 5;
static int C4_2_C4_2_1Reg_N2 = 5;
static int C4_2_C4_2_1Reg_NMid = 1;
static int C4_2_C4_2_1Reg_Mid[][2] = { {9,10} };
static int C4_2_C4_2_1Reg_NPairs = 5;
static int C4_2_C4_2_1Reg_Pairs[][2] = { {0,14}, {1,12}, {2,9}, {3,7}, {4,5} };
static int C4_2_C4_2_1Reg_NHanging = 4;
static int C4_2_C4_2_1Reg_Hanging[] = { 6, 8, 11, 13 };
static HNDesc C4_2_C4_2_1Reg_HangingTypes[] = { HN_C_P4_2D_0, HN_C_P4_2D_1,
                                                 HN_C_P4_2D_2, HN_C_P4_2D_3 };
static int C4_2_C4_2_1Reg_Coupling_0[] = { 0, 1, 2, 3, 4 };
static int C4_2_C4_2_1Reg_Coupling_1[] = { 0, 1, 2, 3, 4 };
static int C4_2_C4_2_1Reg_Coupling_2[] = { 0, 1, 2, 3, 4 };
static int C4_2_C4_2_1Reg_Coupling_3[] = { 0, 1, 2, 3, 4 };
static int *C4_2_C4_2_1Reg_Coupling[] = { C4_2_C4_2_1Reg_Coupling_0,
                                          C4_2_C4_2_1Reg_Coupling_1,
                                          C4_2_C4_2_1Reg_Coupling_2,
      					  C4_2_C4_2_1Reg_Coupling_3 };
static int C4_2_C4_2_1Reg_NFarHanging = 0;
static int *C4_2_C4_2_1Reg_FarHanging = NULL;
static HNDesc *C4_2_C4_2_1Reg_FarHangingTypes = NULL;
static int ****C4_2_C4_2_1Reg_FarCoupling = NULL;
static int C4_2_C4_2_1Reg_NNoopposite = 0;
static int *C4_2_C4_2_1Reg_Nopposite = NULL;
static int C4_2_C4_2_1Reg_NNodes = 15;

TFE2DMapper1Reg *C4_2_C4_2_1Reg = new TFE2DMapper1Reg(
                C4_2_C4_2_1Reg_Name, C4_2_C4_2_1Reg_Desc,
                C4_2_C4_2_1Reg_N0, C4_2_C4_2_1Reg_N1, C4_2_C4_2_1Reg_N2,
                C4_2_C4_2_1Reg_NPairs, (int *)C4_2_C4_2_1Reg_Pairs,
                C4_2_C4_2_1Reg_NMid, (int *)C4_2_C4_2_1Reg_Mid,
                C4_2_C4_2_1Reg_NHanging, C4_2_C4_2_1Reg_Hanging,
                C4_2_C4_2_1Reg_HangingTypes, C4_2_C4_2_1Reg_Coupling,
                C4_2_C4_2_1Reg_NFarHanging, C4_2_C4_2_1Reg_FarHanging,
                C4_2_C4_2_1Reg_FarHangingTypes, C4_2_C4_2_1Reg_FarCoupling,
                C4_2_C4_2_1Reg_NNoopposite, C4_2_C4_2_1Reg_Nopposite,
                C4_2_C4_2_1Reg_NNodes);
