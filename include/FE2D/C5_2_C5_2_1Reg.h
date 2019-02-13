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
   
// P5 and Q5 like behaviour on considered edge
static char C5_2_C5_2_1Reg_Name[] = "C5_2_C5_2_1Reg";
static char C5_2_C5_2_1Reg_Desc[] = "conforming P5 or Q5 element, one regular grid";
static int C5_2_C5_2_1Reg_N0 = 6;
static int C5_2_C5_2_1Reg_N1 = 6;
static int C5_2_C5_2_1Reg_N2 = 6;
static int C5_2_C5_2_1Reg_NMid = 1;
static int C5_2_C5_2_1Reg_Mid[][2] = { {11,12} };
static int C5_2_C5_2_1Reg_NPairs = 6;
static int C5_2_C5_2_1Reg_Pairs[][2] = { {0,17}, {1,15}, {2,13}, {3,10}, {4,8}, {5,6} };
static int C5_2_C5_2_1Reg_NHanging = 5;
static int C5_2_C5_2_1Reg_Hanging[] = { 7, 9, 11, 14, 16 };
static HNDesc C5_2_C5_2_1Reg_HangingTypes[] = { HN_C_P5_2D_0, HN_C_P5_2D_1,
                                                 HN_C_P5_2D_2, HN_C_P5_2D_3, HN_C_P5_2D_4 };
static int C5_2_C5_2_1Reg_Coupling_0[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_1[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_2[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_3[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_4[] = { 0, 1, 2, 3, 4, 5 };
static int *C5_2_C5_2_1Reg_Coupling[] = { C5_2_C5_2_1Reg_Coupling_0,
                                          C5_2_C5_2_1Reg_Coupling_1,
                                          C5_2_C5_2_1Reg_Coupling_2,
					  C5_2_C5_2_1Reg_Coupling_3,
					  C5_2_C5_2_1Reg_Coupling_4 };
static int C5_2_C5_2_1Reg_NFarHanging = 0;
static int *C5_2_C5_2_1Reg_FarHanging = NULL;
static HNDesc *C5_2_C5_2_1Reg_FarHangingTypes = NULL;
static int ****C5_2_C5_2_1Reg_FarCoupling = NULL;
static int C5_2_C5_2_1Reg_NNoopposite = 0;
static int *C5_2_C5_2_1Reg_Nopposite = NULL;
static int C5_2_C5_2_1Reg_NNodes = 18;

TFE2DMapper1Reg *C5_2_C5_2_1Reg = new TFE2DMapper1Reg(
                C5_2_C5_2_1Reg_Name, C5_2_C5_2_1Reg_Desc,
                C5_2_C5_2_1Reg_N0, C5_2_C5_2_1Reg_N1, C5_2_C5_2_1Reg_N2,
                C5_2_C5_2_1Reg_NPairs, (int *)C5_2_C5_2_1Reg_Pairs,
                C5_2_C5_2_1Reg_NMid, (int *)C5_2_C5_2_1Reg_Mid,
                C5_2_C5_2_1Reg_NHanging, C5_2_C5_2_1Reg_Hanging,
                C5_2_C5_2_1Reg_HangingTypes, C5_2_C5_2_1Reg_Coupling,
                C5_2_C5_2_1Reg_NFarHanging, C5_2_C5_2_1Reg_FarHanging,
                C5_2_C5_2_1Reg_FarHangingTypes, C5_2_C5_2_1Reg_FarCoupling,
                C5_2_C5_2_1Reg_NNoopposite, C5_2_C5_2_1Reg_Nopposite,
                C5_2_C5_2_1Reg_NNodes);
