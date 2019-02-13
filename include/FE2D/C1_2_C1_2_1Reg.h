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
   
// P1 and Q1 like behaviour on considered edge
static char C1_2_C1_2_1Reg_Name[] = "C1_2_C1_2_1Reg";
static char C1_2_C1_2_1Reg_Desc[] = "conforming P1 or Q1 element, one regular grid";
static int C1_2_C1_2_1Reg_N0 = 2;
static int C1_2_C1_2_1Reg_N1 = 2;
static int C1_2_C1_2_1Reg_N2 = 2;
static int C1_2_C1_2_1Reg_NMid = 1;
static int C1_2_C1_2_1Reg_Mid[][2] = { {3,4} };
static int C1_2_C1_2_1Reg_NPairs = 2;
static int C1_2_C1_2_1Reg_Pairs[][2] = { {0,5}, {1,2} };
static int C1_2_C1_2_1Reg_NHanging = 1;
static int C1_2_C1_2_1Reg_Hanging[] = { 3 };
static HNDesc C1_2_C1_2_1Reg_HangingTypes[] = { HN_C_P1_2D_0 };
static int C1_2_C1_2_1Reg_Coupling_0[] = { 0, 1 };
static int *C1_2_C1_2_1Reg_Coupling[] = { C1_2_C1_2_1Reg_Coupling_0 };
static int C1_2_C1_2_1Reg_NFarHanging = 0;
static int *C1_2_C1_2_1Reg_FarHanging = NULL;
static HNDesc *C1_2_C1_2_1Reg_FarHangingTypes = NULL;
static int ****C1_2_C1_2_1Reg_FarCoupling = NULL;
static int C1_2_C1_2_1Reg_NNoopposite = 0;
static int *C1_2_C1_2_1Reg_Nopposite = NULL;
static int C1_2_C1_2_1Reg_NNodes = 6;

TFE2DMapper1Reg *C1_2_C1_2_1Reg = new TFE2DMapper1Reg(
                C1_2_C1_2_1Reg_Name, C1_2_C1_2_1Reg_Desc,
                C1_2_C1_2_1Reg_N0, C1_2_C1_2_1Reg_N1, C1_2_C1_2_1Reg_N2,
                C1_2_C1_2_1Reg_NPairs, (int *)C1_2_C1_2_1Reg_Pairs,
                C1_2_C1_2_1Reg_NMid, (int *)C1_2_C1_2_1Reg_Mid,
                C1_2_C1_2_1Reg_NHanging, C1_2_C1_2_1Reg_Hanging,
                C1_2_C1_2_1Reg_HangingTypes, C1_2_C1_2_1Reg_Coupling,
                C1_2_C1_2_1Reg_NFarHanging, C1_2_C1_2_1Reg_FarHanging,
                C1_2_C1_2_1Reg_FarHangingTypes, C1_2_C1_2_1Reg_FarCoupling,
                C1_2_C1_2_1Reg_NNoopposite, C1_2_C1_2_1Reg_Nopposite,
                C1_2_C1_2_1Reg_NNodes);
