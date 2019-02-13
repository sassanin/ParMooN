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
   
// nonconforming P1 and Q1 like behaviour on considered edge
static char N4_2_N4_2_1Reg_Name[] = "N4_2_N4_2_1Reg";
static char N4_2_N4_2_1Reg_Desc[] = "nonconforming P4 or Q4 element, one regular grid";
static int N4_2_N4_2_1Reg_N0 = 4;
static int N4_2_N4_2_1Reg_N1 = 4;
static int N4_2_N4_2_1Reg_N2 = 4;
static int N4_2_N4_2_1Reg_NMid = 0;
static int *N4_2_N4_2_1Reg_Mid = NULL;
static int N4_2_N4_2_1Reg_NPairs = 0;
static int *N4_2_N4_2_1Reg_Pairs = NULL;
static int N4_2_N4_2_1Reg_NHanging = 4;
static int N4_2_N4_2_1Reg_Hanging[] = { 0, 1, 2, 3 };
static HNDesc N4_2_N4_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0, HN_N_P2_2D_0,
                                                HN_N_P3_2D_0, HN_N_P4_2D_0 };
static int N4_2_N4_2_1Reg_Coupling_0[] = { 4, 8 };
static int N4_2_N4_2_1Reg_Coupling_1[] = { 4, 5, 8, 9 };
static int N4_2_N4_2_1Reg_Coupling_2[] = { 4, 5, 6, 8, 9, 10 };
static int N4_2_N4_2_1Reg_Coupling_3[] = { 4, 5, 6, 7, 8, 9, 10, 11  };
static int *N4_2_N4_2_1Reg_Coupling[] = { N4_2_N4_2_1Reg_Coupling_0,
                                          N4_2_N4_2_1Reg_Coupling_1,
                                          N4_2_N4_2_1Reg_Coupling_2,
                                          N4_2_N4_2_1Reg_Coupling_3 };
static int N4_2_N4_2_1Reg_NFarHanging = 0;
static int *N4_2_N4_2_1Reg_FarHanging = NULL;
static HNDesc *N4_2_N4_2_1Reg_FarHangingTypes = NULL;
static int ****N4_2_N4_2_1Reg_FarCoupling = NULL;
static int N4_2_N4_2_1Reg_NNoopposite = 8;
static int N4_2_N4_2_1Reg_Nopposite[] = { 4, 5, 6, 7, 8, 9, 10, 11 };
static int N4_2_N4_2_1Reg_NNodes = 12;

TFE2DMapper1Reg *N4_2_N4_2_1Reg = new TFE2DMapper1Reg(
                N4_2_N4_2_1Reg_Name, N4_2_N4_2_1Reg_Desc,
                N4_2_N4_2_1Reg_N0, N4_2_N4_2_1Reg_N1, N4_2_N4_2_1Reg_N2,
                N4_2_N4_2_1Reg_NPairs, (int *)N4_2_N4_2_1Reg_Pairs,
                N4_2_N4_2_1Reg_NMid, (int *)N4_2_N4_2_1Reg_Mid,
                N4_2_N4_2_1Reg_NHanging, N4_2_N4_2_1Reg_Hanging,
                N4_2_N4_2_1Reg_HangingTypes, N4_2_N4_2_1Reg_Coupling,
                N4_2_N4_2_1Reg_NFarHanging, N4_2_N4_2_1Reg_FarHanging,
                N4_2_N4_2_1Reg_FarHangingTypes, N4_2_N4_2_1Reg_FarCoupling,
                N4_2_N4_2_1Reg_NNoopposite, N4_2_N4_2_1Reg_Nopposite,
                N4_2_N4_2_1Reg_NNodes);
