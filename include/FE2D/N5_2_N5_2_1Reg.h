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
static char N5_2_N5_2_1Reg_Name[] = "N5_2_N5_2_1Reg";
static char N5_2_N5_2_1Reg_Desc[] = "nonconforming P5 or Q5 element, one regular grid";
static int N5_2_N5_2_1Reg_N0 = 5;
static int N5_2_N5_2_1Reg_N1 = 5;
static int N5_2_N5_2_1Reg_N2 = 5;
static int N5_2_N5_2_1Reg_NMid = 0;
static int *N5_2_N5_2_1Reg_Mid = NULL;
static int N5_2_N5_2_1Reg_NPairs = 0;
static int *N5_2_N5_2_1Reg_Pairs = NULL;
static int N5_2_N5_2_1Reg_NHanging = 5;
static int N5_2_N5_2_1Reg_Hanging[] = { 0, 1, 2, 3, 4  };
static HNDesc N5_2_N5_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0, HN_N_P2_2D_0,
                                                HN_N_P3_2D_0, HN_N_P4_2D_0,
                                                HN_N_P5_2D_0 };
static int N5_2_N5_2_1Reg_Coupling_0[] = { 5, 10 };
static int N5_2_N5_2_1Reg_Coupling_1[] = { 5, 6, 10, 11 };
static int N5_2_N5_2_1Reg_Coupling_2[] = { 5, 6, 7, 10, 11, 12 };
static int N5_2_N5_2_1Reg_Coupling_3[] = { 5, 6, 7, 8, 10, 11, 12, 13 };
static int N5_2_N5_2_1Reg_Coupling_4[] = { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
static int *N5_2_N5_2_1Reg_Coupling[] = { N5_2_N5_2_1Reg_Coupling_0,
                                          N5_2_N5_2_1Reg_Coupling_1,
                                          N5_2_N5_2_1Reg_Coupling_2,
                                          N5_2_N5_2_1Reg_Coupling_3,
                                          N5_2_N5_2_1Reg_Coupling_4 };
static int N5_2_N5_2_1Reg_NFarHanging = 0;
static int *N5_2_N5_2_1Reg_FarHanging = NULL;
static HNDesc *N5_2_N5_2_1Reg_FarHangingTypes = NULL;
static int ****N5_2_N5_2_1Reg_FarCoupling = NULL;
static int N5_2_N5_2_1Reg_NNoopposite = 10;
static int N5_2_N5_2_1Reg_Nopposite[] = { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
static int N5_2_N5_2_1Reg_NNodes = 15;

TFE2DMapper1Reg *N5_2_N5_2_1Reg = new TFE2DMapper1Reg(
                N5_2_N5_2_1Reg_Name, N5_2_N5_2_1Reg_Desc,
                N5_2_N5_2_1Reg_N0, N5_2_N5_2_1Reg_N1, N5_2_N5_2_1Reg_N2,
                N5_2_N5_2_1Reg_NPairs, (int *)N5_2_N5_2_1Reg_Pairs,
                N5_2_N5_2_1Reg_NMid, (int *)N5_2_N5_2_1Reg_Mid,
                N5_2_N5_2_1Reg_NHanging, N5_2_N5_2_1Reg_Hanging,
                N5_2_N5_2_1Reg_HangingTypes, N5_2_N5_2_1Reg_Coupling,
                N5_2_N5_2_1Reg_NFarHanging, N5_2_N5_2_1Reg_FarHanging,
                N5_2_N5_2_1Reg_FarHangingTypes, N5_2_N5_2_1Reg_FarCoupling,
                N5_2_N5_2_1Reg_NNoopposite, N5_2_N5_2_1Reg_Nopposite,
                N5_2_N5_2_1Reg_NNodes);
