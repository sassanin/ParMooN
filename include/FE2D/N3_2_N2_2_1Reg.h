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
static char N3_2_N2_2_1Reg_Name[] = "N3_2_N2_2_1Reg";
static char N3_2_N2_2_1Reg_Desc[] = "nonconforming P3 or Q3 element, one regular grid";
static int N3_2_N2_2_1Reg_N0 = 3;
static int N3_2_N2_2_1Reg_N1 = 2;
static int N3_2_N2_2_1Reg_N2 = 2;
static int N3_2_N2_2_1Reg_NMid = 0;
static int *N3_2_N2_2_1Reg_Mid = NULL;
static int N3_2_N2_2_1Reg_NPairs = 0;
static int *N3_2_N2_2_1Reg_Pairs = NULL;
static int N3_2_N2_2_1Reg_NHanging = 2;
static int N3_2_N2_2_1Reg_Hanging[] = { 0, 1 };
static HNDesc N3_2_N2_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0,
                                                HN_N_P2_2D_0 };
static int N3_2_N2_2_1Reg_Coupling_0[] = { 3,    5 };
static int N3_2_N2_2_1Reg_Coupling_1[] = { 3, 4, 5, 6 };
static int *N3_2_N2_2_1Reg_Coupling[] = { N3_2_N2_2_1Reg_Coupling_0,
                                          N3_2_N2_2_1Reg_Coupling_1 };
static int N3_2_N2_2_1Reg_NFarHanging = 0;
static int *N3_2_N2_2_1Reg_FarHanging = NULL;
static HNDesc *N3_2_N2_2_1Reg_FarHangingTypes = NULL;
static int ****N3_2_N2_2_1Reg_FarCoupling = NULL;
static int N3_2_N2_2_1Reg_NNoopposite = 5;
static int N3_2_N2_2_1Reg_Nopposite[] = { 2, 3, 4, 5, 6 };
static int N3_2_N2_2_1Reg_NNodes = 7;

TFE2DMapper1Reg *N3_2_N2_2_1Reg = new TFE2DMapper1Reg(
                N3_2_N2_2_1Reg_Name, N3_2_N2_2_1Reg_Desc,
                N3_2_N2_2_1Reg_N0, N3_2_N2_2_1Reg_N1, N3_2_N2_2_1Reg_N2,
                N3_2_N2_2_1Reg_NPairs, (int *)N3_2_N2_2_1Reg_Pairs,
                N3_2_N2_2_1Reg_NMid, (int *)N3_2_N2_2_1Reg_Mid,
                N3_2_N2_2_1Reg_NHanging, N3_2_N2_2_1Reg_Hanging,
                N3_2_N2_2_1Reg_HangingTypes, N3_2_N2_2_1Reg_Coupling,
                N3_2_N2_2_1Reg_NFarHanging, N3_2_N2_2_1Reg_FarHanging,
                N3_2_N2_2_1Reg_FarHangingTypes, N3_2_N2_2_1Reg_FarCoupling,
                N3_2_N2_2_1Reg_NNoopposite, N3_2_N2_2_1Reg_Nopposite,
                N3_2_N2_2_1Reg_NNodes);
