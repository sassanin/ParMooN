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
static char N2_2_N1_2_1Reg_Name[] = "N2_2_N1_2_1Reg";
static char N2_2_N1_2_1Reg_Desc[] = "nonconforming P2 or Q2 element, one regular grid";
static int N2_2_N1_2_1Reg_N0 = 2;
static int N2_2_N1_2_1Reg_N1 = 1;
static int N2_2_N1_2_1Reg_N2 = 1;
static int N2_2_N1_2_1Reg_NMid = 0;
static int *N2_2_N1_2_1Reg_Mid = NULL;
static int N2_2_N1_2_1Reg_NPairs = 0;
static int *N2_2_N1_2_1Reg_Pairs = NULL;
static int N2_2_N1_2_1Reg_NHanging = 1;
static int N2_2_N1_2_1Reg_Hanging[] = { 0 };
static HNDesc N2_2_N1_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0 };
static int N2_2_N1_2_1Reg_Coupling_0[] = { 2, 3 };
static int *N2_2_N1_2_1Reg_Coupling[] = { N2_2_N1_2_1Reg_Coupling_0 };
static int N2_2_N1_2_1Reg_NFarHanging = 0;
static int *N2_2_N1_2_1Reg_FarHanging = NULL;
static HNDesc *N2_2_N1_2_1Reg_FarHangingTypes = NULL;
static int ****N2_2_N1_2_1Reg_FarCoupling = NULL;
static int N2_2_N1_2_1Reg_NNoopposite = 3;
static int N2_2_N1_2_1Reg_Nopposite[] = { 1, 2, 3};
static int N2_2_N1_2_1Reg_NNodes = 4;

TFE2DMapper1Reg *N2_2_N1_2_1Reg = new TFE2DMapper1Reg(
                N2_2_N1_2_1Reg_Name, N2_2_N1_2_1Reg_Desc,
                N2_2_N1_2_1Reg_N0, N2_2_N1_2_1Reg_N2, N2_2_N1_2_1Reg_N2,
                N2_2_N1_2_1Reg_NPairs, (int *)N2_2_N1_2_1Reg_Pairs,
                N2_2_N1_2_1Reg_NMid, (int *)N2_2_N1_2_1Reg_Mid,
                N2_2_N1_2_1Reg_NHanging, N2_2_N1_2_1Reg_Hanging,
                N2_2_N1_2_1Reg_HangingTypes, N2_2_N1_2_1Reg_Coupling,
                N2_2_N1_2_1Reg_NFarHanging, N2_2_N1_2_1Reg_FarHanging,
                N2_2_N1_2_1Reg_FarHangingTypes, N2_2_N1_2_1Reg_FarCoupling,
                N2_2_N1_2_1Reg_NNoopposite, N2_2_N1_2_1Reg_Nopposite,
                N2_2_N1_2_1Reg_NNodes);
