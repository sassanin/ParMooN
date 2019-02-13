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
   
// P0 and Q0 like behaviour on considered edge
static char C0_2_C0_2_1Reg_Name[] = "C0_2_C0_2_1Reg";
static char C0_2_C0_2_1Reg_Desc[] = "conforming P0 or Q0 element, one regular grid";
static int C0_2_C0_2_1Reg_N0 = 0;
static int C0_2_C0_2_1Reg_N1 = 0;
static int C0_2_C0_2_1Reg_N2 = 0;
static int C0_2_C0_2_1Reg_NMid = 0;
static int *C0_2_C0_2_1Reg_Mid = NULL;
static int C0_2_C0_2_1Reg_NPairs = 0;
static int *C0_2_C0_2_1Reg_Pairs = NULL;
static int C0_2_C0_2_1Reg_NHanging = 0;
static int *C0_2_C0_2_1Reg_Hanging = NULL;
static HNDesc *C0_2_C0_2_1Reg_HangingTypes = NULL;
static int **C0_2_C0_2_1Reg_Coupling = NULL;
static int C0_2_C0_2_1Reg_NFarHanging = 0;
static int *C0_2_C0_2_1Reg_FarHanging = NULL;
static HNDesc *C0_2_C0_2_1Reg_FarHangingTypes = NULL;
static int ****C0_2_C0_2_1Reg_FarCoupling = NULL;
static int C0_2_C0_2_1Reg_NNoopposite = 0;
static int *C0_2_C0_2_1Reg_Nopposite = NULL;
static int C0_2_C0_2_1Reg_NNodes = 0;

TFE2DMapper1Reg *C0_2_C0_2_1Reg = new TFE2DMapper1Reg(
                C0_2_C0_2_1Reg_Name, C0_2_C0_2_1Reg_Desc,
                C0_2_C0_2_1Reg_N0, C0_2_C0_2_1Reg_N1, C0_2_C0_2_1Reg_N2,
                C0_2_C0_2_1Reg_NPairs, (int *)C0_2_C0_2_1Reg_Pairs,
                C0_2_C0_2_1Reg_NMid, (int *)C0_2_C0_2_1Reg_Mid,
                C0_2_C0_2_1Reg_NHanging, C0_2_C0_2_1Reg_Hanging,
                C0_2_C0_2_1Reg_HangingTypes, C0_2_C0_2_1Reg_Coupling,
                C0_2_C0_2_1Reg_NFarHanging, C0_2_C0_2_1Reg_FarHanging,
                C0_2_C0_2_1Reg_FarHangingTypes, C0_2_C0_2_1Reg_FarCoupling,
                C0_2_C0_2_1Reg_NNoopposite, C0_2_C0_2_1Reg_Nopposite,
                C0_2_C0_2_1Reg_NNodes);
