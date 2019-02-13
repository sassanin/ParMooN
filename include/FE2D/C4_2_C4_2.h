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
static char C4_2_C4_2_Name[] = "C4_2_C4_2";
static char C4_2_C4_2_Desc[] = "conforming P4 or Q4 element";
static int C4_2_C4_2_N0 = 5;
static int C4_2_C4_2_N1 = 5;
static int C4_2_C4_2_NPairs = 5;
static int C4_2_C4_2_Pairs[][2] = { {0,9}, {1,8}, {2,7}, {3,6}, {4,5} };
static int C4_2_C4_2_NHanging = 0;
static int *C4_2_C4_2_Hanging = NULL;
static HNDesc *C4_2_C4_2_HangingTypes = NULL;
static int **C4_2_C4_2_Coupling = NULL;
static int C4_2_C4_2_NFarHanging = 0;
static int *C4_2_C4_2_FarHanging = NULL;
static HNDesc *C4_2_C4_2_FarHangingTypes = NULL;
static int ****C4_2_C4_2_FarCoupling = NULL;
static int C4_2_C4_2_NNoopposite = 0;
static int *C4_2_C4_2_Nopposite = NULL;
static int C4_2_C4_2_NNodes = 10;

TFE2DMapper *C4_2_C4_2 = new TFE2DMapper(C4_2_C4_2_Name, C4_2_C4_2_Desc,
                             C4_2_C4_2_N0, C4_2_C4_2_N1,
                             C4_2_C4_2_NPairs, (int *)C4_2_C4_2_Pairs,
                             C4_2_C4_2_NHanging, C4_2_C4_2_Hanging,
                             C4_2_C4_2_HangingTypes, C4_2_C4_2_Coupling,
                             C4_2_C4_2_NFarHanging, C4_2_C4_2_FarHanging,
                             C4_2_C4_2_FarHangingTypes, C4_2_C4_2_FarCoupling,
                             C4_2_C4_2_NNoopposite, C4_2_C4_2_Nopposite,
                             C4_2_C4_2_NNodes);
