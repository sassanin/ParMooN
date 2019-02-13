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
   
// P2 and Q2 like behaviour on considered edge
static char C2_2_C2_2_Name[] = "C2_2_C2_2";
static char C2_2_C2_2_Desc[] = "conforming P2 or Q2 element";
static int C2_2_C2_2_N0 = 3;
static int C2_2_C2_2_N1 = 3;
static int C2_2_C2_2_NPairs = 3;
static int C2_2_C2_2_Pairs[][2] = { {0,5}, {1,4}, {2,3} };
static int C2_2_C2_2_NHanging = 0;
static int *C2_2_C2_2_Hanging = NULL;
static HNDesc *C2_2_C2_2_HangingTypes = NULL;
static int **C2_2_C2_2_Coupling = NULL;
static int C2_2_C2_2_NFarHanging = 0;
static int *C2_2_C2_2_FarHanging = NULL;
static HNDesc *C2_2_C2_2_FarHangingTypes = NULL;
static int ****C2_2_C2_2_FarCoupling = NULL;
static int C2_2_C2_2_NNoopposite = 0;
static int *C2_2_C2_2_Nopposite = NULL;
static int C2_2_C2_2_NNodes = 6;

TFE2DMapper *C2_2_C2_2 = new TFE2DMapper(C2_2_C2_2_Name, C2_2_C2_2_Desc,
                             C2_2_C2_2_N0, C2_2_C2_2_N1,
                             C2_2_C2_2_NPairs, (int *)C2_2_C2_2_Pairs,
                             C2_2_C2_2_NHanging, C2_2_C2_2_Hanging,
                             C2_2_C2_2_HangingTypes, C2_2_C2_2_Coupling,
                             C2_2_C2_2_NFarHanging, C2_2_C2_2_FarHanging,
                             C2_2_C2_2_FarHangingTypes, C2_2_C2_2_FarCoupling,
                             C2_2_C2_2_NNoopposite, C2_2_C2_2_Nopposite,
                             C2_2_C2_2_NNodes);
