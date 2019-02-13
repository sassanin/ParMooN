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
   
// P3 and Q3 like behaviour on considered edge
static char C3_2_C3_2_Name[] = "C3_2_C3_2";
static char C3_2_C3_2_Desc[] = "conforming P3 or Q3 element";
static int C3_2_C3_2_N0 = 4;
static int C3_2_C3_2_N1 = 4;
static int C3_2_C3_2_NPairs = 4;
static int C3_2_C3_2_Pairs[][2] = { {0,7}, {1,6}, {2,5}, {3,4} };
static int C3_2_C3_2_NHanging = 0;
static int *C3_2_C3_2_Hanging = NULL;
static HNDesc *C3_2_C3_2_HangingTypes = NULL;
static int **C3_2_C3_2_Coupling = NULL;
static int C3_2_C3_2_NFarHanging = 0;
static int *C3_2_C3_2_FarHanging = NULL;
static HNDesc *C3_2_C3_2_FarHangingTypes = NULL;
static int ****C3_2_C3_2_FarCoupling = NULL;
static int C3_2_C3_2_NNoopposite = 0;
static int *C3_2_C3_2_Nopposite = NULL;
static int C3_2_C3_2_NNodes = 8;

TFE2DMapper *C3_2_C3_2 = new TFE2DMapper(C3_2_C3_2_Name, C3_2_C3_2_Desc,
                             C3_2_C3_2_N0, C3_2_C3_2_N1,
                             C3_2_C3_2_NPairs, (int *)C3_2_C3_2_Pairs,
                             C3_2_C3_2_NHanging, C3_2_C3_2_Hanging,
                             C3_2_C3_2_HangingTypes, C3_2_C3_2_Coupling,
                             C3_2_C3_2_NFarHanging, C3_2_C3_2_FarHanging,
                             C3_2_C3_2_FarHangingTypes, C3_2_C3_2_FarCoupling,
                             C3_2_C3_2_NNoopposite, C3_2_C3_2_Nopposite,
                             C3_2_C3_2_NNodes);
