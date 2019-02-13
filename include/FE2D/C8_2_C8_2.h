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
   
// P8 and Q8 like behaviour on considered edge
static char C8_2_C8_2_Name[] = "C8_2_C8_2";
static char C8_2_C8_2_Desc[] = "conforming P8 or Q8 element";
static int C8_2_C8_2_N0 = 9;
static int C8_2_C8_2_N1 = 9;
static int C8_2_C8_2_NPairs = 9;
static int C8_2_C8_2_Pairs[][2] = { {0,17}, {1,16}, {2,15}, {3,14},
                                {4,13}, {5,12}, {6,11}, {7,10}, {8,9}};
static int C8_2_C8_2_NHanging = 0;
static int *C8_2_C8_2_Hanging = NULL;
static HNDesc *C8_2_C8_2_HangingTypes = NULL;
static int **C8_2_C8_2_Coupling = NULL;
static int C8_2_C8_2_NFarHanging = 0;
static int *C8_2_C8_2_FarHanging = NULL;
static HNDesc *C8_2_C8_2_FarHangingTypes = NULL;
static int ****C8_2_C8_2_FarCoupling = NULL;
static int C8_2_C8_2_NNoopposite = 0;
static int *C8_2_C8_2_Nopposite = NULL;
static int C8_2_C8_2_NNodes = 18;

TFE2DMapper *C8_2_C8_2 = new TFE2DMapper(C8_2_C8_2_Name, C8_2_C8_2_Desc,
                             C8_2_C8_2_N0, C8_2_C8_2_N1,
                             C8_2_C8_2_NPairs, (int *)C8_2_C8_2_Pairs,
                             C8_2_C8_2_NHanging, C8_2_C8_2_Hanging,
                             C8_2_C8_2_HangingTypes, C8_2_C8_2_Coupling,
                             C8_2_C8_2_NFarHanging, C8_2_C8_2_FarHanging,
                             C8_2_C8_2_FarHangingTypes, C8_2_C8_2_FarCoupling,
                             C8_2_C8_2_NNoopposite, C8_2_C8_2_Nopposite,
                             C8_2_C8_2_NNodes);
