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
   
// P7 and Q7 like behaviour on considered edge
static char C7_2_C7_2_Name[] = "C7_2_C7_2";
static char C7_2_C7_2_Desc[] = "conforming P7 or Q7 element";
static int C7_2_C7_2_N0 = 8;
static int C7_2_C7_2_N1 = 8;
static int C7_2_C7_2_NPairs = 8;
static int C7_2_C7_2_Pairs[][2] = { {0,15}, {1,14}, {2,13}, {3,12},
                                    {4,11}, {5,10}, {6,9}, {7,8}};
static int C7_2_C7_2_NHanging = 0;
static int *C7_2_C7_2_Hanging = NULL;
static HNDesc *C7_2_C7_2_HangingTypes = NULL;
static int **C7_2_C7_2_Coupling = NULL;
static int C7_2_C7_2_NFarHanging = 0;
static int *C7_2_C7_2_FarHanging = NULL;
static HNDesc *C7_2_C7_2_FarHangingTypes = NULL;
static int ****C7_2_C7_2_FarCoupling = NULL;
static int C7_2_C7_2_NNoopposite = 0;
static int *C7_2_C7_2_Nopposite = NULL;
static int C7_2_C7_2_NNodes = 16;

TFE2DMapper *C7_2_C7_2 = new TFE2DMapper(C7_2_C7_2_Name, C7_2_C7_2_Desc,
                             C7_2_C7_2_N0, C7_2_C7_2_N1,
                             C7_2_C7_2_NPairs, (int *)C7_2_C7_2_Pairs,
                             C7_2_C7_2_NHanging, C7_2_C7_2_Hanging,
                             C7_2_C7_2_HangingTypes, C7_2_C7_2_Coupling,
                             C7_2_C7_2_NFarHanging, C7_2_C7_2_FarHanging,
                             C7_2_C7_2_FarHangingTypes, C7_2_C7_2_FarCoupling,
                             C7_2_C7_2_NNoopposite, C7_2_C7_2_Nopposite,
                             C7_2_C7_2_NNodes);
