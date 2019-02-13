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
   
// P5 and Q5 like behaviour on considered edge
static char C5_2_C5_2_Name[] = "C5_2_C5_2";
static char C5_2_C5_2_Desc[] = "conforming P5 or Q5 element";
static int C5_2_C5_2_N0 = 6;
static int C5_2_C5_2_N1 = 6;
static int C5_2_C5_2_NPairs = 6;
static int C5_2_C5_2_Pairs[][2] = { {0,11}, {1,10}, {2,9}, {3,8},
                                    {4,7}, {5,6}};
static int C5_2_C5_2_NHanging = 0;
static int *C5_2_C5_2_Hanging = NULL;
static HNDesc *C5_2_C5_2_HangingTypes = NULL;
static int **C5_2_C5_2_Coupling = NULL;
static int C5_2_C5_2_NFarHanging = 0;
static int *C5_2_C5_2_FarHanging = NULL;
static HNDesc *C5_2_C5_2_FarHangingTypes = NULL;
static int ****C5_2_C5_2_FarCoupling = NULL;
static int C5_2_C5_2_NNoopposite = 0;
static int *C5_2_C5_2_Nopposite = NULL;
static int C5_2_C5_2_NNodes = 12;

TFE2DMapper *C5_2_C5_2 = new TFE2DMapper(C5_2_C5_2_Name, C5_2_C5_2_Desc,
                             C5_2_C5_2_N0, C5_2_C5_2_N1,
                             C5_2_C5_2_NPairs, (int *)C5_2_C5_2_Pairs,
                             C5_2_C5_2_NHanging, C5_2_C5_2_Hanging,
                             C5_2_C5_2_HangingTypes, C5_2_C5_2_Coupling,
                             C5_2_C5_2_NFarHanging, C5_2_C5_2_FarHanging,
                             C5_2_C5_2_FarHangingTypes, C5_2_C5_2_FarCoupling,
                             C5_2_C5_2_NNoopposite, C5_2_C5_2_Nopposite,
                             C5_2_C5_2_NNodes);
