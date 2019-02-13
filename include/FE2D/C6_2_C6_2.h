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
   
// P6 and Q6 like behaviour on considered edge
static char C6_2_C6_2_Name[] = "C6_2_C6_2";
static char C6_2_C6_2_Desc[] = "conforming P6 or Q6 element";
static int C6_2_C6_2_N0 = 7;
static int C6_2_C6_2_N1 = 7;
static int C6_2_C6_2_NPairs = 7;
static int C6_2_C6_2_Pairs[][2] = { {0,13}, {1,12}, {2,11}, {3,10},
                                    {4,9}, {5,8}, {6,7}};
static int C6_2_C6_2_NHanging = 0;
static int *C6_2_C6_2_Hanging = NULL;
static HNDesc *C6_2_C6_2_HangingTypes = NULL;
static int **C6_2_C6_2_Coupling = NULL;
static int C6_2_C6_2_NFarHanging = 0;
static int *C6_2_C6_2_FarHanging = NULL;
static HNDesc *C6_2_C6_2_FarHangingTypes = NULL;
static int ****C6_2_C6_2_FarCoupling = NULL;
static int C6_2_C6_2_NNoopposite = 0;
static int *C6_2_C6_2_Nopposite = NULL;
static int C6_2_C6_2_NNodes = 14;

TFE2DMapper *C6_2_C6_2 = new TFE2DMapper(C6_2_C6_2_Name, C6_2_C6_2_Desc,
                             C6_2_C6_2_N0, C6_2_C6_2_N1,
                             C6_2_C6_2_NPairs, (int *)C6_2_C6_2_Pairs,
                             C6_2_C6_2_NHanging, C6_2_C6_2_Hanging,
                             C6_2_C6_2_HangingTypes, C6_2_C6_2_Coupling,
                             C6_2_C6_2_NFarHanging, C6_2_C6_2_FarHanging,
                             C6_2_C6_2_FarHangingTypes, C6_2_C6_2_FarCoupling,
                             C6_2_C6_2_NNoopposite, C6_2_C6_2_Nopposite,
                             C6_2_C6_2_NNodes);
