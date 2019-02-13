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
   
// nonconfoming Q5 like behaviour on considered edge
static char N5_2_N5_2_Name[] = "N5_2_N5_2";
static char N5_2_N5_2_Desc[] = "nonconforming Q5 element";
static int N5_2_N5_2_N0 = 5;
static int N5_2_N5_2_N1 = 5;
static int N5_2_N5_2_NPairs = 5;
static int N5_2_N5_2_Pairs[][2] = { {0,5}, {1,6}, {2,7}, {3,8}, {4,9} };
static int N5_2_N5_2_NHanging = 0;
static int *N5_2_N5_2_Hanging = NULL;
static HNDesc *N5_2_N5_2_HangingTypes = NULL;
static int **N5_2_N5_2_Coupling = NULL;
static int N5_2_N5_2_NFarHanging = 0;
static int *N5_2_N5_2_FarHanging = NULL;
static HNDesc *N5_2_N5_2_FarHangingTypes = NULL;
static int ****N5_2_N5_2_FarCoupling = NULL;
static int N5_2_N5_2_NNoopposite = 0;
static int *N5_2_N5_2_Nopposite = NULL;
static int N5_2_N5_2_NNodes = 10;

TFE2DMapper *N5_2_N5_2 = new TFE2DMapper(N5_2_N5_2_Name, N5_2_N5_2_Desc,
                             N5_2_N5_2_N0, N5_2_N5_2_N1,
                             N5_2_N5_2_NPairs, (int *)N5_2_N5_2_Pairs,
                             N5_2_N5_2_NHanging, N5_2_N5_2_Hanging,
                             N5_2_N5_2_HangingTypes, N5_2_N5_2_Coupling,
                             N5_2_N5_2_NFarHanging, N5_2_N5_2_FarHanging,
                             N5_2_N5_2_FarHangingTypes, N5_2_N5_2_FarCoupling,
                             N5_2_N5_2_NNoopposite, N5_2_N5_2_Nopposite,
                             N5_2_N5_2_NNodes);
