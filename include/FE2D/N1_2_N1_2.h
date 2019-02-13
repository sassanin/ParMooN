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
   
// nonconfoming P1 and Q1 like behaviour on considered edge
static char N1_2_N1_2_Name[] = "N1_2_N1_2";
static char N1_2_N1_2_Desc[] = "nonconforming P1 or Q1 element";
static int N1_2_N1_2_N0 = 1;
static int N1_2_N1_2_N1 = 1;
static int N1_2_N1_2_NPairs = 1;
static int N1_2_N1_2_Pairs[][2] = { {0,1} };
static int N1_2_N1_2_NHanging = 0;
static int *N1_2_N1_2_Hanging = NULL;
static HNDesc *N1_2_N1_2_HangingTypes = NULL;
static int **N1_2_N1_2_Coupling = NULL;
static int N1_2_N1_2_NFarHanging = 0;
static int *N1_2_N1_2_FarHanging = NULL;
static HNDesc *N1_2_N1_2_FarHangingTypes = NULL;
static int ****N1_2_N1_2_FarCoupling = NULL;
static int N1_2_N1_2_NNoopposite = 0;
static int *N1_2_N1_2_Nopposite = NULL;
static int N1_2_N1_2_NNodes = 2;

TFE2DMapper *N1_2_N1_2 = new TFE2DMapper(N1_2_N1_2_Name, N1_2_N1_2_Desc,
                             N1_2_N1_2_N0, N1_2_N1_2_N1,
                             N1_2_N1_2_NPairs, (int *)N1_2_N1_2_Pairs,
                             N1_2_N1_2_NHanging, N1_2_N1_2_Hanging,
                             N1_2_N1_2_HangingTypes, N1_2_N1_2_Coupling,
                             N1_2_N1_2_NFarHanging, N1_2_N1_2_FarHanging,
                             N1_2_N1_2_FarHangingTypes, N1_2_N1_2_FarCoupling,
                             N1_2_N1_2_NNoopposite, N1_2_N1_2_Nopposite,
                             N1_2_N1_2_NNodes);
