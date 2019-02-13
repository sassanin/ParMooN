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
   
// mapper for hexahedron faces on both sides
static char N3_N3_Name[] = "N3_N3";
static char N3_N3_Desc[] = "nonconforming element of order 3";
static int N3_N3_N0 = 6;
static int N3_N3_N1 = 6;
static int N3_N3_NPairs = 6;
static int N3_N3_Pairs0[][2] = { {0,6}, {1,8}, {2,7}, {3,11}, {4,10}, {5,9} };
static int N3_N3_Pairs1[][2] = { {0,6}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11} };
static int N3_N3_Pairs2[][2] = { {0,6}, {1,8}, {2,7}, {3,11}, {4,10}, {5,9} };
static int N3_N3_Pairs3[][2] = { {0,6}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11} };
static int *N3_N3_Pairs[4] = { (int *)N3_N3_Pairs0, (int *)N3_N3_Pairs1,
                               (int *)N3_N3_Pairs2, (int *)N3_N3_Pairs3 };

static int N3_N3_NNodes = 12;

static int N3_N3_NNoOpposite = 0;
static int **N3_N3_NoOpposite = NULL;

TFE3DMapper *N3_N3 = new TFE3DMapper(N3_N3_Name, N3_N3_Desc,
                             N3_N3_N0, N3_N3_N1,
                             N3_N3_NPairs, N3_N3_Pairs,
                             N3_N3_NNoOpposite, N3_N3_NoOpposite,
                             N3_N3_NNodes);
