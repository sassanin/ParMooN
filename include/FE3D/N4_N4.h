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
static char N4_N4_Name[] = "N4_N4";
static char N4_N4_Desc[] = "nonconforming element of order 4";
static int N4_N4_N0 = 10;
static int N4_N4_N1 = 10;
static int N4_N4_NPairs = 10;
static int N4_N4_Pairs0[][2] = { {0,10}, {1,12}, {2,11}, {3,15}, {4,14},
                                 {5,13}, {6,19}, {7,18}, {8,17}, {9,16} };
static int N4_N4_Pairs1[][2] = { {0,10}, {1,11}, {2,12}, {3,13}, {4,14},
                                 {5,15}, {6,16}, {7,17}, {8,18}, {9,19} };
static int N4_N4_Pairs2[][2] = { {0,10}, {1,12}, {2,11}, {3,15}, {4,14},
                                 {5,13}, {6,19}, {7,18}, {8,17}, {9,16} };
static int N4_N4_Pairs3[][2] = { {0,10}, {1,11}, {2,12}, {3,13}, {4,14},
                                 {5,15}, {6,16}, {7,17}, {8,18}, {9,19} };
static int *N4_N4_Pairs[4] = { (int *)N4_N4_Pairs0, (int *)N4_N4_Pairs1,
                               (int *)N4_N4_Pairs2, (int *)N4_N4_Pairs3 };

static int N4_N4_NNodes = 20;

static int N4_N4_NNoOpposite = 0;
static int **N4_N4_NoOpposite = NULL;

TFE3DMapper *N4_N4 = new TFE3DMapper(N4_N4_Name, N4_N4_Desc,
                             N4_N4_N0, N4_N4_N1,
                             N4_N4_NPairs, N4_N4_Pairs,
                             N4_N4_NNoOpposite, N4_N4_NoOpposite,
                             N4_N4_NNodes);
