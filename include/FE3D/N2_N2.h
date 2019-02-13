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
static char N2_N2_Name[] = "N2_N2";
static char N2_N2_Desc[] = "nonconforming element of order 2";
static int N2_N2_N0 = 3;
static int N2_N2_N1 = 3;
static int N2_N2_NPairs = 3;
static int N2_N2_Pairs0[][2] = { {0,3}, {1,5}, {2,4} };
static int N2_N2_Pairs1[][2] = { {0,3}, {1,4}, {2,5} };
static int N2_N2_Pairs2[][2] = { {0,3}, {1,5}, {2,4} };
static int N2_N2_Pairs3[][2] = { {0,3}, {1,4}, {2,5} };
static int *N2_N2_Pairs[4] = { (int *)N2_N2_Pairs0, (int *)N2_N2_Pairs1,
                               (int *)N2_N2_Pairs2, (int *)N2_N2_Pairs3 };

static int N2_N2_NNodes = 6;

static int N2_N2_NNoOpposite = 0;
static int **N2_N2_NoOpposite = NULL;

TFE3DMapper *N2_N2 = new TFE3DMapper(N2_N2_Name, N2_N2_Desc,
                             N2_N2_N0, N2_N2_N1,
                             N2_N2_NPairs, N2_N2_Pairs,
                             N2_N2_NNoOpposite, N2_N2_NoOpposite,
                             N2_N2_NNodes);
