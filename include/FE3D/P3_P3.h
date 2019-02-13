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
   
// mapper for tetrahedron faces on both sides
static char P3_P3_Name[] = "P3_P3";
static char P3_P3_Desc[] = "conforming P3 element";
static int P3_P3_N0 = 10;
static int P3_P3_N1 = 10;
static int P3_P3_NPairs = 10;
static int P3_P3_Pairs0[][2] = { {0,10}, {1,14}, {2,17}, {3,19},
                                 {4,11}, {5,15}, {6,18}, {7,12},
                                 {8,16}, {9,13} };
static int P3_P3_Pairs1[][2] = { {0,13}, {1,12}, {2,11}, {3,10},
                                 {4,16}, {5,15}, {6,14}, {7,18},
                                 {8,17}, {9,19} };
static int P3_P3_Pairs2[][2] = { {0,19}, {1,18}, {2,16}, {3,13},
                                 {4,17}, {5,15}, {6,12}, {7,14},
                                 {8,11}, {9,10} };
static int *P3_P3_Pairs[3] = { (int *)P3_P3_Pairs0, (int *)P3_P3_Pairs1,
                               (int *)P3_P3_Pairs2 };

static int P3_P3_NNodes = 20;

static int P3_P3_NNoOpposite = 0;
static int **P3_P3_NoOpposite = NULL;

TFE3DMapper *P3_P3 = new TFE3DMapper(P3_P3_Name, P3_P3_Desc,
                             P3_P3_N0, P3_P3_N1,
                             P3_P3_NPairs, P3_P3_Pairs,
                             P3_P3_NNoOpposite, P3_P3_NoOpposite,
                             P3_P3_NNodes);
