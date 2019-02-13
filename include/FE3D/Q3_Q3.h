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
static char Q3_Q3_Name[] = "Q3_Q3";
static char Q3_Q3_Desc[] = "conforming Q3 element";
static int Q3_Q3_N0 = 16;
static int Q3_Q3_N1 = 16;
static int Q3_Q3_NPairs = 16;
static int Q3_Q3_Pairs0[][2] = { {0,16}, {1,20}, {2,24}, {3,28},
                                 {4,17}, {5,21}, {6,25}, {7,29},
                                 {8,18}, {9,22}, {10,26}, {11,30},
                                 {12,19}, {13,23}, {14,27}, {15,31} };
static int Q3_Q3_Pairs1[][2] = { {0,19}, {1,18}, {2,17}, {3,16},
                                 {4,23}, {5,22}, {6,21}, {7,20},
                                 {8,27}, {9,26}, {10,25}, {11,24},
                                 {12,31}, {13,30}, {14,29}, {15,28} };
static int Q3_Q3_Pairs2[][2] = { {0,31}, {1,27}, {2,23}, {3,19},
                                 {4,30}, {5,26}, {6,22}, {7,18},
                                 {8,29}, {9,25}, {10,21}, {11,17},
                                 {12,28}, {13,24}, {14,20}, {15,16} };
static int Q3_Q3_Pairs3[][2] = { {0,28}, {1,29}, {2,30}, {3,31},
                                 {4,24}, {5,25}, {6,26}, {7,27},
                                 {8,20}, {9,21}, {10,22}, {11,23},
                                 {12,16}, {13,17}, {14,18}, {15,19} };
static int *Q3_Q3_Pairs[4] = { (int *)Q3_Q3_Pairs0, (int *)Q3_Q3_Pairs1,
                               (int *)Q3_Q3_Pairs2, (int *)Q3_Q3_Pairs3 };

static int Q3_Q3_NNodes = 32;

static int Q3_Q3_NNoOpposite = 0;
static int **Q3_Q3_NoOpposite = NULL;

TFE3DMapper *Q3_Q3 = new TFE3DMapper(Q3_Q3_Name, Q3_Q3_Desc,
                             Q3_Q3_N0, Q3_Q3_N1,
                             Q3_Q3_NPairs, Q3_Q3_Pairs,
                             Q3_Q3_NNoOpposite, Q3_Q3_NoOpposite,
                             Q3_Q3_NNodes);
