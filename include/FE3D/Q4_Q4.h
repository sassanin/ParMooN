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
static char Q4_Q4_Name[] = "Q4_Q4";
static char Q4_Q4_Desc[] = "conforming Q4 element";
static int Q4_Q4_N0 = 25;
static int Q4_Q4_N1 = 25;
static int Q4_Q4_NPairs = 25;
static int Q4_Q4_Pairs0[][2] = {
                {0,25}, {1,30}, {2,35}, {3,40}, {4,45},
                {5,26}, {6,31}, {7,36}, {8,41}, {9,46},
                {10,27}, {11,32}, {12,37}, {13,42}, {14,47},
                {15,28}, {16,33}, {17,38}, {18,43}, {19,48},
                {20,29}, {21,34}, {22,39}, {23,44}, {24,49}  };
static int Q4_Q4_Pairs1[][2] = {
                {0,29}, {1,28}, {2,27}, {3,26}, {4,25},
                {5,34}, {6,33}, {7,32}, {8,31}, {9,30},
                {10,39}, {11,38}, {12,37}, {13,36}, {14,35},
                {15,44}, {16,43}, {17,42}, {18,41}, {19,40},
                {20,49}, {21,48}, {22,47}, {23,46}, {24,45}  };
static int Q4_Q4_Pairs2[][2] = {
                {0,49}, {1,44}, {2,39}, {3,34}, {4,29},
                {5,48}, {6,43}, {7,38}, {8,33}, {9,28},
                {10,47}, {11,42}, {12,37}, {13,32}, {14,27},
                {15,46}, {16,41}, {17,36}, {18,31}, {19,26},
                {20,45}, {21,40}, {22,35}, {23,30}, {24,25}  };
static int Q4_Q4_Pairs3[][2] = {
                {0,45}, {1,46}, {2,47}, {3,48}, {4,49},
                {5,40}, {6,41}, {7,42}, {8,43}, {9,44},
                {10,35}, {11,36}, {12,37}, {13,38}, {14,39},
                {15,30}, {16,31}, {17,32}, {18,33}, {19,34},
                {20,25}, {21,26}, {22,27}, {23,28}, {24,29}  };
static int *Q4_Q4_Pairs[4] = { (int *)Q4_Q4_Pairs0, (int *)Q4_Q4_Pairs1,
                               (int *)Q4_Q4_Pairs2, (int *)Q4_Q4_Pairs3 };

static int Q4_Q4_NNodes = 50;

static int Q4_Q4_NNoOpposite = 0;
static int **Q4_Q4_NoOpposite = NULL;

TFE3DMapper *Q4_Q4 = new TFE3DMapper(Q4_Q4_Name, Q4_Q4_Desc,
                             Q4_Q4_N0, Q4_Q4_N1,
                             Q4_Q4_NPairs, Q4_Q4_Pairs,
                             Q4_Q4_NNoOpposite, Q4_Q4_NoOpposite,
                             Q4_Q4_NNodes);
