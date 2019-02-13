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
static char Q1_Q1_Name[] = "Q1_Q1";
static char Q1_Q1_Desc[] = "conforming Q1 element";
static int Q1_Q1_N0 = 4;
static int Q1_Q1_N1 = 4;
static int Q1_Q1_NPairs = 4;
static int Q1_Q1_Pairs0[][2] = { {0,4}, {1,6}, {2,5}, {3,7} };
static int Q1_Q1_Pairs1[][2] = { {0,5}, {1,4}, {2,7}, {3,6} };
static int Q1_Q1_Pairs2[][2] = { {0,7}, {1,5}, {2,6}, {3,4} };
static int Q1_Q1_Pairs3[][2] = { {0,6}, {1,7}, {2,4}, {3,5} };
static int *Q1_Q1_Pairs[4] = { (int *)Q1_Q1_Pairs0, (int *)Q1_Q1_Pairs1,
                               (int *)Q1_Q1_Pairs2, (int *)Q1_Q1_Pairs3 };

static int Q1_Q1_NNodes = 8;

static int Q1_Q1_NNoOpposite = 0;
static int **Q1_Q1_NoOpposite = NULL;

TFE3DMapper *Q1_Q1 = new TFE3DMapper(Q1_Q1_Name, Q1_Q1_Desc,
                             Q1_Q1_N0, Q1_Q1_N1,
                             Q1_Q1_NPairs, Q1_Q1_Pairs,
                             Q1_Q1_NNoOpposite, Q1_Q1_NoOpposite,
                             Q1_Q1_NNodes);
