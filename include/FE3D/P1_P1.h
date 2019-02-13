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
static char P1_P1_Name[] = "P1_P1";
static char P1_P1_Desc[] = "conforming P1 element";
static int P1_P1_N0 = 3;
static int P1_P1_N1 = 3;
static int P1_P1_NPairs = 3;
static int P1_P1_Pairs0[][2] = { {0,3}, {1,5}, {2,4} };
static int P1_P1_Pairs1[][2] = { {0,4}, {1,3}, {2,5} };
static int P1_P1_Pairs2[][2] = { {0,5}, {1,4}, {2,3} };
static int *P1_P1_Pairs[3] = { (int *)P1_P1_Pairs0, (int *)P1_P1_Pairs1,
                               (int *)P1_P1_Pairs2 };

static int P1_P1_NNodes = 6;

static int P1_P1_NNoOpposite = 0;
static int **P1_P1_NoOpposite = NULL;

TFE3DMapper *P1_P1 = new TFE3DMapper(P1_P1_Name, P1_P1_Desc,
                             P1_P1_N0, P1_P1_N1,
                             P1_P1_NPairs, P1_P1_Pairs,
                             P1_P1_NNoOpposite, P1_P1_NoOpposite,
                             P1_P1_NNodes);
