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
   
// mapper for tetrahedral faces on both sides
static char NP2_NP2_Name[] = "NP2_NP2";
static char NP2_NP2_Desc[] = "nonconforming P2 element";
static int NP2_NP2_N0 = 3;
static int NP2_NP2_NP2 = 3;
static int NP2_NP2_NPairs = 3;
static int NP2_NP2_Pairs0[][2] = { {0,3}, {1,5}, {2,4} };
static int NP2_NP2_Pairs1[][2] = { {0,4}, {1,3}, {2,5} };
static int NP2_NP2_Pairs2[][2] = { {0,5}, {1,4}, {2,3} };
static int *NP2_NP2_Pairs[3] = { (int *)NP2_NP2_Pairs0,
                                 (int *)NP2_NP2_Pairs1,
                                 (int *)NP2_NP2_Pairs2 };

static int NP2_NP2_NNodes = 6;

static int NP2_NP2_NNoOpposite = 0;
static int **NP2_NP2_NoOpposite = NULL;

TFE3DMapper *NP2_NP2 = new TFE3DMapper(NP2_NP2_Name, NP2_NP2_Desc,
                             NP2_NP2_N0, NP2_NP2_NP2,
                             NP2_NP2_NPairs, NP2_NP2_Pairs,
                             NP2_NP2_NNoOpposite, NP2_NP2_NoOpposite,
                             NP2_NP2_NNodes);
