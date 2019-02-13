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
static char N1_N1_Name[] = "N1_N1";
static char N1_N1_Desc[] = "nonconforming Q1Rot element";
static int N1_N1_N0 = 1;
static int N1_N1_N1 = 1;
static int N1_N1_NPairs = 1;
static int N1_N1_Pairs0[][2] = { {0,1} };
static int N1_N1_Pairs1[][2] = { {0,1} };
static int N1_N1_Pairs2[][2] = { {0,1} };
static int N1_N1_Pairs3[][2] = { {0,1} };
static int *N1_N1_Pairs[4] = { (int *)N1_N1_Pairs0, (int *)N1_N1_Pairs1,
                               (int *)N1_N1_Pairs2, (int *)N1_N1_Pairs3 };

static int N1_N1_NNodes = 2;

static int N1_N1_NNoOpposite = 0;
static int **N1_N1_NoOpposite = NULL;

TFE3DMapper *N1_N1 = new TFE3DMapper(N1_N1_Name, N1_N1_Desc,
                             N1_N1_N0, N1_N1_N1,
                             N1_N1_NPairs, N1_N1_Pairs,
                             N1_N1_NNoOpposite, N1_N1_NoOpposite,
                             N1_N1_NNodes);
