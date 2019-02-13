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
static char NP1_NP1_Name[] = "NP1_NP1";
static char NP1_NP1_Desc[] = "nonconforming P1 element";
static int NP1_NP1_N0 = 1;
static int NP1_NP1_NP1 = 1;
static int NP1_NP1_NPairs = 1;
static int NP1_NP1_Pairs0[][2] = { {0,1} };
static int NP1_NP1_Pairs1[][2] = { {0,1} };
static int NP1_NP1_Pairs2[][2] = { {0,1} };
static int *NP1_NP1_Pairs[3] = { (int *)NP1_NP1_Pairs0,
                                 (int *)NP1_NP1_Pairs1,
                                 (int *)NP1_NP1_Pairs2 };

static int NP1_NP1_NNodes = 2;

static int NP1_NP1_NNoOpposite = 0;
static int **NP1_NP1_NoOpposite = NULL;

TFE3DMapper *NP1_NP1 = new TFE3DMapper(NP1_NP1_Name, NP1_NP1_Desc,
                             NP1_NP1_N0, NP1_NP1_NP1,
                             NP1_NP1_NPairs, NP1_NP1_Pairs,
                             NP1_NP1_NNoOpposite, NP1_NP1_NoOpposite,
                             NP1_NP1_NNodes);
