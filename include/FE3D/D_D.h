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
static char D_D_Name[] = "D_D";
static char D_D_Desc[] = "discontinuous elements";
static int D_D_N0 = 0;
static int D_D_D = 0;
static int D_D_NPairs = 0;
static int *D_D_Pairs0 = NULL;
static int *D_D_Pairs1 = NULL;
static int *D_D_Pairs2 = NULL;
static int *D_D_Pairs3 = NULL;
static int *D_D_Pairs[4] = { (int *)D_D_Pairs0, (int *)D_D_Pairs1,
                               (int *)D_D_Pairs2, (int *)D_D_Pairs3 };

static int D_D_NNodes = 0;

static int D_D_NNoOpposite = 0;
static int **D_D_NoOpposite = NULL;

TFE3DMapper *D_D = new TFE3DMapper(D_D_Name, D_D_Desc,
                             D_D_N0, D_D_D,
                             D_D_NPairs, D_D_Pairs,
                             D_D_NNoOpposite, D_D_NoOpposite,
                             D_D_NNodes);
