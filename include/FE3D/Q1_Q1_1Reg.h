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
   
/*
    TFE3DMapper1Reg(char *name, char *description, int nfine, int ncoarse,
              int n_pairs, int **pairs,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation);
*/
static char Q1_Q1_1Reg_Name[] = "Q1_Q1_1Reg";
static char Q1_Q1_1Reg_Desc[] = "conforming Q1 element, 1-regular";
static int Q1_Q1_1Reg_NFine = 4;
static int Q1_Q1_1Reg_NCoarse = 4;
static int Q1_Q1_1Reg_N_Pairs = 11;
static int Q1_Q1_1Reg_Pairs0[][2] = { {0,16}, { 1,6}, {2,13}, {3,7},
                                      {4,18}, {5,10}, {7,11}, {8,19},
                                      {9,14}, {11,15}, {12,17} };
static int *Q1_Q1_1Reg_Pairs[1] = { (int *)Q1_Q1_1Reg_Pairs0 };

static int Q1_Q1_1Reg_NNodes = 20;

static int Q1_Q1_1Reg_NHanging = 5;
static int Q1_Q1_1Reg_Hanging[5] = { 1, 5, 9, 2, 3 };
static HNDesc Q1_Q1_1Reg_HangingTypes[5] = { HN_C_Q1_3D_E, HN_C_Q1_3D_E,
                                             HN_C_Q1_3D_E, HN_C_Q1_3D_E,
                                             HN_C_Q1_3D_F };
static int Q1_Q1_1Reg_HN0[] = { 0, 4 };
static int Q1_Q1_1Reg_HN1[] = { 4, 8 };
static int Q1_Q1_1Reg_HN2[] = { 8, 12 };
static int Q1_Q1_1Reg_HN3[] = { 0, 12 };
static int Q1_Q1_1Reg_HN4[] = { 0, 4, 8, 12 };

static int *Q1_Q1_1Reg_Coupling[5] = { Q1_Q1_1Reg_HN0, Q1_Q1_1Reg_HN1,
                                       Q1_Q1_1Reg_HN2, Q1_Q1_1Reg_HN3,
                                       Q1_Q1_1Reg_HN4 };

static int Q1_Q1_1Reg_TwistPerm0[] = { 0, 1, 2, 3 };
static int Q1_Q1_1Reg_TwistPerm1[] = { 1, 3, 0, 2 };
static int Q1_Q1_1Reg_TwistPerm2[] = { 3, 2, 1, 0 };
static int Q1_Q1_1Reg_TwistPerm3[] = { 2, 0, 3, 1 };

static int *Q1_Q1_1Reg_TwistPerm[4] = { Q1_Q1_1Reg_TwistPerm0,
                                   Q1_Q1_1Reg_TwistPerm1,
                                   Q1_Q1_1Reg_TwistPerm2,
                                   Q1_Q1_1Reg_TwistPerm3 };

static int Q1_Q1_1Reg_NNoOpposite = 0;
static int **Q1_Q1_1Reg_NoOpposite = NULL;

TFE3DMapper1Reg *Q1_Q1_1Reg = new TFE3DMapper1Reg(
        Q1_Q1_1Reg_Name, Q1_Q1_1Reg_Desc,
        Q1_Q1_1Reg_NFine, Q1_Q1_1Reg_NCoarse,
        Q1_Q1_1Reg_N_Pairs, Q1_Q1_1Reg_Pairs,
        Q1_Q1_1Reg_NNoOpposite, Q1_Q1_1Reg_NoOpposite,
        Q1_Q1_1Reg_NHanging, Q1_Q1_1Reg_Hanging,
        Q1_Q1_1Reg_HangingTypes, Q1_Q1_1Reg_Coupling,
        Q1_Q1_1Reg_NNodes, Q1_Q1_1Reg_TwistPerm);
