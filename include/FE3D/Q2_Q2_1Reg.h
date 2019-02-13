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
static char Q2_Q2_1Reg_Name[] = "Q2_Q2_1Reg";
static char Q2_Q2_1Reg_Desc[] = "conforming Q2 element, 1-regular";
static int Q2_Q2_1Reg_NFine = 9;
static int Q2_Q2_1Reg_NCoarse = 9;
static int Q2_Q2_1Reg_N_Pairs = 20;
static int Q2_Q2_1Reg_Pairs0[][2] = { {0,36}, {2,39}, {5,16}, {6,37},
                                      {7,32}, {8,40}, {9,42}, {11,43},
                                      {14,25}, {2,15}, {8,17}, {18,44},
                                      {20,41}, {23,34}, {11,24}, {17,26},
                                      {27,38}, {6,29}, {20,33}, {26,35} };

static int *Q2_Q2_1Reg_Pairs[1] = { (int *)Q2_Q2_1Reg_Pairs0 };

static int Q2_Q2_1Reg_NNodes = 45;

static int Q2_Q2_1Reg_NHanging = 16;
static int Q2_Q2_1Reg_Hanging[16] = { 1, 3, 10, 12, 19, 21, 28, 30,
                                      5, 14, 23, 7,
                                      4, 13, 22, 31 };
static HNDesc Q2_Q2_1Reg_HangingTypes[16] = { 
        HN_C_Q2_3D_E, HN_C_Q2_3D_E, HN_C_Q2_3D_E, HN_C_Q2_3D_E,
        HN_C_Q2_3D_E, HN_C_Q2_3D_E, HN_C_Q2_3D_E, HN_C_Q2_3D_E,
        HN_C_Q2_3D_E, HN_C_Q2_3D_E, HN_C_Q2_3D_E, HN_C_Q2_3D_E,
        HN_C_Q2_3D_F, HN_C_Q2_3D_F, HN_C_Q2_3D_F, HN_C_Q2_3D_F };
static int Q2_Q2_1Reg_HN0[] = { 0, 2, 9 };
static int Q2_Q2_1Reg_HN1[] = { 0, 6, 27 };
static int Q2_Q2_1Reg_HN2[] = { 9, 11, 18 };
static int Q2_Q2_1Reg_HN3[] = { 9, 2, 0 };
static int Q2_Q2_1Reg_HN4[] = { 18, 20, 27 };
static int Q2_Q2_1Reg_HN5[] = { 18, 11, 9 };
static int Q2_Q2_1Reg_HN6[] = { 27, 29, 0 };
static int Q2_Q2_1Reg_HN7[] = { 27, 20, 18 };

static int Q2_Q2_1Reg_HN8[] = { 2, 8, 20 };
static int Q2_Q2_1Reg_HN9[] = { 11, 8, 6 };
static int Q2_Q2_1Reg_HN10[] = { 20, 8, 2 };
static int Q2_Q2_1Reg_HN11[] = { 6, 8, 11 };

static int Q2_Q2_1Reg_HN12[] = { 0, 6, 27, 2, 8, 20, 9, 11, 18 };
static int Q2_Q2_1Reg_HN13[] = { 9, 2, 0, 11, 8, 6, 18, 20, 27 };
static int Q2_Q2_1Reg_HN14[] = { 18, 11, 9, 20, 8, 2, 27, 6, 0 };
static int Q2_Q2_1Reg_HN15[] = { 27, 20, 18, 6, 8, 11, 0, 2, 9 };

static int *Q2_Q2_1Reg_Coupling[16] = {
     Q2_Q2_1Reg_HN0, Q2_Q2_1Reg_HN1, Q2_Q2_1Reg_HN2, Q2_Q2_1Reg_HN3,
     Q2_Q2_1Reg_HN4, Q2_Q2_1Reg_HN5, Q2_Q2_1Reg_HN6, Q2_Q2_1Reg_HN7,
     Q2_Q2_1Reg_HN8, Q2_Q2_1Reg_HN9, Q2_Q2_1Reg_HN10, Q2_Q2_1Reg_HN11,
     Q2_Q2_1Reg_HN12, Q2_Q2_1Reg_HN13, Q2_Q2_1Reg_HN14, Q2_Q2_1Reg_HN15 };

static int Q2_Q2_1Reg_TwistPerm0[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
static int Q2_Q2_1Reg_TwistPerm1[] = { 2, 5, 8, 1, 4, 7, 0, 3, 6 };
static int Q2_Q2_1Reg_TwistPerm2[] = { 8, 7, 6, 5, 4, 3, 2, 1, 0 };
static int Q2_Q2_1Reg_TwistPerm3[] = { 6, 3, 0, 7, 4, 1, 8, 5, 2 };

static int *Q2_Q2_1Reg_TwistPerm[4] = { Q2_Q2_1Reg_TwistPerm0,
                                   Q2_Q2_1Reg_TwistPerm1,
                                   Q2_Q2_1Reg_TwistPerm2,
                                   Q2_Q2_1Reg_TwistPerm3 };

static int Q2_Q2_1Reg_NNoOpposite = 12;
static int Q2_Q2_1Reg_NoOpposite1[12] = { 4, 13, 22, 31,
        1, 3, 10, 12, 19, 21, 28, 30 };
static int *Q2_Q2_1Reg_NoOpposite[] = { Q2_Q2_1Reg_NoOpposite1 };

TFE3DMapper1Reg *Q2_Q2_1Reg = new TFE3DMapper1Reg(
        Q2_Q2_1Reg_Name, Q2_Q2_1Reg_Desc,
        Q2_Q2_1Reg_NFine, Q2_Q2_1Reg_NCoarse,
        Q2_Q2_1Reg_N_Pairs, Q2_Q2_1Reg_Pairs,
        Q2_Q2_1Reg_NNoOpposite, Q2_Q2_1Reg_NoOpposite,
        Q2_Q2_1Reg_NHanging, Q2_Q2_1Reg_Hanging,
        Q2_Q2_1Reg_HangingTypes, Q2_Q2_1Reg_Coupling,
        Q2_Q2_1Reg_NNodes, Q2_Q2_1Reg_TwistPerm);
