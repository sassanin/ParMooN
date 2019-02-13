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
static char Q3_Q3_1Reg_Name[] = "Q3_Q3_1Reg";
static char Q3_Q3_1Reg_Desc[] = "conforming Q3 element, 1-regular";
static int Q3_Q3_1Reg_NFine = 16;
static int Q3_Q3_1Reg_NCoarse = 16;
static int Q3_Q3_1Reg_N_Pairs = 31;
static int Q3_Q3_1Reg_Pairs0[][2] = { {0,64}, {2,68}, {3,28}, {7,29},
        {8,65}, {10,69}, {11,30}, {12,51}, {13,55}, {14,59}, {15,31},
        {16,76}, {18,77}, {19,44}, {23,45}, {24,72}, {26,73}, {27,46},
        {31,47}, {32,79}, {34,75}, {35,60}, {39,61}, {40,78}, {42,74},
        {43,62}, {47,63}, {48,67}, {50,66}, {56,71}, {58,70} };

static int *Q3_Q3_1Reg_Pairs[1] = { (int *)Q3_Q3_1Reg_Pairs0 };

static int Q3_Q3_1Reg_NNodes = 80;

static int Q3_Q3_1Reg_NHanging = 33;
static int Q3_Q3_1Reg_Hanging[33] = {
        1, 20, 9, 22, 54, 41, 52, 33,
        4, 49, 6, 57, 25, 38, 17, 36,
        3, 11, 43, 35, 12, 14, 27, 19,
        7, 23, 39, 13,
        5, 21, 37, 53,
        15 };

static HNDesc Q3_Q3_1Reg_HangingTypes[33] = { 
        HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1,
        HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1,
        HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1,
        HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1, HN_C_Q3_3D_1,
        HN_C_Q3_3D_2, HN_C_Q3_3D_2, HN_C_Q3_3D_2, HN_C_Q3_3D_2,
        HN_C_Q3_3D_2, HN_C_Q3_3D_2, HN_C_Q3_3D_2, HN_C_Q3_3D_2,
        HN_C_Q3_3D_3, HN_C_Q3_3D_3, HN_C_Q3_3D_3, HN_C_Q3_3D_3,
        HN_C_Q3_3D_4, HN_C_Q3_3D_4, HN_C_Q3_3D_4, HN_C_Q3_3D_4,
        HN_C_Q3_3D_5 };
static int Q3_Q3_1Reg_HN0[] = { 0, 2, 24, 16 };
static int Q3_Q3_1Reg_HN1[] = { 16, 24, 2, 0 };
static int Q3_Q3_1Reg_HN2[] = { 8, 10, 26, 18 };
static int Q3_Q3_1Reg_HN3[] = { 18, 26, 10, 8 };
static int Q3_Q3_1Reg_HN4[] = { 50, 58, 42, 40 };
static int Q3_Q3_1Reg_HN5[] = { 40, 42, 58, 50 };
static int Q3_Q3_1Reg_HN6[] = { 48, 56, 34, 32 };
static int Q3_Q3_1Reg_HN7[] = { 32, 34, 56, 48 };

static int Q3_Q3_1Reg_HN8[] = { 0, 8, 50, 48 };
static int Q3_Q3_1Reg_HN9[] = { 48, 50, 8, 0 };
static int Q3_Q3_1Reg_HN10[] = { 2, 10, 58, 56 };
static int Q3_Q3_1Reg_HN11[] = { 56, 58, 10, 2 };
static int Q3_Q3_1Reg_HN12[] = { 24, 26, 42, 34 };
static int Q3_Q3_1Reg_HN13[] = { 34, 42, 26, 24 };
static int Q3_Q3_1Reg_HN14[] = { 16, 18, 40, 32 };
static int Q3_Q3_1Reg_HN15[] = { 32, 40, 18, 16 };

static int Q3_Q3_1Reg_HN16[] = { 0, 2, 24, 16 };
static int Q3_Q3_1Reg_HN17[] = { 8, 10, 26, 18 };
static int Q3_Q3_1Reg_HN18[] = { 50, 58, 42, 40 };
static int Q3_Q3_1Reg_HN19[] = { 48, 56, 34, 32 };
static int Q3_Q3_1Reg_HN20[] = { 0, 8, 50, 48 };
static int Q3_Q3_1Reg_HN21[] = { 2, 10, 58, 56 };
static int Q3_Q3_1Reg_HN22[] = { 24, 26, 42, 34 };
static int Q3_Q3_1Reg_HN23[] = { 16, 18, 40, 32 };

static int Q3_Q3_1Reg_HN24[] = { 0, 2, 24, 16, 8, 10, 26, 18,
                                 50, 58, 42, 40, 48, 56, 34, 32 };
static int Q3_Q3_1Reg_HN25[] = { 16, 18, 40, 32, 24, 26, 42, 34,
                                 2, 10, 58, 56, 0, 8,  50, 48 };
static int Q3_Q3_1Reg_HN26[] = { 32, 34, 56, 48, 40, 42, 58, 50,
                                 18, 26, 10, 8, 16, 24, 2, 0 };
static int Q3_Q3_1Reg_HN27[] = { 48, 50, 8, 0, 56, 58, 10, 2,
                                 34, 42, 26, 24, 32, 40, 18, 16 };

static int Q3_Q3_1Reg_HN28[] = { 0, 2, 24, 16, 8, 10, 26, 18,
                                 50, 58, 42, 40, 48, 56, 34, 32 };
static int Q3_Q3_1Reg_HN29[] = { 16, 18, 40, 32, 24, 26, 42, 34,
                                 2, 10, 58, 56, 0, 8,  50, 48 };
static int Q3_Q3_1Reg_HN30[] = { 32, 34, 56, 48, 40, 42, 58, 50,
                                 18, 26, 10, 8, 16, 24, 2, 0 };
static int Q3_Q3_1Reg_HN31[] = { 48, 50, 8, 0, 56, 58, 10, 2,
                                 34, 42, 26, 24, 32, 40, 18, 16 };

static int Q3_Q3_1Reg_HN32[] = { 0, 2, 24, 16, 8, 10, 26, 18,
                                 50, 58, 42, 40, 48, 56, 34, 32 };

static int *Q3_Q3_1Reg_Coupling[33] = {
     Q3_Q3_1Reg_HN0, Q3_Q3_1Reg_HN1, Q3_Q3_1Reg_HN2, Q3_Q3_1Reg_HN3,
     Q3_Q3_1Reg_HN4, Q3_Q3_1Reg_HN5, Q3_Q3_1Reg_HN6, Q3_Q3_1Reg_HN7,
     Q3_Q3_1Reg_HN8, Q3_Q3_1Reg_HN9, Q3_Q3_1Reg_HN10, Q3_Q3_1Reg_HN11,
     Q3_Q3_1Reg_HN12, Q3_Q3_1Reg_HN13, Q3_Q3_1Reg_HN14, Q3_Q3_1Reg_HN15,
     Q3_Q3_1Reg_HN16, Q3_Q3_1Reg_HN17, Q3_Q3_1Reg_HN18, Q3_Q3_1Reg_HN19,
     Q3_Q3_1Reg_HN20, Q3_Q3_1Reg_HN21, Q3_Q3_1Reg_HN22, Q3_Q3_1Reg_HN23,
     Q3_Q3_1Reg_HN24, Q3_Q3_1Reg_HN25, Q3_Q3_1Reg_HN26, Q3_Q3_1Reg_HN27,
     Q3_Q3_1Reg_HN28, Q3_Q3_1Reg_HN29, Q3_Q3_1Reg_HN30, Q3_Q3_1Reg_HN31,
     Q3_Q3_1Reg_HN32 };

static int Q3_Q3_1Reg_TwistPerm0[] = {
        0, 1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15 };
static int Q3_Q3_1Reg_TwistPerm1[] = {
        3, 7, 11, 15, 2, 6, 10, 14,
        1, 5, 9, 13, 0, 4, 8, 12 };
static int Q3_Q3_1Reg_TwistPerm2[] = {
        15, 14, 13, 12, 11, 10, 9, 8,
        7, 6, 5, 4, 3, 2, 1, 0 };
static int Q3_Q3_1Reg_TwistPerm3[] = {
        12, 8, 4, 0, 13, 9, 5, 1,
        14, 10, 6, 2, 15, 11, 7, 3 };

static int *Q3_Q3_1Reg_TwistPerm[4] = { Q3_Q3_1Reg_TwistPerm0,
                                   Q3_Q3_1Reg_TwistPerm1,
                                   Q3_Q3_1Reg_TwistPerm2,
                                   Q3_Q3_1Reg_TwistPerm3 };

static int Q3_Q3_1Reg_NNoOpposite = 0;
static int **Q3_Q3_1Reg_NoOpposite = NULL;

TFE3DMapper1Reg *Q3_Q3_1Reg = new TFE3DMapper1Reg(
        Q3_Q3_1Reg_Name, Q3_Q3_1Reg_Desc,
        Q3_Q3_1Reg_NFine, Q3_Q3_1Reg_NCoarse,
        Q3_Q3_1Reg_N_Pairs, Q3_Q3_1Reg_Pairs,
        Q3_Q3_1Reg_NNoOpposite, Q3_Q3_1Reg_NoOpposite,
        Q3_Q3_1Reg_NHanging, Q3_Q3_1Reg_Hanging,
        Q3_Q3_1Reg_HangingTypes, Q3_Q3_1Reg_Coupling,
        Q3_Q3_1Reg_NNodes, Q3_Q3_1Reg_TwistPerm);
