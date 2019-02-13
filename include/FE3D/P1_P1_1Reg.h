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
static char P1_P1_1Reg_Name[] = "P1_P1_1Reg";
static char P1_P1_1Reg_Desc[] = "conforming P1 element, 1-regular";
static int P1_P1_1Reg_NFine = 3;
static int P1_P1_1Reg_NCoarse = 3;
static int P1_P1_1Reg_N_Pairs = 9;
static int P1_P1_1Reg_Pairs0[][2] = { {0,12}, {1,5}, {2,7}, {3,14},
                                      {4,8}, {5,9}, {6,13}, {7,11},
                                      {8,10} };
static int *P1_P1_1Reg_Pairs[1] = { (int *)P1_P1_1Reg_Pairs0 };

static int P1_P1_1Reg_NNodes = 15;

static int P1_P1_1Reg_NHanging = 3;
static int P1_P1_1Reg_Hanging[3] = { 1, 4, 2 };
static HNDesc P1_P1_1Reg_HangingTypes[3] = { HN_C_P1_3D_E, HN_C_P1_3D_E,
                                             HN_C_P1_3D_E };
static int P1_P1_1Reg_HN0[] = { 0, 3 };
static int P1_P1_1Reg_HN1[] = { 3, 6 };
static int P1_P1_1Reg_HN2[] = { 0, 6 };

static int *P1_P1_1Reg_Coupling[3] = { P1_P1_1Reg_HN0, P1_P1_1Reg_HN1,
                                       P1_P1_1Reg_HN2 };

static int P1_P1_1Reg_TwistPerm0[] = { 0, 1, 2 };
static int P1_P1_1Reg_TwistPerm1[] = { 1, 2, 0 };
static int P1_P1_1Reg_TwistPerm2[] = { 2, 0, 1 };

static int *P1_P1_1Reg_TwistPerm[3] = { P1_P1_1Reg_TwistPerm0,
                                   P1_P1_1Reg_TwistPerm1,
                                   P1_P1_1Reg_TwistPerm2 };

static int P1_P1_1Reg_NNoOpposite = 0;
static int **P1_P1_1Reg_NoOpposite = NULL;
        
TFE3DMapper1Reg *P1_P1_1Reg = new TFE3DMapper1Reg(
        P1_P1_1Reg_Name, P1_P1_1Reg_Desc,
        P1_P1_1Reg_NFine, P1_P1_1Reg_NCoarse,
        P1_P1_1Reg_N_Pairs, P1_P1_1Reg_Pairs,
        P1_P1_1Reg_NNoOpposite, P1_P1_1Reg_NoOpposite,
        P1_P1_1Reg_NHanging, P1_P1_1Reg_Hanging,
        P1_P1_1Reg_HangingTypes, P1_P1_1Reg_Coupling,
        P1_P1_1Reg_NNodes, P1_P1_1Reg_TwistPerm);
