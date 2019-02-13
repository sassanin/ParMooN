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
   
// hanging node for conforming Q3 elements in 3D

static int C_Q3_3D_N_Nodes_1 = 4;
static double C_Q3_3D_Coupling_1[4] = { 0.3125, 0.9375, -0.3125, 0.0625 };

THNDesc *HN_C_Q3_3D_1_Obj = new THNDesc(C_Q3_3D_N_Nodes_1,
                                        C_Q3_3D_Coupling_1);

static int C_Q3_3D_N_Nodes_2 = 4;
static double C_Q3_3D_Coupling_2[4] = { -0.0625, 0.5625, 0.5625, -0.0625 };

THNDesc *HN_C_Q3_3D_2_Obj = new THNDesc(C_Q3_3D_N_Nodes_2,
                                        C_Q3_3D_Coupling_2);

static int C_Q3_3D_N_Nodes_3 = 16;
static double C_Q3_3D_Coupling_3[16] =
        { -.01953125000, .1757812500, .1757812500, -.01953125000,
          -.05859375000, .5273437500, .5273437500, -.05859375000,
           .01953125000, -.1757812500, -.1757812500, .01953125000,
        -.003906250000, .03515625000, .03515625000, -.003906250000 };

THNDesc *HN_C_Q3_3D_3_Obj = new THNDesc(C_Q3_3D_N_Nodes_3,
                                        C_Q3_3D_Coupling_3);

static int C_Q3_3D_N_Nodes_4 = 16;
static double C_Q3_3D_Coupling_4[16] =
        { .09765625000, .2929687500, -.09765625000, .01953125000,
          .2929687500, .8789062500, -.2929687500, .05859375000,
         -.09765625000, -.2929687500, .09765625000, -.01953125000,
          .01953125000, .05859375000, -.01953125000, .003906250000 };

THNDesc *HN_C_Q3_3D_4_Obj = new THNDesc(C_Q3_3D_N_Nodes_4,
                                        C_Q3_3D_Coupling_4);

static int C_Q3_3D_N_Nodes_5 = 16;
static double C_Q3_3D_Coupling_5[16] =
        { .003906250000, -.03515625000, -.03515625000, .003906250000,
         -.03515625000, .3164062500, .3164062500, -.03515625000,
         -.03515625000, .3164062500, .3164062500, -.03515625000,
          .003906250000, -.03515625000, -.03515625000, .003906250000 };

THNDesc *HN_C_Q3_3D_5_Obj = new THNDesc(C_Q3_3D_N_Nodes_5,
                                        C_Q3_3D_Coupling_5);
