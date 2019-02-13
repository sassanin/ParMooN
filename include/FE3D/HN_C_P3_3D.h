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
   
// hanging node for conforming P3 elements in 3D

static int C_P3_3D_N_Nodes_E = 4;
static double C_P3_3D_Coupling_E[4] = { 0.3125, 0.9375, -0.3125, 0.0625 };

THNDesc *HN_C_P3_3D_E_Obj = new THNDesc(C_P3_3D_N_Nodes_E,
                                        C_P3_3D_Coupling_E);

static int C_P3_3D_N_Nodes_M = 4;
static double C_P3_3D_Coupling_M[4] = { -0.0625, 0.5625, 0.5625, -0.0625 };

THNDesc *HN_C_P3_3D_M_Obj = new THNDesc(C_P3_3D_N_Nodes_M,
                                        C_P3_3D_Coupling_M);

static int C_P3_3D_N_Nodes_F = 9;
static double C_P3_3D_Coupling_F[9] = 
        { 0.5, 0.5, -0.25, 0.5, -0.25, 0.0625, -0.0625, -0.0625, 0.0625 };

THNDesc *HN_C_P3_3D_F_Obj = new THNDesc(C_P3_3D_N_Nodes_F,
                                        C_P3_3D_Coupling_F);

static int C_P3_3D_N_Nodes_G = 10;
static double C_P3_3D_Coupling_G[10] = 
        { -0.0625, 0.375, 0, 0, 0.1875, 0.75, 0, -0.1875, -0.125, 0.0625 };

THNDesc *HN_C_P3_3D_G_Obj = new THNDesc(C_P3_3D_N_Nodes_G,
                                        C_P3_3D_Coupling_G);
