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
   
// hanging node for conforming P2 elements in 3D

static int N_P2_3D_0_N_Nodes = 12;
static double N_P2_3D_0_Coupling[12] = { 0.25,  0.125, 0.125,
                                         0,     0,     0.125, 
                                         0,     0.125, 0,
                                         0.125, 0,     0.125 };

static int N_P2_3D_1_N_Nodes = 12;
static double N_P2_3D_1_Coupling[12] = { 0,     0,     0.125,
                                         0,     0.125, 0, 
                                         0.25,  0.125, 0.125,
                                         0,     0.125, 0.125 };

static int N_P2_3D_2_N_Nodes = 12;
static double N_P2_3D_2_Coupling[12] = { 0,     0.125, 0,
                                         0.25,  0.125, 0.125,
                                         0,     0,     0.125,
                                         0.125, 0.125, 0 };

THNDesc *HN_N_P2_3D_0_Obj = new THNDesc(N_P2_3D_0_N_Nodes,
                                        N_P2_3D_0_Coupling);

THNDesc *HN_N_P2_3D_1_Obj = new THNDesc(N_P2_3D_1_N_Nodes,
                                        N_P2_3D_1_Coupling);

THNDesc *HN_N_P2_3D_2_Obj = new THNDesc(N_P2_3D_2_N_Nodes,
                                        N_P2_3D_2_Coupling);
