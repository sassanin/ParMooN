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
   
// ***********************************************************************
// Brezzi-Douglas-Duran element of third order on tetrahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_T_BDDF3_3D_NDOF = 60;

// number of dofs on the closure of each joints
static int N_T_BDDF3_3D_JointDOF = 10;

// which local dofs are on the joints
static int N_T_BDDF3_3D_J0[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
static int N_T_BDDF3_3D_J1[10] = {10,11,12,13,14,15,16,17,18,19 };
static int N_T_BDDF3_3D_J2[10] = {20,21,22,23,24,25,26,27,28,29 };
static int N_T_BDDF3_3D_J3[10] = {30,31,32,33,34,35,36,37,38,39 };

static int *N_T_BDDF3_3D_J[4] = { N_T_BDDF3_3D_J0, N_T_BDDF3_3D_J1,
                                N_T_BDDF3_3D_J2, N_T_BDDF3_3D_J3 };

// number of inner dofs
static int N_T_BDDF3_3D_NInner = 20;

// array containing the numbers for the inner dofs (here is no inner dof)
static int N_T_BDDF3_3D_Inner[20] = { 40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59 };

// number of outer dofs
static int N_T_BDDF3_3D_NOuter = 40;

// array containing the numbers for the outer dofs
static int N_T_BDDF3_3D_Outer[40] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 };

static char N_T_BDDF3_3D_String[] = "N_T_BDDF3_3D";

TFEDesc3D *FE_N_T_BDDF3_3D_Obj=new TFEDesc3D(N_T_BDDF3_3D_String, N_T_BDDF3_3D_NDOF,
                                           N_T_BDDF3_3D_JointDOF, N_T_BDDF3_3D_J,
                                           N_T_BDDF3_3D_NInner, N_T_BDDF3_3D_Inner,
                                           N_T_BDDF3_3D_NOuter, N_T_BDDF3_3D_Outer);
