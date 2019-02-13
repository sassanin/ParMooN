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
// P2 element, with face and cell bubble conforming, 3D
// 
// Author:     Sashikumaar Ganesan
//
// ***********************************************************************

// number of degrees of freedom
static int C_T_B2_3D_NDOF = 15;

// number of dofs on the closure of the joints
static int C_T_B2_3D_JointDOF = 7;

// which local dofs are on the joints
static int C_T_B2_3D_J0[7] = {  0,  1,  2,  3,  4,  5, 10 };
static int C_T_B2_3D_J1[7] = {  0,  6,  9,  1,  7,  2, 11 };
static int C_T_B2_3D_J2[7] = {  5,  4,  2,  8,  7,  9, 12 };
static int C_T_B2_3D_J3[7] = {  0,  3,  5,  6,  8,  9, 13 };

static int *C_T_B2_3D_J[4] = { C_T_B2_3D_J0, C_T_B2_3D_J1,
                             C_T_B2_3D_J2, C_T_B2_3D_J3 };

#ifdef _MPI   
// number of dofs on the closure of the edges
static int C_T_B2_3D_EdgeDOF = 3;

// which local dofs are on the joints
static int C_T_B2_3D_E0[3] = { 0, 1, 2 };
static int C_T_B2_3D_E1[3] = { 2, 4, 5 };
static int C_T_B2_3D_E2[3] = { 5, 3, 0 };
static int C_T_B2_3D_E3[3] = { 0, 6, 9 };
static int C_T_B2_3D_E4[3] = { 2, 7, 9 };
static int C_T_B2_3D_E5[3] = { 5, 8, 9 };


static int *C_T_B2_3D_E[6] = { C_T_B2_3D_E0, C_T_B2_3D_E1, C_T_B2_3D_E2, C_T_B2_3D_E3,
                               C_T_B2_3D_E4, C_T_B2_3D_E5};

// number of dofs on the closure of the vertices
static int C_T_B2_3D_VertDOF = 1;

// array containing the numbers for the vertices dofs
static int C_T_B2_3D_Vert[4] =  {0, 2, 5, 9};
#endif



// number of inner dofs
static int C_T_B2_3D_NInner = 1;

// array containing the numbers for the inner dofs (here is no inner dof)
static int C_T_B2_3D_Inner[1] ={ 14 };

static char C_T_B2_3D_String[] = "C_T_B2_3D";

TFEDesc3D *FE_C_T_B2_3D_Obj=new TFEDesc3D(C_T_B2_3D_String, C_T_B2_3D_NDOF, 
                                C_T_B2_3D_JointDOF,
                                C_T_B2_3D_J, C_T_B2_3D_NInner, C_T_B2_3D_Inner
#ifdef _MPI
                                ,C_T_B2_3D_EdgeDOF,  C_T_B2_3D_E, C_T_B2_3D_VertDOF,
                                 C_T_B2_3D_Vert
#endif
                                 );                                
                                
