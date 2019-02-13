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
// Q2 element, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_Q2_3D_NDOF = 27;

// number of dofs on the closure of the joints
static int C_H_Q2_3D_JointDOF = 9;

// which local dofs are on the joints
static int C_H_Q2_3D_J0[9] = {  0,  1,  2,  3,  4,  5,  6,  7,  8 };
static int C_H_Q2_3D_J1[9] = {  0,  9, 18,  1, 10, 19,  2, 11, 20 };
static int C_H_Q2_3D_J2[9] = {  2, 11, 20,  5, 14, 23,  8, 17, 26 };
static int C_H_Q2_3D_J3[9] = {  8, 17, 26,  7, 16, 25,  6, 15, 24 };
static int C_H_Q2_3D_J4[9] = {  0,  3,  6,  9, 12, 15, 18, 21, 24 };
static int C_H_Q2_3D_J5[9] = { 18, 21, 24, 19, 22, 25, 20, 23, 26 };

static int *C_H_Q2_3D_J[6] = { C_H_Q2_3D_J0, C_H_Q2_3D_J1,
                             C_H_Q2_3D_J2, C_H_Q2_3D_J3,
                             C_H_Q2_3D_J4, C_H_Q2_3D_J5};

#ifdef _MPI   
// number of dofs on the closure of the edges
static int C_H_Q2_3D_EdgeDOF = 3;

// which local dofs are on the joints
static int C_H_Q2_3D_E0[3] = {  0,  1,  2 };
static int C_H_Q2_3D_E1[3] = {  2,  5,  8 };
static int C_H_Q2_3D_E2[3] = {  8,  7,  6 };
static int C_H_Q2_3D_E3[3] = {  6,  3,  0 };

static int C_H_Q2_3D_E4[3] = {  0,  9,  18 };
static int C_H_Q2_3D_E5[3] = {  2,  11, 20 };
static int C_H_Q2_3D_E6[3] = {  8,  17, 26 };
static int C_H_Q2_3D_E7[3] = {  6,  15, 24 };

static int C_H_Q2_3D_E8[3] = {  18,  19,  20 };
static int C_H_Q2_3D_E9[3] = {  20,  23,  26 };
static int C_H_Q2_3D_E10[3] = { 26,  25,  24 };
static int C_H_Q2_3D_E11[3] = { 24,  21,  18 };

static int *C_H_Q2_3D_E[12] = { C_H_Q2_3D_E0, C_H_Q2_3D_E1, C_H_Q2_3D_E2, C_H_Q2_3D_E3,
                                C_H_Q2_3D_E4, C_H_Q2_3D_E5, C_H_Q2_3D_E6, C_H_Q2_3D_E7,
                                C_H_Q2_3D_E8, C_H_Q2_3D_E9, C_H_Q2_3D_E10, C_H_Q2_3D_E11};

// number of dofs on the closure of the vertices
static int C_H_Q2_3D_VertDOF = 8;

// array containing the numbers for the vertices dofs
static int C_H_Q2_3D_Vert[8] =  {0, 2, 8, 6, 18, 20, 26, 24};

#endif

// number of inner dofs
static int C_H_Q2_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_H_Q2_3D_Inner[] = { 13 };

// number of outer dof
static int C_H_Q2_3D_NOuter = 26;

// array containing the numbers for the outer dofs
static int C_H_Q2_3D_Outer[26] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,
                                    9, 10, 11, 12,     14, 15, 16, 17,
                                   18, 19, 20, 21, 22, 23, 24, 25, 26 };

static char C_H_Q2_3D_String[] = "C_H_Q2_3D";

TFEDesc3D *FE_C_H_Q2_3D_Obj=new TFEDesc3D(C_H_Q2_3D_String, C_H_Q2_3D_NDOF, 
                                C_H_Q2_3D_JointDOF,
                                C_H_Q2_3D_J,
                                C_H_Q2_3D_NInner, C_H_Q2_3D_Inner,
                                C_H_Q2_3D_NOuter, C_H_Q2_3D_Outer
#ifdef _MPI                                 
                                ,C_H_Q2_3D_EdgeDOF,  C_H_Q2_3D_E, C_H_Q2_3D_VertDOF,
                                 C_H_Q2_3D_Vert
#endif                  
                               );
