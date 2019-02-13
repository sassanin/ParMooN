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
// Q3 element with bubbles, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_UL3_3D_NDOF = 67;

// number of dofs on the closure of the joints
static int C_H_UL3_3D_JointDOF = 16;

// which local dofs are on the joints
static int C_H_UL3_3D_J0[16] = {  0,  1,  2,  3,  4,  5,  6,  7,
                               8,  9, 10, 11, 12, 13, 14, 15 };
static int C_H_UL3_3D_J1[16] = {  0, 16, 28, 40,  1, 17, 29, 41,
                               2, 18, 30, 42,  3, 19, 31, 43 };
static int C_H_UL3_3D_J2[16] = {  3, 19, 31, 43,  7, 20, 32, 47,
                              11, 21, 33, 51, 15, 22, 34, 55 };
static int C_H_UL3_3D_J3[16] = { 15, 22, 34, 55, 14, 23, 35, 54,
                              13, 24, 36, 53, 12, 25, 37, 52 };
static int C_H_UL3_3D_J4[16] = {  0,  4,  8, 12, 16, 27, 26, 25,
                              28, 39, 38, 37, 40, 44, 48, 52 };
static int C_H_UL3_3D_J5[16] = { 40, 44, 48, 52, 41, 45, 49, 53,
                              42, 46, 50, 54, 43, 47, 51, 55 };

static int *C_H_UL3_3D_J[6] = { C_H_UL3_3D_J0, C_H_UL3_3D_J1,
                             C_H_UL3_3D_J2, C_H_UL3_3D_J3,
                             C_H_UL3_3D_J4, C_H_UL3_3D_J5};

// number of inner dofs
static int C_H_UL3_3D_NInner = 11;

// array containing the numbers for the inner dofs
static int C_H_UL3_3D_Inner[11] = { 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66 };

// number of outer dof
static int C_H_UL3_3D_NOuter = 56;

// array containing the numbers for the outer dofs
static int C_H_UL3_3D_Outer[56] = {   0,  1,  2,  3,  4,  5,  6,  7,  8,
                                     9, 10, 11, 12, 13, 14, 15, 16, 17,
									18, 19, 20, 21, 22, 23, 24, 25, 26,
									27, 28, 29, 30, 31, 32, 33, 34, 35,
									36, 37, 38, 39, 40, 41, 42, 43, 44,
									45, 46, 47, 48, 49, 50, 51, 52, 53,
									54, 55 };

static char C_H_UL3_3D_String[] = "C_H_UL3_3D";

TFEDesc3D *FE_C_H_UL3_3D_Obj=new TFEDesc3D(C_H_UL3_3D_String, C_H_UL3_3D_NDOF, 
                                C_H_UL3_3D_JointDOF,
                                C_H_UL3_3D_J,
								C_H_UL3_3D_NInner, C_H_UL3_3D_Inner,
                                C_H_UL3_3D_NOuter, C_H_UL3_3D_Outer);
