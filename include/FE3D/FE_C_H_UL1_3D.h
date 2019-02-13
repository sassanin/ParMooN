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
// Q1 element with bubble, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_UL1_3D_NDOF = 9;

// number of dofs on the closure of the joints
static int C_H_UL1_3D_JointDOF = 4;

// which local dofs are on the joints
static int C_H_UL1_3D_J0[4] = { 0, 1, 2, 3 };
static int C_H_UL1_3D_J1[4] = { 0, 4, 1, 5 };
static int C_H_UL1_3D_J2[4] = { 1, 5, 3, 7 };
static int C_H_UL1_3D_J3[4] = { 3, 7, 2, 6 };
static int C_H_UL1_3D_J4[4] = { 0, 2, 4, 6 };
static int C_H_UL1_3D_J5[4] = { 4, 6, 5, 7 };

static int *C_H_UL1_3D_J[6] = { C_H_UL1_3D_J0, C_H_UL1_3D_J1,
                             C_H_UL1_3D_J2, C_H_UL1_3D_J3,
                             C_H_UL1_3D_J4, C_H_UL1_3D_J5};

// number of inner dofs
static int C_H_UL1_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_H_UL1_3D_Inner[1] = { 8};

// number of outer dof
static int C_H_UL1_3D_NOuter = 8;

// array containing the numbers for the outer dofs
static int C_H_UL1_3D_Outer[8] = {  0,  1,  2,  3,  4,  5,  6,  7};

static char C_H_UL1_3D_String[] = "C_H_UL1_3D";

TFEDesc3D *FE_C_H_UL1_3D_Obj=new TFEDesc3D(C_H_UL1_3D_String, C_H_UL1_3D_NDOF, 
                                C_H_UL1_3D_JointDOF,
								C_H_UL1_3D_J,
								C_H_UL1_3D_NInner, C_H_UL1_3D_Inner,
								C_H_UL1_3D_NOuter, C_H_UL1_3D_Outer);
