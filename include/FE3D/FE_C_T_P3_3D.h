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
// P3 element, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P3_3D_NDOF = 20;

// number of dofs on the closure of the joints
static int C_T_P3_3D_JointDOF = 10;

// which local dofs are on the joints
static int C_T_P3_3D_J0[10] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9 };
static int C_T_P3_3D_J1[10] = { 0, 10, 16, 19,  1, 11, 17,  2, 12,  3 };
static int C_T_P3_3D_J2[10] = { 9,  8,  6,  3, 15, 14, 12, 18, 17, 19 };
static int C_T_P3_3D_J3[10] = { 0,  4,  7,  9, 10, 13, 15, 16, 18, 19 };

static int *C_T_P3_3D_J[4] = { C_T_P3_3D_J0, C_T_P3_3D_J1,
                             C_T_P3_3D_J2, C_T_P3_3D_J3 };

// number of inner dofs
static int C_T_P3_3D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *C_T_P3_3D_Inner = NULL;

static char C_T_P3_3D_String[] = "C_T_P3_3D";

TFEDesc3D *FE_C_T_P3_3D_Obj=new TFEDesc3D(C_T_P3_3D_String, C_T_P3_3D_NDOF, 
                                C_T_P3_3D_JointDOF,
                                C_T_P3_3D_J, C_T_P3_3D_NInner, C_T_P3_3D_Inner);
