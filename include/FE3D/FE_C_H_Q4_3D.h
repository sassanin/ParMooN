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
// Q4 element, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_Q4_3D_NDOF = 125;

// number of dofs on the closure of the joints
static int C_H_Q4_3D_JointDOF = 25;

// which local dofs are on the joints
static int C_H_Q4_3D_J0[25] = {   0,   1,   2,   3,   4,
                                  5,   6,   7,   8,   9,
                                 10,  11,  12,  13,  14,
                                 15,  16,  17,  18,  19,
                                 20,  21,  22,  23,  24 };
static int C_H_Q4_3D_J1[25] = {   0,  25,  50,  75, 100,
                                  1,  26,  51,  76, 101,
                                  2,  27,  52,  77, 102,
                                  3,  28,  53,  78, 103,
                                  4,  29,  54,  79, 104 };
static int C_H_Q4_3D_J2[25] = {   4,  29,  54,  79, 104,
                                  9,  34,  59,  84, 109,
                                 14,  39,  64,  89, 114,
                                 19,  44,  69,  94, 119,
                                 24,  49,  74,  99, 124 };
static int C_H_Q4_3D_J3[25] = {  24,  49,  74,  99, 124,
                                 23,  48,  73,  98, 123,
                                 22,  47,  72,  97, 122,
                                 21,  46,  71,  96, 121,
                                 20,  45,  70,  95, 120 };
static int C_H_Q4_3D_J4[25] = {   0,   5,  10,  15,  20,
                                 25,  30,  35,  40,  45, 
                                 50,  55,  60,  65,  70, 
                                 75,  80,  85,  90,  95, 
                                100, 105, 110, 115, 120 }; 
static int C_H_Q4_3D_J5[25] = { 100, 105, 110, 115, 120,
                                101, 106, 111, 116, 121,
                                102, 107, 112, 117, 122,
                                103, 108, 113, 118, 123,
                                104, 109, 114, 119, 124 };

static int *C_H_Q4_3D_J[6] = { C_H_Q4_3D_J0, C_H_Q4_3D_J1,
                             C_H_Q4_3D_J2, C_H_Q4_3D_J3,
                             C_H_Q4_3D_J4, C_H_Q4_3D_J5};

// number of inner dofs
static int C_H_Q4_3D_NInner = 27;

// array containing the numbers for the inner dofs (here is no inner dof)
static int C_H_Q4_3D_Inner[27] = { 31, 32, 33, 36, 37, 38, 41, 42, 43,
                                   56, 57, 58, 61, 62, 63, 66, 67, 68,
                                   81, 82, 83, 86, 87, 88, 91, 92, 93 };

static char C_H_Q4_3D_String[] = "C_H_Q4_3D";

TFEDesc3D *FE_C_H_Q4_3D_Obj=new TFEDesc3D(
                C_H_Q4_3D_String, C_H_Q4_3D_NDOF, 
                C_H_Q4_3D_JointDOF,
                C_H_Q4_3D_J, C_H_Q4_3D_NInner, C_H_Q4_3D_Inner);
