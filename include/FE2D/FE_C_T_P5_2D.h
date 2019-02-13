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
// conforming P5 element
// ***********************************************************************

// number of degrees of freedom
static int C_T_P5_2D_NDOF = 21;

// number of dofs on the closure of the joints
static int C_T_P5_2D_JointDOF = 6;

// which local dofs are on the joints
static int C_T_P5_2D_J0[6] = {  0,  1,  2,  3,  4,  5 };
static int C_T_P5_2D_J1[6] = {  5, 10, 14, 17, 19, 20 };
static int C_T_P5_2D_J2[6] = { 20, 18, 15, 11,  6,  0 };

static int *C_T_P5_2D_J[3] = { C_T_P5_2D_J0, C_T_P5_2D_J1,
                               C_T_P5_2D_J2 };

// number of inner dofs
static int C_T_P5_2D_NInner = 6;

// array containing the numbers for the inner dofs
static int C_T_P5_2D_Inner[6] = { 7, 8, 9, 12, 13, 16 };

// number of outer dofs
static int C_T_P5_2D_NOuter = 15;

// array containing the numbers for the outer dofs
static int C_T_P5_2D_Outer[15] = { 0,  1,  2,  3,  4, 5, 6, 10, 11, 14,
                                  15, 17, 18, 19, 20 };

static char C_T_P5_2D_String[] = "C_T_P5_2D";

TFEDesc2D *FE_C_T_P5_2D_Obj=new TFEDesc2D(C_T_P5_2D_String, C_T_P5_2D_NDOF,
                              C_T_P5_2D_JointDOF, C_T_P5_2D_J,
                              C_T_P5_2D_NInner, C_T_P5_2D_Inner,
                              C_T_P5_2D_NOuter, C_T_P5_2D_Outer);
