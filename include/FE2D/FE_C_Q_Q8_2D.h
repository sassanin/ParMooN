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
   
//************************************************************
// Q8 element, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_Q_Q8_2D_NDOF = 81;

// number of dofs on the closure of the joints
static int C_Q_Q8_2D_JointDOF = 9;

// which local dofs are on the joints
static int C_Q_Q8_2D_J0[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
static int C_Q_Q8_2D_J1[9] = { 8, 17, 26, 35, 44, 53, 62, 71, 80 };
static int C_Q_Q8_2D_J2[9] = { 80, 79, 78, 77, 76, 75, 74, 73, 72 };
static int C_Q_Q8_2D_J3[9] = { 72, 63, 54, 45, 36, 27, 18, 9, 0 };

static int *C_Q_Q8_2D_J[9] = {  C_Q_Q8_2D_J0,  C_Q_Q8_2D_J1,
                              C_Q_Q8_2D_J2,  C_Q_Q8_2D_J3 };

// number of inner dofs
static int C_Q_Q8_2D_NInner = 49;

// array containing the numbers for the inner dofs
static int C_Q_Q8_2D_Inner[49] = { 10, 11, 12, 13, 14, 15, 16, 19, 20, 21,
                                   22, 23, 24, 25, 28, 29, 30, 31, 32, 33,
                                   34, 37, 38, 39, 40, 41, 42, 43, 46, 47,
                                   48, 49, 50, 51, 52, 55, 56, 57, 58, 59,
                                   60, 61, 64, 65, 66, 67, 68, 69, 70 };

// number of outer dofs
static int C_Q_Q8_2D_NOuter = 32;

// array containing the numbers for the outer dofs
static int C_Q_Q8_2D_Outer[32] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 17, 18,
                                   26, 27, 35, 36, 44, 45, 53, 54, 62, 63,
                                   71, 72, 73, 74, 75, 76, 77, 78, 79, 80 };

static char C_Q_Q8_2D_String[] = "C_Q_Q8_2D";

TFEDesc2D *FE_C_Q_Q8_2D_Obj=new TFEDesc2D(C_Q_Q8_2D_String, C_Q_Q8_2D_NDOF,
                                        C_Q_Q8_2D_JointDOF, C_Q_Q8_2D_J,
                                        C_Q_Q8_2D_NInner, C_Q_Q8_2D_Inner,
                                        C_Q_Q8_2D_NOuter, C_Q_Q8_2D_Outer);
