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
// Q7 element, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_Q_Q7_2D_NDOF = 64;

// number of dofs on the closure of the joints
static int C_Q_Q7_2D_JointDOF = 8;

// which local dofs are on the joints
static int C_Q_Q7_2D_J0[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
static int C_Q_Q7_2D_J1[8] = { 7, 15, 23, 31, 39, 47, 55, 63 };
static int C_Q_Q7_2D_J2[8] = { 63, 62, 61, 60, 59, 58, 57, 56 };
static int C_Q_Q7_2D_J3[8] = { 56, 48, 40, 32, 24, 16, 8, 0 };

static int *C_Q_Q7_2D_J[8] = {  C_Q_Q7_2D_J0,  C_Q_Q7_2D_J1,
                              C_Q_Q7_2D_J2,  C_Q_Q7_2D_J3 };

// number of inner dofs
static int C_Q_Q7_2D_NInner = 36;

// array containing the numbers for the inner dofs
static int C_Q_Q7_2D_Inner[36] = { 9, 10, 11, 12, 13, 14, 17, 18, 19, 20,
                                   21, 22, 25, 26, 27, 28, 29, 30, 33, 34,
                                   35, 36, 37, 38, 41, 42, 43, 44, 45, 46,
                                   49, 50, 51, 52, 53, 54 };

// number of outer dofs
static int C_Q_Q7_2D_NOuter = 28;

// array containing the numbers for the outer dofs
static int C_Q_Q7_2D_Outer[28] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 23,
                                   24, 31, 32, 39, 40, 47, 48, 55, 56,
                                   57, 58, 59, 60, 61, 62, 63 };

static char C_Q_Q7_2D_String[] = "C_Q_Q7_2D";

TFEDesc2D *FE_C_Q_Q7_2D_Obj=new TFEDesc2D(C_Q_Q7_2D_String, C_Q_Q7_2D_NDOF,
                                        C_Q_Q7_2D_JointDOF, C_Q_Q7_2D_J,
                                        C_Q_Q7_2D_NInner, C_Q_Q7_2D_Inner,
                                        C_Q_Q7_2D_NOuter, C_Q_Q7_2D_Outer);
