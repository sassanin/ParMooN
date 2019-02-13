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
// Q5 element with bubbles, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_Q_UL5_2D_NDOF = 38;

// number of dofs on the closure of the joints
static int C_Q_UL5_2D_JointDOF = 6;

// which local dofs are on the joints
static int C_Q_UL5_2D_J0[6] = {  0,  1,  2,  3,  4,  5 };
static int C_Q_UL5_2D_J1[6] = {  5,  6,  7,  8,  9, 10 };
static int C_Q_UL5_2D_J2[6] = { 10, 11, 12, 13, 14, 15 };
static int C_Q_UL5_2D_J3[6] = { 15, 16, 17, 18, 19,  0 };

static int *C_Q_UL5_2D_J[4] = { C_Q_UL5_2D_J0, C_Q_UL5_2D_J1,
                                C_Q_UL5_2D_J2, C_Q_UL5_2D_J3 };

// number of inner dofs
static int C_Q_UL5_2D_NInner = 18;

// array containing the numbers for the inner dofs
static int C_Q_UL5_2D_Inner[18] = { 20, 21, 22, 23, 24, 25, 26, 27, 28,
                                    29, 30, 31, 32, 33, 34, 35, 36, 37 };

// number of outer dofs         
static int C_Q_UL5_2D_NOuter = 20;

// array containing the numbers for the inner dofs
static int C_Q_UL5_2D_Outer[20] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                    12, 13, 14, 15, 16, 17, 18, 19 };

static char C_Q_UL5_2D_String[] = "C_Q_UL5_2D";

TFEDesc2D *FE_C_Q_UL5_2D_Obj=new TFEDesc2D(C_Q_UL5_2D_String, C_Q_UL5_2D_NDOF,
                                C_Q_UL5_2D_JointDOF, C_Q_UL5_2D_J,
                                C_Q_UL5_2D_NInner, C_Q_UL5_2D_Inner,
                                C_Q_UL5_2D_NOuter, C_Q_UL5_2D_Outer);
