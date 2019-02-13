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
// Q4 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_Q4_2D_NDOF = 25;

// number of dofs on the closure of the joints
static int C_Q_Q4_2D_JointDOF = 5;

// which local dofs are on the joints
static int C_Q_Q4_2D_J0[5] = {  0,  1,  2,  3,  4 };
static int C_Q_Q4_2D_J1[5] = {  4,  9, 14, 19, 24 };
static int C_Q_Q4_2D_J2[5] = { 24, 23, 22, 21, 20 };
static int C_Q_Q4_2D_J3[5] = { 20, 15, 10,  5,  0 };

static int *C_Q_Q4_2D_J[4] = { C_Q_Q4_2D_J0, C_Q_Q4_2D_J1, 
                             C_Q_Q4_2D_J2, C_Q_Q4_2D_J3 };

// number of inner dofs
static int C_Q_Q4_2D_NInner = 9;

// array containing the numbers for the inner dofs
static int C_Q_Q4_2D_Inner[9] = { 6, 7, 8, 11, 12, 13, 16, 17, 18 };

// number of outer dofs
static int C_Q_Q4_2D_NOuter = 16;

// array containing the numbers for the outer dofs
static int C_Q_Q4_2D_Outer[16] = { 0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19,
                                   20, 21, 22, 23, 24 };

static char C_Q_Q4_2D_String[] = "C_Q_Q4_2D";

TFEDesc2D *FE_C_Q_Q4_2D_Obj=new TFEDesc2D(C_Q_Q4_2D_String, C_Q_Q4_2D_NDOF,
                                C_Q_Q4_2D_JointDOF, C_Q_Q4_2D_J,
                                C_Q_Q4_2D_NInner, C_Q_Q4_2D_Inner,
                                C_Q_Q4_2D_NOuter, C_Q_Q4_2D_Outer);
