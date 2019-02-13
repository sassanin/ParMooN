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
// Q1 element with bubble, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_UL1_2D_NDOF = 5;

// number of dofs on the closure of the joints
static int C_Q_UL1_2D_JointDOF = 2;

// which local dofs are on the joints
static int C_Q_UL1_2D_J0[2] = { 0, 1 };
static int C_Q_UL1_2D_J1[2] = { 1, 2 };
static int C_Q_UL1_2D_J2[2] = { 2, 3 };
static int C_Q_UL1_2D_J3[2] = { 3, 0 };

static int *C_Q_UL1_2D_J[4] = { C_Q_UL1_2D_J0, C_Q_UL1_2D_J1,
                                C_Q_UL1_2D_J2, C_Q_UL1_2D_J3 };

// number of inner dofs
static int C_Q_UL1_2D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_Q_UL1_2D_Inner[1] = { 4 };

// number of outer dofs
static int C_Q_UL1_2D_NOuter = 4;

// array containing the numbers for the outer dofs
static int C_Q_UL1_2D_Outer[4] = { 0, 1, 2, 3 };

static char C_Q_UL1_2D_String[] = "C_Q_UL1_2D";

TFEDesc2D *FE_C_Q_UL1_2D_Obj=new TFEDesc2D(C_Q_UL1_2D_String, C_Q_UL1_2D_NDOF,
                                C_Q_UL1_2D_JointDOF, C_Q_UL1_2D_J,
                                C_Q_UL1_2D_NInner, C_Q_UL1_2D_Inner,
                                C_Q_UL1_2D_NOuter, C_Q_UL1_2D_Outer);
