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
// P4 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P4_2D_NDOF = 18;

// number of dofs on the closure of the joints
static int N_T_P4_2D_JointDOF = 4;

// which local dofs are on the joints
static int N_T_P4_2D_J0[4] = { 0, 3, 6,  9 };
static int N_T_P4_2D_J1[4] = { 1, 4, 7, 10 };
static int N_T_P4_2D_J2[4] = { 2, 5, 8, 11 };
 
static int *N_T_P4_2D_J[3] = { N_T_P4_2D_J0, N_T_P4_2D_J1, N_T_P4_2D_J2 };

// number of inner dofs
static int N_T_P4_2D_NInner = 6;

// array containing the numbers for the inner dofs
static int N_T_P4_2D_Inner[6] = { 12, 13, 14, 15, 16, 17 };

// number of outer dofs
static int N_T_P4_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int N_T_P4_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

static char N_T_P4_2D_String[] = "N_T_P4_2D";

TFEDesc2D *FE_N_T_P4_2D_Obj=new TFEDesc2D(N_T_P4_2D_String, N_T_P4_2D_NDOF,
                                        N_T_P4_2D_JointDOF, N_T_P4_2D_J,
                                        N_T_P4_2D_NInner, N_T_P4_2D_Inner,
                                        N_T_P4_2D_NOuter, N_T_P4_2D_Outer);
