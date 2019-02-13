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
// P1 Raviart-Thomas vector element, nonconforming , 2D
// History:  17.10.2011 implementation (Ulrich)
// ***********************************************************************

// number of degrees of freedom
static int N_T_RT1_2D_NDOF = 8;

// number of dofs on the closure of the joints
static int N_T_RT1_2D_JointDOF = 2;

// which local dofs are on the joints
static int N_T_RT1_2D_J0[2] = { 0, 1 };
static int N_T_RT1_2D_J1[2] = { 2, 3 };
static int N_T_RT1_2D_J2[2] = { 4, 5 };

 
static int *N_T_RT1_2D_J[3] = { N_T_RT1_2D_J0,
                                N_T_RT1_2D_J1,
                                N_T_RT1_2D_J2
                              };

// number of inner dofs
static int N_T_RT1_2D_NInner = 2;

// array containing the numbers for the inner dofs
static int N_T_RT1_2D_Inner[2] = {6,7};

// number of outer dofs (dofs on edges)
static int N_T_RT1_2D_NOuter = 6;

// array containing the numbers for the outer dofs
static int N_T_RT1_2D_Outer[6] = { 0, 1, 2, 3, 4, 5};

static char N_T_RT1_2D_String[] = "N_T_RT1_2D";

TFEDesc2D *FE_N_T_RT1_2D_Obj=new TFEDesc2D(N_T_RT1_2D_String, N_T_RT1_2D_NDOF,
                                        N_T_RT1_2D_JointDOF, N_T_RT1_2D_J,
                                        N_T_RT1_2D_NInner, N_T_RT1_2D_Inner,
                                        N_T_RT1_2D_NOuter, N_T_RT1_2D_Outer);
