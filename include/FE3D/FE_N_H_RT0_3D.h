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
// Raviart-Thomas element of zero-th order on hexahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_RT0_3D_NDOF = 6;

// number of dofs on the closure of each joints
static int N_H_RT0_3D_JointDOF = 1;

// which local dofs are on the joints
static int N_H_RT0_3D_J0[1] = { 0 };
static int N_H_RT0_3D_J1[1] = { 1 };
static int N_H_RT0_3D_J2[1] = { 2 };
static int N_H_RT0_3D_J3[1] = { 3 };
static int N_H_RT0_3D_J4[1] = { 4 };
static int N_H_RT0_3D_J5[1] = { 5 };

static int *N_H_RT0_3D_J[6] = { N_H_RT0_3D_J0, N_H_RT0_3D_J1,
                                N_H_RT0_3D_J2, N_H_RT0_3D_J3,
                                N_H_RT0_3D_J4, N_H_RT0_3D_J5 };

// number of inner dofs
static int N_H_RT0_3D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *N_H_RT0_3D_Inner = NULL;

// number of outer dofs
static int N_H_RT0_3D_NOuter = 6;

// array containing the numbers for the outer dofs
static int N_H_RT0_3D_Outer[6] = { 0, 1, 2, 3, 4, 5 };

static char N_H_RT0_3D_String[] = "N_H_RT0_3D";

TFEDesc3D *FE_N_H_RT0_3D_Obj=new TFEDesc3D(N_H_RT0_3D_String, N_H_RT0_3D_NDOF,
                                           N_H_RT0_3D_JointDOF, N_H_RT0_3D_J, 
                                           N_H_RT0_3D_NInner, N_H_RT0_3D_Inner,
                                           N_H_RT0_3D_NOuter, N_H_RT0_3D_Outer);
