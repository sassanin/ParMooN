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
// Q3 element, nonconforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_Q3_3D_NDOF = 40;

// number of dofs on the closure of the joints
static int N_H_Q3_3D_JointDOF = 6;

// which local dofs are on the joints
static int N_H_Q3_3D_J0[6] = { 0,  6, 12, 18, 24, 30 };
static int N_H_Q3_3D_J1[6] = { 1,  7, 13, 19, 25, 31 };
static int N_H_Q3_3D_J2[6] = { 2,  8, 14, 20, 26, 32 };
static int N_H_Q3_3D_J3[6] = { 3,  9, 15, 21, 27, 33 };
static int N_H_Q3_3D_J4[6] = { 4, 10, 16, 22, 28, 34 };
static int N_H_Q3_3D_J5[6] = { 5, 11, 17, 23, 29, 35 };

static int *N_H_Q3_3D_J[6] = { N_H_Q3_3D_J0, N_H_Q3_3D_J1,
                               N_H_Q3_3D_J2, N_H_Q3_3D_J3,
                               N_H_Q3_3D_J4, N_H_Q3_3D_J5};

// number of inner dofs
static int N_H_Q3_3D_NInner = 4;

// array containing the numbers for the inner dofs
static int N_H_Q3_3D_Inner[4] = { 36, 37, 38, 39 };

static char N_H_Q3_3D_String[] = "N_H_Q3_3D";

TFEDesc3D *FE_N_H_Q3_3D_Obj=new TFEDesc3D(N_H_Q3_3D_String, N_H_Q3_3D_NDOF, 
                                N_H_Q3_3D_JointDOF,
                                N_H_Q3_3D_J, N_H_Q3_3D_NInner, N_H_Q3_3D_Inner);
