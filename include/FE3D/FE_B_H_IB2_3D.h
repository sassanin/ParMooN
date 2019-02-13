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
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************

// number of degrees of freedom
static int B_H_IB2_3D_NDOF = 1;

// number of dofs on the closure of the joints
static int B_H_IB2_3D_JointDOF = 0;

// which local dofs are on the joints
static int *B_H_IB2_3D_J0 = NULL;
static int *B_H_IB2_3D_J1 = NULL;
static int *B_H_IB2_3D_J2 = NULL;
static int *B_H_IB2_3D_J3 = NULL;
static int *B_H_IB2_3D_J4 = NULL;
static int *B_H_IB2_3D_J5 = NULL;

static int *B_H_IB2_3D_J[6] = { B_H_IB2_3D_J0, B_H_IB2_3D_J1,
                               B_H_IB2_3D_J2, B_H_IB2_3D_J3,
                               B_H_IB2_3D_J4, B_H_IB2_3D_J5};

// number of inner dofs
static int B_H_IB2_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int B_H_IB2_3D_Inner[1] = { 0 };

static char B_H_IB2_3D_String[] = "B_H_IB2_3D";

TFEDesc3D *FE_B_H_IB2_3D_Obj=new TFEDesc3D(B_H_IB2_3D_String, B_H_IB2_3D_NDOF, 
                                B_H_IB2_3D_JointDOF,
                                B_H_IB2_3D_J, B_H_IB2_3D_NInner, B_H_IB2_3D_Inner);
