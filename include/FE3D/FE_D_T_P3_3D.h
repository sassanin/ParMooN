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
// P3 element, discontinuous, 3D
// 
// Author:     Markus Wolff
//
// ***********************************************************************

// number of degrees of freedom
static int D_T_P3_3D_NDOF = 20;

// number of dofs on the closure of the joints
static int D_T_P3_3D_JointDOF = 0;

// which local dofs are on the joints
static int *D_T_P3_3D_J0 = NULL;
static int *D_T_P3_3D_J1 = NULL;
static int *D_T_P3_3D_J2 = NULL;
static int *D_T_P3_3D_J3 = NULL;

static int *D_T_P3_3D_J[4] = { D_T_P3_3D_J0, D_T_P3_3D_J1,
                               D_T_P3_3D_J2, D_T_P3_3D_J3 };

// number of inner dofs
static int D_T_P3_3D_NInner = 20;

// array containing the numbers for the inner dofs 
static int D_T_P3_3D_Inner[] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 };

static char D_T_P3_3D_String[] = "D_T_P3_3D";

TFEDesc3D *FE_D_T_P3_3D_Obj=new TFEDesc3D(D_T_P3_3D_String, D_T_P3_3D_NDOF,
                                D_T_P3_3D_JointDOF,
                                D_T_P3_3D_J, D_T_P3_3D_NInner,
                                D_T_P3_3D_Inner);
