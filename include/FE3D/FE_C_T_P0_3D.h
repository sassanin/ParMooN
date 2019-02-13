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
// P0 element, discontinuous, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P0_3D_NDOF = 1;

// number of dofs on the closure of the joints
static int C_T_P0_3D_JointDOF = 0;

// which local dofs are on the joints
static int *C_T_P0_3D_J0 = NULL;
static int *C_T_P0_3D_J1 = NULL;
static int *C_T_P0_3D_J2 = NULL;
static int *C_T_P0_3D_J3 = NULL;

static int *C_T_P0_3D_J[4] = { C_T_P0_3D_J0, C_T_P0_3D_J1,
                               C_T_P0_3D_J2, C_T_P0_3D_J3 };

// number of inner dofs
static int C_T_P0_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_T_P0_3D_Inner[1] = { 0 };



#ifdef _MPI   
// number of dofs on the closure of the edges
static int  C_T_P0_3D_EdgeDOF = 0;

// which local dofs are on the joints
static int *C_T_P0_3D_E0 = NULL;
static int *C_T_P0_3D_E1 = NULL;
static int *C_T_P0_3D_E2 = NULL;
static int *C_T_P0_3D_E3 = NULL;
static int *C_T_P0_3D_E4 = NULL;
static int *C_T_P0_3D_E5 = NULL;

static int *C_T_P0_3D_E[6] = { C_T_P0_3D_E0, C_T_P0_3D_E1, C_T_P0_3D_E2, C_T_P0_3D_E3,
                               C_T_P0_3D_E4, C_T_P0_3D_E5};
			       
// number of dofs on the closure of the vertices
static int C_T_P0_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *C_T_P0_3D_Vert =  NULL;
#endif


static char C_T_P0_3D_String[] = "C_T_P0_3D";

TFEDesc3D *FE_C_T_P0_3D_Obj=new TFEDesc3D(C_T_P0_3D_String, C_T_P0_3D_NDOF, 
                                C_T_P0_3D_JointDOF,
                                C_T_P0_3D_J, C_T_P0_3D_NInner, C_T_P0_3D_Inner
#ifdef _MPI
                                ,C_T_P0_3D_EdgeDOF,  C_T_P0_3D_E, C_T_P0_3D_VertDOF,
                                 C_T_P0_3D_Vert
#endif
                                 );
