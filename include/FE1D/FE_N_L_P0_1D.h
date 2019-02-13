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
// P0 element, nonconforming, 1D
// ***********************************************************************

// number of degrees of freedom
static int N_L_P0_1D_NDOF = 1;

// number of dofs on the closure of the joints
static int N_L_P0_1D_JointDOF = 0;

// which local dofs are on the joints
static int *N_L_P0_1D_J[2] = {0, 0};

// number of inner dofs
static int N_L_P0_1D_NInner = 1;

// array containing the numbers for the inner dofs
static int N_L_P0_1D_Inner[1] = { 0 };

static char N_L_P0_1D_String[] = "N_L_P0_1D";

TFEDesc1D *FE_N_L_P0_1D_Obj=new TFEDesc1D(N_L_P0_1D_String, N_L_P0_1D_NDOF, 
                              N_L_P0_1D_JointDOF, N_L_P0_1D_J, 
                              N_L_P0_1D_NInner, N_L_P0_1D_Inner);
