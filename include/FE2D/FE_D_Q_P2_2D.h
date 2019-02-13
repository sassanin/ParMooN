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
// P2 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// number of degrees of freedom
static int D_Q_P2_2D_NDOF = 6;

// number of dofs on the closure of the joints
static int D_Q_P2_2D_JointDOF = 0;

// which local dofs are on the joints
static int *D_Q_P2_2D_J0 = NULL;
static int *D_Q_P2_2D_J1 = NULL;
static int *D_Q_P2_2D_J2 = NULL;
static int *D_Q_P2_2D_J3 = NULL;

static int *D_Q_P2_2D_J[4] = { D_Q_P2_2D_J0, D_Q_P2_2D_J1,
                             D_Q_P2_2D_J2, D_Q_P2_2D_J3 };

// number of inner dofs
static int D_Q_P2_2D_NInner = 6;

// array containing the numbers for the inner dofs (here is no inner dof)
static int D_Q_P2_2D_Inner[6] = { 0, 1, 2, 3, 4, 5 };

static char D_Q_P2_2D_String[] = "D_Q_P2_2D";

TFEDesc2D *FE_D_Q_P2_2D_Obj=new TFEDesc2D(D_Q_P2_2D_String, D_Q_P2_2D_NDOF, 
                              D_Q_P2_2D_JointDOF,
                              D_Q_P2_2D_J, D_Q_P2_2D_NInner, D_Q_P2_2D_Inner);
