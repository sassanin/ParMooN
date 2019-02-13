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
   
// =======================================================================
// @(#)FEDesc1D.C        1.1 10/30/98
//
// Class:       TFEDesc1D
// Purpose:     store a finite element descriptor for a 1D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc1D.h>
#include <string.h>

/** constructor, setting all data */
TFEDesc1D::TFEDesc1D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
}
