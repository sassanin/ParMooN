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
// @(#)FEDesc2D.C        1.1 10/30/98
//
// Class:       TFEDesc2D
// Purpose:     store a finite element descriptor for a 2D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc2D.h>
#include <string.h> // in order to use strdup

/** constructor, setting all data with dof on cell boundary */
TFEDesc2D::TFEDesc2D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_outerdof, int *outerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = n_outerdof;
  OuterDOF = outerdof;
}

/** constructor, setting all data with dof on cell boundary */
TFEDesc2D::TFEDesc2D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = 0;
  OuterDOF = NULL;
}

/** return joint on which the i-th local degree of freedom is   
If i is not a dof on an edge, return -1

If i is a dof on two edges (i.e. on a vertex), one of these two edges is 
returned. Don't use this function in this case.

*/
int TFEDesc2D::GetJointOfThisDOF(int localDOF) const
{
  int i,j;
  bool is_DOF_on_edge=false;
 
  
  for (i=0;i<N_OuterDOF;i++)
  {
    if(OuterDOF[i]==localDOF)
    {
      is_DOF_on_edge=true;
      break;
    }
  }
  if(!is_DOF_on_edge)
    return -1;
  //else // continue to find the edge
  i=0;
  while (true)
  {
    // this must terminate, since we already know localDOF is a dof on an edge
    for (j=0;j<N_JointDOF;j++)
    {
      if(JointDOF[i][j] == localDOF)
        return i;
    }
    i++;
  }
}
