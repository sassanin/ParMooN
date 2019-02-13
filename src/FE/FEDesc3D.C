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
// %W% %G%
//
// Class:       TFEDesc3D
// Purpose:     store a finite element descriptor for a 3D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc3D.h>
#include <string.h>
#include <iostream>

/** constructor, setting all data */
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
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
  
#ifdef _MPI    
  N_EdgeDOF = 0;
  EdgeDOF = NULL;
  N_VertDOF = 0;
  VertDOF = NULL;
  EdgeVertData_Filled = 0;  
#endif

}

/** constructor, setting all data */
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
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
  
#ifdef _MPI    
  N_EdgeDOF = 0;
  EdgeDOF = NULL;
  N_VertDOF = 0;
  VertDOF = NULL;
  EdgeVertData_Filled = 0;
#endif
  
}


#ifdef _MPI 
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF )
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = 0;
  OuterDOF = NULL;

  N_EdgeDOF = n_edgeDOF;
  EdgeDOF = edgeDOF;
  N_VertDOF = n_vertDOF;
  VertDOF = vertDOF; 
  EdgeVertData_Filled = 1;
}


TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_outerdof, int *outerdof,
                     int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF )
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = n_outerdof;
  OuterDOF = outerdof;

  N_EdgeDOF = n_edgeDOF;
  EdgeDOF = edgeDOF;
  N_VertDOF = n_vertDOF;
  VertDOF = vertDOF; 
  EdgeVertData_Filled = 1;  
}

#endif // _MPI

/** return face on which the i-th local degree of freedom is   
If i is not a dof on a face, return -1

If i is a dof on two faces (e.g. on a vertex), one of these two faces is 
returned. Don't use this function in this case.

*/
int TFEDesc3D::GetJointOfThisDOF(int localDOF) const
{
  bool is_DOF_on_edge = false;
 
  for (int i = 0; i < N_OuterDOF; i++)
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
  int i = 0;
  while (true)
  {
    // this must terminate, since we already know localDOF is a dof on an edge
    for (int j = 0; j < N_JointDOF; j++)
    {
      if(JointDOF[i][j] == localDOF)
        return i;
    }
    i++;
  }
}
