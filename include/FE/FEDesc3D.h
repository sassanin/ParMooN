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
// History:     Added info of edges and vertices for parallel 
//               (Sashikumaar Ganesan 08.09.2010)
// =======================================================================

#ifndef __FEDESC3D__
#define __FEDESC3D__

#include <Enumerations.h>

/** store a finite element descriptor for a 3D element */
class TFEDesc3D
{
  protected:
    /** description for the object */
    char *Description;

    /** number of degrees of freedom */
    int N_DOF;

    /** number of degrees of freedom on closure of each joint */
    int N_JointDOF;

    /** local numbers of all degrees of freedom on the joints */
    int **JointDOF;

    /** number of inner degrees of freedom */
    int N_InnerDOF;

    /** local numbers of all inner degrees of freedom */
    int *InnerDOF;

    /** number of degrees of freedom on cell boundary */
    int N_OuterDOF;

    /** local numbers of all degrees of freedom on cell boundary */
    int *OuterDOF;
 
#ifdef _MPI      
     /** MPI data set flag */
    int EdgeVertData_Filled;   
    
    /** number of degrees of freedom on closure of each edge */
    int N_EdgeDOF;

    /** local numbers of all degrees of freedom on the edge */
    int **EdgeDOF;

    /** number of degrees of freedom on each vertex */
    int N_VertDOF;

    /** local numbers of all degrees of freedom on the vertices */
    int *VertDOF;
#endif

  public:
    /** constructor, setting all data */
    TFEDesc3D(char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof);
      
    /** constructor, setting all data with dof on cell boundary */
    TFEDesc3D(char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof,
              int n_outerdof, int *outerdof);

#ifdef _MPI 
    /** constructor, setting all data including edge and vertices data*/
    TFEDesc3D(char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof,
              int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF );
      
    /** constructor, setting all data with dof on cell boundary including edge and vertices data* */
    TFEDesc3D(char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof,
              int n_outerdof, int *outerdof,
              int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF );


    /** return number of degrees of freedom per closure of each edge */
    int GetN_EdgeDOF() const
      { return N_EdgeDOF; }
      
    /** return local numbers of degrees of freedom on each edge */
    int **GetEdgeDOF() const
      { return EdgeDOF; }

    /** return local numbers of degrees of freedom on edge i */
    int *GetEdgeDOF(int i) const
      { return EdgeDOF[i]; } 

    /** return number of degrees of freedom per closure of each vertex */
    int GetN_VertDOF() const
      { return N_VertDOF; }

    /** return local numbers of degrees of freedom on vertex i */
    int  GetVertDOF(int i) const
      { return VertDOF[i]; }

    int IsEdgeVertData_Filled() const
      {return EdgeVertData_Filled; }
      
#endif        
      
    /** return description */
    char *GetDescription() const
      { return Description; }

    /** return number of degrees of freedom */
    int GetN_DOF() const
      { return N_DOF; }

    /** return number of degrees of freedom per closure of each joint */
    int GetN_JointDOF() const
      { return N_JointDOF; }

    /** return number of inner degrees of freedom */
    int GetN_InnerDOF() const
      { return N_InnerDOF; }

    /** return local numbers of inner degrees of freedom */
    int *GetInnerDOF() const
      { return InnerDOF; }

    /** return number of degrees of freedom on cell boundary */
    int GetN_OuterDOF() const
      { return N_OuterDOF; }

    /** return local numbers of degrees of freedom on cell boundary */
    int *GetOuterDOF() const
      { return OuterDOF; }

    /** return total number and local numbers of degrees of freedom
        on cell boundary */ 
    void GetOuterDOF(int &n_outerdof, int* &outerdof) const
      { n_outerdof = N_OuterDOF; outerdof = OuterDOF; }

    /** return local numbers of degrees of freedom on each joint */
    int **GetJointDOF() const
      { return JointDOF; }

    /** return local numbers of degrees of freedom on joint i */
    int *GetJointDOF(int i) const
      { return JointDOF[i]; }

    /** return face on which the i-th local degree of freedom is   
    If i is not a dof on a face, return -1.
    If i is a dof on two faces (e.g. on a vertex), one of these two faces is 
    returned. Don't use this function in this case.
    */
    int GetJointOfThisDOF(int localDOF) const;

};

#endif
