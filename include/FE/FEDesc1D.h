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
// @(#)FEDesc1D.h        1.1 10/30/98
//
// Class:       TFEDesc1D
// Purpose:     store a finite element descriptor for a 1D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#ifndef __FEDESC1D__
#define __FEDESC1D__

#include <Enumerations.h>

/** store a finite element descriptor for a 1D element */
class TFEDesc1D
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

  public:
    /** constructor, setting all data */
    TFEDesc1D(char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof);

    /** return description */
    char *GetDescription()
      { return Description; }

    /** return number of degrees of freedom */
    int GetN_DOF()
      { return N_DOF; }

    /** return number of degrees of freedom per closure of each joint */
    int GetN_JointDOF()
      { return N_JointDOF; }

    /** return number of inner degrees of freedom */
    int GetN_InnerDOF()
      { return N_InnerDOF; }

    /** return local numbers of inner degrees of freedom */
    int *GetInnerDOF()
      { return InnerDOF; }

    /** return local numbers of degrees of freedom on each joint */
    int **GetJointDOF()
      { return JointDOF; }

    /** return local numbers of degrees of freedom on joint i */
    int *GetJointDOF(int i)
      { return JointDOF[i]; }

};

#endif
