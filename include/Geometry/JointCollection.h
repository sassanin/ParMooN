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
// @(#)JointCollection.h        4.1
//
// Class:       TJointCollection
// Purpose:     store Joint in an array
//              used by DG matrices assembler
//
// Author:      Sashikumaar Ganesan  03.11.09
//
// History:     03.11.09 Starting implementation
//
// ======================================================================= 

#ifndef __JOINTCOLLECTION__
#define __JOINTCOLLECTION__

#include <Joint.h>

/** store joints in an array, used by DG matrices assembler */
class TJointCollection
{
  protected:
    /** number of joints stored */
    int N_Joints;

    /** array containing the pointers to the cells */
    TJoint **Joints;

  public:
    /** constructor */
    TJointCollection(int n_joints, TJoint **joints);

    /** return number of joints */
    int GetN_Joints()
    { return N_Joints; }

    /** return joint with index i in joint-array */
    TJoint *GetJoint(int i)
    { return Joints[i]; }

    /** destructor: delete arrays */
    ~TJointCollection();
};

#endif
