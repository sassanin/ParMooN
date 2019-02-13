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
// @(#)JointCollection.C        4.1
//
// Class:       TJointCollection
// Purpose:     store Joint in an array
//              used by DG matrices assembler
//
// Author:      Sashikumaa Ganesan  03.11.09
//
// History:     03.11.09 Starting implementation
//
// =======================================================================  

#include <JointCollection.h>
#include <Joint.h>

/** constructor */
TJointCollection::TJointCollection(int n_joints, TJoint **joints)
{
  N_Joints = n_joints;
  Joints = joints;
}


/** destructor: delete arrays */
TJointCollection::~TJointCollection()
{
 int i;

 if(N_Joints)
 {
  for(i=0; i<N_Joints; i++)
   {
    delete Joints[i];
   }
   delete [] Joints;
 }


}
