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
// @(#)PeriodicJoint.h        1.1 04/21/99
// 
// Class:       TPeriodicJoint
// Purpose:     connects two cells on different parts of the domain
//              where we have periodic boundary conditions
//
// Author:      Volker Behns  16.04.99
//
// =======================================================================

#ifndef __PERIODICJOINT__
#define __PERIODICJOINT__

#include <JointEqN.h>

/** connects two cells with periodic boundary conditions */
class TPeriodicJoint : public TJointEqN
{
  public:
    // Constructors
    /** initialize the joint with one neighbour */
    TPeriodicJoint(TBaseCell *neighb0);

    /** initialize the joint with two neighbours */
    TPeriodicJoint(TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();
};

#endif
