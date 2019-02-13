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
// @(#)InterfaceJoint3D.h        1.3 07/19/99
// 
// Class:       TInterfaceJoint3D
// Purpose:     connects two cells on an interface
//
// Author:      Volker Behns  09.03.98
//              Gunar Matthies 05.04.02
//
// =======================================================================

#ifndef __INTERFACEJOINT3D__
#define __INTERFACEJOINT3D__

#include <JointEqN.h>

/** connects two cells on an interface */
class TInterfaceJoint3D : public TJointEqN
{
  protected:
    /** boundary component to which this face belongs */
    TBoundComp3D *BoundComp;

    /** parameters for vertices */
    double Param1[4];
    double Param2[4];

  public:
    // Constructors
    /** initialize the joint with the boundary parameters and one neighbour */
    TInterfaceJoint3D(TBoundComp3D *bdcomp, double *param1, double *param2,
                    TBaseCell *neigh0);
    /** initialize the joint with the boundary parameters and two neighbours */
    TInterfaceJoint3D(TBoundComp3D *bdcomp, double *param1, double *param2,
                    TBaseCell *neigh0, TBaseCell *neigh1);
    /** initialize the joint with the boundary parameters and two neighbours */
    TInterfaceJoint3D(TBoundComp3D *bdcomp, TBaseCell *neighb0);

    // Methods
    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();

    /** return boundary component */
    TBoundComp3D *GetBoundComp()
    { return BoundComp; }

    /** return both parameter vectors */
    void GetParameters(double *param1, double *param2);

    /** set both parameter vectors */
    void SetParameters(double *param1, double *param2);

    /** make sure that Neighb0 is "inside" the domain */
    int CheckOrientation();

    /** check whether Neighb0 is "inside" the domain */
    bool CheckInside(TBaseCell *Me)
    { 
      if (Neighb0 == Me)
        return true;
      else
        return false;
    }
    // Destructor
    virtual ~TInterfaceJoint3D();    
};

#endif
