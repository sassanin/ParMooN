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
// @(#)InterfaceJoint.h        1.3 07/19/99
// 
// Class:       TInterfaceJoint
// Purpose:     connects two cells on an interface
//
// Author:      Volker Behns  09.03.98
//
// =======================================================================

#ifndef __INTERFACEJOINT__
#define __INTERFACEJOINT__

#include <JointEqN.h>
#include <Vertex.h>

/** connects two cells on an interface */
class TInterfaceJoint : public TJointEqN
{
  protected:
    /** boundary component to which this edge belongs */
    TBoundComp2D *BoundComp;

    /** parameter of starting point */
    double T_0;
    /** parameter of end point */
    double T_1;

  public:
    // Constructors
    /** initialize the joint with the boundary parameters and one neighbour */
    TInterfaceJoint(TBoundComp2D *bdcomp, double t_0, double t_1,
                    TBaseCell *neighb0);
    /** initialize the joint with the boundary parameters and two neighbours */
    TInterfaceJoint(TBoundComp2D *bdcomp, double t_0, double t_1,
                    TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();

    /** return start parameter T0 */
    double GetStartParameter()
    { return T_0; }

    /** return end parameter T1 */
    double GetEndParameter()
    { return T_1; }

    /** return parameters */
    void GetParameters(double &t0, double &t1)
    {
      t0 = T_0;
      t1 = T_1;
    }

#ifdef __2D__
    /** update parameters according to the new vertex positions */
    void UpdateParameters(TVertex *Begin, TVertex *End);
#endif

    /** return boundary component */
    TBoundComp2D *GetBoundComp()
    { return BoundComp; }

    /** return the coordinates {X,Y} of parameter value T */
    int GetXYofT(double T, double &X, double &Y);

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
    virtual ~TInterfaceJoint(){};
   
};

#endif
