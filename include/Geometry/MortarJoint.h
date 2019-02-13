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
// @(#)MortarJoint.h        1.2 08/27/99
// 
// Class:       TMortarJoint
// Purpose:     indicates a mortar joint
//
// Author:      Volker Behns  09.03.98
//
// =======================================================================

#ifndef __MORTARJOINT__
#define __MORTARJOINT__

#include <Joint.h>

/** indicates a mortar joint */
class TMortarJoint : public TJoint
{
  protected:
    /** cell number in mortar collection */
    int MEdgeInColl;

  public:
    // Constructors
    TMortarJoint();

    // Methods
    /** checking for matching joints */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                     struct StoreGeom &Tmp);

    #ifdef __MORTAR__
      /** check the refinement pattern on both sides for matching,
          special version for moratr cells */
      virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                    StoreGeomMortar &Tmp);
    #endif

    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me)
    { return new TMortarJoint(); }
    virtual TJoint *NewInst()
    { return new TMortarJoint(); }

    /** set MEdgeInColl */
    void SetMEdgeInColl(int medgeincoll)
    { MEdgeInColl = medgeincoll; }

    /** return MEdgeInColl */
    int GetMEdgeInColl()
    { return MEdgeInColl; }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const
    { return false; }
};

#endif
