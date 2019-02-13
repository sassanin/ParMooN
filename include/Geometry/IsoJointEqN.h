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
// @(#)IsoJointEqN.h        1.2 10/18/99
// 
// Class:       TIsoJointEqN
// Purpose:     connects two cells
//              has additional vertices for isoparametric reference 
//              transformation
//
// Author:      Gunar Matthies 06.08.1999
//
// History:     start of implementation (Gunar Matthies 06.08.99)
//
// =======================================================================

#ifndef __ISOJOINTEQN__
#define __ISOJOINTEQN__

#include <JointEqN.h>
#include <Vertex.h>

/** connects two cells */
class TIsoJointEqN : public TJointEqN
{
  protected:
    /** number of additional vertices */
    int N_Vertices;

    /** array of all additional vertices */
    TVertex **Vertices;

  public:
    // Constructors
    /** constructor with one initial neighbour */
    TIsoJointEqN(TBaseCell *neighb0);
    /** constructor with two initial neighbours */
    TIsoJointEqN(TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me)
    { return new TIsoJointEqN(Me); }
    virtual TJoint *NewInst()
    { return new TIsoJointEqN(NULL); }

    /** return number of additional vertices */
    int GetN_Vertices()
    { return N_Vertices; }

    TVertex **GetVertices()
    { return Vertices; }

    void SetVertices(int n_vertices, TVertex **vertices);
};

#endif
