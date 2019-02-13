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
// @(#)BoundEdge.C        1.1 10/30/98
// 
// Class:       TBoundEdge
// Purpose:     edge on a boundary component
//
// Author:      Volker Behns  02.08.97
//
// =======================================================================

#include <BoundEdge.h>

// Constructors
TBoundEdge::TBoundEdge(TBoundComp2D *bdcomp, double t_0, double t_1)
{
  ID = BoundaryEdge;

  BoundComp = bdcomp;
  T_0 = t_0;
  T_1 = t_1;
}

// Methods
int TBoundEdge::CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#ifdef __2D__
/** update parameters according to the new vertex positions */
void TBoundEdge::UpdateParameters(TVertex *Begin, TVertex *End)
{
  double x1, y1, x2, y2;
  double t1, t2;

#ifdef __2D__
  Begin->GetCoords(x1, y1);
  End->GetCoords(x2, y2);
#else
  double z1, z2;
  Begin->GetCoords(x1, y1, z1);
  End->GetCoords(x2, y2, z2);
#endif

  BoundComp->GetTofXY(x1, y1, t1);
  BoundComp->GetTofXY(x2, y2, t2);

  T_0 = t1;
  T_1 = t2;
}
#endif

#ifdef __MORTAR__

int TBoundEdge::CheckMatchingRef(TBaseCell *Me, int J_i,
                  StoreGeomMortar &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#endif

// create a new instance of this class
TJoint *TBoundEdge::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TBoundEdge(BoundComp, T_0 + newT_0*(T_1 - T_0),
                        T_0 + newT_1*(T_1 - T_0));
}

TJoint *TBoundEdge::NewInst()
{
  return new TBoundEdge(BoundComp, T_0, T_1);
}

/** return the coordinates {X,Y} of parameter value T */
int TBoundEdge::GetXYofT(double T, double &X, double &Y)
{
  BoundComp->GetXYofT(T, X, Y);

  return 0;
}
