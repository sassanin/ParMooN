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
// @(#)InterfaceJoint.C        1.2 10/18/99
// 
// Class:       TInterfaceJoint
// Purpose:     connects two cells on an interface
//
// Author:      Volker Behns  10.03.98
//
// =======================================================================

#include <InterfaceJoint.h>
#include <BoundComp2D.h>
#include <BaseCell.h>

// Constructors
TInterfaceJoint::TInterfaceJoint(TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0) : TJointEqN(neighb0)
{
  ID = InterfaceJoint;

  BoundComp = bdcomp;

  if (t_0 < t_1)
  {
    T_0 = t_0;
    T_1 = t_1;
  }
  else
  {
    T_0 = t_1;
    T_1 = t_0;
  }
}

TInterfaceJoint::TInterfaceJoint(TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0, TBaseCell *neighb1) :
                 TJointEqN(neighb0, neighb1)
{
  ID = InterfaceJoint;

  BoundComp = bdcomp;

  if (t_0 < t_1)
  {
    T_0 = t_0;
    T_1 = t_1;
  }
  else
  {
    T_0 = t_1;
    T_1 = t_0;
  }
}

// Methods
TJoint *TInterfaceJoint::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TInterfaceJoint(BoundComp, T_0 + newT_0*(T_1 - T_0),
                             T_0 + newT_1*(T_1 - T_0), Me);
}

TJoint *TInterfaceJoint::NewInst()
{
  return new TInterfaceJoint(BoundComp, T_0, T_1, NULL);
}

int TInterfaceJoint::GetXYofT(double T, double &X, double &Y)
{
  BoundComp->GetXYofT(T, X, Y);

  return 0;
}

#ifdef __2D__
/** update parameters according to the new vertex positions */
void TInterfaceJoint::UpdateParameters(TVertex *Begin, TVertex *End)
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

int TInterfaceJoint::CheckOrientation()
{
  int i, N_;
  const int *TmpEV;
  TVertex *Vert0;
  double X, Y;
  TBaseCell *Neighb2;

  if (Neighb0)
  {
    N_ = Neighb0->GetRefDesc()->GetN_OrigEdges();
    for (i=0;i<N_;i++)
      if (Neighb0->GetJoint(i) == this) break;

    Neighb0->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV);

    Vert0 = Neighb0->GetVertex(TmpEV[2*i]);
    BoundComp->GetXYofT(T_0, X, Y);

    if (ABS(Vert0->GetX() - X) > 1e-6 || ABS(Vert0->GetY() - Y) > 1e-6)
    {
      Neighb2 = Neighb1;
      Neighb1 = Neighb0;
      Neighb0 = Neighb2;
    }
  }
  else
  {
    cerr << "Erorr in InterfaceJoint: no neighbour given!" << endl;
    return -1;
  }

  return 0;
}
