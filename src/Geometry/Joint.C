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
// @(#)Joint.C        1.7 11/18/99
// 
// Class:       TJoint
// Purpose:     superclass for all joints
//
// Author:      Volker Behns  23.09.98
//
// =======================================================================

#include <BaseCell.h>
#include <Database.h>
#include <Joint.h>
#include <stdlib.h>

// Constructors
TJoint::TJoint()
{
  ID = Joint;

  Neighb0 = NULL;
  Neighb1 = NULL;

  NeibSubDomainLocalJointNo = -1;
#ifdef __3D__
  MapType = -1;
#endif
}

// Methods
int TJoint::SetNeighbour(TBaseCell *Neighb)
{
  switch (this->GetType())
  {
    case JointEqN:
    case MortarBaseJoint:
    case InterfaceJoint:
    case IsoInterfaceJoint:
#ifdef __3D__
    case InterfaceJoint3D:
    case IsoInterfaceJoint3D:
#endif
    case PeriodicJoint:
    case InnerInterfaceJoint:
      if (Neighb0)
        Neighb1 = Neighb;
      else
        Neighb0 = Neighb;

      return 0;

    default:
      return -1;
  }
}

TBaseCell *TJoint::GetNeighbour(TBaseCell *Me) const
{
  if (Neighb0 == Me)
    return Neighb1;
  else
    return Neighb0;
}

int TJoint::SetNeighbour(int i, TBaseCell *Neighb)
{
  switch (this->GetType())
  {
    case JointEqN:
    case MortarBaseJoint:
    case InterfaceJoint:
    case IsoInterfaceJoint:
#ifdef __3D__
    case InterfaceJoint3D:
    case IsoInterfaceJoint3D:
#endif
    case PeriodicJoint:
      if (i)
        Neighb1 = Neighb;
      else
        Neighb0 = Neighb;

      return 0;

    default:
      return -1;
  }
}

TBaseCell *TJoint::GetNeighbour(int i) const
{
  if (i)
    return Neighb1;
  else
    return Neighb0;
}

void TJoint::Delete(TBaseCell *Me)
{
  if (Neighb0 == Me)
    Neighb0 = NULL;
  else
    Neighb1 = NULL;
}

#ifdef __3D__
void TJoint::SetMapType()
{
  int N_, LocJoint0, LocJoint1, MaxLen, aux;
  const int *TmpFV, *TmpLen;
  TVertex *Vert;

  if (Neighb0 && Neighb1)
  {
    N_ = Neighb0->GetN_Faces();
    for (LocJoint0=0;LocJoint0<N_;LocJoint0++)
      if (Neighb0->GetJoint(LocJoint0) == this) break;

    N_ = Neighb1->GetN_Faces();
    for (LocJoint1=0;LocJoint1<N_;LocJoint1++)
      if (Neighb1->GetJoint(LocJoint1) == this) break;

    Neighb0->GetRefDesc()->GetShapeDesc()->
             GetFaceVertex(TmpFV, TmpLen, MaxLen);

    Vert = Neighb0->GetVertex(TmpFV[LocJoint0 * MaxLen]);

    Neighb1->GetRefDesc()->GetShapeDesc()->
             GetFaceVertex(TmpFV, TmpLen, MaxLen);

    N_ = TmpLen[LocJoint1];
    aux = LocJoint1 * MaxLen;
    for (MapType=0;MapType<N_;MapType++)
      if (Neighb1->GetVertex(TmpFV[aux + MapType]) == Vert) break;

    if (MapType == N_)
    {
      /*
      int i;
      N_ = Neighb0->GetN_Vertices();
      for (i=0;i<N_;i++)
        cout << " test 0:" << i << Neighb0->GetVertex(i) << "  " << (int)
                Neighb0->GetVertex(i) << endl;

      N_ = Neighb0->GetN_Vertices();
      for (i=0;i<N_;i++)
        cout << " test 1:" << i << Neighb1->GetVertex(i) << "  " << (int)
                Neighb1->GetVertex(i) << endl;
      */

      cerr << "Error in SetMapType: could not find vertex" << endl;
      exit (-1);
      return ;
    }
  }
}

void TJoint::GetMapperRef(const int *&MapVerts, const int *&MapFaces)
{
  if (MapType != -1)
    switch (Neighb0->GetType())
    {
      case Tetrahedron: TDatabase::MapperDB[MapTriReg0 + MapType]->
                          GetMapperRef(MapVerts, MapFaces);
                        break;

      case Brick:
      case Hexahedron: TDatabase::MapperDB[MapQuadReg0 + MapType]->
                         GetMapperRef(MapVerts, MapFaces);
                       break;
      default:
      break;	
    }
  else
  {
    cerr << "Error in GetMapperRef: wrong MapType " << MapType << endl;
    exit (-1);
  }
}

void TJoint::GetMapperOrig(const int *&MapVerts, const int *&MapEdges)
{
  if (MapType != -1)
    switch (Neighb0->GetType())
    {
      case Tetrahedron: TDatabase::MapperDB[MapTriReg0 + MapType]->
                          GetMapperOrig(MapVerts, MapEdges);
                        break;

      case Brick:
      case Hexahedron: TDatabase::MapperDB[MapQuadReg0 + MapType]->
                         GetMapperOrig(MapVerts, MapEdges);
                       break;
      default:
      break;
    }
  else
  {
    cerr << "Error in GetMapperOrig: wrong MapType " << MapType << endl;
    exit (-1);
  }
}
int TJoint::GetNeighbourEdgeIndex(TBaseCell* me, int LocEdge)
{
  int N_Edges, MaxLen;
  const int *TmpEV;
  TVertex *Vert0, *Vert1;
  TBaseCell* neigh;

  // Set neigh as the neighbour of me
  if(me == Neighb0) neigh=Neighb1;
  else neigh = Neighb0;

  if (me && neigh)
  {
    // Get Vertices of the edge
    me->GetShapeDesc()->GetEdgeVertex(TmpEV);
    Vert0 = me->GetVertex(TmpEV[LocEdge*2]);
    Vert1 = me->GetVertex(TmpEV[LocEdge*2+1]);

    neigh->GetShapeDesc()->GetEdgeVertex(TmpEV);
    N_Edges = neigh->GetShapeDesc()->GetN_Edges();
    // Iterate over all edges of neighbour
    for(int i=0; i<N_Edges; ++i)
    {
        // Check if the vertices of that edge are the ones we searched for
        if((neigh->GetVertex(TmpEV[2*i]) == Vert0 && neigh->GetVertex(TmpEV[2*i+1]) == Vert1) ||
           (neigh->GetVertex(TmpEV[2*i]) == Vert1 && neigh->GetVertex(TmpEV[2*i+1]) == Vert0))
           return i;
    }

    return -1;
  }
  else
  {
    return -1;
  }
}

#endif // __3D__

// Destructor
TJoint::~TJoint()
{
  if(Neighb0)
  { Neighb0 = NULL;}
  
  if(Neighb1)
  { Neighb1 = NULL;}
  
}
