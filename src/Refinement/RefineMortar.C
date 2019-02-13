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
// @(#)RefineMortar.C        1.6 10/18/99
//
// Purpose:     refinement of a mortar cell
//
// Author:      Volker Behns  16.01.98
//
// =======================================================================

#ifndef __MORTAR__
#define __MORTAR__
#endif

#include <GridCell.h>
#include <JointEqN.h>
#include <Database.h>
#include <Quadrangle.h>
#include <RefDesc.h>
#include <Vertex.h>

#ifdef __2D__
  #include <InterfaceJoint.h>
#endif

#include <stdlib.h>

static StoreGeomMortar TmpM[MAXN_ORIGEDGES];

#ifdef __2D__
int TGridCell::RefineMortar(int reflevel)
{
  int i, j, N_1, N_2, MaxLen1, MaxLen2, MaxLen3, index;
  const int *TmpValues1, *TmpValues2, *TmpValues3, *TmpValues4;
  double TmpX, TmpY, T_0, T_1;
  const double *TmpPos;
  
  // allocate needed fields
  TVertex **NewVertices;
  TJoint **NewJoints;

  N_1 = RefDesc->GetN_OrigEdges();
  N_2 = RefDesc->GetN_Children();

  NewVertices = new TVertex*[RefDesc->GetN_Vertices()];
  NewJoints = new TJoint*[RefDesc->GetN_Edges()];

  for (i=0;i<MAXN_ORIGEDGES;i++)
  {
    TmpM[i].Vertices = new TVertex*[N_2+1];
    TmpM[i].Joints = new TJoint*[N_2];
  }

  // check edges of matching and return all allready existing objects
  for (i=0;i<N_1;i++)
//    if (RefDesc->GetEdgeRef(i))
      if (Joints[i]->CheckMatchingRef(this, i, TmpM[i]))
      {
        cerr << "Error: non matching edges (mortar)!!!" << endl;
        exit(-1);
      }

  // create children
  N_1 = RefDesc->GetN_Children();
  Children = (TBaseCell **) new TGridCell*[N_1];

  for (i=0;i<N_1;i++)
  {
    Children[i] = new TGridCell(TDatabase::RefDescDB[RefDesc->
                        GetChildType(i)], reflevel);

    Children[i]->SetParent(this);
  }

  // copy existing vertices
  N_1 = RefDesc->GetN_NewVertEqOldVert();
  RefDesc->GetNewVertEqOldVert(TmpValues1, TmpValues2);

  for (i=0;i<N_1;i++)
    NewVertices[TmpValues1[i]] = Vertices[TmpValues2[i]];

  // create new inner vertices
  N_1 = RefDesc->GetN_InnerVertices();
  RefDesc->GetInnerVerts(TmpValues1, TmpPos, MaxLen1);

  for (i=0;i<N_1;i++)
  {
    TmpX = TmpY = 0;
    for (j=0;j<MaxLen1;j++)
    {
      TmpX += TmpPos[i * MaxLen1 + j] * Vertices[j]->GetX();
      TmpY += TmpPos[i * MaxLen1 + j] * Vertices[j]->GetY();
    }

    NewVertices[TmpValues1[i]] = new TVertex(TmpX, TmpY);
  }

  // copy existing joints
  N_1 = RefDesc->GetN_NewEdgeEqOldEdge();
  RefDesc->GetNewEdgeEqOldEdge(TmpValues1, TmpValues2);
  RefDesc->GetEdgeChild(TmpValues4, TmpValues3, MaxLen3);

  for (i=0;i<N_1;i++)
  {
    index = TmpValues1[i];

    if (TmpM[TmpValues2[i]].Filled)
      NewJoints[index] = TmpM[TmpValues2[i]].Joints[0];
    else
      NewJoints[index] = Joints[TmpValues2[i]]->NewInst();

    NewJoints[index]->SetNeighbour(Children[
                        TmpValues4[index * MaxLen3]]);
  }

  // create new inner edges
  N_1 = RefDesc->GetN_InnerEdges();
  RefDesc->GetInnerEdges(TmpValues1, TmpValues2, MaxLen1);

  for (i=0;i<N_1;i++)
     NewJoints[TmpValues1[i]] = new
          TJointEqN(Children[TmpValues2[TmpValues1[i] * MaxLen1]],
                    Children[TmpValues2[TmpValues1[i] * MaxLen1 + 1]]);

  // refine old edges if necessary
  N_1 = RefDesc->GetN_OrigEdges();

  RefDesc->GetOldEdgeNewVertex(TmpValues1, TmpValues3, MaxLen1);
  RefDesc->GetOldEdgeNewEdge(TmpValues2, TmpValues3, MaxLen2);

  for (i=0;i<N_1;i++)
  {
    if (RefDesc->GetEdgeRef(i))
     {
      if (TmpM[i].Filled)
      {
        N_2 = TmpValues3[i]+1;
        for (j=1;j<N_2;j++)
          NewVertices[TmpValues1[i * MaxLen1 + j]] = TmpM[i].Vertices[N_2-j];

        for (j=0;j<N_2;j++)
        {
          index = TmpValues2[i * MaxLen2 + j];
          NewJoints[index] = TmpM[i].Joints[N_2-1-j];
          NewJoints[index]->SetNeighbour(Children[
                              TmpValues4[index * MaxLen3]]);
        }
      }
      else
      {
        N_2 = TmpValues3[i]+1;
        for (j=1;j<N_2;j++)
          if (LineMidXY(i, j, TmpX, TmpY))
          {
            cerr << "Error: can't generate mid point!!!" << endl;
            return -1;
          }
          else
            NewVertices[TmpValues1[i * MaxLen1 + j]] = new TVertex(TmpX, TmpY);

        for (j=0;j<N_2;j++)
        {
          LineMidT(i, j, T_0, T_1);
          index = TmpValues2[i * MaxLen2 + j];
          NewJoints[index] = Joints[i]->NewInst(T_0, T_1,
                               Children[TmpValues4[index * MaxLen3]]);
        }
      }
     }
  }

  // distribute vertex and joint pointer to children
  N_1 = RefDesc->GetN_Children();
  RefDesc->GetChildVertex(TmpValues1, MaxLen1);
  RefDesc->GetChildEdge(TmpValues2, MaxLen2);

  for (i=0;i<N_1;i++)
  {
    N_2 = TDatabase::ShapeDB[RefDesc->GetChildType(i)]->GetN_Vertices();

    for (j=0;j<N_2;j++)
      Children[i]->SetVertex(j, NewVertices[TmpValues1[i * MaxLen1 + j]]);

    for (j=0;j<N_2;j++)
      Children[i]->SetJoint(j, NewJoints[TmpValues2[i * MaxLen2 + j]]);
  }

  N_1 = RefDesc->GetN_Edges();
  for (i=0;i<N_1;i++)
    if (NewJoints[i]->GetType() == InterfaceJoint ||
        NewJoints[i]->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) NewJoints[i])->CheckOrientation();

  N_1 = RefDesc->GetN_Children();
  for (i=0;i<N_1;i++)
    if (Children[i]->GetType() == Quadrangle)
      Children[i]->SetRefDesc(TDatabase::RefDescDB[((TQuadrangle *)
             TDatabase::RefDescDB[Quadrangle]->GetShapeDesc())->
             CheckQuad(((TGridCell *) Children[i])->GetVertices())]);

  delete NewVertices;
  delete NewJoints;
  
  for (i=0;i<MAXN_ORIGEDGES;i++)
  {
    delete TmpM[i].Vertices;
    delete TmpM[i].Joints;
  }

  return 0;
}
#else

int TGridCell::RefineMortar(int reflevel)
{
  return 0;
}
#endif // __2D__
