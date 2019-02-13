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
// @(#)Refinement.C        1.14 11/26/99
//
// Purpose:     refinement of a cell
//
// Author:      Volker Behns  23.07.97
//
// =======================================================================

#include <DefineParams.h>

#include <GridCell.h>
#include <JointEqN.h>
#include <Database.h>
#include <Quadrangle.h>
#include <RefDesc.h>
#include <Vertex.h>
#include <InterfaceJoint.h>

#ifdef __2D__
  #include <BoundEdge.h>
  #include <InnerInterfaceJoint.h>
#else
  #include <IsoInterfaceJoint3D.h>
  #include <BoundFace.h>
  #include <BoundComp3D.h>
  #include <BdWall.h>
  #include <Hexahedron.h>
#endif

#include <stdlib.h>

static struct StoreGeom Tmp[MAXN_ORIGEDGES];
static TVertex *NewVertices[MAXN_NEWVERTICES];
static TJoint *NewJoints[MAXN_NEWJOINTS];

#ifdef __2D__

int TGridCell::Refine(int reflevel)
{
  int i, j, N_1, N_2, MaxLen1, MaxLen2, MaxLen3, index, auxi;
  const int *TmpValues1, *TmpValues2, *TmpValues3, *TmpValues4;
  double TmpX, TmpY, T_0, T_1, auxd;
  const double *TmpPos;
  bool Inside;
  TVertex *CurrVert;
  
  if (!IsToRefine()) return 1;

  #ifdef __MORTAR__
    if (RefDesc->GetType() >= Mortar)
    {
      RefineMortar(reflevel);
      return 0;
    }
  #endif

  N_1 = RefDesc->GetN_OrigEdges();

  // check edges of matching and return all allready existing objects
  for (i=0;i<N_1;i++)
 //   if (RefDesc->GetEdgeRef(i))
      if (Joints[i]->CheckMatchingRef(this, i, Tmp[i]))
      {
        cerr << "Error: non matching edges!!!" << endl;
//         exit(-1);
      }

  // create children
  N_1 = RefDesc->GetN_Children();
  Children = (TBaseCell **) new TGridCell*[N_1];

  for (i=0;i<N_1;i++)
  {
    Children[i] = new TGridCell(TDatabase::RefDescDB[RefDesc->
                        GetChildType(i)], reflevel);

    Children[i]->SetParent(this);
    // copy physical reference
    Children[i]->SetReference_ID(this->GetReference_ID());
  }

  // copy existing vertices
  N_1 = RefDesc->GetN_NewVertEqOldVert();
  RefDesc->GetNewVertEqOldVert(TmpValues1, TmpValues2);

  for (i=0;i<N_1;i++)
    NewVertices[TmpValues1[i]] = Vertices[TmpValues2[i]];

  // create new inner vertices
  N_1 = RefDesc->GetN_InnerVertices();
  if(N_1)
  {
    RefDesc->GetInnerVerts(TmpValues1, TmpPos, MaxLen1);

    for (i=0;i<N_1;i++)
    {
      TmpX = TmpY = 0;
      auxi = i * MaxLen1;
      for (j=0;j<MaxLen1;j++)
      {
        auxd = TmpPos[auxi++];
        CurrVert = Vertices[j];
        TmpX += auxd * CurrVert->GetX();
        TmpY += auxd * CurrVert->GetY();
      }

      NewVertices[TmpValues1[i]] = new TVertex(TmpX, TmpY);
    }
  }

  // copy existing joints
  RefDesc->GetEdgeChild(TmpValues4, TmpValues3, MaxLen3);
  N_1 = RefDesc->GetN_NewEdgeEqOldEdge();
  if(N_1)
  {
    RefDesc->GetNewEdgeEqOldEdge(TmpValues1, TmpValues2);

    for (i=0;i<N_1;i++)
    {
      index = TmpValues1[i];

      if (Tmp[TmpValues2[i]].Filled)
        NewJoints[index] = Tmp[TmpValues2[i]].Joints[0];
      else
        NewJoints[index] = Joints[TmpValues2[i]]->NewInst();

      NewJoints[index]->SetNeighbour(Children[
                          TmpValues4[index * MaxLen3]]);
    }
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
      if (Tmp[i].Filled)
      {
        N_2 = TmpValues3[i]+1;

        // change for periodic boundary conditions
        if (Joints[i]->GetType() != PeriodicJoint)
        {
          auxi = i * MaxLen1;
          for (j=1;j<N_2;j++)
            NewVertices[TmpValues1[auxi + j]] = Tmp[i].Vertices[N_2-j];
        }
        else
        {
          auxi = i * MaxLen1;
          for (j=1;j<N_2;j++)
            if (LineMidXY(i, j, TmpX, TmpY))
            {
              cerr << "Error in Refinement: can't generate mid point!!!"
                   << endl;
              return -1;
            }
            else
              NewVertices[TmpValues1[auxi + j]] = new TVertex(TmpX, TmpY);
        }

        auxi = i * MaxLen2;
        for (j=0;j<N_2;j++)
        {
          index = TmpValues2[auxi + j];
          NewJoints[index] = Tmp[i].Joints[N_2-1 - j];
          NewJoints[index]->SetNeighbour(Children[
                              TmpValues4[index * MaxLen3]]);
        }
      }
      else
      {
        if (Joints[i]->GetType() == InterfaceJoint ||
            Joints[i]->GetType() == IsoInterfaceJoint)
          Inside = ((TInterfaceJoint *) Joints[i])->CheckInside(this);
        else
          Inside = true;

        N_2 = TmpValues3[i]+1;

        auxi = i * MaxLen1;
        for (j=1;j<N_2;j++)
          if (LineMidXY(i, j, TmpX, TmpY))
          {
            cerr << "Error in Refinement: can't generate mid point!!!"
                 << endl;
            return -1;
          }
          else
            if (Inside)
              NewVertices[TmpValues1[auxi + j]] = new TVertex(TmpX, TmpY);
            else
              NewVertices[TmpValues1[auxi + N_2 - j]] = new
                          TVertex(TmpX, TmpY);

        auxi = i * MaxLen2;
        for (j=0;j<N_2;j++)
        {
          LineMidT(i, j, T_0, T_1);

          if (Inside)
            index = TmpValues2[auxi + j];
          else
            index = TmpValues2[auxi + N_2 -1 - j];

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

    auxi = i * MaxLen1;
    for (j=0;j<N_2;j++)
      Children[i]->SetVertex(j, NewVertices[TmpValues1[auxi + j]]);

    auxi = i * MaxLen2;
    for (j=0;j<N_2;j++)
      Children[i]->SetJoint(j, NewJoints[TmpValues2[auxi + j]]);
  }

  N_1 = RefDesc->GetN_Edges();
  for (i=0;i<N_1;i++)
  {
    if (NewJoints[i]->GetType() == InterfaceJoint ||
        NewJoints[i]->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) NewJoints[i])->CheckOrientation();
    else if (NewJoints[i]->GetType() == InnerInterfaceJoint)
    {
      const int *NewEdgeOnWhichOldEdge;
      RefDesc->GetNewEdgeOldEdge(NewEdgeOnWhichOldEdge);
      TJoint *joint = Joints[NewEdgeOnWhichOldEdge[i]];
      ((TInnerInterfaceJoint *)joint)->SetChild(
                                        (TInnerInterfaceJoint *)NewJoints[i]);
      RefDesc->GetEdgeVertex(TmpValues1); // which Vertices lie on an edge
      ((TInnerInterfaceJoint *)NewJoints[i])->SetParams(
        NewVertices[TmpValues1[2*i]]->GetX(),
        NewVertices[TmpValues1[2*i]]->GetY(),
        NewVertices[TmpValues1[2*i+1]]->GetX()-NewVertices[TmpValues1[2*i]]->GetX(),
        NewVertices[TmpValues1[2*i+1]]->GetY()-NewVertices[TmpValues1[2*i]]->GetY());
      // which Index has this edge in the child(ren) it belongs to
      // we only look at edges which are on the boundary of this cell
      // (belong to only one child)
      RefDesc->GetEdgeChildIndex(TmpValues1,TmpValues2,MaxLen1);
      // which Index have the children to which this edge belongs
      // we only look at edges which are on the boundary of this cell 
      // (belong to only one child)
      RefDesc->GetEdgeChild(TmpValues3,TmpValues4,MaxLen1);
      ((TInnerInterfaceJoint *)NewJoints[i])->SetIndexInNeighbor(
        Children[TmpValues3[i*MaxLen1]],TmpValues1[i*MaxLen1]);
    }
  }

  N_1 = RefDesc->GetN_Children();
  for (i=0;i<N_1;i++)
    if (Children[i]->GetType() == Quadrangle)
      Children[i]->SetRefDesc(TDatabase::RefDescDB[((TQuadrangle *)
             TDatabase::RefDescDB[Quadrangle]->GetShapeDesc())->
             CheckQuad(((TGridCell *) Children[i])->GetVertices())]);

  #ifdef __CORRECT_MP__
    if (RefDesc->GetType() == QuadReg)
    {
      RefDesc->GetInnerVerts(TmpValues1, TmpPos, MaxLen1);
      RefDesc->GetOldEdgeNewVertex(TmpValues2, TmpValues3, MaxLen2);

      for (i=0;i<4;i++)
        if (Joints[i]->GetType() == BoundaryEdge ||
            Joints[i]->GetType() == InterfaceJoint)
        {
          CurrVert = NewVertices[TmpValues2[3*i+1]];
          TmpX = CurrVert->GetX();
          TmpY = CurrVert->GetY();

          CurrVert = NewVertices[TmpValues2[3*((i+2)%4)+1]];
          TmpX = 0.5*(TmpX + CurrVert->GetX());
          TmpY = 0.5*(TmpY + CurrVert->GetY());

          NewVertices[TmpValues1[0]]->SetCoords(TmpX, TmpY);
        }
    }
  #endif // __CORRECT_MP__

  return 0;
}
#else // __3D__

int TGridCell::Refine(int reflevel)
{
  int I, i, j, k, l, m;
  int N_1, N_2, N_3;
  int MaxLen1, MaxLen2, MaxLen3, MaxLen4;
  int index, auxi, auxj;
  const int *TmpValues1, *TmpValues2, *TmpValues3, *TmpValues4;
  const int *TmpIndex, *TmpLen1, *TmpLen2, *TmpLen3, *TmpLen4;
  const int *TmpMV, *TmpME, *TmpMF;
  double TmpX, TmpY, TmpZ, T_0, T_1, auxd;
  double TmpT, TmpS;
  TBoundFace *bdface;
  TBoundComp3D *bdcomp;
  const double *TmpPos;
  bool Inside;
  Refinements FaceRefDescID;
  TRefDesc *FaceRefDesc;
  TShapeDesc *ChildDesc;
  TVertex *CurrVert;
  TBaseCell *CurrCell, *StopCell;
  TJoint *CurrJoint;

  double Param1[100], Param2[100];
  double ParentParam1[4], ParentParam2[4];
  const int *TmpFaceVal1, *TmpFaceVal2;
  int FaceMaxLen;

  const int *MapOrigVert, *MapOrigEdge;
  int LocFace1, LocFace2, LocEdge1, LocEdge2, LocFace1tmp;
  int CurrEdge, LocVertex, CurrCellChild, CurrCellChildVertex;
  int N_Faces, N_Edges, CurrCellVertex;
  Refinements result;

  double X[8], Y[8], Z[8];
  double T[4], S[4], LinComb[4];
  const int *TmpCTI;

  int RefinementOrder[MAXN_JOINTS];

  if (!IsToRefine()) return 1;

  N_1 = RefDesc->GetN_OrigFaces();
  // check faces of matching and return all already existing objects
  for (i=0;i<N_1;i++)
    if (Joints[i]->CheckMatchingRef(this, i, Tmp[i]))
    {
      cerr << "Error: non matching faces!!!" << endl;
      exit(-1);
    }

  // determine the order in which the joint should be refined
  // first: Boundary Faces
  // last: InterfaceJoints ???
  for(i=0;i<N_1;i++)
  {
    if(Joints[i]->GetType() == BoundaryFace ||
        Joints[i]->GetType() == IsoBoundFace)
    {
      bdface = (TBoundFace*)(Joints[i]);
      RefinementOrder[i] = (bdface->GetBoundComp()->GetType())*10000+i;
    }
    else
    {
      if(Joints[i]->GetType() == InterfaceJoint3D ||
          Joints[i]->GetType() == IsoInterfaceJoint3D)
      {
        RefinementOrder[i] = 
          (((TInterfaceJoint3D*)(Joints[i]))->GetBoundComp()
                ->GetType())*5000 + i;
      }
      else
        RefinementOrder[i] = 2000 + i;
    }
  } // endfor i
  // sort array RefinementOrder
  for(i=0;i<N_1-1;i++)
  {
    for(j=i+1;j<N_1;j++)
    {
      if(RefinementOrder[i]<RefinementOrder[j])
      {
        k = RefinementOrder[i];
        RefinementOrder[i] = RefinementOrder[j];
        RefinementOrder[j] = k;
      }
    }
  }
  for(i=0;i<N_1;i++)
    RefinementOrder[i] %= 1000;

  // initialize NewVertices
  for (i=0;i<MAXN_NEWVERTICES;i++)
    NewVertices[i] = NULL;

  // create children
  N_1 = RefDesc->GetN_Children();
  Children = (TBaseCell **) new TGridCell*[N_1];

  for (i=0;i<N_1;i++)
  {
    Children[i] = new TGridCell(TDatabase::RefDescDB[RefDesc->
                        GetChildType(i)], reflevel);

    Children[i]->SetParent(this);
    Children[i]->SetClipBoard(i);
    
    // copy physical reference
    Children[i]->SetReference_ID(this->GetReference_ID());
    Children[i]->SetRegionID(this->GetRegionID());
    Children[i]->SetAsLayerCell(this->IsLayerCell());  
    
  }

  // copy existing vertices
  N_1 = RefDesc->GetN_NewVertEqOldVert();
  RefDesc->GetNewVertEqOldVert(TmpValues1, TmpValues2);

  for (i=0;i<N_1;i++)
    NewVertices[TmpValues1[i]] = Vertices[TmpValues2[i]];

  // create new inner vertices
  if( (N_1 = RefDesc->GetN_InnerVertices()) )
  {
    RefDesc->GetInnerVerts(TmpValues1, TmpPos, MaxLen1);

    for (i=0;i<N_1;i++)
    {
      TmpX = TmpY = TmpZ = 0;
      auxi = i * MaxLen1;
      for (j=0;j<MaxLen1;j++)
      {
        auxd = TmpPos[auxi++];
        CurrVert = Vertices[j];
        TmpX += auxd * CurrVert->GetX();
        TmpY += auxd * CurrVert->GetY();
        TmpZ += auxd * CurrVert->GetZ();
      }

      NewVertices[TmpValues1[i]] = new TVertex(TmpX, TmpY, TmpZ);
    }
  }

  // copy existing joints
  if ((N_1 = RefDesc->GetN_NewFaceEqOldFace()))
  {
    RefDesc->GetNewFaceEqOldFace(TmpValues1, TmpValues2);
    RefDesc->GetFaceChild(TmpValues4, TmpValues3, MaxLen3);

    for (i=0;i<N_1;i++)
    {
      index = TmpValues1[i];

      if (Tmp[TmpValues2[i]].Filled)
        NewJoints[index] = Tmp[TmpValues2[i]].Joints[0];
      else
        NewJoints[index] = Joints[TmpValues2[i]]->NewInst();

      NewJoints[index]->SetNeighbour(Children[
                          TmpValues4[index * MaxLen3]]);
    }
  }

  // create new inner joints
  N_1 = RefDesc->GetN_InnerFaces();
  RefDesc->GetInnerFaces(TmpValues1, TmpValues2, MaxLen1);

  for (i=0;i<N_1;i++)
  {
    auxj = TmpValues1[i];
    auxi = auxj * MaxLen1;
    NewJoints[auxj] = new TJointEqN(Children[TmpValues2[auxi]],
                                    Children[TmpValues2[auxi + 1]]);
  }

  // copy existing vertices on old faces
  N_1 = RefDesc->GetN_OrigFaces();
  RefDesc->GetOldFaceNewVertex(TmpValues1, TmpPos, TmpLen1,
                               MaxLen1, MaxLen2);
  for (i=0;i<N_1;i++)
  {
    if ((FaceRefDescID = RefDesc->GetFaceRef(i)))
      if (Tmp[i].Filled)
      {
        FaceRefDesc = TDatabase::RefDescDB[N_SHAPES + FaceRefDescID];

        N_2 = TmpLen1[i];
        
        auxi = i * MaxLen1;
        for (j=0;j<N_2;j++)
          NewVertices[TmpValues1[auxi++]] = Tmp[i].Vertices[j];
      }
  }

  // check vertices on old edge (diagonally neighboured elements)
  N_1 = RefDesc->GetN_OrigEdges();
  RefDesc->GetOldEdgeNewVertex(TmpValues1, TmpLen1, MaxLen1);
  RefDesc->GetShapeDesc()->GetEdgeFace(TmpValues2, MaxLen2);
  RefDesc->GetShapeDesc()->GetFaceEdge(TmpValues3, TmpLen2, MaxLen3);

  for (i=0;i<N_1;i++)
  {
    N_2 = TmpLen1[i] + 1;
    LocVertex = TmpValues1[i * MaxLen1 + 1];
    auxi = i * MaxLen2;
    // if oldedge i is not refined, then LocVertex == Vert2
    if (!NewVertices[LocVertex])
    {
      LocFace1 = TmpValues2[auxi++];
      
      CurrJoint = Joints[LocFace1];
      if ((CurrCell = CurrJoint->GetNeighbour(this)))
      {
        StopCell = Joints[TmpValues2[auxi]]->GetNeighbour(this);
        LocFace1tmp = TmpValues2[auxi];
      }
      else
      {
        LocFace1 = TmpValues2[auxi];
        CurrJoint = Joints[LocFace1];
        CurrCell = CurrJoint->GetNeighbour(this);
        StopCell = NULL;
      }
      
      // Get Local Index of Edge at Face
      RefDesc->GetShapeDesc()->GetFaceEdge(TmpValues3, TmpLen2, MaxLen3);
      auxj = LocFace1 * MaxLen3;
      N_Edges = TmpLen2[LocFace1];
      for (k=0;k<N_Edges;k++)
        if (TmpValues3[auxj++] == i)
          break;
      LocEdge1 = k;
        
      while ((CurrCell != StopCell) && (CurrCell != NULL))
      {
        CurrJoint->GetMapperOrig(MapOrigVert,MapOrigEdge);
        // Get Local Index of Joint at Cell
        N_Faces = CurrCell->GetN_Faces();
        for(k = 0;k < N_Faces;k++)
          if (CurrCell->GetJoint(k) == CurrJoint)
            break;
        LocFace2 = k;
        CurrCell->GetRefDesc()->GetShapeDesc()->
                    GetFaceEdge(TmpValues4, TmpLen4, MaxLen4);
        LocEdge2 = MapOrigEdge[LocEdge1];
        auxj = LocFace2 * MaxLen4 + LocEdge2;
        CurrEdge = TmpValues4[auxj];

        // At this point, I won't know, how is CurrEdge refined
        // maybe I guess the right vertex
        result = CurrCell->GetRefDesc()->GetEdgeRef(CurrEdge);
        if (result == LineReg)
        {
          CurrCell->GetRefDesc()->
                    GetOldEdgeNewVertex(TmpValues4, TmpLen4, MaxLen4);
          CurrCellVertex = TmpValues4[0 + CurrEdge * MaxLen4 + 1];
          CurrCell->GetRefDesc()->GetVertexChild(TmpValues4, TmpLen4, MaxLen4);
          CurrCellChild = TmpValues4[0 + CurrCellVertex * MaxLen4];
          CurrCell->GetRefDesc()->
                    GetVertexChildIndex(TmpValues4, TmpLen4, MaxLen4);
          CurrCellChildVertex = TmpValues4[0 + CurrCellVertex * MaxLen4];
          NewVertices[LocVertex] = CurrCell->GetChild(CurrCellChild)->
                                     GetVertex(CurrCellChildVertex);
          break;
        }

        CurrCell->GetRefDesc()->GetShapeDesc()->
                  GetEdgeFace(TmpValues3, MaxLen3);
        auxj = CurrEdge * MaxLen3;
        if (TmpValues3[auxj] != LocFace2)
          LocFace1 = TmpValues3[auxj];
        else
          LocFace1 = TmpValues3[auxj + 1];

        CurrCell->GetRefDesc()->GetShapeDesc()->
                  GetFaceEdge(TmpValues3, TmpLen2, MaxLen3);
        auxj = LocFace1 * MaxLen3;
        N_Edges = TmpLen2[LocFace1];
        for (k=0;k<N_Edges;k++)
          if (TmpValues3[auxj++] == CurrEdge)
            break;
        LocEdge1 = k;
          
        CurrJoint = CurrCell->GetJoint(LocFace1);
        CurrCell = CurrJoint->GetNeighbour(CurrCell);
      }
      if ((StopCell != CurrCell) && (!NewVertices[LocVertex])
                                 && (StopCell != NULL) )
      {
          RefDesc->GetShapeDesc()->GetFaceEdge(TmpValues3, TmpLen2, MaxLen3);
          auxj = LocFace1tmp * MaxLen3;
          N_Edges = TmpLen2[LocFace1tmp];
          for (k=0;k<N_Edges;k++)
              if (TmpValues3[auxj++] == i)
                  break;
          LocEdge1 = k;
          CurrCell = StopCell;
          CurrJoint = Joints[LocFace1tmp];
          while ((CurrCell != NULL))
          {
              CurrJoint->GetMapperOrig(MapOrigVert,MapOrigEdge);
              N_Faces = CurrCell->GetN_Faces();
              for(k = 0;k < N_Faces;k++)
                  if (CurrCell->GetJoint(k) == CurrJoint)
                      break;
              LocFace2 = k;
              CurrCell->GetRefDesc()->GetShapeDesc()->
                  GetFaceEdge(TmpValues4, TmpLen4, MaxLen4);
              LocEdge2 = MapOrigEdge[LocEdge1];
              auxj = LocFace2 * MaxLen4 + LocEdge2;
              CurrEdge = TmpValues4[auxj];
              
              // At this point, I won't know, how is CurrEdge refined
              // maybe I guess the right vertex
              result = CurrCell->GetRefDesc()->GetEdgeRef(CurrEdge);
              if (result == LineReg)
              {
                  CurrCell->GetRefDesc()->
                        GetOldEdgeNewVertex(TmpValues4, TmpLen4, MaxLen4);
                  CurrCellVertex = TmpValues4[0 + CurrEdge * MaxLen4 + 1];
                  CurrCell->GetRefDesc()->
                        GetVertexChild(TmpValues4, TmpLen4, MaxLen4);
                  CurrCellChild = TmpValues4[0 + CurrCellVertex * MaxLen4];
                  CurrCell->GetRefDesc()->
                        GetVertexChildIndex(TmpValues4, TmpLen4, MaxLen4);
                  CurrCellChildVertex =
                      TmpValues4[0 + CurrCellVertex * MaxLen4];
                  NewVertices[LocVertex] = CurrCell->GetChild(CurrCellChild)->
                        GetVertex(CurrCellChildVertex);
                  break;
              }
                
              CurrCell->GetRefDesc()->GetShapeDesc()->
                  GetEdgeFace(TmpValues3, MaxLen3);
              auxj = CurrEdge * MaxLen3;
              if (TmpValues3[auxj] != LocFace2)
                  LocFace1 = TmpValues3[auxj];
              else
                  LocFace1 = TmpValues3[auxj + 1];
              
              CurrCell->GetRefDesc()->GetShapeDesc()->
                  GetFaceEdge(TmpValues3, TmpLen2, MaxLen3);
              auxj = LocFace1 * MaxLen3;
              N_Edges = TmpLen2[LocFace1];
              for (k=0;k<N_Edges;k++)
                  if (TmpValues3[auxj++] == CurrEdge)
                      break;
              LocEdge1 = k;
              
              CurrJoint = CurrCell->GetJoint(LocFace1);
              CurrCell = CurrJoint->GetNeighbour(CurrCell);
            }  
        }
    }
  }

  // refine old faces if necessary
  RefDesc->GetFaceChild(TmpValues1, TmpLen1, MaxLen1);
  RefDesc->GetOldFaceNewVertex(TmpValues2, TmpPos, TmpLen2,
                               MaxLen2, MaxLen4);
  RefDesc->GetOldFaceNewFace(TmpValues3, TmpLen3, MaxLen3);
  RefDesc->GetChildTwistIndex(TmpCTI);

  N_1 = RefDesc->GetN_OrigFaces();
  for (I=0;I<N_1;I++)
  {
    i = RefinementOrder[I];
    // cout << "order: " << i << endl;
    if ((FaceRefDescID = RefDesc->GetFaceRef(i)))
    {
      FaceRefDesc = TDatabase::RefDescDB[N_SHAPES + FaceRefDescID];
      if (Tmp[i].Filled)
      {
        // change for periodic boundary conditions
        if (Joints[i]->GetType() == PeriodicJoint)
        {
          // create new vertices on face i (if they do not exist)
          FaceRefDesc->GetNewVertEqOldVert(TmpValues4, TmpIndex);

          N_2 = TmpLen2[i];
          N_3 = FaceRefDesc->GetN_OrigVertices();
          auxi = i * MaxLen2;
          for (j=0;j<N_2;j++)
            if (!NewVertices[TmpValues2[auxi + j]])
            {
              TmpX = TmpY = TmpZ = 0;
              auxj = (auxi + j) * MaxLen4;
              for (k=0;k<N_3;k++)
              {
                auxd = TmpPos[auxj + TmpIndex[k]];
                CurrVert = Vertices[TmpValues2[auxi + TmpValues4[k]]];
                TmpX += auxd * CurrVert->GetX();
                TmpY += auxd * CurrVert->GetY();
                TmpZ += auxd * CurrVert->GetZ();
              }

              NewVertices[TmpValues2[auxi + j]] = new TVertex(TmpX, TmpY, TmpZ);
            }
        }

        N_2 = FaceRefDesc->GetN_Children();
        auxi = i * MaxLen3;
        for (j=0;j<N_2;j++)
        {
          index = TmpValues3[auxi++];
          CurrJoint = NewJoints[index] = Tmp[i].Joints[j];
          CurrJoint->SetNeighbour(Children[TmpValues1[index * MaxLen1]]);
        }
      }
      else
      {
        if (Joints[i]->GetType() == InterfaceJoint3D)
        {
          Inside = ((TInterfaceJoint3D *) Joints[i])->CheckInside(this);
          // cout << "II ";
        }
        else
          Inside = true;

        // cout << i << " Inside: " << Inside << endl;

        // create new vertices on face i (if they do not exist)
        FaceRefDesc->GetNewVertEqOldVert(TmpValues4, TmpIndex);

        N_2 = TmpLen2[i];
        N_3 = FaceRefDesc->GetN_OrigVertices();
        auxi = i * MaxLen2;

        if(Joints[i]->GetType() == BoundaryFace ||
                Joints[i]->GetType() == IsoBoundFace ||
                Joints[i]->GetType() == InterfaceJoint3D ||
                Joints[i]->GetType() == IsoInterfaceJoint3D)
        {
          if(Joints[i]->GetType() == BoundaryFace ||
                Joints[i]->GetType() == IsoBoundFace) 
          {
            bdface = (TBoundFace*)(Joints[i]);
            bdcomp = bdface->GetBoundComp();
            bdface->GetParameters(ParentParam1, ParentParam2);
          }

          if(Joints[i]->GetType() == InterfaceJoint3D ||
                Joints[i]->GetType() == IsoInterfaceJoint3D)
          {
            bdcomp = ((TInterfaceJoint3D*)Joints[i])->GetBoundComp();
            ((TInterfaceJoint3D*)Joints[i])->GetParameters
                        (ParentParam1, ParentParam2);
	    //cout << "ParentParam1: " << (long int)(Joints[i]) << " ";
	    //for(k=0;k<4;k++)
            //   cout << ParentParam1[k] << "  ";
            // cout << endl;
          }

          if(!Inside)
          {
          }

          k = FaceRefDesc->GetN_NewVertEqOldVert();
          for(j=0;j<k;j++)
          {
            Param1[MAXN_nVpoJ*i+TmpValues4[j]] = ParentParam1[TmpIndex[j]];
            Param2[MAXN_nVpoJ*i+TmpValues4[j]] = ParentParam2[TmpIndex[j]];
            //cout << Param1[MAXN_nVpoJ*i+TmpValues4[j]] << " -- ";
            //cout << Param2[MAXN_nVpoJ*i+TmpValues4[j]] << endl;
          }

          for(j=0;j<N_2;j++)
          {
            if (!NewVertices[TmpValues2[auxi + j]] ||
                bdcomp->GetType() == Wall )
            {
              TmpX = TmpY = TmpZ = 0;
              TmpT = TmpS = 0;
              auxj = (auxi + j) * MaxLen4;
              for (k=0;k<N_3;k++)
              {
                LinComb[k] = TmpPos[auxj + TmpIndex[k]];
                T[k] = ParentParam1[k];
                S[k] = ParentParam2[k];
                CurrVert = Vertices[TmpValues2[auxi + TmpValues4[k]]];
                CurrVert->GetCoords(X[k], Y[k], Z[k]);
		//cout <<k << " k " << X[k] <<  " " << Y[k] << " " << Z[k] << endl;
              }
              bdcomp->GetXYZandTS(N_3, LinComb, X, Y, Z, T, S,
                                  TmpX, TmpY, TmpZ, TmpT, TmpS);
              //cout << j << " ";
              //cout << "T: " << TmpT << " " << TmpS << " ";
              //cout << "X: " << TmpX << " " << TmpY << " " << TmpZ << endl;
              Param1[MAXN_nVpoJ*i+j] = TmpT;
              Param2[MAXN_nVpoJ*i+j] = TmpS;

              // bdcomp->GetXYZofTS(TmpT, TmpS, TmpX, TmpY, TmpZ);
              // cout << "T: " << TmpT << " " << TmpS << " ";
              // cout << "X: " << TmpX << " " << TmpY << " " << TmpZ << endl;
  
              if(!NewVertices[TmpValues2[auxi + j]])
                NewVertices[TmpValues2[auxi + j]] =
                                new TVertex(TmpX, TmpY, TmpZ);
              else
                NewVertices[TmpValues2[auxi + j]]->SetCoords(TmpX, TmpY, TmpZ);
            }
            else
            {
              NewVertices[TmpValues2[auxi + j]]->GetCoords(TmpX, TmpY, TmpZ);
              bdcomp->GetTSofXYZ(TmpX, TmpY, TmpZ, TmpT, TmpS);
              Param1[MAXN_nVpoJ*i+j] = TmpT;
              Param2[MAXN_nVpoJ*i+j] = TmpS;
              // cout << "T: " << TmpT << " " << TmpS << " ";
              // cout << "X: " << TmpX << " " << TmpY << " " << TmpZ << endl;
            }
          } // endfor j

          // correct the possibly wrong values in vertices
          for(j=0;j<k;j++)
          {
            Param1[MAXN_nVpoJ*i+TmpValues4[j]] = ParentParam1[TmpIndex[j]];
            Param2[MAXN_nVpoJ*i+TmpValues4[j]] = ParentParam2[TmpIndex[j]];
            // cout << Param1[MAXN_nVpoJ*i+TmpValues4[j]] << " -- ";
            // cout << Param2[MAXN_nVpoJ*i+TmpValues4[j]] << endl;
          }
          // for(j=0;j<N_2;j++)
          // {
          //   TmpT = Param1[MAXN_nVpoJ*i+j];
          //   TmpS = Param2[MAXN_nVpoJ*i+j];
          //   cout << "T: " << TmpT << " " << TmpS << endl;
          // }
        }
        else
        {
          for (j=0;j<N_2;j++)
            if (!NewVertices[TmpValues2[auxi + j]])
            {
              TmpX = TmpY = TmpZ = 0;
              auxj = (auxi + j) * MaxLen4;
              for (k=0;k<N_3;k++)
              {
                auxd = TmpPos[auxj + TmpIndex[k]];
                CurrVert = Vertices[TmpValues2[auxi + TmpValues4[k]]];
                TmpX += auxd * CurrVert->GetX();
                TmpY += auxd * CurrVert->GetY();
                TmpZ += auxd * CurrVert->GetZ();
              }

              NewVertices[TmpValues2[auxi + j]] = new TVertex(TmpX, TmpY, TmpZ);
            }
        }

        N_2 = FaceRefDesc->GetN_Children();
        FaceRefDesc->GetChildVertex(TmpFaceVal1, FaceMaxLen);
        auxi = i * MaxLen3;
        for(j=0;j<N_2;j++)
        {
          index = TmpValues3[auxi++];
          T_0 = 0.0; T_1 = 1.0;
          NewJoints[index] = Joints[i]->NewInst(T_0, T_1,
                      Children[TmpValues1[index * MaxLen1]]);
          if(NewJoints[index]->GetType() == BoundaryFace ||
                NewJoints[index]->GetType() == IsoBoundFace)
          {
            N_3 = FaceRefDesc->GetN_OrigVertices();
            m = TmpCTI[index];
            for(k=0;k<N_3;k++)
            {
              l = TmpFaceVal1[j*FaceMaxLen+k];
              ParentParam1[(m+k)%N_3] = Param1[MAXN_nVpoJ*i+l];
              ParentParam2[(m+k)%N_3] = Param2[MAXN_nVpoJ*i+l];
            }
            ((TBoundFace*)NewJoints[index])->
                SetParameters(ParentParam1, ParentParam2);
          } // endif BoundaryFace

          if(NewJoints[index]->GetType() == InterfaceJoint3D ||
                NewJoints[index]->GetType() == IsoInterfaceJoint3D)
          {
            N_3 = FaceRefDesc->GetN_OrigVertices();
            m = TmpCTI[index];
            for(k=0;k<N_3;k++)
            {
              l = TmpFaceVal1[j*FaceMaxLen+k];
              ParentParam1[(m+k)%N_3] = Param1[MAXN_nVpoJ*i+l];
              ParentParam2[(m+k)%N_3] = Param2[MAXN_nVpoJ*i+l];
            }
            ((TInterfaceJoint3D*)NewJoints[index])->
                SetParameters(ParentParam1, ParentParam2);
            // cout << endl;
            // for(k=0;k<N_3;k++)
            //   cout << ParentParam1[k] << " ";
            // cout << endl;
          } // endif InterfaceJoint3D
        } // endfor j
      } // endif Tmp[i].Filled
    } // endif FaceRefDescID
  } // endfor I

  // distribute vertex and joint pointer to children
  N_1 = RefDesc->GetN_Children();
  RefDesc->GetChildVertex(TmpValues1, MaxLen1);
  RefDesc->GetChildFace(TmpValues2, MaxLen2);

  for (i=0;i<N_1;i++)
  {
    ChildDesc = TDatabase::ShapeDB[RefDesc->GetChildType(i)];
    N_2 = ChildDesc->GetN_Vertices();
    auxi = i * MaxLen1;
    for (j=0;j<N_2;j++)
      Children[i]->SetVertex(j, NewVertices[TmpValues1[auxi + j]]);


    N_2 = ChildDesc->GetN_Faces();
    auxi = i * MaxLen2;
    for (j=0;j<N_2;j++)
      Children[i]->SetJoint(j, NewJoints[TmpValues2[auxi + j]]);
  }

  // set all MapTypes
  N_1 = RefDesc->GetN_Faces();
  for (i=0;i<N_1;i++)
    NewJoints[i]->SetMapType();
//   cout<<"````````````````````````````````````````jfcdjcgsdjcgsdjcgsdkcgsd``````````````````````````"<<endl;
  N_1 = RefDesc->GetN_Faces();
  for (i=0;i<N_1;i++)
    if (NewJoints[i]->GetType() == InterfaceJoint3D)
      ((TInterfaceJoint3D *) NewJoints[i])->CheckOrientation();

  N_1 = RefDesc->GetN_Children();
  for (i=0;i<N_1;i++)
    if (Children[i]->GetType() == Hexahedron)
      Children[i]->SetRefDesc(TDatabase::RefDescDB[((THexahedron *)
             TDatabase::RefDescDB[Hexahedron]->GetShapeDesc())->
             CheckHexa(((TGridCell *) Children[i])->GetVertices())]);

/*
  N_1 = RefDesc->GetN_Children();
  for(i=0;i<N_1;i++)
  {
    ChildDesc = TDatabase::ShapeDB[RefDesc->GetChildType(i)];
    ChildDesc->GetFaceVertex(TmpValues1, TmpValues2, l);
    N_2 = ChildDesc->GetN_Faces();
    for(j=0;j<N_2;j++)
    {
      CurrJoint = Children[i]->GetJoint(j);
      if(CurrJoint->GetType() == BoundaryFace ||
                CurrJoint->GetType() == IsoBoundFace)
      {
        cout << i << " " << j << " child with boundary face" << endl;
        bdface = (TBoundFace*)CurrJoint;
        bdface->GetParameters(Param1, Param2);
        bdcomp = bdface->GetBoundComp();
        for(k=0;k<4;k++)
        {
          // cout << k << " " << Param1[k] << " " << Param2[k] << endl;
          Children[i]->GetVertex(TmpValues1[l*j+k])->
                        GetCoords(TmpX, TmpY, TmpZ);
          // cout << TmpX << " " << TmpY << " " << TmpZ << endl;
          bdcomp->GetTSofXYZ(TmpX, TmpY, TmpZ, TmpT, TmpS);
          // cout << "T, S: " << TmpT << " " << TmpS << endl;
          if( (fabs(Param1[k]-TmpT)>1e-4) || (fabs(Param2[k]-TmpS)>1e-4))
          {
            cout << k << " " << Param1[k] << " " << Param2[k] << endl;
            cout << "T, S: " << TmpT << " " << TmpS << endl;
            cout << "X, Y, Z: " << TmpX << " " << TmpY << " " << TmpZ << endl;
          }
        }
      }
    }
  }
*/

  // cout << "ClipBoard: " << ClipBoard << endl;
  
  // added 25.04.2010 for fixing refinement problem
  CorrectBoundaryVertices(NewVertices, NewJoints);

  return 0;
}
#endif // __2D__

#ifdef __2D__
int TGridCell::Derefine()
{
  int i, j, k, N_1, N_2, MaxLen1, MaxLen2, MaxLen3, MaxLen4;
  const int *TmpIV, *TmpVC, *TmpVCI, *TmpEC, *TmpECI;
  const int *TmpoEnE, *TmpoEnV, *TmpLen;
  const double *TmpPos;
  TBaseCell *Child, *Neighb1, *Neighb2;

  RefDesc->GetInnerVerts(TmpIV, TmpPos, MaxLen1);
  RefDesc->GetVertexChild(TmpVC, TmpLen, MaxLen1);
  RefDesc->GetVertexChildIndex(TmpVCI, TmpLen, MaxLen1);

  if (TmpIV)
  {
    N_1 = RefDesc->GetN_InnerVertices();
    for(i=0;i<N_1;i++)
      delete Children[TmpVC[TmpIV[i]*MaxLen1]]->
               GetVertex(TmpVCI[TmpIV[i]*MaxLen1]);
  }

  RefDesc->GetEdgeChild(TmpEC, TmpLen, MaxLen2);
  RefDesc->GetEdgeChildIndex(TmpECI, TmpLen, MaxLen2);
  RefDesc->GetOldEdgeNewEdge(TmpoEnE, TmpLen, MaxLen3);
  RefDesc->GetOldEdgeNewVertex(TmpoEnV, TmpLen, MaxLen4);

  if (Children)
  {
    N_1 = RefDesc->GetN_OrigEdges();
    
    for (i=0;i<N_1;i++)
    {
      N_2 = TmpLen[i];
      if (N_2)
      {
        k = TmpoEnE[i*MaxLen3]*MaxLen2;
        Child = Children[TmpEC[k]];
        Neighb1 = Child->GetJoint(TmpECI[k])->
                  GetNeighbour(Child);
        for (j=1;j<=N_2;j++)
        {
          k = TmpoEnE[i*MaxLen3 + j]*MaxLen2;
          Child = Children[TmpEC[k]];
          Neighb2 = Child->GetJoint(TmpECI[k])->GetNeighbour(Child);

          if (!(Neighb1 || Neighb2))
	   {
            if (TmpVC[TmpoEnV[i*MaxLen4+j]*MaxLen1] == TmpEC[k])
	    { delete Child->GetVertex(TmpVCI[TmpoEnV[i*MaxLen4+j]*MaxLen1]);}
            else
	    { delete Child->GetVertex(TmpVCI[TmpoEnV[i*MaxLen4+j]*MaxLen1+1]);}
	   }
          Neighb1 = Neighb2;
        }
      }
     }
    N_1 = GetN_Children();
    for (i=0;i<N_1;i++)
    {
      Children[i]->Derefine();
      delete (TGridCell *) Children[i];
    }

    delete [] Children;
    Children = NULL;

    RefDesc = TDatabase::RefDescDB[RefDesc->GetShapeDesc()->GetType()];
  }

  return 0;
}
#else
int TGridCell::Derefine()
{
/*
  int i, j, k, N_1, N_2, MaxLen1, MaxLen2, MaxLen3, MaxLen4;
  const int *TmpIV, *TmpVC, *TmpVCI, *TmpEC, *TmpECI;
  const int *TmpoEnE, *TmpoEnV, *TmpLen;
  const double *TmpPos;
  TBaseCell *Child, *Neighb1, *Neighb2;

  RefDesc->GetInnerVerts(TmpIV, TmpPos, MaxLen1);
  RefDesc->GetVertexChild(TmpVC, TmpLen, MaxLen1);
  RefDesc->GetVertexChildIndex(TmpVCI, TmpLen, MaxLen1);

  if (TmpIV)
  {
    N_1 = RefDesc->GetN_InnerVertices();
    for(i=0;i<N_1;i++)
      delete Children[TmpVC[TmpIV[i]*MaxLen1]]->
               GetVertex(TmpVCI[TmpIV[i]*MaxLen1]);
  }

  RefDesc->GetFaceChild(TmpEC, TmpLen, MaxLen2);
  RefDesc->GetFaceChildIndex(TmpECI, TmpLen, MaxLen2);
  RefDesc->GetOldFaceNewFace(TmpoEnE, TmpLen, MaxLen3);
  RefDesc->GetOldFaceNewVertex(TmpoEnV, TmpLen, MaxLen4);

  if (Children)
  {
    N_1 = RefDesc->GetN_OrigFaces();
    for (i=0;i<N_1;i++)
      if (N_2 = TmpLen[i])
      {
        k = TmpoEnE[i*MaxLen3]*MaxLen2;
        Child = Children[TmpEC[k]];
        Neighb1 = Child->GetJoint(TmpECI[k])->
                  GetNeighbour(Child);
        for (j=1;j<=N_2;j++)
        {
          k = TmpoEnE[i*MaxLen3 + j]*MaxLen2;
          Child = Children[TmpEC[k]];
          Neighb2 = Child->GetJoint(TmpECI[k])->GetNeighbour(Child);

          if (!(Neighb1 || Neighb2))
            if (TmpVC[TmpoEnV[i*MaxLen4+j]*MaxLen1] == TmpEC[k])
              delete Child->GetVertex(TmpVCI[TmpoEnV[i*MaxLen4+j]*MaxLen1]);
            else
              delete Child->GetVertex(TmpVCI[TmpoEnV[i*MaxLen4+j]*MaxLen1+1]);

          Neighb1 = Neighb2;
        }
      }

    N_1 = GetN_Children();
    for (i=0;i<N_1;i++)
    {
      Children[i]->Derefine();
      delete (TGridCell *) Children[i];
    }

    delete Children;
    Children = NULL;

    RefDesc = TDatabase::RefDescDB[RefDesc->GetShapeDesc()->GetType()];
  }
*/
  if(Children)
  {
    int N_1 = GetN_Children();
    for (int i=0;i<N_1;i++)
    {
      Children[i]->Derefine();
      delete (TGridCell *) Children[i];
    }
    delete [] Children;
    Children = NULL;
  }
  return 0;
}
#endif

int TGridCell::RefDeref()
{
  int i, N_;

  if (ExistChildren() && GetClipBoard() == NoRefinement)
  {
    Derefine();
    SetClipBoard(DeRefinement);
  }
  else
    if (!ExistChildren() && GetClipBoard() == Refinement)
    {
      SetRegRefine();
      Refine(0);
      N_ = GetN_Children();
      for(i=0;i<N_;i++) GetChild(i)->SetClipBoard(Refinement);
    }

  return 0;
}

int TGridCell::Gen1RegMarks()
{
  TBaseCell *child, *cell;
  int i, j, k, l, ChildNumber;
//   int ToDerefine=0;
  const int *TmpCE, *TmpnEoE;
  const int *TmpoEnE, *TmpLen1;
  const int *TmpEC, *TmpLen2;
  int MaxLen, LocJointNum;
  int MaxLen1, MaxLen2;
  TRefDesc *RefDesc_tmp;
  int N_ = GetN_Edges();
  TJoint *CurrJoint;

  switch (GetClipBoard())
  {
    case Refinement:
        if (Parent)
        {
          RefDesc_tmp = Parent->GetRefDesc();
          RefDesc_tmp->GetChildEdge(TmpCE, MaxLen);
          RefDesc_tmp->GetNewEdgeOldEdge(TmpnEoE);

          for (i=0;i<N_;i++)
	  {
	    cell = Joints[i]->GetNeighbour(this);
            if(cell)
            {
              if (cell->GetClipBoard() == DeRefinement)
                cell->SetClipBoard(NoRefinement);
            }
            else
            {
              // there is no neighbour on the other side of the joint
              // either joint is boundary joint or 1-regular situation
              ChildNumber = Parent->GetChildNumber(this);

              LocJointNum = TmpnEoE[TmpCE[ChildNumber * MaxLen + i]];
              cell = Parent->GetJoint(LocJointNum)->GetNeighbour(Parent);
              if (cell)
	      { cell->SetClipBoard(Refinement);}
            }
	   }
        }
        break;

    case NoRefinement:
    case DeRefinement:
        for (i=0;i<N_;i++)
	 {
	   cell = Joints[i]->GetNeighbour(this);
           if (cell)
	   {
            if (cell->ExistChildren())
            {
              if (cell->GetClipBoard() == Refinement)
                SetClipBoard(NoRefinement);

              RefDesc_tmp = cell->GetRefDesc();
              j = 0;
              CurrJoint = Joints[i];
              while (cell->GetJoint(j) != CurrJoint) j++;
              // this cell is on local edge j of the neighbour cell

              RefDesc_tmp->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);
              RefDesc_tmp->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);
        
              k = TmpLen1[j]+1;
              for (l=0;l<k;l++) // loop over all new edge on old edge j in cell
              {
                ChildNumber = TmpEC[TmpoEnE[j*MaxLen1+l]*MaxLen2];
                child=cell->GetChild(ChildNumber);
                if (child->GetClipBoard() == Refinement ||
                    child->ExistChildren())
                  SetClipBoard(Refinement);
              }
            }
	   }
	 }
        break;
  }

  return 0;
}
