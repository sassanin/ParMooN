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
// @(#)RefDesc.C        1.3 11/15/99
//
// Class:       TRefDesc
// Purpose:     super class of all refinement descriptors
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#include <RefDesc.h>
#include <Database.h>

// Constructor
TRefDesc::TRefDesc(TShapeDesc *shape)
{
  Shape = shape;

  Type = NoRef;

  N_OrigEdges = Shape->GetN_Edges();
  N_OrigVertices = Shape->GetN_Vertices();

  N_Edges = 0;
  N_Vertices = 0;
  
  N_Children = 0;

  #ifdef __3D__
    N_OrigFaces = Shape->GetN_Faces();
    N_Faces = 0;
  #endif

  N_NewVertEqOldVert = 0;
  N_InnerVertices = 0;
  N_NewEdgeEqOldEdge = 0;
  N_InnerEdges = 0;

  MaxN_VpC = 0;
  MaxN_CpV = 0;
  MaxN_EpC = 0;
  MaxN_CpE = 0;
  MaxN_EpV = 0;
  MaxN_iVpE = 0;
  MaxN_nVpoE = 0;
  MaxN_nEpoE = 0;

  ChildType = NULL;
  EdgeType = NULL;

  ChildVertex = NULL;
  VertexChild = NULL;
  VertexChildIndex = NULL;
  VertexChildLen = NULL;

  ChildEdge = NULL;
  EdgeChild = NULL;
  EdgeChildIndex = NULL;
  EdgeChildLen = NULL;

  EdgeVertex = NULL;
  VertexEdge = NULL;
  VertexEdgeIndex = NULL;
  VertexEdgeLen = NULL;

  NewVertexEqOldVertex = NULL;
  NewVertexEqOldVertexIndex = NULL;

  NewEdgeEqOldEdge = NULL;
  NewEdgeEqOldEdgeIndex = NULL;

  InteriorEdgeOfCell = NULL;
  InteriorVertexOfEdge = NULL;
  InteriorVertexOfEdgeLen = NULL;

  OldEdgeNewVertex = NULL;
  OldEdgeNewVertexLen = NULL;
  OldEdgeNewLocEdge = NULL;
  
  OldEdgeNewEdge = NULL;
  OldEdgeNewEdgeLen = NULL;
  NewEdgeOldEdge = NULL;

  #ifdef __3D__
  
  FaceType = NULL;

  ChildFace = NULL;
  FaceChild = NULL;
  FaceChildIndex = NULL;
  FaceChildLen = NULL;

  FaceVertex = NULL;
  VertexFace = NULL;
  VertexFaceIndex = NULL;
  VertexFaceLen = NULL;

  FaceEdge = NULL;
  EdgeFace = NULL;
  EdgeFaceIndex = NULL;
  EdgeFaceLen = NULL;

  InteriorFaceOfCell = NULL;
  InteriorVertexOfFace = NULL;
  InteriorVertexOfFaceLen = NULL;
  InteriorEdgeOfFace = NULL;
  InteriorEdgeOfFaceLen = NULL;

  NewFaceOldFace = NULL;

  OldFaceNewVertex = NULL;
  OldFaceNewVertexPos = NULL;
  OldFaceNewVertexLen = NULL;
  OldFaceNewEdge = NULL;
  OldFaceNewEdgeLen = NULL;
  OldFaceNewFace = NULL;
  OldFaceNewFaceLen = NULL;

  OldFaceNewLocFace = NULL;
  ChildTwistIndex = NULL;

  #endif
}
 

// Methods
