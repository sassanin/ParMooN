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
   
#ifndef __3D__
#define __3D__
#endif

#include <RefTetraBis15Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, LineReg, NoRef, NoRef, NoRef, LineReg};

static const Refinements DatFaceType [] = {TriBis1, NoRef, TriBis02, TriBis1};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,4,3},{0,4,2,5},{0,4,5,3}};
static const int DatVertexChild[][3] = {{0,1,2},{0},{1},{0,2},{0,1,2},{1,2}};
static const int DatVertexChildIndex[][3] = {{0,0,0},{1},{2},{3,3},{2,1,1},{3,2}};
static const int DatVertexChildLen[] = {3,1,1,2,3,2};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,2,3,8},{1,9,5,6},{9,8,4,7}};
static const int DatFaceChild[][2] = {{0},{1},{0},{0},{2},{1},{1},{2},{0,2},{1,2}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{2},{2},{2},{3},{3},{3,1},{1,0}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,1,8,4,5,9},{8,2,3,11,10,6},{8,10,11,4,9,7}};
static const int DatEdgeChild[][3] = {{0},{0},{1},{1},{0,2},{0},{1},{2},{0,1,2},{0,2},{1,2},{1,2}};
static const int DatEdgeChildIndex[][3] = {{0},{1},{1},{2},{3,3},{4},{5},{5},{2,0,0},{5,4},{4,1},{3,2}};
static const int DatEdgeChildLen[] = {1,1,1,1,2,1,1,1,3,2,2,2};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,4},{4,2},{2,0},{0,3},{1,3},{2,5},{5,3},{4,0},{4,3},{4,5},{0,5}};
static const int DatVertexEdge[][5] = {{0,3,4,8,11},{0,1,5},{2,3,6},{4,5,7,9},{1,2,8,9,10},{6,7,10,11}};
static const int DatVertexEdgeIndex[][5] = {{0,1,0,1,0},{1,0,0},{1,0,0},{1,1,1,1},{1,0,0,0,0},{1,0,1,1}};
static const int DatVertexEdgeLen[] = {5,3,3,4,5,4};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,4},{0,4,2},{0,3,1},{4,1,3},{5,4,3},{2,4,5},{0,2,5},{0,5,3},{0,3,4},{0,5,4}};
static const int DatVertexFace[][7] = {{0,1,2,6,7,8,9},{0,2,3},{1,5,6},{2,3,4,7,8},{0,1,3,4,5,8,9},{4,5,6,7,9}};
static const int DatVertexFaceIndex[][7] = {{0,0,0,0,0,0,0},{1,2,1},{2,0,1},{1,2,2,2,1},{2,1,0,1,1,2,2},{0,2,2,1,1}};
static const int DatVertexFaceLen[] = {7,3,3,5,7,5};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,1,8},{8,2,3},{4,5,0},{1,5,9},{10,9,7},{2,10,6},{3,6,11},{11,7,4},{4,9,8},{11,10,8}};
static const int DatEdgeFace[][4] = {{0,2},{0,3},{1,5},{1,6},{2,7,8},{2,3},{5,6},{4,7},{0,1,8,9},{3,4,8},{4,5,9},{6,7,9}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{1,0},{1,0},{2,0},{0,2,0},{1,1},{2,1},{2,1},{2,0,2,2},{2,1,1},{0,1,1},{2,0,0}};
static const int DatEdgeFaceLen[] = {2,2,2,2,3,2,2,2,4,3,3,3};

static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {2};
static const int DatNewFaceEqOldFaceIndex[] = {1};

static const int DatOldFaceNewVertex[][REFTETRABIS15MAXN_nVpoF] =
  { {0, 1, 2, 4}, {0, 3, 1}, {2, 1, 3, 4, 5}, {0, 2, 3, 5} };
static const int DatOldFaceNewVertexLen[] =
  { 4, 3, 5, 4};
static const double DatOldFaceNewVertexPos[][REFTETRABIS15MAXN_nVpoF][REFTETRABIS15MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} } };

static const int DatInteriorFaceOfCell[] = {8, 9};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRABIS15MAXN_iVpE] =
  { {}, {4}, {}, {}, {}, {5} };
static const int DatInteriorVertexOfEdgeLen[] = {0, 1, 0, 0, 0, 1};

static const int DatInteriorEdgeOfFace[][REFTETRABIS15MAXN_iEpF] =
  { {8}, {}, {9, 10}, {11} };
static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 0, 2, 1};

static const int DatOldEdgeNewVertex[][REFTETRABIS15MAXN_nVpoE] =
  { {0, 1}, {1, 4, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 5, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 3, 2, 2, 2, 3};

static const int DatOldEdgeNewEdge[][REFTETRABIS15MAXN_nEpoE] =
  { {0}, {1, 2}, {3}, {4}, {5}, {6, 7} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 2, 1, 1, 1, 2};

static const int DatOldFaceNewEdge[][REFTETRABIS15MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 5, 0}, {2, 1, 5, 7, 6}, {3, 6, 7, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 3, 5, 4};

static const int DatOldFaceNewFace[][REFTETRABIS15MAXN_nFpoF] =
  { {0, 1}, {2}, {5, 4, 3}, {6, 7} };

static const int DatOldFaceNewFaceLen[] =
  {2, 1, 3, 2};

static const int DatNewEdgeOldEdge[] =
  {0, 1, 1, 2, 3, 4, 5, 5, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 2, 2, 2, 3, 3, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, 2, -1}, {0, -1, 2, 3}, {-1, -1, 2, 3} };

static const int DatChildTwistIndex[] =
  {1, 2, 0, 1, 2, 0, 1, 2, -1, -1};

// Constructor
TRefTetraBis15Desc::TRefTetraBis15Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TetraBis15;

  //set all numbers
  N_Vertices = 6;
  N_Edges = 12;
  N_Faces = 10;
  N_Children = 3;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 1;
  N_InnerEdges = 0;
  N_InnerFaces = 2;

  // initialize all dimension values
  MaxN_VpC = REFTETRABIS15MAXN_VpC;
  MaxN_CpV = REFTETRABIS15MAXN_CpV;
  MaxN_EpC = REFTETRABIS15MAXN_EpC;
  MaxN_CpE = REFTETRABIS15MAXN_CpE;
  MaxN_EpV = REFTETRABIS15MAXN_EpV;
  MaxN_EpF = REFTETRABIS15MAXN_EpF;
  MaxN_FpE = REFTETRABIS15MAXN_FpE;
  MaxN_VpF = REFTETRABIS15MAXN_VpF;
  MaxN_FpV = REFTETRABIS15MAXN_FpV;
  MaxN_FpC = REFTETRABIS15MAXN_FpC;
  MaxN_CpF = REFTETRABIS15MAXN_CpF;
  MaxN_iVpE = REFTETRABIS15MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS15MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS15MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS15MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS15MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS15MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS15MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS15MAXN_nFpoF;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;
  FaceType = (const Refinements *) DatFaceType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildFace = (const int *) DatChildFace;
  FaceChild = (const int *) DatFaceChild;
  FaceChildIndex = (const int *) DatFaceChildIndex;
  FaceChildLen = (const int *) DatFaceChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  FaceVertex = (const int *) DatFaceVertex;
  VertexFace = (const int *) DatVertexFace;
  VertexFaceIndex = (const int *) DatVertexFaceIndex;
  VertexFaceLen = (const int *) DatVertexFaceLen;

  FaceEdge = (const int *) DatFaceEdge;
  EdgeFace = (const int *) DatEdgeFace;
  EdgeFaceIndex = (const int *) DatEdgeFaceIndex;
  EdgeFaceLen = (const int *) DatEdgeFaceLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  NewFaceEqOldFace = (const int *) DatNewFaceEqOldFace;
  NewFaceEqOldFaceIndex = (const int *) DatNewFaceEqOldFaceIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorFaceOfCell = (const int *) DatInteriorFaceOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
  InteriorEdgeOfFace = (const int *) DatInteriorEdgeOfFace;
  InteriorEdgeOfFaceLen = (const int *) DatInteriorEdgeOfFaceLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;
  OldEdgeNewVertexLen = (const int *) DatOldEdgeNewVertexLen;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewEdgeLen = (const int *) DatOldEdgeNewEdgeLen;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
  NewFaceOldFace = (const int *) DatNewFaceOldFace;

  OldFaceNewVertex = (const int *) DatOldFaceNewVertex;
  OldFaceNewVertexPos = (const double *) DatOldFaceNewVertexPos;
  OldFaceNewVertexLen = (const int *) DatOldFaceNewVertexLen;
  OldFaceNewEdge = (const int *) DatOldFaceNewEdge;
  OldFaceNewEdgeLen = (const int *) DatOldFaceNewEdgeLen;
  OldFaceNewFace = (const int *) DatOldFaceNewFace;
  OldFaceNewFaceLen = (const int *) DatOldFaceNewFaceLen;
  OldFaceNewLocFace = (const int *) DatOldFaceNewLocFace;
  ChildTwistIndex = (const int *) DatChildTwistIndex;
}

// Methods
