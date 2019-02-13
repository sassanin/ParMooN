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

#include <RefTetraBis05Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {LineReg, NoRef, NoRef, NoRef, NoRef, LineReg};

static const Refinements DatFaceType [] = {TriBis0, TriBis2, TriBis2, TriBis1};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,4,2,5},{4,1,2,5},{0,4,5,3},{4,1,5,3}};
static const int DatVertexChild[][4] = {{0,2},{1,3},{0,1},{2,3},{0,1,2,3},{0,1,2,3}};
static const int DatVertexChildIndex[][4] = {{0,0},{1,1},{2,2},{3,3},{1,0,1,0},{3,3,2,2}};
static const int DatVertexChildLen[] = {2,2,2,2,4,4};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,10,8,6},{1,11,4,8},{10,2,9,7},{11,3,5,9}};
static const int DatFaceChild[][2] = {{0},{1},{2},{3},{1},{3},{0},{2},{0,1},{2,3},{0,2},{1,3}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{1},{2},{2},{3},{3},{2,3},{2,3},{1,0},{1,0}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,2,2,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,8,3,11,12,6},{1,2,8,12,10,6},{0,12,11,4,9,7},{1,10,12,9,5,7}};
static const int DatEdgeChild[][4] = {{0,2},{1,3},{1},{0},{2},{3},{0,1},{2,3},{0,1},{2,3},{1,3},{0,2},{0,1,2,3}};
static const int DatEdgeChildIndex[][4] = {{0,0},{0,0},{1},{2},{3},{4},{5,5},{5,5},{1,2},{4,3},{4,1},{3,2},{4,3,1,2}};
static const int DatEdgeChildLen[] = {2,2,1,1,1,1,2,2,2,2,2,2,4};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,4},{4,1},{1,2},{2,0},{0,3},{1,3},{2,5},{5,3},{4,2},{4,3},{1,5},{0,5},{4,5}};
static const int DatVertexEdge[][5] = {{0,3,4,11},{1,2,5,10},{2,3,6,8},{4,5,7,9},{0,1,8,9,12},{6,7,10,11,12}};
static const int DatVertexEdgeIndex[][5] = {{0,1,0,0},{1,0,0,0},{1,0,0,1},{1,1,1,1},{1,0,0,0,0},{1,0,1,1,1}};
static const int DatVertexEdgeLen[] = {4,4,4,4,5,5};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,4,2},{4,1,2},{0,3,4},{4,3,1},{2,1,5},{5,1,3},{0,2,5},{0,5,3},{2,4,5},{5,4,3},{0,5,4},{4,5,1}};
static const int DatVertexFace[][8] = {{0,2,6,7,10},{1,3,4,5,11},{0,1,4,6,8},{2,3,5,7,9},{0,1,2,3,8,9,10,11},{4,5,6,7,8,9,10,11}};
static const int DatVertexFaceIndex[][8] = {{0,0,0,0,0},{1,2,1,1,2},{2,2,0,1,0},{1,1,2,2,2},{1,0,2,0,1,1,2,0},{2,0,2,1,2,0,1,1}};
static const int DatVertexFaceLen[] = {5,5,5,5,8,8};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,8,3},{1,2,8},{4,9,0},{9,5,1},{2,10,6},{10,5,7},{3,6,11},{11,7,4},{8,12,6},{12,9,7},{11,12,0},{12,10,1}};
static const int DatEdgeFace[][4] = {{0,2,10},{1,3,11},{1,4},{0,6},{2,7},{3,5},{4,6,8},{5,7,9},{0,1,8},{2,3,9},{4,5,11},{6,7,10},{8,9,10,11}};
static const int DatEdgeFaceIndex[][4] = {{0,2,2},{0,2,2},{1,0},{2,0},{0,2},{1,1},{2,1,2},{2,1,2},{1,2,0},{1,0,1},{1,0,1},{2,0,0},{1,0,1,0}};
static const int DatEdgeFaceLen[] = {3,3,2,2,2,2,3,3,3,3,3,3,4};

static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {0};
static const int DatNewFaceEqOldFaceIndex[] = {0};

static const int DatOldFaceNewVertex[][REFTETRABIS05MAXN_nVpoF] =
  { {0, 1, 2, 4}, {0, 3, 1, 4}, {2, 1, 3, 5}, {0, 2, 3, 5} };
static const int DatOldFaceNewVertexLen[] =
  { 4, 4, 4, 4};
static const double DatOldFaceNewVertexPos[][REFTETRABIS05MAXN_nVpoF][REFTETRABIS05MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} } };

static const int DatInteriorFaceOfCell[] = {8, 9, 10, 11};
static const int DatInteriorEdgeOfCell[] = {12};

static const int DatInteriorVertexOfEdge[][REFTETRABIS05MAXN_iVpE] =
  { {4}, {}, {}, {}, {}, {5} };
static const int DatInteriorVertexOfEdgeLen[] = {1, 0, 0, 0, 0, 1};

static const int DatInteriorEdgeOfFace[][REFTETRABIS05MAXN_iEpF] =
  { {8}, {9}, {10}, {11} };
static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 1, 1, 1};

static const int DatOldEdgeNewVertex[][REFTETRABIS05MAXN_nVpoE] =
  { {0, 4, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 5, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {3, 2, 2, 2, 2, 3};

static const int DatOldEdgeNewEdge[][REFTETRABIS05MAXN_nEpoE] =
  { {0, 1}, {2}, {3}, {4}, {5}, {6, 7} };

static const int DatOldEdgeNewEdgeLen[] =
  {2, 1, 1, 1, 1, 2};

static const int DatOldFaceNewEdge[][REFTETRABIS05MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 5, 1, 0}, {2, 5, 7, 6}, {3, 6, 7, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 4, 4, 4};

static const int DatOldFaceNewFace[][REFTETRABIS05MAXN_nFpoF] =
  { {0, 1}, {3, 2}, {5, 4}, {6, 7} };

static const int DatOldFaceNewFaceLen[] =
  {2, 2, 2, 2};

static const int DatNewEdgeOldEdge[] =
  {0, 0, 1, 2, 3, 4, 5, 5, -1, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 1, 2, 2, 3, 3, -1, -1, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, -1, -1, 3}, {0, -1, 2, -1}, {-1, 1, -1, 3}, {-1, 1, 2, -1} };

static const int DatChildTwistIndex[] =
  {0, 1, 0, 2, 0, 2, 1, 2, -1, -1, -1, -1};

// Constructor
TRefTetraBis05Desc::TRefTetraBis05Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TetraBis05;

  //set all numbers
  N_Vertices = 6;
  N_Edges = 13;
  N_Faces = 12;
  N_Children = 4;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 0;
  N_InnerEdges = 1;
  N_InnerFaces = 4;

  // initialize all dimension values
  MaxN_VpC = REFTETRABIS05MAXN_VpC;
  MaxN_CpV = REFTETRABIS05MAXN_CpV;
  MaxN_EpC = REFTETRABIS05MAXN_EpC;
  MaxN_CpE = REFTETRABIS05MAXN_CpE;
  MaxN_EpV = REFTETRABIS05MAXN_EpV;
  MaxN_EpF = REFTETRABIS05MAXN_EpF;
  MaxN_FpE = REFTETRABIS05MAXN_FpE;
  MaxN_VpF = REFTETRABIS05MAXN_VpF;
  MaxN_FpV = REFTETRABIS05MAXN_FpV;
  MaxN_FpC = REFTETRABIS05MAXN_FpC;
  MaxN_CpF = REFTETRABIS05MAXN_CpF;
  MaxN_iVpE = REFTETRABIS05MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS05MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS05MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS05MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS05MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS05MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS05MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS05MAXN_nFpoF;

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
