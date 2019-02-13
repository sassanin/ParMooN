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

#include <RefTetraBis2Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, NoRef, LineReg, NoRef, NoRef, NoRef};

static const Refinements DatFaceType [] = {TriBis2, NoRef, NoRef, TriBis0};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,4,3},{4,1,2,3}};
static const int DatVertexChild[][2] = {{0},{0,1},{1},{0,1},{0,1}};
static const int DatVertexChildIndex[][2] = {{0},{1,1},{2},{3,3},{2,0}};
static const int DatVertexChildLen[] = {1,2,1,2,2};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,2,6,5},{1,6,3,4}};
static const int DatFaceChild[][2] = {{0},{1},{0},{1},{1},{0},{0,1}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{2},{3},{3},{2,1}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,7,3,4,5,8},{7,1,2,8,5,6}};
static const int DatEdgeChild[][2] = {{0},{1},{1},{0},{0},{0,1},{1},{0,1},{0,1}};
static const int DatEdgeChildIndex[][2] = {{0},{1},{2},{2},{3},{4,4},{5},{1,0},{5,3}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,2,1,2,2};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,2},{2,4},{4,0},{0,3},{1,3},{2,3},{1,4},{4,3}};
static const int DatVertexEdge[][4] = {{0,3,4},{0,1,5,7},{1,2,6},{4,5,6,8},{2,3,7,8}};
static const int DatVertexEdgeIndex[][4] = {{0,1,0},{1,0,0,0},{1,0,0},{1,1,1,1},{1,0,1,0}};
static const int DatVertexEdgeLen[] = {3,4,3,4,4};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,4},{4,1,2},{0,3,1},{2,1,3},{4,2,3},{0,4,3},{4,1,3}};
static const int DatVertexFace[][5] = {{0,2,5},{0,1,2,3,6},{1,3,4},{2,3,4,5,6},{0,1,4,5,6}};
static const int DatVertexFaceIndex[][5] = {{0,0,0},{1,1,2,1,1},{2,0,1},{1,2,2,2,2},{2,0,0,1,0}};
static const int DatVertexFaceLen[] = {3,5,3,5,5};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,7,3},{7,1,2},{4,5,0},{1,5,6},{2,6,8},{3,8,4},{7,5,8}};
static const int DatEdgeFace[][3] = {{0,2},{1,3},{1,4},{0,5},{2,5},{2,3,6},{3,4},{0,1,6},{4,5,6}};
static const int DatEdgeFaceIndex[][3] = {{0,2},{1,0},{2,0},{2,0},{0,2},{1,1,1},{2,1},{1,0,0},{2,1,2}};
static const int DatEdgeFaceLen[] = {2,2,2,2,2,3,2,3,3};

/*
 * New index equal to old index
 */
static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {2, 3};
static const int DatNewFaceEqOldFaceIndex[] = {1, 2};

static const int DatOldFaceNewVertex[][REFTETRABIS2MAXN_nVpoF] =
  { {0, 1, 2, 4}, {0, 3, 1}, {2, 1, 3}, {0, 2, 3, 4} };
static const int DatOldFaceNewVertexLen[] =
  { 4, 3, 3, 4};
static const double DatOldFaceNewVertexPos[][REFTETRABIS2MAXN_nVpoF][REFTETRABIS2MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0} } };

static const int DatInteriorFaceOfCell[] = {6};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRABIS2MAXN_iVpE] =
  { {}, {}, {4}, {}, {}, {} };

static const int DatInteriorVertexOfEdgeLen[] = {0, 0, 1, 0, 0, 0};

static const int DatInteriorEdgeOfFace[][REFTETRABIS2MAXN_iEpF] =
  { {7}, {}, {}, {8} };

static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 0, 0, 1};

/*
 * Old-New Relations
 */
static const int DatOldEdgeNewVertex[][REFTETRABIS2MAXN_nVpoE] =
  { {0, 1}, {1, 2}, {2, 4, 0}, {0, 3}, {1, 3}, {2, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 2, 3, 2, 2, 2};

static const int DatOldEdgeNewEdge[][REFTETRABIS2MAXN_nEpoE] =
  { {0}, {1}, {2, 3}, {4}, {5}, {6} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 1, 2, 1, 1, 1};

static const int DatOldFaceNewEdge[][REFTETRABIS2MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 5, 0}, {1, 5, 6}, {3, 2, 6, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 3, 3, 4};

static const int DatOldFaceNewFace[][REFTETRABIS2MAXN_nFpoF] =
  { {1, 0}, {2}, {3}, {5, 4} };

static const int DatOldFaceNewFaceLen[] =
  {2, 1, 1, 2};

/*
 * New-Old Relations
 */
static const int DatNewEdgeOldEdge[] =
  {0, 1, 2, 2, 3, 4, 5, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 2, 3, 3, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, -1, 3}, {0, -1, 2, 3} };

static const int DatChildTwistIndex[] =
  {0, 2, 0, 0, 1, 0, -1};

// Constructor
TRefTetraBis2Desc::TRefTetraBis2Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TetraBis2;

  //set all numbers
  N_Vertices = 5;
  N_Edges = 9;
  N_Faces = 7;
  N_Children = 2;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 2;
  N_InnerEdges = 0;
  N_InnerFaces = 1;

  // initialize all dimension values
  MaxN_VpC = REFTETRABIS2MAXN_VpC;
  MaxN_CpV = REFTETRABIS2MAXN_CpV;
  MaxN_EpC = REFTETRABIS2MAXN_EpC;
  MaxN_CpE = REFTETRABIS2MAXN_CpE;
  MaxN_EpV = REFTETRABIS2MAXN_EpV;
  MaxN_EpF = REFTETRABIS2MAXN_EpF;
  MaxN_FpE = REFTETRABIS2MAXN_FpE;
  MaxN_VpF = REFTETRABIS2MAXN_VpF;
  MaxN_FpV = REFTETRABIS2MAXN_FpV;
  MaxN_FpC = REFTETRABIS2MAXN_FpC;
  MaxN_CpF = REFTETRABIS2MAXN_CpF;
  MaxN_iVpE = REFTETRABIS2MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS2MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS2MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS2MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS2MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS2MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS2MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS2MAXN_nFpoF;

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
