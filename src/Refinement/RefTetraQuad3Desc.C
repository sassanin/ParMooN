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

#include <RefTetraQuad3Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, NoRef, LineReg, LineReg, NoRef, LineReg};

static const Refinements DatFaceType [] = {TriBis2, TriBis0, TriBis2, TriReg};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,4,6},{4,1,2,5},{6,1,5,3},{4,5,6,1}};
static const int DatVertexChild[][4] = {{0},{0,1,2,3},{1},{2},{0,1,3},{1,2,3},{0,2,3}};
static const int DatVertexChildIndex[][4] = {{0},{1,1,1,3},{2},{3},{2,0,0},{3,2,1},{3,0,2}};
static const int DatVertexChildLen[] = {1,4,1,1,3,3,3};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,2,12,6},{1,10,4,7},{11,3,5,8},{9,10,11,12}};
static const int DatFaceChild[][2] = {{0},{1},{0},{2},{1},{2},{0},{1},{2},{3},{1,3},{2,3},{0,3}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{1},{2},{2},{3},{3},{3},{0},{1,1},{0,2},{2,3}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,1,1,2,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,12,3,4,14,11},{12,1,2,9,13,7},{14,13,10,5,6,8},{9,10,11,12,13,14}};
static const int DatEdgeChild[][3] = {{0},{1},{1},{0},{0},{2},{2},{1},{2},{1,3},{2,3},{0,3},{0,1,3},{1,2,3},{0,2,3}};
static const int DatEdgeChildIndex[][3] = {{0},{1},{2},{2},{3},{3},{4},{5},{5},{3,0},{2,1},{5,2},{1,0,3},{4,1,4},{4,0,5}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,1,1,1,1,2,2,2,3,3,3};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,2},{2,4},{4,0},{0,6},{6,3},{1,3},{2,5},{5,3},{4,5},{5,6},{6,4},{4,1},{5,1},{6,1}};
static const int DatVertexEdge[][6] = {{0,3,4},{0,1,6,12,13,14},{1,2,7},{5,6,8},{2,3,9,11,12},{7,8,9,10,13},{4,5,10,11,14}};
static const int DatVertexEdgeIndex[][6] = {{0,1,0},{1,0,0,1,1,1},{1,0,0},{1,1,1},{1,0,0,1,0},{1,0,1,0,0},{1,0,1,0,0}};
static const int DatVertexEdgeLen[] = {3,6,3,3,5,5,5};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,4},{4,1,2},{0,6,1},{6,3,1},{2,1,5},{5,1,3},{0,4,6},{4,2,5},{6,5,3},{4,5,6},{4,5,1},{5,6,1},{6,4,1}};
static const int DatVertexFace[][9] = {{0,2,6},{0,1,2,3,4,5,10,11,12},{1,4,7},{3,5,8},{0,1,6,7,9,10,12},{4,5,7,8,9,10,11},{2,3,6,8,9,11,12}};
static const int DatVertexFaceIndex[][9] = {{0,0,0},{1,1,2,2,1,1,2,2,2},{2,0,1},{1,2,2},{2,0,1,0,0,0,1},{2,0,2,1,1,1,0},{1,0,2,0,2,1,0}};
static const int DatVertexFaceLen[] = {3,9,3,3,7,7,7};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,12,3},{12,1,2},{4,14,0},{5,6,14},{1,13,7},{13,6,8},{3,11,4},{2,7,9},{10,8,5},{9,10,11},{9,13,12},{10,14,13},{11,12,14}};
static const int DatEdgeFace[][4] = {{0,2},{1,4},{1,7},{0,6},{2,6},{3,8},{3,5},{4,7},{5,8},{7,9,10},{8,9,11},{6,9,12},{0,1,10,12},{4,5,10,11},{2,3,11,12}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{1,0},{2,0},{2,0},{0,2},{0,2},{1,1},{2,1},{2,1},{2,0,0},{0,1,0},{1,2,0},{1,0,2,1},{1,0,1,2},{1,2,1,2}};
static const int DatEdgeFaceLen[] = {2,2,2,2,2,2,2,2,2,3,3,3,4,4,4};

/*
 * New index equal to old index
 */
static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {-1};
static const int DatNewFaceEqOldFaceIndex[] = {-1};

static const int DatInteriorFaceOfCell[] = {10,11,12};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRAQUAD3MAXN_iVpE] =
  { {}, {}, {4}, {6}, {}, {5} };

static const int DatInteriorVertexOfEdgeLen[] = {0, 0, 1, 1, 0, 1};

static const int DatInteriorEdgeOfFace[][REFTETRAQUAD3MAXN_iEpF] =
  { {12}, {14}, {13}, {9, 10, 11} };

static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 1, 1, 3};

/*
 * Old-New Relations
 */
static const int DatOldEdgeNewVertex[][REFTETRAQUAD3MAXN_nVpoE] =
  { {0, 1}, {1, 2}, {2, 4, 0}, {0, 6, 3}, {1, 3}, {2, 5, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 2, 3, 3, 2, 3};

static const int DatOldFaceNewVertex[][REFTETRAQUAD3MAXN_nVpoF] =
  { {0, 1, 2, 4}, {0, 3, 1, 6}, {2, 1, 3, 5}, {0, 2, 3, 4, 5, 6} };

static const int DatOldFaceNewVertexLen[] =
  { 4, 4, 4, 6};

static const double DatOldFaceNewVertexPos[][REFTETRAQUAD3MAXN_nVpoF][REFTETRAQUAD3MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5} }
  };

static const int DatOldEdgeNewEdge[][REFTETRAQUAD3MAXN_nEpoE] =
  { {0}, {1}, {2,3}, {4,5}, {6}, {7,8} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 1, 2, 2, 1, 2};

static const int DatOldFaceNewEdge[][REFTETRAQUAD3MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 5, 6, 0}, {1, 6, 8, 7}, {3, 2, 7, 8, 5, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 4, 4, 6};

static const int DatOldFaceNewFace[][REFTETRAQUAD3MAXN_nFpoF] =
  { {1, 0}, {2, 3}, {5, 4}, {6, 7, 8, 9} };

static const int DatOldFaceNewFaceLen[] =
  {2, 2, 2, 4};

/*
 * New-Old Relations
 */
static const int DatNewEdgeOldEdge[] =
  {0, 1, 2, 2, 3, 3, 4, 5, 5, -1, -1, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 1, 2, 2, 3, 3, 3, 3, -1, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, -1, 3}, {0, -1, 2, 3}, {-1, 1, 2, 3}, {-1, -1, -1, 0} };

static const int DatChildTwistIndex[] =
  {0, 2, 0, 1, 0, 2, 0, 1, 2, 0, -1, -1, -1};


// Constructor
TRefTetraQuad3Desc::TRefTetraQuad3Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TetraQuad3;

  //set all numbers
  N_Vertices = 7;
  N_Edges = 15;
  N_Faces = 13;
  N_Children = 4;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 0;
  N_InnerEdges = 0;
  N_InnerFaces = 3;

  // initialize all dimension values
  MaxN_VpC = REFTETRAQUAD3MAXN_VpC;
  MaxN_CpV = REFTETRAQUAD3MAXN_CpV;
  MaxN_EpC = REFTETRAQUAD3MAXN_EpC;
  MaxN_CpE = REFTETRAQUAD3MAXN_CpE;
  MaxN_EpV = REFTETRAQUAD3MAXN_EpV;
  MaxN_EpF = REFTETRAQUAD3MAXN_EpF;
  MaxN_FpE = REFTETRAQUAD3MAXN_FpE;
  MaxN_VpF = REFTETRAQUAD3MAXN_VpF;
  MaxN_FpV = REFTETRAQUAD3MAXN_FpV;
  MaxN_FpC = REFTETRAQUAD3MAXN_FpC;
  MaxN_CpF = REFTETRAQUAD3MAXN_CpF;
  MaxN_iVpE = REFTETRAQUAD3MAXN_iVpE;
  MaxN_iEpF = REFTETRAQUAD3MAXN_iEpF;
  MaxN_nVpoE = REFTETRAQUAD3MAXN_nVpoE;
  MaxN_nEpoE = REFTETRAQUAD3MAXN_nEpoE;
  MaxN_nVpoF = REFTETRAQUAD3MAXN_nVpoF;
  MaxN_oVpoF = REFTETRAQUAD3MAXN_oVpoF;
  MaxN_nEpoF = REFTETRAQUAD3MAXN_nEpoF;
  MaxN_nFpoF = REFTETRAQUAD3MAXN_nFpoF;

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
