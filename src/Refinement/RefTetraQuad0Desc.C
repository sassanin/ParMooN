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

#include <RefTetraQuad0Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {LineReg, LineReg, LineReg, NoRef, NoRef, NoRef};

static const Refinements DatFaceType [] = {TriReg, TriBis2, TriBis0, TriBis0};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,4,6,3},{4,1,5,3},{6,5,2,3},{4,5,6,3}};
static const int DatVertexChild[][4] = {{0},{1},{2},{0,1,2,3},{0,1,3},{1,2,3},{0,2,3}};
static const int DatVertexChildIndex[][4] = {{0},{1},{2},{3,3,3,3},{1,0,0},{2,1,1},{2,0,2}};
static const int DatVertexChildLen[] = {1,1,1,4,3,3,3};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,4,12,9},{1,5,6,10},{2,11,7,8},{3,10,11,12}};
static const int DatFaceChild[][2] = {{0},{1},{2},{3},{0},{1},{1},{2},{2},{0},{1,3},{2,3},{0,3}};
static const int DatFaceChildIndex[][2] = {{0},{0},{0},{0},{1},{1},{2},{2},{3},{3},{3,1},{1,2},{2,3}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,1,1,2,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,11,5,6,12,14},{1,2,9,12,7,13},{10,3,4,14,13,8},{9,10,11,12,13,14}};
static const int DatEdgeChild[][3] = {{0},{1},{1},{2},{2},{0},{0},{1},{2},{1,3},{2,3},{0,3},{0,1,3},{1,2,3},{0,2,3}};
static const int DatEdgeChildIndex[][3] = {{0},{0},{1},{1},{2},{2},{3},{4},{5},{2,0},{0,1},{1,2},{4,3,3},{5,4,4},{5,3,5}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,1,1,1,1,2,2,2,3,3,3};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,4},{4,1},{1,5},{5,2},{2,6},{6,0},{0,3},{1,3},{2,3},{4,5},{5,6},{6,4},{4,3},{5,3},{6,3}};
static const int DatVertexEdge[][6] = {{0,5,6},{1,2,7},{3,4,8},{6,7,8,12,13,14},{0,1,9,11,12},{2,3,9,10,13},{4,5,10,11,14}};
static const int DatVertexEdgeIndex[][6] = {{0,1,0},{1,0,0},{1,0,0},{1,1,1,1,1,1},{1,0,0,1,0},{1,0,1,0,0},{1,0,1,0,0}};
static const int DatVertexEdgeLen[] = {3,3,3,6,5,5,5};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,4,6},{4,1,5},{6,5,2},{4,5,6},{0,3,4},{4,3,1},{5,1,3},{2,5,3},{6,2,3},{0,6,3},{4,5,3},{5,6,3},{6,4,3}};
static const int DatVertexFace[][9] = {{0,4,9},{1,5,6},{2,7,8},{4,5,6,7,8,9,10,11,12},{0,1,3,4,5,10,12},{1,2,3,6,7,10,11},{0,2,3,8,9,11,12}};
static const int DatVertexFaceIndex[][9] = {{0,0,0},{1,2,1},{2,0,1},{1,1,2,2,2,2,2,2,2},{1,0,0,2,0,0,1},{2,1,1,0,1,1,0},{2,0,2,0,1,1,0}};
static const int DatVertexFaceLen[] = {3,3,3,9,7,7,7};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,11,5},{1,2,9},{10,3,4},{9,10,11},{6,12,0},{12,7,1},{2,7,13},{3,13,8},{4,8,14},{5,14,6},{9,13,12},{10,14,13},{11,12,14}};
static const int DatEdgeFace[][4] = {{0,4},{1,5},{1,6},{2,7},{2,8},{0,9},{4,9},{5,6},{7,8},{1,3,10},{2,3,11},{0,3,12},{4,5,10,12},{6,7,10,11},{8,9,11,12}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{0,2},{1,0},{1,0},{2,0},{2,0},{0,2},{1,1},{2,1},{2,0,0},{0,1,0},{1,2,0},{1,0,2,1},{2,1,1,2},{2,1,1,2}};
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

static const int DatInteriorVertexOfEdge[][REFTETRAQUAD0MAXN_iVpE] =
  { {4}, {5}, {6}, {}, {}, {} };

static const int DatInteriorVertexOfEdgeLen[] = {1, 1, 1, 0, 0, 0};

static const int DatInteriorEdgeOfFace[][REFTETRAQUAD0MAXN_iEpF] =
  { {9,10,11}, {12}, {13}, {14} };

static const int DatInteriorEdgeOfFaceLen[] =
  { 3, 1, 1, 1};

/*
 * Old-New Relations
 */
static const int DatOldEdgeNewVertex[][REFTETRAQUAD0MAXN_nVpoE] =
  { {0, 4, 1}, {1, 5, 2}, {2, 6, 0}, {0, 3}, {1, 3}, {2, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {3, 3, 3, 2, 2, 2};

static const int DatOldFaceNewVertex[][REFTETRAQUAD0MAXN_nVpoF] =
  { {0, 1, 2, 4, 5, 6}, {0, 3, 1, 4}, {2, 1, 3, 5}, {0, 2, 3, 6} };

static const int DatOldFaceNewVertexLen[] =
  { 6, 4, 4, 4};

static const double DatOldFaceNewVertexPos[][REFTETRAQUAD0MAXN_nVpoF][REFTETRAQUAD0MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0} }
  };

static const int DatOldEdgeNewEdge[][REFTETRAQUAD0MAXN_nEpoE] =
  { {0,1}, {2,3}, {4,5}, {6}, {7}, {8} };

static const int DatOldEdgeNewEdgeLen[] =
  {2, 2, 2, 1, 1, 1};

static const int DatOldFaceNewEdge[][REFTETRAQUAD0MAXN_nEpoF] =
  { {0, 1, 2, 3, 4, 5}, {6, 7, 1, 0}, {3, 2, 7, 8}, {5, 4, 8, 6} };

static const int DatOldFaceNewEdgeLen[] =
  {6, 4, 4, 4};

static const int DatOldFaceNewFace[][REFTETRAQUAD0MAXN_nFpoF] =
  { {0, 1, 2, 3}, {5, 4}, {7, 6}, {9, 8} };

static const int DatOldFaceNewFaceLen[] =
  {4, 2, 2, 2};

/*
 * New-Old Relations
 */
static const int DatNewEdgeOldEdge[] =
  {0, 0, 1, 1, 2, 2, 3, 4, 5, -1, -1, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, -1, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, -1, 3}, {0, 1, 2, -1}, {0, -1, 2, 3}, {0, -1, -1, -1} };

static const int DatChildTwistIndex[] =
  {0, 1, 2, 0, 0, 2, 1, 0, 1, 0, -1, -1, -1};

// Constructor
TRefTetraQuad0Desc::TRefTetraQuad0Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TetraQuad0;

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
  MaxN_VpC = REFTETRAQUAD0MAXN_VpC;
  MaxN_CpV = REFTETRAQUAD0MAXN_CpV;
  MaxN_EpC = REFTETRAQUAD0MAXN_EpC;
  MaxN_CpE = REFTETRAQUAD0MAXN_CpE;
  MaxN_EpV = REFTETRAQUAD0MAXN_EpV;
  MaxN_EpF = REFTETRAQUAD0MAXN_EpF;
  MaxN_FpE = REFTETRAQUAD0MAXN_FpE;
  MaxN_VpF = REFTETRAQUAD0MAXN_VpF;
  MaxN_FpV = REFTETRAQUAD0MAXN_FpV;
  MaxN_FpC = REFTETRAQUAD0MAXN_FpC;
  MaxN_CpF = REFTETRAQUAD0MAXN_CpF;
  MaxN_iVpE = REFTETRAQUAD0MAXN_iVpE;
  MaxN_iEpF = REFTETRAQUAD0MAXN_iEpF;
  MaxN_nVpoE = REFTETRAQUAD0MAXN_nVpoE;
  MaxN_nEpoE = REFTETRAQUAD0MAXN_nEpoE;
  MaxN_nVpoF = REFTETRAQUAD0MAXN_nVpoF;
  MaxN_oVpoF = REFTETRAQUAD0MAXN_oVpoF;
  MaxN_nEpoF = REFTETRAQUAD0MAXN_nEpoF;
  MaxN_nFpoF = REFTETRAQUAD0MAXN_nFpoF;

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
