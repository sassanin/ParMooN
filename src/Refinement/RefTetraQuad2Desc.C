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

#include <RefTetraQuad2Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, LineReg, NoRef, NoRef, LineReg, LineReg};

static const Refinements DatFaceType [] = {TriBis1, TriBis1, TriReg, TriBis1};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,4,5},{0,5,6,3},{0,4,2,6},{4,5,6,0}};
static const int DatVertexChild[][4] = {{0,1,2,3},{0},{2},{1},{0,2,3},{0,1,3},{1,2,3}};
static const int DatVertexChildIndex[][4] = {{0,0,0,3},{1},{2},{3},{2,1,0},{3,1,1},{2,3,2}};
static const int DatVertexChildLen[] = {4,1,1,1,3,3,3};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,2,5,10},{11,3,6,9},{1,12,4,8},{7,10,11,12}};
static const int DatFaceChild[][2] = {{0},{2},{0},{1},{2},{0},{1},{3},{2},{1},{0,3},{1,3},{2,3}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{1},{2},{2},{2},{0},{3},{3},{3,1},{0,2},{1,3}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,1,1,2,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,1,12,13,5,9},{13,10,14,4,6,8},{12,2,3,14,11,7},{9,10,11,12,13,14}};
static const int DatEdgeChild[][3] = {{0},{0},{2},{2},{1},{0},{1},{2},{1},{0,3},{1,3},{2,3},{0,2,3},{0,1,3},{1,2,3}};
static const int DatEdgeChildIndex[][3] = {{0},{1},{1},{2},{3},{4},{4},{5},{5},{5,0},{1,1},{4,2},{2,0,3},{3,0,4},{2,3,5}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,1,1,1,1,2,2,2,3,3,3};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,4},{4,2},{2,0},{0,3},{1,5},{5,3},{2,6},{6,3},{4,5},{5,6},{6,4},{4,0},{5,0},{6,0}};
static const int DatVertexEdge[][6] = {{0,3,4,12,13,14},{0,1,5},{2,3,7},{4,6,8},{1,2,9,11,12},{5,6,9,10,13},{7,8,10,11,14}};
static const int DatVertexEdgeIndex[][6] = {{0,1,0,1,1,1},{1,0,0},{1,0,0},{1,1,1},{1,0,0,1,0},{1,0,1,0,0},{1,0,1,0,0}};
static const int DatVertexEdgeLen[] = {6,3,3,3,5,5,5};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,4},{0,4,2},{0,5,1},{0,3,5},{2,4,6},{4,1,5},{6,5,3},{4,5,6},{0,2,6},{0,6,3},{4,5,0},{5,6,0},{6,4,0}};
static const int DatVertexFace[][9] = {{0,1,2,3,8,9,10,11,12},{0,2,5},{1,4,8},{3,6,9},{0,1,4,5,7,10,12},{2,3,5,6,7,10,11},{4,6,7,8,9,11,12}};
static const int DatVertexFaceIndex[][9] = {{0,0,0,0,0,0,2,2,2},{1,2,1},{2,0,1},{1,2,2},{2,1,1,0,0,0,1},{1,2,2,1,1,1,0},{2,0,2,2,1,1,0}};
static const int DatVertexFaceLen[] = {9,3,3,3,7,7,7};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,1,12},{12,2,3},{13,5,0},{4,6,13},{2,11,7},{1,5,9},{10,6,8},{9,10,11},{3,7,14},{14,8,4},{9,13,12},{10,14,13},{11,12,14}};
static const int DatEdgeFace[][4] = {{0,2},{0,5},{1,4},{1,8},{3,9},{2,5},{3,6},{4,8},{6,9},{5,7,10},{6,7,11},{4,7,12},{0,1,10,12},{2,3,10,11},{8,9,11,12}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{1,0},{1,0},{2,0},{0,2},{1,1},{1,1},{2,1},{2,1},{2,0,0},{0,1,0},{1,2,0},{2,0,2,1},{0,2,1,2},{2,0,1,2}};
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

static const int DatInteriorVertexOfEdge[][REFTETRAQUAD2MAXN_iVpE] =
  { {}, {4}, {}, {}, {5}, {6} };

static const int DatInteriorVertexOfEdgeLen[] = {0, 1, 0, 0, 1, 1};

static const int DatInteriorEdgeOfFace[][REFTETRAQUAD2MAXN_iEpF] =
  { {12}, {13}, {9,10,11}, {14} };

static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 1, 3, 1};

/*
 * Old-New Relations
 */
static const int DatOldEdgeNewVertex[][REFTETRAQUAD2MAXN_nVpoE] =
  { {0, 1}, {1, 4, 2}, {2, 0}, {0, 3}, {1, 5, 3}, {2, 6, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 3, 2, 2, 3, 3};

static const int DatOldFaceNewVertex[][REFTETRAQUAD2MAXN_nVpoF] =
  { {0, 1, 2, 4}, {0, 3, 1, 5}, {2, 1, 3, 4, 5, 6}, {0, 2, 3, 6} };

static const int DatOldFaceNewVertexLen[] =
  { 4, 4, 6, 4};

static const double DatOldFaceNewVertexPos[][REFTETRAQUAD2MAXN_nVpoF][REFTETRAQUAD2MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} }
  };

static const int DatOldEdgeNewEdge[][REFTETRAQUAD2MAXN_nEpoE] =
  { {0}, {1, 2}, {3}, {4}, {5, 6}, {7, 8} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 2, 1, 1, 2, 2};

static const int DatOldFaceNewEdge[][REFTETRAQUAD2MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 6, 5, 0}, {2, 1, 5, 6, 9, 7}, {3, 7, 8, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 4, 6, 4};

static const int DatOldFaceNewFace[][REFTETRAQUAD2MAXN_nFpoF] =
  { {0, 1}, {3, 2}, {4, 5, 6, 7}, {8, 9} };

static const int DatOldFaceNewFaceLen[] =
  {2, 2, 4, 2};

/*
 * New-Old Relations
 */
static const int DatNewEdgeOldEdge[] =
  {0, 1, 1, 2, 3, 4, 4, 5, 5, -1, -1, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 1, 2, 2, 2, 2, 3, 3, -1, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, 2, -1}, {-1, 1, 2, 3}, {0, -1, 2, 3}, {-1, -1, 0, -1} };

static const int DatChildTwistIndex[] =
  {1, 2, 2, 1, 0, 1, 2, 0, 1, 2, -1, -1, -1};

// Constructor
TRefTetraQuad2Desc::TRefTetraQuad2Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TetraQuad2;

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
  MaxN_VpC = REFTETRAQUAD2MAXN_VpC;
  MaxN_CpV = REFTETRAQUAD2MAXN_CpV;
  MaxN_EpC = REFTETRAQUAD2MAXN_EpC;
  MaxN_CpE = REFTETRAQUAD2MAXN_CpE;
  MaxN_EpV = REFTETRAQUAD2MAXN_EpV;
  MaxN_EpF = REFTETRAQUAD2MAXN_EpF;
  MaxN_FpE = REFTETRAQUAD2MAXN_FpE;
  MaxN_VpF = REFTETRAQUAD2MAXN_VpF;
  MaxN_FpV = REFTETRAQUAD2MAXN_FpV;
  MaxN_FpC = REFTETRAQUAD2MAXN_FpC;
  MaxN_CpF = REFTETRAQUAD2MAXN_CpF;
  MaxN_iVpE = REFTETRAQUAD2MAXN_iVpE;
  MaxN_iEpF = REFTETRAQUAD2MAXN_iEpF;
  MaxN_nVpoE = REFTETRAQUAD2MAXN_nVpoE;
  MaxN_nEpoE = REFTETRAQUAD2MAXN_nEpoE;
  MaxN_nVpoF = REFTETRAQUAD2MAXN_nVpoF;
  MaxN_oVpoF = REFTETRAQUAD2MAXN_oVpoF;
  MaxN_nEpoF = REFTETRAQUAD2MAXN_nEpoF;
  MaxN_nFpoF = REFTETRAQUAD2MAXN_nFpoF;

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
