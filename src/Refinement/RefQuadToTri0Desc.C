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
// @(#)RefQuadToTri0Desc.C        1.4 10/19/99
//
// Class:       TRefQuadToTri0Desc
// Purpose:     refinement descriptor for refinement of a quadrangle
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//
// =======================================================================

#include <RefQuadToTri0Desc.h>

static const Shapes DatChildType[] = { Triangle, Triangle};

static const Refinements DatEdgeType[] = { NoRef, NoRef, NoRef, NoRef};

static const int DatChildVertex[][QUADTTMAXN_VpC] =
                 { {0, 1, 2},  {2, 3, 0}};
static const int DatVertexChild[][QUADTTMAXN_CpV] =
                 { {0, 1},  {0},  {0, 1},  {1}};
static const int DatVertexChildIndex[][QUADTTMAXN_CpV] =
               { {0, 2},  {1},  {2, 0},  {1}};
static const int DatVertexChildLen[] =
               { 2,  1,  2,  1};

static const int DatChildEdge[][QUADTTMAXN_EpC] =
               { {0, 1, 4},  {2, 3, 4}};
static const int DatEdgeChild[][QUADTTMAXN_CpE] =
               { {0},  {0},  {1},  {1},  {0, 1}};
static const int DatEdgeChildIndex[][QUADTTMAXN_CpE] = 
               { {0},  {1},  {0},  {1},  {2, 2}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  2};

static const int DatEdgeVertex[][2] =
               { {0, 1},  {1, 2},  {2, 3},  {3, 0},  {0, 2}};
static const int DatVertexEdge[][QUADTTMAXN_EpV] = 
               { {0, 3, 4},  {0, 1},  {1, 2, 4},  {2, 3}};
static const int DatVertexEdgeIndex[][QUADTTMAXN_EpV] =
               { {0, 1, 0},  {1, 0},  {1, 0, 1},  {1, 0}};
static const int DatVertexEdgeLen[] =
               { 3,  2,  3,  2};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2, 3};

static const int DatNewEdgeEqOldEdge[] = { 0, 1, 2, 3};
static const int DatNewEdgeEqOldEdgeIndex[] = { 0, 1, 2, 3};

static const int DatInteriorEdgeOfCell[] = { 4};

static const int DatInteriorVertexOfEdgeLen[] = { 0,  0,  0,  0};

static const int DatOldEdgeNewVertex[][QUADTTMAXN_nVpoE] =
               { {0, 1},  {1, 2},  {2, 3},  {3, 0}};

static const int DatOldEdgeNewEdge[][QUADTTMAXN_nEpoE] =
               { {0},  {1},  {2},  {3}};

// static const int DatOldEdgeNewLocEdge[][QUADTTN_E] =
//                { {0, 1, -1}, {2, 3, -1}};

static const int DatOldEdgeNewLocEdge[][4] =
               { {  0,  1, -1, -1 },
                 { -1, -1,  0,  1} };

static const int DatNewEdgeOldEdge[] =
               { 0,  1,  2,  3,  -1};

// Constructor
TRefQuadToTri0Desc::TRefQuadToTri0Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = QuadToTri0;

  // set all numbers
  N_Edges = 5;
  N_Vertices = 4;
  N_Children = 2;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 4;
  N_InnerEdges = 1;

  // initialize all dimension values
  MaxN_VpC = QUADTTMAXN_VpC;
  MaxN_CpV = QUADTTMAXN_CpV;
  MaxN_EpC = QUADTTMAXN_EpC;
  MaxN_CpE = QUADTTMAXN_CpE;
  MaxN_EpV = QUADTTMAXN_EpV;
  MaxN_nVpoE = QUADTTMAXN_nVpoE;
  MaxN_nEpoE = QUADTTMAXN_nEpoE;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;

  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
}

// Methods
