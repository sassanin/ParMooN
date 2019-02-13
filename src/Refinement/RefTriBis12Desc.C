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
// @(#)RefTriBis0Desc.C        1.3 09/13/99
//
// Class:       TRefTriBis0Desc
// Purpose:     bisection of edge 0
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//              Gunar Matthies 03.09.1998
//
// =======================================================================

#include <RefTriBis12Desc.h>

static const Shapes DatChildType[] = { Triangle, Triangle, Triangle};

static const Refinements DatEdgeType[] = { NoRef, LineReg, LineReg };

static const int DatChildVertex[][TRIBI12MAXN_VpC] = {{0,1,3},{2,4,3},{0,3,4}};
static const int DatVertexChild[][TRIBI12MAXN_CpV] = {{0,2},{0},{1},{0,1,2},{1,2}};
static const int DatVertexChildIndex[][TRIBI12MAXN_CpV] = {{0,0},{1},{0},{2,2,1},{1,2}};
static const int DatVertexChildLen[] = {2,1,1,3,2};

static const int DatChildEdge[][TRIBI12MAXN_EpC] = {{0,1,5},{3,6,2},{5,6,4}};
static const int DatEdgeChild[][TRIBI12MAXN_CpE] = {{0},{0},{1},{1},{2},{0,2},{1,2}};
static const int DatEdgeChildIndex[][TRIBI12MAXN_CpE] = {{0},{1},{2},{0},{2},{2,0},{1,1}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,2,2};

static const int DatEdgeVertex[][2] = {{0,1},{1,3},{3,2},{2,4},{4,0},{3,0},{3,4}};
static const int DatVertexEdge[][TRIBI12MAXN_EpV] = {{0,4,5},{0,1},{2,3},{1,2,5,6},{3,4,6}};
static const int DatVertexEdgeIndex[][TRIBI12MAXN_EpV] = {{0,1,1},{1,0},{1,0},{1,0,0,0},{1,0,1}};
static const int DatVertexEdgeLen[] = {3,2,2,4,3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2};

static const int DatNewEdgeEqOldEdge[] = {0};
static const int DatNewEdgeEqOldEdgeIndex[] = {0};

static const int DatInteriorVertexOfEdge[][TRIBI12MAXN_iVpE] = { {}, {3}, {4} };
static const int DatInteriorVertexOfEdgeLen[] = { 0, 1, 1};

static const int DatInteriorEdgeOfCell[] = {5, 6};

static const int DatOldEdgeNewVertex[][TRIBI12MAXN_nVpoE] = { {0, 1}, { 1, 3, 2}, {2, 4, 0}};

static const int DatOldEdgeNewEdge[][TRIBI12MAXN_nEpoE] =
               { {0}, {1, 2}, {3, 4}};

static const int DatOldEdgeNewLocEdge[][TRIBI12N_E] =
               { {0, 1, -1}, {-1, 2, 0}, {-1, -1, 2} };

static const int DatNewEdgeOldEdge[] =
               { 0, 1, 1, 2, 2, -1, -1};

// Constructor
TRefTriBis12Desc::TRefTriBis12Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TriBis12;

  // set all numbers
  N_Edges = 7;
  N_Vertices = 5;
  N_Children = 3;
  N_NewVertEqOldVert = 3;
  N_NewEdgeEqOldEdge = 1;
  N_InnerEdges = 2;

  // initialize all dimension values
  MaxN_VpC = TRIBI12MAXN_VpC;
  MaxN_CpV = TRIBI12MAXN_CpV;
  MaxN_EpC = TRIBI12MAXN_EpC;
  MaxN_CpE = TRIBI12MAXN_CpE;
  MaxN_EpV = TRIBI12MAXN_EpV;
  MaxN_nVpoE = TRIBI12MAXN_nVpoE;
  MaxN_nEpoE = TRIBI12MAXN_nEpoE;

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

  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
}

// Methods
