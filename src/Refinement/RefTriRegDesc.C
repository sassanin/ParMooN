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
// @(#)RefTriRegDesc.C        1.2 09/13/99
//
// Class:       TRefTriRegDesc
// Purpose:     refinement descriptor for regular refinement of a triangle
//
// Author:      Volker Behns  30.07.97
//
// =======================================================================

#include <RefTriRegDesc.h>

static const Shapes DatChildType[] = { Triangle, Triangle, Triangle, Triangle};

static const Refinements DatEdgeType[] = { LineReg, LineReg, LineReg};

static const int DatChildVertex[][TRIRRMAXN_VpC] =
                 { {0, 3, 5},  {1, 4, 3},  {2, 5, 4},  {3, 4, 5}};
static const int DatVertexChild[][TRIRRMAXN_CpV] =
                 { {0}, {1}, {2}, {0, 1, 3},  {1, 2, 3},  {0, 2, 3}};
static const int DatVertexChildIndex[][TRIRRMAXN_CpV] =
                 { {0}, {0}, {0}, {1, 2, 0},  {1, 2, 1},  {2, 1, 2}};
static const int DatVertexChildLen[] = { 1,  1,  1,  3,  3,  3};

static const int DatChildEdge[][TRIRRMAXN_EpC] =
                 { {0, 8, 5},  {2, 6, 1},  {4, 7, 3},  {6, 7, 8}};
static const int DatEdgeChild[][TRIRRMAXN_CpE] =
                 { {0}, {1}, {1}, {2}, {2}, {0}, {1, 3}, {2, 3}, {0, 3}};
static const int DatEdgeChildIndex[][TRIRRMAXN_CpE] = 
                 { {0}, {2}, {0}, {2}, {0}, {2}, {1, 0}, {1, 1}, {1, 2}};
static const int DatEdgeChildLen[] = { 1,  1,  1,  1,  1,  1,  2,  2,  2};

static const int DatEdgeVertex[][2] =
                 { {0, 3},  {3, 1},  {1, 4},  {4, 2},  {2, 5},
                   {5, 0},  {3, 4},  {4, 5},  {5, 3}};
static const int DatVertexEdge[][TRIRRMAXN_EpV] = 
                 { {0, 5},  {1, 2},  {3, 4},  {0, 1, 6, 8},
                   {2, 3, 6, 7},  {4, 5, 7, 8}};
static const int DatVertexEdgeIndex[][TRIRRMAXN_EpV] =
                 { {0, 1},  {1, 0},  {1, 0},  {1, 0, 0, 1},
                   {1, 0, 1, 0},  {1, 0, 1, 0}};
static const int DatVertexEdgeLen[] = { 2,  2,  2,  4,  4,  4};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2};

static const int DatInteriorEdgeOfCell[] = { 6, 7, 8};

static const int DatInteriorVertexOfEdge[][TRIRRMAXN_iVpE] =
                 { {3},  {4},  {5}};
static const int DatInteriorVertexOfEdgeLen[] = { 1,  1,  1};

static const int DatOldEdgeNewVertex[][TRIRRMAXN_nVpoE] =
                 { {0, 3, 1}, {1, 4, 2},  {2, 5, 0}};

static const int DatOldEdgeNewEdge[][TRIRRMAXN_nEpoE] =
                 { {0, 1},  {2, 3},  {4, 5}};

static const int DatOldEdgeNewLocEdge[][TRIRRN_E] =
                 { {0, -1, 2}, {2, 0, -1}, {-1, 2, 0}, {-1, -1, -1} };

static const int DatNewEdgeOldEdge[] =
                 { 0,  0,  1,  1,  2,  2, -1, -1, -1};

// Constructor
TRefTriRegDesc::TRefTriRegDesc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TriReg;

  // set all numbers
  N_Edges = 9;
  N_Vertices = 6;
  N_Children = 4;
  N_NewVertEqOldVert = 3;
  N_InnerEdges = 3;

  // initialize all dimension values
  MaxN_VpC = TRIRRMAXN_VpC;
  MaxN_CpV = TRIRRMAXN_CpV;
  MaxN_EpC = TRIRRMAXN_EpC;
  MaxN_CpE = TRIRRMAXN_CpE;
  MaxN_EpV = TRIRRMAXN_EpV;
  MaxN_iVpE = TRIRRMAXN_iVpE;
  MaxN_nVpoE = TRIRRMAXN_nVpoE;
  MaxN_nEpoE = TRIRRMAXN_nEpoE;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;

  ChildVertex = (const int *) &(DatChildVertex[0][0]);
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

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge= (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
}


// Methods
