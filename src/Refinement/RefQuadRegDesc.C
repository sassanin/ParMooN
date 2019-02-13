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
// @(#)RefQuadRegDesc.C        1.2 09/13/99
//
// Class:       TRefQuadRegDesc
// Purpose:     refinement descriptor for regular refinement of a quadrangle
//
// Author:      Volker Behns  30.07.97
//
// =======================================================================

#include <RefQuadRegDesc.h>

static const Shapes DatChildType[] = { Quadrangle, Quadrangle,
                                       Quadrangle, Quadrangle};

static const Refinements DatEdgeType[] = { LineReg, LineReg, LineReg, LineReg};

static const int DatChildVertex[][QUADRRMAXN_VpC] =
                 { {0, 4, 8, 7},  {1, 5, 8, 4},  {2, 6, 8, 5},  {3, 7, 8, 6}};
static const int DatVertexChild[][QUADRRMAXN_CpV] =
                 { {0},  {1},  {2},  {3},  {0, 1},  {1, 2},  {2, 3},
                   {0, 3},  {0, 1, 2, 3}};
static const int DatVertexChildIndex[][QUADRRMAXN_CpV] =
               { {0},  {0},  {0},  {0},  {1, 3},  {1, 3},  {1, 3},
                 {3, 1},  {2, 2, 2, 2}};
static const int DatVertexChildLen[] =
               { 1,  1,  1,  1,  2,  2,  2,  2,  4};

static const int DatChildEdge[][QUADRRMAXN_EpC] =
               { {0, 8,11, 7},  {2, 9, 8, 1},  {4,10, 9, 3},  {6,11,10, 5}};
static const int DatEdgeChild[][QUADRRMAXN_CpE] =
               { {0},  {1},  {1},  {2},  {2},  {3},  {3},  {0},
                 {0, 1},  {1, 2},  {2, 3},  {0, 3}};
static const int DatEdgeChildIndex[][QUADRRMAXN_CpE] = 
               { {0},  {3},  {0},  {3},  {0},  {3},  {0},  {3},
                 {1, 2},  {1, 2},  {1, 2},  {2, 1}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2};

static const int DatEdgeVertex[][2] =
               { {0, 4},  {4, 1},  {1, 5},  {5, 2},  {2, 6},  {6, 3},
                 {3, 7},  {7, 0},  {4, 8},  {5, 8},  {6, 8},  {7, 8}};
static const int DatVertexEdge[][QUADRRMAXN_EpV] = 
               { {0, 7},  {1, 2},  {3, 4},  {5, 6},  {0, 1, 8},
                 {2, 3, 9},  {4, 5,10},  {6, 7,11},  {8, 9,10,11}};
static const int DatVertexEdgeIndex[][QUADRRMAXN_EpV] =
               { {0, 1},  {1, 0},  {1, 0},  {1,0},  {1, 0, 0},
                 {1, 0, 0},  {1, 0, 0},  {1, 0, 0},  {1, 1, 1, 1}};
static const int DatVertexEdgeLen[] =
               { 2,  2,  2,  2,  3,  3,  3,  3,  4};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2, 3};

static const int DatInteriorVertexOfCell[] = { 8};
static const double DatPositionOfIntVert[][QUADN_V] =
                  { {0.25, 0.25, 0.25, 0.25}};
                                 
static const int DatInteriorEdgeOfCell[] = { 8, 9, 10, 11};
static const int DatInteriorVertexOfEdge[][QUADRRMAXN_iVpE] =
               { {4}, {5}, {6}, {7}};
static const int DatInteriorVertexOfEdgeLen[] = { 1,  1,  1,  1};

static const int DatOldEdgeNewVertex[][QUADRRMAXN_nVpoE] =
               { {0, 4, 1},  {1, 5, 2},  {2, 6, 3},  {3, 7, 0}};

static const int DatOldEdgeNewEdge[][QUADRRMAXN_nEpoE] =
               { {0, 1},  {2, 3},  {4, 5},  {6, 7}};

static const int DatOldEdgeNewLocEdge[][QUADRRN_E] =
               { {0, -1, -1, 3}, {3, 0, -1, -1},
                 {-1, 3, 0, -1}, {-1, -1, 3, 0} };

static const int DatNewEdgeOldEdge[] =
               { 0,  0,  1,  1,  2,  2,  3,  3,  -1,  -1,  -1,  -1};

// Constructor
TRefQuadRegDesc::TRefQuadRegDesc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = QuadReg;

  // set all numbers
  N_Edges = 12;
  N_Vertices = 9;
  N_Children = 4;
  N_NewVertEqOldVert = 4;
  N_InnerVertices = 1;    
  N_InnerEdges = 4;

  // initialize all dimension values
  MaxN_VpC = QUADRRMAXN_VpC;
  MaxN_CpV = QUADRRMAXN_CpV;
  MaxN_EpC = QUADRRMAXN_EpC;
  MaxN_CpE = QUADRRMAXN_CpE;
  MaxN_EpV = QUADRRMAXN_EpV;
  MaxN_iVpE = QUADRRMAXN_iVpE;
  MaxN_nVpoE = QUADRRMAXN_nVpoE;
  MaxN_nEpoE = QUADRRMAXN_nEpoE;

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

  InteriorVertexOfCell = (const int *) DatInteriorVertexOfCell;
  PositionOfIntVert = (const double *) DatPositionOfIntVert;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
}

// Methods
