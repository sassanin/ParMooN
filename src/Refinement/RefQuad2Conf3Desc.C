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
// @(#)RefQuad2Conf3Desc.C        1.2 09/13/99
//
// Class:       TRefQuadConf3Desc
// Purpose:     refinement descriptor for conforming closure of a
//              quadrangle with 2 hanging node
//
// Author:      Matthias Ebeling  20.08.99
//
// =======================================================================

#include <RefQuad2Conf3Desc.h>

static const Shapes DatChildType[] = { Quadrangle, Quadrangle,
                                       Quadrangle };

static const Refinements DatEdgeType[] = { NoRef, NoRef, LineReg, LineReg};

static const int DatChildVertex[][QUAD2ConfMAXN_VpC] =
                 { {0, 4, 6, 5},  {4, 1, 2, 6},  {5, 6, 2, 3}};
static const int DatVertexChild[][QUAD2ConfMAXN_CpV] =
                 { {0},  {1},  {1,2},  {2},  {0, 1},  {0, 2},  {0, 1, 2}};
                 
static const int DatVertexChildIndex[][QUAD2ConfMAXN_CpV] =
               { {0},  {1},  {2, 2},  {3},  {1, 0},  {3, 0},  {2, 3, 1}};

static const int DatVertexChildLen[] =
               { 1,  1,  2,  1,  2,  2,  3};

static const int DatChildEdge[][QUAD2ConfMAXN_EpC] =
               { {0, 6, 8, 5},  {1, 2, 7, 6},  {8, 7, 3, 4}};
static const int DatEdgeChild[][QUAD2ConfMAXN_CpE] =
               { {0},  {1},  {1},  {2},  {2},  {0},  {0, 1},  {1, 2},
                 {0, 2}};
static const int DatEdgeChildIndex[][QUAD2ConfMAXN_CpE] = 
               { {0},  {0},  {1},  {2},  {3},  {3},  {1, 3},  {2, 1},
                 {2, 0}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  1,  1,  2,  2,  2};

static const int DatEdgeVertex[][2] =
               { {0, 4},  {4, 1},  {1, 2},  {2, 3},  {3, 5},  {5, 0},
                 {4, 6},  {2, 6},  {5, 6}};
static const int DatVertexEdge[][QUAD2ConfMAXN_EpV] = 
               { {0, 5},  {1, 2},  {2, 3},  {3, 4},  {0, 1, 6},
                 {5, 8, 4},  {6, 7, 8}};
static const int DatVertexEdgeIndex[][QUAD2ConfMAXN_EpV] =
               { {0, 1},  {1, 0},  {1, 0},  {1,0},  {1, 0, 0},
                 {0, 0, 1},  {1, 1, 1}};
static const int DatVertexEdgeLen[] =
               { 2,  2,  2,  2,  3,  3,  3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 3, 0, 1, 2};

static const int DatInteriorVertexOfCell[] = { 6};
static const double DatPositionOfIntVert[][QUADN_V] =
                  { {0.25, 0.25, 0.25, 0.25}};
                                 
static const int DatInteriorEdgeOfCell[] = { 6, 7, 8};
static const int DatInteriorVertexOfEdge[][QUAD2ConfMAXN_iVpE] =
               { {-1}, {-1}, {5}, {4}};
static const int DatInteriorVertexOfEdgeLen[] = { 0, 0, 1, 1};

static const int DatOldEdgeNewVertex[][QUAD2ConfMAXN_nVpoE] =
               { {1, 2},  {2, 3},  {3, 5, 0}, {0, 4, 1}};

static const int DatOldEdgeNewEdge[][QUAD2ConfMAXN_nEpoE] =
               { {2},  {3},  {4, 5},  {0, 1}};

static const int DatOldEdgeNewLocEdge[][QUADConfN_E] =
               { {-1, -1, 3, 0}, {1, -1, -1, 0},
                 {-1, 2, 3, -1} };

static const int DatNewEdgeOldEdge[] =
               { 3,  3,  0,  1 ,  2,  2,  -1,  -1,  -1};

static const int DatNewEdgeEqOldEdge[] = { 2, 3};
static const int DatNewEdgeEqOldEdgeIndex[] = { 0, 1};

// Constructor
TRefQuad2Conf3Desc::TRefQuad2Conf3Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = Quad2Conf3;

  // set all numbers
  N_Edges = 9;
  N_Vertices = 7;
  N_Children = 3;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 2;
  N_InnerVertices = 1;    
  N_InnerEdges = 3;

  // initialize all dimension values
  MaxN_VpC = QUAD2ConfMAXN_VpC;
  MaxN_CpV = QUAD2ConfMAXN_CpV;
  MaxN_EpC = QUAD2ConfMAXN_EpC;
  MaxN_CpE = QUAD2ConfMAXN_CpE;
  MaxN_EpV = QUAD2ConfMAXN_EpV;
  MaxN_iVpE = QUAD2ConfMAXN_iVpE;
  MaxN_nVpoE = QUAD2ConfMAXN_nVpoE;
  MaxN_nEpoE = QUAD2ConfMAXN_nEpoE;

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

  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;
}

// Methods
