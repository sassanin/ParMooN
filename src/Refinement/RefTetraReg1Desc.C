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
// @(#)RefTetraReg1Desc.C        1.3 11/15/99
//
// Class:       TRefTetraRegDesc
// Purpose:     refinement descriptor for regular refinement of a
//              tetrahedron, the octahedron will be divided by a edge from 
//              vertex 5 to 7
//
// Author:      Volker Behns  31.07.97
//
// =======================================================================

#ifndef __3D__
#define __3D__
#endif

#include <RefTetraReg1Desc.h>

static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron,
                                      Tetrahedron, Tetrahedron, Tetrahedron,
                                      Tetrahedron, Tetrahedron };

static const Refinements DatEdgeType[] = { LineReg, LineReg, LineReg, LineReg, 
                                           LineReg, LineReg};

static const Refinements DatFaceType[] = { TriReg, TriReg, TriReg, TriReg};

static const int DatChildVertex[][REFTETRAREGMAXN_VpC] =
           { { 0, 4, 6, 7}, { 1, 5, 4, 8}, { 2, 6, 5, 9}, { 3, 8, 7, 9}, { 4, 5, 6, 7}, { 4, 8, 5, 7}, { 5, 8, 9, 7}, { 6, 5, 9, 7}};

static const int DatVertexChild[][REFTETRAREGMAXN_CpV] =
           { { 0}, { 1}, { 2}, { 3}, { 0, 1, 4, 5}, { 1, 2, 4, 5, 6, 7}, { 0, 2, 4, 7}, { 0, 3, 4, 5, 6, 7}, { 1, 3, 5, 6}, { 2, 3, 6, 7}};

static const int DatVertexChildIndex[][REFTETRAREGMAXN_CpV] =
           { { 0}, { 0}, { 0}, { 0}, { 1, 2, 0, 0}, { 1, 2, 1, 2, 0, 1}, { 2, 1, 2, 0}, { 3, 2, 3, 3, 3, 3}, { 3, 1, 1, 1}, { 3, 3, 2, 2}};

static const int DatVertexChildLen[] =
           { 1, 1, 1, 1, 4, 6, 4, 6, 4, 4};

static const int DatChildEdge[][REFTETRAREGMAXN_EpC] =
           { { 0, 12, 5, 6, 15, 22}, { 2, 13, 1, 8, 18, 16}, { 4, 14, 3, 10, 21, 19}, { 9, 17, 7, 11, 20, 23}, { 13, 14, 12, 15, 24, 22}, { 16, 18, 13, 15, 17, 24}, { 18, 20, 19, 24, 17, 23}, { 14, 19, 21, 22, 24, 23}};

static const int DatEdgeChild[][REFTETRAREGMAXN_CpE] =
           { { 0}, { 1}, { 1}, { 2}, { 2}, { 0}, { 0}, { 3}, { 1}, { 3}, { 2}, { 3}, { 0, 4}, { 1, 4, 5}, { 2, 4, 7}, { 0, 4, 5}, { 1, 5}, { 3, 5, 6}, { 1, 5, 6}, { 2, 6, 7}, { 3, 6}, { 2, 7}, { 0, 4, 7}, { 3, 6, 7}, { 4, 5, 6, 7}};

static const int DatEdgeChildIndex[][REFTETRAREGMAXN_CpE] =
           { { 0}, { 2}, { 0}, { 2}, { 0}, { 2}, { 3}, { 2}, { 3}, { 0}, { 3}, { 3}, { 1, 2}, { 1, 0, 2}, { 1, 1, 0}, { 4, 3, 3}, { 5, 0}, { 1, 4, 4}, { 4, 1, 0}, { 5, 2, 1}, { 4, 1}, { 4, 2}, { 5, 5, 3}, { 5, 5, 5}, { 4, 5, 3, 4}};

static const int DatEdgeChildLen[] =
           { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 2, 3, 3, 3, 2, 2, 3, 3, 4};

static const int DatChildFace[][REFTETRAREGMAXN_FpC] =
           { { 0, 1, 2, 3}, { 4, 5, 6, 7}, { 8, 9, 10, 11}, { 12, 13, 14, 15}, { 16, 17, 18, 2}, { 6, 19, 20, 17}, { 21, 20, 14, 22}, { 10, 18, 22, 23}};

static const int DatFaceChild[][REFTETRAREGMAXN_CpF] =
           { { 0}, { 0}, { 0, 4}, { 0}, { 1}, { 1}, { 1, 5}, { 1}, { 2}, { 2}, { 2, 7}, { 2}, { 3}, { 3}, { 3, 6}, { 3}, { 4}, { 4, 5}, { 4, 7}, { 5}, { 5, 6}, { 6}, { 6, 7}, { 7}};

static const int DatFaceChildIndex[][REFTETRAREGMAXN_CpF] =
           { { 0}, { 1}, { 2, 3}, { 3}, { 0}, { 1}, { 2, 0}, { 3}, { 0}, { 1}, { 2, 0}, { 3}, { 0}, { 1}, { 2, 2}, { 3}, { 0}, { 1, 3}, { 2, 1}, { 1}, { 2, 1}, { 0}, { 3, 2}, { 3}};

static const int DatFaceChildLen[] =
           { 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, 1};

static const int DatEdgeVertex[][2] =
           { { 0, 4}, { 4, 1}, { 1, 5}, { 5, 2}, { 2, 6}, { 6, 0}, { 0, 7}, { 7, 3}, { 1, 8}, { 8, 3}, { 2, 9}, { 9, 3}, { 6, 4}, { 4, 5}, { 5, 6}, { 7, 4}, { 4, 8}, { 8, 7}, { 8, 5}, { 5, 9}, { 9, 8}, { 9, 6}, { 6, 7}, { 7, 9}, { 5, 7}};

static const int DatVertexEdge[][REFTETRAREGMAXN_EpV] =
           { { 0, 5, 6}, { 1, 2, 8}, { 3, 4, 10}, { 7, 9, 11}, { 0, 1, 12, 13, 15, 16}, { 2, 3, 13, 14, 18, 19, 24}, { 4, 5, 12, 14, 21, 22}, { 6, 7, 15, 17, 22, 23, 24}, { 8, 9, 16, 17, 18, 20}, { 10, 11, 19, 20, 21, 23}};

static const int DatVertexEdgeIndex[][REFTETRAREGMAXN_EpV] =
           { { 0, 1, 0}, { 1, 0, 0}, { 1, 0, 0}, { 1, 1, 1}, { 1, 0, 1, 0, 1, 0}, { 1, 0, 1, 0, 1, 0, 0}, { 1, 0, 0, 1, 1, 0}, { 1, 0, 0, 1, 1, 0, 1}, { 1, 0, 1, 0, 0, 1}, { 1, 0, 1, 0, 0, 1}};

static const int DatVertexEdgeLen[] =
           { 3, 3, 3, 3, 6, 7, 6, 7, 6, 6};

static const int DatFaceVertex[][REFTETRAREGMAXN_VpF] =
           { { 0, 4, 6}, { 0, 7, 4}, { 6, 4, 7}, { 0, 6, 7}, { 1, 5, 4}, { 1, 8, 5}, { 4, 5, 8}, { 1, 4, 8}, { 2, 6, 5}, { 2, 9, 6}, { 5, 6, 9}, { 2, 5, 9}, { 3, 8, 7}, { 3, 9, 8}, { 7, 8, 9}, { 3, 7, 9}, { 4, 5, 6}, { 4, 7, 5}, { 6, 5, 7}, { 4, 7, 8}, { 5, 8, 7}, { 5, 8, 9}, { 5, 9, 7}, { 6, 9, 7} };

static const int DatVertexFace[][REFTETRAREGMAXN_FpV] =
           { { 0, 1, 3}, { 4, 5, 7}, { 8, 9, 11}, { 12, 13, 15}, { 0, 1, 2, 4, 6, 7, 16, 17, 19}, { 4, 5, 6, 8, 10, 11, 16, 17, 18, 20, 21, 22}, { 0, 2, 3, 8, 9, 10, 16, 18, 23}, { 1, 2, 3, 12, 14, 15, 17, 18, 19, 20, 22, 23}, { 5, 6, 7, 12, 13, 14, 19, 20, 21}, { 9, 10, 11, 13, 14, 15, 21, 22, 23}};

static const int DatVertexFaceIndex[][REFTETRAREGMAXN_FpV] =
           { { 0, 0, 2}, { 0, 0, 2}, { 0, 0, 2}, { 0, 0, 2}, { 1, 2, 0, 2, 2, 0, 0, 0, 0}, { 1, 2, 0, 2, 2, 0, 1, 2, 0, 2, 0, 2}, { 2, 2, 0, 1, 2, 0, 2, 2, 2}, { 1, 1, 1, 2, 2, 0, 1, 1, 1, 1, 1, 1}, { 1, 1, 1, 1, 2, 0, 2, 0, 1}, { 1, 1, 1, 1, 1, 1, 2, 0, 0}};

static const int DatVertexFaceLen[] =
           { 3, 3, 3, 3, 9, 12, 9, 12, 9, 9};

static const int DatFaceEdge[][REFTETRAREGMAXN_EpF] =
           { { 0, 12, 5}, { 6, 15, 0}, { 15, 22, 12}, { 22, 6, 5}, { 2, 13, 1}, { 8, 18, 2}, { 18, 16, 13}, { 16, 8, 1}, { 4, 14, 3}, { 10, 21, 4}, { 21, 19, 14}, { 19, 10, 3}, { 9, 17, 7}, { 11, 20, 9}, { 20, 23, 17}, { 23, 11, 7}, { 13, 14, 12}, { 15, 24, 13}, { 24, 22, 14}, { 15, 17, 16}, { 17, 24, 18}, { 18, 20, 19}, { 23, 24, 19}, { 23, 22, 21}};

static const int DatEdgeFace[][REFTETRAREGMAXN_FpE] =
           { { 0, 1}, { 4, 7}, { 4, 5}, { 8, 11}, { 8, 9}, { 0, 3}, { 1, 3}, { 12, 15}, { 5, 7}, { 12, 13}, { 9, 11}, { 13, 15}, { 0, 2, 16}, { 4, 6, 16, 17}, { 8, 10, 16, 18}, { 1, 2, 17, 19}, { 6, 7, 19}, { 12, 14, 19, 20}, { 5, 6, 20, 21}, { 10, 11, 21, 22}, { 13, 14, 21}, { 9, 10, 23}, { 2, 3, 18, 23}, { 14, 15, 22, 23}, { 17, 18, 20, 22}};

static const int DatEdgeFaceIndex[][REFTETRAREGMAXN_FpE] =
           { { 0, 2}, { 2, 2}, { 0, 2}, { 2, 2}, { 0, 2}, { 2, 2}, { 0, 1}, { 2, 2}, { 0, 1}, { 0, 2}, { 0, 1}, { 0, 1}, { 1, 2, 2}, { 1, 2, 0, 2}, { 1, 2, 1, 2}, { 1, 0, 0, 0}, { 1, 0, 2}, { 1, 2, 1, 0}, { 1, 0, 2, 0}, { 1, 0, 2, 2}, { 1, 0, 1}, { 1, 0, 2}, { 1, 0, 1, 1}, { 1, 0, 0, 0}, { 1, 0, 1, 1}};

static const int DatEdgeFaceLen[] =
           { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 3, 4, 4, 4, 3, 3, 4, 4, 4};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2, 3};

static const double DatOldFaceNewVertexPos[][REFTETRAREGMAXN_nVpoF]
                                            [REFTETRAREGMAXN_oVpoF] =
           { { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}},
             { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}},
             { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}},
             { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}}};

static const int DatInteriorEdgeOfCell[] =
           { 24};

static const int DatInteriorFaceOfCell[] =
           { 2, 6, 10, 14, 17, 18, 20, 22};

static const int DatInteriorVertexOfEdge[][REFTETRAREGMAXN_iVpE] =
           {{4}, {5}, {6}, {7}, {8}, {9}};

static const int DatInteriorVertexOfEdgeLen[]=
           { 1, 1, 1, 1, 1, 1};

static const int DatInteriorEdgeOfFace[][REFTETRAREGMAXN_iEpF] =
           { { 12, 13, 14}, { 15, 16, 17}, { 18, 19, 20}, { 21, 22, 23}};

static const int DatInteriorEdgeOfFaceLen[] =
           { 3, 3, 3, 3};

static const int DatOldEdgeNewVertex[][REFTETRAREGMAXN_nVpoE] =
           { { 0, 4, 1}, { 1, 5, 2}, { 2, 6, 0}, { 0, 7, 3}, { 1, 8, 3}, { 2, 9, 3}};

static const int DatOldEdgeNewVertexLen[] =
           { 3, 3, 3, 3, 3, 3};

static const int DatOldEdgeNewEdge[][REFTETRAREGMAXN_nEpoE] =
           { { 0, 1}, { 2, 3}, { 4, 5}, { 6, 7}, { 8, 9}, { 10, 11}};

static const int DatOldEdgeNewEdgeLen[] =
           { 2, 2, 2, 2, 2, 2};

static const int DatOldFaceNewVertex[][REFTETRAREGMAXN_nVpoF] =
           { { 0, 1, 2, 4, 5, 6}, { 0, 3, 1, 7, 8, 4}, 
             { 2, 1, 3, 5, 8, 9}, { 0, 2, 3, 6, 9, 7}};

static const int DatOldFaceNewVertexLen[] =
           { 6, 6, 6, 6};

static const int DatOldFaceNewEdge[][REFTETRAREGMAXN_nEpoF] =
           { { 0, 1, 2, 3, 4, 5, 13, 14, 12}, { 6, 7, 9, 8, 1, 0, 17, 16, 15}, 
             { 3, 2, 8, 9, 11, 10, 18, 20, 19}, { 5, 4, 10, 11, 7, 6, 21, 23, 22}};

static const int DatOldFaceNewEdgeLen[] =
           { 9, 9, 9, 9};

static const int DatOldFaceNewFace[][REFTETRAREGMAXN_nFpoF] =
           { { 0, 4, 8, 16}, { 1, 12, 7, 19}, { 11, 5, 13, 21}, { 3, 9, 15, 23}};

static const int DatOldFaceNewFaceLen[] =
           { 4, 4, 4, 4};

static const int DatNewEdgeOldEdge[] = 
           { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
           { 0, 1, -1, 3, 0, 2, -1, 1, 0, 3, -1, 2, 1, 2, -1, 3, 0, -1, -1, 1, -1, 2, -1, 3};

static const int DatOldFaceNewLocFace[][TETRAN_F] =
           { { 0, 1, -1, 3}, { 0, 3, 1, -1}, { 0, -1, 3, 1}, { -1, 0, 1, 3}, { 0, -1, -1, -1}, { -1, 1, -1, -1}, { -1, -1, 0, -1}, { -1, -1, -1, 3}, };

static const int DatChildTwistIndex[] =
           { 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, -1, 2, -1, 0, -1, 0};

// Constructor
TRefTetraReg1Desc::TRefTetraReg1Desc(TShapeDesc *shape) : TRefDesc(shape)
{

  Type = TetraReg1;

  //set all numbers
  N_Vertices = 10;
  N_Edges = 25;
  N_Faces = 24;
  N_Children = 8;
  N_NewVertEqOldVert = 4;
  N_InnerEdges = 1;
  N_InnerFaces = 8;

  // initialize all dimension values
  MaxN_VpC = REFTETRAREGMAXN_VpC;
  MaxN_CpV = REFTETRAREGMAXN_CpV;
  MaxN_EpC = REFTETRAREGMAXN_EpC;
  MaxN_CpE = REFTETRAREGMAXN_CpE;
  MaxN_EpV = REFTETRAREGMAXN_EpV;
  MaxN_EpF = REFTETRAREGMAXN_EpF;
  MaxN_FpE = REFTETRAREGMAXN_FpE;
  MaxN_VpF = REFTETRAREGMAXN_VpF;
  MaxN_FpV = REFTETRAREGMAXN_FpV;
  MaxN_FpC = REFTETRAREGMAXN_FpC;
  MaxN_CpF = REFTETRAREGMAXN_CpF;
  MaxN_iVpE = REFTETRAREGMAXN_iVpE;
  MaxN_iEpF = REFTETRAREGMAXN_iEpF;
  MaxN_nVpoE = REFTETRAREGMAXN_nVpoE;
  MaxN_nEpoE = REFTETRAREGMAXN_nEpoE;
  MaxN_nVpoF = REFTETRAREGMAXN_nVpoF;
  MaxN_oVpoF = REFTETRAREGMAXN_oVpoF;
  MaxN_nEpoF = REFTETRAREGMAXN_nEpoF;
  MaxN_nFpoF = REFTETRAREGMAXN_nFpoF;

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
