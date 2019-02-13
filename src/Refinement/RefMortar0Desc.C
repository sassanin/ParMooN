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
// @(#)RefMortar0Desc.C        1.2 09/13/99
//
// Class:       TRefMortar0Desc
// Purpose:     refinement descriptor for a mortar cell (base edge 0)
//
// Author:      Volker Behns  30.12.97
//
// =======================================================================

#ifndef __MORTAR__
#define __MORTAR__
#endif

#include <RefMortar0Desc.h>

static Refinements DatEdgeType[] = { MortarLine, NoRef, MortarLine, NoRef};

// Constructor
TRefMortar0Desc::TRefMortar0Desc(TShapeDesc *shape, int Mortar_Ni, int N)
                  : TRefDesc(shape)
{
  int i;

  Shapes *DatChildType;

  int *DatChildVertex;
  int *DatVertexChild;
  int *DatVertexChildIndex;
  int *DatVertexChildLen;

  int *DatChildEdge;
  int *DatEdgeChild;
  int *DatEdgeChildIndex;
  int *DatEdgeChildLen;

  int *DatEdgeVertex;
  int *DatVertexEdge;
  int *DatVertexEdgeIndex;
  int *DatVertexEdgeLen;

  int *DatNewVertexEqOldVertex;
  int *DatNewVertexEqOldVertexIndex;

  int *DatNewEdgeEqOldEdge;
  int *DatNewEdgeEqOldEdgeIndex;

  int *DatInteriorEdgeOfCell;
  int *DatInteriorVertexOfEdge;
  int *DatInteriorVertexOfEdgeLen;

  int *DatOldEdgeNewVertex;

  int *DatOldEdgeNewEdge;
  int *DatOldEdgeNewLocEdge;
  int *DatNewEdgeOldEdge;

  Type = Mortar;

  // set all numbers
  N_Edges = 3*N + 1;
  N_Vertices = 2*N + 2;
  N_Children = N;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 2;
  N_InnerVertices = 0;    
  N_InnerEdges = N - 1;

  // initialize all dimension values
  MaxN_VpC = MORTARRMAXN_VpC;
  MaxN_CpV = MORTARRMAXN_CpV;
  MaxN_EpC = MORTARRMAXN_EpC;
  MaxN_CpE = MORTARRMAXN_CpE;
  MaxN_EpV = MORTARRMAXN_EpV;
  MaxN_iVpE = N - 1;
  MaxN_nVpoE = N + 1;
  MaxN_nEpoE = N;

  // generate data fields
  DatEdgeType[0] = (Refinements) (DatEdgeType[0] + Mortar_Ni);
  DatEdgeType[2] = (Refinements) (DatEdgeType[2] + Mortar_Ni);

  DatChildType = new Shapes[N_Children];
  for (i=0;i<N_Children;i++)
    DatChildType[i] = Quadrangle;

  DatChildVertex = new int[N_Children * MORTARRMAXN_VpC];
  for (i=0;i<N_Children;i++)
  {
    DatChildVertex[i*MORTARRMAXN_VpC    ] = 2*i;
    DatChildVertex[i*MORTARRMAXN_VpC + 1] = 2*i + 2;
    DatChildVertex[i*MORTARRMAXN_VpC + 2] = 2*i + 3;
    DatChildVertex[i*MORTARRMAXN_VpC + 3] = 2*i + 1;
  }

  DatVertexChild = new int[N_Vertices * MORTARRMAXN_CpV];
  DatVertexChildIndex = new int[N_Vertices * MORTARRMAXN_CpV];
  DatVertexChildLen = new int[N_Vertices];

  DatVertexChild[0] = 0;
  DatVertexChildIndex[0] = 0;
  DatVertexChildLen[0] = 1;
  DatVertexChild[MORTARRMAXN_CpV] = 0;
  DatVertexChildIndex[MORTARRMAXN_CpV] = 3;
  DatVertexChildLen[1] = 1;

  for (i=1;i<N_Children;i++)
  {
    DatVertexChild[2*i*MORTARRMAXN_CpV] = i - 1;
    DatVertexChild[2*i*MORTARRMAXN_CpV + 1] = i;
    DatVertexChild[(2*i+1)*MORTARRMAXN_CpV] = i - 1;
    DatVertexChild[(2*i+1)*MORTARRMAXN_CpV + 1] = i;
    DatVertexChildIndex[2*i*MORTARRMAXN_CpV] = 1;
    DatVertexChildIndex[2*i*MORTARRMAXN_CpV + 1] = 0;
    DatVertexChildIndex[(2*i+1)*MORTARRMAXN_CpV] = 2;
    DatVertexChildIndex[(2*i+1)*MORTARRMAXN_CpV + 1] = 3;
    DatVertexChildLen[2*i] = 2;
    DatVertexChildLen[2*i+1] = 2;
  }

  DatVertexChild[(N_Vertices-2)*MORTARRMAXN_CpV] = N_Children-1;
  DatVertexChildIndex[(N_Vertices-2)*MORTARRMAXN_CpV] = 1;
  DatVertexChildLen[N_Vertices-2] = 1;
  DatVertexChild[(N_Vertices-1)*MORTARRMAXN_CpV] = N_Children-1;
  DatVertexChildIndex[(N_Vertices-1)*MORTARRMAXN_CpV] = 2;
  DatVertexChildLen[N_Vertices-1] = 1;

  DatChildEdge = new int[N_Children * MORTARRMAXN_EpC];

  for (i=0;i<N_Children;i++)
  {
    DatChildEdge[i*MORTARRMAXN_EpC    ] = 3*i + 1;
    DatChildEdge[i*MORTARRMAXN_EpC + 1] = 3*i + 2;
    DatChildEdge[i*MORTARRMAXN_EpC + 2] = 3*i + 3;
    DatChildEdge[i*MORTARRMAXN_EpC + 3] = 3*i - 1;
  }
  DatChildEdge[3] = 0;

  DatEdgeChild = new int[N_Edges * MORTARRMAXN_CpE];
  DatEdgeChildIndex = new int[N_Edges * MORTARRMAXN_CpE];
  DatEdgeChildLen = new int[N_Edges];

  DatEdgeChild[0] = 0;
  DatEdgeChildIndex[0] = 3;
  DatEdgeChildLen[0] = 1;

  for (i=0;i<N_Children;i++)
  {
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+1)] = i;
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+2)] = i;
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+2) + 1] = i + 1;
    DatEdgeChild[MORTARRMAXN_CpE*(3*i+3)] = i;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+1)] = 0;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+2)] = 1;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+2) + 1] = 3;
    DatEdgeChildIndex[MORTARRMAXN_CpE*(3*i+3)] = 2;
    DatEdgeChildLen[3*i+1] = 1;
    DatEdgeChildLen[3*i+2] = 2;
    DatEdgeChildLen[3*i+3] = 1;
  }
  DatEdgeChildLen[N_Edges-2] = 1;

  DatEdgeVertex = new int[N_Edges * 2];

  DatEdgeVertex[0] = 0;
  DatEdgeVertex[1] = 1;

  for (i=0;i<N_Children;i++)
  {
    DatEdgeVertex[6*i+2] = 2*i;
    DatEdgeVertex[6*i+3] = 2*i + 2;
    DatEdgeVertex[6*i+4] = 2*i + 2;
    DatEdgeVertex[6*i+5] = 2*i + 3;
    DatEdgeVertex[6*i+6] = 2*i + 3;
    DatEdgeVertex[6*i+7] = 2*i + 1;
  }

  DatVertexEdge = new int[N_Vertices * MORTARRMAXN_EpV];
  DatVertexEdgeIndex = new int[N_Vertices * MORTARRMAXN_EpV];
  DatVertexEdgeLen = new int[N_Vertices];

  DatVertexEdge[0] = 0;
  DatVertexEdge[1] = 1;
  DatVertexEdgeIndex[0] = 3;
  DatVertexEdgeIndex[1] = 0;
  DatVertexEdgeLen[0] = 2;
  DatVertexEdge[MORTARRMAXN_EpV] = 3;
  DatVertexEdge[MORTARRMAXN_EpV + 1] = 0;
  DatVertexEdgeIndex[MORTARRMAXN_EpV] = 2;
  DatVertexEdgeIndex[MORTARRMAXN_EpV + 1] = 3;
  DatVertexEdgeLen[1] = 2;

  for (i=0;i<N_Children;i++)
  {
    DatVertexEdge[(2*i+2)*MORTARRMAXN_EpV    ] = 3*i + 1;
    DatVertexEdge[(2*i+2)*MORTARRMAXN_EpV + 1] = 3*i + 2;
    DatVertexEdge[(2*i+2)*MORTARRMAXN_EpV + 2] = 3*i + 4;
    DatVertexEdge[(2*i+3)*MORTARRMAXN_EpV    ] = 3*i + 2;
    DatVertexEdge[(2*i+3)*MORTARRMAXN_EpV + 1] = 3*i + 3;
    DatVertexEdge[(2*i+3)*MORTARRMAXN_EpV + 2] = 3*i + 5;
    DatVertexEdgeIndex[(2*i+2)*MORTARRMAXN_EpV    ] = 1;
    DatVertexEdgeIndex[(2*i+2)*MORTARRMAXN_EpV + 1] = 0;
    DatVertexEdgeIndex[(2*i+2)*MORTARRMAXN_EpV + 2] = 0;
    DatVertexEdgeIndex[(2*i+3)*MORTARRMAXN_EpV    ] = 1;
    DatVertexEdgeIndex[(2*i+3)*MORTARRMAXN_EpV + 1] = 0;
    DatVertexEdgeIndex[(2*i+3)*MORTARRMAXN_EpV + 2] = 1;
    DatVertexEdgeLen[2*i+2] = 3;
    DatVertexEdgeLen[2*i+3] = 3;
  }
  DatVertexEdgeLen[N_Vertices-2] = 2;
  DatVertexEdgeLen[N_Vertices-1] = 2;

  DatNewVertexEqOldVertex = new int[4];
  DatNewVertexEqOldVertexIndex = new int[4];

  DatNewVertexEqOldVertex[0] = 0;
  DatNewVertexEqOldVertex[1] = N_Vertices - 2;
  DatNewVertexEqOldVertex[2] = N_Vertices - 1;
  DatNewVertexEqOldVertex[3] = 1;
  DatNewVertexEqOldVertexIndex[0] = 0;
  DatNewVertexEqOldVertexIndex[1] = 1;
  DatNewVertexEqOldVertexIndex[2] = 2;
  DatNewVertexEqOldVertexIndex[3] = 3;

  DatNewEdgeEqOldEdge = new int[2];
  DatNewEdgeEqOldEdgeIndex = new int[2];

  DatNewEdgeEqOldEdge[0] = N_Edges - 2;
  DatNewEdgeEqOldEdge[1] = 0;
  DatNewEdgeEqOldEdgeIndex[0] = 1;
  DatNewEdgeEqOldEdgeIndex[1] = 3;

  DatInteriorEdgeOfCell = new int[N_InnerEdges];

  for (i=0;i<N_InnerEdges;i++)
    DatInteriorEdgeOfCell[i] = 3*i + 2;

  DatInteriorVertexOfEdge = new int[N_OrigEdges*MaxN_iVpE];
  DatInteriorVertexOfEdgeLen = new int[N_OrigEdges];

  for (i=0;i<MaxN_iVpE;i++)
  {
    DatInteriorVertexOfEdge[i] = 2*i + 2;
    DatInteriorVertexOfEdge[2*MaxN_iVpE + i] = N_Vertices - 2*i - 3;
  }

  DatInteriorVertexOfEdgeLen[0] = MaxN_iVpE;
  DatInteriorVertexOfEdgeLen[1] = 0;
  DatInteriorVertexOfEdgeLen[2] = MaxN_iVpE;
  DatInteriorVertexOfEdgeLen[3] = 0;

  DatOldEdgeNewVertex = new int[N_OrigEdges*MaxN_nVpoE];

  DatOldEdgeNewVertex[MaxN_nVpoE] = N_Vertices - 2;
  DatOldEdgeNewVertex[MaxN_nVpoE + 1] = N_Vertices - 1;

  for (i=0;i<N_Children+1;i++)
  {
    DatOldEdgeNewVertex[i] = 2*i;
    DatOldEdgeNewVertex[2*MaxN_nVpoE + i] = N_Vertices - 2*i - 1;
  }

  DatOldEdgeNewVertex[3*MaxN_nVpoE    ] = 1;
  DatOldEdgeNewVertex[3*MaxN_nVpoE + 1] = 0;
  
  DatOldEdgeNewEdge = new int[N_OrigEdges*MaxN_nEpoE];

  DatOldEdgeNewEdge[MaxN_nEpoE] = N_Edges - 2;
  DatOldEdgeNewEdge[3*MaxN_nEpoE] = 0;

  for (i=0;i<N_Children;i++)
  {
    DatOldEdgeNewEdge[i] = 3*i + 1;
    DatOldEdgeNewEdge[2*MaxN_nEpoE + i] = N_Edges - 3*i - 1;
  }

  DatOldEdgeNewLocEdge = new int[N_Children * 4];

  for (i=0;i<N_Children;i++)
  {
    DatOldEdgeNewLocEdge[4*i    ] = 0;
    DatOldEdgeNewLocEdge[4*i + 1] = -1;
    DatOldEdgeNewLocEdge[4*i + 2] = 2;
    DatOldEdgeNewLocEdge[4*i + 3] = -1;
  }
  DatOldEdgeNewLocEdge[3] = 3;
  DatOldEdgeNewLocEdge[4*N_Children - 3] = 1;
  
  DatNewEdgeOldEdge = new int[N_Edges];

  for (i=0;i<N_Children;i++)
  {
    DatNewEdgeOldEdge[3*i + 1] = 0;
    DatNewEdgeOldEdge[3*i + 2] = -1;
    DatNewEdgeOldEdge[3*i + 3] = 2;
  }
  DatNewEdgeOldEdge[0] = 3;
  DatNewEdgeOldEdge[N_Edges - 2] = 1;

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
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
}

// Methods
