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
// @(#)RefLineDesc.C        1.2 09/13/99
//
// Class:       TRefLineDesc
// Purpose:     refinement descriptor for line
//
// Author:      Volker Behns  07.08.97
//
// =======================================================================

#include <RefLineDesc.h>

static const Shapes DatChildType[] = { S_Line, S_Line};

static const int DatChildVertex[][LINEMAXN_VpC] = { {0, 2},  {2, 1}};
static const int DatVertexChild[][LINEMAXN_CpV] = { {0},  {1},  {0, 1}};
static const int DatVertexChildIndex[][LINEMAXN_CpV] = { {0},  {1},  {1, 0}};

static const int DatNewVertexEqOldVertex[] = { 0, 1};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1};

static const int DatInteriorVertexOfCell[] = { 2};
static const double DatPositionOfIntVert[][LINEMAXN_V] = { {0.5, 0.5}};

// Constructor
TRefLineDesc::TRefLineDesc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = LineReg;

  N_Edges = 2;
  N_Vertices = 3;
  N_Children = 2;
  N_NewVertEqOldVert = 2;


  MaxN_VpC = LINEMAXN_VpC;
  MaxN_CpV = LINEMAXN_CpV;

  ChildType = (const Shapes *) DatChildType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  InteriorVertexOfCell = (const int *) DatInteriorVertexOfCell;
  PositionOfIntVert = (const double *) DatPositionOfIntVert;
}

// Methods
