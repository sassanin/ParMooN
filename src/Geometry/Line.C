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
// @(#)Line.C        1.1 10/30/98
//
// Class:       TLine
// Purpose:     shape descriptor of a line
//
// Author:      Volker Behns  07.08.97
//
// =======================================================================

#include <Line.h>

static const int DatEdgeVertex[][2] = { {0, 1}};
static const int DatVertexEdge[][LINEMAXN_EpV] = { {0},  {0}};

// Constructor
TLine::TLine()
{
  MaxN_EpV = LINEMAXN_EpV;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  Type = S_Line;
  N_Vertices = 2;
  N_Edges = 2;
  N_Joints = 2;
}

// Methods
double TLine::GetMeasure(TVertex **Verts)
{
  double x1,x2,y1,y2;

  x1 = Verts[0]->GetX();
  y1 = Verts[0]->GetY();
  x2 = Verts[1]->GetX();
  y2 = Verts[1]->GetY();

  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
