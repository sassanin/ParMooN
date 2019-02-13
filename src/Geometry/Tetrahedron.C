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
// @(#)Tetrahedron.C        1.4 11/15/99
//
// Class:       TTetrahedron
// Purpose:     shape descriptor of a tetrahedron
//
// Author:      Volker Behns  26.07.97
//
// =======================================================================

#ifndef __3D__
#define __3D__
#endif

#include <Tetrahedron.h>

static const int DatEdgeVertex[][2] = { {0, 1},  {1, 2},  {2, 0},
                                        {0, 3},  {1, 3},  {2, 3}};
static const int DatVertexEdge[][TETRAMAXN_EpV] =
                 { {2, 0, 3},  {0, 1, 4},  {1, 2, 5},  {3, 4, 5}};

static const int DatFaceVertex[][TETRAMAXN_VpF] =
                   { {0, 1, 2},  {0, 3, 1},  {2, 1, 3},  {0, 2, 3}};

static const int DatFaceVertexLen[] = { 3, 3, 3, 3};

static const int DatVertexFace[][TETRAMAXN_FpV] =
                 { {0, 1, 3},  {0, 2, 1},  {0, 3, 2},  {1, 2, 3}};

static const int DatFaceEdge[][TETRAMAXN_EpF] =
                 { {0, 1, 2},  {3, 4, 0},  {1, 4, 5},  {2, 5, 3}};

static const int DatFaceEdgeLen[] =
                 { 3, 3, 3, 3};

static const int DatEdgeFace[][TETRAMAXN_FpE] =
                 { {1, 0},  {2, 0},  {3, 0}, {3, 1},  {1, 2},  {2, 3}};

static const Shapes DatFaceType[] = { Triangle, Triangle, Triangle, Triangle};

// Constructor
TTetrahedron::TTetrahedron()
{
  MaxN_EpV = TETRAMAXN_EpV;
  MaxN_VpF = TETRAMAXN_VpF;
  MaxN_FpV = TETRAMAXN_FpV;
  MaxN_EpF = TETRAMAXN_EpF;
  MaxN_FpE = TETRAMAXN_FpE;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  FaceVertex = (const int *) DatFaceVertex;
  FaceVertexLen = (const int *) DatFaceVertexLen;

  VertexFace = (const int *) DatVertexFace;

  FaceEdge = (const int *) DatFaceEdge;
  FaceEdgeLen = (const int *) DatFaceEdgeLen;
  EdgeFace = (const int *) DatEdgeFace;

  FaceType = (const Shapes *) DatFaceType;

  Type = Tetrahedron;
  N_Vertices = 4;
  N_Edges = 6;
  N_Faces = 4;
  N_Joints = 4;
}

// Methods
double TTetrahedron::GetDiameter(TVertex **Verts)
{
  int i;
  double len, diam = 0;
  double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[1]->GetCoords(x1, y1, z1);
  Verts[2]->GetCoords(x2, y2, z2);
  Verts[3]->GetCoords(x3, y3, z3);

  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len > diam) diam = len;
  len = (x0-x2)*(x0-x2)+(y0-y2)*(y0-y2)+(z0-z2)*(z0-z2);
  if(len > diam) diam = len;
  len = (x0-x3)*(x0-x3)+(y0-y3)*(y0-y3)+(z0-z3)*(z0-z3);
  if(len > diam) diam = len;
  len = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
  if(len > diam) diam = len;
  len = (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3);
  if(len > diam) diam = len;
  len = (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+(z2-z3)*(z2-z3);
  if(len > diam) diam = len;
  
  return sqrt(diam);
}

double TTetrahedron::GetShortestEdge(TVertex **Verts)
{
  int i;
  double len, diam = 1e10;
  double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[1]->GetCoords(x1, y1, z1);
  Verts[2]->GetCoords(x2, y2, z2);
  Verts[3]->GetCoords(x3, y3, z3);

  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len < diam) diam = len;
  len = (x0-x2)*(x0-x2)+(y0-y2)*(y0-y2)+(z0-z2)*(z0-z2);
  if(len < diam) diam = len;
  len = (x0-x3)*(x0-x3)+(y0-y3)*(y0-y3)+(z0-z3)*(z0-z3);
  if(len < diam) diam = len;
  len = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
  if(len < diam) diam = len;
  len = (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3);
  if(len < diam) diam = len;
  len = (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+(z2-z3)*(z2-z3);
  if(len < diam) diam = len;
  
  return sqrt(diam);
}

double TTetrahedron::GetLengthWithReferenceMap(TVertex **Verts)
{
  double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3;
  double xc0, xc1, xc2, xc3, yc0, yc1, yc2, yc3, zc0, zc1, zc2, zc3;
  double detjk;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[1]->GetCoords(x1, y1, z1);
  Verts[2]->GetCoords(x2, y2, z2);
  Verts[3]->GetCoords(x3, y3, z3);

  xc0=x0;
  xc1=x1-x0;
  xc2=x2-x0;
  xc3=x3-x0;

  yc0=y0;
  yc1=y1-y0;
  yc2=y2-y0;
  yc3=y3-y0;

  zc0=z0;
  zc1=z1-z0;
  zc2=z2-z0;
  zc3=z3-z0;

  detjk= fabs(xc1*yc2*zc3 + xc2*yc3*zc1 + xc3*yc1*zc2
        -xc3*yc2*zc1 - xc2*yc1*zc3 - xc1*yc3*zc2);
	return(pow(detjk,1.0/3.0));
}

double TTetrahedron::GetMeasure(TVertex **Verts)
{
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;

  Verts[0]->GetCoords(x1,y1,z1);
  Verts[1]->GetCoords(x2,y2,z2);
  Verts[2]->GetCoords(x3,y3,z3);
  Verts[3]->GetCoords(x4,y4,z4);

  return fabs(-x2*(y3*z4-y4*z3)+x3*(y2*z4-y4*z2)-x4*(y2*z3-y3*z2)+
	      x1*(y3*z4-y4*z3)-x3*(y1*z4-y4*z1)+x4*(y1*z3-y3*z1)+
	      -x1*(y2*z4-y4*z2)+x2*(y1*z4-y4*z1)-x4*(y1*z2-y2*z1)+
	      x1*(y2*z3-y3*z2)-x2*(y1*z3-y3*z1)+x3*(y1*z2-y2*z1))/6.;
}
