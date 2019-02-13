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
// @(#)Hexahedron.C        1.4 11/15/99
//
// Class:       THexahedron
// Purpose:     shape descriptor of a hexahedron
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __3D__
#define __3D__
#endif

#include <Hexahedron.h>

static const int DatEdgeVertex[][2] = { {0, 1},  {1, 2},  {2, 3},  {3, 0},
                                        {0, 4},  {1, 5},  {2, 6},  {3, 7},
                                       {4, 5},  {5, 6},  {6, 7},  {7, 4}};
static const int DatVertexEdge[][HEXAMAXN_EpV] =
                 { {3, 0, 4},  {0, 1, 5},  {1, 2, 6},  {2, 3, 7},
                   {4, 8,11},  {5, 9, 8},  {6,10, 9},  {7,11,10}};

static const int DatFaceVertex[][HEXAMAXN_VpF] =
                 { {0, 1, 2, 3},  {0, 4, 5, 1},  {1, 5, 6, 2},
                   {2, 6, 7, 3},  {0, 3, 7, 4},  {4, 7, 6, 5}};
static const int DatFaceVertexLen[] = { 4, 4, 4, 4, 4, 4};

static const int DatVertexFace[][HEXAMAXN_FpV] =
                 { {0, 1, 4},   {0, 2, 1},   {0, 3, 2},   {0, 4, 3},
                   {1, 5, 4},   {2, 5, 1},   {3, 5, 2},   {4, 5, 3}};

static const int DatFaceEdge[][HEXAMAXN_EpF] =
                 { {0, 1, 2, 3},  {4, 8, 5, 0},  {5, 9, 6, 1},
                   {6,10, 7, 2},  {3, 7,11, 4},  {11,10, 9,8}};

static const int DatFaceEdgeLen[] =
                 { 4, 4, 4, 4, 4, 4};
static const int DatEdgeFace[][HEXAMAXN_FpE] =
                 { {1, 0},  {2, 0},  {3, 0},  {4, 0},
                   {4, 1},  {1, 2},  {2, 3},  {3, 4},
                   {5, 1},  {5, 2},  {5, 3},  {5, 4}};

static const Shapes DatFaceType[] = { Quadrangle, Quadrangle, Quadrangle,
                                      Quadrangle, Quadrangle, Quadrangle};

// Constructor
THexahedron::THexahedron()
{
  MaxN_EpV = HEXAMAXN_EpV;
  MaxN_VpF = HEXAMAXN_VpF;
  MaxN_FpV = HEXAMAXN_FpV;
  MaxN_EpF = HEXAMAXN_EpF;
  MaxN_FpE = HEXAMAXN_FpE;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  FaceVertex = (const int *) DatFaceVertex;
  FaceVertexLen = (const int *) DatFaceVertexLen;

  VertexFace = (const int *) DatVertexFace;

  FaceEdge = (const int *) DatFaceEdge;
  FaceEdgeLen = (const int *) DatFaceEdgeLen;
  EdgeFace = (const int *) DatEdgeFace;

  FaceType = (const Shapes *) DatFaceType;

  Type = Hexahedron;
  N_Vertices = 8;
  N_Edges = 12;
  N_Faces = 6;
  N_Joints = 6;
}

// Methods
double THexahedron::GetDiameter(TVertex **Verts)
{
  double len, diam = 0;
  double x0, x1, y0, y1, z0, z1;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[6]->GetCoords(x1, y1, z1);
  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len>diam) diam = len;

  Verts[1]->GetCoords(x0, y0, z0);
  Verts[7]->GetCoords(x1, y1, z1);
  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len>diam) diam = len;

  Verts[2]->GetCoords(x0, y0, z0);
  Verts[4]->GetCoords(x1, y1, z1);
  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len>diam) diam = len;

  Verts[3]->GetCoords(x0, y0, z0);
  Verts[5]->GetCoords(x1, y1, z1);
  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len>diam) diam = len;

  return sqrt(diam);
}

// Methods, exact only for parallelepiped
double THexahedron::GetShortestEdge(TVertex **Verts)
{
  double len, diam = 1e10;
  double x0, x1, y0, y1, z0, z1;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[4]->GetCoords(x1, y1, z1);
  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len<diam) diam = len;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[1]->GetCoords(x1, y1, z1);
  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len<diam) diam = len;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[3]->GetCoords(x1, y1, z1);
  len = (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
  if(len<diam) diam = len;

  return sqrt(diam);
}

// Method exact only for parallelepipeds
double THexahedron::GetLengthWithReferenceMap(TVertex **Verts)
{
  double x0, x1, x3, x4, y0, y1, y3, y4, z0, z1, z3, z4;
  double xc0, xc1, xc2, xc3, yc0, yc1, yc2, yc3, zc0, zc1, zc2, zc3;
  double detjk;

  Verts[0]->GetCoords(x0, y0, z0);
  Verts[1]->GetCoords(x1, y1, z1);
  Verts[3]->GetCoords(x3, y3, z3);
  Verts[4]->GetCoords(x4, y4, z4);

  //x0 = 0.5; x1 = 1.4; x3 = 9.8, x4 = -12;
  //y0 = -0.5; y1 = 2.3; y3 = -1.0, y4 = 1;
  //z0 = 0; z1 = -1.4; z3 = 7.2, z4 = -7.2;

  xc0 = (x1 + x3 + x4 - x0) * 0.5;
  xc1 = (x1 - x0) * 0.5;
  xc2 = (x3 - x0) * 0.5;
  xc3 = (x4 - x0) * 0.5;

  yc0 = (y1 + y3 + y4 - y0) * 0.5;
  yc1 = (y1 - y0) * 0.5;
  yc2 = (y3 - y0) * 0.5;
  yc3 = (y4 - y0) * 0.5;

  zc0 = (z1 + z3 + z4 - z0) * 0.5;
  zc1 = (z1 - z0) * 0.5;
  zc2 = (z3 - z0) * 0.5;
  zc3 = (z4 - z0) * 0.5;

  detjk= fabs(xc1*yc2*zc3 + xc2*yc3*zc1 + xc3*yc1*zc2 - xc3*yc2*zc1 - 
	      xc2*yc1*zc3 - xc1*yc3*zc2);

  // factor 2 = 8^(1/3) because of area of unit cube
  return(2*pow(detjk,1.0/3.0));
}

double THexahedron::GetMeasure(TVertex **Verts)
{
  double x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8,z1,z2,z3,z4,z5,z6,z7,z8;

  Verts[0]->GetCoords(x1,y1,z1);
  Verts[1]->GetCoords(x2,y2,z2);
  Verts[2]->GetCoords(x3,y3,z3);
  Verts[3]->GetCoords(x4,y4,z4);
  Verts[4]->GetCoords(x5,y5,z5);
  Verts[5]->GetCoords(x6,y6,z6);
  Verts[6]->GetCoords(x7,y7,z7);
  Verts[7]->GetCoords(x8,y8,z8);

  // hexahedron decomposed into six tetrahedra 
  return fabs(-x2*(y5*z4-y4*z5)+x5*(y2*z4-y4*z2)-x4*(y2*z5-y5*z2)+
	      x1*(y5*z4-y4*z5)-x5*(y1*z4-y4*z1)+x4*(y1*z5-y5*z1)+
	      -x1*(y2*z4-y4*z2)+x2*(y1*z4-y4*z1)-x4*(y1*z2-y2*z1)+
	      x1*(y2*z5-y5*z2)-x2*(y1*z5-y5*z1)+x5*(y1*z2-y2*z1))/6.+
    fabs(-x2*(y3*z4-y4*z3)+x3*(y2*z4-y4*z2)-x4*(y2*z3-y3*z2)+
	 x7*(y3*z4-y4*z3)-x3*(y7*z4-y4*z7)+x4*(y7*z3-y3*z7)+
	 -x7*(y2*z4-y4*z2)+x2*(y7*z4-y4*z7)-x4*(y7*z2-y2*z7)+
	 x7*(y2*z3-y3*z2)-x2*(y7*z3-y3*z7)+x3*(y7*z2-y2*z7))/6.+
    fabs(-x5*(y6*z7-y7*z6)+x6*(y5*z7-y7*z5)-x7*(y5*z6-y6*z5)+
	 x2*(y6*z7-y7*z6)-x6*(y2*z7-y7*z2)+x7*(y2*z6-y6*z2)+
	 -x2*(y5*z7-y7*z5)+x5*(y2*z7-y7*z2)-x7*(y2*z5-y5*z2)+
	 x2*(y5*z6-y6*z5)-x5*(y2*z6-y6*z2)+x6*(y2*z5-y5*z2))/6.+
    fabs(-x5*(y7*z8-y8*z7)+x7*(y5*z8-y8*z5)-x8*(y5*z7-y7*z5)+
	 x4*(y7*z8-y8*z7)-x7*(y4*z8-y8*z4)+x8*(y4*z7-y7*z4)+
	 -x4*(y5*z8-y8*z5)+x5*(y4*z8-y8*z4)-x8*(y4*z5-y5*z4)+
	 x4*(y5*z7-y7*z5)-x5*(y4*z7-y7*z4)+x7*(y4*z5-y5*z4))/6.+
    fabs(-x4*(y5*z7-y7*z5)+x5*(y4*z7-y7*z4)-x7*(y4*z5-y5*z4)+
	 x2*(y5*z7-y7*z5)-x5*(y2*z7-y7*z2)+x7*(y2*z5-y5*z2)+
	 -x2*(y4*z7-y7*z4)+x4*(y2*z7-y7*z2)-x7*(y2*z4-y4*z2)+
	 x2*(y4*z5-y5*z4)-x4*(y2*z5-y5*z2)+x5*(y2*z4-y4*z2))/6.;
}

Shapes THexahedron::CheckHexa(TVertex **Vertices)
{
  double xt1, xt2, xt3, xt4;
  double yt1, yt2, yt3, yt4;
  double zt1, zt2, zt3, zt4;
  double x0, x1, x2, x3, x4, x5, x6, x7;
  double y0, y1, y2, y3, y4, y5, y6, y7;
  double z0, z1, z2, z3, z4, z5, z6, z7;

  Shapes ret = Hexahedron;

  Vertices[0]->GetCoords(x0, y0, z0);
  Vertices[1]->GetCoords(x1, y1, z1);
  Vertices[2]->GetCoords(x2, y2, z2);
  Vertices[3]->GetCoords(x3, y3, z3);
  Vertices[4]->GetCoords(x4, y4, z4);
  Vertices[5]->GetCoords(x5, y5, z5);
  Vertices[6]->GetCoords(x6, y6, z6);
  Vertices[7]->GetCoords(x7, y7, z7);

  xt1 = fabs( x0 - x1 - x3 + x2 + x4 - x5 - x7 + x6);
  xt2 = fabs( x0 - x1 + x3 - x2 - x4 + x5 - x7 + x6);
  xt3 = fabs( x0 + x1 - x3 - x2 - x4 - x5 + x7 + x6);
  xt4 = fabs(-x0 + x1 + x3 - x2 + x4 - x5 - x7 + x6);

  yt1 = fabs( y0 - y1 - y3 + y2 + y4 - y5 - y7 + y6);
  yt2 = fabs( y0 - y1 + y3 - y2 - y4 + y5 - y7 + y6);
  yt3 = fabs( y0 + y1 - y3 - y2 - y4 - y5 + y7 + y6);
  yt4 = fabs(-y0 + y1 + y3 - y2 + y4 - y5 - y7 + y6);

  zt1 = fabs( z0 - z1 - z3 + z2 + z4 - z5 - z7 + z6);
  zt2 = fabs( z0 - z1 + z3 - z2 - z4 + z5 - z7 + z6);
  zt3 = fabs( z0 + z1 - z3 - z2 - z4 - z5 + z7 + z6);
  zt4 = fabs(-z0 + z1 + z3 - z2 + z4 - z5 - z7 + z6);

  if( (xt1+xt2+xt3+xt4 < 1e-8) &&
      (yt1+yt2+yt3+yt4 < 1e-8) &&
      (zt1+zt2+zt3+zt4 < 1e-8) )
    ret = Brick;

  return ret;
}
