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
// @(#)Quadrangle.C        1.1 10/30/98
//
// Class:       TQuadrangle
// Purpose:     shape descriptor of a quadrangle
//
// Author:      Volker Behns  26.07.97
//
// =======================================================================

#include <Quadrangle.h>
#include <Constants.h>
#include <math.h>
#include <stdlib.h>

static const int DatEdgeVertex[][2] = { {0, 1},  {1, 2},  {2, 3},  {3, 0}};
static const int DatVertexEdge[][QUADMAXN_EpV] =
                 { {3, 0},  {0, 1},  {1, 2},  {2, 3}};

// Constructor
TQuadrangle::TQuadrangle()
{
  MaxN_EpV = QUADMAXN_EpV;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  Type = Quadrangle;
  N_Vertices = 4;
  N_Edges = 4;
  N_Joints = 4;
}

// Methods
double TQuadrangle::GetDiameter(TVertex **Verts)
{
  double diffX1 = Verts[0]->GetX() - Verts[2]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[2]->GetY();
  double diffX2 = Verts[1]->GetX() - Verts[3]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[3]->GetY();

  return sqrt(MAX(diffX1*diffX1 + diffY1*diffY1, 
                  diffX2*diffX2 + diffY2*diffY2));
}

double TQuadrangle::GetShortestEdge(TVertex **Verts)
{
  double diffX1 = Verts[0]->GetX() - Verts[1]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[1]->GetY();
  double diffX2 = Verts[1]->GetX() - Verts[2]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[2]->GetY();
  double diffX3 = Verts[2]->GetX() - Verts[3]->GetX();
  double diffY3 = Verts[2]->GetY() - Verts[3]->GetY();
  double diffX4 = Verts[3]->GetX() - Verts[0]->GetX();
  double diffY4 = Verts[3]->GetY() - Verts[0]->GetY();
  double len;

  if (diffX1*diffX1 + diffY1*diffY1<diffX2*diffX2 + diffY2*diffY2)
      len = diffX1*diffX1 + diffY1*diffY1;
  else
      len = diffX2*diffX2 + diffY2*diffY2;
  if (diffX3*diffX3 + diffY3*diffY3<len)
      len = diffX3*diffX3 + diffY3*diffY3;
  if (diffX4*diffX4 + diffY4*diffY4<len)
      len = diffX4*diffX4 + diffY4*diffY4;

  return sqrt(len);
}

// approximation with affine map of vertices 0,1 and 3
double TQuadrangle::GetLengthWithReferenceMap(TVertex **Verts)
{
    double x0, x1, x3, y0, y1, y3;
    double xc1, xc2, yc1, yc2;
    double detjk, rec_detjk;
    double d11, d12, d21, d22;

    x0 = Verts[0]->GetX();
    y0 = Verts[0]->GetY();
    x1 = Verts[1]->GetX();
    y1 = Verts[1]->GetY();
    x3 = Verts[3]->GetX();
    y3 = Verts[3]->GetY();
    

    //double xc0 = (x1 + x3) * 0.5;
    xc1 = (x1 - x0) * 0.5;
    xc2 = (x3 - x0) * 0.5;
    
    //double yc0 = (y1 + y3) * 0.5;
    yc1 = (y1 - y0) * 0.5;
    yc2 = (y3 - y0) * 0.5;
    
    detjk=xc1*yc2-xc2*yc1;
    rec_detjk=1/detjk;
    
    d11 = (y3-y0) * 0.5 * rec_detjk;  //dxi/dx
    d12 = (x0-x3) * 0.5 * rec_detjk;  //dxi/dy
    d21 = (y0-y1) * 0.5 * rec_detjk;  //deta/dx
    d22 = (x1-x0) * 0.5 * rec_detjk;  //deta/dy

    OutPut(d11 << " " << d12 << " " << d21 << " " << d22 << endl);
    OutPut("This gives the same result as GetMeasure()" << endl);
    exit(4711);

    // factor 2 = sqrt(4) because of area of unit square
    return 2*sqrt(fabs(detjk));
}

double TQuadrangle::GetMeasure(TVertex **Verts)
{
  double x1,x2,x3,x4,y1,y2,y3,y4;

  x1 = Verts[0]->GetX();
  y1 = Verts[0]->GetY();
  x2 = Verts[1]->GetX();
  y2 = Verts[1]->GetY();
  x3 = Verts[2]->GetX();
  y3 = Verts[2]->GetY();
  x4 = Verts[3]->GetX();
  y4 = Verts[3]->GetY();

  return 0.5*fabs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)+
    0.5*fabs(x4*y3-x3*y4-x1*y3+x3*y1+x1*y4-x4*y1);
}

Shapes TQuadrangle::CheckQuad(TVertex **Vertices)
{
  double test1, test2;

  test1 = Vertices[0]->GetX() - Vertices[1]->GetX() +
          Vertices[2]->GetX() - Vertices[3]->GetX();
  test2 = Vertices[0]->GetY() - Vertices[1]->GetY() +
          Vertices[2]->GetY() - Vertices[3]->GetY();

  if (ABS(test1) > 1e-8 || ABS(test2) > 1e-8) return Quadrangle;

// check for rectangle
//  test1 = (Vertices[0]->GetX() - Vertices[1]->GetX()) * 
//          (Vertices[2]->GetX() - Vertices[1]->GetX()) +
//          (Vertices[0]->GetY() - Vertices[1]->GetY()) * 
//          (Vertices[2]->GetY() - Vertices[1]->GetY());
//         
//  test2 = (Vertices[0]->GetX() - Vertices[2]->GetX()) * 
//          (Vertices[3]->GetX() - Vertices[1]->GetX()) +
//          (Vertices[0]->GetY() - Vertices[2]->GetY()) * 
//          (Vertices[3]->GetY() - Vertices[1]->GetY());
//         
//  if (ABS(test1) > 1e-8 || ABS(test2) > 1e-8) return Parallelogram;

  test1 = (Vertices[3]->GetX() - Vertices[0]->GetX()) *
          (Vertices[1]->GetX() - Vertices[0]->GetX()) +
          (Vertices[3]->GetY() - Vertices[0]->GetY()) *
          (Vertices[1]->GetY() - Vertices[0]->GetY());

  if (ABS(test1) > 1e-8) return Parallelogram;
  return Rectangle;
}
