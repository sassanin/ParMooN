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
// @(#)Rectangle.C        1.1 10/30/98
//
// Class:       TRectangle
// Purpose:     shape descriptor of a quadrangle, especially a rectangle
//
// Author:      Volker Behns  27.03.98
//
// =======================================================================

#include <Rectangle.h>
#include <math.h>
#include <stdlib.h>

// Constructor
TRectangle::TRectangle() : TQuadrangle()
{
  Type = Rectangle;
}

// Methods
double TRectangle::GetDiameter(TVertex **Verts)
{
  double diffX = Verts[0]->GetX() - Verts[2]->GetX();
  double diffY = Verts[0]->GetY() - Verts[2]->GetY();

  return sqrt(diffX*diffX + diffY*diffY);
}
double TRectangle::GetShortestEdge(TVertex **Verts)
{
  double diffX1 = Verts[0]->GetX() - Verts[1]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[1]->GetY();
  double diffX2 = Verts[1]->GetX() - Verts[2]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[2]->GetY();
  double len;

  if (diffX1*diffX1 + diffY1*diffY1<diffX2*diffX2 + diffY2*diffY2)
      len = sqrt(diffX1*diffX1 + diffY1*diffY1);
  else
      len = sqrt(diffX2*diffX2 + diffY2*diffY2);

  return len;
}

double TRectangle::GetLengthWithReferenceMap(TVertex **Verts)
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
  
  x0 = 3;
  y0 = 8;
  x1 = 9;
  y1 = 12;
  x3 = 11;
  y3 = 10; 

  /*  x0 = -1;
  y0 = -1;
  x1 = 1;
  y1 = -1;
  x3 = -1;
  y3 = 1; */

  // double xc0 = (x1 + x3) * 0.5;
  xc1 = (x1 - x0) * 0.5;
  xc2 = (x3 - x0) * 0.5;
  
  // double yc0 = (y1 + y3) * 0.5;
  yc1 = (y1 - y0) * 0.5;
  yc2 = (y3 - y0) * 0.5;
  
  detjk=xc1*yc2-xc2*yc1;
  rec_detjk=1/detjk;
  
  d11 = (y3-y0) * 0.5 * rec_detjk;  //dxi/dx
  d12 = (x0-x3) * 0.5 * rec_detjk;  //dxi/dy
  d21 = (y0-y1) * 0.5 * rec_detjk;  //deta/dx
  d22 = (x1-x0) * 0.5 * rec_detjk;  //deta/dy

  OutPut(detjk << " " << d11 << " " << d12 << " " << d21 << " " << d22 << 
	 " " <<  2*sqrt(fabs(detjk)) << endl);
  OutPut("This gives the same result as GetMeasure()" << endl);
  exit(4711);
  // factor 2 = sqrt(4) because of area of unit square

  return 2*sqrt(fabs(detjk));
}

// measure of parallelogramm with vector product
double TRectangle::GetMeasure(TVertex **Verts)
{
    double area, x[4], y[4];

    x[0] = Verts[0]->GetX();
    y[0] = Verts[0]->GetY();
    x[1] = Verts[1]->GetX();
    y[1] = Verts[1]->GetY();
    x[3] = Verts[3]->GetX();
    y[3] = Verts[3]->GetY();

    area = fabs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));

    return(area);
}
