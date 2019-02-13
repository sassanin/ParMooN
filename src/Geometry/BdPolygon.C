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
// @(#)BdPolygon.C        1.1 08/12/99
//
// Class:       BdPolygon
// Purpose:     component is polygon
//
// Author:      Gunar Matthies 04.08.1999
//
// =======================================================================

#include <BdPolygon.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
// #include <malloc.h>

// Constructor
TBdPolygon::TBdPolygon(int id, int n_points) : TBoundComp2D(id)
{
  Type = Polygon;

  N_Points = n_points;
  Coords = NULL;
}

// Methods
void TBdPolygon::SetParams (int n_points, double *coords)
{
  N_Points = n_points;
  Coords = coords;
}

int TBdPolygon::GetXYofT(double T, double &X, double &Y)
{
  int N_Intervals, N1, N2;
  double S, D;

  N_Intervals = N_Points-1;
  S = T*N_Intervals;
  N1 = (int)S;
  N2 = N1+1;
  D = S-N1;

  if(N2>N_Intervals)
  {
    // T == 1
    X = Coords[2*N1];
    Y = Coords[2*N1+1];
  }
  else
  {
    X = Coords[2*N1]*(1-D)+Coords[2*N2]*D;
    Y = Coords[2*N1+1]*(1-D)+Coords[2*N2+1]*D;
  }

  return 0;
}

int TBdPolygon::GetTofXY(double X, double Y, double &T)
{
  int i, j; 
  double t1, t2, x1, x2, y1, y2;

  j = -1;
  for(i=1;i<N_Points;i++)
  {
    x1 = Coords[2*i-2];
    y1 = Coords[2*i-1];
    x2 = Coords[2*i];
    y2 = Coords[2*i+1];

    if( fabs(x2-x1) > fabs(y2-y1) )
    {
      t1 = (X-x1)/(x2-x1);
      if( fabs(t1-0.5)<=0.5 && fabs(Y-(y1+t1*(y2-y1))) < (1e-6) )
      {
        t2 = t1;
        j = i;
        break;
      }
    }
    else
    {
      t2 = (Y-y1)/(y2-y1);
      if( fabs(t2-0.5)<=0.5 && fabs(X-(x1+t2*(x2-x1))) < (1e-6) )
      {
        t1 = t2;
        j = i;
        break;
      }
    }
  } // endfor i

  if(j==-1)
  {
    i = -1;
    T = 0;
  }
  else
  {
    i = 0;
    T = t1;
  }

  return i;
}

int TBdPolygon::ReadIn(std::ifstream &dat)
{
  char line[100];
  int k;

  Coords = new double [2 * N_Points];

  for(k=0;k<N_Points;k++)
  {
    dat >> Coords[2*k] >> Coords[2*k+1];
    dat.getline (line, 99);
  }
             
  return 0;
}
