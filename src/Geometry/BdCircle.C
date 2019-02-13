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
// @(#)BdCircle.C        1.2 07/16/99
//   
// Class:       TBdCricle
// Purpose:     a part of a circle as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#include <BdCircle.h>
#include <MooNMD_Io.h>

// Constructor
TBdCircle::TBdCircle(int id) : TBoundComp2D(id)
{
  Type = Circle;
}

// Methods
void TBdCircle::SetParams (double xmid, double ymid, double radius_a,
                           double radius_b, double phi1, double phi2)
{
  Xmid = xmid;
  Ymid = ymid;
  Radius_a = radius_a;
  Radius_b = radius_b;
  Phi1 = phi1;
  Phi2 = phi2;
}

int TBdCircle::ReadIn(std::ifstream &dat)
{
  char line[100];

  dat >> Xmid >> Ymid;
  dat.getline (line, 99);
  dat >> Radius_a >> Radius_b;
  if (Radius_b < 1e-20) Radius_b = Radius_a;
  dat.getline (line, 99);
  dat >> Phi1 >> Phi2;
  dat.getline (line, 99);

  return 0;
}

int TBdCircle::GetXYofT(double T, double &X, double &Y)
{
  double Phi;

  if (T < -0.001 || T > 1.001)
  { 
    cerr << "Error: parameter T out of range " << T << endl;
     return -2;
  }

  Phi = Phi1 + T * (Phi2 - Phi1);
  X = Xmid + Radius_a*cos(Phi);
  Y = Ymid + Radius_b*sin(Phi);

  return 0;
}

int TBdCircle::GetTofXY(double X, double Y, double &T)
{
  T = atan2((Y-Ymid)/Radius_b, (X-Xmid)/Radius_a);
  if (T < MIN(Phi1, Phi2) - 1e-5) T += 2*Pi;
  T = (T - Phi1) / (Phi2 - Phi1);

  if (T < -0.001 || T > 1.001)
    return -1;
  else
    return 0;
}

int TBdCircle::GetN_InitVertsSub(double Phi_a, double Phi_b, int Level)
{
  int j;
  double Phi_c, L_arc, L_sec;
  double Phi, X1, Y1, X2, Y2;
  static double w[]={0.2777777777777778,0.4444444444444444,0.2777777777777778};
  static double x[]={0.11270166537925831149, 0.5, 0.88729833462074168851};
  static int N_InitVerts = 0;
  double R_aa = Radius_a*Radius_a;
  double R_bb = Radius_b*Radius_b;
  
  if (R_aa*sin(Phi_a)*sin(Phi_a) + R_bb*cos(Phi_a)*cos(Phi_a) >
      R_aa*sin(Phi_b)*sin(Phi_b) + R_bb*cos(Phi_b)*cos(Phi_b))
    Phi_c = .6*Phi_a + .4*Phi_b;
  else
    Phi_c = .4*Phi_a + .6*Phi_b;

  if (!Level) N_InitVerts = 2;

  N_InitVerts++;

  for (L_arc=j=0;j<3;j++)
  {
    Phi = Phi_a + x[j]*(Phi_c - Phi_a);
    L_arc += w[j]*sqrt(R_bb*cos(Phi)*cos(Phi) + R_aa*sin(Phi)*sin(Phi));
  }

  L_arc *= ABS(Phi_c - Phi_a);

  X1 = Xmid + Radius_a*cos(Phi_a);
  Y1 = Ymid + Radius_b*sin(Phi_a);
  X2 = Xmid + Radius_a*cos(Phi_c);
  Y2 = Ymid + Radius_b*sin(Phi_c);

  L_sec = sqrt((X1 - X2)*(X1 - X2) + (Y1 - Y2)*(Y1 - Y2));

  if (TOL_SECANT_BOUND * L_sec < L_arc && Level < RECURS_DEPTH)
    GetN_InitVertsSub(Phi_a, Phi_c, Level + 1);

  for (L_arc=j=0;j<3;j++)
  {
    Phi = Phi_c + x[j]*(Phi_b - Phi_c);
    L_arc += w[j]*sqrt(R_bb*cos(Phi)*cos(Phi) + R_aa*sin(Phi)*sin(Phi));
  }

  L_arc *= ABS(Phi_b - Phi_c);

  X1 = Xmid + Radius_a*cos(Phi_c);
  Y1 = Ymid + Radius_b*sin(Phi_c);
  X2 = Xmid + Radius_a*cos(Phi_b);
  Y2 = Ymid + Radius_b*sin(Phi_b);

  L_sec = sqrt((X1 - X2)*(X1 - X2) + (Y1 - Y2)*(Y1 - Y2));

  if (TOL_SECANT_BOUND * L_sec < L_arc && Level < RECURS_DEPTH)
    GetN_InitVertsSub(Phi_c, Phi_b, Level + 1);

  return N_InitVerts;
}

int TBdCircle::GetN_InitVerts()
{
  return GetN_InitVertsSub(Phi1, Phi2, 0);
}

int TBdCircle::GenInitVertsSub(double Phi_a, double Phi_b, int Level,
                               double *&points, int &I_points,
                               int *&edges, int &I_edges)
{
  int j;
  double Phi_c, L_arc, L_sec;
  double Phi, X1, Y1, X2, Y2;
  static double w[]={0.2777777777777778,0.4444444444444444,0.2777777777777778};
  static double x[]={0.11270166537925831149, 0.5, 0.88729833462074168851};
  double R_aa = Radius_a*Radius_a;
  double R_bb = Radius_b*Radius_b;
  
  if (R_aa*sin(Phi_a)*sin(Phi_a) + R_bb*cos(Phi_a)*cos(Phi_a) >
      R_aa*sin(Phi_b)*sin(Phi_b) + R_bb*cos(Phi_b)*cos(Phi_b))
    Phi_c = .6*Phi_a + .4*Phi_b;
  else
    Phi_c = .4*Phi_a + .6*Phi_b;

  for (L_arc=j=0;j<3;j++)
  {
    Phi = Phi_a + x[j]*(Phi_c - Phi_a);
    L_arc += w[j]*sqrt(R_bb*cos(Phi)*cos(Phi) + R_aa*sin(Phi)*sin(Phi));
  }

  L_arc *= ABS(Phi_c - Phi_a);

  X1 = Xmid + Radius_a*cos(Phi_a);
  Y1 = Ymid + Radius_b*sin(Phi_a);
  X2 = Xmid + Radius_a*cos(Phi_c);
  Y2 = Ymid + Radius_b*sin(Phi_c);

  L_sec = sqrt((X1 - X2)*(X1 - X2) + (Y1 - Y2)*(Y1 - Y2));

  if (TOL_SECANT_BOUND * L_sec < L_arc && Level < RECURS_DEPTH)
    GenInitVertsSub(Phi_a, Phi_c, Level + 1, points, I_points,
                    edges, I_edges);
  else
  {
    points[2*I_points    ] = X2;
    points[2*I_points + 1] = Y2;

    edges[2*I_edges    ] = I_points - 1;
    edges[2*I_edges + 1] = I_points++;
    I_edges++;
  }

  for (L_arc=j=0;j<3;j++)
  {
    Phi = Phi_c + x[j]*(Phi_b - Phi_c);
    L_arc += w[j]*sqrt(R_bb*cos(Phi)*cos(Phi) + R_aa*sin(Phi)*sin(Phi));
  }

  L_arc *= ABS(Phi_b - Phi_c);

  X1 = Xmid + Radius_a*cos(Phi_c);
  Y1 = Ymid + Radius_b*sin(Phi_c);
  X2 = Xmid + Radius_a*cos(Phi_b);
  Y2 = Ymid + Radius_b*sin(Phi_b);

  L_sec = sqrt((X1 - X2)*(X1 - X2) + (Y1 - Y2)*(Y1 - Y2));

  if (TOL_SECANT_BOUND * L_sec < L_arc && Level < RECURS_DEPTH)
    GenInitVertsSub(Phi_c, Phi_b, Level + 1, points, I_points,
                    edges, I_edges);
  else
  {
    points[2*I_points    ] = X2;
    points[2*I_points + 1] = Y2;

    edges[2*I_edges    ] = I_points - 1;
    edges[2*I_edges + 1] = I_points++;
    I_edges++;
  }

  return 0;
}

int  TBdCircle::GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges)
{
  points[2*I_points    ] = Xmid + Radius_a*cos(Phi1);
  points[2*I_points + 1] = Ymid + Radius_b*sin(Phi1);
  I_points++;

  return GenInitVertsSub(Phi1, Phi2, 0, points, I_points, edges, I_edges);
}
