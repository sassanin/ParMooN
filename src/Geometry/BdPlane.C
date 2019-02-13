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
// @(#)BdPlane.C        1.2 07/16/99
//
// Class:       TBdPlane
// Purpose:     a plane as a component of a boundary part
//
// Author:      Volker Behns  05.07.99
//
// =======================================================================

#include <BdPlane.h>
#include <MooNMD_Io.h>

// Constructor
TBdPlane::TBdPlane(int id) : TBoundComp3D(id)
{
  Type = Plane;
}

// Methods
void TBdPlane::SetParams (double p_x, double p_y, double p_z,
                         double a_x, double a_y, double a_z,
                         double n_x, double n_y, double n_z)
{
  P_x = p_x;
  P_y = p_y;
  P_z = p_z;
  A_x = a_x;
  A_y = a_y;
  A_z = a_z;
  N_x = n_x;
  N_y = n_y;
  N_z = n_z;

  B_x = A_y * N_z - A_z * N_y;
  B_y = A_z * N_x - A_x * N_z;
  B_z = A_x * N_y - A_y * N_x;
}

int TBdPlane::GetXYZofTS(double T, double S,
                        double &X, double &Y, double &Z)
{
  X = P_x + T * A_x + S * B_x;
  Y = P_y + T * A_y + S * B_y;
  Z = P_z + T * A_z + S * B_z;
//   cout <<" TBdPlane::GetXYZofTS " <<endl;
  return 0;
}

/** return parameters and coordinates of a given linear
    combination of vertices */
int TBdPlane::GetXYZandTS(int N_Points, double *LinComb,
                          double *xp, double *yp, double *zp,
                          double *tp, double *sp,
                          double &X, double &Y, double &Z,
                          double &T, double &S)
{
  int i;
  double t, s, v;

  t = s = 0;
  for(i=0;i<N_Points;i++)
  {
    v = LinComb[i];
    t += v*tp[i];
    s += v*sp[i];
  }
  T = t;
  S = s;

//   X = P_x + T * A_x + S * B_x;
//   Y = P_y + T * A_y + S * B_y;
//   Z = P_z + T * A_z + S * B_z;

  X = Y = Z = 0;
  for (i=0;i<N_Points;++i)
  {
    X += LinComb[i]*xp[i];
    Y += LinComb[i]*yp[i];
    Z += LinComb[i]*zp[i];
  }
//   cout <<" TBdPlane::GetXYZofTS " <<endl;
  return 0;
}

int TBdPlane::GetTSofXYZ(double X, double Y, double Z,
                        double &T, double &S)
{
  double TS_aux;

  TS_aux = A_x * B_y - A_y * B_x;
  if (ABS(TS_aux) > 1e-10)
    S = ((Y - P_y)*A_x - (X - P_x)*A_y)/TS_aux;
  else
  {
    TS_aux = A_x * B_z - A_z * B_x;
    if (ABS(TS_aux) > 1e-10)
      S = ((Z - P_z)*A_x - (X - P_x)*A_z)/TS_aux;
    else
      S = ((Y - P_y)*A_z - (Z - P_z)*A_y)/(A_z * B_y - A_y * B_z);
  }

  if (ABS(A_x) > 1e-10)
    T = (X - P_x - S*B_x)/A_x;
  else
  {
    if (ABS(A_y) > 1e-10)
      T = (Y - P_y - S*B_y)/A_y;
    else
      T = (Z - P_z - S*B_z)/A_z;
  }
//   cout <<" TBdPlane::GetXYZofTS " <<endl;
  return 0;
}

int TBdPlane::ReadIn(std::ifstream &dat)
{
  char line[100];

  dat >> P_x >> P_y >> P_z;
  dat.getline (line, 99);
  dat >> A_x >> A_y >> A_z;
  dat.getline (line, 99);
  dat >> N_x >> N_y >> N_z;
  dat.getline (line, 99);

  B_x = A_y * N_z - A_z * N_y;
  B_y = A_z * N_x - A_x * N_z;
  B_z = A_x * N_y - A_y * N_x;

  return 0;
}

void TBdPlane::GetParams (double &p_x, double &p_y, double &p_z,
                    double &a_x, double &a_y, double &a_z,
                    double &n_x, double &n_y, double &n_z)
{
  p_x = P_x;
  p_y = P_y;
  p_z = P_z;
  
  a_x = A_x;
  a_y = A_y;
  a_z = A_z;
  
  n_x = N_x;
  n_y = N_y;
  n_z = N_z;
}
