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
// @(#)BdSpline.C        1.2 07/16/99
//
// Class:       BdSpline
// Purpose:     splie function as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#include <BdSpline.h>
#include <math.h>

// Constructor
TBdSpline::TBdSpline(int id, int N_Spls) : TBoundComp2D(id)
{
  Type = Spline2D;
  N_Splines = N_Spls;
  Params = new double [10*N_Splines];
}

// Methods
void TBdSpline::SetParams (double *params)
{
  Params = params;
}

int TBdSpline::GetN_Splines ()
{
  return N_Splines;
}

int TBdSpline::GetXYofT(double T, double &X, double &Y)
{
  double phi1, phi2, phi3, phi4;
  int ISpline;

  T *= N_Splines;
  ISpline = (int) T;
  T -= ISpline;
  
  phi1 = (2.*T*T - 3.*T)*T + 1.;
  phi2 = (-2.*T + 3.)*T*T;
  phi3 = (T*T - 2.*T + 1.)*T;
  phi4 = (T - 1)*T*T;

  if (ISpline == N_Splines)
  {
    ISpline -= 1;
    T = 1.;
  }

  ISpline *= 8;

  X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
      Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
  Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
      Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

  return 0;
}

int TBdSpline::GetTofXY(double X, double Y, double &T)
{
  return -1;
}

int TBdSpline::ReadIn(std::ifstream &dat)
{
  char line[100];
  int k, N_Params;

  N_Params = 4 * N_Splines;

  Params = new double [2*N_Params];

  for(k=0;k<N_Params;k++)
  {
    dat >> Params[2*k] >> Params[2*k+1];
    dat.getline (line, 99);
  }
             
  return 0;
}
