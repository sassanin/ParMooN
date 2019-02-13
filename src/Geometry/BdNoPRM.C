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
//
// Class:       TBdNoPRM
// Purpose:     3D domain without PRM file
//
// Author:      Volker John 2008/01/28
//
// =======================================================================

#include <BdNoPRM.h>
#include <MooNMD_Io.h>

// Constructor
TBdNoPRM::TBdNoPRM(int id) : TBoundComp3D(id)
{
  Type = NoPRM;
}

// Methods
// there are no params
void TBdNoPRM::SetParams ()
{
    ;
}

// set dummy values
int TBdNoPRM::GetXYZofTS(double T, double S,
                        double &X, double &Y, double &Z)
{
    X = -4711;
    Y = -4711;
    Z = -4711;

  return 0;
}

/** return parameters and coordinates of a given linear
    combination of vertices */
// set dummy values
int TBdNoPRM::GetXYZandTS(int N_Points, double *LinComb,
                          double *xp, double *yp, double *zp,
                          double *tp, double *sp,
                          double &X, double &Y, double &Z,
                          double &T, double &S)
{
    int i;
    double t;

    X = Y = Z = 0.0;
    for (i=0;i<N_Points;i++)
    {
	t = LinComb[i]; 
	X += t * xp[i];
	Y += t * yp[i];
	Z += t * zp[i];
    }
    T = -4711;
    S = -4711;
    return 0;
}

// set dummu valus
int TBdNoPRM::GetTSofXYZ(double X, double Y, double Z,
                        double &T, double &S)
{
  double TS_aux;

  T = -4711;
  S = -4711;
  return 0;
}

// no readin
int TBdNoPRM::ReadIn(std::ifstream &dat)
{
  return 0;
}
