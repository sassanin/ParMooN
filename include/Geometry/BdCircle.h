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
// @(#)BdCircle.h        1.2 07/16/99
//
// Class:       TBdCircle
// Superclass:  TBoundComp
// Purpose:     a part of a circle as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BDCIRCLE__
#define __BDCIRCLE__

#include <BoundComp2D.h>

/** a part of a circle as a component of a boundary part */
class TBdCircle : public TBoundComp2D
{
  protected:
    /** x coordinate of midpoint */
    double Xmid;
    /** y coordinate of midpoint */
    double Ymid;
    /** radii of the arc */
    double Radius_a, Radius_b;
    /** begin angle */
    double Phi1;
    /** end angle */
    double Phi2;

  public:
    // Constructor
    TBdCircle(int id);

    // Methods
    /** set all parameters to the given values */
    void SetParams (double xmid, double ymid, double radius_a, 
                    double radius_b, double phi1, double phi2);

    /** return the coordinates of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y);

    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T);

    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat);

    /** get number of initial vertices on this component */
    virtual int GetN_InitVerts();
    virtual int GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges);

  protected:
    int GetN_InitVertsSub(double Phi_a, double Phi_b, int Level);
    int GenInitVertsSub(double Phi_a, double Phi_b, int Level,
                        double *&points, int &I_points,
                        int *&edges, int &I_edges);
};

#endif
