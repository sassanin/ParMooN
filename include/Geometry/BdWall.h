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
// @(#)BdWall.h        1.1 07/16/99
//
// Class:       TBdWall
// Superclass:  TBoundComp
// Purpose:     a Wall as a component of a boundary part
//
// Author:      Volker Behns  05.07.99
//
// =======================================================================

#ifndef __BDWALL__
#define __BDWALL__

#include <BoundComp2D.h>
#include <BoundComp3D.h>

/** a Wall as a component of a boundary part */
class TBdWall : public TBoundComp3D
{
  protected:
    /** 2d boundary component */
    TBoundComp2D * BdComp2D;

    /** drift in x-direction */
    double DriftX;
    /** drift in y-direction */
    double DriftY;
    /** drift in z-direction */
    double DriftZ;

  public:
    // Constructor
    TBdWall(int id, TBoundComp2D *bdcomp2d);
    
    // Methods
    int SetParams(double drx, double dry, double drz);

    /** return the coordinates of parameter value T, S */
    virtual int GetXYZofTS(double T, double S,
                           double &X, double &Y, double &Z);

    /** return the parameter value T, S of coordinates */
    virtual int GetTSofXYZ(double X, double Y, double Z,
                           double &T, double &S);

    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S);

    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat);

    /** return BdComp2D */
    TBoundComp2D* GetBdComp2D()
    { return BdComp2D; }
};

#endif
