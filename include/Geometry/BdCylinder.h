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
// @(#)BdSphere.h        1.1 07/16/99
//
// Class:       TBdSphere
// Superclass:  TBoundComp3D
// Purpose:     a Cylinder as a component of a boundary part
//
// Author:      Andreas Hahn 16.04.2010
//
// =======================================================================

#ifndef __BDCYLINDER__
#define __BDCYLINDER__

#include <BoundComp3D.h>

class TBdCylinder : public TBoundComp3D
{
  protected:
    /** radius **/
    double mRadius;
    
    /** coordinates of one point on axis **/
    double mPx, mPy, mPz;
    
    /** direction of axis **/
    double mAx, mAy, mAz;
    
    /** **/
    double mNx, mNy, mNz;
    
    /** **/
    double mBx, mBy, mBz;
    
  protected:
    
  public:
    // CTOR
    TBdCylinder (int id);
    
    virtual ~TBdCylinder () {};
    
    // Methods
    /** return the coordinates {X, Y, Z} of parameter values T and S */
    virtual int GetXYZofTS(double T, double S, double &X, double &Y,
                          double &Z);
    
     /** return the parameter values T and S of coordinates (X, Y, Z) */
    virtual int GetTSofXYZ(double X, double Y, double Z, double &T,
                           double &S);
			   
    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S);
			  
    virtual int ReadIn(std::ifstream &dat);
    
    void SetParams(double r, double px, double py, double pz,
		   double ax, double ay, double az, double nx, double ny, double nz);
};

#endif
