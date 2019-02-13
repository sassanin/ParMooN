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
// @(#)BoundComp3D.h        1.5 08/12/99
//
// Class:       TBoundComp3D
// Purpose:     components of boundary faces
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BOUNDCOMP3D__
#define __BOUNDCOMP3D__

#include <BoundComp.h>

/** components of boundary faces */
class TBoundComp3D : public TBoundComp
{
  public:
    // Constructor
    TBoundComp3D(int id);

    // Methods
    /** return the coordinates {X, Y, Z} of parameter values T and S*/
    virtual int GetXYZofTS(double T, double S, double &X, double &Y,
                          double &Z) = 0;
    /** return the parameter values T and S of coordinates (X, Y, Z) */
    virtual int GetTSofXYZ(double X, double Y, double Z, double &T,
                           double &S) = 0;

    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S) = 0;

};

#endif
