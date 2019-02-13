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
// @(#)BoundComp2D.h        1.5 08/12/99
//
// Class:       TBoundComp2D
// Purpose:     Components of boundary faces
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BOUNDCOMP2D__
#define __BOUNDCOMP2D__

#include <BoundComp.h>

/** Components of boundary faces */
class TBoundComp2D : public TBoundComp
{
  public:
    // Constructor
    TBoundComp2D(int id);

    // Methods
    /** return the coordinates {X,Y} of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y) = 0;
    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T) = 0;

    /** get number of initial vertices on a Comp2Donent */
    virtual int GetN_InitVerts() = 0;
    virtual int GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges) = 0;

};

#endif
