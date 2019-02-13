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
// @(#)NodalFunctional1D.h        1.1 11/20/98
// 
// Class:       TNodalFunctional1D
// Purpose:     realize nodal functionals in 1D
//
// Author:      Gunar Matthies (08.10.98)
//
// History:     start of implementation 08.10.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __NODALFUNCTIONAL1D__
#define __NODALFUNCTIONAL1D__

#include <Enumerations.h>
#include <Constants.h>

/**  realize nodal functionals in 1D */
class TNodalFunctional1D
{
  protected:
    /** number of all functionals */
    int N_AllFunctionals;

    /** number of functionals on one edge */
    int N_EdgeFunctionals;

    /** number of points needed for all nodal functionals */
    int N_PointsAll;

    /** values at points with following xi-coordinates are needed
        to evaluate all nodal functionals */
    double *Xi;

    /** values at points with following eta-coordinates are needed
        to evaluate all nodal functionals */
    double *Eta;

    /** routine for evaluating all functionals */
    DoubleFunctVect *EvalAll;

    /** Number of points needed for edge nodal functionals */
    int N_PointsEdge;

    /** values at edge points with following edge parameters in [-1,1]
        are needed to evaluate functional on edge */
    double *T;

    /** routine for evaluating the edge functionals */
    DoubleFunctVect *EvalEdge;

    /** ID for this set of nodal functionals */
    NodalFunctional1D ID;

  public:
    /** constructor */
    TNodalFunctional1D(NodalFunctional1D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);

    /** return information for points for all functionals */
    void GetPointsForAll(int &n_points, double* &xi, double* &eta)
    { n_points = N_PointsAll; xi = Xi; eta = Eta; }

    /** return information for points for edge functionals */
    void GetPointsForEdge(int &n_points, double* &t)
    { n_points = N_PointsEdge; t = T; }

    /** return values for all nodal functionals */
    void GetAllFunctionals(double *PointValues, double *Functionals)
    { EvalAll(PointValues, Functionals); }

    /** return values for edge nodal functional */
    void GetEdgeFunctionals(double *PointValues, double *Functionals)
    { EvalEdge(PointValues, Functionals); }

    /** return ID for this set */
    NodalFunctional1D GetID()
    { return ID; }

};

#endif
