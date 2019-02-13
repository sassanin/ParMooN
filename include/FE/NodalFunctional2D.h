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
// @(#)NodalFunctional2D.h        1.1 10/30/98
// 
// Class:       TNodalFunctional2D
// Purpose:     realize nodal functionals in 2D
//
// Author:      Gunar Matthies (02.09.98)
//
// History:     start of implementation 02.09.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __NODALFUNCTIONAL2D__
#define __NODALFUNCTIONAL2D__

#include <Enumerations.h>
#include <Constants.h>

/**  realize nodal functionals in 2D */
class TNodalFunctional2D
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
    EvalAllNF *EvalAll;

    /** Number of points needed for edge nodal functionals */
    int N_PointsEdge;

    /** values at edge points with following edge parameters in [-1,1]
        are needed to evaluate functional on edge */
    double *T;

    /** routine for evaluating the edge functionals */
    EvalJointNF *EvalEdge;

    /** ID for this set of nodal functionals */
    NodalFunctional2D ID;

  public:
    /** constructor */
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       EvalAllNF *evalall,
                       EvalJointNF *evaledge);

    /** return information for points for all functionals */
    void GetPointsForAll(int &n_points, double* &xi, double* &eta)
    { n_points = N_PointsAll; xi = Xi; eta = Eta; }

    /** return information for points for edge functionals */
    void GetPointsForEdge(int &n_points, double* &t)
    { n_points = N_PointsEdge; t = T; }

    /** return values for all nodal functionals */
    void GetAllFunctionals(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                           double *Functionals)
    { EvalAll(Coll, Cell, PointValues, Functionals); }

    /** return values for edge nodal functional */
    void GetEdgeFunctionals(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                            double *Functionals)
    { EvalEdge(Coll, Cell, Joint, PointValues, Functionals); }

    /** return ID for this set */
    NodalFunctional2D GetID()
    { return ID; }

};

#endif
