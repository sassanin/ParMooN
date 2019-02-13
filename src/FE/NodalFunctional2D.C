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
// @(#)NodalFunctional2D.C        1.1 10/30/98
// 
// Class:       TNodalFunctional2D
// Purpose:     realize nodal functionals in 2D
//
// Author:      Gunar Matthies (02.09.98)
//
// History:     start of implementation 02.09.98 (Gunar Matthies)
//
// =======================================================================

#include <NodalFunctional2D.h>

/** constructor */
TNodalFunctional2D::TNodalFunctional2D(NodalFunctional2D id,
                int n_allfunctionals, 
                int n_edgefunctionals,
                int n_pointsall, int n_pointsedge,
                double *xi, double *eta, double *t,
                EvalAllNF *evalall,
                EvalJointNF *evaledge)
{
  ID = id;
  N_AllFunctionals = n_allfunctionals;
  N_EdgeFunctionals = n_edgefunctionals;
  N_PointsAll = n_pointsall;
  N_PointsEdge = n_pointsedge;
  Xi = xi;
  Eta = eta;
  T = t;
  EvalAll = evalall;
  EvalEdge = evaledge;
}
