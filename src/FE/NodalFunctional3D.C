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
// %W% %G% 
// 
// Class:       TNodalFunctional3D
// Purpose:     realize nodal functionals in 3D
//
// Author:      Gunar Matthies (30.05.2000)
//
// History:     start of implementation 30.05.00 (Gunar Matthies)
//
// =======================================================================

#include <NodalFunctional3D.h>

/** constructor */
TNodalFunctional3D::TNodalFunctional3D(NodalFunctional3D id,
                   int n_allfunctionals, int *n_facefunctionals,
                   int n_pointsall, int *n_pointsface,
                   double *xi, double *eta, double *zeta,
                   double **xiarray, double **etaarray,
                   double **zetaarray,
                   double *t, double *s,
                   EvalAllNF *evalall,
                   EvalJointNF *evalface)
{
  ID = id;

  N_AllFunctionals = n_allfunctionals;
  N_FaceFunctionals = n_facefunctionals;
  N_PointsAll = n_pointsall;
  N_PointsFace = n_pointsface;

  Xi = xi;
  Eta = eta;
  Zeta = zeta;

  XiArray = xiarray;
  EtaArray = etaarray;
  ZetaArray = zetaarray;

  T = t;
  S = s;

  EvalAll = evalall;
  EvalFace = evalface;
}

/** return information for points for all functionals */
void TNodalFunctional3D::GetPointsForAll(int &n_points, double* &xi,
                        double* &eta, double* &zeta)
{ 
  n_points = N_PointsAll;
  xi = Xi; 
  eta = Eta;
  zeta = Zeta;
}

/** return information for points for face functionals 
    on joint j */
void TNodalFunctional3D::GetPointsForFace(int j, int &n_points,
                        double* &xi, double* &eta, double* &zeta)
{ 
  n_points = N_PointsFace[j];
  xi   = XiArray[j];
  eta  = EtaArray[j];
  zeta = ZetaArray[j];
}

/** return information for points for face functionals */
void TNodalFunctional3D::GetPointsForFace(int &n_points, double* &t,
                                          double* &s)
{
  n_points = N_PointsFace[0],
  t = T;
  s = S;
}
