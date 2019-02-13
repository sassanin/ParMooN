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
   
#ifndef __CG__
#define __CG__

#include <ItMethod.h>

/** iteration method */
class TCg : public TItMethod
{
  protected:

    /** arrays for cg depending on number of dof */
    double *r;
    double *z;
    double *p;
    double *Ap;
    //TODO eingefuegt fuer Polak-Ribiere
    double *r_last;
    //TODO eingefuegt fuer Polak-Ribiere
    //TODO eingefuegt fer flexible Paper
    double **d;
    //TODO eingefuegt fer flexible Paper

    /** matrices for flexible gmres depending on restart */

  public:
    /** constructor */
    TCg(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
      int n_aux, int N_Unknowns, int scalar);

    /** destructor */
    ~TCg();

    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol,
      double *rhs);
};
#endif
