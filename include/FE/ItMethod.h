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
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TItMethod
// Purpose:     defines an iteration method
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================

#ifndef __ITMETHOD__
#define __ITMETHOD__

#include <Constants.h>

/** iteration method */
class TItMethod
{
  protected:
    /** routine for matrix vector product */
    MatVecProc *matvec;

    /** routine for defect */
    DefectProc *matvecdefect;

    /** preconditioner */
    TItMethod *prec;

    /** number of degrees of freedom  */
    int N_DOF;

    /** id of the system type  */
    int system_id;

    /** maximal number of preconditioner iterations */
    int prec_maxit;

    /** absolute tolerance for stopping */
    double res_norm_min;

    /** relative tolerance for stopping */
    double red_factor;

    /** limit divergence */
    double div_factor;

    /** maximal number of iterations */
    int maxit;

    /** minimal number of iterations */
    int minit;

    /** restart */
    int restart;

    /** array for defect */
    double *defect;

    /** auxiliary arrays */
    double **AuxArray; 

    /** number of auxiliary arrays */
    int N_Aux;

  public:
    /** constructor */
    TItMethod(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
              int n_aux, int N_Unknowns);

    /** destructor */
    virtual ~TItMethod();

    /** iterate */
    /** input **A, **B */
    /** output (solution) *sol */
    /** input (rhs) *rhs */
    virtual int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                        double *rhs) = 0;

    /** return system id */
    int GetSystemId() const
    { return system_id; };

    /** return maximal number of preconditioner iterations */
    int GetPrecMaxit() const
    { return prec_maxit; };
    
    /** set maximal number of preconditioner iterations */
    void SetPrecMaxit(int Prec_Maxit)
    { prec_maxit = Prec_Maxit; }

    /** return absolute tolerance for stopping */
    double GetResNormMin() const
    { return res_norm_min; };
    
    /** set absolute tolerance for stopping */
    void SetResNormMin(double Res_Norm_Min)
   { res_norm_min = Res_Norm_Min; }

    /** return relative tolerance for stopping */
    double GetRedFactor() const
    { return red_factor; };
    
    /** set relative tolerance for stopping */
    void SetRedFactor(double Red_Factor)
    { red_factor = Red_Factor; }

    /** return tolerance for divergence */
    double GetDivFactor() const
    { return div_factor; };
    
    /** set tolerance for divergence */
    void SetDivFactor(double Div_Factor)
    { div_factor = Div_Factor; }

    /** return maximal number of iterations */
    int GetMaxit() const
    { return maxit; };
    
    /** set maximal number of iterations */
    void SetMaxit(int Maxit)
    { maxit = Maxit; }

    /** return minimal number of iterations */
    int GetMinit() const
    { return minit; };

    /** set minimal number of iterations */
    void SetMinit(int Minit)
    { minit = Minit; }

    /** return restart */
    int GetRestart() const
    { return restart; };
    
    /** set restart */
    void SetRestart(int Restart)
    { restart = Restart; };
};

#endif
