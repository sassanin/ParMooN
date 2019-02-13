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
// Class:       TFixedPointIte
// Purpose:     defines the fixed point iteration
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================
#ifndef __FIXEDPOINT__
#define __FIXEDPOINT__

#include <ItMethod.h>
#include <NSE_MultiGrid.h>
#ifdef _MPI   
   #ifdef __3D__
    #include <ParFECommunicator3D.h>
   #else
    #include <ParFECommunicator2D.h>
   #endif
#endif 

/** iteration method */
class TFixedPointIte : public TItMethod
{
  protected: 
#ifdef _MPI   
   #ifdef __3D__
      TParFECommunicator3D *ParComm;
   #else
      TParFECommunicator2D *ParComm;
   #endif
#endif   
  
  public:
   
    /** constructor */
    TFixedPointIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
                   int n_aux, int N_Unknowns, int scalar);

#ifdef _MPI   
    TFixedPointIte(MatVecProc *MatVec, 
                               DefectProc *Defect, 
                               TItMethod *Prec,
                               int n_aux, int n_dof,
                               int scalar, 
  #ifdef  __3D__			       
			       TParFECommunicator3D *parcomm
  #else			       
                               TParFECommunicator2D *parcomm
  #endif
                                );
#endif     
    
    /** destructor */
    virtual ~TFixedPointIte();
    
    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                double *rhs);    
};
#endif
