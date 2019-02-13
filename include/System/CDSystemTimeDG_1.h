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
   
/** ************************************************************************ 
*
* @class     TCDSystemTimeDG_1
* @brief     stores the information of a 2D scalar dG(1) in time discretization
* @author    Sashikumaar Ganesan, 
* @date      24.12.15
* @History    
 ************************************************************************  */


#ifndef __CDSYSTEMTIMEDG_1__
#define __CDSYSTEMTIMEDG_1__

#include <SquareMatrix2D.h>
#include <CDSystemTimeDG.h>

/**class for 2D scalar system  dG(1) in time discretization */
 
class TCDSystemTimeDG_1  : public TCDSystemTimeDG
{
  protected:
    /** dG type */
    int Type;    
    
  public:
    /** constructor */
     TCDSystemTimeDG_1(TSquareMatrix2D *mat_m, TSquareMatrix2D *mat_A);

    /** destrcutor */
    ~TCDSystemTimeDG_1();
    
    /** assemble the system matrix */
    virtual void AssembleSysMat(double *Mu_old, double *Rhs);
    
    /** assemble the system matrix */   
    virtual void AssembleALESysMat_Qp1(double *Mu_old, double *Rhs);
    
    
    /** solve dG system and return the Sol at end t^n system matrix */    
    virtual void SoveTimedG(double *Sol); 
    
};

#endif
