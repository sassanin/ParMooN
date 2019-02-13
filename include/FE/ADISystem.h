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
// @(#)ADISystem.h        4.1 13.11.09
// 
// Class:       TADISYSTEM
// Purpose:     general super class for all ADI/operator-splitting System
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (13.11.09)
//
// History:     start of implementation 13.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __ADISYSTEM__
#define __ADISYSTEM__

#include <BaseCell.h>
#include <Collection.h>

/** general super class for ADI/operator splitting scheme, special spaces are
    implemented in subclasses */
class TADISystem
{
  protected:
// =======================================================================
// // information of cell and the dimension of the internal domain
// =======================================================================
    /** all internal domains use the same coll */
    TCollection *Collection_Intl;

    /** number of degrees of freedom  (same for all QuadPt) */
    int N_V;

    /** solution vector: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *sol;

    /** solution vector: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *oldsol;

    /** rhs vector  for each QuadPt in internal direction [N_V]*/
    double *rhs;

    /** working rhs: if we have, we no need to allocate & deallocate for each Intl poin in internal direction */
    double *B;

    /** tmp array, if we have, we no need to allocate & deallocate for each Intl point*/
    double *defect;

    /** No. M mat value  */
    int N_MMatValues;

    /** DoubleFunct2D to calculate growth rate*/
    DoubleFunctND *GetGrowthAndNuc;

    /** store l values */
    double *IntlPosL;

    /** No. A mat value  */
//     int N_AMatValues;

    /** store the SUPG S-mat value for next time step rhs calculation, for other than Euler schems, time-dep. growth */
//     double *SMatValues;

    /** No.  SUPG S-mat values*/
//     int N_SMatValues;

    /** store the SUPG K-mat value for next time step rhs calculation, for other than Euler schems, time-dep. growth */
//     double *KMatValues;

    /** No. SUPG K-mat values*/
//     int N_KMatValues;

  public:
    /** constructor */
    TADISystem(double *Sol, double *OldSol, double *b, double *Defect, DoubleFunctND *growthfunct);

    /** destrcutor */
    ~TADISystem();

};

#endif
