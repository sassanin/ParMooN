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
// @(#)ADICell.h        4.1 07.11.09
// 
// Class:       TADICell
// Purpose:     general super class for all ADICells
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (07.11.09)
//
// History:     start of implementation 07.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __ADICELL__
#define __ADICELL__

#include <BaseCell.h>
#include <Collection.h>

/** general super class for all ADICells, special spaces are
    implemented in subclasses */
class TADICell
{
  protected:
// =======================================================================
// // information of cell and the dimension of the internal domain
// =======================================================================
    /** cells for which the info of internal direction is constructed*/
    TBaseCell *Cell;

    /** all internal domains use the same coll */
    TCollection *Collection_Internal;

    /** number of quadrature points in this cell, where the internal domains have to be constructed*/
    int N_QuadPts;

    /** number of degrees of freedom  (same for all QuadPt) */
    int N_V;

    /** solution vector  for each QuadPt in internal direction [N_QuadPts][N_V]*/
    double *sol;

    /** solution vector in internal direction for each QuadPt [N_V][N_QuadPts]*/
    double *solT;

    /** rhs vector  for each QuadPt in internal direction [N_QuadPts][N_V]*/
    double *rhs;

    /** rhs vector in internal direction for each QuadPt [N_V][N_QuadPts] */
    double *rhsT;

  public:
    /** constructor */
    TADICell(TBaseCell *cell, int N_quadPts);

    /** destrcutor */
    ~TADICell();

};

#endif
