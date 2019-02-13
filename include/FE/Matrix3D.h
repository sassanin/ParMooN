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
// @(#)Matrix3D.h        1.2 11/20/98
// 
// Class:       TMatrix3D
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#ifndef __MATRIX3D__
#define __MATRIX3D__

#include <Structure3D.h>
#include <Matrix.h>

class TMatrix3D : public TMatrix
{
  protected:
    /** matrix structure */
    TStructure3D *structure;

  public:
    /** generate the matrix */
    TMatrix3D(TStructure3D *structure);

    /** destructor: free Entries array */
    ~TMatrix3D();

    TStructure3D *GetStructure()
    { return structure; }
    
    /** @brief set all Dirichlet rows to zero. That means all rows where the 
     * test space has nonactive degrees of freedom. 
     */
    void resetNonActive();
};

#endif
