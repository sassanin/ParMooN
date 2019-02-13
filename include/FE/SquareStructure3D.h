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
// @(#)SquareStructure3D.h        1.6 09/14/99
// 
// Class:       TSquareStructure3D
//
// Purpose:     build and store a structure for a square matrix in 3D
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE3D__
#define __SQUARESTRUCTURE3D__

#include <FESpace3D.h>
#include <SquareStructure.h>

class TSquareStructure3D : public TSquareStructure
{
  protected:
    /** FE space */
    TFESpace3D *FESpace;

  public:
    /** dummy constructor, needed only for derived classes */
    TSquareStructure3D();
  
    /** generate the matrix structure, only one space needed */
    TSquareStructure3D(TFESpace3D *space);

    /** generate the matrix structure, all arrays are already defined */
    TSquareStructure3D(int n, int N_entries, int *col_ptr,
      int *row_ptr);
    
    /** Generates an empty n*n Structure for a Zero-Matrix */
    explicit TSquareStructure3D(int n);

    /** destructor: free all used arrays */
    ~TSquareStructure3D();

    /** return FESpace */
    TFESpace3D *GetFESpace()
    { return FESpace; }

};

#endif
