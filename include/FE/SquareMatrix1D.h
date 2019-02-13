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
// @(#)SquareMatrix2D.h        
// 
// Class:       TSquareMatrix2D
//
// Purpose:     store a square matrix (ansatz = test space) in 1d
//
// Author:      Sashikumaar Ganesan
//
// History:     17.05.2007 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIX1D__
#define __SQUAREMATRIX1D__

#include <SquareMatrix.h>
#include <SquareStructure1D.h>

class TSquareMatrix1D : public TSquareMatrix
{
  protected:
    /** matrix strcuture */
    TSquareStructure1D *structure;

  public:
    /** generate the matrix */
    TSquareMatrix1D(TSquareStructure1D *squarestructure);

    /** destructor: free Entries array */
    ~TSquareMatrix1D();

    /** return FESpace */
    TFESpace1D *GetFESpace() const
    { return structure->GetFESpace(); }

    /** return used matrix structure */
    TSquareStructure1D *GetMatrixStructure() const
    { return structure; }

};

#endif
