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
// @(#)SquareStructure.h        1.6 09/14/99
// 
// Class:       TSquareStructure
//
// Purpose:     build and store a structure for a square matrix
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE__
#define __SQUARESTRUCTURE__

#include <Structure.h>

class TSquareStructure : public TStructure
{
  protected:
    /** number of active rows */
    int ActiveBound;

    /** ordering of the column entries */
    /** 0 - no special ordering */
    /** 1 - increasing ordering (like used for umfpack) */
    /** 2 - diagonal entry first, then increasing ordering */
    int ColOrder;

    /** sort an integer array */
    void IntSort(int *BeginPtr, int *AfterEndPtr);

  public:
    /** generate the matrix structure, only one space needed */
    TSquareStructure();

    /** destructor: free all used arrays */
    ~TSquareStructure();

    /** generate the matrix structure, all arrays are already defined */
    TSquareStructure(int n, int N_entries, int *col_ptr,
      int *row_ptr);
    
    /** Generates an empty n*n Structure for a Zero-Matrix */
    explicit TSquareStructure(int n);

    /** return ActiveBound */
    int GetActiveBound() const
    { return ActiveBound; }

    /** return ordering of columns */
    int GetColOrder() const
    { return ColOrder;}

    /** sort column numbers in each row, increasing indices */
    void Sort();

    void SetColOrder(int n)  
    {ColOrder = n;}    
    
    /** sort column numbers: diag is first element, other numbers are
    increasing */
    void SortDiagFirst();
    
    /** @brief Comparision Operators */
    friend bool operator==(const TSquareStructure &lhs, 
                           const TSquareStructure &rhs);
    friend bool operator!=(const TSquareStructure &lhs,
                           const TSquareStructure &rhs);
};

#endif
