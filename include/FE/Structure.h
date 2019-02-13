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
// @(#)Structure.h        1.3 09/14/99
// 
// Class:       TStructure
//
// Purpose:     build and store a matrix structure
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE__
#define __STRUCTURE__

class TStructure
{
  protected:
    /** number of rows */
    int N_Rows;

    /** number columns */
    int N_Columns;

    /** number of matrix entries */
    int N_Entries;

    /** number of matrix entries in hanging nodes part */
    int HangingN_Entries;

    /** in which column is the current entry */
    int *KCol;

    /** in which column is the current entry (hanging nodes part */
    int *HangingKCol;

    /** index in KCol where each row starts */
    int *RowPtr;

    /** index in HangingKCol where each row starts */
    int *HangingRowPtr;

  public:
    /** generate the matrix structure, both space with 2D collection */
    TStructure();
    
    /** generate a (square) matrix structure, all arrays are already defined */
    TStructure(int n, int N_entries, int *col_ptr, int *row_ptr);

    /** generate the matrix structure, all arrays are already defined */
    TStructure(int nRows, int nCols, int N_entries, int *col_ptr, int *row_ptr);
    
    /** Generates an empty nRows*nCols Structure for a Zero-Matrix */
    TStructure(int nRows, int nCols);

    /** destructor: free all used arrays */
    ~TStructure();

    /** return number of rows */
    int GetN_Rows() const
    { return N_Rows; }

    /** return number of columns */
    int GetN_Columns() const
    { return N_Columns; }
    
    /** return number of matrix entries */
    int GetN_Entries() const
    { return N_Entries; }

    /** return number of matrix entries (hanging nodes part) */
    int GetHangingN_Entries() const
    { return HangingN_Entries; }

    /** return array KCol */
    int *GetKCol() const
    { return KCol; }

    /** return array HangingKCol */
    int *GetHangingKCol() const
    { return HangingKCol; }

    /** return array RowPtr */
    int *GetRowPtr() const
    { return RowPtr; }

    /** return array HangingRowPtr */
    int *GetHangingRowPtr() const
    { return HangingRowPtr; }
    
    /** @brief set member variables. Careful, this can produce inconsistencies! */
    void setN_Rows(int n) { N_Rows = n; }
    void setN_Columns(int n) { N_Columns = n; }
    void setN_Entries(int n) { N_Entries = n; }
    void setKCol(int * p) { KCol = p; }
    void setRowPtr(int * p) { RowPtr = p; }

    /** sort one row */
    void SortRow(int *BeginPtr, int *AfterEndPtr);
    
    /** sort rows */
    void Sort();
    
    
    /**
     * @brief find the index of a given entry
     * 
     * If the (i,j)-th entry is not in the sparsity pattern, -1 is returned. 
     * This is how this function can be used to check whether an entry is in the
     * sparsity pattern.
     * 
     * @param i row of entry to check
     * @param j column of entry to check
     */ 
    int index_of_entry(const int i, const int j) const;
    
    /** return a new structure for a transposed matrix 
     * If this is an object of a derived class (e.g. TStructure2D, 
     * TSquareStructure), then the number of active degrees of freedom is not 
     * taken into account. The returned TMatrix is really the algebraic 
     * transposed matrix.
     * */
    TStructure* GetTransposed();
    
    /**
     * @brief return a structure for the matrix-matrix-product A*B
     * 
     * if A and B are matrices with structures 'strucA' and 'strucB', this
     * function computes a structure for the product C = A*B
     * 
     * @param strucA structure of left factor
     * @param strucB structure of right factor
     */
    friend TStructure* get_product_structure(TStructure const * const strucA,
                                             TStructure const * const strucB);
    
    /** @brief Comparision Operator 
     * 
     * It is not explicitly checked if the arrays are the same, only the 
     * integers are compared.
     */
    friend bool operator==(const TStructure &lhs, const TStructure &rhs);
    friend bool operator!=(const TStructure &lhs, const TStructure &rhs);
};

#endif
