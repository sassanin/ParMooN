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
// @(#)Structure.C        1.10 11/24/99
// 
// Class:       TStrcuture
//
// Purpose:     build and store a matrix structure
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation (Gunar Matthies)
//
//              04.08.1998 start reimplementation (Gunar Matthies)
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================

#include <Structure.h>
#include <Constants.h>
#include <MooNMD_Io.h>
#include <cstring>
#include <vector>
#include <stdlib.h>


/** generate the matrix structure, both spaces are 2D */
TStructure::TStructure() 
 : N_Rows(0), N_Columns(0), N_Entries(0), HangingN_Entries(0), KCol(NULL),
   HangingKCol(NULL), RowPtr(NULL), HangingRowPtr(NULL)
{
}

TStructure::TStructure(int n, int N_entries, int *col_ptr, int *row_ptr)
 : N_Rows(n), N_Columns(n), N_Entries(N_entries), HangingN_Entries(0), 
   KCol(col_ptr), HangingKCol(NULL), RowPtr(row_ptr), HangingRowPtr(NULL)
{
}

/** generate the matrix structure, all arrays are already defined */
TStructure::TStructure(int nRows, int nCols, int N_entries, int *col_ptr, 
                       int *row_ptr)
 : N_Rows(nRows), N_Columns(nCols), N_Entries(N_entries), HangingN_Entries(0), 
   KCol(col_ptr), HangingKCol(NULL), RowPtr(row_ptr), HangingRowPtr(NULL)
{
}

TStructure::TStructure(int nRows, int nCols)
 : N_Rows(nRows), N_Columns(nCols), N_Entries(0), HangingN_Entries(0), 
   KCol(NULL), HangingKCol(NULL), RowPtr(new int[nRows+1]), HangingRowPtr(NULL)
{
  memset(RowPtr, 0, (N_Rows+1)*SizeOfInt);
}



/** sort one row [BeginPtr, AfterEndPtr) */
void TStructure::SortRow(int *BeginPtr, int *AfterEndPtr)
{
  int *IPtr, *JPtr, T;

  for(IPtr=BeginPtr;IPtr<AfterEndPtr;IPtr++)
  {
    for(JPtr=IPtr+1;JPtr<AfterEndPtr;JPtr++)
    {
      if( *IPtr > *JPtr )
      {
        T = *IPtr;
        *IPtr = *JPtr;
        *JPtr = T;
      }
    } // endfor JPtr
  } // endfor IPtr
}

/** sort numbers within each row */
void TStructure::Sort()
{
  int end, begin;

  end = 0;
  for(int i=0; i<N_Rows; i++)
  {
    begin = end;
    end = RowPtr[i+1];
    SortRow(KCol+begin, KCol+end);
  } // endfor i
}

/** destructor */
TStructure::~TStructure()
{
  delete [] KCol;
  delete [] HangingKCol;
  delete [] RowPtr;
  delete [] HangingRowPtr;
}

int TStructure::index_of_entry(const int i, const int j) const
{
  if(i < 0 || i >= this->GetN_Rows())
  {
    ErrMsg("row index is too large or too small");
    exit(1);
  }
  if(j < 0 || j > this->GetN_Columns())
  {
    ErrMsg("column index is too large or too small");
    exit(1);
  }
  
  for (int m=RowPtr[i];m < RowPtr[i+1]; m++) 
  {
    if (KCol[m]== j) 
    {
      // index found in sparsity pattern
      return m;
    }
  }
  // index not in the sparsity pattern
  return -1;
}


/** return a new structure for a transposed matrix */
TStructure* TStructure::GetTransposed()
{
  if(HangingN_Entries!=0)
  {
    Error("TStructure::GetTransposed(): Hanging entries not supported! Exit\n");
    exit(0);
  }
  // variables for transposed structure:
  int nRowsT = N_Columns;
  int nColsT = N_Rows;
  // number of entries does not change
  int *rowsT = new int[nRowsT+1];  memset(rowsT, 0, (nRowsT+1)*SizeOfInt);
  int *colsT = new int[N_Entries]; memset(colsT, 0, N_Entries *SizeOfInt);
  
  int *ColB_count = new int[N_Columns]; 
  memset(ColB_count, 0, N_Columns*SizeOfInt);
  
  // count number of entries per column (will be number of entries in each row)
  for(int i=0;i<RowPtr[N_Rows];i++)
    rowsT[KCol[i]]++;
  // change to increasing numbering as in RowPtr
  for(int i=0,k=0;i<=nRowsT;i++)
  {
    int j = rowsT[i];
    rowsT[i] = k;
    k += j;
  }
  
  // fill 'colsT'
  // loop over (non-transposed) rows
  for(int i=0; i<N_Rows; i++)
  {
    // loop over all entries in this row
    for(int j=RowPtr[i]; j<RowPtr[i+1]; j++)
    {
      int col = KCol[j]; // (non-transposed) column = transposed row
      int l  = rowsT[col];
      int offset = ColB_count[col];
      colsT[l+offset] = i;
      ColB_count[col]++;
    }
  }
  delete []ColB_count;
  
  TStructure* structureT = new TStructure(nRowsT, nColsT, N_Entries, colsT, rowsT); 
  return structureT;
}

TStructure* get_product_structure(TStructure const * const strucA,
                                  TStructure const * const strucB)
{
  const int n_A_rows = strucA->GetN_Rows();   // = n_C_rows
  const int n_A_cols = strucA->GetN_Columns();
  const int n_B_rows = strucB->GetN_Rows();
  const int n_B_cols = strucB->GetN_Columns();   // = n_C_cols
  
  if(n_A_cols != n_B_rows)
  {
    ErrMsg("dimension mismatch during matrix-matrix multiplication");
    exit(1);
  }
  const int * const a_rows = strucA->GetRowPtr();
  const int * const a_cols = strucA->GetKCol();
  
  // everything needed to call the constructor of TStructure later on:
  int n_c_entries = 0; // number of entries in product structure C
  int * c_rows = new int[n_A_rows + 1]; // row pointer
  memset(c_rows,0.0, (n_A_rows + 1)*SizeOfInt);
  int * c_cols; // columns pointer
  std::vector<std::vector<int> > dofs(n_A_rows);
  
  // loop over all rows of C
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all columns of C
    for(int col = 0; col < n_B_cols; col++)
    {
      // check whether 'this row of A' x 'this column of B' would give an entry
      int n_a_entries_in_row = a_rows[row+1] - a_rows[row];
      // loop over all entries in this row in A
      for(int i = 0; i < n_a_entries_in_row; i++)
      {
        if(strucB->index_of_entry(a_cols[i+a_rows[row]], col) != -1)
        {
          dofs[row].push_back(col);
          break;
        }
      }
    }
    n_c_entries += dofs[row].size();
    c_rows[row+1] = n_c_entries; // c_rows[0] has been set to 0 already
  }
  
  // now fill the array c_cols
  c_cols = new int[n_c_entries];
  // loop over all rows of C
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all columns of C
    int n_entries_in_this_row = c_rows[row+1] - c_rows[row];
    for(int col = 0; col < n_entries_in_this_row; col++)
    {
      c_cols[c_rows[row] + col] = dofs[row].at(col);
    }
  }  
  return new TStructure(n_A_rows, n_B_cols, n_c_entries, c_cols, c_rows);
}


bool operator==(const TStructure &lhs, const TStructure &rhs)
{
  return lhs.N_Rows == rhs.N_Rows 
      && lhs.N_Columns == rhs.N_Columns
      && lhs.N_Entries == rhs.N_Entries
      && lhs.HangingN_Entries == rhs.HangingN_Entries;
}

bool operator!=(const TStructure &lhs, const TStructure &rhs)
{
  return !(rhs == lhs);
}
