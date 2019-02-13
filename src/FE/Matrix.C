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
// @(#)Matrix.C        1.2 11/20/98
// 
// Class:       TMatrix
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Matrix.h>
#include <string.h>
#include <Constants.h>
#include <LinAlg.h>
#include <stdlib.h>

#include <MooNMD_Io.h>

TMatrix::TMatrix(TStructure *structure)
{
  this->structure = structure;
  Entries = new double[structure->GetN_Entries()];
  memset(Entries, 0, structure->GetN_Entries()*SizeOfDouble);
}

TMatrix::TMatrix(TStructure *structure, double* Entries)
{
  this->structure = structure;
  this->Entries = Entries;
}

TMatrix::TMatrix(int nRows, int nCols)
{
  structure = new TStructure(nRows,nCols);
  Entries = NULL;
}

void TMatrix::SetStructure(TStructure *structure)
{
  this->structure = structure;

  if (Entries) delete Entries;
  Entries = new double[structure->GetN_Entries()];
  memset(Entries, 0, structure->GetN_Entries()*SizeOfDouble);
}

void TMatrix::Reset()
{
  memset(Entries, 0, structure->GetN_Entries()*SizeOfDouble);
}

TMatrix::~TMatrix()
{
  delete [] Entries;
}


void TMatrix::Print(const char *name) const
{
  int* RowPtr = structure->GetRowPtr();
  int* KCol = structure->GetKCol();
  int begin, end, pos=0;
  
  for (int i=0; i<structure->GetN_Rows(); ++i)
  {
    begin = RowPtr[i];
    end   = RowPtr[i+1];
    
    for (int j=begin; j<end; ++j)
    {
      OutPut(name<< "(" << i+1 << "," << KCol[pos]+1 << ") = " << Entries[pos] 
             << ";\n");
      ++pos;
    }
  }
}


void TMatrix::PrintFull(std::string name, int fieldWidth) const
{
  int* rowPtr = structure->GetRowPtr();
  int* KCol = structure->GetKCol();

  cout << endl << name << " = " << endl;
  for (int curRow = 0; curRow < structure->GetN_Rows(); curRow++) 
  {
    int rowEnd = rowPtr[curRow+1];
    int posKCol = rowPtr[curRow];
    for (int curCol = 0; curCol < structure->GetN_Columns(); curCol++) 
    {
      if(curCol == KCol[posKCol] && posKCol < rowEnd) 
      {
        cout << setw(fieldWidth) << Entries[posKCol] << ", ";
        posKCol++;
      }
      else 
      {
        cout << setw(fieldWidth) << 0.0 << ", ";
      }
    }
    cout << endl;
  }
  cout << endl;
}

// add val to a matrix element
// return an error if the entry is not in the sparse structure
void TMatrix::add(int i,int j, double val)
{
  if(val != 0.0)
    this->get(i, j) += val;
}

void TMatrix::add(int i, std::map<int,double> vals, double factor)
{
  if(i < 0 || i > this->GetN_Rows())
  {
    ErrMsg("This matrix does not have a row " << i << 
          ".\nThe dimension of this matrix is " << this->GetN_Rows() <<
          " x " << this->GetN_Columns());
    exit(0);
  }
  int* RowPtr = structure->GetRowPtr();
  int* KCol = structure->GetKCol();
  std::map<int,double>::iterator it = vals.begin();
  for (int m=RowPtr[i];m < RowPtr[i+1] && it != vals.end(); m++) 
  {
    if (KCol[m] == it->first) 
    {
      Entries[m] += factor*it->second;
      ++it;
    }
  }
  if(it != vals.end())
  {
    Error("Error in TMatrix::add. There are entries in 'vals' which are not "
          << "in the sparse structure. row " << i << ", column " << it->first
          << ".\nExit\n");
    exit(0);
  }
}

void TMatrix::add(std::map<int, std::map<int,double> > vals, double factor)
{
  // add every row to the matrix
  std::map<int,std::map<int,double> >::iterator it;
  for(it = vals.begin(); it != vals.end(); ++it)
    add(it->first, it->second, factor);
}
  
  

// set val of a matrix element
// return an error if the entry is not in the sparse structure
void TMatrix::set(int i,int j, double val)
{
  this->Entries[this->structure->index_of_entry(i, j)] = val;
}

// get val of a matrix element
// return an error if the entry is not in the sparse structure
const double& TMatrix::get(int i,int j) const
{
  // index of the entry (i,j) within the vector Entries
  int index = this->structure->index_of_entry(i, j);
  if(index >= 0 )
    return this->Entries[index];
  ErrMsg("could not find the entry (" << i << "," << j 
         << ") in the sparsity structure");
  exit(1);
}

double& TMatrix::get(int i,int j)
{
  // index of the entry (i,j) within the vector Entries
  int index = this->structure->index_of_entry(i, j);
  if(index >= 0 )
    return this->Entries[index];
  ErrMsg("could not find the entry (" << i << "," << j 
         << ") in the sparsity structure");
  exit(1);
}

double & TMatrix::operator()(const int i, const int j)
{
  return this->get(i, j);
}
//was ist der unterschied ausser const?
const double & TMatrix::operator()(const int i, const int j) const
{
  return this->get(i, j);
}

double TMatrix::GetNorm(int p) const
{
  double result = 0.0;
  switch(p)
  {
    case -2:
      result = Dnorm(this->GetN_Entries(), Entries);
      break;
    case -1:
    {
      int *rows = this->GetRowPtr();
      for(int row=0; row<this->GetN_Rows(); row++)
      {
        double row_sum = 0.0;
//#pragma omp parallel for
        for(int i=rows[row]; i<rows[row+1]; i++)
        {
          row_sum += fabs(Entries[i]);
        }
        if(row_sum>result && row_sum!=1.0)
          result = row_sum;
      }
      break;
    }
    case 0:
      for(int i=0; i<this->GetN_Entries(); i++)
      {
        double a = fabs(Entries[i]);
        if(a > result)
          result = a;
      }
      break;
    case 1:
      Error("spectral norm of a matrix not yet implemented!\n");
      exit(0);
      break;
    case 2:
      Error(" maximum absolute column sum norm of a matrix not yet "
            << "implemented!\n");
      exit(0);
      break;
    default:
      Error("undefined norm of a matrix!\n");
      exit(0);
      break;
  }
  return result;
}


double* operator*(const TMatrix & A,const double* x)
{
  double *AEntries = A.GetEntries();
  int *ARowPtr = A.GetRowPtr();
  int *AColIndex = A.GetKCol();

  int nrows = A.GetN_Rows();
  
  double *y = new double[nrows];
  double value;
  int index;

  for(int i=0;i<nrows;i++) 
  {
    value = 0;
//#pragma omp parallel for
    for (int j=ARowPtr[i]; j<ARowPtr[i+1]; j++) 
    {

      index = AColIndex[j];
      value += AEntries[j] * x[index];
    }
    y[i] = value;
  }
  return y;
}

TMatrix & TMatrix::operator+=(const TMatrix* A)
{
  if(this->structure != A->GetStructure() // compare pointers
     && (*(this->structure)) != (*(A->GetStructure()))) // compare objects
  {
    ErrMsg("TMatrix::operator+= : the two matrices do not match.");
    exit(1);
  }
  
  int n_entries = this->GetN_Entries();
  double *AEntries = A->GetEntries();
  for(int i = 0; i < n_entries; ++i)
  {
    this->Entries[i] += AEntries[i];
  }
  return *this;
}



TMatrix & TMatrix::operator-=(const TMatrix* A)
{
  if(this->structure != A->GetStructure() // compare pointers
     && (*(this->structure)) != (*(A->GetStructure()))) // compare objects
  {
    ErrMsg("TMatrix::operator-= : the two matrices do not match.");
    exit(1);
  }
  
  int n_entries = this->GetN_Entries();
  double *AEntries = A->GetEntries();
  for(int i = 0; i < n_entries; ++i)
  {
    this->Entries[i] -= AEntries[i];
  }
  return *this;
}


TMatrix & TMatrix::operator=(const TMatrix& A)
{
  // compare structures (first pointers, then structures themselves)
  if(this->structure != A.GetStructure() // compare pointers
     && ((*(this->structure)) != (*(A.GetStructure())))) // compare objects
  {
    OutPut("WARNING: TMatrix& operator= Matrices have different structures.\n"
           << "         The structure in this matrix (on left hand side) is "
           << "         reset to use the same structure as the right hand "
           << "side\n         Are you sure this is what you want?\n");
    this->SetStructure(A.GetStructure());
  }
  
  // copy matrix entries.
  int n_entries = this->GetN_Entries();
  double *AEntries = A.GetEntries();
  for(int i = 0; i < n_entries; i++)
    this->Entries[i] = AEntries[i];
  return *this;
}

void TMatrix::multiply(const double * const x, double *y, double a) const
{
  if(a == 0.0)
    return;
  
  int *rowPtr = GetRowPtr();
  int *colIndex = GetKCol();

  int nrows = GetN_Rows();
  
  double value;
  int i, j, end;
  for(i = 0; i < nrows; i++) 
  {
    value = 0;
    end = rowPtr[i+1];
    for (j = rowPtr[i]; j < end; j++) 
    {
      value += Entries[j] * x[colIndex[j]];
    }
    y[i] += a * value;
  }
}

TMatrix* TMatrix::multiply(const TMatrix * const B, double a) const
{
  const int n_A_rows = this->GetN_Rows();   // = n_C_rows
  const int n_A_cols = this->GetN_Columns();
  const int n_B_rows = B->GetN_Rows();
  
  if(n_A_cols != n_B_rows)
  {
    ErrMsg("dimention mismatch during matrix-matrix multiplication");
    exit(1);
  }
  const int * const a_rows = this->GetRowPtr();
  const int * const a_cols = this->GetKCol();
  
  TStructure * strucB = B->GetStructure();
  const double * const b_entries = B->GetEntries();
  
  TStructure * struc_c = get_product_structure(this->GetStructure(), strucB);
  int * c_rows = struc_c->GetRowPtr();
  int * c_cols = struc_c->GetKCol();
  double * c_entries = new double[struc_c->GetN_Entries()];
  memset(c_entries, 0.0, struc_c->GetN_Entries() * SizeOfDouble);
  
  // fill the entries
  // loop over all rows in C
//#pragma omp parallel for
  for(int row = 0; row < n_A_rows; row++)
  {
    
    // loop over all entries in this row in C
    for(int col = c_rows[row]; col < c_rows[row + 1]; col++)
    {
      // multiply 'this row of A' x 'this column of B'
      // loop over all entries in this row in A
      for(int i = a_rows[row]; i < a_rows[row+1]; i++)
      {
        int ib = strucB->index_of_entry(a_cols[i], c_cols[col]);
        if(ib != -1)
        {
          c_entries[col] += Entries[i] * b_entries[ib];
        }
      }
    }
  }
  return new TMatrix(struc_c, c_entries);
}


void TMatrix::remove_zeros(double tol)
{
  if(tol < 0)
    tol = this->GetNorm(0) * 1e-15; // largest (in magnitude) entry
  // we want to call this->changeRows(new_entries, true)
  std::map<int,std::map<int,double> > new_entries;
  int *rows = structure->GetRowPtr();
  int *cols = structure->GetKCol();
  int n_rows = structure->GetN_Rows();
  int n_removed = 0; // number of entries to be removed
  // loop over all rows of this matrix
  for(int row=0; row<n_rows; row++)
  {
    new_entries[row]; // empty row
    for(int col = rows[row]; col < rows[row + 1]; col++)
    {
      if(fabs(Entries[col]) > tol)
        (new_entries[row])[cols[col]] = Entries[col];
    }
    n_removed += rows[row + 1] - rows[row] - new_entries[row].size();
  }
  if(n_removed != 0)
  {
    OutPut("TMatrix::remove_zeros: tol " << tol << "\tn_removed " << n_removed
            << "\tratio " << (double)n_removed/(rows[n_rows]) << endl);
    this->changeRows(new_entries, true);
  }
  else
    OutPut("TMatrix::remove_zeros: no removable entries\n");
}


void TMatrix::scale(const double * const factor, bool from_left)
{
  int *rowPtr = GetRowPtr();
  int *colIndex = GetKCol();
  
  if(from_left)
  {
    for(int i = 0, nrows = GetN_Rows(); i < nrows; i++) 
    {
      int end = rowPtr[i+1];
      // scale entire row with the same factor
      for(int j = rowPtr[i]; j < end; j++) 
      {
        Entries[j] *= factor[i];
      }
    }
  }
  else
  {
    for(int i = 0, nrows = GetN_Rows(); i < nrows; i++) 
    {
      int end = rowPtr[i+1];
      for(int j = rowPtr[i]; j < end; j++) 
      {
        // scale columnwise
        Entries[j] *= factor[colIndex[j]];
      }
    }
  }
}

TMatrix & TMatrix::operator*=(const double a)
{
  Dscal(this->GetN_Entries(), a, Entries);
  return *this;
}


void TMatrix::changeRows(std::map<int,std::map<int,double> > entries,
                         bool deleteOldArrays)
{
  if(entries.size() == 0)
    return; // nothing needs to be done
  
  int *oldRows = structure->GetRowPtr();
  int *oldCols = structure->GetKCol();
  
  // find out how many entries there are after all changes are applied, i.e. 
  // how many entries are deleted/created
  int offset = 0;
  for(std::map<int,std::map<int,double> >::iterator it=entries.begin(); 
       it!=entries.end(); ++it)
  {
    int row = it->first;
    offset -= oldRows[row+1]-oldRows[row];// number of entries in old structure
    offset += (it->second).size();        // number of entries in new structure
  }
  
  int n_rows = structure->GetN_Rows();// new number of rows = old number of rows
  int n_entries = structure->GetN_Entries() + offset; // new number of entries
  int *columns = new int[n_entries];  // new pointer to columns
  int *rows = new int[n_rows+1];      // new row pointer
  rows[0] = 0;
  
  // create new array to store the entries
  double * new_entries = new double[n_entries];
  
  // fill the arrays 'rows', 'columns' and 'new_entries'
  for(int row=0; row<n_rows; row++)
  {
    std::map<int,std::map<int,double> >::iterator it = entries.find(row);
    if(it == entries.end())
    {
      // this row stays unchanged
      // copy pointer to columns in this row
      memcpy(columns+rows[row],oldCols+oldRows[row],
             (oldRows[row+1]-oldRows[row])*SizeOfInt);
      // update row pointer
      rows[row+1] = rows[row] + oldRows[row+1] - oldRows[row];
      // copy entries
      memcpy(new_entries+rows[row],Entries+oldRows[row],
             (oldRows[row+1]-oldRows[row])*SizeOfDouble);
    }
    else
    {
      // this row will be replaced
      std::map<int,double> newRow = it->second;
      // loop over all new entries in this row
      int columnIndex=0;
      for(std::map<int,double>::iterator it2 = newRow.begin(); 
          it2 != newRow.end(); ++it2)
      {
        int colInd = it2->first;
        double entry = it2->second;
        columns[columnIndex+rows[row]] = colInd;
        new_entries[columnIndex+rows[row]] = entry;
        columnIndex++;
      }
      rows[row+1] = rows[row] + newRow.size();
      //if(newRow.size() != columnIndex)
      //  OutPut("ERROR: wrong number of columns in this row "<< newRow.size()
      //      << "\t" << columnIndex << "\t" << row << endl);
    }
  }
  
  // change Structure of this matrix
  structure->setN_Entries(n_entries);
  structure->setKCol(columns);
  structure->setRowPtr(rows);
  
  if(deleteOldArrays)
  {
    delete [] oldRows;
    delete [] oldCols;
    delete [] Entries; // maybe we should always delete this
  }
  Entries = new_entries;
}
