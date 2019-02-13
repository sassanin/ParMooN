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
// @(#)Matrix.h        1.2 11/20/98
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

#ifndef __MATRIX__
#define __MATRIX__

#include <Structure.h>
#include <string>
#include <map>

class TMatrix
{
  protected:
    /** Sparse structure of the matrix */
    TStructure *structure;
    
    /** matrix elements in an array */
    double *Entries;
    
  public:
    /** generate the matrix, intialize entries with zeros */
    TMatrix(TStructure *structure);
    
    /** generate the matrix with given entries */
    TMatrix(TStructure *structure, double* Entries);
    
    /** @brief reset the structure, this may mean that the entries need to be 
     *         reallocated */
    void SetStructure(TStructure *structure);
    
    /** create a nRows*nCols zero matrix */
    TMatrix(int nRows, int nCols);

    /** destructor: free Entries array */
    ~TMatrix();

    /** reset matrix entries to zero */
    void Reset();

    /** return number of rows */
    int GetN_Rows() const
    { return structure->GetN_Rows(); }

    /** return number of columns */
    int GetN_Columns() const
    { return structure->GetN_Columns(); }
    
    /** return number of matrix entries */
    int GetN_Entries() const
    { return structure->GetN_Entries(); }

    /** return number of matrix entries for hanging node data */
    int GetHangingN_Entries() const
    { return structure->GetHangingN_Entries(); }

    /** return array KCol */
    int *GetKCol() const
    { return structure->GetKCol(); }

    /** return array HangingKCol */
    int *GetHangingKCol() const
    { return structure->GetHangingKCol(); }

    /** return array HangingRowPtr */
    int *GetHangingRowPtr() const
    { return structure->GetHangingRowPtr(); }

    /** return array RowPtr */
    int *GetRowPtr() const
    { return structure->GetRowPtr(); }

    /** return structure */
    TStructure* GetStructure() const
    { return structure; }
    
    /** return matrix entries */
    double *GetEntries() const
    { return Entries; }

    /** return the norm of the matrix. p is 
     * -2 for Frobenius norm
     * -1 for maximum absolute row sum
     *  0 for maximum entry
     *  1 for maximum absolute column sum
     *  2 for euclidean norm, 
     * 
     */
    double GetNorm(int p=-1) const;

    /** write matrix into file */
    int Write(const char *filename) const;
    
    /** @brief Print matrix into the shell */
    void Print(const char *name = "a") const;
    
    /** @brief print the full matrix, including all zeros
     * 
     * This is only meaningful for very small matrices.
     */
    void PrintFull(std::string name="", int fieldWidth=4) const;

    // add a value at selected entry
    void add(int i,int j, double val);
    // add values in row 'i' given by the map 'vals', multiplied by 'factor'
    // this should be faster than adding all values in 'vals' individually
    void add(int i, std::map<int,double> vals, double factor = 1.0);
    // add values 'vals[i][j]' to this matrix at the positions (i,j) for all
    // i,j defined in the map 'vals'
    void add(std::map<int, std::map<int,double> > vals, double factor = 1.0);
    // set a value at selected entry
    void set(int i, int j, double val);
    // get a value at selected entry
    const double& get(int i,int j) const;
    // get a value at selected entry (you may change that value (similar to set)
    double& get(int i,int j);
    
    // set the whole elements array
    void setEntries(double* entries) {
      this->Entries = entries;
    }
    
    /** @brief return a new TMatrix which is the transposed of this matrix 
     * 
     * If this is an object of a derived class (e.g. TMatrix2D, TSquareMatrix),
     * then the number of active degrees of freedom is not taken into account. 
     * The returned TMatrix is really the algebraic transposed matrix.
     * */
    TMatrix* GetTransposed() const;
    
    /** 
     * @brief replace several rows in the matrix with new entries.
     * 
     * Replace rows by new ones. This changes the sparsity pattern of the 
     * matrix. Therefore reallocation is necessary. The old arrays are no 
     * longer needed (for this matrix) and will be deleted if the flag 
     * 'deleteOldArrays' is set to true.  
     * 
     * If there are no rows to change, i.e. if entries.size()==0, nothing is 
     * done.
     * 
     * @param entries for every row a map of columns-to-entries map
     * @param deleteOldArrays remove old arrays in matrix and structure
     */
    void changeRows(std::map<int,std::map<int,double> > entries,
                    bool deleteOldArrays=false);
    
    
    /** @brief compute y = A*x   (Matrix-Vector-Multiplication)
     *
     * Note that 'y' is created here and it is up to the user to delete it. 
     */
    friend double* operator*(const TMatrix & A, const double* x);

    /** @brief add another matrix to this one
     * 
     * This is of course only possible if the corresponding structures are the
     * same. 
     */
    virtual TMatrix & operator+=(const TMatrix * A);
    /** @brief substract another matrix to this one
     * 
     * This is of course only possible if the corresponding structures are the
     * same. 
     */
    virtual TMatrix & operator-=(const TMatrix * A);
    /** @brief add another matrix to this one
     * 
     * This is of course only possible if the corresponding structures are the
     * same. This method exists only for convenience and uses the same method 
     * with a pointer to A instead of the reference.
     */
    virtual TMatrix & operator+=(const TMatrix & A) 
    { *this += &A; return *this; }
    /** @brief scale matrix by a factor */
    virtual TMatrix & operator*=(const double a);
    /** @brief copy entries from A to this 
     * 
     * This is of course only possible if the corresponding structures are the
     * same. 
     */
    TMatrix & operator=(const TMatrix& A);
    
    /** @brief compute y += a * A*x
     *
     * 'A' is this TMatrix and 'x', and 'y' are given vectors. The scalar 'a'
     * is a scaling factor
     * 
     * @param x array representing the vector which is multiplied by this matrix
     * @param y array representing the vector to which a * A*x is added
     * @param a scaling factor, default is 1.0
     */
    void multiply(const double * const x, double *y, double a = 1.0) const;
    
    /**
     * @brief compute matrix-matrix product C = a*A*B, 
     * 
     * 'A' is this matrix, 'a' is a scalar factor, 'B' is given. Then matrix
     * 'C' is created during this function and the user is responsible to 
     * delete C.
     * 
     * Note that this is rather slow.
     * 
     * @param B matrix to be multiplied (from right) to this matrix
     * @param a scaling factor, default is 1.0
     */ 
    TMatrix* multiply(const TMatrix * const B, double a = 1.0) const;
    
    /**
     * @brief scale a matrix using a vector
     * 
     * think of this as multipling this matrix with a diagonal matrix from the 
     * left (if the second argument is true). The parameter 'factor' are the 
     * entries of that diagonal matrix.
     * 
     * The i-th row (colum if from_left is false) is scaled by the i-th entry in
     * 'factor'.
     * 
     * The array 'factor' must be of size this->GetN_Rows() if 'from_left' is 
     * true or of size this->GetN_Columns() otherwise.
     * 
     * If all entries in 'factor' are the same, you can use operator*= as well.
     * 
     * @param factor array of scaling factors
     * @param from_left scale rows (true) or columns (false)
     */
    void scale(const double * const factor, bool from_left = true);
    
    /**
     * @brief remove all entries from sparsity structure where a zero is stored
     * 
     * This changes the sparsity structure of this matrix. Afterwards all stored
     * entries are nonzero. This can help if a lot of zeros are explicitly 
     * stored in the matrix. Afterwards matrix operations should be faster.
     * 
     * if tol is greater or equal to zero, all entries with magnitude smaller 
     * than tol are removed. If tol is smaller than zero, tol is set such that 
     * the ratio of the largest over the smallest entry in the resulting matrix
     * is smaller than 10^{15}.
     */
    void remove_zeros(double tol = 0.0);
    
    /** @brief get/set a specific matrix entry 
     * 
     * This will give an error if that entry is not in the sparsity structure
     */
    double & operator()(const int i, const int j);
    /** @brief get a specific matrix entry 
     * 
     * This will give an error if that entry is not in the sparsity structure
     */
    const double & operator()(const int i, const int j) const;
};

#endif
