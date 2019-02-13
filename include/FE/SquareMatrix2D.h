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
// @(#)SquareMatrix2D.h        1.3 11/20/98
// 
// Class:       TSquareMatrix2D
//
// Purpose:     store a square matrix (ansatz = test space) in 2d
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIX2D__
#define __SQUAREMATRIX2D__

#include <SquareMatrix.h>
#include <SquareStructure2D.h>

class TSquareMatrix2D : public TSquareMatrix
{
  protected:
    /** matrix strcuture */
    TSquareStructure2D *structure;

  public:
    /** generate the matrix */
    TSquareMatrix2D(TSquareStructure2D *squarestructure);

//     /** generate an empty nxn matrix */
//     explicit TSquareMatrix2D(int n);
    
    /** fill empty matrix, 
     you can either call the constructor TSquareMatrix2D(TSquareStructure2D*);
     or use TSquareMatrix2D(); and then  void SetStructure(TSquareStructure2D*);
    */
    void SetStructure(TSquareStructure2D *squarestructure);
    
    /** destructor: free Entries array */
    ~TSquareMatrix2D();

    /** return FESpace */
    TFESpace2D *GetFESpace() const
    { return structure->GetFESpace(); }

    /** return used matrix structure */
    TSquareStructure2D *GetMatrixStructure() const
    { return structure; }
    TSquareStructure2D *GetStructure() const
    { return structure; }

    
    /** @brief scale matrix by scalar (only active entries) */
    TSquareMatrix2D& operator*=(double alpha);
    /** @brief add another matrix to this one (only active entries) */
    TSquareMatrix2D& operator+=(TSquareMatrix2D& rhsMat);
    /** @brief add another matrix to this one (only active entries) */
    TSquareMatrix2D& operator+=(TSquareMatrix2D* rhsMat)
    { *this += *rhsMat; return *this; }
    
    /** @brief copy matrix 'rhs' to this */
    TSquareMatrix2D& operator=(const TSquareMatrix2D& rhs);

    /** @brief add to matrices A and B
     * 
     * note: only active DOF are added
     * note: only works for matrices with the same sparsity pattern
    */
    friend TSquareMatrix2D& operator+(const TSquareMatrix2D & A, 
                                      const TSquareMatrix2D & B);
    
    /**  @brief C= A*alpha 
     * 
     * note: only active DOF are multiplied, others are just copied
     * note: the user is responsible to delete the newly created matrix C
     */
    friend TSquareMatrix2D& operator*(const TSquareMatrix2D & A,
                                      const double alpha);
    /**  @brief C= A*alpha 
     * 
     * Same as TSquareMatrix2D& operator+(const TSquareMatrix2D & A, 
     *                                    const TSquareMatrix2D & B);
     * 
     * note: only active DOF are multiplied, others are just copied
     * note: the user is responsible to delete the newly created matrix C
     */
    friend TSquareMatrix2D& operator*(const double alpha,
                                      const TSquareMatrix2D & A);

    /** @brief y= A*x
     * 
     * note: only active DOF are multiplied, others are just copied from x
     */
    friend double* operator*(const TSquareMatrix2D & A, const double* x);
};

#endif
