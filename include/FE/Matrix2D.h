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
// @(#)Matrix2D.h        1.2 11/20/98
// 
// Class:       TMatrix2D
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#ifndef __MATRIX2D__
#define __MATRIX2D__

#include <Structure2D.h>
#include <Matrix.h>

class TMatrix2D : public TMatrix
{
  protected:
    /** matrix structure */
    TStructure2D *structure;

  public:
    /** generate the matrix */
    TMatrix2D(TStructure2D *structure);
    
    /** @brief generate empty matrix */
    TMatrix2D();
    
    /** @brief fill empty matrix, 
     * 
     * you can either call the constructor TMatrix2D(TStructure2D*);
     * or use TMatrix2D(); and then  void SetStructure(TStructure2D*);
     */
    void SetStructure(TStructure2D *structure);

    /** destructor: free Entries array */
    ~TMatrix2D();

    TStructure2D *GetStructure() const
    { return structure; }
    
    /** @brief scale all active rows */
    TMatrix2D& operator*=(double alpha);

    /** @brief set all Dirichlet rows to zero. That means all rows where the test 
     *         space has nonactive degrees of freedom. */
    void resetNonActive();

    /** @brief add two matrices A and B
     * 
     * note: only active DOF are added
     * note: only works for matrices with the same sparsity pattern 
     */
    friend TMatrix2D& operator+(const TMatrix2D & A, const TMatrix2D & B);
    /** @brief substract matrices A and B, i.e. C = A - B
     * 
     * note: only active DOF are substracted
     * note: only works for matrices with the same sparsity pattern
     */
    friend TMatrix2D& operator-(const TMatrix2D & A, const TMatrix2D & B);
    
    /** @brief C = A*alpha 
     * 
     * C will consist of the active entries of A scaled by alpha. C will have
     * zero entries in nonactive rows.
     */
    friend TMatrix2D& operator*(const TMatrix2D & A, const double alpha);
    /** @brief C = alpha*A 
     * 
     * Same as TMatrix2D& operator*(const TMatrix2D & A, const double alpha);
     */
    friend TMatrix2D& operator*(const double alpha, const TMatrix2D & A);

    /** @brief y = A*x
     * 
     * note: only active DOF are multiplied, others are just copied from x
     * note: the user has to delete y
     */
    friend double* operator*(const TMatrix2D & A, const double* x);
};

void AllocateMatricesNSE_2D(int mg_level,
			    TFESpace2D *velocity_space, 
			    TFESpace2D *pressure_space,
			    TSquareStructure2D *&sqstructureA, 
			    TSquareStructure2D *&sqstructureC, 
			    TStructure2D *&structureB, 
			    TStructure2D *&structureBT,
			    TSquareMatrix2D *&sqmatrixA,
			    TSquareMatrix2D *&sqmatrixA11,
			    TSquareMatrix2D *&sqmatrixA12,
			    TSquareMatrix2D *&sqmatrixA21,
			    TSquareMatrix2D *&sqmatrixA22,
			    TSquareMatrix2D *&sqmatrixC,
			    TMatrix2D *&matrixB1,
			    TMatrix2D *&matrixB2,
			    TMatrix2D *&matrixB1T,
			    TMatrix2D *&matrixB2T,
			    TSquareMatrix2D **MatricesA,
			    TSquareMatrix2D **MatricesA11,
			    TSquareMatrix2D **MatricesA12,
			    TSquareMatrix2D **MatricesA21,
			    TSquareMatrix2D **MatricesA22,
			    TSquareMatrix2D **MatricesC,
			    TMatrix2D **MatricesB1,			    
			    TMatrix2D **MatricesB2,			    
			    TMatrix2D **MatricesB1T,			    
			    TMatrix2D **MatricesB2T);


#endif // __MATRIX2D__
