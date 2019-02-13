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
// @(#)SquareMatrixNSE3D.h        1.3 11/20/98
// 
// Class:       TSquareMatrixNSE3D
//
// Purpose:     store a square matrix (ansatz = test space) in 3D
//
// Author:      Gunar Matthies
//
// History:     29.07.05 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIXNSE3D__
#define __SQUAREMATRIXNSE3D__

#include <SquareMatrix3D.h>
#include <SquareStructure3D.h>

class TSquareMatrixNSE3D : public TSquareMatrix3D
{
  protected:
    /** begin array for jb array */
    int *BeginJb;

    /** local number of special edge dof */
    int *jb;

    /** number of dof on each joint */
    int N_DOFperJoint;

    /** coefficient array due to special edge dof */
    double *Alpha;

    /** begin array for bubble correction */
    int *BeginC;

    /** bubble correction coefficients */
    double *C;

    /** begin array for nonconstant pressure */
    int *BeginP;

    /** nonconstant pressure parts */
    double *P;

  public:
    /** generate the matrix */
    TSquareMatrixNSE3D(TSquareStructure3D *squarestructure);

    /** destructor: free Entries array */
    ~TSquareMatrixNSE3D();

    int *GetBeginJb()
    { return BeginJb; }

    int *GetJb()
    { return jb; }

    int GetN_DOFperJoint()
    { return N_DOFperJoint; }

    double *GetAlpha()
    { return Alpha; }

    int *GetBeginC()
    { return BeginC; }

    void SetBeginC(int *beginc)
    { BeginC = beginc; }

    double *GetC()
    { return C; }

    void SetC(double *c)
    { C = c; }

    int *GetBeginP()
    { return BeginP; }

    void SetBeginP(int *beginp)
    { BeginP = beginp; }

    double *GetP()
    { return P; }

    void SetP(double *p)
    { P = p; }
};

#endif
