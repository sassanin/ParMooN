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
// @(#)QuadFormulaTria.h        1.3 12/08/99
//
// Class:      TQuadFormulaTria
// Superclass: TQuadFormula2D
//
// Purpose:    quadrature formula for a 1D integral
// Author:     Gunar Matthies
//
// History:    29.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_TRIA__
#define __QUAD_FORMULA_TRIA__

#include <QuadFormula2D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 2D integral */
class TQuadFormulaTria : public TQuadFormula2D
{
  public:
    /** constructor */
    TQuadFormulaTria();
    /** constructor */
    TQuadFormulaTria(int n_points, double* weights, double* xi, 
                      double* eta, int acc);

    /** BaryCenter */
    void BaryCenter();
    /** MidPoint */
    void MidPoint();
    /** SevenPoint */
    void SevenPoint();
    /** high order formula, nearly a Gauss3 formula */ 
    void Gauss3();
    /** Vertex */
    void Vertex();
    /** formula of degree 8 */
    void Degree8();
      /** formula of degree 9 */  
    void Degree9();
    /** formula of degree 11 */
    void Degree11();
    /** formula of degree 19 */
    void Degree19();
    /** Gauss-like formula, composed on three sub-triangles */
    void CompGauss3();
    /** Gauss-like formula, composed on four sub-triangles */
    void CompGauss4();
    /** Gauss formula of degree 8 */
    void Gauss_Degree8();
};

#endif
