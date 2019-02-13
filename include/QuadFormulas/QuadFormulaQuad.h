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
// @(#)QuadFormulaQuad.h        1.3 05/04/99
//
// Class:      TQuadFormulaQuad
// Superclass: TQuadFormula2D
//
// Purpose:    quadrature formula for a 2D integral on the unit sqare
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_QUAD__
#define __QUAD_FORMULA_QUAD__

#include <QuadFormula2D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 2D integral */
class TQuadFormulaQuad : public TQuadFormula2D
{
  public:
    /** constructor */
    TQuadFormulaQuad();
    /** constructor */
    TQuadFormulaQuad(int n_points, double* weights, double* xi, 
                      double* eta, int acc);

    /** Gauss9x9 */
    void Gauss9();

    /** Gauss8x8 */
    void Gauss8();

    /** Gauss7x7 */
    void Gauss7();

    /** Gauss6x6 */
    void Gauss6();

    /** Gauss5x5 */
    void Gauss5();

    /** Gauss4x4 */
    void Gauss4();

    /** Gauss3x3 */
    void Gauss3();

    /** Gauss2x2 */
    void Gauss2();

    /** Vertex */
    void Vertex();

    /** Simpson */
    void Simpson();

};

#endif
