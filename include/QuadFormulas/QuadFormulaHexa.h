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
// @(#)QuadFormulaHexa.h        1.2 05/04/99
//
// Class:      TQuadFormulaHexa
// Superclass: TQuadFormula3D
//
// Purpose:    quadrature formula for a 3D integral on a hexahedron
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_HEXA__
#define __QUAD_FORMULA_HEXA__

#include <QuadFormula3D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 3D integral on a hexahedron */
class TQuadFormulaHexa : public TQuadFormula3D
{
  public:
    /** constructor */
    TQuadFormulaHexa();
    /** constructor */
    TQuadFormulaHexa(int n_points, double* weights, 
      double* xi, double* eta, double* zeta, int acc);

    /** Vertex */
    void Vertex();
    /** Gauss2x2x2 */
    void Gauss2();

    /** Gauss 3x3x3 */
    void Gauss3();

    /** Gauss 4x4x4 */
    void Gauss4();

    /** Gauss 5x5x5 */
    void Gauss5();

    /** Gauss 6x6x6 */
    void Gauss6();

    /** Gauss 7x7x7 */
    void Gauss7();

    /** Gauss 8x8x8 */
    void Gauss8();

    /** Gauss 9x9x9 */
    void Gauss9();

    /** quad rule with vertices and origin */
    void VerticesAndOrigin();

    /** quad rule with vertices and origin, 15 quad points */
    void VerticesAndOrigin15();
 
    /** quad rule with vertices and origin, 57 quad points */
    void VerticesAndOrigin57();

    /** quad rule for polynomials of degree 7 with 38 quad points */
    void Degree7_Points38();
};

#endif
