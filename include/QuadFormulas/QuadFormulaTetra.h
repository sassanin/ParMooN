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
// @(#)QuadFormulaTetra.h        1.2 05/04/99
//
// Class:      TQuadFormulaTetra
// Superclass: TQuadFormula3D
//
// Purpose:    quadrature formula for a 3D integral on a tetrahedron
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_TETRA__
#define __QUAD_FORMULA_TETRA__

#include <QuadFormula3D.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 3D integral on a tetrahedron */
class TQuadFormulaTetra : public TQuadFormula3D
{
  public:
    /** constructor */
    TQuadFormulaTetra();
    /** constructor */
    TQuadFormulaTetra(int n_points, double* weights, 
      double* xi, double* eta, double* zeta, int acc);

    /** BaryCenter */
    void BaryCenter();

    /** Vertex */
    void Vertex();

    /** P2 exact formula */
    void P2Exact();

    /** P4 exact formula */
    void P4Exact();

    /** P5 exact formula */
    void P5Exact();

    /** P8 exact formula */
    void P8Exact();
};

#endif
