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
// @(#)QuadFormula.h        1.3 05/04/99
//
// Class:    TQuadFormula
//
// Purpose:  base class for quadrature formulas
// Author:   Gunar Matthies
//
// History:  29.08.1997 start implementation
//           07.04.1999 add Accuracy
// 
// =======================================================================

#ifndef __QUAD_FORMULA__
#define __QUAD_FORMULA__

#include <Constants.h>
#include <MooNMD_Io.h>

/** base class for quadrature formulas */
class TQuadFormula
{
  protected:
    /** number of quadrature points */
    int N_QuadPoints;
    /** weights for the formula */
    double *Weights;

    /** accuracy of this formula */
    int Accuracy;

  protected:
    /** constructor */
    TQuadFormula();

  public:
    /** return number of quadrature points */
    int GetN_QuadPoints()
    { return N_QuadPoints; }

    /** return weights of the formula */
    double *GetWeights()
    { return Weights; }
    /** return coordinates of the formula */
    virtual double *GetCoords(int i);

    /** print information on this formula */
    friend std::ostream & operator << (std::ostream &s, TQuadFormula *qf);
};

#endif
