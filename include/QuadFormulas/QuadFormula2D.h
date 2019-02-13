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
// @(#)QuadFormula2D.h        1.2 05/04/99
//
// Class:      TQuadFormula2D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 2D integral
// Author:     Gunar Matthies
//
// History:    29.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_2D__
#define __QUAD_FORMULA_2D__

#include <Enumerations.h>
#include <QuadFormula.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 2D integral */
class TQuadFormula2D : public TQuadFormula
{
  protected:
    /** first coordinate in [0,1]x[0,1] for the formula */
    double *Xi;
    /** second coordinate in [0,1]x[0,1] for the formula */
    double *Eta;

  protected:
    /** This is a private method for initializing the data structure */
    void InitObject(int n, double* w, double* xi, double* eta, int acc);

  public:
    /** constructor */
    TQuadFormula2D();
    /** constructor */
    TQuadFormula2D(int n_points, double* weights, double* xi, 
                    double* eta, int acc);

    /** return coordinates of the formula */
    virtual double *GetCoords(int i);
    /** return all data of the quadrature formula */
    void GetFormulaData(int &n_points, double* &weights, 
                        double* &xi, double* &eta);

// #ifdef __2D__
    /** return a quadrature formula which can be used for 
        all given elements */
    static void FindQuadFormula2D(FE2D *UsedElements,
        QuadFormula1D &qf1, QuadFormula2D &qf2);

    /** return a quadrature formula which can be used for
        the given elements */
    static void FindQF_2D(FE2D CurrentElement,
        QuadFormula1D &qf1, QuadFormula2D &qf2);

    /** find a quadrature formula for all given elements */
    static void FindLocalQuadFormula2D(int N_LocalUsedElements,
        FE2D *LocalUsedElements,
        QuadFormula1D &qf1, QuadFormula2D &qf2);
// #endif // __2D__

    /** print all information of this formula */
    friend std::ostream & operator << (std::ostream &s, TQuadFormula2D *qf);
};

#endif
