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
// @(#)QuadFormula3D.h        1.2 05/04/99
//
// Class:      TQuadFormula3D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 3D integral
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_3D__
#define __QUAD_FORMULA_3D__

#include <Enumerations.h>
#include <QuadFormula.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 3D integral */
class TQuadFormula3D : public TQuadFormula
{
  protected:
    /** first coordinate for the formula */
    double *Xi;
    /** second coordinate for the formula */
    double *Eta;
    /** third coordinate for the formula */
    double *Zeta;

  protected:
    /** This is a private method for initializing the data structure */
    void InitObject(int n, double* w, double* xi, 
                    double* eta, double* zeta, int acc);

  public:
    /** constructor */
    TQuadFormula3D();
    /** constructor */
    TQuadFormula3D(int n_points, double* weights, 
      double* xi, double* eta, double* zeta, int acc);

    /** return coordinates of the formula */
    virtual double *GetCoords(int i);
    /** return all data of the quadrature formula */
    void GetFormulaData(int &n_points, double* &weights, 
                        double* &xi, double* &eta, double* &zeta);

#ifdef __3D__
    /** find a quadrature formula for all given elements */
    static void FindLocalQuadFormula3D(int N_LocalUsedElements,
        FE3D *LocalUsedElements,
        QuadFormula2D &qf1, QuadFormula3D &qf2);
#endif // __3D__

    /** print all information of this formula */
    friend std::ostream & operator << (std::ostream &s, TQuadFormula3D *qf);
};

#endif
