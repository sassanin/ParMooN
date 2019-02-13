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
// @(#)QuadFormula1D.h        1.2 05/04/99
//
// Class:      TQuadFormula1D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 1D integral
// Author:     Gunar Matthies
//
// History:    29.08.1997 start implementation
// 
// =======================================================================

#ifndef __QUAD_FORMULA_1D__
#define __QUAD_FORMULA_1D__

#include <QuadFormula.h>
#include <MooNMD_Io.h>

/** quadrature formula for a 1D integral */
class TQuadFormula1D : public TQuadFormula
{
  protected:
    /** coordinates in [0,1] for the formula */
    double *Xi;

  protected:
    /** This is a private method for initializing the data structure */
    void InitObject(int n, double* w, double* xi, int acc);

  public:
    /** constructor */
    TQuadFormula1D();
    /** constructor */
    TQuadFormula1D(int n_points, double* weights, double* xi, int acc);

    /** return coordinates of the formula */
    virtual double *GetCoords(int i);
    /** return all data of the quadrature formula */
    void GetFormulaData(int &n_points, double* &weights, double* &xi);

    /** return accuracy of this formula */
    int GetAccuracy()
    { return Accuracy; }

    /** 1-Point-Gauss */
    void Gauss1();
    /** 2-Points-Gauss */
    void Gauss2();
    /** 3-Points-Gauss */
    void Gauss3();
    /** 4-Points-Gauss */
    void Gauss4();
    /** 5-Points-Gauss */
    void Gauss5();
    /** 6-Points-Gauss */
    void Gauss6();
    /** 7-Points-Gauss */
    void Gauss7();
    /** 8-Points-Gauss */
    void Gauss8();
    /** 9-Points-Gauss */
    void Gauss9();
    /** 10-Points-Gauss */
    void Gauss10();
    /** 11-Points-Gauss */
    void Gauss11();
    /** 12-Points-Gauss */
    void Gauss12();
    /** 2-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
    void Gauss2W1();
    /** 4-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
    void Gauss4W1();
    /** 6-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
    void Gauss6W1();
    /** 8-Points-Gauss with ln(1/x) as the weighting function on the interval [0;1] */
    void Gauss8W1();
    /** 16-Points-Gauss for strongly singular integral of order O(1/x) on the interval [-1;1] */
    void Gauss16W2();
    
    /** print information on this formula */
    friend std::ostream & operator << (std::ostream &s, TQuadFormula1D *qf);
    
    /** destructor */    
    ~TQuadFormula1D();
};

#endif
