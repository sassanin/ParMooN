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
// @(#)BdNonUniformSpline.h        1.2 07/16/99
//
// Class:       TBdNonUniformSpline
// Superclass:  TBoundComp
// Purpose:     spline function as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BDNONUNIFORMSPLINE__
#define __BDNONUNIFORMSPLINE__

#include <BoundComp2D.h>

/** spline function as a component of a boundary part */
/*    determined by 10 parameters;
    every subspline is presented as X(t), Y(t), t=0..1,
    X(t) = Params[0]*phi1(t) + Params[2]*phi2(t) +
           Params[4]*phi3(t) + Params[6]*phi4(t)
    Y(t) = Params[1]*phi1(t) + Params[3]*phi2(t) +
           Params[5]*phi3(t) + Params[7]*phi4(t); 
    phi1(0) = 1, phi2(0) = 0, phi3(0) = 0, phi4(0) = 0;
    phi1(1) = 0, phi2(1) = 1, phi3(1) = 0, phi4(1) = 0;
    d_phi1(0) = 0, d_phi2(0) = 0, d_phi3(0) = 1, d_phi4(0) = 0;
    d_phi1(1) = 0, d_phi2(1) = 0, d_phi3(1) = 0, d_phi4(1) = 1; 
 Params[8] = T[i], T=0..1 through the whole spline;
 Params[9] is not used */

class TBdNonUniformSpline : public TBoundComp2D
{
  protected:
    /** number of subsplines */
    int N_Splines;
    /** array for all parameters */
    double *Params;
    /** array for the values of parameter T at the i-th boundary points */
    double *Param9;

  public:
    // Constuctor
    /** constructor initializes the parameter array */
    TBdNonUniformSpline (int id, int N_Splines);
    // Destructor
    ~TBdNonUniformSpline () {delete Params; delete Param9;};

    // Methods
    /** set all parameters */
    void SetParams (double *params);
    /** get number of splines */
    int GetN_Splines();
    /** return new inner points for a good boundary approximation */
    //double *GetBoundPoints(QuadFormula *formula, int &N_NewPt);

    /** return the coordinates of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y);

    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T);

    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat);

    /** get number of initial vertices on this component */
    virtual int GetN_InitVerts()
    { return 4; }
    virtual int GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges)
    { return -1; }
    
    /** return the X-coordinate of parameter value T from [0;1] on the ISpline-th subspline */
    double GetLocalXofT(int ISpline, double LocalT);
    
    /** return the Y-coordinate of parameter value T from [0;1] on the ISpline-th subspline */
    double GetLocalYofT(int ISpline, double LocalT);
    
    /** Generate all parameters according to the info about the pisition of boundary points */
    /** and the values of derivatives at the ends of boundary; the parameter is the cumulative */
    /** length of the piecewise line that joints the boundary points */ 
    void GenerateParams1(double *x, double *y, double dx0, double dy0, double dx1, double dy1);

    /** Solution of a tridiagonal matrix */
    /** with diagonal elements 'a', subdiagonal elemets 'b' and updiagonal elements 'c' */
    /** a[i], i=0..N_Splines-1, b[i], i=1..N_Splines-1, c[i], i=0..N_Splines-2 */
    void Solver_3diag(double *a, double *b, double *c, double *rhs, double *sol);

    /** for BEM */
    /** (x-x_middle)/t */
    double AdXofT(int ISpline, double LocalT);
    /** (y-y_middle)/t */
    double AdYofT(int ISpline, double LocalT);
    /** x'(t) */
    double GetLocalDXofT(int ISpline, double LocalT);
    /** y'(t) */
    double GetLocalDYofT(int ISpline, double LocalT);
};

#endif
