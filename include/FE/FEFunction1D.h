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
// @(#)FEFunction1D.h
// 
// Class:       TFEFunction1D
// Purpose:     a function from a finite element space in 1D
// 
// Author:      Sashikumaar Ganesan (17.05.2007)
//
// History:     start of implementation 17.05.2007
//
// =======================================================================

#ifndef __FEFUNCTION1D__
#define __FEFUNCTION1D__

#include <AllClasses.h>
#include <FESpace1D.h>
#include <Constants.h>

/** a function from a finite element space */
class TFEFunction1D
{
  protected:
    /** name of the function */
    char *Name;

    /** some more words describing the function */
    char *Description;

    /** space to which this function belongs to */
    TFESpace1D *FESpace1D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /** constructor with vector initialization */
    TFEFunction1D(TFESpace1D *fespace1D, char *name, char *description,
                  double *values, int length);

    /** destructor */
    ~TFEFunction1D();

    /** return name */
    char *GetName()
    { return Name; }

    /** return description */
    char *GetDescription()
    { return Description; }

    /** return fe space */
    TFESpace1D *GetFESpace1D()
    { return FESpace1D; }

    /** return length */
    int GetLength()
    { return Length; }

    /** return vector of data */
    double *GetValues()
    { return Values; }


    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct2D *Exact);

    /** calculate the interpolation of an exact function */
    void Interpolate(int ConstCoord, double x, DoubleFunct2D *Exact);

   /** calculate for 1d function for given 2D , PBE*/
    void Interpolate(double x, double y, DoubleFunct3D *Exact);

    void InterpolateNodalPts(int N_Coord, double *Coords, DoubleFunctND *Exact, double *val);

    /** convert current grid to vector-values FE function */
    void GridToData();

};

#endif
