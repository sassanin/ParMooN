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
// @(#)FEVectFunct2D.h        
// 
// Class:       TFEVectFunct1D
// Purpose:     a function from a finite element space in 1D
//
// Author:      Sashikumaar Ganesan (26.09.09)
//
// History:     start of implementation 26.09.09 (Sashikumaar Ganesan)
//
//
// =======================================================================

#ifndef __FEVECTFUNCT1D__
#define __FEVECTFUNCT1D__

#include <FEFunction1D.h>

/** a function from a finite element space */
class TFEVectFunct1D : public TFEFunction1D
{
  protected:
    /** number of components */
    int N_Components;

  public:
    /** constructor with vector initialization */
    TFEVectFunct1D(TFESpace1D *fespace1D, char *name, char *description,
                  double *values, int length, int n_components);

    /** return number of components */
    int GetN_Components()
    { return N_Components; }

    /** return i-th component as FEFunction2D */
    TFEFunction1D *GetComponent(int i)
    {
      return new TFEFunction1D(FESpace1D, Name, Description,
                               Values+i*Length, Length);
    }

    /** convert current grid to vector-values FE function */
//     void GridToData();

    /** use current data for grid replacement */
//     void DataToGrid();



};

#endif
