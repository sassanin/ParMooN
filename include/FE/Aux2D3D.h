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
// %W% %G%
//
// Class:       TAux2D3D
// Purpose:     calculate parameters from 3d fefuntion for the use
//              in 2d
//
// Author:      Gunar Matthies (15.05.01)
//
// History:     start of implementation 15.05.01 (Gunar Matthies)
// =======================================================================

#ifndef __AUX2D3D__
#define __AUX2D3D__

#include <FEFunction3D.h>

/** calculate parameters from 3d fefuntion for the use in 2d */
class TAux2D3D
{
  protected:
    /** collection of 3d cells */
    TCollection *Coll;

    /** 3d cell number */
    int *CellNumbers;

    /** joint numbers for 3d cells */
    int *JointNumbers;

    /** 3d finite element function */
    TFEFunction3D *FEFunct;

    /** 3d finite element space */
    TFESpace3D *FESpace;

    /** BeginIndex */
    int *BeginIndex;

    /** GlobalNumbers */
    int *GlobalNumbers;

    /** value vector of FEFunct */
    double *Values;

    /** shift in parameter list */
    int Shift;

  public:
    /** constructor */
    TAux2D3D(int *cellnumbers, int *jointnumbers, TFEFunction3D *fefunct,
             int shift);

    /** calculate gradient for local coordinates (xi,eta) on face
        JointNumbers[num] of cell CellNumbers[num] */
    void GetGradient(int num, int N_Points, double *xi, double *eta,
                     double **Param);
};

#endif
