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
// @(#)StructureNSE2D.h        1.6 09/14/99
// 
// Class:       TStructureStokes2D
//
// Purpose:     build and store matrix for (Navier-)Stokes equations in 2d
//
// Author:      Gunar Matthies
//
// History:     10.10.2003 start implementation
//
// =======================================================================

#ifndef __STRUCTURENSE2D__
#define __STRUCTURENSE2D__

#include <SquareStructure2D.h>

class TStructureNSE2D : public TSquareStructure2D
{
  protected:
    /** begin array for jb array */
    int *BeginJb;

    /** local number of special edge dof */
    int *jb;

    /** number of dof on each joint */
    int N_DOFperJoint;

    /** coefficient array due to special edge dof */
    double *Alpha;

    void GenerateAlpha();

  public:
    /** generate the matrix structure, only one space needed */
    TStructureNSE2D(TFESpace2D *space);

    /** destructor: free all used arrays */
    ~TStructureNSE2D();

    int *GetBeginJb()
    { return BeginJb; }

    int *GetJb()
    { return jb; }

    int GetN_DOFperJoint()
    { return N_DOFperJoint; }

    double *GetAlpha()
    { return Alpha; }

};

#endif
