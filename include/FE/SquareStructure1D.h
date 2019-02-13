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
// @(#)SquareStructure1D.C        
// 
// Class:       TSquareStructure1D
//
// Purpose:     build and store a structure for a square matrix in 1d
//
// Author:      Sashikumaar Ganesan
//
// History:     17.05.2007 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE1D__
#define __SQUARESTRUCTURE1D__

#include <FESpace1D.h>
#include <SquareStructure.h>

class TSquareStructure1D : public TSquareStructure
{
  protected:
    /** FE space */
    TFESpace1D *FESpace;

  public:
    /** dummy constructor, needed only derives classes */
    TSquareStructure1D();

    /** generate the matrix structure, only one space needed */
    TSquareStructure1D(TFESpace1D *space);

    /** destructor: free all used arrays */
    ~TSquareStructure1D();

    /** return FESpace */
    TFESpace1D *GetFESpace()
    { return FESpace; }

};

#endif

 
