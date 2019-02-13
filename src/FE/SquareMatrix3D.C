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
// @(#)SquareMatrix3D.C        1.2 11/20/98
// 
// Class:       TSquareMatrix3D
//
// Purpose:     store a square matrix (ansatz = test space) in 3d
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#include <SquareMatrix3D.h>
#include <string.h>

TSquareMatrix3D::TSquareMatrix3D(TSquareStructure3D *squarestructure)
  : TSquareMatrix(squarestructure), structure(squarestructure)
{
}

// not clear! what is the initial value of structure while calling contructor/
// commented by Sashi 
// TSquareMatrix3D::TSquareMatrix3D(int n) 
//  : structure(new TSquareStructure3D(n)), TSquareMatrix(structure)
// {
// }

TSquareMatrix3D::~TSquareMatrix3D()
{
}
