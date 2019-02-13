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
// @(#)Matrix3D.C        1.2 11/20/98
// 
// Class:       TMatrix3D
//
// Purpose:     store a  matrix3D (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Matrix3D.h>
#include <string.h>

TMatrix3D::TMatrix3D(TStructure3D *structure)
 : TMatrix(structure)
{
  this->structure = structure;
}

TMatrix3D::~TMatrix3D()
{
}


void TMatrix3D::resetNonActive()
{
  int n_active = this->structure->GetTestSpace()->GetN_DegreesOfFreedom()
                -this->structure->GetTestSpace()->GetN_Dirichlet();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active];
  int n_nonactive_entries = rowPtr[this->structure->GetN_Rows()]
                           - index_nonactive;
  memset(Entries + index_nonactive, 0.0, n_nonactive_entries * SizeOfDouble);
}
