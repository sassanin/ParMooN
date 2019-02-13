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
// @(#)ADICell.C        4.1 07.11.09
// 
// Class:       src for TADICell
// Purpose:     general super class for all ADICells
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (07.11.09)
//
// History:     start of implementation 07.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#include <ADICell.h>
#include <Database.h>

TADICell::TADICell(TBaseCell *cell, int N_quadPts)
{
  Cell = cell;
  N_QuadPts = N_quadPts;

  sol = NULL;
  solT = NULL;
  rhs = NULL;
  rhsT = NULL;
}

TADICell::~TADICell()
{

}
