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
// @(#)ADISystem.C        4.1 13.11.09
// 
// Class:       TADISYSTEM
// Purpose:     general super class for all ADI/operator-splitting System
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (13.11.09)
//
// History:     start of implementation 13.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#include <ADISystem.h>

TADISystem::TADISystem(double *Sol, double *OldSol, double *b, double *Defect, DoubleFunctND *growthAndNuc)
{
  sol = Sol;
  oldsol = OldSol;
  rhs = NULL;
  B = b;
  defect = Defect;
  GetGrowthAndNuc = growthAndNuc;  
}

TADISystem::~TADISystem()
{

}
