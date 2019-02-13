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
// @(#)It_Finest.C        1.1 10/30/98
// 
// Class:       TIt_Finest
// Purpose:     iterator to produce a series of cells which lie
//              on top of cell tree

//
// Author:      Volker Behns  08.08.97
//
// =======================================================================

#include <It_Finest.h>

// Methods
int TIt_Finest::Init(int level)
{
  TIt_Search::Init(level);

  Level = MAX_ItLevel;

  return 0;
}

int TIt_Finest::GetMaxLevel()
{
  int MaxLevel = 0;

  do
  {
    if (ActiveLevel)
      do
      {
        ActiveCell = ActiveCell->GetParent();
        Status[ActiveLevel].N_Children = 0;
        Status[ActiveLevel--].CurrentChild = 0;
      } while (ActiveLevel &&
          Status[ActiveLevel].N_Children == Status[ActiveLevel].CurrentChild);

    if (!ActiveLevel &&
         Status[ActiveLevel].N_Children == Status[ActiveLevel].CurrentChild)
    {
      if (ActiveRootCell < N_RootCells)
      {
        ActiveCell = CellTree[ActiveRootCell++];
        Status[ActiveLevel].N_Children = ActiveCell->GetN_Children();
        Status[ActiveLevel].CurrentChild = 0;
      }
      else
        break;
    }

    while (Status[ActiveLevel].N_Children)
    {
      ActiveCell = ActiveCell->GetChild(Status[ActiveLevel].CurrentChild++);
      Status[++ActiveLevel].N_Children = ActiveCell->GetN_Children();

      if (ActiveLevel > MaxLevel) MaxLevel = ActiveLevel;
    }

  } while (TRUE);

  return MaxLevel;
}
