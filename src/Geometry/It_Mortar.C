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
// @(#)It_Mortar.C        1.3 07/28/99
// 
// Class:       TIt_Mortar
// Purpose:     iterator which gets all cells on a given mortar face
//
// Author:      Volker Behns  20.03.98
//
// =======================================================================

#ifndef __MORTAR__
#define __MORTAR__
#endif

#include <It_Mortar.h>
#include <MortarBaseJoint.h>
#include <RefDesc.h>

// Methods
int TIt_Mortar::Init(int level)
{
  TMortarFace *Tmp;
  const int *TmpoEnE, *TmpLen;
  int MaxLen, MortarFaceID;

  MortarSide = (bool) (level < 0);
  if (MortarSide) level = - level;

  TIt_Search::Init(level);

  MortarFaceID = Level >> 8;
  Level -= MortarFaceID << 8;

  Tmp = Domain->GetMortarFace(MortarFaceID);
  ActiveCell = Tmp->Cell;

  if (MortarSide)
  {
    ActiveCell = ((TMortarBaseJoint *) ActiveCell->
                 GetJoint(Tmp->LocFaceNumber[0]))->GetNeighbour(ActiveCell);
    LocMortarFace[0] = Tmp->LocFaceNumber[1];
  }
  else
    LocMortarFace[0] = Tmp->LocFaceNumber[0];

  ActiveCell->GetRefDesc()->GetOldEdgeNewEdge(TmpoEnE, TmpLen, MaxLen);
  if (TmpLen)
    Status[0].N_Children = TmpLen[LocMortarFace[0]] + 1;
  else
    // set to 1 in order to prevent a stop in first step
    Status[0].N_Children = 1;

  return 0;
}

int TIt_Mortar::GetMaxLevel()
{
  return -1;
}

TBaseCell *TIt_Mortar::Next(int &info)
{
  const int *TmpoEnE, *TmpEC, *TmpECI, *TmpLen1, *TmpLen2;
  int MaxLen1, MaxLen2, LocEdge;

  if (ActiveLevel)
    do
    {
      ActiveCell = ActiveCell->GetParent();
      Status[ActiveLevel].N_Children = 0;
      Status[ActiveLevel--].CurrentChild = 0;
    } while (ActiveLevel &&
        Status[ActiveLevel].N_Children == Status[ActiveLevel].CurrentChild);

  if (!ActiveLevel)
  {
    if (Status[0].N_Children == Status[0].CurrentChild)
      return NULL;
    else
    {
      if (!ActiveCell->GetN_Children())
      {
        Status[0].CurrentChild++;
        info = LocMortarFace[0];
        return ActiveCell;
      }
    }
  }

  ActiveCell->GetRefDesc()->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);

  while (Status[ActiveLevel].N_Children && ActiveLevel != Level)
  {
    ActiveCell->GetRefDesc()->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);
    ActiveCell->GetRefDesc()->GetEdgeChildIndex(TmpECI, TmpLen2, MaxLen2);

    if (MortarSide)
      LocEdge = TmpoEnE[MaxLen1 * LocMortarFace[ActiveLevel] +
                  Status[ActiveLevel].N_Children
                  - ++Status[ActiveLevel].CurrentChild];
    else
      LocEdge = TmpoEnE[MaxLen1 * LocMortarFace[ActiveLevel] +
                  Status[ActiveLevel].CurrentChild++];

    ActiveCell = ActiveCell->GetChild(TmpEC[MaxLen2 * LocEdge]);

    if (ActiveCell->GetN_Children())
    {
      ActiveCell->GetRefDesc()->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);

      LocMortarFace[++ActiveLevel] = TmpECI[MaxLen2 * LocEdge];

      Status[ActiveLevel].N_Children =
               TmpLen1[LocMortarFace[ActiveLevel]] + 1;
    }
    else
      Status[++ActiveLevel].N_Children = 0;
  }
  
  info = TmpECI[MaxLen2 * LocEdge];
  return ActiveCell;
}
