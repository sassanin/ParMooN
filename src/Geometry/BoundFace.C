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
// @(#)BoundFace.C        1.2 10/18/99
// 
// Class:       TBoundFace
// Purpose:     face on a boundary component
//
// Author:      Volker Behns  28.08.97
//
// =======================================================================

#ifndef __3D__
#define __3D__
#endif

#include <BoundFace.h>

// Constructors
TBoundFace::TBoundFace(TBoundComp3D *bdcomp, double *param1, double *param2)
 : TJoint()
{
  register int i;
  ID = BoundaryFace;

  BoundComp = bdcomp;
  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }
}

TBoundFace::TBoundFace(TBoundComp3D *bdcomp) : TJoint()
{
  ID = BoundaryFace;

  BoundComp = bdcomp;
  Param1[0] = 0.0;
  Param2[0] = 0.0;
  Param1[1] = 1.0;
  Param2[1] = 0.0;
  Param1[2] = 1.0;
  Param2[2] = 1.0;
  Param1[3] = 0.0;
  Param2[3] = 1.0;
}

// Methods
int TBoundFace::CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

// create a new instance of this class
TJoint *TBoundFace::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TBoundFace(BoundComp);
}

TJoint *TBoundFace::NewInst()
{
  return new TBoundFace(BoundComp);
}

void TBoundFace::SetParameters(double *param1, double *param2)
{
  register int i;

  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }
}

void TBoundFace::GetParameters(double *param1, double *param2)
{
  register int i;

  for(i=0;i<4;i++)
  {
    param1[i] = Param1[i];
    param2[i] = Param2[i];
  }
}
