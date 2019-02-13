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
// @(#)Structure3D.h        1.3 09/14/99
// 
// Class:       TStructure3D
//
// Purpose:     build and store a matrix Structure3D
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE3D__
#define __STRUCTURE3D__

#include <FESpace2D.h>
#include <FESpace3D.h>
#include <Structure.h>

class TStructure3D : public TStructure
{
  protected:
    /** Ansatzspace */
    TFESpace2D *AnsatzSpace2D;
    TFESpace3D *AnsatzSpace3D;

    /** Testspace */
    TFESpace2D *TestSpace2D;
    TFESpace3D *TestSpace3D;

  public:
    /** generate the matrix Structure3D, both space with 3D collection */
    TStructure3D(TFESpace3D *testspace, TFESpace3D *ansatzspace);

    /** return AnsatzSpace */
    TFESpace *GetAnsatzSpace()
    {
      if (AnsatzSpace2D)
        return AnsatzSpace2D;
      else
        return AnsatzSpace3D;
    }

    /** return TestSpace */
    TFESpace *GetTestSpace()
    {
      if (TestSpace2D)
        return TestSpace2D;
      else
        return TestSpace3D;
    }

};

#endif
