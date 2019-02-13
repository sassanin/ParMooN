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
// @(#)Structure2D.h        1.3 09/14/99
// 
// Class:       TStructure2D
//
// Purpose:     build and store a matrix Structure2D
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE2D__
#define __STRUCTURE2D__

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <Structure.h>

class TStructure2D : public TStructure
{
  protected:
    /** Ansatzspace */
    TFESpace1D *AnsatzSpace1D;
    TFESpace2D *AnsatzSpace2D;

    /** Testspace */
    TFESpace1D *TestSpace1D;
    TFESpace2D *TestSpace2D;

    int *AnsatzMortarSpaceGlobNo;
    int *TestMortarSpaceGlobNo;

    int *AnsatzNonMortarSpaceGlobNo;
    int *TestNonMortarSpaceGlobNo;

  public:
    /** generate the matrix Structure2D, both space with 2D collection */
    TStructure2D(TFESpace2D *testspace, TFESpace2D *ansatzspace);

    /** destructor: free all used arrays */
    ~TStructure2D();

    /** generate the matrix structure, both spaces are 2D */
    /** both spaces are defined on different grids */
    TStructure2D(TFESpace2D *testspace, int test_level, 
                 TFESpace2D *ansatzspace, int ansatz_level);
    #ifdef __MORTAR__
    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
    TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace);
    #endif

    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
     TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace, int **ansatzcelljoints);
     
    /** generate the matrix Structure2D, one space with 1D and the other
        with 2D collection */
     TStructure2D(TFESpace1D *testspace, TFESpace2D *ansatzspace, TNonMortarData *NonMortarFEData);
     
     TStructure2D(TFESpace2D *testspace, TFESpace1D *ansatzspace, TNonMortarData *NonMortarFEData);

    /** return AnsatzSpace */
    TFESpace2D *GetAnsatzSpace2D() const
    { return AnsatzSpace2D; }
    
    /** return AnsatzSpace */
    TFESpace *GetAnsatzSpace()
    {
      if (AnsatzSpace1D)
        return AnsatzSpace1D;
      else
        return AnsatzSpace2D;
    }
    
    /** return TestSpace */
    TFESpace2D *GetTestSpace2D() const
    { return TestSpace2D; }
    
    /** return TestSpace */
    TFESpace *GetTestSpace()
    {
      if (TestSpace1D)
        return TestSpace1D;
      else
        return TestSpace2D;
    }

};

#endif
