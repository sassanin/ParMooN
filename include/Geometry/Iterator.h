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
// @(#)Iterator.h        1.2 02/08/99
// 
// Class:       TIterator
// Purpose:     iterator to produce a series of cells which lye
//              exactly on a special level
//
// Author:      Volker Behns  04.08.97
//
// =======================================================================

#ifndef __ITERATOR__
#define __ITERATOR__

#define MAX_ItLevel  50
#define N_ITERATORS   9

enum Iterators {It_EQ, It_LE, It_Finest, It_EQLevel, It_LELevel,
                It_Between, It_OCAF, It_Mortar1, It_Mortar2};

#include <Domain.h>

/** iterator to produce a series of cells with some
    special properties */
class TIterator
{
  protected:
    /** current level */
    int Level;

    /** domain on which the iterator works */
    TDomain *Domain;
    /** a copy of tree of cells */
    TBaseCell **CellTree;
    /** number of cells on trees root */
    int N_RootCells;

  public:
    // Constructors

    // Methods
    /** set all parameters to the given values */
    int SetParam(TDomain *domain);

    /** return the next cell */
    virtual TBaseCell *Next(int &info) = 0;
    /** return the previous cell */
    virtual TBaseCell *Prev() = 0;

    /** Initialize at level "level" */
    virtual int Init(int level) = 0;

    /** return the maximum level */
    virtual int GetMaxLevel() = 0;
};

#endif
