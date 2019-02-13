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
// 
// Class:       TVector
//
// Purpose:     emulates a dynamic data structure for storing data
//              of arbitrary type (template)
//
// Author:      Gunar Matthies
// Version:     1.0
//
// =======================================================================

#ifndef __VECTOR__
#define __VECTOR__

#ifdef __AIX__
#pragma implementation("../../src/General/Vector.C")
#endif

#include <Constants.h>
#include <AllClasses.h>

template <class Type>
class TVector
{
  protected:
    /** number of all elements */
    int N_Elements;

    /** initial size */
    int init_size;
    
    /** increment */
    int init_incr;

    /** number of lists used */
    int NumberOfLists;

    /** number of free lists now */
    int NumberOfFreeLists;

    /** number of free entries in current list */
    int NumberOfFreeEntries;

    /** list for next entry */
    int ListForNext;
    
    /** index for next entry */
    int IndexForNext;

    /** array of arrays containing the data */
    Type **Lists;

    /** if there is no space left in the array Lists, create a new array
        which is IncrForListNumbers entries longer than the old */
    int IncrForListNumbers;

  public:
    /** constructor */
    TVector(int i_size, int i_incr);

    /** destructor = free all allocated memory */
    ~TVector();

    /** return the number of elements strored */
    int GetN_Elements()
    { return N_Elements; }

    /** return the element i */
    Type GetElement(int i);

    /** set the value of the already existing element i to value */
    void SetElement(int i, Type value);

    /** add a new element at the end */
    int AddElement(Type value);

};


#endif
