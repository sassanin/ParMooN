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
// @(#)Vector.C        1.2 07/28/99
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

#include <MooNMD_Io.h>

#ifndef __AIX__
#include <Vector.h>
#endif

#ifdef __VECTOR__
template <class Type>
TVector<Type>::TVector(int i_size, int i_incr)
{
  N_Elements=0;
  IncrForListNumbers=5;

  init_size=i_size;
  init_incr=i_incr;

  Lists=new Type*[IncrForListNumbers];

  NumberOfLists=1;
  NumberOfFreeLists=4;
  NumberOfFreeEntries=init_size;
  ListForNext=0;
  IndexForNext=0;

  Lists[ListForNext]=new Type[NumberOfFreeEntries];

}

template <class Type>
Type TVector<Type>::GetElement(int i)
{
  int h,index,list;
  Type ret;
  h=i-init_size;

  if(i>=N_Elements)
  {
    ret=(Type)0;
  }
  else
  {
    if(h<0)
    {
      ret=Lists[0][i];
    }
    else
    {
      list=(h / init_incr);
      index=h - list*init_incr;
  
      ret=Lists[list+1][index];
    }
  }

  return ret;
}

template <class Type>
void TVector<Type>::SetElement(int i, Type value)
{
  int h,index,list;
  h=i-init_size;

  if(i>=N_Elements)
  {
    Error("Error in SetElement:" << endl);
  }
  else
  {
    if(h<0)
    {
      Lists[0][i]=value;
    }
    else
    {
      list=(h / init_incr);
      index=h - list*init_incr;
  
      Lists[list+1][index]=value;
    }
  }
}

template <class Type>
int TVector<Type>::AddElement(Type value)
{
  if(NumberOfFreeEntries==0)
  {
    // there is no space left on the end of list ListForNext
    if(NumberOfFreeLists==0)
    {
      int i;
      Type **NewLists;

      NewLists=new Type*[NumberOfLists+IncrForListNumbers];
      for(i=0;i<NumberOfLists;i++) NewLists[i]=Lists[i];
      delete Lists;
      Lists=NewLists;
      NumberOfFreeLists=IncrForListNumbers;
    }

    ListForNext++;
    NumberOfLists++;
    NumberOfFreeLists--;
    Lists[ListForNext]=new Type[init_incr];
    NumberOfFreeEntries=init_incr;
    IndexForNext=0;
  }

  // now there is space in all cases
  Lists[ListForNext][IndexForNext]=value;
  NumberOfFreeEntries--;
  IndexForNext++;
  N_Elements++;
  return (N_Elements-1);
}

template <class Type>
TVector<Type>::~TVector()
{
  int i;

  for(i=0;i<NumberOfLists;i++)
    delete Lists[i];
  delete Lists;
}

// #ifndef __IRIX__
#ifndef __AIX__
template class TVector<int>;
template class TVector<THangingNode *>;
#endif
// #endif

#endif
