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
// @(#)Vertex.C        1.1 10/30/98
// 
// Class:       TVertex
// Purpose:     a vertex in a grid 
//
// Author:      Volker Behns  09.07.97
//              Sashikumaar Ganesan 05.11.09 (added parallel methods)
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Vertex.h>
#include <string.h>
#include <Database.h>

// Constructors

#ifdef __3D__
TVertex::TVertex(double initX, double initY, double initZ)
{
  X = initX;
  Y = initY;
  Z = initZ;

#ifdef _MPI
  N_Cells = 0;
  SubDomainVert = FALSE;
  N_SubDomains = 0;
  SubDomain_Ranks = NULL;
  SubDomainGlobalCellNo = NULL;
  SubDomainLocVertNo = NULL;
  N_CrossNeibCells = 0;
  Cells = NULL;
  Bd_id = -1; // necessary for gmesh implementation
#endif

  BoundVert=FALSE; 
  ClipBoard = 0;
}
#else
TVertex::TVertex(double initX, double initY)
{
  X = initX;
  Y = initY;
  
  BoundVert=FALSE; 
  ClipBoard = 0;  
}
#endif

// Methods

// set coordinates
#ifdef __3D__
void TVertex::SetCoords(double initX, double initY, double initZ)
{
  X = initX;
  Y = initY;
  Z = initZ;
}
#else
void TVertex::SetCoords(double initX, double initY)
{
  X = initX;
  Y = initY;
}
#endif

#ifndef __3D__
std::ostream& operator << (std::ostream& s, TVertex *v)
{
  return s << " X= " << setw(12) << v->X << ", Y= " << setw(12) << v->Y;
}
#else
std::ostream& operator << (std::ostream& s, TVertex *v)
{
  return s << " X= " << setw(12) << v->X << ", Y= " << setw(12) << v->Y
           << ", Z= " << setw(12) << v->Z;
}
#endif

#ifdef _MPI

void TVertex::SetVertexCells(int n_Cells, TBaseCell **cells)
{

 int i;

  N_Cells = n_Cells;
  Cells = new TBaseCell*[N_Cells];

  for(i=0;i<N_Cells;i++)
   Cells[i] = cells[i];
}


#ifdef __3D__
/** add only one cell (lowest index cell) from each subdomain as a Hallo cell for this vertex */
void TVertex::SetSubDomainInfo(int n_SubDomains, int *subDomain_Ranks, int *subDomainGlobalCellNo, 
                               int *subDomainLocVertNo)
{
 int i;

   if(N_SubDomains)
    delete [] SubDomain_Ranks;

   N_SubDomains = n_SubDomains;
   SubDomain_Ranks = new int[N_SubDomains];
   SubDomainGlobalCellNo = new int[N_SubDomains];
   SubDomainLocVertNo = new int[N_SubDomains];


   for(i=0;i<N_SubDomains;i++)
    {
     SubDomain_Ranks[i] = subDomain_Ranks[i];
     SubDomainGlobalCellNo[i] = subDomainGlobalCellNo[i];
     SubDomainLocVertNo[i] = subDomainLocVertNo[i];
    }
}

void TVertex::AddCrossNeib(int Neib_ID)
{
 int i, tmp;
 
  CrossVert = TRUE;
 
    for(i=N_CrossNeibCells;i<N_SubDomains;i++)
     if(SubDomain_Ranks[i] == Neib_ID) // swap, cross ID will be at the begning
     {
      tmp = SubDomain_Ranks[N_CrossNeibCells];
      SubDomain_Ranks[N_CrossNeibCells] = SubDomain_Ranks[i];
      SubDomain_Ranks[i] = tmp;

      tmp = SubDomainGlobalCellNo[N_CrossNeibCells];
      SubDomainGlobalCellNo[N_CrossNeibCells] = SubDomainGlobalCellNo[i];
      SubDomainGlobalCellNo[i] = tmp;

      tmp = SubDomainLocVertNo[N_CrossNeibCells];
      SubDomainLocVertNo[N_CrossNeibCells] = SubDomainLocVertNo[i];
      SubDomainLocVertNo[i] = tmp;

      N_CrossNeibCells++;
     }
    
}

#endif

#endif

TVertex::~TVertex()
{

#ifdef _MPI
//  if(N_Cells)
//    delete [] Cells; // cells itself will be deleted seperately
// 
//  if(N_SubDomains)
//   {
//    delete [] SubDomain_Ranks;
//    delete [] SubDomainGlobalCellNo;
//    delete [] SubDomainLocVertNo;
//   }
#endif

}
