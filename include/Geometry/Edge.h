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
// @(#)Edge.h  
// 
// Class:       TEdge
// Purpose:     superclass for edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#ifndef __EDGE__
#define __EDGE__

#include <Vertex.h>
#include <BaseCell.h>
#include <MooNMD_Io.h>
#include <Constants.h>

/** an edge in a 3D grid */
class TEdge
{
  protected:
    JointType EdgeID;

   /** Number of 3D cells containing this edge **/
    int N_Neibs;

    /** neighbours */
    TBaseCell **Neighbs;

    /** an integer for storing clipboard information*/
    int ClipBoard;

#ifdef  _MPI
  /** a bool to check this edge belongs to two or more SubDomains **/
    bool SubDomainEdge;

  /** a bool to check this SubDomainEdge is a cross edge **/
  /** i.e, this edge in their cell is the only connection between two subdomains */
    bool SubDomainCrossEdge;

   /** Number of 3D SubDomain cells containing this edge **/
    int N_SubDomainNeibs;

    /** SubDomain neighbours' rank */
    int *SubDomainNeibsRank;

    /** global cell number of the neibs' cell, which contains this edge */
    int *CrossNeibsGlobalNo;

   /** Number of 3D SubDomain cells containing this edge as cross egde**/
    int N_CrossNeibs;

    /** SubDomain neighbours' rank */
    int *CrossNeibsRank;

    /** local edge number of this joint in the neib cell */
    int *CrossNeibsLocalEdgeNo;

    /** local edge number of this joint in the neib cell */
    int *CrossEdgeMaptype;
        
    int Bd_id;

#endif

  public:
    // Constructors
    TEdge(int n_Neibs, TBaseCell **neighbs);

 // methods
    /** return type */
    JointType GetType()
    { return EdgeID; }

    /** set value in ClipBoard */
    void SetClipBoard(int value)
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard()
    { return ClipBoard; }

#ifdef  _MPI

  void InitSubDomainInfo(int rank);

  void SetAsNotCrossEdgeFor(int rank, int Neib_ID);

  void SetCrossNeibInfo(TVertex *Vert_a);

  bool IsSubDomainEdge()
   {return SubDomainEdge;}

//   void SetAsSubDomainCrossEdge()
//   {SubDomainCrossEdge = TRUE;}

  bool IsSubDomainCrossEdge()
   {return SubDomainCrossEdge;}

  void GetNeibs(int &n_Neibs, TBaseCell **&neighbs)
   {
    n_Neibs = N_Neibs;
    neighbs = Neighbs;
   }

  int GetN_CrossNeibs()
   {return N_CrossNeibs; }

  void GetNeibSubDomainRanks(int &n_SubDomainNeibs, int *&subDomainNeibsRank)
   {
     n_SubDomainNeibs = N_SubDomainNeibs;
     subDomainNeibsRank = SubDomainNeibsRank;
   }

  void GetCrossEdgeNeibs(int &n_CrossNeibs, int *&crossNeibsRank)
   {
    n_CrossNeibs = N_CrossNeibs;
    crossNeibsRank = CrossNeibsRank;
   }

  void GetCrossEdgeNeibs(int &n_CrossNeibs, int *&crossNeibsRank, int *&crossNeibsGlobalNo,
                         int *&crossNeibsLocalEdgeNo, int *&crossEdgeMaptype);

#endif

  // Destructor
  ~TEdge();

};

#endif
