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
// @(#)Vertex.h        1.1 10/30/98
// 
// Class:       TVertex
// Purpose:     a vertex in a grid
//
// Author:      Volker Behns  09.07.97
//              Sashikumaar Ganesan 05.11.09 (added parallel methods)
//              Sashikumaar Ganesan 08.09.2010 (added 3D parallel methods)
// =======================================================================

#ifndef __VERTEX__
#define __VERTEX__

#include <MooNMD_Io.h>
#include <Constants.h>

/** a vertex in a grid */
class TVertex
{
  protected:
    /** first coordinate */
    double X;
    /** second coordinate */
    double Y;
#ifdef __3D__
      /** third coordinate (3D) */
      double Z;
#endif

    /** an integer for storing clipboard information*/
    int ClipBoard;

#ifdef _MPI

    /** Number of 3D cells containing this cells **/
    /** Note !this info only set for dependent cells !!!!!!! */
    int N_Cells;

    /** cells */
    /** Note ! this info only set for dependent cells !!!!!!!!!!*/
    TBaseCell **Cells;

    /** marking this vertex as subdomain vertex */
    bool SubDomainVert;

    /** marking this vertex as cross vertex */
    bool CrossVert;

    /** an integer which stores the number of ranks (SubDomains) contain this vertex */
    int N_SubDomains;

    /** an integer which stores the rank of SubDomains, which contain this vertex */
    int *SubDomain_Ranks;

    /** list of neib cell Globalnumbers, which incident this vertex */
    int *SubDomainGlobalCellNo;

    /** list of neib cell local vert no */
    int *SubDomainLocVertNo;

    /** an integer which stores the number of Cross neib cells, which incident this vertex */
    int N_CrossNeibCells;  
    
    int Bd_id; 
    
#endif

    
    /** marking this vertex as Bound vertex */
    bool BoundVert;   
    
  public:
    // Constructors

#ifdef __3D__
      /** 3D vertex */
      TVertex(double initX, double initY, double initZ);
#else
      /** 2D vertex */ 
      TVertex(double initX, double initY);
#endif

    // Destructor
    ~TVertex();

    // Methods

    // set coordinates
#ifdef __3D__
      /** set the coordinates in 3D */
      void SetCoords(double initX, double initY, double initZ);
#else
      /** set the coordinate in 2D */
      void SetCoords(double initX, double initY);
#endif

    /** return the x coordinate */
    double GetX() const
    { return X; }
    /** return the y coordinate */
    double GetY() const
    { return Y; }
#ifdef __3D__
      /** return the z coordinate (3D) */
      double GetZ() const
      { return Z; }
      /** return all three coordinates */
      void GetCoords(double& x, double& y, double& z) const
      {
        x = X;
        y = Y;
        z = Z;
      }
#else
      /** return all two coordinates */
      void GetCoords(double& x, double& y) const
      {
        x = X;
        y = Y;
      }
#endif

    /** write some information of the vertex in stream s */
    friend std::ostream& operator << (std::ostream& s, TVertex *v);

    /** set value in ClipBoard */
    void SetClipBoard(int value)
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard() const
    { return ClipBoard; }

     void SetAsBoundVert()
      { BoundVert=TRUE; }    
      
     bool IsBoundVert() const
     { return BoundVert; }     
    
#ifdef _MPI

    /** Note ! this info only set for dependent cells !!!!!! */
    void SetVertexCells(int n_Cells, TBaseCell **cells);

    void SetSubDomainInfo(int n_SubDomains, int *subDomain_Ranks, int *subDomainGlobalCellNo, 
                          int *subDomainLocVertNo);

    void AddCrossNeib(int Neib_ID);

    void SetAsSubDomainVert()
     { SubDomainVert = TRUE; }

     bool IsSubDomainVert()
     { return SubDomainVert; }

     void SetAsCrossVert()
     {  CrossVert = TRUE; }
     
     bool IsCrossVert()
     { return CrossVert; }

     void GetCrossNeibs(int &n_VertCrossNeibs, int *&vertCrossNeibs) 
      {
       n_VertCrossNeibs = N_CrossNeibCells;
       vertCrossNeibs = SubDomain_Ranks;
      }

     void GetCrossNeibsInfo(int &N_NeibCells, int *&NeibCellRank, 
                            int *&GlobalNo, int *&LocVertNo)
      {
       N_NeibCells = N_CrossNeibCells;
       NeibCellRank = SubDomain_Ranks;
       GlobalNo = SubDomainGlobalCellNo;
       LocVertNo = SubDomainLocVertNo;
      }

     void GetNeibs(int &n_Neibs, TBaseCell **&neighbs)
      {
       n_Neibs = N_Cells;
       neighbs = Cells;
      }
      
      int GetNNeibs()
      {  return N_Cells; } 
#ifdef _MPI
      void set_Bd_id(int key)
      {      
	Bd_id = key;
      }
      
      int get_Bd_id()
      {      
	return Bd_id;
      }
#endif    
#endif
};

#endif
