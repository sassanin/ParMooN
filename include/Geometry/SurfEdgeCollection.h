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
   
#ifndef __SURFEDGECOLLECTION__
#define __SURFEDGECOLLECTION__

#include <Collection.h>

// #define iFLIP(x,y) if ( x > y ) { int ivar; ivar = a; a = b; b = ivar; }

class TSurfEdgeCollection
{
  protected:
    /** number of vertices */
    int mN_Vertices;
    
    /** number of edges */
    int mN_Edges;
    
    /** hash table for edges */
    int **mEdgeHash;
  
  protected:
    void Init(TCollection *Coll, int N_SurfaceJoints,
	      int *CellNumbers, int *JointNumbers);
	      
    inline int Flip(int &a, int &b)
    {
      if ( a > b )
      { 
	int ivar;
	ivar = a; a = b; b = ivar;
      }
      
      return a+b;
    }
    
    void IncreaseBucket(int a, int b, int &counter);
    int FindEdge(int *Bucket, int a, int b);
    
    
  public:
    TSurfEdgeCollection (TCollection *Coll, int N_SurfaceJoints,
			 int *CellNumbers, int *JointNumbers);
			 
    ~TSurfEdgeCollection ();
    
    int GetN_Edges()
    { return mN_Edges; }
    
    int GetN_Vertices()
    { return mN_Vertices; }
    
    int GetEdge(int a, int b)
    { 
      if ( a < 0 || b < 0 ) return -1;
      else return FindEdge(mEdgeHash[Flip(a,b)], a, b);
    }
};
#endif
