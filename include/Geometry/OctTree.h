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
   
#ifndef __OCTREE__
#define __OCTREE__

#include <Collection.h>
#include <BaseCell.h>

class TOctTree
{
  public:
    typedef struct {
      double StartX;
      double StartY;
      double StartZ;
      double BoundX;
      double BoundY;
      double BoundZ;
    } TBoundBox;
  
  protected:   
    class TNode
    {
      protected:
      	/// parameters describing box
	TBoundBox Box;
	
	/// parent
	TNode *Parent;
	
	/// sub boxes
	TNode *Childs[8];
	
	///
	int Depth;
	
	int N_Cells;
	TBaseCell **Cells;

      protected:
	void MakeSubBoxes();
	bool Intersect(TBaseCell *Cell);
	bool PointInBox(double x, double y, double z);
		
      public:
	TNode (TNode *Parent, TBoundBox *box, int n_cells, TBaseCell **cells);	
	~TNode ();
	
	TNode *GetLeaf(double x, double y, double z);
	void GetCells(int &n_cells, TBaseCell **&cells);
    };
 
  protected:
    
    /// Collection from which the octtree is build
    TCollection *Coll;
    
    /// head of the tree
    TNode *Head;
    
    /// bounding box for whole domain
    TBoundBox Bound;
    
  protected:
    void BuildTree();
    
  public:
    TOctTree(TCollection *coll, TBoundBox *bounds); 
    ~TOctTree();
    
    void GetCells(double x, double y, double z, int &N_Cells, TBaseCell **&Cells);
};

#endif
