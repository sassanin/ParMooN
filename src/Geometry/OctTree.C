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
   
// #include <OctTree.h>
// 
// #include <stdlib.h>
// 
// #define MAXDEPTH 5
// #define EPS 1e-8
// 
// // CTOR
// TOctTree::TOctTree(TCollection *coll, TBoundBox *bound)
// {  
//   Coll = coll;
//   
//   Bound.StartX = bound->StartX;
//   Bound.StartY = bound->StartY;
//   Bound.StartZ = bound->StartZ;
//   Bound.BoundX = bound->BoundX;
//   Bound.BoundY = bound->BoundY;
//   Bound.BoundZ = bound->BoundZ;
//   
//   Head = new TNode (NULL, &Bound, Coll->GetN_Cells(), Coll->GetCells());
// }
// 
// TOctTree::TNode::TNode (TNode *parent, TBoundBox *box, int n_cells, TBaseCell **cells)
// {
//   int counter;
//   TBaseCell *Cell;
//   
//   Box.StartX = box->StartX - EPS;
//   Box.StartY = box->StartY - EPS;
//   Box.StartZ = box->StartZ - EPS;
//   Box.BoundX = box->BoundX + 2*EPS;
//   Box.BoundY = box->BoundY + 2*EPS;
//   Box.BoundZ = box->BoundZ + 2*EPS;
//   
//   N_Cells = 0;
//   Cells = NULL;
//   Depth = 0;
//   Parent = parent;
//   
//   if ( parent || 1)
//   {  
//     for (int i=0;i<n_cells;++i)
//     {
//       Cell = cells[i];
//       
//       if ( Intersect(Cell) )
//       {
// 	Cell->SetClipBoard(1);
// 	++N_Cells;
//       }
//       else
//       {
// 	Cell->SetClipBoard(0);
//       }
//       
//       if ( Parent == NULL )
// 	Cell->SetCellIndex(i);
//     }
//     
//     Cells = new TBaseCell* [N_Cells];
//     
//     counter = 0;
//     for (int i=0;i<n_cells;++i)
//     {
//       Cell = cells[i];
//       
//       if ( Cell->GetClipBoard() == 1 )
//       {
// 	Cells[counter++] = Cell;
//       }      
//     }
//   }
//   
//   if ( parent ) Depth = parent->Depth + 1;
//   
//   if ( Depth < MAXDEPTH && N_Cells > 1)
//   {
//       MakeSubBoxes();
//   }
// //   else
// //   {
// //     OutPut("Depth: " << Depth);
// //     OutPut("  |  N_Cells: " << N_Cells << endl);
// //   }
//   
//   Box.BoundX = Box.StartX + Box.BoundX;
//   Box.BoundY = Box.StartY + Box.BoundY;
//   Box.BoundZ = Box.StartZ + Box.BoundZ;
// // 
// }
// 
// // DTOR
// TOctTree::~TOctTree ()
// {
//   delete Head;
// }
// 
// TOctTree::TNode::~TNode ()
// {
//   if ( Cells) delete [] Cells;
//   else
//   {
//     for (int i=0;i<8;++i)
//       delete Childs[i];
//   }
// }
// 
// /// Methods
// 
// void TOctTree::TNode::MakeSubBoxes()
// {
//   TBoundBox box;
//   
//   // upper 1
//   box.StartX = Box.StartX;
//   box.BoundX = 0.5*Box.BoundX;
//   box.StartY = Box.StartY;
//   box.BoundY = 0.5*Box.BoundY;
//   box.StartZ = Box.StartZ + 0.5*Box.BoundZ;
//   box.BoundZ = 0.5*Box.BoundZ;
//   Childs[0] = new TNode (this, &box, N_Cells, Cells);
//   
//   // upper 2
//   box.StartX = Box.StartX + box.BoundX;
//   box.StartY = Box.StartY;
//   Childs[1] = new TNode (this, &box, N_Cells, Cells);
//   
//   // upper 3
//   box.StartY = Box.StartY + box.BoundY;
//   Childs[2] = new TNode (this, &box, N_Cells, Cells);
//   
//   // upper 4
//   box.StartX = Box.StartX;
//   Childs[3] = new TNode (this, &box, N_Cells, Cells);
//   
//   // lower 1
//   box.StartY = Box.StartY;
//   box.StartZ = Box.StartZ;
//   Childs[4] = new TNode (this, &box, N_Cells, Cells);
//   
//   // lower 2
//   box.StartX = Box.StartX + box.BoundX;
//   Childs[5] = new TNode (this, &box, N_Cells, Cells);
//   
//   // lower 3
//   box.StartY = Box.StartY + box.BoundY;
//   Childs[6] = new TNode (this, &box, N_Cells, Cells);
//   
//   // lower 4
//   box.StartX = Box.StartX;
//   Childs[7] = new TNode (this, &box, N_Cells, Cells);
//   
//   delete [] Cells;
//   Cells = NULL;
//   N_Cells = 0;
// }
// 
// bool TOctTree::TNode::Intersect(TBaseCell *Cell)
// {
//   int N_Vertices;
//   double x, y, z, xmax, xmin, ymax, ymin, zmax, zmin;
//   double left, right, top, bottom, front, back;
//   bool status;
//   
//   left = Box.StartX;
//   right = left + Box.BoundX;
//   bottom = Box.StartZ;
//   top = bottom + Box.BoundZ;
//   front = Box.StartY;
//   back = front + Box.BoundY;
//   
//   N_Vertices = Cell->GetN_Vertices();
//   
//   // get bounding box
//   for (int i=0;i<N_Vertices;++i)
//   {
//     Cell->GetVertex(i)->GetCoords(x, y, z);
//     
//     if ( i )
//     {
//       if ( x > xmax ) xmax = x;
//       if ( x < xmin ) xmin = x;
//       
//       if ( y > ymax ) ymax = y;
//       if ( y < ymin ) ymin = y;
//       
//       if ( z > zmax ) zmax = z;
//       if ( z < zmin ) zmin = z;
//     }
//     else
//     {
//       xmax = x;
//       xmin = x;
//       ymax = y;
//       ymin = y;
//       zmax = z;
//       zmin = z;
//     }
//   }
//   
//   status = ( xmax > left ) && ( xmin < right ) && 
// 	   ( ymax > front ) && ( ymin < back ) &&
// 	   ( zmax > bottom ) && ( zmin < top ) ;
//   
//   return status;
// }
// 
// 
// void TOctTree::GetCells(double x, double y, double z, int &N_Cells, TBaseCell **&Cells)
// {
//   TNode *Node;
//   
//   Node = Head->GetLeaf(x, y, z);
//   
//   if (Node)  Node->GetCells(N_Cells, Cells);  
//   else 
//   {
//     N_Cells = 0;
//     Cells = NULL;
//   }
// }
// 
// TOctTree::TNode* TOctTree::TNode::GetLeaf(double x, double y, double z)
// {
//   TNode *Node;
//   
//   if ( PointInBox(x, y, z) )
//   {
//     if ( Cells )
//       return this;
//     else
//     {
//       for (int i=0;i<8;++i)
//       {
// 	Node = Childs[i]->GetLeaf(x, y, z);
// 	
// 	if ( Node ) 
// 	  return Node;
//       }
//       
//       cerr << __FILE__ << ":" << __LINE__ << ": could not find leaf for point (";
//       cerr << x << "," << y << "," << z << ") " << Depth << endl;
//       
//       cerr << "BBox: " << Box.StartX << " - " << Box.BoundX << endl;
//       cerr << "      " << Box.StartY << " - " << Box.BoundY << endl;
//       cerr << "      " << Box.StartZ << " - " << Box.BoundZ << endl;
//       
//       exit(0);
//     }
//   }
//   
//   else return NULL;
// }
// 
// void TOctTree::TNode::GetCells(int &n_cells, TBaseCell **&cells)
// {
//   n_cells = N_Cells;
//   cells = Cells;
// }
// 
// bool TOctTree::TNode::PointInBox(double x, double y, double z)
// {
//   return ( x >= Box.StartX && x <= Box.BoundX ) &&
// 	 ( y >= Box.StartY && y <= Box.BoundY ) &&
// 	 ( z >= Box.StartZ && z <= Box.BoundZ ) ;
// }
// 
