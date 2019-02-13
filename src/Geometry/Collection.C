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
// @(#)Collection.C        1.2 08/12/99
//
// Class:       TCollection
// Purpose:     store cells in an array
//              used by cell iterators
//
// Author:      Gunar Matthies  14.10.97
//
// History:     14.10.97 Starting implementation
// =======================================================================

#include <algorithm>
#include <Collection.h>
#include <BaseCell.h>
#include <string.h>
#include <JointCollection.h>   
#include <IsoBoundEdge.h>
#include <BoundComp.h>

/** constructor */
TCollection::TCollection(int n_cells, TBaseCell **cells)
{
  N_Cells = n_cells;
  Cells = cells;

  SortedCells = NULL;
  Index = NULL;

  #ifdef  _MPI
  N_OwnCells = 0;
  GlobalIndex = new int[N_Cells];

  for(int i=0; i<N_Cells; i++)
   GlobalIndex[i] = Cells[i]->GetGlobalCellNo();
  #endif
}

void TCollection::GenerateSortedArrays()
{
  int i;

  if(!SortedCells)
  {
    SortedCells = new TBaseCell*[N_Cells];
    Index = new int[N_Cells];

    memcpy(SortedCells, Cells, N_Cells*SizeOfPointer);
    std::sort(SortedCells, SortedCells+N_Cells);

    for(i=0;i<N_Cells;i++)
      Index[GetSortedIndex(Cells[i])] = i;
  }
}

/** destructor: delete arrays */
TCollection::~TCollection()
{
  if(Index) delete [] Index;
  if(SortedCells) delete [] SortedCells;
  if(Cells) delete [] Cells;
}

/** get maximal and minimal diameter */
int TCollection::GetHminHmax(double *hmin, double *hmax)
{
  int i;
  double h_min = 1e10 , h_max= 0, h;
  TBaseCell *cell;
    
  for (i=0;i<N_Cells;i++)
  {
    cell = GetCell(i);
    h = cell->GetDiameter();
    if (h<h_min)
      h_min = h;
    if (h>h_max)
      h_max = h;
  }
  *hmin = h_min;
  *hmax = h_max;
  return(0);
}

/** return Index of cell in SortedCells-array */
int TCollection::GetSortedIndex(TBaseCell *cell)
{
  int Left = 0;
  int Mid;
  int Right = N_Cells - 1;
  int ret = -1;
  TBaseCell *b;

  while (Left <= Right )
  {
    Mid = Left + ((Right - Left) / 2);

    b = SortedCells[Mid];

    if (b == cell)
    {
      ret = Mid;
      break;
    }

    if(b > cell)
      Right = Mid - 1;
    else
      Left = Mid + 1;
  }
  return ret;
}

/** return Index of cell in SortedCells-array */
int TCollection::GetIndex(TBaseCell *cell)
{
  int ret, gsi;

  GenerateSortedArrays();

  gsi = GetSortedIndex(cell);
  ret = (gsi==-1)?(N_Cells+1):Index[gsi];

  return ret;
}

 /** return the Index of the vertex in the sorted array */
int TCollection::GetIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l=0, r=Length, m=(r+l)/2;
  TVertex *Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element)
    {
      l=m;
    }
    else
    {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }

  return m;
}

/** mark the vertices that are on the boundary */
int TCollection::MarkBoundaryVertices()
{
  int i, j, N_, comp;
  double x0, y0, x2, y2, t0, t1, eps=1e-8;
  const int *TmpEdVer;  
  TBaseCell *cell;
  TVertex *vertex0, *vertex1;
  TRefDesc *refdesc;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  TJoint *joint;

  // initialization 
  // loop over all mesh cells
  for (i=0;i<N_Cells;i++)
  {
    cell = GetCell(i);
    N_ = cell->GetN_Vertices();
    // loop over the vertices
    for (j=0;j<N_;j++)
    {
     vertex1 = cell->GetVertex(j);
     vertex1->SetClipBoard(-1);
    }
  }

  // set ClipBoard for vertices on the boundary
  for (i=0;i<N_Cells;i++)
  {
    cell = GetCell(i);
    // get refinement descriptor
    refdesc=cell->GetRefDesc();                   
    // get information to compute vertices from edges
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
    // number of edges
    N_ = cell->GetN_Edges();
    for (j=0;j<N_;j++)
    {
      joint = cell->GetJoint(j);
	// if boundary edge
	if (joint->GetType() == BoundaryEdge ||
	    joint->GetType() == IsoBoundEdge)
	{
	  boundedge=(TBoundEdge *)joint;
	  BoundComp=boundedge->GetBoundComp();      // get boundary component
	  comp=BoundComp->GetID();                  // boundary id
	  boundedge->GetParameters(t0, t1);         // parameter interval
	  vertex0 = cell->GetVertex(TmpEdVer[2*j]);
	  x0 = vertex0->GetX();
	  y0 = vertex0->GetY();
	  vertex1 = cell->GetVertex(TmpEdVer[2*j+1]);
	  boundedge->GetXYofT(t0,x2,y2);
	  if ((fabs(x0-x2)<eps)&&(fabs(y0-y2)<eps))
	    {
	       vertex0->SetClipBoard(comp+t0);
	       vertex1->SetClipBoard(comp+t1);
	    }
	  else
	    {
	       vertex0->SetClipBoard(comp+t1);
	       vertex1->SetClipBoard(comp+t0);
	    }
	  //OutPut(t0 << " " << vertex0->GetClipBoard() << " " << vertex1->GetClipBoard() << endl);
	}
    }
  }
  
  return(0);
  
}

// methods for TJointCollection, 03.11.09  (Sashi)
static void Sort(TJoint **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, *rr, len;
  TJoint *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}

static void Sort(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, *rr, len;
  TVertex *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}

/** return Index of joints in Cells-array */
TJointCollection *TCollection::GetJointCollection()
{
 int i, j, N_Joints, N, N_RootJoints;
 TBaseCell *Me;
 TJoint **joints, **RootJoints, *Last, *Current;
 TJointCollection *JointColl;


 #ifdef __3D__ 
 joints = new TJoint*[6*N_Cells];
 #else 
 joints = new TJoint*[4*N_Cells];
 #endif

 N=0;
 for(i=0; i<N_Cells; i++)
  {
   Me = Cells[i];
   Me->SetCellIndex(i); // needed for DG matrices assembling
   N_Joints = Me->GetN_Joints();
   
   for(j=0; j<N_Joints; j++) 
    {
     joints[N] = Me->GetJoint(j);
     N++;
    } //for(j=0; j<N_Joints; j++) 
  } //for(i=0; i<N_Cells; i++)


  N--;
  // sort the Vertices array
  Sort(joints, N);
  N++;

  Last=NULL;
  N_RootJoints=0;
  for(i=0;i<N;i++)
   {
    Current=joints[i];
    if(Current!=Last)
    {
      N_RootJoints++;
      Last=Current;
    }
   }

  RootJoints =  new TJoint*[N_RootJoints];
  Last=NULL;
  N_RootJoints=0;
  for(i=0;i<N;i++)
   {
    Current=joints[i];
    if(Current!=Last)
    {
      RootJoints[N_RootJoints] = Current;
      Last = Current;
      N_RootJoints++;
    }
   }

  JointColl = new TJointCollection(N_RootJoints, RootJoints);

 return JointColl;
}

// for operator-split nodal point collection, 14.07.2010 (Sashi)
void TCollection::GenerateCellVertNeibs()
{
 int i, j, k, m, N, N_VertInCell, N_LocVertices;
 int Max_N_VertInCell, N_RootVertices, *NumberVertex, *VertexNumbers;
 TVertex *Current, *Last, **Vertices;

   Max_N_VertInCell = 0;
   for(i=0;i<N_Cells;i++)
    if(Max_N_VertInCell<Cells[i]->GetN_Vertices())
      Max_N_VertInCell=Cells[i]->GetN_Vertices();

   Vertices=new TVertex*[Max_N_VertInCell*N_Cells];

   N=0;
   for(i=0;i<N_Cells;i++)
    {
     N_VertInCell = Cells[i]->GetN_Vertices();

     for(j=0;j<N_VertInCell;j++)
      {
       Vertices[N]=Cells[i]->GetVertex(j);
       N++;
      } // j
    } // i
  N_LocVertices = N;
  N--;
  NumberVertex =new int[N_LocVertices];
  VertexNumbers= new int[N_LocVertices];
  // sort the Vertices array based on vertices pointer values
  Sort(Vertices, N);

  Last=NULL;
  N_RootVertices=-1;
  for(i=0;i<N_LocVertices;i++)
   {
    if((Current=Vertices[i])!=Last)
    {
      N_RootVertices++;
      Last=Current;
    }
    NumberVertex[i]=N_RootVertices;
   }
  N_RootVertices++;

  m=0;
  for(i=0;i<N_Cells;i++)
   {
    k=Cells[i]->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=Cells[i]->GetVertex(j);
      N=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[N];
      m++;
    } // endfor j
  } //endfor i



}
