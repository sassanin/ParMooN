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
   

#ifdef __2D__

#ifdef _MPI
#  include "mpi.h"
#endif
#include <math.h>
#include <Database.h>
#include <Domain.h>
#include <JointEqN.h>
#include <MacroCell.h>
#include <Vector.h>
#include <MooNMD_Io.h>
#include <BdSphere.h>
#include <string.h>
#include <Quadrangle.h>
#include <IsoInterfaceJoint.h>

#include <sstream>
#include <iostream> 
#include <stdlib.h>

#include <DefineParams.h>

#include <GridCell.h>
#include <JointEqN.h>
#include <Database.h>
#include <Quadrangle.h>
#include <RefDesc.h>
#include <Vertex.h>
#include <InterfaceJoint.h>
#include <IsoBoundEdge.h>

#include <Joint.h>


static int GetIndex(TVertex **Array, int Length, TVertex *Element)
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

static void Sort(TVertex **Array, int length) /// Sorting the vertex list
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

int TDomain::AdaptRefineAll() 
{
  
  int i, j, N, MaxCpV, N_Cells, N_VertInCell, N_JointsInCell, N_AllLocVert, *VertexNumbers, *PointNeighb, *vncn;
  int N_RootVertices, *NumberVertex;
  
  TBaseCell *cell, *neib_cell, *vcell;
  TCollection *coll;
  TVertex **Vertices, *Last, *CurrVert;
  TJoint *Joint; 
  //TIsoInterfaceJoint *isojoint, *isojoint2;
  //TIsoBoundEdge *isoboundedge, *isoboundedge2;
  
        double nrml[3];
      double tvect[3];
      double D[100000];
      double dprod[4];
      double xcoords[4];
      double ycoords[4];
      double rangeMax = 1.0;
      double stepSize = 0.01;
      double novx[1001];
      double novy[1001];
      
  
  MaxCpV = 0;
  coll = this->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  cell = coll->GetCell(0);
  N_VertInCell = cell->GetN_Vertices();
  N_JointsInCell = cell->GetN_Joints();
  N_AllLocVert = N_VertInCell*N_Cells;
  
  NumberVertex=new int[N_AllLocVert];
  VertexNumbers= new int[N_AllLocVert];
  Vertices=new TVertex*[N_AllLocVert];
  
//   GlobalIndex = new int[N_Cells];
//   Cells=cells;
//   for(i=0; i<N_Cells; i++)
//   {
//     GlobalIndex[i] = Cells[i]->GetGlobalCellNo();
//   }
  
  /** Set Golbal Cell Number */
  for(i=0;i<N_Cells;i++)
  {
    coll->GetCell(i)->SetGlobalCellNo(i);
  }
  
  N = 0;    
     for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      for(j=0;j<N_VertInCell;j++)
       {
        Vertices[N]=cell->GetVertex(j);
        N++;
       }
     } 
   
    // sort the Vertices array
    Sort(Vertices, N);
    
    Last=NULL;
    N_RootVertices=-1;
    for(i=0;i<N_AllLocVert;i++)
     {
       if((CurrVert=Vertices[i])!=Last)
        {
         N_RootVertices++;
	
         Last=CurrVert;
        }
       NumberVertex[i]=N_RootVertices;
      }
     N_RootVertices++;
    
    int m;
    m=0;
     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
 
       for(j=0;j<N_VertInCell;j++)
        {
         CurrVert=cell->GetVertex(j);
         N=GetIndex(Vertices, N_AllLocVert, CurrVert);
         VertexNumbers[m]=NumberVertex[N];
         m++;
        } // endfor j
      } //endfor i
   
// cout<< "N_RootVertices  " << N_RootVertices << endl;  
    
   
   //collect the info for subdomain mesh manupulation
   PointNeighb = new int[N_RootVertices];
   memset(PointNeighb, 0, N_RootVertices*SizeOfInt);
 
     // find cells per vertex 
     for (i=0;i<N_AllLocVert;i++)
      PointNeighb[VertexNumbers[i]]++;
   
     // find maximum cells per vertex ( maxCpV )
     for (i=0;i<N_RootVertices;i++)
      if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];
  
     printf("Max number of cells per vertex %d \n", MaxCpV);
  
    /** PointNeighb's first column contains number of neib cells associated with each vertex*/
    /** further columns contain the cell numbers associated with this vertex*/
     MaxCpV++;
 
     int M;
     delete [] PointNeighb;
     PointNeighb = new int[N_RootVertices*MaxCpV];
     memset(PointNeighb, 0, (N_RootVertices*MaxCpV)*SizeOfInt);
  
     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
  
       for(j=0;j<N_VertInCell;j++)
        {
         M = VertexNumbers[i*N_VertInCell + j] *MaxCpV ;
         PointNeighb[M]++;
         PointNeighb[M + PointNeighb[M]  ] = i;
        } // for(j=0;j<k;j++) 
      } //  for(i=0;i<N_Cells;i++) 
   
     
      cout<< "AAdaptRefineAll N_Cells " << N_Cells << endl;
//    exit(0);   
    
    RefLevel++;
    TDatabase::IteratorDB[It_Finest]->Init(0);
  
  // loop over all cells
//    while (CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info))
//    {
//     CurrCell->SetRegRefine();
//     CurrCell->Refine(RefLevel);
/// POINT REFINEMENT AND LAYER REFINEMENT (FOR STRAIGHT LINE SHAPES)   
/** Point Refinement START */    

    //coll = this->GetCollection(It_Finest, 0);
    //N_Cells = coll->GetN_Cells();
    int jj, l, vn, kk;
#ifdef __2D__
    double x, y;
#else
    double x, y, z;
#endif
    //vncn = (int*) calloc(qq, sizeof(int));
    //vertexneibcellno = (int*) malloc(vn, sizeof(int));
    //int *vcn = new int[1000]; /// vertex cell number
    vncn = new int[1000]; /// vertex cell number
    // int *vertexneibcellno = new int [1000]; /// vertex neibghor's cell number 
    cout<< "Number of cells " << N_Cells << endl;
    //int mm=0; 
     vn=0;
    //int vn1=0;
     x=0.75;
     y=0.75;
#ifdef __3D__
     z=0.5;
#endif
     /// comment here if u dont want to run this routine
     /*
     for(i=0;i<N_Cells;i++)
     {
       //if current cell is the one we want to refine (write condition)
       //if (mm==39)
       cell = coll->GetCell(i);
#ifdef __2D__
       if(cell->PointInCell(x,y)) /// Note that PointInCell command in this context is used to mark the cell which contains the singularity
#else
       if(cell->PointInCell(x,y,z))
#endif       
       {
	 cell = coll->GetCell(i);
	 cell->SetRegRefine();
	 cell->Refine(RefLevel);
	 //break;
	 for(j=0;j<N_VertInCell;j++)
	 {
	   CurrVert = cell->GetVertex(j);
	   N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;
	   N_CellsInThisVert =  PointNeighb[N];
	   for(ii=1;ii<=N_CellsInThisVert;ii++)
	   {
	     VertexCellNo = PointNeighb[N + ii];
	     Vertexcell = coll->GetCell(VertexCellNo);
	     if (!Vertexcell->ExistChildren())
	     {
	       vertexneibcellno[vn] = VertexCellNo;
	       cout<< "Vertex cell numbers " << VertexCellNo << endl;
	       Vertexcell->SetRegRefine();
	       Vertexcell->Refine(RefLevel);
	       vn++;
	     } //end if (!Vertexcell->ExistChildren()) 
	   } //end for(ii=1;ii<=N_CellsInThisVert;ii++) 
	 } //end for(j=0;j<N_VertInCell;j++) 

       } //end if (mm==1) 
       mm++;
     } //end for(i=0;i<N_Cells;i++)
     /// This is written to prevent hanging nodes
      for(l=0;l<vn;l++)
      {
        cout<< "neibnos " << l << endl;  
        vcell = coll->GetCell(vertexneibcellno[l]);
        N_JointsInCell = vcell->GetN_Joints();
        cout<< "njcell " <<  vertexneibcellno[l]<< endl;  
        for(jj=0;jj<N_JointsInCell;jj++)
        {
	  Joint = vcell->GetJoint(jj);
 	  JointType jointtype = Joint->GetType();
 	  if(Joint->InnerJoint())
 	  {
 	    neib_cell = Joint->GetNeighbour(vcell);
 	    if (!neib_cell->ExistChildren())
 	    {
 	      for(kk=0;kk<N_JointsInCell;kk++)
 	      {
    	        if( neib_cell->GetJoint(kk) ==  Joint)
 	        {
 		  neib_cell->Set1Refine(kk);
    		  neib_cell->Refine(RefLevel);
		} //end if( neib_cell->GetJoint(kk) ==  Joint)
 	      } //end for(kk=0;kk<N_JointsInCell;kk++)
 	    } //end if (!neib_cell->ExistChildren())
	  } //end if(Joint->InnerJoint())
	} //end for(jj=0;jj<N_JointsInCell;jj++)
      } //end for(l=0;l<vn;l++) */

/** Point refinement END */

/** Refinement Level 2 (OPTIONAL FOR POINT REFINEMENT)*/

/*   vncn = vertexneibcellno;
     vertexneibcellno = (int*) calloc(m, sizeof(int));
     for(l=0;l<vn;l++)
     {
       cout<< "neibnos " << l << endl;  
       vcell = coll->GetCell(vncn[l]);
       for(j=0;j<N_VertInCell;j++)
	 {
	   CurrVert = cell->GetVertex(j);
	   N = VertexNumbers[vncn[l]*N_VertInCell + j]*MaxCpV ;
	   N_CellsInThisVert =  PointNeighb[N];
	   for(ii=1;ii<=N_CellsInThisVert;ii++)
	   {
	     VertexCellNo = PointNeighb[N + ii];
	     Vertexcell = coll->GetCell(VertexCellNo);
	     if (!Vertexcell->ExistChildren())
	     {
	       vertexneibcellno[vn1] = VertexCellNo;
	       cout<< "Vertex cell numbers " << VertexCellNo << endl;
	       Vertexcell->SetRegRefine();
	       Vertexcell->Refine(RefLevel);
	       vn1++;
	     } //end if (!Vertexcell->ExistChildren()) 
	   } //end for(ii=1;ii<=N_CellsInThisVert;ii++) 
	 } //end for(j=0;j<N_VertInCell;j++)
     } //end for(l=0;l<vn;l++)
     cout<< "vn1 " << vn1 << endl;
     
     for(l=0;l<vn1;l++)
     {
       cout<< "neibnos " << l << endl;  
       vcell = coll->GetCell(vertexneibcellno[l]);
       N_JointsInCell = vcell->GetN_Joints();
       cout<< "njcell " <<  vertexneibcellno[l]<< endl;  
       for(jj=0;jj<N_JointsInCell;jj++)
       {
	 Joint = vcell->GetJoint(jj);
  	 JointType jointtype = Joint->GetType();
  	 if(Joint->InnerJoint())
  	 {
  	   neib_cell = Joint->GetNeighbour(vcell);
  	   if (!neib_cell->ExistChildren())
  	   {
  	     for(kk=0;kk<N_JointsInCell;kk++)
  	     {
     	       if( neib_cell->GetJoint(kk) ==  Joint)
  	       {
  		 neib_cell->Set1Refine(kk);
     		 neib_cell->Refine(RefLevel);
	       } //end if( neib_cell->GetJoint(kk) ==  Joint)
  	     } //end for(kk=0;kk<N_JointsInCell;kk++)
  	   } //end if (!neib_cell->ExistChildren())
 	 } //end if(Joint->InnerJoint())
       } //end for(jj=0;jj<N_JointsInCell;jj++)
     } //end for(l=0;l<vn;l++) */
     
/** Refinement Level 2 END */

/** Refinement in layer for the 2D box Two Boundary Layer Problem */

/// Type 1 (PRESENTLY WORKING): Using pointincell command to mark and refine initial cells containing the discontinuity 
      /*   
      x=1.0;
      //y=0.99;
      y=1.0;
      
      double rangeMin = 0.0;
      double rangeMax = 1.0;
      double stepSize = 0.001;
      
      for(i=0;i<N_Cells;i++)
      {
	cell = coll->GetCell(i);
        for (double xx = rangeMin; xx <= rangeMax; xx+= stepSize)
        {
	  if(cell->PointInCell(xx,y))
	  {
	    if (!cell->ExistChildren())
	    {
	      cell = coll->GetCell(i);
	      cell->SetRegRefine();
	      cell->Refine(RefLevel);
	      cout<<"cellno"<<i<<endl;
	      vncn[vn]=i;
	      vn++;
	    } 
	  }
 	}
      }
      
      x=1.0;
      y=1.0;
      
      for(i=0;i<N_Cells;i++)
      {
	cell = coll->GetCell(i);
        for (double yy = rangeMin; yy <= rangeMax; yy+= stepSize)
        {
	  if(cell->PointInCell(x,yy))
	  {
	    if (!cell->ExistChildren())
	    {
	      cell = coll->GetCell(i);
	      cell->SetRegRefine();
	      cell->Refine(RefLevel);
	      cout<<"cellno"<<i<<endl;
	      vncn[vn]=i;
	      vn++;
	    } 
	  }
 	}
      } // */
      
/// Type 1 END
      
/// Type 2 (PRESENTLY WORKING): Using product of normal and tangent vectors algorithm (but specific to this problem)

      double small, xcord, ycord, xx, yy, xt, yt;
      int k, ff, r, s;
      x=1.0;
      y=1.0;
      yy=1.0;
      xx=1.0;
      

      
      for(r=0;r<100;r++)
      {
	novx[r]=rangeMax;
	rangeMax=rangeMax-stepSize;
      }
      for(r=0;r<100;r++)
      {
	novy[r]=1.0;
      }
      
      
      for(i=0;i<N_Cells;i++)
      {
	cell = coll->GetCell(i);
	for(k=0;k<N_VertInCell;k++)
	{
	  cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k]);
	  xt = xcoords[k];
	  yt = ycoords[k];
	  ff=0;
	  for (s=0;s<100;s++)
	  {
	    D[ff] = sqrt(pow(novx[s]-xt, 2) + pow(novy[s]-yt, 2));
	    ff++;
	  }
	  small=D[0];
	  for (kk=1; kk<ff;kk++)
	  {
	    if (small < D[kk])
	    {
	      small = D[kk];
	    }
	  }
	  
	  xcord = sqrt(pow(small, 2) - pow(y-yt,2));
	  nrml[1] = y-yy;
	  nrml[2] = xcord-stepSize;
	  tvect[1] = xcord-xt;
	  tvect[2] = y-yt;
	  dprod[k] = nrml[1]*tvect[1] + nrml[2]*tvect[2];
	}
	for (l=0;l<N_VertInCell;l++)
	{
	  if (dprod[l] <= 0)
	  {
	    vncn[vn]=i;
	    cell->SetRegRefine();
	    cell->Refine(RefLevel);
	    vn++;
	  }
	}
      }
      
       for(r=0;r<100;r++)
       {
 	novy[r]=rangeMax;
 	rangeMax=rangeMax-stepSize;
       }
       for(r=0;r<100;r++)
       {
 	novx[r]=rangeMax;
       }
      
      for(i=0;i<N_Cells;i++)
      {
	cell = coll->GetCell(i);
	for(k=0;k<N_VertInCell;k++)
	{
	  cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k]);
	  xt = xcoords[k];
	  yt = ycoords[k];
	  ff=0;
// 	  for (ys = rangeMax; ys >= rangeMin; ys-= stepSize)
// 	  {
	    for(s=0;s<100;s++)
	    {
	    
	    D[ff] = sqrt(pow(novx[s]-xt, 2) + pow(novy[s]-yt, 2));
	    ff++;
	  }
	  small=D[0];
	  for (kk=1; kk<ff;kk++)
	  {
	    if (small < D[kk])
	    {
	      small = D[kk];
	    }
	  }
	  
	  ycord = sqrt(pow(small, 2) - pow(x-xt,2));
	  nrml[1] = ycord-stepSize;
	  nrml[2] = x-xx;
	  tvect[1] = x-xt;
	  tvect[2] = ycord-yt;
	  dprod[k] = nrml[1]*tvect[1] + nrml[2]*tvect[2];
	}
	for (l=0;l<N_VertInCell;l++)
	{
	  if (dprod[l] <= 0)
	  {
	    vncn[vn]=i;
	    cell->SetRegRefine();
	    cell->Refine(RefLevel);
	    vn++;
	  }
	}
      } //*/


/// Type 2 End
     
/// This part below is written to prevent hanging nodes

/// Part 1: uncomment this part 1 and comment next part 2 to refine one layer up, i.e, after regular initial refinement again do regular refinement for the vertex neighbours and then do 1refine
/*
     for(i=0;i<vn;i++)
     {
 	 cell = coll->GetCell(vncn[i]);
	 for(j=0;j<N_VertInCell;j++)
	 {
	   CurrVert = cell->GetVertex(j);
	   N = VertexNumbers[vncn[i]*N_VertInCell + j]*MaxCpV ;
	   N_CellsInThisVert =  PointNeighb[N];
	   for(ii=1;ii<=N_CellsInThisVert;ii++)
	   {
	     VertexCellNo = PointNeighb[N + ii];
	     Vertexcell = coll->GetCell(VertexCellNo);
	     if (!Vertexcell->ExistChildren())
	     {
	       vertexneibcellno[vn1] = VertexCellNo;
	       Vertexcell->SetRegRefine();
	       Vertexcell->Refine(RefLevel);
	       vn1++;
	     } //end if (!Vertexcell->ExistChildren()) 
	   } //end for(ii=1;ii<=N_CellsInThisVert;ii++) 
	 } //end for(j=0;j<N_VertInCell;j++) 
     } //end for(i=0;i<N_Cells;i++) 
          
     
     for(l=0;l<vn1;l++)
      {
        vcell = coll->GetCell(vertexneibcellno[l]);
        N_JointsInCell = vcell->GetN_Joints();
        for(jj=0;jj<N_JointsInCell;jj++)
        {
	  
	  Joint = vcell->GetJoint(jj);
 	  JointType jointtype = Joint->GetType();
 	  if(Joint->InnerJoint())
 	  {
 	    neib_cell = Joint->GetNeighbour(vcell);
 	    if (!neib_cell->ExistChildren())
 	    {
 	      for(kk=0;kk<N_JointsInCell;kk++)
 	      {
    	        if( neib_cell->GetJoint(kk) ==  Joint)
 	        {
 		  neib_cell->Set1Refine(kk);
		  cout<< "neibnosss " << l << endl;
    		  neib_cell->Refine(RefLevel);
		} //end if( neib_cell->GetJoint(kk) ==  Joint)
 	      } //end for(kk=0;kk<N_JointsInCell;kk++)
 	    } //end if (!neib_cell->ExistChildren())
	  } //end if(Joint->InnerJoint())
	} //end for(jj=0;jj<N_JointsInCell;jj++)
      } //end for(l=0;l<vn1;l++) */
///Part 2: 

     for(l=0;l<vn;l++)
      {
        vcell = coll->GetCell(vncn[l]);
        N_JointsInCell = vcell->GetN_Joints();
        for(jj=0;jj<N_JointsInCell;jj++)
        {
	  
	  Joint = vcell->GetJoint(jj);
 	  //JointType jointtype = Joint->GetType();
 	  if(Joint->InnerJoint())
 	  {
 	    neib_cell = Joint->GetNeighbour(vcell);
 	    if (!neib_cell->ExistChildren())
 	    {
 	      for(kk=0;kk<N_JointsInCell;kk++)
 	      {
    	        if( neib_cell->GetJoint(kk) ==  Joint)
 	        {
 		  neib_cell->Set1Refine(kk);
		  neib_cell->Refine(RefLevel);
		} //end if( neib_cell->GetJoint(kk) ==  Joint)
 	      } //end for(kk=0;kk<N_JointsInCell;kk++)
 	    } //end if (!neib_cell->ExistChildren())
	  } //end if(Joint->InnerJoint())
	} //end for(jj=0;jj<N_JointsInCell;jj++)
      } //end for(l=0;l<vn1;l++) */
     

/// END     
      
/// edgewise refinement for layer start (Optional)

/*        
	  
	  rr=0;
	  int aaa;
          for(l=0;l<qq;l++)
          {
  
	  vcell = coll->GetCell(vncn[l]);
	  aaa = vcell->GetGlobalCellNo();
	  cout<<"aaa"<<aaa<<endl;
	  
	  N_JointsInCell = vcell->GetN_Joints();
	  cout<<"aaa"<<aaa<<endl;
	      for(jj=0;jj<N_JointsInCell;jj++)
              {
		Joint = vcell->GetJoint(jj);
		JointType jointtype = Joint->GetType();
		if(Joint->InnerJoint())
		{
		  neib_cell = Joint->GetNeighbour(vcell);
		  if (!neib_cell->ExistChildren())
		  {
		    for(kk=0;kk<N_JointsInCell;kk++)
		    {
		      if( neib_cell->GetJoint(kk) ==  Joint)
		      {
			
			vncn[rr] = neib_cell->GetGlobalCellNo();
			
			//cout<<"vncn"<<vncn[rr]<<endl;
			//cout<<"rr"<<rr<<endl;
			neib_cell->SetRegRefine();
			cout<<"aaa"<<aaa<<endl;
			
			neib_cell->Refine(RefLevel);
			//cout<<"aaa"<<aaa<<endl;
			rr++;
		      } //end if( neib_cell->GetJoint(kk) ==  Joint)
		    } //end for(kk=0;kk<N_JointsInCell;kk++)
		  } //end if (!neib_cell->ExistChildren())
		}
	      }
	  }
	  
// 	  cout<<"rr"<<rr<<endl;
  	  int oo, tt=0;
   	  for(l=0;l<rr;l++)
             {
    
  	      vvcell = coll->GetCell(vncn[l]);
 	      N_JointsInCell2 = vvcell->GetN_Joints();
	      cout<<"N_JointsInCell"<<N_JointsInCell2<<endl;
  	      for(oo=0;oo<N_JointsInCell;oo++)
                {
  		Joint2 = vvcell->GetJoint(oo);
  		jointtype2 = Joint2->GetType();
 		
 
  		if(Joint2->InnerJoint())
//  		//if (jointtype2 != BoundaryEdge || jointtype2 != InterfaceJoint)
   		{
		  //cout<<"rr"<<rr<<endl;
   		  neib_cell2 = Joint2->GetNeighbour(vvcell);
   		  if (!neib_cell2->ExistChildren())
   		  {
   		    for(kk=0;kk<N_JointsInCell;kk++)
   		    {
   		      if( neib_cell2->GetJoint(kk) ==  Joint2)
   		      {
   			vertexneibcellno[tt] = neib_cell2->GetGlobalCellNo();
   			//vncn[tt] = neib_cell2->GetGlobalCellNo();
			cout<<"tt"<<tt<<endl;
   			neib_cell2->SetRegRefine();
   			neib_cell2->Refine(RefLevel);
   			tt++;
   		      } //end if( neib_cell->GetJoint(kk) ==  Joint)
   		    } //end for(kk=0;kk<N_JointsInCell;kk++)
   		  } //end if (!neib_cell->ExistChildren())
   		}
  	      }
   	  }
	  cout<<"tttt"<<tt<<endl;
	  for(l=0;l<tt;l++)
             {
    
  	      vvcell = coll->GetCell(vertexneibcellno[l]);
 	      N_JointsInCell2 = vvcell->GetN_Joints();
	      //cout<<"N_JointsInCell"<<N_JointsInCell2<<endl;
  	      for(oo=0;oo<N_JointsInCell;oo++)
                {
  		Joint2 = vvcell->GetJoint(oo);
  		jointtype2 = Joint2->GetType();
 		
 
  		if(Joint2->InnerJoint())
//  		//if (jointtype2 != BoundaryEdge || jointtype2 != InterfaceJoint)
   		{
		  //cout<<"rr"<<rr<<endl;
   		  neib_cell2 = Joint2->GetNeighbour(vvcell);
   		  if (!neib_cell2->ExistChildren())
   		  {
   		    for(kk=0;kk<N_JointsInCell;kk++)
   		    {
   		      if( neib_cell2->GetJoint(kk) ==  Joint2)
   		      {
   			//vertexneibcellno[tt] = neib_cell2->GetGlobalCellNo();
   			//vncn[tt] = neib_cell2->GetGlobalCellNo();
			cout<<"tttt"<<tt<<endl;
   			neib_cell2->Set1Refine(kk);
   			neib_cell2->Refine(RefLevel);
   			//tt++;
   		      } //end if( neib_cell->GetJoint(kk) ==  Joint)
   		    } //end for(kk=0;kk<N_JointsInCell;kk++)
   		  } //end if (!neib_cell->ExistChildren())
   		}
  	      }
   	  } // */
/// edgewise refinement for layer end 

/// Type 3 (TO MODIFY): Generalized product of normal and tangent vector algorithm
      /*
     
      /// for cylinder temperature problem
      m=0;
      for(i=0;i<N_Cells;i++)
      {
	cell = coll->GetCell(i);
        N_JointsInCell = cell->GetN_Joints();
	
        for(jj=0;jj<N_JointsInCell;jj++)
        {
	  Joint = cell->GetJoint(jj);
 	  JointType jointtype = Joint->GetType();
	  
  	  if(jointtype == IsoBoundEdge) // && jointtype != InterfaceJoint)
 	  {
	   
 	    //cout<<"aaaa"<<x<<endl;
 	    isojoint = (TIsoBoundEdge *)Joint;
 	    //Vertices = ((TIsoBoundEdge *)joint)->GetVertices()
 	    //cout<<"h"<<h<<endl;
 	    IsoVertices = isojoint->GetVertices();
// 	    h=2;
//  	   for(k=0;k<h;k++)
// 	   {
 	      cell->GetVertex(jj)->GetCoords(xc[2], yc[2]);
 	      //IsoVertices->GetCoords(xc[2], yc[2]);
  	      cout<<"xc"<<xc[2]<<endl;
  	      cout<<"yc"<<yc[2]<<endl;
 	      //rad = sqrt(pow((xc-0.2),2) + pow((yc-0.2),2));
  	      //cout<<"rad"<<rad<<endl;
  	      novx[m] = xc[1];
 	      novy[m] = yc[1];
 	      vcn[m] = i;
 	      m=m+1;
 	      novx[m] = xc[2];
 	      novy[m] = yc[2];
  	      vcn[m] = i;
	      m++;
 	    //}
	    
	  }
	}
      }
      
      cout<<"meq"<<m<<endl;
      //exit(0);
      
     
      
      for(i=0;i<N_Cells;i++)
      {
 	cout<<"i"<<x<<endl;
	
	cell = coll->GetCell(i);
	for(k=0;k<N_VertInCell;k++)
	{
	  cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k]);
	  xt = xcoords[k];
	  yt = ycoords[k];
	  ff=0;
	  for (r=0;r<m;r++)
	  {
	    D[r] = sqrt(pow(novx[r]-xt, 2) + pow(novy[r]-yt, 2));
 	    //cout<<"D"<<D[r]<<endl;
	    ff++;
	    //exit(0);
	  }
	  //exit(0);
	  //cout<<"ff"<<ff<<endl;
	  //exit(0);

	  small=D[0];
	  for (kk=1; kk<m;kk++)
	  {
	    if (small < D[kk])
	    {
	      small = D[kk];
	      //isocell = coll->GetCell(kk);
  // 	      xcord=novx[kk];
  // 	      ycord=novy[kk];
	    }
	  }
	  //cout<<"ll"<<ll<<endl;
	  for (ll=0;ll<m;ll++)
	  {
	    if (small==D[ll])
	    {
	      
	      isocell = coll->GetCell(vcn[ll]);
	      xcord=novx[ll];
	      ycord=novy[ll];
	      xxcord=novx[ll+1];
	      yycord=novy[ll+1];
	      cout<<"xcord"<<xcord<<endl;
	      cout<<"ycord"<<ycord<<endl;
	      cout<<"xxcord"<<xxcord<<endl;
	      cout<<"yycord"<<yycord<<endl;
	    }
	  }
// 	  cout<<"ll"<<ll<<endl;
// 	  cout<<"xcord"<<xcord<<endl;
// 	  cout<<"iiiiiii"<<x<<endl;
	
	  //exit(0);
	  
	  /// get the other cordinate in the iso edge to compute normal vector n tangent vectors
	  for(s=0;s<N_JointsInCell;s++)
	  {
	    Joint2 = isocell->GetJoint(s);
	    jointtype2 = Joint2->GetType();
// 	    cout<<"iiiiiii"<<x<<endl;
	    if(jointtype2 == IsoBoundEdge)
	    {
	    
	      isojoint2 = (TIsoBoundEdge *)Joint2;
	      //h = isojoint2->GetN_Vertices();
	      IsoVertices2 = isojoint2->GetVertices();
	      h=2;
	       
	       isocell->GetVertex(s)->GetCoords(xxc[2], yyc[2]);
		//IsoVertices2->GetCoords(novxx[h], novyy[h]);
 		novxx[1] = xxc[1];
 		novyy[1] = yyc[1];
		novxx[2] = xxc[2];
 		novyy[2] = yyc[2];
	      
	      for(q=0;q<h;q++)
	      {
		if(novxx[h]!=xcord && novyy[h]!=ycord)
		{
		  xxcord=novxx[h];
		  yycord=novyy[h];
		}
 	      }
	    }
	  }
	  nrml[1] = ycord-yycord;
	  nrml[2] = xcord-xxcord;
	  tvect[1] = xcord-xt;
	  tvect[2] = ycord-yt;
	  dprod[k] = nrml[1]*tvect[1] + nrml[2]*tvect[2];
	  //cout<<"dprod"<<dprod[k]<<endl;
	  if (dprod[k] < 0.0)
	    cout<<"dprod"<<dprod[1]<<endl;
	    
	}
	//exit(0);
  	for (l=0;l<N_VertInCell;l++)
  	{
	  //DP = dprod[1]*dprod[2]*dprod[3];
 	  if (dprod[l] <= 0.0)
 	  {
 	    vncn[vn]=i;
 	    cell->SetRegRefine();
 	    cell->Refine(RefLevel);
 	    vn++;
 	  }
  	}
      }*/



     return 0;
     
     
     
     
     
     
     

}


#endif // #ifdef __2D__
