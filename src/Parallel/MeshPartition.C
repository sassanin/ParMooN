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
// @(#)MeshPartition.C
// 
// Purpose:     partition the domain into "npart" parts for parallel computing
// 
// Author:      Sashikumaar Ganesan
// History:      start of implementation  07/09/09 (Sashikumaar Ganesan)
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <Database.h>
#include <Domain.h>
#include <Output2D.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
// #include <malloc.h>
#include <Joint.h>
#include <SubDomainJoint.h>
#include <SubDomainHaloJoint.h>
#include <Vertex.h>
#include <BaseCell.h>
#include <Quadrangle.h>
#include <MacroCell.h>
#include <Edge.h>

#ifdef __2D__
  #include <FEDatabase2D.h>
#endif

#ifdef _MPI
extern "C"
{
  #include <metis.h>
  #include <parmetis.h>
}
#endif


static void Sort(int *Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  int Mid, Temp;
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
  int i, j, k, *rr, len, s;
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



static void Sort(TEdge **Array, int length)
{
  int n=0, l=0, r=length-1, m, max=0;
  int i, j, k, *rr, len, s;
  TEdge *Mid, *Temp;
  double lend = length;

  len=(int)(4*log(lend)/log(2.0) +2);
  

  
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
	if(max<n) max=n;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

    
//   printf("%d :Time taken for Domain Decomposition is %d\n", max, len);
  
  delete [] rr;

}



static void Sort(TBaseCell **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  TBaseCell *Mid, *Temp;
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

#ifdef _MPI

#ifdef __2D__


#else // 3D

int Partition_Mesh3D(MPI_Comm comm, TDomain *Domain, int &MaxRankPerV)
{
 idx_t *Cell_Rank, *Vert_Rank, ne, nn, *eptr, nparts, objval, edgecut=0, ncommon=3; 

 
 int i, j, k, ll, m, M, MM, n, rank, size, N_Cells, N, NN, ID, Neib_ID, out_rank= TDatabase::ParamDB->Par_P0;
 int *GlobalCellIndex, N_RootVertices, N_VertInCell, N_JointsInCell, N_AllLocVert;
 int *VertexNumbers, *PointNeighb, maptype, MaxCpV;
 int MaxLen, N_Edges, N_FaceEdges, ii, jj, GblCellNr;
 int N_EdgeDel=0, N_VertexDel=0, N_CellDel=0, N_CellsInHaloVert;
 const int *TmpFV, *TmpLen, *TmpVE, *EdgeVertex, *NeibEdgeVertex;
 
 int VertexCellNo, VertexcellID, N_CellsInThisVert, *VertNeibRanks;
 int N_LocalCells, N_OwnCells, N_HalloCells, N_OwnIncidentCells, N_NeibIncidentCells;
 int a, b, N_SubDomInThisVert, *HaloCellIndex;
 int M1, M2, N1, N2, N_CellsIn_a, N_CellsIn_b;
 int N_CrossNeibs, *CrossNeibsRank, *HaloCellGlobalNo, *HaloCellLocVertNo;
 int N_SubDomIn_a, N_SubDomIn_b, *Temp, kk, GlobCellNo, test_b, EdgeCellID;

 bool UPDATE, UPDATE_1;

 TVertex *Vert_a, **VertexDel, *Last;
 TCollection *coll;
 TBaseCell *cell, *neib_cell, *Vertexcell, **SubDomainCells, *cell_a, *cell_b, **IncCells, **IncHalloCells, *HalloCell; 
 TBaseCell **CellDel, *LastCell;
 TJoint *Joint, *NewJoint;
 TEdge *edge, **EdgeDel, *LastEdge;
 TShapeDesc *ShapeDesc, *NeibShapeDesc;
 TVertex *CurrVert, *NeibCurrVert, *NeibVert_a; 
 
 MPI_Status status, status1, status2, status3;
 MPI_Request request, request1;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // root will not take part in computations
  MaxCpV = 0;
  MaxRankPerV = -1;

//   if(rank==0)
//    MaxRankPerV = 1;
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  cell = coll->GetCell(0);
  N_VertInCell = cell->GetN_Vertices();
  N_JointsInCell = cell->GetN_Joints();
  N_AllLocVert = N_VertInCell*N_Cells;
  VertexNumbers= new int[N_AllLocVert];
  Vert_Rank = new idx_t[N_AllLocVert];  
  Cell_Rank = new idx_t[N_Cells];  
  eptr = new idx_t[N_Cells+1];  
  
 // if(rank==0)
 //  printf("Number of  ranks: %d\n",  size);

  int global_flag,flag = 0;
  if(N_Cells<size)
   {
//      printf("Number of cells less than number of processors !!!, %d\n",  N_Cells);
     flag = 1;
//      MPI_Finalize();
//      exit(0);
   }

   MPI_Allreduce(&flag, &global_flag, 1, MPI_INT, MPI_SUM, comm);
     
   if(global_flag != 0)
     {
       if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("Warning :: some ranks have Number of cells less than number of processors !!!. \n");
	  printf("metis type being changed from %d to %d\n",TDatabase::ParamDB->Par_P2,abs(TDatabase::ParamDB->Par_P2-1));
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
       
       TDatabase::ParamDB->Par_P2 = abs(TDatabase::ParamDB->Par_P2 - 1);
       return 1;  		//partinioning unsuccessful
     }
   
 if(size==1)
   {
    cout <<  "Total number of process should be grater than 1 (or 2 if root is not involved in computation), but you have " << size<<endl;
     MPI_Finalize();
     exit(0);
   }

 //check, all cell have same number of vertices  
  for(i=1;i<N_Cells;i++)
   if(((coll->GetCell(i))->GetN_Vertices() != N_VertInCell) || ( (coll->GetCell(i))->GetN_Joints() != N_JointsInCell) )
    {
     cout << N_JointsInCell << "Mesh partition for heterogeneous cells are not yet implemented " << N_VertInCell<<endl;
     MPI_Finalize();
     exit(0);
    }


// partition the mesh in the root and send info to all processors
 if(rank==0)
  {
  // cout <<  "Total Cells " << N_Cells<<endl;

   int  m1,  *NumberVertex;
//    int numflag=0; // c-style numbering (array starting with 0)
   int type=TDatabase::ParamDB->Par_P2;
   double t1, t2;
   
   
   idx_t *MetisVertexNumbers; 
   idx_t options[METIS_NOPTIONS];
 
   METIS_SetDefaultOptions(options);
   
   options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; 
   options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; 
   options[METIS_OPTION_NUMBERING] = 0; 
   options[METIS_OPTION_CONTIG] = 1;

   TVertex **Vertices;    

   
    if(!(N_VertInCell==4 || N_VertInCell ==8))
     {
      cout<<" Error only  Tetra or Hexa mesh can be partitioned !!" <<endl;
      MPI_Finalize();
      exit(0);
     }
//    cout<<"N_VertInCell::"<<N_VertInCell<<endl;
//    exit(0);
   eptr[0] = 0;  
   for(i=1;i<=N_Cells;i++)
    eptr[i] = eptr[i-1] + N_VertInCell;
   

    MetisVertexNumbers= new idx_t[N_AllLocVert];
    NumberVertex=new int[N_AllLocVert];
    Vertices=new TVertex*[N_AllLocVert];
   
  /** *********************************************/
  /** STEP 1 : STORE VERTEX POINTERS CELL-WISE */
  /** *********************************************/
    N = 0;    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      for(j=0;j<N_VertInCell;j++)
       {
        Vertices[N]=cell->GetVertex(j);
        N++;
       }
     } // for i=0

  /** *********************************************/
  /** STEP 2 : SORT THE VERTICES ARRAY */ 
  /** *********************************************/
    Sort(Vertices, N);
    
  /** ***************************************************/
  /**STEP 3: STORE THE SORTED POINTER ARRAY AS INDICES */
  /** ***************************************************/
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


  /** *********************************************/   
  /** STEP 4 : STORE THE INDICES CELL-WISE */
  /** *********************************************/ 
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

   memcpy(MetisVertexNumbers, VertexNumbers, N_AllLocVert*SizeOfInt);

 //  cout << "Total Vertices " <<N_RootVertices<<endl;
   //cout << "Vertices N " <<N_AllLocVert<<endl;
   
  /** 1. SEND INFO TO ALL PROCESSORS */
  for(i=1;i<size;i++)
   MPI_Send(&N_RootVertices, 1, MPI_INT, i, 75, comm);
  
  //   for (i=0; i<N_LocVertices; i++) 
 //    cout<< i <<"  " << "elements " << MetisVertexNumbers[i] <<endl;
     
  /** 2. SEND INFO TO ALL PROCESSORS */    
  for(i=1;i<size;i++)
    MPI_Send(VertexNumbers, N_AllLocVert, MPI_INT, i, 80, comm);  
  
  /** *********************************************/
  /** STEP 5 : MESH PARTITION */
  /** *********************************************/
  
  ne = (idx_t)N_Cells;
  nn =  (idx_t)N_RootVertices;
  nparts = (idx_t)size;
  
  t1 = MPI_Wtime();
   if(type == 0)
    {
     METIS_PartMeshNodal(&ne, &nn, eptr, MetisVertexNumbers, NULL, NULL, &nparts, NULL, options, &edgecut, Cell_Rank, Vert_Rank);
    }
   else if(type == 1)
    {
     METIS_PartMeshDual(&ne, &nn,  eptr, MetisVertexNumbers,  NULL, NULL,  &ncommon,
                       &nparts, NULL, options, &edgecut, Cell_Rank, Vert_Rank);
    }
   else
    {
     cout<<" Error METIS_PartMesh implemented for Par_P2 = 0 or 1 !!" <<endl;
     MPI_Abort(comm, 0);
    }

  t2 = MPI_Wtime();
  OutPut( "Time taken for METIS mesh partinioning "<< t2-t1<< " sec"<<endl); 


  /** *********************************************/
  /** STEP 6 : TO FIND MAXIMUM CELLS PER VERTEX */
  /** *********************************************/ 
   PointNeighb = new int[N_RootVertices];
   memset(PointNeighb, 0, N_RootVertices*SizeOfInt);

  
  for (i=0;i<N_AllLocVert;i++)
    PointNeighb[VertexNumbers[i]]++;
 
   // find maximum cells per vertex ( maxCpV )
   for (i=0;i<N_RootVertices;i++)
    if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];

//       printf("Max number of cells per vertex %d \n", MaxCpV);
   MaxCpV++;
  for(i=1;i<size;i++)
   MPI_Send(&MaxCpV, 1, MPI_INT, i, 85, comm);

   delete [] PointNeighb;
  
  /** ********************************************************************************************/
  /** STEP 7 : CREATE AN ARRAY CONTAINING INFO REGARDING THE NEIGHBOURS OF VERTICES (CELL INDEX) */
  /** PointNeighb's first column contains number of neib cells associated with each vertex*/
  /** further columns contain the cell numbers associated with this vertex*/
  /** ********************************************************************************************/
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

   for(i=1;i<size;i++)
    MPI_Send(PointNeighb, N_RootVertices*MaxCpV, MPI_INT, i, 90, comm); 

   delete [] MetisVertexNumbers;
   delete [] NumberVertex;
   delete [] Vertices;
   
  }
 else
  {
    MPI_Recv(&N_RootVertices, 1, MPI_INT, 0, 75, comm, &status);
    // printf("%d SubDomain_N_Cells in rank test 1, %d  \n", rank, N_RootVertices );

    MPI_Recv(VertexNumbers, N_AllLocVert, MPI_INT, 0, 80, comm, &status);

    MPI_Recv(&MaxCpV, 1, MPI_INT, 0, 85, comm, &status);

    PointNeighb = new int[N_RootVertices*MaxCpV];

    MPI_Recv(PointNeighb, N_RootVertices*MaxCpV, MPI_INT, 0, 90, comm, &status);     
  } //else if(rank==0)


    /**Metis partition done!! code for all processors */  
    MPI_Bcast(Cell_Rank, N_Cells, MPI_INT, 0, comm);
    MPI_Bcast(Vert_Rank, N_RootVertices, MPI_INT, 0, comm);
  
  /** ********************************************************/
  /** STEP 8.1 : SETTING SUBDOMAIN NUMBER AND GLOBAL CELL NO */
  /** ********************************************************/
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      cell->SetSubDomainNo(Cell_Rank[i]);
      cell->SetGlobalCellNo(i); 
     }

  /** *********************************************/
  /** STEP 8.2 : SET ALL VERTICES in OWN cells TO -1 */
  /** *********************************************/ 
    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      if(cell->GetSubDomainNo()==rank) 
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

  /** *********************************************/
  /** STEP 9 : Fill the VertNeibRanks info        */
  /** first column contains how many ranks contain this vertex further columns contain the rank ID of the subdomains */
  /** *********************************************/ 
    VertNeibRanks = new int[N_RootVertices*MaxCpV];
    memset(VertNeibRanks, 0, N_RootVertices*MaxCpV*SizeOfInt);

    HaloCellIndex = new int[size];
    for(i=0;i<size;i++)
      HaloCellIndex[i] = -1;

    HaloCellGlobalNo = new int[MaxCpV];
    HaloCellLocVertNo = new int[MaxCpV];

    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      ID = cell->GetSubDomainNo();

      /**run only through own cells */
      if(ID==rank) 
       {
        cell->SetAsOwnCell();

        /** set SubDomainVert if any vert in this cell is so needed for moving meshes and setting cross vertex */
        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);
	  
	  if(CurrVert->GetClipBoard() != -1)
           continue;

          CurrVert->SetClipBoard(5);
          N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;

          N_CellsInThisVert =  PointNeighb[N];
          N_OwnIncidentCells = 0;
          N_NeibIncidentCells = 0;

          // check! any subdomain cell containg this vert has to be added
          // Setting all halo cells
          UPDATE = TRUE;
          for(ii=1;ii<=N_CellsInThisVert;ii++)
           {
            VertexCellNo = PointNeighb[N + ii];
            Vertexcell = coll->GetCell(VertexCellNo);
            VertexcellID = Vertexcell->GetSubDomainNo();

            /** Hallo Cell */    
            if(VertexcellID != rank) 
             {
              //vertex belong to diff subdomain 
              if(UPDATE)
               {
                cell->SetAsDependentCell();
                CurrVert->SetAsSubDomainVert();
                HaloCellIndex[VertexcellID]= VertexCellNo;
                UPDATE = FALSE;
               }
              Vertexcell->SetAsHaloCell();

               /** changed on 4 Feb 2012 - by Sashi */
               if(HaloCellIndex[VertexcellID] == -1)
                HaloCellIndex[VertexcellID]= VertexCellNo;
    
             }
            } //  for(ii=1;ii< 

           if(CurrVert->IsSubDomainVert())
             for(ii=1;ii<=N_CellsInThisVert;ii++)
              {
               VertexCellNo = PointNeighb[N + ii];
               Vertexcell = coll->GetCell(VertexCellNo);
               VertexcellID = Vertexcell->GetSubDomainNo();

               // all cells associated with this vert are DepCells
               Vertexcell->SetAsDependentCell();

               if(VertexcellID==rank)
                continue;

               N_CellsIn_b = VertNeibRanks[N];
               UPDATE = TRUE;

               for(jj=1;jj<=N_CellsIn_b;jj++)
                if( VertNeibRanks[N + jj]==VertexcellID)
                 {
                  UPDATE = FALSE;
                  break;
                 }

               if(UPDATE)
                {
                 N2 = HaloCellIndex[VertexcellID];
                 HaloCellGlobalNo[VertNeibRanks[N]] = N2;

                 //find the local index of this vertex in the neib cell
                 jj=0;
                 while(CurrVert != (coll->GetCell(N2))->GetVertex(jj)) jj++;
                 HaloCellLocVertNo[VertNeibRanks[N]] = jj;

                 VertNeibRanks[N]++;
                 VertNeibRanks[N + VertNeibRanks[N]] = VertexcellID;

                 if(MaxRankPerV<VertNeibRanks[N])
                   MaxRankPerV=VertNeibRanks[N];
                } //   if(UPDATE)
              } //  for(ii=1;ii<=N_CellsInT

           N_SubDomIn_a = VertNeibRanks[N];
           Temp = VertNeibRanks+(N+1);
 
    
           /** set only one cell (lowest index cell) from each subdomain as the Halo cell for this vertex */
           CurrVert->SetSubDomainInfo(N_SubDomIn_a,  Temp, HaloCellGlobalNo, HaloCellLocVertNo);

            //reset
            for(jj=0;jj<size;jj++)
             HaloCellIndex[jj] = -1;

           } // for(j=0;j<N_V
       } // if(ID==ra
     } //  for(i=0;i<N_Cel 
    
    MaxRankPerV++;    // including own rank
     
    delete [] HaloCellIndex;
    delete [] HaloCellGlobalNo;
    delete [] HaloCellLocVertNo;

  /** *************************************************/
  /** STEP 10 : Get own SubDomain collection of cells */
  /** *************************************************/ 
     N_LocalCells = 0;
     N_OwnCells = 0;

     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       if(ID==rank || cell->IsHaloCell() )
        N_LocalCells++;

       if(ID==rank)
        N_OwnCells++;
      } //  for (i=0;i<N_Cells;i++)

      N_HalloCells = N_LocalCells - N_OwnCells;
      
   //  printf("Rank: %d N_own Cells %d\n ", rank,  N_OwnCells);   
      //      int aa;
//      for(aa=0;aa<size;aa++){
//      if(rank==aa)
//       printf("Rank:: %d N_own Cells:: %d  N_cells:: %d\n ", rank,  N_OwnCells, N_LocalCells);
//      //MPI_Barrier(comm);
//      }
//      MPI_Barrier(comm);
     //exit(0);
 
     global_flag=0;
     flag = 0;
     if(N_OwnCells == 0){
       flag = 1;
     }
     
     MPI_Allreduce(&flag, &global_flag, 1, MPI_INT, MPI_SUM, comm);
     
     if(global_flag != 0)
     {
       if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("Warning :: some ranks have not been allocated any cells. \n");
	  printf("metis type being changed from %d to %d\n",TDatabase::ParamDB->Par_P2,abs(TDatabase::ParamDB->Par_P2-1));
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
       
       TDatabase::ParamDB->Par_P2 = abs(TDatabase::ParamDB->Par_P2 - 1);
       return 1;  		//partinioning unsuccessful
     }
      
     if(N_LocalCells)
      {
       // collect the own collection
       SubDomainCells = new TBaseCell*[N_LocalCells];
       GlobalCellIndex = new int[N_LocalCells];
      }
     else
      {
       SubDomainCells = NULL;
       GlobalCellIndex = NULL; 
      }


      
  /** *************************************************************/
  /** STEP 11 : Fill subdomain own cells and their global numbers */
  /** *************************************************************/ 
  VertexDel = new TVertex*[N_AllLocVert];
  EdgeDel = new TEdge*[N_Cells*12]; // max 12 in Hexa
  CellDel = new TBaseCell*[N_Cells];
     N=0;
     M=0;
     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       // in the collection, first own cells then Halo cells
       if(ID==rank)
        {
         SubDomainCells[N] =  cell;
         GlobalCellIndex[N] =  i;
         N++;
        }
      else if(cell->IsHaloCell() )
        {
         SubDomainCells[N_OwnCells + M] =  cell;
         GlobalCellIndex[N_OwnCells + M] =  i;
         M++;
        }
      else //cell no longer needed for this rank, delete!!!!!
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
          EdgeDel[N_EdgeDel++] = cell->GetEdge(j);

        for(j=0;j<N_VertInCell;j++)
         VertexDel[N_VertexDel++] = cell->GetVertex(j);

         CellDel[N_CellDel++] = cell;
       }
      } //  for (i=0;i<N_Cells;i++)

  /** ***********************************************************************************/
  /** STEP 12 : Delete other processors cells from own collection, excluding halo cells */
  /** ***********************************************************************************/ 
    // 
/*    if(N_VertexDel)
     Sort(VertexDel, N_VertexDel); 
    if(N_EdgeDel)
     Sort(EdgeDel, N_EdgeDel);*/
//     if(N_CellDel)
//      Sort(CellDel, N_CellDel);

  /*  Last=NULL;
    for(i=0;i<N_VertexDel;i++)
     {
      if( VertexDel[i] != Last)
       {
        delete VertexDel[i];
        Last = VertexDel[i];
       }
     }
   delete [] VertexDel;


    LastEdge=NULL;
    for(i=0;i<N_EdgeDel;i++)
     {
      if( EdgeDel[i] != LastEdge)
       {
        delete EdgeDel[i];
        LastEdge = EdgeDel[i];
       }
     }
   delete [] EdgeDel;

    LastCell=NULL;
    for(i=0;i<N_CellDel;i++)
     {
      if( CellDel[i] != LastCell)
       {
        delete CellDel[i];
        LastCell = CellDel[i];
       }
     }
   delete [] CellDel;

    if((N+M)!=N_LocalCells)
     {
      printf("Error in mesh partition %d \n", rank);
      MPI_Abort(comm, 0);
     } 
     
   // ********************************************************************************** */
   // STEP 13 : Edges will be flaged as SubDomainEdges and subdomain infos are filled     */
   // *********************************************************************************** */ 
  
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);
          if(edge->GetClipBoard()==-1)
           {
            edge->SetClipBoard(5);
            edge->InitSubDomainInfo(rank);
           } // if(edge->GetCli
         } // for(j=0;j<N_Edges;
       } // if(cell->IsDependentCell
     } // for(i=0;i<N_Cells;i++

              
  
 /** change all subdomain joint from JointEqN into TSubDomainJoint */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      /**run only through Dependent cells */
      if(cell->IsDependentCell()) 
       {
        GblCellNr = GlobalCellIndex[i];
        ID = cell->GetSubDomainNo();

	//if( GblCellNr == 111615)
          // printf("rank %d local cell no %d\n", rank, i);
	
        /** change the faces, also edges in 3D */
        for(j=0;j<N_JointsInCell;j++)
         {
          Joint = cell->GetJoint(j);

          if(Joint->GetType() == JointEqN)
           {
            neib_cell = Joint->GetNeighbour(cell);
            Neib_ID = neib_cell->GetSubDomainNo();


            if(ID!=Neib_ID)
             {
              GlobCellNo =  neib_cell->GetGlobalCellNo();

              // this joint belongs to two SubDomains
              cell->SetAsSubDomainInterfaceCell();

              // find neib cell local face number
              m=0;
              while(Joint!=neib_cell->GetJoint(m)) m++;

              delete Joint;
              NewJoint = new TSubDomainJoint(cell, neib_cell, Neib_ID, GlobCellNo, m);
              cell->SetJoint(j, NewJoint);
              neib_cell->SetJoint(m, NewJoint);
              NewJoint->SetMapType();

              // set all edges in this face as SubDomainEdges
              ShapeDesc = cell->GetShapeDesc();
              ShapeDesc->GetFaceEdge(TmpFV, TmpLen, MaxLen);
              N_FaceEdges = TmpLen[j];
              for(n=0;n<N_FaceEdges;n++)
               {
                edge = cell->GetEdge(TmpFV[j*MaxLen+n]);
                edge->SetAsNotCrossEdgeFor(rank, Neib_ID);
               } // for(n=0;n<N_Edge

              // remove the Neib_ID from VertNeibRanks
              // these face vertices cannot be cross vertices for processors ID and Neib_ID
              ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
              for(ii=0;ii<TmpLen[j];ii++)
              {
               jj = TmpFV[j*MaxLen + ii];
               M = VertexNumbers[GblCellNr*N_VertInCell + jj] *MaxCpV;
               N = VertNeibRanks[M];

               for(jj=1;jj<=N;jj++)
                if(VertNeibRanks[M + jj] == Neib_ID )
                 {
                  VertNeibRanks[M + jj] = VertNeibRanks[M + N];
                  VertNeibRanks[M]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
              } // for(ii=0;ii<TmpL
             }//if(ID!=Neib_ID
           } // if(Joint->GetT
         }// for (j=0;j< ;
       } // if(cell->IsDependentCell() 
     }// for(i=0;i<N_Cel

     
    /** find cross edges (if any) */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];
        N_Edges=cell->GetN_Edges();
        ShapeDesc= cell->GetShapeDesc();
        ShapeDesc->GetEdgeVertex(EdgeVertex);

        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);

          if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
           {
            edge->SetClipBoard(5);
            N_CrossNeibs = edge->GetN_CrossNeibs();

            if(N_CrossNeibs==0)
             continue;

            cell->SetAsCrossEdgeCell();

            a = EdgeVertex[2*j];
            b = EdgeVertex[2*j+ 1];

            Vert_a = cell->GetVertex(a);
            edge->SetCrossNeibInfo(Vert_a);
            edge->GetCrossEdgeNeibs(N_CrossNeibs, CrossNeibsRank);

            M1 = VertexNumbers[M*N_VertInCell + a]*MaxCpV;
            N1 = VertNeibRanks[M1];

            M2 = VertexNumbers[M*N_VertInCell + b]*MaxCpV;
            N2 = VertNeibRanks[M2];

            for(jj=0;jj<N_CrossNeibs;jj++)
             {
              Neib_ID = CrossNeibsRank[jj];

               for(kk=1;kk<=N1;kk++)
                if(VertNeibRanks[M1 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M1 + kk] = VertNeibRanks[M1 + N1];
                  VertNeibRanks[M1]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)

               for(kk=1;kk<=N2;kk++)
                if(VertNeibRanks[M2 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M2 + kk] = VertNeibRanks[M2 + N2];
                  VertNeibRanks[M2]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
             } // for(jj=0;jj<N_CrossNei
           } //   if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
          } // for(j=0;j<N_Edge
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel
         
    /** find cross vertex (if any) */
    // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    for(i=0;i<N_HalloCells;i++)
     {
      SubDomainCells[N_OwnCells + i]->SetClipBoard(-1); 
     }
     
    /** set incident cell list for all vertices */
    IncCells = new TBaseCell*[MaxCpV];
    IncHalloCells = new TBaseCell*[MaxCpV];
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      M = GlobalCellIndex[i];    
     
      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        {
         // set the vertexcells for all vertices in dep. cells
         CurrVert = cell->GetVertex(j);
         if(CurrVert->GetClipBoard() != -1)
          continue;

         CurrVert->SetClipBoard(5);

         N = VertexNumbers[M*N_VertInCell + j] *MaxCpV;
         N_CellsInThisVert =  PointNeighb[N];

         for(ii=1;ii<=N_CellsInThisVert;ii++)
           IncCells[ii-1] = coll->GetCell(PointNeighb[N + ii]);

         CurrVert->SetVertexCells(N_CellsInThisVert, IncCells);


         /** added, hallo cells vertCells neib info on 11 Mar 2015 - Sashi */  
         for(ii=0;ii<N_CellsInThisVert;ii++) 
          {
           /** fill only for Hallo cells, since own Dependentcells' vertices got filled above */    
	   if(IncCells[ii]->IsHaloCell())
	    {
	      /** this Hallo cell is already handled */
	      if(IncCells[ii]->GetClipBoard() != -1)
               continue;

              IncCells[ii]->SetClipBoard(5);

//       if(rank==0 && IncCells[ii]->GetGlobalCellNo()==32)
//       printf("N_CellsInThisVert :%d, Vert_index: %d,   GlobalCellNo: %d , Neib_ID: %d \n",  N_CellsInThisVert, j, IncCells[ii]->GetGlobalCellNo(), IncCells[ii]->GetSubDomainNo()  );  

              for(jj=0;jj<N_VertInCell;jj++)
               {
                ll = 0; 
                NeibCurrVert = IncCells[ii]->GetVertex(jj);

                /** even in Hallo cells, set neib info only for unset vertices */
                if( NeibCurrVert->GetNNeibs()!= 0)
                 continue;

                MM = IncCells[ii]->GetGlobalCellNo();
                NN = VertexNumbers[MM*N_VertInCell + jj] *MaxCpV;
                N_CellsInHaloVert =  PointNeighb[NN];

                for(kk=1;kk<=N_CellsInHaloVert;kk++)
                 {
                  HalloCell = coll->GetCell(PointNeighb[NN + kk]);
 
                  /** include only owna and Halo in the neib list */
                  if(HalloCell->GetSubDomainNo() == rank || HalloCell->IsHaloCell())
                   {    
                    IncHalloCells[ll] = HalloCell;
                    ll++;
                   } // if(HalloCell->GetSubDomainNo(
                  } // for(kk=1;kk<=N_CellsInHaloVert;kk++)
                  
                 NeibCurrVert->SetVertexCells(ll, IncHalloCells);    
                  
                } //   for(jj=0;jj<N_VertInCell;jj++)   
              }  // if(IncCells[ii]->IsHalo
          } //  for(jj=0;jj<N_VertInCell;jj++)
 
  
        } //  for (j=0;j<N_VertInCell;
     } //  for(i=0;i<N_OwnCell
    delete [] IncCells;
    delete [] IncHalloCells;
    
     // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];

        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);

          // continue if the vert is already handled or not a subdomain vert
          if( (CurrVert->GetClipBoard() != -1) || !(CurrVert->IsSubDomainVert()))
           continue;


          CurrVert->SetClipBoard(5);

          M1 = VertexNumbers[M*N_VertInCell + j] *MaxCpV ;
          N1 = VertNeibRanks[M1];

          if(N1!=0)
           cell->SetAsCrossVertexCell();

          for(ii=1;ii<=N1;ii++)
            CurrVert->AddCrossNeib(VertNeibRanks[M1+ii]);

//            TDatabase::ParamDB->Par_P5 =0;
          } // for(j=0;j<N_VertInCell;j
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel


    if(N_LocalCells)
     { Domain->ReplaceTreeInfo(N_LocalCells, SubDomainCells, GlobalCellIndex, N_OwnCells); }
    else
     { Domain->SetN_OwnCells(N_OwnCells); }

        // initialize iterators
    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   delete [] VertNeibRanks;

//    delete [] Vert_Rank; // for more no. processor showing double free or corruption 
   delete [] VertexNumbers;
   delete [] Cell_Rank;
   delete [] PointNeighb;

   
   //Barrier is needed, before calling FECommunicator, since neib process must have all info
   MPI_Barrier(comm);
   
//   if(rank==out_rank)
//      printf("Mesh Partition end rank %d, N_HalloCells %d\n" , rank, N_HalloCells);
// 
//      MPI_Finalize();
//     exit(0);  
   return 0;		//partinioning successful
}
#endif // 3D

void Domain_Crop(MPI_Comm comm, TDomain *Domain)
{
  
#if 1
  //Variable list
  
 idx_t *Cell_Rank, *Vert_Rank;
 int i, j, k, m, ll, M, MM, n, rank, size, N_Cells, N, NN, ID, Neib_ID, out_rank= TDatabase::ParamDB->Par_P0;
 int *GlobalCellIndex, N_RootVertices, N_VertInCell, N_JointsInCell, N_AllLocVert;
 int *VertexNumbers, *PointNeighb, maptype, MaxCpV;
 int MaxLen, N_Edges, N_FaceEdges, ii, jj, GblCellNr;
 int N_EdgeDel=0, N_VertexDel=0, N_CellDel=0, N_CellsInHaloVert;
 const int *TmpFV, *TmpLen, *TmpVE, *EdgeVertex, *NeibEdgeVertex;

 int MaxRankPerV,VertexCellNo, VertexcellID, N_CellsInThisVert, *VertNeibRanks, *NumberVertex;
 int N_LocalCells, N_OwnCells, N_HalloCells, N_OwnIncidentCells, N_NeibIncidentCells;
 int a, b, N_SubDomInThisVert, *HaloCellIndex;
 int M1, M2, N1, N2, N_CellsIn_a, N_CellsIn_b;
 int N_CrossNeibs, *CrossNeibsRank, *HaloCellGlobalNo, *HaloCellLocVertNo;
 int N_SubDomIn_a, N_SubDomIn_b, *Temp, kk, GlobCellNo, test_b, EdgeCellID;

 bool UPDATE, UPDATE_1;

 TVertex *Vert_a, **VertexDel, *Last;
 TCollection *coll;
 TBaseCell *cell, *neib_cell, *Vertexcell, **SubDomainCells,**OwnCells, *cell_a, *cell_b, **IncCells, **IncHalloCells, *HalloCell;
 TBaseCell **CellDel, *LastCell;
 TJoint *Joint, *NewJoint;
 TEdge *edge, **EdgeDel, *LastEdge;
 TShapeDesc *ShapeDesc, *NeibShapeDesc;
 TVertex *CurrVert, *NeibVert_a, *NeibCurrVert;
  
 MPI_Status status, status1, status2, status3;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
 
 // root will not take part in computations
  MaxCpV = 0;
  MaxRankPerV = -1;

  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  cell = coll->GetCell(0);
  N_VertInCell = cell->GetN_Vertices();
  N_JointsInCell = cell->GetN_Joints();
  N_AllLocVert = N_VertInCell*N_Cells;
  VertexNumbers= new int[N_AllLocVert];
  Vert_Rank = new idx_t[N_AllLocVert];
  Cell_Rank = new idx_t[N_Cells];

  VertexDel = new TVertex*[N_AllLocVert];
  EdgeDel = new TEdge*[N_Cells*12]; // max 12 in Hexa
  CellDel = new TBaseCell*[N_Cells];

   TVertex **Vertices;
   
#endif
   
  
   
   int  m1;
   int etype, numflag=0; // c-style numbering (array starting with 0)
   int type=TDatabase::ParamDB->Par_P2;

   double t1, t2;

    NumberVertex=new int[N_AllLocVert];
    Vertices=new TVertex*[N_AllLocVert];

    N = 0;
   /** *********************************************/
    /** STEP 1 : STORE VERTEX POINTERS CELL-WISE */
    /** *********************************************/
    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      for(j=0;j<N_VertInCell;j++)
       {
        Vertices[N]=cell->GetVertex(j);
        N++;
       }
     } 
  /** *********************************************/
   /** STEP 2 : SORT THE VERTICES ARRAY */ 
  /** *********************************************/
    Sort(Vertices, N);

  /** ***************************************************/
    /**STEP 3: STORE THE SORTED POINTER ARRAY AS INDICES */
  /** ***************************************************/
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

 /** *********************************************/   
  /** STEP 4 : STORE THE INDICES CELL-WISE */
 /** *********************************************/ 
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
       } 
     } 


  /** *********************************************/
   /** STEP 6 : TO FIND MAXIMUM CELLS PER VERTEX */
  /** *********************************************/ 
  
  PointNeighb = new int[N_RootVertices];
  memset(PointNeighb, 0, N_RootVertices*SizeOfInt);
  
  for (i=0;i<N_AllLocVert;i++)
    PointNeighb[VertexNumbers[i]]++;
  
  for (i=0;i<N_RootVertices;i++)
    if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];
  
   MaxRankPerV = MaxCpV;
   MaxCpV++; // ACCOUNTING FOR AN EXTRA COLUMN GIVING INFO OF NUMBER OF CELLS SHARING THE VERTEX
   
   delete [] PointNeighb;
  
  /** ********************************************************************************************/
  /** STEP 7 : CREATE AN ARRAY CONTAINING INFO REGARDING THE NEIGHBOURS OF VERTICES (CELL INDEX) */
  /** ********************************************************************************************/
   
   PointNeighb = new int[N_RootVertices*MaxCpV] ;  
   memset(PointNeighb, 0, (N_RootVertices*MaxCpV)*SizeOfInt);

   for(i=0;i<N_Cells;i++)
    {
     for(j=0;j<N_VertInCell;j++)
      {
       M = VertexNumbers[i*N_VertInCell + j] *MaxCpV ;
       PointNeighb[M]++;
       PointNeighb[M + PointNeighb[M]] = i;
      }  
    }
    
   delete [] NumberVertex;
   delete [] Vertices;
   


    /**Metis partition done!! code for all processors */
    //MPI_Bcast(Cell_Rank, N_Cells, MPI_INT, 0, comm);
    

/********** VARIABLES FOR MULTIGRID ****************/	

  int Nchildren,childn,parentglobalno;

 
  Nchildren = coll->GetCell(0)->GetParent()->GetN_Children();

  /** copy the parent MPI_ID and global cell_no to children */
  for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);	
      cell->SetSubDomainNo((cell->GetParent())->GetSubDomainNo());
      childn = (cell->GetParent())->GetChildNumber(cell);	
      parentglobalno = (cell->GetParent())->GetGlobalCellNo();	
      cell->SetGlobalCellNo(parentglobalno*Nchildren+childn);
      cell->SetLocalCellNo(i);
      Cell_Rank[i] = cell->GetSubDomainNo();
     }
  
  /** *********************************************/
  /** STEP 8.2 : SET ALL VERTICES TO -1 */
  /** *********************************************/ 
    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      if(cell->GetSubDomainNo()==rank) 
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

  /** *********************************************/
  /** STEP 9 : Fill the VertNeibRanks info        */
  /** first column contains how many ranks contain this vertex further columns contain the rank ID of the subdomains */
  /** *********************************************/ 
    VertNeibRanks = new int[N_RootVertices*MaxCpV];
    memset(VertNeibRanks, 0, N_RootVertices*MaxCpV*SizeOfInt);

    HaloCellIndex = new int[size];
    for(i=0;i<size;i++)
      HaloCellIndex[i] = -1;

    HaloCellGlobalNo = new int[MaxCpV];
    HaloCellLocVertNo = new int[MaxCpV];
   
   
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      ID = cell->GetSubDomainNo();

      /**run only through own cells */
      if(ID==rank) 
       {
        cell->SetAsOwnCell();
      //set SubDomainVert if any vert in this cell is so needed for moving meshes and setting cross vertex
        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);
	  
	  if(CurrVert->GetClipBoard() != -1)
           continue;

          CurrVert->SetClipBoard(5);
          N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;

          N_CellsInThisVert =  PointNeighb[N];
          N_OwnIncidentCells = 0;
          N_NeibIncidentCells = 0;

          // check! any subdomain cell containg this vert has to be added
          // Setting all halo cells
          UPDATE = TRUE;
          for(ii=1;ii<=N_CellsInThisVert;ii++)
           {
            VertexCellNo = PointNeighb[N + ii];
            Vertexcell = coll->GetCell(VertexCellNo);
            VertexcellID = Vertexcell->GetSubDomainNo();

            if(VertexcellID != rank)
             {
              //vertex belong to diff subdomain 
              if(UPDATE)
               {
                cell->SetAsDependentCell();
		
                CurrVert->SetAsSubDomainVert();
                HaloCellIndex[VertexcellID]= VertexCellNo;
                UPDATE = FALSE;
               }
              Vertexcell->SetAsHaloCell();

               /** changed on 4 Feb 2012 - by Sashi */
               if(HaloCellIndex[VertexcellID] == -1)
                HaloCellIndex[VertexcellID]= VertexCellNo;

             }
            } //  for(ii=1;ii< 

           if(CurrVert->IsSubDomainVert())
             for(ii=1;ii<=N_CellsInThisVert;ii++)
              {
               VertexCellNo = PointNeighb[N + ii];
               Vertexcell = coll->GetCell(VertexCellNo);
               VertexcellID = Vertexcell->GetSubDomainNo();

               // all cells associated with this vert are DepCells
               Vertexcell->SetAsDependentCell();

               if(VertexcellID==rank)
                continue;

               N_CellsIn_b = VertNeibRanks[N];
               UPDATE = TRUE;

               for(jj=1;jj<=N_CellsIn_b;jj++)
                if( VertNeibRanks[N + jj]==VertexcellID)
                 {
                  UPDATE = FALSE;
                  break;
                 }

               if(UPDATE)
                {
                 N2 = HaloCellIndex[VertexcellID];
                 HaloCellGlobalNo[VertNeibRanks[N]] = N2;

                 //find the local index of this vertex in the neib cell
                 jj=0;
                 while(CurrVert != (coll->GetCell(N2))->GetVertex(jj)) jj++;
                 HaloCellLocVertNo[VertNeibRanks[N]] = jj;

                 VertNeibRanks[N]++;
                 VertNeibRanks[N + VertNeibRanks[N]] = VertexcellID;

                 if(MaxRankPerV<VertNeibRanks[N])
                   MaxRankPerV=VertNeibRanks[N];
                } //   if(UPDATE)
              } //  for(ii=1;ii<=N_CellsInT

           N_SubDomIn_a = VertNeibRanks[N];
           Temp = VertNeibRanks+(N+1);
 
    
           /** set only one cell (lowest index cell) from each subdomain as the Halo cell for this vertex */
           CurrVert->SetSubDomainInfo(N_SubDomIn_a,  Temp, HaloCellGlobalNo, HaloCellLocVertNo);

            //reset
            for(jj=0;jj<size;jj++)
             HaloCellIndex[jj] = -1;

           } // for(j=0;j<N_V
       } // if(ID==ra
     } //  for(i=0;i<N_Cel 
    
    MaxRankPerV++;    // including own rank

    delete [] HaloCellIndex;
    delete [] HaloCellGlobalNo;
    delete [] HaloCellLocVertNo;

     
   
  /** *************************************************/
  /** STEP 10 : Get own SubDomain collection of cells */
  /** *************************************************/ 
     N_LocalCells = 0;
     N_OwnCells = 0;

     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       if(ID==rank || cell->IsHaloCell() )
        N_LocalCells++;

       if(ID==rank)
        N_OwnCells++;
      } //  for (i=0;i<N_Cells;i++)

    N_HalloCells = N_LocalCells - N_OwnCells;      

//     printf("Rank: %d N_own Cells %d\n ", rank,  N_OwnCells);   
      

     if(N_LocalCells)
      {
       // collect the own collection
       SubDomainCells = new TBaseCell*[N_LocalCells];
       GlobalCellIndex = new int[N_LocalCells];
      }
     else
      {
       SubDomainCells = NULL;
       GlobalCellIndex = NULL;
      }


      
  /** *************************************************************/
  /** STEP 11 : Fill subdomain own cells and their global numbers */
  /** *************************************************************/ 
  VertexDel = new TVertex*[N_AllLocVert];
  EdgeDel = new TEdge*[N_Cells*12]; // max 12 in Hexa
  CellDel = new TBaseCell*[N_Cells];
     N=0;
     M=0;
     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       // in the collection, first own cells then Halo cells
       if(ID==rank)
        {
         SubDomainCells[N] =  cell;
         GlobalCellIndex[N] =  i;
         N++;
        }
      else if(cell->IsHaloCell() )
        {
         SubDomainCells[N_OwnCells + M] =  cell;
         GlobalCellIndex[N_OwnCells + M] =  i;
         M++;
        }
      else //cell no longer needed for this rank, delete!!!!!
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
          EdgeDel[N_EdgeDel++] = cell->GetEdge(j);

        for(j=0;j<N_VertInCell;j++)
         VertexDel[N_VertexDel++] = cell->GetVertex(j);

         CellDel[N_CellDel++] = cell;
       }
      } //  for (i=0;i<N_Cells;i++)


  /** ***********************************************************************************/
  /** STEP 12 : Delete other processors cells from own collection, excluding halo cells */
  /** ***********************************************************************************/ 
    // 
/*    if(N_VertexDel)
     Sort(VertexDel, N_VertexDel); 
    if(N_EdgeDel)
     Sort(EdgeDel, N_EdgeDel);*/
//     if(N_CellDel)
//      Sort(CellDel, N_CellDel);

  /*  Last=NULL;
    for(i=0;i<N_VertexDel;i++)
     {
      if( VertexDel[i] != Last)
       {
        delete VertexDel[i];
        Last = VertexDel[i];
       }
     }
   delete [] VertexDel;


    LastEdge=NULL;
    for(i=0;i<N_EdgeDel;i++)
     {
      if( EdgeDel[i] != LastEdge)
       {
        delete EdgeDel[i];
        LastEdge = EdgeDel[i];
       }
     }
   delete [] EdgeDel;

    LastCell=NULL;
    for(i=0;i<N_CellDel;i++)
     {
      if( CellDel[i] != LastCell)
       {
        delete CellDel[i];
        LastCell = CellDel[i];
       }
     }
   delete [] CellDel;

    if((N+M)!=N_LocalCells)
     {
      printf("Error in mesh partition %d \n", rank);
      MPI_Abort(comm, 0);
     } 
     
   // ***********************************************************************************/
   //** STEP 13 : Edges will be flaged as SubDomainEdges and subdomain infos are filled    */
   //** ************************************************************************************/ 
  
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);
          if(edge->GetClipBoard()==-1)
           {
            edge->SetClipBoard(5);
            edge->InitSubDomainInfo(rank);
           } // if(edge->GetCli
         } // for(j=0;j<N_Edges;
       } // if(cell->IsDependentCell
     } // for(i=0;i<N_Cells;i++

              
  
 /** change all subdomain joint from JointEqN into TSubDomainJoint */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      /**run only through Dependent cells */
      if(cell->IsDependentCell()) 
       {
        GblCellNr = GlobalCellIndex[i];
        ID = cell->GetSubDomainNo();

	//if( GblCellNr == 111615)
          // printf("rank %d local cell no %d\n", rank, i);
	
        /** change the faces, also edges in 3D */
        for(j=0;j<N_JointsInCell;j++)
         {
          Joint = cell->GetJoint(j);

          if(Joint->GetType() == JointEqN)
           {
            neib_cell = Joint->GetNeighbour(cell);
            Neib_ID = neib_cell->GetSubDomainNo();


            if(ID!=Neib_ID)
             {
              GlobCellNo =  neib_cell->GetGlobalCellNo();

              // this joint belongs to two SubDomains
              cell->SetAsSubDomainInterfaceCell();

              // find neib cell local face number
              m=0;
              while(Joint!=neib_cell->GetJoint(m)) m++;

              delete Joint;
              NewJoint = new TSubDomainJoint(cell, neib_cell, Neib_ID, GlobCellNo, m);
              cell->SetJoint(j, NewJoint);
              neib_cell->SetJoint(m, NewJoint);
              NewJoint->SetMapType();

              // set all edges in this face as SubDomainEdges
              ShapeDesc = cell->GetShapeDesc();
              ShapeDesc->GetFaceEdge(TmpFV, TmpLen, MaxLen);
              N_FaceEdges = TmpLen[j];
              for(n=0;n<N_FaceEdges;n++)
               {
                edge = cell->GetEdge(TmpFV[j*MaxLen+n]);
                edge->SetAsNotCrossEdgeFor(rank, Neib_ID);
               } // for(n=0;n<N_Edge

              // remove the Neib_ID from VertNeibRanks
              // these face vertices cannot be cross vertices for processors ID and Neib_ID
              ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
              for(ii=0;ii<TmpLen[j];ii++)
              {
               jj = TmpFV[j*MaxLen + ii];
               M = VertexNumbers[GblCellNr*N_VertInCell + jj] *MaxCpV;
               N = VertNeibRanks[M];

               for(jj=1;jj<=N;jj++)
                if(VertNeibRanks[M + jj] == Neib_ID )
                 {
                  VertNeibRanks[M + jj] = VertNeibRanks[M + N];
                  VertNeibRanks[M]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
              } // for(ii=0;ii<TmpL
             }//if(ID!=Neib_ID
           } // if(Joint->GetT
         }// for (j=0;j< ;
       } // if(cell->IsDependentCell() 
     }// for(i=0;i<N_Cel

     //      ===========================================================================================   
    /** set the Hallo BD face as BD face: 07 April 2015 - by Sashi  */
    for(i=0;i<N_Cells;i++)
     {     
      cell = coll->GetCell(i);
      
       // cell in this collection
       if(cell->GetSubDomainNo()==rank || cell->IsHaloCell() )
        { cell->SetClipBoard(i); }
       else
        {cell->SetClipBoard(-1); }
     } // 
     
    for(i=0;i<N_HalloCells;i++)
     {     
      cell =  SubDomainCells[N_OwnCells + i]; 
      for(j=0;j<N_JointsInCell;j++)
       {
        Joint = cell->GetJoint(j);

        if(Joint->GetType() == JointEqN)
	//if(Joint->GetType() == SubDomainHaloJoint)  
         {
	   //exit(0);
          neib_cell = Joint->GetNeighbour(cell);
          //j is a Halo BD face
          if(neib_cell->GetClipBoard() == -1)
	   {  
             Joint-> ChangeType(SubDomainHaloJoint);//BoundaryFace);
	   } // if(neib_cell->GetCl
       } //   if(Joint->GetTy  
      } //for(j=0;j<N_JointsInCel
     } //   for(i=0;i<N_HalloCells;i+
//      ===========================================================================================
     
    /** find cross edges (if any) */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];
        N_Edges=cell->GetN_Edges();
        ShapeDesc= cell->GetShapeDesc();
        ShapeDesc->GetEdgeVertex(EdgeVertex);

        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);

          if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
           {
            edge->SetClipBoard(5);
            N_CrossNeibs = edge->GetN_CrossNeibs();

            if(N_CrossNeibs==0)
             continue;

            cell->SetAsCrossEdgeCell();

            a = EdgeVertex[2*j];
            b = EdgeVertex[2*j+ 1];

            Vert_a = cell->GetVertex(a);
            edge->SetCrossNeibInfo(Vert_a);
            edge->GetCrossEdgeNeibs(N_CrossNeibs, CrossNeibsRank);

            M1 = VertexNumbers[M*N_VertInCell + a]*MaxCpV;
            N1 = VertNeibRanks[M1];

            M2 = VertexNumbers[M*N_VertInCell + b]*MaxCpV;
            N2 = VertNeibRanks[M2];

            for(jj=0;jj<N_CrossNeibs;jj++)
             {
              Neib_ID = CrossNeibsRank[jj];

               for(kk=1;kk<=N1;kk++)
                if(VertNeibRanks[M1 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M1 + kk] = VertNeibRanks[M1 + N1];
                  VertNeibRanks[M1]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)

               for(kk=1;kk<=N2;kk++)
                if(VertNeibRanks[M2 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M2 + kk] = VertNeibRanks[M2 + N2];
                  VertNeibRanks[M2]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
             } // for(jj=0;jj<N_CrossNei
           } //   if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
          } // for(j=0;j<N_Edge
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel
      
    /** find cross vertex (if any) */
    // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    for(i=0;i<N_HalloCells;i++)
     {
      SubDomainCells[N_OwnCells + i]->SetClipBoard(-1); 
      
      // store GlobalCellIndex
      SubDomainCells[N_OwnCells + i]->SetClipBoard_Par(GlobalCellIndex[N_OwnCells + i]);
     }
     
    /** set incident cell list for all vertices */
    IncCells = new TBaseCell*[MaxCpV];
    IncHalloCells = new TBaseCell*[MaxCpV];
	
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      M = GlobalCellIndex[i];    
     
      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        {
         // set the vertexcells for all vertices in dep. cells
         CurrVert = cell->GetVertex(j);
         if(CurrVert->GetClipBoard() != -1)
          continue;

         CurrVert->SetClipBoard(5);

         N = VertexNumbers[M*N_VertInCell + j] *MaxCpV;
         N_CellsInThisVert =  PointNeighb[N];

         for(ii=1;ii<=N_CellsInThisVert;ii++)
          IncCells[ii-1] = coll->GetCell(PointNeighb[N + ii]);

         CurrVert->SetVertexCells(N_CellsInThisVert, IncCells);

         /** added, hallo cells vertCells neib info on 11 Mar 2015 - Sashi */  
         for(ii=0;ii<N_CellsInThisVert;ii++) 
          {
           /** fill only for Hallo cells, since own Dependentcells' vertices got filled above */    
	   if(IncCells[ii]->IsHaloCell())
	    {
	      /** this Hallo cell is already handled */
	      if(IncCells[ii]->GetClipBoard() != -1)
               continue;

              IncCells[ii]->SetClipBoard(5);
 
              for(jj=0;jj<N_VertInCell;jj++)
               {
                ll = 0; 
                NeibCurrVert = IncCells[ii]->GetVertex(jj);

                /** even in Hallo cells, set neib info only for unset vertices */
                if( NeibCurrVert->GetNNeibs()!= 0)
                 continue;

                MM = IncCells[ii]->GetClipBoard_Par();
                NN = VertexNumbers[MM*N_VertInCell + jj] *MaxCpV;
                N_CellsInHaloVert =  PointNeighb[NN];

                for(kk=1;kk<=N_CellsInHaloVert;kk++)
                 {
                  HalloCell = coll->GetCell(PointNeighb[NN + kk]);
 
                  /** include only owna and Halo in the neib list */
                  if(HalloCell->GetSubDomainNo() == rank || HalloCell->IsHaloCell())
                   {    
                    IncHalloCells[ll] = HalloCell;
                    ll++;
                   } // if(HalloCell->GetSubDomainNo(
                  } // for(kk=1;kk<=N_CellsInHaloVert;kk++)
         
                 NeibCurrVert->SetVertexCells(ll, IncHalloCells);    
     
                } //   for(jj=0;jj<N_VertInCell;jj++)   
              }  // if(IncCells[ii]->IsHalo
          } //  for(jj=0;jj<N_VertInCell;jj++)
        } //  for (j=0;j<N_VertInCell;
     } //  for(i=0;i<N_OwnCell
    delete [] IncCells;
    delete [] IncHalloCells;
    
//    printf("Rank %d : Domain_Crop \n\n",rank);   
//    MPI_Finalize();  
//    exit(0);
    
        
     // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];

        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);

          // continue if the vert is already handled or not a subdomain vert
          if( (CurrVert->GetClipBoard() != -1) || !(CurrVert->IsSubDomainVert()))
           continue;


          CurrVert->SetClipBoard(5);

          M1 = VertexNumbers[M*N_VertInCell + j] *MaxCpV ;
          N1 = VertNeibRanks[M1];

          if(N1!=0)
           cell->SetAsCrossVertexCell();

          for(ii=1;ii<=N1;ii++)
            CurrVert->AddCrossNeib(VertNeibRanks[M1+ii]);

//            TDatabase::ParamDB->Par_P5 =0;
          } // for(j=0;j<N_VertInCell;j
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    if(N_LocalCells)
     { Domain->ReplaceTreeInfo(N_LocalCells, SubDomainCells, GlobalCellIndex, N_OwnCells); }
    else
     { Domain->SetN_OwnCells(N_OwnCells); }

        // initialize iterators
    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   delete [] VertNeibRanks;

 coll = Domain->GetCollection(It_Finest, 0);
 for(i=0;i<N_LocalCells;i++)
   {  
     cell = coll->GetCell(i);
     cell->SetCellIndex(i);
   }
   
 
  delete [] VertexNumbers;
  //delete [] Vert_Rank;
  delete [] Cell_Rank;
  delete [] PointNeighb;
//    printf("Rank %d : Domain_Crop \n\n",rank);   
//    MPI_Finalize();  
//    exit(0);
   //Barrier is needed, before calling FECommunicator, since neib process must have all info
   MPI_Barrier(comm);
} // Partition_Mesh2D()

#endif //mpi
