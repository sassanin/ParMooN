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
// @(#)SquareStructure1D.C
//
// Class:       TSquareStructure1D
//
// Purpose:     build and store a structure for a square matrix in 1d
//
// Author:      Sashikumaar Ganesan
//
// History:     17.05.2007 start implementation
//
// =======================================================================

// #include <DefineParams.h>

#include <FEDatabase2D.h>
#include <SquareStructure1D.h>
#include <string.h>
#include <Database.h>
#include <stdlib.h>
/** dummy constructor, needed only derives classes */
TSquareStructure1D::TSquareStructure1D()
: TSquareStructure()
{
}


TSquareStructure1D::~TSquareStructure1D()
{

}


/** generate the matrix structure, both space are 1D */
TSquareStructure1D::TSquareStructure1D(TFESpace1D *Space)
: TSquareStructure()
{
  int i,j,k,l,n,N_, n1,n2, m;
  int *GlobalNumbers;
  int *BeginIndex;
  int *Numbers, *nieb_Numbers;
  int N_Dirichlet, end, N_Hanging, Offset;
  int N_Inner, NE, nieb_i, nieb_e, nieb_n1;
  int *AuxPtr, *AuxPtr_delete, *HangingAuxPtr;
  int *KColAux, *HangingKColAux;
  int index, oldindex, EdgeIntegrals=Space->IsDGSpace();

  TCollection *Coll;

  TBaseCell *cell, *neigh;
  FE1D CurrentElement, CurrentNeighbour;

  FESpace = Space;

  // all dof are treated as unknowns !!!
  // no boundary description is used so far!!!
  N_Inner=Space->GetN_Inner();
  ActiveBound = N_Inner;
  N_Rows = N_Inner;
  N_Columns = N_Rows;
  N_Hanging = 0;
  HangingN_Entries=0;
  // AuxPtr[i] will contain an upper bound for the number of ...
  // matrix entries in row i
  l=N_Rows+1;
  AuxPtr=new int[l];
  memset(AuxPtr, 0, l*sizeof(int));

  GlobalNumbers=Space->GetGlobalNumbers();
  BeginIndex=Space->GetBeginIndex();

  // loop over all elements
  N_=Space->GetN_Cells();
  Coll = FESpace->GetCollection();
  
  // associate each cell with her number in the collection
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = Space->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();
    for(j=0;j<n1;j++)
    {
      k=Numbers[j];
      AuxPtr[k]+=n1;

     //DG additional entries necessary if integrals over the edges(points in 1D)
     // appear in the discretization
     if(EdgeIntegrals)
     {
      l=cell->GetN_Edges();                   // # edges
      for(m=0;m<l;m++)                        // for all edges
      {
       neigh = (cell->GetJoint(m))->GetNeighbour(cell);

       if(neigh)
       {
        n = neigh->GetClipBoard();
        CurrentNeighbour = Space->GetFE1D(n, neigh);
        n2 = TFEDatabase2D::GetFE1D(CurrentNeighbour)->GetN_DOF();
        AuxPtr[k]+=n2;
       }//if(neigh)
      }
     }//if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)

    }                                             // for(j=0;j<n1;
  }                                               //  for(i=0;i<N

  N_Entries = 0;

//   for(i=0;i<=N_;i++)
//     cout << i << "   " << AuxPtr[i] << endl;
//   cout << endl;
// exit(0);
  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=N_Rows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  //   cout << l << "  " << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    //     cout << AuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "AuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << AuxPtr[i] << endl;
  cout << endl;
  */
// exit(0);

  if(TDatabase::ParamDB->SC_VERBOSE)
    cout << "Upper bound: " << AuxPtr[N_Rows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows];                               // upper bound for number of matrix entries
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  N_=Space->GetN_Cells();
  Coll = FESpace->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = Space->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)                             // test space
    {
      n=Numbers[j];
      for(k=0;k<n1;k++)                           // ansatz space
      {
        m=Numbers[k];
        // this  node is a real node (inner or Neumann)
        index=AuxPtr[n];
        l=KColAux[index];
        // check whether this column is already in this row
        while(l!=-1 && l!=m)
        {
          index++;
          l=KColAux[index];
        }
        if(l==-1)
        {
          // this is a new column for this row
          KColAux[index]=m;
          N_Entries++;
        }
      }                                           // endfor k
      
     // DG part
     if(EdgeIntegrals)
     {
      NE=cell->GetN_Edges();                   // # edges
      for(nieb_e=0;nieb_e<NE;nieb_e++)                        // for all edges
      {
       neigh = (cell->GetJoint(nieb_e))->GetNeighbour(cell);

       if(neigh)
       {
        nieb_i = neigh->GetClipBoard();
        CurrentNeighbour = Space->GetFE1D(nieb_i, neigh);
        nieb_Numbers=GlobalNumbers+BeginIndex[nieb_i];
        nieb_n1 = TFEDatabase2D::GetFE1D(CurrentNeighbour)->GetN_DOF();
 
        for(k=0;k<nieb_n1;k++)                           // ansatz space
         {
          m=nieb_Numbers[k];
          // this  node is a real node (inner or Neumann)
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
           index++;
           l=KColAux[index];
          }
         if(l==-1)
          {
           // this is a new column for this row
           KColAux[index]=m;
           N_Entries++;
          }   
         } // for(k=0;k<nieb_n1;k++)	
       }  //if(neigh)
       
      } //  for(nieb_e=0;nieb_e<NE;nieb_e++)
     }//if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)

      
    }                                             // endfor j
  }                                               // for(i=0;

  
  /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=ActiveBound;
  N_=N_Rows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index;j++) //  && KColAux[j]!=-1
  cout << setw(4) << KColAux[j];
  cout << endl;
  }
  exit(0);
  */

  //   cout << "Number of matrix entries: ";
  //   cout << N_Entries << endl;
  //   cout << endl;

  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=ActiveBound;
  KCol=new int[N_Entries];
  RowPtr=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      KCol[index]=KColAux[j];
      index++;
    }                                             // endfor j
    RowPtr[i]=oldindex;
    //     cout << setw(4) << i << " RowPtr[i]: " << RowPtr[i] << endl;
  }                                               // endfor i

  // cout << "index: " << index << endl;
  // cout << "RowPtr[N_]: " << RowPtr[N_] << endl;
  Offset=index-RowPtr[N_];
  for(i=0,j=ActiveBound;i<=N_Hanging;i++,j++)
  {
    // cout << "HangingRowPtr[i]: " << HangingRowPtr[i] << endl;
    RowPtr[j]+=Offset;
    // cout << setw(4) << j << " RowPtr[j]: " << RowPtr[j] << endl;
  }

  /*
    // print out the whole matrix structure
    cout << endl;
    N_=N_Rows;
    for(i=0;i<N_;i++)
    {
      cout << RowPtr[i] << "---" << RowPtr[i+1]-1 << endl;
      cout << "Rows: " << setw(4) << i << ": ";
      end=RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
        cout << setw(4) << KCol[j];
  cout << endl;
  }
//   */

  // free KColAux
  delete [] KColAux;
  if(TDatabase::ParamDB->SC_VERBOSE)
  {
    cout << "Information on the stored matrix structure" << endl;
    cout << "Number of rows: " << N_Rows << endl;
    cout << "Number of columns: " << N_Columns << endl;
    cout << "Number of matrix entries: " << N_Entries << endl;
  }
}
