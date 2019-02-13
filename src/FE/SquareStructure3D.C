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
// @(#)SquareStructure3D.C        1.6 09/17/99
// 
// Class:       TSquareStructure3D
//
// Purpose:     build and store a structure for a square matrix in 3d
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================

#include <DefineParams.h>
#include <FEDatabase3D.h>
#include <SquareStructure3D.h>
#include <HNDesc.h>
#include <HangingNode.h>
#include <string.h>
#include <Database.h>

/** dummy constructor, needed only for derived classes */
TSquareStructure3D::TSquareStructure3D()
  : TSquareStructure()
{
}

/** generate the matrix structure, both spaces are 3D */
TSquareStructure3D::TSquareStructure3D(TFESpace3D *Space)
: TSquareStructure()
{
  int i,j,k,l,n,N_, n1, m; 
  int *Numbers;

  int N_Dirichlet;
  int N_Inner;
  int N_Hanging;
  int *BoundNodeBounds, N_NonDiri, N_BoundaryNodeTypes;
  int HangingBound;

  int *AuxPtr, *HangingAuxPtr;
  int *KColAux, *HangingKColAux;
  int index, oldindex;

  int Offset, *DOF, end, begin;
  THangingNode **HangingNodes;
  THangingNode *hn;

  TCollection *Coll;

  TBaseCell *cell;
  FE3D CurrentElement;

  FESpace = Space;

  ActiveBound=Space->GetN_ActiveDegrees();
  HangingBound=Space->GetHangingBound();

  N_BoundaryNodeTypes=Space->GetN_DiffBoundaryNodeTypes();
  BoundNodeBounds=Space->GetN_BoundaryNodes();
  N_NonDiri = 0;
  for(i=0;i<N_BoundaryNodeTypes;i++)
    N_NonDiri += BoundNodeBounds[i];

  N_Dirichlet=Space->GetN_Dirichlet();
  N_Inner=Space->GetN_Inner();
  N_Hanging=Space->GetN_Hanging();

  N_Rows = ActiveBound+N_Hanging+N_Dirichlet;
  N_Columns = N_Rows;
  // assembles matrices without shorter rows for non-active dof
  if (TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
  {
    N_NonDiri += N_Dirichlet;
    N_Dirichlet = 0;
    ActiveBound = N_Inner + N_NonDiri;
  }
  ColOrder = 0;
  // AuxPtr[i] will contain an upper bound for the number of 
  // matrix entries in row i
  l=N_Rows+1;
  AuxPtr=new int[l];
  memset(AuxPtr, 0, l*sizeof(int));

  l=N_Hanging+1;
  HangingAuxPtr=new int[l];
  memset(HangingAuxPtr, 0, l*sizeof(int));

  Offset=ActiveBound;

  // loop over all elements 
  N_=Space->GetN_Cells();
  Coll = FESpace->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=Space->GetGlobalDOF(i);

    CurrentElement = Space->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      k=Numbers[j];
      if(k<ActiveBound) 
      {
        AuxPtr[k]+=n1;
      }
      else
      {
        if(k<HangingBound) HangingAuxPtr[k-Offset]+=n1;
      } // endif
    } // endfor j
  } // endfor i

  // add rows for Dirichlet nodes in  space
  N_Entries=N_Dirichlet;
  Offset=ActiveBound+N_Hanging;
  for(i=0,j=Offset;i<N_Dirichlet;i++,j++)
  {
    AuxPtr[j]=1;
  }

// #ifdef __3D__
  // add couplings for hanging nodes of  space
  HangingNodes=Space->GetHangingNodes();
  Offset=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    // there is the additional entry in diagonal
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes() + 1;
    AuxPtr[j]=n;
    N_Entries+=n;
    // cout << "AuxPtr[j]: " << AuxPtr[j] << endl;
  }

  // additional space for storing the columns caused by the 
  // hanging nodes of  space => some new columns
  Offset=ActiveBound;
  HangingNodes=Space->GetHangingNodes();
  // cout << "N_Hanging: " << N_Hanging << endl;
  for(i=0;i<N_Hanging;i++)
  {
    hn=HangingNodes[i];
    m=HangingAuxPtr[i];
    // cout << "HangingAuxPtr[i]: " << m << endl;
    DOF=hn->GetDOF();
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes();
    for(j=0;j<n;j++)
    {
      k=DOF[j];
      if(k<ActiveBound)
        AuxPtr[k] += m;
    }
  }
// #endif

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=N_Rows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  // cout << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    // cout << AuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "AuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << AuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound: " << AuxPtr[N_Rows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[N_Rows]; // upper bound for number of matrix entries 
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

// #ifdef __3D__
  // sum up the array HangingAuxPtr, HangingAuxPtr[i] will now contain 
  // the index for Hanging KColAux array, the column numbers for row i 
  // are in the intervall // [ HangingAuxPtr[i], HangingAuxPtr[i+1] )
  N_=N_Hanging;
  l=HangingAuxPtr[0];
  HangingAuxPtr[0]=0;
  // cout << HangingAuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=HangingAuxPtr[i+1];
    HangingAuxPtr[i+1]=HangingAuxPtr[i]+l;
    l=k;
    // cout << HangingAuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "HangingAuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << HangingAuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound for hanging nodes: ";
  // cout << HangingAuxPtr[N_Hanging] << endl;

  // get memory for HangingKColAux array, initialize it with -1
  l=HangingAuxPtr[N_Hanging]; //upper bound for number of matrix entries 
  // cout << "l= " << l << endl;
  HangingKColAux=new int[l];
  memset(HangingKColAux, -1, sizeof(int)*l);
// #endif
  HangingN_Entries=0;

  N_=Space->GetN_Cells();
  Coll = FESpace->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=Space->GetGlobalDOF(i);

    CurrentElement = Space->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      for(k=0;k<n1;k++)
      {
        m=Numbers[k];
        n=Numbers[j];
        if(n<ActiveBound)
        {
          // this  node is a real node (inner or Neumann)
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            N_Entries++;
          }
        }
        else
        {
// #ifdef __3D__
          if(n<HangingBound)
          {
            // this node is a hanging node in  space
            index=HangingAuxPtr[n-ActiveBound];
            l=HangingKColAux[index];
            // check whether this column is already in this row
            while(l!=-1 && l!=m)
            {
              index++; l=HangingKColAux[index];
            }
            if(l==-1)
            {
              // this is a new column for this row
              HangingKColAux[index]=m;
              HangingN_Entries++;
            }
          } // endif
// #endif
        } // endif
      } // endfor k
    } // endfor j
  } // endfor i

// #ifdef __3D__
  // check hanging node data
/*
  cout << endl;
  cout << "check hanging node data" << endl;
  N_=N_Hanging;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<index && HangingKColAux[j]!=-1;j++) 
      cout << setw(4) << HangingKColAux[j];
    cout << endl;
  }
*/
  
  // cout << "Number of matrix entries (hanging nodes): ";
  // cout << HangingN_Entries << endl;
  // cout << endl;

  // compress HangingKColAux array to HangingKCol by deleting all -1's
  // build the HangingRowPtr array
  N_=N_Hanging;
  HangingKCol=new int[HangingN_Entries];
  HangingRowPtr=HangingAuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<m && HangingKColAux[j]!=-1;j++)
    {
      HangingKCol[index]=HangingKColAux[j];
      index++;
    } // endfor j
    HangingRowPtr[i]=oldindex;
    // cout << HangingRowPtr[i] << endl;
  } // endfor i
  HangingRowPtr[N_]=index;

  // free HangingKColAux
  delete HangingKColAux;

  // add the additional columns from hanging nodes to other nodes
  Offset=ActiveBound;
  HangingNodes=Space->GetHangingNodes();
  for(i=0;i<N_Hanging;i++)
  {
    // cout << "hanging node: " << i << endl;
    hn=HangingNodes[i];
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();

    for(j=0;j<n;j++) // loop over all nodes in coupling
    {
      k=DOF[j];
      // cout << "k= " << k << endl;
      if(k<ActiveBound)
      {
        // node is either inner or Neumann node

        end=HangingAuxPtr[i+1];
        for(oldindex=HangingAuxPtr[i];oldindex<end;oldindex++)
        {
          m=HangingKCol[oldindex];
          // cout << "m= " << m << endl;
          index=AuxPtr[k];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            N_Entries++;
          }
        } // endfor
      } // endif
    } // endfor j
  } // endfor i
// #endif

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
    for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++) 
      cout << setw(4) << KColAux[j];
    cout << endl;
  }
  */
  
  // cout << "Number of matrix entries: ";
  // cout << N_Entries << endl;
  // cout << endl;

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
    } // endfor j
    RowPtr[i]=oldindex;
    // cout << setw(4) << i << " RowPtr[i]: " << RowPtr[i] << endl;
  } // endfor i

  // cout << "index: " << index << endl;
  // cout << "RowPtr[N_]: " << RowPtr[N_] << endl;
  Offset=index-RowPtr[N_];
  for(i=0,j=ActiveBound;i<=N_Hanging;i++,j++)
  {
    // cout << "HangingRowPtr[i]: " << HangingRowPtr[i] << endl;
    RowPtr[j]+=Offset;
    // cout << setw(4) << j << " RowPtr[j]: " << RowPtr[j] << endl;
  }

  j=ActiveBound+N_Hanging;
  Offset=RowPtr[j];
  for(i=0;i<N_Dirichlet;i++,j++)
  {
    RowPtr[j+1]=RowPtr[j]+1;
    // cout << setw(4) << j+1 << " RowPtr[j+1]: " << RowPtr[j+1] << endl;
  }

// #ifdef __3D__
  // add information for hanging and Dirichlet nodes into matrix
  HangingNodes=Space->GetHangingNodes();
  Offset=ActiveBound;
  m=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();
    index=AuxPtr[j];
    KCol[index]=m;
    index++;
    m++;
    for(k=0;k<n;k++)
    {
      // cout << "index: " << index << " DOF[k]" << DOF[k] << endl;
      KCol[index]=DOF[k];
      index++;
    }
  }
// #endif

  // add Dirichlet rows
  j=HangingBound;
  index=RowPtr[ActiveBound+N_Hanging];
  for(i=0;i<N_Dirichlet;i++)
  {
    // cout << "index: " << index << endl;
    KCol[index]=j;
    j++;
    index++;
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
*/

  // free KColAux
  delete KColAux;


#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif
  {
   if(TDatabase::ParamDB->SC_VERBOSE>1)
   OutPut("Information on the stored matrix structure" << endl);
   OutPut("Number of rows: " << N_Rows << endl);
   OutPut("Number of columns: " << N_Columns << endl);
   OutPut("Number of matrix entries: " << N_Entries << endl);
  }

} 

/** generate the matrix structure, all arrays are already defined */
TSquareStructure3D::TSquareStructure3D(int n, int N_entries, int *col_ptr,
int *row_ptr)
: TSquareStructure()
{
  N_Rows = n;
  N_Columns = n;
  ActiveBound = n;
  N_Entries = N_entries;
  HangingN_Entries = 0;
  KCol = col_ptr;
  HangingKCol = NULL;
  RowPtr = row_ptr;
  HangingRowPtr = NULL;
  FESpace = NULL;
  ColOrder = 0;
}

TSquareStructure3D::TSquareStructure3D(int n)
 : TSquareStructure(n), FESpace(NULL)
{
}

TSquareStructure3D::~TSquareStructure3D()
{
}
