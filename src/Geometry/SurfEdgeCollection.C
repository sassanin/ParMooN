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
   
#include <SurfEdgeCollection.h>
#include <MooNMD_Io.h>

#include <stdlib.h>
#include <string.h>

// CTOR
TSurfEdgeCollection::TSurfEdgeCollection (TCollection *Coll, int N_SurfaceJoints,
					  int *CellNumbers, int *JointNumbers)
{
  Init (Coll, N_SurfaceJoints, CellNumbers, JointNumbers);
  
  OutPut("N_Vertices: " << mN_Vertices << endl);
  OutPut("N_Faces: " << N_SurfaceJoints << endl);
  OutPut("N_Edges: " << mN_Edges << endl);
}
					  
// DTOR
TSurfEdgeCollection::~TSurfEdgeCollection ()
{  
  for (int i=0;i<2*mN_Vertices;++i)
  {
    if ( mEdgeHash[i] ) delete [] mEdgeHash[i];
  }    
  
  delete [] mEdgeHash;
}

/// Methods
void TSurfEdgeCollection::Init(TCollection *Coll, int N_SurfaceJoints,
			      int *CellNumbers, int *JointNumbers)
{
  int CellNr, JointNr, counter;
  const int *TmpFV, *TmpLen, *indices;
  int MaxLen, len, a, b, hash;
  TBaseCell *Cell;
  TVertex *Vertex;
  
//   // reset clipboard
//   for (int i=0;i<N_SurfaceJoints;++i)
//   {
//     CellNr = CellNumbers[i];
//     JointNr = JointNumbers[i];
//     
//     Cell = Coll->GetCell(CellNr);
//     Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);    
//     
//     len = TmpLen[JointNr];
//     indices = TmpFV+JointNr*MaxLen;
//     
//     for (int j=0;j<len;++j)
//     {
//       Vertex = Cell->GetVertex(indices[j]);
//       Vertex->SetClipBoard(-1);
//     }
//   } // end for i
  
  // count vertices on surface
  counter = 0;
  for (int i=0;i<N_SurfaceJoints;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    
    Cell = Coll->GetCell(CellNr);
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);    
    
    len = TmpLen[JointNr];
    indices = TmpFV+JointNr*MaxLen;
    
    for (int j=0;j<len;++j)
    {
      Vertex = Cell->GetVertex(indices[j]);
      if ( Vertex->GetClipBoard() == -1 )
      {
	Vertex->SetClipBoard(counter);
	++counter;
      }
    }
  } // end for i
  
  mN_Vertices = counter;
  mEdgeHash = new int* [2*mN_Vertices];
  memset(mEdgeHash, 0, 2*mN_Vertices*sizeof(mEdgeHash[0]));
  
  // find edges
  counter = 0;
  for (int i=0;i<N_SurfaceJoints;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    
    Cell = Coll->GetCell(CellNr);
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);    
    
    len = TmpLen[JointNr];
    indices = TmpFV+JointNr*MaxLen;
    
    for (int j=0;j<len;++j)
    {
      a = Cell->GetVertex(indices[j])->GetClipBoard();
      b = Cell->GetVertex(indices[(j+1)%len])->GetClipBoard();
      
      IncreaseBucket(a, b, counter);
    }
  } // end for i
  
  mN_Edges = counter;
}

void TSurfEdgeCollection::IncreaseBucket(int a, int b, int &counter)
{
  int hash, *Bucket, len;
    
  hash = Flip(a,b);
  
  if ( mEdgeHash[hash] == NULL )
  {
    mEdgeHash[hash] = new int [4];
    Bucket = mEdgeHash[hash];
    Bucket[0] = 1;
    Bucket[1] = a;
    Bucket[2] = b;
    Bucket[3] = counter;
    
    ++counter;
  }
  else 
  {
    Bucket = mEdgeHash[hash];
    len = Bucket[0];
    
    if ( FindEdge(Bucket, a, b) == -1 )
    {
      mEdgeHash[hash] = new int[3*len+4];
      memcpy(mEdgeHash[hash], Bucket, (3*len+1)*sizeof(Bucket[0]));
      
      delete [] Bucket;
      
      Bucket = mEdgeHash[hash];
      Bucket[0] = len+1;
      Bucket[3*len+1] = a;
      Bucket[3*len+2] = b;
      Bucket[3*len+3] = counter;
      
      ++counter;
    }
  }
}

int TSurfEdgeCollection::FindEdge(int *Bucket, int a, int b)
{
  int len;
  
  if ( Bucket == NULL) return -1;
  
  len = Bucket[0];
  
  Flip(a, b);
  
  for (int i=0;i<len;++i)
  {
    if ( a == Bucket[3*i  +1] && b == Bucket[3*i+1+1] )
      return Bucket[3*i+2+1];
  }
  
  return -1;
}
