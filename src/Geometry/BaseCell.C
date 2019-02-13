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
   

#ifdef _MPI
#  include "mpi.h"
#endif

#include <BaseCell.h>
#include <Joint.h>
#include <stdlib.h>

#include <BoundComp3D.h>
#include <BoundFace.h>
#include <Edge.h>


// Constructor
TBaseCell::TBaseCell(TRefDesc *refdesc)
{
 int i, N_;

  RefDesc = refdesc;

  N_ = RefDesc->GetShapeDesc()->GetN_Joints();  
  Joints = new TJoint*[N_];

  for (i=0;i<N_;i++)
   Joints[i] = NULL;

  ClipBoard = 0;
  Reference_ID = 0;
  CellIndex = -1;
  region = 0;
  LayerCell = 0;  
  
  normalOrientation = NULL;

#ifdef __3D__   
  N_ = RefDesc->GetN_OrigEdges();
  Edges = new TEdge*[N_];
  for (i=0;i<N_;i++)
   Edges[i] = NULL;
#endif

#ifdef  _MPI
  ClipBoard_Par  = -1;
  SubDomainNumber = 0;
  GlobalCellNo = -1;
  SubDomainLocalCellNo = -1;

  OwnCell=FALSE;
  HaloCell=FALSE;
  SubDomainInterfaceCell=FALSE;
  DependentCell=FALSE;
  N_NeibProcesses = 0;
  NeibProcessesIds = NULL;
#endif
}

// Destructor
TBaseCell::~TBaseCell()
{
  int i, N_;
  TJoint *CurrJoint;

#ifdef __2D__
  N_ = RefDesc->GetN_OrigEdges();
#else
  N_ = RefDesc->GetN_OrigFaces();
#endif
  for (i=0;i<N_;i++)
  {
    CurrJoint = Joints[i];
    switch (CurrJoint->GetType())
    {
      case JointEqN:
      case MortarBaseJoint:
      case InterfaceJoint:
      case IsoInterfaceJoint:
#ifdef __3D__
      case InterfaceJoint3D:
      case IsoInterfaceJoint3D:
#endif
        if (CurrJoint->GetNeighbour(this))
	{ CurrJoint->Delete(this);}
        else
	{ delete CurrJoint;}
        break;

      default:
        delete CurrJoint;
        break;
    }
  }

  delete Joints;
}

double TBaseCell::Get_hK(int cell_measure)
{
  switch (cell_measure)
  {
    case 0:                                     // diameter
      return this->GetDiameter();
      //case 1: // with reference map
      //OutPut("cell measure " << endl);
      //return this->GetLengthWithReferenceMap();
    case 2:                                     // shortest edge
      return this->GetShortestEdge();
      break;
    case 1:                                     // with reference map
    case 3:                                     // measure
      return sqrt(this->GetMeasure());
      break;
    case 4:                                     // mesh size in convection direction, this is just a dummy
      return this->GetDiameter();
      break;
    case 5:                                     // take value from an array
      // this is in general not the diameter but a pw constant value
      // which is needed for some reasons
      return this->GetDiameter();
      break;
    default:                                    // diameter
      OutPut("CELL_MEASURE " << cell_measure << " not available!!!" << endl);
      return this->GetDiameter();
      break;
  }
}

// added 25.04.2010 for fixing refinement problem
#ifdef __3D__
void TBaseCell::CorrectBoundaryVertices(TVertex **NewVertices, TJoint **NewJoints)
{
  int N_NewFaces = RefDesc->GetN_Faces();
  const int *FaceVertex;
  int MaxN_VpF;
  
  TBoundComp3D *BoundComp;
  TVertex *BoundVertices[3], *Vertex;
  double x, y, z, t, s;

  for (int i=0;i<N_NewFaces;++i)
  {
    if ( NewJoints[i]->GetType() == BoundaryFace )
    {
      BoundComp = ((TBoundFace*) NewJoints[i])->GetBoundComp();
      
      if ( BoundComp->GetType() == Plane ) continue;
      
      RefDesc->GetFaceVertex(FaceVertex, MaxN_VpF);
      
      for (int j=0;j<MaxN_VpF;++j)
      {
	Vertex = NewVertices[FaceVertex[MaxN_VpF*i+j]];
	
	Vertex->GetCoords(x ,y, z);
	
	BoundComp->GetTSofXYZ(x, y, z, t, s);
	BoundComp->GetXYZofTS(t, s, x, y, z);
	
	Vertex->SetCoords(x, y, z);
      }     
    }
    else continue;
  }
  
}
#endif

// on each joint, decide whether the 
// global normal is outgoing (+1) or ingoing (-1)
void TBaseCell::SetNormalOrientation()
{
  TBaseCell *neighbCell;
  int nEdges = RefDesc->GetN_OrigEdges();
  #ifdef __3D__
  nEdges = RefDesc->GetN_OrigFaces();
  #endif
  if(normalOrientation != NULL) // nothing more to do
    return;
  normalOrientation = new int[nEdges];
  for (int i=0; i<nEdges;i++)
  {
    normalOrientation[i] = 1;
    TJoint *joint = Joints[i];
    if(joint->InnerJoint())
    {
      neighbCell = joint->GetNeighbour(this);
      if(neighbCell->GetCellIndex() < GetCellIndex()&& 
         Reference_ID == neighbCell->GetReference_ID()) 
        normalOrientation[i] = -1;
    }
  }
}



// Methods
#ifdef  _MPI
void TBaseCell::SetNeibProcessesIds(int *Neiblist)
{
 int i;

 if(N_NeibProcesses==0)
  {
   printf(" Set the N_NeibProcesses for this cell first !!!! \n");
   MPI_Finalize();
   exit(0);
  }

 if(NeibProcessesIds) delete [] NeibProcessesIds;

 NeibProcessesIds = new int[N_NeibProcesses];

 for(i=0;i<N_NeibProcesses;i++)
  NeibProcessesIds[i] = Neiblist[i];

}
#endif
