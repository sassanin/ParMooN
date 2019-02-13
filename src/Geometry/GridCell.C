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
   

#include <Database.h>
#include <GridCell.h>
#include <Vertex.h>
#include <stdlib.h>
#include <MooNMD_Io.h>

#ifdef __2D__
  #include <IsoBoundEdge.h>
  #include <IsoInterfaceJoint.h>
#else
  #include <BoundFace.h>
  #include <FEDatabase3D.h>
  #include <TetraAffin.h>
#endif

#ifdef __COMPAQ__
  #include <Utilities.h>
#endif

// Constructors
TGridCell::TGridCell(TRefDesc *refdesc, int reflevel) : TBaseCell(refdesc)
{
  Vertices = new TVertex*[RefDesc->GetN_OrigVertices()];

  Children = NULL;
  Parent = NULL;

  RefLevel = reflevel;
}

// Destructor
TGridCell::~TGridCell()
{
  if (Children)
  {
    cerr << "Error: Deleting Cell, but there are children!!!" << endl;
    exit (-1);
  }

  delete [] Vertices;
}

// Methods
int TGridCell::SetVertex(int Vert_i, TVertex *Vert)
{
  Vertices[Vert_i] = Vert;
  return 0;
}

TVertex *TGridCell::GetVertex(int Vert_i)
{
  return Vertices[Vert_i];
}

int TGridCell::GetN_Children()
{
  return Children ? RefDesc->GetN_Children() : 0;
}

int TGridCell::GetN_Parents()
{
  if ( Parent == NULL)
    return 0;
  else
    return 1;
}

TBaseCell *TGridCell::GetChild(int C_i)
{
  return Children[C_i];
}

TBaseCell *TGridCell::GetParent()
{
  return Parent;
}

int TGridCell::SetParent(TBaseCell *parent)
{
  Parent = parent;
  return 0;
}

int TGridCell::GetChildNumber(TBaseCell *Me)
{
  int i, N_ = GetN_Children();
  
  for (i=0;i<N_;i++)
    if (Children[i] == Me)
      return i;

  cerr << "Error: could not find child" << endl;
  return -1;
}

// #define __GRIDCELL_WITH_NUMBERS__

int TGridCell::PS(std::ofstream &dat, double scale, double StartX,
                  double StartY)
{
  int i;
  double text_x1,text_x2,text_y1,text_y2;
#ifdef __GRIDCELL_WITH_NUMBERS__
  double x, y;
#endif
  dat<<"newpath"<<endl;
  dat << (30 + (Vertices[0]->GetX() - StartX) * scale + .5) << " " <<
         (30 + (Vertices[0]->GetY() - StartY) * scale + .5) <<
         " M" << endl;
  text_x1 = Vertices[0]->GetX();
  text_y1 = Vertices[0]->GetY();

#ifdef __GRIDCELL_WITH_NUMBERS__
  x = Vertices[0]->GetX();
  y = Vertices[0]->GetY();
#endif
//if(!IsHaloCell() && IsDependentCell())

  for (i=1;i<RefDesc->GetN_OrigVertices();i++)
  {
    dat << (30 + (Vertices[i]->GetX() - StartX) * scale + .5) << " " <<
           (30 + (Vertices[i]->GetY() - StartY) * scale + .5) <<
           " L" << endl;
    if(Vertices[i]->GetX() != text_x1)
	text_x2 = Vertices[i]->GetX();
    if(Vertices[i]->GetY() != text_y1)
	text_y2 = Vertices[i]->GetY();

#ifdef __GRIDCELL_WITH_NUMBERS__
    x += Vertices[i]->GetX();
    y += Vertices[i]->GetY(); 
#endif
  }

text_x1 = (text_x1 + text_x2)/2;
text_y1 = (text_y1 + text_y2)/2;

#ifdef __GRIDCELL_WITH_NUMBERS__
  x /= RefDesc->GetN_OrigVertices();
  y /= RefDesc->GetN_OrigVertices();
#endif

  dat << (30 + (Vertices[0]->GetX() - StartX) * scale + .5) << " " <<
         (30 + (Vertices[0]->GetY() - StartY) * scale + .5) <<
         " L" << endl;
dat<<"closepath"<<endl;
#ifdef __GRIDCELL_WITH_NUMBERS__
  dat << (30 + (x - StartX) * scale + .5) << " " <<
         (30 + (y - StartY) * scale + .5);
  dat << " moveto" << endl;
  dat << "(" << ClipBoard << ") show" << endl;
#endif
  
#ifdef _MPI
if(IsHaloCell())
{
dat<<"gsave"<<endl;
//dat<< r<<" setgray"<<endl;
dat<<"0 1 0 setrgbcolor"<<endl;
dat<<"fill"<<endl;
dat<<"grestore"<<endl;
}
else if(IsDependentCell())
{
dat<<"gsave"<<endl;
dat<<"1 0 0 setrgbcolor"<<endl;
dat<<"fill"<<endl;
dat<<"grestore"<<endl;
}
else
{
dat<<"stroke"<<endl;
}
#else
 dat<<"stroke"<<endl;
#endif

#ifdef _MPI
for (i=0;i<RefDesc->GetN_OrigVertices();i++)
  {
    if(Vertices[i]->IsCrossVert())
    {
     //printf(" I'M HERE\n");
    dat << "newpath" <<endl;
    dat << (30 + (Vertices[i]->GetX() - StartX) * scale + .5) << " " <<
           (30 + (Vertices[i]->GetY() - StartY) * scale + .5) <<" "<<"M"<< endl;
    dat << (30 + (Vertices[i]->GetX() - StartX) * scale + .5) << " " <<
           (30 + (Vertices[i]->GetY() - StartY) * scale + .5) <<" "<<"5 0 360 arc fill"<< endl;
    dat<<"closepath"<<endl;
    }
  }
#endif
dat << (30 + (text_x1 - StartX) * scale + .5) << " " <<
         (30 + (text_y1 - StartY) * scale + .5) <<
         " M" << endl;
#ifdef _MPI
//dat <<"("<<GlobalCellNo<<") show"<<endl;
//dat <<"("<<CellIndex<<") show"<<endl;
dat <<"("<<GetSubDomainNo()<<") show"<<endl;
#endif
  return 0;
}


int TGridCell::Draw(std::ofstream &dat, double scale, double StartX,
                  double StartY)
{
  int i,j,n;

  n=RefDesc->GetN_OrigVertices();
  for (i=0;i<n;i++)
  {
    j=Joints[i]->GetType();
    if(j == BoundaryEdge || j == InterfaceJoint || 
       j == IsoBoundEdge || j == IsoInterfaceJoint)
    {
    dat << (int) (300 + (Vertices[i]->GetX() - StartX) * scale + .5) << " " <<
           (int) (300 + (Vertices[i]->GetY() - StartY) * scale + .5) <<
           " M" << endl;
    j=(i+1) % n;
    dat << (int) (300 + (Vertices[j]->GetX() - StartX) * scale + .5) << " " <<
           (int) (300 + (Vertices[j]->GetY() - StartY) * scale + .5) <<
           " L" << endl;
    }
  }

  return 0;
}

int TGridCell::MD_raw(std::ofstream &dat)
{
  int i, N_;
  float OutValues[8];

  N_ = GetN_Vertices();
  for (i=0;i<N_;i++)
  {
    OutValues[2*i] = Vertices[i]->GetX();
    OutValues[2*i + 1] = Vertices[i]->GetY();
  }

#ifdef __COMPAQ__
  SwapFloatArray(OutValues, 2*N_); 
#endif
  dat.write((const char *) OutValues, 2*N_*SizeOfFloat);
  
  return 0;
}

#ifdef __2D__
int TGridCell::Gen1RegGrid()
{
  int i, N_, ChildNumber, LocJointNum;
  const int *TmpCE, *TmpnEoE;
  int MaxLen;
  TRefDesc *ParentRefDesc = Parent->GetRefDesc();
  TBaseCell *RefCell;

  ParentRefDesc->GetChildEdge(TmpCE, MaxLen);
  ParentRefDesc->GetNewEdgeOldEdge(TmpnEoE);

  N_ = GetN_Edges();

  for (i=0;i<N_;i++)
  {
    if (Joints[i]->GetType() != BoundaryEdge &&
        Joints[i]->GetType() != IsoBoundEdge)

      if (!Joints[i]->GetNeighbour(this))
      {
        ChildNumber = Parent->GetChildNumber(this);

        LocJointNum = TmpnEoE[TmpCE[ChildNumber * MaxLen + i]];

        if (!Parent->GetJoint(LocJointNum)->GetNeighbour(Parent))
          Parent->Ref1Reg(LocJointNum, RefCell);
      }

    if (Joints[i]->GetType() == InterfaceJoint ||
        Joints[i]->GetType() == IsoInterfaceJoint)
      if (!Joints[i]->GetNeighbour(this))
        Ref1Reg(i, RefCell);
  }

  return 0;
}

#else

int TGridCell::Gen1RegGrid()
{
  if(!Parent)
    return 0;

  int i, N_, CurrLocEdge;
  const int *TmpCE, *TmpnEoE;
  int MaxLen;
  TRefDesc *ParentRefDesc = Parent->GetRefDesc(), *GrandParentRefDesc;
  TBaseCell *RefCell;
  TBaseCell* Grandfather, *CurrCell;
  TJoint* LastJoint, *Joint;
  const int *TmpEF, *TmpEF2;
  int TmpEFMaxLen, TmpEF2MaxLen;

  if(RefLevel > 1 && (Grandfather = Parent->GetParent()))
  {
    GrandParentRefDesc = Grandfather->GetRefDesc();
    Grandfather->GetShapeDesc()->GetEdgeFace(TmpEF, TmpEFMaxLen);

    for(int i=0; i<Grandfather->GetN_Edges(); ++i)
    {
      for(int j=0; j<TmpEFMaxLen; ++j)
      {
        LastJoint = Grandfather->GetJoint(TmpEF[2*i+j]);

        if(!(CurrCell = LastJoint->GetNeighbour(Grandfather)))
          continue;

        CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(Grandfather, i);

        while(CurrCell != Grandfather && CurrCell)
        {
          if(CurrLocEdge == -1)
          {
            std::cout << "Error!\n"; exit(-1);
          }
          
          // If neighbour is not refined, refine it regular
          //if(CurrCell->GetRefDesc()->GetType() == NoRef)
          if(CurrCell->GetN_Children() == 0)
          {
            CurrCell->SetRegRefine();
            CurrCell->Refine(RefLevel-1);
            
            if(RefLevel>2)
            {
              for(int k=0; k<CurrCell->GetN_Children(); ++k)
                CurrCell->GetChild(k)->Gen1RegGrid();
            }
          }
          
          // Get next joint which contains this edge
          CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
          if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
            Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
          else
            Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

          // Get new element and the index of our edge in this element
          CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
          CurrCell = Joint->GetNeighbour(CurrCell);
          
          LastJoint = Joint;
        }
      } // j=0..N_JointsPerEdge
    }
  }

  return 0;
}
#endif


#ifdef __2D__
int TGridCell::Ref1Reg(int LocJointNum, TBaseCell *&RefCell)
{
  int ChildNumber;
  TBaseCell *LocCell;
  const int *TmpCE, *TmpnEoE;
  int MaxLen;
  TRefDesc *ParentRefDesc = Parent->GetRefDesc();

  ParentRefDesc->GetChildEdge(TmpCE, MaxLen);
  ParentRefDesc->GetNewEdgeOldEdge(TmpnEoE);

  ChildNumber = Parent->GetChildNumber(this);

  LocCell = Parent->GetJoint(TmpnEoE[TmpCE[ChildNumber * MaxLen +
              LocJointNum]])->GetNeighbour(Parent);

  if (LocCell)
  {
    LocCell->SetRegRefine();
    LocCell->Refine(RefLevel);
  }
  else
  {
    cerr << "Error: generation of 1-regular grid impossible" << endl;
    return -1;
  }

  RefCell = LocCell;
  return 0;
}

#else

int TGridCell::Ref1Reg(int LocJointNum, TBaseCell *&RefCell)
{
  int ChildNumber;
  TBaseCell *LocCell;
  const int *TmpCF, *TmpnFoF;
  int MaxLen;
  TRefDesc *ParentRefDesc = Parent->GetRefDesc();

  ParentRefDesc->GetChildFace(TmpCF, MaxLen);
  ParentRefDesc->GetNewFaceOldFace(TmpnFoF);

  ChildNumber = Parent->GetChildNumber(this);

  LocCell = Parent->GetJoint(TmpnFoF[TmpCF[ChildNumber * MaxLen +
              LocJointNum]])->GetNeighbour(Parent);

  if (LocCell)
  {
    LocCell->SetRegRefine();
    LocCell->Refine(RefLevel);
  }
  else
  {
    cerr << "Error: generation of 1-regular grid impossible" << endl;
    return -1;
  }

  RefCell = LocCell;
  return 0;
}
#endif


#ifdef __2D__
int TGridCell::MakeConfClosure()
{
  int i, clip, N_;
  TBaseCell *Neighb;

  N_ = GetN_Edges();
  // set triangle/quadrangle bit(7/8)
  clip = 1 << (N_ + 4); // mac64 warning bracket added
  // set edge bits(0-3) and number of edges (bits 4-6)
  for (i=0;i<N_;i++)
      // there is a neighbour cell
    if ((Neighb = Joints[i]->GetNeighbour(this)))
	// this neighbour cell is refined regularly or 
	// it has children 
	if (Neighb->ExistChildren() || Neighb->GetClipBoard() > 512)
         {
	  if (Neighb->ExistChildren() ) OutPut("1a");
	  if (Neighb->GetClipBoard()>512) OutPut("2b");
	  OutPut(endl);
        clip += (1 << i) + 16;
      }
  // set regular refinement bit(9)
  if ((clip > 160 && clip < 256) || clip > 304)
    clip |= 512;

  ClipBoard = clip;

  N_ = clip & 496;
  if (N_ == 160 || N_ == 304)
  {
    if (clip & 128)
      i = (clip & 7) ^ 7;
    else
      i = (clip & 15) ^ 15;

    i--;
    if (i == 3) i--;
    if (i == 7) i = 3;

    if (Joints[i]->InnerJoint())
    {
      if ((Neighb = Joints[i]->GetNeighbour(this)))
        Neighb->MakeConfClosure();
      else
      {
        Ref1Reg(i, Neighb);
        Neighb->Check1Reg();
        Joints[i]->GetNeighbour(this)->MakeConfClosure();
      }
    }
  }

  return 0;
}

#else // 3D
int TGridCell::MakeConfClosure()
{
  int i, j, clip, N_Edges, LocEdge, LocFace, LocFace2, CurrLocEdge, N_ToRefine=0, CurrClip;
  TBaseCell *Neigh, *CurrCell;
  const int *TmpEF, *TmpEF2;
  int TmpEFMaxLen, TmpEF2MaxLen;
  TJoint *Joint, *LastJoint;
  bool EdgeRefined;

  if(Children)
    {
      Derefine();
      std::cerr << "Element has already been refined\n";
      //exit(-1);
    }

  clip = 0; //GetClipBoard(); 

  GetShapeDesc()->GetEdgeFace(TmpEF, TmpEFMaxLen);
  N_Edges = GetN_Edges();

  for (i=0;i<N_Edges;i++)
    {
      EdgeRefined = false;
      LocEdge = i;
      for(j=0; j<TmpEFMaxLen; ++j)
        {
    LastJoint = GetJoint(TmpEF[2*LocEdge+j]);

    if(!(CurrCell = LastJoint->GetNeighbour(this)))
      continue;

    CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(this, LocEdge);

    while(CurrCell != this && CurrCell)
            {
        if(CurrLocEdge == -1)
                {
      std::cout << "Error!\n"; exit(-1);
                }

        // Check if Edge is refined in CurrCell
        CurrClip = CurrCell->GetClipBoard();
                
        if(CurrClip != -1 && ((CurrClip >> CurrLocEdge) & 1))
                {
      // Check if Cell is already marked for refinement
      if((clip >> LocEdge) & 1)
        std::cout << "Should not happen\n";

      clip |= 1 << LocEdge;
      EdgeRefined = true;
      break;
                }

        // Get next joint which contains this edge
        CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
        if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
    Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
        else
    Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

        // Get new element and the index of our edge in this element
        CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
        CurrCell = Joint->GetNeighbour(CurrCell);

        LastJoint = Joint;
            }

    if(EdgeRefined)
      break;
        } // j=0..N_JointsPerEdge

      //if(EdgeRefined)
      //    std::cout << "Edge " << i << " is refined by an other element\n";
    } // i=0--N_edges

  // Check if a RefDesc exists for this setting
  switch(clip)
    {
      // NoRef
    case 0:
      break;
      std::cout << "Unexpected behaviour!\n";
      exit(-1);
      // Bisections
    case 1:
    case 2:
    case 4:
    case 8:
    case 16:
    case 32:
      // Bis0X
    case 3:
    case 5:
    case 9:
    case 17:
    case 33:
      // Bis 1X
    case 6:
    case 10:
    case 18:
    case 34:
      // Bis2X
    case 12:
    case 20:
    case 36:
      // Bis 3X
    case 24:
    case 40:
      // Bis 4X
    case 48:
      // QuadX
    case 7:
    case 25:
    case 50:
    case 44:
      break;
      // Refinement Rule was found
      // Reg
    case 63:
    default:
      // No Suitable Refinement Rule was found

      N_Edges = GetN_Edges();
      for(i=0; i<N_Edges; ++i)
        {
    // Check if edge is marked for refinement
    if( clip >> i) continue;
//     if( (clip >> i) & 1 == 1) continue;
    
    LocEdge = i;
    for(j=0; j<TmpEFMaxLen; ++j)
            {
        LastJoint = GetJoint(TmpEF[2*LocEdge+j]);

        if(!(CurrCell = LastJoint->GetNeighbour(this)))
    continue;

        CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(this, LocEdge);

        while(CurrCell != this && CurrCell)
                {
      // Add this cell to Gk
      if(CurrCell->GetClipBoard() > 0 && ((CurrCell->GetClipBoard() >> CurrLocEdge) & 1) == 1)
        std::cout << "Edge has been marked for refinement already before\n";

      if(CurrCell->GetClipBoard() != 0)
        {
          N_ToRefine++;         
          CurrCell->SetClipBoard(0);
        }
      
      // Get next joint which contains this edge
      CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
      if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
        Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
      else
        Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

      // Get new element and the index of our edge in this element
      CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
      CurrCell = Joint->GetNeighbour(CurrCell);
      LastJoint = Joint;
                }
            } // j=0..N_JointsPerEdge
        }
      clip = 63;
    }

  if(clip == 0)
    std::cerr << "Fehler!\n";

  SetClipBoard(clip);

  return N_ToRefine;
}
#endif


#ifdef __2D__
int TGridCell::Check1Reg()
{
  int i, N_;
  TBaseCell *Neighb;

  N_ = GetN_Edges();
  for (i=0;i<N_;i++)
    if (Joints[i]->InnerJoint())
      if (!(Neighb = Joints[i]->GetNeighbour(this)))
      {
        Ref1Reg(i, Neighb);
        Neighb->Check1Reg();
      }

  return 0;
}
#else
int TGridCell::Check1Reg()
{
  int i, N_;
  TBaseCell *Neighb;

  N_ = GetN_Faces();
  for (i=0;i<N_;i++)
    if (Joints[i]->InnerJoint())
      if (!(Neighb = Joints[i]->GetNeighbour(this)))
      {
        Ref1Reg(i, Neighb);
        Neighb->Check1Reg();
      }

  return 0;
}
#endif


#ifdef __2D__
int TGridCell::SetRegRefine()
{
  if (ExistChildren()) return -1;

  switch (RefDesc->GetShapeDesc()->GetType())
  {
    case Triangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriReg];
           break;
    case Parallelogram:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + ParallReg];
           break;
    case Quadrangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + QuadReg];
           break;
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + RectReg];
           break;
     default:
          cerr << "Only 2D cells allowed" << endl;
     exit(-1);
     break;
  }

  return 0;
}

int TGridCell::Set1Refine(int i)
{
  if (ExistChildren()) return -1;

  switch (RefDesc->GetShapeDesc()->GetType())
  {
    case Triangle:
      
        switch (i)
          {
           case 0:  
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriBis0];
           break;
           
           case 1:  
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriBis1];
           break;
           
           case 2:  
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriBis2];
           break;                     

          }   
	   
	   
           break;
    case Parallelogram:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + ParallReg];
           break;
    case Quadrangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + QuadReg];
           break;
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + RectReg];
           break;
     default:
          cerr << "Only 2D cells allowed" << endl;
     exit(-1);	    
  }

  return 0;
}



#endif

#ifdef __3D__
int TGridCell::Set1Refine(int i)
{
 
    cerr << "Error: Not yet implemented !!!" << endl;
    exit (-1);
 
}

int TGridCell::SetRegRefine()
{
/*
  double x0, y0, z0;
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double x5, y5, z5;
  double x6, y6, z6;
  double x7, y7, z7;
  double x8, y8, z8;
  double x9, y9, z9;
  double l0, l1, l2;
  double min;
  int index;
*/

  if (ExistChildren()) return -1;

  switch (RefDesc->GetShapeDesc()->GetType())
  {
    case Triangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriReg];
           break;
    case Parallelogram:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + ParallReg];
           break;
    case Quadrangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + QuadReg];
           break;
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + RectReg];
           break;
    case Tetrahedron:
         /*
           Vertices[0]->GetCoords(x0, y0, z0);
           Vertices[1]->GetCoords(x1, y1, z1);
           Vertices[2]->GetCoords(x2, y2, z2);
           Vertices[3]->GetCoords(x3, y3, z3);
           x4 = 0.5*(x0+x1);
           y4 = 0.5*(y0+y1);
           z4 = 0.5*(z0+z1);
           x5 = 0.5*(x1+x2);
           y5 = 0.5*(y1+y2);
           z5 = 0.5*(z1+z2);
           x6 = 0.5*(x0+x2);
           y6 = 0.5*(y0+y2);
           z6 = 0.5*(z0+z2);
           x7 = 0.5*(x0+x3);
           y7 = 0.5*(y0+y3);
           z7 = 0.5*(z0+z3);
           x8 = 0.5*(x1+x3);
           y8 = 0.5*(y1+y3);
           z8 = 0.5*(z1+z3);
           x9 = 0.5*(x2+x3);
           y9 = 0.5*(y2+y3);
           z9 = 0.5*(z2+z3);

           l0 = (x6-x8)*(x6-x8) + (y6-y8)*(y6-y8) +(z6-z8)*(z6-z8);
           l1 = (x5-x7)*(x5-x7) + (y5-y7)*(y5-y7) +(z5-z7)*(z5-z7);
           l2 = (x4-x9)*(x4-x9) + (y4-y9)*(y4-y9) +(z4-z9)*(z4-z9);

           min = l0;
           index = 0;
           if(l1 < min)
           {
             min = l1;
             index = 1;
           }
           if(l2 < min)
           {
             min = l2;
             index = 2;
           }

           RefDesc = TDatabase::RefDescDB[N_SHAPES + TetraReg0 + index];
          */
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TetraReg];
           break;
    case Hexahedron:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + HexaReg];
           break;
    case Brick:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + BrickReg];
           break;
    default:
     break;
	   
  }

  return 0;
}
#endif

int TGridCell::SetNoRefinement()
{
  switch (GetType())
  {
    case Triangle:
           RefDesc = TDatabase::RefDescDB[Triangle];
           break;
    case Quadrangle:
    case Parallelogram:
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[Quadrangle];
           break;
     default:
          cerr << "Only 2D cells allowed" << endl;
     exit(-1);
     break;    
  }

  return 0;
}

int TGridCell::IsToRefine()
{
  return RefDesc->IsToRefine() ? Children == NULL ? TRUE : FALSE : FALSE;
}

#ifdef __2D__
// return coordinates of mid point P_j on edge J_i
int TGridCell::LineMidXY(int J_i, int P_j, double &X, double &Y)
{
  const int *TmpVerts, *TmpEV;
  const double *TmpPos;
  int MaxLen;
  double T_0, T;

  TDatabase::RefDescDB[N_SHAPES + RefDesc->GetEdgeRef(J_i)]->
                  GetInnerVerts(TmpVerts, TmpPos, MaxLen);

  RefDesc->GetShapeDesc()->GetEdgeVertex(TmpEV);

  --P_j;

  if (Joints[J_i]->GetType() == BoundaryEdge ||
      Joints[J_i]->GetType() == InterfaceJoint ||
      Joints[J_i]->GetType() == IsoBoundEdge ||
      Joints[J_i]->GetType() == IsoInterfaceJoint)
  {
    T_0 = TmpPos[P_j * MaxLen + 1];

    switch(Joints[J_i]->GetType())
    {
      case IsoBoundEdge:
      T = T_0*((TIsoBoundEdge *) Joints[J_i])->GetEndParameter()
          +(1.-T_0)*((TIsoBoundEdge *) Joints[J_i])->GetStartParameter();

      ((TIsoBoundEdge *) Joints[J_i])->GetXYofT(T, X, Y);
      break;

      case IsoInterfaceJoint:
      T = T_0*((TIsoInterfaceJoint *) Joints[J_i])->GetEndParameter()
        +(1.-T_0)*((TIsoInterfaceJoint *) Joints[J_i])->GetStartParameter();

      ((TIsoInterfaceJoint *) Joints[J_i])->GetXYofT(T, X, Y);
      break;

      case InterfaceJoint:
      T = T_0*((TInterfaceJoint *) Joints[J_i])->GetEndParameter()
          +(1.-T_0)*((TInterfaceJoint *) Joints[J_i])->GetStartParameter();

      ((TInterfaceJoint *) Joints[J_i])->GetXYofT(T, X, Y);
      break;

      case BoundaryEdge:
      T = T_0 * ((TBoundEdge *) Joints[J_i])->GetEndParameter() 
          + (1. - T_0) * ((TBoundEdge *) Joints[J_i])->GetStartParameter();

      ((TBoundEdge *) Joints[J_i])->GetXYofT(T, X, Y);
    break;
     default:
          cerr << "Only edges" << endl;
     exit(-1);
      break;  
    }
  }
  else
  {
    X = Vertices[TmpEV[J_i * 2]]->GetX();
    Y = Vertices[TmpEV[J_i * 2]]->GetY();

    X += TmpPos[P_j * MaxLen + 1] * (Vertices[TmpEV[J_i * 2 + 1]]->GetX() - X);
    Y += TmpPos[P_j * MaxLen + 1] * (Vertices[TmpEV[J_i * 2 + 1]]->GetY() - Y);
  }

  return 0;
}

// return parameters on boundary of subedge SJ_j on edge J_i
int TGridCell::LineMidT(int J_i, int SJ_j, double &T_0, double &T_1)
{
  const int *TmpVerts;
  const double *TmpPos;
  int MaxLen;

  TDatabase::RefDescDB[N_SHAPES + RefDesc->GetEdgeRef(J_i)]->
                  GetInnerVerts(TmpVerts, TmpPos, MaxLen);

  if (SJ_j)
    T_0 = TmpPos[(SJ_j - 1) * MaxLen + 1];
  else
    T_0 = 0.;

  if (SJ_j + 1 != TDatabase::RefDescDB[N_SHAPES + RefDesc->
                    GetEdgeRef(J_i)]->GetN_Edges())
    T_1 = TmpPos[SJ_j * MaxLen + 1];
  else
    T_1 = 1.;

  return 0;
}
#endif

bool TGridCell::PointInCell(double X, double Y)
{
  int i, j, N_ = RefDesc->GetN_OrigEdges();
  bool test = TRUE;
  double NX, NY, DX, DY;
//   double len;


  for (i=0;i<N_;i++)
  {
    j = (i + 1) % N_;

    DX = Vertices[i]->GetX();
    DY = Vertices[i]->GetY();

    NX = DY - Vertices[j]->GetY();
    NY = Vertices[j]->GetX() - DX;

//     len = sqrt(NX*NX+NY*NY);
//     NX /= len;
//     NY /= len;

    DX = X - DX;
    DY = Y - DY;

    // check if given point is identical to vertex
    if(DX*DX + DY*DY < 1e-8)
    {
      return true;
    }

//     len = sqrt(DX*DX+DY*DY);
//     DX /= len;
//     DY /= len;

    test = (bool) (test && ((NX * DX + NY * DY) > -1e-8)); 
  }

  return test;
}

#ifdef __3D__
bool TGridCell::PointInCell(double X, double Y, double Z)
{
  double xmin = 1e+8,  ymin = 1e+8, zmin = 1e+8;
  double xmax = -1e+8,  ymax = -1e+8, zmax = -1e+8;
  double x, y, z;
  int i;
  bool ret = FALSE;
  TTetraAffin *rt;
  double xi, eta, zeta;

  rt = (TTetraAffin *)TFEDatabase3D::GetRefTrans3D(TetraAffin);
  
  // nur fuer kantenparallele Hexaeder !!!
  if (GetType()==Hexahedron || GetType()==Brick) 
  {
    for (i=0; i<8; i++)
    {
      Vertices[i]->GetCoords(x,y,z);
      if (x > xmax) xmax=x;
      if (y > ymax) ymax=y;
      if (z > zmax) zmax=z;
      if (x < xmin) xmin=x;
      if (y < ymin) ymin=y;
      if (z < zmin) zmin=z;
    } // endfor i
    if (xmin <= X && X <= xmax && 
        ymin <= Y && Y <= ymax &&
        zmin <= Z && Z <= zmax)
    { 
      ret = TRUE;
    }
  }
  else
  {
    // Error("PointInCell works only for hexahedrons !" << endl);
    rt->SetCell(this); 
    rt->GetRefFromOrig(X, Y, Z, xi, eta, zeta);
    if(-1e-4 < xi && xi < 1.0001 &&
       -1e-4 < eta && eta < 1.0001 &&
       -1e-4 < zeta && zeta < 1.0001 &&
       xi + eta + zeta < 1.0001)
    {
      ret = TRUE;
    }
    // cout << "ref: " << xi << "  " << eta << "  " << zeta << endl;
    rt->GetRefFromOrig(xi, eta, zeta, X, Y, Z);
    // cout << "orig: " << X << "  " << Y << "  " << Z << endl;
  }
  return ret;
}
#endif // __2D__

int TGridCell::GetGeoLevel()
{
  int i = 0;
  TBaseCell *parent = this;

  while ((parent = parent->GetParent()))
    i++;

  return i;
}

int TGridCell::GetSubGridID()
{
  if (Parent)
    return Parent->GetSubGridID();
  else
    return -1;
}
/** compute number of edges at the boundary */
/** call: ((TGridCell *)cell)->GetN_BoundaryEdges(); */
int TGridCell::GetN_BoundaryEdges()
{
  int bdry_edges = 0, i, N_;
    
  N_ = GetN_Edges();

  for (i=0;i<N_;i++)
  {
    if (Joints[i]->GetType() == BoundaryEdge ||
        Joints[i]->GetType() == IsoBoundEdge)
    {
	bdry_edges++;
    }
  }
   
  return bdry_edges;    
}
#ifdef __3D__
/** compute number of faces at the boundary */
int TGridCell::GetN_BoundaryFaces()
{
  int bdry_faces = 0, i, N_;
    
  N_ = GetN_Faces();

  for (i=0;i<N_;i++)
  {
    if (Joints[i]->GetType() ==  BoundaryFace  ||
        Joints[i]->GetType() == IsoBoundFace)
    {
	bdry_faces++;
    }
  }
   
  return bdry_faces;    
}
#endif 

/** compute number of vertices at the boundary */
/** BEFORE THIS ROUTINE CAN BE CALLED */
/**    collection->MarkBoundaryVertices(); */
/** has to be called */
/** call: ((TGridCell *)cell)->GetN_BoundaryVertices(); */
int TGridCell::GetN_BoundaryVertices()
{
  int bdry_vertices = 0, i, N_;
  TVertex *vertex;
    
  N_ = GetN_Vertices();

  // loop over the edges
  for (i=0;i<N_;i++)
  {
      vertex = GetVertex(i);
      if (vertex->GetClipBoard())
	  bdry_vertices++;
  }
   
  return bdry_vertices;    
}
