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
// %W% %G%
// 
// Class:       TFESpace3D
// Purpose:     class for all 3D finite element spaces
//
// Author:      Gunar Matthies (22.11.97)
//
// History:     start of implementation 22.11.97 (Gunar Matthies)
//
// =======================================================================

#include <DefineParams.h>

#include <Constants.h>
#include <FESpace3D.h>
#include <Joint.h>
#include <BoundFace.h>
#include <FEDatabase3D.h>

#include <Database.h>
#include <IsoInterfaceJoint3D.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>

#include <HangingNode.h>
#include <HNDesc.h>

#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>

#include <Edge.h>

/** Constructor */
TFESpace3D::TFESpace3D(TCollection *coll, char *name, char *description) :
     TFESpace(coll, name, description)
{
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  AllElements = NULL;
  ElementForShape = NULL;
}

// =====================================================================
// rules: if a neighbour element is not in the collection, its clipboard
//        must be set to -1
// =====================================================================

/** constructor for building a space with elements of order k */
TFESpace3D::TFESpace3D(TCollection *coll, char *name, char *description, 
                       BoundCondFunct3D *BoundaryCondition, int ord) :
     TFESpace(coll, name, description)
{
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  AllElements = NULL;
 
  ElementForShape = new FE3D[N_SHAPES];

  // build ElementForShape array
  switch(ord)
  {
     case 0: ElementForShape[Tetrahedron] = C_P0_3D_T_A;
             ElementForShape[Hexahedron] = C_Q0_3D_H_M;
             ElementForShape[Brick] = C_Q0_3D_H_A;
             break;

    case 1:  ElementForShape[Tetrahedron] = C_P1_3D_T_A;
             ElementForShape[Brick] = C_Q1_3D_H_A;
             ElementForShape[Hexahedron] = C_Q1_3D_H_M;
             break;
    case 12:
    case 2:  ElementForShape[Tetrahedron] = C_P2_3D_T_A;
             ElementForShape[Brick] = C_Q2_3D_H_A;
             ElementForShape[Hexahedron] = C_Q2_3D_H_M;
             break;
    case 13:
    case 3:  ElementForShape[Tetrahedron] = C_P3_3D_T_A;
             ElementForShape[Brick] = C_Q3_3D_H_A;
             ElementForShape[Hexahedron] = C_Q3_3D_H_M;
             break;
    case 14:
    case 4:  ElementForShape[Tetrahedron] = C_P3_3D_T_A;
             Error("No elements of order 4 for Tetrahedron!" << endl);
             ElementForShape[Brick] = C_Q4_3D_H_A;
             ElementForShape[Hexahedron] = C_Q4_3D_H_M;
             break;
    case -1: ElementForShape[Hexahedron] = N_Q1_3D_H_M;
             ElementForShape[Brick] = N_Q1_3D_H_A;
             ElementForShape[Tetrahedron] = N_P1_3D_T_A;
             break;

    case -2: ElementForShape[Hexahedron] = N_Q2_3D_H_M;
             ElementForShape[Brick] = N_Q2_3D_H_A;
             ElementForShape[Tetrahedron] = N_P2_3D_T_A;
             break;

    case -3: ElementForShape[Hexahedron] = N_Q3_3D_H_M;
             ElementForShape[Brick] = N_Q3_3D_H_A;
             ElementForShape[Tetrahedron] = N_P3_3D_T_A;
             Error("No nonconforming elements of order 3 for Tetrahedron!" << endl);
             break;

    case -4: ElementForShape[Hexahedron] = N_Q4_3D_H_M;
             ElementForShape[Brick] = N_Q4_3D_H_A;
             ElementForShape[Tetrahedron] = N_P4_3D_T_A;
             Error("No nonconforming elements of order 4 for Tetrahedron!" << endl);
             break;
    case -11: 
             ElementForShape[Hexahedron] = D_Q1_3D_H_M;
             ElementForShape[Brick] = D_Q1_3D_H_A;
             ElementForShape[Tetrahedron] = D_P1_3D_T_A;
             break;
    case -12:
             ElementForShape[Hexahedron] = D_Q2_3D_H_M;
             ElementForShape[Brick] = D_Q2_3D_H_A;
             ElementForShape[Tetrahedron] = D_P2_3D_T_A;
             break;
    case -13:
             //ElementForShape[Hexahedron] = D_Q2_3D_H_M;
             //ElementForShape[Brick] = D_Q2_3D_H_A;
             ElementForShape[Tetrahedron] = D_P3_3D_T_A;
             break;
    case -110:
             ElementForShape[Hexahedron] = D_P1_3D_H_M;
             ElementForShape[Brick] = D_P1_3D_H_A;
             ElementForShape[Tetrahedron] = D_P1_3D_T_A;
             break;
    case -120:
             ElementForShape[Hexahedron] = D_P2_3D_H_M;
             ElementForShape[Brick] = D_P2_3D_H_A;
             ElementForShape[Tetrahedron] = D_P2_3D_T_A;
             break;

    //========LOCALPROJECTION=============
    // Q1+bubble*P0
    case 100:
             ElementForShape[Hexahedron] = C_UL1_3D_H_M;
             ElementForShape[Brick] = C_UL1_3D_H_A;
             // ElementForShape[Tetrahedron] = C_UL1_3D_T_A;
             break;
    // Q2+bubble*P1
    case 201:
             ElementForShape[Hexahedron] = C_UL2_3D_H_M;
             ElementForShape[Brick] = C_UL2_3D_H_A;
             // ElementForShape[Tetrahedron] = C_UL2_3D_T_A;
             break;
    // Q3+bubble*P2
    case 302:
             ElementForShape[Hexahedron] = C_UL3_3D_H_M;
             ElementForShape[Brick] = C_UL3_3D_H_A;
             // ElementForShape[Tetrahedron] = C_UL3_3D_T_A;
             break;
    //====================================
    //Raviar-Thomas-0
    case 1000:
      //ElementForShape[Hexahedron] = N_RT0_3D_H_M;
      ElementForShape[Brick] = N_RT0_3D_H_A;
      ElementForShape[Tetrahedron] = N_RT0_3D_T_A;
      break;
   //Raviar-Thomas-1
    case 1001:
      //ElementForShape[Hexahedron] = N_RT1_3D_H_M;
      ElementForShape[Brick] = N_RT1_3D_H_A;
      ElementForShape[Tetrahedron] = N_RT1_3D_T_A;
      break;
    //Raviar-Thomas-2
    case 1002:
      ElementForShape[Brick] = N_RT2_3D_H_A;
      ElementForShape[Tetrahedron] = N_RT2_3D_T_A;
      break;
    //Raviar-Thomas-3
    case 1003:
      ElementForShape[Tetrahedron] = N_RT3_3D_T_A;
      break;


    //Brezzi-Douglas-Duran-Fortin-1
    case 1011:
        ElementForShape[Brick] = N_BDDF1_3D_H_A;
        ElementForShape[Tetrahedron] = N_BDDF1_3D_T_A;
        break;

    //Brezzi-Douglas-Duran-Fortin-2
    case 1012:
        ElementForShape[Brick] = N_BDDF2_3D_H_A;
        ElementForShape[Tetrahedron] = N_BDDF2_3D_T_A;
        break;
    //Brezzi-Douglas-Duran-Fortin-3
    case 1013:
        ElementForShape[Brick] = N_BDDF3_3D_H_A;
        ElementForShape[Tetrahedron] = N_BDDF3_3D_T_A;
        break;

    default: Error("unknown order" << endl);
             exit(-1);
             break;
  } // endswitch

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace(BoundaryCondition);

}

/** constructor for building a space with the given elements */
TFESpace3D::TFESpace3D(TCollection *coll, char *name, char *description,
               BoundCondFunct3D *BoundaryCondition,
               FE3D *fes) : TFESpace(coll, name, description)
{
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  ElementForShape = NULL;

  AllElements = fes;

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace(BoundaryCondition);
}

/** constructor for building a space with elements of order k */
TFESpace3D::TFESpace3D(TCollection *coll, char *name, char *description, 
                       BoundCondFunct3D *BoundaryCondition, SpaceType type,
                       int ord) :
     TFESpace(coll, name, description)
{
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  AllElements = NULL;
 
  ElementForShape = new FE3D[N_SHAPES];

  // build ElementForShape array
  switch(type)
  {
    // find velo space for discontinuous pressure
    case DiscP_USpace:
      switch(ord)
      {
        case 1:
          Error("This makes no sense." << endl);
          ElementForShape[Tetrahedron] = C_P1_3D_T_A;
          ElementForShape[Hexahedron] = C_Q1_3D_H_M;
          ElementForShape[Brick] = C_Q1_3D_H_A;
        break;

        case 2:
          ElementForShape[Tetrahedron] = C_B2_3D_T_A;
          ElementForShape[Hexahedron] = C_Q2_3D_H_M;
          ElementForShape[Brick] = C_Q2_3D_H_A;
        break;

        case 3:
          // ElementForShape[Tetrahedron] = C_B3_3D_T_A;
          ElementForShape[Hexahedron] = C_Q3_3D_H_M;
          ElementForShape[Brick] = C_Q3_3D_H_A;
        break;

        case 4:
          // ElementForShape[Tetrahedron] = C_B4_3D_T_A;
          ElementForShape[Hexahedron] = C_Q4_3D_H_M;
          ElementForShape[Brick] = C_Q4_3D_H_A;
        break;

        default:
          Error("Space is not available" << endl);
          exit(-1);
      }
    break;
    
    // find pressure space for discontinuous pressure
    case DiscP_PSpace:
      switch(ord)
      {
        case 0:
          ElementForShape[Tetrahedron] = C_P0_3D_T_A;
          ElementForShape[Hexahedron] = C_Q0_3D_H_M;
          ElementForShape[Brick] = C_Q0_3D_H_A;
        break;

        case 1:
	    //Error("This will not work on Tetrahedrons!!" << endl);
          ElementForShape[Tetrahedron] = D_P1_3D_T_A;
          ElementForShape[Hexahedron] = D_P1_3D_H_M;
          ElementForShape[Brick] = D_P1_3D_H_A;
        break;

        case 2:
          Error("This will not work on Tetrahedrons!!" << endl);
          // ElementForShape[Tetrahedron] = D_P2_3D_T_A;
          ElementForShape[Hexahedron] = D_P2_3D_H_M;
          ElementForShape[Brick] = D_P2_3D_H_A;
        break;

        case 3:
          Error("This will not work on Tetrahedrons!!" << endl);
          // ElementForShape[Tetrahedron] = D_P3_3D_T_A;
          ElementForShape[Hexahedron] = D_P3_3D_H_M;
          ElementForShape[Brick] = D_P3_3D_H_A;
        break;

        default:
          Error("Space is not available" << endl);
          exit(-1);
      }
    break;
    
    // find pressure and velo space for continuous pressure
    case ContP_USpace:
    case ContP_PSpace:
      switch(ord)
      {
        case 1:
          ElementForShape[Tetrahedron] = C_P1_3D_T_A;
          ElementForShape[Hexahedron] = C_Q1_3D_H_M;
          ElementForShape[Brick] = C_Q1_3D_H_A;
        break;

        case 2:
          ElementForShape[Tetrahedron] = C_P2_3D_T_A;
          ElementForShape[Hexahedron] = C_Q2_3D_H_M;
          ElementForShape[Brick] = C_Q2_3D_H_A;
        break;

        case 3:
          ElementForShape[Tetrahedron] = C_P3_3D_T_A;
          ElementForShape[Hexahedron] = C_Q3_3D_H_M;
          ElementForShape[Brick] = C_Q3_3D_H_A;
        break;

        case 4:
          ElementForShape[Tetrahedron] = C_P3_3D_T_A;
          Error("No elements of order 4 for Tetrahedron!" << endl);
          ElementForShape[Hexahedron] = C_Q4_3D_H_M;
          ElementForShape[Brick] = C_Q4_3D_H_A;
        break;

        case -2:
          ElementForShape[Tetrahedron] = C_P3_3D_T_A;
          Error("No elements of order 4 for Tetrahedron!" << endl);
          ElementForShape[Hexahedron] = B_IB2_3D_H_M;
          ElementForShape[Brick] = B_IB2_3D_H_A;
        break;
        default:
          ErrMsg("Space is not available.");
          exit(-1);
      } // endswitch
    break;
    
    // find velo space for nonconforming fe
    case Non_USpace:
      switch(ord)
      {
        case 1:
          ElementForShape[Tetrahedron] = N_P1_3D_T_A;
          ElementForShape[Hexahedron] = N_Q1_3D_H_M;
          ElementForShape[Brick] = N_Q1_3D_H_A;
        break;

        case 2:
          ElementForShape[Hexahedron] = N_Q2_3D_H_M;
          ElementForShape[Brick] = N_Q2_3D_H_A;
          ElementForShape[Tetrahedron] = N_P2_3D_T_A;
        break;

        case 3:
          ElementForShape[Hexahedron] = N_Q3_3D_H_M;
          ElementForShape[Brick] = N_Q3_3D_H_A;
          ElementForShape[Tetrahedron] = N_P3_3D_T_A;
          Error("No nonconforming elements of order 3 for Tetrahedron!" << endl);
        break;

        case 4:
          ElementForShape[Hexahedron] = N_Q4_3D_H_M;
          ElementForShape[Brick] = N_Q4_3D_H_A;
          ElementForShape[Tetrahedron] = N_P4_3D_T_A;
          Error("No nonconforming elements of order 4 for Tetrahedron!" << endl);
        break;

        default:
          ErrMsg("Only first order nonconforming spaces available");
          exit(-1);
      }
    break;

    default:
      Error("Wrong space type" << endl);
      exit(-1);
  }

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace(BoundaryCondition);
}

/** return the FE Id for element i, corresponding to cell */
FE3D TFESpace3D::GetFE3D(int i, TBaseCell *cell)
{
  FE3D ret;

  if(AllElements)
    ret=AllElements[i];
  else
    ret=ElementForShape[cell->GetType()];

  return ret;
}

void TFESpace3D::FindUsedElements()
{
  TBaseCell *cell;
  int i, j, N_;
  int Used[N_FEs3D];

  memset(Used,0,N_FEs3D*SizeOfInt);

  N_ = N_Cells;
  for(i=0;i<N_;i++)
  {
    cell = Collection->GetCell(i);
    Used[GetFE3D(i, cell)] = 1;
  }

  for(i=0;i<N_FEs3D;i++)
    if(Used[i])
      N_UsedElements++;

  UsedElements = new FE3D[N_UsedElements];
  j=0;
  for(i=0;i<N_FEs3D;i++)
    if(Used[i])
    {
      UsedElements[j]=(FE3D)i;
      j++;
    }

/*
  OutPut(endl << "N_UsedElements: " << N_UsedElements << endl);
  for(i=0;i<N_UsedElements;i++)
    OutPut("UsedElement[" << i << "]: " << UsedElements[i] << endl);
*/
}

void TFESpace3D::ConstructSpace(BoundCondFunct3D *BoundaryCondition)
{
  int i, j, k, l, m, m2, n, comp, N_Faces, NFaces;
  int *v, N_FaceEdges;
  TBaseCell *cell, *neigh, *child1, *child2, *child3, *child4;
  TJoint *joint;
  TBoundComp3D *BoundComp;
  TBoundFace *BoundFace;
  double t0,t1;
  BoundCond Cond0, Cond1;

  TFE3DMapper *mapper;
  TFE3DMapper1Reg *mapper1reg;

  const int *TmpoFnF, *TmpLen1, *TmpFC, *TmpLen2, *TmpoFnlF;
  const int *TmpCTI;
  int MaxLen1, MaxLen2;
  TRefDesc *refdesc;

  int SumLocDOF;
  int count, *BoundaryUpperBound, DirichletUpperBound;
  int DirichletCounter, *BoundCounter, Counter;
  int SlaveMark, *BoundMark, InnerMark, DirichletMark;
  int DirichletOffset, InnerOffset, SlaveOffset;
  int *BoundOffset;
  int N_Slave, EMaxLen;

  FE3D FEType0, FEType1, FEType2;
  TFE3D *FE0, *FE1, *FE2;
  TFEDesc3D *FEDesc0_Obj, *FEDesc1_Obj, *FEDesc2_Obj;
  FEDesc3D FEDesc0, FEDesc1, FEDesc2;

  int I_K0, I_K1, I_K2, I_K3, I_K4, Bd_ID;
  int *J_K0, *J_K1, *J_K2, *J_K3, *J_K4;
  int *Indices0, *Indices1, *Indices2, *Indices3, *Indices4;
  int c1, c2, c3, c4, f1, f2, f3, f4;
  int chnum1, chnum2, chnum3, chnum4;
  int Twist1, Twist2, Twist3, Twist4;

  double X, Y, Z, T, S;
  double xp[4], yp[4], zp[4], tp[4], sp[4];
  double LinComb[4];
  const int *TmpFV, *TmpLen, *TmpFE, *ETmpLen;
  int MaxLen, N_Points;
  int N_Edges, N_EdgeNeibs, N_EdgeDOF, *EdgeDof, *NeibEdgeDof;
  int N_VertDof, N_VertInCell, N_VertNeibs;
  int owndof, neibdof, maptype, w, w0, w1, v0, v1, e;
  TEdge *edge;
  TBaseCell **EdgeNeibs, **VertNeibs;
  TVertex *Vert;
  const int *EdgeVertex, *NeibEdgeVertex;

  TInterfaceJoint3D *InterFace;

  THangingNode *hn;
  TVector<THangingNode *> *VHN = new TVector<THangingNode *>(20,20);
  TVector<int> *HNNumbers = new TVector<int>(20,20);

#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  N_EdgeDOF = -1;
  N_VertDof = -1;
#endif
 
  N_ActiveDegrees=0;
  N_SlaveDegrees=0;

  N_DiffBoundNodeTypes=N_BOUNDCOND-1; // not only Neumann nodes possible
  BoundaryUpperBound=new int[N_DiffBoundNodeTypes];
  BoundaryNodeTypes=new BoundCond[N_DiffBoundNodeTypes];
  BoundCounter=new int[N_DiffBoundNodeTypes];
  BoundMark=new int[N_DiffBoundNodeTypes];
  N_BoundaryNodes=new int[N_DiffBoundNodeTypes];
  BoundOffset=new int[N_DiffBoundNodeTypes];
  BoundaryNodesBound=new int[N_DiffBoundNodeTypes];

  BoundaryNodeTypes[0] = NEUMANN;
  BoundaryNodeTypes[1] = ROBIN;
  BoundaryNodeTypes[2] = SLIP;
  BoundaryNodeTypes[3] = FREESURF;
  BoundaryNodeTypes[4] = SLIP_FRICTION_PENETRATION_RESISTANCE;
  BoundaryNodeTypes[5] = INTERFACE;
  BoundaryNodeTypes[6] = SUBDOMAIN_INTERFACE;
  BoundaryNodeTypes[7] = SUBDOMAIN_HALOBOUND;
  BoundaryNodeTypes[8] = DIRICHLET_WEAK;
  N_BoundaryNodes[0] = 0;
  N_BoundaryNodes[1] = 0;
  N_BoundaryNodes[2] = 0;
  N_BoundaryNodes[3] = 0;
  N_BoundaryNodes[4] = 0;
  N_BoundaryNodes[5] = 0;
  N_BoundaryNodes[6] = 0;
  N_BoundaryNodes[7] = 0;
  N_BoundaryNodes[8] = 0;

//   OutPut("Number of Cells: " << N_Cells << endl);
     
  // reset clipboards to -1
  for(i=0;i<N_Cells;i++)
   {
    cell=Collection->GetCell(i);
    k=cell->GetN_Joints();
    for(j=0;j<k;j++)
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);
      if(neigh) neigh->SetClipBoard(-1);
    }
    cell->SetClipBoard(-1);
    
    /** set clipboard to all bound/interface edges, since only an edge from a cell be on the bound (in particaul  Dirchlet) */
    /** Finding BD (in particaul Dirichlet) Dof only may miss the edge BD DOF */
    /** unfortunately if this processor is the master of the missed BD  (in particaul Dirichlet) dof, it will give wrong result */ 
    /** only vert in a cell may be on the BD, so BD vertices also included 30.06.2012, sashi */
#ifdef _MPI   

    cell->SetClipBoard_Par(-1);
    
    FEType0 = GetFE3D(i, cell);
    FE0 = TFEDatabase3D::GetFE3D(FEType0);
    FEDesc0_Obj = FE0->GetFEDesc3D();    

    if(!(FEDesc0_Obj->IsEdgeVertData_Filled()) )
     {
       printf("Rank %d, FESpace3D Error! Edge and vertes data are not set in FEdesc3D for this FE \n", rank);
       MPI_Abort(MPI_COMM_WORLD,  0);  
     }

     N_Edges=cell->GetN_Edges();      
     for(j=0;j<N_Edges;j++)
      (cell->GetEdge(j))->SetClipBoard(-1);
     
      N_VertInCell = cell->GetN_Vertices();
      for(j=0;j<N_VertInCell;j++)
       (cell->GetVertex(j))->SetClipBoard(-1);    
      
      /** set dept. vert neib cells to -1, especiall to use own_FEspace and FEspace with (own + Hallo) cells */
      if(cell->IsDependentCell())
       {   
        for(j=0;j<N_VertInCell;j++)
         {
          Vert=cell->GetVertex(j);
          Vert->GetNeibs(N_VertNeibs, VertNeibs);

          for(k=0;k<N_VertNeibs;k++)
           VertNeibs[k]->SetClipBoard_Par(-1);  
          
         } //for(j=0;j<N_ 
        } // f(cell->IsDependentCell())
#endif   
   } // endfor i

  for(i=0;i<N_Cells;i++)
    Collection->GetCell(i)->SetClipBoard(i);
   
  /** if Hallo cells in this coll, the it will set to "i", otherwise -1 */
#ifdef _MPI  
  for(i=0;i<N_Cells;i++)
    Collection->GetCell(i)->SetClipBoard_Par(i);
#endif  
  
  // set number i into clipboard, count the number of local degrees of
  // freedom
  count=0;
  DirichletUpperBound=0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
    BoundaryUpperBound[i]=0;

  BeginIndex=new int[N_Cells+1];
  BeginIndex[0]=0;
 
  
  //==========================================================
  // ROutine to set: Boundary upper bound node numbers
  // 
  //==========================================================
  for(i=0;i<N_Cells;i++)
  {
    cell=Collection->GetCell(i);
    N_Faces = cell->GetN_Joints();
    N_Edges=cell->GetN_Edges();
    N_VertInCell = cell->GetN_Vertices();
    
    FEType0 = GetFE3D(i, cell);
    FE0 = TFEDatabase3D::GetFE3D(FEType0);
    count += FE0->GetSize();

    BeginIndex[i+1]=count;

    // find number of boundary nodes in this element
    FEDesc0_Obj = FE0->GetFEDesc3D();
    l = FEDesc0_Obj->GetN_JointDOF();

    cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
    for(j=0;j<N_Faces;j++) // loop over all Faces of cell
    {
      joint=cell->GetJoint(j);
      
      // Get count for Dirichlet nodes 
      // Set also in boundary node types and increment
      // Set vertex clipboards to 'i'      
      if(joint->GetType() == BoundaryFace ||
          joint->GetType() == IsoBoundFace)
      {
        // boundary joint
        BoundFace=(TBoundFace *)joint;
        BoundComp=BoundFace->GetBoundComp();
        Bd_ID = BoundComp->GetID();
	  
        N_Points = TmpLen[j];
//         for(k=0;k<N_Points;k++)
//         {
//           LinComb[k] = 1.0/N_Points;
//           cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
//         }
//         BoundFace->GetParameters(tp, sp);
// 
//         BoundComp->GetXYZandTS(N_Points, LinComb,
//                                xp, yp, zp, tp, sp,
//                                X, Y, Z, T, S);
       /** modified by sashi, since BD is not regular and T, S are not defined when tetgen is used */
        X = Y = Z = 0;
        for(k=0;k<N_Points;k++)
        {
          LinComb[k] = 1.0/N_Points;
          cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
          X += xp[k];
          Y += yp[k];
          Z += zp[k];
         }
        X /= N_Points;
        Y /= N_Points;
        Z /= N_Points;

        BoundaryCondition(Bd_ID, X, Y, Z, Cond0);
        switch(Cond0)
        {
          case DIRICHLET: // boundary nodes are Dirichlet nodes
               DirichletUpperBound += l;
            break;
          default:        // non Dirichlet boundary node
               BoundaryUpperBound[Cond0-1] += l;
  
        } // endswitch l
#ifdef _MPI         
        /** edges on this face are already set, so no need of checking in BD edges on this face */
        N_FaceEdges = ETmpLen[j];
        for(n=0;n<N_FaceEdges;n++)
         {
          edge = cell->GetEdge(TmpFE[j*EMaxLen+n]);
          edge->SetClipBoard(i);
         } // fo   
         
         /** vertices on this face are already set, so no need of checking in BD vert on this face, 30.06.12, sashi */        
        for(n=0;n<N_Points;n++)
          cell->GetVertex(TmpFV[j*MaxLen+n])->SetClipBoard(i);         
#endif        
      } // endif
      else
      {
        neigh = joint->GetNeighbour(cell);
        if(neigh)
        {
	    if(neigh->GetClipBoard() == -1 &&
		neigh->GetRefDesc()->GetType() == NoRef )
          {
	      if(joint->GetType() == InterfaceJoint3D ||
               joint->GetType() == IsoInterfaceJoint3D)
            {
              InterFace=(TInterfaceJoint3D *)joint;
              BoundComp=InterFace->GetBoundComp();
              Bd_ID = BoundComp->GetID();

              N_Points = TmpLen[j];
              for(k=0;k<N_Points;k++)
              {
                LinComb[k] = 1.0/N_Points;
                cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
              }
              InterFace->GetParameters(tp, sp);

              BoundComp->GetXYZandTS(N_Points, LinComb,
                                     xp, yp, zp, tp, sp,
                                     X, Y, Z, T, S);

              BoundaryCondition(Bd_ID, X, Y, Z, Cond0);

              switch(Cond0)
              {
                case DIRICHLET: 
                  // boundary nodes are Dirichlet nodes
                  DirichletUpperBound += l;
                break;

                default:
                  // non Dirichlet boundary node
                  BoundaryUpperBound[Cond0-1] += l;
              } // endswitch Cond0
            } // endif InterfaceJoint3D
#ifdef _MPI
	    else if(joint->GetType() == SubDomainHaloJoint)
            {
             // the neighbour is not a member of current collection
             Cond0 = SUBDOMAIN_HALOBOUND;
//              exit(0);
             BoundaryUpperBound[Cond0 - 1] += l; //  SUBDOMAIN_INTERFACE == 8, see constants.h
            }
            else if(joint->GetType() == SubDomainJoint)
            { 
             // the neighbour is not a member of current collection
             Cond0 = SUBDOMAIN_INTERFACE;
             BoundaryUpperBound[Cond0 - 1] += l; //  SUBDOMAIN_INTERFACE == 8, see constants.h
            }
#endif
            else
            {
              X = Y = Z = 0;
              N_Points = TmpLen[j];
              for(k=0;k<N_Points;k++)
              {
                cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
                X += xp[k];
                Y += yp[k];
                Z += zp[k];
              }
              X /= N_Points;
              Y /= N_Points;
              Z /= N_Points;

              BoundaryCondition(-1, X, Y, Z, Cond0);

              switch(Cond0)
              {
                case DIRICHLET: 
                  // boundary nodes are Dirichlet nodes
                  DirichletUpperBound += l;
                break;

                default:
                  // non Dirichlet boundary node
                  BoundaryUpperBound[Cond0-1] += l;
              } // endswitch Cond0
            } // end else
          } // end == -1
        } // endif neigh
      } // end no boundary
    } // endfor j  
    
#ifdef _MPI 
   /** identify BD edges but not part of BD faces */
  if( (l = FEDesc0_Obj->GetN_EdgeDOF()) > 0) //conforming FE
  {
   cell->GetShapeDesc()->GetEdgeVertex(EdgeVertex);   
    
   for(j=0;j<N_Edges;j++)
   {
    edge = cell->GetEdge(j);  
    
    if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)
    {
     if(edge->GetClipBoard()!=-1)
     continue;

     edge->SetClipBoard(1);

     cell->GetVertex(EdgeVertex[2*j])->GetCoords(xp[0], yp[0], zp[0]);
     cell->GetVertex(EdgeVertex[2*j+1])->GetCoords(xp[1], yp[1], zp[1]);
     X = (xp[0] + xp[1])/2.;
     Y = (yp[0] + yp[1])/2.;
     Z = (zp[0] + zp[1])/2.;   

//      if(rank==1)
//      {
//        cout<<rank << " xp[0] " << xp[0] << " yp[0] " << yp[0] << " zp[0] " << zp[0] <<endl;
//        cout<<rank << " xp[1] " << xp[1] << " yp[1] " << yp[1] << " zp[1] " << zp[1] <<endl;   
//      }
     
     BoundaryCondition(0,X, Y, Z, Cond0);

     switch(Cond0)
      {
       case DIRICHLET: 
        // boundary nodes are Dirichlet nodes
        DirichletUpperBound += l;
       break;

       default:
        // non Dirichlet boundary node
        BoundaryUpperBound[Cond0-1] += l;
       } // endswitch Cond0     

     /** vertices on this edge are already set, so no need of checking in BD vert on this edge, 30.06.12, sashi */             
     cell->GetVertex(EdgeVertex[2*j])->SetClipBoard(i);
     cell->GetVertex(EdgeVertex[2*j+1])->SetClipBoard(i); 
       
    } //  if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)    
   }//  for(j=0;j<N_Edges;j  
  }
   
  /** identify BD vert but not part of BD faces/edge */
  /** I hope, it is enough to check BD vertices, sice other vertices are anyway inner, so no prob with Dirichlet BD */
  if( (l=FEDesc0_Obj->GetN_VertDOF()) > 0) //conforming FE
   {
    for(j=0;j<N_VertInCell;j++)
     {
      Vert = cell->GetVertex(j);
      
      if(Vert->IsBoundVert() && Vert->GetClipBoard()!=-1 )
       {
        // BD vert is not yet set 
        Vert->SetClipBoard(i);

        Vert->GetCoords(X, Y, Z);       
        BoundaryCondition(0,X, Y, Z, Cond0);

        switch(Cond0)
         {
          case DIRICHLET: 
           // boundary nodes are Dirichlet nodes
           DirichletUpperBound += l;
          break;

          default:
           // non Dirichlet boundary node
          BoundaryUpperBound[Cond0-1] += l;
         } // endswitch Cond0     
      }     
    } //   for(j=0;j<N_VertInCell;j++)
   } //     if( FEDesc0_Obj->GetN_VertDOF() > 0) 
#endif   
  } // endfor i
   
   
 /*
  OutPut("DirichletUpperBound: " << DirichletUpperBound << endl);
  OutPut("N_DiffBoundNodeTypes: " << N_DiffBoundNodeTypes << endl);
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    OutPut("type[" << i << "]: " << BoundaryNodeTypes[i]);
    OutPut(" number of: " << BoundaryUpperBound[i] << endl);
  }
 */
  GlobalNumbers = new int[count];
  memset(GlobalNumbers, -1, SizeOfInt*count);
  SumLocDOF = count;

  // start DOF manager
  DirichletCounter=FIRSTMARK+1;
  l=-DirichletUpperBound + FIRSTMARK +1;


  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundCounter[i] = l;
    l -=  BoundaryUpperBound[i];
  }
  Counter = l;
  DirichletBound = Counter;
  
//   OutPut("DirichletCounter: " << DirichletCounter << endl);
//   for(i=0;i<N_DiffBoundNodeTypes;i++)
//     OutPut(i << "   " << BoundCounter[i] << endl);
//   OutPut("Counter: " << Counter << endl);

#ifdef _MPI
 for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i); 
    FEType0 = GetFE3D(i, cell);
    FE0 = TFEDatabase3D::GetFE3D(FEType0);
    FEDesc0_Obj = FE0->GetFEDesc3D();    
    
    if(cell->IsDependentCell())
     {       
      if(N_EdgeDOF < FEDesc0_Obj->GetN_EdgeDOF())
       N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF();

      if(N_VertDof < FEDesc0_Obj->GetN_VertDOF())
       N_VertDof = FEDesc0_Obj->GetN_VertDOF();

      N_Edges=cell->GetN_Edges();
      for(j=0;j<N_Edges;j++)
       (cell->GetEdge(j))->SetClipBoard(-1);
      
      N_VertInCell = cell->GetN_Vertices();
      for(j=0;j<N_VertInCell;j++)
       (cell->GetVertex(j))->SetClipBoard(-1);          
     }     
  }// for(i=0;i<N_Cells;i++) 
#endif 

  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    //OutPut("cell: " << i << endl);   
    
    N_Faces=cell->GetN_Joints();
    N_Edges=cell->GetN_Edges();
      
    FEType0 = GetFE3D(i, cell);
    FE0 = TFEDatabase3D::GetFE3D(FEType0);
    FE0->GetFEDesc3D(FEDesc0, FEDesc0_Obj);
    // FEDesc0 = FE0->GetFEDesc3D_ID();
    // FEDesc0_Obj = FE0->GetFEDesc3D();
    // FEDesc0_Obj = TFEDatabase3D::GetFEDesc3D(FEDesc0);
    I_K0 = BeginIndex[i];
    J_K0 = GlobalNumbers + BeginIndex[i];

    cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
    for(j=0;j<N_Faces;j++)
    {
      joint = cell->GetJoint(j);
      Indices0 = FEDesc0_Obj->GetJointDOF(j);

      if(joint->GetType() == BoundaryFace ||
         joint->GetType() == IsoBoundFace)
      {
        // boundary joint
        BoundFace=(TBoundFace *)joint;
        BoundComp=BoundFace->GetBoundComp();
        Bd_ID = BoundComp->GetID();
 
        N_Points = TmpLen[j];
//         for(k=0;k<N_Points;k++)
//         {
//           LinComb[k] = 1.0/N_Points;
//           cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
//         }
// 
//         BoundFace->GetParameters(tp, sp);
// 
//         BoundComp->GetXYZandTS(N_Points, LinComb,
//                                xp, yp, zp, tp, sp,
//                                X, Y, Z, T, S);
	
       /** modified by sashi, since BD is not regular and T, S are not properly defined when using tetgen */
        X = Y = Z = 0;
        for(k=0;k<N_Points;k++)
        {
          LinComb[k] = 1.0/N_Points;
          cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
          X += xp[k];
          Y += yp[k];
          Z += zp[k];
         }
        X /= N_Points;
        Y /= N_Points;
        Z /= N_Points;
 
        BoundaryCondition(Bd_ID, X, Y, Z, Cond0);

        switch(Cond0)
        {
          case DIRICHLET:
            // boundary nodes are Dirichlet nodes
            // OutPut("Dirichlet Boundary Face" << endl);

            mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
            mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                             DirichletCounter,
                             VHN, HNNumbers);
          break;

          default:
            // non Dirichlet boundary node
	      //OutPut("non Dirichlet Boundary Face" << endl);
            mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
            mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                             BoundCounter[Cond0-1],
                             VHN, HNNumbers);
	    //OutPut("done"<<endl);
        } // endswitch Cond0
#ifdef _MPI         
        /** edges on this face are already set, so no need to check in BD edges on this face */
        N_FaceEdges = ETmpLen[j];
        for(n=0;n<N_FaceEdges;n++)
         {
          edge = cell->GetEdge(TmpFE[j*EMaxLen+n]);
          edge->SetClipBoard(i);
         } // fo 
         
         /** vertices on this face are already set, so no need to check in BD vert on this face, 30.06.12, sashi */        
        for(n=0;n<N_Points;n++)
          cell->GetVertex(TmpFV[j*MaxLen+n])->SetClipBoard(i);
#endif         
      } // boundary joint
      else
      {
        // no boundary joint
         neigh = joint->GetNeighbour(cell);
        if (!neigh || joint->GetType() == MortarJoint ||
            joint->GetType() == MortarBaseJoint)
        {
          // there is no neighbour
          // => either mortar joint
          //    or finer cell in 1 regular grid
          //    will be handle from coarser cell
          if(joint->GetType() == MortarJoint ||
             joint->GetType() == MortarBaseJoint)
          {
            // do mortar mapping
            mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
            mapper->MapBound(GlobalNumbers, I_K0, Indices0, Counter,
                             VHN, HNNumbers);
          }
        } // !neigh
        else
        {
          // there is a neighbour on same level
          n = neigh->GetClipBoard();
          if( n == -1 )
          {
            // the neighbour is no member of current collection
            if(joint->GetType() == InterfaceJoint3D ||
               joint->GetType() == IsoInterfaceJoint3D)
            {
              InterFace=(TInterfaceJoint3D *)joint;
              BoundComp=InterFace->GetBoundComp();
              Bd_ID = BoundComp->GetID();

              N_Points = TmpLen[j];
              for(k=0;k<N_Points;k++)
              {
                LinComb[k] = 1.0/N_Points;
                cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
              }
              InterFace->GetParameters(tp, sp);

              BoundComp->GetXYZandTS(N_Points, LinComb,
                                     xp, yp, zp, tp, sp,
                                     X, Y, Z, T, S);

              BoundaryCondition(Bd_ID, X, Y, Z, Cond0);

              switch(Cond0)
              {
                case DIRICHLET:
                  // boundary nodes are Dirichlet nodes
                  // OutPut("Dirichlet Boundary Face" << endl);
                  mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
                  mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                   DirichletCounter,
                                   VHN, HNNumbers);
                break;

                default:
                  // non Dirichlet boundary node
                  // OutPut("non Dirichlet Boundary Face" << endl);
                  mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
                  mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                   BoundCounter[Cond0-1],
                                   VHN, HNNumbers);
              } // endswitch Cond0
            } // endif InterfaceJoint3D
#ifdef _MPI
	    else if(joint->GetType() == SubDomainHaloJoint)
            {
             Cond0 = SUBDOMAIN_HALOBOUND;
//              OutPut("Halo Boundary Edge" << endl);
//              exit(0);
             mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
             mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                              BoundCounter[Cond0-1],
                              VHN, HNNumbers);
            }
            else if(joint->GetType() == SubDomainJoint)
            {
             Cond0 = SUBDOMAIN_INTERFACE;
//              OutPut("SUBDOMAIN_INTERFACE Edge" << endl);
//             printf("SUBDOMAIN_INTERFACE  i %d j %d %d\n", i, j, neigh->GetClipBoard());
             mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
             mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                              BoundCounter[Cond0-1],
                              VHN, HNNumbers);
            }
#endif
            else
            {
              // find local edge of neigh on which cell is
              l=0;
              while(neigh->GetJoint(l) != joint) l++;
              // OutPut("l= " << l << endl);

              // check for children of neigh on its face l
              refdesc=neigh->GetRefDesc();
              if(refdesc->GetType() != NoRef)
              {
                refdesc->GetOldFaceNewFace(TmpoFnF, TmpLen1, MaxLen1);
                refdesc->GetFaceChild(TmpFC, TmpLen2, MaxLen2);
                refdesc->GetOldFaceNewLocFace(TmpoFnlF);
                refdesc->GetChildTwistIndex(TmpCTI);
                NFaces=refdesc->GetShapeDesc()->GetN_Joints();
  
                // OutPut("NFaces: " << NFaces << endl);
                // OutPut("TmpLen1[l]: " << TmpLen1[l] << endl);
                // OutPut("MaxLen1: " << MaxLen1 << endl);
  
                if(TmpLen1[l] == 4)
                {
                  f1=TmpoFnF[l*MaxLen1+0];
                  // OutPut("face1: " << f1 << endl);
                  chnum1=TmpFC[f1*MaxLen2];
                  child1=neigh->GetChild(chnum1);
                  c1=child1->GetClipBoard();
                  // OutPut("child: " << 1 << " " << c1 << endl);
                  FEType1 = GetFE3D(c1, child1);
                  FE1 = TFEDatabase3D::GetFE3D(FEType1);
                  FE1->GetFEDesc3D(FEDesc1, FEDesc1_Obj);
                  I_K1 = BeginIndex[c1];
                  J_K1 = GlobalNumbers + BeginIndex[c1];
                  m=TmpoFnlF[chnum1*NFaces+l];
                  Indices1 = FEDesc1_Obj->GetJointDOF(m);
                  Twist1 = TmpCTI[f1];
                  // OutPut("Twist1: " << Twist1 << endl);
  
                  f2=TmpoFnF[l*MaxLen1+1];
                  // OutPut("face2: " << f2 << endl);
                  chnum2=TmpFC[f2*MaxLen2];
                  child2=neigh->GetChild(chnum2);
                  c2=child2->GetClipBoard();
                  // OutPut("child: " << 2 << " " << c2 << endl);
                  I_K2 = BeginIndex[c2];
                  J_K2 = GlobalNumbers + BeginIndex[c2];
                  m=TmpoFnlF[chnum2*NFaces+l];
                  Indices2 = FEDesc1_Obj->GetJointDOF(m);
                  Twist2 = TmpCTI[f2];
                  // OutPut("Twist2: " << Twist2 << endl);
  
                  f3=TmpoFnF[l*MaxLen1+2];
                  // OutPut("face3: " << f3 << endl);
                  chnum3=TmpFC[f3*MaxLen2];
                  child3=neigh->GetChild(chnum3);
                  c3=child3->GetClipBoard();
                  // OutPut("child: " << 3 << " " << c3 << endl);
                  I_K3 = BeginIndex[c3];
                  J_K3 = GlobalNumbers + BeginIndex[c3];
                  m=TmpoFnlF[chnum3*NFaces+l];
                  Indices3 = FEDesc1_Obj->GetJointDOF(m);
                  Twist3 = TmpCTI[f3];
                  // OutPut("Twist3: " << Twist3 << endl);
  
                  f4=TmpoFnF[l*MaxLen1+3];
                  // OutPut("face4: " << f4 << endl);
                  chnum4=TmpFC[f4*MaxLen2];
                  child4=neigh->GetChild(chnum4);
                  c4=child4->GetClipBoard();
                  // OutPut("child: " << 4 << " " << c4 << endl);
                  I_K4 = BeginIndex[c4];
                  J_K4 = GlobalNumbers + BeginIndex[c4];
                  m=TmpoFnlF[chnum4*NFaces+l];
                  Indices4 = FEDesc1_Obj->GetJointDOF(m);
                  Twist4 = TmpCTI[f4];
                  // OutPut("Twist4: " << Twist4 << endl);
  
                  // OutPut("big: " << i << endl);
                  // OutPut("children: " << c1 << " " << c2 << " " << c3);
                  // OutPut(" " << c4 << endl);
  
                  mapper1reg = TFEDatabase3D::GetFE3DMapper1Reg(FEDesc0, FEDesc1);
                  mapper1reg->Map(GlobalNumbers, I_K1, I_K2, I_K3, I_K4,
                                  I_K0, Indices1, Indices2, Indices3,
                                  Indices4, Indices0,
                                  Twist1, Twist2, Twist3, Twist4,
                                  joint->GetMapType(),
                                  DirichletBound,
                                  Counter, VHN, HNNumbers);
                }
                // There is exactly one children on the neighbour sharing the current face
                else if(TmpLen1[l] == 1)
                {
                  // Get Child of Neighbour that contains this face
                  f1 = TmpoFnF[l * MaxLen1];
                  chnum1 = TmpFC[f1 * MaxLen2];
                  child1 = neigh->GetChild(chnum1);
                  c1 = child1->GetClipBoard();
                  
                  FEType1 = GetFE3D(c1, child1);
                  FE1 = TFEDatabase3D::GetFE3D(FEType1);
                  FE1->GetFEDesc3D(FEDesc1, FEDesc1_Obj);
                  
                  I_K1 = BeginIndex[c1];
                  J_K1 = GlobalNumbers + BeginIndex[c1];
                  m = TmpoFnlF[chnum1 * NFaces + l];
                  if(m == -1)
                  {
                    std::cerr << "Error!\n";
                    exit(-1);
                  }
                  Indices1 = FEDesc1_Obj->GetJointDOF(m);
                  
                  mapper = TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc1);
                  mapper->Map(joint->GetMapType(), GlobalNumbers, I_K0, I_K1,
                              Indices0, Indices1, 0, 0, FEDesc0_Obj,
                              FEDesc1_Obj, Counter, VHN, HNNumbers);
                }
                else
                {
                  Error("Only regular refinement are allowed!!" << endl);
                  exit(-1);
                }
              }
              else
              {
                // neighbour is not refined
                N_Points = TmpLen[j];
                X = Y = Z = 0.0;
                for(k=0;k<N_Points;k++)
                {
                  cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
                  X += xp[k];
                  Y += yp[k];
                  Z += zp[k];
                }
                X /= N_Points;
                Y /= N_Points;
                Z /= N_Points;
  
                BoundaryCondition(0, X, Y, Z, Cond0);
  
                switch(Cond0)
                {
                  case DIRICHLET:
                    // nodes are Dirichlet nodes
                    // OutPut("Dirichlet inner" << endl);
                    mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
                    mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                     DirichletCounter,
                                     VHN, HNNumbers);
                  break;

                  default:
                    // non Dirichlet node
                    // OutPut("non Dirichlet inner" << endl);
                    mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
                    mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                     BoundCounter[Cond0-1],
                                     VHN, HNNumbers);
                } // endswitch Cond0
              }
            }
          } // n == -1
          else
          {
	      // neighbour is member of this collection
            // => using mappers
            if (n>i)
            {
              // this joint was not handled until now

              FEType1 = GetFE3D(n, neigh);
              FE1 = TFEDatabase3D::GetFE3D(FEType1);
              FE1->GetFEDesc3D(FEDesc1, FEDesc1_Obj);
              // FEDesc1 = FE1->GetFEDesc3D_ID();
              // FEDesc1_Obj = FE1->GetFEDesc3D();
              // FEDesc1_Obj = TFEDatabase3D::GetFEDesc3D(FEDesc1);
              I_K1 = BeginIndex[n];
              J_K1 = GlobalNumbers + BeginIndex[n];
              // Indices1 = FEDesc1_Obj->GetJointDOF(l);

              // find the local edge of neigh on which cell is
              l=0;
              // while(neigh->GetJoint(l)->GetNeighbour(neigh)!=cell) l++;
              while(neigh->GetJoint(l) != joint) l++;

              Indices1 = FEDesc1_Obj->GetJointDOF(l);

              // OutPut("i= " << i << endl);
              // OutPut("j= " << j << endl);
              // OutPut("n= " << n << endl);
              // OutPut("l= " << l << endl);
              // OutPut("-------------" << endl);

              //OutPut("FEDesc0: " << FEDesc0 << endl);
              //OutPut("FEDesc1: " << FEDesc1 << endl);
              mapper = TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc1);
              //OutPut("MapType: " << joint->GetMapType() << endl);
              mapper->Map(joint->GetMapType(),
                          GlobalNumbers, I_K0, I_K1, Indices0, Indices1,
                          j, l, FEDesc0_Obj, FEDesc1_Obj, Counter,
                          VHN, HNNumbers);

	    } // n>i

          } // n != -1
        } // end neigh
      } // no boundary joint
//   if(rank==7){
//     for(jj=0;jj<BeginIndex[i+1];jj++)
//      printf("FESpace3D I_K0 i %d j %d global %d\n",  i, jj, GlobalNumbers[jj]);
//     printf("\n");
// }

    } // endfor j
 
 
#ifdef _MPI   
  if( (N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF())>0)
  {
   cell->GetShapeDesc()->GetEdgeVertex(EdgeVertex); 

   for(j=0;j<N_Edges;j++)
   {
    edge = cell->GetEdge(j);  
    
    if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)
    {
     if(edge->GetClipBoard()!=-1)
     continue;
    
     Indices0 = FEDesc0_Obj->GetEdgeDOF(j);   
     edge->SetClipBoard(i);

     cell->GetVertex(EdgeVertex[2*j])->GetCoords(xp[0], yp[0], zp[0]);
     cell->GetVertex(EdgeVertex[2*j+1])->GetCoords(xp[1], yp[1], zp[1]);
     X = (xp[0] + xp[1])/2.;
     Y = (yp[0] + yp[1])/2.;
     Z = (zp[0] + zp[1])/2.;   
    
     BoundaryCondition(0,X, Y, Z, Cond0);

     switch(Cond0)
      {
       case DIRICHLET:
         // boundary nodes are Dirichlet nodes
         // OutPut("Dirichlet Boundary Face" << endl);

         mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
         mapper->MapBoundEdge(N_EdgeDOF, GlobalNumbers, I_K0, Indices0,
                              DirichletCounter,
                              VHN, HNNumbers);
          break;

          default:
            // non Dirichlet boundary node
            //OutPut("non Dirichlet Boundary Face" << endl);
            mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);
            mapper->MapBoundEdge(N_EdgeDOF, GlobalNumbers, I_K0, Indices0,
                                 BoundCounter[Cond0-1],
                                 VHN, HNNumbers);
          //OutPut("done"<<endl);
        } // endswitch Cond0  
        
     /** vertices on this edge are already set, so no need of checking in BD vert on this edge, 30.06.12, sashi */             
     cell->GetVertex(EdgeVertex[2*j])->SetClipBoard(i);
     cell->GetVertex(EdgeVertex[2*j+1])->SetClipBoard(i); 
    }    // if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)
   }//  for(j=0;j<N_Edges;j 
  } // if(N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF()>0)
  
  /** identify BD vert but not part of BD faces/edge */
  /** I hope, it is enough to check BD vertices, since other vertices are anyway inner, so no prob with Dirichlet BD */
  if(FEDesc0_Obj->GetN_VertDOF() > 0)
   {    
    for(j=0;j<N_VertInCell;j++)
     {
      Vert = cell->GetVertex(j);
      
      if(Vert->IsBoundVert() && Vert->GetClipBoard()==-1 )
       {   
        // BD vert is not yet set 
        Vert->SetClipBoard(i);

        Vert->GetCoords(X, Y, Z);       
        BoundaryCondition(0,X, Y, Z, Cond0);
        mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0); 

        switch(Cond0)
         {
          case DIRICHLET: 
           // boundary nodes are Dirichlet nodes
           mapper->MapBoundVert(GlobalNumbers, I_K0, FEDesc0_Obj->GetVertDOF(j),
                                DirichletCounter, VHN, HNNumbers);    
    
          break;

          default:
           // non Dirichlet boundary node
           mapper->MapBoundVert(GlobalNumbers, I_K0, FEDesc0_Obj->GetVertDOF(j),
                                BoundCounter[Cond0-1], VHN, HNNumbers);    
         } // endswitch Cond0       
       }  //   if(Vert->IsBoundVert() && Vert->GetClipBoard()==-1 )
     } //   for(j=0;j<N_VertInCell;j++)
   } //    if(cell->IsBoundCell(
#endif    
 
    // handle inner degrees of freedom
    k = FEDesc0_Obj->GetN_InnerDOF();
    Indices0 = FEDesc0_Obj->GetInnerDOF();
    for(j=0;j<k;j++)
    {
     Counter--;
     J_K0[Indices0[j]] = Counter;
    } // endfor j 
  } // endfor i

// ===============================================================================================
// When the domain is partioned, cells in same subdomain may have only an edge(s) or vertex(ies) as common 
// between other cells in the same subcoll, therefore just joint map will not give a unique 
// globaldof on such edges or vertex
//  - sashikumaar Ganesan
// ===============================================================================================


#ifdef _MPI
 for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i); 
    FEType0 = GetFE3D(i, cell);
    FE0 = TFEDatabase3D::GetFE3D(FEType0);
    FEDesc0_Obj = FE0->GetFEDesc3D();    
    
    if(cell->IsDependentCell())
     {       
      if(N_EdgeDOF < FEDesc0_Obj->GetN_EdgeDOF())
       N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF();

      if(N_VertDof < FEDesc0_Obj->GetN_VertDOF())
       N_VertDof = FEDesc0_Obj->GetN_VertDOF();

      N_Edges=cell->GetN_Edges();

      for(j=0;j<N_Edges;j++)
       (cell->GetEdge(j))->SetClipBoard(-1);
     }     
  }// for(i=0;i<N_Cells;i++) 
#endif 

#ifdef _MPI
  if(N_EdgeDOF>0) // only for cont. FESpace
   {
    int *w_array = new int[N_EdgeDOF];

    for(i=0;i<N_Cells;i++)
    {
     cell = Collection->GetCell(i);

     if(cell->IsDependentCell())
     {

      N_VertInCell = cell->GetN_Vertices();

      for(j=0;j<N_VertInCell;j++)
       (cell->GetVertex(j))->SetClipBoard(-1);

      I_K0 = BeginIndex[i];
      FEType0 = GetFE3D(i, cell);
      FE0 = TFEDatabase3D::GetFE3D(FEType0);
      FE0->GetFEDesc3D(FEDesc0, FEDesc0_Obj);
      N_Edges=cell->GetN_Edges();

      N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF();
      (cell->GetShapeDesc())->GetEdgeVertex(EdgeVertex);

      for(j=0;j<N_Edges;j++)
      {
       Vert = cell->GetVertex(EdgeVertex[2*j]); // start vertex of the edge
       edge = cell->GetEdge(j);

       if(edge->GetClipBoard()!=-1)
          continue;

        edge->SetClipBoard(1);
        EdgeDof = FEDesc0_Obj->GetEdgeDOF(j);
     
        /** fill the already assigned global dof */
        for(k=0;k<N_EdgeDOF;k++)
         {
          w_array[k]= I_K0 +  EdgeDof[k];  
          while((v0 = GlobalNumbers[w_array[k]])>-1 )
          { w_array[k] = v0; }
         } // for(k=0;k<N_EdgeDOF;k+

        edge->GetNeibs(N_EdgeNeibs, EdgeNeibs);

        for(k=0;k<N_EdgeNeibs;k++)
         {
          neigh = EdgeNeibs[k];

          /** hallo cell included, as multigrid fespace will contain Halo cells - Sashi:10-03-15 */
          if( (rank!=(neigh->GetSubDomainNo())) && !(neigh->IsHaloCell()) ) // if neigh is not in this domain and it is also not a halo cell
           continue;
  
          /** own cell or neib cell (i.e. Hallo cells) is not in this collection */
          if(neigh==cell || neigh->GetClipBoard_Par()==-1 )
           continue;
  
         // either a halo cell or a cell in this domain(owncell??) and in the collection
         l=0; // find the edge local index in the neib cell
         while(edge!= neigh->GetEdge(l) ) l++;

         (neigh->GetShapeDesc())->GetEdgeVertex(NeibEdgeVertex);

         if(neigh->GetVertex(NeibEdgeVertex[2*l])==Vert)
           { maptype =  1; }
         else if(neigh->GetVertex(NeibEdgeVertex[2*l+1])==Vert)
          { maptype =  -1; }
         else
          {
           printf("FESpace3D Error in finding cross edge maptype Edge \n");
           MPI_Abort(MPI_COMM_WORLD, 0);
          }

         n = neigh->GetClipBoard();
         FEType1 = GetFE3D(n, neigh);
         FE1 = TFEDatabase3D::GetFE3D(FEType1);
         FE1->GetFEDesc3D(FEDesc1, FEDesc1_Obj);
         NeibEdgeDof = FEDesc1_Obj->GetEdgeDOF(l);

         for(m=0;m<N_EdgeDOF;m++)
          {
           if(maptype==1)
            {neibdof = BeginIndex[n] + NeibEdgeDof[m]; }
           else
            {neibdof = BeginIndex[n] + NeibEdgeDof[N_EdgeDOF-1 - m]; }

           w1 = neibdof;
           while( (v1=GlobalNumbers[w1])>-1 )
            { w1=v1; }

           w0 = w_array[m];

/*      if(neigh->IsHaloCell()  )
     {
      #ifdef _MPI        
      printf("Test IISc  rank %d:, count %d , w0 %d, w1 %d\n", rank, m,k, w1); 
      MPI_Finalize();
      #endif      
      exit(0);   
     }  */ 
   
           if(GlobalNumbers[w0] != GlobalNumbers[w1])
            {
//       if(TDatabase::ParamDB->PBE_P9>0)
//      {
//       #ifdef _MPI        
//       printf("Test IISc  rank %d:, count %d , w0 %d, w1 %d\n", rank, m,k, w1); 
//       MPI_Finalize();
//       #endif      
//       exit(0);   
//      }   

// if(rank== TDatabase::ParamDB->Par_P5)
// if(rank== 7)
// {     
//   cout << endl;
//   cout <<"Rank: " << rank << " Cell: " << i <<" Edge: " << j<< " NeibCell: " << n <<" NeibEdge: " << l << endl;
//   cout << "owndof: " << owndof << endl;
//   cout << "Global[owndof]: " << GlobalNumbers[owndof] << endl;
//   cout << "neibdof: " << neibdof << endl;
//   cout << "Global[neibdof]: " << GlobalNumbers[neibdof] << endl;
//   cout << "w0: " << w0 << endl;
//   cout << "Global[w0]: " << GlobalNumbers[w0] << endl;
//   cout << "w1: " << w1 << endl;
//   cout << "Global[w1]: " << GlobalNumbers[w1] << endl;
//   }

            e = MAX( GlobalNumbers[w0], GlobalNumbers[w1]);
            w_array[m] = MIN(w0,w1);

            GlobalNumbers[w0] =  w_array[m]; 
            GlobalNumbers[w1] =  w_array[m];
            GlobalNumbers[w_array[m]] = e;
           }

         } // for(m=0;m<N_EdgeDOF;m++
        } //  for(k=0;k<N_EdgeNeibs;k++)
      } //  for(j=0;j<N_Edges;j+
     } // if(cell->IsDependentCell())
   } // for(i=0;i<N_Cells;i++)


     
   delete [] w_array;
 }//  if(N_EdgeDOF>0)  only for cont. FESpace


// #ifdef _MPI      
//    if(rank==0)
//     printf("FESpace3D    %d \n",  rank );  
// //     MPI_Finalize(); 
// #endif    
// //   exit(0); 
  
  if(N_VertDof>0)// only for cont. FESpace
  {
    for(i=0;i<N_Cells;i++)
    {
     cell = Collection->GetCell(i);

     if(cell->IsDependentCell())
     {
      I_K0 = BeginIndex[i];
      FEType0 = GetFE3D(i, cell);
      FE0 = TFEDatabase3D::GetFE3D(FEType0);
      FE0->GetFEDesc3D(FEDesc0, FEDesc0_Obj);
      N_VertInCell = cell->GetN_Vertices();

      N_VertDof = FEDesc0_Obj->GetN_VertDOF();

      for(j=0;j<N_VertInCell;j++)
       {
        Vert=cell->GetVertex(j);

        if( Vert->GetClipBoard()!=-1) 
         continue;
 
        Vert->SetClipBoard(5);
        w0 = I_K0+FEDesc0_Obj->GetVertDOF(j);

        while((v0 = GlobalNumbers[w0])>-1 )
          { w0 = v0; }

        Vert->GetNeibs(N_VertNeibs, VertNeibs);
   
        for(k=0;k<N_VertNeibs;k++)
         {
          neigh = VertNeibs[k];

          /** hallo cells are also need to be considered, so modified 11 Mar 2015 by Sashi */
          if( (rank!=neigh->GetSubDomainNo() ) && !(neigh->IsHaloCell()) )
           {  continue; }

          /** own cell or neib cell (i.e. Hallo cells) is not in this collection */
          if(neigh==cell || neigh->GetClipBoard_Par()==-1 )
           continue; 

          n = neigh->GetClipBoard();
          FEType1 = GetFE3D(n, neigh);
          FE1 = TFEDatabase3D::GetFE3D(FEType1);
          FE1->GetFEDesc3D(FEDesc1, FEDesc1_Obj);

          l=0; // find the vert local index in the neib cell
          while(Vert!= neigh->GetVertex(l) ) l++;

          w1 = BeginIndex[n] + FEDesc1_Obj->GetVertDOF(l);

          while((v1=GlobalNumbers[w1])>-1)
           { w1=v1; }

         if( GlobalNumbers[w0] != GlobalNumbers[w1])
          {

// if(rank== TDatabase::ParamDB->Par_P5)
// {     
//   cout << endl;
//   cout <<"Rank: " << rank << " Cell: " << i <<" Edge: " << j<< " NeibCell: " << n <<" NeibEdge: " << l << endl;
//   cout << "w0: " << w0 << endl;
//   cout << "Global[w0]: " << GlobalNumbers[w0] << endl;
//   cout << "w1: " << w1 << endl;
//   cout << "Global[w1]: " << GlobalNumbers[w1] << endl;
//   }
            e = MAX( GlobalNumbers[w0], GlobalNumbers[w1]);
            w = MIN(w0,w1);

            GlobalNumbers[w0] =  w; 
            GlobalNumbers[w1] =  w;
            GlobalNumbers[w] = e;
            w0 = w;
           }

         } //for(k=0;k<N_VertNeibs;k++)
      } //  for(j=0;j<N_VertInCell;j++)
     } // if(cell->IsDependentCell())
   } // for(i=0;i<N_Cells;i++)
  }//  if(N_VertDof>0)

#endif
//  printf("Rank %d : MPI_Finalize test :: %d\n\n",rank, N_Cells);   
//   MPI_Finalize();  
//   exit(0);
//   if(rank==7){
//     for(jj=0;jj<BeginIndex[1];jj++)
//      printf("FESpace3D I_K0 i %d j %d global %d\n",  i, jj, GlobalNumbers[jj]);
//     printf("\n");
// }

//   if(rank==7)
//    printf("FESpace3D !!! %d\n", FIRSTMARK);
//  MPI_Finalize();
//  exit(0);

  // find global numbers
  l=0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
    l += BoundaryUpperBound[i];

  m = -SumLocDOF - DirichletUpperBound - l +FIRSTMARK;
  n = VHN->GetN_Elements();
  l = 0;
  //OutPut("HN" << endl);
  for(i=0;i<n;i++)
  {
    j=HNNumbers->GetElement(i);
    // OutPut(j << endl);
    if(j<0) continue;

    l++;

    while( (k=GlobalNumbers[j]) > -1 )
      j=k;

    GlobalNumbers[j] = m;
    m--;
  }
  // OutPut("number of hanging nodes: " << l << endl);

  l = FIRSTMARK - DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundMark[i] = l;
    l -= BoundaryUpperBound[i];
  }
  InnerMark = l;
  SlaveMark = l - SumLocDOF;

  //for(i=0;i<N_DiffBoundNodeTypes;i++)
  //    OutPut(i << "   " << BoundMark[i] << endl);
  // OutPut("InnerMark: " << InnerMark << endl);
  // OutPut("SlaveMark: " << SlaveMark << endl);

  DirichletCounter=0;
  l=DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundCounter[i] = l;
    l += BoundaryUpperBound[i];
  }
  count = l;

  N_Dirichlet = 0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
    N_BoundaryNodes[i]=0;
  N_Inner = 0;
  N_Slave = 0;

  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    // OutPut(endl << "cell: " << i << endl);

    FEType0 = GetFE3D(i, cell);
    FE0 = TFEDatabase3D::GetFE3D(FEType0);
    FE0->GetFEDesc3D(FEDesc0, FEDesc0_Obj);
    // FEDesc0 = FE0->GetFEDesc3D_ID();
    // FEDesc0_Obj = FE0->GetFEDesc3D();
    // FEDesc0_Obj = TFEDatabase3D::GetFEDesc3D(FEDesc0);
    J_K0 = GlobalNumbers + BeginIndex[i];

    k=FEDesc0_Obj->GetN_DOF();
    for(j=0;j<k;j++)
    {
      l = J_K0[j];
      if (l < -1)
      {
        // OutPut(endl << "new node" << endl);
        if(l<=SlaveMark)
        {
          // hanging node
          J_K0[j] = -l + SumLocDOF;
          N_Slave++;
          // OutPut("slave" << endl);
          continue;
        }

        if(l<=InnerMark)
        {
          // inner node
          J_K0[j] = count;
          count++;
          N_Inner++;
          // OutPut("inner" << endl);
          continue;
        }

        for(m=N_DiffBoundNodeTypes-1;m>=0;m--)
        {
          if(l<=BoundMark[m])
          {
            J_K0[j] = BoundCounter[m];
            BoundCounter[m]++;
            N_BoundaryNodes[m]++;
            // OutPut("type " << m << endl);
            m = -2;
            break;
          }
        }

        if(m!=-2 && l<=FIRSTMARK) // no match in loop above
        {
          // Dirichlet nodes
          J_K0[j] = DirichletCounter;
          DirichletCounter++;
          N_Dirichlet++;
          // OutPut("Dirichlet" << endl);
          continue;
        }

      } // l < -1
      else
      {
        if (l >= 0)
        {
          J_K0[j] = GlobalNumbers[l];
        }
        else
        {
          OutPut("J_K0[j]==-1 at locdof: " << j << " in ele: " << i << endl);
        }
      } // l >= -1
    } // endfor j

  } // endfor i

/*
  OutPut("N_Inner: " << N_Inner << endl);
  OutPut("N_Slave: " << N_Slave << endl);
  OutPut("N_Dirichlet: " << N_Dirichlet << endl);
  for(i=0;i<N_DiffBoundNodeTypes;i++)
    OutPut(i << " N_BoundaryNodes: " << N_BoundaryNodes[i] << endl);
*/

  // create real numbers
  l = 0; m = 0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += N_BoundaryNodes[i];
    m += BoundaryUpperBound[i];
  }

  DirichletOffset = N_Inner + N_Slave + l - 0;
  // OutPut("DirichletOffset: " << DirichletOffset << endl);
  SlaveOffset = N_Inner + l - 2*SumLocDOF - DirichletUpperBound 
                - m + FIRSTMARK;
  InnerOffset = 0 - DirichletUpperBound - m;

  l = N_Inner; m = DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundOffset[i] = l-m;
    l += N_BoundaryNodes[i];
    m += BoundaryUpperBound[i];
  }

  DirichletMark = DirichletUpperBound;
  InnerMark = 2*SumLocDOF-FIRSTMARK;
  l = DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += BoundaryUpperBound[i];
    BoundMark[i] = l;
  }

  InnerMark = l + N_Inner+1; // Change made by ravi
  // OutPut("DirichletMark: " << DirichletMark << endl);
  // OutPut("InnerMark: " << InnerMark << endl);

  for(i=0;i<SumLocDOF;i++) //Looping over all dofs in all cells in the collection only
  {
    n=GlobalNumbers[i];
    // OutPut(i << "  " << n << endl);
    if(n<DirichletMark)
    {
      // Dirichlet node
      GlobalNumbers[i] += DirichletOffset;
      // OutPut("Diri" << GlobalNumbers[i] << endl);
    }
    else
    {
      if(n<BoundMark[N_DiffBoundNodeTypes-1])
      {
        // non Dirichlet boundary type
        for(m=0;m<N_DiffBoundNodeTypes;m++)
        {
          if(n<BoundMark[m])
          {
            // node of type m
            GlobalNumbers[i] += BoundOffset[m];
            // OutPut("type: " << m << endl);
            break;
          }
        } // endfor m
      }
      else
      {
        // no boundary node
        if(n<InnerMark)
        {
          // inner node
          GlobalNumbers[i] += InnerOffset;
          // OutPut("inner" << endl);
        }
        else
        {
          // slave node
          GlobalNumbers[i] += SlaveOffset;
          // OutPut("slave: " << endl);
        }
      }
    } // non Dirichlet node
  } // endfor i

/*
  // print for all elements for global numbers of their local dofs
  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    OutPut("cell number: " << i << endl);

    J_K0 = GlobalNumbers + BeginIndex[i];
    FEType0 = GetFE3D(i, cell);
    FE0 = TFEDatabase3D::GetFE3D(FEType0);
    k = FE0->GetSize();
    for(j=0;j<k;j++)
    {
      OutPut(j << ": " << " number: " << J_K0[j] << endl);
    }

    OutPut(endl);

  } // endfor i
*/

  // fill in information for hanging nodes
  HangingNodeArray = new THangingNode*[N_Slave];

  m = VHN->GetN_Elements();
  n = 0;
  for(i=0;i<m;i++)
  {
    hn=VHN->GetElement(i);
    l=HNNumbers->GetElement(i);
    if(l<0) continue;
    // OutPut("number: " << l << endl);
    k=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes();
    v=hn->GetDOF();
    // OutPut("HN: " << i << " ");
    // OutPut("number of nodes in coupling: " << k << endl);
    t0 = 0;
    for(j=0;j<k;j++)
    {
      // OutPut(j << ": ");
      v[j] = GlobalNumbers[v[j]];
      // OutPut(v[j] << endl);
      t0 += v[j];
    }
    // OutPut("sum: " << t0 << endl);

    HangingNodeArray[n] = hn;
    n++;
  }
  delete VHN;
  delete HNNumbers;

  InnerBound = N_Inner;
  l = N_Inner;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += N_BoundaryNodes[i];
    BoundaryNodesBound[i] = l;
  }
  HangingBound = l+N_Slave;
  DirichletBound = HangingBound + N_Dirichlet;

  N_ActiveDegrees = l;
  N_SlaveDegrees = N_Slave;
  N_DegreesOfFreedom = N_ActiveDegrees + N_Slave + N_Dirichlet;

/*
  OutPut("N_ActiveDegrees: " << N_ActiveDegrees << endl);
  OutPut("N_SlaveDegrees: " << N_SlaveDegrees << endl);
  OutPut("N_DegreesOfFreedom: " << N_DegreesOfFreedom << endl << endl);
*/

  ActiveBound = N_ActiveDegrees;

  delete BoundaryUpperBound;
  delete BoundCounter;
  delete BoundMark;
  delete BoundOffset;
  
//   #ifdef _MPI
//   int rank, size;
//   MPI_Comm comm = MPI_COMM_WORLD;
  
//   MPI_Comm_size(comm, &size);
//   MPI_Comm_rank(comm, &rank);

//   BEGIN_SEQ
//   OutPut("SUBDOMAIN_INTERFACE: " << N_BoundaryNodes[SUBDOMAIN_INTERFACE-1] << endl);
//   OutPut("SUBDOMAIN_HALOBOUND: " << N_BoundaryNodes[SUBDOMAIN_HALOBOUND-1] << endl);
//   END_SEQ
//   #endif
}

TFESpace3D::~TFESpace3D()
{
  delete BoundaryNodesBound;

  if (UsedElements)
    delete UsedElements;

  if(AllElements)
    delete AllElements;

  if(ElementForShape)
    delete ElementForShape;
}

/** return position of all dofs */
void TFESpace3D::GetDOFPosition(double *x, double *y, double *z)
{
  int i,j,k;
  TBaseCell *cell;
  int N_Joints;
  TJoint *joint;
  JointType jointtype;
  FE3D FEid;
  int *DOF;
  TNodalFunctional3D *nf;
  double *xi, *eta, *zeta;
  int N_Points;
  RefTrans3D RefTrans, *RefTransArray;
  int IsIsoparametric;
  BoundTypes bdtype;
  BF3DRefElements RefElement;
  TRefTrans3D *rt;
  double X[MaxN_BaseFunctions3D], Y[MaxN_BaseFunctions3D];
  double Z[MaxN_BaseFunctions3D], absdetjk[MaxN_BaseFunctions3D];

  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();

  for(i=0;i<N_Cells;i++)
  {
    DOF = GlobalNumbers + BeginIndex[i];

    cell  = Collection->GetCell(i);
    FEid = GetFE3D(i, cell); 
    RefTrans = RefTransArray[FEid];

    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEid);

    nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(FEid);
    nf->GetPointsForAll(N_Points, xi, eta, zeta);

    N_Joints = cell->GetN_Joints();

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Joints;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryFace)
        {
          bdtype = ((TBoundFace *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Plane)
            IsIsoparametric = TRUE;
        }
        if(jointtype == InterfaceJoint3D)
        {
          bdtype = ((TInterfaceJoint3D *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Plane)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint3D)
          IsIsoparametric = TRUE;
  
        if(jointtype == IsoJointEqN)
          IsIsoparametric = TRUE;
  
        if(jointtype == IsoBoundFace)
          IsIsoparametric = TRUE;
      }
    }// endif

    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          RefTrans = HexaIsoparametric;
        break;
  
        case BFUnitTetrahedron:
          RefTrans = TetraIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    rt = TFEDatabase3D::GetRefTrans3D(RefTrans);
    switch(RefTrans)
    {
      case TetraAffin:
       ((TTetraAffin *)rt)->SetCell(cell);
       ((TTetraAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                                X, Y, Z, absdetjk);
      break;

      case TetraIsoparametric:
       ((TTetraIsoparametric *)rt)->SetCell(cell);
       ((TTetraIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                                X, Y, Z, absdetjk);
      break;

      case HexaAffin:
       ((THexaAffin *)rt)->SetCell(cell);
       ((THexaAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                                X, Y, Z, absdetjk);
      break;

      case HexaTrilinear:
       ((THexaTrilinear *)rt)->SetCell(cell);
       ((THexaTrilinear *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                                X, Y, Z, absdetjk);
      break;

      case HexaIsoparametric:
       ((THexaIsoparametric *)rt)->SetCell(cell);
       ((THexaIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                                X, Y, Z, absdetjk);
      break;
    } // endswitch RefTrans

    for(j=0;j<N_Points;j++)
    {
      k = DOF[j];
      x[k] = X[j];
      y[k] = Y[j];
      z[k] = Z[j];
    }

  } // endfor i
} // end GetDOFPosition

/** return position of all dofs */
void TFESpace3D::GetDOFPosition(int dof, double &x, double &y, double &z)
{
  int i,j,k;
  TBaseCell *cell;
  int N_Joints;
  TJoint *joint;
  JointType jointtype;
  FE3D FEid;
  int *DOF;
  TNodalFunctional3D *nf;
  double *xi, *eta, *zeta;
  int N_Points;
  RefTrans3D RefTrans, *RefTransArray;
  int IsIsoparametric;
  BoundTypes bdtype;
  BF3DRefElements RefElement;
  TRefTrans3D *rt;
  double absdetjk[1];
  int DOFFound;

  if(dof > N_DegreesOfFreedom)
  {
    Error(dof << " dof number is larger than total number of degrees of freedom" << endl);
    x = -1; y = -1; z = -1;
  }
  
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();

  for(i=0;i<N_Cells;i++)
  {
    DOF = GlobalNumbers + BeginIndex[i];
    k = BeginIndex[i+1] - BeginIndex[i];

    DOFFound = -1;
    for(j=0;j<k;j++)
    {
      if(DOF[j] == dof) 
      {
        DOFFound = j;
        break;
      } // endif
    } // endfor


    if(DOFFound>-1) // i.e. dof was found
    {
      //cout << "dof " << dof << " found in cell: " << i << endl;
      cell  = Collection->GetCell(i);
      FEid = GetFE3D(i, cell); 
      RefTrans = RefTransArray[FEid];
  
      RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEid);
  
      nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(FEid);
      nf->GetPointsForAll(N_Points, xi, eta, zeta);
  
      N_Joints = cell->GetN_Joints();
  
      IsIsoparametric = FALSE;
      if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
      {
        for(j=0;j<N_Joints;j++)
        {
          joint = cell->GetJoint(j);
          jointtype = joint->GetType();
          if(jointtype == BoundaryFace)
          {
            bdtype = ((TBoundFace *)(joint))->GetBoundComp()->GetType();
            if(bdtype != Plane)
              IsIsoparametric = TRUE;
          }
          if(jointtype == InterfaceJoint3D)
          {
            bdtype = ((TInterfaceJoint3D *)(joint))->GetBoundComp()->GetType();
            if(bdtype != Plane)
              IsIsoparametric = TRUE;
          }
          if(jointtype == IsoInterfaceJoint3D)
            IsIsoparametric = TRUE;
    
          if(jointtype == IsoJointEqN)
            IsIsoparametric = TRUE;
    
          if(jointtype == IsoBoundFace)
            IsIsoparametric = TRUE;
        }
      }// endif
  
      if(IsIsoparametric)
      {
        switch(RefElement)
        {
          case BFUnitHexahedron:
            RefTrans = HexaIsoparametric;
          break;
    
          case BFUnitTetrahedron:
            RefTrans = TetraIsoparametric;
          break;
        }
      } // endif IsIsoparametric
  
      rt = TFEDatabase3D::GetRefTrans3D(RefTrans);
      switch(RefTrans)
      {
        case TetraAffin:
         ((TTetraAffin *)rt)->SetCell(cell);
         ((TTetraAffin *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound, zeta+DOFFound,
                                                &x, &y, &z, absdetjk);
        break;
  
        case TetraIsoparametric:
         ((TTetraIsoparametric *)rt)->SetCell(cell);
         ((TTetraIsoparametric *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound, zeta+DOFFound,
                                                  &x, &y, &z, absdetjk);
        break;
  
        case HexaAffin:
         ((THexaAffin *)rt)->SetCell(cell);
         ((THexaAffin *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound, zeta+DOFFound,
                                                  &x, &y, &z, absdetjk);
        break;
  
        case HexaTrilinear:
         ((THexaTrilinear *)rt)->SetCell(cell);
         ((THexaTrilinear *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound, zeta+DOFFound,
                                                  &x, &y, &z, absdetjk);
        break;
  
        case HexaIsoparametric:
         ((THexaIsoparametric *)rt)->SetCell(cell);
         ((THexaIsoparametric *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound, zeta+DOFFound,
                                                  &x, &y, &z, absdetjk);
        break;
      } // endswitch RefTrans
  
      break;
    } // endif DOFFound > -1
  } // endfor i
} // end GetDOFPosition
