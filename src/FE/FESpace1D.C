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
// @(#)FESpace1D.C        1.9 06/13/00
// 
// Class:       TFESpace1D
// Purpose:     class for all 1D finite element spaces
//
// Author:      Gunar Matthies (04.11.97)
//
// History:     start of implementation 09.10.98 (Gunar Matthies)
//        :     Reimplemented 04.06.2011 (Sashikumaar Ganesan)
// =======================================================================

#include <Constants.h>
#include <FESpace1D.h>
#include <Joint.h>
#include <Database.h>
#include <FEDatabase2D.h>

#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>

#include "Vector.h"

/** Constructor */
TFESpace1D::TFESpace1D(TCollection *coll, char *name, char *description) :
     TFESpace(coll, name, description)
{
  UsedElements = NULL;
  AllElements = NULL;
  ElementForShape = NULL;
  DGSpace=0;
}

// =====================================================================
// rules: if a neighbour element is not in the collection, its clipboard
//        must be set to -1
// =====================================================================

/** constructor for building a space with elements of order k */
TFESpace1D::TFESpace1D(TCollection *coll, char *name, char *description, 
                       int ord) : TFESpace(coll, name, description)
{
  UsedElements = NULL;
  AllElements = NULL;

  ElementForShape = new FE1D[N_SHAPES];

  // build ElementForShape array
  switch(ord)
  {
    case -1: ElementForShape[S_Line] = N_P0_1D_L_A;
//              ElementForShape[1] = N_P0_1D_L_A;
             break;
    case 1:  ElementForShape[S_Line] = C_P1_1D_L_A;
//              ElementForShape[1] = C_P0_1D_L_A;
             break;
    case 2:  ElementForShape[S_Line] = C_P2_1D_L_A;
//              ElementForShape[1] = C_P1_1D_L_A;
             break;
    case 3:  ElementForShape[S_Line] = C_P3_1D_L_A;
//              ElementForShape[1] = C_P2_1D_L_A;
             break;
    // for NS
    case -2: ElementForShape[S_Line] = C_P2_1D_L_A;
//              ElementForShape[1] = C_P1_1D_L_A;
             break;
    case -3: ElementForShape[S_Line] = C_P3_1D_L_A;
//              ElementForShape[1] = C_P2_1D_L_A;
             break;
	     
              // discontinous elements
    case -11: ElementForShape[S_Line] = D_P1_1D_L_A;
//               ElementForShape[1] = N_P0_1D_L_A;
              break;
    case -12: ElementForShape[S_Line] = D_P2_1D_L_A;
//               ElementForShape[1] = D_P1_1D_L_A;
              break;	     

    default: cerr << "unknown order" << endl;
             exit(-1);
             break;
  } // endswitch

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace();
  
  DGSpace=0;
}

/** constructor for building a space with the given elements */
TFESpace1D::TFESpace1D(TCollection *coll, char *name, char *description,
               FE1D *fes) : TFESpace(coll, name, description)
{
  UsedElements = NULL;
  ElementForShape = NULL;

  AllElements = fes;

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace();
  
  DGSpace=0;
}

/** return the FE Id for element i, corresponding to cell */
FE1D TFESpace1D::GetFE1D(int i, TBaseCell *cell)
{
  FE1D ret;

  if(AllElements)
    ret=AllElements[i];
  else
  {
    // cout << "cell: " << (int)cell << endl;
    // cout << "j0: " << (int)(cell->GetJoint(0)) << endl;
    // cout << "n0: " << (int)(cell->GetJoint(0)->GetNeighbour(cell)) << endl;
    // cout << "j1: " << (int)(cell->GetJoint(1)) << endl;
    // cout << "n1: " << (int)(cell->GetJoint(1)->GetNeighbour(cell)) << endl;
//     if(cell->GetJoint(0)->GetNeighbour(cell) &&
//        cell->GetJoint(1)->GetNeighbour(cell))
//           ret=ElementForShape[0];
//      else
//           ret=ElementForShape[1];
    ret=ElementForShape[cell->GetType()];
  }

  return ret;
}

void TFESpace1D::FindUsedElements()
{
  TBaseCell *cell;
  int i, j, N_;
  int Used[N_FEs1D];

  memset(Used,0,N_FEs1D*SizeOfInt);

  N_ = N_Cells;
  for(i=0;i<N_;i++)
  {
    cell=Collection->GetCell(i);
    Used[GetFE1D(i, cell)] = 1;
  }

  for(i=0;i<N_FEs1D;i++)
    if(Used[i])
      N_UsedElements++;

  UsedElements = new FE1D[N_UsedElements];
  j=0;
  for(i=0;i<N_FEs1D;i++)
    if(Used[i])
    {
      UsedElements[j]=(FE1D)i;
      j++;
    }

  if(TDatabase::ParamDB->SC_VERBOSE>3)
   {
    cout << "N_UsedElements: " << N_UsedElements << endl;
    for(i=0;i<N_UsedElements;i++)
     cout << "UsedElement[" << i << "]: " << UsedElements[i] << endl;
   }
}

void TFESpace1D::ConstructSpace()
{
  int i, j, k, l, m, m2, n, comp, N_Joints, NEdges;
  int *v;
  TBaseCell *cell, *neigh, *child1, *child2;
  TJoint *joint;
  double t0,t1;

  TFE2DMapper *mapper;
  TFE2DMapper1Reg *mapper1reg;

  const int *TmpoEnE, *TmpLen1, *TmpEC, *TmpLen2, *TmpoEnlE;
  int MaxLen1, MaxLen2;
  TRefDesc *refdesc;

  int SumLocDOF, N_JointDOF;
  int count, *BoundaryUpperBound;
  int *BoundCounter, Counter;
  int *BoundMark, DirichletMark;
  int DirichletOffset, InnerOffset, SlaveOffset;
  int *BoundOffset;

  FE1D FEType0, FEType1, FEType2;
  TFE1D *FE0, *FE1, *FE2;
  TFEDesc1D *FEDesc0_Obj, *FEDesc1_Obj, *FEDesc2_Obj;
  FEDesc1D FEDesc0, FEDesc1, FEDesc2;

  int I_K0, I_K1, I_K2;
  int *J_K0, *J_K1, *J_K2;
  int *Indices0, *Indices1, *Indices2;
  int c1, c2, e1, e2, chnum1, chnum2;

  // boundary structures (no boundary values used)
  N_DiffBoundNodeTypes=N_BOUNDCOND-1; // not only Neumann nodes possible
  BoundaryNodeTypes=new BoundCond[N_DiffBoundNodeTypes];
  N_BoundaryNodes=new int[N_DiffBoundNodeTypes];
  BoundaryNodesBound=new int[N_DiffBoundNodeTypes];

  BoundaryNodeTypes[0]=NEUMANN;
  BoundaryNodeTypes[1]=ROBIN;
  N_BoundaryNodes[0]=0;
  N_BoundaryNodes[1]=0;

  // reset clipboards to -1
  for(i=0;i<N_Cells;i++)
  {
    cell=Collection->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);
      if(neigh) neigh->SetClipBoard(-1);
    }
    cell->SetClipBoard(-1);
  } // endfor i

  // set number i into clipboard, count the number of local degrees of
  // freedom
  count=0;

  BeginIndex=new int[N_Cells+1];
  BeginIndex[0]=0;
  for(i=0;i<N_Cells;i++)
  {
    cell=Collection->GetCell(i);
    cell->SetClipBoard(i);

    FEType0 = GetFE1D(i, cell);
    FE0 = TFEDatabase2D::GetFE1D(FEType0);
    count += FE0->GetSize();

    BeginIndex[i+1]=count;
  } // endfor i

  // cout << "count: " << count << endl;
  SumLocDOF = count;
  GlobalNumbers = new int[SumLocDOF];
  memset(GlobalNumbers, -1, SizeOfInt*SumLocDOF);

  // start DOF manager
  l = FIRSTMARK + 1;
  Counter = l;

  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);

    N_Joints = cell->GetN_Vertices();

    FEType0 = GetFE1D(i, cell);
    FE0 = TFEDatabase2D::GetFE1D(FEType0);
    FE0->GetFEDesc1D(FEDesc0, FEDesc0_Obj);
    I_K0 = BeginIndex[i];
    J_K0 = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      Indices0 = FEDesc0_Obj->GetJointDOF(j);
      N_JointDOF = FEDesc0_Obj->GetN_JointDOF();
      
      // if (N_JointDOF)
      if(joint)
       {
        if(!(neigh = (joint->GetNeighbour(cell))))
        {
          // boundary joint
          for(k=0;k<N_JointDOF;k++)
          {
            if(J_K0[Indices0[k]]==-1)
            {
              // not handled yet
              Counter--;
              J_K0[Indices0[k]] = Counter;
            }
          }
        } // boundary joint
       else
        {
          // no boundary joint
          n = neigh->GetClipBoard();
          if (n>i && N_JointDOF!=0)
          {
            // this joint was not handled until now

            FEType1 = GetFE1D(n, neigh);
            FE1 = TFEDatabase2D::GetFE1D(FEType1);
            FE1->GetFEDesc1D(FEDesc1, FEDesc1_Obj);
            I_K1 = BeginIndex[n];
            J_K1 = GlobalNumbers + BeginIndex[n];
            Indices1 = FEDesc1_Obj->GetJointDOF(l);

            // find the local edge of neigh on which cell is
            l=0;
            // while(neigh->GetJoint(l)->GetNeighbour(neigh)!=cell) l++;
            while(neigh->GetJoint(l) != joint) l++;

            Indices1 = FEDesc1_Obj->GetJointDOF(l);

            // cout << "i= " << i << endl;
            // cout << "j= " << j << endl;
            // cout << "n= " << n << endl;
            // cout << "l= " << l << endl;

            // cout << "FEDesc0: " << FEDesc0 << endl;
            // cout << "FEDesc1: " << FEDesc1 << endl;

            k = J_K0[Indices0[0]];
            if(k<=FIRSTMARK)
            {
              // DOF already handled
              J_K1[Indices1[0]] = I_K0+Indices0[0];
            }
            else
            {
              // DOF not handle yet
              J_K1[Indices1[0]] = I_K0+Indices0[0];
              Counter--;
              J_K0[Indices0[0]] = Counter;
            }
          } // end neigh
        } // no boundary joint
       }
    } // endfor j

    // handle inner degrees of freedom
    k = FEDesc0_Obj->GetN_InnerDOF();
    Indices0 = FEDesc0_Obj->GetInnerDOF();
    for(j=0;j<k;j++)
    {
      Counter--;
      J_K0[Indices0[j]] = Counter;
    } // endfor j
  } // endfor i

  // find global numbers
  l=0;
  count = l;

  N_Inner = 0;

  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);

    FEType0 = GetFE1D(i, cell);
    FE0 = TFEDatabase2D::GetFE1D(FEType0);
    FE0->GetFEDesc1D(FEDesc0, FEDesc0_Obj);
    J_K0 = GlobalNumbers + BeginIndex[i];

    k=FEDesc0_Obj->GetN_DOF();
    for(j=0;j<k;j++)
    {
      l = J_K0[j];
      if (l < -1)
      {
        // node to assign here
        J_K0[j] = count;
        count++;
        N_Inner++;
      } // l < -1
      else
      {
        if (l >= 0)
        {
          J_K0[j] = GlobalNumbers[l];
        }
        else
        {
          cerr << "J_K0[j]==-1" << endl;
        }
      } // l >= -1
    } // endfor j
  } // endfor i

/*
  // print for all elements for global numbers of their local dofs
  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    cout << "cell number: " << i << endl;

    J_K0 = GlobalNumbers + BeginIndex[i];
    FEType0 = GetFE1D(i, cell);
    FE0 = TFEDatabase2D::GetFE1D(FEType0);
    k = FE0->GetSize();
    for(j=0;j<k;j++)
    {
      cout << j << ": " << " number: " << J_K0[j] << endl;
    }
    cout << endl;
  } // endfor i
*/

  InnerBound = N_Inner;
  l = N_Inner;
  DirichletBound = l;
  ActiveBound = N_Inner;
  N_DegreesOfFreedom = l;
}


TFESpace1D::~TFESpace1D()
{
  delete [] BoundaryNodesBound;

  if (UsedElements)
    delete [] UsedElements;

  if(AllElements)
    delete [] AllElements;

  if(ElementForShape)
    delete [] ElementForShape;
}
