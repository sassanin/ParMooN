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
// @(#)InterfaceJoint.C        1.2 01/26/12
// 
// Class:       TInnerInterfaceJoint
// Purpose:     connects two cells of different collection
//              (not treated as a boundary joint)
//
//
// =======================================================================


#include <InnerInterfaceJoint.h>
#include <BoundComp2D.h>
#include <BaseCell.h>
#include <stdlib.h>

// Constructors

TInnerInterfaceJoint::TInnerInterfaceJoint(
    TBaseCell *neighb0, TBaseCell *neighb1) :
  TJointEqN(neighb0, neighb1)
{
  ID = InnerInterfaceJoint;
  NeibGlobalCellNo[0]=neighb0->GetCellIndex();
  NeibGlobalCellNo[1]=neighb1->GetCellIndex();
  SubDomainsID[0]=neighb0->GetReference_ID();
  SubDomainsID[1]=neighb1->GetReference_ID();
  children[0] = NULL;
  children[1] = NULL;
  
  // find coordinates of edge end points 
  // it does not matter which neighbor we take here 
  // first find which edge of neighb0 is this joint
  
  Xstart = -4711;
  Ystart = -4711;
  delX = -4711;
  delY = -4711;
  
}

TInnerInterfaceJoint::TInnerInterfaceJoint(TBaseCell *neighb0): 
  TJointEqN(neighb0)
{
  ID = InnerInterfaceJoint;
  NeibGlobalCellNo[0]=neighb0->GetCellIndex();
  NeibGlobalCellNo[1]=-1;
  SubDomainsID[0]=neighb0->GetReference_ID();
  SubDomainsID[1]=-1;
  children[0] = NULL;
  children[1] = NULL;
  
  Xstart = -4711;
  Ystart = -4711;
  delX = -4711;
  delY = -4711;
}


// Methods
TJoint *TInnerInterfaceJoint::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TInnerInterfaceJoint(Me);
}

TJoint *TInnerInterfaceJoint::NewInst()
{
  return new TInnerInterfaceJoint(NULL);
}

/** Remember which joints are children of this joint */
void TInnerInterfaceJoint::SetChild(TInnerInterfaceJoint *child)
{
  if(children[0] == NULL)
    children[0] = child;
  else if(children[1] == NULL)
    children[1] = child;
  // if two children are already set, then nothing happens,
  if(children[0]!=child && children[1]!=child)
  {
    ErrMsg("ERROR, TInnerInterfaceJoint::SetChild");
  }
}

/** return one of the two children of this edge */
TInnerInterfaceJoint *TInnerInterfaceJoint::GetChild(int child) const
{
  if(child == 0 || child == 1)
    return children[child];
  else
  {
    ErrMsg("ERROR, TInnerInterfaceJoint::GetChild, no such child available");
    exit(-4711);
  }
}

/** Get either one of the two neighbors */
TBaseCell *TInnerInterfaceJoint::GetNeighbour(int i) const
{
  if (i==1)
    return Neighb1;
  else if(i==0)
    return Neighb0;
  else
  {
    ErrMsg("TInnerInterfaceJoint::GetNeighbour: only two neighbors " <<i);
    exit(0);
  }
}

/** set the coordinates of the start point and the vector pointing to the     
    second point */
void TInnerInterfaceJoint::SetParams(double xstart, double ystart, 
                                     double delx, double dely)
{
  Xstart = xstart;
  Ystart = ystart;
  delX = delx;
  delY = dely;
}
/** get the coordinates of the start point and the vector pointing to the     
    second point */
void TInnerInterfaceJoint::GetParams(double &xstart, double &ystart, 
                                     double &delx, double &dely) const
{
  xstart = Xstart;
  ystart = Ystart;
  delx = delX;
  dely = delY;
}

double TInnerInterfaceJoint::GetLength() const
{
  return sqrt(delX*delX + delY*delY);
}

/** the unit normal of this edge */
void TInnerInterfaceJoint::GetNormal(double &nx, double &ny) const
{
  double length = sqrt(delX*delX + delY*delY);
  nx =  delY/length;
  ny = -delX/length;
}
  
/** the unit tangential of this edge */
void TInnerInterfaceJoint::GetTangent(double &tx, double &ty) const
{
  double length = sqrt(delX*delX + delY*delY);
  tx = delX/length;
  ty = delY/length;
}

/** get the index of this joint in given neighbor */
int TInnerInterfaceJoint::GetIndexInNeighbor(TBaseCell const * const neigh)
{
  if(neigh == Neighb0)
    return IndexInNeighbor[0];
  else if(neigh == Neighb1)
    return IndexInNeighbor[1];
  else
  {
    ErrMsg("ERROR, TInnerInterfaceJoint::GetIndexInNeighbor !!!!!!!!");
    exit(-4711);
  }
}

/** set the index of this joint in given neighbor */
void TInnerInterfaceJoint::SetIndexInNeighbor(TBaseCell *neigh, int index)
{
  if(neigh == Neighb0)
    IndexInNeighbor[0] = index;
  else if(neigh == Neighb1)
    IndexInNeighbor[1] = index;
  else
  {
    ErrMsg("ERROR, TInnerInterfaceJoint::SetIndexInNeighbor !!!!!!!!");
    exit(-4711);
  }
}
