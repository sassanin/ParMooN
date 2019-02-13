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
// @(#)IsoInterfaceJoint.C        1.2 09/13/99
// 
// Class:       TIsoInterfaceJoint
// Purpose:     connects two cells on an interface with additional
//              vertices for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//
// =======================================================================

#include <IsoInterfaceJoint.h>
#include <BoundComp2D.h>
#include <BaseCell.h>

// Constructors
TIsoInterfaceJoint::TIsoInterfaceJoint(TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0) 
 : TInterfaceJoint(bdcomp, t_0, t_1, neighb0) 
{
  ID = IsoInterfaceJoint;

  N_Vertices = 0;
  Vertices = NULL;
}

TIsoInterfaceJoint::TIsoInterfaceJoint(TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0, TBaseCell *neighb1)
 : TInterfaceJoint(bdcomp, t_0, t_1, neighb0, neighb1) 
{
  ID = IsoInterfaceJoint;

  N_Vertices = 0;
  Vertices = NULL;
}

// Methods
TJoint *TIsoInterfaceJoint::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TIsoInterfaceJoint(BoundComp, T_0 + newT_0*(T_1 - T_0),
                             T_0 + newT_1*(T_1 - T_0), Me);
}

TJoint *TIsoInterfaceJoint::NewInst()
{
  return new TIsoInterfaceJoint(BoundComp, T_0, T_1, NULL);
}

void TIsoInterfaceJoint::SetVertices(int n_vertices, TVertex **vertices)
{
  delete Vertices;

  N_Vertices = n_vertices;
  Vertices = vertices;
}

void TIsoInterfaceJoint::GenerateVertices(int n_vertices)
{
  int i;
  TVertex *vertex;
  double t, x, y;
  
  if(N_Vertices != n_vertices)
  {
    if(Vertices)
      delete Vertices;

    N_Vertices = n_vertices;
    Vertices = new TVertex*[N_Vertices];
  }

#ifdef __2D__
  for(i=0;i<N_Vertices;i++)
  {
    t = T_0 + ((i+1)*(T_1-T_0))/(N_Vertices+1);
    
    BoundComp->GetXYofT(t, x, y);

    vertex = new TVertex(x, y);

    Vertices[i] = vertex;
  } // endfor i
#endif // __2D__
}

void TIsoInterfaceJoint::GeneratemidVert(int n_vertices, double*X, double*Y)
{
  int i;
  TVertex *vertex;
  double  x, y;

   if(Vertices)
     delete Vertices;

     N_Vertices = n_vertices;
     Vertices = new TVertex*[N_Vertices];

#ifdef __2D__
  for(i=0;i<N_Vertices;i++)
  {
   x = (X[0]+X[1])/2.0;
   y = (Y[0]+Y[1])/2.0;

   vertex = new TVertex(x, y);
   Vertices[i] = vertex;
//     cout <<N_Vertices << "test ("<< x <<", "<< y << ")"<<endl;
  } // endfor i
#endif // __2D__
}
