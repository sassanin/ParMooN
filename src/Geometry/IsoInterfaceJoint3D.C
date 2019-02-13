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
// @(#)IsoInterfaceJoint3D.C        1.2 09/13/99
// 
// Class:       TIsoInterfaceJoint3D
// Purpose:     connects two cells on an interface with additional
//              vertices for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//              Gunar Matthies  05.04.02
//
// =======================================================================

#include <IsoInterfaceJoint3D.h>
#include <BoundComp2D.h>
#include <BaseCell.h>

// Constructors
/** initialize the joint with the boundary parameters and one neighbour */
TIsoInterfaceJoint3D::TIsoInterfaceJoint3D(TBoundComp3D *bdcomp,
                double *param1, double *param2, TBaseCell *neigh0)
   : TInterfaceJoint3D(bdcomp, param1, param2, neigh0)
{
  ID = IsoInterfaceJoint3D;

  N_Vertices = 0;
  Vertices = NULL;
}

/** initialize the joint with the boundary parameters and two neighbours */
TIsoInterfaceJoint3D::TIsoInterfaceJoint3D(TBoundComp3D *bdcomp,
                double *param1, double *param2,
                TBaseCell *neigh0, TBaseCell *neigh1)
   : TInterfaceJoint3D(bdcomp, param1, param2, neigh0, neigh1)
{
  ID = IsoInterfaceJoint3D;

  N_Vertices = 0;
  Vertices = NULL;
}

/** initialize the joint with the boundary parameters and two neighbours */
TIsoInterfaceJoint3D::TIsoInterfaceJoint3D(TBoundComp3D *bdcomp,
                TBaseCell *neighb0)
   : TInterfaceJoint3D(bdcomp, neighb0)
{
  ID = IsoInterfaceJoint3D;

  N_Vertices = 0;
  Vertices = NULL;
}


// Destructor
TIsoInterfaceJoint3D::~TIsoInterfaceJoint3D()
   {
    if(Neighb0)
     { Neighb0 = NULL;}
  
    if(Neighb1)
     { Neighb1 = NULL;}   
   }
   

// Methods
TJoint *TIsoInterfaceJoint3D::NewInst(double newtT_0, double newT_1, TBaseCell *Me)
{
  return new TIsoInterfaceJoint3D(BoundComp, Me);
}

TJoint *TIsoInterfaceJoint3D::NewInst()
{
  return new TIsoInterfaceJoint3D(BoundComp, NULL);
}

void TIsoInterfaceJoint3D::SetVertices(int n_vertices, TVertex **vertices)
{
  delete Vertices;

  N_Vertices = n_vertices;
  Vertices = vertices;
}

void TIsoInterfaceJoint3D::GenerateVertices(int n_vertices)
{
}
