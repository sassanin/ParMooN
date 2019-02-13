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
// @(#)IsoEdge3D
// 
// Class:       TIsoEdge3D
// Purpose:     class for iso edges in 3D
//
// Author:      Sashikumaar Ganesan  06.09.2010
//
// History:
//
// ======================================================================= 

#ifndef __ISOEDGE3D__
#define __ISOEDGE3D__

#include <Vertex.h>
#include <Edge.h>


class TIsoEdge3D : public TEdge
{
  protected:
    int N_Vertices;

    TVertex **Vertices;

    int *LocDof;

    int LocVert[2];

  public:
    /** constructor with neighbours */
    TIsoEdge3D(int n_Neibs, TBaseCell **neighbs);

    /** methods */
    void SetIsoVertInfo(int n_vertices, TVertex **vertices, int *locdof, int *locvert);

    ~TIsoEdge3D();

    int GetN_Vertices()
    { return N_Vertices; }

    TVertex *GetVertex(int i)
    { return Vertices[i]; }

    int GetLocalDof (int i)
    { return LocDof[i]; }

    void GetLocalVertices(int &a, int &b)
    { a = LocVert[0]; b = LocVert[1]; }
};

#endif
