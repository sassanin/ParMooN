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
// @(#)IsoBoundFace.C        1.2 09/13/99
// 
// Class:       TIsoBoundFace
// Purpose:     face on a boundary component with additional vertices
//              for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
// 
// =======================================================================

#include <BoundComp2D.h>
#include <IsoBoundFace.h>
#include <Vertex.h>

#include <MooNMD_Io.h>

#include <stdlib.h>

// Constructors
TIsoBoundFace::TIsoBoundFace(TBoundComp3D *bdcomp, double *param1,
                             double *param2)
 : TBoundFace(bdcomp, param1, param2)
{
  ID = IsoBoundFace;

  Vertices = NULL;
  N_Vertices = 0;
}

TIsoBoundFace::TIsoBoundFace(TBoundComp3D *bdcomp)
 : TBoundFace(bdcomp)
{
  ID = IsoBoundFace;

  Vertices = NULL;
  N_Vertices = 0;
}

// Methods
int TIsoBoundFace::CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

// create a new instance of this class
TJoint *TIsoBoundFace::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TIsoBoundFace(BoundComp);
}

TJoint *TIsoBoundFace::NewInst()
{
  return new TIsoBoundFace(BoundComp);
}

void TIsoBoundFace::SetVertices(int n_vertices, TVertex **vertices)
{
  if(Vertices)
    delete Vertices;

  N_Vertices = n_vertices;
  Vertices = vertices;

}

void TIsoBoundFace::GenVert(int N_NewVert,const int N_V, double **LinComb,
                            double *X, double *Y, double *Z )
{
    cerr << __FILE__ << ":" << __LINE__ << ": GentVert() is currently not working properly" << endl;
    exit(0);
  
//   int i, n_vertices;
//   TVertex *vertex;
//   double x, y, z, r, T, S;
// 
//     if(Vertices)
//       delete Vertices;
// 
//     N_Vertices = N_NewVert; 
//     Vertices = new TVertex*[N_Vertices];
//     RefVertices = new TVertex*[N_Vertices];
//     IsoParam1 = new double[N_Vertices];
//     IsoParam2 = new double[N_Vertices];
// 
//     for(i=0; i<N_NewVert; i++)
//      {
//       BoundComp->GetXYZ(N_V, LinComb[i],
//                         X, Y, Z,
//                         x, y, z);
// 
//       vertex = new TVertex(x, y, z);
//       Vertices[i] = vertex;
// 
//       vertex = new TVertex(x, y, z);
//       RefVertices[i] = vertex;
// 
// // cout<< " X: " << x<<" Y: " << y<<" Z: " << z<<endl;
// // needed for evolving ellipsoid problem
//       r = sqrt(x*x + y*y + z*z);
//       IsoParam1[i] = atan2(y, x);
//       IsoParam2[i] = acos(z/r);
//      }

}

void TIsoBoundFace::GenHexaVert(int N_NewVert,const int N_V, double **LinComb,
                            double *X, double *Y, double *Z )
{
  cerr << __FILE__ << ":" << __LINE__ << ": GentHexaVert() is currently not working properly" << endl;
  exit(0);
  
//   int i, n_vertices;
//   TVertex *vertex;
//   double x, y, z, T, S;
// 
// 
// //implemented for second order for other order one have to change
//     if(Vertices)
//       delete Vertices;
// 
//     N_Vertices = N_NewVert;
//     Vertices = new TVertex*[N_Vertices];
// 
//     for(i=0; i<N_Vertices; i++)
//      {
//       BoundComp->GetXYZ(N_V, LinComb[i],
//                         X, Y, Z,
//                         x, y, z);
// // cout<< " X: " << x<<" Y: " << y<<" Z: " << z<<endl;
//       vertex = new TVertex(x, y, z);
//       Vertices[i] = vertex;
//      }
}

void TIsoBoundFace::GenerateVertices(int n_vertices)
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

  Error("TIsoBoundFace GenerateVertices cannot be called directly!" << endl);
}



void TIsoBoundFace::GetParameters(double *OutIsoParam1, double *OutIsoParam2)
{
  int i;
// cout<< " X: " <<N_Vertices << endl;
  for(i=0;i<N_Vertices;i++)
  {
    OutIsoParam1[i] = IsoParam1[i];
    OutIsoParam2[i] = IsoParam2[i];

// cout<< " X: " << OutIsoParam1[i] <<" Y: " <<  OutIsoParam2[i] <<endl;
  }
}
