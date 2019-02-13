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
// @(#)ShapeDesc.h        1.5 11/15/99
//
// Class:       TShapeDesc
// Purpose:     super class of all shape descriptors
//
// Author:      Volker Behns  09.07.97
//
// =======================================================================

#ifndef __SHAPEDESC__
#define __SHAPEDESC__

#include <Constants.h>
#include <Vertex.h>

#define N_SHAPES  8
enum Shapes {S_Line, Triangle, Quadrangle, Parallelogram, Rectangle,
             Tetrahedron, Hexahedron, Brick};

#ifndef __3D__
  #define MAXN_JOINTS   4
#else
  #define MAXN_JOINTS   6
  #define MAXN_EDGES   12 
#endif

/** super class of all shape descriptors */
class TShapeDesc
{
  protected:
    /** type of shape */
    Shapes Type;

    /** number of vertices */
    int N_Vertices;
    /** number of edges */
    int N_Edges;
    /** maximum number of edges per vertex */
    int MaxN_EpV;

    #ifdef __3D__
      /** number of faces (3D) */
      int N_Faces;
      /** maximum number of vertices per face */
      int MaxN_VpF;
      /** maximum number of faces per vertex */
      int MaxN_FpV;
      /** maximum number of edges per face */
      int MaxN_EpF;
      /** maximum number of faces per edge */
      int MaxN_FpE;
    #endif

    /** number of joints */
    int N_Joints;

    /** which vertices belong to one edge */
    const int *EdgeVertex;
    /** which edges meet at a vertex */
    const int *VertexEdge;

    #ifdef __3D__
      /** which vertices are on one face */
      const int *FaceVertex;
      /** number of  vertices on one face */
      const int *FaceVertexLen;
      /** which edges are on one face */
      const int *FaceEdge;
      /** number of edges on one face */
      const int *FaceEdgeLen;
      /** which shapes have the faces got */
      const Shapes *FaceType;
      /** which faces meet at a vertex */
      const int *VertexFace;
      /** which faces meet a one edge */
      const int *EdgeFace;
    #endif

  public:
    // Constructor

    // Methods
    /** return the number of vertices */
    int GetN_Vertices()
    { return N_Vertices; }
    /** return the number of edges */
    int GetN_Edges()
    { return N_Edges; }
    /** return the number of joints */
    int GetN_Joints()
    { return N_Joints; }

    #ifdef __3D__
      /** return the number of faces */
      int GetN_Faces()
      { return N_Faces; }
    #endif

    /** return the shape */
    Shapes GetType()
    { return Type; }

    /** return the EdgeVertex array */
    int GetEdgeVertex(const int *&TmpEV)
    {
      TmpEV = EdgeVertex;
      return 0;
    }

    /** return the MaxN_EpV */
    int GetMaxN_EpV()
    { return MaxN_EpV; }

    /** return the EdgeVertex array */
    int GetVertexEdge(const int *&TmpVE)
    {
      TmpVE = VertexEdge;
      return 0;
    }

    #ifdef __3D__
      /** return the FaceVertex array */
      int GetFaceVertex(const int *&TmpFV, const int *&TmpLen, int &MaxLen)
      {
        TmpFV = FaceVertex;
        TmpLen = FaceVertexLen;
        MaxLen = MaxN_VpF;
        return 0;
      }
      /** return the FaceEdge array */
      int GetFaceEdge(const int *&TmpFV, const int *&TmpLen, int &MaxLen)
      {
        TmpFV = FaceEdge;
        TmpLen = FaceEdgeLen;
        MaxLen = MaxN_EpF;
        return 0;
      }

      /** return the EdgeFace array */
      int GetEdgeFace(const int *&TmpEF, int &MaxLen)
      {
        TmpEF = EdgeFace;
        MaxLen = MaxN_FpE;
        return 0;
      }

      /** return the FaceType array */
      int GetFaceType(const Shapes *&TmpFT)
      {
        TmpFT = FaceType;
        return 0;
      }
    #endif
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts) = 0;

    /** return shortest of a cell */
    virtual double GetShortestEdge(TVertex **Verts) = 0;

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts) = 0;

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts) = 0;
};

#endif
