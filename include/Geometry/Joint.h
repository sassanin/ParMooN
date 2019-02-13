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
// @(#)Joint.h        1.7 11/15/99
// 
// Class:       TJoint
// Purpose:     superclass for edges and faces
//
// Author:      Volker Behns  23.07.97
//
// History:     Add method for getting a certain neighbour
//              (Gunar Matthies 17.10.97)
//
// =======================================================================

#ifndef __JOINT__
#define __JOINT__

#include <AllClasses.h>
#include <Constants.h>
#include <Mapper.h>

#ifndef __3D__
  #define MAXN_nVpoJ  3
  #define MAXN_nJpoJ  2
#else
  #define MAXN_nVpoJ  9
  #define MAXN_nJpoJ  4
#endif

struct StoreGeom
{
  TVertex *Vertices[MAXN_nVpoJ];
  TJoint *Joints[MAXN_nJpoJ];
  bool Filled;
};

#ifdef __MORTAR__
  struct StoreGeomMortarStruct
  {
    TVertex **Vertices;
    TJoint **Joints;
    bool Filled;
  };

  typedef struct StoreGeomMortarStruct StoreGeomMortar;
#endif

/** supercall for edges and faces */
class TJoint
{
  protected:
    JointType ID;

    /** first neighbour */
    TBaseCell *Neighb0;
    /** second neighbour */
    TBaseCell *Neighb1;

    /** value in ClipBoard (dG method)*/
    int ClipBoard;
    
    /** */
    int NeibSubDomainLocalJointNo;
    
    
#ifdef __3D__
    int MapType;
#endif

  public:
    // Constructors
    TJoint();

    // Methods
    /** return type */
    JointType GetType()
    { return ID; }  
    
    void ChangeType(JointType New_ID)
     {ID = New_ID;}

#ifdef __3D__
    /** set mapper type automatically */
    void SetMapType();

    /** set mapper type */
    void SetMapType(int maptype)
    { MapType = maptype; }

    /** return mapper type */
    int GetMapType() const
    { return MapType; }
    
    /** Function is used to get local edge index on neighboured element */
    int GetNeighbourEdgeIndex(TBaseCell*, int);
#endif

    /** check the refinement pattern on both sides for matching,
        return already existing object on the joint in Tmp */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp) = 0;

    #ifdef __MORTAR__
      /** check the refinement pattern on both sides for matching,
          special version for moratr cells */
      virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                   StoreGeomMortar &Tmp) = 0;
    #endif

    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me) = 0;
    virtual TJoint *NewInst() = 0;

    /** set the neighbour to Neighb */
    int SetNeighbour(TBaseCell *Neighb);
    /** return the neighbour of this joint which is not equal to Me */
    TBaseCell *GetNeighbour(TBaseCell *Me) const;

    /** set neighbour i to Neighb */
    int SetNeighbour(int i, TBaseCell *Neighb);
    /** return neighbour with number i */
    TBaseCell *GetNeighbour(int i) const;

    /** remove a neighbour */
    void Delete(TBaseCell *Neighb);

    /** function for debug purpose only */
    TBaseCell *GetNeighb(int i) const
    { return(i ? Neighb1 : Neighb0); }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const = 0;

    #ifdef __3D__
      /** return mapper of refined vertices and faces */
      void GetMapperRef(const int *&MapVerts, const int *&MapFaces);

      /** return mapper of original vertices and edges */
      void GetMapperOrig(const int *&MapVerts, const int *&MapEdges);
    #endif
            
    /** set value in ClipBoard */
    void SetClipBoard(int value)
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard()
    { return ClipBoard; }
    
    void SetNeibSubDomainLocalJointNo(int value)
    { NeibSubDomainLocalJointNo = value; }
    
    int GetNeibSubDomainLocalJointNo()
    { return NeibSubDomainLocalJointNo; }  
    
    // Destructor
    virtual ~TJoint();

};

#endif
