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
// @(#)Domain.C        1.23 05/05/00
// 
// Class:       TTetGenMeshLoader
// Purpose:     creates domain with tetgen from a .smesh file
//
// Author:      Andreas Hahn 26.04.2010
//
// History:
//
// =======================================================================

#ifndef __TETGENMESHLOADER__
#define __TETGENMESHLOADER__

#define __TETGEN_14X__ 

#include <tetgen.h>

#include <Constants.h>
#include <MooNMD_Io.h>

class TBoundPart;
class TBoundComp3D;
class TVertex;
class TBaseCell;
class TJoint;
class TDomain;

class TTetGenMeshLoader
{
  protected:
    /** Filename of .smesh file without postfix **/
    char *mFileName;
    
    /** filename of .bound file with postfix **/
    char *mBoundFile;
    
    /** list of all boundary components **/
    TBoundComp3D **mAllBoundComps;
    
    /** list of all vertices **/
    TVertex **mAllVertices;
    
    /** list of all joints **/
    TJoint **mAllJoints;
    
    /** */
    int **mTrifaceHash;
    
    /** tetgen io objects **/
    tetgenio mTetIn;
    tetgenio mTetOut;
#ifdef __TETGEN_14X__
    tetgenbehavior mTetBeh;
    tetgenio mTetAddIn;
#endif
    
    /** options */
    bool plc;
    bool reconstruct;
    bool insertpoints;
    
  protected:
    int ReadBoundFile(int N_BoundComps);
    int ReadBoundComp(std::ifstream *in, int &CurrBoundComp, int &Range);
    int ReadBoundType(std::ifstream *in, char *buff);
    int ReadBoundParams(std::ifstream *in, char *type, int CurrBoundComp, int Range);
    int ReadBoundParamsCylinder(std::ifstream *in, double &radius,
				double &px, double &py, double &pz,
				double &ax, double &ay, double &az,
				double &nx, double &ny, double &nz);
    int ReadBoundParamsSphere(std::ifstream *in, double &radius, 
			      double &mx, double &my, double &mz);
				
    int ExtractValue(const char *buff, int increment, double &val);
    int ExtractVector(const char *buff, int increment,
		      double &val1, double &val2, double &val3);
		      
    int SetBdPlaneParams(int N_BoundComps);
    
    void Normalize3D(double *vec);
    void Cross3D(double *vec1, double *vec2, double *res);
    
    int AllocVertices(double &StartX, double &StartY, double &StartZ,
		      double &BoundX, double &BoundY, double &BoundZ);
    int AllocRootCells(TBaseCell **&CellTree, int &N_RootCells);
    int AllocJoints(TBaseCell **CellTree);
    int DistributeJoints(TBaseCell **CellTree);
    
    int  FindTriface(int a, int b, int c);
    void HashTrifaces();
    int  CreateAdjacency();
    
    void MakeOptionString(std::string &opts);
    
    void MakeBoundaryLayer_smesh();
    void MakeBoundaryLayer_msh();
    void FindBoundaryPoints_msh(int *pointlist, double *normallist);
    
    int Load_msh(const char *filename);    
    
    /** build boundary from .smesh file **/
    int BuildBoundary  (TBoundPart **&BdParts, int &N_BoundParts,
		        int &N_BoundComps, int *&StartBdCompID, int *&Interfaces);
		      
    /** build boundary from trifaces */
    int BuildBoundary2 (TBoundPart **&BdParts, int &N_BoundParts,
		        int &N_BoundComps, int *&StartBdCompID, int *&Interfaces);
		  
    /** triangulation with tegen **/
    int Tetgen();
    
    /** build rootcells */
    int BuildMooNMDMesh(TBaseCell **&CellTree, int &N_RootCells,
			double &StartX, double &StartY, double &StartZ,
			double &BoundX, double &BoundY, double &BoundZ);
  
  public:
    typedef struct
    {
      TBoundPart **BdParts;
      int N_BoundParts;
      int N_BoundComps;
      int *StartBdCompID;
      int *Interfaces;
      TBaseCell **CellTree;
      int N_RootCells;
      double StartX;
      double StartY;
      double StartZ;
      double BoundX;
      double BoundY;
      double BoundZ;
    } TDummyDomain;
			
  public:
    // CTOR
//     TTetGenMeshLoader() : mFileName(0), mAllBoundComps(0), mBoundFile(0),
// 			  mAllVertices(0), mAllJoints(0) {};
    TTetGenMeshLoader();
    TTetGenMeshLoader(const char *FileName);
    
    // DTOR
    ~TTetGenMeshLoader();
    
    /** generate mesh */
    int Generate(TDomain *Domain);
    int Generate(int N_Points, double *Points, int N_Facets, int *Facets,
		 int N_Regions, double *Regions, TDummyDomain *Domain);
};
#endif
