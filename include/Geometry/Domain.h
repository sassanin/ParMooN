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
   
/** ************************************************************************ 
*
* @class TDomain  
* @date  09.07.97
* @brief  contains the boundary description, the virtual cell tree and macro grid
* @author Volker Behns
* @History: collection methods (Gunar Matthies 14.10.97), 
            methods for Refine/Derefine algorithm  (Gunar Matthies 17.10.97)
            mesh class, edge generation, parallel methods and documentation (Sashikumaar Ganesan 08.08.2014)
   
************************************************************************  */

#ifndef __DOMAIN__
#define __DOMAIN__

#include <BaseCell.h>
#include <BoundPart.h>
#include <Collection.h>
#include <Iterator.h>
#include <Vertex.h>

class TDatabase;

#ifdef __MORTAR__
struct TMortarFaceStruct
       {
         TBaseCell *Cell;
         int LocFaceNumber[2];
       };

typedef struct TMortarFaceStruct TMortarFace;
#endif

/** contains the boundary description, the virtual cell
    tree and macro grid */
class TDomain
{
  protected:
    
    /** @brief number of boundary parts */
    int N_BoundParts;    
    
    /** @brief boundary parts of domain */
    TBoundPart **BdParts;

    /** @brief number of all boundary components */
    int N_BoundComps;
    
    /** @brief start id of boundary component on each boundary part */
    int *StartBdCompID;

    /** @brief boundary part id's of all Interfaces */
    int *Interfaces;

    /** @brief number of holes */
    int N_Holes;
    
    /** @brief point in each hole */
    double *PointInHole;

    /** @brief number of regions */
    int N_Regions;
    
    /** @brief point in each region */
    double *PointInRegion;

    /** @brief array of all root cells of cell tree */
    TBaseCell **CellTree;
    
    /** @brief number of all root cells of cell tree */
    int N_RootCells;

    /** @brief number of virtuell cells on initial level */
    int N_InitVCells;

    /** @brief x coordinate of the start point (2D) */
    double StartX;
    /** @brief y coordinate of the start point (2D) */
    double StartY;
    /** @brief x length of bounding box */
    double BoundX;
    /** @brief y length of bounding box */
    double BoundY;

#ifdef __3D__
      /** @brief third coordinate of start point (3D) */
      double StartZ;
      /** @brief return number of cell in Y direction */
      int N_InitVCellsY;
      /** @brief z length of the bounding box */
      double BoundZ;
#endif

    /** @brief current refinment level */
    int RefLevel;

#ifdef __MORTAR__
      /** @brief number of mortar faces */
      int N_MortarFaces;
      /** @brief structur for mortar faces */
      TMortarFace *MortarFaces;

      /** @brief begin of each mortar face on coll */
      int *BeginMFace;
#endif
    
    friend class TTetGenMeshLoader;

#ifdef  _MPI
      /** @brief array contains the global cell number of local cells (including Halo cells) */
      int *GlobalCellIndex;

      /** @brief Number of own cells (excluding Halo cells) */
      int N_OwnCells;
#endif

  public:
    // Constructors
    TDomain();

    TDomain(char *ParamFile);
    
    // Methods
    /** @brief read geometry file */
    int ReadGeo(char *GeoFile);

    /** @brief read Gmsh mesh */
    int GmshGen(char *GeoFile);
    
#ifdef __3D__
    /** @brief read sandwich geometry */
    int ReadSandwichGeo(char *GeoFile);

    /** @brief make boundary parameter consistent */
    void MakeBdParamsConsistent(TCollection *coll);
    
    int CloseGrid(int level);
    
    /** @brief read TetGen mesh */
    int TetrameshGen(char *GeoFile);
    
#endif

#ifdef __2D__
    /** @brief set boundary parameters and cell shape according to
        possibly moving vertices */
    void CorrectParametersAndShapes();
    
    /** @brief mesh genration using Triangle for give IN */
    void TriMeshGen(struct triangulateio *In);
#endif

    /** @brief read parameter file */
    int ReadParam(char *ParamFile);
    /** @brief read boundary parameterization */
    int ReadBdParam(char *ParamFile, int &Flag);
    /** @brief read mapping and mortar information */
    int ReadMapFile(char *MapFile, TDatabase *database);

    /** @brief get boundary part of BdCompID */
    int GetBdPartID(int BdCompID);
    /** @brief get local number of boundary component */
    int GetLocalBdCompID(int BdCompID);
    /** @brief get local number of last boundary component on part BdPartID*/
    int GetLastLocalComp(int BdPartID)
    { return StartBdCompID[BdPartID+1] - StartBdCompID[BdPartID] - 1; }
    /** @brief set start BdCompID on boundary part i */
    void SetStartBdCompID(int BdCompID, int i)
    { StartBdCompID[i] = BdCompID; }

    /** @brief get i-th boundary part */
    TBoundPart* GetBdPart(int i)
    { return BdParts[i]; }

    /** @brief get tree of cells */
    void GetTreeInfo(TBaseCell **&celltree, int &N_rootcells)
    { 
      celltree = CellTree;
      N_rootcells = N_RootCells;
    }

    /** @brief set tree of cells */
    void SetTreeInfo(TBaseCell **celltree, int N_rootcells)
    {
      CellTree = celltree;
      N_RootCells = N_rootcells;

      #ifdef  _MPI
      N_OwnCells = 0;
      #endif 
    }

    #ifdef __MORTAR__
      /** @brief set subgrid ID's on all MacroCells and generate mortar structurs */
      int SetSubGridIDs(IntFunct2D *TestFunc);
      /** @brief generate mortar structurs */
      int GenMortarStructs();
      /** @brief return number of mortar face structs */
      int GetN_MortarFace()
      { return N_MortarFaces; }
      /** @brief return mortar face struct */
      TMortarFace *GetMortarFace(int i)
      { return &(MortarFaces[i]); }

      /** @brief get a collection of all mortar cells */
      TCollection *GetMortarColl(Iterators it, int level);

      /** @brief initialize all mortar joints with FE-information */
      int InitMortarJoints(Iterators it, int level, TCollection *coll);
    #endif

    /** @brief generate initial grid using external mesh generator */
    int GenInitGrid();
    #ifdef __2D__
      /** @brief make initial 2D grid */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int N_Vertices,
                   int NVE);
      /** @brief make initial 2D grid, extended version which sets ReferenceID
       *         in cells */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                   int N_Vertices, int NVE);
    #else
      /** @brief make initial 3D grid */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                   int N_Vertices, int NVE, int *BoundFaces, int *FaceParam,
                   int NBF, int NVpF,
                   int *Interfaceparam, int N_Interfaces);
      /** @brief make initial sandwich grid */
      int MakeSandwichGrid(double *DCORVG, int *KVERT, int *KNPR,
                           int N_Vertices, int NVE,
                           double DriftX, double DriftY, double DriftZ,
                           int N_Layers, double *Lambda);
    #endif

    /** @brief Init process for current domain */
    void Init(char *PRM, char *GEO);

    /** @brief write domain boundary  into a postscript file */
    int Draw(char *name, Iterators iterator, int arg);
    /** @brief write mesh into a postscript file */
    int PS(const char *name, Iterators iterator, int arg);
    /** @brief write collection into a postscript file */
    int PS(const char *name, TCollection *Coll);
    /** @brief write files for MD-Out format */
    int MD_raw(const char *name, Iterators iterator, int arg);

    /** @brief refine the grid according the cell refinement descriptors */
    int Refine();
    /** @brief refine all cells regular */
    int RegRefineAll();
    /** @brief refine all cells in subgrid ID regular */
    int RegRefineSub(int ID);
    /** @brief refine only in one direction */
    int RefineallxDirection();
    /** @brief generate a 1-regular grid */
    int Gen1RegGrid();
    /** @brief refine the finest grid according the given indicator function */
    int RefineByIndicator(DoubleFunct2D *Indicator);
    /** @brief refine the finest grid if necessary in order to get a 
        grid with conforming closures */
    int MakeConfClosure();

    /** @brief refine the finest grid according a given error estimate */
    int RefineByErrorEstimator(TCollection *Collection,double *eta_K,
                               double eta_max,double tolerance,
                               bool ConfClosure);

    /** @brief refine/derefine algorithm for a 1-regular grid, geolevel of all
        cells on the finest grid is between MinLevel and MaxLevel */
    void Refine1Reg(int MinLevel, int MaxLevel);

    /** @brief derefinemnt */
    void DeRefine();

    /** @brief convert all finest quadrangles into two triangles */
    int ConvertQuadToTri(int type);

    /** @brief produce a collection with all cells returned by iterator it */
    TCollection *GetCollection(Iterators it, int level);
    
    /** @brief produce a collection with all cells of a given Collection which
     *         have the given reference as ReferenceID 
     * 
     * This will give you a subcollection.
     */
    TCollection *GetCollection(TCollection *coll, int reference);

#ifdef  _MPI 
    /** @brief produce a own collection with all cells returned by iterator it */
    TCollection *GetOwnCollection(Iterators it, int level, int ID);
#endif

    /** @brief produce a collection with all cells in the finest grid, sort 
        they according to their geometry level and return in Indices 
        the indices where each level starts */
    void GetSortedCollection(TCollection* &Coll, int* &Indices);

    /** @brief get bounding box parameters */
    void GetBoundBox(double &startx, double &starty,
                     double &boundx, double &boundy)
    {
      startx = StartX;
      starty = StartY;
      boundx = BoundX;
      boundy = BoundY;
    }
    
#ifdef __3D__
    /** @brief get bounding box parameters */
    void GetBoundBox(double &startx, double &starty, double &startz,
                     double &boundx, double &boundy, double &boundz)
    {
      startx = StartX;
      starty = StartY;
      startz = StartZ;
      boundx = BoundX;
      boundy = BoundY;
      boundz = BoundZ;
    }
    
    void SetBoundBox(double startx, double starty, double startz,
                     double boundx, double boundy, double boundz)
    {
      StartX = startx;
      StartY = starty;
      StartZ = startz;
      BoundX = boundx;
      BoundY = boundy;
      BoundZ = boundz;
    }
#endif
    
    // test
    #ifndef __3D__
      void TestGrid1();
      void TestGrid2();
      void TestGrid3();
      void TestMortar();
      void TestShishkin();
      void TriangleShishkin();
      void UnitSquare();
      void UnitSquareRef();
      void TwoTriangles();
      void TwoTrianglesRef();
      void SquareInSquare();
      void SquareInSquareRef();
      void SetBoundBox(double boundx, double boundy);
      void SetBoundBoxstart(double startx , double starty);
      void RefOnMortarEdge();
      void RefCardioide(double A);
      void PeriodicSquares();
      void PeriodicSquaresLarge();
      void PeriodicRectangle_2_4();
      void PeriodicTrianglesLarge();
      void QuadShishkin(double tau1, double tau2);
      void Rectangular(int dimx, int dimy);
      void TestTriaConf();
      void TestTriaConf2();
      void CheckCells();

    #else
      void TestGrid3D();
      void SetBoundBox(double boundx, double boundy, double boundz);
    #endif

    #ifdef __3D__
      int Grape(const char *name, TCollection *coll);

    #endif

    void TetrameshGen();


    #ifdef  _MPI
      void ReplaceTreeInfo(int n_cells, TBaseCell **cells, int *GLOB_cellIndex, int n_OwnCells)
       {
        if(CellTree) delete CellTree;
        N_RootCells = n_cells;
        CellTree = cells;
        GlobalCellIndex = GLOB_cellIndex;
        N_OwnCells = n_OwnCells;
       }

      void SetN_OwnCells(int n_OwnCells)
       { N_OwnCells = n_OwnCells; }

      int GetN_OwnCells()
       { return N_OwnCells; }

      int GetN_HaloCells()
       { return (N_RootCells - N_OwnCells); }
     #endif
    
#ifdef __3D__
  public: 

//      int Tetgen(const char*);

     /** @brief generate edge info in 3D mesh **/
     int GenerateEdgeInfo();

#endif

  /** @brief adaptive refine  */
  int AdaptRefineAll();    
     
};

#endif
