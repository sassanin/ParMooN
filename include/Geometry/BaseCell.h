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
* @class      TBaseCell
* @brief      prototype of a cell
*

* @author     Volker Behns (97) & Gunar Matthies & Sashikumaar Ganesan
* @date       09.07.97
* @History    Making some methods inline (Gunar Matthies,  10.10.97)
              Adding ClipBoard (Gunar Matthies, 10.10.97)
              Twophase flows (Sashikumaar Ganesan. 20.09.09)
              Added edge info in 3D (Sashikumaar Ganesan, 04.09.10)
              MPI methods (Sashikumaar Ganesan, 08.08.14)
************************************************************************  */

#ifndef __BASECELL__
#define __BASECELL__

#include <Edge.h>
#include <Joint.h>
#include <RefDesc.h>
#include <fstream>

 /**  @brief information for finite element data structure */
class TBaseCell
{
  protected:
    /**  @brief current property of refinement (including shape descriptor) */
    TRefDesc *RefDesc;

    /**  @brief array of all joints */
    TJoint **Joints;

    /**  @brief an integer for storing clipboard information*/
    int ClipBoard;

    /** @brief an integer for storing physical reference of cell (e.g. material properties) **/
    int Reference_ID;
    
    /**  @brief an integer for storing boundary part (surface meshes) */
    int Bd_Part;

    /**  @brief cell index value in the collection */
    int CellIndex;

  /**  @brief an integer for storing the global cell number **/
    int GlobalCellNo;
    
    /** @brief normal orientation of joints (for each joint this is 1 or -1) */
    int *normalOrientation;
    
  /**  @brief an integer for storing the region of this cell number **/
    int region;

  /**  @brief an integer for indicating layer cells **/
    int LayerCell;
    
#ifdef __3D__
    /**  @brief array of all Edges in 3D */
    TEdge **Edges;
#endif

#ifdef  _MPI
  /**  @brief an integer for storing clipboard information in parallel FEspace mapping*/ 
  int  ClipBoard_Par; 
    
  /**  @brief an integer for storing which subdomain contains this cell **/
    int SubDomainNumber;

  /**  @brief an integer for storing the global cell number **/
    int SubDomainLocalCellNo;

  /**  @brief a bool to check this cell is own cell of the SubDomain or not **/
    bool OwnCell;

  /**  @brief a bool to check this cell is a halo cell for the SubDomain or not **/
    bool HaloCell;

  /**  @brief a bool to check this cell is neibs' halo cell or not **/
    bool DependentCell;
    
  /**  @brief a bool to check this cell contains SubDomainInterface(s) or not **/
    bool SubDomainInterfaceCell;

  /**  @brief a bool to check this cell contains  cross edges or not **/
    bool CrossEdgeCell;
    
  /**  @brief a bool to check this cell contains cross vertices or not **/
    bool CrossVertexCell;

  /**  @brief Number of neib processes for this cell **/
    int N_NeibProcesses;

  /**  @brief if(N_NeibProcesses), the rank ID of neib processes **/
    int *NeibProcessesIds;
 
#endif

  /**  @brief an integer for storing the local cell number (2Phase flows) **/
    int LocalCellNo;   
    

  public:
    // Constructor
    TBaseCell(TRefDesc *refdesc);

    // Destructor
    virtual ~TBaseCell();

    // Methods
    /**  @brief set refinement descriptor to newrefdesc */
    int SetRefDesc(TRefDesc *newrefdesc)
    {
      if (newrefdesc->GetShapeDesc() == GetShapeDesc())
      {
        RefDesc = newrefdesc;
        return 0;
      }
      else
        if (newrefdesc->GetShapeDesc()->GetN_Vertices() == 
              GetShapeDesc()->GetN_Vertices())
        {
          RefDesc = newrefdesc;
          return 0;
        }
        else
          return -1;
    }

    /**  @brief return refinement descriptor */
    TRefDesc *GetRefDesc() const
    { return RefDesc; }
    /**  @brief return shape descriptor of refinement descriptor */
    TShapeDesc *GetShapeDesc() const
    { return RefDesc->GetShapeDesc(); }
    /**  @brief return shape type of refinement descriptor */
    Shapes GetType() const
    { return RefDesc->GetShapeDesc()->GetType(); }

    /**  @brief return refinement descriptor of edge i */
    Refinements GetEdgeRef(int i) const
    { return RefDesc->GetEdgeRef(i); }

#ifdef __3D__
    /**  @brief return refinement descriptor of face i */
    Refinements GetFaceRef(int i)
    { return RefDesc->GetFaceRef(i); }

    /**  @brief set the pointer to edge E_i to E */
    int SetEdge(int E_i, TEdge *E)
    {
      Edges[E_i] = E;
      return 0;
    }

    /**  @brief return the pointer to edge with number E_i */
    TEdge *GetEdge(int E_i)
    { return Edges[E_i]; }

#endif

    /**  @brief set the pointer of vertex Vert\_i to Vert */
    virtual int SetVertex(int Vert_i, TVertex *Vert) = 0;
    /**  @brief return the pointer to vertex with number i */
    virtual TVertex *GetVertex(int Vert_i) = 0;

    /**  @brief set the pointer to face J\_i to J */
    int SetJoint(int J_i, TJoint *J)
    {
      Joints[J_i] = J;
      return 0;
    }

    /**  @brief return the pointer to face with number i */
    TJoint *GetJoint(int J_i)
    { return Joints[J_i]; }

    /**  @brief return the number of vertices of the cell */
    int GetN_Vertices()
    { return RefDesc->GetN_OrigVertices(); }
    /**  @brief return the number of edges of the cell */
    int GetN_Edges()
    {  return RefDesc->GetN_OrigEdges(); }
    /**  @brief return the number of joints */
    int GetN_Joints()
    {  return RefDesc->GetShapeDesc()->GetN_Joints(); }

    #ifdef __3D__
      /**  @brief return the number of faces of the cell */
      int GetN_Faces()
      { return RefDesc->GetN_OrigFaces(); }
    #endif

    /**  @brief return the number of children of the cell */
    virtual int GetN_Children() = 0;
    /**  @brief return the number of parents of the cell */
    virtual int GetN_Parents() = 0;

    /**  @brief return the child with the number C\_i */
    virtual TBaseCell *GetChild(int C_i) =  0;
    /**  @brief return the parent cell */
    virtual TBaseCell *GetParent() =  0;
    /**  @brief set the parent to parent */
    virtual int SetParent(TBaseCell *parent) =  0;
    /**  @brief return the child number of cell Me */
    virtual int GetChildNumber(TBaseCell *Me) =  0;

    /**  @brief write boundary and interface joints to stream dat */
    virtual int Draw(std::ofstream &dat, double scale, double StartX,
                   double StartY) = 0;
    /**  @brief write the postscript cell data to stream dat */
    virtual int PS(std::ofstream &dat, double scale, double StartX,
                   double StartY) = 0;
    /**  @brief write cell data according to MD format in stream dat */
    virtual int MD_raw(std::ofstream &dat) = 0;

    /**  @brief refine the current cell on level RefLevel according actual
        refinement descriptor */
    virtual int Refine(int RefLevel) = 0;

    #ifdef __MORTAR__
      /**  @brief refine a mortar cell */
      virtual int RefineMortar(int RefLevel) = 0;
    #endif

    /**  @brief derefine the current cell, remove the children */
    virtual int Derefine() = 0;

    /**  @brief make refinement or derefinement according to cell's clipboard */
    virtual int RefDeref() = 0;
    /**  @brief set marks in neighbour cell in order to maintain 1 regularity */
    virtual int Gen1RegMarks() = 0;

    /**  @brief generate such refinement information in the neighbouring cells
        so that the traingulation will become 1-regular */
    virtual int Gen1RegGrid() = 0;
    /**  @brief regular refinement of joint LocJointNum */
    virtual int Ref1Reg(int LocJointNum, TBaseCell *&RefCell) = 0;
    /**  @brief check whether the surroundings of cell is 1-regular */
    virtual int Check1Reg() = 0;
    /**  @brief set refinement descriptor to no refinement */
    virtual int SetNoRefinement() = 0;
    /**  @brief set refinement descriptor to regular refinement */
    virtual int SetRegRefine() = 0;
    /**  @brief set refinement descriptor to adaptive refinement */
    virtual int Set1Refine(int i)= 0;    
    /**  @brief is the cell to refine */
    virtual int IsToRefine() = 0;
    /**  @brief are there any children of this cell */
    virtual int ExistChildren() = 0;
    /**  @brief generate conforming closures */
    virtual int MakeConfClosure() = 0;

    #ifdef __2D__
      /**  @brief return (x,y) coordinates */
      virtual int LineMidXY(int J_i, int P_j, double &X, double &Y) = 0;
      /**  @brief return parameter values */
      virtual int LineMidT(int J_i, int SJ_j, double &T_0, double &T_1) = 0;
    #else
    #endif

    /**  @brief set value in ClipBoard */
    void SetClipBoard(int value)
    { ClipBoard=value; }
    /**  @brief get value from ClipBoard */
    int GetClipBoard()
    { return ClipBoard; }

    /**  @brief get diameter of a cell */
    virtual double GetDiameter() = 0;

    /**  @brief return shortest edge of a cell */
    virtual double GetShortestEdge() = 0;

    /**  @brief return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap() = 0;

     /**  @brief get measure of a cell */
    virtual double GetMeasure() = 0;
    
    /** @brief get the value of hK 
     * 
     * This function calls either this->GetDiameter(), thid->GetShortestEdge, or
     * this->GetMeasure(), depending on cell_measure.
     * 
     * Typically you should set cell_measure = TDatabase::ParamDB->CELL_MEASURE.
     */
    double Get_hK(int cell_measure);

    /**  @brief return whether a point is inside a cell */
    virtual bool PointInCell(double X, double Y) = 0;

#ifdef __3D__
    virtual bool PointInCell(double X, double Y, double Z) = 0;

     // added 25.04.2010 for fixing refinement problem
     void CorrectBoundaryVertices(TVertex **NewVertices, TJoint **NewJoints);
#endif

    /**  @brief get geometry level */
    virtual int GetGeoLevel() = 0;

    /**  @brief get subgrid ID */
    virtual int GetSubGridID() = 0;

    /** @brief set reference number to this cell   */
    void SetReference_ID(int val)
    { Reference_ID = val;}
    
    /** @brief get reference number of this cell   */
    int GetReference_ID() const
    {return Reference_ID;}
    
    /**  @brief set phase number to this cell   */
    void SetBd_Part(int val)
       {Bd_Part = val;}

    /**  @brief get phase number to this cell   */
    int GetBd_Part() const
       {return Bd_Part;}

    /**  @brief set subdomain number to this cell   */
    void SetCellIndex(int val)
       {CellIndex = val;}

    /**  @brief set subdomain number to this cell   */
    int GetCellIndex() const
       {return CellIndex;}

    /**  @brief set subdomain number to this cell   */
    void SetGlobalCellNo(int val)
       {GlobalCellNo = val;}

    /**  @brief set subdomain number to this cell   */
    int GetGlobalCellNo() const
       {return GlobalCellNo;}

    /**  @brief set subdomain number to this cell   */
    void SetLocalCellNo(int val)
       {LocalCellNo = val;}

    /**  @brief set LocalCellNo number to this cell   */
    int GetLocalCellNo() const 
       {return LocalCellNo;}

    /**  @brief set region number to this cell   */
    void SetRegionID(int val)
       {region = val;}

    /**  @brief set region number to this cell   */
    int GetRegionID() const
       {return region;}       

     /**  @brief set as LayerCell cell   */
    void SetAsLayerCell(int val)
       {LayerCell = val;}

    /**  @brief get LayerCell info  */
    int IsLayerCell() const
       {return LayerCell;}          

    /** @brief compute normal orientation w.r.t cell */
    void SetNormalOrientation();
    
    /** @brief get normal orientation w.r.t cell at i-th joint */
    int GetNormalOrientation(int i) const
    { return normalOrientation[i]; }

       
#ifdef  _MPI

    /**  @brief set value in ClipBoard */
    void SetClipBoard_Par(int value)
    { ClipBoard_Par=value; }
    
    /**  @brief get value from ClipBoard */
    int GetClipBoard_Par()
    { return ClipBoard_Par; }
    
    /**  @brief set subdomain number to this cell   */
    void SetSubDomainNo(int val)
       {SubDomainNumber = val;}

    /**  @brief set subdomain number to this cell   */
    int GetSubDomainNo() const
       {return SubDomainNumber;}

    void SetAsOwnCell()
     { OwnCell=TRUE; }

    void SetAsDependentCell()
     { DependentCell=TRUE; }

    void SetAsHaloCell()
     { 
      OwnCell=FALSE;
      HaloCell=TRUE;
     }

    bool IsHaloCell() const
     {  return HaloCell; }

    void SetAsSubDomainInterfaceCell()
     { SubDomainInterfaceCell=TRUE; }

    bool IsSubDomainInterfaceCell() const
     {  return SubDomainInterfaceCell; }

    void SetAsCrossEdgeCell()
     { CrossEdgeCell=TRUE; }

    bool IsCrossEdgeCell() const
     {  return CrossEdgeCell; }


    void SetAsCrossVertexCell()
     { CrossVertexCell=TRUE; }

    bool IsCrossVertexCell() const
     {  return CrossVertexCell; }

    bool IsDependentCell() const
     {  return DependentCell; }

    void SetN_NeibProcesses(int n)
     { N_NeibProcesses = n; }

    int GetN_NeibProcesses() const
     { return N_NeibProcesses; }

    void SetNeibProcessesIds(int *Neiblist);

    int *GetNeibProcessesIds() const
     { return NeibProcessesIds; } 
     
#endif
};

#endif
