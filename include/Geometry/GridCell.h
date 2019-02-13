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
* @class     TGridCell
* @brief     represent geometric information of the cell
* @author    Volker Behns  09.07.97
* @History 
 ************************************************************************  */

#ifndef __GRIDCELL__
#define __GRIDCELL__

#include <BaseCell.h>

/**  @brief represent geometric information of the cell */
class TGridCell : public TBaseCell
{
  protected:
    /**  @brief field of pointer to children */
    TBaseCell **Children;
    /**  @brief pointer to parent cell */
    TBaseCell *Parent;

    /**  @brief field of all vertices */
    TVertex  **Vertices;

    /**  @brief grid level on with this cell was generated */
    int RefLevel;

  public:
    // Constructor
    TGridCell(TRefDesc *refdesc, int reflevel);

    // Destructor
    ~TGridCell();

    // Methods
    /**  @brief set the pointer to vertex with number i */
    virtual int SetVertex(int Vert_i, TVertex *Vert);
    /**  @brief return the pointer to vertex with number i */
    virtual TVertex *GetVertex(int Vert_i);
    /**  @brief return field of pointers to all vertices */
    TVertex **GetVertices()
    { return Vertices; }

    /**  @brief return number of children */
    virtual int GetN_Children();
    /**  @brief return number of parents */
    virtual int GetN_Parents();

    /**  @brief return pointer to child cell with number C\_i */
    virtual TBaseCell *GetChild(int C_i);
    /**  @brief return pointer to parent cell */
    virtual TBaseCell *GetParent();
    /**  @brief set parent */
    virtual int SetParent(TBaseCell *parent);
    /**  @brief return local number of child Me */
    virtual int GetChildNumber(TBaseCell *Me);

    /**  @brief put boundary and interface joints to postscript file */
    virtual int Draw(std::ofstream &dat, double scale, double StartX,
                   double StartY);
    /**  @brief put out postscript data to a file */
    virtual int PS(std::ofstream &dat, double scale, double StartX,
                   double StartY);
    /**  @brief put out raw data for MD-Format */
    virtual int MD_raw(std::ofstream &dat);

    /**  @brief derefine the cell */
    virtual int Derefine();
    /**  @brief refine or derefine the cell according to cell's clipboard */
    virtual int RefDeref();
    /**  @brief set marks in neighbour cells in order to maintain 1-regularity */
    virtual int Gen1RegMarks();
    /**  @brief generate conforming closures */
    virtual int MakeConfClosure();
    /**  @brief refine a cell */
    virtual int Refine(int RefLevel);

    #ifdef __MORTAR__
      /**  @brief refine a mortar cell */
      virtual int RefineMortar(int RefLevel);
    #endif

    /**  @brief generate a 2-regular grid */
    virtual int Gen1RegGrid();
    /**  @brief set refinement for the neighbour of your parent on joint LocJointNum */
    virtual int Ref1Reg(int LocJointNum, TBaseCell *&RefCell);
    /**  @brief check whether the surroundings of cell is 1-regular */
    virtual int Check1Reg();
    /**  @brief set RefDesc to no refinement */
    virtual int SetNoRefinement();
    /**  @brief set RefDesc to regular refinement */
    virtual int SetRegRefine();
    /**  @brief set RefDesc to adaptive refinement */
    virtual int Set1Refine(int i);      
    /**  @brief check whether a cell should be refined */
    virtual int IsToRefine();
    /**  @brief check whether exist some children */
    virtual int ExistChildren()
    { return Children == NULL ? FALSE : TRUE; }

#ifdef __2D__
    /**  @brief return coordinates of mid point P\_j on edge J\_i */
    virtual int LineMidXY(int J_i, int P_j, double &X, double &Y);
    /**  @brief return parameters on boundary of subedge SJ\_j on edge J\_i */
    virtual int LineMidT(int J_i, int SJ_j, double &T_0, double &T_1);
#else
    /**  @brief return whether a point is inside a cell */
    virtual bool PointInCell(double X, double Y, double Z);
#endif

    /**  @brief return whether a point is inside a cell */
    virtual bool PointInCell(double X, double Y);

    /**  @brief get diameter of a cell */
    virtual double GetDiameter()
    { return RefDesc->GetShapeDesc()->GetDiameter(Vertices); }

    /**  @brief get shortest edge of a cell */
    virtual double GetShortestEdge()
    { return RefDesc->GetShapeDesc()->GetShortestEdge(Vertices); }

    /**  @brief return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap()
    { return RefDesc->GetShapeDesc()->GetLengthWithReferenceMap(Vertices); }

     /**  @brief get measure of a cell */
    virtual double GetMeasure()
    { return RefDesc->GetShapeDesc()->GetMeasure(Vertices); }

    /**  @brief get geometry level */
    virtual int GetGeoLevel();

    /**  @brief return subgrid ID */
    virtual int GetSubGridID();

    /**  @brief compute number of edges at the boundary */
    virtual int GetN_BoundaryEdges();
#ifdef __3D__
    /**  @brief compute number of faces at the boundary */
    virtual int GetN_BoundaryFaces();    
#endif 
    /**  @brief compute number of vertices at the boundary */
    virtual int GetN_BoundaryVertices();
};

#endif
