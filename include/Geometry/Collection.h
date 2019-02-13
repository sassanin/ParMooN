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
* @class TCollection 
* @date  14.10.97
* @brief store cells in an array
* @author Gunar Matthies & Sashikumaar Ganesan
* @History: MPI methods (Sashikumaar Ganesan, 08.08.14)
   
****************************************************************************/

#ifndef __COLLECTION__
#define __COLLECTION__

#include <BaseCell.h>
#include <JointCollection.h>

/** @brief store cells in an array, used by cell iterators */
class TCollection
{
  protected:
    /** @brief number of cells stored */
    int N_Cells;

    /** @brief array containing the pointers to the cells */
    TBaseCell **Cells;

    /** @brief array with all cells sorted by pointer */
    TBaseCell **SortedCells;

    /** @brief array with index of SortedCells in Cells */
    int *Index;

    #ifdef  _MPI
    /** @brief Number of own cells (excluding Halo cells) */
    int N_OwnCells;

    /** @brief array for Globalcell number in Cells */
    int *GlobalIndex;
    #endif

  public:
    /** @brief constructor */
    TCollection(int n_cells, TBaseCell **cells);

    /** @brief return number of cells */
    int GetN_Cells() const
    { return N_Cells; }

    /** @brief return Cell with index i in Cells-array */
    TBaseCell *GetCell(int i) const
    { return Cells[i]; }

    /** @brief destructor: delete arrays */
    ~TCollection();

    /** @brief get maximal and minimal diameter */
    int GetHminHmax(double *hmin, double *hmax);

    /** @brief return Index of cell in Cells-array */
    int GetIndex(TBaseCell *cell);

    /** @brief mark the vertices that are on the boundary */
    int MarkBoundaryVertices();

    /** @brief return Index of joints in Cells-array */
    TJointCollection  *GetJointCollection();

    /** @brief Generate   Vertex Neibs for all cells in the collection */
    void GenerateCellVertNeibs();

    /** @brief return the Index of the vertex in the sorted array */
    int GetIndex(TVertex **Array, int Length, TVertex *Element);
    
#ifdef  _MPI
    void SetN_OwnCells(int n_OwnCells)
     { N_OwnCells = n_OwnCells; }

    int GetN_OwnCells()
     { return N_OwnCells; }

    int GetN_HaloCells()
     { return (N_Cells - N_OwnCells); }

    int *GetGlobalIndex()
     {
      return GlobalIndex;
     }
#endif

   void Replace_Coll(int n_cells, TBaseCell **cells)
     {
      N_Cells = n_cells;
      Cells = cells;
     }

  private:
    /** @brief provide additional arrays */
    void GenerateSortedArrays();

    /** @brief return Index of cell in SortedCells-array */
    int GetSortedIndex(TBaseCell *cell);

};

#endif
