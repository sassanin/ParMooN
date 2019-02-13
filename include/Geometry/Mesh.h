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
* @class     TMesh
* @brief     stores the information of a mesh 
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History 
 ************************************************************************  */

#ifndef __MESH__
#define __MESH__

#include <BaseCell.h>
#include <Vertex.h>
#include <Joint.h>

/**  @brief mesh details */
class TMesh : 
{
  protected:
    /**  @brief number of vertices in the mesh */    
    int N_RootVertices;
    
    /**  @brief number of joints (2D:edge, 3D:Face) in the mesh */    
    int N_Joints;    
    
    /**  @brief number of cells in the mesh */    
    int N_Cells;   

     /**  @brief cell-vetrtices index */    
    int *CellVertices;
    
    /**  @brief cell-joints (2D:edge, 3D:Face) in the mesh */    
    int *CellJoints;      
    
    /**  @brief array of pointers to vertices in the mesh  */
    TVertex  **Vertices;
 
    /**  @brief array of pointers to joints in the mesh  */
    TJoint **Joints;
    
    /**  @brief array of pointers to cells in the mesh */
    TBaseCell **CellTree;
 
    /**  @brief grid level on with this cell was generated */
    int RefLevel;

  public:
    // Constructor
    TMesh();

    TMesh(int N_RootVertices, int N_Joints, int N_Cells, int *CellVertices, int *CellJoints, TVertex  **Vertices, TJoint **Joints, TBaseCell **CellTree);
	
    // Destructor
    ~TMesh();

    // Methods
 

#ifdef __2D__
 
#else
 
 
#endif 
 
};

#endif
