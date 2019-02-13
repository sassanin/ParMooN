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
// @(#)Tetrahedron.h        1.2 10/18/99
//
// Class:       TTetrahedron
// Purpose:     shape descriptor of a tetrahedron
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __TETRAHEDRON__
#define __TETRAHEDRON__

#include <ShapeDesc.h>

#define TETRAMAXN_EpV   3
#define TETRAMAXN_VpF   3
#define TETRAMAXN_FpV   3
#define TETRAMAXN_EpF   3
#define TETRAMAXN_FpE   2

/** shape descriptor of a tetrahedron */
class TTetrahedron : public TShapeDesc
{
  public:
    // Constructor
    /** build the shape descriptor for a tetrahedron */
    TTetrahedron();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts);

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts);

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts);

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts);
};

#endif
