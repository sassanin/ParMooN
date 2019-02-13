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
// @(#)Hexahedron.h        1.2 10/18/99
//
// Class:       THexahedron
// Purpose:     shape descriptor of a hexahedron
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __HEXAHEDRON__
#define __HEXAHEDRON__

#include <ShapeDesc.h>

#define HEXAMAXN_EpV   3
#define HEXAMAXN_VpF   4
#define HEXAMAXN_FpV   3
#define HEXAMAXN_EpF   4
#define HEXAMAXN_FpE   2

/** shape descriptor of a hexahedron */
class THexahedron : public TShapeDesc
{
  public:
    // Constructor
    THexahedron();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts);

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts);

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts);

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts);

    /** check whether this hexahedron is a brick */
    Shapes CheckHexa(TVertex **Vertices);
};

#endif
