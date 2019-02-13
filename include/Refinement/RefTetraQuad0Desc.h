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
   
#ifndef __REFTETRAQUAD0DESC__
#define __REFTETRAQUAD0DESC__

#include <RefDesc.h>

#define REFTETRAQUAD0MAXN_VpC   4
#define REFTETRAQUAD0MAXN_CpV   4
#define REFTETRAQUAD0MAXN_EpC   6
#define REFTETRAQUAD0MAXN_CpE   3
#define REFTETRAQUAD0MAXN_FpC   4
#define REFTETRAQUAD0MAXN_CpF   2
#define REFTETRAQUAD0MAXN_EpV   6
#define REFTETRAQUAD0MAXN_VpF   3
#define REFTETRAQUAD0MAXN_FpV   9
#define REFTETRAQUAD0MAXN_EpF   3
#define REFTETRAQUAD0MAXN_FpE   4
#define REFTETRAQUAD0MAXN_nVpoF 6
#define REFTETRAQUAD0MAXN_oVpoF 3
#define REFTETRAQUAD0MAXN_iVpE  1
#define REFTETRAQUAD0MAXN_iEpF  3
#define REFTETRAQUAD0MAXN_nEpoE 2
#define REFTETRAQUAD0MAXN_nVpoE 3
#define REFTETRAQUAD0MAXN_nEpoF 6
#define REFTETRAQUAD0MAXN_nFpoF 4

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraQuad0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraQuad0Desc(TShapeDesc *shape);

    // Methods
};

#endif
