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
   
#ifndef __REFTETRABIS3DESC__
#define __REFTETRABIS3DESC__

#include <RefDesc.h>

#define REFTETRABIS3MAXN_VpC   4
#define REFTETRABIS3MAXN_CpV   2
#define REFTETRABIS3MAXN_EpC   6
#define REFTETRABIS3MAXN_CpE   2
#define REFTETRABIS3MAXN_FpC   4
#define REFTETRABIS3MAXN_CpF   2
#define REFTETRABIS3MAXN_EpV   4
#define REFTETRABIS3MAXN_VpF   3
#define REFTETRABIS3MAXN_FpV   5
#define REFTETRABIS3MAXN_EpF   3
#define REFTETRABIS3MAXN_FpE   3
#define REFTETRABIS3MAXN_nVpoF 4
#define REFTETRABIS3MAXN_oVpoF 3
#define REFTETRABIS3MAXN_iVpE  1
#define REFTETRABIS3MAXN_iEpF  1
#define REFTETRABIS3MAXN_nEpoE 2
#define REFTETRABIS3MAXN_nVpoE 3
#define REFTETRABIS3MAXN_nEpoF 4
#define REFTETRABIS3MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis3Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis3Desc(TShapeDesc *shape);

    // Methods
};

#endif
