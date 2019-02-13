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
// @(#)FE1D.C        1.1 10/30/98
//
// Class:       TFE1D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  08.10.98
//
// =======================================================================

#include <FE1D.h>
#include <FEDatabase2D.h>

/** constructor */
TFE1D::TFE1D()
{
}

/** constructor with data */
TFE1D::TFE1D(BaseFunct1D basefunct_id, NodalFunctional1D nodalfunctional_id,
         RefTrans1D reftransid, FEDesc1D fedesc_id, int n_info)
{
  BaseFunct_ID = basefunct_id;
  BaseFunct = TFEDatabase2D::GetBaseFunct1D(BaseFunct_ID);

  NodalFunctional_ID = nodalfunctional_id;
  NodalFunctional  = TFEDatabase2D::GetNodalFunctional1D(NodalFunctional_ID);

  RefTransID = reftransid;

  FEDesc_ID = fedesc_id;
  FEDesc = TFEDatabase2D::GetFEDesc1D(FEDesc_ID);

  N_Info = n_info;
  N_DOF = BaseFunct->GetDimension();

  Size = N_Info + N_DOF;
}
