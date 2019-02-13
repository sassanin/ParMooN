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
// @(#)HNDesc.h        1.1 10/30/98
//
// Class:       THNDesc
// Purpose:     describe fixed information for a hanging node
//
// Author:      Gunar Matthies  18.11.97
//
// =======================================================================

#ifndef __HNDESC__
#define __HNDESC__

/** describe fixed information for a hanging node */
class THNDesc
{
  protected:
    /** number of degrees in coupling */
    int N_Nodes;

    /** coefficient of other nodes in this coupling */
    double *Coeff;

  public:
    /** constructor, filling all data */
    THNDesc(int n_nodes, double *coeff);

    // Methods
    /** return number of nodes in coupling */
    int GetN_Nodes() 
    { return N_Nodes; }

    /** return coefficients of degrees of freedom in coupling */
    double *GetCoeff()
    { return Coeff; }

};

#endif
