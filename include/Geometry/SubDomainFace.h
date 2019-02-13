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
// @(#)Joint.h        
// 
// Class:       TSubDomainFace
// Purpose:     superclass for faces between two subdomains
//
// Author:      Sashikumaar Ganesan  07.09.09
//
// History:     
//
// =======================================================================

#ifdef __3D__

#ifndef __SUBDOMAINFACE__
#define __SUBDOMAINFACE__

#include <stdlib.h>
#include <JointEqN.h>

/** for setting the neib cell information while refinement */
struct SubDomainMap
{
    int ParentLocalCell_No;
    int NeibParentLocalCell_No;
    int ParentLocalEdge_No;
    int NeibParentLocalEdge_No;
    int NewEdge;
    int NeibNewEdge;
    bool MapFilled;
};

/** connects two cells, which are in two/(more in 3D) different sub domains */
class TSubDomainFace : public TJointEqN
{

#endif
#endif
