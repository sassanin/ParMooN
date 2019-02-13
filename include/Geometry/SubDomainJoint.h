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
// @(#)SubDomainJoint.h        
// 
// Class:       TSubDomainJoint
// Purpose:     superclass for edges in two/more subdomains
//
// Author:      Sashikumaar Ganesan  07.09.09
//
// History:     
//
// =======================================================================

#ifndef __SUBDOMAINJOINT__
#define __SUBDOMAINJOINT__

#include <stdlib.h>
#include <JointEqN.h>

/** edges which are in interface between two different sub domains */
class TSubDomainJoint : public TJointEqN
{
  protected:
    /** number of subdomains contain this joint */
    int N_NeibSubDomains;

    /** subdomain number of the neib cell, which contains this joint */
    int NeibSubDomainRank;

    /** global cell number of the neibs' cell, which contains this joint */
    int NeibSubDomainGlobalCellNo;

    /** local joint number of this joint in the neib cell */
    int NeibSubDomainLocalJointNo;

  public:
    // Constructors
    TSubDomainJoint(TBaseCell *neighb0);

    /** constructor with one initial neighbour */
    TSubDomainJoint(TBaseCell *neighb0, TBaseCell *neighb1);

     /** constructor with neighbour info */
    TSubDomainJoint(TBaseCell *neighb0, TBaseCell *neighb1, int neibID, 
                    int neibSubDomainGlobalCellNo,  int neibSubDomainLocalJointNo);

    // Methods
    /** return whether this is an interior joint */
    virtual bool InnerJoint()
    { return false; }

    /** return the subdomain number of the neib cell joint */
    int GetNeibRank()
    { return NeibSubDomainRank; }

    /** return the subdomain number of the neib cell joint */
    int GetNeibGlobalCellNo()
    { return NeibSubDomainGlobalCellNo; }

    /** return the subdomain number of the neib cell joint */
    int GetNeibLocalJointNo()
    { return NeibSubDomainLocalJointNo; }


  };
#endif
