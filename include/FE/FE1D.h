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
// @(#)FE1D.h        1.1 10/30/98
//
// Class:       TFE1D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  08.10.98
//
// =======================================================================

#ifndef __FE1D__
#define __FE1D__

#include <AllClasses.h>
#include <Constants.h>
#include <Enumerations.h>
#include <BaseFunct1D.h>
#include <FEDesc1D.h>

/** store all information for one finite element class */
class TFE1D
{
  protected:
   /** ID for set of basis functions */
   BaseFunct1D BaseFunct_ID;

   /** set of basis function */
   TBaseFunct1D *BaseFunct;

   /** ID for set of nodal functional */
   NodalFunctional1D NodalFunctional_ID;

   /** set of nodal functional */
   TNodalFunctional1D *NodalFunctional;

   /** ID for reference transformation */
   RefTrans1D RefTransID;

   /** ID for element description */
   FEDesc1D FEDesc_ID;

   /** element description */
   TFEDesc1D *FEDesc;

   /** number of needed integer entries (numbers + infos) */
   int Size;

   /** number of degrees of freedom */
   int N_DOF;

   /** number of info blocks */
   int N_Info;

  public:
    /** constructor */
    TFE1D();

    /** constructor with data */
    TFE1D(BaseFunct1D basefunct_id, NodalFunctional1D nodalfunctional_id,
        RefTrans1D reftransid, FEDesc1D fedesc_id,
        int info);

    /** return BaseFunct1D_ID */
    BaseFunct1D GetBaseFunct1D_ID()
      { return BaseFunct_ID; };

    /** return BaseFunct1D */
    TBaseFunct1D *GetBaseFunct1D()
      { return BaseFunct; };

    /** return BaseFunct1D_ID and BaseFunct1D */
    void GetBaseFunct1D(BaseFunct1D &ID,
                        TBaseFunct1D* &Obj)
      { ID = BaseFunct_ID; Obj = BaseFunct; };

    /** return NodalFunctional1D_ID */
    NodalFunctional1D GetNodalFunctional1D_ID()
      { return NodalFunctional_ID; };

    /** return NodalFunctional1D */
    TNodalFunctional1D *GetNodalFunctional1D()
      { return NodalFunctional; };

    /** return NodalFunctional1D_ID and NodalFunctional1D */
    void GetNodalFunctional1D(NodalFunctional1D &ID,
                              TNodalFunctional1D* &Obj)
      { ID = NodalFunctional_ID; Obj = NodalFunctional; };

    /** return RefTransID */
    RefTrans1D GetRefTransID()
      { return RefTransID; };

    /** return FEDesc1D_ID */
    FEDesc1D GetFEDesc1D_ID()
      { return FEDesc_ID; };

    /** return FEDesc1D */
    TFEDesc1D *GetFEDesc1D()
      { return FEDesc; };

    /** return FEDesc1D_ID and FEDesc1D */
    void GetFEDesc1D(FEDesc1D &ID, TFEDesc1D* &Obj)
      { ID = FEDesc_ID; Obj = FEDesc; };

    /** return size */
    int GetSize()
      { return Size; };

    /** return number of degrees of freedom */
    int GetN_DOF()
      { return N_DOF; };

    /** return number of info blocks */
    int GetN_Info()
      { return N_Info; };
};

#endif
