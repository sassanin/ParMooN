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
// %W% %G%
//
// Class:       TFE3D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  19.11.99
//
// =======================================================================

#ifndef __FE3D__
#define __FE3D__

#include <AllClasses.h>
#include <Constants.h>
#include <Enumerations.h>
#include <BaseFunct3D.h>
#include <FEDesc3D.h>
#include <NodalFunctional3D.h>

/** store all information for one finite element class */
class TFE3D
{
  protected:
   /** ID for set of basis functions */
   BaseFunct3D BaseFunct_ID;

   /** set of basis function */
   TBaseFunct3D *BaseFunct;

   /** ID for set of nodal functional */
   NodalFunctional3D NodalFunctional_ID;

   /** set of nodal functional */
   TNodalFunctional3D *NodalFunctional;

   /** ID for reference transformation */
   RefTrans3D RefTransID;

   /** ID for element description */
   FEDesc3D FEDesc_ID;

   /** element description */
   TFEDesc3D *FEDesc;

   /** number of needed integer entries (numbers + infos) */
   int Size;

   /** number of degrees of freedom */
   int N_DOF;

   /** number of info blocks */
   int N_Info;

  public:
    /** constructor */
    TFE3D();

    /** constructor with data */
    TFE3D(BaseFunct3D basefunct_id, NodalFunctional3D nodalfunctional_id,
        RefTrans3D reftransid, FEDesc3D fedesc_id,
        int info);

    /** return BaseFunct3D_ID */
    BaseFunct3D GetBaseFunct3D_ID()
      { return BaseFunct_ID; };

    /** return BaseFunct3D */
    TBaseFunct3D *GetBaseFunct3D()
      { return BaseFunct; };

    /** return BaseFunct3D_ID and BaseFunct3D */
    void GetBaseFunct3D(BaseFunct3D &ID,
                        TBaseFunct3D* &Obj)
      { ID = BaseFunct_ID; Obj = BaseFunct; };

    /** return NodalFunctional3D_ID */
    NodalFunctional3D GetNodalFunctional3D_ID()
      { return NodalFunctional_ID; };

    /** return NodalFunctional3D */
    TNodalFunctional3D *GetNodalFunctional3D()
      { return NodalFunctional; };

    /** return NodalFunctional3D_ID and NodalFunctional3D */
    void GetNodalFunctional3D(NodalFunctional3D &ID,
                              TNodalFunctional3D* &Obj)
      { ID = NodalFunctional_ID; Obj = NodalFunctional; };

    /** return RefTransID */
    RefTrans3D GetRefTransID()
      { return RefTransID; };

    /** return FEDesc3D_ID */
    FEDesc3D GetFEDesc3D_ID()
      { return FEDesc_ID; };

    /** return FEDesc3D */
    TFEDesc3D *GetFEDesc3D()
      { return FEDesc; };

    /** return FEDesc3D_ID and FEDesc3D */
    void GetFEDesc3D(FEDesc3D &ID, TFEDesc3D* &Obj)
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

    /** check N[i](b[j]) = delta[ij] */
    void CheckNFandBF();
};

#endif
