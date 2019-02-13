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
// @(#)FE2D.h        1.2 05/04/99
//
// Class:       TFE2D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  09.07.98
//
// =======================================================================

#ifndef __FE2D__
#define __FE2D__

#include <Enumerations.h>
#include <BaseFunct2D.h>
#include <NodalFunctional2D.h>
#include <FEDesc2D.h>

/** store all information for one finite element class */
class TFE2D
{
  protected:
   /** ID for set of basis functions */
   BaseFunct2D BaseFunct_ID;

   /** set of basis function */
   TBaseFunct2D *BaseFunct;

   /** ID for set of nodal functional */
   NodalFunctional2D NodalFunctional_ID;

   /** set of nodal functional */
   TNodalFunctional2D *NodalFunctional;

   /** ID for reference transformation */
   RefTrans2D RefTransID;

   /** ID for element description */
   FEDesc2D FEDesc_ID;

   /** element description */
   TFEDesc2D *FEDesc;

   /** number of needed integer entries (numbers + infos) */
   int Size;

   /** number of degrees of freedom */
   int N_DOF;

   /** number of info blocks */
   int N_Info;

  public:
    /** constructor */
    TFE2D();

    /** constructor with data */
    TFE2D(BaseFunct2D basefunct_id, NodalFunctional2D nodalfunctional_id,
        RefTrans2D reftransid, FEDesc2D fedesc_id,
        int info);

    /** return BaseFunct2D_ID */
    BaseFunct2D GetBaseFunct2D_ID() const
      { return BaseFunct_ID; };

    /** return BaseFunct2D */
    TBaseFunct2D *GetBaseFunct2D() const
      { return BaseFunct; };

    /** return BaseFunct2D_ID and BaseFunct2D */
    void GetBaseFunct2D(BaseFunct2D &ID,
                        TBaseFunct2D* &Obj) const
      { ID = BaseFunct_ID; Obj = BaseFunct; };

    /** return NodalFunctional2D_ID */
    NodalFunctional2D GetNodalFunctional2D_ID() const
      { return NodalFunctional_ID; };

    /** return NodalFunctional2D */
    TNodalFunctional2D *GetNodalFunctional2D() const
      { return NodalFunctional; };

    /** return NodalFunctional2D_ID and NodalFunctional2D */
    void GetNodalFunctional2D(NodalFunctional2D &ID,
                              TNodalFunctional2D* &Obj) const
      { ID = NodalFunctional_ID; Obj = NodalFunctional; };

    /** return RefTransID */
    RefTrans2D GetRefTransID() const
      { return RefTransID; };

    /** return FEDesc2D_ID */
    FEDesc2D GetFEDesc2D_ID() const
      { return FEDesc_ID; };

    /** return FEDesc2D */
    TFEDesc2D *GetFEDesc2D() const
      { return FEDesc; };

    /** return FEDesc2D_ID and FEDesc2D */
    void GetFEDesc2D(FEDesc2D &ID, TFEDesc2D* &Obj) const
      { ID = FEDesc_ID; Obj = FEDesc; };

    /** return size */
    int GetSize() const
      { return Size; };

    /** return number of degrees of freedom */
    int GetN_DOF() const
      { return N_DOF; };

    /** return number of info blocks */
    int GetN_Info() const
      { return N_Info; };

    /** check N[i](b[j]) = delta[ij] */
    void CheckNFandBF();
};

#endif
