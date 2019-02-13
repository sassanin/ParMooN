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
// Class:      TRefTrans3D
//
// Purpose:    reference transformations for 3D geometric objects
//
// Author:     Gunar Matthies
//
// History:    01.02.00 start implementation
// 
// =======================================================================

#ifndef __REFTRANS3D__
#define __REFTRANS3D__

#include <Constants.h>
#include <Enumerations.h>
#include <BaseCell.h>

/** reference transformations for 3D geometric objects */
class TRefTrans3D
{
  protected:
    TBaseCell *Cell;

  public:
    /** constuctor */
    TRefTrans3D();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double xi, double eta, double zeta,
                        double &x, double &y, double &z);

    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z,
                        double &xi, double &eta, double &zeta);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(TBaseCell *cell);

    /** set original element to cell */
    virtual void SetCell(TBaseCell *cell)
    {  Cell = cell; }

    /** return outer normal unit vector */
    void GetOuterNormal(int j, double s, double t,
                        double &n1, double &n2, double &n3);

    /** return two tangent vectors */
    void GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23);
    
    virtual void PiolaMapOrigFromRef(int N_Functs, double *refD00, 
                                     double *origD00)
    { ErrMsg(" Piola Map not defined for this element\n"); };

};

#endif
