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
// @(#)BoundComp.h        1.5 08/12/99
//
// Class:       TBoundComp
// Purpose:     components of boundary faces
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BOUNDCOMP__
#define __BOUNDCOMP__

#include <Constants.h>
#include <fstream>

#define TOL_SECANT_BOUND 1.01
#define MAXVERTS 20
#define RECURS_DEPTH 8

enum BoundTypes {Line = 1, Circle, Spline2D, Polygon, NonUniformSpline2D,
                 Plane = 10, Sphere, Cylinder, 
                 Wall = 100,
                 NoPRM = 4711 };

/** components of boundary faces */
class TBoundComp
{
  protected:
    /** component identifier */
    int ID;
    /** type of component */
    BoundTypes Type;
    /** TRUE if component is on free boundary */
    bool FreeBoundaryStatus;
    
    /** reference identifier */
    int refID;
    
  public:
    // Constructor
    TBoundComp(int id, int ref = -1);

    // Methods
    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat) = 0;

    /** return ID */
    int GetID() const
    { return ID; }

    /** get type of component */
    BoundTypes GetType() const
    { return Type; }

    /** get free boundary status */
    bool IsFreeBoundary() const
    { return FreeBoundaryStatus; }

    /** set free boundary status */
    void SetFreeBoundaryStatus(bool status)
    { FreeBoundaryStatus = status; }
    
    void ChangeType(BoundTypes New_Type)
     { Type = New_Type; }
     
     /** set reference */
    void SetRefID(int _ref)
    { refID=_ref; }

    /** return reference */
    int GetRefID() const
    { return refID; }

};

#endif
