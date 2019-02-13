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
// @(#)FESpace2D.h        1.5 06/13/00
// 
// Class:       TFESpace2D
// Purpose:     class for all 2D finite element spaces
//
// Author:      Gunar Matthies (03.11.97)
//
// History:     start of implementation 03.11.97 (Gunar Matthies)
//
//              split FESpace into TFESpacexD
//              15.04.1998 (Volker Behns)
//
//              start of reimplementation
//              30.07.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __FESPACE2D__
#define __FESPACE2D__

#include <FESpace.h>
#include <FE2D.h>

class THangingNode;

/** class for all 2D finite element spaces */
class TFESpace2D : public TFESpace
{
  protected:
    /** number of active degrees of freedom */
    int N_ActiveDegrees;

    /** number of slave degrees of freedom */
    int N_SlaveDegrees;

    /** NeumannBound <= i < HangingBound for all hanging nodes i */
    /** => HangingBound <= i < DirichletBound for all Dirichlet node i */
    int HangingBound;

    /** array containing the used elements */
    FE2D *UsedElements; 

    /** array with an element for each shape */
    FE2D *ElementForShape;

    /** array of hanging nodes */
    THangingNode **HangingNodeArray;

    /** array storing the fe for each element, if necessary */
    FE2D *AllElements;
    
    /** boundary condition used to create this space */
   BoundCondFunct2D *BoundCondition;

#ifdef __MORTAR__
    /** 1D collection of mortar cells */
    TCollection *MortarColl;
#endif

   /** no. points in internal coordinate direction --- ADI scheme*/
   int N_IntlPts;

   /** internal points' X coordinate --- ADI scheme*/
   double *X_Intl;

   /** internal points' Y coordinate --- ADI scheme*/
   double *Y_Intl;

   /** indices for mapping between Nodalfunctional/Nodal-interpolation point 
      in operator-splitting methods ---  Sashikumaar Ganesan */
   int *IntlPtIndexOfPts;


  public:
    /** constructor */
    TFESpace2D(TCollection *coll, char *name, char *description);

    /** constructor for building a space with elements of order k */
    // if no mortar is used, the last argument can be set to NULL
    TFESpace2D(TCollection *coll, char *name, char *description, 
               BoundCondFunct2D *BoundaryCondition, int k,
               TCollection *mortarcoll);

    TFESpace2D(TCollection *coll, char *name, char *description, 
               BoundCondFunct2D *BoundaryCondition, SpaceType type,
               int k, TCollection *mortarcoll);

    /** constructor for building a space with the given elements */
    TFESpace2D(TCollection *coll, char *name, char *description,
               BoundCondFunct2D *BoundaryCondition,
               FE2D *fes, TCollection *mortarcoll);

    /** destructor */
    ~TFESpace2D();

    /** find used elements */
    void FindUsedElements();

    /** construct space */
    void ConstructSpace(BoundCondFunct2D *BoundaryCondition);

    /** return number of active degrees of freedom */
    int GetN_ActiveDegrees() const
    { return N_ActiveDegrees; }

    /** return number of slave degrees of freedom */
    int GetN_SlaveDegrees() const
    { return N_SlaveDegrees; }

    /** return HangingBound */
    int GetHangingBound() const
    { return HangingBound; }

    /** return N_Hanging=N_SlaveDegrees */
    int GetN_Hanging() const
    { return N_SlaveDegrees; }

    /** return identifiers of used elements */
    FE2D *GetUsedElements() const
    { return UsedElements; }

    /** return array with all hanging nodes */
    THangingNode **GetHangingNodes() const
    { return HangingNodeArray; }

    /** return the FE Id for element i, corresponding to cell */
    FE2D GetFE2D(int i, TBaseCell *cell);

    /** return position of one given DOF */
    void GetDOFPosition(int dof, double &x, double &y);

    /** return position of all dofs */
    void GetDOFPosition(double *x, double *y);

    void SetIntlPtIndexOfPts(int *intlPtIndexOfPts)
     { IntlPtIndexOfPts = intlPtIndexOfPts; }

    int *GetIntlPtIndexOfPts() const
     { return IntlPtIndexOfPts; }
     
     FE2D *GetAllElements() const
     { return AllElements; }
     
     /** return boundary condition */
    BoundCondFunct2D *GetBoundCondition() const
    { return BoundCondition; }
    
    friend  bool operator== (const TFESpace2D &lhs, const TFESpace2D &rhs);
    friend  bool operator!= (const TFESpace2D &lhs, const TFESpace2D &rhs);
};

#ifdef __MORTAR__
  double GetLambda(double , double , TVertex *, double , double );
#endif // __MORTAR__

#endif
