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
// @(#)FEDatabase2D.h        1.5 02/08/00
//
// Class:       TFEDatabase
// Purpose:     store all used information for a FEM
//
// Author:      Gunar Matthies  09.07.98
//
// =======================================================================

#ifndef __FEDATABASE2D__
#define __FEDATABASE2D__

#include <Constants.h>
#include <Enumerations.h>

#include <BaseFunct1D.h>
#include <FEDesc1D.h>
#include <FE1D.h>

#include <BaseFunct2D.h>
#include <FEDesc2D.h>
#include <FE2D.h>
#include <FE2DMapper.h>
#include <FE2DMapper1Reg.h>
#include <HNDesc.h>
#include <RefTrans2D.h>

#include <QuadFormula1D.h>
#include <QuadFormulaQuad.h>
#include <QuadFormulaTria.h>


/** store all used information for a FEM */
class TFEDatabase2D
{
  protected:
//======================================================================
//      1D data
//======================================================================
    /** 1D (line) quadrature formulas */
    static TQuadFormula1D *QuadFormulas1D[N_QuadFormulas_1D];

    /** all sets of basis functions in 1D */
    static TBaseFunct1D *BaseFuncts1D[N_BaseFuncts1D];

    /** all sets of nodal functional */
    static TNodalFunctional1D *NodalFunctionals1D[N_NodalFunctionals1D];

    /** all descriptors for 1D finite elements */
    static TFEDesc1D *FEDescs1D[N_FEDescs1D];

    /** all 1D finite elements */
    static TFE1D *FEs1D[N_FEs1D];

    /** values of FE functions and their derivatives on the
        corresponding reference element */
    static double **RefElementValues1D[N_BaseFuncts1D][N_QuadFormulas_1D]
    [N_MultiIndices1D];

    /** values of FE functions and their derivatives on the
        current element */
    static double **OrigElementValues1D[N_BaseFuncts1D][N_MultiIndices1D];

    /** get line quadrature formula for given acuracy */
    static QuadFormula1D QFLineFromDegree[MAXDEGREE];

    /** highest accuracy for which a quadrature formula is available */
    static int HighestAccuracyLine;

    /** reference transformations */
    static TRefTrans1D *ReferenceTrans1D[1];

//======================================================================
//      2D data
//======================================================================
    /** all 2D finite elements */
    static TFE2D *FEs2D[N_FEs2D];

    /** all descriptors for 2D finite elements */
    static TFEDesc2D *FEDescs2D[N_FEDescs2D];

    /** all sets of basis functions in 2D */
    static TBaseFunct2D *BaseFuncts2D[N_BaseFuncts2D];

    /** all sets of nodal functional */
    static TNodalFunctional2D *NodalFunctionals2D[N_NodalFunctionals2D];

    /** quadrature formulas for 2D */
    static TQuadFormula2D *QuadFormulas2D[N_QuadFormulas_2D];

    /** 2D mapper for regular triangulation */
    static TFE2DMapper *FE2DMapper[N_FEDescs2D][N_FEDescs2D];

    /** 2D mapper for one regular triangulation */
    static TFE2DMapper1Reg *FE2DMapper1Reg[N_FEDescs2D][N_FEDescs2D];

    /** 2D hanging node descriptors */
    static THNDesc *HNDescs2D[N_HNDescs];

    /** values of FE functions and their derivatives on the
        corresponding reference element */
    static double **RefElementValues2D[N_BaseFuncts2D][N_QuadFormulas_2D]
    [N_MultiIndices2D];

    /** values of FE base functions on the edges of the corresponding
        reference element */
    static double **JointValues2D[N_BaseFuncts2D][N_QuadFormulas_1D]
    [MAXN_JOINTS];

    /** derivatives of FE base functions on the edges of the corresponding
        reference element */
    static double **JointDerivatives2D[N_BaseFuncts2D][N_QuadFormulas_1D]
    [MAXN_JOINTS][N_MultiIndices2D];

    /** values of FE functions and their derivatives on the
        current element */
    static double **OrigElementValues2D[N_BaseFuncts2D][N_MultiIndices2D];

    /** reference transformations */
    static TRefTrans2D *ReferenceTrans2D[N_RefTrans2D];

    /** prolongation matrix storage */
    static double *ProlongationMatrix2D[MaxN_BaseFunctions2D]
        [N_REFDESC][MaxN_BaseFunctions2D][MAXN_CHILDREN];

    /** function restriction matrix storage */
    static double *RestrictionMatrix2D[MaxN_BaseFunctions2D]
        [N_REFDESC][MaxN_BaseFunctions2D][MAXN_CHILDREN];

    /** get triangle quadrature formula for given acuracy */
    static QuadFormula2D QFTriaFromDegree[MAXDEGREE];

    /** highest accuracy for which a quadrature formula is available */
    static int HighestAccuracyTria;

    /** get quadrilateral quadrature formula for given acuracy */
    static QuadFormula2D QFQuadFromDegree[MAXDEGREE];

    /** highest accuracy for which a quadrature formula is available */
    static int HighestAccuracyQuad;

//======================================================================
//      2D arrays for easier access of information
//======================================================================
     /** Id of FEDesc from FE Id */ 
     static FEDesc2D FEDesc2D_IDFromFE2D[N_FEs2D];

     /** TFEDesc2D object from FE Id */ 
     static TFEDesc2D *FEDesc2DFromFE2D[N_FEs2D];

     /** Id of BaseFunct2D from FE Id */ 
     static BaseFunct2D BaseFunct2D_IDFromFE2D[N_FEs2D];

     /** number of base functions from FE Id */ 
     static int N_BaseFunctFromFE2D[N_FEs2D];

     /** polynomial degree of base functions from FE Id */ 
     static int PolynomialDegreeFromFE2D[N_FEs2D];

     /** accuracy of base functions from FE Id */ 
     static int AccuracyFromFE2D[N_FEs2D];

     /** TBaseFunct2DFEDesc object from FE Id */ 
     static TBaseFunct2D *BaseFunct2DFromFE2D[N_FEs2D];

     /** Id of NodalFunctional2D from FE Id */ 
     static NodalFunctional2D NodalFunctional2D_IDFromFE2D[N_FEs2D];

     /** TNodalFunctional2D object from FE Id */ 
     static TNodalFunctional2D *NodalFunctional2DFromFE2D[N_FEs2D];

     /** Id of RefTrans2D from FE Id */ 
     static RefTrans2D RefTrans2D_IDFromFE2D[N_FEs2D];

     /** reference element from FE Id */
     static BF2DRefElements RefElementFromFE2D[N_FEs2D];

//======================================================================
//      3D data
//======================================================================
    /** quadrature formulas for 3D */
    static TQuadFormula3D *QuadFormulas3D[N_QuadFormulas_3D];

//======================================================================
//      method only used by constructor
//======================================================================
  protected:
    /** register all known quadrature formulas into database */
    static void RegisterAllQuadFormulas();

    /** register all FE descriptors */
    static void RegisterAllFEDescs();

    /** register all base functions */
    static void RegisterAllBaseFunctions();

    /** register all nodal functionals */
    static void RegisterAllNodalFunctionals();

    /** register all finite elements */
    static void RegisterAllFEs();

    /** register all FE mappers */
    static void RegisterAllFEMappers();

    /** register all hanging node descriptors */
    static void RegisterAllHangingNodes();

    /** register all reference tranformations */
    static void RegisterAllRefTrans();

    /** generate some arrays form registered information */
    static void GenerateArrays();

//======================================================================
//      constructor
//======================================================================
  public:
    /** initialize the database */
    TFEDatabase2D();

//======================================================================
//      FE1D
//======================================================================
    /** return FE for given element */
    static TFE1D *GetFE1D (FE1D FE)
      { return FEs1D[FE]; };

    /** register FE1D for given element */
    static void RegisterFE1D(FE1D FE, TFE1D *element)
     { FEs1D[FE] = element; };

//======================================================================
//      FEDesc1D
//======================================================================
    /** return FEDesc1D for given element */
    static TFEDesc1D *GetFEDesc1D(FEDesc1D FEDesc)
      { return FEDescs1D[FEDesc]; };

    /** register FEDesc1D for given element */
    static void RegisterFEDesc1D(FEDesc1D FEDesc, TFEDesc1D *FEDesc1D)
      { FEDescs1D[FEDesc] = FEDesc1D; };

//======================================================================
//      BaseFunct1D
//======================================================================
    /** return BaseFunctions for given element */
    static TBaseFunct1D *GetBaseFunct1D(BaseFunct1D BaseFunct)
      { return BaseFuncts1D[BaseFunct]; };

    /** register BaseFunct1D for given element */
    static void RegisterBaseFunct1D(BaseFunct1D BaseFunct, 
                             TBaseFunct1D *BaseFunct1D)
      { BaseFuncts1D[BaseFunct] = BaseFunct1D; };

//======================================================================
//      NodalFunctional1D
//======================================================================
    /** return NodalFunctionals for given element */
    static TNodalFunctional1D *GetNodalFunctional1D
        (NodalFunctional1D NodalFunctional)
      { return NodalFunctionals1D[NodalFunctional]; };

    /** register NodalFunctional1D for given element */
    static void RegisterNodalFunctional1D
                (NodalFunctional1D NodalFunctional,
                 TNodalFunctional1D *NodalFunctional1D)
     { NodalFunctionals1D[NodalFunctional] = NodalFunctional1D; };

//======================================================================
//      FE function values and derivatives 1D
//======================================================================
    /** register FE function values or derivatives on ref element */
    static void RegisterRefElementValues
        (BaseFunct1D BaseFunct, QuadFormula1D QuadFormula,
         MultiIndex1D MultiIndex, double **Values)
    {
      RefElementValues1D[BaseFunct][QuadFormula][MultiIndex]
        = Values;
    }

    /** return requested FE function values or derivatives */
    static double **GetRefElementValues
        (BaseFunct1D BaseFunct, QuadFormula1D QuadFormula,
         MultiIndex1D MultiIndex)
    {
      return RefElementValues1D[BaseFunct][QuadFormula][MultiIndex];
    }

    /** register FE function values or derivatives on current element */
    static void RegisterOrigElementValues
        (BaseFunct1D BaseFunct, MultiIndex1D MultiIndex, double **Values)
    {
      OrigElementValues1D[BaseFunct][MultiIndex] = Values;
    }

    /** return requested FE function values or derivatives */
    static double **GetOrigElementValues
        (BaseFunct1D BaseFunct, MultiIndex1D MultiIndex)
    {
      return OrigElementValues1D[BaseFunct][MultiIndex];
    }

     /** return reference transformation */
     static TRefTrans1D *GetRefTrans1D(RefTrans1D reftrans)
     {
       return ReferenceTrans1D[reftrans];
     }


//======================================================================
//      FEDesc2D
//======================================================================
    /** return FEDesc2D for given element */
    static TFEDesc2D *GetFEDesc2D(FEDesc2D FEDesc)
      { return FEDescs2D[FEDesc]; };

    /** register FEDesc2D for given element */
    static void RegisterFEDesc2D(FEDesc2D FEDesc, TFEDesc2D *FEDesc2D)
      { FEDescs2D[FEDesc] = FEDesc2D; };

//======================================================================
//      BaseFunct2D
//======================================================================
    /** return BaseFunctions for given element */
    static TBaseFunct2D *GetBaseFunct2D(BaseFunct2D BaseFunct)
      { return BaseFuncts2D[BaseFunct]; };

    /** register BaseFunct2D for given element */
    static void RegisterBaseFunct2D(BaseFunct2D BaseFunct, 
                             TBaseFunct2D *BaseFunct2D)
      { BaseFuncts2D[BaseFunct] = BaseFunct2D; };

//======================================================================
//      NodalFunctional2D
//======================================================================
    /** return NodalFunctionals for given element */
    static TNodalFunctional2D *GetNodalFunctional2D
        (NodalFunctional2D NodalFunctional)
      { return NodalFunctionals2D[NodalFunctional]; };

    /** register NodalFunctional2D for given element */
    static void RegisterNodalFunctional2D
                (NodalFunctional2D NodalFunctional,
                 TNodalFunctional2D *NodalFunctional2D)
     { NodalFunctionals2D[NodalFunctional] = NodalFunctional2D; };

//======================================================================
//      FE2D
//======================================================================
    /** return FE for given element */
    static TFE2D *GetFE2D(FE2D FE)
      { return FEs2D[FE]; };

    /** register FE2D for given element */
    static void RegisterFE2D(FE2D FE, TFE2D *element)
     { FEs2D[FE] = element; };

//======================================================================
//      access to some arrays
//======================================================================
    /** return Id of FEDesc from FE Id */
    static FEDesc2D GetFEDesc2D_IDFromFE2D(FE2D ele)
      { return FEDesc2D_IDFromFE2D[ele]; };

    /** return TFEDesc2D object from FE Id */ 
    static TFEDesc2D *GetFEDesc2DFromFE2D(FE2D ele)
      { return FEDesc2DFromFE2D[ele]; };

    /** return Id of BaseFunct2D from FE Id */ 
    static BaseFunct2D GetBaseFunct2D_IDFromFE2D(FE2D ele)
      { return BaseFunct2D_IDFromFE2D[ele]; };

    /** return number of base functions from FE Id */
    static int GetN_BaseFunctFromFE2D(FE2D ele)
      { return N_BaseFunctFromFE2D[ele]; };

    /** return polynomial degree of base functions from FE Id */
    static int GetPolynomialDegreeFromFE2D(FE2D ele)
      { return PolynomialDegreeFromFE2D[ele]; };

    /** return accuracy of base functions from FE Id */
    static int GetAccuracyFromFE2D(FE2D ele)
      { return AccuracyFromFE2D[ele]; };

    /** return TBaseFunct2DFEDesc object from FE Id */
    static TBaseFunct2D *GetBaseFunct2DFromFE2D(FE2D ele)
      { return BaseFunct2DFromFE2D[ele]; };

    /** return Id of NodalFunctional2D from FE Id */
    static NodalFunctional2D GetNodalFunctional2D_IDFromFE2D(FE2D ele)
      { return NodalFunctional2D_IDFromFE2D[ele]; };

    /** return TNodalFunctional2D object from FE Id */
    static TNodalFunctional2D *GetNodalFunctional2DFromFE2D(FE2D ele)
      { return NodalFunctional2DFromFE2D[ele]; };

    /** return Id of RefTrans2D from FE Id */
    static RefTrans2D GetRefTrans2D_IDFromFE2D(FE2D ele)
      { return RefTrans2D_IDFromFE2D[ele]; };

    /** return reference element from FE Id */
    static BF2DRefElements GetRefElementFromFE2D(FE2D ele)
      { return RefElementFromFE2D[ele]; };


    /** return array of Id of FEDesc from FE Id */
    static FEDesc2D *GetFEDesc2D_IDFromFE2D()
      { return FEDesc2D_IDFromFE2D; };

    /** return array of TFEDesc2D object from FE Id */ 
    static TFEDesc2D **GetFEDesc2DFromFE2D()
      { return FEDesc2DFromFE2D; };

    /** return array of Id of BaseFunct2D from FE Id */ 
    static BaseFunct2D *GetBaseFunct2D_IDFromFE2D()
      { return BaseFunct2D_IDFromFE2D; };

    /** return array of number of base functions from FE Id */
    static int *GetN_BaseFunctFromFE2D()
      { return N_BaseFunctFromFE2D; };

    /** return array of polynomial degree of base functions from FE Id */
    static int *GetPolynomialDegreeFromFE2D()
      { return PolynomialDegreeFromFE2D; };

    /** return array of accuracy of base functions from FE Id */
    static int *GetAccuracyFromFE2D()
      { return AccuracyFromFE2D; };

    /** return array of TBaseFunct2DFEDesc object from FE Id */
    static TBaseFunct2D **GetBaseFunct2DFromFE2D()
      { return BaseFunct2DFromFE2D; };

    /** return array of Id of NodalFunctional2D from FE Id */
    static NodalFunctional2D *GetNodalFunctional2D_IDFromFE2D()
      { return NodalFunctional2D_IDFromFE2D; };

    /** return array of TNodalFunctional2D object from FE Id */
    static TNodalFunctional2D **GetNodalFunctional2DFromFE2D()
      { return NodalFunctional2DFromFE2D; };

    /** return array of Id of RefTrans2D from FE Id */
    static RefTrans2D *GetRefTrans2D_IDFromFE2D()
      { return RefTrans2D_IDFromFE2D; };

    /** return array of reference elements from FE Id */
    static BF2DRefElements *GetRefElementFromFE2D()
      { return RefElementFromFE2D; };

//======================================================================
//      prolongation matrices in 2D
//======================================================================
    /** return prolongation matrix for given situation */
    static double *GetProlongationMatrix2D (FE2D parent, 
        Refinements refine, FE2D child, int childnumber);

    /** register prolongation matrix for given situation */
    static void RegisterProlongationMatrix2D(BaseFunct2D parent, 
        Refinements refine, BaseFunct2D child, int childnumber, 
        double *matrix)
     { ProlongationMatrix2D[parent][refine][child][childnumber] = 
                matrix; };

//======================================================================
//      function restriction matrices in 2D
//======================================================================
    /** return function restriction matrix for given situation */
    static double *GetRestrictionMatrix2D (FE2D parent, 
        Refinements refine, FE2D child, int childnumber);

    /** register function restriction matrix for given situation */
    static void RegisterRestrictionMatrix2D(BaseFunct2D parent, 
        Refinements refine, BaseFunct2D child, int childnumber, 
        double *matrix)
     { RestrictionMatrix2D[parent][refine][child][childnumber] = 
                matrix; };

//======================================================================
//      FE mapper, regular
//======================================================================
    /** return FE mapper */
    static TFE2DMapper *GetFE2DMapper(FEDesc2D FE1, FEDesc2D FE2)
      { return FE2DMapper[FE1][FE2]; };

    /** register FE2D for given element */
    static void RegisterFE2DMapper(FEDesc2D FE1, FEDesc2D FE2, 
                                   TFE2DMapper *mapper)
     { FE2DMapper[FE1][FE2] = mapper; };

//======================================================================
//      FE mapper, one regular
//======================================================================
    /** return FE mapper */
    static TFE2DMapper1Reg *GetFE2DMapper1Reg (FEDesc2D FE1, FEDesc2D FE2)
      { return FE2DMapper1Reg[FE1][FE2]; };

    /** register FE2D for given element */
    static void RegisterFE2DMapper1Reg(FEDesc2D FE1, FEDesc2D FE2, 
                                   TFE2DMapper1Reg *mapper)
     { FE2DMapper1Reg[FE1][FE2] = mapper; };

//======================================================================
//      HNDesc
//======================================================================
    /** return HNDesc for given element */
    static THNDesc *GetHNDesc2D(HNDesc Desc)
      { return HNDescs2D[Desc]; };

    /** register HNDesc for given element */
    static void RegisterHNDesc2D(HNDesc Desc, THNDesc *HNDesc_Obj)
      { HNDescs2D[Desc] = HNDesc_Obj; };

//======================================================================
//      FE function values and derivatives 2D
//======================================================================
    /** register FE function values or derivatives on ref element */
    static void RegisterRefElementValues
        (BaseFunct2D BaseFunct, QuadFormula2D QuadFormula,
         MultiIndex2D MultiIndex, double **Values)
    {
      RefElementValues2D[BaseFunct][QuadFormula][MultiIndex]
        = Values;
    }

    /** return requested FE function values or derivatives */
    static double **GetRefElementValues
        (BaseFunct2D BaseFunct, QuadFormula2D QuadFormula,
         MultiIndex2D MultiIndex)
    {
      return RefElementValues2D[BaseFunct][QuadFormula][MultiIndex];
    }

    /** register joint values of FE base function on ref element */
    static void RegisterJointValues2D
        (BaseFunct2D BaseFunct, QuadFormula1D formula, int joint, 
         double **Values)
    {
      JointValues2D[BaseFunct][formula][joint] = Values;
    }

    /** return requested joint values of FE base functions */
    static double **GetJointValues2D
        (BaseFunct2D BaseFunct, QuadFormula1D formula, int joint)
    {
      return JointValues2D[BaseFunct][formula][joint];
    }

    /** register joint derivatives of FE base function on ref element */
    static void RegisterJointDerivatives2D
        (BaseFunct2D BaseFunct, QuadFormula1D formula, int joint, 
         MultiIndex2D MultiIndex, double **Values)
    {
      JointDerivatives2D[BaseFunct][formula][joint][MultiIndex] = Values;
    }

    /** return requested joint values of FE base functions */
    static double **GetJointDerivatives2D
        (BaseFunct2D BaseFunct, QuadFormula1D formula, int joint,
        MultiIndex2D MultiIndex)
    {
      return JointDerivatives2D[BaseFunct][formula][joint][MultiIndex];
    }

    /** register FE function values or derivatives on current element */
    static void RegisterOrigElementValues
        (BaseFunct2D BaseFunct, MultiIndex2D MultiIndex, double **Values)
    {
      OrigElementValues2D[BaseFunct][MultiIndex] = Values;
    }

    /** return requested FE function values or derivatives */
    static double **GetOrigElementValues
        (BaseFunct2D BaseFunct, MultiIndex2D MultiIndex)
    {
      return OrigElementValues2D[BaseFunct][MultiIndex];
    }

    /** calculate the values of base functions and their derivatives
        on the original element */
    /*
    static RefTrans2D GetOrigValues(TBaseCell *cell, TFE2D *element, 
                                    int N_Points, 
                                    double *xi, double *eta,
                                    int N_Functs,
                                    BaseFunct2D BaseFunct,
                                    QuadFormula2D QuadFormula);
    */

    /** calculate functions and derivatives from reference element
        to original element */
    static void GetOrigValues(RefTrans2D RefTrans,
                double xi, double eta, TBaseFunct2D *bf,
                TCollection *Coll, TGridCell *cell,
                double *uref, double *uxiref, double *uetaref,
                double *uorig, double *uxorig, double *uyorig);
    
    static void GetOrigValues(RefTrans2D RefTrans,
                              double zeta, TBaseFunct2D *bf, int edgeNumber,
                              double *uref, double *uxiref, double *uetaref,
                              double *uorig, double *uxorig, double *uyorig);

//======================================================================
//      reference transformation 2D
//======================================================================
     /** return reference transformation */
     static TRefTrans2D *GetRefTrans2D(RefTrans2D reftrans)
     {
       return ReferenceTrans2D[reftrans];
     }

     /** calculate points on original element */
     static void GetOrigFromRef(RefTrans2D RefTrans, int n_points, 
                        double *xi, double *eta,
                        double *X, double *Y, double *absdetjk);

     /** calculate base functions with derivatives and coordinates
         from reference to original element */
     static RefTrans2D GetOrig(int N_LocalUsedElements, FE2D *LocalUsedElements,
                         TCollection *Coll,
                         TBaseCell *cell, bool *Needs2ndDer,
                         int &N_Points, double* &xi, double* &eta, 
                         double* &weights, double* X, double* Y,
                         double* absdetjk);

     /** calculate points on reference element */
     static void GetRefFromOrig(RefTrans2D RefTrans,
                        double X, double Y,
                        double &xi, double &eta);

     /** set cell for reference transformation */
     static void SetCellForRefTrans(TBaseCell *cell,
                                    RefTrans2D reftrans);

//======================================================================
//      QuadFormula1D
//======================================================================
     /** return QuadFormula1D */
     static TQuadFormula1D *GetQuadFormula1D(QuadFormula1D QF)
       { return QuadFormulas1D[QF]; };

     /** register QuadFormula1D */
     static void RegisterQuadFormula1D(QuadFormula1D QF,
                                 TQuadFormula1D *QuadForm)
       { QuadFormulas1D[QF] = QuadForm; };

    /** get line quadrature formula for given acuracy */
    static QuadFormula1D GetQFLineFromDegree(int accuracy)
    {
      if(accuracy<=HighestAccuracyLine)
         return QFLineFromDegree[accuracy];
      else
        return QFLineFromDegree[HighestAccuracyLine];
    };

//======================================================================
//      QuadFormula2D
//======================================================================
     /** return QuadFormula2D */
     static TQuadFormula2D *GetQuadFormula2D(QuadFormula2D QF)
       { return QuadFormulas2D[QF]; };

     /** register QuadFormula2D */
     static void RegisterQuadFormula2D(QuadFormula2D QF,
                                 TQuadFormula2D *QuadForm)
       { QuadFormulas2D[QF] = QuadForm; };

    /** get triangle quadrature formula for given acuracy */
    static QuadFormula2D GetQFTriaFromDegree(int accuracy)
    {
      if(accuracy<=HighestAccuracyTria)
         return QFTriaFromDegree[accuracy];
      else
        return QFTriaFromDegree[HighestAccuracyTria];
    };
    /** get quadrilateral quadrature formula for given acuracy */
    static QuadFormula2D GetQFQuadFromDegree(int accuracy)
    {
      if(accuracy<=HighestAccuracyQuad)
         return QFQuadFromDegree[accuracy];
      else
        return QFQuadFromDegree[HighestAccuracyQuad];
    };
    
    /** get quadrature formula on given reference element */
    static QuadFormula2D GetQFFromDegree(int accuracy, BF2DRefElements RefElem);


//======================================================================
//      QuadFormula3D
//======================================================================
     /** return QuadFormula3D */
     static TQuadFormula3D *GetQuadFormula3D(QuadFormula3D QF)
       { return QuadFormulas3D[QF]; };

     /** register QuadFormula3D */
     static void RegisterQuadFormula3D(QuadFormula3D QF,
                                 TQuadFormula3D *QuadForm)
       { QuadFormulas3D[QF] = QuadForm; };
  
};

#endif
