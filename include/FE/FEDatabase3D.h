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
// Class:       TFEDatabase3D
// Purpose:     store all used information for a FEM for 3D
//
// Author:      Gunar Matthies  26.11.99
//
// =======================================================================

#ifndef __FEDATABASE3D__
#define __FEDATABASE3D__

#include <Constants.h>
#include <Enumerations.h>

#include <BaseFunct3D.h>
#include <FEDesc3D.h>
#include <FE3D.h>
#include <FE3DMapper1Reg.h>

#include <QuadFormula1D.h>
#include <QuadFormulaHexa.h>
#include <QuadFormulaTetra.h>
#include <QuadFormulaTria.h>
#include <QuadFormulaQuad.h>

/** store all used information for a FEM */
class TFEDatabase3D
{
  protected:
//======================================================================
//      1D data
//======================================================================
    /** 1D (line) quadrature formulas */
    static TQuadFormula1D *QuadFormulas1D[N_QuadFormulas_1D];

    /** get line quadrature formula for given acuracy */
    static QuadFormula1D QFLineFromDegree[MAXDEGREE];

    /** highest accuracy for which a quadrature formula is available */
    static int HighestAccuracyLine;

//======================================================================
//      2D arrays for easier access of information
//======================================================================
    /** quadrature formulas for 2D */
    static TQuadFormula2D *QuadFormulas2D[N_QuadFormulas_2D];

    /** get triangle quadrature formula for given acuracy */
    static QuadFormula2D QFTriaFromDegree[MAXDEGREE];

    /** get quadrilateral quadrature formula for given acuracy */
    static QuadFormula2D QFQuadFromDegree[MAXDEGREE];

//======================================================================
//      3D arrays for easier access of information
//======================================================================
     /** get tetrahedron quadrature formula for given acuracy */
     static QuadFormula3D QFTetraFromDegree[MAXDEGREE];

     /** get hexahedron quadrature formula for given acuracy */
     static QuadFormula3D QFHexaFromDegree[MAXDEGREE];

     /** get hexahedron quadrature formula for given acuracy */
     static QuadFormula3D QFConvolutionHexaFromDegree[MAXDEGREE];

     /** Id of FEDesc from FE Id */ 
     static FEDesc3D FEDesc3D_IDFromFE3D[N_FEs3D];

     /** TFEDesc3D object from FE Id */ 
     static TFEDesc3D *FEDesc3DFromFE3D[N_FEs3D];

     /** Id of BaseFunct3D from FE Id */ 
     static BaseFunct3D BaseFunct3D_IDFromFE3D[N_FEs3D];

     /** number of base functions from FE Id */ 
     static int N_BaseFunctFromFE3D[N_FEs3D];

     /** polynomial degree of base functions from FE Id */ 
     static int PolynomialDegreeFromFE3D[N_FEs3D];

     /** accuracy of base functions from FE Id */ 
     static int AccuracyFromFE3D[N_FEs3D];

     /** TBaseFunct3DFEDesc object from FE Id */ 
     static TBaseFunct3D *BaseFunct3DFromFE3D[N_FEs3D];

     /** Id of NodalFunctional3D from FE Id */ 
     static NodalFunctional3D NodalFunctional3D_IDFromFE3D[N_FEs3D];

     /** TNodalFunctional3D object from FE Id */ 
     static TNodalFunctional3D *NodalFunctional3DFromFE3D[N_FEs3D];

     /** Id of RefTrans3D from FE Id */ 
     static RefTrans3D RefTrans3D_IDFromFE3D[N_FEs3D];

     /** reference element from FE Id */
     static BF3DRefElements RefElementFromFE3D[N_FEs3D];

//======================================================================
//      3D data
//======================================================================
    /** quadrature formulas for 3D */
    static TQuadFormula3D *QuadFormulas3D[N_QuadFormulas_3D];

    /** all 3D finite elements */
    static TFE3D *FEs3D[N_FEs3D];

    /** all descriptors for 3D finite elements */
    static TFEDesc3D *FEDescs3D[N_FEDescs3D];

    /** all sets of basis functions in 3D */
    static TBaseFunct3D *BaseFuncts3D[N_BaseFuncts3D];

    /** all sets of nodal functional */
    static TNodalFunctional3D *NodalFunctionals3D[N_NodalFunctionals3D];

    /** 3D mapper for regular triangulation */
    static TFE3DMapper *FE3DMapper[N_FEDescs3D][N_FEDescs3D];

    /** 3D mapper for regular triangulation */
    static TFE3DMapper1Reg *FE3DMapper1Reg[N_FEDescs3D][N_FEDescs3D];

    /** 3D hanging node descriptors */
    static THNDesc *HNDescs3D[N_HNDescs];

    /** reference transformations */
    static TRefTrans3D *ReferenceTrans3D[N_RefTrans3D];

    /** prolongation matrix storage */
    static double *ProlongationMatrix3D[MaxN_BaseFunctions3D]
        [N_REFDESC][MaxN_BaseFunctions3D][MAXN_CHILDREN];

    /** function restriction matrix storage */
    static double *RestrictionMatrix3D[MaxN_BaseFunctions3D]
        [N_REFDESC][MaxN_BaseFunctions3D][MAXN_CHILDREN];

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

    /** values of FE functions and their derivatives on the
        corresponding reference element */
    static double **RefElementValues3D[N_BaseFuncts3D][N_QuadFormulas_3D]
      [N_MultiIndices3D];
    
    /** values of FE base functions on the edges of the corresponding
        reference element */
    static double **JointValues3D[N_BaseFuncts3D][N_QuadFormulas_2D]
      [MAXN_JOINTS];

    /** derivatives of FE base functions on the edges of the corresponding
        reference element */
    static double **JointDerivatives3D[N_BaseFuncts3D][N_QuadFormulas_2D]
    [MAXN_JOINTS][N_MultiIndices3D];

    /** values of FE functions and their derivatives on the
        current element */
    static double **OrigElementValues3D[N_BaseFuncts3D][N_MultiIndices3D];


//======================================================================
//      constructor
//======================================================================
  public:
    /** initialize the database */
    TFEDatabase3D();
//======================================================================
//      QuadFormula1D
//======================================================================
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
       { return QFTriaFromDegree[accuracy]; };

    /** get quadrilateral quadrature formula for given acuracy */
    static QuadFormula2D GetQFQuadFromDegree(int accuracy)
       { return QFQuadFromDegree[accuracy]; };

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

     /** get tetrahedron quadrature formula for given acuracy */
     static QuadFormula3D GetQFTetraFromDegree(int accuracy)
       { return QFTetraFromDegree[accuracy]; };
     
     /** get hexahedron quadrature formula for given acuracy */
     static QuadFormula3D GetQFHexaFromDegree(int accuracy)
       { return QFHexaFromDegree[accuracy]; };   

     /** get hexahedron quadrature formula for convolution */
     static QuadFormula3D GetQFConvolutionHexaFromDegree(int accuracy)
       { return QFConvolutionHexaFromDegree[accuracy]; };   

//======================================================================
//      FEDesc3D
//======================================================================
    /** return FEDesc3D for given element */
    static TFEDesc3D *GetFEDesc3D(FEDesc3D FEDesc)
      { return FEDescs3D[FEDesc]; };

    /** register FEDesc3D for given element */
    static void RegisterFEDesc3D(FEDesc3D FEDesc, TFEDesc3D *FEDesc3D)
      { FEDescs3D[FEDesc] = FEDesc3D; };

//======================================================================
//      BaseFunct3D
//======================================================================
    /** return BaseFunctions for given element */
    static TBaseFunct3D *GetBaseFunct3D(BaseFunct3D BaseFunct)
      { return BaseFuncts3D[BaseFunct]; };

    /** register BaseFunct3D for given element */
    static void RegisterBaseFunct3D(BaseFunct3D BaseFunct, 
                             TBaseFunct3D *BaseFunct3D)
      { BaseFuncts3D[BaseFunct] = BaseFunct3D; };

//======================================================================
//      NodalFunctional3D
//======================================================================
    /** return NodalFunctionals for given element */
    static TNodalFunctional3D *GetNodalFunctional3D
        (NodalFunctional3D NodalFunctional)
      { return NodalFunctionals3D[NodalFunctional]; };

    /** register NodalFunctional3D for given element */
    static void RegisterNodalFunctional3D
                (NodalFunctional3D NodalFunctional,
                 TNodalFunctional3D *NodalFunctional3D)
     { NodalFunctionals3D[NodalFunctional] = NodalFunctional3D; };

//======================================================================
//      FE3D
//======================================================================
    /** return FE for given element */
    static TFE3D *GetFE3D (FE3D FE)
      { return FEs3D[FE]; };

    /** register FE3D for given element */
    static void RegisterFE3D(FE3D FE, TFE3D *element)
     { FEs3D[FE] = element; };

//======================================================================
//      access to some arrays
//======================================================================
    /** return Id of FEDesc from FE Id */
    static FEDesc3D GetFEDesc3D_IDFromFE3D(FE3D ele)
      { return FEDesc3D_IDFromFE3D[ele]; };

    /** return TFEDesc3D object from FE Id */ 
    static TFEDesc3D *GetFEDesc3DFromFE3D(FE3D ele)
      { return FEDesc3DFromFE3D[ele]; };

    /** return Id of BaseFunct3D from FE Id */ 
    static BaseFunct3D GetBaseFunct3D_IDFromFE3D(FE3D ele)
      { return BaseFunct3D_IDFromFE3D[ele]; };

    /** return number of base functions from FE Id */
    static int GetN_BaseFunctFromFE3D(FE3D ele)
      { return N_BaseFunctFromFE3D[ele]; };

    /** return polynomial degree of base functions from FE Id */
    static int GetPolynomialDegreeFromFE3D(FE3D ele)
      { return PolynomialDegreeFromFE3D[ele]; };

    /** return accuracy of base functions from FE Id */
    static int GetAccuracyFromFE3D(FE3D ele)
      { return AccuracyFromFE3D[ele]; };

    /** return TBaseFunct3DFEDesc object from FE Id */
    static TBaseFunct3D *GetBaseFunct3DFromFE3D(FE3D ele)
      { return BaseFunct3DFromFE3D[ele]; };

    /** return Id of NodalFunctional3D from FE Id */
    static NodalFunctional3D GetNodalFunctional3D_IDFromFE3D(FE3D ele)
      { return NodalFunctional3D_IDFromFE3D[ele]; };

    /** return TNodalFunctional3D object from FE Id */
    static TNodalFunctional3D *GetNodalFunctional3DFromFE3D(FE3D ele)
      { return NodalFunctional3DFromFE3D[ele]; };

    /** return Id of RefTrans3D from FE Id */
    static RefTrans3D GetRefTrans3D_IDFromFE3D(FE3D ele)
      { return RefTrans3D_IDFromFE3D[ele]; };

    /** return reference element from FE Id */
    static BF3DRefElements GetRefElementFromFE3D(FE3D ele)
      { return RefElementFromFE3D[ele]; };


    /** return array of Id of FEDesc from FE Id */
    static FEDesc3D *GetFEDesc3D_IDFromFE3D()
      { return FEDesc3D_IDFromFE3D; };

    /** return array of TFEDesc3D object from FE Id */ 
    static TFEDesc3D **GetFEDesc3DFromFE3D()
      { return FEDesc3DFromFE3D; };

    /** return array of Id of BaseFunct3D from FE Id */ 
    static BaseFunct3D *GetBaseFunct3D_IDFromFE3D()
      { return BaseFunct3D_IDFromFE3D; };

    /** return array of number of base functions from FE Id */
    static int *GetN_BaseFunctFromFE3D()
      { return N_BaseFunctFromFE3D; };

    /** return array of polynomial degree of base functions from FE Id */
    static int *GetPolynomialDegreeFromFE3D()
      { return PolynomialDegreeFromFE3D; };

    /** return array of accuracy of base functions from FE Id */
    static int *GetAccuracyFromFE3D()
      { return AccuracyFromFE3D; };

    /** return array of TBaseFunct3DFEDesc object from FE Id */
    static TBaseFunct3D **GetBaseFunct3DFromFE3D()
      { return BaseFunct3DFromFE3D; };

    /** return array of Id of NodalFunctional3D from FE Id */
    static NodalFunctional3D *GetNodalFunctional3D_IDFromFE3D()
      { return NodalFunctional3D_IDFromFE3D; };

    /** return array of TNodalFunctional3D object from FE Id */
    static TNodalFunctional3D **GetNodalFunctional3DFromFE3D()
      { return NodalFunctional3DFromFE3D; };

    /** return array of Id of RefTrans3D from FE Id */
    static RefTrans3D *GetRefTrans3D_IDFromFE3D()
      { return RefTrans3D_IDFromFE3D; };

    /** return array of reference elements from FE Id */
    static BF3DRefElements *GetRefElementFromFE3D()
      { return RefElementFromFE3D; };

//======================================================================
//      prolongation matrices in 3D
//======================================================================
    /** return prolongation matrix for given situation */
    static double *GetProlongationMatrix3D (FE3D parent, 
        Refinements refine, FE3D child, int childnumber);

    /** register prolongation matrix for given situation */
    static void RegisterProlongationMatrix3D(BaseFunct3D parent, 
        Refinements refine, BaseFunct3D child, int childnumber, 
        double *matrix)
     { ProlongationMatrix3D[parent][refine][child][childnumber] = 
                matrix; };

//======================================================================
//      function restriction matrices in 3D
//======================================================================
    /** return function restriction matrix for given situation */
    static double *GetRestrictionMatrix3D (FE3D parent, 
        Refinements refine, FE3D child, int childnumber);

    /** register function restriction matrix for given situation */
    static void RegisterRestrictionMatrix3D(BaseFunct3D parent, 
        Refinements refine, BaseFunct3D child, int childnumber, 
        double *matrix)
     { RestrictionMatrix3D[parent][refine][child][childnumber] = 
                matrix; };

//======================================================================
//      FE mapper 3D, regular
//======================================================================
    /** return FE mapper */
    static TFE3DMapper *GetFE3DMapper (FEDesc3D FE1, FEDesc3D FE2)
      { return FE3DMapper[FE1][FE2]; };

    /** register FE3DMapper for given element */
    static void RegisterFE3DMapper(FEDesc3D FE1, FEDesc3D FE2, 
                                   TFE3DMapper *mapper)
     { FE3DMapper[FE1][FE2] = mapper; };

//======================================================================
//      FE mapper 3D, 1-regular
//======================================================================
    /** return FE mapper */
    static TFE3DMapper1Reg *GetFE3DMapper1Reg (FEDesc3D FE1, FEDesc3D FE2)
      { return FE3DMapper1Reg[FE1][FE2]; };

    /** register FE3DMapper1Reg for given element */
    static void RegisterFE3DMapper1Reg(FEDesc3D FE1, FEDesc3D FE2, 
                                       TFE3DMapper1Reg *mapper)
     { FE3DMapper1Reg[FE1][FE2] = mapper; };

    static void SetCellForRefTrans(TBaseCell *cell, RefTrans3D reftrans);

//======================================================================
//      HNDesc
//======================================================================
    /** return HNDesc for given element */
    static THNDesc *GetHNDesc3D(HNDesc Desc)
      { return HNDescs3D[Desc]; }

    /** register HNDesc for given element */
    static void RegisterHNDesc3D(HNDesc Desc, THNDesc *HNDesc_Obj)
      { HNDescs3D[Desc] = HNDesc_Obj; }

//======================================================================
//      reference transformation 3D
//======================================================================
     /** return reference transformation */
     static TRefTrans3D *GetRefTrans3D(RefTrans3D reftrans)
     {
       return ReferenceTrans3D[reftrans];
     }

     /** calculate points on original element */
     static void GetOrigFromRef(RefTrans3D RefTrans, int n_points, 
                        double *xi, double *eta, double *zeta,
                        double *X, double *Y, double *Z, double *absdetjk);

     /** calculate base functions with derivatives and coordinates
         from reference to original element */
     static RefTrans3D GetOrig(int N_LocalUsedElements, FE3D *LocalUsedElements,
                         TCollection *Coll, TBaseCell *cell,
                         bool *Needs2ndDer,
                         int &N_Points, double* &xi, double* &eta, double* &zeta, 
                         double* &weights, double* X, double* Y, double* Z,
                         double* absdetjk);

     /** calculate points on reference element */
     static void GetRefFromOrig(RefTrans3D RefTrans,
                        double X, double Y, double Z,
                        double &xi, double &eta, double &zeta);

//======================================================================
//      FE function values and derivatives 3D
//======================================================================
    /** register FE function values or derivatives on ref element */
    static void RegisterRefElementValues
        (BaseFunct3D BaseFunct, QuadFormula3D QuadFormula,
         MultiIndex3D MultiIndex, double **Values)
    {
      RefElementValues3D[BaseFunct][QuadFormula][MultiIndex]
        = Values;
    }

    /** return requested FE function values or derivatives */
    static double **GetRefElementValues
        (BaseFunct3D BaseFunct, QuadFormula3D QuadFormula,
         MultiIndex3D MultiIndex)
    {
      return RefElementValues3D[BaseFunct][QuadFormula][MultiIndex];
    }

    /** register joint values of FE base function on ref element */
    static void RegisterJointValues3D
        (BaseFunct3D BaseFunct, QuadFormula2D formula, int joint, 
         double **Values)
    {
      JointValues3D[BaseFunct][formula][joint] = Values;
    }

    /** return requested joint values of FE base functions */
    static double **GetJointValues3D
        (BaseFunct3D BaseFunct, QuadFormula2D formula, int joint)
    {
      return JointValues3D[BaseFunct][formula][joint];
    }

    /** register joint derivatives of FE base function on ref element */
    static void RegisterJointDerivatives3D
        (BaseFunct3D BaseFunct, QuadFormula2D formula, int joint, 
         MultiIndex3D MultiIndex, double **Values)
    {
      JointDerivatives3D[BaseFunct][formula][joint][MultiIndex] = Values;
    }

    /** return requested joint values of FE base functions */
    static double **GetJointDerivatives3D
        (BaseFunct3D BaseFunct, QuadFormula2D formula, int joint,
        MultiIndex3D MultiIndex)
    {
      return JointDerivatives3D[BaseFunct][formula][joint][MultiIndex];
    }

    /** register FE function values or derivatives on current element */
    static void RegisterOrigElementValues
        (BaseFunct3D BaseFunct, MultiIndex3D MultiIndex, double **Values)
    {
      OrigElementValues3D[BaseFunct][MultiIndex] = Values;
    }

    /** return requested FE function values or derivatives */
    static double **GetOrigElementValues
        (BaseFunct3D BaseFunct, MultiIndex3D MultiIndex)
    {
      return OrigElementValues3D[BaseFunct][MultiIndex];
    }

    /** calculate the values of base functions and their derivatives
        on the original element */
    static RefTrans3D GetOrigValues(TBaseCell *cell, TFE3D *element, 
                                    int N_Points, 
                                    double *xi, double *eta, double *zeta,
                                    int N_Functs,
                                    BaseFunct3D BaseFunct,
                                    QuadFormula3D QuadFormula);

    /** calculate functions and derivatives from reference element
        to original element */
    static void GetOrigValues(RefTrans3D RefTrans,
                double xi, double eta, double zeta,
                TBaseFunct3D *bf, TCollection *Coll, TBaseCell *cell,
                double *uref, double *uxiref, double *uetaref, double *uzetaref,
                double *uorig, double *uxorig, double *uyorig, double *uzorig);

};

#endif
