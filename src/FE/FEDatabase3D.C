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
// Purpose:     store all used information for a FEM in 3D
//
// Author:      Gunar Matthies  26.11.99
//
// History:     start reimplementation 26.11.99 (GM)
//
// =======================================================================

#include <FEDatabase3D.h>
#include <RefTrans3D.h>
#include <AllRefTrans3D.h>
#include <Database.h>
#include <string.h>

#include <BoundFace.h>
#include <IsoInterfaceJoint3D.h>

#include <AllFEDescs3D.h>
#include <AllNodalFunctionals3D.h>
#include <AllBaseFunctions3D.h>
#include <AllFE3DMappers.h>
#include <AllHNDescs3D.h>

#include <LinAlg.h>
#include <MooNMD_Io.h>

#include <stdlib.h>

// =======================================================================
// initialize static members
// =======================================================================
TQuadFormula1D *TFEDatabase3D::QuadFormulas1D[N_QuadFormulas_1D] =  { NULL };
TQuadFormula2D *TFEDatabase3D::QuadFormulas2D[N_QuadFormulas_2D] =  { NULL };
TQuadFormula3D *TFEDatabase3D::QuadFormulas3D[N_QuadFormulas_3D] = { NULL };

QuadFormula1D TFEDatabase3D::QFLineFromDegree[MAXDEGREE] = { Gauss1Line };
int TFEDatabase3D::HighestAccuracyLine = 0;

QuadFormula2D TFEDatabase3D::QFTriaFromDegree[MAXDEGREE] = { BaryCenterTria };
QuadFormula2D TFEDatabase3D::QFQuadFromDegree[MAXDEGREE] = { VertexQuad };

TFE3D *TFEDatabase3D::FEs3D[N_FEs3D] = { NULL };
TFEDesc3D *TFEDatabase3D::FEDescs3D[N_FEDescs3D] = { NULL };
TBaseFunct3D *TFEDatabase3D::BaseFuncts3D[N_BaseFuncts3D] = { NULL };
TNodalFunctional3D *TFEDatabase3D::NodalFunctionals3D[N_NodalFunctionals3D] 
   = { NULL };

TFE3DMapper *TFEDatabase3D::FE3DMapper[N_FEDescs3D][N_FEDescs3D] = { NULL };
TFE3DMapper1Reg *TFEDatabase3D::FE3DMapper1Reg[N_FEDescs3D][N_FEDescs3D]
   = { NULL };

THNDesc *TFEDatabase3D::HNDescs3D[N_HNDescs] = { NULL };

TRefTrans3D *TFEDatabase3D::ReferenceTrans3D[N_RefTrans3D] = { NULL };

double **TFEDatabase3D::RefElementValues3D[N_BaseFuncts3D][N_QuadFormulas_3D]
        [N_MultiIndices3D] = { NULL };
double **TFEDatabase3D::OrigElementValues3D[N_BaseFuncts3D]
        [N_MultiIndices3D] = { NULL };
double **TFEDatabase3D::JointValues3D[N_BaseFuncts3D][N_QuadFormulas_2D]
        [MAXN_JOINTS] = { NULL };
double **TFEDatabase3D::JointDerivatives3D[N_BaseFuncts3D][N_QuadFormulas_2D]
        [MAXN_JOINTS][N_MultiIndices3D] = { NULL };

double *TFEDatabase3D::ProlongationMatrix3D[MaxN_BaseFunctions3D]
        [N_REFDESC][MaxN_BaseFunctions3D][MAXN_CHILDREN] = { NULL };

double *TFEDatabase3D::RestrictionMatrix3D[MaxN_BaseFunctions3D]
        [N_REFDESC][MaxN_BaseFunctions3D][MAXN_CHILDREN] = { NULL };

/** Id of FEDesc from FE Id */ 
FEDesc3D TFEDatabase3D::FEDesc3D_IDFromFE3D[N_FEs3D] = { FE_C_T_P1_3D };

/** TFEDesc3D object from FE Id */ 
TFEDesc3D *TFEDatabase3D::FEDesc3DFromFE3D[N_FEs3D] = { NULL };

/** Id of BaseFunct3D from FE Id */ 
BaseFunct3D TFEDatabase3D::BaseFunct3D_IDFromFE3D[N_FEs3D] = { BF_C_T_P1_3D };

/** number of base functions from FE Id */ 
int TFEDatabase3D::N_BaseFunctFromFE3D[N_FEs3D] = { 0 };

/** polynomial degree of base functions from FE Id */ 
int TFEDatabase3D::PolynomialDegreeFromFE3D[N_FEs3D] = { 0 };

/** accuracy of base functions from FE Id */ 
int TFEDatabase3D::AccuracyFromFE3D[N_FEs3D] = { 0 };

/** TBaseFunct3DFEDesc object from FE Id */ 
TBaseFunct3D *TFEDatabase3D::BaseFunct3DFromFE3D[N_FEs3D] = { NULL };

/** Id of NodalFunctional3D from FE Id */ 
NodalFunctional3D TFEDatabase3D::NodalFunctional3D_IDFromFE3D[N_FEs3D] = {
NF_C_T_P1_3D };

/** TNodalFunctional3D object from FE Id */ 
TNodalFunctional3D *TFEDatabase3D::NodalFunctional3DFromFE3D[N_FEs3D] = { NULL };

/** Id of RefTrans3D from FE Id */ 
RefTrans3D TFEDatabase3D::RefTrans3D_IDFromFE3D[N_FEs3D] = { TetraAffin };

/** reference element from FE Id */
BF3DRefElements TFEDatabase3D::RefElementFromFE3D[N_FEs3D] = { BFUnitTetrahedron };

/** get tetrahedron quadrature formula for given acuracy */
QuadFormula3D TFEDatabase3D::QFTetraFromDegree[MAXDEGREE] = { BaryCenterTetra };

/** get hexahedron quadrature formula for given acuracy */
QuadFormula3D TFEDatabase3D::QFHexaFromDegree[MAXDEGREE] = { VertexHexa };

/** get hexahedron quadrature formula for convolution */
QuadFormula3D TFEDatabase3D::QFConvolutionHexaFromDegree[MAXDEGREE] = { VerticesAndOrigin };

TFEDatabase3D::TFEDatabase3D()
{
  RegisterAllQuadFormulas();
  RegisterAllFEDescs();
  RegisterAllBaseFunctions();
  RegisterAllNodalFunctionals();
  RegisterAllFEs();
  RegisterAllHangingNodes();
  RegisterAllFEMappers();
  RegisterAllRefTrans();

  GenerateArrays();
}

void TFEDatabase3D::RegisterAllQuadFormulas()
{
  TQuadFormula1D *qf1d;
  TQuadFormulaTetra *qftetra;
  TQuadFormulaHexa *qfhexa;
  TQuadFormulaTria *qftria;
  TQuadFormulaQuad *qfquad;
  
  // ===================================================================
  // register 1d quadrature formulas
  // ===================================================================
  qf1d = new TQuadFormula1D();
  qf1d->Gauss1();
  RegisterQuadFormula1D(Gauss1Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss2();
  RegisterQuadFormula1D(Gauss2Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss3();
  RegisterQuadFormula1D(Gauss3Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss4();
  RegisterQuadFormula1D(Gauss4Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss5();
  RegisterQuadFormula1D(Gauss5Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss6();
  RegisterQuadFormula1D(Gauss6Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss7();
  RegisterQuadFormula1D(Gauss7Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss8();
  RegisterQuadFormula1D(Gauss8Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss9();
  RegisterQuadFormula1D(Gauss9Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss10();
  RegisterQuadFormula1D(Gauss10Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss11();
  RegisterQuadFormula1D(Gauss11Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss12();
  RegisterQuadFormula1D(Gauss12Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss2W1();
  RegisterQuadFormula1D(Gauss2W1Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss4W1();
  RegisterQuadFormula1D(Gauss4W1Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss6W1();
  RegisterQuadFormula1D(Gauss6W1Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss8W1();
  RegisterQuadFormula1D(Gauss8W1Line, qf1d);

  qf1d = new TQuadFormula1D();
  qf1d->Gauss16W2();
  RegisterQuadFormula1D(Gauss16W2Line, qf1d);

  // =====================================================================
  // register triangle quadrature formulas
  // =====================================================================
  qftria = new TQuadFormulaTria();
  qftria->BaryCenter();
  RegisterQuadFormula2D(BaryCenterTria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->MidPoint();
  RegisterQuadFormula2D(MidPointTria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->SevenPoint();
  RegisterQuadFormula2D(SevenPointTria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->Gauss3();
  RegisterQuadFormula2D(Gauss3Tria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->Vertex();
  RegisterQuadFormula2D(VertexTria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->Degree8();
  RegisterQuadFormula2D(Degree8Tria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->Degree11();
  RegisterQuadFormula2D(Degree11Tria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->Degree19();
  RegisterQuadFormula2D(Degree19Tria, qftria);

  // =====================================================================
  // register quadrangle quadrature formulas
  // =====================================================================
  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss2();
  RegisterQuadFormula2D(Gauss2Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss3();
  RegisterQuadFormula2D(Gauss3Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss4();
  RegisterQuadFormula2D(Gauss4Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss5();
  RegisterQuadFormula2D(Gauss5Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss6();
  RegisterQuadFormula2D(Gauss6Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss7();
  RegisterQuadFormula2D(Gauss7Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss8();
  RegisterQuadFormula2D(Gauss8Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Gauss9();
  RegisterQuadFormula2D(Gauss9Quad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Vertex();
  RegisterQuadFormula2D(VertexQuad, qfquad);

  qfquad = new TQuadFormulaQuad();
  qfquad->Simpson();
  RegisterQuadFormula2D(SimpsonQuad, qfquad);

  // =====================================================================
  // register hexahedron quadrature formulas
  // =====================================================================
  qfhexa = new TQuadFormulaHexa();
  qfhexa->Vertex();
  RegisterQuadFormula3D(VertexHexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss2();
  RegisterQuadFormula3D(Gauss2Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss3();
  RegisterQuadFormula3D(Gauss3Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss4();
  RegisterQuadFormula3D(Gauss4Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss5();
  RegisterQuadFormula3D(Gauss5Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss6();
  RegisterQuadFormula3D(Gauss6Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss7();
  RegisterQuadFormula3D(Gauss7Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss8();
  RegisterQuadFormula3D(Gauss8Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Gauss9();
  RegisterQuadFormula3D(Gauss9Hexa, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->VerticesAndOrigin();
  RegisterQuadFormula3D(VerticesAndOrigin, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->VerticesAndOrigin15();
  RegisterQuadFormula3D(VerticesAndOrigin15, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->VerticesAndOrigin57();
  RegisterQuadFormula3D(VerticesAndOrigin57, qfhexa);

  qfhexa = new TQuadFormulaHexa();
  qfhexa->Degree7_Points38();
  RegisterQuadFormula3D(Degree7_Points38, qfhexa);
  // =====================================================================
  // register tetrahedron quadrature formulas
  // =====================================================================
  qftetra = new TQuadFormulaTetra();
  qftetra->BaryCenter();
  RegisterQuadFormula3D(BaryCenterTetra, qftetra);

  qftetra = new TQuadFormulaTetra();
  qftetra->Vertex();
  RegisterQuadFormula3D(VertexTetra, qftetra);

  qftetra = new TQuadFormulaTetra();
  qftetra->P2Exact();
  RegisterQuadFormula3D(P2Tetra, qftetra);

  qftetra = new TQuadFormulaTetra();
  qftetra->P4Exact();
  RegisterQuadFormula3D(P4Tetra, qftetra);

  qftetra = new TQuadFormulaTetra();
  qftetra->P5Exact();
  RegisterQuadFormula3D(P5Tetra, qftetra);

  qftetra = new TQuadFormulaTetra();
  qftetra->P8Exact();
  RegisterQuadFormula3D(P8Tetra, qftetra);
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "quadrature formulas registered" << endl;
}

void TFEDatabase3D::RegisterAllFEDescs()
{
  RegisterFEDesc3D(FE_C_T_P00_3D, FE_C_T_P00_3D_Obj);
  RegisterFEDesc3D(FE_C_T_P0_3D, FE_C_T_P0_3D_Obj);
  RegisterFEDesc3D(FE_C_T_P1_3D, FE_C_T_P1_3D_Obj);
  RegisterFEDesc3D(FE_C_T_P2_3D, FE_C_T_P2_3D_Obj);
  RegisterFEDesc3D(FE_C_T_P3_3D, FE_C_T_P3_3D_Obj);

  RegisterFEDesc3D(FE_C_T_B2_3D, FE_C_T_B2_3D_Obj);

  RegisterFEDesc3D(FE_C_H_Q00_3D, FE_C_H_Q00_3D_Obj);
  RegisterFEDesc3D(FE_C_H_Q0_3D, FE_C_H_Q0_3D_Obj);
  RegisterFEDesc3D(FE_C_H_Q1_3D, FE_C_H_Q1_3D_Obj);
  RegisterFEDesc3D(FE_C_H_Q2_3D, FE_C_H_Q2_3D_Obj);
  RegisterFEDesc3D(FE_C_H_Q3_3D, FE_C_H_Q3_3D_Obj);
  RegisterFEDesc3D(FE_C_H_Q4_3D, FE_C_H_Q4_3D_Obj);

  RegisterFEDesc3D(FE_N_H_Q1_3D, FE_N_H_Q1_3D_Obj);
  RegisterFEDesc3D(FE_N_T_P1_3D, FE_N_T_P1_3D_Obj);

  RegisterFEDesc3D(FE_D_T_P1_3D, FE_D_T_P1_3D_Obj);
  RegisterFEDesc3D(FE_D_T_P2_3D, FE_D_T_P2_3D_Obj);
  RegisterFEDesc3D(FE_D_T_P3_3D, FE_D_T_P3_3D_Obj);

  RegisterFEDesc3D(FE_D_H_P1_3D, FE_D_H_P1_3D_Obj);
  RegisterFEDesc3D(FE_D_H_P2_3D, FE_D_H_P2_3D_Obj);
  RegisterFEDesc3D(FE_D_H_P3_3D, FE_D_H_P3_3D_Obj);
  
  RegisterFEDesc3D(FE_D_H_Q1_3D, FE_D_H_Q1_3D_Obj);
  RegisterFEDesc3D(FE_D_H_Q2_3D, FE_D_H_Q2_3D_Obj);

  RegisterFEDesc3D(FE_B_H_IB2_3D, FE_B_H_IB2_3D_Obj);

  RegisterFEDesc3D(FE_N_T_P2_3D, FE_N_T_P2_3D_Obj);

  RegisterFEDesc3D(FE_N_H_Q2_3D, FE_N_H_Q2_3D_Obj);
  RegisterFEDesc3D(FE_N_H_Q3_3D, FE_N_H_Q3_3D_Obj);
  RegisterFEDesc3D(FE_N_H_Q4_3D, FE_N_H_Q4_3D_Obj);
  
  RegisterFEDesc3D(FE_C_H_UL1_3D, FE_C_H_UL1_3D_Obj);  
  RegisterFEDesc3D(FE_C_H_UL2_3D, FE_C_H_UL2_3D_Obj);
  RegisterFEDesc3D(FE_C_H_UL3_3D, FE_C_H_UL3_3D_Obj);
  
  RegisterFEDesc3D(FE_N_T_RT0_3D, FE_N_T_RT0_3D_Obj);
  RegisterFEDesc3D(FE_N_T_RT1_3D, FE_N_T_RT1_3D_Obj);
  RegisterFEDesc3D(FE_N_T_RT2_3D, FE_N_T_RT2_3D_Obj);
  RegisterFEDesc3D(FE_N_T_RT3_3D, FE_N_T_RT3_3D_Obj);

  RegisterFEDesc3D(FE_N_T_BDDF1_3D, FE_N_T_BDDF1_3D_Obj);
  RegisterFEDesc3D(FE_N_T_BDDF2_3D, FE_N_T_BDDF2_3D_Obj);
  RegisterFEDesc3D(FE_N_T_BDDF3_3D, FE_N_T_BDDF3_3D_Obj);

  RegisterFEDesc3D(FE_N_H_RT0_3D, FE_N_H_RT0_3D_Obj);
  RegisterFEDesc3D(FE_N_H_RT1_3D, FE_N_H_RT1_3D_Obj);
  RegisterFEDesc3D(FE_N_H_RT2_3D, FE_N_H_RT2_3D_Obj);
  
  RegisterFEDesc3D(FE_N_H_BDDF1_3D, FE_N_H_BDDF1_3D_Obj);
  RegisterFEDesc3D(FE_N_H_BDDF2_3D, FE_N_H_BDDF2_3D_Obj);
  RegisterFEDesc3D(FE_N_H_BDDF3_3D, FE_N_H_BDDF3_3D_Obj);
  
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "FE descriptors registered" << endl;  
}

void TFEDatabase3D::RegisterAllBaseFunctions()
{
  RegisterBaseFunct3D(BF_C_T_P00_3D, BF_C_T_P00_3D_Obj);
  RegisterBaseFunct3D(BF_C_T_P0_3D, BF_C_T_P0_3D_Obj);
  RegisterBaseFunct3D(BF_C_T_P1_3D, BF_C_T_P1_3D_Obj);
  RegisterBaseFunct3D(BF_C_T_P2_3D, BF_C_T_P2_3D_Obj);
  RegisterBaseFunct3D(BF_C_T_P3_3D, BF_C_T_P3_3D_Obj);

  RegisterBaseFunct3D(BF_C_T_B2_3D, BF_C_T_B2_3D_Obj);

  RegisterBaseFunct3D(BF_C_H_Q00_3D, BF_C_H_Q00_3D_Obj);
  RegisterBaseFunct3D(BF_C_H_Q0_3D, BF_C_H_Q0_3D_Obj);
  RegisterBaseFunct3D(BF_C_H_Q1_3D, BF_C_H_Q1_3D_Obj);
  RegisterBaseFunct3D(BF_C_H_Q2_3D, BF_C_H_Q2_3D_Obj);
  RegisterBaseFunct3D(BF_C_H_Q3_3D, BF_C_H_Q3_3D_Obj);
  RegisterBaseFunct3D(BF_C_H_Q4_3D, BF_C_H_Q4_3D_Obj);

  RegisterBaseFunct3D(BF_N_H_Q1_3D, BF_N_H_Q1_3D_Obj);
  RegisterBaseFunct3D(BF_N_T_P1_3D, BF_N_T_P1_3D_Obj);

  RegisterBaseFunct3D(BF_D_T_P1_3D, BF_D_T_P1_3D_Obj);
  RegisterBaseFunct3D(BF_D_T_P2_3D, BF_D_T_P2_3D_Obj);
  RegisterBaseFunct3D(BF_D_T_P3_3D, BF_D_T_P3_3D_Obj);

  RegisterBaseFunct3D(BF_D_H_P1_3D, BF_D_H_P1_3D_Obj);
  RegisterBaseFunct3D(BF_D_H_P2_3D, BF_D_H_P2_3D_Obj);
  RegisterBaseFunct3D(BF_D_H_P3_3D, BF_D_H_P3_3D_Obj);
  
  RegisterBaseFunct3D(BF_D_H_Q1_3D, BF_D_H_Q1_3D_Obj);
  RegisterBaseFunct3D(BF_D_H_Q2_3D, BF_D_H_Q2_3D_Obj);

  RegisterBaseFunct3D(BF_B_H_IB2_3D, BF_B_H_IB2_3D_Obj);

  RegisterBaseFunct3D(BF_N_T_P2_3D, BF_N_T_P2_3D_Obj);

  RegisterBaseFunct3D(BF_N_H_Q2_3D, BF_N_H_Q2_3D_Obj);
  RegisterBaseFunct3D(BF_N_H_Q3_3D, BF_N_H_Q3_3D_Obj);
  RegisterBaseFunct3D(BF_N_H_Q4_3D, BF_N_H_Q4_3D_Obj);
  
  RegisterBaseFunct3D(BF_C_H_UL1_3D, BF_C_H_UL1_3D_Obj);
  RegisterBaseFunct3D(BF_C_H_UL2_3D, BF_C_H_UL2_3D_Obj);
  RegisterBaseFunct3D(BF_C_H_UL3_3D, BF_C_H_UL3_3D_Obj);
  
  RegisterBaseFunct3D(BF_N_T_RT0_3D, BF_N_T_RT0_3D_Obj);
  RegisterBaseFunct3D(BF_N_T_RT1_3D, BF_N_T_RT1_3D_Obj);
  RegisterBaseFunct3D(BF_N_T_RT2_3D, BF_N_T_RT2_3D_Obj);
  RegisterBaseFunct3D(BF_N_T_RT3_3D, BF_N_T_RT3_3D_Obj);

  RegisterBaseFunct3D(BF_N_T_BDDF1_3D, BF_N_T_BDDF1_3D_Obj);
  RegisterBaseFunct3D(BF_N_T_BDDF2_3D, BF_N_T_BDDF2_3D_Obj);
  RegisterBaseFunct3D(BF_N_T_BDDF3_3D, BF_N_T_BDDF3_3D_Obj);

  RegisterBaseFunct3D(BF_N_H_RT0_3D, BF_N_H_RT0_3D_Obj);
  RegisterBaseFunct3D(BF_N_H_RT1_3D, BF_N_H_RT1_3D_Obj);
  RegisterBaseFunct3D(BF_N_H_RT2_3D, BF_N_H_RT2_3D_Obj);

  RegisterBaseFunct3D(BF_N_H_BDDF1_3D, BF_N_H_BDDF1_3D_Obj);
  RegisterBaseFunct3D(BF_N_H_BDDF2_3D, BF_N_H_BDDF2_3D_Obj);
  RegisterBaseFunct3D(BF_N_H_BDDF3_3D, BF_N_H_BDDF3_3D_Obj);

#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "base functions registered" << endl;  
}

void TFEDatabase3D::RegisterAllNodalFunctionals()
{
  RegisterNodalFunctional3D(NF_C_T_P00_3D, NF_C_T_P00_3D_Obj);
  RegisterNodalFunctional3D(NF_C_T_P0_3D, NF_C_T_P0_3D_Obj);
  RegisterNodalFunctional3D(NF_C_T_P1_3D, NF_C_T_P1_3D_Obj);
  RegisterNodalFunctional3D(NF_C_T_P2_3D, NF_C_T_P2_3D_Obj);
  RegisterNodalFunctional3D(NF_C_T_P3_3D, NF_C_T_P3_3D_Obj);

  RegisterNodalFunctional3D(NF_C_T_B2_3D, NF_C_T_B2_3D_Obj);

  RegisterNodalFunctional3D(NF_C_H_Q00_3D, NF_C_H_Q00_3D_Obj);
  RegisterNodalFunctional3D(NF_C_H_Q0_3D, NF_C_H_Q0_3D_Obj);
  RegisterNodalFunctional3D(NF_C_H_Q1_3D, NF_C_H_Q1_3D_Obj);
  RegisterNodalFunctional3D(NF_C_H_Q2_3D, NF_C_H_Q2_3D_Obj);
  RegisterNodalFunctional3D(NF_C_H_Q3_3D, NF_C_H_Q3_3D_Obj);
  RegisterNodalFunctional3D(NF_C_H_Q4_3D, NF_C_H_Q4_3D_Obj);

  RegisterNodalFunctional3D(NF_N_H_Q1_3D, NF_N_H_Q1_3D_Obj);
  RegisterNodalFunctional3D(NF_N_T_P1_3D, NF_N_T_P1_3D_Obj);

  RegisterNodalFunctional3D(NF_D_T_P1_3D, NF_D_T_P1_3D_Obj);
  RegisterNodalFunctional3D(NF_D_T_P2_3D, NF_D_T_P2_3D_Obj);
  RegisterNodalFunctional3D(NF_D_T_P3_3D, NF_D_T_P3_3D_Obj);

  RegisterNodalFunctional3D(NF_D_H_P1_3D, NF_D_H_P1_3D_Obj);
  RegisterNodalFunctional3D(NF_D_H_P2_3D, NF_D_H_P2_3D_Obj);
  RegisterNodalFunctional3D(NF_D_H_P3_3D, NF_D_H_P3_3D_Obj);
  
  RegisterNodalFunctional3D(NF_D_H_Q1_3D, NF_D_H_Q1_3D_Obj);
  RegisterNodalFunctional3D(NF_D_H_Q2_3D, NF_D_H_Q2_3D_Obj);

  RegisterNodalFunctional3D(NF_B_H_IB2_3D, NF_B_H_IB2_3D_Obj);

  RegisterNodalFunctional3D(NF_S_H_Q2_3D, NF_S_H_Q2_3D_Obj);

  RegisterNodalFunctional3D(NF_N_T_P2_3D, NF_N_T_P2_3D_Obj);

  RegisterNodalFunctional3D(NF_N_H_Q2_3D, NF_N_H_Q2_3D_Obj);
  RegisterNodalFunctional3D(NF_N_H_Q3_3D, NF_N_H_Q3_3D_Obj);
  RegisterNodalFunctional3D(NF_N_H_Q4_3D, NF_N_H_Q4_3D_Obj);
  
  RegisterNodalFunctional3D(NF_C_H_UL1_3D, NF_C_H_UL1_3D_Obj);
  RegisterNodalFunctional3D(NF_C_H_UL2_3D, NF_C_H_UL2_3D_Obj);
  RegisterNodalFunctional3D(NF_C_H_UL3_3D, NF_C_H_UL3_3D_Obj);
  
  RegisterNodalFunctional3D(NF_N_T_RT0_3D, NF_N_T_RT0_3D_Obj);
  RegisterNodalFunctional3D(NF_N_T_RT1_3D, NF_N_T_RT1_3D_Obj);
  RegisterNodalFunctional3D(NF_N_T_RT2_3D, NF_N_T_RT2_3D_Obj);
  RegisterNodalFunctional3D(NF_N_T_RT3_3D, NF_N_T_RT3_3D_Obj);

  RegisterNodalFunctional3D(NF_N_T_BDDF1_3D, NF_N_T_BDDF1_3D_Obj);
  RegisterNodalFunctional3D(NF_N_T_BDDF2_3D, NF_N_T_BDDF2_3D_Obj);
  RegisterNodalFunctional3D(NF_N_T_BDDF3_3D, NF_N_T_BDDF3_3D_Obj);

  RegisterNodalFunctional3D(NF_N_H_RT0_3D, NF_N_H_RT0_3D_Obj);
  RegisterNodalFunctional3D(NF_N_H_RT1_3D, NF_N_H_RT1_3D_Obj);
  RegisterNodalFunctional3D(NF_N_H_RT2_3D, NF_N_H_RT2_3D_Obj);

  RegisterNodalFunctional3D(NF_N_H_BDDF1_3D, NF_N_H_BDDF1_3D_Obj);
  RegisterNodalFunctional3D(NF_N_H_BDDF2_3D, NF_N_H_BDDF2_3D_Obj);
  RegisterNodalFunctional3D(NF_N_H_BDDF3_3D, NF_N_H_BDDF3_3D_Obj);

#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "nodal functionals registered" << endl;
}

void TFEDatabase3D::RegisterAllFEs()
{
  TFE3D *ele3D;

  // ======================================================================
  // elements on tetrahedrons
  // ======================================================================
  ele3D = new TFE3D(BF_C_T_P00_3D, NF_C_T_P00_3D, TetraAffin, FE_C_T_P00_3D, 0);
  RegisterFE3D(C_P00_3D_T_A, ele3D);
  ele3D = new TFE3D(BF_C_T_P0_3D, NF_C_T_P0_3D, TetraAffin, FE_C_T_P0_3D, 0);
  RegisterFE3D(C_P0_3D_T_A, ele3D);
  ele3D = new TFE3D(BF_C_T_P1_3D, NF_C_T_P1_3D, TetraAffin, FE_C_T_P1_3D, 0);
  RegisterFE3D(C_P1_3D_T_A, ele3D);
  ele3D = new TFE3D(BF_C_T_P2_3D, NF_C_T_P2_3D, TetraAffin, FE_C_T_P2_3D, 0);
  RegisterFE3D(C_P2_3D_T_A, ele3D);
  ele3D = new TFE3D(BF_C_T_P3_3D, NF_C_T_P3_3D, TetraAffin, FE_C_T_P3_3D, 0);
  RegisterFE3D(C_P3_3D_T_A, ele3D);

  ele3D = new TFE3D(BF_C_T_B2_3D, NF_C_T_B2_3D, TetraAffin, FE_C_T_B2_3D, 0);
  RegisterFE3D(C_B2_3D_T_A, ele3D);

  ele3D = new TFE3D(BF_N_T_P1_3D, NF_N_T_P1_3D, TetraAffin, FE_N_T_P1_3D, 0);
  RegisterFE3D(N_P1_3D_T_A, ele3D);

  ele3D = new TFE3D(BF_N_T_P2_3D, NF_N_T_P2_3D, TetraAffin, FE_N_T_P2_3D, 0);
  RegisterFE3D(N_P2_3D_T_A, ele3D);

  ele3D = new TFE3D(BF_D_T_P1_3D, NF_D_T_P1_3D, TetraAffin, FE_D_T_P1_3D, 0);
  RegisterFE3D(D_P1_3D_T_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_D_T_P2_3D, NF_D_T_P2_3D, TetraAffin, FE_D_T_P2_3D, 0);
  RegisterFE3D(D_P2_3D_T_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_D_T_P3_3D, NF_D_T_P3_3D, TetraAffin, FE_D_T_P3_3D, 0);
  RegisterFE3D(D_P3_3D_T_A, ele3D);
  //ele3D->CheckNFandBF();


  // ======================================================================
  // elements on arbitrary hexahedrons
  // ======================================================================
  ele3D = new TFE3D(BF_C_H_Q00_3D, NF_C_H_Q00_3D, HexaTrilinear, FE_C_H_Q00_3D, 0);
  RegisterFE3D(C_Q00_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_C_H_Q0_3D, NF_C_H_Q0_3D, HexaTrilinear, FE_C_H_Q0_3D, 0);
  RegisterFE3D(C_Q0_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_C_H_Q1_3D, NF_C_H_Q1_3D, HexaTrilinear, FE_C_H_Q1_3D, 0);
  RegisterFE3D(C_Q1_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_C_H_Q2_3D, NF_C_H_Q2_3D, HexaTrilinear, FE_C_H_Q2_3D, 0);
  RegisterFE3D(C_Q2_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_C_H_Q3_3D, NF_C_H_Q3_3D, HexaTrilinear, FE_C_H_Q3_3D, 0);
  RegisterFE3D(C_Q3_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_C_H_Q4_3D, NF_C_H_Q4_3D, HexaTrilinear, FE_C_H_Q4_3D, 0);
  RegisterFE3D(C_Q4_3D_H_M, ele3D);
  
  ele3D = new TFE3D(BF_N_H_Q1_3D, NF_N_H_Q1_3D, HexaTrilinear, FE_N_H_Q1_3D, 0);
  RegisterFE3D(N_Q1_3D_H_M, ele3D);
  
  ele3D = new TFE3D(BF_D_H_P1_3D, NF_D_H_P1_3D, HexaTrilinear, FE_D_H_P1_3D, 0);
  RegisterFE3D(D_P1_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_D_H_P2_3D, NF_D_H_P2_3D, HexaTrilinear, FE_D_H_P2_3D, 0);
  RegisterFE3D(D_P2_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_D_H_P3_3D, NF_D_H_P3_3D, HexaTrilinear, FE_D_H_P3_3D, 0);
  RegisterFE3D(D_P3_3D_H_M, ele3D);
  
  ele3D = new TFE3D(BF_D_H_Q1_3D, NF_D_H_Q1_3D, HexaTrilinear, FE_D_H_Q1_3D, 0);
  RegisterFE3D(D_Q1_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_D_H_Q2_3D, NF_D_H_Q2_3D, HexaTrilinear, FE_D_H_Q2_3D, 0);
  RegisterFE3D(D_Q2_3D_H_M, ele3D);

  ele3D = new TFE3D(BF_B_H_IB2_3D, NF_B_H_IB2_3D, HexaTrilinear, FE_B_H_IB2_3D, 0);
  RegisterFE3D(B_IB2_3D_H_M, ele3D);

  ele3D = new TFE3D(BF_N_H_Q2_3D, NF_N_H_Q2_3D, HexaTrilinear, FE_N_H_Q2_3D, 0);
  RegisterFE3D(N_Q2_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_N_H_Q3_3D, NF_N_H_Q3_3D, HexaTrilinear, FE_N_H_Q3_3D, 0);
  RegisterFE3D(N_Q3_3D_H_M, ele3D);
  ele3D = new TFE3D(BF_N_H_Q4_3D, NF_N_H_Q4_3D, HexaTrilinear, FE_N_H_Q4_3D, 0);
  RegisterFE3D(N_Q4_3D_H_M, ele3D);
  ele3D->CheckNFandBF();
  
  //========LOCALPROJECTION==============
  ele3D = new TFE3D(BF_C_H_UL1_3D, NF_C_H_UL1_3D, HexaTrilinear, FE_C_H_UL1_3D, 0);
  RegisterFE3D(C_UL1_3D_H_M, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_C_H_UL2_3D, NF_C_H_UL2_3D, HexaTrilinear, FE_C_H_UL2_3D, 0);
  RegisterFE3D(C_UL2_3D_H_M, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_C_H_UL3_3D, NF_C_H_UL3_3D, HexaTrilinear, FE_C_H_UL3_3D, 0);
  RegisterFE3D(C_UL3_3D_H_M, ele3D);
  ele3D->CheckNFandBF();
  //=====================================
  
  // ======================================================================
  // elements on affine hexahedrons (bricks)
  // ======================================================================
  ele3D = new TFE3D(BF_C_H_Q00_3D, NF_C_H_Q00_3D, HexaAffin, FE_C_H_Q00_3D, 0);
  RegisterFE3D(C_Q00_3D_H_A, ele3D);
  ele3D = new TFE3D(BF_C_H_Q0_3D, NF_C_H_Q0_3D, HexaAffin, FE_C_H_Q0_3D, 0);
  RegisterFE3D(C_Q0_3D_H_A, ele3D);
  ele3D = new TFE3D(BF_C_H_Q1_3D, NF_C_H_Q1_3D, HexaAffin, FE_C_H_Q1_3D, 0);
  RegisterFE3D(C_Q1_3D_H_A, ele3D);
  ele3D = new TFE3D(BF_C_H_Q2_3D, NF_C_H_Q2_3D, HexaAffin, FE_C_H_Q2_3D, 0);
  RegisterFE3D(C_Q2_3D_H_A, ele3D);
  ele3D = new TFE3D(BF_C_H_Q3_3D, NF_C_H_Q3_3D, HexaAffin, FE_C_H_Q3_3D, 0);
  RegisterFE3D(C_Q3_3D_H_A, ele3D);
  ele3D = new TFE3D(BF_C_H_Q4_3D, NF_C_H_Q4_3D, HexaAffin, FE_C_H_Q4_3D, 0);
  RegisterFE3D(C_Q4_3D_H_A, ele3D);

  ele3D = new TFE3D(BF_N_H_Q1_3D, NF_N_H_Q1_3D, HexaAffin, FE_N_H_Q1_3D, 0);
  RegisterFE3D(N_Q1_3D_H_A, ele3D);

  ele3D = new TFE3D(BF_D_H_P1_3D, NF_D_H_P1_3D, HexaAffin, FE_D_H_P1_3D, 0);
  RegisterFE3D(D_P1_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_D_H_P2_3D, NF_D_H_P2_3D, HexaAffin, FE_D_H_P2_3D, 0);
  RegisterFE3D(D_P2_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_D_H_P3_3D, NF_D_H_P3_3D, HexaAffin, FE_D_H_P3_3D, 0);
  RegisterFE3D(D_P3_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  
  ele3D = new TFE3D(BF_D_H_Q1_3D, NF_D_H_Q1_3D, HexaAffin, FE_D_H_Q1_3D, 0);
  RegisterFE3D(D_Q1_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_D_H_Q2_3D, NF_D_H_Q2_3D, HexaAffin, FE_D_H_Q2_3D, 0);
  RegisterFE3D(D_Q2_3D_H_A, ele3D);
  ele3D->CheckNFandBF();

  ele3D = new TFE3D(BF_B_H_IB2_3D, NF_B_H_IB2_3D, HexaAffin, FE_B_H_IB2_3D, 0);
  RegisterFE3D(B_IB2_3D_H_A, ele3D);

  ele3D = new TFE3D(BF_N_H_Q2_3D, NF_N_H_Q2_3D, HexaAffin, FE_N_H_Q2_3D, 0);
  RegisterFE3D(N_Q2_3D_H_A, ele3D);
  ele3D = new TFE3D(BF_N_H_Q3_3D, NF_N_H_Q3_3D, HexaAffin, FE_N_H_Q3_3D, 0);
  RegisterFE3D(N_Q3_3D_H_A, ele3D);
  ele3D = new TFE3D(BF_N_H_Q4_3D, NF_N_H_Q4_3D, HexaAffin, FE_N_H_Q4_3D, 0);
  RegisterFE3D(N_Q4_3D_H_A, ele3D);
  ele3D->CheckNFandBF();

  ele3D = new TFE3D(BF_C_H_UL1_3D, NF_C_H_UL1_3D, HexaTrilinear, FE_C_H_UL1_3D, 0);
  RegisterFE3D(C_UL1_3D_H_M, ele3D);
  
  //========LOCALPROJECTION==============
  ele3D = new TFE3D(BF_C_H_UL1_3D, NF_C_H_UL1_3D, HexaAffin, FE_C_H_UL1_3D, 0);
  RegisterFE3D(C_UL1_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_C_H_UL2_3D, NF_C_H_UL2_3D, HexaAffin, FE_C_H_UL2_3D, 0);
  RegisterFE3D(C_UL2_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_C_H_UL3_3D, NF_C_H_UL3_3D, HexaAffin, FE_C_H_UL3_3D, 0);
  RegisterFE3D(C_UL3_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  //=====================================
  
  // for mixed problems
  ele3D = new TFE3D(BF_N_T_RT0_3D, NF_N_T_RT0_3D, TetraAffin, FE_N_T_RT0_3D, 0);
  RegisterFE3D(N_RT0_3D_T_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_T_RT1_3D, NF_N_T_RT1_3D, TetraAffin, FE_N_T_RT1_3D, 0);
  RegisterFE3D(N_RT1_3D_T_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_T_RT2_3D, NF_N_T_RT2_3D, TetraAffin, FE_N_T_RT2_3D, 0);
  RegisterFE3D(N_RT2_3D_T_A, ele3D);
  //ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_T_RT3_3D, NF_N_T_RT3_3D, TetraAffin, FE_N_T_RT3_3D, 0);
  RegisterFE3D(N_RT3_3D_T_A, ele3D);
  //ele3D->CheckNFandBF();

  ele3D = new TFE3D(BF_N_T_BDDF1_3D, NF_N_T_BDDF1_3D, TetraAffin, FE_N_T_BDDF1_3D, 0);
  RegisterFE3D(N_BDDF1_3D_T_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_T_BDDF2_3D, NF_N_T_BDDF2_3D, TetraAffin, FE_N_T_BDDF2_3D, 0);
  RegisterFE3D(N_BDDF2_3D_T_A, ele3D);
  //ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_T_BDDF3_3D, NF_N_T_BDDF3_3D, TetraAffin, FE_N_T_BDDF3_3D, 0);
  RegisterFE3D(N_BDDF3_3D_T_A, ele3D);
  //ele3D->CheckNFandBF();

  ele3D = new TFE3D(BF_N_H_RT0_3D, NF_N_H_RT0_3D, HexaAffin, FE_N_H_RT0_3D, 0);
  RegisterFE3D(N_RT0_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_H_RT1_3D, NF_N_H_RT1_3D, HexaAffin, FE_N_H_RT1_3D, 0);
  RegisterFE3D(N_RT1_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_H_RT2_3D, NF_N_H_RT2_3D, HexaAffin, FE_N_H_RT2_3D, 0);
  RegisterFE3D(N_RT2_3D_H_A, ele3D);
  
  ele3D = new TFE3D(BF_N_H_BDDF1_3D, NF_N_H_BDDF1_3D, HexaAffin, FE_N_H_BDDF1_3D, 0);
  RegisterFE3D(N_BDDF1_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_H_BDDF2_3D, NF_N_H_BDDF2_3D, HexaAffin, FE_N_H_BDDF2_3D, 0);
  RegisterFE3D(N_BDDF2_3D_H_A, ele3D);
  ele3D->CheckNFandBF();
  ele3D = new TFE3D(BF_N_H_BDDF3_3D, NF_N_H_BDDF3_3D, HexaAffin, FE_N_H_BDDF3_3D, 0);
  RegisterFE3D(N_BDDF3_3D_H_A, ele3D);
  ele3D->CheckNFandBF();

  
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "finite element registered" << endl;
}

void TFEDatabase3D::RegisterAllFEMappers()
{
  // ========================================================================
  // regular grid, same pattern on both sides
  // ========================================================================
  RegisterFE3DMapper(FE_C_T_P00_3D, FE_C_T_P00_3D, D_D);
  RegisterFE3DMapper(FE_C_T_P0_3D, FE_C_T_P0_3D, D_D);
  RegisterFE3DMapper(FE_C_T_P1_3D, FE_C_T_P1_3D, P1_P1);
  RegisterFE3DMapper(FE_C_T_P2_3D, FE_C_T_P2_3D, P2_P2);
  RegisterFE3DMapper(FE_C_T_P3_3D, FE_C_T_P3_3D, P3_P3);

  RegisterFE3DMapper(FE_C_T_B2_3D, FE_C_T_B2_3D, P2B_P2B);

  RegisterFE3DMapper(FE_C_H_Q00_3D, FE_C_H_Q00_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q0_3D, FE_C_H_Q0_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q1_3D, FE_C_H_Q1_3D, Q1_Q1);
  RegisterFE3DMapper(FE_C_H_Q2_3D, FE_C_H_Q2_3D, Q2_Q2);
  RegisterFE3DMapper(FE_C_H_Q3_3D, FE_C_H_Q3_3D, Q3_Q3);
  RegisterFE3DMapper(FE_C_H_Q4_3D, FE_C_H_Q4_3D, Q4_Q4);

  RegisterFE3DMapper(FE_N_H_Q1_3D, FE_N_H_Q1_3D, N1_N1);
  RegisterFE3DMapper(FE_N_T_P1_3D, FE_N_T_P1_3D, NP1_NP1);

  RegisterFE3DMapper(FE_D_T_P1_3D, FE_D_T_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_T_P2_3D, FE_D_T_P2_3D, D_D);
  RegisterFE3DMapper(FE_D_T_P3_3D, FE_D_T_P3_3D, D_D);

  RegisterFE3DMapper(FE_D_H_P1_3D, FE_D_H_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P2_3D, FE_D_H_P2_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P3_3D, FE_D_H_P3_3D, D_D);
  
  RegisterFE3DMapper(FE_D_H_Q1_3D, FE_D_H_Q1_3D, D_D);
  RegisterFE3DMapper(FE_D_H_Q2_3D, FE_D_H_Q2_3D, D_D);

  RegisterFE3DMapper(FE_B_H_IB2_3D, FE_B_H_IB2_3D, D_D);

  RegisterFE3DMapper(FE_N_T_P2_3D, FE_N_T_P2_3D, P1_P1);
  
  RegisterFE3DMapper(FE_N_H_Q2_3D, FE_N_H_Q2_3D, N2_N2);
  RegisterFE3DMapper(FE_N_H_Q3_3D, FE_N_H_Q3_3D, N3_N3);
  RegisterFE3DMapper(FE_N_H_Q4_3D, FE_N_H_Q4_3D, N4_N4);
  
  //========LOCALPROJECTION==============
  RegisterFE3DMapper(FE_C_H_UL1_3D, FE_C_H_Q1_3D, Q1_Q1);
  RegisterFE3DMapper(FE_C_H_UL2_3D, FE_C_H_Q2_3D, Q2_Q2);
  RegisterFE3DMapper(FE_C_H_Q1_3D, FE_C_H_UL1_3D, Q1_Q1);
  RegisterFE3DMapper(FE_C_H_Q2_3D, FE_C_H_UL2_3D, Q2_Q2);
  RegisterFE3DMapper(FE_C_H_UL1_3D, FE_C_H_UL1_3D, Q1_Q1);
  RegisterFE3DMapper(FE_C_H_UL2_3D, FE_C_H_UL2_3D, Q2_Q2);

  RegisterFE3DMapper(FE_C_H_UL3_3D, FE_C_H_Q3_3D, Q3_Q3);
  RegisterFE3DMapper(FE_C_H_Q3_3D, FE_C_H_UL3_3D, Q3_Q3);
  RegisterFE3DMapper(FE_C_H_UL3_3D, FE_C_H_UL3_3D, Q3_Q3);
  //=====================================

  // ========================================================================
  // discontinuous finite element spaces
  // ========================================================================
  RegisterFE3DMapper(FE_C_H_Q0_3D, FE_D_H_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P1_3D, FE_C_H_Q0_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q0_3D, FE_D_H_P2_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P2_3D, FE_C_H_Q0_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q0_3D, FE_D_H_P3_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P3_3D, FE_C_H_Q0_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P1_3D, FE_D_H_P2_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P2_3D, FE_D_H_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P1_3D, FE_D_H_P3_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P3_3D, FE_D_H_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P2_3D, FE_D_H_P3_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q0_3D, FE_D_H_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P3_3D, FE_D_H_P2_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q00_3D, FE_D_H_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P1_3D, FE_C_H_Q00_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q00_3D, FE_D_H_P2_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P2_3D, FE_C_H_Q00_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q00_3D, FE_D_H_P3_3D, D_D);
  RegisterFE3DMapper(FE_D_H_P3_3D, FE_C_H_Q00_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q00_3D, FE_C_H_Q0_3D, D_D);
  RegisterFE3DMapper(FE_C_H_Q0_3D, FE_C_H_Q00_3D, D_D);
  RegisterFE3DMapper(FE_C_T_P0_3D, FE_C_T_P00_3D, D_D);
  RegisterFE3DMapper(FE_C_T_P00_3D, FE_C_T_P0_3D, D_D);
  RegisterFE3DMapper(FE_C_T_P00_3D, FE_D_T_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_T_P1_3D, FE_C_T_P00_3D, D_D);
  RegisterFE3DMapper(FE_C_T_P0_3D, FE_D_T_P1_3D, D_D);
  RegisterFE3DMapper(FE_D_T_P1_3D, FE_C_T_P0_3D, D_D);
  
  RegisterFE3DMapper(FE_N_T_RT0_3D, FE_N_T_RT0_3D, NP1_NP1);
  RegisterFE3DMapper(FE_N_T_RT1_3D, FE_N_T_RT1_3D, NP2_NP2);
  RegisterFE3DMapper(FE_N_T_RT2_3D, FE_N_T_RT2_3D, P2_P2);
  RegisterFE3DMapper(FE_N_T_RT3_3D, FE_N_T_RT3_3D, P3_P3);
  RegisterFE3DMapper(FE_N_T_BDDF1_3D, FE_N_T_BDDF1_3D, NP2_NP2);
  RegisterFE3DMapper(FE_N_T_BDDF2_3D, FE_N_T_BDDF2_3D, P2_P2);
  RegisterFE3DMapper(FE_N_T_BDDF3_3D, FE_N_T_BDDF3_3D, P3_P3);

  RegisterFE3DMapper(FE_N_H_RT0_3D, FE_N_H_RT0_3D, N1_N1);
  RegisterFE3DMapper(FE_N_H_RT1_3D, FE_N_H_RT1_3D, Q1_Q1);
  RegisterFE3DMapper(FE_N_H_RT2_3D, FE_N_H_RT2_3D, Q2_Q2);
  
  RegisterFE3DMapper(FE_N_H_BDDF1_3D, FE_N_H_BDDF1_3D, N2_N2);
  RegisterFE3DMapper(FE_N_H_BDDF2_3D, FE_N_H_BDDF2_3D, N3_N3);
  RegisterFE3DMapper(FE_N_H_BDDF3_3D, FE_N_H_BDDF3_3D, N4_N4);

  // ========================================================================
  // 1-regular grid, same pattern on both sides
  // ========================================================================
  RegisterFE3DMapper1Reg(FE_C_T_P1_3D, FE_C_T_P1_3D, P1_P1_1Reg);
  RegisterFE3DMapper1Reg(FE_C_T_P2_3D, FE_C_T_P2_3D, P2_P2_1Reg);
  RegisterFE3DMapper1Reg(FE_C_T_P3_3D, FE_C_T_P3_3D, P3_P3_1Reg);
  RegisterFE3DMapper1Reg(FE_C_H_Q1_3D, FE_C_H_Q1_3D, Q1_Q1_1Reg);
  RegisterFE3DMapper1Reg(FE_C_H_Q2_3D, FE_C_H_Q2_3D, Q2_Q2_1Reg);
  RegisterFE3DMapper1Reg(FE_C_H_Q3_3D, FE_C_H_Q3_3D, Q3_Q3_1Reg);

  RegisterFE3DMapper1Reg(FE_N_T_P1_3D, FE_N_T_P1_3D, NP1_NP1_1Reg);
  RegisterFE3DMapper1Reg(FE_N_T_P2_3D, FE_N_T_P2_3D, NP2_NP2_1Reg);
  
  RegisterFE3DMapper(FE_C_H_UL1_3D, FE_C_H_UL1_3D, Q1_Q1);  
  RegisterFE3DMapper(FE_C_H_Q1_3D, FE_C_H_UL1_3D, Q1_Q1);
  RegisterFE3DMapper(FE_C_H_UL1_3D, FE_C_H_Q1_3D, Q1_Q1);  
  
  RegisterFE3DMapper1Reg(FE_C_H_UL1_3D, FE_C_H_UL1_3D, Q1_Q1_1Reg);
  RegisterFE3DMapper1Reg(FE_C_H_UL1_3D, FE_C_H_Q1_3D, Q1_Q1_1Reg);
  RegisterFE3DMapper1Reg(FE_C_H_Q1_3D, FE_C_H_UL1_3D, Q1_Q1_1Reg);
      
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "FE mapper registered" << endl;
}

void TFEDatabase3D::RegisterAllHangingNodes()
{
  RegisterHNDesc3D(HN_C_P1_3D_E, HN_C_P1_3D_E_Obj);

  RegisterHNDesc3D(HN_C_P2_3D_E, HN_C_P2_3D_E_Obj);
  RegisterHNDesc3D(HN_C_P2_3D_F, HN_C_P2_3D_F_Obj);
  
  RegisterHNDesc3D(HN_C_P3_3D_E, HN_C_P3_3D_E_Obj);
  RegisterHNDesc3D(HN_C_P3_3D_M, HN_C_P3_3D_M_Obj);
  RegisterHNDesc3D(HN_C_P3_3D_F, HN_C_P3_3D_F_Obj);
  RegisterHNDesc3D(HN_C_P3_3D_G, HN_C_P3_3D_G_Obj);
  
  RegisterHNDesc3D(HN_C_Q1_3D_E, HN_C_Q1_3D_E_Obj);
  RegisterHNDesc3D(HN_C_Q1_3D_F, HN_C_Q1_3D_F_Obj);
  
  RegisterHNDesc3D(HN_C_Q2_3D_E, HN_C_Q2_3D_E_Obj);
  RegisterHNDesc3D(HN_C_Q2_3D_F, HN_C_Q2_3D_F_Obj);
  
  RegisterHNDesc3D(HN_C_Q3_3D_1, HN_C_Q3_3D_1_Obj);
  RegisterHNDesc3D(HN_C_Q3_3D_2, HN_C_Q3_3D_2_Obj);
  RegisterHNDesc3D(HN_C_Q3_3D_3, HN_C_Q3_3D_3_Obj);
  RegisterHNDesc3D(HN_C_Q3_3D_4, HN_C_Q3_3D_4_Obj);
  RegisterHNDesc3D(HN_C_Q3_3D_5, HN_C_Q3_3D_5_Obj);
  
  RegisterHNDesc3D(HN_N_P1_3D_E, HN_N_P1_3D_E_Obj);

  RegisterHNDesc3D(HN_N_P2_3D_0, HN_N_P2_3D_0_Obj);
  RegisterHNDesc3D(HN_N_P2_3D_1, HN_N_P2_3D_1_Obj);
  RegisterHNDesc3D(HN_N_P2_3D_2, HN_N_P2_3D_2_Obj);
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "Hanging nodes registered" << endl;
}

void TFEDatabase3D::RegisterAllRefTrans()
{
  ReferenceTrans3D[TetraAffin] = new TTetraAffin();
  ReferenceTrans3D[TetraIsoparametric] = new TTetraIsoparametric();
  ReferenceTrans3D[HexaAffin] = new THexaAffin();
  ReferenceTrans3D[HexaTrilinear] = new THexaTrilinear();
  ReferenceTrans3D[HexaIsoparametric] = new THexaIsoparametric();
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << "Reference Transformations registered" << endl;
}

void TFEDatabase3D::GenerateArrays()
{
  int i;
  TFE3D *ele;
  TBaseFunct3D *bf;

  for(i=0;i<N_FEs3D;i++)
  {
    ele = FEs3D[i];
    if(ele)
    {
      ele->GetFEDesc3D(FEDesc3D_IDFromFE3D[i],
                       FEDesc3DFromFE3D[i]);

      ele->GetBaseFunct3D(BaseFunct3D_IDFromFE3D[i],
                          bf);
      BaseFunct3DFromFE3D[i] = bf; 
      N_BaseFunctFromFE3D[i] = bf->GetDimension();
      PolynomialDegreeFromFE3D[i] = bf->GetPolynomialDegree();
      AccuracyFromFE3D[i] = bf->GetAccuracy();
      RefElementFromFE3D[i] = bf->GetRefElement();

      ele->GetNodalFunctional3D(NodalFunctional3D_IDFromFE3D[i],
                                NodalFunctional3DFromFE3D[i]);

      RefTrans3D_IDFromFE3D[i] = ele->GetRefTransID(); 
    } // endif
  } // endfor i

  QFHexaFromDegree[0] = Gauss2Hexa;
  QFHexaFromDegree[1] = Gauss2Hexa;
  QFHexaFromDegree[2] = Gauss2Hexa;
  QFHexaFromDegree[3] = Gauss3Hexa;
  QFHexaFromDegree[4] = Gauss3Hexa;
  QFHexaFromDegree[5] = Gauss3Hexa;
  QFHexaFromDegree[6] = Gauss4Hexa;
//  QFHexaFromDegree[6] = Degree7_Points38;
  QFHexaFromDegree[7] = Gauss4Hexa;             
  QFHexaFromDegree[8] = Gauss5Hexa;
  QFHexaFromDegree[9] = Gauss5Hexa;
  QFHexaFromDegree[10] = Gauss6Hexa;
  QFHexaFromDegree[11] = Gauss6Hexa;
  QFHexaFromDegree[12] = Gauss7Hexa;
  QFHexaFromDegree[13] = Gauss7Hexa;
  QFHexaFromDegree[14] = Gauss8Hexa;
  QFHexaFromDegree[15] = Gauss8Hexa;
  QFHexaFromDegree[16] = Gauss9Hexa;
  QFHexaFromDegree[17] = Gauss9Hexa;
  
  QFConvolutionHexaFromDegree[0] = VerticesAndOrigin;
  QFConvolutionHexaFromDegree[1] = VerticesAndOrigin15;
  QFConvolutionHexaFromDegree[2] = VerticesAndOrigin57;

  QFTetraFromDegree[0] = P2Tetra;
  QFTetraFromDegree[1] = P2Tetra;
  QFTetraFromDegree[2] = P2Tetra;
  QFTetraFromDegree[3] = P4Tetra;
  QFTetraFromDegree[4] = P4Tetra;
  QFTetraFromDegree[5] = P5Tetra;
  QFTetraFromDegree[6] = P5Tetra;
  QFTetraFromDegree[7] = P8Tetra;
  QFTetraFromDegree[8] = P8Tetra;
  QFTetraFromDegree[9] = P8Tetra;
  QFTetraFromDegree[10] = P8Tetra;

  QFQuadFromDegree[0] = Gauss2Quad;
  QFQuadFromDegree[1] = Gauss2Quad;
  QFQuadFromDegree[2] = Gauss2Quad;
  QFQuadFromDegree[3] = Gauss2Quad;
  QFQuadFromDegree[4] = Gauss3Quad;
  QFQuadFromDegree[5] = Gauss3Quad;
  QFQuadFromDegree[6] = Gauss4Quad;
  QFQuadFromDegree[7] = Gauss4Quad;
  QFQuadFromDegree[8] = Gauss5Quad;
  QFQuadFromDegree[9] = Gauss5Quad;
  QFQuadFromDegree[10] = Gauss6Quad;
  QFQuadFromDegree[11] = Gauss6Quad;
  QFQuadFromDegree[12] = Gauss7Quad;
  QFQuadFromDegree[13] = Gauss7Quad;
  QFQuadFromDegree[14] = Gauss8Quad;
  QFQuadFromDegree[15] = Gauss8Quad;
  QFQuadFromDegree[16] = Gauss9Quad;
  QFQuadFromDegree[17] = Gauss9Quad;

  QFTriaFromDegree[0] = MidPointTria;
  QFTriaFromDegree[1] = MidPointTria;
  QFTriaFromDegree[2] = MidPointTria;
  QFTriaFromDegree[3] = SevenPointTria;
  QFTriaFromDegree[4] = Gauss3Tria;
  QFTriaFromDegree[5] = Gauss3Tria;
  QFTriaFromDegree[6] = Degree8Tria;
  QFTriaFromDegree[7] = Degree8Tria;
  QFTriaFromDegree[8] = Degree8Tria;
  QFTriaFromDegree[9] = Degree19Tria;
  QFTriaFromDegree[10] = Degree19Tria;
  QFTriaFromDegree[11] = Degree19Tria;
  QFTriaFromDegree[12] = Degree19Tria;
  QFTriaFromDegree[13] = Degree19Tria;
  QFTriaFromDegree[14] = Degree19Tria;
  QFTriaFromDegree[15] = Degree19Tria;
  QFTriaFromDegree[16] = Degree19Tria;
  QFTriaFromDegree[17] = Degree19Tria;
  QFTriaFromDegree[18] = Degree19Tria;
  QFTriaFromDegree[19] = Degree19Tria;
  
  QFLineFromDegree[0] = Gauss1Line;
  QFLineFromDegree[1] = Gauss1Line;
  QFLineFromDegree[2] = Gauss2Line;
  QFLineFromDegree[3] = Gauss2Line;
  QFLineFromDegree[4] = Gauss3Line;
  QFLineFromDegree[5] = Gauss3Line;
  QFLineFromDegree[6] = Gauss4Line;
  QFLineFromDegree[7] = Gauss4Line;
  QFLineFromDegree[8] = Gauss5Line;
  QFLineFromDegree[9] = Gauss5Line;
  QFLineFromDegree[10] = Gauss6Line;
  QFLineFromDegree[11] = Gauss6Line;
  QFLineFromDegree[12] = Gauss7Line;
  QFLineFromDegree[13] = Gauss7Line;
  QFLineFromDegree[14] = Gauss8Line;
  QFLineFromDegree[15] = Gauss8Line;
  QFLineFromDegree[16] = Gauss9Line;
  QFLineFromDegree[17] = Gauss9Line;
  QFLineFromDegree[18] = Gauss10Line;
  QFLineFromDegree[19] = Gauss10Line;
  QFLineFromDegree[20] = Gauss11Line;
  QFLineFromDegree[21] = Gauss11Line;
  QFLineFromDegree[22] = Gauss12Line;
  QFLineFromDegree[23] = Gauss12Line;
  
  HighestAccuracyLine = 23;
} // end TFEDatabase3D::GenerateArrays


/** calculate base functions with derivatives and coordinates
    from reference to original element */
RefTrans3D TFEDatabase3D::GetOrig(int N_LocalUsedElements, 
                          FE3D *LocalUsedElements,
                          TCollection *Coll,
                          TBaseCell *cell, bool *Needs2ndDer,
                          int &N_Points, double* &xi, double* &eta,
                          double* &zeta, double* &weights,
                          double* X, double* Y, double* Z, double* absdetjk)
{
  int i,j, MaxPolynomialDegree, PolynomialDegree, N_Faces, N_terms;
  BF3DRefElements RefElement;
  QuadFormula3D QuadFormula;
  TQuadFormula3D *qf2;
  RefTrans3D RefTrans, *RefTransArray, CurrentRefTrans;
  TRefTrans3D *rt;
  BaseFunct3D BaseFuncts[N_FEs3D];
  JointType jointtype;
  TJoint *joint;
  bool IsIsoparametric;
  BoundTypes bdtype;
  int MaxApproxOrder, ApproxOrder;
  int vertex = 0;

  BaseFunct3D BaseFunct;
  double **origvaluesD000, **origvaluesD100, **origvaluesD010, **origvaluesD001;
  double **origvaluesD200, **origvaluesD110, **origvaluesD101, **origvaluesD020;
  double **origvaluesD011, **origvaluesD002;
  int N_Functs;

#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
#endif 

  // find adequate quadrature formula for all elements
  // and find needed reference transformation
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();
  RefTrans = RefTransArray[LocalUsedElements[0]];
  MaxPolynomialDegree = 0;
  MaxApproxOrder = 0;
  RefElement = 
        TFEDatabase3D::GetRefElementFromFE3D(LocalUsedElements[0]);
  for(i=0;i<N_LocalUsedElements;i++)
  {
    BaseFuncts[i] = 
        TFEDatabase3D::GetBaseFunct3D_IDFromFE3D(LocalUsedElements[i]);
    PolynomialDegree = PolynomialDegreeFromFE3D[LocalUsedElements[i]];
    if(PolynomialDegree > MaxPolynomialDegree) 
      MaxPolynomialDegree = PolynomialDegree;
    ApproxOrder = AccuracyFromFE3D[LocalUsedElements[i]];
    if(ApproxOrder > MaxApproxOrder)
      MaxApproxOrder = ApproxOrder;
    CurrentRefTrans = RefTransArray[LocalUsedElements[i]];
    if(CurrentRefTrans > RefTrans)
      RefTrans = CurrentRefTrans;
  }
/*
  cout << "MaxPolynomialDegree: " << MaxPolynomialDegree << endl;
  cout << "RefElement: " << RefElement << endl;
  cout << "RefTrans: " << RefTrans << endl;
*/
  
  if (TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR)
    N_terms = 2;
  else
    N_terms = 3;

  switch(RefElement)
  {
    case BFUnitHexahedron:
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE==0)
      {
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(N_terms*MaxPolynomialDegree);
	
        if (TDatabase::ParamDB->INTERNAL_QUAD_HEXA<N_terms*MaxPolynomialDegree)
        {
#ifdef _MPI
        if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif
         {
          switch(N_terms*MaxPolynomialDegree)
          {
            case 0:
              OutPut("Quadrature formula for hexahedra is Gauss2" << endl); 
              break;
            case 3:
              OutPut("Quadrature formula for hexahedra is Gauss3" << endl); 
              break;
            case 6:
               OutPut("Quadrature formula for hexahedra is Gauss4" << endl); 
               //OutPut("Quadrature formula for hexahedra is Degree7_Points38" << endl); 
             break;
            case 9:
              OutPut("Quadrature formula for hexahedra is Gauss5" << endl); 
              break;
            case 12:
              OutPut("Quadrature formula for hexahedra is Gauss7" << endl); 
              break;
            case 15:
              OutPut("Quadrature formula for hexahedra is Gauss8" << endl); 
              break;
          }
         }
          TDatabase::ParamDB->INTERNAL_QUAD_HEXA = N_terms*MaxPolynomialDegree;
        }
      }
      else
      {      
        if (TDatabase::ParamDB->INTERNAL_QUAD_RULE==1)
        {
          // only for convolution
           QuadFormula = TFEDatabase3D::GetQFConvolutionHexaFromDegree(2);
          //OutPut("Quadrature formula for convolution on hexahedra is VerticesAndOrigin57" << endl); 
        }
        else
        {
          if (TDatabase::ParamDB->INTERNAL_QUAD_RULE==2)
          {
          // only adaptive VMS
           QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(3);
          //OutPut("Quadrature formula for convolution on hexahedra is VerticesAndOrigin57" << endl); 
          }
        }
       }
      N_Faces = 6;
    break;

    case BFUnitTetrahedron:
      QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(N_terms*MaxPolynomialDegree);
      
      if (TDatabase::ParamDB->INTERNAL_QUAD_TETRA<N_terms*MaxPolynomialDegree)
      {
#ifdef _MPI
      if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif
       {
        switch(N_terms*MaxPolynomialDegree)
        {
          case 0:
            OutPut("Quadrature formula for tetrahedra is P2Tetra" << endl); 
            break;
          case 3:
            OutPut("Quadrature formula for tetrahedra is P4Tetra" << endl); 
            break;
          case 6:
            OutPut("Quadrature formula for tetrahedra is P5Tetra" << endl); 
            break;
          case 9:
            OutPut("Quadrature formula for tetrahedra is P8Tetra" << endl); 
            break;
         }
       }
        TDatabase::ParamDB->INTERNAL_QUAD_TETRA = N_terms*MaxPolynomialDegree;
      }
      else
      {
          if (TDatabase::ParamDB->INTERNAL_QUAD_RULE==2)
          {
          // only adaptive VMS
           QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(3);
	   //OutPut("Quadrature formula for resolved small scales set" << endl); 
          }
      }
      N_Faces = 4;
    break;

    default:
      N_Faces = 0;
  } // endswitch

// =================================================================  
//   // only for terahertz problem (sashi)
//   if(cell->IsLayerCell())
//    {
//     switch(RefElement)
//      { 
//       case BFUnitTetrahedron: 
//        QuadFormula = VertexTetra;
//        //OutPut("LayerCell Quadrature formula for tetrahedra is P8Tetra" << endl); 
//       break;
//       default:
// 	
//       break;	
//       
//      } // endswitch 
//    } //if(cell->IsLayerCell())
// =================================================================

  IsIsoparametric = FALSE;
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
  {
    for(i=0;i<N_Faces;i++)
    {
      joint = cell->GetJoint(i);
      jointtype = joint->GetType();
      if(jointtype == BoundaryFace)
      {
        bdtype = ((TBoundFace *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Plane)
          IsIsoparametric = TRUE;
      }
      if(jointtype == InterfaceJoint3D)
      {
        bdtype = ((TInterfaceJoint3D *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Plane)
          IsIsoparametric = TRUE;
      }
      if(jointtype == IsoInterfaceJoint3D)
        IsIsoparametric = TRUE;

      if(jointtype == IsoJointEqN)
        IsIsoparametric = TRUE;

      if(jointtype == IsoBoundFace)
        IsIsoparametric = TRUE;
    } // endfor
    
    if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY==1)
    // for 3d benchmark on tetrahedral meshes 
    // set all cells isoparametric which has two vertices on the cylinder
    // - outcomment before that all cells become isoparametric which
    //   have a face on the boundary
    {
      if (N_Faces==4)
      {
        double xp,yp,zp;
        for(i=0;i<4;i++)
        {
          cell->GetVertex(i)->GetCoords(xp, yp, zp);
          if(fabs(sqrt((xp-0.5)*(xp-0.5)+(yp-0.2)*(yp-0.2))-0.05)<1e-6)
            vertex++;
        }
        if (vertex>=2)
          IsIsoparametric = TRUE;
      }
    } // endif  
        
  }// endif 

  if(IsIsoparametric)
  {
    switch(RefElement)
    {
      case BFUnitHexahedron:
        RefTrans = HexaIsoparametric;
      break;

      case BFUnitTetrahedron:
        RefTrans = TetraIsoparametric;
      break;
    }
  } // endif IsIsoparametric

  // cout << "QuadFormula: " << QuadFormula << endl;
  qf2 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
  qf2->GetFormulaData(N_Points, weights, xi, eta, zeta);

  // calculate the values of base functions and their derivatives
  // on the original element
  rt = TFEDatabase3D::ReferenceTrans3D[RefTrans];
  switch(RefTrans)
  {
    case TetraAffin:
      // cout << "TetraAffin" << endl;
      ((TTetraAffin *)rt)->SetCell(cell);
      ((TTetraAffin *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta, zeta,
                               QuadFormula,
                               Needs2ndDer);
      ((TTetraAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case TetraIsoparametric:
      // cout << "TetraIsoparametric" << endl;
      ((TTetraIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
      ((TTetraIsoparametric *)rt)->SetCell(cell);
      ((TTetraIsoparametric *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta, zeta,
                               QuadFormula,
                               Needs2ndDer);
      ((TTetraIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case HexaAffin:
      // cout << "HexaAffin" << endl;
      ((THexaAffin *)rt)->SetCell(cell);
      ((THexaAffin *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta, zeta,
                               QuadFormula,
                               Needs2ndDer);
      ((THexaAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case HexaTrilinear:
      // cout << "HexaTrilinear" << endl;
      ((THexaTrilinear *)rt)->SetCell(cell);
      ((THexaTrilinear *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta, zeta,
                               QuadFormula,
                               Needs2ndDer);
      ((THexaTrilinear *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case HexaIsoparametric:
      // cout << "HexaIsoparametric" << endl;
      ((THexaIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
      ((THexaIsoparametric *)rt)->SetCell(cell);
      ((THexaIsoparametric *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta, zeta,
                               QuadFormula,
                               Needs2ndDer);
      ((THexaIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
  } // endswitch

  for(i=0;i<N_LocalUsedElements;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
    origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);

    BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                      origvaluesD000);
    BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                      origvaluesD100);
    BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                      origvaluesD010);
    BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                      origvaluesD001);

    if(Needs2ndDer[i])
    {
      origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
      origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
      origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
      origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
      origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
      origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
      BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD200);
      BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD110);
      BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD101);
      BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD020);
      BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD011);
      BaseFuncts3D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD002);
    } // endif Needs2ndDer[i]
  } // endfor i

  return RefTrans;
}

/** calculate points on original element */
void TFEDatabase3D::GetOrigFromRef(RefTrans3D RefTrans, int n_points,
                   double *xi, double *eta, double *zeta,
                   double *X, double *Y, double *Z, double *absdetjk)
{
  TRefTrans3D *rt;

  rt = TFEDatabase3D::ReferenceTrans3D[RefTrans];
  switch(RefTrans)
  {
    case TetraAffin:
      // cout << "TetraAffin" << endl;
      ((TTetraAffin *)rt)->GetOrigFromRef(n_points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case TetraIsoparametric:
      // cout << "TetraIsoparametric" << endl;
      ((TTetraIsoparametric *)rt)->GetOrigFromRef(n_points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case HexaAffin:
      // cout << "HexaAffin" << endl;
      ((THexaAffin *)rt)->GetOrigFromRef(n_points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case HexaTrilinear:
      // cout << "HexaTrilinear" << endl;
      ((THexaTrilinear *)rt)->GetOrigFromRef(n_points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
    case HexaIsoparametric:
      // cout << "HexaIsoparametric" << endl;
      ((THexaIsoparametric *)rt)->GetOrigFromRef(n_points, xi, eta, zeta,
                                         X, Y, Z, absdetjk);
    break;
  } // endswitch
}


/** calculate the values of base functions or their derivatives
    on the original element */
void TFEDatabase3D::GetOrigValues(RefTrans3D RefTrans,
                double xi, double eta, double zeta,
                TBaseFunct3D *bf, TCollection *Coll, TBaseCell *cell,
                double *uref, double *uxiref, double *uetaref, double *uzetaref,
                double *uorig, double *uxorig, double *uyorig, double *uzorig)
{
  TRefTrans3D *rt;
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();
  rt = ReferenceTrans3D[RefTrans];
  switch(RefTrans)
  {
    case TetraAffin: 
      ((TTetraAffin *)rt)->GetOrigValues(xi, eta, zeta, N_BaseFunct, 
          uref, uxiref, uetaref, uzetaref, uorig, uxorig, uyorig, uzorig, 
          BaseVectDim);
      break;
    case TetraIsoparametric:
      ((TTetraIsoparametric *)rt)->GetOrigValues(xi, eta, zeta, N_BaseFunct, 
          uref, uxiref, uetaref, uzetaref, uorig, uxorig, uyorig, uzorig);
      break;
    case HexaAffin:
      ((THexaAffin *)rt)->GetOrigValues(xi, eta, zeta, N_BaseFunct, 
          uref, uxiref, uetaref, uzetaref, uorig, uxorig, uyorig, uzorig, 
          BaseVectDim);
      break;
    case HexaTrilinear:
      ((THexaTrilinear *)rt)->GetOrigValues(xi, eta, zeta, N_BaseFunct, 
          uref, uxiref, uetaref, uzetaref, uorig, uxorig, uyorig, uzorig,
          BaseVectDim);
      break;
    case HexaIsoparametric:
      ((THexaIsoparametric *)rt)->GetOrigValues(xi, eta, zeta, N_BaseFunct, 
          uref, uxiref, uetaref, uzetaref, uorig, uxorig, uyorig, uzorig);
      break;
    default:
      cout << "wrong reference transformation identifier" << endl;
      break;
  } // endswitch

  bf->ChangeBF(Coll, cell, uorig);
  bf->ChangeBF(Coll, cell, uxorig);
  bf->ChangeBF(Coll, cell, uyorig);
  bf->ChangeBF(Coll, cell, uzorig);
}

/** set cell for reference transformation */
void TFEDatabase3D::SetCellForRefTrans(TBaseCell *cell, 
                                     RefTrans3D reftrans)
{
  TRefTrans3D *rt;

  rt = ReferenceTrans3D[reftrans];

  switch(reftrans)
  {
    case TetraAffin: 
      ((TTetraAffin *)rt)->SetCell(cell);
      break;
    case TetraIsoparametric:
      ((TTetraIsoparametric *)rt)->SetCell(cell);
      break;
    case HexaAffin:
      ((THexaAffin *)rt)->SetCell(cell);
      break;
    case HexaTrilinear:
      ((THexaTrilinear *)rt)->SetCell(cell);
      break;
    case HexaIsoparametric:
      ((THexaIsoparametric *)rt)->SetCell(cell);
      break;
    default:
      cout << "wrong reference transformation identifier" << endl;
      break;
  } // endswitch
}
/** calculate points on reference element */
void TFEDatabase3D::GetRefFromOrig (RefTrans3D RefTrans,
                                 double X, double Y, double Z,
                                 double &xi, double &eta, double &zeta)
{
  TRefTrans3D *rt;

  rt = ReferenceTrans3D[RefTrans];
  switch(RefTrans)
  {
    case TetraAffin: 
      ((TTetraAffin *)rt)->GetRefFromOrig(X, Y, Z,  xi, eta, zeta);
      break;
    case TetraIsoparametric: 
      ((TTetraIsoparametric *)rt)->GetRefFromOrig(X, Y, Z,  xi, eta, zeta);
      break;
    case HexaAffin: 
      ((THexaAffin *)rt)->GetRefFromOrig(X, Y, Z,  xi, eta, zeta);
      break;
    case HexaTrilinear: 
      ((THexaTrilinear *)rt)->GetRefFromOrig(X, Y, Z,  xi, eta, zeta);
      break;
    case HexaIsoparametric: 
      ((THexaIsoparametric *)rt)->GetRefFromOrig(X, Y, Z,  xi, eta, zeta);
      break;
   default:
     cout << "wrong reference transformation identifier" << RefTrans << endl;
      break;
  } // endswitch
}

double *TFEDatabase3D::GetProlongationMatrix3D (FE3D parent, 
    Refinements refine, FE3D child, int childnumber)
{ 
  double *ret, *ret2;
  int i,j,k,l;
  int N_Coarse, N_Fine, N_Points, N_Children;
  double *xi, *eta, *zeta;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double AllPointValues[MaxN_PointsForNodal3D][MaxN_BaseFunctions3D];
  double PointValues[MaxN_PointsForNodal3D];
  TFE3D *CoarseElement, *FineElement;
  TRefDesc *RefDesc;
  TBaseFunct3D *BaseFunctions;
  BaseFunct3D Coarse, Fine;
  TGridCell *RefCell, *cell;
  TNodalFunctional3D *nf;
  BF3DRefElements RefElement;
  RefTrans3D F_K;
  TRefTrans3D *rt;

  CoarseElement = TFEDatabase3D::GetFE3D(parent);
  Coarse = CoarseElement->GetBaseFunct3D_ID();
  FineElement = TFEDatabase3D::GetFE3D(child);
  Fine = FineElement->GetBaseFunct3D_ID();

  ret = ProlongationMatrix3D[Coarse][refine][Fine][childnumber]; 

  if(ret == NULL)
  {
    // cerr << "ret == NULL" << endl;
    // prolongation matrix was not generated yet

    BaseFunctions = CoarseElement->GetBaseFunct3D();
    N_Coarse = BaseFunctions->GetDimension();

    if(refine == NoRef)
    {
      ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
      ret2 = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];

      RefCell = BaseFunctions->GenerateRefElement();
      FineElement = TFEDatabase3D::GetFE3D(child);
      N_Fine = FineElement->GetBaseFunct3D()->GetDimension();

      nf = FineElement->GetNodalFunctional3D();
      nf->GetPointsForAll(N_Points, xi, eta, zeta);

      RefElement = BaseFunctions->GetRefElement();

      switch(RefElement)
      {
        case BFUnitHexahedron:
          rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
          ((THexaAffin *)rt)->SetCell(RefCell);
          F_K = HexaAffin;
          break;
        case BFUnitTetrahedron:
          rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
          ((TTetraAffin *)rt)->SetCell(RefCell);
          F_K = TetraAffin;
          break;
        default:
          F_K = TetraAffin;
      }
      TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                  X, Y, Z, AbsDetjk);

      for(k=0;k<N_Points;k++)
        BaseFunctions->GetDerivatives(D000, X[k], Y[k], Z[k],
                                      AllPointValues[k]);

      for(k=0;k<N_Coarse;k++)
      {
        for(l=0;l<N_Points;l++)
          PointValues[l] = AllPointValues[l][k];

        nf->GetAllFunctionals(NULL, NULL, PointValues, ret2+k*MaxN_BaseFunctions3D);
      }

      for(k=0;k<MaxN_BaseFunctions3D;k++)
        for(l=0;l<MaxN_BaseFunctions3D;l++)
          ret[k*MaxN_BaseFunctions3D+l] = 
                        ret2[l*MaxN_BaseFunctions3D+k];
  
      RegisterProlongationMatrix3D(Coarse, refine, Fine, 
                                   childnumber, ret);
    }
    else
    {
      RefDesc = TDatabase::RefDescDB[refine+N_SHAPES];
      N_Children = RefDesc->GetN_Children();
  
      RefCell = BaseFunctions->GenerateRefElement();
      RefCell->SetRefDesc(RefDesc);
      RefCell->Refine(1);
  
      for(j=0;j<N_Children;j++)
      {
        ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
        ret2 = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  
        cell = (TGridCell *)RefCell->GetChild(j);
        FineElement = TFEDatabase3D::GetFE3D(child);
        Fine = FineElement->GetBaseFunct3D_ID();
        N_Fine = FineElement->GetBaseFunct3D()->GetDimension();
  
        nf = FineElement->GetNodalFunctional3D();
        nf->GetPointsForAll(N_Points, xi, eta, zeta);
  
        switch(cell->GetType())
        {
          case Hexahedron:
          case Brick:
            rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
            ((THexaAffin *)rt)->SetCell(cell);
            F_K = HexaAffin;
          break;
          case Tetrahedron:
            rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
            ((TTetraAffin *)rt)->SetCell(cell);
            F_K = TetraAffin;
          break;
	  default:
	    
	  break;
	  
        }
        TFEDatabase3D::GetOrigFromRef(F_K ,N_Points, xi, eta, zeta,
                                    X, Y, Z, AbsDetjk);
  
        for(k=0;k<N_Points;k++)
          BaseFunctions->GetDerivatives(D000, X[k], Y[k], Z[k],
                        AllPointValues[k]);
  
        for(k=0;k<N_Coarse;k++)
        {
          for(l=0;l<N_Points;l++)
            PointValues[l] = AllPointValues[l][k];
  
          nf->GetAllFunctionals(NULL, NULL, PointValues, ret2+k*MaxN_BaseFunctions3D);
        }

        for(k=0;k<MaxN_BaseFunctions3D;k++)
          for(l=0;l<MaxN_BaseFunctions3D;l++)
            ret[k*MaxN_BaseFunctions3D+l] = ret2[l*MaxN_BaseFunctions3D+k];
  
        RegisterProlongationMatrix3D(Coarse, refine, Fine, j, ret);
      } // endfor j
    }

    ret = ProlongationMatrix3D[Coarse][refine][Fine][childnumber]; 

    // RefCell->Derefine();
    // delete (TGridCell *)RefCell;
  }

  return ret;
}

double *TFEDatabase3D::GetRestrictionMatrix3D (FE3D parent, 
    Refinements refine, FE3D child, int childnumber)
{ 
  double *ret, *ret2;
  int i,j,k,l, l1, l2;
  int N_Coarse, N_Fine, N_Points, N_Children;
  double AllPointValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];
  double PointValues[MaxN_PointsForNodal3D];
  TFE3D *CoarseElement, *FineElement;
  TRefDesc *RefDesc;
  TBaseFunct3D *BaseFunctions, *FineBF;
  BaseFunct3D Coarse, Fine;
  TBaseCell *RefCell, *cell;
  TNodalFunctional3D *nf;
  RefTrans3D F_K;
  TRefTrans3D *rt;

  double G[MaxN_BaseFunctions3D][MaxN_BaseFunctions3D];
  double Gret[MaxN_BaseFunctions3D][MaxN_BaseFunctions3D];
  double R[MaxN_BaseFunctions3D][MaxN_BaseFunctions3D];
  TQuadFormula3D *Formula;
  QuadFormula3D QuadFormula;
  QuadFormula2D LineQuadFormula;
  double **CoarseBFData, **FineBFData, *PointData;
  double *FinePointData;
  int N_QuadPoints;
  double *xi, *eta, *zeta, *weights, sum, w;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  int LDA=MaxN_BaseFunctions3D;

  CoarseElement = TFEDatabase3D::GetFE3D(parent);
  Coarse = CoarseElement->GetBaseFunct3D_ID();
  FineElement = TFEDatabase3D::GetFE3D(child);
  Fine = FineElement->GetBaseFunct3D_ID();

  ret = RestrictionMatrix3D[Coarse][refine][Fine][childnumber]; 

  if(ret == NULL)
  {
    // restriction matrix was not generated yet

    BaseFunctions = CoarseElement->GetBaseFunct3D();
    N_Coarse = BaseFunctions->GetDimension();

    memset(G, 0, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D*
                 SizeOfDouble);

    // build matrix G, gij = (uiH, ujH)
    TQuadFormula3D::FindLocalQuadFormula3D
            (1, &parent, LineQuadFormula, QuadFormula);
    Formula = GetQuadFormula3D(QuadFormula);
    Formula->GetFormulaData(N_QuadPoints, weights, xi, eta, zeta);

    TFEDatabase3D::GetBaseFunct3D(Coarse)->MakeRefElementData(QuadFormula);
    CoarseBFData = GetRefElementValues(Coarse, QuadFormula, D000);
    for(k=0;k<N_QuadPoints;k++)
    {
      PointData = CoarseBFData[k];
      w = weights[k];
      for(i=0;i<N_Coarse;i++)
      {
        for(j=0;j<N_Coarse;j++)
        {
          // G is symmetric
          G[j][i] += w*PointData[i]*PointData[j];
        }
      }
    } // endfor k

    memcpy(Gret, G, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D
                    *SizeOfDouble);

    if(refine == NoRef)
    {
      memset(R, 0, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D*
                   SizeOfDouble);

      BaseFunctions = FineElement->GetBaseFunct3D();
      N_Fine = BaseFunctions->GetDimension();

      // build matrix R, rij = (uiH, ujh)
      TQuadFormula3D::FindLocalQuadFormula3D
              (1, &child, LineQuadFormula, QuadFormula);
      Formula = GetQuadFormula3D(QuadFormula);
      Formula->GetFormulaData(N_QuadPoints, weights, xi, eta, zeta);
      
      TFEDatabase3D::GetBaseFunct3D(Fine)->MakeRefElementData(QuadFormula);
      FineBFData = GetRefElementValues(Fine, QuadFormula, D000);
      TFEDatabase3D::GetBaseFunct3D(Coarse)->MakeRefElementData(QuadFormula);
      CoarseBFData = GetRefElementValues(Coarse, QuadFormula, D000);
      for(k=0;k<N_QuadPoints;k++)
      {
        FinePointData = FineBFData[k];
        PointData = CoarseBFData[k];
        w = weights[k];
        for(l1=0;l1<N_Coarse;l1++)
        {
          for(l2=0;l2<N_Fine;l2++)
          {
            // = r(l1,l2) if row stored
            R[l2][l1] += w*PointData[l1]*FinePointData[l2];
          }
        }
      } // endfor k

      memcpy(G, Gret, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D
                      *SizeOfDouble);

      // determine inv(G)*R
      SolveMultipleSystems((double *)G, (double *)R, 
                           N_Coarse, LDA, LDA, N_Fine);

      ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
      for(l1=0;l1<N_Coarse;l1++)
        for(l2=0;l2<N_Fine;l2++)
          ret[l1*MaxN_BaseFunctions3D+l2] = R[l2][l1]; 
  
      RegisterRestrictionMatrix3D(Coarse, refine, Fine, childnumber, ret);
    }
    else
    {
      RefDesc = TDatabase::RefDescDB[refine+N_SHAPES];
      N_Children = RefDesc->GetN_Children();
  
      RefCell = BaseFunctions->GenerateRefElement();
      RefCell->SetRefDesc(RefDesc);
      RefCell->Refine(1);
  
      for(j=0;j<N_Children;j++)
      {
        // cout << "child: " << j << endl;
        memset(R, 0, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D*
                     SizeOfDouble);
  
        cell = RefCell->GetChild(j);
        FineElement = TFEDatabase3D::GetFE3D(child);
        Fine = FineElement->GetBaseFunct3D_ID();
        FineBF = FineElement->GetBaseFunct3D();
        N_Fine = FineBF->GetDimension();
  
        TQuadFormula3D::FindLocalQuadFormula3D
                (1, &child, LineQuadFormula, QuadFormula);
        Formula = GetQuadFormula3D(QuadFormula);
        Formula->GetFormulaData(N_QuadPoints, weights, xi, eta, zeta);

        TFEDatabase3D::GetBaseFunct3D(Fine)->MakeRefElementData(QuadFormula);
        FineBFData = GetRefElementValues(Fine, QuadFormula, D000);
  
        F_K = FineElement->GetRefTransID();
  
        switch(F_K)
        {
          case HexaAffin:
          case HexaTrilinear:
          case HexaIsoparametric:
            rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
            ((THexaAffin *)rt)->SetCell(RefCell->GetChild(j));
            F_K = HexaAffin;
            break;
          case TetraAffin:
          case TetraIsoparametric:
            rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
            ((TTetraAffin *)rt)->SetCell(RefCell->GetChild(j));
            break;
        }
        TFEDatabase3D::GetOrigFromRef(F_K ,N_QuadPoints, xi, eta, zeta,
                                    X, Y, Z, AbsDetjk);
  
        for(k=0;k<N_QuadPoints;k++)
          BaseFunctions->GetDerivatives(D000, X[k], Y[k], Z[k], AllPointValues[k]);

        for(k=0;k<N_QuadPoints;k++)
        {
          FinePointData = FineBFData[k];
          PointData = AllPointValues[k];
          w = weights[k]*AbsDetjk[k];
          for(l1=0;l1<N_Coarse;l1++)
          {
            for(l2=0;l2<N_Fine;l2++)
            {
              R[l2][l1] += w*PointData[l1]*FinePointData[l2];
            }
          }
        } // endfor k
       
        memcpy(G, Gret, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D
                        *SizeOfDouble);

        // determine inv(G)*R
        SolveMultipleSystems((double *)G, (double *)R, N_Coarse, 
                             LDA, LDA, N_Fine);
  
        ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
        // change from column to row storage
        for(l1=0;l1<N_Coarse;l1++)
          for(l2=0;l2<N_Fine;l2++)
            ret[l1*MaxN_BaseFunctions3D+l2] = R[l2][l1]; 
        RegisterRestrictionMatrix3D(Coarse, refine, Fine, j, ret);
      } // endfor j

      // RefCell->Derefine();
      // delete (TGridCell *)RefCell;

    }

    ret = RestrictionMatrix3D[Coarse][refine][Fine][childnumber]; 

  }

  return ret;
}
