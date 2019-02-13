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
// @(#)FEDatabase.C        1.21 05/05/00
//
// Class:       TFEDatabase2D
// Purpose:     store all used information for a FEM
//
// Author:      Gunar Matthies  09.07.98
//
// History:     start reimplementation 09.07.98 (GM)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <FEDatabase2D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <string.h>

#include <LinAlg.h>
#include <BoundEdge.h>
#include <InterfaceJoint.h>

#include <AllFEDescs1D.h>
#include <AllBaseFunctions1D.h>
#include <AllNodalFunctionals1D.h>
#include <AllFEDescs2D.h>
#include <AllBaseFunctions2D.h>
#include <AllFE2DMappers.h>
#include <AllHNDescs2D.h>
#include <AllRefTrans.h>
#include <AllNodalFunctionals2D.h>


#include <stdlib.h>
// #define __SCOTT_VOGELIUS__
// #define __ULTRA_LOCAL_PROJECTION__

// =======================================================================
// initialize static members
// =======================================================================
TQuadFormula1D *TFEDatabase2D::QuadFormulas1D[N_QuadFormulas_1D] =  { NULL };
TQuadFormula2D *TFEDatabase2D::QuadFormulas2D[N_QuadFormulas_2D] =  { NULL };
TQuadFormula3D *TFEDatabase2D::QuadFormulas3D[N_QuadFormulas_3D] = { NULL };

QuadFormula1D TFEDatabase2D::QFLineFromDegree[MAXDEGREE] = { Gauss1Line };
int TFEDatabase2D::HighestAccuracyLine = 0;

QuadFormula2D TFEDatabase2D::QFTriaFromDegree[MAXDEGREE] = { BaryCenterTria };
int TFEDatabase2D::HighestAccuracyTria = 0;

QuadFormula2D TFEDatabase2D::QFQuadFromDegree[MAXDEGREE] = { VertexQuad };
int TFEDatabase2D::HighestAccuracyQuad = 0;

TFE1D *TFEDatabase2D::FEs1D[N_FEs1D] = { NULL };
TFE2D *TFEDatabase2D::FEs2D[N_FEs2D] = { NULL };
TFEDesc1D *TFEDatabase2D::FEDescs1D[N_FEDescs1D] = { NULL };
TFEDesc2D *TFEDatabase2D::FEDescs2D[N_FEDescs2D] = { NULL };
TBaseFunct2D *TFEDatabase2D::BaseFuncts2D[N_BaseFuncts2D] = { NULL };
TBaseFunct1D *TFEDatabase2D::BaseFuncts1D[N_BaseFuncts1D] = { NULL };
TNodalFunctional2D *TFEDatabase2D::NodalFunctionals2D[N_NodalFunctionals2D] 
   = { NULL };

TNodalFunctional1D *TFEDatabase2D::NodalFunctionals1D[N_NodalFunctionals1D] 
   = { NULL };

TFE2DMapper *TFEDatabase2D::FE2DMapper[N_FEDescs2D][N_FEDescs2D] = { { NULL } };
TFE2DMapper1Reg *TFEDatabase2D::FE2DMapper1Reg[N_FEDescs2D][N_FEDescs2D] = { { NULL } };

THNDesc *TFEDatabase2D::HNDescs2D[N_HNDescs] = { NULL };

double **TFEDatabase2D::RefElementValues1D[N_BaseFuncts1D][N_QuadFormulas_1D]
        [N_MultiIndices1D] ={ { { NULL } } };
double **TFEDatabase2D::OrigElementValues1D[N_BaseFuncts1D]
        [N_MultiIndices1D] = { { NULL } };

double **TFEDatabase2D::RefElementValues2D[N_BaseFuncts2D][N_QuadFormulas_2D]
        [N_MultiIndices2D] = { { { NULL } } };
double **TFEDatabase2D::OrigElementValues2D[N_BaseFuncts2D]
        [N_MultiIndices2D] = { { NULL } };
double **TFEDatabase2D::JointValues2D[N_BaseFuncts2D][N_QuadFormulas_1D]
        [MAXN_JOINTS] = { { { NULL } } };
double **TFEDatabase2D::JointDerivatives2D[N_BaseFuncts2D][N_QuadFormulas_1D]
        [MAXN_JOINTS][N_MultiIndices2D] = { { { { NULL } } } };

TRefTrans2D *TFEDatabase2D::ReferenceTrans2D[N_RefTrans2D] = { NULL };

TRefTrans1D *TFEDatabase2D::ReferenceTrans1D[N_RefTrans1D] = { NULL };

double *TFEDatabase2D::ProlongationMatrix2D[MaxN_BaseFunctions2D]
        [N_REFDESC][MaxN_BaseFunctions2D][MAXN_CHILDREN] = { { { { NULL } } } };

double *TFEDatabase2D::RestrictionMatrix2D[MaxN_BaseFunctions2D]
        [N_REFDESC][MaxN_BaseFunctions2D][MAXN_CHILDREN] = { { { { NULL } } } };

FEDesc2D TFEDatabase2D::
           FEDesc2D_IDFromFE2D[N_FEs2D] = { FE_C_T_P0_2D };
TFEDesc2D *TFEDatabase2D::
           FEDesc2DFromFE2D[N_FEs2D] = { NULL };
BaseFunct2D TFEDatabase2D::
           BaseFunct2D_IDFromFE2D[N_FEs2D] = { BF_C_T_P0_2D };
int TFEDatabase2D::
           N_BaseFunctFromFE2D[N_FEs2D] = { 0 };
int TFEDatabase2D::
           PolynomialDegreeFromFE2D[N_FEs2D] = { 0 };
int TFEDatabase2D::
           AccuracyFromFE2D[N_FEs2D] = { 0 };
TBaseFunct2D *TFEDatabase2D::
           BaseFunct2DFromFE2D[N_FEs2D] = { NULL };
NodalFunctional2D TFEDatabase2D::
           NodalFunctional2D_IDFromFE2D[N_FEs2D] = { NF_C_T_P0_2D };
TNodalFunctional2D *TFEDatabase2D::
           NodalFunctional2DFromFE2D[N_FEs2D] = { NULL };
RefTrans2D TFEDatabase2D::
           RefTrans2D_IDFromFE2D[N_FEs2D] = { TriaAffin };
BF2DRefElements TFEDatabase2D::
           RefElementFromFE2D[N_FEs2D] = { BFUnitTriangle };

TFEDatabase2D::TFEDatabase2D()
{
  RegisterAllQuadFormulas();
  RegisterAllFEDescs();
  RegisterAllBaseFunctions();
  RegisterAllNodalFunctionals();
  RegisterAllFEs();
  RegisterAllFEMappers();
  RegisterAllHangingNodes();
  RegisterAllRefTrans();

  GenerateArrays();
}

void TFEDatabase2D::RegisterAllQuadFormulas()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  TQuadFormula1D *qf1d;
  TQuadFormulaTria *qftria;
  TQuadFormulaQuad *qfquad;
  TQuadFormula3D *qf3d;

  // =====================================================================
  // register 1d quadrature formulas
  // =====================================================================
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
  // register triangle quadratur formulas
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
  qftria->Degree9();
  RegisterQuadFormula2D(Degree9Tria, qftria);
  
  qftria = new TQuadFormulaTria();
  qftria->Degree11();
  RegisterQuadFormula2D(Degree11Tria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->Degree19();
  RegisterQuadFormula2D(Degree19Tria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->CompGauss3();
  RegisterQuadFormula2D(CompGauss3Tria, qftria);

  qftria = new TQuadFormulaTria();
  qftria->Gauss_Degree8();
  RegisterQuadFormula2D(Gauss_Degree8Tria, qftria);

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

#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "quadrature formulas registered" << endl;
}

void TFEDatabase2D::RegisterAllFEDescs()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  RegisterFEDesc1D(FE_C_L_P0_1D, FE_C_L_P0_1D_Obj);
  RegisterFEDesc1D(FE_C_L_P1_1D, FE_C_L_P1_1D_Obj);
  RegisterFEDesc1D(FE_C_L_P2_1D, FE_C_L_P2_1D_Obj);
  RegisterFEDesc1D(FE_C_L_P3_1D, FE_C_L_P3_1D_Obj);

  RegisterFEDesc1D(FE_N_L_P0_1D, FE_N_L_P0_1D_Obj);

  RegisterFEDesc1D(FE_D_L_P1_1D, FE_D_L_P1_1D_Obj);
  RegisterFEDesc1D(FE_D_L_P2_1D, FE_D_L_P2_1D_Obj);
  
  RegisterFEDesc2D(FE_C_T_P00_2D, FE_C_T_P00_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P0_2D, FE_C_T_P0_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P1_2D, FE_C_T_P1_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P2_2D, FE_C_T_P2_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P3_2D, FE_C_T_P3_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P4_2D, FE_C_T_P4_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P5_2D, FE_C_T_P5_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P6_2D, FE_C_T_P6_2D_Obj);
  
  
  RegisterFEDesc2D(FE_C_Q_Q00_2D, FE_C_Q_Q00_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q0_2D, FE_C_Q_Q0_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q1_2D, FE_C_Q_Q1_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q2_2D, FE_C_Q_Q2_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q3_2D, FE_C_Q_Q3_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q4_2D, FE_C_Q_Q4_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q5_2D, FE_C_Q_Q5_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q6_2D, FE_C_Q_Q6_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q7_2D, FE_C_Q_Q7_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q8_2D, FE_C_Q_Q8_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_Q9_2D, FE_C_Q_Q9_2D_Obj);

  RegisterFEDesc2D(FE_N_T_P1_2D, FE_N_T_P1_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_Q1_2D, FE_N_Q_Q1_2D_Obj);
  
  RegisterFEDesc2D(FE_N_Q_RT0_2D, FE_N_Q_RT0_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_RT1_2D, FE_N_Q_RT1_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_RT2_2D, FE_N_Q_RT2_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_RT3_2D, FE_N_Q_RT3_2D_Obj);
  RegisterFEDesc2D(FE_N_T_RT0_2D, FE_N_T_RT0_2D_Obj);
  RegisterFEDesc2D(FE_N_T_RT1_2D, FE_N_T_RT1_2D_Obj);
  RegisterFEDesc2D(FE_N_T_RT2_2D, FE_N_T_RT2_2D_Obj);
  RegisterFEDesc2D(FE_N_T_RT3_2D, FE_N_T_RT3_2D_Obj);

  RegisterFEDesc2D(FE_N_Q_BDM1_2D, FE_N_Q_BDM1_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_BDM2_2D, FE_N_Q_BDM2_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_BDM3_2D, FE_N_Q_BDM3_2D_Obj);
  RegisterFEDesc2D(FE_N_T_BDM1_2D, FE_N_T_BDM1_2D_Obj);
  RegisterFEDesc2D(FE_N_T_BDM2_2D, FE_N_T_BDM2_2D_Obj);
  RegisterFEDesc2D(FE_N_T_BDM3_2D, FE_N_T_BDM3_2D_Obj);

  
  RegisterFEDesc2D(FE_D_Q_P1_2D, FE_D_Q_P1_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_P2_2D, FE_D_Q_P2_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_P3_2D, FE_D_Q_P3_2D_Obj);

  RegisterFEDesc2D(FE_C_T_B2_2D, FE_C_T_B2_2D_Obj);
  RegisterFEDesc2D(FE_C_T_B3_2D, FE_C_T_B3_2D_Obj);
  RegisterFEDesc2D(FE_C_T_B4_2D, FE_C_T_B4_2D_Obj);

  RegisterFEDesc2D(FE_D_T_P1_2D, FE_D_T_P1_2D_Obj);
  RegisterFEDesc2D(FE_D_T_P2_2D, FE_D_T_P2_2D_Obj);

  RegisterFEDesc2D(FE_N_Q_Q2_2D, FE_N_Q_Q2_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_Q3_2D, FE_N_Q_Q3_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_Q4_2D, FE_N_Q_Q4_2D_Obj);
  RegisterFEDesc2D(FE_N_Q_Q5_2D, FE_N_Q_Q5_2D_Obj);

  RegisterFEDesc2D(FE_D_Q_P4_2D, FE_D_Q_P4_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_P5_2D, FE_D_Q_P5_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_P6_2D, FE_D_Q_P6_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_P7_2D, FE_D_Q_P7_2D_Obj);

  RegisterFEDesc2D(FE_N_T_P1MOD_2D, FE_N_T_P1MOD_2D_Obj);
  RegisterFEDesc2D(FE_C_T_P1MINI_2D, FE_C_T_P1MINI_2D_Obj);

  RegisterFEDesc2D(FE_N_T_P2_2D, FE_N_T_P2_2D_Obj);
  RegisterFEDesc2D(FE_N_T_P3_2D, FE_N_T_P3_2D_Obj);
  RegisterFEDesc2D(FE_N_T_P4_2D, FE_N_T_P4_2D_Obj);
  RegisterFEDesc2D(FE_N_T_P5_2D, FE_N_T_P5_2D_Obj);

  RegisterFEDesc2D(FE_D_T_P3_2D, FE_D_T_P3_2D_Obj);
  RegisterFEDesc2D(FE_D_T_P4_2D, FE_D_T_P4_2D_Obj);

  RegisterFEDesc2D(FE_B_Q_IB2_2D, FE_B_Q_IB2_2D_Obj);

  RegisterFEDesc2D(FE_D_Q_Q1_2D, FE_D_Q_Q1_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_Q2_2D, FE_D_Q_Q2_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_Q3_2D, FE_D_Q_Q3_2D_Obj);
  RegisterFEDesc2D(FE_D_Q_Q4_2D, FE_D_Q_Q4_2D_Obj);

  RegisterFEDesc2D(FE_D_Q_D2_2D, FE_D_Q_D2_2D_Obj);
  
  RegisterFEDesc2D(FE_C_Q_UL1_2D, FE_C_Q_UL1_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL2_2D, FE_C_Q_UL2_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL3_2D, FE_C_Q_UL3_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL4_2D, FE_C_Q_UL4_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL5_2D, FE_C_Q_UL5_2D_Obj);

  RegisterFEDesc2D(FE_C_T_UL1_2D, FE_C_T_UL1_2D_Obj);
  RegisterFEDesc2D(FE_C_T_UL2_2D, FE_C_T_UL2_2D_Obj);
  RegisterFEDesc2D(FE_C_T_UL3_2D, FE_C_T_UL3_2D_Obj);
  RegisterFEDesc2D(FE_C_T_UL4_2D, FE_C_T_UL4_2D_Obj);
  RegisterFEDesc2D(FE_C_T_UL5_2D, FE_C_T_UL5_2D_Obj);

  RegisterFEDesc2D(FE_C_Q_UL2S_2D, FE_C_Q_UL2S_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL3S_2D, FE_C_Q_UL3S_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL4S_2D, FE_C_Q_UL4S_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL5S_2D, FE_C_Q_UL5S_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL6S_2D, FE_C_Q_UL6S_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL7S_2D, FE_C_Q_UL7S_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL8S_2D, FE_C_Q_UL8S_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL9S_2D, FE_C_Q_UL9S_2D_Obj);

  RegisterFEDesc2D(FE_C_Q_UL2SE_2D, FE_C_Q_UL2SE_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL3SE_2D, FE_C_Q_UL3SE_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL4SE_2D, FE_C_Q_UL4SE_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL5SE_2D, FE_C_Q_UL5SE_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL6SE_2D, FE_C_Q_UL6SE_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL7SE_2D, FE_C_Q_UL7SE_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL8SE_2D, FE_C_Q_UL8SE_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_UL9SE_2D, FE_C_Q_UL9SE_2D_Obj);

  RegisterFEDesc2D(FE_C_Q_M2_2D, FE_C_Q_M2_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_M3_2D, FE_C_Q_M3_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_M4_2D, FE_C_Q_M4_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_M5_2D, FE_C_Q_M5_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_M6_2D, FE_C_Q_M6_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_M7_2D, FE_C_Q_M7_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_M8_2D, FE_C_Q_M8_2D_Obj);
  RegisterFEDesc2D(FE_C_Q_M9_2D, FE_C_Q_M9_2D_Obj);

  RegisterFEDesc2D(FE_C_Q_EL1_2D, FE_C_Q_EL1_2D_Obj);
  
  RegisterFEDesc2D(FE_C_T_SV2_2D, FE_C_T_SV2_2D_Obj);
  RegisterFEDesc2D(FE_D_T_SV1_2D, FE_D_T_SV1_2D_Obj);
  
#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "FE descriptors registered" << endl;
}

void TFEDatabase2D::RegisterAllBaseFunctions()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  RegisterBaseFunct1D(BF_C_L_P0_1D, BF_C_L_P0_1D_Obj);
  RegisterBaseFunct1D(BF_C_L_P1_1D, BF_C_L_P1_1D_Obj);
  RegisterBaseFunct1D(BF_C_L_P2_1D, BF_C_L_P2_1D_Obj);
  RegisterBaseFunct1D(BF_C_L_P3_1D, BF_C_L_P3_1D_Obj);

  RegisterBaseFunct1D(BF_D_L_P1_1D, BF_D_L_P1_1D_Obj);
  RegisterBaseFunct1D(BF_D_L_P2_1D, BF_D_L_P2_1D_Obj);
  
  RegisterBaseFunct2D(BF_C_T_P00_2D, BF_C_T_P00_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P0_2D, BF_C_T_P0_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P1_2D, BF_C_T_P1_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P2_2D, BF_C_T_P2_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P3_2D, BF_C_T_P3_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P4_2D, BF_C_T_P4_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P5_2D, BF_C_T_P5_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P6_2D, BF_C_T_P6_2D_Obj);
  
  RegisterBaseFunct2D(BF_C_Q_Q00_2D, BF_C_Q_Q00_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q0_2D, BF_C_Q_Q0_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q1_2D, BF_C_Q_Q1_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q2_2D, BF_C_Q_Q2_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q3_2D, BF_C_Q_Q3_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q4_2D, BF_C_Q_Q4_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q5_2D, BF_C_Q_Q5_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q6_2D, BF_C_Q_Q6_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q7_2D, BF_C_Q_Q7_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q8_2D, BF_C_Q_Q8_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_Q9_2D, BF_C_Q_Q9_2D_Obj);
  
  RegisterBaseFunct2D(BF_N_Q_RT0_2D, BF_N_Q_RT0_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_RT1_2D, BF_N_Q_RT1_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_RT2_2D, BF_N_Q_RT2_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_RT3_2D, BF_N_Q_RT3_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_RT0_2D, BF_N_T_RT0_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_RT1_2D, BF_N_T_RT1_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_RT2_2D, BF_N_T_RT2_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_RT3_2D, BF_N_T_RT3_2D_Obj);

  RegisterBaseFunct2D(BF_N_Q_BDM1_2D, BF_N_Q_BDM1_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_BDM2_2D, BF_N_Q_BDM2_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_BDM3_2D, BF_N_Q_BDM3_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_BDM1_2D, BF_N_T_BDM1_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_BDM2_2D, BF_N_T_BDM2_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_BDM3_2D, BF_N_T_BDM3_2D_Obj);

  RegisterBaseFunct2D(BF_N_T_P1_2D, BF_N_T_P1_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_Q1_2D, BF_N_Q_Q1_2D_Obj);
  
  RegisterBaseFunct2D(BF_D_Q_P1_2D, BF_D_Q_P1_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_P2_2D, BF_D_Q_P2_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_P3_2D, BF_D_Q_P3_2D_Obj);

  RegisterBaseFunct2D(BF_C_T_B2_2D, BF_C_T_B2_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_B3_2D, BF_C_T_B3_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_B4_2D, BF_C_T_B4_2D_Obj);

  RegisterBaseFunct2D(BF_D_T_P1_2D, BF_D_T_P1_2D_Obj);
  RegisterBaseFunct2D(BF_D_T_P2_2D, BF_D_T_P2_2D_Obj);

  RegisterBaseFunct2D(BF_N_Q_Q2_2D, BF_N_Q_Q2_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_Q3_2D, BF_N_Q_Q3_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_Q4_2D, BF_N_Q_Q4_2D_Obj);
  RegisterBaseFunct2D(BF_N_Q_Q5_2D, BF_N_Q_Q5_2D_Obj);

  RegisterBaseFunct2D(BF_D_Q_P4_2D, BF_D_Q_P4_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_P5_2D, BF_D_Q_P5_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_P6_2D, BF_D_Q_P6_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_P7_2D, BF_D_Q_P7_2D_Obj);

  RegisterBaseFunct2D(BF_N_T_P1MOD_2D, BF_N_T_P1MOD_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_P1MINI_2D, BF_C_T_P1MINI_2D_Obj);

  RegisterBaseFunct2D(BF_N_T_P2_2D, BF_N_T_P2_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_P3_2D, BF_N_T_P3_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_P4_2D, BF_N_T_P4_2D_Obj);
  RegisterBaseFunct2D(BF_N_T_P5_2D, BF_N_T_P5_2D_Obj);

  RegisterBaseFunct2D(BF_D_T_P3_2D, BF_D_T_P3_2D_Obj);
  RegisterBaseFunct2D(BF_D_T_P4_2D, BF_D_T_P4_2D_Obj);

  RegisterBaseFunct2D(BF_B_Q_IB2_2D, BF_B_Q_IB2_2D_Obj);

  RegisterBaseFunct2D(BF_D_Q_Q1_2D, BF_D_Q_Q1_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_Q2_2D, BF_D_Q_Q2_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_Q3_2D, BF_D_Q_Q3_2D_Obj);
  RegisterBaseFunct2D(BF_D_Q_Q4_2D, BF_D_Q_Q4_2D_Obj);

  RegisterBaseFunct2D(BF_D_Q_D2_2D, BF_D_Q_D2_2D_Obj);

  RegisterBaseFunct2D(BF_C_Q_UL1_2D, BF_C_Q_UL1_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL2_2D, BF_C_Q_UL2_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL3_2D, BF_C_Q_UL3_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL4_2D, BF_C_Q_UL4_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL5_2D, BF_C_Q_UL5_2D_Obj);
  
  RegisterBaseFunct2D(BF_C_T_UL1_2D, BF_C_T_UL1_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_UL2_2D, BF_C_T_UL2_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_UL3_2D, BF_C_T_UL3_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_UL4_2D, BF_C_T_UL4_2D_Obj);
  RegisterBaseFunct2D(BF_C_T_UL5_2D, BF_C_T_UL5_2D_Obj);
  
  RegisterBaseFunct2D(BF_C_Q_UL2S_2D, BF_C_Q_UL2S_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL3S_2D, BF_C_Q_UL3S_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL4S_2D, BF_C_Q_UL4S_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL5S_2D, BF_C_Q_UL5S_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL6S_2D, BF_C_Q_UL6S_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL7S_2D, BF_C_Q_UL7S_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL8S_2D, BF_C_Q_UL8S_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL9S_2D, BF_C_Q_UL9S_2D_Obj);

  RegisterBaseFunct2D(BF_C_Q_UL2SE_2D, BF_C_Q_UL2SE_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL3SE_2D, BF_C_Q_UL3SE_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL4SE_2D, BF_C_Q_UL4SE_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL5SE_2D, BF_C_Q_UL5SE_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL6SE_2D, BF_C_Q_UL6SE_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL7SE_2D, BF_C_Q_UL7SE_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL8SE_2D, BF_C_Q_UL8SE_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_UL9SE_2D, BF_C_Q_UL9SE_2D_Obj);

  RegisterBaseFunct2D(BF_C_Q_M2_2D, BF_C_Q_M2_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_M3_2D, BF_C_Q_M3_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_M4_2D, BF_C_Q_M4_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_M5_2D, BF_C_Q_M5_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_M6_2D, BF_C_Q_M6_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_M7_2D, BF_C_Q_M7_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_M8_2D, BF_C_Q_M8_2D_Obj);
  RegisterBaseFunct2D(BF_C_Q_M9_2D, BF_C_Q_M9_2D_Obj);

  RegisterBaseFunct2D(BF_C_Q_EL1_2D, BF_C_Q_EL1_2D_Obj);  
  
  RegisterBaseFunct2D(BF_C_T_SV2_2D, BF_C_T_SV2_2D_Obj);
  RegisterBaseFunct2D(BF_D_T_SV1_2D, BF_D_T_SV1_2D_Obj);
  
#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "base functions registered" << endl;
}

void TFEDatabase2D::RegisterAllNodalFunctionals()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  RegisterNodalFunctional1D(NF_C_L_P0_1D, NF_C_L_P0_1D_Obj);
  RegisterNodalFunctional1D(NF_C_L_P1_1D, NF_C_L_P1_1D_Obj);
  RegisterNodalFunctional1D(NF_C_L_P2_1D, NF_C_L_P2_1D_Obj);
  RegisterNodalFunctional1D(NF_C_L_P3_1D, NF_C_L_P3_1D_Obj);

  RegisterNodalFunctional1D(NF_D_L_P1_1D, NF_D_L_P1_1D_Obj);
  RegisterNodalFunctional1D(NF_D_L_P2_1D, NF_D_L_P2_1D_Obj); 

  RegisterNodalFunctional2D(NF_C_T_P00_2D, NF_C_T_P00_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P0_2D, NF_C_T_P0_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P1_2D, NF_C_T_P1_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P2_2D, NF_C_T_P2_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P3_2D, NF_C_T_P3_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P4_2D, NF_C_T_P4_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P5_2D, NF_C_T_P5_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P6_2D, NF_C_T_P6_2D_Obj);
  
  RegisterNodalFunctional2D(NF_C_Q_Q00_2D, NF_C_Q_Q00_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q0_2D, NF_C_Q_Q0_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q1_2D, NF_C_Q_Q1_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q2_2D, NF_C_Q_Q2_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q3_2D, NF_C_Q_Q3_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q4_2D, NF_C_Q_Q4_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q5_2D, NF_C_Q_Q5_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q6_2D, NF_C_Q_Q6_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q7_2D, NF_C_Q_Q7_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q8_2D, NF_C_Q_Q8_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_Q9_2D, NF_C_Q_Q9_2D_Obj);
  
  RegisterNodalFunctional2D(NF_N_Q_RT0_2D, NF_N_Q_RT0_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_RT1_2D, NF_N_Q_RT1_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_RT2_2D, NF_N_Q_RT2_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_RT3_2D, NF_N_Q_RT3_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_RT0_2D, NF_N_T_RT0_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_RT1_2D, NF_N_T_RT1_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_RT2_2D, NF_N_T_RT2_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_RT3_2D, NF_N_T_RT3_2D_Obj);

  RegisterNodalFunctional2D(NF_N_Q_BDM1_2D, NF_N_Q_BDM1_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_BDM2_2D, NF_N_Q_BDM2_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_BDM3_2D, NF_N_Q_BDM3_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_BDM1_2D, NF_N_T_BDM1_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_BDM2_2D, NF_N_T_BDM2_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_BDM3_2D, NF_N_T_BDM3_2D_Obj);

  RegisterNodalFunctional2D(NF_N_T_P1_2D, NF_N_T_P1_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_Q1_2D, NF_N_Q_Q1_2D_Obj);
  
  RegisterNodalFunctional2D(NF_D_Q_P1_2D, NF_D_Q_P1_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_P2_2D, NF_D_Q_P2_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_P3_2D, NF_D_Q_P3_2D_Obj);

  RegisterNodalFunctional2D(NF_C_T_B2_2D, NF_C_T_B2_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_B3_2D, NF_C_T_B3_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_B4_2D, NF_C_T_B4_2D_Obj);

  RegisterNodalFunctional2D(NF_D_T_P1_2D, NF_D_T_P1_2D_Obj);
  RegisterNodalFunctional2D(NF_D_T_P2_2D, NF_D_T_P2_2D_Obj);

  RegisterNodalFunctional2D(NF_N_Q_Q2_2D, NF_N_Q_Q2_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_Q3_2D, NF_N_Q_Q3_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_Q4_2D, NF_N_Q_Q4_2D_Obj);
  RegisterNodalFunctional2D(NF_N_Q_Q5_2D, NF_N_Q_Q5_2D_Obj);

  RegisterNodalFunctional2D(NF_D_Q_P4_2D, NF_D_Q_P4_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_P5_2D, NF_D_Q_P5_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_P6_2D, NF_D_Q_P6_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_P7_2D, NF_D_Q_P7_2D_Obj);

  RegisterNodalFunctional2D(NF_N_T_P1MOD_2D, NF_N_T_P1MOD_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_P1MINI_2D, NF_C_T_P1MINI_2D_Obj);

  RegisterNodalFunctional2D(NF_N_T_P2_2D, NF_N_T_P2_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_P3_2D, NF_N_T_P3_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_P4_2D, NF_N_T_P4_2D_Obj);
  RegisterNodalFunctional2D(NF_N_T_P5_2D, NF_N_T_P5_2D_Obj);

  RegisterNodalFunctional2D(NF_D_T_P3_2D, NF_D_T_P3_2D_Obj);
  RegisterNodalFunctional2D(NF_D_T_P4_2D, NF_D_T_P4_2D_Obj);

  RegisterNodalFunctional2D(NF_B_Q_IB2_2D, NF_B_Q_IB2_2D_Obj);

  RegisterNodalFunctional2D(NF_S_Q_Q2_2D, NF_S_Q_Q2_2D_Obj);

  RegisterNodalFunctional2D(NF_D_Q_Q1_2D, NF_D_Q_Q1_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_Q2_2D, NF_D_Q_Q2_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_Q3_2D, NF_D_Q_Q3_2D_Obj);
  RegisterNodalFunctional2D(NF_D_Q_Q4_2D, NF_D_Q_Q4_2D_Obj);

  RegisterNodalFunctional2D(NF_D_Q_D2_2D, NF_D_Q_D2_2D_Obj);
  
  RegisterNodalFunctional2D(NF_C_Q_UL1_2D, NF_C_Q_UL1_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL2_2D, NF_C_Q_UL2_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL3_2D, NF_C_Q_UL3_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL4_2D, NF_C_Q_UL4_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL5_2D, NF_C_Q_UL5_2D_Obj);

  RegisterNodalFunctional2D(NF_C_T_UL1_2D, NF_C_T_UL1_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_UL2_2D, NF_C_T_UL2_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_UL3_2D, NF_C_T_UL3_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_UL4_2D, NF_C_T_UL4_2D_Obj);
  RegisterNodalFunctional2D(NF_C_T_UL5_2D, NF_C_T_UL5_2D_Obj);

  RegisterNodalFunctional2D(NF_C_Q_UL2S_2D, NF_C_Q_UL2S_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL3S_2D, NF_C_Q_UL3S_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL4S_2D, NF_C_Q_UL4S_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL5S_2D, NF_C_Q_UL5S_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL6S_2D, NF_C_Q_UL6S_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL7S_2D, NF_C_Q_UL7S_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL8S_2D, NF_C_Q_UL8S_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL9S_2D, NF_C_Q_UL9S_2D_Obj);

  RegisterNodalFunctional2D(NF_C_Q_UL2SE_2D, NF_C_Q_UL2SE_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL3SE_2D, NF_C_Q_UL3SE_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL4SE_2D, NF_C_Q_UL4SE_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL5SE_2D, NF_C_Q_UL5SE_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL6SE_2D, NF_C_Q_UL6SE_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL7SE_2D, NF_C_Q_UL7SE_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL8SE_2D, NF_C_Q_UL8SE_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_UL9SE_2D, NF_C_Q_UL9SE_2D_Obj);

  RegisterNodalFunctional2D(NF_C_Q_M2_2D, NF_C_Q_M2_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_M3_2D, NF_C_Q_M3_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_M4_2D, NF_C_Q_M4_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_M5_2D, NF_C_Q_M5_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_M6_2D, NF_C_Q_M6_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_M7_2D, NF_C_Q_M7_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_M8_2D, NF_C_Q_M8_2D_Obj);
  RegisterNodalFunctional2D(NF_C_Q_M9_2D, NF_C_Q_M9_2D_Obj);
 
  RegisterNodalFunctional2D(NF_C_Q_EL1_2D, NF_C_Q_EL1_2D_Obj);
  
  RegisterNodalFunctional2D(NF_C_T_SV2_2D, NF_C_T_SV2_2D_Obj);
  RegisterNodalFunctional2D(NF_D_T_SV1_2D, NF_D_T_SV1_2D_Obj);
  
#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "nodal functionals registered" << endl;
}

void TFEDatabase2D::RegisterAllFEs()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  TFE1D *ele1D;
  TFE2D *ele2D;

  // ======================================================================
  // elements on lines
  // ======================================================================
  ele1D = new TFE1D(BF_C_L_P0_1D, NF_C_L_P0_1D, LineAffin, FE_C_L_P0_1D, 0);
  RegisterFE1D(C_P0_1D_L_A, ele1D);
  ele1D = new TFE1D(BF_C_L_P1_1D, NF_C_L_P1_1D, LineAffin, FE_C_L_P1_1D, 0);
  RegisterFE1D(C_P1_1D_L_A, ele1D);
  ele1D = new TFE1D(BF_C_L_P2_1D, NF_C_L_P2_1D, LineAffin, FE_C_L_P2_1D, 0);
  RegisterFE1D(C_P2_1D_L_A, ele1D);
  ele1D = new TFE1D(BF_C_L_P3_1D, NF_C_L_P3_1D, LineAffin, FE_C_L_P3_1D, 0);
  RegisterFE1D(C_P3_1D_L_A, ele1D);

  ele1D = new TFE1D(BF_C_L_P0_1D, NF_C_L_P0_1D, LineAffin, FE_N_L_P0_1D, 0);
  RegisterFE1D(N_P0_1D_L_A, ele1D);

  ele1D = new TFE1D(BF_D_L_P1_1D, NF_D_L_P1_1D, LineAffin, FE_D_L_P1_1D, 0);
  RegisterFE1D(D_P1_1D_L_A, ele1D);
  ele1D = new TFE1D(BF_D_L_P2_1D, NF_D_L_P2_1D, LineAffin, FE_D_L_P2_1D, 0);
  RegisterFE1D(D_P2_1D_L_A, ele1D);  
  // ======================================================================
  // elements on triangles
  // ======================================================================
  ele2D = new TFE2D(BF_C_T_P00_2D, NF_C_T_P00_2D, TriaAffin, FE_C_T_P00_2D, 0);
  RegisterFE2D(C_P00_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_P0_2D, NF_C_T_P0_2D, TriaAffin, FE_C_T_P0_2D, 0);
  RegisterFE2D(C_P0_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_P1_2D, NF_C_T_P1_2D, TriaAffin, FE_C_T_P1_2D, 0);
  RegisterFE2D(C_P1_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_P2_2D, NF_C_T_P2_2D, TriaAffin, FE_C_T_P2_2D, 0);
  RegisterFE2D(C_P2_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_P3_2D, NF_C_T_P3_2D, TriaAffin, FE_C_T_P3_2D, 0);
  RegisterFE2D(C_P3_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_P4_2D, NF_C_T_P4_2D, TriaAffin, FE_C_T_P4_2D, 0);
  RegisterFE2D(C_P4_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_P5_2D, NF_C_T_P5_2D, TriaAffin, FE_C_T_P5_2D, 0);
  RegisterFE2D(C_P5_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_P6_2D, NF_C_T_P6_2D, TriaAffin, FE_C_T_P6_2D, 0);
  RegisterFE2D(C_P6_2D_T_A, ele2D);
  
  ele2D = new TFE2D(BF_N_T_P1_2D, NF_N_T_P1_2D, TriaAffin, FE_N_T_P1_2D, 0);
  RegisterFE2D(N_P1_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_B2_2D, NF_C_T_B2_2D, TriaAffin, FE_C_T_B2_2D, 0);
  RegisterFE2D(C_B2_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_SV2_2D, NF_C_T_SV2_2D, TriaAffin, FE_C_T_SV2_2D, 0);
  RegisterFE2D(C_SV2_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_B3_2D, NF_C_T_B3_2D, TriaAffin, FE_C_T_B3_2D, 0);
  RegisterFE2D(C_B3_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_B4_2D, NF_C_T_B4_2D, TriaAffin, FE_C_T_B4_2D, 0);
  RegisterFE2D(C_B4_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_D_T_P1_2D, NF_D_T_P1_2D, TriaAffin, FE_D_T_P1_2D, 0);
  RegisterFE2D(D_P1_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_D_T_SV1_2D, NF_D_T_SV1_2D, TriaAffin, FE_D_T_SV1_2D, 0);
  RegisterFE2D(D_SV1_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_D_T_P2_2D, NF_D_T_P2_2D, TriaAffin, FE_D_T_P2_2D, 0);
  RegisterFE2D(D_P2_2D_T_A, ele2D);

  ele2D = new TFE2D(BF_N_T_P1MOD_2D, NF_N_T_P1MOD_2D, TriaAffin,
                    FE_N_T_P1MOD_2D, 0);
  RegisterFE2D(N_P1MOD_2D_T_A, ele2D);

  ele2D = new TFE2D(BF_C_T_P1MINI_2D, NF_C_T_P1MINI_2D, TriaAffin,
                    FE_C_T_P1MINI_2D, 0);
  RegisterFE2D(C_P1MINI_2D_T_A, ele2D);

  ele2D = new TFE2D(BF_N_T_P2_2D, NF_N_T_P2_2D, TriaAffin, FE_N_T_P2_2D, 0);
  RegisterFE2D(N_P2_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_P3_2D, NF_N_T_P3_2D, TriaAffin, FE_N_T_P3_2D, 0);
  RegisterFE2D(N_P3_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_P4_2D, NF_N_T_P4_2D, TriaAffin, FE_N_T_P4_2D, 0);
  RegisterFE2D(N_P4_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_P5_2D, NF_N_T_P5_2D, TriaAffin, FE_N_T_P5_2D, 0);
  RegisterFE2D(N_P5_2D_T_A, ele2D);

  ele2D = new TFE2D(BF_D_T_P3_2D, NF_D_T_P3_2D, TriaAffin, FE_D_T_P3_2D, 0);
  RegisterFE2D(D_P3_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_D_T_P4_2D, NF_D_T_P4_2D, TriaAffin, FE_D_T_P4_2D, 0);
  RegisterFE2D(D_P4_2D_T_A, ele2D);
  
  ele2D = new TFE2D(BF_N_T_RT0_2D, NF_N_T_RT0_2D, TriaAffin, FE_N_T_RT0_2D, 0);
  RegisterFE2D(N_RT0_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_RT1_2D, NF_N_T_RT1_2D, TriaAffin, FE_N_T_RT1_2D, 0);
  RegisterFE2D(N_RT1_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_RT2_2D, NF_N_T_RT2_2D, TriaAffin, FE_N_T_RT2_2D, 0);
  RegisterFE2D(N_RT2_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_RT3_2D, NF_N_T_RT3_2D, TriaAffin, FE_N_T_RT3_2D, 0);
  RegisterFE2D(N_RT3_2D_T_A, ele2D);

  ele2D = new TFE2D(BF_N_T_BDM1_2D, NF_N_T_BDM1_2D, TriaAffin, FE_N_T_BDM1_2D, 0);
  RegisterFE2D(N_BDM1_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_BDM2_2D, NF_N_T_BDM2_2D, TriaAffin, FE_N_T_BDM2_2D, 0);
  RegisterFE2D(N_BDM2_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_N_T_BDM3_2D, NF_N_T_BDM3_2D, TriaAffin, FE_N_T_BDM3_2D, 0);
  RegisterFE2D(N_BDM3_2D_T_A, ele2D);

  //========LOCALPROJECTION==============
  ele2D = new TFE2D(BF_C_T_UL1_2D,NF_C_T_UL1_2D, TriaAffin, FE_C_T_UL1_2D, 0);
  RegisterFE2D(C_UL1_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_UL2_2D,NF_C_T_UL2_2D, TriaAffin, FE_C_T_UL2_2D, 0);
  RegisterFE2D(C_UL2_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_UL3_2D,NF_C_T_UL3_2D, TriaAffin, FE_C_T_UL3_2D, 0);
  RegisterFE2D(C_UL3_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_UL4_2D,NF_C_T_UL4_2D, TriaAffin, FE_C_T_UL4_2D, 0);
  RegisterFE2D(C_UL4_2D_T_A, ele2D);
  ele2D = new TFE2D(BF_C_T_UL5_2D,NF_C_T_UL5_2D, TriaAffin, FE_C_T_UL5_2D, 0);
  RegisterFE2D(C_UL5_2D_T_A, ele2D);
  // ele2D->CheckNFandBF();
  //=====================================

  // ======================================================================
  // elements on parallelograms
  // ======================================================================
  ele2D = new TFE2D(BF_C_Q_Q00_2D, NF_C_Q_Q00_2D, QuadAffin, FE_C_Q_Q00_2D, 0);
  RegisterFE2D(C_Q00_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q0_2D, NF_C_Q_Q0_2D, QuadAffin, FE_C_Q_Q0_2D, 0);
  RegisterFE2D(C_Q0_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q1_2D, NF_C_Q_Q1_2D, QuadAffin, FE_C_Q_Q1_2D, 0);
  RegisterFE2D(C_Q1_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q2_2D, NF_C_Q_Q2_2D, QuadAffin, FE_C_Q_Q2_2D, 0);
  RegisterFE2D(C_Q2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q3_2D, NF_C_Q_Q3_2D, QuadAffin, FE_C_Q_Q3_2D, 0);
  RegisterFE2D(C_Q3_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q4_2D, NF_C_Q_Q4_2D, QuadAffin, FE_C_Q_Q4_2D, 0);
  RegisterFE2D(C_Q4_2D_Q_A, ele2D);
  // ele2D->CheckNFandBF();
  ele2D = new TFE2D(BF_C_Q_Q5_2D, NF_C_Q_Q5_2D, QuadAffin, FE_C_Q_Q5_2D, 0);
  RegisterFE2D(C_Q5_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q6_2D, NF_C_Q_Q6_2D, QuadAffin, FE_C_Q_Q6_2D, 0);
  RegisterFE2D(C_Q6_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q7_2D, NF_C_Q_Q7_2D, QuadAffin, FE_C_Q_Q7_2D, 0);
  RegisterFE2D(C_Q7_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q8_2D, NF_C_Q_Q8_2D, QuadAffin, FE_C_Q_Q8_2D, 0);
  RegisterFE2D(C_Q8_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q9_2D, NF_C_Q_Q9_2D, QuadAffin, FE_C_Q_Q9_2D, 0);
  RegisterFE2D(C_Q9_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_N_Q_RT0_2D, NF_N_Q_RT0_2D, QuadAffin, FE_N_Q_RT0_2D, 0);
  RegisterFE2D(N_RT0_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_RT1_2D, NF_N_Q_RT1_2D, QuadAffin, FE_N_Q_RT1_2D, 0);
  RegisterFE2D(N_RT1_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_RT2_2D, NF_N_Q_RT2_2D, QuadAffin, FE_N_Q_RT2_2D, 0);
  RegisterFE2D(N_RT2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_RT3_2D, NF_N_Q_RT3_2D, QuadAffin, FE_N_Q_RT3_2D, 0);
  RegisterFE2D(N_RT3_2D_Q_A, ele2D);
  
  ele2D = new TFE2D(BF_N_Q_BDM1_2D, NF_N_Q_BDM1_2D, QuadAffin, FE_N_Q_BDM1_2D, 0);
  RegisterFE2D(N_BDM1_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_BDM2_2D, NF_N_Q_BDM2_2D, QuadAffin, FE_N_Q_BDM2_2D, 0);
  RegisterFE2D(N_BDM2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_BDM3_2D, NF_N_Q_BDM3_2D, QuadAffin, FE_N_Q_BDM3_2D, 0);
  RegisterFE2D(N_BDM3_2D_Q_A, ele2D);
  
  ele2D = new TFE2D(BF_N_Q_Q1_2D, NF_N_Q_Q1_2D, QuadAffin, FE_N_Q_Q1_2D, 0);
  RegisterFE2D(N_Q1_2D_Q_A, ele2D);
  
  ele2D = new TFE2D(BF_D_Q_P1_2D, NF_D_Q_P1_2D, QuadAffin, FE_D_Q_P1_2D, 0);
  RegisterFE2D(D_P1_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_P2_2D, NF_D_Q_P2_2D, QuadAffin, FE_D_Q_P2_2D, 0);
  RegisterFE2D(D_P2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_P3_2D, NF_D_Q_P3_2D, QuadAffin, FE_D_Q_P3_2D, 0);
  RegisterFE2D(D_P3_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_N_Q_Q2_2D, NF_N_Q_Q2_2D, QuadAffin, FE_N_Q_Q2_2D, 0);
  RegisterFE2D(N_Q2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_Q3_2D, NF_N_Q_Q3_2D, QuadAffin, FE_N_Q_Q3_2D, 0);
  RegisterFE2D(N_Q3_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_Q4_2D, NF_N_Q_Q4_2D, QuadAffin, FE_N_Q_Q4_2D, 0);
  RegisterFE2D(N_Q4_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_N_Q_Q5_2D, NF_N_Q_Q5_2D, QuadAffin, FE_N_Q_Q5_2D, 0);
  RegisterFE2D(N_Q5_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_D_Q_P4_2D, NF_D_Q_P4_2D, QuadAffin, FE_D_Q_P4_2D, 0);
  RegisterFE2D(D_P4_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_P5_2D, NF_D_Q_P5_2D, QuadAffin, FE_D_Q_P5_2D, 0);
  RegisterFE2D(D_P5_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_P6_2D, NF_D_Q_P6_2D, QuadAffin, FE_D_Q_P6_2D, 0);
  RegisterFE2D(D_P6_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_P7_2D, NF_D_Q_P7_2D, QuadAffin, FE_D_Q_P7_2D, 0);
  RegisterFE2D(D_P7_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_B_Q_IB2_2D, NF_B_Q_IB2_2D, QuadAffin, FE_B_Q_IB2_2D, 0);
  RegisterFE2D(B_IB2_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_D_Q_Q1_2D, NF_D_Q_Q1_2D, QuadAffin, FE_D_Q_Q1_2D, 0);
  RegisterFE2D(D_Q1_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_Q2_2D, NF_D_Q_Q2_2D, QuadAffin, FE_D_Q_Q2_2D, 0);
  RegisterFE2D(D_Q2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_Q3_2D, NF_D_Q_Q3_2D, QuadAffin, FE_D_Q_Q3_2D, 0);
  RegisterFE2D(D_Q3_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_D_Q_Q4_2D, NF_D_Q_Q4_2D, QuadAffin, FE_D_Q_Q4_2D, 0);
  RegisterFE2D(D_Q4_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_D_Q_D2_2D, NF_D_Q_D2_2D, QuadAffin, FE_D_Q_D2_2D, 0);
  RegisterFE2D(D_D2_2D_Q_A, ele2D);
  
    //========LOCALPROJECTION==============
  ele2D = new TFE2D(BF_C_Q_UL1_2D,NF_C_Q_UL1_2D, QuadAffin, FE_C_Q_UL1_2D, 0);
  RegisterFE2D(C_UL1_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL2_2D,NF_C_Q_UL2_2D, QuadAffin, FE_C_Q_UL2_2D, 0);
  RegisterFE2D(C_UL2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL3_2D,NF_C_Q_UL3_2D, QuadAffin, FE_C_Q_UL3_2D, 0);
  RegisterFE2D(C_UL3_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL4_2D,NF_C_Q_UL4_2D, QuadAffin, FE_C_Q_UL4_2D, 0);
  RegisterFE2D(C_UL4_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL5_2D,NF_C_Q_UL5_2D, QuadAffin, FE_C_Q_UL5_2D, 0);
  RegisterFE2D(C_UL5_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_C_Q_UL2S_2D,NF_C_Q_UL2S_2D, QuadAffin, FE_C_Q_UL2S_2D, 0);
  RegisterFE2D(C_UL2S_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL3S_2D,NF_C_Q_UL3S_2D, QuadAffin, FE_C_Q_UL3S_2D, 0);
  RegisterFE2D(C_UL3S_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL4S_2D,NF_C_Q_UL4S_2D, QuadAffin, FE_C_Q_UL4S_2D, 0);
  RegisterFE2D(C_UL4S_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL5S_2D,NF_C_Q_UL5S_2D, QuadAffin , FE_C_Q_UL5S_2D, 0);
  RegisterFE2D(C_UL5S_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL6S_2D,NF_C_Q_UL6S_2D, QuadAffin, FE_C_Q_UL6S_2D, 0);
  RegisterFE2D(C_UL6S_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL7S_2D,NF_C_Q_UL7S_2D, QuadAffin, FE_C_Q_UL7S_2D, 0);
  RegisterFE2D(C_UL7S_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL8S_2D,NF_C_Q_UL8S_2D, QuadAffin, FE_C_Q_UL8S_2D, 0);
  RegisterFE2D(C_UL8S_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL9S_2D,NF_C_Q_UL9S_2D, QuadAffin , FE_C_Q_UL9S_2D, 0);
  RegisterFE2D(C_UL9S_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_C_Q_UL2SE_2D,NF_C_Q_UL2SE_2D, QuadAffin, FE_C_Q_UL2SE_2D, 0);
  RegisterFE2D(C_UL2SE_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL3SE_2D,NF_C_Q_UL3SE_2D, QuadAffin, FE_C_Q_UL3SE_2D, 0);
  RegisterFE2D(C_UL3SE_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL4SE_2D,NF_C_Q_UL4SE_2D, QuadAffin, FE_C_Q_UL4SE_2D, 0);
  RegisterFE2D(C_UL4SE_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL5SE_2D,NF_C_Q_UL5SE_2D, QuadAffin , FE_C_Q_UL5SE_2D, 0);
  RegisterFE2D(C_UL5SE_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL6SE_2D,NF_C_Q_UL6SE_2D, QuadAffin, FE_C_Q_UL6SE_2D, 0);
  RegisterFE2D(C_UL6SE_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL7SE_2D,NF_C_Q_UL7SE_2D, QuadAffin, FE_C_Q_UL7SE_2D, 0);
  RegisterFE2D(C_UL7SE_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL8SE_2D,NF_C_Q_UL8SE_2D, QuadAffin, FE_C_Q_UL8SE_2D, 0);
  RegisterFE2D(C_UL8SE_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL9SE_2D,NF_C_Q_UL9SE_2D, QuadAffin , FE_C_Q_UL9SE_2D, 0);
  RegisterFE2D(C_UL9SE_2D_Q_A, ele2D);

  ele2D = new TFE2D(BF_C_Q_M2_2D,NF_C_Q_M2_2D, QuadAffin, FE_C_Q_M2_2D, 0);
  RegisterFE2D(C_M2_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_M3_2D,NF_C_Q_M3_2D, QuadAffin, FE_C_Q_M3_2D, 0);
  RegisterFE2D(C_M3_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_M4_2D,NF_C_Q_M4_2D, QuadAffin, FE_C_Q_M4_2D, 0);
  RegisterFE2D(C_M4_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_M5_2D,NF_C_Q_M5_2D, QuadAffin , FE_C_Q_M5_2D, 0);
  RegisterFE2D(C_M5_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_M6_2D,NF_C_Q_M6_2D, QuadAffin, FE_C_Q_M6_2D, 0);
  RegisterFE2D(C_M6_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_M7_2D,NF_C_Q_M7_2D, QuadAffin, FE_C_Q_M7_2D, 0);
  RegisterFE2D(C_M7_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_M8_2D,NF_C_Q_M8_2D, QuadAffin, FE_C_Q_M8_2D, 0);
  RegisterFE2D(C_M8_2D_Q_A, ele2D);
  ele2D = new TFE2D(BF_C_Q_M9_2D,NF_C_Q_M9_2D, QuadAffin , FE_C_Q_M9_2D, 0);
  RegisterFE2D(C_M9_2D_Q_A, ele2D);
  //=====================================

  // ======================================================================
  // elements on arbitrary quadrangles
  // ======================================================================
  ele2D = new TFE2D(BF_C_Q_Q00_2D, NF_C_Q_Q00_2D, QuadBilinear, FE_C_Q_Q00_2D, 0);
  RegisterFE2D(C_Q00_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q0_2D, NF_C_Q_Q0_2D, QuadBilinear, FE_C_Q_Q0_2D, 0);
  RegisterFE2D(C_Q0_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q1_2D, NF_C_Q_Q1_2D, QuadBilinear, FE_C_Q_Q1_2D, 0);
  RegisterFE2D(C_Q1_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q2_2D, NF_C_Q_Q2_2D, QuadBilinear, FE_C_Q_Q2_2D, 0);
  RegisterFE2D(C_Q2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q3_2D, NF_C_Q_Q3_2D, QuadBilinear, FE_C_Q_Q3_2D, 0);
  RegisterFE2D(C_Q3_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q4_2D, NF_C_Q_Q4_2D, QuadBilinear, FE_C_Q_Q4_2D, 0);
  RegisterFE2D(C_Q4_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q5_2D, NF_C_Q_Q5_2D, QuadBilinear, FE_C_Q_Q5_2D, 0);
  RegisterFE2D(C_Q5_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q6_2D, NF_C_Q_Q6_2D, QuadBilinear, FE_C_Q_Q6_2D, 0);
  RegisterFE2D(C_Q6_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q7_2D, NF_C_Q_Q7_2D, QuadBilinear, FE_C_Q_Q7_2D, 0);
  RegisterFE2D(C_Q7_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q8_2D, NF_C_Q_Q8_2D, QuadBilinear, FE_C_Q_Q8_2D, 0);
  RegisterFE2D(C_Q8_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_Q9_2D, NF_C_Q_Q9_2D, QuadBilinear, FE_C_Q_Q9_2D, 0);
  RegisterFE2D(C_Q9_2D_Q_M, ele2D);
  
  // Raviart-Thomas elements (RT)
  ele2D = new TFE2D(BF_N_Q_RT0_2D, NF_N_Q_RT0_2D, QuadBilinear, FE_N_Q_RT0_2D, 0);
  RegisterFE2D(N_RT0_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_RT1_2D, NF_N_Q_RT1_2D, QuadBilinear, FE_N_Q_RT1_2D, 0);
  RegisterFE2D(N_RT1_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_RT2_2D, NF_N_Q_RT2_2D, QuadBilinear, FE_N_Q_RT2_2D, 0);
  RegisterFE2D(N_RT2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_RT3_2D, NF_N_Q_RT3_2D, QuadBilinear, FE_N_Q_RT3_2D, 0);
  RegisterFE2D(N_RT3_2D_Q_M, ele2D);
  // Brezzi-Douglas-Marini elements (BDM)
  ele2D = new TFE2D(BF_N_Q_BDM1_2D, NF_N_Q_BDM1_2D, QuadBilinear, FE_N_Q_BDM1_2D, 0);
  RegisterFE2D(N_BDM1_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_BDM2_2D, NF_N_Q_BDM2_2D, QuadBilinear, FE_N_Q_BDM2_2D, 0);
  RegisterFE2D(N_BDM2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_BDM3_2D, NF_N_Q_BDM3_2D, QuadBilinear, FE_N_Q_BDM3_2D, 0);
  RegisterFE2D(N_BDM3_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_N_Q_Q1_2D, NF_N_Q_Q1_2D, QuadBilinear, FE_N_Q_Q1_2D, 0);
  RegisterFE2D(N_Q1_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_D_Q_P1_2D, NF_D_Q_P1_2D, QuadBilinear, FE_D_Q_P1_2D, 0);
  RegisterFE2D(D_P1_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_D_Q_P2_2D, NF_D_Q_P2_2D, QuadBilinear, FE_D_Q_P2_2D, 0);
  RegisterFE2D(D_P2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_D_Q_P3_2D, NF_D_Q_P3_2D, QuadBilinear, FE_D_Q_P3_2D, 0);
  RegisterFE2D(D_P3_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_N_Q_Q2_2D, NF_N_Q_Q2_2D, QuadBilinear, FE_N_Q_Q2_2D, 0);
  RegisterFE2D(N_Q2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_Q3_2D, NF_N_Q_Q3_2D, QuadBilinear, FE_N_Q_Q3_2D, 0);
  RegisterFE2D(N_Q3_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_Q4_2D, NF_N_Q_Q4_2D, QuadBilinear, FE_N_Q_Q4_2D, 0);
  RegisterFE2D(N_Q4_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_N_Q_Q5_2D, NF_N_Q_Q5_2D, QuadBilinear, FE_N_Q_Q5_2D, 0);
  RegisterFE2D(N_Q5_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_D_Q_P4_2D, NF_D_Q_P4_2D, QuadBilinear, FE_D_Q_P4_2D, 0);
  RegisterFE2D(D_P4_2D_Q_M, ele2D);
  // ele2D->CheckNFandBF();
  ele2D = new TFE2D(BF_D_Q_P5_2D, NF_D_Q_P5_2D, QuadBilinear, FE_D_Q_P5_2D, 0);
  RegisterFE2D(D_P5_2D_Q_M, ele2D);
  // ele2D->CheckNFandBF();
  ele2D = new TFE2D(BF_D_Q_P6_2D, NF_D_Q_P6_2D, QuadBilinear, FE_D_Q_P6_2D, 0);
  RegisterFE2D(D_P6_2D_Q_M, ele2D);
  // ele2D->CheckNFandBF();
  ele2D = new TFE2D(BF_D_Q_P7_2D, NF_D_Q_P7_2D, QuadBilinear, FE_D_Q_P7_2D, 0);
  RegisterFE2D(D_P7_2D_Q_M, ele2D);
  // ele2D->CheckNFandBF();

  ele2D = new TFE2D(BF_B_Q_IB2_2D, NF_B_Q_IB2_2D, QuadBilinear, FE_B_Q_IB2_2D, 0);
  RegisterFE2D(B_IB2_2D_Q_M, ele2D);
 
  ele2D = new TFE2D(BF_D_Q_Q1_2D, NF_D_Q_Q1_2D, QuadBilinear, FE_D_Q_Q1_2D, 0);
  RegisterFE2D(D_Q1_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_D_Q_Q2_2D, NF_D_Q_Q2_2D, QuadBilinear, FE_D_Q_Q2_2D, 0);
  RegisterFE2D(D_Q2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_D_Q_Q3_2D, NF_D_Q_Q3_2D, QuadBilinear, FE_D_Q_Q3_2D, 0);
  RegisterFE2D(D_Q3_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_D_Q_Q4_2D, NF_D_Q_Q4_2D, QuadBilinear, FE_D_Q_Q4_2D, 0);
  RegisterFE2D(D_Q4_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_D_Q_D2_2D, NF_D_Q_D2_2D, QuadBilinear, FE_D_Q_D2_2D, 0);
  RegisterFE2D(D_D2_2D_Q_M, ele2D);
  
  //========LOCALPROJECTION==============
  ele2D = new TFE2D(BF_C_Q_UL1_2D,NF_C_Q_UL1_2D, QuadBilinear, FE_C_Q_UL1_2D, 0);
  RegisterFE2D(C_UL1_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL2_2D,NF_C_Q_UL2_2D, QuadBilinear, FE_C_Q_UL2_2D, 0);
  RegisterFE2D(C_UL2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL3_2D,NF_C_Q_UL3_2D, QuadBilinear, FE_C_Q_UL3_2D, 0);
  RegisterFE2D(C_UL3_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL4_2D,NF_C_Q_UL4_2D, QuadBilinear, FE_C_Q_UL4_2D, 0);
  RegisterFE2D(C_UL4_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL5_2D,NF_C_Q_UL5_2D, QuadBilinear, FE_C_Q_UL5_2D, 0);
  RegisterFE2D(C_UL5_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_C_Q_UL2S_2D,NF_C_Q_UL2S_2D, QuadBilinear, FE_C_Q_UL2S_2D, 0);
  RegisterFE2D(C_UL2S_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL3S_2D,NF_C_Q_UL3S_2D, QuadBilinear, FE_C_Q_UL3S_2D, 0);
  RegisterFE2D(C_UL3S_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL4S_2D,NF_C_Q_UL4S_2D, QuadBilinear, FE_C_Q_UL4S_2D, 0);
  RegisterFE2D(C_UL4S_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL5S_2D,NF_C_Q_UL5S_2D, QuadBilinear, FE_C_Q_UL5S_2D, 0);
  RegisterFE2D(C_UL5S_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL6S_2D,NF_C_Q_UL6S_2D, QuadBilinear, FE_C_Q_UL6S_2D, 0);
  RegisterFE2D(C_UL6S_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL7S_2D,NF_C_Q_UL7S_2D, QuadBilinear, FE_C_Q_UL7S_2D, 0);
  RegisterFE2D(C_UL7S_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL8S_2D,NF_C_Q_UL8S_2D, QuadBilinear, FE_C_Q_UL8S_2D, 0);
  RegisterFE2D(C_UL8S_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL9S_2D,NF_C_Q_UL9S_2D, QuadBilinear, FE_C_Q_UL9S_2D, 0);
  RegisterFE2D(C_UL9S_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_C_Q_UL2SE_2D,NF_C_Q_UL2SE_2D, QuadBilinear, FE_C_Q_UL2SE_2D, 0);
  RegisterFE2D(C_UL2SE_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL3SE_2D,NF_C_Q_UL3SE_2D, QuadBilinear, FE_C_Q_UL3SE_2D, 0);
  RegisterFE2D(C_UL3SE_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL4SE_2D,NF_C_Q_UL4SE_2D, QuadBilinear, FE_C_Q_UL4SE_2D, 0);
  RegisterFE2D(C_UL4SE_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL5SE_2D,NF_C_Q_UL5SE_2D, QuadBilinear, FE_C_Q_UL5SE_2D, 0);
  RegisterFE2D(C_UL5SE_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL6SE_2D,NF_C_Q_UL6SE_2D, QuadBilinear, FE_C_Q_UL6SE_2D, 0);
  RegisterFE2D(C_UL6SE_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL7SE_2D,NF_C_Q_UL7SE_2D, QuadBilinear, FE_C_Q_UL7SE_2D, 0);
  RegisterFE2D(C_UL7SE_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL8SE_2D,NF_C_Q_UL8SE_2D, QuadBilinear, FE_C_Q_UL8SE_2D, 0);
  RegisterFE2D(C_UL8SE_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_UL9SE_2D,NF_C_Q_UL9SE_2D, QuadBilinear, FE_C_Q_UL9SE_2D, 0);
  RegisterFE2D(C_UL9SE_2D_Q_M, ele2D);

  ele2D = new TFE2D(BF_C_Q_M2_2D,NF_C_Q_M2_2D, QuadBilinear, FE_C_Q_M2_2D, 0);
  RegisterFE2D(C_M2_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_M3_2D,NF_C_Q_M3_2D, QuadBilinear, FE_C_Q_M3_2D, 0);
  RegisterFE2D(C_M3_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_M4_2D,NF_C_Q_M4_2D, QuadBilinear, FE_C_Q_M4_2D, 0);
  RegisterFE2D(C_M4_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_M5_2D,NF_C_Q_M5_2D, QuadBilinear, FE_C_Q_M5_2D, 0);
  RegisterFE2D(C_M5_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_M6_2D,NF_C_Q_M6_2D, QuadBilinear, FE_C_Q_M6_2D, 0);
  RegisterFE2D(C_M6_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_M7_2D,NF_C_Q_M7_2D, QuadBilinear, FE_C_Q_M7_2D, 0);
  RegisterFE2D(C_M7_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_M8_2D,NF_C_Q_M8_2D, QuadBilinear, FE_C_Q_M8_2D, 0);
  RegisterFE2D(C_M8_2D_Q_M, ele2D);
  ele2D = new TFE2D(BF_C_Q_M9_2D,NF_C_Q_M9_2D, QuadBilinear, FE_C_Q_M9_2D, 0);
  RegisterFE2D(C_M9_2D_Q_M, ele2D);
  //=====================================

  ele2D = new TFE2D(BF_C_Q_EL1_2D,NF_C_Q_EL1_2D, QuadAffin, FE_C_Q_EL1_2D, 0);
  RegisterFE2D(C_EL1_2D_Q_A, ele2D);
  
  ele2D = new TFE2D(BF_C_Q_EL1_2D,NF_C_Q_EL1_2D, QuadBilinear, FE_C_Q_EL1_2D, 0);
  RegisterFE2D(C_EL1_2D_Q_M, ele2D);
  
#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "finite element registered" << endl;
}

void TFEDatabase2D::RegisterAllFEMappers()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
// ========================================================================
// regular grid, same pattern on both sides
// ========================================================================
  RegisterFE2DMapper(FE_C_T_P00_2D,FE_C_T_P00_2D,C0_2_C0_2);
  RegisterFE2DMapper(FE_C_T_P00_2D,FE_C_Q_Q00_2D,C0_2_C0_2);
  RegisterFE2DMapper(FE_C_Q_Q00_2D,FE_C_T_P00_2D,C0_2_C0_2);
  RegisterFE2DMapper(FE_C_Q_Q00_2D,FE_C_Q_Q00_2D,C0_2_C0_2);

  RegisterFE2DMapper(FE_C_T_P0_2D,FE_C_T_P0_2D,C0_2_C0_2);
  RegisterFE2DMapper(FE_C_T_P0_2D,FE_C_Q_Q0_2D,C0_2_C0_2);
  RegisterFE2DMapper(FE_C_Q_Q0_2D,FE_C_T_P0_2D,C0_2_C0_2);
  RegisterFE2DMapper(FE_C_Q_Q0_2D,FE_C_Q_Q0_2D,C0_2_C0_2);

  RegisterFE2DMapper(FE_C_T_P1_2D,FE_C_T_P1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_T_P1_2D,FE_C_Q_Q1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_Q_Q1_2D,FE_C_T_P1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_Q_Q1_2D,FE_C_Q_Q1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_T_P1_2D,FE_C_T_UL1_2D,C1_2_C1_2);  
  RegisterFE2DMapper(FE_C_T_UL1_2D,FE_C_T_P1_2D,C1_2_C1_2); 
  
  RegisterFE2DMapper(FE_C_T_P2_2D,FE_C_T_P2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_P2_2D,FE_C_Q_Q2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_Q2_2D,FE_C_T_P2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_Q2_2D,FE_C_Q_Q2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_P2_2D,FE_C_T_UL2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_UL2_2D,FE_C_T_P2_2D,C2_2_C2_2);
  
  RegisterFE2DMapper(FE_C_T_B2_2D,FE_C_T_P2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_B2_2D,FE_C_Q_Q2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_P2_2D,FE_C_T_B2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_Q2_2D,FE_C_T_B2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_B2_2D,FE_C_T_B2_2D,C2_2_C2_2);
  
  RegisterFE2DMapper(FE_C_T_SV2_2D,FE_C_T_P2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_SV2_2D,FE_C_Q_Q2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_P2_2D,FE_C_T_SV2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_Q2_2D,FE_C_T_SV2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_SV2_2D,FE_C_T_SV2_2D,C2_2_C2_2);

  RegisterFE2DMapper(FE_C_T_P3_2D,FE_C_T_P3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_P3_2D,FE_C_Q_Q3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_Q3_2D,FE_C_T_P3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_Q3_2D,FE_C_Q_Q3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_P3_2D,FE_C_T_UL3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_UL3_2D,FE_C_T_P3_2D,C3_2_C3_2);  
  
  RegisterFE2DMapper(FE_C_T_B3_2D,FE_C_T_P3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_B3_2D,FE_C_Q_Q3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_P3_2D,FE_C_T_B3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_Q3_2D,FE_C_T_B3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_B3_2D,FE_C_T_B3_2D,C3_2_C3_2);

  RegisterFE2DMapper(FE_C_T_P4_2D,FE_C_T_P4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_P4_2D,FE_C_Q_Q4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_Q4_2D,FE_C_T_P4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_Q4_2D,FE_C_Q_Q4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_P4_2D,FE_C_T_UL4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_UL4_2D,FE_C_T_P4_2D,C4_2_C4_2);
  
  RegisterFE2DMapper(FE_C_T_B4_2D,FE_C_T_P4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_B4_2D,FE_C_Q_Q4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_P4_2D,FE_C_T_B4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_Q4_2D,FE_C_T_B4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_B4_2D,FE_C_T_B4_2D,C4_2_C4_2);

  RegisterFE2DMapper(FE_C_T_P5_2D,FE_C_T_P5_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_T_P5_2D,FE_C_Q_Q5_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_Q_Q5_2D,FE_C_T_P5_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_Q_Q5_2D,FE_C_Q_Q5_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_T_P5_2D,FE_C_T_UL5_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_T_UL5_2D,FE_C_T_P5_2D,C5_2_C5_2);  
  
  RegisterFE2DMapper(FE_C_T_P6_2D,FE_C_T_P6_2D,C6_2_C6_2);
  RegisterFE2DMapper(FE_C_T_P6_2D,FE_C_Q_Q6_2D,C6_2_C6_2);
  RegisterFE2DMapper(FE_C_Q_Q6_2D,FE_C_T_P6_2D,C6_2_C6_2);
  RegisterFE2DMapper(FE_C_Q_Q6_2D,FE_C_Q_Q6_2D,C6_2_C6_2);

  RegisterFE2DMapper(FE_C_T_P7_2D,FE_C_T_P7_2D,C7_2_C7_2);
  RegisterFE2DMapper(FE_C_T_P7_2D,FE_C_Q_Q7_2D,C7_2_C7_2);
  RegisterFE2DMapper(FE_C_Q_Q7_2D,FE_C_T_P7_2D,C7_2_C7_2);
  RegisterFE2DMapper(FE_C_Q_Q7_2D,FE_C_Q_Q7_2D,C7_2_C7_2);

  RegisterFE2DMapper(FE_C_T_P8_2D,FE_C_T_P8_2D,C8_2_C8_2);
  RegisterFE2DMapper(FE_C_T_P8_2D,FE_C_Q_Q8_2D,C8_2_C8_2);
  RegisterFE2DMapper(FE_C_Q_Q8_2D,FE_C_T_P8_2D,C8_2_C8_2);
  RegisterFE2DMapper(FE_C_Q_Q8_2D,FE_C_Q_Q8_2D,C8_2_C8_2);

  RegisterFE2DMapper(FE_C_T_P9_2D,FE_C_T_P9_2D,C9_2_C9_2);
  RegisterFE2DMapper(FE_C_T_P9_2D,FE_C_Q_Q9_2D,C9_2_C9_2);
  RegisterFE2DMapper(FE_C_Q_Q9_2D,FE_C_T_P9_2D,C9_2_C9_2);
  RegisterFE2DMapper(FE_C_Q_Q9_2D,FE_C_Q_Q9_2D,C9_2_C9_2);

  RegisterFE2DMapper(FE_N_T_RT0_2D,FE_N_T_RT0_2D,N1_2_N1_2);
  RegisterFE2DMapper(FE_N_T_RT0_2D,FE_N_Q_RT0_2D,N1_2_N1_2);
  RegisterFE2DMapper(FE_N_T_RT1_2D,FE_N_T_RT1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_N_T_RT2_2D,FE_N_T_RT2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_N_T_RT3_2D,FE_N_T_RT3_2D,C3_2_C3_2);
  
  RegisterFE2DMapper(FE_N_T_BDM1_2D,FE_N_T_BDM1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_N_T_BDM2_2D,FE_N_T_BDM2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_N_T_BDM3_2D,FE_N_T_BDM3_2D,C3_2_C3_2);

  RegisterFE2DMapper(FE_N_Q_RT0_2D,FE_N_Q_RT0_2D,N1_2_N1_2);
  RegisterFE2DMapper(FE_N_Q_Q1_2D,FE_N_T_RT0_2D,N1_2_N1_2);
  RegisterFE2DMapper(FE_N_Q_RT1_2D,FE_N_Q_RT1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_N_Q_RT2_2D,FE_N_Q_RT2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_N_Q_RT3_2D,FE_N_Q_RT3_2D,C3_2_C3_2);
  
  RegisterFE2DMapper(FE_N_Q_BDM1_2D,FE_N_Q_BDM1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_N_Q_BDM2_2D,FE_N_Q_BDM2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_N_Q_BDM3_2D,FE_N_Q_BDM3_2D,C3_2_C3_2);
  
  RegisterFE2DMapper(FE_N_T_P1_2D,FE_N_T_P1_2D,N1_2_N1_2);
  RegisterFE2DMapper(FE_N_T_P1_2D,FE_N_Q_Q1_2D,N1_2_N1_2);
  RegisterFE2DMapper(FE_N_Q_Q1_2D,FE_N_T_P1_2D,N1_2_N1_2);
  RegisterFE2DMapper(FE_N_Q_Q1_2D,FE_N_Q_Q1_2D,N1_2_N1_2);
  
  RegisterFE2DMapper(FE_D_Q_P1_2D,FE_D_Q_P1_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P1_2D,FE_D_T_P1_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P1_2D,FE_D_Q_P1_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_P1_2D,FE_D_T_P1_2D,C0_2_C0_2); // ==== //

  RegisterFE2DMapper(FE_D_Q_P2_2D,FE_D_Q_P2_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P2_2D,FE_D_T_P2_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P2_2D,FE_D_Q_P2_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_P2_2D,FE_D_T_P2_2D,C0_2_C0_2); // ==== //

  RegisterFE2DMapper(FE_D_Q_P3_2D,FE_D_Q_P3_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P3_2D,FE_D_T_P3_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P3_2D,FE_D_Q_P3_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_P3_2D,FE_D_T_P3_2D,C0_2_C0_2); // ==== //

  RegisterFE2DMapper(FE_D_Q_P4_2D,FE_D_Q_P4_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P4_2D,FE_D_T_P4_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_T_P4_2D,FE_D_Q_P4_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_P4_2D,FE_D_T_P4_2D,C0_2_C0_2); // ==== //

  RegisterFE2DMapper(FE_D_Q_P5_2D,FE_D_Q_P5_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_P6_2D,FE_D_Q_P6_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_P7_2D,FE_D_Q_P7_2D,C0_2_C0_2); // ==== //

  RegisterFE2DMapper(FE_D_Q_Q1_2D,FE_D_Q_Q1_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_Q2_2D,FE_D_Q_Q2_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_Q3_2D,FE_D_Q_Q3_2D,C0_2_C0_2); // ==== //
  RegisterFE2DMapper(FE_D_Q_Q4_2D,FE_D_Q_Q4_2D,C0_2_C0_2); // ==== //

  RegisterFE2DMapper(FE_N_T_P1MOD_2D,FE_N_T_P1MOD_2D,N2_2_N2_2);

  RegisterFE2DMapper(FE_C_T_P1MINI_2D,FE_C_T_P1MINI_2D,C1_2_C1_2);

  RegisterFE2DMapper(FE_N_Q_Q2_2D,FE_N_Q_Q2_2D,N2_2_N2_2);
  RegisterFE2DMapper(FE_N_Q_Q2_2D,FE_N_T_P2_2D,N2_2_N2_2);
  RegisterFE2DMapper(FE_N_T_P2_2D,FE_N_Q_Q2_2D,N2_2_N2_2);
  RegisterFE2DMapper(FE_N_T_P2_2D,FE_N_T_P2_2D,N2_2_N2_2);

  RegisterFE2DMapper(FE_N_Q_Q3_2D,FE_N_Q_Q3_2D,N3_2_N3_2);
  RegisterFE2DMapper(FE_N_Q_Q3_2D,FE_N_T_P3_2D,N3_2_N3_2);
  RegisterFE2DMapper(FE_N_T_P3_2D,FE_N_Q_Q3_2D,N3_2_N3_2);
  RegisterFE2DMapper(FE_N_T_P3_2D,FE_N_T_P3_2D,N3_2_N3_2);

  RegisterFE2DMapper(FE_N_Q_Q4_2D,FE_N_Q_Q4_2D,N4_2_N4_2);
  RegisterFE2DMapper(FE_N_Q_Q4_2D,FE_N_T_P4_2D,N4_2_N4_2);
  RegisterFE2DMapper(FE_N_T_P4_2D,FE_N_Q_Q4_2D,N4_2_N4_2);
  RegisterFE2DMapper(FE_N_T_P4_2D,FE_N_T_P4_2D,N4_2_N4_2);

  RegisterFE2DMapper(FE_N_Q_Q5_2D,FE_N_Q_Q5_2D,N5_2_N5_2);
  RegisterFE2DMapper(FE_N_Q_Q5_2D,FE_N_T_P5_2D,N5_2_N5_2);
  RegisterFE2DMapper(FE_N_T_P5_2D,FE_N_Q_Q5_2D,N5_2_N5_2);
  RegisterFE2DMapper(FE_N_T_P5_2D,FE_N_T_P5_2D,N5_2_N5_2);


  RegisterFE2DMapper(FE_B_Q_IB2_2D,FE_B_Q_IB2_2D,C0_2_C0_2);

  RegisterFE2DMapper(FE_D_Q_D2_2D,FE_D_Q_D2_2D,C0_2_C0_2);
  
  //========LOCALPROJECTION=============
  RegisterFE2DMapper(FE_C_Q_UL1_2D,FE_C_Q_Q1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_Q_UL2_2D,FE_C_Q_Q2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_UL3_2D,FE_C_Q_Q3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_UL4_2D,FE_C_Q_Q4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_UL5_2D,FE_C_Q_Q5_2D,C5_2_C5_2);

  RegisterFE2DMapper(FE_C_Q_Q1_2D,FE_C_Q_UL1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_Q_Q2_2D,FE_C_Q_UL2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_Q3_2D,FE_C_Q_UL3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_Q4_2D,FE_C_Q_UL4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_Q5_2D,FE_C_Q_UL5_2D,C5_2_C5_2);

  RegisterFE2DMapper(FE_C_Q_UL1_2D,FE_C_Q_UL1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_Q_UL2_2D,FE_C_Q_UL2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_UL3_2D,FE_C_Q_UL3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_UL4_2D,FE_C_Q_UL4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_UL5_2D,FE_C_Q_UL5_2D,C5_2_C5_2);

  RegisterFE2DMapper(FE_C_T_UL1_2D,FE_C_Q_UL1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_T_UL2_2D,FE_C_Q_UL2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_UL3_2D,FE_C_Q_UL3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_UL4_2D,FE_C_Q_UL4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_UL5_2D,FE_C_Q_UL5_2D,C5_2_C5_2);

  RegisterFE2DMapper(FE_C_Q_UL1_2D,FE_C_T_UL1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_Q_UL2_2D,FE_C_T_UL2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_UL3_2D,FE_C_T_UL3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_UL4_2D,FE_C_T_UL4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_UL5_2D,FE_C_T_UL5_2D,C5_2_C5_2);

  RegisterFE2DMapper(FE_C_T_UL1_2D,FE_C_T_UL1_2D,C1_2_C1_2);
  RegisterFE2DMapper(FE_C_T_UL2_2D,FE_C_T_UL2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_UL3_2D,FE_C_T_UL3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_UL4_2D,FE_C_T_UL4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_UL5_2D,FE_C_T_UL5_2D,C5_2_C5_2);

  RegisterFE2DMapper(FE_C_Q_UL2S_2D,FE_C_Q_UL2S_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_UL3S_2D,FE_C_Q_UL3S_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_UL4S_2D,FE_C_Q_UL4S_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_UL5S_2D,FE_C_Q_UL5S_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_Q_UL6S_2D,FE_C_Q_UL6S_2D,C6_2_C6_2);
  RegisterFE2DMapper(FE_C_Q_UL7S_2D,FE_C_Q_UL7S_2D,C7_2_C7_2);
  RegisterFE2DMapper(FE_C_Q_UL8S_2D,FE_C_Q_UL8S_2D,C8_2_C8_2);
  RegisterFE2DMapper(FE_C_Q_UL9S_2D,FE_C_Q_UL9S_2D,C9_2_C9_2);

  RegisterFE2DMapper(FE_C_Q_UL2SE_2D,FE_C_Q_UL2SE_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_UL3SE_2D,FE_C_Q_UL3SE_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_UL4SE_2D,FE_C_Q_UL4SE_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_UL5SE_2D,FE_C_Q_UL5SE_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_Q_UL6SE_2D,FE_C_Q_UL6SE_2D,C6_2_C6_2);
  RegisterFE2DMapper(FE_C_Q_UL7SE_2D,FE_C_Q_UL7SE_2D,C7_2_C7_2);
  RegisterFE2DMapper(FE_C_Q_UL8SE_2D,FE_C_Q_UL8SE_2D,C8_2_C8_2);
  RegisterFE2DMapper(FE_C_Q_UL9SE_2D,FE_C_Q_UL9SE_2D,C9_2_C9_2);

  RegisterFE2DMapper(FE_C_Q_M2_2D,FE_C_Q_M2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_M3_2D,FE_C_Q_M3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_M4_2D,FE_C_Q_M4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_M5_2D,FE_C_Q_M5_2D,C5_2_C5_2);
  RegisterFE2DMapper(FE_C_Q_M6_2D,FE_C_Q_M6_2D,C6_2_C6_2);
  RegisterFE2DMapper(FE_C_Q_M7_2D,FE_C_Q_M7_2D,C7_2_C7_2);
  RegisterFE2DMapper(FE_C_Q_M8_2D,FE_C_Q_M8_2D,C8_2_C8_2);
  RegisterFE2DMapper(FE_C_Q_M9_2D,FE_C_Q_M9_2D,C9_2_C9_2);

  RegisterFE2DMapper(FE_C_T_UL2_2D,FE_C_Q_UL2S_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_T_UL3_2D,FE_C_Q_UL3S_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_T_UL4_2D,FE_C_Q_UL4S_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_T_UL5_2D,FE_C_Q_UL5S_2D,C5_2_C5_2);

  RegisterFE2DMapper(FE_C_Q_UL2S_2D,FE_C_T_UL2_2D,C2_2_C2_2);
  RegisterFE2DMapper(FE_C_Q_UL3S_2D,FE_C_T_UL3_2D,C3_2_C3_2);
  RegisterFE2DMapper(FE_C_Q_UL4S_2D,FE_C_T_UL4_2D,C4_2_C4_2);
  RegisterFE2DMapper(FE_C_Q_UL5S_2D,FE_C_T_UL5_2D,C5_2_C5_2);

  //====================================

// ========================================================================
// ONE regular grid, same pattern on both sides
// ========================================================================
  RegisterFE2DMapper1Reg(FE_C_T_P0_2D,FE_C_T_P0_2D,C0_2_C0_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P0_2D,FE_C_Q_Q0_2D,C0_2_C0_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q0_2D,FE_C_T_P0_2D,C0_2_C0_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q0_2D,FE_C_Q_Q0_2D,C0_2_C0_2_1Reg);

  RegisterFE2DMapper1Reg(FE_C_T_P1_2D,FE_C_T_P1_2D,C1_2_C1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P1_2D,FE_C_Q_Q1_2D,C1_2_C1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q1_2D,FE_C_T_P1_2D,C1_2_C1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q1_2D,FE_C_Q_Q1_2D,C1_2_C1_2_1Reg);

  RegisterFE2DMapper1Reg(FE_C_T_P2_2D,FE_C_T_P2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P2_2D,FE_C_Q_Q2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q2_2D,FE_C_T_P2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q2_2D,FE_C_Q_Q2_2D,C2_2_C2_2_1Reg);

  RegisterFE2DMapper1Reg(FE_C_T_B2_2D,FE_C_T_P2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_B2_2D,FE_C_Q_Q2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P2_2D,FE_C_T_B2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q2_2D,FE_C_T_B2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_B2_2D,FE_C_T_B2_2D,C2_2_C2_2_1Reg);
  
  RegisterFE2DMapper1Reg(FE_C_T_SV2_2D,FE_C_T_P2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_SV2_2D,FE_C_Q_Q2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P2_2D,FE_C_T_SV2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q2_2D,FE_C_T_SV2_2D,C2_2_C2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_SV2_2D,FE_C_T_SV2_2D,C2_2_C2_2_1Reg);

  RegisterFE2DMapper1Reg(FE_C_T_P3_2D,FE_C_T_P3_2D,C3_2_C3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P3_2D,FE_C_Q_Q3_2D,C3_2_C3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q3_2D,FE_C_T_P3_2D,C3_2_C3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q3_2D,FE_C_Q_Q3_2D,C3_2_C3_2_1Reg);

  RegisterFE2DMapper1Reg(FE_C_T_B3_2D,FE_C_T_P3_2D,C3_2_C3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_B3_2D,FE_C_Q_Q3_2D,C3_2_C3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P3_2D,FE_C_T_B3_2D,C3_2_C3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q3_2D,FE_C_T_B3_2D,C3_2_C3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_B3_2D,FE_C_T_B3_2D,C3_2_C3_2_1Reg);

//fourth order
  RegisterFE2DMapper1Reg(FE_C_T_P4_2D,FE_C_T_P4_2D,C4_2_C4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P4_2D,FE_C_Q_Q4_2D,C4_2_C4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q4_2D,FE_C_T_P4_2D,C4_2_C4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q4_2D,FE_C_Q_Q4_2D,C4_2_C4_2_1Reg);

  RegisterFE2DMapper1Reg(FE_C_T_B4_2D,FE_C_T_P4_2D,C4_2_C4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_B4_2D,FE_C_Q_Q4_2D,C4_2_C4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P4_2D,FE_C_T_B4_2D,C4_2_C4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q4_2D,FE_C_T_B4_2D,C4_2_C4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_B4_2D,FE_C_T_B4_2D,C4_2_C4_2_1Reg);

//fifth order
  RegisterFE2DMapper1Reg(FE_C_T_P5_2D,FE_C_T_P5_2D,C5_2_C5_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_T_P5_2D,FE_C_Q_Q5_2D,C5_2_C5_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q5_2D,FE_C_T_P5_2D,C5_2_C5_2_1Reg);
  RegisterFE2DMapper1Reg(FE_C_Q_Q5_2D,FE_C_Q_Q5_2D,C5_2_C5_2_1Reg);

  RegisterFE2DMapper1Reg(FE_N_T_P1_2D,FE_N_T_P1_2D,N1_2_N1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P1_2D,FE_N_Q_Q1_2D,N1_2_N1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q1_2D,FE_N_T_P1_2D,N1_2_N1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q1_2D,FE_N_Q_Q1_2D,N1_2_N1_2_1Reg);

// second order nonconforming
  RegisterFE2DMapper1Reg(FE_N_Q_Q2_2D,FE_N_Q_Q2_2D,N2_2_N2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P2_2D,FE_N_T_P2_2D,N2_2_N2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q2_2D,FE_N_T_P2_2D,N2_2_N2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P2_2D,FE_N_Q_Q2_2D,N2_2_N2_2_1Reg);

// third order nonconforming
  RegisterFE2DMapper1Reg(FE_N_Q_Q3_2D,FE_N_Q_Q3_2D,N3_2_N3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P3_2D,FE_N_T_P3_2D,N3_2_N3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q3_2D,FE_N_T_P3_2D,N3_2_N3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P3_2D,FE_N_Q_Q3_2D,N3_2_N3_2_1Reg);

// fourth order nonconforming
  RegisterFE2DMapper1Reg(FE_N_Q_Q4_2D,FE_N_Q_Q4_2D,N4_2_N4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P4_2D,FE_N_T_P4_2D,N4_2_N4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q4_2D,FE_N_T_P4_2D,N4_2_N4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P4_2D,FE_N_Q_Q4_2D,N4_2_N4_2_1Reg);

// fifth order nonconforming
  RegisterFE2DMapper1Reg(FE_N_Q_Q5_2D,FE_N_Q_Q5_2D,N5_2_N5_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P5_2D,FE_N_T_P5_2D,N5_2_N5_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q5_2D,FE_N_T_P5_2D,N5_2_N5_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P5_2D,FE_N_Q_Q5_2D,N5_2_N5_2_1Reg);

  RegisterFE2DMapper1Reg(FE_D_Q_P1_2D,FE_D_Q_P1_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_T_P1_2D,FE_D_T_P1_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_T_P1_2D,FE_D_Q_P1_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_Q_P1_2D,FE_D_T_P1_2D,C0_2_C0_2_1Reg); // ==== //
  
  RegisterFE2DMapper1Reg(FE_D_T_SV1_2D,FE_D_T_SV1_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_T_SV1_2D,FE_D_Q_P1_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_Q_P1_2D,FE_D_T_SV1_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_Q_P2_2D,FE_D_Q_P2_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_T_P2_2D,FE_D_T_P2_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_T_P2_2D,FE_D_Q_P2_2D,C0_2_C0_2_1Reg); // ==== //
  RegisterFE2DMapper1Reg(FE_D_Q_P2_2D,FE_D_T_P2_2D,C0_2_C0_2_1Reg); // ==== //

// ========================================================================
// ONE regular grid, different pattern
// ========================================================================
  // coarse: second order, fine: first order
  RegisterFE2DMapper1Reg(FE_N_Q_Q2_2D,FE_N_Q_Q1_2D,N2_2_N1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P2_2D,FE_N_T_P1_2D,N2_2_N1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q2_2D,FE_N_T_P1_2D,N2_2_N1_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P2_2D,FE_N_Q_Q1_2D,N2_2_N1_2_1Reg);

  // coarse: third order, fine: second order
  RegisterFE2DMapper1Reg(FE_N_Q_Q3_2D,FE_N_Q_Q2_2D,N3_2_N2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P3_2D,FE_N_T_P2_2D,N3_2_N2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q3_2D,FE_N_T_P2_2D,N3_2_N2_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P3_2D,FE_N_Q_Q2_2D,N3_2_N2_2_1Reg);

  // coarse: fourth order, fine: third order
  RegisterFE2DMapper1Reg(FE_N_Q_Q4_2D,FE_N_Q_Q3_2D,N4_2_N3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P4_2D,FE_N_T_P3_2D,N4_2_N3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q4_2D,FE_N_T_P3_2D,N4_2_N3_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P4_2D,FE_N_Q_Q3_2D,N4_2_N3_2_1Reg);

  // coarse: fifth order, fine: fourth order
  RegisterFE2DMapper1Reg(FE_N_Q_Q5_2D,FE_N_Q_Q4_2D,N5_2_N4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P5_2D,FE_N_T_P4_2D,N5_2_N4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_Q_Q5_2D,FE_N_T_P4_2D,N5_2_N4_2_1Reg);
  RegisterFE2DMapper1Reg(FE_N_T_P5_2D,FE_N_Q_Q4_2D,N5_2_N4_2_1Reg);

  RegisterFE2DMapper(FE_C_Q_EL1_2D,FE_C_Q_EL1_2D,C1_2_C1_2);   
  RegisterFE2DMapper(FE_C_Q_UL1_2D,FE_C_Q_EL1_2D,C1_2_C1_2); 
  RegisterFE2DMapper(FE_C_Q_EL1_2D,FE_C_Q_UL1_2D,C1_2_C1_2);   
  RegisterFE2DMapper(FE_C_Q_EL1_2D,FE_C_Q_Q1_2D,C1_2_C1_2); 
  RegisterFE2DMapper(FE_C_Q_Q1_2D,FE_C_Q_EL1_2D,C1_2_C1_2);   
  
#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "FE mapper registered" << endl;
}

void TFEDatabase2D::RegisterAllHangingNodes()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  RegisterHNDesc2D(HN_C_P1_2D_0, HN_C_P1_2D_0_Obj);

  RegisterHNDesc2D(HN_C_P2_2D_0, HN_C_P2_2D_0_Obj);
  RegisterHNDesc2D(HN_C_P2_2D_1, HN_C_P2_2D_1_Obj);

  RegisterHNDesc2D(HN_C_P3_2D_0, HN_C_P3_2D_0_Obj);
  RegisterHNDesc2D(HN_C_P3_2D_1, HN_C_P3_2D_1_Obj);
  RegisterHNDesc2D(HN_C_P3_2D_2, HN_C_P3_2D_2_Obj);

//fourth order
  RegisterHNDesc2D(HN_C_P4_2D_0, HN_C_P4_2D_0_Obj);
  RegisterHNDesc2D(HN_C_P4_2D_1, HN_C_P4_2D_1_Obj);
  RegisterHNDesc2D(HN_C_P4_2D_2, HN_C_P4_2D_2_Obj);
  RegisterHNDesc2D(HN_C_P4_2D_3, HN_C_P4_2D_3_Obj);

//fifth order
  RegisterHNDesc2D(HN_C_P5_2D_0, HN_C_P5_2D_0_Obj);
  RegisterHNDesc2D(HN_C_P5_2D_1, HN_C_P5_2D_1_Obj);
  RegisterHNDesc2D(HN_C_P5_2D_2, HN_C_P5_2D_2_Obj);
  RegisterHNDesc2D(HN_C_P5_2D_3, HN_C_P5_2D_3_Obj);
  RegisterHNDesc2D(HN_C_P5_2D_4, HN_C_P5_2D_4_Obj);

  RegisterHNDesc2D(HN_N_P1_2D_0, HN_N_P1_2D_0_Obj);
  RegisterHNDesc2D(HN_N_P2_2D_0, HN_N_P2_2D_0_Obj);
  RegisterHNDesc2D(HN_N_P3_2D_0, HN_N_P3_2D_0_Obj);
  RegisterHNDesc2D(HN_N_P4_2D_0, HN_N_P4_2D_0_Obj);
  RegisterHNDesc2D(HN_N_P5_2D_0, HN_N_P5_2D_0_Obj);


#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "Hanging nodes registered" << endl;
}

void TFEDatabase2D::RegisterAllRefTrans()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  ReferenceTrans1D[LineAffin] = new TLineAffin();
  ReferenceTrans2D[TriaAffin] = new TTriaAffin();
  ReferenceTrans2D[QuadAffin] = new TQuadAffin();
  ReferenceTrans2D[QuadBilinear] = new TQuadBilinear();
  ReferenceTrans2D[TriaIsoparametric] = new TTriaIsoparametric();
  ReferenceTrans2D[QuadIsoparametric] = new TQuadIsoparametric();


#ifdef _MPI
  if(rank==out_rank)
#endif
  cout << "Reference Transformations registered" << endl;
}

void TFEDatabase2D::GenerateArrays()
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  int i;
  TFE2D *ele;
  TBaseFunct2D *bf;

  for(i=0;i<N_FEs2D;i++)
  {
    ele = FEs2D[i];
    if(ele)
    {
      ele->GetFEDesc2D(FEDesc2D_IDFromFE2D[i],
                       FEDesc2DFromFE2D[i]);

      ele->GetBaseFunct2D(BaseFunct2D_IDFromFE2D[i],
                          bf);
      BaseFunct2DFromFE2D[i] = bf; 
      N_BaseFunctFromFE2D[i] = bf->GetDimension();
      PolynomialDegreeFromFE2D[i] = bf->GetPolynomialDegree();
      AccuracyFromFE2D[i] = bf->GetAccuracy();
      RefElementFromFE2D[i] = bf->GetRefElement();

      ele->GetNodalFunctional2D(NodalFunctional2D_IDFromFE2D[i],
                                NodalFunctional2DFromFE2D[i]);

      RefTrans2D_IDFromFE2D[i] = ele->GetRefTransID(); 
    } // endif
  } // endfor i

  QFQuadFromDegree[0] = Gauss2Quad;
  QFQuadFromDegree[1] = Gauss2Quad;
  QFQuadFromDegree[2] = Gauss2Quad;
  QFQuadFromDegree[3] = Gauss3Quad;
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
  QFQuadFromDegree[18] = Gauss9Quad;

  HighestAccuracyQuad = 18;

//   QFTriaFromDegree[0] = MidPointTria; // points on boundary problem in Axial3D assembling
//   QFTriaFromDegree[1] = MidPointTria;
//   QFTriaFromDegree[2] = MidPointTria;  
//   QFTriaFromDegree[3] = SevenPointTria;  
  QFTriaFromDegree[0] = Gauss3Tria; 
  QFTriaFromDegree[1] = Gauss3Tria;
  QFTriaFromDegree[2] = Gauss3Tria;  
  QFTriaFromDegree[3] = Gauss3Tria;  
  QFTriaFromDegree[4] = Gauss3Tria;
  QFTriaFromDegree[5] = Gauss3Tria;
  QFTriaFromDegree[6] = Degree9Tria;
  QFTriaFromDegree[7] = Degree9Tria;
  QFTriaFromDegree[8] = Degree9Tria;
  QFTriaFromDegree[9] = Degree9Tria;
  QFTriaFromDegree[10] = Degree19Tria;
  QFTriaFromDegree[11] = Degree11Tria;
  QFTriaFromDegree[12] = Degree19Tria;
  QFTriaFromDegree[13] = Degree19Tria;
  QFTriaFromDegree[14] = Degree19Tria;
  QFTriaFromDegree[15] = Degree19Tria;
  QFTriaFromDegree[16] = Degree19Tria;
  QFTriaFromDegree[17] = Degree19Tria;
  QFTriaFromDegree[18] = Degree19Tria;
  QFTriaFromDegree[19] = Degree19Tria;

  QFTriaFromDegree[20] = Gauss_Degree8Tria;
  // Scott-Vogelius
  QFTriaFromDegree[21] = CompGauss3Tria;
  HighestAccuracyTria = 21;

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

#ifdef __SCOTT_VOGELIUS__
  for(i=0;i<HighestAccuracyTria;i++)
    QFTriaFromDegree[i] = CompGauss3Tria;
#endif

} // end TFEDatabase2D::GenerateArrays

/** calculate base functions with derivatives and coordinates
    from reference to original element */
RefTrans2D TFEDatabase2D::GetOrig(int N_LocalUsedElements, 
                          FE2D *LocalUsedElements,
                          TCollection *Coll,
                          TBaseCell *cell, bool *Needs2ndDer,
                          int &N_Points, double* &xi, double* &eta, 
                          double* &weights, double* X, double* Y,
                          double* absdetjk)
{
  int i,j, MaxPolynomialDegree, PolynomialDegree, N_Edges, N_terms;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  RefTrans2D RefTrans, *RefTransArray, CurrentRefTrans;
  TRefTrans2D *rt;
  BaseFunct2D BaseFuncts[N_FEs2D];
  JointType jointtype;
  TJoint *joint;
  bool IsIsoparametric;
  BoundTypes bdtype;
  int MaxApproxOrder, ApproxOrder;

  BaseFunct2D BaseFunct;
  double **origvaluesD00, **origvaluesD10, **origvaluesD01;
  double **origvaluesD20, **origvaluesD11, **origvaluesD02;
  int N_Functs;
  int BaseVectDim = 1;

#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  // find adequate quadrature formula for all elements
  // and find needed reference transformation
  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();
  RefTrans = RefTransArray[LocalUsedElements[0]];
  MaxPolynomialDegree = 0;
  MaxApproxOrder = 0;
  RefElement = TFEDatabase2D::GetRefElementFromFE2D(LocalUsedElements[0]);
  
  for(i=0;i<N_LocalUsedElements;i++)
  {
    BaseFuncts[i] = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(LocalUsedElements[i]);
    PolynomialDegree = PolynomialDegreeFromFE2D[LocalUsedElements[i]];
    if(PolynomialDegree > MaxPolynomialDegree) 
      MaxPolynomialDegree = PolynomialDegree;
    ApproxOrder = AccuracyFromFE2D[LocalUsedElements[i]];
    if(ApproxOrder > MaxApproxOrder)
      MaxApproxOrder = ApproxOrder;

    CurrentRefTrans = RefTransArray[LocalUsedElements[i]];
    if(CurrentRefTrans > RefTrans)
      RefTrans = CurrentRefTrans;
  }
//   cout << "MaxPolynomialDegree: " << MaxPolynomialDegree << endl;
  // cout << "MaxApproxOrder: " << MaxApproxOrder << endl;
  // cout << "RefElement: " << RefElement << endl;
  // cout << "RefTrans: " << RefTrans << endl;

  // for isoparametric pressure error
  if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
  {
    if(MaxApproxOrder < (ABS(TDatabase::ParamDB->VELOCITY_SPACE))%10)
      MaxApproxOrder = (ABS(TDatabase::ParamDB->VELOCITY_SPACE))%10;
    // cout << "MaxApproxOrder: " << MaxApproxOrder << endl;
  }
  
  // number of terms in products that need to be assembled
  if (TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR)
    N_terms = 2;
  else
    N_terms = 3;

  switch(RefElement)
  {
    case BFUnitSquare:
      QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(N_terms*MaxPolynomialDegree);
      
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 95)
      {
        QuadFormula = GetQFQuadFromDegree(3);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 97)
      {
	  QuadFormula = GetQFQuadFromDegree(18);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 96)
      {
	  QuadFormula = GetQFQuadFromDegree(9);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 99)
      {
	  QuadFormula = GetQFQuadFromDegree(18);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_QUAD<N_terms*MaxPolynomialDegree)
      {
       #ifdef _MPI
       if(rank==out_rank)
       #endif
       {
        switch(3*MaxPolynomialDegree)
        {
          case 0:
            OutPut("Quadrature formula for quadrilateral is Gauss2" << endl); 
            break;
          case 3:
            OutPut("Quadrature formula for quadrilateral is Gauss3" << endl); 
            break;
          case 6:
            OutPut("Quadrature formula for quadrilateral is Gauss4" << endl); 
            break;
          case 9:
            OutPut("Quadrature formula for quadrilateral is Gauss5" << endl); 
            break;
          case 12:
            OutPut("Quadrature formula for quadrilateral is Gauss7" << endl); 
            break;
          case 15:
            OutPut("Quadrature formula for quadrilateral is Gauss8" << endl); 
            break;
          case 18:
            OutPut("Quadrature formula for quadrilateral is Gauss9" << endl); 
            break;
        }
        if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 95)
      {
        OutPut("Quadrature formula for quadrilateral is Gauss3" << endl); 
      }
	if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 97)
	{
	    OutPut("Quadrature formula for quadrilateral is set to Gauss9" << endl); 
	}
	if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 96)
	{
	    OutPut("Quadrature formula for quadrilateral is set to Gauss5" << endl); 
	}
	if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 99)
	{
	    OutPut("Quadrature formula for quadrilateral is set to Gauss9" << endl); 
	}
        }


        TDatabase::ParamDB->INTERNAL_QUAD_QUAD = N_terms*MaxPolynomialDegree;
	TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE = 4;
      }
      N_Edges = 4;
    break;

    case BFUnitTriangle:
      if (MaxPolynomialDegree>0)
	  QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(N_terms*MaxPolynomialDegree-1);
      else
          QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(0);
      // more accurate quad rule for JohnMaubachTobiska1997.h
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 97)
      {
         QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(11);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 98)
      {
         QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(5);
      }      
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 99)
      {
	  QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(20);
      }
      
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 96)
      {
	  QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(20);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == -22)
      {
    QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(21);
      }
            if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 90)
      {
    QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
      }
      
      
      //   QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(5);
      if (TDatabase::ParamDB->INTERNAL_QUAD_TRIA<N_terms*MaxPolynomialDegree)
      {
       #ifdef _MPI
       if(rank==0)
       #endif
       {
        switch(N_terms*MaxPolynomialDegree)
        {
          case 3:
            OutPut("Quadrature formula for triangles is MidPointTria" << endl); 
            break;
          case 6:
            OutPut("Quadrature formula for triangles is Gauss3Tria" << endl); 
            break;
          case 9:
            OutPut("Quadrature formula for triangles is  Degree9Tria"<< endl); 
            break;
          case 12:
            OutPut("Quadrature formula for triangles is Degree19Tria" << endl); 
            break;
          case 15:
            OutPut("Quadrature formula for triangles is Degree19Tria" << endl); 
            break;
          case 18:
            OutPut("Quadrature formula for triangles is Degree19Tria" << endl); 
            break;
         }
        }
       
        
        TDatabase::ParamDB->INTERNAL_QUAD_TRIA = N_terms*MaxPolynomialDegree;
	TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE = 3;

        #ifdef _MPI
        if(rank==out_rank)
        #endif
        {
	if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 97)
	{
	    OutPut("Quadrature formula for triangles is set to Degree19Tria" << endl); 
	}
	if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 98)
	{
	    OutPut("Quadrature formula for triangles is set to Gauss3Tria" << endl); 
 
	}
	if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 99)
	{
	    OutPut("Quadrature formula for triangles is set to Gauss_Degree8Tria" << endl); 
 
	}
	if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 96)
	{
	    OutPut("Quadrature formula for triangles is set to Gauss_Degree8Tria" << endl); 
 
	}
		if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 90)
	{
	    OutPut("Quadrature formula for triangles is set to Degree9Tria" << endl); 
 
	}
	
       }
      }
      N_Edges = 3;
    break;

    default:
      N_Edges = 0;
  } // endswitch

  IsIsoparametric = FALSE;
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
  {
    for(i=0;i<N_Edges;i++)
    {
      joint = cell->GetJoint(i);
      jointtype = joint->GetType();
      if(jointtype == BoundaryEdge)
      {
        bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Line)
          IsIsoparametric = TRUE;
      }
      if(jointtype == InterfaceJoint)
      {
        bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Line)
          IsIsoparametric = TRUE;
      }
      if(jointtype == IsoInterfaceJoint ||
         jointtype == IsoBoundEdge)
        IsIsoparametric = TRUE;
    }
  }// endif 

  if(IsIsoparametric)
  {
    switch(RefElement)
    {
      case BFUnitSquare:
        RefTrans = QuadIsoparametric;
      break;

      case BFUnitTriangle:
        RefTrans = TriaIsoparametric;
      break;
    }
  } // endif IsIsoparametric

  //cout << "QuadFormula: " << QuadFormula << endl;
  qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
  qf2->GetFormulaData(N_Points, weights, xi, eta);

  // calculate the values of base functions and their derivatives
  // on the original element
  rt = TFEDatabase2D::ReferenceTrans2D[RefTrans];
  switch(RefTrans)
  {
    case TriaAffin:
	//cout << "TriaAffin" << endl;
      ((TTriaAffin *)rt)->SetCell(cell);
      ((TTriaAffin *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta,
                               QuadFormula,
                               Needs2ndDer);

      ((TTriaAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, 
                                         X, Y, absdetjk);
      break;
    case QuadAffin:
	// cout << "QuadAffin" << endl;
      ((TQuadAffin *)rt)->SetCell(cell);
      ((TQuadAffin *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta,
                               QuadFormula,
                               Needs2ndDer);
      ((TQuadAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, 
                                         X, Y, absdetjk);
    break;
    case QuadBilinear:
      // cout << "QuadBilinear" << endl;
      ((TQuadBilinear *)rt)->SetCell(cell);
      ((TQuadBilinear *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta,
                               QuadFormula,
                               Needs2ndDer);
      ((TQuadBilinear *)rt)->GetOrigFromRef(N_Points, xi, eta, 
                                            X, Y, absdetjk);
    break;
    case TriaIsoparametric:
      // cout << "TriaIsoparametric" << endl;
      ((TTriaIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
      ((TTriaIsoparametric *)rt)->SetCell(cell);
      ((TTriaIsoparametric *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta,
                               QuadFormula,
                               Needs2ndDer);
      ((TTriaIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, 
                                         X, Y, absdetjk);
    break;
    case QuadIsoparametric:
      // cout << "QuadIsoparametric" << endl;
      ((TQuadIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
      ((TQuadIsoparametric *)rt)->SetCell(cell);
      ((TQuadIsoparametric *)rt)->GetOrigValues(N_LocalUsedElements, 
                               BaseFuncts,
                               N_Points, xi, eta,
                               QuadFormula,
                               Needs2ndDer);
      ((TQuadIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, 
                                         X, Y, absdetjk);
    break;
  } // endswitch

  for(i=0;i<N_LocalUsedElements;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();
    origvaluesD00=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);

    BaseFuncts2D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                      origvaluesD00);
    BaseFuncts2D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                      origvaluesD10);
    BaseFuncts2D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                      origvaluesD01);
    if(Needs2ndDer[i])
    {
      origvaluesD20=TFEDatabase2D::GetOrigElementValues(BaseFunct, D20);
      origvaluesD11=TFEDatabase2D::GetOrigElementValues(BaseFunct, D11);
      origvaluesD02=TFEDatabase2D::GetOrigElementValues(BaseFunct, D02);
      BaseFuncts2D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD20);
      BaseFuncts2D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD11);
      BaseFuncts2D[BaseFunct]->ChangeBF(Coll, cell, N_Points,
                                        origvaluesD02);
    } // endif Needs2ndDer[i]
  } // endfor i

  return RefTrans;
}


/** calculate the values of base functions or their derivatives
    on the original element */
void TFEDatabase2D::GetOrigValues(RefTrans2D RefTrans,
                double xi, double eta, TBaseFunct2D *bf,
                TCollection *Coll, TGridCell *cell,
                double *uref, double *uxiref, double *uetaref,
                double *uorig, double *uxorig, double *uyorig)
{
  TRefTrans2D *rt;
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();
  rt = ReferenceTrans2D[RefTrans];
  switch(RefTrans)
  {
    case TriaAffin: 
      ((TTriaAffin *)rt)->GetOrigValues(xi, eta, N_BaseFunct, 
                  uref, uxiref, uetaref, uorig, uxorig, uyorig, BaseVectDim);
      break;
    case QuadAffin:
       ((TQuadAffin *)rt)->GetOrigValues(xi, eta, N_BaseFunct, 
                  uref, uxiref, uetaref, uorig, uxorig, uyorig, BaseVectDim);
      break;
    case QuadBilinear:
      ((TQuadBilinear *)rt)->GetOrigValues(xi, eta, N_BaseFunct, 
                  uref, uxiref, uetaref, uorig, uxorig, uyorig, BaseVectDim);
      break;
    case TriaIsoparametric: 
      ((TTriaIsoparametric *)rt)->GetOrigValues(xi, eta, N_BaseFunct, 
                  uref, uxiref, uetaref, uorig, uxorig, uyorig);
      break;
    case QuadIsoparametric:
      ((TQuadIsoparametric *)rt)->GetOrigValues(xi, eta, N_BaseFunct, 
                  uref, uxiref, uetaref, uorig, uxorig, uyorig);
      break;
    default:
      cout << "wrong reference transformation identifier" << endl;
      break;
  } // endswitch

  bf->ChangeBF(Coll, cell, uorig);
  bf->ChangeBF(Coll, cell, uxorig);
  bf->ChangeBF(Coll, cell, uyorig);
}

void TFEDatabase2D::GetOrigValues(RefTrans2D RefTrans,
                                  double zeta, TBaseFunct2D *bf, int edgeNumber,
                                  double *uref, double *uxiref, double *uetaref,
                                  double *uorig, double *uxorig, double *uyorig)
{
  TRefTrans2D *rt;
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();
  rt = ReferenceTrans2D[RefTrans];
  switch(RefTrans)
  {
    case TriaAffin:
      ((TTriaAffin *)rt)->GetOrigValues(edgeNumber, zeta, N_BaseFunct,
          uref, uxiref, uetaref,
          uorig, uxorig, uyorig,
          BaseVectDim);
      break;
    case QuadAffin:
      ((TQuadAffin *)rt)->GetOrigValues(edgeNumber, zeta, N_BaseFunct,
          uref, uxiref, uetaref,
          uorig, uxorig, uyorig,
          BaseVectDim);
      break;
    case QuadBilinear:
      ((TQuadBilinear *)rt)->GetOrigValues(edgeNumber, zeta, N_BaseFunct,
          uref, uxiref, uetaref, uorig, uxorig, uyorig, BaseVectDim);
      break;
    case TriaIsoparametric:
      ((TTriaIsoparametric *)rt)->GetOrigValues(edgeNumber, zeta, N_BaseFunct,
          uref, uxiref, uetaref, uorig, uxorig, uyorig);
      break;
    case QuadIsoparametric:
      ((TQuadIsoparametric *)rt)->GetOrigValues(edgeNumber, zeta, N_BaseFunct,
          uref, uxiref, uetaref, uorig, uxorig, uyorig);
      break;
    default:
      cout << "wrong reference transformation identifier" << endl;
      break;
    }                                               // endswitch
}

void TFEDatabase2D::GetOrigFromRef(RefTrans2D RefTrans, int n_points, 
                        double *xi, double *eta,
                        double *X, double *Y, double *absdetjk)
{
  TRefTrans2D *rt;

  rt = ReferenceTrans2D[RefTrans];
  switch(RefTrans)
  {
    case TriaAffin: 
      ((TTriaAffin *)rt)->GetOrigFromRef(n_points, xi, eta, 
                                         X, Y, absdetjk);
      break;
    case QuadAffin:
      ((TQuadAffin *)rt)->GetOrigFromRef(n_points, xi, eta, 
                                         X, Y, absdetjk);
      break;
    case QuadBilinear:
      ((TQuadBilinear *)rt)->GetOrigFromRef(n_points, xi, eta, 
                                            X, Y, absdetjk);
      break;
    case TriaIsoparametric: 
      ((TTriaIsoparametric *)rt)->GetOrigFromRef(n_points, xi, eta, 
                                         X, Y, absdetjk);
      break;
    case QuadIsoparametric:
      ((TQuadIsoparametric *)rt)->GetOrigFromRef(n_points, xi, eta, 
                                         X, Y, absdetjk);
      break;
    default:
      cout << "wrong reference transformation identifier" << endl;
      break;
  } // endswitch
}

/** calculate points on reference element */
void TFEDatabase2D::GetRefFromOrig(RefTrans2D RefTrans,
                                 double X, double Y,
                                 double &xi, double &eta)
{
  TRefTrans2D *rt;

  rt = ReferenceTrans2D[RefTrans];
  switch(RefTrans)
  {
    case TriaAffin: 
      ((TTriaAffin *)rt)->GetRefFromOrig(X, Y, xi, eta);
      break;
    case QuadAffin:
      ((TQuadAffin *)rt)->GetRefFromOrig(X, Y, xi, eta);
      break;
    case QuadBilinear:
      ((TQuadBilinear *)rt)->GetRefFromOrig(X, Y, xi, eta);
      break;
    case TriaIsoparametric: 
      ((TTriaIsoparametric *)rt)->GetRefFromOrig(X, Y, xi, eta);
      break;
    case QuadIsoparametric:
      ((TQuadIsoparametric *)rt)->GetRefFromOrig(X, Y, xi, eta);
      break;
    default:
      cout << "wrong reference transformation identifier" << endl;
      break;
  } // endswitch
}

/** set cell for reference transformation */
void TFEDatabase2D::SetCellForRefTrans(TBaseCell *cell, 
                                     RefTrans2D reftrans)
{
  TRefTrans2D *rt;

  rt = ReferenceTrans2D[reftrans];

  switch(reftrans)
  {
    case TriaAffin: 
      ((TTriaAffin *)rt)->SetCell(cell);
      break;
    case QuadAffin:
      ((TQuadAffin *)rt)->SetCell(cell);
      break;
    case QuadBilinear:
      ((TQuadBilinear *)rt)->SetCell(cell);
      break;
    case TriaIsoparametric: 
      ((TTriaIsoparametric *)rt)->SetCell(cell);
      break;
    case QuadIsoparametric:
      ((TQuadIsoparametric *)rt)->SetCell(cell);
      break;
    default:
      cout << "wrong reference transformation identifier" << endl;
      break;
  } // endswitch
}

double *TFEDatabase2D::GetProlongationMatrix2D (FE2D parent, 
    Refinements refine, FE2D child, int childnumber)
{ 
  double *ret, *ret2;
  int i,j,k,l;
  int N_Coarse, N_Fine, N_Points, N_Children;
  double *xi, *eta;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double AllPointValues[MaxN_PointsForNodal2D][MaxN_BaseFunctions2D];
  double PointValues[MaxN_PointsForNodal2D];
  TFE2D *CoarseElement, *FineElement;
  TRefDesc *RefDesc;
  TBaseFunct2D *BaseFunctions;
  BaseFunct2D Coarse, Fine;
  TGridCell *RefCell, *cell;
  TNodalFunctional2D *nf;
  BF2DRefElements RefElement;
  RefTrans2D F_K;
  TRefTrans2D *rt;

  CoarseElement = TFEDatabase2D::GetFE2D(parent);
  Coarse = CoarseElement->GetBaseFunct2D_ID();
  FineElement = TFEDatabase2D::GetFE2D(child);
  Fine = FineElement->GetBaseFunct2D_ID();

  ret = ProlongationMatrix2D[Coarse][refine][Fine][childnumber]; 

  if(ret == NULL)
  {
    // cerr << "ret == NULL" << endl;
    // prolongation matrix was not generated yet

    BaseFunctions = CoarseElement->GetBaseFunct2D();
    N_Coarse = BaseFunctions->GetDimension();

    ret2 = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];

    if(refine == NoRef)
    {
      ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];

      RefCell = BaseFunctions->GenerateRefElement();
      FineElement = TFEDatabase2D::GetFE2D(child);
      N_Fine = FineElement->GetBaseFunct2D()->GetDimension();

      nf = FineElement->GetNodalFunctional2D();
      nf->GetPointsForAll(N_Points, xi, eta);

      RefElement = BaseFunctions->GetRefElement();

      switch(RefElement)
      {
        case BFUnitSquare:
          rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
          ((TQuadAffin *)rt)->SetCell(RefCell);
          F_K = QuadAffin;
          break;
        case BFUnitTriangle:
          rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
          ((TTriaAffin *)rt)->SetCell(RefCell);
          F_K = TriaAffin;
          break;
        default:
          F_K = TriaAffin;
          break;
      }
      TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta,
                                  X, Y, AbsDetjk);

      for(k=0;k<N_Points;k++)
        BaseFunctions->GetDerivatives(D00, X[k], Y[k], 
                                      AllPointValues[k]);

      for(k=0;k<N_Coarse;k++)
      {
        for(l=0;l<N_Points;l++)
          PointValues[l] = AllPointValues[l][k];

        nf->GetAllFunctionals(NULL, NULL, PointValues,
                              ret2+k*MaxN_BaseFunctions2D);
      }

      for(k=0;k<MaxN_BaseFunctions2D;k++)
        for(l=0;l<MaxN_BaseFunctions2D;l++)
          ret[k*MaxN_BaseFunctions2D+l] = 
                        ret2[l*MaxN_BaseFunctions2D+k];
  
      RegisterProlongationMatrix2D(Coarse, refine, Fine, 
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
        ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  
        cell = (TGridCell *)RefCell->GetChild(j);
        FineElement = TFEDatabase2D::GetFE2D(child);
        Fine = FineElement->GetBaseFunct2D_ID();
        N_Fine = FineElement->GetBaseFunct2D()->GetDimension();
  
        nf = FineElement->GetNodalFunctional2D();
        nf->GetPointsForAll(N_Points, xi, eta);

        switch(cell->GetType())
        {
          case Quadrangle:
          case Parallelogram:
          case Rectangle:
            rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
            ((TQuadAffin *)rt)->SetCell(cell);
            F_K = QuadAffin;
          break;
          case Triangle:
            rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
            ((TTriaAffin *)rt)->SetCell(cell);
            F_K = TriaAffin;
          break;
         default:
           cout << "Not a 2D cell type" << endl;
         break;          
        }
        TFEDatabase2D::GetOrigFromRef(F_K ,N_Points, xi, eta,
                                      X, Y, AbsDetjk);
  
        for(k=0;k<N_Points;k++)
          BaseFunctions->GetDerivatives(D00, X[k], Y[k], 
                        AllPointValues[k]);
  
        for(k=0;k<N_Coarse;k++)
        {
          for(l=0;l<N_Points;l++)
            PointValues[l] = AllPointValues[l][k];
  
          nf->GetAllFunctionals(NULL, NULL, PointValues,
                                ret2+k*MaxN_BaseFunctions2D);
        }

        for(k=0;k<MaxN_BaseFunctions2D;k++)
          for(l=0;l<MaxN_BaseFunctions2D;l++)
            ret[k*MaxN_BaseFunctions2D+l] = ret2[l*MaxN_BaseFunctions2D+k];
  
        RegisterProlongationMatrix2D(Coarse, refine, Fine, j, ret);
      } // endfor j
    }

    ret = ProlongationMatrix2D[Coarse][refine][Fine][childnumber]; 

    RefCell->Derefine();
    delete (TGridCell *)RefCell;
    delete [] ret2;
  }

  return ret;
}

double *TFEDatabase2D::GetRestrictionMatrix2D (FE2D parent, 
    Refinements refine, FE2D child, int childnumber)
{ 
  double *ret, *ret2;
  int i,j,k,l, l1, l2;
  int N_Coarse, N_Fine, N_Points, N_Children;
  double AllPointValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double PointValues[MaxN_PointsForNodal2D];
  TFE2D *CoarseElement, *FineElement;
  TRefDesc *RefDesc;
  TBaseFunct2D *BaseFunctions, *FineBF;
  BaseFunct2D Coarse, Fine;
  TBaseCell *RefCell, *cell;
  TNodalFunctional2D *nf;
  RefTrans2D F_K;
  TRefTrans2D *rt;

  double G[MaxN_BaseFunctions2D][MaxN_BaseFunctions2D];
  double Gret[MaxN_BaseFunctions2D][MaxN_BaseFunctions2D];
  double R[MaxN_BaseFunctions2D][MaxN_BaseFunctions2D];
  TQuadFormula2D *Formula;
  QuadFormula2D QuadFormula;
  QuadFormula1D LineQuadFormula;
  double **CoarseBFData, **FineBFData, *PointData;
  double *FinePointData;
  int N_QuadPoints;
  double *xi, *eta, *weights, sum, w;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  int LDA=MaxN_BaseFunctions2D;

  CoarseElement = TFEDatabase2D::GetFE2D(parent);
  Coarse = CoarseElement->GetBaseFunct2D_ID();
  FineElement = TFEDatabase2D::GetFE2D(child);
  Fine = FineElement->GetBaseFunct2D_ID();

  ret = RestrictionMatrix2D[Coarse][refine][Fine][childnumber]; 

  if(ret == NULL)
  {
    // restriction matrix was not generated yet

    BaseFunctions = CoarseElement->GetBaseFunct2D();
    N_Coarse = BaseFunctions->GetDimension();

    memset(G, 0, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D*
                 SizeOfDouble);

    // build matrix G, gij = (uiH, ujH)
    TQuadFormula2D::FindLocalQuadFormula2D
            (1, &parent, LineQuadFormula, QuadFormula);
    Formula = GetQuadFormula2D(QuadFormula);
    Formula->GetFormulaData(N_QuadPoints, weights, xi, eta);

    TFEDatabase2D::GetBaseFunct2D(Coarse)->MakeRefElementData(QuadFormula);
    CoarseBFData = GetRefElementValues(Coarse, QuadFormula, D00);
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

    memcpy(Gret, G, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D
                    *SizeOfDouble);

    if(refine == NoRef)
    {
      memset(R, 0, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D*
                   SizeOfDouble);

      BaseFunctions = FineElement->GetBaseFunct2D();
      N_Fine = BaseFunctions->GetDimension();

      // build matrix R, rij = (uiH, ujh)
      TQuadFormula2D::FindLocalQuadFormula2D
              (1, &child, LineQuadFormula, QuadFormula);
      Formula = GetQuadFormula2D(QuadFormula);
      Formula->GetFormulaData(N_QuadPoints, weights, xi, eta);
      
      TFEDatabase2D::GetBaseFunct2D(Fine)->MakeRefElementData(QuadFormula);
      FineBFData = GetRefElementValues(Fine, QuadFormula, D00);
      TFEDatabase2D::GetBaseFunct2D(Coarse)->MakeRefElementData(QuadFormula);
      CoarseBFData = GetRefElementValues(Coarse, QuadFormula, D00);
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

      memcpy(G, Gret, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D
                      *SizeOfDouble);

      // determine inv(G)*R
      SolveMultipleSystems((double *)G, (double *)R, 
                           N_Coarse, LDA, LDA, N_Fine);

      ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
      for(l1=0;l1<N_Coarse;l1++)
        for(l2=0;l2<N_Fine;l2++)
          ret[l1*MaxN_BaseFunctions2D+l2] = R[l2][l1]; 
  
      RegisterRestrictionMatrix2D(Coarse, refine, Fine, childnumber, ret);
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
        memset(R, 0, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D*
                     SizeOfDouble);
  
        cell = RefCell->GetChild(j);
        FineElement = TFEDatabase2D::GetFE2D(child);
        Fine = FineElement->GetBaseFunct2D_ID();
        FineBF = FineElement->GetBaseFunct2D();
        N_Fine = FineBF->GetDimension();
  
        TQuadFormula2D::FindLocalQuadFormula2D
                (1, &child, LineQuadFormula, QuadFormula);
        Formula = GetQuadFormula2D(QuadFormula);
        Formula->GetFormulaData(N_QuadPoints, weights, xi, eta);

        TFEDatabase2D::GetBaseFunct2D(Fine)->MakeRefElementData(QuadFormula);
        FineBFData = GetRefElementValues(Fine, QuadFormula, D00);
  
        F_K = FineElement->GetRefTransID();
  
        switch(F_K)
        {
          case QuadAffin:
          case QuadBilinear:
          case QuadIsoparametric:
            rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
            ((TQuadAffin *)rt)->SetCell(RefCell->GetChild(j));
            F_K = QuadAffin;
            break;
          case TriaAffin:
          case TriaIsoparametric:
            rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
            ((TTriaAffin *)rt)->SetCell(RefCell->GetChild(j));
            break;
        }
        TFEDatabase2D::GetOrigFromRef(F_K ,N_QuadPoints, xi, eta,
                                    X, Y, AbsDetjk);
  
        for(k=0;k<N_QuadPoints;k++)
          BaseFunctions->GetDerivatives(D00, X[k], Y[k], AllPointValues[k]);

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
       
        memcpy(G, Gret, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D
                        *SizeOfDouble);

        // determine inv(G)*R
        SolveMultipleSystems((double *)G, (double *)R, N_Coarse, 
                             LDA, LDA, N_Fine);
  
        ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
        // change from column to row storage
        for(l1=0;l1<N_Coarse;l1++)
          for(l2=0;l2<N_Fine;l2++)
            ret[l1*MaxN_BaseFunctions2D+l2] = R[l2][l1]; 
        RegisterRestrictionMatrix2D(Coarse, refine, Fine, j, ret);
      } // endfor j

      RefCell->Derefine();
      delete (TGridCell *)RefCell;

    }

    ret = RestrictionMatrix2D[Coarse][refine][Fine][childnumber]; 

  }

  return ret;
}

QuadFormula2D TFEDatabase2D::GetQFFromDegree(int accuracy, 
                                             BF2DRefElements RefElem)
{
  switch(RefElem)
  {
    case BFUnitSquare:
      return GetQFQuadFromDegree(accuracy);
      break;
    case BFUnitTriangle:
      return GetQFTriaFromDegree(accuracy);
      break;
    default:
      ErrMsg("unknown reference element");
      exit(0);
  }
};
