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
   
// this file contains a lot of enumerations used by object classes

#ifndef __ENUMERATIONS__
#define __ENUMERATIONS__

// ======================================================================
// general definitions
// ======================================================================
#define MAXDEGREE 24

// ======================================================================
// definitions for 1D objects
// ======================================================================

#define N_MultiIndices1D 3
enum MultiIndex1D { D0, D1, D2 };

#define N_FEs1D 7
enum FE1D { C_P0_1D_L_A, C_P1_1D_L_A, C_P2_1D_L_A, C_P3_1D_L_A,
            N_P0_1D_L_A, D_P1_1D_L_A, D_P2_1D_L_A };

#define N_FEDescs1D 7
enum FEDesc1D { FE_C_L_P0_1D, FE_C_L_P1_1D, FE_C_L_P2_1D, FE_C_L_P3_1D,
                FE_N_L_P0_1D, FE_D_L_P1_1D, FE_D_L_P2_1D };

#define MaxN_BaseFunctions1D 4
#define N_BaseFuncts1D 6
enum BaseFunct1D { BF_C_L_P0_1D, BF_C_L_P1_1D, BF_C_L_P2_1D, BF_C_L_P3_1D,
                   BF_D_L_P1_1D, BF_D_L_P2_1D};

#define MaxN_PointsForNodal1D 4
#define N_NodalFunctionals1D 6
enum NodalFunctional1D { NF_C_L_P0_1D, NF_C_L_P1_1D, NF_C_L_P2_1D,
                         NF_C_L_P3_1D, NF_D_L_P1_1D, NF_D_L_P2_1D };

#define N_RefTrans1D 1
enum RefTrans1D { LineAffin };

#define MaxN_QuadPoints_1D 16
#define N_QuadFormulas_1D 18
enum QuadFormula1D { 
                      Gauss1Line, Gauss2Line, Gauss3Line,
                      Gauss4Line, Gauss5Line, Gauss6Line,
                      Gauss7Line, Gauss8Line, Gauss9Line,
                      Gauss10Line, Gauss11Line, Gauss12Line,
                      Gauss2W1Line, Gauss4W1Line, Gauss6W1Line,
                      Gauss8W1Line,
                      Gauss16W2Line, Dummy
                    };

// ======================================================================
// definitions for 2D objects
// ======================================================================

#define N_MultiIndices2D 6
enum MultiIndex2D { D00, 
                    D10, D01, 
                    D20, D11, D02 };

#define N_FEs2D 171
enum FE2D { C_P00_2D_T_A, C_P0_2D_T_A, C_P1_2D_T_A, C_P2_2D_T_A, C_P3_2D_T_A,
            C_P4_2D_T_A, C_P5_2D_T_A, C_P6_2D_T_A, C_P7_2D_T_A,
            C_P8_2D_T_A, C_P9_2D_T_A,
            N_P1_2D_T_A, // 12
            C_Q00_2D_Q_A, C_Q0_2D_Q_A, C_Q1_2D_Q_A, C_Q2_2D_Q_A, C_Q3_2D_Q_A,
            C_Q4_2D_Q_A, C_Q5_2D_Q_A, C_Q6_2D_Q_A, C_Q7_2D_Q_A,
            C_Q8_2D_Q_A, C_Q9_2D_Q_A,
            N_Q1_2D_Q_A, // 24
            C_Q00_2D_Q_M, C_Q0_2D_Q_M, C_Q1_2D_Q_M, C_Q2_2D_Q_M, C_Q3_2D_Q_M,
            C_Q4_2D_Q_M, C_Q5_2D_Q_M, C_Q6_2D_Q_M, C_Q7_2D_Q_M,
            C_Q8_2D_Q_M, C_Q9_2D_Q_M,
            N_Q1_2D_Q_M, // 36
            D_P1_2D_Q_A, D_P2_2D_Q_A, D_P3_2D_Q_A,
            D_P1_2D_Q_M, D_P2_2D_Q_M, D_P3_2D_Q_M,
            C_B2_2D_T_A, C_B3_2D_T_A, C_SV2_2D_T_A,
            D_P1_2D_T_A, D_P2_2D_T_A, D_SV1_2D_T_A, // 48

            N_Q2_2D_Q_A, N_Q3_2D_Q_A, N_Q4_2D_Q_A, N_Q5_2D_Q_A,
            N_Q2_2D_Q_M, N_Q3_2D_Q_M, N_Q4_2D_Q_M, N_Q5_2D_Q_M,

            D_P4_2D_Q_A, D_P5_2D_Q_A, D_P6_2D_Q_A, D_P7_2D_Q_A,
            D_P4_2D_Q_M, D_P5_2D_Q_M, D_P6_2D_Q_M, D_P7_2D_Q_M, // 64

            N_P1MOD_2D_T_A,

            N_P2_2D_T_A, N_P3_2D_T_A, N_P4_2D_T_A, N_P5_2D_T_A,

            D_P3_2D_T_A, D_P4_2D_T_A,

            C_P1MINI_2D_T_A,

            B_IB2_2D_Q_A, B_IB2_2D_Q_M, // 74

            D_Q1_2D_Q_A, D_Q2_2D_Q_A, D_Q3_2D_Q_A, D_Q4_2D_Q_A,
            D_Q1_2D_Q_M, D_Q2_2D_Q_M, D_Q3_2D_Q_M, D_Q4_2D_Q_M,

            C_B4_2D_T_A,

            D_D2_2D_Q_A, D_D2_2D_Q_M, // 85

            C_UL1_2D_Q_A, C_UL2_2D_Q_A, C_UL3_2D_Q_A, C_UL4_2D_Q_A, C_UL5_2D_Q_A,
            C_UL1_2D_Q_M, C_UL2_2D_Q_M, C_UL3_2D_Q_M, C_UL4_2D_Q_M, C_UL5_2D_Q_M,

            C_UL1_2D_T_A, C_UL2_2D_T_A, C_UL3_2D_T_A, C_UL4_2D_T_A, C_UL5_2D_T_A, // 100

            C_UL2S_2D_Q_A, C_UL3S_2D_Q_A, C_UL4S_2D_Q_A, C_UL5S_2D_Q_A,
            C_UL6S_2D_Q_A, C_UL7S_2D_Q_A, C_UL8S_2D_Q_A, C_UL9S_2D_Q_A,
            C_UL2S_2D_Q_M, C_UL3S_2D_Q_M, C_UL4S_2D_Q_M, C_UL5S_2D_Q_M,
            C_UL6S_2D_Q_M, C_UL7S_2D_Q_M, C_UL8S_2D_Q_M, C_UL9S_2D_Q_M, // 116

            C_UL2SE_2D_Q_A, C_UL3SE_2D_Q_A, C_UL4SE_2D_Q_A, C_UL5SE_2D_Q_A,
            C_UL6SE_2D_Q_A, C_UL7SE_2D_Q_A, C_UL8SE_2D_Q_A, C_UL9SE_2D_Q_A,
            C_UL2SE_2D_Q_M, C_UL3SE_2D_Q_M, C_UL4SE_2D_Q_M, C_UL5SE_2D_Q_M,
            C_UL6SE_2D_Q_M, C_UL7SE_2D_Q_M, C_UL8SE_2D_Q_M, C_UL9SE_2D_Q_M, // 132

            C_M2_2D_Q_A, C_M3_2D_Q_A, C_M4_2D_Q_A, C_M5_2D_Q_A,
            C_M6_2D_Q_A, C_M7_2D_Q_A, C_M8_2D_Q_A, C_M9_2D_Q_A,
            C_M2_2D_Q_M, C_M3_2D_Q_M, C_M4_2D_Q_M, C_M5_2D_Q_M,
            C_M6_2D_Q_M, C_M7_2D_Q_M, C_M8_2D_Q_M, C_M9_2D_Q_M,
            
            C_EL1_2D_Q_A, C_EL1_2D_Q_M, // 150
            
            N_RT0_2D_Q_A, N_RT0_2D_Q_M, N_RT1_2D_Q_A, N_RT1_2D_Q_M,
            N_RT2_2D_Q_A, N_RT2_2D_Q_M, N_RT3_2D_Q_A, N_RT3_2D_Q_M,
            N_RT0_2D_T_A, N_RT1_2D_T_A, N_RT2_2D_T_A, N_RT3_2D_T_A, // 162
            
            N_BDM1_2D_Q_A, N_BDM2_2D_Q_A, N_BDM3_2D_Q_A,
            N_BDM1_2D_Q_M, N_BDM2_2D_Q_M, N_BDM3_2D_Q_M,
            N_BDM1_2D_T_A, N_BDM2_2D_T_A, N_BDM3_2D_T_A // 171
          };

#define N_FEDescs2D 105
enum FEDesc2D { FE_C_T_P00_2D, FE_C_T_P0_2D, FE_C_T_P1_2D, FE_C_T_P2_2D, 
                FE_C_T_P3_2D, FE_C_T_P4_2D, FE_C_T_P5_2D, FE_C_T_P6_2D,
                FE_C_T_P7_2D, FE_C_T_P8_2D, FE_C_T_P9_2D, FE_N_T_P1_2D,
                FE_C_Q_Q00_2D, FE_C_Q_Q0_2D, FE_C_Q_Q1_2D, FE_C_Q_Q2_2D,
                FE_C_Q_Q3_2D, FE_C_Q_Q4_2D, FE_C_Q_Q5_2D, FE_C_Q_Q6_2D, 
                FE_C_Q_Q7_2D, FE_C_Q_Q8_2D, FE_C_Q_Q9_2D, FE_N_Q_Q1_2D, // 24
                FE_D_Q_P1_2D, FE_D_Q_P2_2D, FE_D_Q_P3_2D,
                FE_C_T_B2_2D, FE_C_T_B3_2D, FE_C_T_SV2_2D,
                FE_D_T_P1_2D, FE_D_T_P2_2D, FE_D_T_SV1_2D, // 33

                FE_N_Q_Q2_2D, FE_N_Q_Q3_2D, FE_N_Q_Q4_2D, FE_N_Q_Q5_2D,
                FE_D_Q_P4_2D, FE_D_Q_P5_2D, FE_D_Q_P6_2D, FE_D_Q_P7_2D,
                
                FE_N_T_P1MOD_2D, FE_N_T_P2_2D, FE_N_T_P3_2D, FE_N_T_P4_2D, 
                FE_N_T_P5_2D, FE_D_T_P3_2D, FE_D_T_P4_2D, FE_C_T_P1MINI_2D, //49

                FE_B_Q_IB2_2D, 
                FE_D_Q_Q1_2D, FE_D_Q_Q2_2D, FE_D_Q_Q3_2D, FE_D_Q_Q4_2D, 
                FE_C_T_B4_2D, FE_D_Q_D2_2D, // 56

                FE_C_T_UL1_2D, FE_C_T_UL2_2D, FE_C_T_UL3_2D, FE_C_T_UL4_2D, 
                FE_C_T_UL5_2D, FE_C_Q_UL1_2D, FE_C_Q_UL2_2D, FE_C_Q_UL3_2D, 
                FE_C_Q_UL4_2D, FE_C_Q_UL5_2D, // 66

                FE_C_Q_UL2S_2D, FE_C_Q_UL3S_2D, FE_C_Q_UL4S_2D, FE_C_Q_UL5S_2D,
                FE_C_Q_UL6S_2D, FE_C_Q_UL7S_2D, FE_C_Q_UL8S_2D, FE_C_Q_UL9S_2D,

                FE_C_Q_UL2SE_2D, FE_C_Q_UL3SE_2D, FE_C_Q_UL4SE_2D, FE_C_Q_UL5SE_2D,
                FE_C_Q_UL6SE_2D, FE_C_Q_UL7SE_2D, FE_C_Q_UL8SE_2D, FE_C_Q_UL9SE_2D, // 82

                FE_C_Q_M2_2D, FE_C_Q_M3_2D, FE_C_Q_M4_2D, FE_C_Q_M5_2D,
                FE_C_Q_M6_2D, FE_C_Q_M7_2D, FE_C_Q_M8_2D, FE_C_Q_M9_2D,
                
                FE_C_Q_EL1_2D, // 91
                
                FE_N_Q_RT0_2D, FE_N_Q_RT1_2D, FE_N_Q_RT2_2D, FE_N_Q_RT3_2D,
                FE_N_T_RT0_2D, FE_N_T_RT1_2D, FE_N_T_RT2_2D, FE_N_T_RT3_2D,
                FE_N_Q_BDM1_2D, FE_N_Q_BDM2_2D, FE_N_Q_BDM3_2D,
                FE_N_T_BDM1_2D, FE_N_T_BDM2_2D, FE_N_T_BDM3_2D //105
              };

#define MaxN_BaseFunctions2D 100

#define N_BaseFuncts2D 105
enum BaseFunct2D { BF_C_T_P00_2D, BF_C_T_P0_2D, BF_C_T_P1_2D, BF_C_T_P2_2D, BF_C_T_P3_2D,
                   BF_N_T_P1_2D,
                   BF_C_Q_Q00_2D, BF_C_Q_Q0_2D, BF_C_Q_Q1_2D, BF_C_Q_Q2_2D, BF_C_Q_Q3_2D,
                   BF_N_Q_Q1_2D,
                   BF_D_Q_P1_2D, BF_D_Q_P2_2D, BF_D_Q_P3_2D, // 15
                   BF_C_T_P4_2D, BF_C_T_P5_2D, BF_C_T_P6_2D, BF_C_T_P7_2D,
                   BF_C_T_P8_2D, BF_C_T_P9_2D,
                   BF_C_Q_Q4_2D, BF_C_Q_Q5_2D, BF_C_Q_Q6_2D, BF_C_Q_Q7_2D,
                   BF_C_Q_Q8_2D, BF_C_Q_Q9_2D,
                   BF_C_T_B2_2D, BF_C_T_B3_2D, BF_C_T_SV2_2D,
                   BF_D_T_P1_2D, BF_D_T_P2_2D, BF_D_T_SV1_2D, // 33

                   BF_N_Q_Q2_2D, BF_N_Q_Q3_2D, BF_N_Q_Q4_2D, BF_N_Q_Q5_2D,

                   BF_D_Q_P4_2D, BF_D_Q_P5_2D, BF_D_Q_P6_2D, BF_D_Q_P7_2D,

                   BF_N_T_P1MOD_2D, // 42

                   BF_N_T_P2_2D, BF_N_T_P3_2D, BF_N_T_P4_2D, BF_N_T_P5_2D,

                   BF_D_T_P3_2D, BF_D_T_P4_2D,

                   BF_C_T_P1MINI_2D,

                   BF_B_Q_IB2_2D,

                   BF_D_Q_Q1_2D, BF_D_Q_Q2_2D, BF_D_Q_Q3_2D, BF_D_Q_Q4_2D,

                   BF_C_T_B4_2D,
                   
                   BF_D_Q_D2_2D, // 56

                   BF_C_Q_UL1_2D, BF_C_Q_UL2_2D, BF_C_Q_UL3_2D, BF_C_Q_UL4_2D, 
                   BF_C_Q_UL5_2D,

                   BF_C_T_UL1_2D, BF_C_T_UL2_2D, BF_C_T_UL3_2D, BF_C_T_UL4_2D, 
                   BF_C_T_UL5_2D, // 66

                   BF_C_Q_UL2S_2D, BF_C_Q_UL3S_2D, BF_C_Q_UL4S_2D, BF_C_Q_UL5S_2D,
                   BF_C_Q_UL6S_2D, BF_C_Q_UL7S_2D, BF_C_Q_UL8S_2D, BF_C_Q_UL9S_2D,

                   BF_C_Q_UL2SE_2D, BF_C_Q_UL3SE_2D, BF_C_Q_UL4SE_2D, BF_C_Q_UL5SE_2D,
                   BF_C_Q_UL6SE_2D, BF_C_Q_UL7SE_2D, BF_C_Q_UL8SE_2D, BF_C_Q_UL9SE_2D,

                   BF_C_Q_M2_2D, BF_C_Q_M3_2D, BF_C_Q_M4_2D, BF_C_Q_M5_2D,
                   BF_C_Q_M6_2D, BF_C_Q_M7_2D, BF_C_Q_M8_2D, BF_C_Q_M9_2D, // 90
                   
                   BF_C_Q_EL1_2D,
                   
                   BF_N_Q_RT0_2D, BF_N_Q_RT1_2D, BF_N_Q_RT2_2D, BF_N_Q_RT3_2D,
                   BF_N_T_RT0_2D, BF_N_T_RT1_2D, BF_N_T_RT2_2D, BF_N_T_RT3_2D,
                   
                   BF_N_Q_BDM1_2D, BF_N_Q_BDM2_2D, BF_N_Q_BDM3_2D,
                   BF_N_T_BDM1_2D, BF_N_T_BDM2_2D, BF_N_T_BDM3_2D //105
                 };

enum BF2DRefElements { BFUnitTriangle, BFUnitSquare };

#define MaxN_PointsForNodal2D 100
#define N_NodalFunctionals2D 106
enum NodalFunctional2D {
                NF_C_T_P00_2D, NF_C_T_P0_2D, NF_C_T_P1_2D, NF_C_T_P2_2D, NF_C_T_P3_2D,
                NF_C_T_P4_2D, NF_C_T_P5_2D, NF_C_T_P6_2D, NF_C_T_P7_2D,
                NF_C_T_P8_2D, NF_C_T_P9_2D,
                NF_N_T_P1_2D, // 12
                NF_C_Q_Q00_2D, NF_C_Q_Q0_2D, NF_C_Q_Q1_2D, NF_C_Q_Q2_2D, NF_C_Q_Q3_2D,
                NF_C_Q_Q4_2D, NF_C_Q_Q5_2D, NF_C_Q_Q6_2D, NF_C_Q_Q7_2D,
                NF_C_Q_Q8_2D, NF_C_Q_Q9_2D,
                NF_N_Q_Q1_2D, // 24
                NF_D_Q_P1_2D, NF_D_Q_P2_2D, NF_D_Q_P3_2D,
                NF_C_T_B2_2D, NF_C_T_B3_2D, NF_C_T_SV2_2D,
                NF_D_T_P1_2D, NF_D_T_P2_2D, NF_D_T_SV1_2D, // 33

                NF_N_Q_Q2_2D, NF_N_Q_Q3_2D, NF_N_Q_Q4_2D, NF_N_Q_Q5_2D,

                NF_D_Q_P4_2D, NF_D_Q_P5_2D, NF_D_Q_P6_2D, NF_D_Q_P7_2D,

                NF_N_T_P1MOD_2D, // 42

                NF_N_T_P2_2D, NF_N_T_P3_2D, NF_N_T_P4_2D, NF_N_T_P5_2D,

                NF_D_T_P3_2D, NF_D_T_P4_2D,

                NF_C_T_P1MINI_2D,

                NF_B_Q_IB2_2D,

                NF_S_Q_Q2_2D,

                NF_D_Q_Q1_2D, NF_D_Q_Q2_2D, NF_D_Q_Q3_2D, NF_D_Q_Q4_2D,

                NF_C_T_B4_2D,

                NF_D_Q_D2_2D, // 57

                NF_C_Q_UL1_2D, NF_C_Q_UL2_2D, NF_C_Q_UL3_2D, NF_C_Q_UL4_2D, 
                NF_C_Q_UL5_2D,

                NF_C_T_UL1_2D, NF_C_T_UL2_2D, NF_C_T_UL3_2D, NF_C_T_UL4_2D, 
                NF_C_T_UL5_2D, // 67

                NF_C_Q_UL2S_2D, NF_C_Q_UL3S_2D, NF_C_Q_UL4S_2D, NF_C_Q_UL5S_2D,
                NF_C_Q_UL6S_2D, NF_C_Q_UL7S_2D, NF_C_Q_UL8S_2D, NF_C_Q_UL9S_2D,

                NF_C_Q_UL2SE_2D, NF_C_Q_UL3SE_2D, NF_C_Q_UL4SE_2D, NF_C_Q_UL5SE_2D,
                NF_C_Q_UL6SE_2D, NF_C_Q_UL7SE_2D, NF_C_Q_UL8SE_2D, NF_C_Q_UL9SE_2D,

                NF_C_Q_M2_2D, NF_C_Q_M3_2D, NF_C_Q_M4_2D, NF_C_Q_M5_2D,
                NF_C_Q_M6_2D, NF_C_Q_M7_2D, NF_C_Q_M8_2D, NF_C_Q_M9_2D, // 91
                
                NF_C_Q_EL1_2D,
                
                NF_N_Q_RT0_2D, NF_N_Q_RT1_2D, NF_N_Q_RT2_2D, NF_N_Q_RT3_2D,
                NF_N_T_RT0_2D, NF_N_T_RT1_2D, NF_N_T_RT2_2D, NF_N_T_RT3_2D,
                
                NF_N_Q_BDM1_2D, NF_N_Q_BDM2_2D, NF_N_Q_BDM3_2D,
                NF_N_T_BDM1_2D, NF_N_T_BDM2_2D, NF_N_T_BDM3_2D // 106
              };

#define N_HNDescs 70
enum HNDesc { HN_C_P1_2D_0, 
              HN_C_P2_2D_0, HN_C_P2_2D_1,
              HN_C_P3_2D_0, HN_C_P3_2D_1, HN_C_P3_2D_2,
              HN_C_P4_2D_0, HN_C_P4_2D_1, HN_C_P4_2D_2, HN_C_P4_2D_3,
              HN_C_P5_2D_0, HN_C_P5_2D_1, HN_C_P5_2D_2, HN_C_P5_2D_3,
              HN_C_P5_2D_4,
              HN_C_P6_2D_0, HN_C_P6_2D_1, HN_C_P6_2D_2, HN_C_P6_2D_3,
              HN_C_P6_2D_4, HN_C_P6_2D_5,
              HN_C_P7_2D_0, HN_C_P7_2D_1, HN_C_P7_2D_2, HN_C_P7_2D_3,
              HN_C_P7_2D_4, HN_C_P7_2D_5, HN_C_P7_2D_6,
              HN_C_P8_2D_0, HN_C_P8_2D_1, HN_C_P8_2D_2, HN_C_P8_2D_3,
              HN_C_P8_2D_4, HN_C_P8_2D_5, HN_C_P8_2D_6, HN_C_P8_2D_7,
              HN_C_P9_2D_0, HN_C_P9_2D_1, HN_C_P9_2D_2, HN_C_P9_2D_3,
              HN_C_P9_2D_4, HN_C_P9_2D_5, HN_C_P9_2D_6, HN_C_P9_2D_7,
              HN_C_P9_2D_8,
              HN_N_P1_2D_0, 
              HN_N_P2_2D_0, HN_N_P3_2D_0, 
              HN_N_P4_2D_0, HN_N_P5_2D_0, 
              HN_C_Q1_3D_E, HN_C_Q1_3D_F,
              HN_C_Q2_3D_E, HN_C_Q2_3D_F,
              HN_C_Q3_3D_1, HN_C_Q3_3D_2, HN_C_Q3_3D_3, HN_C_Q3_3D_4,
              HN_C_Q3_3D_5,
              HN_C_P1_3D_E,
              HN_C_P2_3D_E, HN_C_P2_3D_F,
              HN_C_P3_3D_E, HN_C_P3_3D_M,
              HN_C_P3_3D_F, HN_C_P3_3D_G,
              HN_N_P1_3D_E,
              HN_N_P2_3D_0, HN_N_P2_3D_1, HN_N_P2_3D_2
            };

#define N_RefTrans2D 5
enum RefTrans2D { TriaAffin, QuadAffin, QuadBilinear, 
                  TriaIsoparametric, QuadIsoparametric };

#define MaxN_QuadPoints_2D 81
#define N_QuadFormulas_2D 22
enum QuadFormula2D { BaryCenterTria, MidPointTria, SevenPointTria,
                     Gauss3Tria, VertexTria, Degree8Tria,  Degree9Tria,
                     Degree11Tria, Degree19Tria,
                     Gauss2Quad, Gauss3Quad, Gauss4Quad, 
                     Gauss5Quad, Gauss6Quad, Gauss7Quad,
                     Gauss8Quad, Gauss9Quad, VertexQuad,
                     SimpsonQuad,
                     CompGauss3Tria, CompGauss4Tria,
                     Gauss_Degree8Tria};

enum SpaceType { DiscP_USpace, DiscP_PSpace, ContP_USpace, ContP_PSpace,
                 Non_USpace };

// ======================================================================
// definitions for 3D objects
// ======================================================================

#define N_MultiIndices3D 10
enum MultiIndex3D { D000, 
                    D100, D010, D001, 
                    D200, D110, D101, D020, D011, D002 };

#define MaxN_QuadPoints_3D 729
#define N_QuadFormulas_3D 19
enum QuadFormula3D { BaryCenterTetra, VertexTetra, P2Tetra, P4Tetra,
                     P5Tetra, P8Tetra,
                     VertexHexa, Gauss2Hexa, Gauss3Hexa, Gauss4Hexa,
                     Gauss5Hexa, Gauss6Hexa, Gauss7Hexa, Gauss8Hexa,
                     Gauss9Hexa, VerticesAndOrigin, VerticesAndOrigin15,
                     VerticesAndOrigin57, Degree7_Points38};

#define N_FEs3D 76
enum FE3D { C_P00_3D_T_A, C_P0_3D_T_A, C_P1_3D_T_A, C_P2_3D_T_A, C_P3_3D_T_A,
            N_P1_3D_T_A,
            C_Q00_3D_H_A, C_Q0_3D_H_A, C_Q1_3D_H_A, C_Q2_3D_H_A, C_Q3_3D_H_A,
            C_Q4_3D_H_A,
            C_Q00_3D_H_M, C_Q0_3D_H_M, C_Q1_3D_H_M, C_Q2_3D_H_M, C_Q3_3D_H_M,
            C_Q4_3D_H_M, // 18
            N_Q1_3D_H_A,
            N_Q1_3D_H_M,
            D_P1_3D_H_A, D_P2_3D_H_A, D_P3_3D_H_A,
            D_P1_3D_H_M, D_P2_3D_H_M, D_P3_3D_H_M,
            
            D_Q1_3D_H_A, D_Q1_3D_H_M, D_Q2_3D_H_A, D_Q2_3D_H_M,
            
            B_IB2_3D_H_A, B_IB2_3D_H_M, // 32

            N_P2_3D_T_A, N_P3_3D_T_A, N_P4_3D_T_A, N_P5_3D_T_A,
            N_Q2_3D_H_A, N_Q3_3D_H_A, N_Q4_3D_H_A, N_Q5_3D_H_A,
            N_Q2_3D_H_M, N_Q3_3D_H_M, N_Q4_3D_H_M, N_Q5_3D_H_M,

            C_B2_3D_T_A, D_P1_3D_T_A, D_P2_3D_T_A, D_P3_3D_T_A, // 48
            
            C_UL1_3D_H_A, C_UL2_3D_H_A, C_UL3_3D_H_A, C_UL4_3D_H_A, C_UL5_3D_H_A,
            C_UL1_3D_H_M, C_UL2_3D_H_M, C_UL3_3D_H_M, C_UL4_3D_H_M, C_UL5_3D_H_M,
            C_UL1_3D_T_A, C_UL2_3D_T_A, C_UL3_3D_T_A, C_UL4_3D_T_A, C_UL5_3D_T_A, // 63
            
            N_RT0_3D_T_A, N_RT1_3D_T_A, N_RT2_3D_T_A, N_RT3_3D_T_A,
            N_RT0_3D_H_A, N_RT1_3D_H_A, N_RT2_3D_H_A,
            N_BDDF1_3D_T_A, N_BDDF2_3D_T_A, N_BDDF3_3D_T_A,
            N_BDDF1_3D_H_A, N_BDDF2_3D_H_A, N_BDDF3_3D_H_A // 76
          };

#define N_FEDescs3D 54
enum FEDesc3D { FE_C_T_P00_3D, FE_C_T_P0_3D, FE_C_T_P1_3D, FE_C_T_P2_3D, FE_C_T_P3_3D,
                FE_N_T_P1_3D,
                FE_C_H_Q00_3D, FE_C_H_Q0_3D, FE_C_H_Q1_3D, FE_C_H_Q2_3D, FE_C_H_Q3_3D,
                FE_C_H_Q4_3D,
                FE_N_H_Q1_3D,
                FE_D_H_P1_3D, FE_D_H_P2_3D, FE_D_H_P3_3D, // 16
                
                FE_D_H_Q1_3D,  FE_D_H_Q2_3D,

                FE_B_H_IB2_3D,

                FE_N_T_P2_3D, FE_N_T_P3_3D, FE_N_T_P4_3D, FE_N_T_P5_3D,
                FE_N_H_Q2_3D, FE_N_H_Q3_3D, FE_N_H_Q4_3D, FE_N_H_Q5_3D,

                FE_C_T_B2_3D, FE_D_T_P1_3D, FE_D_T_P2_3D, FE_D_T_P3_3D, // 31
                
                FE_C_H_UL1_3D, FE_C_H_UL2_3D, FE_C_H_UL3_3D, FE_C_H_UL4_3D, 
                FE_C_H_UL5_3D, 
                FE_C_T_UL1_3D, FE_C_T_UL2_3D, FE_C_T_UL3_3D, FE_C_T_UL4_3D, 
                FE_C_T_UL5_3D, // 41
                
                FE_N_T_RT0_3D, FE_N_T_RT1_3D, FE_N_T_RT2_3D, FE_N_T_RT3_3D,
                FE_N_H_RT0_3D, FE_N_H_RT1_3D, FE_N_H_RT2_3D,
                FE_N_T_BDDF1_3D, FE_N_T_BDDF2_3D, FE_N_T_BDDF3_3D,
                FE_N_H_BDDF1_3D, FE_N_H_BDDF2_3D, FE_N_H_BDDF3_3D // 54
              };

#define MaxN_BaseFunctions3D 125
#define N_BaseFuncts3D 54
enum BaseFunct3D { BF_C_T_P00_3D, BF_C_T_P0_3D, BF_C_T_P1_3D, BF_C_T_P2_3D, BF_C_T_P3_3D,
                   BF_N_T_P1_3D,
                   BF_C_H_Q00_3D, BF_C_H_Q0_3D, BF_C_H_Q1_3D, BF_C_H_Q2_3D, BF_C_H_Q3_3D,
                   BF_C_H_Q4_3D,
                   BF_N_H_Q1_3D,
                   BF_D_H_P1_3D, BF_D_H_P2_3D, BF_D_H_P3_3D, // 16
                   
                   BF_D_H_Q1_3D, BF_D_H_Q2_3D,
                   
                   BF_B_H_IB2_3D,

                   BF_N_T_P2_3D, BF_N_T_P3_3D, BF_N_T_P4_3D, BF_N_T_P5_3D,
                   BF_N_H_Q2_3D, BF_N_H_Q3_3D, BF_N_H_Q4_3D, BF_N_H_Q5_3D,

                   BF_C_T_B2_3D, BF_D_T_P1_3D, BF_D_T_P2_3D, BF_D_T_P3_3D, // 31
                   
                   BF_C_H_UL1_3D, BF_C_H_UL2_3D, BF_C_H_UL3_3D, BF_C_H_UL4_3D,
                   BF_C_H_UL5_3D, 
                   BF_C_T_UL1_3D, BF_C_T_UL2_3D, BF_C_T_UL3_3D, BF_C_T_UL4_3D,
                   BF_C_T_UL5_3D, // 41
                   
                   BF_N_T_RT0_3D, BF_N_T_RT1_3D, BF_N_T_RT2_3D, BF_N_T_RT3_3D,
                   BF_N_H_RT0_3D, BF_N_H_RT1_3D, BF_N_H_RT2_3D,
                   BF_N_T_BDDF1_3D, BF_N_T_BDDF2_3D, BF_N_T_BDDF3_3D,
                   BF_N_H_BDDF1_3D, BF_N_H_BDDF2_3D, BF_N_H_BDDF3_3D // 54
                 };

enum BF3DRefElements { BFUnitTetrahedron, BFUnitHexahedron };

#define MaxN_PointsForNodal3D 275
#define N_NodalFunctionals3D 55
enum NodalFunctional3D {
                NF_C_T_P00_3D, NF_C_T_P0_3D, NF_C_T_P1_3D, NF_C_T_P2_3D, NF_C_T_P3_3D,
                NF_N_T_P1_3D,
                NF_C_H_Q00_3D, NF_C_H_Q0_3D, NF_C_H_Q1_3D, NF_C_H_Q2_3D, NF_C_H_Q3_3D,
                NF_C_H_Q4_3D,
                NF_N_H_Q1_3D,
                NF_D_H_P1_3D, NF_D_H_P2_3D, NF_D_H_P3_3D, // 16
                
                NF_D_H_Q1_3D, NF_D_H_Q2_3D,

                NF_B_H_IB2_3D,

                NF_S_H_Q2_3D,

                NF_N_T_P2_3D, NF_N_T_P3_3D, NF_N_T_P4_3D, NF_N_T_P5_3D,
                NF_N_H_Q2_3D, NF_N_H_Q3_3D, NF_N_H_Q4_3D, NF_N_H_Q5_3D,

                NF_C_T_B2_3D, NF_D_T_P1_3D, NF_D_T_P2_3D, NF_D_T_P3_3D, // 32
                
                NF_C_H_UL1_3D, NF_C_H_UL2_3D, NF_C_H_UL3_3D, NF_C_H_UL4_3D, 
                NF_C_H_UL5_3D, 
                NF_C_T_UL1_3D, NF_C_T_UL2_3D, NF_C_T_UL3_3D, NF_C_T_UL4_3D,
                NF_C_T_UL5_3D, // 42
                
                NF_N_T_RT0_3D,  NF_N_T_RT1_3D, NF_N_T_RT2_3D, NF_N_T_RT3_3D,
                NF_N_H_RT0_3D, NF_N_H_RT1_3D, NF_N_H_RT2_3D,
                NF_N_T_BDDF1_3D, NF_N_T_BDDF2_3D, NF_N_T_BDDF3_3D,
                NF_N_H_BDDF1_3D, NF_N_H_BDDF2_3D, NF_N_H_BDDF3_3D // 55
              };

#define N_RefTrans3D 5
enum RefTrans3D { TetraAffin, TetraIsoparametric, 
                  HexaAffin, HexaTrilinear, HexaIsoparametric };

#endif
