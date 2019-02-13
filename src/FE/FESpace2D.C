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
// @(#)FESpace2D.C        1.18 06/27/00
// 
// Class:       TFESpace2D
// Purpose:     class for all 2D finite element spaces
//
// Author:      Gunar Matthies (04.11.97)
//
// History:     start of implementation 04.11.97 (Gunar Matthies)
//
//              split TFESpace into TFESpacexD (15.04.1998) Volker Behns
//              parallel methods  (Sashikumaar Ganesan) 09.09.09
// =======================================================================

#ifdef _MPI
#  include "mpi.h"
#endif

#include <DefineParams.h>
#include <Constants.h>
#include <FESpace2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp2D.h>
#include <FEDatabase2D.h>

#ifdef __MORTAR__
  #include <MortarBaseJoint.h>
#endif

#include <Database.h>
#include <InterfaceJoint.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <NodalFunctional2D.h>

#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>

/** Constructor */
TFESpace2D::TFESpace2D(TCollection *coll, char *name, char *description) :
     TFESpace(coll, name, description)
{
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  AllElements = NULL;
  ElementForShape = NULL;
  DGSpace = 0;
}

// =====================================================================
// rules: if a neighbour element is not in the collection, its clipboard
//        must be set to -1
// =====================================================================

/** constructor for building a space with elements of order k */
TFESpace2D::TFESpace2D(TCollection *coll, char *name, char *description, 
                       BoundCondFunct2D *BoundaryCondition, int ord,
                       TCollection *mortarcoll) :
     TFESpace(coll, name, description)
{
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  AllElements = NULL;
  DGSpace = 0;
#ifdef __MORTAR__
  MortarColl = mortarcoll;
#endif

  ElementForShape = new FE2D[N_SHAPES];

  // build ElementForShape array
  switch(ord)
  {
    case 0:  ElementForShape[Triangle] = C_P0_2D_T_A;
             ElementForShape[Quadrangle] = C_Q0_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q0_2D_Q_A;
             ElementForShape[Rectangle] = C_Q0_2D_Q_A;
             break;
    case 1:  ElementForShape[Triangle] = C_P1_2D_T_A;
             ElementForShape[Quadrangle] = C_Q1_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
             ElementForShape[Rectangle] = C_Q1_2D_Q_A;
             break;
    case 2:  ElementForShape[Triangle] = C_P2_2D_T_A;
             // ElementForShape[Triangle] = C_B2_2D_T_A;
             ElementForShape[Quadrangle] = C_Q2_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
             ElementForShape[Rectangle] = C_Q2_2D_Q_A;
             break;
    case 222:
      switch((int) (TDatabase::ParamDB->REACTOR_P11 + 0.5))
      {
        case 1:
        case 2:
          ElementForShape[Triangle] = C_B2_2D_T_A;
          break;
          
        default:
          ElementForShape[Triangle] = C_P2_2D_T_A;
      }
      ElementForShape[Quadrangle] = C_Q2_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
      ElementForShape[Rectangle] = C_Q2_2D_Q_A;
      break;
    case 3:  ElementForShape[Triangle] = C_P3_2D_T_A;
             // ElementForShape[Triangle] = C_B3_2D_T_A;
             ElementForShape[Quadrangle] = C_Q3_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q3_2D_Q_A;
             ElementForShape[Rectangle] = C_Q3_2D_Q_A;
             break;
    case 4:  ElementForShape[Triangle] = C_P4_2D_T_A;
             ElementForShape[Quadrangle] = C_Q4_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q4_2D_Q_A;
             ElementForShape[Rectangle] = C_Q4_2D_Q_A;
             break;
    case 5:  ElementForShape[Triangle] = C_P5_2D_T_A;
             ElementForShape[Quadrangle] = C_Q5_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q5_2D_Q_A;
             ElementForShape[Rectangle] = C_Q5_2D_Q_A;
             break;
    case 6:  ElementForShape[Triangle] = C_P6_2D_T_A;
             ElementForShape[Quadrangle] = C_Q6_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q6_2D_Q_A;
             ElementForShape[Rectangle] = C_Q6_2D_Q_A;
             break;
    case 7:  ElementForShape[Triangle] = C_P7_2D_T_A;
             ElementForShape[Quadrangle] = C_Q7_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q7_2D_Q_A;
             ElementForShape[Rectangle] = C_Q7_2D_Q_A;
             break;
    case 8:  ElementForShape[Triangle] = C_P8_2D_T_A;
             ElementForShape[Quadrangle] = C_Q8_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q8_2D_Q_A;
             ElementForShape[Rectangle] = C_Q8_2D_Q_A;
             break;
    case 9:  ElementForShape[Triangle] = C_P9_2D_T_A;
             ElementForShape[Quadrangle] = C_Q9_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q9_2D_Q_A;
             ElementForShape[Rectangle] = C_Q9_2D_Q_A;
             break;
    case -1: ElementForShape[Triangle] = N_P1_2D_T_A;
             ElementForShape[Quadrangle] = N_Q1_2D_Q_M;
             ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
             ElementForShape[Rectangle] = N_Q1_2D_Q_A;
             break;

    // P2/Q2 nonconforming
    case -2: ElementForShape[Triangle] = N_P2_2D_T_A;
             ElementForShape[Quadrangle] = N_Q2_2D_Q_M;
             ElementForShape[Parallelogram] = N_Q2_2D_Q_A;
             ElementForShape[Rectangle] = N_Q2_2D_Q_A;
             break;

    // P3/Q3 nonconforming
    case -3: ElementForShape[Triangle] = N_P3_2D_T_A;
             ElementForShape[Quadrangle] = N_Q3_2D_Q_M;
             ElementForShape[Parallelogram] = N_Q3_2D_Q_A;
             ElementForShape[Rectangle] = N_Q3_2D_Q_A;
             break;
    // P4/Q4 nonconforming
    case -4: ElementForShape[Triangle] = N_P4_2D_T_A;
             ElementForShape[Quadrangle] = N_Q4_2D_Q_M;
             ElementForShape[Parallelogram] = N_Q4_2D_Q_A;
             ElementForShape[Rectangle] = N_Q4_2D_Q_A;
             break;
    // P5/Q5 nonconforming
    case -5: ElementForShape[Triangle] = N_P5_2D_T_A;
             ElementForShape[Quadrangle] = N_Q5_2D_Q_M;
             ElementForShape[Parallelogram] = N_Q5_2D_Q_A;
             ElementForShape[Rectangle] = N_Q5_2D_Q_A;
             break;

    // P1mod
    case -101: 
             ElementForShape[Triangle] = N_P1MOD_2D_T_A;
             ElementForShape[Quadrangle] = N_Q1_2D_Q_M;
             ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
             ElementForShape[Rectangle] = N_Q1_2D_Q_A;
             cout << "P1MOD works only on triangles" << endl;
             break;

    // P1mini
    case 101: 
             ElementForShape[Triangle] = C_P1MINI_2D_T_A;
             ElementForShape[Quadrangle] = C_Q1_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
             ElementForShape[Rectangle] = C_Q1_2D_Q_A;
             cout << "P1MINI works only on triangles" << endl;
             break;
            
    //========LOCALPROJECTION=============
    // Q1+bubble*P0
    case 100: 
             ElementForShape[Triangle] = C_UL1_2D_T_A;
             ElementForShape[Quadrangle] = C_UL1_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL1_2D_Q_A;
             ElementForShape[Rectangle] = C_UL1_2D_Q_A;
             break;
    // Q2+bubble*P1
    case 201:  
             ElementForShape[Triangle] = C_UL2_2D_T_A;
             ElementForShape[Quadrangle] = C_UL2_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL2_2D_Q_A;
             ElementForShape[Rectangle] = C_UL2_2D_Q_A;
             break;
    // Q3+bubble*P2
    case 302: 
             ElementForShape[Triangle] = C_UL3_2D_T_A;
             ElementForShape[Quadrangle] = C_UL3_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL3_2D_Q_A;
             ElementForShape[Rectangle] = C_UL3_2D_Q_A;
             break;
    // Q4+bubble*P3
    case 403:  
             ElementForShape[Triangle] = C_UL4_2D_T_A;
             ElementForShape[Quadrangle] = C_UL4_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL4_2D_Q_A;
             ElementForShape[Rectangle] = C_UL4_2D_Q_A;
             break;
    // Q5+bubble*P4
    case 504:  
             ElementForShape[Triangle] = C_UL5_2D_T_A;
             ElementForShape[Quadrangle] = C_UL5_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL5_2D_Q_A;
             ElementForShape[Rectangle] = C_UL5_2D_Q_A;
             break; 

    case 200: // enriched Q_1-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL2_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL2S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL2S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL2S_2D_Q_A;
             break;

    case 211: // enriched Q_1-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL2_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL2SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL2SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL2SE_2D_Q_A;
             break;

    case 221: // enriched P_2-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL2_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M2_2D_Q_M;
             ElementForShape[Parallelogram] = C_M2_2D_Q_A;
             ElementForShape[Rectangle] = C_M2_2D_Q_A;
             break;

    case 301: // enriched Q_2-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL3_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL3S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL3S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL3S_2D_Q_A;
             break;

    case 312: // enriched Q_2-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL3_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL3SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL3SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL3SE_2D_Q_A;
             break;

    case 322: // enriched P_3-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL3_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M3_2D_Q_M;
             ElementForShape[Parallelogram] = C_M3_2D_Q_A;
             ElementForShape[Rectangle] = C_M3_2D_Q_A;
             break;

    case 402: // enriched Q_3-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL4_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL4S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL4S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL4S_2D_Q_A;
             break;

    case 413: // enriched Q_3-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL4_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL4SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL4SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL4SE_2D_Q_A;
             break;

    case 423: // enriched P_4-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL4_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M4_2D_Q_M;
             ElementForShape[Parallelogram] = C_M4_2D_Q_A;
             ElementForShape[Rectangle] = C_M4_2D_Q_A;
             break;
    
    case 503: // enriched Q_4-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL5_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL5S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL5S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL5S_2D_Q_A;
             break;

    case 514: // enriched Q_4-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL5_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL5SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL5SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL5SE_2D_Q_A;
             break;
    
    case 524: // enriched P_5-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_UL5_2D_T_A;
             OutPut("Using usual local projection element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M5_2D_Q_M;
             ElementForShape[Parallelogram] = C_M5_2D_Q_A;
             ElementForShape[Rectangle] = C_M5_2D_Q_A;
             break;

    case 604: // enriched Q_5-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P6_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL6S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL6S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL6S_2D_Q_A;
             break;

    case 615: // enriched Q_5-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P6_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL6SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL6SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL6SE_2D_Q_A;
             break;
    
    case 625: // enriched P_6-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P6_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M6_2D_Q_M;
             ElementForShape[Parallelogram] = C_M6_2D_Q_A;
             ElementForShape[Rectangle] = C_M6_2D_Q_A;
             break;

    case 705: // enriched Q_6-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P7_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL7S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL7S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL7S_2D_Q_A;
             break;

    case 716: // enriched Q_6-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P7_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL7SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL7SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL7SE_2D_Q_A;
             break;
    
    case 726: // enriched P_8-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P8_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M8_2D_Q_M;
             ElementForShape[Parallelogram] = C_M8_2D_Q_A;
             ElementForShape[Rectangle] = C_M8_2D_Q_A;
             break;

    case 806: // enriched Q_7-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P8_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL8S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL8S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL8S_2D_Q_A;
             break;

    case 817: // enriched Q_7-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P8_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL8SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL8SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL8SE_2D_Q_A;
             break;
    
    case 827: // enriched P_8-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P8_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M8_2D_Q_M;
             ElementForShape[Parallelogram] = C_M8_2D_Q_A;
             ElementForShape[Rectangle] = C_M8_2D_Q_A;
             break;

    case 907: // enriched Q_8-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P9_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL9S_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL9S_2D_Q_A;
             ElementForShape[Rectangle] = C_UL9S_2D_Q_A;
             break;

    case 918: // enriched Q_8-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P9_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_UL9SE_2D_Q_M;
             ElementForShape[Parallelogram] = C_UL9SE_2D_Q_A;
             ElementForShape[Rectangle] = C_UL9SE_2D_Q_A;
             break;
    
    case 928: // enriched P_9-element, used on Shishkin meshes
             ElementForShape[Triangle] = C_P9_2D_T_A;
             OutPut("Using usual element on triangles" << endl);
             ElementForShape[Quadrangle] = C_M9_2D_Q_M;
             ElementForShape[Parallelogram] = C_M9_2D_Q_A;
             ElementForShape[Rectangle] = C_M9_2D_Q_A;
             break;

	     // only function zero	     
    case 102:  ElementForShape[Triangle] = C_P00_2D_T_A;
             ElementForShape[Quadrangle] = C_Q00_2D_Q_M;
             ElementForShape[Parallelogram] = C_Q00_2D_Q_A;
             ElementForShape[Rectangle] = C_Q00_2D_Q_A;
             break;

             // discontinous elements
    case -11:
             ElementForShape[Triangle] = D_P1_2D_T_A;
             ElementForShape[Quadrangle] = D_Q1_2D_Q_M;
             ElementForShape[Parallelogram] = D_Q1_2D_Q_A;
             ElementForShape[Rectangle] = D_Q1_2D_Q_A;
             break;

    case -12:
             ElementForShape[Triangle] = D_P2_2D_T_A;
             ElementForShape[Quadrangle] = D_Q2_2D_Q_M;
             ElementForShape[Parallelogram] = D_Q2_2D_Q_A;
             ElementForShape[Rectangle] = D_Q2_2D_Q_A;
             break;

    case -13:
             ElementForShape[Triangle] = D_P3_2D_T_A;
             ElementForShape[Quadrangle] = D_Q3_2D_Q_M;
             ElementForShape[Parallelogram] = D_Q3_2D_Q_A;
             ElementForShape[Rectangle] = D_Q3_2D_Q_A;
             break;

    case -14:
             ElementForShape[Triangle] = D_P4_2D_T_A;
             ElementForShape[Quadrangle] = D_Q4_2D_Q_M;
             ElementForShape[Parallelogram] = D_Q4_2D_Q_A;
             ElementForShape[Rectangle] = D_Q4_2D_Q_A;
             break;

              // discontionuous P-elements on quadrangles
    case -110:
             ElementForShape[Triangle] = D_P1_2D_T_A;
             ElementForShape[Quadrangle] = D_P1_2D_Q_M;
             ElementForShape[Parallelogram] = D_P1_2D_Q_A;
             ElementForShape[Rectangle] = D_P1_2D_Q_A;
             break;

    case -120:
             ElementForShape[Triangle] = D_P2_2D_T_A;
             ElementForShape[Quadrangle] = D_P2_2D_Q_M;
             ElementForShape[Parallelogram] = D_P2_2D_Q_A;
             ElementForShape[Rectangle] = D_P2_2D_Q_A;
             break;

    case -130:
             ElementForShape[Triangle] = D_P3_2D_T_A;
             ElementForShape[Quadrangle] = D_P3_2D_Q_M;
             ElementForShape[Parallelogram] = D_P3_2D_Q_A;
             ElementForShape[Rectangle] = D_P3_2D_Q_A;
             break;

    case -140:
             ElementForShape[Triangle] = D_P4_2D_T_A;
             ElementForShape[Quadrangle] = D_P4_2D_Q_M;
             ElementForShape[Parallelogram] = D_P4_2D_Q_A;
             ElementForShape[Rectangle] = D_P4_2D_Q_A;
             break;

    //====================================
    //========Vector basis Raviart-Thomas  element=============
    case 1000:
      ElementForShape[Triangle] = N_RT0_2D_T_A;
      ElementForShape[Quadrangle] = N_RT0_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT0_2D_Q_A;
      ElementForShape[Rectangle] = N_RT0_2D_Q_A;
      break;
    case 1001:
      ElementForShape[Triangle] = N_RT1_2D_T_A;
      ElementForShape[Quadrangle] = N_RT1_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT1_2D_Q_A;
      ElementForShape[Rectangle] = N_RT1_2D_Q_A;
      break;
    case 1002:
      ElementForShape[Triangle] = N_RT2_2D_T_A;
      ElementForShape[Quadrangle] = N_RT2_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT2_2D_Q_A;
      ElementForShape[Rectangle] = N_RT2_2D_Q_A;
      break;
    case 1003:
      ElementForShape[Triangle] = N_RT3_2D_T_A;
      ElementForShape[Quadrangle] = N_RT3_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT3_2D_Q_A;
      ElementForShape[Rectangle] = N_RT3_2D_Q_A;
      break;
      //========Vector basis BDM  element=============
    case 1011:
      ElementForShape[Triangle] = N_BDM1_2D_T_A;
      ElementForShape[Quadrangle] = N_BDM1_2D_Q_M;
      ElementForShape[Parallelogram] = N_BDM1_2D_Q_A;
      ElementForShape[Rectangle] = N_BDM1_2D_Q_A;
      break;
    case 1012:
      ElementForShape[Triangle] = N_BDM2_2D_T_A;
      ElementForShape[Quadrangle] = N_BDM2_2D_Q_M;
      ElementForShape[Parallelogram] = N_BDM2_2D_Q_A;
      ElementForShape[Rectangle] = N_BDM2_2D_Q_A;
      break;
    case 1013:
      ElementForShape[Triangle] = N_BDM3_2D_T_A;
      ElementForShape[Quadrangle] = N_BDM3_2D_Q_M;
      ElementForShape[Parallelogram] = N_BDM3_2D_Q_A;
      ElementForShape[Rectangle] = N_BDM3_2D_Q_A;
      break;

    //========LOCALPROJECTION WITH EXP BUBBLE=============
    // Q1+bubble*P0
    case 1100: 
             ElementForShape[Triangle] = C_UL1_2D_T_A;
             OutPut("Using usual LPS bubble element on triangles" << endl);     
             ElementForShape[Quadrangle] = C_EL1_2D_Q_M;
             ElementForShape[Parallelogram] = C_EL1_2D_Q_A;
             ElementForShape[Rectangle] = C_EL1_2D_Q_A;
             break;     
  
    default: cerr << "unknown order" << endl;
             exit(-1);
             break;
  } // endswitch

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace(BoundaryCondition);
}

/** constructor for building a space with elements of order k */
TFESpace2D::TFESpace2D(TCollection *coll, char *name, char *description, 
                       BoundCondFunct2D *BoundaryCondition, SpaceType type,
                       int ord, TCollection *mortarcoll) :
     TFESpace(coll, name, description)
{
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  AllElements = NULL;
  DGSpace = 0;
#ifdef __MORTAR__
  MortarColl = mortarcoll;
#endif

  ElementForShape = new FE2D[N_SHAPES];

  // build ElementForShape array
  switch(type)
  {
    // find velo space for discontinuous pressure
    case DiscP_USpace:
      switch(ord)
      {
        case 1:
          Error("This makes no sense." << endl);
          ElementForShape[Triangle] = C_P1_2D_T_A;
          ElementForShape[Quadrangle] = C_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
          ElementForShape[Rectangle] = C_Q1_2D_Q_A;
        break;

        case 2:
          ElementForShape[Triangle] = C_B2_2D_T_A;
          ElementForShape[Quadrangle] = C_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
          ElementForShape[Rectangle] = C_Q2_2D_Q_A;
        break;

        case 3:
          ElementForShape[Triangle] = C_B3_2D_T_A;
          ElementForShape[Quadrangle] = C_Q3_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q3_2D_Q_A;
          ElementForShape[Rectangle] = C_Q3_2D_Q_A;
        break;

        case 4:
          ElementForShape[Triangle] = C_B4_2D_T_A;
          ElementForShape[Quadrangle] = C_Q4_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q4_2D_Q_A;
          ElementForShape[Rectangle] = C_Q4_2D_Q_A;
        break;

        // Scott-Vogelius
        case 222:
          ElementForShape[Triangle] = C_SV2_2D_T_A;
          ElementForShape[Quadrangle] = C_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
          ElementForShape[Rectangle] = C_Q2_2D_Q_A;
          break;

        default:
          Error("Space is not available" << endl);
          exit(-1);
      }
    break;
    
    // find pressure space for discontinuous pressure
    case DiscP_PSpace:
      switch(ord)
      {
        case 0:
          ElementForShape[Triangle] = C_P0_2D_T_A;
          ElementForShape[Quadrangle] = C_Q0_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q0_2D_Q_A;
          ElementForShape[Rectangle] = C_Q0_2D_Q_A;
        break;

        case 1:
          ElementForShape[Triangle] = D_P1_2D_T_A;
          ElementForShape[Quadrangle] = D_P1_2D_Q_M;
          ElementForShape[Parallelogram] = D_P1_2D_Q_A;
          ElementForShape[Rectangle] = D_P1_2D_Q_A;
        break;

        case 2:
          ElementForShape[Triangle] = D_P2_2D_T_A;
          ElementForShape[Quadrangle] = D_P2_2D_Q_M;
          ElementForShape[Parallelogram] = D_P2_2D_Q_A;
          ElementForShape[Rectangle] = D_P2_2D_Q_A;
        break;

        case 3:
          ElementForShape[Triangle] = D_P3_2D_T_A;
          ElementForShape[Quadrangle] = D_P3_2D_Q_M;
          ElementForShape[Parallelogram] = D_P3_2D_Q_A;
          ElementForShape[Rectangle] = D_P3_2D_Q_A;
        break;

        case 4:
          ElementForShape[Triangle] = D_P4_2D_T_A;
          ElementForShape[Quadrangle] = D_P4_2D_Q_M;
          ElementForShape[Parallelogram] = D_P4_2D_Q_A;
          ElementForShape[Rectangle] = D_P4_2D_Q_A;
        break;

        case 5:
          ElementForShape[Triangle] = D_P4_2D_T_A;
          OutPut("Using P4 on triangles" << endl);
          ElementForShape[Quadrangle] = D_P5_2D_Q_M;
          ElementForShape[Parallelogram] = D_P5_2D_Q_A;
          ElementForShape[Rectangle] = D_P5_2D_Q_A;
        break;

        case 6:
          ElementForShape[Triangle] = D_P4_2D_T_A;
          OutPut("Using P4 on triangles" << endl);
          ElementForShape[Quadrangle] = D_P6_2D_Q_M;
          ElementForShape[Parallelogram] = D_P6_2D_Q_A;
          ElementForShape[Rectangle] = D_P6_2D_Q_A;
        break;

        case 7:
          ElementForShape[Triangle] = D_P4_2D_T_A;
          OutPut("Using P4 on triangles" << endl);
          ElementForShape[Quadrangle] = D_P7_2D_Q_M;
          ElementForShape[Parallelogram] = D_P7_2D_Q_A;
          ElementForShape[Rectangle] = D_P7_2D_Q_A;
        break;

        case 92:
          ElementForShape[Triangle] = D_P2_2D_T_A;
          ElementForShape[Quadrangle] = D_D2_2D_Q_M;
          ElementForShape[Parallelogram] = D_D2_2D_Q_A;
          ElementForShape[Rectangle] = D_D2_2D_Q_A;
        break;
        
        // only function zero
        case 102:
          ElementForShape[Triangle] = C_P00_2D_T_A;
          ElementForShape[Quadrangle] = C_Q00_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q00_2D_Q_A;
          ElementForShape[Rectangle] = C_Q00_2D_Q_A;
          break;
        // Scott-Vogelius
        case 222:
          ElementForShape[Triangle] = D_SV1_2D_T_A;
          ElementForShape[Quadrangle] = D_P1_2D_Q_M;
          ElementForShape[Parallelogram] = D_P1_2D_Q_A;
          ElementForShape[Rectangle] = D_P1_2D_Q_A;
          break;      
	  default:
          Error("Space is not available" << endl);
          exit(-1);
      }
    break;
    
    // find pressure and velo space for continuous pressure
    case ContP_USpace:
    case ContP_PSpace:
      switch(ord)
      {
        case 1:
          ElementForShape[Triangle] = C_P1_2D_T_A;
          ElementForShape[Quadrangle] = C_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
          ElementForShape[Rectangle] = C_Q1_2D_Q_A;
        break;

        case 2:
          ElementForShape[Triangle] = C_P2_2D_T_A;
          ElementForShape[Quadrangle] = C_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
          ElementForShape[Rectangle] = C_Q2_2D_Q_A;
        break;

        case 3:
          ElementForShape[Triangle] = C_P3_2D_T_A;
          ElementForShape[Quadrangle] = C_Q3_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q3_2D_Q_A;
          ElementForShape[Rectangle] = C_Q3_2D_Q_A;
        break;

        case 4:
          ElementForShape[Triangle] = C_P4_2D_T_A;
          ElementForShape[Quadrangle] = C_Q4_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q4_2D_Q_A;
          ElementForShape[Rectangle] = C_Q4_2D_Q_A;
        break;

        case 5:
          ElementForShape[Triangle] = C_P5_2D_T_A;
          ElementForShape[Quadrangle] = C_Q5_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q5_2D_Q_A;
          ElementForShape[Rectangle] = C_Q5_2D_Q_A;
        break;

        case 6:
          ElementForShape[Triangle] = C_P6_2D_T_A;
          ElementForShape[Quadrangle] = C_Q6_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q6_2D_Q_A;
          ElementForShape[Rectangle] = C_Q6_2D_Q_A;
        break;

        case 7:
          ElementForShape[Triangle] = C_P7_2D_T_A;
          ElementForShape[Quadrangle] = C_Q7_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q7_2D_Q_A;
          ElementForShape[Rectangle] = C_Q7_2D_Q_A;
        break;
        
        case 8:
          ElementForShape[Triangle] = C_P8_2D_T_A;
          ElementForShape[Quadrangle] = C_Q8_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q8_2D_Q_A;
          ElementForShape[Rectangle] = C_Q8_2D_Q_A;
        break;

        case 9:
          ElementForShape[Triangle] = C_P9_2D_T_A;
          ElementForShape[Quadrangle] = C_Q9_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q9_2D_Q_A;
          ElementForShape[Rectangle] = C_Q9_2D_Q_A;
        break;

        // P1mini
        case 101: 
          ElementForShape[Triangle] = C_P1MINI_2D_T_A;
          ElementForShape[Quadrangle] = C_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
          ElementForShape[Rectangle] = C_Q1_2D_Q_A;
          cout << "P1MINI works only on triangles" << endl;
        break;

        case -2:
          ElementForShape[Triangle] = C_P9_2D_T_A;
          ElementForShape[Quadrangle] = B_IB2_2D_Q_M;
          ElementForShape[Parallelogram] = B_IB2_2D_Q_A;
          ElementForShape[Rectangle] = B_IB2_2D_Q_A;
          OutPut("Space for triangles not yet implemented !!!!!!"<< endl)
        break;

        //========LOCALPROJECTION=============
        // Q1+bubble*P0
        case 100: 
                 ElementForShape[Triangle] = C_UL1_2D_T_A;
                 ElementForShape[Quadrangle] = C_UL1_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL1_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL1_2D_Q_A;
                 break;
        // Q2+bubble*P1
        case 201:  
                 ElementForShape[Triangle] = C_UL2_2D_T_A;
                 ElementForShape[Quadrangle] = C_UL2_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL2_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL2_2D_Q_A;
                 break;
        // Q3+bubble*P2
        case 302: 
                 ElementForShape[Triangle] = C_UL3_2D_T_A;
                 ElementForShape[Quadrangle] = C_UL3_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL3_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL3_2D_Q_A;
                 break;
        // Q4+bubble*P3
        case 403:  
                 ElementForShape[Triangle] = C_UL4_2D_T_A;
                 ElementForShape[Quadrangle] = C_UL4_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL4_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL4_2D_Q_A;
                 break;
        // Q5+bubble*P4
        case 504:  
                 ElementForShape[Triangle] = C_UL5_2D_T_A;
                 ElementForShape[Quadrangle] = C_UL5_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL5_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL5_2D_Q_A;
                 break; 

        case 200: // enriched Q_1-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL2_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL2S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL2S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL2S_2D_Q_A;
                 break;
    
        case 211: // enriched Q_1-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL2_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL2SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL2SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL2SE_2D_Q_A;
                 break;
    
        case 221: // enriched P_2-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL2_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M2_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M2_2D_Q_A;
                 ElementForShape[Rectangle] = C_M2_2D_Q_A;
                 break;

        case 301: // enriched Q_2-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL3_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL3S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL3S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL3S_2D_Q_A;
                 break;
    
        case 312: // enriched Q_2-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL3_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL3SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL3SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL3SE_2D_Q_A;
                 break;
    
        case 322: // enriched P_3-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL3_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M3_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M3_2D_Q_A;
                 ElementForShape[Rectangle] = C_M3_2D_Q_A;
                 break;

        case 402: // enriched Q_3-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL4_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL4S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL4S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL4S_2D_Q_A;
                 break;
    
        case 413: // enriched Q_3-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL4_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL4SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL4SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL4SE_2D_Q_A;
                 break;
    
        case 423: // enriched P_4-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL4_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M4_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M4_2D_Q_A;
                 ElementForShape[Rectangle] = C_M4_2D_Q_A;
                 break;
    
        case 503: // enriched Q_4-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL5_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL5S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL5S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL5S_2D_Q_A;
                 break;
    
        case 514: // enriched Q_4-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL5_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL5SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL5SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL5SE_2D_Q_A;
                 break;
    
        case 524: // enriched P_5-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_UL5_2D_T_A;
                 OutPut("Using usual local projection element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M5_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M5_2D_Q_A;
                 ElementForShape[Rectangle] = C_M5_2D_Q_A;
                 break;

        case 604: // enriched Q_5-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P6_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL6S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL6S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL6S_2D_Q_A;
                 break;
    
        case 615: // enriched Q_5-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P6_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL6SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL6SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL6SE_2D_Q_A;
                 break;
        
        case 625: // enriched P_6-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P6_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M6_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M6_2D_Q_A;
                 ElementForShape[Rectangle] = C_M6_2D_Q_A;
                 break;
    
        case 705: // enriched Q_6-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P7_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL7S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL7S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL7S_2D_Q_A;
                 break;
    
        case 716: // enriched Q_6-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P7_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL7SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL7SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL7SE_2D_Q_A;
                 break;
        
        case 726: // enriched P_8-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P8_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M8_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M8_2D_Q_A;
                 ElementForShape[Rectangle] = C_M8_2D_Q_A;
                 break;
    
        case 806: // enriched Q_7-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P8_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL8S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL8S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL8S_2D_Q_A;
                 break;
    
        case 817: // enriched Q_7-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P8_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL8SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL8SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL8SE_2D_Q_A;
                 break;
        
        case 827: // enriched P_8-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P8_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M8_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M8_2D_Q_A;
                 ElementForShape[Rectangle] = C_M8_2D_Q_A;
                 break;
    
        case 907: // enriched Q_8-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P9_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL9S_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL9S_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL9S_2D_Q_A;
                 break;
    
        case 918: // enriched Q_8-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P9_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_UL9SE_2D_Q_M;
                 ElementForShape[Parallelogram] = C_UL9SE_2D_Q_A;
                 ElementForShape[Rectangle] = C_UL9SE_2D_Q_A;
                 break;
        
        case 928: // enriched P_9-element, used on Shishkin meshes
                 ElementForShape[Triangle] = C_P9_2D_T_A;
                 OutPut("Using usual element on triangles" << endl);
                 ElementForShape[Quadrangle] = C_M9_2D_Q_M;
                 ElementForShape[Parallelogram] = C_M9_2D_Q_A;
                 ElementForShape[Rectangle] = C_M9_2D_Q_A;
                 break;
    
        //====================================
        //====================================
    
        default:
          Error("Space is not available." << endl);
          exit(-1);
      } // endswitch
    break;
    
    // find velo space for nonconforming fe
    case Non_USpace:
      switch(ord)
      {
        case 1:
          ElementForShape[Triangle] = N_P1_2D_T_A;
          ElementForShape[Quadrangle] = N_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
          ElementForShape[Rectangle] = N_Q1_2D_Q_A;
        break;

        case 2:
          ElementForShape[Triangle] = N_P2_2D_T_A;
          ElementForShape[Quadrangle] = N_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q2_2D_Q_A;
          ElementForShape[Rectangle] = N_Q2_2D_Q_A;
        break;

        case 3:
          ElementForShape[Triangle] = N_P3_2D_T_A;
          ElementForShape[Quadrangle] = N_Q3_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q3_2D_Q_A;
          ElementForShape[Rectangle] = N_Q3_2D_Q_A;
        break;

        case 4:
          ElementForShape[Triangle] = N_P4_2D_T_A;
          ElementForShape[Quadrangle] = N_Q4_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q4_2D_Q_A;
          ElementForShape[Rectangle] = N_Q4_2D_Q_A;
        break;

        case 5:
          ElementForShape[Triangle] = N_P5_2D_T_A;
          ElementForShape[Quadrangle] = N_Q5_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q5_2D_Q_A;
          ElementForShape[Rectangle] = N_Q5_2D_Q_A;
        break;

        case 101: 
          ElementForShape[Triangle] = N_P1MOD_2D_T_A;
          ElementForShape[Quadrangle] = N_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
          ElementForShape[Rectangle] = N_Q1_2D_Q_A;
          cout << "P1MOD works only on triangles" << endl;
        break;

        //========Vector basis Raviart-Thomas  element=============
        case 1000:
          ElementForShape[Triangle] = N_RT0_2D_T_A;
          ElementForShape[Quadrangle] = N_RT0_2D_Q_M;
          ElementForShape[Parallelogram] = N_RT0_2D_Q_A;
          ElementForShape[Rectangle] = N_RT0_2D_Q_A;
          break;
          
        case 1001:
          ElementForShape[Triangle] = N_RT1_2D_T_A;
          ElementForShape[Quadrangle] = N_RT1_2D_Q_M;
          ElementForShape[Parallelogram] = N_RT1_2D_Q_A;
          ElementForShape[Rectangle] = N_RT1_2D_Q_A;
          break;
        case 1002:
          ElementForShape[Triangle] = N_RT2_2D_T_A;
          ElementForShape[Quadrangle] = N_RT2_2D_Q_M;
          ElementForShape[Parallelogram] = N_RT2_2D_Q_A;
          ElementForShape[Rectangle] = N_RT2_2D_Q_A;
          break;
        case 1003:
          ElementForShape[Triangle] = N_RT3_2D_T_A;
          ElementForShape[Quadrangle] = N_RT3_2D_Q_M;
          ElementForShape[Parallelogram] = N_RT3_2D_Q_A;
          ElementForShape[Rectangle] = N_RT3_2D_Q_A;
          break;
          //========Vector basis BDM  element=============
        case 1011:
          ElementForShape[Triangle] = N_BDM1_2D_T_A;
          ElementForShape[Quadrangle] = N_BDM1_2D_Q_M;
          ElementForShape[Parallelogram] = N_BDM1_2D_Q_A;
          ElementForShape[Rectangle] = N_BDM1_2D_Q_A;
          break;
        case 1012:
          ElementForShape[Triangle] = N_BDM2_2D_T_A;
          ElementForShape[Quadrangle] = N_BDM2_2D_Q_M;
          ElementForShape[Parallelogram] = N_BDM2_2D_Q_A;
          ElementForShape[Rectangle] = N_BDM2_2D_Q_A;
          break;
        case 1013:
          ElementForShape[Triangle] = N_BDM3_2D_T_A;
          ElementForShape[Quadrangle] = N_BDM3_2D_Q_M;
          ElementForShape[Parallelogram] = N_BDM3_2D_Q_A;
          ElementForShape[Rectangle] = N_BDM3_2D_Q_A;
          break;
    
        default:
          Error("This nonconforming space (order " << ord << ")");
          Error("is not available." << endl);
          exit(-1);
      }
    break;

    default:
      Error("Wrong space type" << endl);
      exit(-1);
  }

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace(BoundaryCondition);
}

/** constructor for building a space with the given elements */
TFESpace2D::TFESpace2D(TCollection *coll, char *name, char *description,
               BoundCondFunct2D *BoundaryCondition,
               FE2D *fes, TCollection *mortarcoll) :
    TFESpace(coll, name, description)
{
  int i, N_;
  
  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;
  UsedElements = NULL;
  ElementForShape = NULL;
  DGSpace = 0;
  //AllElements = fes;
  N_ = coll->GetN_Cells();
  AllElements = new FE2D[N_];
  
  for(i=0;i<N_;i++)
    AllElements[i] =  fes[i];  
  
#ifdef __MORTAR__
  MortarColl = mortarcoll;
#endif

  // find out all used elements
  FindUsedElements();

  // construct space
  ConstructSpace(BoundaryCondition);
}

/** return the FE Id for element i, corresponding to cell */
FE2D TFESpace2D::GetFE2D(int i, TBaseCell *cell)
{
  FE2D ret;

  if(AllElements)
    ret=AllElements[i];
  else
    ret=ElementForShape[cell->GetType()];

  return ret;
}

void TFESpace2D::FindUsedElements()
{
  TBaseCell *cell;
  int i, j, N_;
  int Used[N_FEs2D];

  memset(Used,0,N_FEs2D*SizeOfInt);

  N_ = N_Cells;
  for(i=0;i<N_;i++)
  {
    cell = Collection->GetCell(i);
    Used[GetFE2D(i, cell)] = 1;
  }

  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
      N_UsedElements++;

  UsedElements = new FE2D[N_UsedElements];
  j=0;
  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
    {
      UsedElements[j]=(FE2D)i;
      j++;
    }

  // OutPut(endl << "N_UsedElements: " << N_UsedElements << endl);
  // for(i=0;i<N_UsedElements;i++)
  //   OutPut("UsedElement[" << i << "]: " << UsedElements[i] << endl);
}

void TFESpace2D::ConstructSpace(BoundCondFunct2D *BoundaryCondition)
{
  int i, j, k, l, m, m2, n, comp, N_Edges, NEdges;
  int *v;
  TBaseCell *cell, *neigh, *child1, *child2;
  TJoint *joint;
  TBoundComp2D *BoundComp;
  TBoundEdge *BoundEdge;
  double t0,t1;
  BoundCond Cond0, Cond1;

  TFE2DMapper *mapper;
  TFE2DMapper1Reg *mapper1reg;

  const int *TmpoEnE, *TmpLen1, *TmpEC, *TmpLen2, *TmpoEnlE;
  int MaxLen1, MaxLen2;
  TRefDesc *refdesc;

  int SumLocDOF;
  int count, *BoundaryUpperBound, DirichletUpperBound;
  int DirichletCounter, *BoundCounter, Counter;
  int SlaveMark, *BoundMark, InnerMark, DirichletMark;
  int DirichletOffset, InnerOffset, SlaveOffset;
  int *BoundOffset;
  int N_Slave;

  FE2D FEType0, FEType1, FEType2;
  TFE2D *FE0, *FE1, *FE2;
  TFEDesc2D *FEDesc0_Obj, *FEDesc1_Obj, *FEDesc2_Obj;
  FEDesc2D FEDesc0, FEDesc1, FEDesc2;

  int I_K0, I_K1, I_K2;
  int *J_K0, *J_K1, *J_K2;
  int *Indices0, *Indices1, *Indices2;
  int c1, c2, e1, e2, chnum1, chnum2;
  double eps=1e-6;

  THangingNode *hn;
  TVector<THangingNode *> *VHN = new TVector<THangingNode *>(20,20);
  TVector<int> *HNNumbers = new TVector<int>(20,20);

  N_ActiveDegrees = 0;
  N_SlaveDegrees = 0;

  N_DiffBoundNodeTypes = N_BOUNDCOND-1; // not only Neumann nodes possible
  BoundaryUpperBound = new int[N_DiffBoundNodeTypes];
  BoundaryNodeTypes = new BoundCond[N_DiffBoundNodeTypes];
  BoundCounter = new int[N_DiffBoundNodeTypes];
  BoundMark = new int[N_DiffBoundNodeTypes];
  N_BoundaryNodes = new int[N_DiffBoundNodeTypes];
  BoundOffset = new int[N_DiffBoundNodeTypes];
  BoundaryNodesBound = new int[N_DiffBoundNodeTypes];

  BoundaryNodeTypes[0] = NEUMANN;
  BoundaryNodeTypes[1] = ROBIN;
  BoundaryNodeTypes[2] = SLIP;
  BoundaryNodeTypes[3] = FREESURF;
  BoundaryNodeTypes[4] = SLIP_FRICTION_PENETRATION_RESISTANCE;
  BoundaryNodeTypes[5] = INTERFACE;
  BoundaryNodeTypes[6] = SUBDOMAIN_INTERFACE;
  BoundaryNodeTypes[7] = SUBDOMAIN_HALOBOUND;
  N_BoundaryNodes[0] = 0;
  N_BoundaryNodes[1] = 0;
  N_BoundaryNodes[2] = 0;
  N_BoundaryNodes[3] = 0;
  N_BoundaryNodes[4] = 0;
  N_BoundaryNodes[5] = 0;
  N_BoundaryNodes[6] = 0;
  N_BoundaryNodes[7] = 0;

  // reset clipboards to -1
  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    k=cell->GetN_Edges();
    for(j=0;j<k;j++)
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);
      if(neigh)
      { 
        m = neigh->GetN_Children();
        for(l=0;l<m;l++)
        {
          child1 = neigh->GetChild(l);
          if(child1) child1->SetClipBoard(-1);
        }
        neigh->SetClipBoard(-1);
      }
    }
    cell->SetClipBoard(-1);
  } // endfor i

  for(i=0;i<N_Cells;i++)
    Collection->GetCell(i)->SetClipBoard(i);

  // set number i into clipboard, count the number of local degrees of
  // freedom
  count=0;
  DirichletUpperBound=0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
    BoundaryUpperBound[i]=0;

  BeginIndex=new int[N_Cells+1];
  BeginIndex[0]=0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);

    FEType0 = GetFE2D(i, cell);
    FE0 = TFEDatabase2D::GetFE2D(FEType0);
    count += FE0->GetSize();

    BeginIndex[i+1]=count;

    // find number of boundary nodes in this element
    FEDesc0_Obj = FE0->GetFEDesc2D();
    l = FEDesc0_Obj->GetN_JointDOF();

    k = cell->GetN_Edges();
    for(j=0;j<k;j++) // loop over all edges of cell
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == BoundaryEdge ||
         joint->GetType() == IsoBoundEdge)
      {
        BoundEdge = (TBoundEdge *)joint;
        BoundComp = BoundEdge->GetBoundComp();
        BoundEdge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        if (t0 < t1)
        {
          BoundaryCondition(comp, t0+eps, Cond0);
          BoundaryCondition(comp, t1-eps, Cond1);
        }
        else
        {
          BoundaryCondition(comp, t0-eps, Cond0);
          BoundaryCondition(comp, t1+eps, Cond1);
        }

        if(Cond0==Cond1)
        {
          switch(Cond0)
          {
            case DIRICHLET:
              // boundary nodes are Dirichlet nodes
              DirichletUpperBound += l;
            break;

            default:
              // non Dirichlet boundary node
              BoundaryUpperBound[Cond0-1] += l;
          } // endswitch l
        } // endif
        else
        {
          OutPut("different boundary condition on one ");
          OutPut("edge are not allowed" << endl);
          exit(4711);
        }
      } // endif
      else
      {
        neigh = joint->GetNeighbour(cell);
        if(neigh)
        {
          if(neigh->GetClipBoard() == -1 &&
             neigh->GetRefDesc()->GetType() == NoRef)
          {
            if(joint->GetType() == InterfaceJoint ||
               joint->GetType() == IsoInterfaceJoint)
            {
              // the neighbour is not a member of current collection
              BoundEdge = (TBoundEdge *)joint;
              BoundComp = BoundEdge->GetBoundComp();
              BoundEdge->GetParameters(t0, t1);
              comp=BoundComp->GetID();
              t0 = 0.5*(t0+t1);
              BoundaryCondition(comp, t0, Cond0);
              switch(Cond0)
              {
                case DIRICHLET: 
                  // boundary nodes are Dirichlet nodes
                  DirichletUpperBound += l;
                break;

                default:
                  // non Dirichlet boundary node
                  BoundaryUpperBound[Cond0-1] += l;
              } // endswitch Cond0
            } // endif InterfaceJoint
#ifdef _MPI
           else 
// 	     if(joint->GetType() == SubDomainHaloJoint || joint->GetType() == JointEqN)
            { 
             // the neighbour is not a member of current collection
             Cond0 = SUBDOMAIN_HALOBOUND;
             BoundaryUpperBound[Cond0 - 1] += l; //  SUBDOMAIN_INTERFACE == 8, see constants.h
            }
#else
            else
              DirichletUpperBound += l;
#endif	    
          } // end == -1
        } // if(neigh)
      } // end no boundary
    } // endfor j

  } // endfor i

/*
  OutPut("DirichletUpperBound: " << DirichletUpperBound << endl);
  OutPut("N_DiffBoundNodeTypes: " << N_DiffBoundNodeTypes << endl);
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    OutPut("type[" << i << "]: " << BoundaryNodeTypes[i]);
    OutPut(" number of: " << BoundaryUpperBound[i] << endl);
  }
*/
// #ifdef _MPI
//    MPI_Finalize();
// exit(0);
// #endif
  GlobalNumbers = new int[count];
  memset(GlobalNumbers, -1, SizeOfInt*count);
  SumLocDOF = count;

  // start DOF manager
  DirichletCounter=FIRSTMARK+1;
  l=-DirichletUpperBound + FIRSTMARK +1;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundCounter[i] = l;
    l -=  BoundaryUpperBound[i];
  }
  Counter = l;

  // OutPut("DirichletCounter: " << DirichletCounter << endl);
  // for(i=0;i<N_DiffBoundNodeTypes;i++)
  //   OutPut(i << "   " << BoundCounter[i] << endl);
  // OutPut("Counter: " << Counter << endl);

//   OutPut("Number of cells: " << N_Cells << endl);
  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    N_Edges=cell->GetN_Edges();

    FEType0 = GetFE2D(i, cell);
    FE0 = TFEDatabase2D::GetFE2D(FEType0);
    FE0->GetFEDesc2D(FEDesc0, FEDesc0_Obj);
    // FEDesc0 = FE0->GetFEDesc2D_ID();
    // FEDesc0_Obj = FE0->GetFEDesc2D();
    // FEDesc0_Obj = TFEDatabase2D::GetFEDesc2D(FEDesc0);
    I_K0 = BeginIndex[i];
    J_K0 = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      Indices0 = FEDesc0_Obj->GetJointDOF(j);
      if(joint->GetType() == BoundaryEdge ||
         joint->GetType() == IsoBoundEdge)
      {
        // boundary joint
        
        BoundEdge=(TBoundEdge *)joint;
        BoundComp=BoundEdge->GetBoundComp();
        BoundEdge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        if (t0 < t1)
        {
          BoundaryCondition(comp, t0+eps, Cond0);
          BoundaryCondition(comp, t1-eps, Cond1);
        }
        else
        {
          BoundaryCondition(comp, t0-eps, Cond0);
          BoundaryCondition(comp, t1+eps, Cond1);
        }

        if(Cond0==Cond1)
        {
          switch(Cond0)
          {
            case DIRICHLET: 
              // boundary nodes are Dirichlet nodes
              // OutPut("Dirichlet Boundary Edge" << endl);
              mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc0);
              mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                               DirichletCounter);
             break;

             default:
               // non Dirichlet boundary node
               // OutPut("non Dirichlet Boundary Edge" << endl);
               mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc0);
               mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                BoundCounter[Cond0-1]);
          } // endswitch l
        } // endif
        else
        {
          OutPut("different boundary condition on one edge are not allowed" << endl);
          exit(4711);
        }
      } // boundary joint
      else
      {
        // no boundary joint
        neigh = joint->GetNeighbour(cell);
        if (!neigh || joint->GetType() == MortarJoint ||
            joint->GetType() == MortarBaseJoint)
        {
          // there is no neighbour
          // => either mortar joint
          //    or finer cell in 1 regular grid
          //    will be handle from coarser cell
          if(joint->GetType() == MortarJoint ||
             joint->GetType() == MortarBaseJoint)
          {
            // do mortar mapping
            mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc0);
            mapper->MapBound(GlobalNumbers, I_K0, Indices0, Counter);
          }
        } // !neigh
        else
        {
          // there is a neighbour on same level
          n = neigh->GetClipBoard();
          if(n == -1)
          {
            // the neighbour is not a member of current collection
            if(joint->GetType() == InterfaceJoint ||
               joint->GetType() == IsoInterfaceJoint)
            {
              BoundEdge = (TBoundEdge *)joint;
              BoundComp = BoundEdge->GetBoundComp();
              BoundEdge->GetParameters(t0, t1);
              comp=BoundComp->GetID();
              t0 = 0.5*(t0+t1);
              BoundaryCondition(comp, t0, Cond0);
              switch(Cond0)
              {
                case DIRICHLET:
                  // boundary nodes are Dirichlet nodes
                  // OutPut("Dirichlet Boundary Edge" << endl);
                  mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc0);
                  mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                   DirichletCounter);
                break;
                default:
                  // non Dirichlet boundary node
                  // OutPut("non Dirichlet Boundary Edge" << endl);
                  mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc0);
                  mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                   BoundCounter[Cond0-1]);
              } // endswitch Cond0
             } // endif InterfaceJoint
#ifdef _MPI
           else 
// 	     if(joint->GetType() == SubDomainHaloJoint || joint->GetType() == JointEqN)
            {
             Cond0 = SUBDOMAIN_HALOBOUND;
             // OutPut("Halo Boundary Edge" << endl);
             mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc0);
             mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                              BoundCounter[Cond0-1]);
            }
#else
            else
             {
              // find local edge of neigh on which cell is
              l=0;
              // while(neigh->GetJoint(l)->GetNeighbour(neigh)!=cell) l++;
              while(neigh->GetJoint(l) != joint) l++;
  
              // check for children of neigh on its face l
              refdesc=neigh->GetRefDesc();
              if(refdesc->GetType() != NoRef)
              {
                refdesc->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);
                refdesc->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);
                refdesc->GetOldEdgeNewLocEdge(TmpoEnlE);
                NEdges=refdesc->GetShapeDesc()->GetN_Edges();
  
                if(TmpLen1[l]==0)
                {
                  // there is only ONE child on this edge
                  e1=TmpoEnE[l*MaxLen1+0];
                  chnum1=TmpEC[e1*MaxLen2];
                  child1=neigh->GetChild(chnum1);
                  c1=child1->GetClipBoard();
          
                  FEType1 = GetFE2D(c1, child1);
                  FE1 = TFEDatabase2D::GetFE2D(FEType1);
                  FE1->GetFEDesc2D(FEDesc1, FEDesc1_Obj);
                  // FEDesc1 = FE1->GetFEDesc2D_ID();
                  // FEDesc1_Obj = FE1->GetFEDesc2D();
                  // FEDesc1_Obj = TFEDatabase2D::GetFEDesc2D(FEDesc1);
                  I_K1 = BeginIndex[c1];
                  J_K1 = GlobalNumbers + BeginIndex[c1];
  
                  m=TmpoEnlE[chnum1*NEdges+l];
                  Indices1 = FEDesc1_Obj->GetJointDOF(m);
  
                  // OutPut("arg " << chnum1*NEdges+l << endl);
                  // OutPut("m= " << m << endl);
  
                  mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc1);
                  mapper->Map(GlobalNumbers, I_K0, I_K1, Indices0, Indices1,
                              j, m, FEDesc0_Obj, FEDesc1_Obj, Counter, 
                              VHN, HNNumbers);
                } // TmpLen1[l]==0
                else
                {
                  // there are at least two children on this edge
                  if(TmpLen1[l]==1)
                  {
                    // there are exactly two children
                    e1=TmpoEnE[l*MaxLen1+0];
                    chnum1=TmpEC[e1*MaxLen2];
                    child1=neigh->GetChild(chnum1);
                    c1=child1->GetClipBoard();
              
                    e2=TmpoEnE[l*MaxLen1+1];
                    chnum2=TmpEC[e2*MaxLen2];
                    child2=neigh->GetChild(chnum2);
                    c2=child2->GetClipBoard();
                    
                    if(c1 == -1 || c2 == -1)
                      break;
        
                    // OutPut("large cell: " << i << " small ");
                    // OutPut("cell numbers: " << c1 << "      " << c2 << endl);
        
                    FEType1 = GetFE2D(c1, child1);
                    FE1 = TFEDatabase2D::GetFE2D(FEType1);
                    FE1->GetFEDesc2D(FEDesc1, FEDesc1_Obj);
                    // FEDesc1 = FE1->GetFEDesc2D_ID();
                    // FEDesc1_Obj = FE1->GetFEDesc2D();
                    // FEDesc1_Obj = TFEDatabase2D::GetFEDesc2D(FEDesc1);
                    I_K1 = BeginIndex[c1];
                    J_K1 = GlobalNumbers + BeginIndex[c1];
                    m=TmpoEnlE[chnum1*NEdges+l];
                    
                    Indices1 = FEDesc1_Obj->GetJointDOF(m);
     
                    FEType2 = GetFE2D(c2, child2);
                    FE2 = TFEDatabase2D::GetFE2D(FEType2);
                    FE2->GetFEDesc2D(FEDesc2, FEDesc2_Obj);
                    // FEDesc2 = FE2->GetFEDesc2D_ID();
                    // FEDesc2_Obj = FE2->GetFEDesc2D();
                    // FEDesc2_Obj = TFEDatabase2D::GetFEDesc2D(FEDesc2);
                    I_K2 = BeginIndex[c2];
                    J_K2 = GlobalNumbers + BeginIndex[c2];
                    m2=TmpoEnlE[chnum2*NEdges+l];
    
                    Indices2 = FEDesc2_Obj->GetJointDOF(m2);
    
                    mapper1reg = TFEDatabase2D::GetFE2DMapper1Reg(FEDesc0, FEDesc1);
                    if(i<c1)
                    {
                      // // OutPut(i << "i<c1" << endl);
                      mapper1reg->MapCoarseFine(GlobalNumbers, I_K0, I_K1, I_K2,
                            Indices0, Indices1, Indices2,
                            j, m, m2, FEDesc0_Obj, FEDesc1_Obj, FEDesc2_Obj,
                            Counter, (c1<c2), VHN, HNNumbers);
                    }
                    else
                    {
                      // // OutPut(i << "i>c1" << endl);
                      mapper1reg->MapFineCoarse(GlobalNumbers, I_K0, I_K1, I_K2,
                            Indices0, Indices1, Indices2,
                            j, m, m2, FEDesc0_Obj, FEDesc1_Obj, FEDesc2_Obj,
                            Counter, (c1<c2), VHN, HNNumbers);
                    }
                  } // TmpLen1[l]==1
                  else
                  {
                    cerr << "more then two children on one edge ";
                    cerr << "are not allowed" << endl;
                  } // TmpLen1[l]!=1
                } // TmpLen1[l]!=0
              }
              else
              {
                // neighbour is not refined
                mapper=TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc0);
                mapper->MapBound(GlobalNumbers, I_K0, Indices0,
                                 DirichletCounter);
              }
            } // end no interfacejoint
#endif                
            
          } // n == -1
          else
          {
            // neighbour is member of this collection
            // => using mappers
            if (n>i)
            {
              // this joint was not handled until now

              FEType1 = GetFE2D(n, neigh);
              FE1 = TFEDatabase2D::GetFE2D(FEType1);
              FE1->GetFEDesc2D(FEDesc1, FEDesc1_Obj);
              // FEDesc1 = FE1->GetFEDesc2D_ID();
              // FEDesc1_Obj = FE1->GetFEDesc2D();
              // FEDesc1_Obj = TFEDatabase2D::GetFEDesc2D(FEDesc1);
              I_K1 = BeginIndex[n];
              J_K1 = GlobalNumbers + BeginIndex[n];

              // find the local edge of neigh on which cell is
              l=0;
              // while(neigh->GetJoint(l)->GetNeighbour(neigh)!=cell) l++;
              while(neigh->GetJoint(l) != joint) l++;

              Indices1 = FEDesc1_Obj->GetJointDOF(l);

              // OutPut("i= " << i << endl);
              // OutPut("j= " << j << endl);
              // OutPut("n= " << n << endl);
              // OutPut("l= " << l << endl);
              // OutPut("-------------" << endl);

              // OutPut("FEDesc0: " << FEDesc0 << endl);
              // OutPut("FEDesc1: " << FEDesc1 << endl);
              mapper = TFEDatabase2D::GetFE2DMapper(FEDesc0, FEDesc1);
              mapper->Map(GlobalNumbers, I_K0, I_K1, Indices0, Indices1,
                          j, l, FEDesc0_Obj, FEDesc1_Obj, Counter, 
                          VHN, HNNumbers);

            } // n>i
            
          } // n != -1
        } // end neigh
      } // no boundary joint
    } // endfor j

    // handle inner degrees of freedom
    k = FEDesc0_Obj->GetN_InnerDOF();
    Indices0 = FEDesc0_Obj->GetInnerDOF();
    for(j=0;j<k;j++)
    {
      Counter--;
      J_K0[Indices0[j]] = Counter;
    } // endfor j
  } // endfor i

#ifdef __MORTAR__
#ifdef __CONNECT_CROSSPOINTS__
  bool Connect_Crosspoints = false;
  const int *TmpCE, *TmpEV, *TmpnEoE, *TmpVE, *TmpnVEqoV, *TmpIndex;
  int N_, N2_, C_loc, J_loc, V_loc, V_inc, clip, v0, w0, v1, w1;
  int N_JointDOFs, **JointDOFs;
  TBaseCell *cell1D, *cell2, *cell3;
  TJoint *CurrJoint;
  TVertex *CurrVert;
  FE2D CurrElementID;
  TFE2D *CurrElement;
  TFEDesc2D *CurrDesc;
  JointType type;

  for (i=0;i<N_UsedElements;i++)
  {
    CurrElementID = UsedElements[i];
    if (CurrElementID != N_P1_2D_T_A && CurrElementID != N_Q1_2D_Q_A &&
        CurrElementID != N_Q1_2D_Q_M && CurrElementID != C_P0_2D_T_A &&
        CurrElementID < D_P1_2D_Q_A)
    {
      Connect_Crosspoints = true;
      break;
    }
  }

  if (Connect_Crosspoints)
  {
    // connect DOFs on crosspoints of mortar edges
    for (i=0;i<N_Cells;i++)
    {
      cell = Collection->GetCell(i);
      N_ = cell->GetN_Edges();
      for (j=0;j<N_;j++)
      {
        CurrJoint = cell->GetJoint(j);
        type = CurrJoint->GetType();
        if (type == MortarJoint || type == MortarBaseJoint)
        {
          if (type == MortarJoint)
            k = ((TMortarJoint *) CurrJoint)->GetMEdgeInColl();
          else
            k = ((TMortarBaseJoint *) CurrJoint)->GetMEdgeInColl();

          // consider only cells on the non-mortar side
          if (k < 0) continue;
          cell1D = MortarColl->GetCell(k);

          // bound cell of mortar edge
          if (cell1D->GetJoint(0)->GetNeighbour((int) 1) &&
              cell1D->GetJoint(1)->GetNeighbour((int) 1))
            continue;

          cell2 = cell;
          J_loc = j;
          while (cell3 = cell2->GetParent())
          {
            N2_ = cell3->GetN_Edges();
            for (k=0;k<N2_;k++)
              if (cell3->GetChild(k) == cell2) break;

            cell3->GetRefDesc()->GetChildEdge(TmpCE, MaxLen1);
            cell3->GetRefDesc()->GetNewEdgeOldEdge(TmpnEoE);

            J_loc = TmpnEoE[TmpCE[MaxLen1*k + J_loc]];

            cell2 = cell3;
          }

          // switch to mortar side
          CurrJoint = cell2->GetJoint(J_loc);
          cell3 = CurrJoint->GetNeighbour(cell2);

          N2_ = cell3->GetN_Edges();
          for (J_loc=0;J_loc<N2_;J_loc++)
            if (cell3->GetJoint(J_loc) == CurrJoint) break;

          // check whether cell3 belongs to the collection
          clip = cell3->GetClipBoard();

         if (clip >= 0 && clip < N_Cells)
            clip = Collection->GetCell(clip) != cell3 ? 1 : 0;
          else
            clip = 1;

          if (clip)
          {
            // select the right vertex
            if (cell1D->GetJoint(1)->GetNeighbour((int) 1))
              CurrVert = cell1D->GetVertex(0);
            else
              CurrVert = cell1D->GetVertex(1);

            cell3->GetShapeDesc()->GetEdgeVertex(TmpEV);
            V_inc = cell3->GetVertex(TmpEV[2*J_loc]) == CurrVert ? 0 : 1;

            while (cell3->ExistChildren())
            {
              cell3->GetShapeDesc()->GetEdgeVertex(TmpEV);
              cell3->GetRefDesc()->GetEdgeChild(TmpEC, TmpLen1, MaxLen2);
              cell3->GetRefDesc()->GetNewVertEqOldVert(TmpnVEqoV, TmpIndex);
              cell3->GetRefDesc()->GetVertexEdge(TmpVE, TmpLen1, MaxLen1);
              cell3->GetRefDesc()->GetNewEdgeOldEdge(TmpnEoE);
              cell3->GetRefDesc()->GetOldEdgeNewLocEdge(TmpoEnlE);

              V_loc = TmpnVEqoV[TmpEV[2*J_loc + V_inc]];
              N2_ = TmpLen1[V_loc];
              for (k=0;k<N2_;k++)
                if (TmpnEoE[TmpVE[MaxLen1*V_loc + k]] == J_loc) break;

              C_loc = TmpEC[MaxLen2*TmpVE[MaxLen1*V_loc + k]];
              N2_ = cell3->GetN_Edges();
              J_loc = TmpoEnlE[C_loc*N2_ + J_loc];

              cell3 = cell3->GetChild(C_loc);

              clip = cell3->GetClipBoard();

              if (clip >= 0 && clip < N_Cells)
                if (Collection->GetCell(clip) == cell3) break;
            }
          }
          else
            clip = cell3->GetClipBoard();

          // consistency check
          if (Collection->GetCell(clip) != cell3)
          {
            cerr << "ERROR in TFESpace2D mortar part: wrong cell"<< endl;
            exit (-1);
          }

          // get parameter of non-mortar cell
          CurrElementID = GetFE2D(i, cell);
          CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
          CurrDesc = CurrElement->GetFEDesc2D();

          JointDOFs = CurrDesc->GetJointDOF();
          N_JointDOFs = CurrDesc->GetN_JointDOF();

          w0 = BeginIndex[i] + JointDOFs[j][(!V_inc)*(N_JointDOFs-1)];
          while ((v0 = GlobalNumbers[w0]) > -1) w0 = v0;

          // get parameter of mortar cell
          CurrElementID = GetFE2D(clip, cell3);
          CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
          CurrDesc = CurrElement->GetFEDesc2D();

          JointDOFs = CurrDesc->GetJointDOF();
          N_JointDOFs = CurrDesc->GetN_JointDOF();

          w1 = BeginIndex[clip] + JointDOFs[J_loc][V_inc*(N_JointDOFs-1)];
          while ((v1 = GlobalNumbers[w1]) > -1) w1 = v1;

          // connect both DOFs
          if (v0 != v1)
            if (v0 > v1)
              GlobalNumbers[w1] = w0;
            else
              GlobalNumbers[w0] = w1;
        }
      }
    }
  }
#endif // __CONNECT_CROSSPOINTS__
#endif // __MORTAR__

  // find global numbers

  // do something with hanging nodes LATER
  l=0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
    l += BoundaryUpperBound[i];

  m = -SumLocDOF - DirichletUpperBound - l +FIRSTMARK;
  n = VHN->GetN_Elements();
  // OutPut("number of hanging nodes: " << n << endl);
  for(i=0;i<n;i++)
  {
    j=HNNumbers->GetElement(i);

    while( (k=GlobalNumbers[j]) > -1 )
      j=k;

    GlobalNumbers[j] = m;
    m--;
  }

  l = FIRSTMARK - DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundMark[i] = l;
    l -= BoundaryUpperBound[i];
  }
  InnerMark = l;
  SlaveMark = l - SumLocDOF;

  // for(i=0;i<N_DiffBoundNodeTypes;i++)
  //   OutPut(i << "   " << BoundMark[i] << endl);
  // OutPut("InnerMark: " << InnerMark << endl);
  // OutPut("SlaveMark: " << SlaveMark << endl);

  DirichletCounter=0;
  l=DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundCounter[i] = l;
    l += BoundaryUpperBound[i];
  }
  count = l;

  N_Dirichlet = 0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
    N_BoundaryNodes[i]=0;
  N_Inner = 0;
  N_Slave = 0;

  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);

    FEType0 = GetFE2D(i, cell);
    FE0 = TFEDatabase2D::GetFE2D(FEType0);
    FE0->GetFEDesc2D(FEDesc0, FEDesc0_Obj);
    // FEDesc0 = FE0->GetFEDesc2D_ID();
    // FEDesc0_Obj = FE0->GetFEDesc2D();
    // FEDesc0_Obj = TFEDatabase2D::GetFEDesc2D(FEDesc0);
    J_K0 = GlobalNumbers + BeginIndex[i];
    
    k=FEDesc0_Obj->GetN_DOF();
    for(j=0;j<k;j++)
    {
      l = J_K0[j];
      if (l < -1)
      {
        // OutPut(endl << "new node" << endl);
        if(l<=SlaveMark)
        {
          // hanging node
          J_K0[j] = -l;
          N_Slave++;
          // OutPut("slave" << endl);
          continue;
        }

        if(l<=InnerMark)
        {
          // inner node
          J_K0[j] = count;
          count++;
          N_Inner++;
          // OutPut("inner" << endl);
          continue;
        }

        for(m=N_DiffBoundNodeTypes-1;m>=0;m--)
        {
          if(l<=BoundMark[m])
          {
            J_K0[j] = BoundCounter[m];
            BoundCounter[m]++;
            N_BoundaryNodes[m]++;
            // OutPut("type " << m << endl);
            m = -2;
            break;
          }
        }
        
        if(m!=-2 && l<=FIRSTMARK) // no match in loop above
        {
          // Dirichlet nodes
          J_K0[j] = DirichletCounter;
          DirichletCounter++;
          N_Dirichlet++;
          // OutPut("Dirichlet" << endl);
          continue;
        }
        
      } // l < -1
      else
      {
        if (l >= 0)
        {
          J_K0[j] = GlobalNumbers[l];
        }
        else
        {
          // OutPut("J_K0[j]==-1" << endl);
        }
      } // l >= -1
    } // endfor j

  } // endfor i

  // OutPut("N_Inner: " << N_Inner << endl);
  // OutPut("N_Slave: " << N_Slave << endl);
  // OutPut("N_Dirichlet: " << N_Dirichlet << endl);
  // for(i=0;i<N_DiffBoundNodeTypes;i++)
    // OutPut(i << " N_BoundaryNodes: " << N_BoundaryNodes[i] << endl);

  // create real numbers
  l = 0; m = 0;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += N_BoundaryNodes[i];
    m += BoundaryUpperBound[i];
  }
  
  DirichletOffset = N_Inner + N_Slave + l - 0;
  SlaveOffset = N_Inner + l - SumLocDOF - DirichletUpperBound 
                - m + FIRSTMARK;
  InnerOffset = 0 - DirichletUpperBound - m;

  l = N_Inner; m = DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    BoundOffset[i] = l-m;
    l += N_BoundaryNodes[i];
    m += BoundaryUpperBound[i];
  }

  DirichletMark = DirichletUpperBound;
  InnerMark = SumLocDOF-FIRSTMARK;
  l = DirichletUpperBound;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += BoundaryUpperBound[i];
    BoundMark[i] = l;
  }

  // OutPut("DirichletMark: " << DirichletMark << endl);

  for(i=0;i<SumLocDOF;i++)
  {
    n=GlobalNumbers[i];
    // OutPut(i << "  " << n << endl);
    if(n<DirichletMark)
    {
      // Dirichlet node
      GlobalNumbers[i] += DirichletOffset;
      // OutPut("Diri" << endl);
    }
    else
    {
      if(n<BoundMark[N_DiffBoundNodeTypes-1])
      {
        // non Dirichlet boundary type
        for(m=0;m<N_DiffBoundNodeTypes;m++)
        {
          if(n<BoundMark[m])
          {
            // node of type m
            GlobalNumbers[i] += BoundOffset[m];
            // OutPut("type: " << m << endl);
            break;
          }
        } // endfor m
      }
      else
      {
        // no boundary node
        if(n<InnerMark)
        {
          // inner node
          GlobalNumbers[i] += InnerOffset;
          // OutPut("inner" << endl);
        }
        else
        {
          // slave node
          GlobalNumbers[i] += SlaveOffset;
          // OutPut("slave: " << endl);
        }
      }
    } // non Dirichlet node
  } // endfor i

/*
  // print for all elements for global numbers of their local dofs
  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    OutPut("cell number: " << i << endl);

    J_K0 = GlobalNumbers + BeginIndex[i];
    FEType0 = GetFE2D(i, cell);
    FE0 = TFEDatabase2D::GetFE2D(FEType0);
    k = FE0->GetSize();
    for(j=0;j<k;j++)
    {
      OutPut(j << ": " << " number: " << J_K0[j] << endl);
    }

    OutPut(endl);

  } // endfor i
*/

  // fill in information for hanging nodes
  HangingNodeArray = new THangingNode*[N_Slave];

  for(i=0;i<N_Slave;i++)
  {
    hn=VHN->GetElement(i);
    k=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    v=hn->GetDOF();
    for(j=0;j<k;j++)
      v[j] = GlobalNumbers[v[j]];

    HangingNodeArray[i] = hn;

    // OutPut(hn << endl);
  }
  delete VHN;
  delete HNNumbers;

  InnerBound = N_Inner;
  l = N_Inner;
  for(i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += N_BoundaryNodes[i];
    BoundaryNodesBound[i] = l;
  }
  HangingBound = l+N_Slave;
  DirichletBound = HangingBound + N_Dirichlet;

  N_ActiveDegrees = l;
  N_SlaveDegrees = N_Slave;
  N_DegreesOfFreedom = N_ActiveDegrees + N_Slave + N_Dirichlet;

  ActiveBound = N_ActiveDegrees;

  delete [] BoundaryUpperBound;
  delete [] BoundCounter;
  delete [] BoundMark;
  delete [] BoundOffset;
}

TFESpace2D::~TFESpace2D()
{
  delete [] BoundaryNodesBound;

  if (UsedElements)
    delete [] UsedElements;

  delete [] HangingNodeArray;

  if(AllElements)
    delete [] AllElements;

  if(ElementForShape)
    delete [] ElementForShape;
}

/** return position of all dofs */
void TFESpace2D::GetDOFPosition(double *x, double *y)
{
  int i,j,k;
  TBaseCell *cell;
  int N_Joints;
  TJoint *joint;
  JointType jointtype;
  FE2D FEid;
  int *DOF;
  TNodalFunctional2D *nf;
  double *xi, *eta;
  int N_Points;
  RefTrans2D RefTrans, *RefTransArray;
  int IsIsoparametric;
  BoundTypes bdtype;
  BF2DRefElements RefElement;
  TRefTrans2D *rt;
  double X[MaxN_BaseFunctions2D], Y[MaxN_BaseFunctions2D];
  double absdetjk[MaxN_BaseFunctions2D];

  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

  for(i=0;i<N_Cells;i++)
  {
    DOF = GlobalNumbers + BeginIndex[i];

    cell  = Collection->GetCell(i);
    FEid = GetFE2D(i, cell); 
    RefTrans = RefTransArray[FEid];

    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEid);

    nf = TFEDatabase2D::GetNodalFunctional2DFromFE2D(FEid);
    nf->GetPointsForAll(N_Points, xi, eta);

    N_Joints = cell->GetN_Joints();

    IsIsoparametric = false;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Joints;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == IsoInterfaceJoint ||
           jointtype == IsoBoundEdge)
          IsIsoparametric = true;
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

    rt = TFEDatabase2D::GetRefTrans2D(RefTrans);
    switch(RefTrans)
    {
      case TriaAffin:
       ((TTriaAffin *)rt)->SetCell(cell);
       ((TTriaAffin *)rt)->GetOrigFromRef(N_Points, xi, eta,
                                                X, Y, absdetjk);
      break;

      case TriaIsoparametric:
       ((TTriaIsoparametric *)rt)->SetCell(cell);
       ((TTriaIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta,
                                                X, Y, absdetjk);
      break;

      case QuadAffin:
       ((TQuadAffin *)rt)->SetCell(cell);
       ((TQuadAffin *)rt)->GetOrigFromRef(N_Points, xi, eta,
                                                X, Y, absdetjk);
      break;

      case QuadBilinear:
       ((TQuadBilinear *)rt)->SetCell(cell);
       ((TQuadBilinear *)rt)->GetOrigFromRef(N_Points, xi, eta,
                                                X, Y, absdetjk);
      break;

      case QuadIsoparametric:
       ((TQuadIsoparametric *)rt)->SetCell(cell);
       ((TQuadIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta,
                                                X, Y, absdetjk);
      break;
    } // endswitch RefTrans

    for(j=0;j<N_Points;j++)
    {
      k = DOF[j];
      x[k] = X[j];
      y[k] = Y[j];
    }

  } // endfor i
} // end GetDOFPosition

/** return position of one given dof */
void TFESpace2D::GetDOFPosition(int dof, double &x, double &y)
{
  int i,j,k;
  TBaseCell *cell;
  int N_Joints;
  TJoint *joint;
  JointType jointtype;
  FE2D FEid;
  int *DOF;
  TNodalFunctional2D *nf;
  double *xi, *eta;
  int N_Points;
  RefTrans2D RefTrans, *RefTransArray;
  int IsIsoparametric;
  BoundTypes bdtype;
  BF2DRefElements RefElement;
  TRefTrans2D *rt;
  int DOFFound;
  double absdetjk[1];

  if(dof > N_DegreesOfFreedom)
  {
    Error("dof number is larger than total number of degrees of freedom" << endl);
    x = -1; y = -1;
  }

  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

  for(i=0;i<N_Cells;i++)
  {
    DOF = GlobalNumbers + BeginIndex[i];
    k = BeginIndex[i+1] - BeginIndex[i];

    DOFFound = -1;
    for(j=0;j<k;j++)
    {
      if(DOF[j] == dof) 
      {
        DOFFound = j;
        break;
      } // endif
    } // endfor

    if(DOFFound>-1) // i.e. dof was found
    {
      //cout << "dof " << dof << " found in cell: " << i << endl;
      cell  = Collection->GetCell(i);
      FEid = GetFE2D(i, cell); 
      RefTrans = RefTransArray[FEid];
  
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEid);
  
      nf = TFEDatabase2D::GetNodalFunctional2DFromFE2D(FEid);
      nf->GetPointsForAll(N_Points, xi, eta);
  
      N_Joints = cell->GetN_Joints();
  
      IsIsoparametric = false;
      if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
      {
        for(j=0;j<N_Joints;j++)
        {
          joint = cell->GetJoint(j);
          jointtype = joint->GetType();
          if(jointtype == BoundaryEdge)
          {
            bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
            if(bdtype != Line)
              IsIsoparametric = true;
          }
          if(jointtype == InterfaceJoint)
          {
            bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
            if(bdtype != Line)
              IsIsoparametric = true;
          }
          if(jointtype == IsoInterfaceJoint ||
             jointtype == IsoBoundEdge)
            IsIsoparametric = true;
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
  
      rt = TFEDatabase2D::GetRefTrans2D(RefTrans);
      switch(RefTrans)
      {
        case TriaAffin:
         ((TTriaAffin *)rt)->SetCell(cell);
         ((TTriaAffin *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound,
                                                  &x, &y, absdetjk);
        break;
  
        case TriaIsoparametric:
         ((TTriaIsoparametric *)rt)->SetCell(cell);
         ((TTriaIsoparametric *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound,
                                                  &x, &y, absdetjk);
        break;
  
        case QuadAffin:
         ((TQuadAffin *)rt)->SetCell(cell);
         ((TQuadAffin *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound,
                                                  &x, &y, absdetjk);
        break;
  
        case QuadBilinear:
         ((TQuadBilinear *)rt)->SetCell(cell);
         ((TQuadBilinear *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound,
                                                  &x, &y, absdetjk);
        break;
  
        case QuadIsoparametric:
         ((TQuadIsoparametric *)rt)->SetCell(cell);
         ((TQuadIsoparametric *)rt)->GetOrigFromRef(1, xi+DOFFound, eta+DOFFound,
                                                  &x, &y, absdetjk);
        break;
      } // endswitch RefTrans

      break;
    } // endif DOFFound > -1
  } // endfor i
} // end GetDOFPosition


/** check if FE spaces lhs_space and rhs_space are equal*/
bool operator==(const TFESpace2D &lhs_space, const TFESpace2D &rhs_space)
{
  if(&lhs_space == &rhs_space) // compare pointers
    return true;
  if((lhs_space.N_DegreesOfFreedom == rhs_space.N_DegreesOfFreedom)
     && (lhs_space.N_UsedElements == rhs_space.N_UsedElements)
     && (lhs_space.BoundCondition == rhs_space.BoundCondition)
     && (lhs_space.N_ActiveDegrees == rhs_space.N_ActiveDegrees)
     //&& (lhs_space.OrderOfSpace == rhs_space.OrderOfSpace)
     && (lhs_space.Collection == rhs_space.Collection))
  {
    return true;
  }
  
  return false;
}

/** check if FE spaces lhs_space and rhs_space are not equal*/
bool operator!=(const TFESpace2D &lhs_space, const TFESpace2D &rhs_space)

{
  return !(lhs_space == rhs_space);
}
