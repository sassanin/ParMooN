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
   
// ======================================================================
// TimeDiscRout.C       07/12/18
//
// Volker John and Joachim Rang
//
// routines which are used for the different temporal discretizations
//
// ======================================================================

#include <Database.h>
#include <MooNMD_Io.h>
#ifdef __2D__
  #include <SquareMatrix2D.h>
  #include <Matrix2D.h>
  #include <AuxParam2D.h>
  #include <Assemble2D.h>
#endif
#include <TimeDiscRout.h>

#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <math.h>

/******************************************************************************/
/*                                                                            */
/* GENERAL                                                                    */
/*                                                                            */
/******************************************************************************/

void SetTimeDiscParameters(int increase_count)
{
  static int count=0;
  static int rb_type = TDatabase::TimeDB->RB_TYPE;

  static double theta=1-sqrt(double(2.0))/2;
  static double thetap=1-2*theta;
  static double alpha=thetap/(1-theta);
  static double beta=1-alpha;
  double gamma;
  double b4, c3;
  double c2, a22, a31, a32, a33;
  
  if (increase_count)
    count++;

  switch(TDatabase::TimeDB->TIME_DISC)
  {
    // case 0 - 4: theta-schemes without pressure correction for TNSE
    // case 5    : Rosenbrock schemes (p0 is needed)
    // case 6    : DIRK (p0 is needed)
    case 0: // FORWARD_EULER
      TDatabase::TimeDB->THETA1 = 0;
      TDatabase::TimeDB->THETA2 = 1;
      TDatabase::TimeDB->THETA3 = 1;
      TDatabase::TimeDB->THETA4 = 0;

      TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
        TDatabase::TimeDB->TIMESTEPLENGTH;

      if(count==1) count=0;
    break;

    case 1: // BACKWARD_EULER
      if(count==1) count=0;
      TDatabase::TimeDB->THETA1 = 1;
      TDatabase::TimeDB->THETA2 = 0;
      TDatabase::TimeDB->THETA3 = 0;
      TDatabase::TimeDB->THETA4 = 1;


      TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = 
        TDatabase::TimeDB->TIMESTEPLENGTH;

      if(count==1) count=0;
    break;

    case 2: // CRANK_NICOLSON
      TDatabase::TimeDB->THETA1 = 0.5;
      TDatabase::TimeDB->THETA2 = 0.5;
      TDatabase::TimeDB->THETA3 = 0.5;
      TDatabase::TimeDB->THETA4 = 0.5;

      TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = 
        TDatabase::TimeDB->TIMESTEPLENGTH;

      if(count==1) count=0;
    break;

    case 3: // FRACTIONAL_STEP WITH CANONICAL RHS
      switch(count)
      {
        case 1:
          TDatabase::TimeDB->THETA1 = alpha;
          TDatabase::TimeDB->THETA2 = beta;
          TDatabase::TimeDB->THETA3 = beta;
          TDatabase::TimeDB->THETA4 = alpha;

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = 
            theta * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
        case 2:
          TDatabase::TimeDB->THETA1 = beta;
          TDatabase::TimeDB->THETA2 = alpha;
          TDatabase::TimeDB->THETA3 = alpha;
          TDatabase::TimeDB->THETA4 = beta;

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = 
            thetap * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
        case 3:
          TDatabase::TimeDB->THETA1 = alpha;
          TDatabase::TimeDB->THETA2 = beta;
          TDatabase::TimeDB->THETA3 = beta;
          TDatabase::TimeDB->THETA4 = alpha;

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = 
            theta * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
      default:
          OutPut("Pass 1 as the argument in SetTimeDiscParameters() for FRACTIONAL_STEP scheme; !!!"<<endl);
          exit(4711);
      }
      if(count==3) count=0;
    break; // FRACTIONAL_STEP

    case 4: // FRACTIONAL_STEP WITH SIMPLIFED RHS
      switch(count)
      {
        case 1:
          TDatabase::TimeDB->THETA1 = alpha;
          TDatabase::TimeDB->THETA2 = beta;
          TDatabase::TimeDB->THETA3 = 1;
          TDatabase::TimeDB->THETA4 = 0;

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = 
            theta * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
        case 2:
          TDatabase::TimeDB->THETA1 = beta;
          TDatabase::TimeDB->THETA2 = alpha;
          TDatabase::TimeDB->THETA3 = 0;
          TDatabase::TimeDB->THETA4 = 1;

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
            thetap * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
        case 3:
          TDatabase::TimeDB->THETA1 = alpha;
          TDatabase::TimeDB->THETA2 = beta;
          TDatabase::TimeDB->THETA3 = 1;
          TDatabase::TimeDB->THETA4 = 0;

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
            theta * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
      }
      if (count==3) count=0;
    break; // FRACTIONAL_STEP
    
    case 9: // FRACTIONAL_STEP WITH 2 STEPS
      switch(count)
      {
        case 1:
          TDatabase::TimeDB->THETA2 = 0.1;
          TDatabase::TimeDB->THETA3 = 0.1;
          TDatabase::TimeDB->THETA1 = 1-(TDatabase::TimeDB->THETA2);
          TDatabase::TimeDB->THETA4 = (TDatabase::TimeDB->THETA1);

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = 
            1.0/3 * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
        case 2:
          TDatabase::TimeDB->THETA2 = 3.0/5;
          TDatabase::TimeDB->THETA3 = (TDatabase::TimeDB->THETA2);
          TDatabase::TimeDB->THETA1 = 1-(TDatabase::TimeDB->THETA2);
          TDatabase::TimeDB->THETA4 = (TDatabase::TimeDB->THETA1);

          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
            2.0/3 * TDatabase::TimeDB->TIMESTEPLENGTH;

        break;
      }
      if (count==2) count=0;
    break; // FRACTIONAL_STEP WITH 2 STEPS

    case 5: // Rosenbrock methods
      switch(rb_type)
      {
        case 0: // linear impl. Euler
          TDatabase::TimeDB->RB_GAMMA_I  =     5.00000000000000000000e-01;
          TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
          TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
          TDatabase::TimeDB->RB_M_I      =     2.00000000000000000000e+00;
          TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
          count = 0;
          break; // linear impl. Euler
        case 1: // ROS2
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.70710678118654701763e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     1.70710678118654701763e+00;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     8.78679656440357836900e-01;
              TDatabase::TimeDB->RB_MS_I     =     5.85786437626905076570e-01;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.70710678118654790580e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     1.70710678118654701763e+00;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     5.85786437626905076570e-01;
              TDatabase::TimeDB->RB_M_I      =     2.92893218813452538285e-01;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     5.85786437626905076570e-01;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.17157287525381059723e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.43145750507620028724e-01;
            break;
          }
          if (count==2)  count = 0;
          break; // ROS2
        case 2: // ROWDA3
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.23672704529659016615e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.05935616764594087158e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     6.04455284065558817730e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.99999999999999955591e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.60599625219532926579e+00;
              TDatabase::TimeDB->RB_M_I      =     2.25006773096964485248e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.69401431934652568767e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.60599625219532926579e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -8.87404441065782090270e-01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.68460566009349044236e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     6.37978879934487963510e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.99999999999999955591e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.60599625219532926579e+00;
              TDatabase::TimeDB->RB_M_I      =    -2.09251404439032034910e-01;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.60599625219532926579e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -2.39874797163503465924e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.68460566009349044236e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -5.26372237156212907649e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==3)  count = 0;
          break; // ROWDA3
        case 3: // ROS3P
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.11324865405187134471e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -2.11324865405187106715e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.26794919243112280682e+00;
              TDatabase::TimeDB->RB_M_I      =     5.77350269189625731059e-01;
              TDatabase::TimeDB->RB_MS_I     =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.26794919243112280682e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.60769515458673617481e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.60769515458673617481e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.07735026918962573106e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.26794919243112280682e+00;
              TDatabase::TimeDB->RB_M_I      =     4.22649730810374213430e-01;
              TDatabase::TimeDB->RB_MS_I     =     4.22649730810374213430e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.26794919243112280682e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     3.46410161513775483044e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.60769515458673617481e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     1.73205080756887741522e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==3)  count = 0;
          break; // ROS3P
        case 4: // GRK4A
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     3.95000000000000017764e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     3.95000000000000017764e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.84568324040708864331e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.89400194217594530777e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -3.72672395483999996380e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     3.95000000000000017764e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     4.38000000000000000444e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.10886075949367080007e+00;
              TDatabase::TimeDB->RB_M_I      =     1.36979689436154350446e-01;
              TDatabase::TimeDB->RB_MS_I     =    -5.10131175670240022413e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.10886075949367080007e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.92018840239705212980e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     2.80724242909790078215e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     6.62919654459999951879e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     3.95000000000000017764e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.69999999999999995559e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.35694476970006894234e+00;
              TDatabase::TimeDB->RB_M_I      =     5.40602212188314390495e-01;
              TDatabase::TimeDB->RB_MS_I     =     9.31597444379746786325e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.37708526198442093857e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.05558868604935085500e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     6.40885667370853084890e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.85011498891139242184e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -3.35181726766864285239e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.68383541496555033667e-01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     3.55389155740999995725e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     3.95000000000000017764e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.69999999999999995559e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.35694476970006894234e+00;
              TDatabase::TimeDB->RB_M_I      =     8.05218958546835450463e-01;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.37708526198442093857e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -3.24956722543420051252e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     6.40885667370853084890e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.85011498891139242184e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -3.41099762766929570645e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.68383541496555033667e-01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     1.69967830598942470921e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // GRK4A
        case 5: // GRK4T
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.31000000000000010880e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.31000000000000010880e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     3.95750374664077231301e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.59011493022284389198e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -3.96296677524429971640e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     2.31000000000000010880e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     4.62000000000000021760e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     4.62489238836330684990e+00;
              TDatabase::TimeDB->RB_MS_I     =    -1.39340789777533657912e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     5.07167533877631626638e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     8.65800865800865793176e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     5.50778939578913995234e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.31000000000000010880e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.80208333333333370341e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     3.61179418775468608072e-01;
              TDatabase::TimeDB->RB_M_I      =     6.17477263875011006533e-01;
              TDatabase::TimeDB->RB_MS_I     =     9.79647124616345066350e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.52470820737311463233e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -6.02015272865080941500e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.65567939706882327755e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     4.16352878859764619079e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -1.59750684672702480960e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.80239341497733605024e+01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -5.53509845705278877293e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     2.31000000000000010880e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.80208333333333370341e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     3.61179418775468608072e-01;
              TDatabase::TimeDB->RB_M_I      =     1.28261294526903757429e+00;
              TDatabase::TimeDB->RB_MS_I     =    -6.84615728561236536187e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.52470820737311463233e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.85634361868609865098e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.65567939706882327755e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     4.16352878859764619079e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     8.50538085817980338277e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.80239341497733605024e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     2.08407513602318772428e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // GRK4T
        case 6: // Shampine 1982
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.11111519999999996955e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.79629629629591991424e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.50000000000000000000e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     4.99998400000000009502e-01;
              TDatabase::TimeDB->RB_MS_I     =     3.05555555555599989148e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     8.00000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.00000000000000000000e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.41999999999999992895e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.99999999999999977796e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.67999999999999993783e+00;
              TDatabase::TimeDB->RB_M_I      =     2.31479999999999991322e-01;
              TDatabase::TimeDB->RB_MS_I     =     2.31481481482000012173e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.91999999999999992895e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.48800000000000007816e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.79999999999999982236e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.39999999999999991118e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -2.39999999999999991118e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.79999999999999982236e-01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.16000000000000005884e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.99999999999999977796e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.67999999999999993783e+00;
              TDatabase::TimeDB->RB_M_I      =     1.15739999999999998437e+00;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.91999999999999992895e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     8.96000000000000018652e-01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.79999999999999982236e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.39999999999999991118e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     4.31999999999999995115e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.79999999999999982236e-01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     4.00000000000000022204e-01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // Shampine 1982
        case 7: // Velhuizen1
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.25709999999999993969e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.25709999999999993969e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     4.28919943191591368503e+00;
              TDatabase::TimeDB->RB_MS_I     =    -4.67730849944162851983e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -4.59900000000000030997e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     2.25709999999999993969e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     4.51419999999999987939e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     5.03603558423743269401e+00;
              TDatabase::TimeDB->RB_MS_I     =     4.69074035867010952217e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     5.33320204404498188211e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     8.86092773913428821686e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     5.17759999999999998010e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.25709999999999993969e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.75600000000000044942e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.34009371686359324283e-01;
              TDatabase::TimeDB->RB_M_I      =     6.08555550199671113631e-01;
              TDatabase::TimeDB->RB_MS_I     =     7.91468196423286496355e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.81214060202617588402e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -6.10037092228335797728e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     5.16962389579347014745e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     4.57813123033981650423e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -1.80469118855169541327e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     2.02832450061575322309e+01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -3.80500000000000004885e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     2.25709999999999993969e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.75600000000000044942e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.34009371686359324283e-01;
              TDatabase::TimeDB->RB_M_I      =     1.35594346728102421729e+00;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.81214060202617588402e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.54029031068434063556e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     5.16962389579347014745e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     4.57813123033981650423e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     9.44364979019760752976e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     2.02832450061575322309e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     1.98841872308338873943e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // Velhuizen1
        case 8: // L-stable method
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     5.72815999999999991843e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.72815999999999991843e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.25556617004763682033e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.97402207064192558583e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.76917681090700007474e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     5.72815999999999991843e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.14563199999999998369e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.87055150561301375411e-01;
              TDatabase::TimeDB->RB_MS_I     =     2.14289450758942900954e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     7.13765047493353321784e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.49152258316806785032e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     7.59294430611999970893e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.72815999999999991843e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.55214947223999977233e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.63350264059105132475e+00;
              TDatabase::TimeDB->RB_M_I      =     4.35311872314439129994e-01;
              TDatabase::TimeDB->RB_MS_I     =     3.27099068126588621297e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.86794821820585266181e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -2.58092694955929014000e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.11581025366953578981e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.34445577614801253796e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -6.51630418482569195859e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.09286014382980345427e-01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.04894507469999995197e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.72815999999999991843e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.55214947223999977233e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.63350264059105132475e+00;
              TDatabase::TimeDB->RB_M_I      =     1.09350773869619555256e+00;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.86794821820585266181e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.13711475285956975512e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.11581025366953578981e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.34445577614801253796e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     3.21469570357977130204e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.09286014382980345427e-01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     6.94965924392651257513e-01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // L-stable method
        case 9: // ROWDAIND2
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.99999999999999988898e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.99999999999999988898e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.83076923076923070433e+00;
              TDatabase::TimeDB->RB_MS_I     =    -8.62563780642876665183e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.87820512820512819374e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.99999999999999988898e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.66666666666666674068e+00;
              TDatabase::TimeDB->RB_M_I      =     2.39999999999999991118e+00;
              TDatabase::TimeDB->RB_MS_I     =     3.70998116760828651195e+01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.66666666666666674068e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.24643874643874630337e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     5.55555555555555535818e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     2.99999999999999988898e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.30769230769230782041e-01;
              TDatabase::TimeDB->RB_M_I      =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     3.38041431261770242145e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.83076923076923070433e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.22678062678062680391e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -4.23931623931623935420e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.39999999999999991118e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     4.26666666666666642982e+01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     8.00000000000000000000e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.99999999999999982944e-18;
              TDatabase::TimeDB->RB_GAMMA_II =     2.99999999999999988898e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.30769230769230782041e-01;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =    -1.54237288135593217930e+01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.83076923076923070433e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     5.82462804685026439011e-02;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -4.23931623931623935420e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.39999999999999991118e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     3.25925925925925907833e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     8.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -3.70370370370370349811e-01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // ROWDAIND2
        case 10: // DAE34
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =    -5.58287037393940721586e+00;
              TDatabase::TimeDB->RB_MS_I     =    -5.54305555779703684749e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.50000000000000000000e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.34027777916749979426e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.29675926035499999855e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -4.00000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     0.00000000000000000000e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.40740740499999994473e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     5.00000000000000041668e-32;
              TDatabase::TimeDB->RB_M_I      =     7.00000000950000034194e-01;
              TDatabase::TimeDB->RB_MS_I     =     9.00000000699999969100e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =    -5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -3.46296296200000019994e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -3.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     1.50000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.00000000000000000000e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.56250000000000000000e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     7.50000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =    -1.04166666125000001442e-01;
              TDatabase::TimeDB->RB_M_I      =     4.35555555639999969486e+00;
              TDatabase::TimeDB->RB_MS_I     =     4.17777777840000030807e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =    -2.32291666612499980715e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -9.08333333224999961431e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -9.54166666449999922861e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.09375000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     3.31250000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     2.75000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[2]  =     1.12500000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     2.25000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     2.25000000000000000000e+00;
            break;
            case 5: // 5. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.11111111149999997050e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     5.00000000030000002482e+00;
              TDatabase::TimeDB->RB_M_I      =     4.00000000000000022204e-01;
              TDatabase::TimeDB->RB_MS_I     =     4.00000000000000022204e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.55787038450000003831e-01;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.24822530883802471635e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.01157407749999994628e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -2.14120370400000009647e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     4.70138888927499998260e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     2.07175925919999981417e+00;
              TDatabase::TimeDB->RB_A_IJ[2]  =     2.75000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -1.16666666650000006022e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     5.50000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[3]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[3]  =     1.58518518519999993543e+01;
              TDatabase::TimeDB->RB_S_IJ[3]  =     4.00000000000000000000e+00;
            break;
          }
          if (count==5)  count = 0;
          break; // DAE34
        case 11: // ros3w
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.12548994752803022124e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.05820860602268629336e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     7.99373358398527078528e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.66666666666666629659e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.52952024018602772415e+00;
              TDatabase::TimeDB->RB_M_I      =     9.78349767388509317101e-01;
              TDatabase::TimeDB->RB_MS_I     =     3.74919173845941255951e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.52952024018602772415e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.91339906955403749045e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.50914824770808619903e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -6.17619939953493068963e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.66666666666666629659e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.52952024018602772415e+00;
              TDatabase::TimeDB->RB_M_I      =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_MS_I     =     3.18650050038755794368e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.52952024018602772415e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.06053924969355861663e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.50914824770808619903e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     8.09559354637497841090e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==3)  count = 0;
          break; // ros3w
        case 12: // ros3Dw
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.91931393562181629164e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.96207579734685988448e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     9.11187903127909315515e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.71733043016918007773e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.07431947473051159214e+00;
              TDatabase::TimeDB->RB_MS_I     =     7.92844107656903385184e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -2.50195979011212088494e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.58856072055808361654e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -6.48565537178490969517e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     3.82132943717632289626e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     8.76720107786980284992e-01;
              TDatabase::TimeDB->RB_M_I      =     1.25734648470373233664e+00;
              TDatabase::TimeDB->RB_MS_I     =     8.90631832214905938550e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     8.76720107786980284992e-01;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.52731248285624232608e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     2.01144172475739324568e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     5.64853010944563860285e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==3)  count = 0;
          break; // ros3Dw
        case 13: // ros3Pw
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.63397459621556140341e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.99444650053486505215e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.57735026918962573106e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.94228634059947813384e-01;
              TDatabase::TimeDB->RB_MS_I     =     6.54700538379251573140e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.53589838486224561365e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     2.53589838486224561365e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -5.28312163512967766787e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     7.88675134594812865529e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     6.33974596215561403412e-01;
              TDatabase::TimeDB->RB_M_I      =     1.07179676972449078320e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.07179676972449078320e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     6.33974596215561403412e-01;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.62740473580835520728e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     8.03847577293368087403e-01;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     2.74519052838329002952e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==3)  count = 0;
          break; // ros3Pw
        case 14: // ros34PW1a
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     6.15383214653102150749e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.15835114880704193041e+01;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.78292094614482721227e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     2.21878746765328616064e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     5.09052051067020450148e+00;
              TDatabase::TimeDB->RB_M_I      =    -8.36423375973235905256e-01;
              TDatabase::TimeDB->RB_MS_I     =     4.90957657168138439374e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     5.09052051067020450148e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.16790812312282881180e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.16790812312282881180e+01;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     3.33333333333333314830e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =    -8.61479212095767610258e-01;
              TDatabase::TimeDB->RB_MS_I     =    -8.61479212095767610258e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     7.10095263654306196877e-01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     4.16546077167549849696e-02;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.25807049614762522793e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.78370379319140659469e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     4.36216969948677313340e+00;
              TDatabase::TimeDB->RB_M_I      =     2.29428036027904180827e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.29428036027904180827e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.00517369636786391141e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.19795577622266034012e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.00035701597476247571e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.93164702379441555191e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     4.80544005238949689662e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.90957657168138439374e-01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -1.43504930216552795130e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     2.63186118578106453825e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // ros34PW1a
        case 15: // ros34PW1b
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     5.22582761233093862074e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.03942797401713331595e+01;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.78292094614482721227e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     2.21878746765328616064e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     5.09052051067020450148e+00;
              TDatabase::TimeDB->RB_M_I      =    -5.56971148154164819033e-01;
              TDatabase::TimeDB->RB_MS_I     =     7.06548277884253672632e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     5.09052051067020450148e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.16790812312282881180e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.16790812312282881180e+01;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -2.46541900496934207609e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     2.21878746765328616064e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     5.09052051067020450148e+00;
              TDatabase::TimeDB->RB_M_I      =     3.57979469353645218810e-01;
              TDatabase::TimeDB->RB_MS_I     =     3.57979469353645218810e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     5.09052051067020450148e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.64057326467366806355e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.16790812312282881180e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     2.77268164715849529944e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -8.05529997906369588101e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.55392337535788360725e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     3.92438391154034249553e+00;
              TDatabase::TimeDB->RB_M_I      =     1.72337398521064066870e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.72337398521064066870e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.97628111010787321788e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     8.38103960500475864137e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.25014103693966855957e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.77268164715849495250e-02;
              TDatabase::TimeDB->RB_C_IJ[1]  =     8.48328409199343047575e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.27226180967637575447e-01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     2.29428036027904180827e-01;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -2.87009860433105612465e-01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     5.26372237156212952058e-01;
            break;
          }
          if (count==4)  count = 0;
          break; // ros34PW1b
        case 16: // ros34PW2
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     4.18476048231916042397e+00;
              TDatabase::TimeDB->RB_MS_I     =     3.90701053467119185925e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.71733043016918007773e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =    -2.85192017355495874842e-01;
              TDatabase::TimeDB->RB_MS_I     =     1.11804787782050318867e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.58856072055808361654e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.58856072055808361654e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -4.13333376233886440332e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     7.31579957788852430767e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.67844952912873446316e+00;
              TDatabase::TimeDB->RB_M_I      =     2.29428036027904180827e+00;
              TDatabase::TimeDB->RB_MS_I     =     5.21650232611490682899e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.41921731745576473749e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.18476048231916042397e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.25608241840666678968e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -2.59232211672969725669e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -2.85192017355495874842e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =    -5.94751371992993771443e-01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.00000000000000007154e-17;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.91339906955403771249e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.18476048231916042397e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     6.36817920012835703147e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.17316616301130984823e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -2.85192017355495874842e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     6.79562094446683584437e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     5.59055033583923874363e-02;
              TDatabase::TimeDB->RB_A_IJ[2]  =     2.29428036027904180827e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -2.87009860433105590261e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     5.26372237156212907649e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // ros34PW2
        case 17: // ros34PW3
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.06857902130162885079e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     1.06857902130162885079e+00;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.54687580769067034581e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.34054870555389404885e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.44696658076125284076e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     1.06857902130162885079e+00;
              TDatabase::TimeDB->RB_ALPHA_I  =     2.51554560206288169155e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.35410348876090891324e+00;
              TDatabase::TimeDB->RB_M_I      =     5.52780970443755181876e-01;
              TDatabase::TimeDB->RB_MS_I     =     5.50176738302150081239e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.35410348876090891324e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.20302237067446027297e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     2.20302237067446027297e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -7.71476248531343222758e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     1.06857902130162885079e+00;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.25777280103144084578e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.17705174438045423457e+00;
              TDatabase::TimeDB->RB_M_I      =     9.20521963049926061906e-01;
              TDatabase::TimeDB->RB_MS_I     =     9.13798568129854826836e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.12745185174323347965e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.75006011499346758598e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.99091673084859821508e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     7.01866670643065848623e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     8.40856963108111976624e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     6.56822431146109275701e-01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -2.93711722610809189415e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     1.06857902130162885079e+00;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.28886400515720422888e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.00688480056840878873e+00;
              TDatabase::TimeDB->RB_M_I      =     7.20167498325635424550e-01;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.65735413669071252052e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     3.00778719731558963346e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     2.13382405288825172107e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     3.79983651191293847482e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     7.07077462561825020870e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     4.96890356906715180418e-01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     7.67753793376751203503e-01;
              TDatabase::TimeDB->RB_C_IJ[2]  =     1.18743627492743542007e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     7.18481065107899663502e-01;
            break;
          }
          if (count==4)  count = 0;
          break; // ros34PW3
        case 18: // RODAS3
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.00000000000000000000e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.50000000000000000000e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -4.00000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     0.00000000000000000000e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     0.00000000000000000000e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[2]  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     2.66666666666666651864e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =     2.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // RODAS3
        case 19: // RODAS
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.22122450922664116391e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.22122450922664094186e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.04300000000000003819e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     3.86000000000000009770e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.54400000000000003908e+00;
              TDatabase::TimeDB->RB_M_I      =     6.01913448128862960118e+00;
              TDatabase::TimeDB->RB_MS_I     =     6.01913448128862960118e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.54400000000000003908e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     5.66880000000000006111e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     6.17600000000000015632e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.03499999999999980904e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     2.09999999999999992228e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     8.07577091656891954408e-01;
              TDatabase::TimeDB->RB_M_I      =     1.25370833293208701065e+01;
              TDatabase::TimeDB->RB_MS_I     =     1.25370833293208701065e+01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     9.46678528081582593146e-01;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.43009335683387517335e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.65702247895389831456e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.55701169898328417585e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     2.06359915709191488187e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.02280467959331367034e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -3.62000000000001834199e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.30000000000000115463e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.93149530386442469521e+00;
              TDatabase::TimeDB->RB_M_I      =    -6.87886036105895048998e-01;
              TDatabase::TimeDB->RB_MS_I     =    -6.87886036105895048998e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     3.31482518706852102852e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.07352905815137886214e-01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.05651238005194265668e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.89612401597220081584e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     9.59456225102335480415e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.07691601022100975626e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     9.98641913997781682788e-01;
              TDatabase::TimeDB->RB_C_IJ[2]  =     2.04702861480961608720e+01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     3.99456765599112673115e+00;
            break;
            case 5: // 5. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -5.00142462167563195408e-16;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     9.99999999999999444888e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     9.99999999999999333866e-01;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.22122450922664049777e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -7.49644331396764762587e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -6.35636504793053269680e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     6.01913448128862871300e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     1.02468043146435210389e+01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.46486913464229253634e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     1.25370833293208701065e+01;
              TDatabase::TimeDB->RB_C_IJ[2]  =     3.39999035281990487078e+01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     3.88149166317527516412e+01;
              TDatabase::TimeDB->RB_A_IJ[3]  =    -6.87886036105895048998e-01;
              TDatabase::TimeDB->RB_C_IJ[3]  =    -1.17089089320616004386e+01;
              TDatabase::TimeDB->RB_S_IJ[3]  =    -2.75154414442358019599e+00;
            break;
            case 6: // 6. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.12973866372989562024e-16;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     9.99999999999999000799e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.00000000000000022204e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.22122450922664094186e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -8.08324679592152151031e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -7.49644331396764584952e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     6.01913448128862960118e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     7.98113298806489446235e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.02468043146435228152e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =     1.25370833293208701065e+01;
              TDatabase::TimeDB->RB_C_IJ[2]  =     3.15215943287437099229e+01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     3.39999035281990487078e+01;
              TDatabase::TimeDB->RB_A_IJ[3]  =    -6.87886036105895048998e-01;
              TDatabase::TimeDB->RB_C_IJ[3]  =    -1.63193054312313599041e+01;
              TDatabase::TimeDB->RB_S_IJ[3]  =    -1.17089089320616004386e+01;
              TDatabase::TimeDB->RB_A_IJ[4]  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[4]  =     6.05881823883405346010e+00;
              TDatabase::TimeDB->RB_S_IJ[4]  =     4.00000000000000000000e+00;
            break;
          }
          if (count==6)  count = 0;
          break; // RODAS
        case 20: // RODASP
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =    -7.17045496242302604628e+00;
              TDatabase::TimeDB->RB_MS_I     =    -7.17045496242302515810e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     7.50000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     3.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =    -4.74163667148178546995e+00;
              TDatabase::TimeDB->RB_MS_I     =    -4.74163667148178546995e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     3.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.20000000000000000000e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.20000000000000000000e+01;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -2.35039999999999416058e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     2.09999999999999992228e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     8.39999999999999968914e-01;
              TDatabase::TimeDB->RB_M_I      =    -1.63100263133097094226e+01;
              TDatabase::TimeDB->RB_MS_I     =    -1.63100263133097094226e+01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.83103679348675907335e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     8.79179517394703502475e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     7.32414717394703629338e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     4.95518396743379496705e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     2.20786558697351820157e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.98207358697351820886e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -3.62000000000003499534e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     6.29999999999999893419e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.22107428571428400232e+00;
              TDatabase::TimeDB->RB_M_I      =    -1.06200404411140092442e+00;
              TDatabase::TimeDB->RB_MS_I     =    -1.06200404411140092442e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.30437658269266920641e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -1.08179305685715299035e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =     7.49037998156431861219e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -5.24927524574300141680e-02;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -6.78027061142826603657e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =    -4.75682755861467576608e-01;
              TDatabase::TimeDB->RB_A_IJ[2]  =    -1.17679876183278198098e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -1.95348594464241003266e+01;
              TDatabase::TimeDB->RB_S_IJ[2]  =    -4.70719504733112792394e+00;
            break;
            case 5: // 5. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.15006745446399882837e-15;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.00000000000000222045e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =    -7.17045496242302515810e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -3.41909500674967645750e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -3.31756975032767584821e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -4.74163667148178546995e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -1.54967115372596335732e+01;
              TDatabase::TimeDB->RB_S_IJ[1]  =    -1.59537223481944021586e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =    -1.63100263133097094226e+01;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -5.47476087596413023562e+01;
              TDatabase::TimeDB->RB_S_IJ[2]  =    -4.94930656966754582982e+01;
              TDatabase::TimeDB->RB_A_IJ[3]  =    -1.06200404411140092442e+00;
              TDatabase::TimeDB->RB_C_IJ[3]  =    -1.41600539214853391456e+01;
              TDatabase::TimeDB->RB_S_IJ[3]  =    -4.24801617644560369769e+00;
            break;
            case 6: // 6. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.07005333413456860114e-15;
              TDatabase::TimeDB->RB_GAMMA_II =     2.50000000000000000000e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     9.99999999999998889777e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.00000000000000022204e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =    -7.17045496242302515810e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -3.46260583093053284642e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -3.41909500674967645750e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -4.74163667148178546995e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -1.53008497611447342734e+01;
              TDatabase::TimeDB->RB_S_IJ[1]  =    -1.54967115372596317968e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =    -1.63100263133097094226e+01;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -5.69995557866266722158e+01;
              TDatabase::TimeDB->RB_S_IJ[2]  =    -5.47476087596413023562e+01;
              TDatabase::TimeDB->RB_A_IJ[3]  =    -1.06200404411140092442e+00;
              TDatabase::TimeDB->RB_C_IJ[3]  =    -1.84080700979309490606e+01;
              TDatabase::TimeDB->RB_S_IJ[3]  =    -1.41600539214853391456e+01;
              TDatabase::TimeDB->RB_A_IJ[4]  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[4]  =     5.71428571428571707713e+00;
              TDatabase::TimeDB->RB_S_IJ[4]  =     4.00000000000000000000e+00;
            break;
          }
          if (count==6)  count = 0;
          break; // RODASP
        case 21: // ROSI2P1
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.83337514882783159109e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.74778579810360490399e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -6.41334784915409961137e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_M_I      =     3.95341788699960128284e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.77038079363523226384e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.63186118578106453825e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     2.63186118578106453825e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.45563307177156459060e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     7.50000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.72071027020928135620e+00;
              TDatabase::TimeDB->RB_M_I      =    -1.21522771421847153306e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.57316038155499193785e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.78576458718195851816e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     3.01131047554100383934e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.09704462045337791665e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     4.42124760965982632754e-01;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -3.34203214637756484962e-01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.01435815587731981147e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -1.35847884055847717422e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.82552237199314371097e+00;
              TDatabase::TimeDB->RB_M_I      =     1.16541744745930642146e+00;
              TDatabase::TimeDB->RB_MS_I     =     3.43556220548094648493e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.50623951095167241121e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     3.74359059430178309213e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     7.24139738721911996322e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     4.55821087656518209030e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     1.08994123815715848735e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     1.23102185539136641523e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =    -1.37361554490644821591e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     1.71836543021444199120e+00;
              TDatabase::TimeDB->RB_S_IJ[2]  =    -3.15145916725285824000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // ROSI2P1
        case 22: // ROSI2P2
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     2.80734818821136977718e+00;
              TDatabase::TimeDB->RB_MS_I     =     4.20084258522924747226e-01;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -6.41334784915409961137e-02;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     5.00000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_M_I      =     3.48693217206767203109e+00;
              TDatabase::TimeDB->RB_MS_I     =    -5.94329934171131846199e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.63186118578106453825e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     2.63186118578106453825e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     1.20849664917601007375e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.29428036027904180827e+00;
              TDatabase::TimeDB->RB_M_I      =    -1.00000000000000005597e-32;
              TDatabase::TimeDB->RB_MS_I     =     3.60559439940373371858e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.80734818821136977718e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -4.97638997727638887625e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     6.44084381267829630957e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     3.48693217206767203109e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -6.18104102134040900296e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     8.00000000000000000000e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     3.00000000000000006055e-17;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.29428036027904180827e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =    -3.34373891242489040820e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.80734818821136977718e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     1.76105018434538296290e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     6.44084381267829630957e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     3.48693217206767203109e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     6.54597265243972703530e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     8.00000000000000000000e+00;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[2]  =     5.39706236424998597734e-01;
              TDatabase::TimeDB->RB_S_IJ[2]  =     0.00000000000000000000e+00;
            break;
          }
          if (count==4)  count = 0;
          break; // ROSI2P2
        case 23: // ROSI2Pw
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     3.28685750853367641966e+00;
              TDatabase::TimeDB->RB_MS_I     =     3.39043774753558357915e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.71733043016918007773e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.54503745896544728566e+01;
              TDatabase::TimeDB->RB_MS_I     =     5.86507841261477924633e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.58856072055808361654e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.58856072055808361654e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -4.18867127163060515294e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     7.50000000000000000000e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     1.72071027020928135620e+00;
              TDatabase::TimeDB->RB_M_I      =    -1.50445558288074323627e+01;
              TDatabase::TimeDB->RB_MS_I     =    -4.96246395067967771553e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.63034046718533942588e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.56739138878289008261e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     3.74045811443148146225e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -9.03698030239419025644e-02;
              TDatabase::TimeDB->RB_C_IJ[1]  =     6.83107605436873155380e-02;
              TDatabase::TimeDB->RB_S_IJ[1]  =    -2.07333664240115456145e-01;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.99999999999999989549e-17;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =    -1.32075118456237405873e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.60825116305745285938e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     3.28685750853367641966e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.99296365068222680605e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -4.89987118381377939613e+00;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.54503745896544728566e+01;
              TDatabase::TimeDB->RB_C_IJ[1]  =     4.26472578877039012468e+01;
              TDatabase::TimeDB->RB_S_IJ[1]  =     3.13005430424391839495e+01;
              TDatabase::TimeDB->RB_A_IJ[2]  =    -1.50445558288074323627e+01;
              TDatabase::TimeDB->RB_C_IJ[2]  =    -4.36510246478373389323e+01;
              TDatabase::TimeDB->RB_S_IJ[2]  =    -3.45164289671544750604e+01;
            break;
          }
          if (count==4)  count = 0;
          break; // ROSI2Pw
        case 24: // ROSI2PW
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     0.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     4.33856072055808361654e+00;
              TDatabase::TimeDB->RB_MS_I     =     3.31383144448515931657e+00;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     8.71733043016918007773e-01;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_M_I      =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_MS_I     =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_A_IJ[0]  =     2.00000000000000000000e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     4.58856072055808361654e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     4.58856072055808361654e+00;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_GAMMA_I  =     6.56544000523295512295e+00;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =    -1.59874671679705415706e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =    -3.66797319340808058996e+00;
              TDatabase::TimeDB->RB_M_I      =    -5.95593546374978410896e-02;
              TDatabase::TimeDB->RB_MS_I     =     8.47038665840773273563e-03;
              TDatabase::TimeDB->RB_A_IJ[0]  =    -5.50195979011212088494e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =    -4.83965596116219671785e+01;
              TDatabase::TimeDB->RB_S_IJ[0]  =    -1.26230382894992363418e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =    -1.83398659670404029498e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =    -1.61321865372073247613e+01;
              TDatabase::TimeDB->RB_S_IJ[1]  =    -4.20767942983307907667e+00;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_GAMMA_I  =    -4.00000000000000028617e-18;
              TDatabase::TimeDB->RB_GAMMA_II =     4.35866521508459003886e-01;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_SIGMA_I  =     2.91339906955403726840e+00;
              TDatabase::TimeDB->RB_M_I      =     1.00000000000000000000e+00;
              TDatabase::TimeDB->RB_MS_I     =     2.60825116305745285938e-01;
              TDatabase::TimeDB->RB_A_IJ[0]  =     4.33856072055808361654e+00;
              TDatabase::TimeDB->RB_C_IJ[0]  =     2.31911621201704987172e+00;
              TDatabase::TimeDB->RB_S_IJ[0]  =     1.20845224961108961281e+01;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_C_IJ[1]  =     1.14714018013952090413e+00;
              TDatabase::TimeDB->RB_S_IJ[1]  =     3.34207713346653090269e+00;
              TDatabase::TimeDB->RB_A_IJ[2]  =    -5.95593546374978410896e-02;
              TDatabase::TimeDB->RB_C_IJ[2]  =     7.45075552140237323817e-02;
              TDatabase::TimeDB->RB_S_IJ[2]  =    -1.36645857615705773602e-01;
            break;
          }
          if (count==4)  count = 0;
          break; // ROSI2Pw
          default:
            OutPut("This Rosenbrock-scheme is not implemented" << endl);
            exit(4711);
          break;
      }
      TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
        TDatabase::TimeDB->TIMESTEPLENGTH;

      break;
      
      case 6: // DIRK-methods (p0 is needed)
      switch(rb_type)
      {
        case 0: // impl. Euler
          TDatabase::TimeDB->RB_A_IJ[0]  =     1.0;
          TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
          TDatabase::TimeDB->RB_M_I      =     1.0;
          TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
          if (count==1) count = 0;
          break; // impl. Euler
        case 1: // GL1, mitpoint-rule
          TDatabase::TimeDB->RB_A_IJ[0]  =     0.5;
          TDatabase::TimeDB->RB_ALPHA_I  =     0.5;
          TDatabase::TimeDB->RB_M_I      =     1.0;
          TDatabase::TimeDB->RB_MS_I     =     0.00000000000000000000e+00;
          if (count==1) count = 0;
          break; // GL1, mitpoint-rule
        case 2: // Crank-Nicolson (pressure corrected)
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     0.5;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.5;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.5;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     0.5;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==2)  count = 0;
          break; // Crank-Nicolson (pressure corrected)
        case 3: // fractional-step-theta-scheme (pressure corrected)
        
          b4 = (4*sqrt(2.0) - 5)/(sqrt(2.0)+1);
          b4 = (4*sqrt(2.0) - 5)/3/(sqrt(2.0)+1);
        
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     theta*beta;
              TDatabase::TimeDB->RB_MS_I     =     1 - b4 - (2 - 5*sqrt(2.0)/4 + b4*sqrt(2.0)/2) - (2 - sqrt(2.0) - b4);
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     theta*beta;
              TDatabase::TimeDB->RB_A_IJ[1]  =     theta*alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     theta;
              TDatabase::TimeDB->RB_M_I      =     (1-theta)*alpha;
              TDatabase::TimeDB->RB_MS_I     =     2 - sqrt(2.0) - b4;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     theta*beta;
              TDatabase::TimeDB->RB_A_IJ[1]  =     (1-theta)*alpha;
              TDatabase::TimeDB->RB_A_IJ[2]  =     thetap*beta;
              TDatabase::TimeDB->RB_ALPHA_I  =     1-theta;
              TDatabase::TimeDB->RB_M_I      =     (1-theta)*beta;
              TDatabase::TimeDB->RB_MS_I     =     2 - 5*sqrt(2.0)/4 + b4*sqrt(2.0)/2;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     theta*beta;
              TDatabase::TimeDB->RB_A_IJ[1]  =     (1-theta)*alpha;
              TDatabase::TimeDB->RB_A_IJ[2]  =     (1-theta)*beta;
              TDatabase::TimeDB->RB_A_IJ[3]  =     theta*alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     theta*alpha;
              TDatabase::TimeDB->RB_MS_I     =     b4;
            break;
          }
          if (count==4)  count = 0;
          break; // fractional-step-theta-scheme (pressure corrected)
        case 4: // DIRK3
          gamma = 0.5 + sqrt(3.0)/6;
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     1.0 + 1.0/12/gamma/(2*gamma-1) - (3*gamma-1)/3/(2*gamma-1);
              TDatabase::TimeDB->RB_MS_I     =     5.0/12 + sqrt(3.0)/12;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     gamma;
              TDatabase::TimeDB->RB_A_IJ[1]  =     gamma;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*gamma;
              TDatabase::TimeDB->RB_M_I      =     -1.0/12/gamma/(2*gamma-1);
              TDatabase::TimeDB->RB_MS_I     =     3.0/4 + sqrt(3.0)/12;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0 + 1.0/12/gamma/(2*gamma-1) - (3*gamma-1)/3/(2*gamma-1);
              TDatabase::TimeDB->RB_A_IJ[1]  =     -1.0/12/gamma/(2*gamma-1);
              TDatabase::TimeDB->RB_A_IJ[2]  =     (3*gamma-1)/3/(2*gamma-1);
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     (3*gamma-1)/3/(2*gamma-1);
              TDatabase::TimeDB->RB_MS_I     =     -1.0/6 - sqrt(3.0)/6;
            break;
          }
          if (count==3)  count = 0;
          break; // DIRK3
        case 5: // DIRK3L
          gamma = 1.0 - sqrt(2.0)/2;
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     1.0 - 1.0/4/(1-gamma) - (1-2*gamma)/2/(1-gamma);
              TDatabase::TimeDB->RB_MS_I     =     0.5 - sqrt(2.0)/8;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     gamma;
              TDatabase::TimeDB->RB_A_IJ[1]  =     gamma;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*gamma;
              TDatabase::TimeDB->RB_M_I      =     1.0/4/(1-gamma);
              TDatabase::TimeDB->RB_MS_I     =     0.5 - sqrt(2.0)/8;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0 - 1.0/4/(1-gamma) - (1-2*gamma)/2/(1-gamma);
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.0/4/(1-gamma);
              TDatabase::TimeDB->RB_A_IJ[2]  =     (1-2*gamma)/2/(1-gamma);
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     (1-2*gamma)/2/(1-gamma);
              TDatabase::TimeDB->RB_MS_I     =     sqrt(2.0)/4;
            break;
          }
          if (count==3)  count = 0;
          break; // DIRK3L
        case 6: // embedded method of DIRK34L
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     1.0 - 1.0724862707343697992166039413780 - 0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_M_I      =     1.0724862707343697992166039413780;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0 - 1.0724862707343697992166039413780 - 0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.0724862707343697992166039413780;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==3)  count = 0;
          break; // embedded method of DIRK34L
        case 7: // DIRK34L
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     1.0 - 0.76852982927695365400080128104719 - 0.096648360979159732288805255170684 - 0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_MS_I     =     1.0 - 1.0724862707343697992166039413780 - 0.15898389998867654678259477757479;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_M_I      =     0.76852982927695365400080128104719;
              TDatabase::TimeDB->RB_MS_I     =     1.0724862707343697992166039413780;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0 - 1.0724862707343697992166039413780
                                                       - 0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.0724862707343697992166039413780;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     0.096648360979159732288805255170684;
              TDatabase::TimeDB->RB_MS_I     =     0.15898389998867654678259477757479;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0 - 0.76852982927695365400080128104719 - 0.096648360979159732288805255170684 - 0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.76852982927695365400080128104719;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.096648360979159732288805255170684;
              TDatabase::TimeDB->RB_A_IJ[3]  =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     0.15898389998867654678259477757479;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==4)  count = 0;
          break; // DIRK34L
        case 8: // embedded method of DIRK44L
          gamma = 5.0/6;
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     1.0 - (3*gamma-2)/3/(3*gamma-1) - (6*gamma-1)/6/(3*gamma-1);
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.25;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.25;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.5;
              TDatabase::TimeDB->RB_M_I      =     (3*gamma-2)/3/(3*gamma-1);
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0 - (3*gamma-2)/3/(3*gamma-1) - (6*gamma-1)/6/(3*gamma-1);
              TDatabase::TimeDB->RB_A_IJ[1]  =     (3*gamma-2)/3/(3*gamma-1);
              TDatabase::TimeDB->RB_A_IJ[2]  =     (6*gamma-1)/6/(3*gamma-1);
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     (6*gamma-1)/6/(3*gamma-1);
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==3)  count = 0;
          break; // embedded method of DIRK44L
        case 9: // DIRK44L
          gamma = 5.0/6;
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     1.0/6;
              TDatabase::TimeDB->RB_MS_I     =     1.0 - (3*gamma-2)/3/(3*gamma-1) - (6*gamma-1)/6/(3*gamma-1);
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.25;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.25;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.5;
              TDatabase::TimeDB->RB_M_I      =     2.0/3;
              TDatabase::TimeDB->RB_MS_I     =     (3*gamma-2)/3/(3*gamma-1);
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0 - (3*gamma-2)/3/(3*gamma-1) - (6*gamma-1)/6/(3*gamma-1);
              TDatabase::TimeDB->RB_A_IJ[1]  =     (3*gamma-2)/3/(3*gamma-1);
              TDatabase::TimeDB->RB_A_IJ[2]  =     (6*gamma-1)/6/(3*gamma-1);
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     1.0/6 -gamma;
              TDatabase::TimeDB->RB_MS_I     =     (6*gamma-1)/6/(3*gamma-1);
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     1.0/6;
              TDatabase::TimeDB->RB_A_IJ[1]  =     2.0/3;
              TDatabase::TimeDB->RB_A_IJ[2]  =     1.0/6 -gamma;
              TDatabase::TimeDB->RB_A_IJ[3]  =     gamma;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     gamma;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==4)  count = 0;
          break; // DIRK44L
        case 10: // DIRK 2 (Williams et al)
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     7.0/18;
              TDatabase::TimeDB->RB_MS_I     =     0.5;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.5;
              TDatabase::TimeDB->RB_A_IJ[1]  =     0.5;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     1.0/3;
              TDatabase::TimeDB->RB_MS_I     =     0.5;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     5.0/8;
              TDatabase::TimeDB->RB_A_IJ[1]  =     3.0/8;
              TDatabase::TimeDB->RB_A_IJ[2]  =     0.5;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.5;
              TDatabase::TimeDB->RB_M_I      =     -2.0/9;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 4: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     7.0/18;
              TDatabase::TimeDB->RB_A_IJ[1]  =     1.0/3;
              TDatabase::TimeDB->RB_A_IJ[2]  =     -2.0/9;
              TDatabase::TimeDB->RB_A_IJ[3]  =     0.5;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     0.5;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==4)  count = 0;
          break; // DIRK 2 (Williams et al)
        case 11: // Alexander: SDIRK(3,2) - order 3
          alpha = 1 - sqrt(2.0)/2*cos(1.0/3*atan(sqrt(2.0)/4)) + sqrt(6.0)/2*sin(1.0/3*atan(sqrt(2.0)/4));
          c3 = 18.0/13*alpha*alpha - 2*alpha + 14.0/13;
        
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     (18*alpha*c3 - 12*alpha*alpha*c3 - 3*c3 - 12*alpha + 12*alpha*alpha +2)/12/alpha/c3;
              TDatabase::TimeDB->RB_MS_I     =     (18*alpha*c3 - 12*alpha*alpha*c3 - 3*c3 - 12*alpha + 12*alpha*alpha +2)/12/alpha/c3 + (24*c3-71+94*alpha)/500/alpha;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*alpha;
              TDatabase::TimeDB->RB_M_I      =     (2-3*c3+6*alpha*c3-6*alpha)/12/alpha/(2*alpha-c3);
              TDatabase::TimeDB->RB_MS_I     =     (2-3*c3+6*alpha*c3-6*alpha)/12/alpha/(2*alpha-c3) + (24*c3-71)/500/alpha;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (6*alpha*c3 - c3*c3 - 4*alpha*alpha)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     -c3*(2*alpha-c3)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[2]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     c3;
              TDatabase::TimeDB->RB_M_I      =     -(6*alpha*alpha + 1 - 6*alpha)/3/c3/(2*alpha-c3);
              TDatabase::TimeDB->RB_MS_I     =     -(6*alpha*alpha + 1 - 6*alpha)/3/c3/(2*alpha-c3) + 12.0/125;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (18*alpha*c3 - 12*alpha*alpha*c3 - 3*c3 - 12*alpha + 12*alpha*alpha +2)/12/alpha/c3;
              TDatabase::TimeDB->RB_A_IJ[1]  =     (2-3*c3+6*alpha*c3-6*alpha)/12/alpha/(2*alpha-c3);
              TDatabase::TimeDB->RB_A_IJ[2]  =     -(6*alpha*alpha + 1 - 6*alpha)/3/c3/(2*alpha-c3);
              TDatabase::TimeDB->RB_A_IJ[3]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     alpha;
              TDatabase::TimeDB->RB_MS_I     =     alpha - 71.0/250;
            break;
          }
          if (count==4)  count = 0;
          break;
        case 12: // Alexander: SDIRK(3,4) - method 4
        
          alpha = 1 - sqrt(2.0)/2*cos(1.0/3*atan(sqrt(2.0)/4)) + sqrt(6.0)/2*sin(1.0/3*atan(sqrt(2.0)/4));
          c3 = 13.0/18*alpha*alpha - 2*alpha + 14.0/13;
        
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     (12*alpha*c3 - 4*alpha - 2*c3 + 1)/24/alpha/c3;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*alpha;
              TDatabase::TimeDB->RB_M_I      =     (2*c3-1)/24/alpha/(2*alpha-c3)*(2*alpha-1);
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (6*alpha*c3 - c3*c3 - 4*alpha*alpha)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     -c3*(2*alpha-c3)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[2]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     c3;
              TDatabase::TimeDB->RB_M_I      =     (4*alpha-1)/12/c3/(2*alpha-c3)/(c3-1);
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (18*alpha*c3 - 12*alpha*alpha*c3 - 3*c3 - 12*alpha + 12*alpha*alpha +2)/12/alpha/c3;
              TDatabase::TimeDB->RB_A_IJ[1]  =     (2-3*c3+6*alpha*c3-6*alpha)/12/alpha/(2*alpha-c3);
              TDatabase::TimeDB->RB_A_IJ[2]  =     -(6*alpha*alpha + 1 - 6*alpha)/3/c3/(2*alpha-c3);
              TDatabase::TimeDB->RB_A_IJ[3]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     (3+12*alpha*c3-4*c3-8*alpha)/12/(2*alpha-1)/(c3-1);
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==4)  count = 0;
          break;
        case 13: // Alexander: SDIRK(3,4) - order 3
          alpha = 1 - sqrt(2.0)/2*cos(1.0/3*atan(sqrt(2.0)/4)) + sqrt(6.0)/2*sin(1.0/3*atan(sqrt(2.0)/4));
          c3 = 13.0/18*alpha*alpha - 2*alpha + 14.0/13;

          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     (18*alpha*c3 - 12*alpha*alpha*c3 - 3*c3 - 12*alpha + 12*alpha*alpha +2)/12/alpha/c3;
              TDatabase::TimeDB->RB_MS_I     =     (12*alpha*c3 - 4*alpha - 2*c3 + 1)/24/alpha/c3;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*alpha;
              TDatabase::TimeDB->RB_M_I      =     (2-3*c3+6*alpha*c3-6*alpha)/12/alpha/(2*alpha-c3);
              TDatabase::TimeDB->RB_MS_I     =     (2*c3-1)/24/alpha/(2*alpha-c3)*(2*alpha-1);
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (6*alpha*c3 - c3*c3 - 4*alpha*alpha)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     -c3*(2*alpha-c3)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[2]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     c3;
              TDatabase::TimeDB->RB_M_I      =     -(6*alpha*alpha + 1 - 6*alpha)/3/c3/(2*alpha-c3);
              TDatabase::TimeDB->RB_MS_I     =     (4*alpha-1)/12/c3/(2*alpha-c3)/(c3-1);
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (18*alpha*c3 - 12*alpha*alpha*c3 - 3*c3 - 12*alpha + 12*alpha*alpha +2)/12/alpha/c3;
              TDatabase::TimeDB->RB_A_IJ[1]  =     (2-3*c3+6*alpha*c3-6*alpha)/12/alpha/(2*alpha-c3);
              TDatabase::TimeDB->RB_A_IJ[2]  =     -(6*alpha*alpha + 1 - 6*alpha)/3/c3/(2*alpha-c3);
              TDatabase::TimeDB->RB_A_IJ[3]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     alpha;
              TDatabase::TimeDB->RB_MS_I     =     (3+12*alpha*c3-4*c3-8*alpha)/12/(2*alpha-1)/(c3-1);
            break;
          }
          if (count==4)  count = 0;
          break;
        case 14: // ESDIRK 3/2 (a) with 4 stages and gamma = 0.43 (Kvaerno)
          alpha = 0.4358665215;
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     (6*alpha-1)/12/alpha;
              TDatabase::TimeDB->RB_MS_I     =     (-4*alpha*alpha+6*alpha-1)/4/alpha;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     2*alpha;
              TDatabase::TimeDB->RB_M_I      =     -1/(24*alpha-12)/alpha;
              TDatabase::TimeDB->RB_MS_I     =     (1-2*alpha)/4/alpha;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (-4*alpha*alpha+6*alpha-1)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     (1-2*alpha)/4/alpha;
              TDatabase::TimeDB->RB_A_IJ[2]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     1;
              TDatabase::TimeDB->RB_M_I      =     (-6*alpha*alpha+6*alpha-1)/(6*alpha-3);
              TDatabase::TimeDB->RB_MS_I     =     alpha;
            break;
            case 4: // 4. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     (6*alpha-1)/12/alpha;
              TDatabase::TimeDB->RB_A_IJ[1]  =     -1/(24*alpha-12)/alpha;
              TDatabase::TimeDB->RB_A_IJ[2]  =     (-6*alpha*alpha+6*alpha-1)/(6*alpha-3);
              TDatabase::TimeDB->RB_A_IJ[3]  =     alpha;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     alpha;
              TDatabase::TimeDB->RB_MS_I     =     0;
            break;
          }
          if (count==4)  count = 0;
          break;
        case 15: // This method is only A-stable but not strongly A-stable
        
          c2 = 0.5 + sqrt(1.0/6);
          a22 = 0.5*c2;
          a32 = 0.5;
          a33 = 0.5 - a32*c2;
          a31 = 1.0 - a32 - a33;
        
          switch(count)
          {
            case 1: // 1. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     0.0;
              TDatabase::TimeDB->RB_ALPHA_I  =     0.0;
              TDatabase::TimeDB->RB_M_I      =     a31;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 2: // 2. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     a22;
              TDatabase::TimeDB->RB_A_IJ[1]  =     a22;
              TDatabase::TimeDB->RB_ALPHA_I  =     c2;
              TDatabase::TimeDB->RB_M_I      =     a32;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
            case 3: // 3. substep
              TDatabase::TimeDB->RB_A_IJ[0]  =     a31;
              TDatabase::TimeDB->RB_A_IJ[1]  =     a32;
              TDatabase::TimeDB->RB_A_IJ[2]  =     a33;
              TDatabase::TimeDB->RB_ALPHA_I  =     1.0;
              TDatabase::TimeDB->RB_M_I      =     a33;
              TDatabase::TimeDB->RB_MS_I     =     0.0;
            break;
          }
          if (count==3)  count = 0;
          break;
      }
      TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
          TDatabase::TimeDB->TIMESTEPLENGTH;
      break;
    case 7:
    case 8:
      // Extrapolation Euler and CN
        TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
          TDatabase::TimeDB->TIMESTEPLENGTH;
      break;
   }
}


int GetN_SubSteps()
{
  int ret;

  switch(TDatabase::TimeDB->TIME_DISC)
  {
    case 0: // BACKWARD_EULER
    case 1: // FORWARD_EULER
    case 2: // CRANK_NICOLSON
      ret = 1;
    break;

    case 3: // FRACTIONAL_STEP
    case 4: // FRACTIONAL_STEP
      ret = 3;
    break;
    case 5:
      switch(TDatabase::TimeDB->RB_TYPE)
      {
        case 0:
          ret = 1;
          break;
        case 1:
          ret = 2;
          break;
        case 2:
          ret = 3;
          break;
        case 3:
          ret = 3;
          break;
        case 4:
          ret = 4;
          break;
        case 5:
          ret = 4;
          break;
        case 6:
          ret = 4;
          break;
        case 7:
          ret = 4;
          break;
        case 8:
          ret = 4;
          break;
        case 9:
          ret = 4;
          break;
        case 10:
          ret = 5;
          break;
        case 11:
          ret = 3;
          break;
        case 12:
          ret = 3;
          break;
        case 13:
          ret = 3;
          break;
        case 14:
          ret = 4;
          break;
        case 15:
          ret = 4;
          break;
        case 16:
          ret = 4;
          break;
        case 17:
          ret = 4;
          break;
        case 18:
          ret = 4;
          break;
        case 19:
          ret = 6;
          break;
        case 20:
          ret = 6;
          break;
        case 21:
          ret = 4;
          break;
        case 22:
          ret = 4;
          break;
        case 23:
          ret = 4;
          break;
        case 24:
          ret = 4;
          break;
        default:
          cout << "This RB-TYPE does not exist!" << endl;
          exit(4711);
      }
    break;
    case 6:
    switch(TDatabase::TimeDB->RB_TYPE)
      {
        case 0:
          ret = 1;
          break;
        case 1:
          ret = 1;
          break;
        case 2:
          ret = 2;
          break;
        case 3:
          ret = 4;
          break;
        case 4:
          ret = 3;
          break;
        case 5:
          ret = 3;
          break;
        case 6:
          ret = 3;
          break;
        case 7:
          ret = 4;
          break;
        case 8:
          ret = 3;
          break;
        case 9:
          ret = 4;
          break;
        case 10:
          ret = 4;
          break;
        case 11:
          ret = 4;
          break;
        case 12:
          ret = 4;
          break;
        case 13:
          ret = 4;
          break;
        case 14:
          ret = 4;
          break;
        case 15:
          ret = 3;
          break;
        default:
          cout << "This DIRK-TYPE (" << (TDatabase::TimeDB->RB_TYPE) << ") does not exist!" << endl;
          exit(4711);
        }
        break;
    case 7:
    case 8:
      ret = 1;
      break;
      default:
          OutPut("Time discretization not implemented !!!"<<endl);
          exit(4711);
  } // endswitch

  return ret;
}

/******************************************************************************/
/*                                                                            */
/* ROSENBROCK                                                                 */
/*                                                                            */
/******************************************************************************/

void AllocateAuxiliaryVectorsRB(int &rb_order, int &RB_s,
                                double* &RB_m, double* &RB_ms,
                                double* &RB_alpha, double* &RB_gamma,
                                double* &RB_sigma, double* &RB_RHS_YN,
                                double* &old_sol_rbU, double* &rb_mein,
                                double* &rb_diff, double* &sol_tilde,                                
                                double* &B1, double* &B2, double* &dq,                                
                                double* RB_A, double* RB_C, double *RB_S,
                                int N_Unknowns)
{
  int i, i1, i2, il, j;
  double RB_gamma_ii, val;

  switch(TDatabase::TimeDB->RB_TYPE)
  {
    case 0:
    case 1:
      rb_order = 2;
     break;
    case 2:
    case 3:
    case 11:
    case 12:
    case 13:
      rb_order = 3;
     break;
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
    case 21:
    case 22:
    case 23:
    case 24:
      rb_order = 4;
     break;
    case 10:
      rb_order = 5;
     break;
    case 19:
    case 20:
      rb_order = 6;
      break;
    default: 
      OutPut("RB_TYPE " << TDatabase::TimeDB->RB_TYPE << " not implemented!\n");
      exit(4711);
      break;
  }

  RB_s = rb_order;
  RB_m  = new double[rb_order+2];
  RB_ms = new double[rb_order+2];
  RB_alpha = new double[rb_order+2];
  RB_gamma = new double[rb_order+2];
  RB_sigma = new double[rb_order+2];
  RB_RHS_YN = new double[rb_order+2];
  old_sol_rbU = new double[(rb_order+1)*N_Unknowns];
  rb_mein = new double[(rb_order+1)*N_Unknowns];
  
  memset(RB_A, 0, 100*SizeOfDouble); 
  memset(RB_C, 0, 100*SizeOfDouble); 
  memset(RB_S, 0, 100*SizeOfDouble); 

  for (il=1;il<=rb_order;il++) 
  {
    SetTimeDiscParameters();
    old_sol_rbU[il] =  old_sol_rbU[0] + il*N_Unknowns;
    rb_mein[il] =  rb_mein[0] + il*N_Unknowns;
    
    RB_m[il-1]  = TDatabase::TimeDB->RB_M_I;
    RB_ms[il-1] = TDatabase::TimeDB->RB_MS_I;
    
    for (i2=1;i2<il;++i2)
    {
      RB_C[(il-1)*10+i2-1] = TDatabase::TimeDB->RB_C_IJ[i2-1];
      RB_A[(il-1)*10+i2-1] = TDatabase::TimeDB->RB_A_IJ[i2-1];
      RB_S[(il-1)*10+i2-1] = TDatabase::TimeDB->RB_S_IJ[i2-1];
    }
    RB_alpha[il-1] = TDatabase::TimeDB->RB_ALPHA_I;
    RB_gamma[il-1] = TDatabase::TimeDB->RB_GAMMA_I;
    RB_sigma[il-1] = TDatabase::TimeDB->RB_SIGMA_I;
  }
  
  RB_gamma_ii = TDatabase::TimeDB->RB_GAMMA_II;
  
  for (i2=0;i2<rb_order;++i2)
  {
    RB_C[(il-1)*10+il-1] = 1/RB_gamma_ii;
  }
  
  // rechte Seite auswerten
  RB_RHS_YN[0] = 1;
  for (i1=1;i1<rb_order;++i1) 
  {
    val = 0;
    for (i2=0;i2<rb_order;++i2)
    {
      val += (RB_A[(i1-1)*10+i2] - RB_A[i1*10+i2])*(RB_A[(i1-1)*10+i2] - RB_A[i1*10+i2]);
    }
    if (val < 1.0e-10)
    {
       RB_RHS_YN[i1] = 0;
    } else {
       RB_RHS_YN[i1] = 1;
    }
    cout << "AUSW: " << i1 << " : " <<RB_RHS_YN[i1] << endl;
  }
  cout << "N_Unknowns: " << N_Unknowns << endl;
  rb_diff = new double[N_Unknowns];
  sol_tilde = new double[N_Unknowns];
  B1 = new double[N_Unknowns];
  B2 = new double[N_Unknowns];
  dq = new double[N_Unknowns];

  printf("Hallo\n");
  printf("    case %d: // \n", (TDatabase::TimeDB->RB_TYPE));
  printf("      *s = %d;\n", rb_order);
  printf("      *p = %d;\n", rb_order);
  for (i=0;i<RB_s;i++)
  {
    for (j=0;j<RB_s;j++)
    {
      printf("      C[%d][%d] = %30.20e;\n", i, j, RB_C[i*10+j]);
    }
  }
  for (i=0;i<RB_s;i++)
  {
    for (j=0;j<RB_s;j++)
    {
      printf("      A[%d][%d] = %30.20e;\n", i, j, RB_A[i*10+j]);
    }
  }
  for (i=0;i<RB_s;i++)
  {
    for (j=0;j<RB_s;j++)
    {
      printf("      S[%d][%d] = %30.20e;\n", i, j, RB_S[i*10+j]);
    }
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      gamma[%d] = %30.20e;\n", i, RB_gamma[i]);
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      alpha[%d] = %30.20e;\n", i, RB_alpha[i]);
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      sigma[%d] = %30.20e;\n", i, RB_sigma[i]);
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      m[%d] = %30.20e;\n", i, RB_m[i]);
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      ms[%d] = %30.20e;\n", i, RB_ms[i]);
  }
  
  printf("    break;\n");
}

#ifdef __2D__
void AssembleRHS_RB_DIRK(TFESpace2D **fesp, TFEFunction2D **fefct, TFESpace2D **ferhs,
TFESpace2D **USpaces, TFEFunction2D **U1Array, TFEFunction2D **U2Array,
TDiscreteForm2D *DiscreteForm, TDiscreteForm2D *DiscreteFormRHS,
double **RHSs, double *rhs, TAuxParam2D *aux,
BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **my_BoundValues,
int mg_type, int mg_level, int N_U, int N_Unknowns)
{
  int N_FESpaces, N_Rhs;

  fesp[0] = USpaces[mg_level-1];

  fefct[0] = U1Array[mg_level-1];
  fefct[1] = U2Array[mg_level-1];

  // current rhs
  ferhs[0] = USpaces[mg_level-1];
  ferhs[1] = USpaces[mg_level-1];

  DiscreteForm = DiscreteFormRHS;

  aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
  // initialize array
  N_Rhs = 2;
  memset(RHSs[0], 0, N_Unknowns*SizeOfDouble);

  RHSs[0] = rhs;
  RHSs[1] = rhs + N_U;
  RHSs[2] = rhs + 2*N_U;

  N_FESpaces = 1;
  Assemble2D(N_FESpaces, fesp,
    0, NULL,
    0, NULL,
    N_Rhs, RHSs, ferhs,
    DiscreteForm,
    BoundaryConditions,
    my_BoundValues,
    aux);

  delete aux;
}



#endif  


/******************************************************************************/
/*                                                                            */
/* DIRK                                                                       */
/*                                                                            */
/******************************************************************************/

void AllocateAuxiliaryVectorsDIRK(int &rb_order, int &RB_s,
                                  int &stiff_acc1, int &stiff_acc2,
                                  double* &RB_m, double* &RB_ms,
                                  double* &RB_alpha, 
                                  double* &RB_RHS_YN,
                                  double* &old_sol_rbU, double* &rb_mein,
                                  double* &old_sol_rbK, double* &old_rhs,
                                  double* &rb_diff, double* &sol_tilde,
                                  double* &B1, double* &B2,
                                  double* RB_A,
                                  int N_Unknowns)
{
  int i, i1, i2, il, j;
  double val;

  OutPut("DIRK-Methods" << endl);
  // set order of the methods for stepsize-selection
  switch(TDatabase::TimeDB->RB_TYPE)
  {
    case 0:
    case 1:
      rb_order = 1;
      break;
    case 2:
    case 3:
    case 6:
    case 8:
    case 15:
      rb_order = 2;
      break;
    case 4:
    case 5:
    case 7:
    case 9:
    case 10:
    case 14:
      rb_order = 3;
      break;
    case 11:
    case 12:
    case 13:
      rb_order = 4;
      break;
    default: 
      OutPut("RB_TYPE " << TDatabase::TimeDB->RB_TYPE << " not implemented!\n");
      exit(4711);
  }

  // allocate vectors and matrices
  RB_s = GetN_SubSteps();
  RB_m  = new double[RB_s+2];
  RB_ms = new double[RB_s+2];
  RB_alpha = new double[RB_s+2];
  RB_RHS_YN = new double[RB_s+2];
  // stores the values k_i for stepsize selction, if embedded is not stiffly accurate
  old_sol_rbK = new double[(RB_s+1)*N_Unknowns];
  // stores the values U_i
  old_sol_rbU = new double[(RB_s+1)*N_Unknowns];
  // stores the old rhs
  old_rhs = new double[(RB_s+1)*N_Unknowns];
  rb_mein = new double[(RB_s+1)*N_Unknowns];
  B2 = new double[N_Unknowns];

  memset(RB_A, 0, 100*SizeOfDouble); 
  for (il=1;il<=RB_s;il++)
  {
    SetTimeDiscParameters();
    old_sol_rbK[il] =  old_sol_rbK[0] + il*N_Unknowns;
    old_sol_rbU[il] =  old_sol_rbU[0] + il*N_Unknowns;
    old_rhs[il] =  old_rhs[0] + il*N_Unknowns;
    rb_mein[il] =  rb_mein[0] + il*N_Unknowns;

    RB_m[il-1]  = TDatabase::TimeDB->RB_M_I;
    RB_ms[il-1] = TDatabase::TimeDB->RB_MS_I;

    for (i2=1;i2<=il;++i2)
    {
      RB_A[(il-1)*10+(i2-1)] = TDatabase::TimeDB->RB_A_IJ[i2-1];
    }
    RB_alpha[il-1] = TDatabase::TimeDB->RB_ALPHA_I;
  }

  // are the method1 and method2 stiffly accurate?
  //  if yes : stiff_acc=1 else =0
  stiff_acc1 = 1;
  stiff_acc2 = 1;
  
  for (il=1;il<=RB_s;il++)
  {
    if (fabs(RB_m[il-1]  - RB_A[(RB_s-1)*10+il-1])>1.0E-6) stiff_acc1 = 0;
    if (fabs(RB_ms[il-1] - RB_A[(RB_s-2)*10+il-1])>1.0E-6) stiff_acc2 = 0;
  }
  
  OutPut("DIRK(m):  stiffly-accurate : " << stiff_acc1 << endl);
  OutPut("DIRK(ms): stiffly-accurate : " << stiff_acc2 << endl);
  
  OutPut("DIRK-coefficients: " << endl);
  for (il=1;il<=RB_s;il++)
  {
    for (i2=1;i2<=RB_s;i2++)
    {
      OutPut("DIRK: a[" << il << "," << i2 << "] = " << RB_A[(il-1)*10+i2-1]
             << endl);
    }
    OutPut("DIRK: alpha[" << il << "] = " << RB_alpha[il-1] << endl);        
    OutPut("DIRK: b[" << il << "] = " << RB_m[il-1] << endl);        
    OutPut("DIRK: bs[" << il << "] = " << RB_ms[il-1] << endl);        
  }
  OutPut("END OF DIRK-coefficients: " << endl);
  
  // rechte Seite auswerten
  RB_RHS_YN[0] = 1;
  for (i1=1;i1<RB_s;++i1) 
  {
    val = 0;
    for (i2=0;i2<RB_s;++i2)
    {
      val += (RB_A[(i1-1)+i2] - RB_A[i1*10+i2])*(RB_A[(i1-1)*10+i2] - RB_A[i1*10+i2]);
    }
    if (val < 1.0e-10)
    {
       RB_RHS_YN[i1] = 0;
    } else
    {
       RB_RHS_YN[i1] = 1;
    }
  }
  
  rb_diff = new double[N_Unknowns];
  sol_tilde = new double[N_Unknowns];
  B1 = new double[N_Unknowns];
  
  printf("DIRK\n");
  printf("    case %d: // \n", (TDatabase::TimeDB->RB_TYPE));
  printf("      *s = %d;\n", RB_s);
  printf("      *p = %d;\n", rb_order);
  for (i=0;i<RB_s;i++)
  {
    for (j=0;j<RB_s;j++)
    {
      printf("      A[%d][%d] = %30.20e;\n", i, j, RB_A[i*10+j]);
    }
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      alpha[%d] = %30.20e;\n", i, RB_alpha[i]);
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      m[%d] = %30.20e;\n", i, RB_m[i]);
  }
  for (i=0;i<RB_s;i++)
  {
    printf("      ms[%d] = %30.20e;\n", i, RB_ms[i]);
  }
  
  printf("    break;\n");
}
