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
   
static double NF_N_T_P5_2D_Xi[45] = 
        { 
          0.0337652428984239860938492227505,
          0.169395306766867743169300202487,
          0.380690406958401545684749139159,
          0.619309593041598454315250860839,
          0.830604693233132256830699797512,
          0.966234757101576013906150777249,

          0.966234757101576013906150777249,
          0.830604693233132256830699797512,
          0.619309593041598454315250860839,
          0.380690406958401545684749139159,
          0.169395306766867743169300202487,
          0.0337652428984239860938492227505,

          0, 0, 0, 0, 0, 0,

          0.3236494811127589E-1, 0.3236494811127589E-1, 0.9352701037774482,
          0.1193509122825813, 0.1193509122825813, 0.7612981754348374,
          0.5346110482707583, 0.5346110482707583, -0.6922209654151662E-1,
          0.2033099004312825, 0.2033099004312825, 0.5933801991374351,
          0.3989693029658552, 0.3989693029658552, 0.2020613940682896,
          0.5932012134282128, 0.5932012134282128, 0.5017813831049467E-1,
          0.5017813831049467E-1, 0.3566206482612926, 0.3566206482612926,
          0.8074890031597922, 0.8074890031597922, 0.210220165361663E-1,
          0.210220165361663E-1, 0.1714889803040415, 0.1714889803040415
        };

static double NF_N_T_P5_2D_Eta[45] = 
        { 
          0, 0, 0, 0, 0, 0,

          0.0337652428984239860938492227505,
          0.169395306766867743169300202487,
          0.380690406958401545684749139159,
          0.619309593041598454315250860839,
          0.830604693233132256830699797512,
          0.966234757101576013906150777249,

          0.966234757101576013906150777249,
          0.830604693233132256830699797512,
          0.619309593041598454315250860839,
          0.380690406958401545684749139159,
          0.169395306766867743169300202487,
          0.0337652428984239860938492227505,

          0.3236494811127589E-1, 0.9352701037774482, 0.3236494811127589E-1,
          0.1193509122825813, 0.7612981754348374, 0.1193509122825813,
          0.5346110482707583, -0.6922209654151662E-1, 0.5346110482707583,
          0.2033099004312825, 0.5933801991374351, 0.2033099004312825,
          0.3989693029658552, 0.2020613940682896, 0.3989693029658552,
          0.5017813831049467E-1, 0.3566206482612926, 0.5932012134282128,
          0.3566206482612926, 0.5932012134282128, 0.5017813831049467E-1,
          0.210220165361663E-1, 0.1714889803040415, 0.8074890031597922,
          0.1714889803040415, 0.8074890031597922, 0.210220165361663E-1
        };

static double NF_N_T_P5_2D_T[6] = 
        { 
          -0.932469514203152027812301554495,
          -0.661209386466264513661399595021,
          -0.238619186083196908630501721681, 
           0.238619186083196908630501721681,
           0.661209386466264513661399595021, 
           0.932469514203152027812301554495
        };

static double NF_N_T_P5_2D_EdgeWeight0[6] = {
        0.0856622461895851725201480710875,
        0.180380786524069303784916756918,
        0.233956967286345523694935171995,
        0.233956967286345523694935171995,
        0.180380786524069303784916756918,
        0.0856622461895851725201480710875 };

static double NF_N_T_P5_2D_EdgeWeight1[6] = {
       -0.239632299269859890105518370648,
       -0.357808407563646294838734602300,      
       -0.167479863337082682628861535781,
        0.167479863337082682628861535781,
        0.357808407563646294838734602300,
        0.239632299269859890105518370648 };

static double NF_N_T_P5_2D_EdgeWeight2[6] = {
        0.344468918744913553951585996483,
        0.140513727783900954708003276278,
       -0.484982646528814508659589272762,
       -0.484982646528814508659589272762,
        0.140513727783900954708003276278,
        0.344468918744913553951585996483 };

static double NF_N_T_P5_2D_EdgeWeight3[6] = {
       -0.376721097993911994360121892502,
        0.339803199487927794882556709270,
        0.530551948742438674219510675340,
       -0.530551948742438674219510675340,
       -0.339803199487927794882556709270,
        0.376721097993911994360121892502 };

static double NF_N_T_P5_2D_EdgeWeight4[6] = {
        0.32534907297640428362836948428805172,
       -0.69522592887681074559282921930495380,
        0.36987685590040646196445973503819405,
        0.36987685590040646196445973503819405,  
       -0.69522592887681074559282921930495380,
        0.32534907297640428362836948428805172 };

static double NF_N_T_P5_2D_CellWeight0[27] = {
        0.1117514851689078471460708,
        0.4630915880614403340099810e-5,
        0.4630915880614403340099810e-5,
        0.1596568030840561983457594,
        0.0006151772848242272951167625,
        0.0006151772848242272951167625,
       -0.3074800516466622838992547e-5,
        0.001416434473412530323717482,
        0.001416434473412530323717482,
        0.1239430245205130961597144,
        0.004985382936139870119727217,
        0.004985382936139870119727217,
        0.006364779273706852291409020,
        0.04899499958107480071432628,
        0.04899499958107480071432628,
        0.02373718471524705173174519,
        0.00006612312750449697159795492,
        0.02373718471524705173174519,
        0.1092485914182767059820645,
        0.00006612312750449697159795492,
        0.1092485914182767059820645,
        0.001044334612342460313333734,
        0.1923774386767252901075908e-5,
        0.001044334612342460313333734,
        0.1090287085375766676977504,
        0.1923774386767252901075908e-5,
        0.1090287085375766676977504 };

static double NF_N_T_P5_2D_CellWeight1[27] = {
        0.4630915880614394755041929e-5,
        0.4630915880614394755041929e-5,
        0.1117514851689078399769235,
        0.0006151772848242272951167625,
        0.0006151772848242272951167625,
        0.1596568030840561983457594,
        0.001416434473412530164749446,
        0.001416434473412530164749446,
       -0.3074800516466625504153532e-5,
        0.004985382936139877476058009,
        0.004985382936139877476058009,
        0.1239430245205131588225862,
        0.04899499958107480071432628,
        0.04899499958107480071432628,
        0.006364779273706852291409020,
        0.1092485914182767446573128,
        0.1092485914182767446573128,
        0.00006612312750449724832916185,
        0.00006612312750449724832916185,
        0.02373718471524706570964808,
        0.02373718471524706570964808,
        0.1090287085375766676977504,
        0.1090287085375766676977504,
        0.1923774386767252901075908e-5,
        0.1923774386767252901075908e-5,
        0.001044334612342460313333734,
        0.001044334612342460313333734 };

static double NF_N_T_P5_2D_CellWeight2[27] = {
        0.4630915880614394755041929e-5,
        0.1117514851689078399769235,
        0.4630915880614394755041929e-5,
        0.0006151772848242272951167625,
        0.1596568030840561983457594,
        0.0006151772848242272951167625,
        0.001416434473412530164749446,
       -0.3074800516466625504153532e-5,
        0.001416434473412530164749446,
        0.004985382936139877476058009,
        0.1239430245205131588225862,
        0.004985382936139877476058009,
        0.04899499958107480071432628,
        0.006364779273706852291409020,
        0.04899499958107480071432628,
        0.00006612312750449724832916185,
        0.02373718471524706570964808,
        0.1092485914182767446573128,
        0.02373718471524706570964808,
        0.1092485914182767446573128,
        0.00006612312750449724832916185,
        0.1923774386767252901075908e-5,
        0.001044334612342460313333734,
        0.1090287085375766676977504,
        0.001044334612342460313333734,
        0.1090287085375766676977504,
        0.1923774386767252901075908e-5 };


static double NF_N_T_P5_2D_CellWeight3[27] = {
        0.01160145396792356276700694,
        0.00001389274764184320143524155,
        0.0004014673987446835280592390,
        0.07508957350115266701528283,
        0.001845531854472681885350287,
        0.01177200916730636231568522,
        0.00007124122539451513582433278,
        0.004249303420237590812184409,
       -0.0005502050369915114418360994,
        0.1273998222946151501560467,
        0.01495614880841961771551244,
        0.04365101030221956034589423,
        0.03770168312563656790901463,
        0.1469849987432244021429788,
        0.07444180174363820154398787,
        0.1184529850846203691495780,
        0.002345104110596251010004287,
        0.01001978777911980044258594,
        0.02772359937734555959380499,
        0.001409829462970686862074220,
        0.1970341729481329266677232,
        0.01475235400415651230017852,
        0.0002216856778514819455227546,
        0.0003840600040486093626694198,
        0.008515298554500189756935485,
        0.00004708008492251608533997387,
        0.06946430964809912667775407 };

static double NF_N_T_P5_2D_CellWeight4[27] = {
        0.01160145396792356276700694,
        0.0004014673987446835280592390,
        0.00001389274764184320143524155,
        0.07508957350115266701528283,
        0.01177200916730636231568522,
        0.001845531854472681885350287,
        0.00007124122539451513582433278,
       -0.0005502050369915114418360994,
        0.004249303420237590812184409,
        0.1273998222946151501560467,
        0.04365101030221956034589423,
        0.01495614880841961771551244,
        0.03770168312563656790901463,
        0.07444180174363820154398787,
        0.1469849987432244021429788,
        0.01001978777911980044258594,
        0.001409829462970686862074220,
        0.1184529850846203691495780,
        0.1970341729481329266677232,
        0.002345104110596251010004287,
        0.02772359937734555959380499,
        0.0003840600040486093626694198,
        0.00004708008492251608533997387,
        0.01475235400415651230017852,
        0.06946430964809912667775407,
        0.0002216856778514819455227546,
        0.008515298554500189756935485 };

static double NF_N_T_P5_2D_CellWeight5[27] = {
        0.0004014673987446830404687462,
        0.00001389274764184319285018367,
        0.01160145396792356943997865,
        0.01177200916730636231568522,
        0.001845531854472681885350287,
        0.07508957350115266701528283,
       -0.0005502050369915112417013094,
        0.004249303420237590653216373,
        0.00007124122539451517965624771,
        0.04365101030221959592993202,
        0.01495614880841962507184323,
        0.1273998222946151304335435,
        0.07444180174363820154398787,
        0.1469849987432244021429788,
        0.03770168312563656790901463,
        0.1970341729481329344940390,
        0.02772359937734552746153713,
        0.001409829462970690518851339,
        0.002345104110596257276253480,
        0.01001978777911979039819138,
        0.1184529850846204016732391,
        0.06946430964809912667775407,
        0.008515298554500189756935485,
        0.00004708008492251608533997387,
        0.0002216856778514819455227546,
        0.0003840600040486093626694198,
        0.01475235400415651230017852 };

static double NF_N_T_P5_2D_CellWeight6[27] = {
        0.00001389274764184318426512579,
        0.0004014673987446830318836883,
        0.01160145396792356227083139,
        0.001845531854472681885350287,
        0.01177200916730636231568522,
        0.07508957350115266701528283,
        0.004249303420237590494248337,
       -0.0005502050369915114006693454,
        0.00007124122539451517699108672,
        0.01495614880841963242817403,
        0.04365101030221960328626281,
        0.1273998222946151930964153,
        0.1469849987432244021429788,
        0.07444180174363820154398787,
        0.03770168312563656790901463,
        0.02772359937734556613678539,
        0.1970341729481329731692872,
        0.002345104110596257552984687,
        0.001409829462970690795582546,
        0.1184529850846204156511420,
        0.01001978777911980437609427,
        0.008515298554500189756935485,
        0.06946430964809912667775407,
        0.0002216856778514819455227546,
        0.00004708008492251608533997387,
        0.01475235400415651230017852,
        0.0003840600040486093626694198 };

static double NF_N_T_P5_2D_CellWeight7[27] = {
        0.0004014673987446830404687462,
        0.01160145396792356943997865,
        0.00001389274764184319285018367,
        0.01177200916730636231568522,
        0.07508957350115266701528283,
        0.001845531854472681885350287,
       -0.0005502050369915112417013094,
        0.00007124122539451517965624771,
        0.004249303420237590653216373,
        0.04365101030221959592993202,
        0.1273998222946151304335435,
        0.01495614880841962507184323,
        0.07444180174363820154398787,
        0.03770168312563656790901463,
        0.1469849987432244021429788,
        0.001409829462970690518851339,
        0.01001978777911979039819138,
        0.1970341729481329344940390,
        0.1184529850846204016732391,
        0.02772359937734552746153713,
        0.002345104110596257276253480,
        0.00004708008492251608533997387,
        0.0003840600040486093626694198,
        0.06946430964809912667775407,
        0.01475235400415651230017852,
        0.008515298554500189756935485,
        0.0002216856778514819455227546 };

static double NF_N_T_P5_2D_CellWeight8[27] = {
        0.00001389274764184318426512579,
        0.01160145396792356227083139,
        0.0004014673987446830318836883,
        0.001845531854472681885350287,
        0.07508957350115266701528283,
        0.01177200916730636231568522,
        0.004249303420237590494248337,
        0.00007124122539451517699108672,
       -0.0005502050369915114006693454,
        0.01495614880841963242817403,
        0.1273998222946151930964153,
        0.04365101030221960328626281,
        0.1469849987432244021429788,
        0.03770168312563656790901463,
        0.07444180174363820154398787,
        0.002345104110596257552984687,
        0.1184529850846204156511420,
        0.02772359937734556613678539,
        0.01001978777911980437609427,
        0.1970341729481329731692872,
        0.001409829462970690795582546,
        0.0002216856778514819455227546,
        0.01475235400415651230017852,
        0.008515298554500189756935485,
        0.0003840600040486093626694198,
        0.06946430964809912667775407,
        0.00004708008492251608533997387 };

static double NF_N_T_P5_2D_CellWeight9[27] = {
        0.0008029347974893660809374924,
        0.0008029347974893665599429273,
        0.0008029347974893665599429273,
        0.02354401833461272463137044,
        0.02354401833461272463137044,
        0.02354401833461272463137044,
       -0.001100410073983022483402619,
       -0.001100410073983022842505445,
       -0.001100410073983022842505445,
        0.08730202060443919185986404,
        0.08730202060443916363215705,
        0.08730202060443916363215705,
        0.1488836034872764030879757,
        0.1488836034872764030879757,
        0.1488836034872764030879757,
        0.03333374159822688806778858,
        0.03333374159822684810920494,
        0.03333374159822688806778858,
        0.03333374159822689067726065,
        0.03333374159822684810920494,
        0.03333374159822689067726065,
        0.003616841493522468607614618,
        0.003616841493522468607614618,
        0.003616841493522468607614618,
        0.003616841493522468607614618,
        0.003616841493522468607614618,
        0.003616841493522468607614618 };

void NF_N_T_P5_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] =( NF_N_T_P5_2D_EdgeWeight0[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight0[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight0[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight0[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight0[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight0[5]*PointValues[5]);
  Functionals[1] =( NF_N_T_P5_2D_EdgeWeight0[0]*PointValues[6]
                   +NF_N_T_P5_2D_EdgeWeight0[1]*PointValues[7]
                   +NF_N_T_P5_2D_EdgeWeight0[2]*PointValues[8]
                   +NF_N_T_P5_2D_EdgeWeight0[3]*PointValues[9]
                   +NF_N_T_P5_2D_EdgeWeight0[4]*PointValues[10]
                   +NF_N_T_P5_2D_EdgeWeight0[5]*PointValues[11]);
  Functionals[2] =( NF_N_T_P5_2D_EdgeWeight0[0]*PointValues[12]
                   +NF_N_T_P5_2D_EdgeWeight0[1]*PointValues[13]
                   +NF_N_T_P5_2D_EdgeWeight0[2]*PointValues[14]
                   +NF_N_T_P5_2D_EdgeWeight0[3]*PointValues[15]
                   +NF_N_T_P5_2D_EdgeWeight0[4]*PointValues[16]
                   +NF_N_T_P5_2D_EdgeWeight0[5]*PointValues[17]);

  Functionals[3] =( NF_N_T_P5_2D_EdgeWeight1[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight1[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight1[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight1[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight1[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight1[5]*PointValues[5]);
  Functionals[4] =( NF_N_T_P5_2D_EdgeWeight1[0]*PointValues[6]
                   +NF_N_T_P5_2D_EdgeWeight1[1]*PointValues[7]
                   +NF_N_T_P5_2D_EdgeWeight1[2]*PointValues[8]
                   +NF_N_T_P5_2D_EdgeWeight1[3]*PointValues[9]
                   +NF_N_T_P5_2D_EdgeWeight1[4]*PointValues[10]
                   +NF_N_T_P5_2D_EdgeWeight1[5]*PointValues[11]);
  Functionals[5] =( NF_N_T_P5_2D_EdgeWeight1[0]*PointValues[12]
                   +NF_N_T_P5_2D_EdgeWeight1[1]*PointValues[13]
                   +NF_N_T_P5_2D_EdgeWeight1[2]*PointValues[14]
                   +NF_N_T_P5_2D_EdgeWeight1[3]*PointValues[15]
                   +NF_N_T_P5_2D_EdgeWeight1[4]*PointValues[16]
                   +NF_N_T_P5_2D_EdgeWeight1[5]*PointValues[17]);

  Functionals[6] =( NF_N_T_P5_2D_EdgeWeight2[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight2[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight2[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight2[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight2[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight2[5]*PointValues[5]);
  Functionals[7] =( NF_N_T_P5_2D_EdgeWeight2[0]*PointValues[6]
                   +NF_N_T_P5_2D_EdgeWeight2[1]*PointValues[7]
                   +NF_N_T_P5_2D_EdgeWeight2[2]*PointValues[8]
                   +NF_N_T_P5_2D_EdgeWeight2[3]*PointValues[9]
                   +NF_N_T_P5_2D_EdgeWeight2[4]*PointValues[10]
                   +NF_N_T_P5_2D_EdgeWeight2[5]*PointValues[11]);
  Functionals[8] =( NF_N_T_P5_2D_EdgeWeight2[0]*PointValues[12]
                   +NF_N_T_P5_2D_EdgeWeight2[1]*PointValues[13]
                   +NF_N_T_P5_2D_EdgeWeight2[2]*PointValues[14]
                   +NF_N_T_P5_2D_EdgeWeight2[3]*PointValues[15]
                   +NF_N_T_P5_2D_EdgeWeight2[4]*PointValues[16]
                   +NF_N_T_P5_2D_EdgeWeight2[5]*PointValues[17]);

  Functionals[9] =( NF_N_T_P5_2D_EdgeWeight3[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight3[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight3[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight3[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight3[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight3[5]*PointValues[5]);
  Functionals[10]=( NF_N_T_P5_2D_EdgeWeight3[0]*PointValues[6]
                   +NF_N_T_P5_2D_EdgeWeight3[1]*PointValues[7]
                   +NF_N_T_P5_2D_EdgeWeight3[2]*PointValues[8]
                   +NF_N_T_P5_2D_EdgeWeight3[3]*PointValues[9]
                   +NF_N_T_P5_2D_EdgeWeight3[4]*PointValues[10]
                   +NF_N_T_P5_2D_EdgeWeight3[5]*PointValues[11]);
  Functionals[11]=( NF_N_T_P5_2D_EdgeWeight3[0]*PointValues[12]
                   +NF_N_T_P5_2D_EdgeWeight3[1]*PointValues[13]
                   +NF_N_T_P5_2D_EdgeWeight3[2]*PointValues[14]
                   +NF_N_T_P5_2D_EdgeWeight3[3]*PointValues[15]
                   +NF_N_T_P5_2D_EdgeWeight3[4]*PointValues[16]
                   +NF_N_T_P5_2D_EdgeWeight3[5]*PointValues[17]);

  Functionals[12]=( NF_N_T_P5_2D_EdgeWeight4[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight4[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight4[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight4[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight4[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight4[5]*PointValues[5]);
  Functionals[13]=( NF_N_T_P5_2D_EdgeWeight4[0]*PointValues[6]
                   +NF_N_T_P5_2D_EdgeWeight4[1]*PointValues[7]
                   +NF_N_T_P5_2D_EdgeWeight4[2]*PointValues[8]
                   +NF_N_T_P5_2D_EdgeWeight4[3]*PointValues[9]
                   +NF_N_T_P5_2D_EdgeWeight4[4]*PointValues[10]
                   +NF_N_T_P5_2D_EdgeWeight4[5]*PointValues[11]);
  Functionals[14]=( NF_N_T_P5_2D_EdgeWeight4[0]*PointValues[12]
                   +NF_N_T_P5_2D_EdgeWeight4[1]*PointValues[13]
                   +NF_N_T_P5_2D_EdgeWeight4[2]*PointValues[14]
                   +NF_N_T_P5_2D_EdgeWeight4[3]*PointValues[15]
                   +NF_N_T_P5_2D_EdgeWeight4[4]*PointValues[16]
                   +NF_N_T_P5_2D_EdgeWeight4[5]*PointValues[17]);


  Functionals[15] =( NF_N_T_P5_2D_CellWeight0[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight0[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight0[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight0[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight0[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight0[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight0[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight0[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight0[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight0[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight0[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight0[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight0[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight0[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight0[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight0[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight0[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight0[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight0[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight0[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight0[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight0[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight0[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight0[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight0[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight0[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight0[26]*PointValues[44] );
  Functionals[16] =( NF_N_T_P5_2D_CellWeight1[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight1[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight1[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight1[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight1[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight1[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight1[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight1[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight1[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight1[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight1[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight1[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight1[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight1[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight1[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight1[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight1[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight1[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight1[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight1[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight1[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight1[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight1[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight1[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight1[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight1[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight1[26]*PointValues[44] );
  Functionals[17] =( NF_N_T_P5_2D_CellWeight2[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight2[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight2[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight2[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight2[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight2[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight2[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight2[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight2[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight2[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight2[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight2[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight2[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight2[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight2[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight2[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight2[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight2[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight2[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight2[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight2[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight2[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight2[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight2[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight2[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight2[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight2[26]*PointValues[44] );
  Functionals[18] =( NF_N_T_P5_2D_CellWeight3[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight3[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight3[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight3[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight3[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight3[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight3[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight3[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight3[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight3[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight3[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight3[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight3[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight3[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight3[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight3[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight3[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight3[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight3[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight3[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight3[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight3[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight3[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight3[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight3[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight3[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight3[26]*PointValues[44] );
  Functionals[19] =( NF_N_T_P5_2D_CellWeight4[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight4[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight4[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight4[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight4[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight4[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight4[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight4[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight4[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight4[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight4[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight4[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight4[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight4[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight4[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight4[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight4[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight4[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight4[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight4[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight4[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight4[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight4[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight4[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight4[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight4[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight4[26]*PointValues[44] );
  Functionals[20] =( NF_N_T_P5_2D_CellWeight5[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight5[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight5[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight5[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight5[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight5[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight5[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight5[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight5[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight5[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight5[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight5[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight5[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight5[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight5[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight5[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight5[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight5[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight5[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight5[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight5[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight5[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight5[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight5[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight5[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight5[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight5[26]*PointValues[44] );
  Functionals[21] =( NF_N_T_P5_2D_CellWeight6[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight6[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight6[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight6[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight6[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight6[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight6[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight6[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight6[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight6[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight6[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight6[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight6[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight6[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight6[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight6[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight6[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight6[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight6[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight6[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight6[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight6[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight6[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight6[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight6[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight6[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight6[26]*PointValues[44] );
  Functionals[22] =( NF_N_T_P5_2D_CellWeight7[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight7[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight7[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight7[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight7[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight7[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight7[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight7[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight7[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight7[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight7[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight7[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight7[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight7[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight7[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight7[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight7[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight7[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight7[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight7[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight7[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight7[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight7[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight7[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight7[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight7[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight7[26]*PointValues[44] );
  Functionals[23] =( NF_N_T_P5_2D_CellWeight8[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight8[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight8[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight8[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight8[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight8[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight8[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight8[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight8[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight8[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight8[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight8[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight8[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight8[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight8[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight8[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight8[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight8[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight8[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight8[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight8[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight8[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight8[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight8[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight8[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight8[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight8[26]*PointValues[44] );
  Functionals[24] =( NF_N_T_P5_2D_CellWeight9[ 0]*PointValues[18]
                    +NF_N_T_P5_2D_CellWeight9[ 1]*PointValues[19]
                    +NF_N_T_P5_2D_CellWeight9[ 2]*PointValues[20]
                    +NF_N_T_P5_2D_CellWeight9[ 3]*PointValues[21]
                    +NF_N_T_P5_2D_CellWeight9[ 4]*PointValues[22]
                    +NF_N_T_P5_2D_CellWeight9[ 5]*PointValues[23]
                    +NF_N_T_P5_2D_CellWeight9[ 6]*PointValues[24]
                    +NF_N_T_P5_2D_CellWeight9[ 7]*PointValues[25]
                    +NF_N_T_P5_2D_CellWeight9[ 8]*PointValues[26]
                    +NF_N_T_P5_2D_CellWeight9[ 9]*PointValues[27]
                    +NF_N_T_P5_2D_CellWeight9[10]*PointValues[28]
                    +NF_N_T_P5_2D_CellWeight9[11]*PointValues[29]
                    +NF_N_T_P5_2D_CellWeight9[12]*PointValues[30]
                    +NF_N_T_P5_2D_CellWeight9[13]*PointValues[31]
                    +NF_N_T_P5_2D_CellWeight9[14]*PointValues[32]
                    +NF_N_T_P5_2D_CellWeight9[15]*PointValues[33]
                    +NF_N_T_P5_2D_CellWeight9[16]*PointValues[34]
                    +NF_N_T_P5_2D_CellWeight9[17]*PointValues[35]
                    +NF_N_T_P5_2D_CellWeight9[18]*PointValues[36]
                    +NF_N_T_P5_2D_CellWeight9[19]*PointValues[37]
                    +NF_N_T_P5_2D_CellWeight9[20]*PointValues[38]
                    +NF_N_T_P5_2D_CellWeight9[21]*PointValues[39]
                    +NF_N_T_P5_2D_CellWeight9[22]*PointValues[40]
                    +NF_N_T_P5_2D_CellWeight9[23]*PointValues[41]
                    +NF_N_T_P5_2D_CellWeight9[24]*PointValues[42]
                    +NF_N_T_P5_2D_CellWeight9[25]*PointValues[43]
                    +NF_N_T_P5_2D_CellWeight9[26]*PointValues[44] );

  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
    {
      Functionals[3] = -Functionals[3];
      Functionals[9] = -Functionals[9];
    }
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
    {
      Functionals[4] = -Functionals[4];
      Functionals[10] = -Functionals[10];
    }
    if(Cell->GetVertex(2) > Cell->GetVertex(0))
    {
      Functionals[5] = -Functionals[5];
      Functionals[11] = -Functionals[11];
    }
  }
  */

  if(Cell)
  {
    OwnNum = Coll->GetIndex(Cell);

    neigh = Cell->GetJoint(0)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[3] = -Functionals[3];
        Functionals[9] = -Functionals[9];
      }
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[ 4] = -Functionals[ 4];
        Functionals[10] = -Functionals[10];
      }
    } // endif neigh

    neigh = Cell->GetJoint(2)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[ 5] = -Functionals[ 5];
        Functionals[11] = -Functionals[11];
      }
    } // endif neigh
  } // endif Cell
}

void NF_N_T_P5_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                           double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] =( NF_N_T_P5_2D_EdgeWeight0[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight0[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight0[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight0[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight0[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight0[5]*PointValues[5]);

  Functionals[1] =( NF_N_T_P5_2D_EdgeWeight1[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight1[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight1[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight1[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight1[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight1[5]*PointValues[5]);

  Functionals[2] =( NF_N_T_P5_2D_EdgeWeight2[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight2[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight2[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight2[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight2[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight2[5]*PointValues[5]);

  Functionals[3] =( NF_N_T_P5_2D_EdgeWeight3[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight3[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight3[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight3[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight3[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight3[5]*PointValues[5]);

  Functionals[4] =( NF_N_T_P5_2D_EdgeWeight4[0]*PointValues[0]
                   +NF_N_T_P5_2D_EdgeWeight4[1]*PointValues[1]
                   +NF_N_T_P5_2D_EdgeWeight4[2]*PointValues[2]
                   +NF_N_T_P5_2D_EdgeWeight4[3]*PointValues[3]
                   +NF_N_T_P5_2D_EdgeWeight4[4]*PointValues[4]
                   +NF_N_T_P5_2D_EdgeWeight4[5]*PointValues[5]);

  if(Joint != -1)
  {
    // if(Cell->GetVertex(Joint) > Cell->GetVertex((Joint+1)%3))
    // {
    //   Functionals[1] = -Functionals[1];
    //   Functionals[3] = -Functionals[3];
    // }
    neigh = Cell->GetJoint(Joint)->GetNeighbour(Cell);
    if(neigh)
    {
      OwnNum = Coll->GetIndex(Cell);
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[1] = -Functionals[1];
        Functionals[3] = -Functionals[3];
      }
    } // endif neigh
    // */
  }
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_N_T_P5_2D_Obj = new TNodalFunctional2D
        (NF_N_T_P5_2D, 25, 5, 45, 6, NF_N_T_P5_2D_Xi, NF_N_T_P5_2D_Eta,
         NF_N_T_P5_2D_T, NF_N_T_P5_2D_EvalAll, NF_N_T_P5_2D_EvalEdge);
