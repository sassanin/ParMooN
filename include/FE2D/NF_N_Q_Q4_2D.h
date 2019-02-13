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
   

static double NF_N_Q_Q4_2D_Xi[45] = 
        {
          -0.906179845938663992797626878,
          -0.538469310105683091036314421,
           0,
           0.538469310105683091036314421,
           0.906179845938663992797626878,
           1, 1, 1, 1, 1,
           0.906179845938663992797626878,
           0.538469310105683091036314421,
           0,
          -0.538469310105683091036314421,
          -0.906179845938663992797626878,
          -1, -1, -1, -1, -1,

          -0.906179845938663992797626878,
          -0.538469310105683091036314421,
           0,
           0.538469310105683091036314421,
           0.906179845938663992797626878,

          -0.906179845938663992797626878,
          -0.538469310105683091036314421,
           0,
           0.538469310105683091036314421,
           0.906179845938663992797626878,

          -0.906179845938663992797626878,
          -0.538469310105683091036314421,
           0,
           0.538469310105683091036314421,
           0.906179845938663992797626878,

          -0.906179845938663992797626878,
          -0.538469310105683091036314421,
           0,
           0.538469310105683091036314421,
           0.906179845938663992797626878,

          -0.906179845938663992797626878,
          -0.538469310105683091036314421,
           0,
           0.538469310105683091036314421,
           0.906179845938663992797626878
        };

static double NF_N_Q_Q4_2D_Eta[45] = 
        { -1, -1, -1, -1, -1,
          -0.906179845938663992797626878,
          -0.538469310105683091036314421,
           0,
           0.538469310105683091036314421,
           0.906179845938663992797626878,
           1, 1, 1, 1, 1,
           0.906179845938663992797626878,
           0.538469310105683091036314421,
           0,
          -0.538469310105683091036314421,
          -0.906179845938663992797626878,

          -0.906179845938663992797626878,
          -0.906179845938663992797626878,
          -0.906179845938663992797626878,
          -0.906179845938663992797626878,
          -0.906179845938663992797626878,

          -0.538469310105683091036314421,
          -0.538469310105683091036314421,
          -0.538469310105683091036314421,
          -0.538469310105683091036314421,
          -0.538469310105683091036314421,

           0, 0, 0, 0, 0,

           0.538469310105683091036314421,
           0.538469310105683091036314421,
           0.538469310105683091036314421,
           0.538469310105683091036314421,
           0.538469310105683091036314421,

           0.906179845938663992797626878,
           0.906179845938663992797626878,
           0.906179845938663992797626878,
           0.906179845938663992797626878,
           0.906179845938663992797626878
        };

static double NF_N_Q_Q4_2D_T[5] = 
        {
               -0.906179845938663992797626878,
               -0.538469310105683091036314421,
                0,
                0.538469310105683091036314421,
                0.906179845938663992797626878
        };

static double NF_N_Q_Q4_2D_EdgeWeight0[5] = {
        0.1184634425280945437571320205,
        0.2393143352496832340206457575,
        0.2844444444444444444444444444,
        0.2393143352496832340206457575,
        0.1184634425280945437571320205 };

static double NF_N_Q_Q4_2D_EdgeWeight1[5] = {
       -0.3220475522984174693134711033,
       -0.3865902750008912622604589675,
        0.0,
        0.3865902750008912622604589675,
        0.3220475522984174693134711033 };

static double NF_N_Q_Q4_2D_EdgeWeight2[5] = {
        0.4334238969965230841044414528,
       -0.0778683414409675285488858967,
       -0.7111111111111111111111111111,
       -0.0778683414409675285488858967,
        0.4334238969965230841044414528 };

static double NF_N_Q_Q4_2D_EdgeWeight3[5] = {
       -0.4154771413508325868510801940,
        0.6991986448892333762714923941,
        0.0,
       -0.6991986448892333762714923941,
        0.4154771413508325868510801940 };

static double NF_N_Q_Q4_2D_CellWeight0[25] = {
        .0140335872156071589886627895397, .02835,
        .0336962680968802257798064413932, .02835,
        .0140335872156071589886627895397, .02835,
        .0572713510559977792829421488044, .0680716331376876754547614599244,
        .0572713510559977792829421488044, .02835,
        .0336962680968802257798064413932, .0680716331376876754547614599244,
        .0809086419753086419753086420069, .0680716331376876754547614599244,
        .0336962680968802257798064413932, .02835,
        .0572713510559977792829421488044, .0680716331376876754547614599244,
        .0572713510559977792829421488044, .02835,
        .0140335872156071589886627895397, .02835,
        .0336962680968802257798064413932, .02835,
        .0140335872156071589886627895397 };

static double NF_N_Q_Q4_2D_CellWeight1[25] = {
       -.0381508617030170997669370509012, -.0457968148244883468926385415758, 0.,
        .0457968148244883468926385415758, .0381508617030170997669370509012,
        -.0770705958970833725874381660912, -.0925165946758305269383256669817,
        0., .0925165946758305269383256669817, .0770705958970833725874381660912,
        -.0916046370982165246047206693712, -.109963456000253514598530550856,
        0., .109963456000253514598530550856, .0916046370982165246047206693712,
        -.0770705958970833725874381660912, -.0925165946758305269383256669817,
        0., .0925165946758305269383256669817, .0770705958970833725874381660912,
        -.0381508617030170997669370509012, -.0457968148244883468926385415758,
        0., .0457968148244883468926385415758, .0381508617030170997669370509012 };

static double NF_N_Q_Q4_2D_CellWeight2[25] = {
        -.0381508617030170997669370509012, -.0770705958970833725874381660912,
        -.0916046370982165246047206693712, -.0770705958970833725874381660912,
        -.0381508617030170997669370509012, -.0457968148244883468926385415758,
        -.0925165946758305269383256669817, -.109963456000253514598530550856,
        -.0925165946758305269383256669817, -.0457968148244883468926385415758,
        0., 0., 0., 0., 0., .0457968148244883468926385415758,
        .0925165946758305269383256669817, .109963456000253514598530550856,
        .0925165946758305269383256669817, .0457968148244883468926385415758,
        .0381508617030170997669370509012, .0770705958970833725874381660912,
        .0916046370982165246047206693712, .0770705958970833725874381660912,
        .0381508617030170997669370509012 };

static double NF_N_Q_Q4_2D_CellWeight3[25] = {
        .0513448869121503817132846634880, -.00922455179105009948852661168580,
        -.0842406702422005644495161034830, -.00922455179105009948852661168580,
        .0513448869121503817132846634880, .103724551791050099488526611782,
        -.0186350103689405051700747867548, -.170179082844219188636903649811,
        -.0186350103689405051700747867548, .103724551791050099488526611782,
        .123285019590122121700818902026, -.0221492171209863192316830993719,
        -.202271604938271604938271605017, -.0221492171209863192316830993719,
        .123285019590122121700818902026, .103724551791050099488526611782,
        -.0186350103689405051700747867548, -.170179082844219188636903649811,
        -.0186350103689405051700747867548, .103724551791050099488526611782,
        .0513448869121503817132846634880, -.00922455179105009948852661168580,
        -.0842406702422005644495161034830, -.00922455179105009948852661168580,
        .0513448869121503817132846634880 };

static double NF_N_Q_Q4_2D_CellWeight4[25] = {
        .103714625941401935021929964806, .124500451806409120356135620186, -0.,
        -.124500451806409120356135620186, -.103714625941401935021929964806,
        .124500451806409120356135620186, .149452040725264731644736702307, -0.,
        -.149452040725264731644736702307, -.124500451806409120356135620186,
        -0., -0., 0., 0., 0., -.124500451806409120356135620186,
        -.149452040725264731644736702307, 0., .149452040725264731644736702307,
        .124500451806409120356135620186, -.103714625941401935021929964806,
        -.124500451806409120356135620186, 0., .124500451806409120356135620186,
        .103714625941401935021929964806 };

static double NF_N_Q_Q4_2D_CellWeight5[25] = {
        .0513448869121503817132846634880, .103724551791050099488526611782,
        .123285019590122121700818902026, .103724551791050099488526611782,
        .0513448869121503817132846634880, -.00922455179105009948852661168580,
        -.0186350103689405051700747867548, -.0221492171209863192316830993719,
        -.0186350103689405051700747867548, -.00922455179105009948852661168580,
        -.0842406702422005644495161034830, -.170179082844219188636903649811,
        -.202271604938271604938271605017, -.170179082844219188636903649811,
        -.0842406702422005644495161034830, -.00922455179105009948852661168580,
        -.0186350103689405051700747867548, -.0221492171209863192316830993719,
        -.0186350103689405051700747867548, -.00922455179105009948852661168580,
        .0513448869121503817132846634880, .103724551791050099488526611782,
        .123285019590122121700818902026, .103724551791050099488526611782,
        .0513448869121503817132846634880 };

void NF_N_Q_Q4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] =( NF_N_Q_Q4_2D_EdgeWeight0[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight0[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight0[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight0[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight0[4]*PointValues[4]);
  Functionals[1] =( NF_N_Q_Q4_2D_EdgeWeight0[0]*PointValues[5]
                   +NF_N_Q_Q4_2D_EdgeWeight0[1]*PointValues[6]
                   +NF_N_Q_Q4_2D_EdgeWeight0[2]*PointValues[7]
                   +NF_N_Q_Q4_2D_EdgeWeight0[3]*PointValues[8]
                   +NF_N_Q_Q4_2D_EdgeWeight0[4]*PointValues[9]);
  Functionals[2] =( NF_N_Q_Q4_2D_EdgeWeight0[0]*PointValues[10]
                   +NF_N_Q_Q4_2D_EdgeWeight0[1]*PointValues[11]
                   +NF_N_Q_Q4_2D_EdgeWeight0[2]*PointValues[12]
                   +NF_N_Q_Q4_2D_EdgeWeight0[3]*PointValues[13]
                   +NF_N_Q_Q4_2D_EdgeWeight0[4]*PointValues[14]);
  Functionals[3] =( NF_N_Q_Q4_2D_EdgeWeight0[0]*PointValues[15]
                   +NF_N_Q_Q4_2D_EdgeWeight0[1]*PointValues[16]
                   +NF_N_Q_Q4_2D_EdgeWeight0[2]*PointValues[17]
                   +NF_N_Q_Q4_2D_EdgeWeight0[3]*PointValues[18]
                   +NF_N_Q_Q4_2D_EdgeWeight0[4]*PointValues[19]);

  Functionals[4] =( NF_N_Q_Q4_2D_EdgeWeight1[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight1[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight1[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight1[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight1[4]*PointValues[4]);
  Functionals[5] =( NF_N_Q_Q4_2D_EdgeWeight1[0]*PointValues[5]
                   +NF_N_Q_Q4_2D_EdgeWeight1[1]*PointValues[6]
                   +NF_N_Q_Q4_2D_EdgeWeight1[2]*PointValues[7]
                   +NF_N_Q_Q4_2D_EdgeWeight1[3]*PointValues[8]
                   +NF_N_Q_Q4_2D_EdgeWeight1[4]*PointValues[9]);
  Functionals[6] =( NF_N_Q_Q4_2D_EdgeWeight1[0]*PointValues[10]
                   +NF_N_Q_Q4_2D_EdgeWeight1[1]*PointValues[11]
                   +NF_N_Q_Q4_2D_EdgeWeight1[2]*PointValues[12]
                   +NF_N_Q_Q4_2D_EdgeWeight1[3]*PointValues[13]
                   +NF_N_Q_Q4_2D_EdgeWeight1[4]*PointValues[14]);
  Functionals[7] =( NF_N_Q_Q4_2D_EdgeWeight1[0]*PointValues[15]
                   +NF_N_Q_Q4_2D_EdgeWeight1[1]*PointValues[16]
                   +NF_N_Q_Q4_2D_EdgeWeight1[2]*PointValues[17]
                   +NF_N_Q_Q4_2D_EdgeWeight1[3]*PointValues[18]
                   +NF_N_Q_Q4_2D_EdgeWeight1[4]*PointValues[19]);

  Functionals[8] =( NF_N_Q_Q4_2D_EdgeWeight2[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight2[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight2[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight2[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight2[4]*PointValues[4]);
  Functionals[9] =( NF_N_Q_Q4_2D_EdgeWeight2[0]*PointValues[5]
                   +NF_N_Q_Q4_2D_EdgeWeight2[1]*PointValues[6]
                   +NF_N_Q_Q4_2D_EdgeWeight2[2]*PointValues[7]
                   +NF_N_Q_Q4_2D_EdgeWeight2[3]*PointValues[8]
                   +NF_N_Q_Q4_2D_EdgeWeight2[4]*PointValues[9]);
  Functionals[10]=( NF_N_Q_Q4_2D_EdgeWeight2[0]*PointValues[10]
                   +NF_N_Q_Q4_2D_EdgeWeight2[1]*PointValues[11]
                   +NF_N_Q_Q4_2D_EdgeWeight2[2]*PointValues[12]
                   +NF_N_Q_Q4_2D_EdgeWeight2[3]*PointValues[13]
                   +NF_N_Q_Q4_2D_EdgeWeight2[4]*PointValues[14]);
  Functionals[11]=( NF_N_Q_Q4_2D_EdgeWeight2[0]*PointValues[15]
                   +NF_N_Q_Q4_2D_EdgeWeight2[1]*PointValues[16]
                   +NF_N_Q_Q4_2D_EdgeWeight2[2]*PointValues[17]
                   +NF_N_Q_Q4_2D_EdgeWeight2[3]*PointValues[18]
                   +NF_N_Q_Q4_2D_EdgeWeight2[4]*PointValues[19]);

  Functionals[12]=( NF_N_Q_Q4_2D_EdgeWeight3[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight3[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight3[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight3[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight3[4]*PointValues[4]);
  Functionals[13]=( NF_N_Q_Q4_2D_EdgeWeight3[0]*PointValues[5]
                   +NF_N_Q_Q4_2D_EdgeWeight3[1]*PointValues[6]
                   +NF_N_Q_Q4_2D_EdgeWeight3[2]*PointValues[7]
                   +NF_N_Q_Q4_2D_EdgeWeight3[3]*PointValues[8]
                   +NF_N_Q_Q4_2D_EdgeWeight3[4]*PointValues[9]);
  Functionals[14]=( NF_N_Q_Q4_2D_EdgeWeight3[0]*PointValues[10]
                   +NF_N_Q_Q4_2D_EdgeWeight3[1]*PointValues[11]
                   +NF_N_Q_Q4_2D_EdgeWeight3[2]*PointValues[12]
                   +NF_N_Q_Q4_2D_EdgeWeight3[3]*PointValues[13]
                   +NF_N_Q_Q4_2D_EdgeWeight3[4]*PointValues[14]);
  Functionals[15]=( NF_N_Q_Q4_2D_EdgeWeight3[0]*PointValues[15]
                   +NF_N_Q_Q4_2D_EdgeWeight3[1]*PointValues[16]
                   +NF_N_Q_Q4_2D_EdgeWeight3[2]*PointValues[17]
                   +NF_N_Q_Q4_2D_EdgeWeight3[3]*PointValues[18]
                   +NF_N_Q_Q4_2D_EdgeWeight3[4]*PointValues[19]);

  Functionals[16] = NF_N_Q_Q4_2D_CellWeight0[0]*PointValues[20]
                   +NF_N_Q_Q4_2D_CellWeight0[1]*PointValues[21]
                   +NF_N_Q_Q4_2D_CellWeight0[2]*PointValues[22]
                   +NF_N_Q_Q4_2D_CellWeight0[3]*PointValues[23]
                   +NF_N_Q_Q4_2D_CellWeight0[4]*PointValues[24]
                   +NF_N_Q_Q4_2D_CellWeight0[5]*PointValues[25]
                   +NF_N_Q_Q4_2D_CellWeight0[6]*PointValues[26]
                   +NF_N_Q_Q4_2D_CellWeight0[7]*PointValues[27]
                   +NF_N_Q_Q4_2D_CellWeight0[8]*PointValues[28]
                   +NF_N_Q_Q4_2D_CellWeight0[9]*PointValues[29]
                   +NF_N_Q_Q4_2D_CellWeight0[10]*PointValues[30]
                   +NF_N_Q_Q4_2D_CellWeight0[11]*PointValues[31]
                   +NF_N_Q_Q4_2D_CellWeight0[12]*PointValues[32]
                   +NF_N_Q_Q4_2D_CellWeight0[13]*PointValues[33]
                   +NF_N_Q_Q4_2D_CellWeight0[14]*PointValues[34]
                   +NF_N_Q_Q4_2D_CellWeight0[15]*PointValues[35]
                   +NF_N_Q_Q4_2D_CellWeight0[16]*PointValues[36]
                   +NF_N_Q_Q4_2D_CellWeight0[17]*PointValues[37]
                   +NF_N_Q_Q4_2D_CellWeight0[18]*PointValues[38]
                   +NF_N_Q_Q4_2D_CellWeight0[19]*PointValues[39]
                   +NF_N_Q_Q4_2D_CellWeight0[20]*PointValues[40]
                   +NF_N_Q_Q4_2D_CellWeight0[21]*PointValues[41]
                   +NF_N_Q_Q4_2D_CellWeight0[22]*PointValues[42]
                   +NF_N_Q_Q4_2D_CellWeight0[23]*PointValues[43]
                   +NF_N_Q_Q4_2D_CellWeight0[24]*PointValues[44];
  Functionals[17] = NF_N_Q_Q4_2D_CellWeight1[0]*PointValues[20]
                   +NF_N_Q_Q4_2D_CellWeight1[1]*PointValues[21]
                   +NF_N_Q_Q4_2D_CellWeight1[2]*PointValues[22]
                   +NF_N_Q_Q4_2D_CellWeight1[3]*PointValues[23]
                   +NF_N_Q_Q4_2D_CellWeight1[4]*PointValues[24]
                   +NF_N_Q_Q4_2D_CellWeight1[5]*PointValues[25]
                   +NF_N_Q_Q4_2D_CellWeight1[6]*PointValues[26]
                   +NF_N_Q_Q4_2D_CellWeight1[7]*PointValues[27]
                   +NF_N_Q_Q4_2D_CellWeight1[8]*PointValues[28]
                   +NF_N_Q_Q4_2D_CellWeight1[9]*PointValues[29]
                   +NF_N_Q_Q4_2D_CellWeight1[10]*PointValues[30]
                   +NF_N_Q_Q4_2D_CellWeight1[11]*PointValues[31]
                   +NF_N_Q_Q4_2D_CellWeight1[12]*PointValues[32]
                   +NF_N_Q_Q4_2D_CellWeight1[13]*PointValues[33]
                   +NF_N_Q_Q4_2D_CellWeight1[14]*PointValues[34]
                   +NF_N_Q_Q4_2D_CellWeight1[15]*PointValues[35]
                   +NF_N_Q_Q4_2D_CellWeight1[16]*PointValues[36]
                   +NF_N_Q_Q4_2D_CellWeight1[17]*PointValues[37]
                   +NF_N_Q_Q4_2D_CellWeight1[18]*PointValues[38]
                   +NF_N_Q_Q4_2D_CellWeight1[19]*PointValues[39]
                   +NF_N_Q_Q4_2D_CellWeight1[20]*PointValues[40]
                   +NF_N_Q_Q4_2D_CellWeight1[21]*PointValues[41]
                   +NF_N_Q_Q4_2D_CellWeight1[22]*PointValues[42]
                   +NF_N_Q_Q4_2D_CellWeight1[23]*PointValues[43]
                   +NF_N_Q_Q4_2D_CellWeight1[24]*PointValues[44];
  Functionals[18] = NF_N_Q_Q4_2D_CellWeight2[0]*PointValues[20]
                   +NF_N_Q_Q4_2D_CellWeight2[1]*PointValues[21]
                   +NF_N_Q_Q4_2D_CellWeight2[2]*PointValues[22]
                   +NF_N_Q_Q4_2D_CellWeight2[3]*PointValues[23]
                   +NF_N_Q_Q4_2D_CellWeight2[4]*PointValues[24]
                   +NF_N_Q_Q4_2D_CellWeight2[5]*PointValues[25]
                   +NF_N_Q_Q4_2D_CellWeight2[6]*PointValues[26]
                   +NF_N_Q_Q4_2D_CellWeight2[7]*PointValues[27]
                   +NF_N_Q_Q4_2D_CellWeight2[8]*PointValues[28]
                   +NF_N_Q_Q4_2D_CellWeight2[9]*PointValues[29]
                   +NF_N_Q_Q4_2D_CellWeight2[10]*PointValues[30]
                   +NF_N_Q_Q4_2D_CellWeight2[11]*PointValues[31]
                   +NF_N_Q_Q4_2D_CellWeight2[12]*PointValues[32]
                   +NF_N_Q_Q4_2D_CellWeight2[13]*PointValues[33]
                   +NF_N_Q_Q4_2D_CellWeight2[14]*PointValues[34]
                   +NF_N_Q_Q4_2D_CellWeight2[15]*PointValues[35]
                   +NF_N_Q_Q4_2D_CellWeight2[16]*PointValues[36]
                   +NF_N_Q_Q4_2D_CellWeight2[17]*PointValues[37]
                   +NF_N_Q_Q4_2D_CellWeight2[18]*PointValues[38]
                   +NF_N_Q_Q4_2D_CellWeight2[19]*PointValues[39]
                   +NF_N_Q_Q4_2D_CellWeight2[20]*PointValues[40]
                   +NF_N_Q_Q4_2D_CellWeight2[21]*PointValues[41]
                   +NF_N_Q_Q4_2D_CellWeight2[22]*PointValues[42]
                   +NF_N_Q_Q4_2D_CellWeight2[23]*PointValues[43]
                   +NF_N_Q_Q4_2D_CellWeight2[24]*PointValues[44];
  Functionals[19] = NF_N_Q_Q4_2D_CellWeight3[0]*PointValues[20]
                   +NF_N_Q_Q4_2D_CellWeight3[1]*PointValues[21]
                   +NF_N_Q_Q4_2D_CellWeight3[2]*PointValues[22]
                   +NF_N_Q_Q4_2D_CellWeight3[3]*PointValues[23]
                   +NF_N_Q_Q4_2D_CellWeight3[4]*PointValues[24]
                   +NF_N_Q_Q4_2D_CellWeight3[5]*PointValues[25]
                   +NF_N_Q_Q4_2D_CellWeight3[6]*PointValues[26]
                   +NF_N_Q_Q4_2D_CellWeight3[7]*PointValues[27]
                   +NF_N_Q_Q4_2D_CellWeight3[8]*PointValues[28]
                   +NF_N_Q_Q4_2D_CellWeight3[9]*PointValues[29]
                   +NF_N_Q_Q4_2D_CellWeight3[10]*PointValues[30]
                   +NF_N_Q_Q4_2D_CellWeight3[11]*PointValues[31]
                   +NF_N_Q_Q4_2D_CellWeight3[12]*PointValues[32]
                   +NF_N_Q_Q4_2D_CellWeight3[13]*PointValues[33]
                   +NF_N_Q_Q4_2D_CellWeight3[14]*PointValues[34]
                   +NF_N_Q_Q4_2D_CellWeight3[15]*PointValues[35]
                   +NF_N_Q_Q4_2D_CellWeight3[16]*PointValues[36]
                   +NF_N_Q_Q4_2D_CellWeight3[17]*PointValues[37]
                   +NF_N_Q_Q4_2D_CellWeight3[18]*PointValues[38]
                   +NF_N_Q_Q4_2D_CellWeight3[19]*PointValues[39]
                   +NF_N_Q_Q4_2D_CellWeight3[20]*PointValues[40]
                   +NF_N_Q_Q4_2D_CellWeight3[21]*PointValues[41]
                   +NF_N_Q_Q4_2D_CellWeight3[22]*PointValues[42]
                   +NF_N_Q_Q4_2D_CellWeight3[23]*PointValues[43]
                   +NF_N_Q_Q4_2D_CellWeight3[24]*PointValues[44];
  Functionals[20] = NF_N_Q_Q4_2D_CellWeight4[0]*PointValues[20]
                   +NF_N_Q_Q4_2D_CellWeight4[1]*PointValues[21]
                   +NF_N_Q_Q4_2D_CellWeight4[2]*PointValues[22]
                   +NF_N_Q_Q4_2D_CellWeight4[3]*PointValues[23]
                   +NF_N_Q_Q4_2D_CellWeight4[4]*PointValues[24]
                   +NF_N_Q_Q4_2D_CellWeight4[5]*PointValues[25]
                   +NF_N_Q_Q4_2D_CellWeight4[6]*PointValues[26]
                   +NF_N_Q_Q4_2D_CellWeight4[7]*PointValues[27]
                   +NF_N_Q_Q4_2D_CellWeight4[8]*PointValues[28]
                   +NF_N_Q_Q4_2D_CellWeight4[9]*PointValues[29]
                   +NF_N_Q_Q4_2D_CellWeight4[10]*PointValues[30]
                   +NF_N_Q_Q4_2D_CellWeight4[11]*PointValues[31]
                   +NF_N_Q_Q4_2D_CellWeight4[12]*PointValues[32]
                   +NF_N_Q_Q4_2D_CellWeight4[13]*PointValues[33]
                   +NF_N_Q_Q4_2D_CellWeight4[14]*PointValues[34]
                   +NF_N_Q_Q4_2D_CellWeight4[15]*PointValues[35]
                   +NF_N_Q_Q4_2D_CellWeight4[16]*PointValues[36]
                   +NF_N_Q_Q4_2D_CellWeight4[17]*PointValues[37]
                   +NF_N_Q_Q4_2D_CellWeight4[18]*PointValues[38]
                   +NF_N_Q_Q4_2D_CellWeight4[19]*PointValues[39]
                   +NF_N_Q_Q4_2D_CellWeight4[20]*PointValues[40]
                   +NF_N_Q_Q4_2D_CellWeight4[21]*PointValues[41]
                   +NF_N_Q_Q4_2D_CellWeight4[22]*PointValues[42]
                   +NF_N_Q_Q4_2D_CellWeight4[23]*PointValues[43]
                   +NF_N_Q_Q4_2D_CellWeight4[24]*PointValues[44];
  Functionals[21] = NF_N_Q_Q4_2D_CellWeight5[0]*PointValues[20]
                   +NF_N_Q_Q4_2D_CellWeight5[1]*PointValues[21]
                   +NF_N_Q_Q4_2D_CellWeight5[2]*PointValues[22]
                   +NF_N_Q_Q4_2D_CellWeight5[3]*PointValues[23]
                   +NF_N_Q_Q4_2D_CellWeight5[4]*PointValues[24]
                   +NF_N_Q_Q4_2D_CellWeight5[5]*PointValues[25]
                   +NF_N_Q_Q4_2D_CellWeight5[6]*PointValues[26]
                   +NF_N_Q_Q4_2D_CellWeight5[7]*PointValues[27]
                   +NF_N_Q_Q4_2D_CellWeight5[8]*PointValues[28]
                   +NF_N_Q_Q4_2D_CellWeight5[9]*PointValues[29]
                   +NF_N_Q_Q4_2D_CellWeight5[10]*PointValues[30]
                   +NF_N_Q_Q4_2D_CellWeight5[11]*PointValues[31]
                   +NF_N_Q_Q4_2D_CellWeight5[12]*PointValues[32]
                   +NF_N_Q_Q4_2D_CellWeight5[13]*PointValues[33]
                   +NF_N_Q_Q4_2D_CellWeight5[14]*PointValues[34]
                   +NF_N_Q_Q4_2D_CellWeight5[15]*PointValues[35]
                   +NF_N_Q_Q4_2D_CellWeight5[16]*PointValues[36]
                   +NF_N_Q_Q4_2D_CellWeight5[17]*PointValues[37]
                   +NF_N_Q_Q4_2D_CellWeight5[18]*PointValues[38]
                   +NF_N_Q_Q4_2D_CellWeight5[19]*PointValues[39]
                   +NF_N_Q_Q4_2D_CellWeight5[20]*PointValues[40]
                   +NF_N_Q_Q4_2D_CellWeight5[21]*PointValues[41]
                   +NF_N_Q_Q4_2D_CellWeight5[22]*PointValues[42]
                   +NF_N_Q_Q4_2D_CellWeight5[23]*PointValues[43]
                   +NF_N_Q_Q4_2D_CellWeight5[24]*PointValues[44];

  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
    {
      Functionals[4] = -Functionals[4];
      Functionals[12] = -Functionals[12];
    }
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
    {
      Functionals[5] = -Functionals[5];
      Functionals[13] = -Functionals[13];
    }
    if(Cell->GetVertex(2) > Cell->GetVertex(3))
    {
      Functionals[6] = -Functionals[6];
      Functionals[14] = -Functionals[14];
    }
    if(Cell->GetVertex(3) > Cell->GetVertex(0))
    {
      Functionals[7] = -Functionals[7];
      Functionals[15] = -Functionals[15];
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
        Functionals[ 4] = -Functionals[ 4];
        Functionals[12] = -Functionals[12];
      }
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[ 5] = -Functionals[ 5];
        Functionals[13] = -Functionals[13];
      }
    } // endif neigh

    neigh = Cell->GetJoint(2)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[ 6] = -Functionals[ 6];
        Functionals[14] = -Functionals[14];
      }
    } // endif neigh

    neigh = Cell->GetJoint(3)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[ 7] = -Functionals[ 7];
        Functionals[15] = -Functionals[15];
      }
    } // endif neigh
  } // endif Cell
}

void NF_N_Q_Q4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                           double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] =( NF_N_Q_Q4_2D_EdgeWeight0[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight0[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight0[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight0[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight0[4]*PointValues[4]);

  Functionals[1] =( NF_N_Q_Q4_2D_EdgeWeight1[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight1[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight1[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight1[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight1[4]*PointValues[4]);

  Functionals[2] =( NF_N_Q_Q4_2D_EdgeWeight2[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight2[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight2[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight2[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight2[4]*PointValues[4]);

  Functionals[3] =( NF_N_Q_Q4_2D_EdgeWeight3[0]*PointValues[0]
                   +NF_N_Q_Q4_2D_EdgeWeight3[1]*PointValues[1]
                   +NF_N_Q_Q4_2D_EdgeWeight3[2]*PointValues[2]
                   +NF_N_Q_Q4_2D_EdgeWeight3[3]*PointValues[3]
                   +NF_N_Q_Q4_2D_EdgeWeight3[4]*PointValues[4]);

  if(Joint != -1)
  {
    // if(Cell->GetVertex(Joint) > Cell->GetVertex((Joint+1)%4))
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

TNodalFunctional2D *NF_N_Q_Q4_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_Q4_2D, 22, 4, 45, 5, NF_N_Q_Q4_2D_Xi, NF_N_Q_Q4_2D_Eta,
         NF_N_Q_Q4_2D_T, NF_N_Q_Q4_2D_EvalAll, NF_N_Q_Q4_2D_EvalEdge);
