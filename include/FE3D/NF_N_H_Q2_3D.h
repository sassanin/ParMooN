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
   
/* for all functionals */
static double NF_N_H_Q2_3D_Xi[] = {
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0,
          0,
          0,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531
};

static double NF_N_H_Q2_3D_Eta[] = {
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531
};

static double NF_N_H_Q2_3D_Zeta[] = {
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531
};

/* face 0 */
static double NF_N_H_Q2_3D_F0_Xi[]   = {
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F0_Eta[]   = {
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F0_Zeta[]   = {
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1
       };

/* face 1 */
static double NF_N_H_Q2_3D_F1_Xi[]   = {
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F1_Eta[]   = {
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1
       };
static double NF_N_H_Q2_3D_F1_Zeta[]   = {
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531
       };

/* face 2 */
static double NF_N_H_Q2_3D_F2_Xi[]   = {
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1
       };
static double NF_N_H_Q2_3D_F2_Eta[]   = {
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F2_Zeta[]   = {
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531
       };

/* face 3 */
static double NF_N_H_Q2_3D_F3_Xi[]   = {
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0,
          0,
          0,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F3_Eta[]   = {
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1
       };
static double NF_N_H_Q2_3D_F3_Zeta[]   = {
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531
       };

/* face 4 */
static double NF_N_H_Q2_3D_F4_Xi[]   = {
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1,
          -1
       };
static double NF_N_H_Q2_3D_F4_Eta[]   = {
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F4_Zeta[]   = {
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531
       };

/* face 5 */
static double NF_N_H_Q2_3D_F5_Xi[]   = {
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0,
          0,
          0.7745966692414833770358531,
          0.7745966692414833770358531,
          0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F5_Eta[]   = {
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531,
          -0.7745966692414833770358531,
          0,
          0.7745966692414833770358531
       };
static double NF_N_H_Q2_3D_F5_Zeta[]   = {
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1
       };

static double NF_N_H_Q2_3D_T[]   = {
          0.1127016653792583114820735,
          0.5,
          0.8872983346207416885179265,
          0.1127016653792583114820735,
          0.5,
          0.8872983346207416885179265,
          0.1127016653792583114820735,
          0.5,
          0.8872983346207416885179265
       };
static double NF_N_H_Q2_3D_S[]   = {
          0.1127016653792583114820735,
          0.1127016653792583114820735,
          0.1127016653792583114820735,
          0.5,
          0.5,
          0.5,
          0.8872983346207416885179265,
          0.8872983346207416885179265,
          0.8872983346207416885179265
       };

static double *NF_N_H_Q2_3D_XiArray[6] = {
                        NF_N_H_Q2_3D_F0_Xi,
                        NF_N_H_Q2_3D_F1_Xi,
                        NF_N_H_Q2_3D_F2_Xi,
                        NF_N_H_Q2_3D_F3_Xi,
                        NF_N_H_Q2_3D_F4_Xi,
                        NF_N_H_Q2_3D_F5_Xi };

static double *NF_N_H_Q2_3D_EtaArray[6] = {
                        NF_N_H_Q2_3D_F0_Eta,
                        NF_N_H_Q2_3D_F1_Eta,
                        NF_N_H_Q2_3D_F2_Eta,
                        NF_N_H_Q2_3D_F3_Eta,
                        NF_N_H_Q2_3D_F4_Eta,
                        NF_N_H_Q2_3D_F5_Eta };

static double *NF_N_H_Q2_3D_ZetaArray[6] = {
                        NF_N_H_Q2_3D_F0_Zeta,
                        NF_N_H_Q2_3D_F1_Zeta,
                        NF_N_H_Q2_3D_F2_Zeta,
                        NF_N_H_Q2_3D_F3_Zeta,
                        NF_N_H_Q2_3D_F4_Zeta,
                        NF_N_H_Q2_3D_F5_Zeta };

static double NF_N_H_Q2_3D_FaceWeight0[] = {
          0.07716049382716049382716049,
          0.1234567901234567901234568,
          0.07716049382716049382716049,
          0.1234567901234567901234568,
          0.1975308641975308641975309,
          0.1234567901234567901234568,
          0.07716049382716049382716049,
          0.1234567901234567901234568,
          0.07716049382716049382716049
      };

static double NF_N_H_Q2_3D_FaceWeight1[] = {
          -0.1793047845466396706101512,
          0,
          0.1793047845466396706101512,
          -0.2868876552746234729762419,
          0,
          0.2868876552746234729762419,
          -0.1793047845466396706101512,
          0,
          0.1793047845466396706101512
      };

static double NF_N_H_Q2_3D_FaceWeight2[] = {
          -0.1793047845466396706101512,
          -0.2868876552746234729762419,
          -0.1793047845466396706101512,
          0,
          0,
          0,
          0.1793047845466396706101512,
          0.2868876552746234729762419,
          0.1793047845466396706101512
      };

static double NF_N_H_Q2_3D_CellWeight0[] = {
          0.02143347050754458161865569,
          0.03429355281207133058984911,
          0.02143347050754458161865569,
          0.03429355281207133058984911,
          0.05486968449931412894375857,
          0.03429355281207133058984911,
          0.02143347050754458161865569,
          0.03429355281207133058984911,
          0.02143347050754458161865569,
          0.03429355281207133058984911,
          0.05486968449931412894375857,
          0.03429355281207133058984911,
          0.05486968449931412894375857,
          0.08779149519890260631001372,
          0.05486968449931412894375857,
          0.03429355281207133058984911,
          0.05486968449931412894375857,
          0.03429355281207133058984911,
          0.02143347050754458161865569,
          0.03429355281207133058984911,
          0.02143347050754458161865569,
          0.03429355281207133058984911,
          0.05486968449931412894375857,
          0.03429355281207133058984911,
          0.02143347050754458161865569,
          0.03429355281207133058984911,
          0.02143347050754458161865569
      };

void NF_N_H_Q2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  TJoint *joint;
  TBaseCell *neigh;
  int maptype;

  Functionals[0] = 
          +NF_N_H_Q2_3D_FaceWeight0[0]*PointValues[0]
          +NF_N_H_Q2_3D_FaceWeight0[1]*PointValues[1]
          +NF_N_H_Q2_3D_FaceWeight0[2]*PointValues[2]
          +NF_N_H_Q2_3D_FaceWeight0[3]*PointValues[3]
          +NF_N_H_Q2_3D_FaceWeight0[4]*PointValues[4]
          +NF_N_H_Q2_3D_FaceWeight0[5]*PointValues[5]
          +NF_N_H_Q2_3D_FaceWeight0[6]*PointValues[6]
          +NF_N_H_Q2_3D_FaceWeight0[7]*PointValues[7]
          +NF_N_H_Q2_3D_FaceWeight0[8]*PointValues[8];
  Functionals[1] = 
          +NF_N_H_Q2_3D_FaceWeight0[0]*PointValues[9]
          +NF_N_H_Q2_3D_FaceWeight0[1]*PointValues[10]
          +NF_N_H_Q2_3D_FaceWeight0[2]*PointValues[11]
          +NF_N_H_Q2_3D_FaceWeight0[3]*PointValues[12]
          +NF_N_H_Q2_3D_FaceWeight0[4]*PointValues[13]
          +NF_N_H_Q2_3D_FaceWeight0[5]*PointValues[14]
          +NF_N_H_Q2_3D_FaceWeight0[6]*PointValues[15]
          +NF_N_H_Q2_3D_FaceWeight0[7]*PointValues[16]
          +NF_N_H_Q2_3D_FaceWeight0[8]*PointValues[17];
  Functionals[2] = 
          +NF_N_H_Q2_3D_FaceWeight0[0]*PointValues[18]
          +NF_N_H_Q2_3D_FaceWeight0[1]*PointValues[19]
          +NF_N_H_Q2_3D_FaceWeight0[2]*PointValues[20]
          +NF_N_H_Q2_3D_FaceWeight0[3]*PointValues[21]
          +NF_N_H_Q2_3D_FaceWeight0[4]*PointValues[22]
          +NF_N_H_Q2_3D_FaceWeight0[5]*PointValues[23]
          +NF_N_H_Q2_3D_FaceWeight0[6]*PointValues[24]
          +NF_N_H_Q2_3D_FaceWeight0[7]*PointValues[25]
          +NF_N_H_Q2_3D_FaceWeight0[8]*PointValues[26];
  Functionals[3] = 
          +NF_N_H_Q2_3D_FaceWeight0[0]*PointValues[27]
          +NF_N_H_Q2_3D_FaceWeight0[1]*PointValues[28]
          +NF_N_H_Q2_3D_FaceWeight0[2]*PointValues[29]
          +NF_N_H_Q2_3D_FaceWeight0[3]*PointValues[30]
          +NF_N_H_Q2_3D_FaceWeight0[4]*PointValues[31]
          +NF_N_H_Q2_3D_FaceWeight0[5]*PointValues[32]
          +NF_N_H_Q2_3D_FaceWeight0[6]*PointValues[33]
          +NF_N_H_Q2_3D_FaceWeight0[7]*PointValues[34]
          +NF_N_H_Q2_3D_FaceWeight0[8]*PointValues[35];
  Functionals[4] = 
          +NF_N_H_Q2_3D_FaceWeight0[0]*PointValues[36]
          +NF_N_H_Q2_3D_FaceWeight0[1]*PointValues[37]
          +NF_N_H_Q2_3D_FaceWeight0[2]*PointValues[38]
          +NF_N_H_Q2_3D_FaceWeight0[3]*PointValues[39]
          +NF_N_H_Q2_3D_FaceWeight0[4]*PointValues[40]
          +NF_N_H_Q2_3D_FaceWeight0[5]*PointValues[41]
          +NF_N_H_Q2_3D_FaceWeight0[6]*PointValues[42]
          +NF_N_H_Q2_3D_FaceWeight0[7]*PointValues[43]
          +NF_N_H_Q2_3D_FaceWeight0[8]*PointValues[44];
  Functionals[5] = 
          +NF_N_H_Q2_3D_FaceWeight0[0]*PointValues[45]
          +NF_N_H_Q2_3D_FaceWeight0[1]*PointValues[46]
          +NF_N_H_Q2_3D_FaceWeight0[2]*PointValues[47]
          +NF_N_H_Q2_3D_FaceWeight0[3]*PointValues[48]
          +NF_N_H_Q2_3D_FaceWeight0[4]*PointValues[49]
          +NF_N_H_Q2_3D_FaceWeight0[5]*PointValues[50]
          +NF_N_H_Q2_3D_FaceWeight0[6]*PointValues[51]
          +NF_N_H_Q2_3D_FaceWeight0[7]*PointValues[52]
          +NF_N_H_Q2_3D_FaceWeight0[8]*PointValues[53];

  Functionals[6] = 
          +NF_N_H_Q2_3D_FaceWeight1[0]*PointValues[0]
          +NF_N_H_Q2_3D_FaceWeight1[1]*PointValues[1]
          +NF_N_H_Q2_3D_FaceWeight1[2]*PointValues[2]
          +NF_N_H_Q2_3D_FaceWeight1[3]*PointValues[3]
          +NF_N_H_Q2_3D_FaceWeight1[4]*PointValues[4]
          +NF_N_H_Q2_3D_FaceWeight1[5]*PointValues[5]
          +NF_N_H_Q2_3D_FaceWeight1[6]*PointValues[6]
          +NF_N_H_Q2_3D_FaceWeight1[7]*PointValues[7]
          +NF_N_H_Q2_3D_FaceWeight1[8]*PointValues[8];
  Functionals[7] = 
          +NF_N_H_Q2_3D_FaceWeight1[0]*PointValues[9]
          +NF_N_H_Q2_3D_FaceWeight1[1]*PointValues[10]
          +NF_N_H_Q2_3D_FaceWeight1[2]*PointValues[11]
          +NF_N_H_Q2_3D_FaceWeight1[3]*PointValues[12]
          +NF_N_H_Q2_3D_FaceWeight1[4]*PointValues[13]
          +NF_N_H_Q2_3D_FaceWeight1[5]*PointValues[14]
          +NF_N_H_Q2_3D_FaceWeight1[6]*PointValues[15]
          +NF_N_H_Q2_3D_FaceWeight1[7]*PointValues[16]
          +NF_N_H_Q2_3D_FaceWeight1[8]*PointValues[17];
  Functionals[8] = 
          +NF_N_H_Q2_3D_FaceWeight1[0]*PointValues[18]
          +NF_N_H_Q2_3D_FaceWeight1[1]*PointValues[19]
          +NF_N_H_Q2_3D_FaceWeight1[2]*PointValues[20]
          +NF_N_H_Q2_3D_FaceWeight1[3]*PointValues[21]
          +NF_N_H_Q2_3D_FaceWeight1[4]*PointValues[22]
          +NF_N_H_Q2_3D_FaceWeight1[5]*PointValues[23]
          +NF_N_H_Q2_3D_FaceWeight1[6]*PointValues[24]
          +NF_N_H_Q2_3D_FaceWeight1[7]*PointValues[25]
          +NF_N_H_Q2_3D_FaceWeight1[8]*PointValues[26];
  Functionals[9] = 
          +NF_N_H_Q2_3D_FaceWeight1[0]*PointValues[27]
          +NF_N_H_Q2_3D_FaceWeight1[1]*PointValues[28]
          +NF_N_H_Q2_3D_FaceWeight1[2]*PointValues[29]
          +NF_N_H_Q2_3D_FaceWeight1[3]*PointValues[30]
          +NF_N_H_Q2_3D_FaceWeight1[4]*PointValues[31]
          +NF_N_H_Q2_3D_FaceWeight1[5]*PointValues[32]
          +NF_N_H_Q2_3D_FaceWeight1[6]*PointValues[33]
          +NF_N_H_Q2_3D_FaceWeight1[7]*PointValues[34]
          +NF_N_H_Q2_3D_FaceWeight1[8]*PointValues[35];
  Functionals[10] = 
          +NF_N_H_Q2_3D_FaceWeight1[0]*PointValues[36]
          +NF_N_H_Q2_3D_FaceWeight1[1]*PointValues[37]
          +NF_N_H_Q2_3D_FaceWeight1[2]*PointValues[38]
          +NF_N_H_Q2_3D_FaceWeight1[3]*PointValues[39]
          +NF_N_H_Q2_3D_FaceWeight1[4]*PointValues[40]
          +NF_N_H_Q2_3D_FaceWeight1[5]*PointValues[41]
          +NF_N_H_Q2_3D_FaceWeight1[6]*PointValues[42]
          +NF_N_H_Q2_3D_FaceWeight1[7]*PointValues[43]
          +NF_N_H_Q2_3D_FaceWeight1[8]*PointValues[44];
  Functionals[11] = 
          +NF_N_H_Q2_3D_FaceWeight1[0]*PointValues[45]
          +NF_N_H_Q2_3D_FaceWeight1[1]*PointValues[46]
          +NF_N_H_Q2_3D_FaceWeight1[2]*PointValues[47]
          +NF_N_H_Q2_3D_FaceWeight1[3]*PointValues[48]
          +NF_N_H_Q2_3D_FaceWeight1[4]*PointValues[49]
          +NF_N_H_Q2_3D_FaceWeight1[5]*PointValues[50]
          +NF_N_H_Q2_3D_FaceWeight1[6]*PointValues[51]
          +NF_N_H_Q2_3D_FaceWeight1[7]*PointValues[52]
          +NF_N_H_Q2_3D_FaceWeight1[8]*PointValues[53];

  Functionals[12] = 
          +NF_N_H_Q2_3D_FaceWeight2[0]*PointValues[0]
          +NF_N_H_Q2_3D_FaceWeight2[1]*PointValues[1]
          +NF_N_H_Q2_3D_FaceWeight2[2]*PointValues[2]
          +NF_N_H_Q2_3D_FaceWeight2[3]*PointValues[3]
          +NF_N_H_Q2_3D_FaceWeight2[4]*PointValues[4]
          +NF_N_H_Q2_3D_FaceWeight2[5]*PointValues[5]
          +NF_N_H_Q2_3D_FaceWeight2[6]*PointValues[6]
          +NF_N_H_Q2_3D_FaceWeight2[7]*PointValues[7]
          +NF_N_H_Q2_3D_FaceWeight2[8]*PointValues[8];
  Functionals[13] = 
          +NF_N_H_Q2_3D_FaceWeight2[0]*PointValues[9]
          +NF_N_H_Q2_3D_FaceWeight2[1]*PointValues[10]
          +NF_N_H_Q2_3D_FaceWeight2[2]*PointValues[11]
          +NF_N_H_Q2_3D_FaceWeight2[3]*PointValues[12]
          +NF_N_H_Q2_3D_FaceWeight2[4]*PointValues[13]
          +NF_N_H_Q2_3D_FaceWeight2[5]*PointValues[14]
          +NF_N_H_Q2_3D_FaceWeight2[6]*PointValues[15]
          +NF_N_H_Q2_3D_FaceWeight2[7]*PointValues[16]
          +NF_N_H_Q2_3D_FaceWeight2[8]*PointValues[17];
  Functionals[14] = 
          +NF_N_H_Q2_3D_FaceWeight2[0]*PointValues[18]
          +NF_N_H_Q2_3D_FaceWeight2[1]*PointValues[19]
          +NF_N_H_Q2_3D_FaceWeight2[2]*PointValues[20]
          +NF_N_H_Q2_3D_FaceWeight2[3]*PointValues[21]
          +NF_N_H_Q2_3D_FaceWeight2[4]*PointValues[22]
          +NF_N_H_Q2_3D_FaceWeight2[5]*PointValues[23]
          +NF_N_H_Q2_3D_FaceWeight2[6]*PointValues[24]
          +NF_N_H_Q2_3D_FaceWeight2[7]*PointValues[25]
          +NF_N_H_Q2_3D_FaceWeight2[8]*PointValues[26];
  Functionals[15] = 
          +NF_N_H_Q2_3D_FaceWeight2[0]*PointValues[27]
          +NF_N_H_Q2_3D_FaceWeight2[1]*PointValues[28]
          +NF_N_H_Q2_3D_FaceWeight2[2]*PointValues[29]
          +NF_N_H_Q2_3D_FaceWeight2[3]*PointValues[30]
          +NF_N_H_Q2_3D_FaceWeight2[4]*PointValues[31]
          +NF_N_H_Q2_3D_FaceWeight2[5]*PointValues[32]
          +NF_N_H_Q2_3D_FaceWeight2[6]*PointValues[33]
          +NF_N_H_Q2_3D_FaceWeight2[7]*PointValues[34]
          +NF_N_H_Q2_3D_FaceWeight2[8]*PointValues[35];
  Functionals[16] = 
          +NF_N_H_Q2_3D_FaceWeight2[0]*PointValues[36]
          +NF_N_H_Q2_3D_FaceWeight2[1]*PointValues[37]
          +NF_N_H_Q2_3D_FaceWeight2[2]*PointValues[38]
          +NF_N_H_Q2_3D_FaceWeight2[3]*PointValues[39]
          +NF_N_H_Q2_3D_FaceWeight2[4]*PointValues[40]
          +NF_N_H_Q2_3D_FaceWeight2[5]*PointValues[41]
          +NF_N_H_Q2_3D_FaceWeight2[6]*PointValues[42]
          +NF_N_H_Q2_3D_FaceWeight2[7]*PointValues[43]
          +NF_N_H_Q2_3D_FaceWeight2[8]*PointValues[44];
  Functionals[17] = 
          +NF_N_H_Q2_3D_FaceWeight2[0]*PointValues[45]
          +NF_N_H_Q2_3D_FaceWeight2[1]*PointValues[46]
          +NF_N_H_Q2_3D_FaceWeight2[2]*PointValues[47]
          +NF_N_H_Q2_3D_FaceWeight2[3]*PointValues[48]
          +NF_N_H_Q2_3D_FaceWeight2[4]*PointValues[49]
          +NF_N_H_Q2_3D_FaceWeight2[5]*PointValues[50]
          +NF_N_H_Q2_3D_FaceWeight2[6]*PointValues[51]
          +NF_N_H_Q2_3D_FaceWeight2[7]*PointValues[52]
          +NF_N_H_Q2_3D_FaceWeight2[8]*PointValues[53];

  Functionals[18] = 
          +NF_N_H_Q2_3D_CellWeight0[0]*PointValues[54]
          +NF_N_H_Q2_3D_CellWeight0[1]*PointValues[55]
          +NF_N_H_Q2_3D_CellWeight0[2]*PointValues[56]
          +NF_N_H_Q2_3D_CellWeight0[3]*PointValues[57]
          +NF_N_H_Q2_3D_CellWeight0[4]*PointValues[58]
          +NF_N_H_Q2_3D_CellWeight0[5]*PointValues[59]
          +NF_N_H_Q2_3D_CellWeight0[6]*PointValues[60]
          +NF_N_H_Q2_3D_CellWeight0[7]*PointValues[61]
          +NF_N_H_Q2_3D_CellWeight0[8]*PointValues[62]
          +NF_N_H_Q2_3D_CellWeight0[9]*PointValues[63]
          +NF_N_H_Q2_3D_CellWeight0[10]*PointValues[64]
          +NF_N_H_Q2_3D_CellWeight0[11]*PointValues[65]
          +NF_N_H_Q2_3D_CellWeight0[12]*PointValues[66]
          +NF_N_H_Q2_3D_CellWeight0[13]*PointValues[67]
          +NF_N_H_Q2_3D_CellWeight0[14]*PointValues[68]
          +NF_N_H_Q2_3D_CellWeight0[15]*PointValues[69]
          +NF_N_H_Q2_3D_CellWeight0[16]*PointValues[70]
          +NF_N_H_Q2_3D_CellWeight0[17]*PointValues[71]
          +NF_N_H_Q2_3D_CellWeight0[18]*PointValues[72]
          +NF_N_H_Q2_3D_CellWeight0[19]*PointValues[73]
          +NF_N_H_Q2_3D_CellWeight0[20]*PointValues[74]
          +NF_N_H_Q2_3D_CellWeight0[21]*PointValues[75]
          +NF_N_H_Q2_3D_CellWeight0[22]*PointValues[76]
          +NF_N_H_Q2_3D_CellWeight0[23]*PointValues[77]
          +NF_N_H_Q2_3D_CellWeight0[24]*PointValues[78]
          +NF_N_H_Q2_3D_CellWeight0[25]*PointValues[79]
          +NF_N_H_Q2_3D_CellWeight0[26]*PointValues[80];

  if(Cell)
  {
    joint = Cell->GetJoint(0);
    maptype = joint->GetMapType();
    if(maptype == 1 || maptype == 2)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[6] = -Functionals[6];
      }
    }
    if(maptype == 2 || maptype == 3)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[12] = -Functionals[12];
      }
    }
  
    joint = Cell->GetJoint(1);
    maptype = joint->GetMapType();
    if(maptype == 1 || maptype == 2)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[7] = -Functionals[7];
      }
    }
    if(maptype == 2 || maptype == 3)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[13] = -Functionals[13];
      }
    }
  
    joint = Cell->GetJoint(2);
    maptype = joint->GetMapType();
    if(maptype == 1 || maptype == 2)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[8] = -Functionals[8];
      }
    }
    if(maptype == 2 || maptype == 3)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[14] = -Functionals[14];
      }
    }
  
    joint = Cell->GetJoint(3);
    maptype = joint->GetMapType();
    if(maptype == 1 || maptype == 2)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[9] = -Functionals[9];
      }
    }
    if(maptype == 2 || maptype == 3)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[15] = -Functionals[15];
      }
    }
  
    joint = Cell->GetJoint(4);
    maptype = joint->GetMapType();
    if(maptype == 1 || maptype == 2)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[10] = -Functionals[10];
      }
    }
    if(maptype == 2 || maptype == 3)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[16] = -Functionals[16];
      }
    }
  
    joint = Cell->GetJoint(5);
    maptype = joint->GetMapType();
    if(maptype == 1 || maptype == 2)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[11] = -Functionals[11];
      }
    }
    if(maptype == 2 || maptype == 3)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[17] = -Functionals[17];
      }
    }
  }
}

void NF_N_H_Q2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                           double *PointValues, double *Functionals)
{
  TJoint *joint;
  TBaseCell *neigh;
  int maptype;

  Functionals[0] = 
          +NF_N_H_Q2_3D_FaceWeight0[0]*PointValues[0]
          +NF_N_H_Q2_3D_FaceWeight0[1]*PointValues[1]
          +NF_N_H_Q2_3D_FaceWeight0[2]*PointValues[2]
          +NF_N_H_Q2_3D_FaceWeight0[3]*PointValues[3]
          +NF_N_H_Q2_3D_FaceWeight0[4]*PointValues[4]
          +NF_N_H_Q2_3D_FaceWeight0[5]*PointValues[5]
          +NF_N_H_Q2_3D_FaceWeight0[6]*PointValues[6]
          +NF_N_H_Q2_3D_FaceWeight0[7]*PointValues[7]
          +NF_N_H_Q2_3D_FaceWeight0[8]*PointValues[8];
  Functionals[1] = 
          +NF_N_H_Q2_3D_FaceWeight1[0]*PointValues[0]
          +NF_N_H_Q2_3D_FaceWeight1[1]*PointValues[1]
          +NF_N_H_Q2_3D_FaceWeight1[2]*PointValues[2]
          +NF_N_H_Q2_3D_FaceWeight1[3]*PointValues[3]
          +NF_N_H_Q2_3D_FaceWeight1[4]*PointValues[4]
          +NF_N_H_Q2_3D_FaceWeight1[5]*PointValues[5]
          +NF_N_H_Q2_3D_FaceWeight1[6]*PointValues[6]
          +NF_N_H_Q2_3D_FaceWeight1[7]*PointValues[7]
          +NF_N_H_Q2_3D_FaceWeight1[8]*PointValues[8];
  Functionals[2] = 
          +NF_N_H_Q2_3D_FaceWeight2[0]*PointValues[0]
          +NF_N_H_Q2_3D_FaceWeight2[1]*PointValues[1]
          +NF_N_H_Q2_3D_FaceWeight2[2]*PointValues[2]
          +NF_N_H_Q2_3D_FaceWeight2[3]*PointValues[3]
          +NF_N_H_Q2_3D_FaceWeight2[4]*PointValues[4]
          +NF_N_H_Q2_3D_FaceWeight2[5]*PointValues[5]
          +NF_N_H_Q2_3D_FaceWeight2[6]*PointValues[6]
          +NF_N_H_Q2_3D_FaceWeight2[7]*PointValues[7]
          +NF_N_H_Q2_3D_FaceWeight2[8]*PointValues[8];

  if(Cell)
  {
    joint = Cell->GetJoint(Joint);
    maptype = joint->GetMapType();
    if(maptype == 1 || maptype == 2)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[1] = -Functionals[1];
      }
    }
    if(maptype == 2 || maptype == 3)
    {
      neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
      if(neigh != NULL && neigh-Cell < 0)
      {
        Functionals[2] = -Functionals[2];
      }
    }
  }
}

static int NF_N_H_Q2_3D_N_AllFunctionals = 19;
static int NF_N_H_Q2_3D_N_PointsAll = 81;
static int NF_N_H_Q2_3D_N_FaceFunctionals[] = { 3, 3, 3, 3, 3, 3 };
static int NF_N_H_Q2_3D_N_PointsFace[] = { 9, 9, 9, 9, 9, 9 };

TNodalFunctional3D *NF_N_H_Q2_3D_Obj = new TNodalFunctional3D
        (NF_N_H_Q2_3D, NF_N_H_Q2_3D_N_AllFunctionals,
         NF_N_H_Q2_3D_N_FaceFunctionals, NF_N_H_Q2_3D_N_PointsAll,
         NF_N_H_Q2_3D_N_PointsFace,
         NF_N_H_Q2_3D_Xi, NF_N_H_Q2_3D_Eta, NF_N_H_Q2_3D_Zeta,
         NF_N_H_Q2_3D_XiArray, NF_N_H_Q2_3D_EtaArray,
         NF_N_H_Q2_3D_ZetaArray,
         NF_N_H_Q2_3D_T, NF_N_H_Q2_3D_S,
         NF_N_H_Q2_3D_EvalAll, NF_N_H_Q2_3D_EvalFace);
