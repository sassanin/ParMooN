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
   
// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of third order on hexahedra, 3D
// ***********************************************************************

//tschebyscheff point
static double BDDF3_tp_a = -0.923879532511287;
static double BDDF3_tp_b = -0.382683432365090;

static double NF_N_H_BDDF3_3D_Xi[] = {BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,
    BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a,
    1,1,1,1,1,1,1,1,1,1,
    BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a,
    -1, -1, -1, -1, -1, -1, -1, -1 ,-1, -1,
    BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a,
    // innere dofs
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526
};
static double NF_N_H_BDDF3_3D_Eta[]  = {BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a,
    -1, -1, -1, -1, -1, -1, -1, -1 ,-1, -1,
    BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a,
    1,1,1,1,1,1,1,1,1,1,
    BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,
    BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,
    //innere dofs
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526
};
static double NF_N_H_BDDF3_3D_Zeta[] = {-1, -1, -1, -1, -1, -1, -1, -1 ,-1, -1,
    BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,
    BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,
    BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,
    BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a,
    1,1,1,1,1,1,1,1,1,1,
    //innere dofs
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526
};

// face 0
static double NF_N_H_BDDF3_3D_F0_Xi[]   = {BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F0_Eta[]  = {BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F0_Zeta[] = { -1, -1, -1, -1, -1, -1, -1, -1 ,-1, -1};

// face 1                               1
static double NF_N_H_BDDF3_3D_F1_Xi[]   = {BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F1_Eta[]  = { -1, -1, -1, -1, -1, -1, -1, -1 ,-1, -1};
static double NF_N_H_BDDF3_3D_F1_Zeta[] = {BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a};

// face 2                               2
static double NF_N_H_BDDF3_3D_F2_Xi[]   = { 1,1,1,1,1,1,1,1,1,1 };
static double NF_N_H_BDDF3_3D_F2_Eta[]  = {BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F2_Zeta[] = {BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a};

 //face 3                               3
static double NF_N_H_BDDF3_3D_F3_Xi[]   = {BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F3_Eta[]  = { 1,1,1,1,1,1,1,1,1,1 };
static double NF_N_H_BDDF3_3D_F3_Zeta[] = {BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a};

// face 4                               4
static double NF_N_H_BDDF3_3D_F4_Xi[]   = { -1, -1, -1, -1, -1, -1, -1, -1 ,-1, -1};
static double NF_N_H_BDDF3_3D_F4_Eta[]  = {BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F4_Zeta[] = {BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a};

// face 5                               5
static double NF_N_H_BDDF3_3D_F5_Xi[]   = {BDDF3_tp_a,BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,BDDF3_tp_a,BDDF3_tp_b,-BDDF3_tp_b,-BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F5_Eta[]  = {BDDF3_tp_a,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a,-BDDF3_tp_a,-BDDF3_tp_b,BDDF3_tp_b,BDDF3_tp_a};
static double NF_N_H_BDDF3_3D_F5_Zeta[] = { 1,1,1,1,1,1,1,1,1,1 };

static double *NF_N_H_BDDF3_3D_XiArray[6] = {
                        NF_N_H_BDDF3_3D_F0_Xi,
                        NF_N_H_BDDF3_3D_F1_Xi,
                        NF_N_H_BDDF3_3D_F2_Xi,
                        NF_N_H_BDDF3_3D_F3_Xi,
                        NF_N_H_BDDF3_3D_F4_Xi,
                        NF_N_H_BDDF3_3D_F5_Xi };

static double *NF_N_H_BDDF3_3D_EtaArray[6] = {
                        NF_N_H_BDDF3_3D_F0_Eta,
                        NF_N_H_BDDF3_3D_F1_Eta,
                        NF_N_H_BDDF3_3D_F2_Eta,
                        NF_N_H_BDDF3_3D_F3_Eta,
                        NF_N_H_BDDF3_3D_F4_Eta,
                        NF_N_H_BDDF3_3D_F5_Eta };

static double *NF_N_H_BDDF3_3D_ZetaArray[6] = {
                        NF_N_H_BDDF3_3D_F0_Zeta,
                        NF_N_H_BDDF3_3D_F1_Zeta,
                        NF_N_H_BDDF3_3D_F2_Zeta,
                        NF_N_H_BDDF3_3D_F3_Zeta,
                        NF_N_H_BDDF3_3D_F4_Zeta,
                        NF_N_H_BDDF3_3D_F5_Zeta };

static double NF_N_H_BDDF3_3D_T[] = {-100};//???
static double NF_N_H_BDDF3_3D_S[] = {-100};//???



void NF_N_H_BDDF3_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  //face 0
  Functionals[0]  = -PointValues[248] * 4.0;
  Functionals[1]  = -PointValues[249] * 4.0;
  Functionals[2]  = -PointValues[250] * 4.0;
  Functionals[3]  = -PointValues[251] * 4.0;
  Functionals[4]  = -PointValues[252] * 4.0;
  Functionals[5]  = -PointValues[253] * 4.0;
  Functionals[6]  = -PointValues[254] * 4.0;
  Functionals[7]  = -PointValues[255] * 4.0;
  Functionals[8]  = -PointValues[256] * 4.0;
  Functionals[9]  = -PointValues[257] * 4.0;
  //face 1
  Functionals[10] = -PointValues[134] * 4.0;
  Functionals[11] = -PointValues[135] * 4.0;
  Functionals[12] = -PointValues[136] * 4.0;
  Functionals[13] = -PointValues[137] * 4.0;
  Functionals[14] = -PointValues[138] * 4.0;
  Functionals[15] = -PointValues[139] * 4.0;
  Functionals[16] = -PointValues[140] * 4.0;
  Functionals[17] = -PointValues[141] * 4.0;
  Functionals[18] = -PointValues[142] * 4.0;
  Functionals[19] = -PointValues[143] * 4.0;
  //face 2
  Functionals[20] =  PointValues[20]  * 4.0;
  Functionals[21] =  PointValues[21]  * 4.0;
  Functionals[22] =  PointValues[22]  * 4.0;
  Functionals[23] =  PointValues[23]  * 4.0;
  Functionals[24] =  PointValues[24]  * 4.0;
  Functionals[25] =  PointValues[25]  * 4.0;
  Functionals[26] =  PointValues[26]  * 4.0;
  Functionals[27] =  PointValues[27]  * 4.0;
  Functionals[28] =  PointValues[28]  * 4.0;
  Functionals[29] =  PointValues[29]  * 4.0;
  //face 3
  Functionals[30] =  PointValues[154] * 4.0;
  Functionals[31] =  PointValues[155] * 4.0;
  Functionals[32] =  PointValues[156] * 4.0;
  Functionals[33] =  PointValues[157] * 4.0;
  Functionals[34] =  PointValues[158] * 4.0;
  Functionals[35] =  PointValues[159] * 4.0;
  Functionals[36] =  PointValues[160] * 4.0;
  Functionals[37] =  PointValues[161] * 4.0;
  Functionals[38] =  PointValues[162] * 4.0;
  Functionals[39] =  PointValues[163] * 4.0;
  //face 4
  Functionals[40] = -PointValues[40]  * 4.0;
  Functionals[41] = -PointValues[41]  * 4.0;
  Functionals[42] = -PointValues[42]  * 4.0;
  Functionals[43] = -PointValues[43]  * 4.0;
  Functionals[44] = -PointValues[44]  * 4.0;
  Functionals[45] = -PointValues[45]  * 4.0;
  Functionals[46] = -PointValues[46]  * 4.0;
  Functionals[47] = -PointValues[47]  * 4.0;
  Functionals[48] = -PointValues[48]  * 4.0;
  Functionals[49] = -PointValues[49]  * 4.0;
  //face 5
  Functionals[50] =  PointValues[298] * 4.0;
  Functionals[51] =  PointValues[299] * 4.0;
  Functionals[52] =  PointValues[300] * 4.0;
  Functionals[53] =  PointValues[301] * 4.0;
  Functionals[54] =  PointValues[302] * 4.0;
  Functionals[55] =  PointValues[303] * 4.0;
  Functionals[56] =  PointValues[304] * 4.0;
  Functionals[57] =  PointValues[305] * 4.0;
  Functionals[58] =  PointValues[306] * 4.0;
  Functionals[59] =  PointValues[307] * 4.0;
  //inner dofs
  int i;
  double s;

  //x-component
  s = 0;
  //constant 1
  for(i=0;i<64;i++)
    s += PointValues[i+60] * NF_D_H_P3_3D_Weights[i];
  Functionals[60] = s;

  s = 0;
  //x
  for(i=0;i<64;i++)
    s += PointValues[i+60] * NF_D_H_P3_3D_Array2[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[61] = s;

  s = 0;
  //y
  for(i=0;i<64;i++)
    s += PointValues[i+60] * NF_D_H_P3_3D_Array3[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[62] = s;

  s = 0;
  ///z
  for(i=0;i<64;i++)
    s += PointValues[i+60] * NF_D_H_P3_3D_Array4[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[63] = s;

  //y-component
  s = 0;
  //constant 1
  for(i=0;i<64;i++)
    s += PointValues[i+184] * NF_D_H_P3_3D_Weights[i];
  Functionals[64] = s;

  s = 0;
  //x
  for(i=0;i<64;i++)
    s += PointValues[i+184] * NF_D_H_P3_3D_Array2[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[65] = s;

  s = 0;
  //y
  for(i=0;i<64;i++)
    s += PointValues[i+184] * NF_D_H_P3_3D_Array3[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[66] = s;

  s = 0;
  ///z
  for(i=0;i<64;i++)
    s += PointValues[i+184] * NF_D_H_P3_3D_Array4[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[67] = s;

  //z-component
  s = 0;
  //constant 1
  for(i=0;i<64;i++)
    s += PointValues[i+308] * NF_D_H_P3_3D_Weights[i];
  Functionals[68] = s;

  s = 0;
  //x
  for(i=0;i<64;i++)
    s += PointValues[i+308] * NF_D_H_P3_3D_Array2[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[69] = s;

  s = 0;
  //y
  for(i=0;i<64;i++)
    s += PointValues[i+308] * NF_D_H_P3_3D_Array3[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[70] = s;

  s = 0;
  ///z
  for(i=0;i<64;i++)
    s += PointValues[i+308] * NF_D_H_P3_3D_Array4[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[71] = s;
}

void NF_N_H_BDDF3_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
                           double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  #ifdef __3D__
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 4, length == {4,4,4,4}
  Cell->GetVertex(faceVertex[4*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[4*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[4*face + 2])->GetCoords(x2,y2,z2);
  #endif
  // compute measure of this face
  s = sqrt( POW((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + POW((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + POW((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  Functionals[0] = PointValues[0]*s;
  Functionals[1] = PointValues[1]*s;
  Functionals[2] = PointValues[2]*s;
  Functionals[3] = PointValues[3]*s;
  Functionals[4] = PointValues[4]*s;
  Functionals[5] = PointValues[5]*s;
  Functionals[6] = PointValues[6]*s;
  Functionals[7] = PointValues[7]*s;
  Functionals[8] = PointValues[8]*s;
  Functionals[9] = PointValues[9]*s;
}

static int NF_N_H_BDDF3_3D_N_AllFunctionals = 72;
static int NF_N_H_BDDF3_3D_N_PointsAll = 124;
static int NF_N_H_BDDF3_3D_N_FaceFunctionals[] = { 10, 10, 10, 10, 10, 10 };
static int NF_N_H_BDDF3_3D_N_PointsFace[] = { 10, 10, 10, 10, 10, 10 };

TNodalFunctional3D *NF_N_H_BDDF3_3D_Obj = new TNodalFunctional3D
        (NF_N_H_BDDF3_3D, NF_N_H_BDDF3_3D_N_AllFunctionals,
         NF_N_H_BDDF3_3D_N_FaceFunctionals, NF_N_H_BDDF3_3D_N_PointsAll,
         NF_N_H_BDDF3_3D_N_PointsFace,
         NF_N_H_BDDF3_3D_Xi, NF_N_H_BDDF3_3D_Eta, NF_N_H_BDDF3_3D_Zeta,
         NF_N_H_BDDF3_3D_XiArray, NF_N_H_BDDF3_3D_EtaArray,
         NF_N_H_BDDF3_3D_ZetaArray,
         NF_N_H_BDDF3_3D_T, NF_N_H_BDDF3_3D_S,
         NF_N_H_BDDF3_3D_EvalAll, NF_N_H_BDDF3_3D_EvalFace);
