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
// Raviart-Thomas element of third order on hexahedra, 3D
// ***********************************************************************

static double RT3T_tp = 0.16666666666666666667;

/* for all functionals */
static double NF_N_T_RT3_3D_Xi[]  = {
    RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp,
    RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp,
    RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp,
    0,0,0,0,0,0,0,0,0,0,
    //inner dof
    0.25,
    0.206829931610673204083980900024961,
    0.206829931610673204083980900024961,
    0.379510205167980387748057299925117,
    0.206829931610673204083980900024961,
    0.821035883105467230906058078714215e-1,
    0.821035883105467230906058078714215e-1,
    0.753689235068359830728182576385735,
    0.821035883105467230906058078714215e-1,
    0.578195050519799725317663886414270e-2,
    0.578195050519799725317663886414270e-2,
    0.982654148484406008240470083407572,
    0.578195050519799725317663886414270e-2,
    0.50532740018894224425624528557907e-1,
    0.449467259981105775574375471442092,
    0.50532740018894224425624528557907e-1,
    0.449467259981105775574375471442092,
    0.50532740018894224425624528557907e-1,
    0.449467259981105775574375471442092,
    0.229066536116811139600408854554753,
    0.356395827885340437169173969506114e-1,
    0.229066536116811139600408854554753,
    0.356395827885340437169173969506114e-1,
    0.229066536116811139600408854554753,
    0.229066536116811139600408854554753,
    0.506227344977843677082264893939883,
    0.506227344977843677082264893939883,
    0.506227344977843677082264893939883,
    0.229066536116811139600408854554753,
    0.229066536116811139600408854554753,
    0.356395827885340437169173969506114e-1,
    0.366077495531974236787738546327104e-1,
    0.190486041934633455699433285315099,
    0.366077495531974236787738546327104e-1,
    0.366077495531974236787738546327104e-1,
    0.366077495531974236787738546327104e-1,
    0.190486041934633455699433285315099,
    0.736298458958971696943019005419480,
    0.736298458958971696943019005419480,
    0.736298458958971696943019005419480,
    0.366077495531974236787738546327104e-1,
    0.366077495531974236787738546327104e-1,
    0.190486041934633455699433285315099
};
static double NF_N_T_RT3_3D_Eta[] = {
    RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp,
    0,0,0,0,0,0,0,0,0,0,
    4*RT3T_tp,0.5,2*RT3T_tp,RT3T_tp,0.5,2*RT3T_tp,RT3T_tp,2*RT3T_tp,RT3T_tp,RT3T_tp,
    RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp,
    //inner dof
    0.25,
    0.206829931610673204083980900024961,
    0.379510205167980387748057299925117,
    0.206829931610673204083980900024961,
    0.206829931610673204083980900024961,
    0.821035883105467230906058078714215e-1,
    0.753689235068359830728182576385735,
    0.821035883105467230906058078714215e-1,
    0.821035883105467230906058078714215e-1,
    0.578195050519799725317663886414270e-2,
    0.982654148484406008240470083407572,
    0.578195050519799725317663886414270e-2,
    0.578195050519799725317663886414270e-2,
    0.449467259981105775574375471442092,
    0.50532740018894224425624528557907e-1,
    0.50532740018894224425624528557907e-1,
    0.449467259981105775574375471442092,
    0.449467259981105775574375471442092,
    0.50532740018894224425624528557907e-1,
    0.356395827885340437169173969506114e-1,
    0.229066536116811139600408854554753,
    0.229066536116811139600408854554753,
    0.506227344977843677082264893939883,
    0.506227344977843677082264893939883,
    0.506227344977843677082264893939883,
    0.229066536116811139600408854554753,
    0.356395827885340437169173969506114e-1,
    0.229066536116811139600408854554753,
    0.229066536116811139600408854554753,
    0.356395827885340437169173969506114e-1,
    0.229066536116811139600408854554753,
    0.190486041934633455699433285315099,
    0.366077495531974236787738546327104e-1,
    0.366077495531974236787738546327104e-1,
    0.736298458958971696943019005419480,
    0.736298458958971696943019005419480,
    0.736298458958971696943019005419480,
    0.366077495531974236787738546327104e-1,
    0.190486041934633455699433285315099,
    0.366077495531974236787738546327104e-1,
    0.190486041934633455699433285315099,
    0.366077495531974236787738546327104e-1,
    0.366077495531974236787738546327104e-1
};
static double NF_N_T_RT3_3D_Zeta[]= {
    0,0,0,0,0,0,0,0,0,0,
    RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp,
    RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp,
    RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp,
    //inner dof
    0.25,
    0.379510205167980387748057299925117,
    0.206829931610673204083980900024961,
    0.206829931610673204083980900024961,
    0.206829931610673204083980900024961,
    0.753689235068359830728182576385737,
    0.82103588310546723090605807871423e-1,
    0.821035883105467230906058078714225e-1,
    0.821035883105467230906058078714225e-1,
    0.982654148484406008240470083407571,
    0.5781950505197997253176638864142e-2,
    0.578195050519799725317663886414230e-2,
    0.578195050519799725317663886414260e-2,
    0.449467259981105775574375471442094,
    0.449467259981105775574375471442094,
    0.449467259981105775574375471442094,
    0.50532740018894224425624528557909e-1,
    0.50532740018894224425624528557909e-1,
    0.50532740018894224425624528557909e-1,
    0.506227344977843677082264893939883,
    0.506227344977843677082264893939883,
    0.506227344977843677082264893939883,
    0.229066536116811139600408854554753,
    0.229066536116811139600408854554753,
    0.35639582788534043716917396950611e-1,
    0.35639582788534043716917396950611e-1,
    0.229066536116811139600408854554753,
    0.229066536116811139600408854554753,
    0.35639582788534043716917396950611e-1,
    0.229066536116811139600408854554753,
    0.229066536116811139600408854554753,
    0.736298458958971696943019005419481,
    0.736298458958971696943019005419481,
    0.736298458958971696943019005419481,
    0.190486041934633455699433285315100,
    0.36607749553197423678773854632711e-1,
    0.36607749553197423678773854632711e-1,
    0.190486041934633455699433285315100,
    0.36607749553197423678773854632711e-1,
    0.366077495531974236787738546327106e-1,
    0.36607749553197423678773854632711e-1,
    0.190486041934633455699433285315100,
    0.366077495531974236787738546327106e-1
};

static double NF_N_T_RT3_3D_Weights[]= {
    -0.205001886586399158405865177642941e-1,
     0.142503058228669012484397415358704e-1,
     0.142503058228669012484397415358704e-1,
     0.142503058228669012484397415358704e-1,
     0.142503058228669012484397415358704e-1,
     0.196703331313390098756280342445466e-2,
     0.196703331313390098756280342445466e-2,
     0.196703331313390098756280342445466e-2,
     0.196703331313390098756280342445466e-2,
     0.169834109092887379837744566704016e-3,
     0.169834109092887379837744566704016e-3,
     0.169834109092887379837744566704016e-3,
     0.169834109092887379837744566704016e-3,
     0.457968382446728180074351446297276e-2,
     0.457968382446728180074351446297276e-2,
     0.457968382446728180074351446297276e-2,
     0.457968382446728180074351446297276e-2,
     0.457968382446728180074351446297276e-2,
     0.457968382446728180074351446297276e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.570448580868191850680255862783040e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2,
     0.214051914116209259648335300092023e-2
};

/* face 0                               0 */
static double NF_N_T_RT3_3D_F0_Xi[]   = {RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp};
static double NF_N_T_RT3_3D_F0_Eta[]  = {RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp};
static double NF_N_T_RT3_3D_F0_Zeta[] = {0,0,0,0,0,0,0,0,0,0};

/* face 1                               1 */
static double NF_N_T_RT3_3D_F1_Xi[]   = {RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp};
static double NF_N_T_RT3_3D_F1_Eta[]  = {0,0,0,0,0,0,0,0,0,0};
static double NF_N_T_RT3_3D_F1_Zeta[] = {RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp};

/* face 2                               2 */
static double NF_N_T_RT3_3D_F2_Xi[]   = {RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp};
static double NF_N_T_RT3_3D_F2_Eta[]  = {4*RT3T_tp,0.5,2*RT3T_tp,RT3T_tp,0.5,2*RT3T_tp,RT3T_tp,2*RT3T_tp,RT3T_tp,RT3T_tp};
static double NF_N_T_RT3_3D_F2_Zeta[] = {RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp};

/* face 3                               3 */
static double NF_N_T_RT3_3D_F3_Xi[]   = {0,0,0,0,0,0,0,0,0,0};
static double NF_N_T_RT3_3D_F3_Eta[]  = {RT3T_tp,2*RT3T_tp,0.5,4*RT3T_tp,RT3T_tp,2*RT3T_tp,0.5,RT3T_tp,2*RT3T_tp,RT3T_tp};
static double NF_N_T_RT3_3D_F3_Zeta[] = {RT3T_tp,RT3T_tp,RT3T_tp,RT3T_tp,2*RT3T_tp,2*RT3T_tp,2*RT3T_tp,0.5,0.5,4*RT3T_tp};

static double *NF_N_T_RT3_3D_XiArray[4] = {
                        NF_N_T_RT3_3D_F0_Xi,
                        NF_N_T_RT3_3D_F1_Xi,
                        NF_N_T_RT3_3D_F2_Xi,
                        NF_N_T_RT3_3D_F3_Xi };

static double *NF_N_T_RT3_3D_EtaArray[4] = {
                        NF_N_T_RT3_3D_F0_Eta,
                        NF_N_T_RT3_3D_F1_Eta,
                        NF_N_T_RT3_3D_F2_Eta,
                        NF_N_T_RT3_3D_F3_Eta };

static double *NF_N_T_RT3_3D_ZetaArray[4] = {
                        NF_N_T_RT3_3D_F0_Zeta,
                        NF_N_T_RT3_3D_F1_Zeta,
                        NF_N_T_RT3_3D_F2_Zeta,
                        NF_N_T_RT3_3D_F3_Zeta };

static double NF_N_T_RT3_3D_T[1] = {};// ???
static double NF_N_T_RT3_3D_S[1] = {};// ???

void NF_N_T_RT3_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  // PointValues[4*i + j] means i-th component (i=0 for x, i=1 for y, i=2 for z)
  // at j-th evaluation point (see NF_N_T_RT3_3D_Xi, ...Eta, ...Zeta)
  //face 0
  Functionals[0]  = -PointValues[166];
  Functionals[1]  = -PointValues[167];
  Functionals[2]  = -PointValues[168];
  Functionals[3]  = -PointValues[169];
  Functionals[4]  = -PointValues[170];
  Functionals[5]  = -PointValues[171];
  Functionals[6]  = -PointValues[172];
  Functionals[7]  = -PointValues[173];
  Functionals[8]  = -PointValues[174];
  Functionals[9]  = -PointValues[175];
  //face 1
  Functionals[10] = -PointValues[93];
  Functionals[11] = -PointValues[94];
  Functionals[12] = -PointValues[95];
  Functionals[13] = -PointValues[96];
  Functionals[14] = -PointValues[97];
  Functionals[15] = -PointValues[98];
  Functionals[16] = -PointValues[99];
  Functionals[17] = -PointValues[100];
  Functionals[18] = -PointValues[101];
  Functionals[19] = -PointValues[102];
  //face 2
  Functionals[20] =  PointValues[20]+PointValues[103]+PointValues[186];
  Functionals[21] =  PointValues[21]+PointValues[104]+PointValues[187];
  Functionals[22] =  PointValues[22]+PointValues[105]+PointValues[188];
  Functionals[23] =  PointValues[23]+PointValues[106]+PointValues[189];
  Functionals[24] =  PointValues[24]+PointValues[107]+PointValues[190];
  Functionals[25] =  PointValues[25]+PointValues[108]+PointValues[191];
  Functionals[26] =  PointValues[26]+PointValues[109]+PointValues[192];
  Functionals[27] =  PointValues[27]+PointValues[110]+PointValues[193];
  Functionals[28] =  PointValues[28]+PointValues[111]+PointValues[194];
  Functionals[29] =  PointValues[29]+PointValues[112]+PointValues[195];
  //face 3
  Functionals[30] = -PointValues[30];
  Functionals[31] = -PointValues[31];
  Functionals[32] = -PointValues[32];
  Functionals[33] = -PointValues[33];
  Functionals[34] = -PointValues[34];
  Functionals[35] = -PointValues[35];
  Functionals[36] = -PointValues[36];
  Functionals[37] = -PointValues[37];
  Functionals[38] = -PointValues[38];
  Functionals[39] = -PointValues[39];

  //inner dofs
  int i;
  double s;

  //x-component
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[40] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Xi[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[41] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[42] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[43] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Xi[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[44] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[45] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[46] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Eta[i+40]*NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[47] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Eta[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[48] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+40] * NF_N_T_RT3_3D_Zeta[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[49] = s;

  //y-component
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Weights[i];
  Functionals[50] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Xi[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[51] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[52] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[53] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Xi[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[54] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[55] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[56] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Eta[i+40]*NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[57] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Eta[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[58] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+123] * NF_N_T_RT3_3D_Zeta[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[59] = s;

  //z-component
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Weights[i];
  Functionals[60] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Xi[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[61] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[62] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[63] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Xi[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[64] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[65] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Xi[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[66] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Eta[i+40]*NF_N_T_RT3_3D_Eta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[67] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Eta[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[68] = s;
  s = 0;
  for(i=0;i<43;i++)
    s += PointValues[i+206] * NF_N_T_RT3_3D_Zeta[i+40]*NF_N_T_RT3_3D_Zeta[i+40] * NF_N_T_RT3_3D_Weights[i];
  Functionals[69] = s;
}

void NF_N_T_RT3_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
                           double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  #ifdef __3D__
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 3, length == {3,3,3,3}
  Cell->GetVertex(faceVertex[3*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[3*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[3*face + 2])->GetCoords(x2,y2,z2);
  #endif
  // compute measure of this face
  s = sqrt( POW((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + POW((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + POW((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  for(int i=0; i<10; i++)
    Functionals[i] = PointValues[i]*s;
}

static int NF_N_T_RT3_3D_N_AllFunctionals = 70;
static int NF_N_T_RT3_3D_N_PointsAll = 83;
static int NF_N_T_RT3_3D_N_FaceFunctionals[] = {10,10,10,10};
static int NF_N_T_RT3_3D_N_PointsFace[] = {10,10,10,10};

TNodalFunctional3D *NF_N_T_RT3_3D_Obj = new TNodalFunctional3D
        (NF_N_T_RT3_3D, NF_N_T_RT3_3D_N_AllFunctionals,
         NF_N_T_RT3_3D_N_FaceFunctionals, NF_N_T_RT3_3D_N_PointsAll,
         NF_N_T_RT3_3D_N_PointsFace,
         NF_N_T_RT3_3D_Xi, NF_N_T_RT3_3D_Eta, NF_N_T_RT3_3D_Zeta,
         NF_N_T_RT3_3D_XiArray, NF_N_T_RT3_3D_EtaArray,
         NF_N_T_RT3_3D_ZetaArray,
         NF_N_T_RT3_3D_T, NF_N_T_RT3_3D_S,
         NF_N_T_RT3_3D_EvalAll, NF_N_T_RT3_3D_EvalFace);
