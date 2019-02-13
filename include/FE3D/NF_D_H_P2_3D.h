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
   
/*
    TNodalFunctional3D(NodalFunctional3D id,
                       int n_allfunctionals, int *n_facefunctionals,
                       int n_pointsall, int *n_pointsface,
                       double *xi, double *eta, double *zeta,
                       double **xiarray, double **etaarray,
                       double **zetaarray,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evalface);
*/

/* for all functionals */
static double NF_D_H_P2_3D_Xi[]   = {
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530,
        -0.7745966692414833770358530, 0,
         0.7745966692414833770358530 };
static double NF_D_H_P2_3D_Eta[]  = {
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0,
         0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0,
         0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0,
         0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530 };
static double NF_D_H_P2_3D_Zeta[] = {
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, -0.7745966692414833770358530,
        -0.7745966692414833770358530, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530, 0.7745966692414833770358530,
         0.7745966692414833770358530 };

static double NF_D_H_P2_3D_Weight[] = {
        0.171467764060356652949245541838,
        0.274348422496570644718792866941,
        0.171467764060356652949245541838,
        0.274348422496570644718792866941,
        0.438957475994513031550068587106,
        0.274348422496570644718792866941,
        0.171467764060356652949245541838,
        0.274348422496570644718792866941,
        0.171467764060356652949245541838,
        0.274348422496570644718792866941,
        0.438957475994513031550068587106,
        0.274348422496570644718792866941,
        0.438957475994513031550068587106,
        0.702331961591220850480109739369,
        0.438957475994513031550068587106,
        0.274348422496570644718792866941,
        0.438957475994513031550068587106,
        0.274348422496570644718792866941,
        0.171467764060356652949245541838,
        0.274348422496570644718792866941,
        0.171467764060356652949245541838,
        0.274348422496570644718792866941,
        0.438957475994513031550068587106,
        0.274348422496570644718792866941,
        0.171467764060356652949245541838,
        0.274348422496570644718792866941,
        0.171467764060356652949245541838 };

static double NF_D_H_P2_3D_XiXi[] = {
        0.8, -1, 0.8,
        0.8, -1, 0.8,
        0.8, -1, 0.8,
        0.8, -1, 0.8,
        0.8, -1, 0.8,
        0.8, -1, 0.8,
        0.8, -1, 0.8,
        0.8, -1, 0.8,
        0.8, -1, 0.8 };

static double NF_D_H_P2_3D_XiEta[] = {
        0.6, 0, -0.6,
        0, 0, 0,
       -0.6, 0, 0.6,
        0.6, 0, -0.6,
        0, 0, 0,
       -0.6, 0, 0.6,
        0.6, 0, -0.6,
        0, 0, 0,
       -0.6, 0, 0.6 };


static double NF_D_H_P2_3D_XiZeta[] = {
        0.6, 0, -0.6,
        0.6, 0, -0.6,
        0.6, 0, -0.6,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
       -0.6, 0, 0.6,
       -0.6, 0, 0.6,
       -0.6, 0, 0.6 };

static double NF_D_H_P2_3D_EtaEta[] = {
        0.8, 0.8, 0.8,
       -1, -1, -1,
        0.8, 0.8, 0.8,
        0.8, 0.8, 0.8,
       -1, -1, -1,
        0.8, 0.8, 0.8,
        0.8, 0.8, 0.8,
       -1, -1, -1,
        0.8, 0.8, 0.8 };

static double NF_D_H_P2_3D_EtaZeta[] = {
        0.6, 0.6, 0.6, 0, 0, 0, -0.6, -0.6, -0.6,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
       -0.6, -0.6, -0.6, 0, 0, 0, 0.6, 0.6, 0.6 };

static double NF_D_H_P2_3D_ZetaZeta[] = {
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
       -1, -1, -1, -1, -1, -1, -1, -1, -1,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 };

/* face 0                               0 */
static double *NF_D_H_P2_3D_F0_Xi = NULL;
static double *NF_D_H_P2_3D_F0_Eta = NULL;
static double *NF_D_H_P2_3D_F0_Zeta = NULL;

/* face 1                               1 */
static double *NF_D_H_P2_3D_F1_Xi = NULL;
static double *NF_D_H_P2_3D_F1_Eta = NULL;
static double *NF_D_H_P2_3D_F1_Zeta = NULL;

/* face 2                               2 */
static double *NF_D_H_P2_3D_F2_Xi = NULL;
static double *NF_D_H_P2_3D_F2_Eta = NULL;
static double *NF_D_H_P2_3D_F2_Zeta = NULL;

/* face 3                               3 */
static double *NF_D_H_P2_3D_F3_Xi = NULL;
static double *NF_D_H_P2_3D_F3_Eta = NULL;
static double *NF_D_H_P2_3D_F3_Zeta = NULL;

/* face 4                               4 */
static double *NF_D_H_P2_3D_F4_Xi = NULL;
static double *NF_D_H_P2_3D_F4_Eta = NULL;
static double *NF_D_H_P2_3D_F4_Zeta = NULL;

/* face 5                               5 */
static double *NF_D_H_P2_3D_F5_Xi = NULL;
static double *NF_D_H_P2_3D_F5_Eta = NULL;
static double *NF_D_H_P2_3D_F5_Zeta = NULL;

static double *NF_D_H_P2_3D_XiArray[6] = { 
                        NF_D_H_P2_3D_F0_Xi,
                        NF_D_H_P2_3D_F1_Xi,
                        NF_D_H_P2_3D_F2_Xi,
                        NF_D_H_P2_3D_F3_Xi,
                        NF_D_H_P2_3D_F4_Xi,
                        NF_D_H_P2_3D_F5_Xi };

static double *NF_D_H_P2_3D_EtaArray[6] = { 
                        NF_D_H_P2_3D_F0_Eta,
                        NF_D_H_P2_3D_F1_Eta,
                        NF_D_H_P2_3D_F2_Eta,
                        NF_D_H_P2_3D_F3_Eta,
                        NF_D_H_P2_3D_F4_Eta,
                        NF_D_H_P2_3D_F5_Eta };

static double *NF_D_H_P2_3D_ZetaArray[6] = { 
                        NF_D_H_P2_3D_F0_Zeta,
                        NF_D_H_P2_3D_F1_Zeta,
                        NF_D_H_P2_3D_F2_Zeta,
                        NF_D_H_P2_3D_F3_Zeta,
                        NF_D_H_P2_3D_F4_Zeta,
                        NF_D_H_P2_3D_F5_Zeta };

static double *NF_D_H_P2_3D_T = NULL;
static double *NF_D_H_P2_3D_S = NULL;

void NF_D_H_P2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  Functionals[0] = 
    ( +PointValues[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_Weight[26] ) * 0.125;
  Functionals[1] = 
    ( +PointValues[0]*NF_D_H_P2_3D_Xi[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_Xi[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_Xi[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_Xi[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_Xi[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_Xi[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_Xi[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_Xi[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_Xi[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_Xi[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_Xi[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_Xi[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_Xi[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_Xi[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_Xi[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_Xi[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_Xi[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_Xi[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_Xi[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_Xi[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_Xi[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_Xi[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_Xi[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_Xi[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_Xi[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_Xi[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_Xi[26]*NF_D_H_P2_3D_Weight[26] ) * 0.375;
  Functionals[2] = 
    ( +PointValues[0]*NF_D_H_P2_3D_Eta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_Eta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_Eta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_Eta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_Eta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_Eta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_Eta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_Eta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_Eta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_Eta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_Eta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_Eta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_Eta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_Eta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_Eta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_Eta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_Eta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_Eta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_Eta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_Eta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_Eta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_Eta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_Eta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_Eta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_Eta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_Eta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_Eta[26]*NF_D_H_P2_3D_Weight[26] ) * 0.375;
  Functionals[3] = 
    ( +PointValues[0]*NF_D_H_P2_3D_Zeta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_Zeta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_Zeta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_Zeta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_Zeta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_Zeta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_Zeta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_Zeta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_Zeta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_Zeta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_Zeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_Zeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_Zeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_Zeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_Zeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_Zeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_Zeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_Zeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_Zeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_Zeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_Zeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_Zeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_Zeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_Zeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_Zeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_Zeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_Zeta[26]*NF_D_H_P2_3D_Weight[26] ) * 0.375;
  Functionals[4] = 
    ( +PointValues[0]*NF_D_H_P2_3D_XiXi[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_XiXi[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_XiXi[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_XiXi[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_XiXi[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_XiXi[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_XiXi[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_XiXi[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_XiXi[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_XiXi[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_XiXi[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_XiXi[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_XiXi[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_XiXi[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_XiXi[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_XiXi[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_XiXi[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_XiXi[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_XiXi[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_XiXi[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_XiXi[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_XiXi[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_XiXi[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_XiXi[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_XiXi[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_XiXi[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_XiXi[26]*NF_D_H_P2_3D_Weight[26] ) * 0.15625;
  Functionals[5] = 
    ( +PointValues[0]*NF_D_H_P2_3D_XiEta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_XiEta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_XiEta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_XiEta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_XiEta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_XiEta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_XiEta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_XiEta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_XiEta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_XiEta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_XiEta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_XiEta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_XiEta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_XiEta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_XiEta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_XiEta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_XiEta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_XiEta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_XiEta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_XiEta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_XiEta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_XiEta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_XiEta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_XiEta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_XiEta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_XiEta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_XiEta[26]*NF_D_H_P2_3D_Weight[26] ) * 1.125;
  Functionals[6] = 
    ( +PointValues[0]*NF_D_H_P2_3D_XiZeta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_XiZeta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_XiZeta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_XiZeta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_XiZeta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_XiZeta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_XiZeta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_XiZeta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_XiZeta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_XiZeta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_XiZeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_XiZeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_XiZeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_XiZeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_XiZeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_XiZeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_XiZeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_XiZeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_XiZeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_XiZeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_XiZeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_XiZeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_XiZeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_XiZeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_XiZeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_XiZeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_XiZeta[26]*NF_D_H_P2_3D_Weight[26] ) * 1.125;
  Functionals[7] = 
    ( +PointValues[0]*NF_D_H_P2_3D_EtaEta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_EtaEta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_EtaEta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_EtaEta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_EtaEta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_EtaEta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_EtaEta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_EtaEta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_EtaEta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_EtaEta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_EtaEta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_EtaEta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_EtaEta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_EtaEta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_EtaEta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_EtaEta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_EtaEta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_EtaEta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_EtaEta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_EtaEta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_EtaEta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_EtaEta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_EtaEta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_EtaEta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_EtaEta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_EtaEta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_EtaEta[26]*NF_D_H_P2_3D_Weight[26] ) * 0.15625;
  Functionals[8] = 
    ( +PointValues[0]*NF_D_H_P2_3D_EtaZeta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_EtaZeta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_EtaZeta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_EtaZeta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_EtaZeta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_EtaZeta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_EtaZeta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_EtaZeta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_EtaZeta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_EtaZeta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_EtaZeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_EtaZeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_EtaZeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_EtaZeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_EtaZeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_EtaZeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_EtaZeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_EtaZeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_EtaZeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_EtaZeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_EtaZeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_EtaZeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_EtaZeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_EtaZeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_EtaZeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_EtaZeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_EtaZeta[26]*NF_D_H_P2_3D_Weight[26] ) * 1.125;
  Functionals[9] = 
    ( +PointValues[0]*NF_D_H_P2_3D_ZetaZeta[0]*NF_D_H_P2_3D_Weight[0]
      +PointValues[1]*NF_D_H_P2_3D_ZetaZeta[1]*NF_D_H_P2_3D_Weight[1]
      +PointValues[2]*NF_D_H_P2_3D_ZetaZeta[2]*NF_D_H_P2_3D_Weight[2]
      +PointValues[3]*NF_D_H_P2_3D_ZetaZeta[3]*NF_D_H_P2_3D_Weight[3]
      +PointValues[4]*NF_D_H_P2_3D_ZetaZeta[4]*NF_D_H_P2_3D_Weight[4]
      +PointValues[5]*NF_D_H_P2_3D_ZetaZeta[5]*NF_D_H_P2_3D_Weight[5]
      +PointValues[6]*NF_D_H_P2_3D_ZetaZeta[6]*NF_D_H_P2_3D_Weight[6]
      +PointValues[7]*NF_D_H_P2_3D_ZetaZeta[7]*NF_D_H_P2_3D_Weight[7]
      +PointValues[8]*NF_D_H_P2_3D_ZetaZeta[8]*NF_D_H_P2_3D_Weight[8]
      +PointValues[9]*NF_D_H_P2_3D_ZetaZeta[9]*NF_D_H_P2_3D_Weight[9]
      +PointValues[10]*NF_D_H_P2_3D_ZetaZeta[10]*NF_D_H_P2_3D_Weight[10]
      +PointValues[11]*NF_D_H_P2_3D_ZetaZeta[11]*NF_D_H_P2_3D_Weight[11]
      +PointValues[12]*NF_D_H_P2_3D_ZetaZeta[12]*NF_D_H_P2_3D_Weight[12]
      +PointValues[13]*NF_D_H_P2_3D_ZetaZeta[13]*NF_D_H_P2_3D_Weight[13]
      +PointValues[14]*NF_D_H_P2_3D_ZetaZeta[14]*NF_D_H_P2_3D_Weight[14]
      +PointValues[15]*NF_D_H_P2_3D_ZetaZeta[15]*NF_D_H_P2_3D_Weight[15]
      +PointValues[16]*NF_D_H_P2_3D_ZetaZeta[16]*NF_D_H_P2_3D_Weight[16]
      +PointValues[17]*NF_D_H_P2_3D_ZetaZeta[17]*NF_D_H_P2_3D_Weight[17]
      +PointValues[18]*NF_D_H_P2_3D_ZetaZeta[18]*NF_D_H_P2_3D_Weight[18]
      +PointValues[19]*NF_D_H_P2_3D_ZetaZeta[19]*NF_D_H_P2_3D_Weight[19]
      +PointValues[20]*NF_D_H_P2_3D_ZetaZeta[20]*NF_D_H_P2_3D_Weight[20]
      +PointValues[21]*NF_D_H_P2_3D_ZetaZeta[21]*NF_D_H_P2_3D_Weight[21]
      +PointValues[22]*NF_D_H_P2_3D_ZetaZeta[22]*NF_D_H_P2_3D_Weight[22]
      +PointValues[23]*NF_D_H_P2_3D_ZetaZeta[23]*NF_D_H_P2_3D_Weight[23]
      +PointValues[24]*NF_D_H_P2_3D_ZetaZeta[24]*NF_D_H_P2_3D_Weight[24]
      +PointValues[25]*NF_D_H_P2_3D_ZetaZeta[25]*NF_D_H_P2_3D_Weight[25]
      +PointValues[26]*NF_D_H_P2_3D_ZetaZeta[26]*NF_D_H_P2_3D_Weight[26] ) * 0.15625;
}

void NF_D_H_P2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                           double *PointValues, double *Functionals)
{
}

static int NF_D_H_P2_3D_N_AllFunctionals = 10;
static int NF_D_H_P2_3D_N_PointsAll = 27;
static int NF_D_H_P2_3D_N_FaceFunctionals[] = { 0, 0, 0, 0, 0, 0 };
static int NF_D_H_P2_3D_N_PointsFace[] = { 0, 0, 0, 0, 0, 0 };

TNodalFunctional3D *NF_D_H_P2_3D_Obj = new TNodalFunctional3D
        (NF_D_H_P2_3D, NF_D_H_P2_3D_N_AllFunctionals,
         NF_D_H_P2_3D_N_FaceFunctionals, NF_D_H_P2_3D_N_PointsAll,
         NF_D_H_P2_3D_N_PointsFace,
         NF_D_H_P2_3D_Xi, NF_D_H_P2_3D_Eta, NF_D_H_P2_3D_Zeta,
         NF_D_H_P2_3D_XiArray, NF_D_H_P2_3D_EtaArray,
         NF_D_H_P2_3D_ZetaArray,
         NF_D_H_P2_3D_T, NF_D_H_P2_3D_S,
         NF_D_H_P2_3D_EvalAll, NF_D_H_P2_3D_EvalFace);
