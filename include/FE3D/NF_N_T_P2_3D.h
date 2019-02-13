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
                       int B_T_Pointsall, int *B_T_Pointsface,
                       double *xi, double *eta, double *zeta,
                       double **xiarray, double **etaarray,
                       double **zetaarray,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evalface);
*/

/* for all functionals */
static double NF_N_T_P2_3D_Xi[39]   = { 
 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.05971587178976982, 0.4701420641051151,
 0.4701420641051151,

 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.05971587178976982,
 0.4701420641051151,

 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.05971587178976982, 0.4701420641051151,
 0.4701420641051151,

 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0,
 0.0,

 0.25,
 0.07142857142857143, 0.7857142857142857,
 0.07142857142857143, 0.07142857142857143,
 0.3994035761667992, 0.3994035761667992,
 0.1005964238332008, 0.1005964238332008,
 0.1005964238332008, 0.3994035761667992 };

static double NF_N_T_P2_3D_Eta[39]  = {
 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.5971587178976982E-1,
 0.4701420641051151,

 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0,
 0.0,

 0.3333333333333333, 0.1012865073234563, 0.1012865073234563,
 0.7974269853530873, 0.4701420641051151, 0.4701420641051151,
 0.5971587178976982E-1,

 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.5971587178976982E-1, 0.4701420641051151,
 0.4701420641051151,

 0.25,
 0.7142857142857143E-1, 0.7142857142857143E-1,
 0.7857142857142857, 0.7142857142857143E-1,
 0.3994035761667992, 0.1005964238332008,
 0.3994035761667992, 0.1005964238332008,
 0.3994035761667992, 0.1005964238332008 };

static double NF_N_T_P2_3D_Zeta[39] = {
 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0,
 0.0,

 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.5971587178976982E-1, 0.4701420641051151,
 0.4701420641051151,

 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.5971587178976982E-1,
 0.4701420641051151,

 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.5971587178976982E-1,
 0.4701420641051151,

 0.25,
 0.7142857142857143E-1, 0.7142857142857143E-1,
 0.7142857142857143E-1, 0.7857142857142857,
 0.1005964238332008, 0.3994035761667992,
 0.3994035761667992, 0.3994035761667992,
 0.1005964238332008, 0.1005964238332008 };

/* face 0                               0 */
static double NF_N_T_P2_3D_F0_Xi[7]   = {
 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.05971587178976982, 0.4701420641051151,
 0.4701420641051151 };
static double NF_N_T_P2_3D_F0_Eta[7]  = {
 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.5971587178976982E-1,
 0.4701420641051151 };
static double NF_N_T_P2_3D_F0_Zeta[7] = {
 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0,
 0.0 };

/* face 1                               1 */
static double NF_N_T_P2_3D_F1_Xi[7]   = {
 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.05971587178976982,
 0.4701420641051151 };
static double NF_N_T_P2_3D_F1_Eta[7]  = {
 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0,
 0.0 };
static double NF_N_T_P2_3D_F1_Zeta[7] = {
 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.5971587178976982E-1, 0.4701420641051151,
 0.4701420641051151 };

/* face 2                               2 */
static double NF_N_T_P2_3D_F2_Xi[7]   = {
 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.05971587178976982, 0.4701420641051151,
 0.4701420641051151 };
static double NF_N_T_P2_3D_F2_Eta[7]  = { 
 0.3333333333333333, 0.1012865073234563, 0.1012865073234563,
 0.7974269853530873, 0.4701420641051151, 0.4701420641051151,
 0.5971587178976982E-1 };
static double NF_N_T_P2_3D_F2_Zeta[7] = {
 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.5971587178976982E-1,
 0.4701420641051151 };

/* face 3                               3 */
static double NF_N_T_P2_3D_F3_Xi[7]   = {
 0.0, 0.0, 0.0,
 0.0, 0.0, 0.0,
 0.0 };
static double NF_N_T_P2_3D_F3_Eta[7]  = {
 0.3333333333333333, 0.7974269853530873, 0.1012865073234563,
 0.1012865073234563, 0.5971587178976982E-1, 0.4701420641051151,
 0.4701420641051151 };
static double NF_N_T_P2_3D_F3_Zeta[7] = {
 0.3333333333333333, 0.1012865073234563, 0.7974269853530873,
 0.1012865073234563, 0.4701420641051151, 0.5971587178976982E-1,
 0.4701420641051151 };

static double *NF_N_T_P2_3D_XiArray[4] = { 
                        NF_N_T_P2_3D_F0_Xi,
                        NF_N_T_P2_3D_F1_Xi,
                        NF_N_T_P2_3D_F2_Xi,
                        NF_N_T_P2_3D_F3_Xi };

static double *NF_N_T_P2_3D_EtaArray[4] = { 
                        NF_N_T_P2_3D_F0_Eta,
                        NF_N_T_P2_3D_F1_Eta,
                        NF_N_T_P2_3D_F2_Eta,
                        NF_N_T_P2_3D_F3_Eta };

static double *NF_N_T_P2_3D_ZetaArray[4] = { 
                        NF_N_T_P2_3D_F0_Zeta,
                        NF_N_T_P2_3D_F1_Zeta,
                        NF_N_T_P2_3D_F2_Zeta,
                        NF_N_T_P2_3D_F3_Zeta };


static double NF_N_T_P2_3D_T[7] =
                { 0.333333333333333333333333333333333,
                  0.797426985353087322398025276169754,
                  0.101286507323456338800987361915123,
                  0.101286507323456338800987361915123,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447, 
                  0.470142064105115089770441209513447 };

static double NF_N_T_P2_3D_S[7] =
                { 0.333333333333333333333333333333333, 
                  0.101286507323456338800987361915123,
                  0.797426985353087322398025276169754, 
                  0.101286507323456338800987361915123, 
                  0.470142064105115089770441209513447,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447 }; 

static double NF_N_T_P2_3D_W0[7] = {
      0.375E-1,
      0.6377969866281863E-2,
      0.6377969866281863E-2,
      0.5021365053984985E-1,
      0.3112203013371814E-1,
      0.3112203013371814E-1,
      0.3953016126816816E-2 };

static double NF_N_T_P2_3D_W1[7] = {
      0.375E-1,
      0.5021365053984985E-1,
      0.6377969866281863E-2,
      0.6377969866281863E-2,
      0.3953016126816816E-2,
      0.3112203013371814E-1,
      0.3112203013371814E-1 };

static double NF_N_T_P2_3D_W2[7] = {
      0.375E-1,
      0.6377969866281863E-2,
      0.5021365053984985E-1,
      0.6377969866281863E-2,
      0.3112203013371814E-1,
      0.3953016126816816E-2,
      0.3112203013371814E-1 };

static double NF_N_T_P2_3D_W[11] = {
                   -.0131555555555555555555555555556,
                    .00762222222222222222222222222222,
                    .00762222222222222222222222222222,
                    .00762222222222222222222222222222,
                    .00762222222222222222222222222222,
                    .0248888888888888888888888888889,
                    .0248888888888888888888888888889,
                    .0248888888888888888888888888889,
                    .0248888888888888888888888888889,
                    .0248888888888888888888888888889,
                    .0248888888888888888888888888889 };

void NF_N_T_P2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  Functionals[0] = ( NF_N_T_P2_3D_W0[0]*PointValues[0]
                    +NF_N_T_P2_3D_W0[1]*PointValues[1]
                    +NF_N_T_P2_3D_W0[2]*PointValues[2]
                    +NF_N_T_P2_3D_W0[3]*PointValues[3]
                    +NF_N_T_P2_3D_W0[4]*PointValues[4]
                    +NF_N_T_P2_3D_W0[5]*PointValues[5]
                    +NF_N_T_P2_3D_W0[6]*PointValues[6])*2;
  Functionals[1] = ( NF_N_T_P2_3D_W1[0]*PointValues[0]
                    +NF_N_T_P2_3D_W1[1]*PointValues[1]
                    +NF_N_T_P2_3D_W1[2]*PointValues[2]
                    +NF_N_T_P2_3D_W1[3]*PointValues[3]
                    +NF_N_T_P2_3D_W1[4]*PointValues[4]
                    +NF_N_T_P2_3D_W1[5]*PointValues[5]
                    +NF_N_T_P2_3D_W1[6]*PointValues[6])*2;
  Functionals[2] = ( NF_N_T_P2_3D_W2[0]*PointValues[0]
                    +NF_N_T_P2_3D_W2[1]*PointValues[1]
                    +NF_N_T_P2_3D_W2[2]*PointValues[2]
                    +NF_N_T_P2_3D_W2[3]*PointValues[3]
                    +NF_N_T_P2_3D_W2[4]*PointValues[4]
                    +NF_N_T_P2_3D_W2[5]*PointValues[5]
                    +NF_N_T_P2_3D_W2[6]*PointValues[6])*2;

  Functionals[3] = ( NF_N_T_P2_3D_W0[0]*PointValues[7]
                    +NF_N_T_P2_3D_W0[1]*PointValues[8]
                    +NF_N_T_P2_3D_W0[2]*PointValues[9]
                    +NF_N_T_P2_3D_W0[3]*PointValues[10]
                    +NF_N_T_P2_3D_W0[4]*PointValues[11]
                    +NF_N_T_P2_3D_W0[5]*PointValues[12]
                    +NF_N_T_P2_3D_W0[6]*PointValues[13])*2;
  Functionals[4] = ( NF_N_T_P2_3D_W1[0]*PointValues[7]
                    +NF_N_T_P2_3D_W1[1]*PointValues[8]
                    +NF_N_T_P2_3D_W1[2]*PointValues[9]
                    +NF_N_T_P2_3D_W1[3]*PointValues[10]
                    +NF_N_T_P2_3D_W1[4]*PointValues[11]
                    +NF_N_T_P2_3D_W1[5]*PointValues[12]
                    +NF_N_T_P2_3D_W1[6]*PointValues[13])*2;
  Functionals[5] = ( NF_N_T_P2_3D_W2[0]*PointValues[7]
                    +NF_N_T_P2_3D_W2[1]*PointValues[8]
                    +NF_N_T_P2_3D_W2[2]*PointValues[9]
                    +NF_N_T_P2_3D_W2[3]*PointValues[10]
                    +NF_N_T_P2_3D_W2[4]*PointValues[11]
                    +NF_N_T_P2_3D_W2[5]*PointValues[12]
                    +NF_N_T_P2_3D_W2[6]*PointValues[13])*2;

  Functionals[6] = ( NF_N_T_P2_3D_W0[0]*PointValues[14]
                    +NF_N_T_P2_3D_W0[1]*PointValues[15]
                    +NF_N_T_P2_3D_W0[2]*PointValues[16]
                    +NF_N_T_P2_3D_W0[3]*PointValues[17]
                    +NF_N_T_P2_3D_W0[4]*PointValues[18]
                    +NF_N_T_P2_3D_W0[5]*PointValues[19]
                    +NF_N_T_P2_3D_W0[6]*PointValues[20])*2;
  Functionals[7] = ( NF_N_T_P2_3D_W1[0]*PointValues[14]
                    +NF_N_T_P2_3D_W1[1]*PointValues[15]
                    +NF_N_T_P2_3D_W1[2]*PointValues[16]
                    +NF_N_T_P2_3D_W1[3]*PointValues[17]
                    +NF_N_T_P2_3D_W1[4]*PointValues[18]
                    +NF_N_T_P2_3D_W1[5]*PointValues[19]
                    +NF_N_T_P2_3D_W1[6]*PointValues[20])*2;
  Functionals[8] = ( NF_N_T_P2_3D_W2[0]*PointValues[14]
                    +NF_N_T_P2_3D_W2[1]*PointValues[15]
                    +NF_N_T_P2_3D_W2[2]*PointValues[16]
                    +NF_N_T_P2_3D_W2[3]*PointValues[17]
                    +NF_N_T_P2_3D_W2[4]*PointValues[18]
                    +NF_N_T_P2_3D_W2[5]*PointValues[19]
                    +NF_N_T_P2_3D_W2[6]*PointValues[20])*2;

  Functionals[9] = ( NF_N_T_P2_3D_W0[0]*PointValues[21]
                    +NF_N_T_P2_3D_W0[1]*PointValues[22]
                    +NF_N_T_P2_3D_W0[2]*PointValues[23]
                    +NF_N_T_P2_3D_W0[3]*PointValues[24]
                    +NF_N_T_P2_3D_W0[4]*PointValues[25]
                    +NF_N_T_P2_3D_W0[5]*PointValues[26]
                    +NF_N_T_P2_3D_W0[6]*PointValues[27])*2;
  Functionals[10] = ( NF_N_T_P2_3D_W1[0]*PointValues[21]
                     +NF_N_T_P2_3D_W1[1]*PointValues[22]
                     +NF_N_T_P2_3D_W1[2]*PointValues[23]
                     +NF_N_T_P2_3D_W1[3]*PointValues[24]
                     +NF_N_T_P2_3D_W1[4]*PointValues[25]
                     +NF_N_T_P2_3D_W1[5]*PointValues[26]
                     +NF_N_T_P2_3D_W1[6]*PointValues[27])*2;
  Functionals[11] = ( NF_N_T_P2_3D_W2[0]*PointValues[21]
                     +NF_N_T_P2_3D_W2[1]*PointValues[22]
                     +NF_N_T_P2_3D_W2[2]*PointValues[23]
                     +NF_N_T_P2_3D_W2[3]*PointValues[24]
                     +NF_N_T_P2_3D_W2[4]*PointValues[25]
                     +NF_N_T_P2_3D_W2[5]*PointValues[26]
                     +NF_N_T_P2_3D_W2[6]*PointValues[27])*2;

  Functionals[12] = ( NF_N_T_P2_3D_W[ 0]*PointValues[28]
                     +NF_N_T_P2_3D_W[ 1]*PointValues[29]
                     +NF_N_T_P2_3D_W[ 2]*PointValues[30]
                     +NF_N_T_P2_3D_W[ 3]*PointValues[31]
                     +NF_N_T_P2_3D_W[ 4]*PointValues[32]
                     +NF_N_T_P2_3D_W[ 5]*PointValues[33]
                     +NF_N_T_P2_3D_W[ 6]*PointValues[34]
                     +NF_N_T_P2_3D_W[ 7]*PointValues[35]
                     +NF_N_T_P2_3D_W[ 8]*PointValues[36]
                     +NF_N_T_P2_3D_W[ 9]*PointValues[37]
                     +NF_N_T_P2_3D_W[10]*PointValues[38])*6;
}

void NF_N_T_P2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                           double *PointValues, double *Functionals)
{
  Functionals[0] = ( NF_N_T_P2_3D_W0[0]*PointValues[0]
                    +NF_N_T_P2_3D_W0[1]*PointValues[1]
                    +NF_N_T_P2_3D_W0[2]*PointValues[2]
                    +NF_N_T_P2_3D_W0[3]*PointValues[3]
                    +NF_N_T_P2_3D_W0[4]*PointValues[4]
                    +NF_N_T_P2_3D_W0[5]*PointValues[5]
                    +NF_N_T_P2_3D_W0[6]*PointValues[6])*2;
  Functionals[1] = ( NF_N_T_P2_3D_W1[0]*PointValues[0]
                    +NF_N_T_P2_3D_W1[1]*PointValues[1]
                    +NF_N_T_P2_3D_W1[2]*PointValues[2]
                    +NF_N_T_P2_3D_W1[3]*PointValues[3]
                    +NF_N_T_P2_3D_W1[4]*PointValues[4]
                    +NF_N_T_P2_3D_W1[5]*PointValues[5]
                    +NF_N_T_P2_3D_W1[6]*PointValues[6])*2;
  Functionals[2] = ( NF_N_T_P2_3D_W2[0]*PointValues[0]
                    +NF_N_T_P2_3D_W2[1]*PointValues[1]
                    +NF_N_T_P2_3D_W2[2]*PointValues[2]
                    +NF_N_T_P2_3D_W2[3]*PointValues[3]
                    +NF_N_T_P2_3D_W2[4]*PointValues[4]
                    +NF_N_T_P2_3D_W2[5]*PointValues[5]
                    +NF_N_T_P2_3D_W2[6]*PointValues[6])*2;
}

static int NF_N_T_P2_3D_N_AllFunctionals = 13;
static int NF_N_T_P2_3D_N_T_PointsAll = 39;
static int NF_N_T_P2_3D_N_FaceFunctionals[] = { 3, 3, 3, 3 };
static int NF_N_T_P2_3D_N_T_PointsFace[] = { 7, 7, 7, 7 };

TNodalFunctional3D *NF_N_T_P2_3D_Obj = new TNodalFunctional3D
        (NF_N_T_P2_3D, NF_N_T_P2_3D_N_AllFunctionals,
         NF_N_T_P2_3D_N_FaceFunctionals, NF_N_T_P2_3D_N_T_PointsAll,
         NF_N_T_P2_3D_N_T_PointsFace,
         NF_N_T_P2_3D_Xi, NF_N_T_P2_3D_Eta, NF_N_T_P2_3D_Zeta,
         NF_N_T_P2_3D_XiArray, NF_N_T_P2_3D_EtaArray,
         NF_N_T_P2_3D_ZetaArray,
         NF_N_T_P2_3D_T, NF_N_T_P2_3D_S,
         NF_N_T_P2_3D_EvalAll, NF_N_T_P2_3D_EvalFace);
