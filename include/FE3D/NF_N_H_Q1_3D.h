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

#define __POINTVALUES__

#ifdef __POINTVALUES__
/* for all functionals */
static double NF_N_H_Q1_3D_Xi[]   = {  0,  0,  1,  0, -1,  0 };
static double NF_N_H_Q1_3D_Eta[]  = {  0, -1,  0,  1,  0,  0 };
static double NF_N_H_Q1_3D_Zeta[] = { -1,  0,  0,  0,  0,  1 };

/* face 0                               0 */
static double NF_N_H_Q1_3D_F0_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F0_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F0_Zeta[] = { -1 };

/* face 1                               1 */
static double NF_N_H_Q1_3D_F1_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F1_Eta[]  = { -1 };
static double NF_N_H_Q1_3D_F1_Zeta[] = {  0 };

/* face 2                               2 */
static double NF_N_H_Q1_3D_F2_Xi[]   = {  1 };
static double NF_N_H_Q1_3D_F2_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F2_Zeta[] = {  0 };

/* face 3                               3 */
static double NF_N_H_Q1_3D_F3_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F3_Eta[]  = {  1 };
static double NF_N_H_Q1_3D_F3_Zeta[] = {  0 };

/* face 4                               4 */
static double NF_N_H_Q1_3D_F4_Xi[]   = { -1 };
static double NF_N_H_Q1_3D_F4_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F4_Zeta[] = {  0 };

/* face 5                               5 */
static double NF_N_H_Q1_3D_F5_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F5_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F5_Zeta[] = {  1 };

static double *NF_N_H_Q1_3D_XiArray[6] = { 
                        NF_N_H_Q1_3D_F0_Xi,
                        NF_N_H_Q1_3D_F1_Xi,
                        NF_N_H_Q1_3D_F2_Xi,
                        NF_N_H_Q1_3D_F3_Xi,
                        NF_N_H_Q1_3D_F4_Xi,
                        NF_N_H_Q1_3D_F5_Xi };

static double *NF_N_H_Q1_3D_EtaArray[6] = { 
                        NF_N_H_Q1_3D_F0_Eta,
                        NF_N_H_Q1_3D_F1_Eta,
                        NF_N_H_Q1_3D_F2_Eta,
                        NF_N_H_Q1_3D_F3_Eta,
                        NF_N_H_Q1_3D_F4_Eta,
                        NF_N_H_Q1_3D_F5_Eta };

static double *NF_N_H_Q1_3D_ZetaArray[6] = { 
                        NF_N_H_Q1_3D_F0_Zeta,
                        NF_N_H_Q1_3D_F1_Zeta,
                        NF_N_H_Q1_3D_F2_Zeta,
                        NF_N_H_Q1_3D_F3_Zeta,
                        NF_N_H_Q1_3D_F4_Zeta,
                        NF_N_H_Q1_3D_F5_Zeta };

static double NF_N_H_Q1_3D_T[1] = { 0.5 };
static double NF_N_H_Q1_3D_S[1] = { 0.5 };

void NF_N_H_Q1_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
}

void NF_N_H_Q1_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                        double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

static int NF_N_H_Q1_3D_N_AllFunctionals = 6;
static int NF_N_H_Q1_3D_N_PointsAll = 6;
static int NF_N_H_Q1_3D_N_FaceFunctionals[] = { 1, 1, 1, 1, 1, 1 };
static int NF_N_H_Q1_3D_N_PointsFace[] = { 1, 1, 1, 1, 1, 1 };

#else

/* for all functionals */
static double NF_N_H_Q1_3D_Xi[]   = { 
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
 -0.0000000000000000e+00,
 -0.0000000000000000e+00,
 -0.0000000000000000e+00,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  7.7459666924148340e-01 };

static double NF_N_H_Q1_3D_Eta[]  = {
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01 };

static double NF_N_H_Q1_3D_Zeta[] = { 
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -1.0000000000000000e+00,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
 -7.7459666924148340e-01,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  0.0000000000000000e+00,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  7.7459666924148340e-01,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00,
  1.0000000000000000e+00};

/* face 0                               0 */
static double NF_N_H_Q1_3D_F0_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F0_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F0_Zeta[] = { -1 };

/* face 1                               1 */
static double NF_N_H_Q1_3D_F1_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F1_Eta[]  = { -1 };
static double NF_N_H_Q1_3D_F1_Zeta[] = {  0 };

/* face 2                               2 */
static double NF_N_H_Q1_3D_F2_Xi[]   = {  1 };
static double NF_N_H_Q1_3D_F2_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F2_Zeta[] = {  0 };

/* face 3                               3 */
static double NF_N_H_Q1_3D_F3_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F3_Eta[]  = {  1 };
static double NF_N_H_Q1_3D_F3_Zeta[] = {  0 };

/* face 4                               4 */
static double NF_N_H_Q1_3D_F4_Xi[]   = { -1 };
static double NF_N_H_Q1_3D_F4_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F4_Zeta[] = {  0 };

/* face 5                               5 */
static double NF_N_H_Q1_3D_F5_Xi[]   = {  0 };
static double NF_N_H_Q1_3D_F5_Eta[]  = {  0 };
static double NF_N_H_Q1_3D_F5_Zeta[] = {  1 };

static double *NF_N_H_Q1_3D_XiArray[6] = { 
                        NF_N_H_Q1_3D_F0_Xi,
                        NF_N_H_Q1_3D_F1_Xi,
                        NF_N_H_Q1_3D_F2_Xi,
                        NF_N_H_Q1_3D_F3_Xi,
                        NF_N_H_Q1_3D_F4_Xi,
                        NF_N_H_Q1_3D_F5_Xi };

static double *NF_N_H_Q1_3D_EtaArray[6] = { 
                        NF_N_H_Q1_3D_F0_Eta,
                        NF_N_H_Q1_3D_F1_Eta,
                        NF_N_H_Q1_3D_F2_Eta,
                        NF_N_H_Q1_3D_F3_Eta,
                        NF_N_H_Q1_3D_F4_Eta,
                        NF_N_H_Q1_3D_F5_Eta };

static double *NF_N_H_Q1_3D_ZetaArray[6] = { 
                        NF_N_H_Q1_3D_F0_Zeta,
                        NF_N_H_Q1_3D_F1_Zeta,
                        NF_N_H_Q1_3D_F2_Zeta,
                        NF_N_H_Q1_3D_F3_Zeta,
                        NF_N_H_Q1_3D_F4_Zeta,
                        NF_N_H_Q1_3D_F5_Zeta };

static double NF_N_H_Q1_3D_T[9] = { 0.11270166537925831149, 0.5,
                                    0.88729833462074168851,
                                    0.11270166537925831149, 0.5,
                                    0.88729833462074168851,
                                    0.11270166537925831149, 0.5,
                                    0.88729833462074168851 };
static double NF_N_H_Q1_3D_S[9] = { 0.11270166537925831149,
                                    0.11270166537925831149,
                                    0.11270166537925831149,
                                    0.5, 0.5, 0.5,
                                    0.88729833462074168851,
                                    0.88729833462074168851,
                                    0.88729833462074168851 };

void NF_N_H_Q1_3D_EvalAll(TCollection *Coll, TBaseCell *cell,
                           double *PointValues, double *Functionals)
{
  Functionals[0] = ( 25*PointValues[0]
                    +40*PointValues[1]
                    +25*PointValues[2]
                    +40*PointValues[3]
                    +64*PointValues[4]
                    +40*PointValues[5]
                    +25*PointValues[6]
                    +40*PointValues[7]
                    +25*PointValues[8])/324;
  Functionals[1] = ( 25*PointValues[9]
                    +40*PointValues[10]
                    +25*PointValues[11]
                    +40*PointValues[12]
                    +64*PointValues[13]
                    +40*PointValues[14]
                    +25*PointValues[15]
                    +40*PointValues[16]
                    +25*PointValues[17])/324;
  Functionals[2] = ( 25*PointValues[18]
                    +40*PointValues[19]
                    +25*PointValues[20]
                    +40*PointValues[21]
                    +64*PointValues[22]
                    +40*PointValues[23]
                    +25*PointValues[24]
                    +40*PointValues[25]
                    +25*PointValues[26])/324;
  Functionals[3] = ( 25*PointValues[27]
                    +40*PointValues[28]
                    +25*PointValues[29]
                    +40*PointValues[30]
                    +64*PointValues[31]
                    +40*PointValues[32]
                    +25*PointValues[33]
                    +40*PointValues[34]
                    +25*PointValues[35])/324;
  Functionals[4] = ( 25*PointValues[36]
                    +40*PointValues[37]
                    +25*PointValues[38]
                    +40*PointValues[39]
                    +64*PointValues[40]
                    +40*PointValues[41]
                    +25*PointValues[42]
                    +40*PointValues[43]
                    +25*PointValues[44])/324;
  Functionals[5] = ( 25*PointValues[45]
                    +40*PointValues[46]
                    +25*PointValues[47]
                    +40*PointValues[48]
                    +64*PointValues[49]
                    +40*PointValues[50]
                    +25*PointValues[51]
                    +40*PointValues[52]
                    +25*PointValues[53])/324;
}

void NF_N_H_Q1_3D_EvalFace(TCollection *Coll, TBaseCell *cell, int Joint,
                           double *PointValues, double *Functionals)
{
  Functionals[0] = ( 25*PointValues[0]
                    +40*PointValues[1]
                    +25*PointValues[2]
                    +40*PointValues[3]
                    +64*PointValues[4]
                    +40*PointValues[5]
                    +25*PointValues[6]
                    +40*PointValues[7]
                    +25*PointValues[8])/324;
}

static int NF_N_H_Q1_3D_N_AllFunctionals = 6;
static int NF_N_H_Q1_3D_N_PointsAll = 54;
static int NF_N_H_Q1_3D_N_FaceFunctionals[] = { 1, 1, 1, 1, 1, 1 };
static int NF_N_H_Q1_3D_N_PointsFace[] = { 9, 9, 9, 9, 9, 9 };
#endif // __POINTVALUES__

TNodalFunctional3D *NF_N_H_Q1_3D_Obj = new TNodalFunctional3D
        (NF_N_H_Q1_3D, NF_N_H_Q1_3D_N_AllFunctionals,
         NF_N_H_Q1_3D_N_FaceFunctionals, NF_N_H_Q1_3D_N_PointsAll,
         NF_N_H_Q1_3D_N_PointsFace,
         NF_N_H_Q1_3D_Xi, NF_N_H_Q1_3D_Eta, NF_N_H_Q1_3D_Zeta,
         NF_N_H_Q1_3D_XiArray, NF_N_H_Q1_3D_EtaArray,
         NF_N_H_Q1_3D_ZetaArray,
         NF_N_H_Q1_3D_T, NF_N_H_Q1_3D_S,
         NF_N_H_Q1_3D_EvalAll, NF_N_H_Q1_3D_EvalFace);
