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
// P2 element, discontinuous, 3D
// 
// Author:     Markus Wolff
//
// ***********************************************************************

/* for all functionals */
static double NF_D_T_P2_3D_Xi[]   = {
  0.25,
  0.0714285714285714285714285714286,
  0.785714285714285714285714285714,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.399403576166799204996102147462,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.100596423833200795003897852538,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462};
static double NF_D_T_P2_3D_Eta[]  = {
  0.25,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.785714285714285714285714285714,
  0.0714285714285714285714285714286,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538};
static double NF_D_T_P2_3D_Zeta[] = {
  0.25,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.0714285714285714285714285714286,
  0.785714285714285714285714285714,
  0.100596423833200795003897852538,
  0.399403576166799204996102147462,
  0.399403576166799204996102147462,
  0.399403576166799204996102147462,
  0.100596423833200795003897852538,
  0.100596423833200795003897852538 };

static double NF_D_T_P2_3D_Weights[]= {
   -0.0131555555555555555555555555556,
    0.00762222222222222222222222222222,
    0.00762222222222222222222222222222,
    0.00762222222222222222222222222222,
    0.00762222222222222222222222222222,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889,
    0.0248888888888888888888888888889
};

/* face 0                               0 */
static double *NF_D_T_P2_3D_F0_Xi = NULL;
static double *NF_D_T_P2_3D_F0_Eta = NULL;
static double *NF_D_T_P2_3D_F0_Zeta = NULL;

/* face 1                               1 */
static double *NF_D_T_P2_3D_F1_Xi = NULL;
static double *NF_D_T_P2_3D_F1_Eta = NULL;
static double *NF_D_T_P2_3D_F1_Zeta = NULL;

/* face 2                               2 */
static double *NF_D_T_P2_3D_F2_Xi = NULL;
static double *NF_D_T_P2_3D_F2_Eta = NULL;
static double *NF_D_T_P2_3D_F2_Zeta = NULL;

/* face 3                               3 */
static double *NF_D_T_P2_3D_F3_Xi = NULL;
static double *NF_D_T_P2_3D_F3_Eta = NULL;
static double *NF_D_T_P2_3D_F3_Zeta = NULL;

static double *NF_D_T_P2_3D_XiArray[4] = {
                        NF_D_T_P2_3D_F0_Xi,
                        NF_D_T_P2_3D_F1_Xi,
                        NF_D_T_P2_3D_F2_Xi,
                        NF_D_T_P2_3D_F3_Xi };

static double *NF_D_T_P2_3D_EtaArray[4] = {
                        NF_D_T_P2_3D_F0_Eta,
                        NF_D_T_P2_3D_F1_Eta,
                        NF_D_T_P2_3D_F2_Eta,
                        NF_D_T_P2_3D_F3_Eta };

static double *NF_D_T_P2_3D_ZetaArray[4] = {
                        NF_D_T_P2_3D_F0_Zeta,
                        NF_D_T_P2_3D_F1_Zeta,
                        NF_D_T_P2_3D_F2_Zeta,
                        NF_D_T_P2_3D_F3_Zeta };

static double *NF_D_T_P2_3D_T = NULL;
static double *NF_D_T_P2_3D_S = NULL;

void NF_D_T_P2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues, double *Functionals)
{
  int i;
  double s;

  s = 0;
  //constant
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[0] = s;

  s = 0;
  //x
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[1] = s;

  s = 0;
  //y
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[2] = s;

  s = 0;
  //z
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[3] = s;

  s = 0;
  //x*x
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[4] = s;

  s = 0;
  //x*y
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[5] = s;

  s = 0;
  //x*z
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Xi[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[6] = s;

  s = 0;
  //yy
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[7] = s;

  s = 0;
  //y*z
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Eta[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[8] = s;

  s = 0;
  //zz
  for(i=0;i<11;i++)
    s += PointValues[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Zeta[i] * NF_D_T_P2_3D_Weights[i];
  Functionals[9] = s;
}

void NF_D_T_P2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint,
			   double *PointValues, double *Functionals)
{
}

static int NF_D_T_P2_3D_N_AllFunctionals = 10;
static int NF_D_T_P2_3D_N_PointsAll = 11;
static int NF_D_T_P2_3D_N_FaceFunctionals[] = { 0, 0, 0, 0 };
static int NF_D_T_P2_3D_N_PointsFace[] = { 0, 0, 0, 0 };

TNodalFunctional3D *NF_D_T_P2_3D_Obj = new TNodalFunctional3D
        (NF_D_T_P2_3D, NF_D_T_P2_3D_N_AllFunctionals,
         NF_D_T_P2_3D_N_FaceFunctionals, NF_D_T_P2_3D_N_PointsAll,
         NF_D_T_P2_3D_N_PointsFace,
         NF_D_T_P2_3D_Xi, NF_D_T_P2_3D_Eta, NF_D_T_P2_3D_Zeta,
         NF_D_T_P2_3D_XiArray, NF_D_T_P2_3D_EtaArray,
         NF_D_T_P2_3D_ZetaArray,
         NF_D_T_P2_3D_T, NF_D_T_P2_3D_S,
         NF_D_T_P2_3D_EvalAll, NF_D_T_P2_3D_EvalFace);
