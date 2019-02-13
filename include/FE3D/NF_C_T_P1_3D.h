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
static double NF_C_T_P1_3D_Xi[]   = { 0, 1, 0, 0 };
static double NF_C_T_P1_3D_Eta[]  = { 0, 0, 1, 0 };
static double NF_C_T_P1_3D_Zeta[] = { 0, 0, 0, 1 };

/* face 0                               0   1   2 */
static double NF_C_T_P1_3D_F0_Xi[]   = {  0,  1,  0 };
static double NF_C_T_P1_3D_F0_Eta[]  = {  0,  0,  1 };
static double NF_C_T_P1_3D_F0_Zeta[] = {  0,  0,  0 };

/* face 1                               0   3   1 */
static double NF_C_T_P1_3D_F1_Xi[]   = {  0,  0,  1 };
static double NF_C_T_P1_3D_F1_Eta[]  = {  0,  0,  0 };
static double NF_C_T_P1_3D_F1_Zeta[] = {  0,  1,  0 };

/* face 2                               2   1   3 */
static double NF_C_T_P1_3D_F2_Xi[]   = {  0,  1,  0 };
static double NF_C_T_P1_3D_F2_Eta[]  = {  1,  0,  0 };
static double NF_C_T_P1_3D_F2_Zeta[] = {  0,  0,  1 };

/* face 3                               0   2   3 */
static double NF_C_T_P1_3D_F3_Xi[]   = {  0,  0,  0 };
static double NF_C_T_P1_3D_F3_Eta[]  = {  0,  1,  0 };
static double NF_C_T_P1_3D_F3_Zeta[] = {  0,  0,  1 };

static double *NF_C_T_P1_3D_XiArray[4] = { 
                        NF_C_T_P1_3D_F0_Xi,
                        NF_C_T_P1_3D_F1_Xi,
                        NF_C_T_P1_3D_F2_Xi,
                        NF_C_T_P1_3D_F3_Xi };

static double *NF_C_T_P1_3D_EtaArray[4] = { 
                        NF_C_T_P1_3D_F0_Eta,
                        NF_C_T_P1_3D_F1_Eta,
                        NF_C_T_P1_3D_F2_Eta,
                        NF_C_T_P1_3D_F3_Eta };

static double *NF_C_T_P1_3D_ZetaArray[4] = { 
                        NF_C_T_P1_3D_F0_Zeta,
                        NF_C_T_P1_3D_F1_Zeta,
                        NF_C_T_P1_3D_F2_Zeta,
                        NF_C_T_P1_3D_F3_Zeta };

static double NF_C_T_P1_3D_T[3] = { 0, 1, 0 };
static double NF_C_T_P1_3D_S[3] = { 0, 0, 1 };

void NF_C_T_P1_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
}

void NF_C_T_P1_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                           double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}

static int NF_C_T_P1_3D_N_AllFunctionals = 4;
static int NF_C_T_P1_3D_N_PointsAll = 4;
static int NF_C_T_P1_3D_N_FaceFunctionals[] = { 3, 3, 3, 3 };
static int NF_C_T_P1_3D_N_PointsFace[] = { 3, 3, 3, 3 };

TNodalFunctional3D *NF_C_T_P1_3D_Obj = new TNodalFunctional3D
        (NF_C_T_P1_3D, NF_C_T_P1_3D_N_AllFunctionals,
         NF_C_T_P1_3D_N_FaceFunctionals, NF_C_T_P1_3D_N_PointsAll,
         NF_C_T_P1_3D_N_PointsFace,
         NF_C_T_P1_3D_Xi, NF_C_T_P1_3D_Eta, NF_C_T_P1_3D_Zeta,
         NF_C_T_P1_3D_XiArray, NF_C_T_P1_3D_EtaArray,
         NF_C_T_P1_3D_ZetaArray,
         NF_C_T_P1_3D_T, NF_C_T_P1_3D_S,
         NF_C_T_P1_3D_EvalAll, NF_C_T_P1_3D_EvalFace);
