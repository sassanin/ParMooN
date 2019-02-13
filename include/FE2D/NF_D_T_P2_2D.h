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
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_D_T_P2_2D_Xi[] = { 
                0.3333333333333333, 
                0.7974269853530873, 0.1012865073234563, 
                0.1012865073234563, 
                0.05971587178976982, 0.4701420641051151, 
                0.4701420641051151 };
static double NF_D_T_P2_2D_Eta[] = { 
                0.3333333333333333, 
                0.1012865073234563, 0.7974269853530873, 
                0.1012865073234563, 
                0.4701420641051151, 0.4701420641051151, 
                0.05971587178976982 };

/* weight : 1 */
static double NF_D_T_P2_2D_W1[] = {
                1.0000000000000000e+00,
                1.0000000000000000e+00,
                1.0000000000000000e+00,
                1.0000000000000000e+00,
                1.0000000000000000e+00,
                1.0000000000000000e+00,
                1.0000000000000000e+00 };
/* weight : 2*xi+eta-1 */
static double NF_D_T_P2_2D_W2[] = {
                0.0000000000000000e+00,
                6.9614047802963097e-01,
                0.0000000000000000e+00,
               -6.9614047802963108e-01,
               -4.1042619231534527e-01,
                4.1042619231534538e-01,
                0.0000000000000000e+00 };
/* weight : xi+2*eta-1 */
static double NF_D_T_P2_2D_W3[] = {
                0.0000000000000000e+00,
                0.0000000000000000e+00,
                6.9614047802963097e-01,
               -6.9614047802963108e-01,
                0.0000000000000000e+00,
                4.1042619231534538e-01,
               -4.1042619231534527e-01 };
/* weight : 180*xi^2-144*xi+18 */
static double NF_D_T_P2_2D_W4[] = {
               -1.0000000000000000e+01,
                1.7630677563631750e+01,
                5.2613551272635135e+00,
                5.2613551272635135e+00,
                1.0042791824123347e+01,
               -9.9144163517533173e+00,
               -9.9144163517533173e+00 };
/* weight : 360*xi*eta-72*xi-72*eta+18 */
static double NF_D_T_P2_2D_W5[] = {
                1.0000000000000000e+01,
               -1.7630677563631764e+01,
               -1.7630677563631764e+01,
                7.1079673091047368e+00,
               -1.0042791824123348e+01,
                2.9871624527629962e+01,
               -1.0042791824123348e+01 };
/* weight : 180*eta^2-144*eta+18 */
static double NF_D_T_P2_2D_W6[] = {
               -1.0000000000000000e+01,
                5.2613551272635135e+00,
                1.7630677563631750e+01,
                5.2613551272635135e+00,
               -9.9144163517533173e+00,
               -9.9144163517533173e+00,
                1.0042791824123347e+01 };

static double NF_D_T_P2_2D_T_P[] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

static double NF_D_T_P2_2D_W[7] = {
                0.1125, 
                0.0629695902724136, 0.0629695902724136,
                0.0629695902724136, 
                0.0661970763942531, 0.0661970763942531, 
                0.0661970763942531 };

void NF_D_T_P2_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0]*NF_D_T_P2_2D_W1[0]*NF_D_T_P2_2D_W[0]
                  +PointValues[1]*NF_D_T_P2_2D_W1[1]*NF_D_T_P2_2D_W[1]
                  +PointValues[2]*NF_D_T_P2_2D_W1[2]*NF_D_T_P2_2D_W[2]
                  +PointValues[3]*NF_D_T_P2_2D_W1[3]*NF_D_T_P2_2D_W[3]
                  +PointValues[4]*NF_D_T_P2_2D_W1[4]*NF_D_T_P2_2D_W[4]
                  +PointValues[5]*NF_D_T_P2_2D_W1[5]*NF_D_T_P2_2D_W[5]
                  +PointValues[6]*NF_D_T_P2_2D_W1[6]*NF_D_T_P2_2D_W[6];
  Functionals[0] *= 2;
  Functionals[1] = PointValues[0]*NF_D_T_P2_2D_W2[0]*NF_D_T_P2_2D_W[0]
                  +PointValues[1]*NF_D_T_P2_2D_W2[1]*NF_D_T_P2_2D_W[1]
                  +PointValues[2]*NF_D_T_P2_2D_W2[2]*NF_D_T_P2_2D_W[2]
                  +PointValues[3]*NF_D_T_P2_2D_W2[3]*NF_D_T_P2_2D_W[3]
                  +PointValues[4]*NF_D_T_P2_2D_W2[4]*NF_D_T_P2_2D_W[4]
                  +PointValues[5]*NF_D_T_P2_2D_W2[5]*NF_D_T_P2_2D_W[5]
                  +PointValues[6]*NF_D_T_P2_2D_W2[6]*NF_D_T_P2_2D_W[6];
  Functionals[2] = PointValues[0]*NF_D_T_P2_2D_W3[0]*NF_D_T_P2_2D_W[0]
                  +PointValues[1]*NF_D_T_P2_2D_W3[1]*NF_D_T_P2_2D_W[1]
                  +PointValues[2]*NF_D_T_P2_2D_W3[2]*NF_D_T_P2_2D_W[2]
                  +PointValues[3]*NF_D_T_P2_2D_W3[3]*NF_D_T_P2_2D_W[3]
                  +PointValues[4]*NF_D_T_P2_2D_W3[4]*NF_D_T_P2_2D_W[4]
                  +PointValues[5]*NF_D_T_P2_2D_W3[5]*NF_D_T_P2_2D_W[5]
                  +PointValues[6]*NF_D_T_P2_2D_W3[6]*NF_D_T_P2_2D_W[6];
  Functionals[3] = PointValues[0]*NF_D_T_P2_2D_W4[0]*NF_D_T_P2_2D_W[0]
                  +PointValues[1]*NF_D_T_P2_2D_W4[1]*NF_D_T_P2_2D_W[1]
                  +PointValues[2]*NF_D_T_P2_2D_W4[2]*NF_D_T_P2_2D_W[2]
                  +PointValues[3]*NF_D_T_P2_2D_W4[3]*NF_D_T_P2_2D_W[3]
                  +PointValues[4]*NF_D_T_P2_2D_W4[4]*NF_D_T_P2_2D_W[4]
                  +PointValues[5]*NF_D_T_P2_2D_W4[5]*NF_D_T_P2_2D_W[5]
                  +PointValues[6]*NF_D_T_P2_2D_W4[6]*NF_D_T_P2_2D_W[6];
  Functionals[4] = PointValues[0]*NF_D_T_P2_2D_W5[0]*NF_D_T_P2_2D_W[0]
                  +PointValues[1]*NF_D_T_P2_2D_W5[1]*NF_D_T_P2_2D_W[1]
                  +PointValues[2]*NF_D_T_P2_2D_W5[2]*NF_D_T_P2_2D_W[2]
                  +PointValues[3]*NF_D_T_P2_2D_W5[3]*NF_D_T_P2_2D_W[3]
                  +PointValues[4]*NF_D_T_P2_2D_W5[4]*NF_D_T_P2_2D_W[4]
                  +PointValues[5]*NF_D_T_P2_2D_W5[5]*NF_D_T_P2_2D_W[5]
                  +PointValues[6]*NF_D_T_P2_2D_W5[6]*NF_D_T_P2_2D_W[6];
  Functionals[5] = PointValues[0]*NF_D_T_P2_2D_W6[0]*NF_D_T_P2_2D_W[0]
                  +PointValues[1]*NF_D_T_P2_2D_W6[1]*NF_D_T_P2_2D_W[1]
                  +PointValues[2]*NF_D_T_P2_2D_W6[2]*NF_D_T_P2_2D_W[2]
                  +PointValues[3]*NF_D_T_P2_2D_W6[3]*NF_D_T_P2_2D_W[3]
                  +PointValues[4]*NF_D_T_P2_2D_W6[4]*NF_D_T_P2_2D_W[4]
                  +PointValues[5]*NF_D_T_P2_2D_W6[5]*NF_D_T_P2_2D_W[5]
                  +PointValues[6]*NF_D_T_P2_2D_W6[6]*NF_D_T_P2_2D_W[6];
}

void NF_D_T_P2_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

TNodalFunctional2D *NF_D_T_P2_2D_Obj = new TNodalFunctional2D
        (NF_D_T_P2_2D, 6, 0, 7, 0, NF_D_T_P2_2D_Xi, NF_D_T_P2_2D_Eta,
         NF_D_T_P2_2D_T_P, NF_D_T_P2_2D_EvalAll, NF_D_T_P2_2D_EvalEdge);
