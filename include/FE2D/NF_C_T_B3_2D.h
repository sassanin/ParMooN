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

static double NF_C_T_B3_2D_Xi[] = {
                0.0000000000000000, 0.3333333333333333,
                0.6666666666666667, 1.0000000000000000,
                0.0000000000000000, 0.6666666666666667,
                0.0000000000000000, 0.3333333333333333,
                0.0000000000000000,
                0.3333333333333333, 
                0.7974269853530873, 0.1012865073234563, 
                0.1012865073234563, 
                0.05971587178976982, 0.4701420641051151, 
                0.4701420641051151 };
static double NF_C_T_B3_2D_Eta[] = {
                0.0000000000000000, 0.0000000000000000,
                0.0000000000000000, 0.0000000000000000,
                0.3333333333333333, 0.3333333333333333,
                0.6666666666666667, 0.6666666666666667,
                1.0000000000000000,
                0.3333333333333333, 
                0.1012865073234563, 0.7974269853530873, 
                0.1012865073234563, 
                0.4701420641051151, 0.4701420641051151, 
                0.05971587178976982 };
static double NF_C_T_B3_2D_T[] = {
               -1.0000000000000000, -0.3333333333333333,
                0.3333333333333333,  1.0000000000000000 };

static double NF_C_T_B3_2D_W[] = {
                0.1125, 
                0.0629695902724136, 0.0629695902724136,
                0.0629695902724136, 
                0.0661970763942531, 0.0661970763942531, 
                0.0661970763942531 };

void NF_C_T_B3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = PointValues[6];
  Functionals[7] = PointValues[7];
  Functionals[8] = PointValues[8];
  Functionals[9] =(NF_C_T_B3_2D_W[0]*PointValues[9]
                  +NF_C_T_B3_2D_W[1]*PointValues[10]
                  +NF_C_T_B3_2D_W[2]*PointValues[11]
                  +NF_C_T_B3_2D_W[3]*PointValues[12]
                  +NF_C_T_B3_2D_W[4]*PointValues[13]
                  +NF_C_T_B3_2D_W[5]*PointValues[14]
                  +NF_C_T_B3_2D_W[6]*PointValues[15])*360;
  Functionals[10] =(NF_C_T_B3_2D_W[0]*PointValues[9]*NF_C_T_B3_2D_Xi[9]
                   +NF_C_T_B3_2D_W[1]*PointValues[10]*NF_C_T_B3_2D_Xi[10]
                   +NF_C_T_B3_2D_W[2]*PointValues[11]*NF_C_T_B3_2D_Xi[11]
                   +NF_C_T_B3_2D_W[3]*PointValues[12]*NF_C_T_B3_2D_Xi[12]
                   +NF_C_T_B3_2D_W[4]*PointValues[13]*NF_C_T_B3_2D_Xi[13]
                   +NF_C_T_B3_2D_W[5]*PointValues[14]*NF_C_T_B3_2D_Xi[14]
                   +NF_C_T_B3_2D_W[6]*PointValues[15]*NF_C_T_B3_2D_Xi[15])*360;
  Functionals[11] =(NF_C_T_B3_2D_W[0]*PointValues[9]*NF_C_T_B3_2D_Eta[9]
                   +NF_C_T_B3_2D_W[1]*PointValues[10]*NF_C_T_B3_2D_Eta[10]
                   +NF_C_T_B3_2D_W[2]*PointValues[11]*NF_C_T_B3_2D_Eta[11]
                   +NF_C_T_B3_2D_W[3]*PointValues[12]*NF_C_T_B3_2D_Eta[12]
                   +NF_C_T_B3_2D_W[4]*PointValues[13]*NF_C_T_B3_2D_Eta[13]
                   +NF_C_T_B3_2D_W[5]*PointValues[14]*NF_C_T_B3_2D_Eta[14]
                   +NF_C_T_B3_2D_W[6]*PointValues[15]*NF_C_T_B3_2D_Eta[15])*360;
}

void NF_C_T_B3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
}

TNodalFunctional2D *NF_C_T_B3_2D_Obj = new TNodalFunctional2D
        (NF_C_T_B3_2D, 12, 4, 16, 4, NF_C_T_B3_2D_Xi, NF_C_T_B3_2D_Eta,
         NF_C_T_B3_2D_T, NF_C_T_B3_2D_EvalAll, NF_C_T_B3_2D_EvalEdge);
