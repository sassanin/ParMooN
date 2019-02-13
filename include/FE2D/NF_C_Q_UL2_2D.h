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
   
static double NF_C_Q_UL2_2D_Xi[] = {
-1, 0, 1, 1, 1, 0, -1, -1, -.7745966692414833770358530, 0., .7745966692414833770358530,
-.7745966692414833770358530, 0., .7745966692414833770358530, -.7745966692414833770358530, 0.,
.7745966692414833770358530
};
static double NF_C_Q_UL2_2D_Eta[] = {
-1, -1, -1, 0, 1, 1, 1, 0, -.7745966692414833770358530, -.7745966692414833770358530,
-.7745966692414833770358530, 0., 0., 0., .7745966692414833770358530, .7745966692414833770358530,
.7745966692414833770358530
};
static double NF_C_Q_UL2_2D_T[] = { -1, 0, 1 };

static double NF_C_Q_UL2_2D_W8[] = {
.3086419753086419753086424, .4938271604938271604938270, .3086419753086419753086424,
.4938271604938271604938270, .7901234567901234567901219, .4938271604938271604938270,
.3086419753086419753086424, .4938271604938271604938270, .3086419753086419753086424
};

static double NF_C_Q_UL2_2D_W9[] = {
-.2390730460621862274802019, -.3825168736994979639683223, -.2390730460621862274802019, 0., 0.,
0., .2390730460621862274802019, .3825168736994979639683223, .2390730460621862274802019
};

static double NF_C_Q_UL2_2D_W10[] = {
-.2390730460621862274802019, 0., .2390730460621862274802019, -.3825168736994979639683223, 0.,
.3825168736994979639683223, -.2390730460621862274802019, 0., .2390730460621862274802019
};

void NF_C_Q_UL2_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
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

  Functionals[8] = NF_C_Q_UL2_2D_W8[0]*PointValues[ 8]
                  +NF_C_Q_UL2_2D_W8[1]*PointValues[ 9]
                  +NF_C_Q_UL2_2D_W8[2]*PointValues[10]
                  +NF_C_Q_UL2_2D_W8[3]*PointValues[11]
                  +NF_C_Q_UL2_2D_W8[4]*PointValues[12]
                  +NF_C_Q_UL2_2D_W8[5]*PointValues[13]
                  +NF_C_Q_UL2_2D_W8[6]*PointValues[14]
                  +NF_C_Q_UL2_2D_W8[7]*PointValues[15]
                  +NF_C_Q_UL2_2D_W8[8]*PointValues[16];

  Functionals[9] = NF_C_Q_UL2_2D_W9[0]*PointValues[ 8]
                  +NF_C_Q_UL2_2D_W9[1]*PointValues[ 9]
                  +NF_C_Q_UL2_2D_W9[2]*PointValues[10]
                  +NF_C_Q_UL2_2D_W9[3]*PointValues[11]
                  +NF_C_Q_UL2_2D_W9[4]*PointValues[12]
                  +NF_C_Q_UL2_2D_W9[5]*PointValues[13]
                  +NF_C_Q_UL2_2D_W9[6]*PointValues[14]
                  +NF_C_Q_UL2_2D_W9[7]*PointValues[15]
                  +NF_C_Q_UL2_2D_W9[8]*PointValues[16];

  Functionals[10] = NF_C_Q_UL2_2D_W10[0]*PointValues[ 8]
                   +NF_C_Q_UL2_2D_W10[1]*PointValues[ 9]
                   +NF_C_Q_UL2_2D_W10[2]*PointValues[10]
                   +NF_C_Q_UL2_2D_W10[3]*PointValues[11]
                   +NF_C_Q_UL2_2D_W10[4]*PointValues[12]
                   +NF_C_Q_UL2_2D_W10[5]*PointValues[13]
                   +NF_C_Q_UL2_2D_W10[6]*PointValues[14]
                   +NF_C_Q_UL2_2D_W10[7]*PointValues[15]
                   +NF_C_Q_UL2_2D_W10[8]*PointValues[16];
}

void NF_C_Q_UL2_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_C_Q_UL2_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL2_2D, 11, 3, 17, 3, NF_C_Q_UL2_2D_Xi, NF_C_Q_UL2_2D_Eta,
         NF_C_Q_UL2_2D_T, NF_C_Q_UL2_2D_EvalAll, NF_C_Q_UL2_2D_EvalEdge);
