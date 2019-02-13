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
   
static double NF_C_Q_UL3_2D_Xi[] = {
-1., -.3333333333333333333333333, .3333333333333333333333333, 1., 1., 1., 1.,
.3333333333333333333333333, -.3333333333333333333333333, -1., -1., -1., -.8611363115940525752239466,
-.3399810435848562648026657, .3399810435848562648026657, .8611363115940525752239466,
-.8611363115940525752239466, -.3399810435848562648026657, .3399810435848562648026657,
.8611363115940525752239466, -.8611363115940525752239466, -.3399810435848562648026657,
.3399810435848562648026657, .8611363115940525752239466, -.8611363115940525752239466,
-.3399810435848562648026657, .3399810435848562648026657, .8611363115940525752239466
};
                                  
static double NF_C_Q_UL3_2D_Eta[] = {
-1., -1., -1., -1., -.3333333333333333333333333, .3333333333333333333333333, 1., 1., 1., 1.,
.3333333333333333333333333, -.3333333333333333333333333, -.8611363115940525752239466,
-.8611363115940525752239466, -.8611363115940525752239466, -.8611363115940525752239466,
-.3399810435848562648026657, -.3399810435848562648026657, -.3399810435848562648026657,
-.3399810435848562648026657, .3399810435848562648026657, .3399810435848562648026657,
.3399810435848562648026657, .3399810435848562648026657, .8611363115940525752239466,
.8611363115940525752239466, .8611363115940525752239466, .8611363115940525752239466
};

static double NF_C_Q_UL3_2D_T[] = { -1, -0.33333333333333333333,
                                 0.33333333333333333333, 1 };

static double NF_C_Q_UL3_2D_W12[] = {
.1210029932856020055212120, .2268518518518518518518517, .2268518518518518518518520,
.1210029932856020055212121, .2268518518518518518518517, .4252933030106942907750839,
.4252933030106942907750845, .2268518518518518518518518, .2268518518518518518518520,
.4252933030106942907750845, .4252933030106942907750852, .2268518518518518518518521,
.1210029932856020055212121, .2268518518518518518518518, .2268518518518518518518521,
.1210029932856020055212121
};

static double NF_C_Q_UL3_2D_W13[] = {
-.1042000713298032202198559, -.1953503669819841490091359, -.1953503669819841490091362,
-.1042000713298032202198560, -0.7712532933174980081171578e-1, -.1445916609872263377762813,
-.1445916609872263377762815, -0.7712532933174980081171581e-1, 0.7712532933174980081171588e-1,
.1445916609872263377762815, .1445916609872263377762817, 0.7712532933174980081171591e-1,
.1042000713298032202198560, .1953503669819841490091360, .1953503669819841490091362,
.1042000713298032202198560
};

static double NF_C_Q_UL3_2D_W14[] = {
-.1042000713298032202198559, -0.7712532933174980081171578e-1, 0.7712532933174980081171588e-1,
.1042000713298032202198560, -.1953503669819841490091359, -.1445916609872263377762813,
.1445916609872263377762815, .1953503669819841490091360, -.1953503669819841490091362,
-.1445916609872263377762815, .1445916609872263377762817, .1953503669819841490091362,
-.1042000713298032202198560, -0.7712532933174980081171581e-1, 0.7712532933174980081171591e-1,
.1042000713298032202198560
};

static double NF_C_Q_UL3_2D_W15[] = {
0.8973046509278393012803587e-1, 0.6641542163121961915050669e-1,
-0.6641542163121961915050678e-1, -0.8973046509278393012803594e-1, 0.6641542163121961915050669e-1,
0.4915842379610495876085292e-1, -0.4915842379610495876085299e-1, -0.6641542163121961915050672e-1,
-0.6641542163121961915050678e-1, -0.4915842379610495876085299e-1, 0.4915842379610495876085307e-1,
0.6641542163121961915050681e-1, -0.8973046509278393012803594e-1, -0.6641542163121961915050672e-1,
0.6641542163121961915050681e-1, 0.8973046509278393012803594e-1
};

static double NF_C_Q_UL3_2D_W16[] = {
0.7409420099637489243144783e-1, .1389090158111897072462627, .1389090158111897072462628,
0.7409420099637489243144789e-1, -0.7409420099637489243144779e-1, -.1389090158111897072462626,
-.1389090158111897072462628, -0.7409420099637489243144782e-1, -0.7409420099637489243144789e-1,
-.1389090158111897072462628, -.1389090158111897072462630, -0.7409420099637489243144792e-1,
0.7409420099637489243144789e-1, .1389090158111897072462627, .1389090158111897072462629,
0.7409420099637489243144789e-1
};

static double NF_C_Q_UL3_2D_W17[] = {
0.7409420099637489243144783e-1, -0.7409420099637489243144779e-1,
-0.7409420099637489243144789e-1, 0.7409420099637489243144789e-1, .1389090158111897072462627,
-.1389090158111897072462626, -.1389090158111897072462628, .1389090158111897072462627,
.1389090158111897072462628, -.1389090158111897072462628, -.1389090158111897072462630,
.1389090158111897072462629, 0.7409420099637489243144789e-1, -0.7409420099637489243144782e-1,
-0.7409420099637489243144792e-1, 0.7409420099637489243144789e-1
};

void NF_C_Q_UL3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
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
  Functionals[9] = PointValues[9];
  Functionals[10] = PointValues[10];
  Functionals[11] = PointValues[11];
  Functionals[12] = NF_C_Q_UL3_2D_W12[0]*PointValues[12]
                   +NF_C_Q_UL3_2D_W12[1]*PointValues[13]
                   +NF_C_Q_UL3_2D_W12[2]*PointValues[14]
                   +NF_C_Q_UL3_2D_W12[3]*PointValues[15]
                   +NF_C_Q_UL3_2D_W12[4]*PointValues[16]
                   +NF_C_Q_UL3_2D_W12[5]*PointValues[17]
                   +NF_C_Q_UL3_2D_W12[6]*PointValues[18]
                   +NF_C_Q_UL3_2D_W12[7]*PointValues[19]
                   +NF_C_Q_UL3_2D_W12[8]*PointValues[20]
                   +NF_C_Q_UL3_2D_W12[9]*PointValues[21]
                   +NF_C_Q_UL3_2D_W12[10]*PointValues[22]
                   +NF_C_Q_UL3_2D_W12[11]*PointValues[23]
                   +NF_C_Q_UL3_2D_W12[12]*PointValues[24]
                   +NF_C_Q_UL3_2D_W12[13]*PointValues[25]
                   +NF_C_Q_UL3_2D_W12[14]*PointValues[26]
                   +NF_C_Q_UL3_2D_W12[15]*PointValues[27];
  Functionals[13] = NF_C_Q_UL3_2D_W13[0]*PointValues[12]
                   +NF_C_Q_UL3_2D_W13[1]*PointValues[13]
                   +NF_C_Q_UL3_2D_W13[2]*PointValues[14]
                   +NF_C_Q_UL3_2D_W13[3]*PointValues[15]
                   +NF_C_Q_UL3_2D_W13[4]*PointValues[16]
                   +NF_C_Q_UL3_2D_W13[5]*PointValues[17]
                   +NF_C_Q_UL3_2D_W13[6]*PointValues[18]
                   +NF_C_Q_UL3_2D_W13[7]*PointValues[19]
                   +NF_C_Q_UL3_2D_W13[8]*PointValues[20]
                   +NF_C_Q_UL3_2D_W13[9]*PointValues[21]
                   +NF_C_Q_UL3_2D_W13[10]*PointValues[22]
                   +NF_C_Q_UL3_2D_W13[11]*PointValues[23]
                   +NF_C_Q_UL3_2D_W13[12]*PointValues[24]
                   +NF_C_Q_UL3_2D_W13[13]*PointValues[25]
                   +NF_C_Q_UL3_2D_W13[14]*PointValues[26]
                   +NF_C_Q_UL3_2D_W13[15]*PointValues[27];
  Functionals[14] = NF_C_Q_UL3_2D_W14[0]*PointValues[12]
                   +NF_C_Q_UL3_2D_W14[1]*PointValues[13]
                   +NF_C_Q_UL3_2D_W14[2]*PointValues[14]
                   +NF_C_Q_UL3_2D_W14[3]*PointValues[15]
                   +NF_C_Q_UL3_2D_W14[4]*PointValues[16]
                   +NF_C_Q_UL3_2D_W14[5]*PointValues[17]
                   +NF_C_Q_UL3_2D_W14[6]*PointValues[18]
                   +NF_C_Q_UL3_2D_W14[7]*PointValues[19]
                   +NF_C_Q_UL3_2D_W14[8]*PointValues[20]
                   +NF_C_Q_UL3_2D_W14[9]*PointValues[21]
                   +NF_C_Q_UL3_2D_W14[10]*PointValues[22]
                   +NF_C_Q_UL3_2D_W14[11]*PointValues[23]
                   +NF_C_Q_UL3_2D_W14[12]*PointValues[24]
                   +NF_C_Q_UL3_2D_W14[13]*PointValues[25]
                   +NF_C_Q_UL3_2D_W14[14]*PointValues[26]
                   +NF_C_Q_UL3_2D_W14[15]*PointValues[27];
  Functionals[15] = NF_C_Q_UL3_2D_W15[0]*PointValues[12]
                   +NF_C_Q_UL3_2D_W15[1]*PointValues[13]
                   +NF_C_Q_UL3_2D_W15[2]*PointValues[14]
                   +NF_C_Q_UL3_2D_W15[3]*PointValues[15]
                   +NF_C_Q_UL3_2D_W15[4]*PointValues[16]
                   +NF_C_Q_UL3_2D_W15[5]*PointValues[17]
                   +NF_C_Q_UL3_2D_W15[6]*PointValues[18]
                   +NF_C_Q_UL3_2D_W15[7]*PointValues[19]
                   +NF_C_Q_UL3_2D_W15[8]*PointValues[20]
                   +NF_C_Q_UL3_2D_W15[9]*PointValues[21]
                   +NF_C_Q_UL3_2D_W15[10]*PointValues[22]
                   +NF_C_Q_UL3_2D_W15[11]*PointValues[23]
                   +NF_C_Q_UL3_2D_W15[12]*PointValues[24]
                   +NF_C_Q_UL3_2D_W15[13]*PointValues[25]
                   +NF_C_Q_UL3_2D_W15[14]*PointValues[26]
                   +NF_C_Q_UL3_2D_W15[15]*PointValues[27];
  Functionals[16] = NF_C_Q_UL3_2D_W16[0]*PointValues[12]
                   +NF_C_Q_UL3_2D_W16[1]*PointValues[13]
                   +NF_C_Q_UL3_2D_W16[2]*PointValues[14]
                   +NF_C_Q_UL3_2D_W16[3]*PointValues[15]
                   +NF_C_Q_UL3_2D_W16[4]*PointValues[16]
                   +NF_C_Q_UL3_2D_W16[5]*PointValues[17]
                   +NF_C_Q_UL3_2D_W16[6]*PointValues[18]
                   +NF_C_Q_UL3_2D_W16[7]*PointValues[19]
                   +NF_C_Q_UL3_2D_W16[8]*PointValues[20]
                   +NF_C_Q_UL3_2D_W16[9]*PointValues[21]
                   +NF_C_Q_UL3_2D_W16[10]*PointValues[22]
                   +NF_C_Q_UL3_2D_W16[11]*PointValues[23]
                   +NF_C_Q_UL3_2D_W16[12]*PointValues[24]
                   +NF_C_Q_UL3_2D_W16[13]*PointValues[25]
                   +NF_C_Q_UL3_2D_W16[14]*PointValues[26]
                   +NF_C_Q_UL3_2D_W16[15]*PointValues[27];
  Functionals[17] = NF_C_Q_UL3_2D_W17[0]*PointValues[12]
                   +NF_C_Q_UL3_2D_W17[1]*PointValues[13]
                   +NF_C_Q_UL3_2D_W17[2]*PointValues[14]
                   +NF_C_Q_UL3_2D_W17[3]*PointValues[15]
                   +NF_C_Q_UL3_2D_W17[4]*PointValues[16]
                   +NF_C_Q_UL3_2D_W17[5]*PointValues[17]
                   +NF_C_Q_UL3_2D_W17[6]*PointValues[18]
                   +NF_C_Q_UL3_2D_W17[7]*PointValues[19]
                   +NF_C_Q_UL3_2D_W17[8]*PointValues[20]
                   +NF_C_Q_UL3_2D_W17[9]*PointValues[21]
                   +NF_C_Q_UL3_2D_W17[10]*PointValues[22]
                   +NF_C_Q_UL3_2D_W17[11]*PointValues[23]
                   +NF_C_Q_UL3_2D_W17[12]*PointValues[24]
                   +NF_C_Q_UL3_2D_W17[13]*PointValues[25]
                   +NF_C_Q_UL3_2D_W17[14]*PointValues[26]
                   +NF_C_Q_UL3_2D_W17[15]*PointValues[27];
}

void NF_C_Q_UL3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  memcpy(Functionals, PointValues, 4*SizeOfDouble);
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_C_Q_UL3_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL3_2D, 18, 4, 28, 4, NF_C_Q_UL3_2D_Xi, NF_C_Q_UL3_2D_Eta,
         NF_C_Q_UL3_2D_T, NF_C_Q_UL3_2D_EvalAll, NF_C_Q_UL3_2D_EvalEdge);
