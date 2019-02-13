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
   
static double NF_D_Q_D2_2D_W0[16] = {
0.302507483214005E-1, 0.5671296296296296E-1, 0.5671296296296296E-1, 0.302507483214005E-1,
0.5671296296296296E-1, 0.1063233257526736, 0.1063233257526736, 0.5671296296296296E-1,
0.5671296296296296E-1, 0.1063233257526736, 0.1063233257526736, 0.5671296296296296E-1,
0.302507483214005E-1, 0.5671296296296296E-1, 0.5671296296296296E-1, 0.302507483214005E-1 };

static double NF_D_Q_D2_2D_W1[16] = {
-0.7815005349735242E-1, -0.5784399699881235E-1, 0.5784399699881235E-1, 0.7815005349735242E-1,
-0.1465127752364881, -0.1084437457404198, 0.1084437457404198, 0.1465127752364881,
-0.1465127752364881, -0.1084437457404198, 0.1084437457404198, 0.1465127752364881,
-0.7815005349735242E-1, -0.5784399699881235E-1, 0.5784399699881235E-1, 0.7815005349735242E-1 };

static double NF_D_Q_D2_2D_W2[16] = {
-0.7815005349735242E-1, -0.1465127752364881, -0.1465127752364881, -0.7815005349735242E-1,
-0.5784399699881235E-1, -0.1084437457404198, -0.1084437457404198, -0.5784399699881235E-1,
0.5784399699881235E-1, 0.1084437457404198, 0.1084437457404198, 0.5784399699881235E-1,
0.7815005349735242E-1, 0.1465127752364881, 0.1465127752364881, 0.7815005349735242E-1 };

static double NF_D_Q_D2_2D_W3[16] = {
0.9261775124546862E-1, -0.9261775124546862E-1, -0.9261775124546862E-1, 0.9261775124546862E-1,
0.1736362697639871, -0.1736362697639871, -0.1736362697639871, 0.1736362697639871,
0.1736362697639871, -0.1736362697639871, -0.1736362697639871, 0.1736362697639871,
0.9261775124546862E-1, -0.9261775124546862E-1, -0.9261775124546862E-1, 0.9261775124546862E-1 };

static double NF_D_Q_D2_2D_W4[16] = {
0.2018935464587638, 0.1494346986702441, -0.1494346986702441, -0.2018935464587638,
0.1494346986702441, 0.1106064535412362, -0.1106064535412362, -0.1494346986702441,
-0.1494346986702441, -0.1106064535412362, 0.1106064535412362, 0.1494346986702441,
-0.2018935464587638, -0.1494346986702441, 0.1494346986702441, 0.2018935464587638 };

static double NF_D_Q_D2_2D_W5[16] = {
0.9261775124546862E-1, 0.1736362697639871, 0.1736362697639871, 0.9261775124546862E-1,
-0.9261775124546862E-1, -0.1736362697639871, -0.1736362697639871, -0.9261775124546862E-1,
-0.9261775124546862E-1, -0.1736362697639871, -0.1736362697639871, -0.9261775124546862E-1,
0.9261775124546862E-1, 0.1736362697639871, 0.1736362697639871, 0.9261775124546862E-1 };

static double NF_D_Q_D2_2D_W6[16] = {
-0.2392695260869749, 0.2392695260869749, 0.2392695260869749, -0.2392695260869749,
-0.1770991205956259, 0.1770991205956259, 0.1770991205956259, -0.1770991205956259,
0.1770991205956259, -0.1770991205956259, -0.1770991205956259, 0.1770991205956259,
0.2392695260869749, -0.2392695260869749, -0.2392695260869749, 0.2392695260869749 };

static double NF_D_Q_D2_2D_W7[16] = {
-0.2392695260869749, -0.1770991205956259, 0.1770991205956259, 0.2392695260869749,
0.2392695260869749, 0.1770991205956259, -0.1770991205956259, -0.2392695260869749,
0.2392695260869749, 0.1770991205956259, -0.1770991205956259, -0.2392695260869749,
-0.2392695260869749, -0.1770991205956259, 0.1770991205956259, 0.2392695260869749 };

static double NF_D_Q_D2_2D_W8[16] = {
0.0, -0.2728291844522666, 0.2728291844522666, 0.0,
0.2728291844522666, 0.0, 0.0, -0.2728291844522666,
-0.2728291844522666, 0.0, 0.0, 0.2728291844522666,
0.0, 0.2728291844522666, -0.2728291844522666, 0.0 };

static double NF_D_Q_D2_2D_Xi[16] = {
-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526,
-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526,
-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526,
-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 };

static double NF_D_Q_D2_2D_Eta[16] = {
-0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526,
-0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563,
0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563,
0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526 };

static double *NF_D_Q_D2_2D_t = NULL;

void NF_D_Q_D2_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] =  NF_D_Q_D2_2D_W0[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W0[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W0[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W0[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W0[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W0[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W0[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W0[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W0[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W0[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W0[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W0[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W0[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W0[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W0[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W0[15]*PointValues[15];
  Functionals[1] =  NF_D_Q_D2_2D_W1[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W1[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W1[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W1[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W1[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W1[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W1[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W1[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W1[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W1[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W1[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W1[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W1[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W1[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W1[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W1[15]*PointValues[15];
  Functionals[2] =  NF_D_Q_D2_2D_W2[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W2[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W2[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W2[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W2[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W2[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W2[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W2[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W2[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W2[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W2[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W2[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W2[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W2[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W2[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W2[15]*PointValues[15];
  Functionals[3] =  NF_D_Q_D2_2D_W3[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W3[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W3[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W3[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W3[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W3[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W3[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W3[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W3[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W3[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W3[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W3[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W3[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W3[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W3[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W3[15]*PointValues[15];
  Functionals[4] =  NF_D_Q_D2_2D_W4[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W4[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W4[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W4[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W4[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W4[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W4[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W4[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W4[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W4[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W4[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W4[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W4[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W4[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W4[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W4[15]*PointValues[15];
  Functionals[5] =  NF_D_Q_D2_2D_W5[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W5[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W5[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W5[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W5[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W5[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W5[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W5[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W5[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W5[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W5[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W5[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W5[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W5[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W5[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W5[15]*PointValues[15];
  Functionals[6] =  NF_D_Q_D2_2D_W6[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W6[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W6[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W6[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W6[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W6[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W6[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W6[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W6[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W6[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W6[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W6[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W6[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W6[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W6[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W6[15]*PointValues[15];
  Functionals[7] =  NF_D_Q_D2_2D_W7[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W7[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W7[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W7[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W7[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W7[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W7[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W7[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W7[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W7[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W7[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W7[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W7[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W7[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W7[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W7[15]*PointValues[15];
  Functionals[8] =  NF_D_Q_D2_2D_W8[ 0]*PointValues[ 0]
                   +NF_D_Q_D2_2D_W8[ 1]*PointValues[ 1]
                   +NF_D_Q_D2_2D_W8[ 2]*PointValues[ 2]
                   +NF_D_Q_D2_2D_W8[ 3]*PointValues[ 3]
                   +NF_D_Q_D2_2D_W8[ 4]*PointValues[ 4]
                   +NF_D_Q_D2_2D_W8[ 5]*PointValues[ 5]
                   +NF_D_Q_D2_2D_W8[ 6]*PointValues[ 6]
                   +NF_D_Q_D2_2D_W8[ 7]*PointValues[ 7]
                   +NF_D_Q_D2_2D_W8[ 8]*PointValues[ 8]
                   +NF_D_Q_D2_2D_W8[ 9]*PointValues[ 9]
                   +NF_D_Q_D2_2D_W8[10]*PointValues[10]
                   +NF_D_Q_D2_2D_W8[11]*PointValues[11]
                   +NF_D_Q_D2_2D_W8[12]*PointValues[12]
                   +NF_D_Q_D2_2D_W8[13]*PointValues[13]
                   +NF_D_Q_D2_2D_W8[14]*PointValues[14]
                   +NF_D_Q_D2_2D_W8[15]*PointValues[15];
}

void NF_D_Q_D2_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_D_Q_D2_2D_Obj = new TNodalFunctional2D
        (NF_D_Q_D2_2D, 9, 0, 16, 0, NF_D_Q_D2_2D_Xi, NF_D_Q_D2_2D_Eta,
         NF_D_Q_D2_2D_t, NF_D_Q_D2_2D_EvalAll, NULL);
