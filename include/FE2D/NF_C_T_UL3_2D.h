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

static double NF_C_T_UL3_2D_W0[] = {
                  0.1338815352798323483836615,
                  0.08781911358244156480372712,
                  0.05857143528033178431852058,
                  0.05306124886956075869379972,
                  3.211668771336243759680046e-10,
                  0.1338815352798323483836615,
                  0.08781911358244156480372712,
                  0.05857143528033178431852058,
                  0.05306124886956075869379972,
                  3.211668771336243759680046e-10,
                  0.1338815352798323483836615,
                  0.08781911358244156480372712,
                  0.05857143528033178431852058,
                  0.05306124886956075869379972,
                  3.211668771336243759680046e-10 };

static double NF_C_T_UL3_2D_W1[] = {
                  0.02060002306022122324914420,
                  0.01129019454258610471158756,
                  0.01807012227527386947012354,
                0.000009792147884325017198945303,
                5.322943133234498398641960e-10,
                 -0.01555072009943916871059616,
                  0.01393558585908052317821126,
                 -0.009907282607036290058979157,
                  0.02135882696742077274497469,
                5.262923413432595370639873e-10,
                 -0.005049302960782054538548041,
                 -0.02522578040166662788979882,
                 -0.008162839668237579411144377,
                 -0.02136861911530509776217363,
                -1.058586654666709376928183e-9 };

static double NF_C_T_UL3_2D_W2[] = {
                 0.005049302960782054538548041,
                  0.02522578040166662788979882,
                 0.008162839668237579411144374,
                  0.02136861911530509776217364,
                1.058586654666709376928183e-9,
                 -0.02060002306022122324914420,
                 -0.01129019454258610471158757,
                 -0.01807012227527386947012353,
               -0.000009792147884325017198945303,
               -5.322943133234498398641958e-10,
                  0.01555072009943916871059616,
                 -0.01393558585908052317821126,
                 0.009907282607036290058979162,
                 -0.02135882696742077274497469,
               -5.262923413432595370639871e-10 };

static double NF_C_T_UL3_2D_W3[] = {
                  -.5681561240792843260361997,
                  -.4147449330059048571735955,
                  0.01787438064203683735652336,
                  0.2490748470886491633963465,
                -1.649363563393114180399687e-9,
                  -.5524478787902309208087563,
                 -0.05385400553858738142636380,
                  -.2768234313569137721564739,
                  0.7691212129434811378559316,
                2.985545465286077101560408e-7,
                  -.2676548426456185647840447,
                  0.4603790251744214261310966,
                  0.3869750551988903332977606,
                  0.2502560711293055021705686,
                3.263345724569626099174574e-7 };

static double NF_C_T_UL3_2D_W4[] = {
                  0.8529491602238966820609113,
                  0.9289779637189136647310559,
                  0.6459241059137672680977100,
                  -.7679399889028247990817089,
                2.942938949174801394181639e-8,
                  0.2519465973565651595565995,
                  -.8212699526417389018783287,
                 -0.09227724319993972378476302,
                  -.7703024369841374766301528,
                -6.265384825489634342538976e-7,
                  0.2833630879346719700114886,
                 -0.09948809770710395038386405,
                  -.6816728671978409428107582,
                  0.2697902947255264722890158,
                -2.613066236496178558101700e-8 };

static double NF_C_T_UL3_2D_W5[] = {
                  -.5524478787902309208087563,
                 -0.05385400553858738142636380,
                  -.2768234313569137721564739,
                  0.7691212129434811378559316,
                2.985545465286077101560408e-7,
                  -.2676548426456185647840447,
                  0.4603790251744214261310966,
                  0.3869750551988903332977606,
                  0.2502560711293055021705686,
                3.263345724569626099174574e-7,
                  -.5681561240792843260361997,
                  -.4147449330059048571735955,
                  0.01787438064203683735652336,
                  0.2490748470886491633963465,
                -1.649363563393114180399687e-9 };

static double NF_C_T_UL3_2D_Xi[24] = {
                   0.0, 0.33333333333333333333, 0.66666666666666666667,
                   1.0, 0.66666666666666666667, 0.33333333333333333333,
                   0.0, 0.0, 0.0,
                   0.513346920639454149493588964941289,
                   0.313251210671725306955957434773977,
                   0.651775303648795707537235440171623,
                   0.0651019934589391663281047918600097,
                   0.345792011168269028821402703032657,
                   0.281041247315110390572737910356218,
                   0.630621434318956140102957434914989,
                   0.313477887523733007173578085409864,
                   0.870165101563563060777479924328995,
                   3.62316822156926166667912853628633,
                   0.205611832045435459933673124702493,
                   0.056127355009318552941085130311034,
                   0.034746808827471285289186474418513,
                   0.0647329049774977728944152838109953,
                  -2.968960232737530695500531239318987 };

static double NF_C_T_UL3_2D_Eta[24] = {
                   0.0, 0.0, 0.0, 0.0,
                   0.33333333333333333333, 0.66666666666666666667,
                   1.0, 0.66666666666666666667, 0.33333333333333333333,
                   0.281041247315110390572737910356218,
                   0.630621434318956140102957434914989,
                   0.313477887523733007173578085409864,
                   0.870165101563563060777479924328995,
                   3.62316822156926166667912853628633,
                   0.205611832045435459933673124702493,
                   0.056127355009318552941085130311034,
                   0.034746808827471285289186474418513,
                   0.0647329049774977728944152838109953,
                  -2.968960232737530695500531239318987,
                   0.513346920639454149493588964941289,
                   0.313251210671725306955957434773977,
                   0.651775303648795707537235440171623,
                   0.0651019934589391663281047918600097,
                   0.345792011168269028821402703032657 };

static double NF_C_T_UL3_2D_T[4] = { -1, -0.33333333333333333333,
                                     0.33333333333333333333, 1 };

void NF_C_T_UL3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
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

  Functionals[9] =      PointValues[ 9]*NF_C_T_UL3_2D_W0[ 0]
                      + PointValues[10]*NF_C_T_UL3_2D_W0[ 1]
                      + PointValues[11]*NF_C_T_UL3_2D_W0[ 2]
                      + PointValues[12]*NF_C_T_UL3_2D_W0[ 3]
                      + PointValues[13]*NF_C_T_UL3_2D_W0[ 4]
                      + PointValues[14]*NF_C_T_UL3_2D_W0[ 5]
                      + PointValues[15]*NF_C_T_UL3_2D_W0[ 6]
                      + PointValues[16]*NF_C_T_UL3_2D_W0[ 7]
                      + PointValues[17]*NF_C_T_UL3_2D_W0[ 8]
                      + PointValues[18]*NF_C_T_UL3_2D_W0[ 9]
                      + PointValues[19]*NF_C_T_UL3_2D_W0[10]
                      + PointValues[20]*NF_C_T_UL3_2D_W0[11]
                      + PointValues[21]*NF_C_T_UL3_2D_W0[12]
                      + PointValues[22]*NF_C_T_UL3_2D_W0[13]
                      + PointValues[23]*NF_C_T_UL3_2D_W0[14];

  Functionals[10] =     PointValues[ 9]*NF_C_T_UL3_2D_W1[ 0]
                      + PointValues[10]*NF_C_T_UL3_2D_W1[ 1]
                      + PointValues[11]*NF_C_T_UL3_2D_W1[ 2]
                      + PointValues[12]*NF_C_T_UL3_2D_W1[ 3]
                      + PointValues[13]*NF_C_T_UL3_2D_W1[ 4]
                      + PointValues[14]*NF_C_T_UL3_2D_W1[ 5]
                      + PointValues[15]*NF_C_T_UL3_2D_W1[ 6]
                      + PointValues[16]*NF_C_T_UL3_2D_W1[ 7]
                      + PointValues[17]*NF_C_T_UL3_2D_W1[ 8]
                      + PointValues[18]*NF_C_T_UL3_2D_W1[ 9]
                      + PointValues[19]*NF_C_T_UL3_2D_W1[10]
                      + PointValues[20]*NF_C_T_UL3_2D_W1[11]
                      + PointValues[21]*NF_C_T_UL3_2D_W1[12]
                      + PointValues[22]*NF_C_T_UL3_2D_W1[13]
                      + PointValues[23]*NF_C_T_UL3_2D_W1[14];

  Functionals[11] =     PointValues[ 9]*NF_C_T_UL3_2D_W2[ 0]
                      + PointValues[10]*NF_C_T_UL3_2D_W2[ 1]
                      + PointValues[11]*NF_C_T_UL3_2D_W2[ 2]
                      + PointValues[12]*NF_C_T_UL3_2D_W2[ 3]
                      + PointValues[13]*NF_C_T_UL3_2D_W2[ 4]
                      + PointValues[14]*NF_C_T_UL3_2D_W2[ 5]
                      + PointValues[15]*NF_C_T_UL3_2D_W2[ 6]
                      + PointValues[16]*NF_C_T_UL3_2D_W2[ 7]
                      + PointValues[17]*NF_C_T_UL3_2D_W2[ 8]
                      + PointValues[18]*NF_C_T_UL3_2D_W2[ 9]
                      + PointValues[19]*NF_C_T_UL3_2D_W2[10]
                      + PointValues[20]*NF_C_T_UL3_2D_W2[11]
                      + PointValues[21]*NF_C_T_UL3_2D_W2[12]
                      + PointValues[22]*NF_C_T_UL3_2D_W2[13]
                      + PointValues[23]*NF_C_T_UL3_2D_W2[14];

  Functionals[12] =     PointValues[ 9]*NF_C_T_UL3_2D_W3[ 0]
                      + PointValues[10]*NF_C_T_UL3_2D_W3[ 1]
                      + PointValues[11]*NF_C_T_UL3_2D_W3[ 2]
                      + PointValues[12]*NF_C_T_UL3_2D_W3[ 3]
                      + PointValues[13]*NF_C_T_UL3_2D_W3[ 4]
                      + PointValues[14]*NF_C_T_UL3_2D_W3[ 5]
                      + PointValues[15]*NF_C_T_UL3_2D_W3[ 6]
                      + PointValues[16]*NF_C_T_UL3_2D_W3[ 7]
                      + PointValues[17]*NF_C_T_UL3_2D_W3[ 8]
                      + PointValues[18]*NF_C_T_UL3_2D_W3[ 9]
                      + PointValues[19]*NF_C_T_UL3_2D_W3[10]
                      + PointValues[20]*NF_C_T_UL3_2D_W3[11]
                      + PointValues[21]*NF_C_T_UL3_2D_W3[12]
                      + PointValues[22]*NF_C_T_UL3_2D_W3[13]
                      + PointValues[23]*NF_C_T_UL3_2D_W3[14];

  Functionals[13] =     PointValues[ 9]*NF_C_T_UL3_2D_W4[ 0]
                      + PointValues[10]*NF_C_T_UL3_2D_W4[ 1]
                      + PointValues[11]*NF_C_T_UL3_2D_W4[ 2]
                      + PointValues[12]*NF_C_T_UL3_2D_W4[ 3]
                      + PointValues[13]*NF_C_T_UL3_2D_W4[ 4]
                      + PointValues[14]*NF_C_T_UL3_2D_W4[ 5]
                      + PointValues[15]*NF_C_T_UL3_2D_W4[ 6]
                      + PointValues[16]*NF_C_T_UL3_2D_W4[ 7]
                      + PointValues[17]*NF_C_T_UL3_2D_W4[ 8]
                      + PointValues[18]*NF_C_T_UL3_2D_W4[ 9]
                      + PointValues[19]*NF_C_T_UL3_2D_W4[10]
                      + PointValues[20]*NF_C_T_UL3_2D_W4[11]
                      + PointValues[21]*NF_C_T_UL3_2D_W4[12]
                      + PointValues[22]*NF_C_T_UL3_2D_W4[13]
                      + PointValues[23]*NF_C_T_UL3_2D_W4[14];

  Functionals[14] =     PointValues[ 9]*NF_C_T_UL3_2D_W5[ 0]
                      + PointValues[10]*NF_C_T_UL3_2D_W5[ 1]
                      + PointValues[11]*NF_C_T_UL3_2D_W5[ 2]
                      + PointValues[12]*NF_C_T_UL3_2D_W5[ 3]
                      + PointValues[13]*NF_C_T_UL3_2D_W5[ 4]
                      + PointValues[14]*NF_C_T_UL3_2D_W5[ 5]
                      + PointValues[15]*NF_C_T_UL3_2D_W5[ 6]
                      + PointValues[16]*NF_C_T_UL3_2D_W5[ 7]
                      + PointValues[17]*NF_C_T_UL3_2D_W5[ 8]
                      + PointValues[18]*NF_C_T_UL3_2D_W5[ 9]
                      + PointValues[19]*NF_C_T_UL3_2D_W5[10]
                      + PointValues[20]*NF_C_T_UL3_2D_W5[11]
                      + PointValues[21]*NF_C_T_UL3_2D_W5[12]
                      + PointValues[22]*NF_C_T_UL3_2D_W5[13]
                      + PointValues[23]*NF_C_T_UL3_2D_W5[14];
}

void NF_C_T_UL3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
}

TNodalFunctional2D *NF_C_T_UL3_2D_Obj = new TNodalFunctional2D
        (NF_C_T_UL3_2D, 15, 4, 24, 4, NF_C_T_UL3_2D_Xi, NF_C_T_UL3_2D_Eta,
         NF_C_T_UL3_2D_T, NF_C_T_UL3_2D_EvalAll, NF_C_T_UL3_2D_EvalEdge);
