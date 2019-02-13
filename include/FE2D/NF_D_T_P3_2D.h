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

static double NF_D_T_P3_2D_Xi[15] = {
    0.5133469206394541, 0.3132512106717253, 0.6517753036487957,
    0.6510199345893917E-1, 0.345792011168269, 0.2810412473151104,
    0.6306214343189561, 0.313477887523733, 0.8701651015635631,
    0.3623168221569262E1, 0.2056118320454355, 0.5612735500931855E-1,
    0.3474680882747129E-1, 0.6473290497749777E-1, -0.2968960232737531E1 };

static double NF_D_T_P3_2D_Eta[15] = {
    0.2810412473151104, 0.6306214343189561, 0.313477887523733,
    0.8701651015635631, 0.3623168221569262E1, 0.2056118320454355,
    0.5612735500931855E-1, 0.3474680882747129E-1, 0.6473290497749777E-1,
    -0.2968960232737531E1, 0.5133469206394541, 0.3132512106717253,
    0.6517753036487957, 0.6510199345893917E-1, 0.345792011168269 };


void NF_D_T_P3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // w = 1
  static double w0[15] = {
        .6694076763991617e-1, .4390955679122078e-1, .2928571764016589e-1,
        .2653062443478038e-1, .1605834385668122e-9, .6694076763991617e-1,
        .4390955679122078e-1, .2928571764016589e-1, .2653062443478038e-1,
        .1605834385668122e-9, .6694076763991617e-1, .4390955679122078e-1,
        .2928571764016589e-1, .2653062443478038e-1, .1605834385668122e-9 };

  // w = 2*xi+eta-1
  static double w1[15] = {
        .206000230602212159639612674e-1,
        .112901945425861017222420458e-1, .180701222752738674862483255e-1,
        .97921478843262528732922309e-5, .532294313323449923963583090e-9,
        -.155507200994391637926030882e-1, .139355858590805187649833588e-1,
        -.990728260703628961071419105e-2, .213588269674207752751667694e-1,
        .526292341343259634598983724e-9, -.504930296078205217135817915e-2,
        -.252257804016666270736589233e-1, -.816283966823757875410566365e-2,
        -.213686191153050983443651294e-1, -.105858665466670955856256681e-8 };

  // w = xi+2*eta-1
  static double w2[15] = {
        .504930296078205217135817915e-2,
        .252257804016666226827032442e-1, .816283966823757816839131085e-2,
        .213686191153051004668150842e-1, .105858665466670955856256681e-8,
        -.206000230602212159639612674e-1, -.112901945425861061131977249e-1,
        -.180701222752738680719626783e-1, -.97921478843241304233374485e-5,
        -.532294313323449923963583090e-9, .155507200994391637926030882e-1,
        -.139355858590805231559390379e-1, .990728260703628902499983824e-2,
        -.213588269674207731527168146e-1, -.526292341343259634598983724e-9 };

  // w = 180*xi^2-144*xi+18
  static double w3[15] = {
        -.568156124079284425650535186, -.414744933005904824948405918,
        .178743806420368160315568071e-1, .249074847088649157782787817,
        -.164936356339311421349730024e-8, -.552447878790230913240060076,
        -.538540055385875246772541652e-1, -.276823431356913745202383217,
        .769121212943481332920287488, .298554546528607794604190635e-6,
        -.267654842645618735714108454, .460379025174421416935243448,
        .386975055198890286625307968, .250256071129305517599417835,
        .326334572456962693640367494e-6 };

  // w = 360*xi*eta-72*xi-72*eta+18
  static double w4[15] = {
        .852949160223896603176486808, .928977963718913494773585717,
        .645924105913767177288007812, -.767939988902824845017883949,
        .294293894917480132496741587e-7, .251946597356565223303633344,
        -.821269952641738785770211250, -.922772431999396988470543191e-1,
        -.770302436984137564792150915, -.626538482548963602458055430e-6,
        .283363087934672248124583564, -.994880977071039343873622809e-1,
        -.681672867197840856981151753, .269790294725526477915740096,
        -.261306623649617848226795582e-7 };

  // w = 180*eta^2-144*eta+18
  static double w5[15] = {
        -.552447878790230913240060076,
        -.538540055385875246772541652e-1, -.276823431356913745202383217,
        .769121212943481332920287488, .298554546528607794604190635e-6,
        -.267654842645618735714108454, .460379025174421416935243448,
        .386975055198890286625307968, .250256071129305517599417835,
        .326334572456962693640367494e-6, -.568156124079284425650535186,
        -.414744933005904824948405918, .178743806420368160315568071e-1,
        .249074847088649157782787817, -.164936356339311421349730024e-8 };

  // 35*xi^3-45*xi^2+15*xi-1
  static double w6[15] = {
        -.283582300812142314864875229e-1,
        .157600324085359098022255315e-1, -.190070260589904781867463408e-1,
        -.542644219965077966003735613e-2, .406730337639505305471483802e-10,
        .293370559238441772914955737e-1, -.289303678459884980822835164e-1,
        .104920404810938679338774661e-1, .275859250926960655799067950e-1,
        .181026916911544332889429182e-6, .325321795238610297791925570e-1,
        -.128946421695611868701754462e-1, -.155700419125490718558154235e-1,
        -.552044613055773574787347754e-2, -.218099108980700210189162121e-6 };

  // 105*xi^2*eta-60*xi*eta-15*xi^2+10*xi+5*eta-1
  static double w7[15] = {
        .472567169306760663658070139e-1,
        -.676818513328929399517743523e-1, .713599033825593937884471142e-1,
        .245817333396178085819438091e-1, -.175075129056328751103916923e-8,
        -.724344547436548992447909943e-2, -.695572540194121786964284835e-2,
        .157987544488321124145029297e-1, -.415370051759828139164174601e-1,
        -.581863021060588212948512507e-6, -.713289786794576036778561628e-1,
        .566322349972855693713038826e-2, .384238526960336452592664918e-1,
        -.833663002182206502225356168e-2, .354027868461204908497880664e-7 };

  // 105*xi*eta^2-60*xi*eta-158eta^2=5*xi+10*eta-1
  static double w8[15] = {
        -.807677222971670419500076218e-1,
        .937468289399066925339941225e-1, -.472748758921137151939194844e-1,
        -.412207701021053485512960996e-1, .387822703259552142802249597e-7,
        -.262675598921254856597215084e-1, .330207030089549993633697615e-1,
        .828627304161356941731930473e-2, .248979684134952739171531599e-1,
        .618894540095980139717698298e-6, .378179733129666280936555549e-1,
        .204017541072852182112868982e-1, -.143388252055879623028106437e-1,
        -.830240674066546796871400745e-2, .162873218927143591939772409e-8 };

  // 35*eta^3-45*eta^2+15*eta-1
  static double w9[15] = {
        .293370559238441772914955737e-1, -.289303678459884980822835164e-1,
        .104920404810938679338774661e-1, .275859250926960655799067950e-1,
        .181026916911544332889429182e-6, .325321795238610297791925570e-1,
        -.128946421695611868701754462e-1, -.155700419125490718558154235e-1,
        -.552044613055773574787347754e-2, -.218099108980700210189162121e-6,
        -.283582300812142314864875229e-1, .157600324085359098022255315e-1,
        -.190070260589904781867463408e-1, -.542644219965077966003735613e-2,
        .406730337639505305471483802e-10 };

  Functionals[0] =  PointValues[ 0]*w0[ 0] +PointValues[ 1]*w0[ 1]
                   +PointValues[ 2]*w0[ 2] +PointValues[ 3]*w0[ 3]
                   +PointValues[ 4]*w0[ 4] +PointValues[ 5]*w0[ 5]
                   +PointValues[ 6]*w0[ 6] +PointValues[ 7]*w0[ 7]
                   +PointValues[ 8]*w0[ 8] +PointValues[ 9]*w0[ 9]
                   +PointValues[10]*w0[10] +PointValues[11]*w0[11]
                   +PointValues[12]*w0[12] +PointValues[13]*w0[13]
                   +PointValues[14]*w0[14];
  Functionals[0] *= 2;
  Functionals[1] =  PointValues[ 0]*w1[ 0] +PointValues[ 1]*w1[ 1]
                   +PointValues[ 2]*w1[ 2] +PointValues[ 3]*w1[ 3]
                   +PointValues[ 4]*w1[ 4] +PointValues[ 5]*w1[ 5]
                   +PointValues[ 6]*w1[ 6] +PointValues[ 7]*w1[ 7]
                   +PointValues[ 8]*w1[ 8] +PointValues[ 9]*w1[ 9]
                   +PointValues[10]*w1[10] +PointValues[11]*w1[11]
                   +PointValues[12]*w1[12] +PointValues[13]*w1[13]
                   +PointValues[14]*w1[14];
  Functionals[2] =  PointValues[ 0]*w2[ 0] +PointValues[ 1]*w2[ 1]
                   +PointValues[ 2]*w2[ 2] +PointValues[ 3]*w2[ 3]
                   +PointValues[ 4]*w2[ 4] +PointValues[ 5]*w2[ 5]
                   +PointValues[ 6]*w2[ 6] +PointValues[ 7]*w2[ 7]
                   +PointValues[ 8]*w2[ 8] +PointValues[ 9]*w2[ 9]
                   +PointValues[10]*w2[10] +PointValues[11]*w2[11]
                   +PointValues[12]*w2[12] +PointValues[13]*w2[13]
                   +PointValues[14]*w2[14];
  Functionals[3] =  PointValues[ 0]*w3[ 0] +PointValues[ 1]*w3[ 1]
                   +PointValues[ 2]*w3[ 2] +PointValues[ 3]*w3[ 3]
                   +PointValues[ 4]*w3[ 4] +PointValues[ 5]*w3[ 5]
                   +PointValues[ 6]*w3[ 6] +PointValues[ 7]*w3[ 7]
                   +PointValues[ 8]*w3[ 8] +PointValues[ 9]*w3[ 9]
                   +PointValues[10]*w3[10] +PointValues[11]*w3[11]
                   +PointValues[12]*w3[12] +PointValues[13]*w3[13]
                   +PointValues[14]*w3[14];
  Functionals[4] =  PointValues[ 0]*w4[ 0] +PointValues[ 1]*w4[ 1]
                   +PointValues[ 2]*w4[ 2] +PointValues[ 3]*w4[ 3]
                   +PointValues[ 4]*w4[ 4] +PointValues[ 5]*w4[ 5]
                   +PointValues[ 6]*w4[ 6] +PointValues[ 7]*w4[ 7]
                   +PointValues[ 8]*w4[ 8] +PointValues[ 9]*w4[ 9]
                   +PointValues[10]*w4[10] +PointValues[11]*w4[11]
                   +PointValues[12]*w4[12] +PointValues[13]*w4[13]
                   +PointValues[14]*w4[14];
  Functionals[5] =  PointValues[ 0]*w5[ 0] +PointValues[ 1]*w5[ 1]
                   +PointValues[ 2]*w5[ 2] +PointValues[ 3]*w5[ 3]
                   +PointValues[ 4]*w5[ 4] +PointValues[ 5]*w5[ 5]
                   +PointValues[ 6]*w5[ 6] +PointValues[ 7]*w5[ 7]
                   +PointValues[ 8]*w5[ 8] +PointValues[ 9]*w5[ 9]
                   +PointValues[10]*w5[10] +PointValues[11]*w5[11]
                   +PointValues[12]*w5[12] +PointValues[13]*w5[13]
                   +PointValues[14]*w5[14];
  Functionals[6] =  PointValues[ 0]*w6[ 0] +PointValues[ 1]*w6[ 1]
                   +PointValues[ 2]*w6[ 2] +PointValues[ 3]*w6[ 3]
                   +PointValues[ 4]*w6[ 4] +PointValues[ 5]*w6[ 5]
                   +PointValues[ 6]*w6[ 6] +PointValues[ 7]*w6[ 7]
                   +PointValues[ 8]*w6[ 8] +PointValues[ 9]*w6[ 9]
                   +PointValues[10]*w6[10] +PointValues[11]*w6[11]
                   +PointValues[12]*w6[12] +PointValues[13]*w6[13]
                   +PointValues[14]*w6[14];
  Functionals[7] =  PointValues[ 0]*w7[ 0] +PointValues[ 1]*w7[ 1]
                   +PointValues[ 2]*w7[ 2] +PointValues[ 3]*w7[ 3]
                   +PointValues[ 4]*w7[ 4] +PointValues[ 5]*w7[ 5]
                   +PointValues[ 6]*w7[ 6] +PointValues[ 7]*w7[ 7]
                   +PointValues[ 8]*w7[ 8] +PointValues[ 9]*w7[ 9]
                   +PointValues[10]*w7[10] +PointValues[11]*w7[11]
                   +PointValues[12]*w7[12] +PointValues[13]*w7[13]
                   +PointValues[14]*w7[14];
  Functionals[8] =  PointValues[ 0]*w8[ 0] +PointValues[ 1]*w8[ 1]
                   +PointValues[ 2]*w8[ 2] +PointValues[ 3]*w8[ 3]
                   +PointValues[ 4]*w8[ 4] +PointValues[ 5]*w8[ 5]
                   +PointValues[ 6]*w8[ 6] +PointValues[ 7]*w8[ 7]
                   +PointValues[ 8]*w8[ 8] +PointValues[ 9]*w8[ 9]
                   +PointValues[10]*w8[10] +PointValues[11]*w8[11]
                   +PointValues[12]*w8[12] +PointValues[13]*w8[13]
                   +PointValues[14]*w8[14];
  Functionals[9] =  PointValues[ 0]*w9[ 0] +PointValues[ 1]*w9[ 1]
                   +PointValues[ 2]*w9[ 2] +PointValues[ 3]*w9[ 3]
                   +PointValues[ 4]*w9[ 4] +PointValues[ 5]*w9[ 5]
                   +PointValues[ 6]*w9[ 6] +PointValues[ 7]*w9[ 7]
                   +PointValues[ 8]*w9[ 8] +PointValues[ 9]*w9[ 9]
                   +PointValues[10]*w9[10] +PointValues[11]*w9[11]
                   +PointValues[12]*w9[12] +PointValues[13]*w9[13]
                   +PointValues[14]*w9[14];
}

void NF_D_T_P3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

TNodalFunctional2D *NF_D_T_P3_2D_Obj = new TNodalFunctional2D
        (NF_D_T_P3_2D, 10, 0, 15, 0, NF_D_T_P3_2D_Xi, NF_D_T_P3_2D_Eta,
         NULL, NF_D_T_P3_2D_EvalAll, NF_D_T_P3_2D_EvalEdge);
