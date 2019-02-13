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

static double NF_D_T_P4_2D_Xi[15] = {
    0.5133469206394541, 0.3132512106717253, 0.6517753036487957,
    0.6510199345893917E-1, 0.345792011168269, 0.2810412473151104,
    0.6306214343189561, 0.313477887523733, 0.8701651015635631,
    0.3623168221569262E1, 0.2056118320454355, 0.5612735500931855E-1,
    0.3474680882747129E-1, 0.6473290497749777E-1, -0.2968960232737531E1 };

static double NF_D_T_P4_2D_Eta[15] = {
    0.2810412473151104, 0.6306214343189561, 0.313477887523733,
    0.8701651015635631, 0.3623168221569262E1, 0.2056118320454355,
    0.5612735500931855E-1, 0.3474680882747129E-1, 0.6473290497749777E-1,
    -0.2968960232737531E1, 0.5133469206394541, 0.3132512106717253,
    0.6517753036487957, 0.6510199345893917E-1, 0.345792011168269 };


void NF_D_T_P4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
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

  // 126*xi^4-224*xi^3+126*xi^2-24*xi+1:
  static double w10[15] = {
        .221767308627192790731057838e-1, 
        .762966838390765433457705185e-2, -.116943725176883175552288990e-1,
        -.233386148664526976200066638e-2, .492695122875488653688224976e-10,
        .138956063914820570481510846e-2, -.121226872459650998707182898e-1,
        .511844273204934938973536974e-2, .458667993997063190145207201e-2,
        .202773605370327817929277680e-5, -.220771067553461410490535003e-1,
        .505936135607316659040987976e-3, .904893850824439171228863726e-2,
        -.223266044184371696482408914e-2, .270346051855358960039121004e-5 };

  // 504*xi^3*eta+1-18*xi-6*eta+63*xi^2+126*xi*eta-56*xi^3-504*xi^2*eta:
  static double w11[15] = {
        -.593850102393492470686138462e-1,
        -.214896763642209307535780894e-1, .297848223036231999841425693e-1,
        .686684127724466548920153410e-2, -.108070797989861472219013971e-8,
        -.109250599679565181918467558e-1, .275034022478316083366964602e-1,
        -.568908912708172747139696201e-2, -.943820270386623476471406851e-2,
        -.879445675195838933072176304e-5, .673317007142630778027258182e-1,
        .196043956928960775238882782e-2, -.290417506217523237324422229e-1,
        .253104540365830776994193708e-2, -.666954223600022711654758158e-6 };

  // 756*xi^2*eta^2+1-12*xi-12*eta+21*xi^2+21*eta^2+168*xi*eta-336*xi^2*eta-336*xi*eta^2:
  static double w12[15] = {
        .837479457077192586905313168e-2,
        .550128462585687444490579899e-2, .624727791375825138104382341e-3,
        -.244380311639749815591620249e-2, -.675663382735412659540126963e-8,
        .328126106397522462157027483e-1, -.427279864184728747971790117e-2,
        -.130186612196750903788218604e-1, -.164927464462249928130641159e-2,
        .142101972998281447898132270e-4, -.367198509709601408981537041e-1,
        -.131897341633600061349066994e-1, .198129595961155305009788046e-1,
        .415355179546487742060410306e-2, -.970314069332467756975875379e-8 };

  // 504*xi*eta^3+1-18*eta-6*xi+63*eta^2+126*xi*eta-56*eta^3-504*xi*eta^2:
  static double w13[15] = {
        .536681741136369537258632193e-2,
        .209873467360287924132811255e-1, -.147846818011156699401349118e-1,
        -.890851705601623543636072237e-2, .683512537145276613550655826e-6,
        .209767263071214863934881832e-1, -.398418411171887345200876298e-2,
        -.715400341122524246971515244e-2, .639959636371656028156548574e-2,
        -.101468878506143356899100820e-4, -.293219132115278692238092890e-1,
        -.902899717140968550023859389e-2, .169926677671300697778305007e-1,
        .246860466933641377722977765e-2, .883629930748419260714849723e-9 };

  // 126*eta^4-224*eta^3+126*eta^2-24*eta+1:
  static double w14[15] = {
        .138956063914820570481510846e-2,
        -.121226872459650998707182898e-1, .511844273204934938973536974e-2,
        .458667993997063190145207201e-2, .202773605370327817929277680e-5,
        -.220771067553461410490535003e-1, .505936135607316659040987976e-3,
        .904893850824439171228863726e-2, -.223266044184371696482408914e-2,
        .270346051855358960039121004e-5, .221767308627192790731057838e-1,
        .762966838390765433457705185e-2, -.116943725176883175552288990e-1,
        -.233386148664526976200066638e-2, .492695122875488653688224976e-10 };

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
  Functionals[10]=  PointValues[ 0]*w10[ 0] +PointValues[ 1]*w10[ 1]
                   +PointValues[ 2]*w10[ 2] +PointValues[ 3]*w10[ 3]
                   +PointValues[ 4]*w10[ 4] +PointValues[ 5]*w10[ 5]
                   +PointValues[ 6]*w10[ 6] +PointValues[ 7]*w10[ 7]
                   +PointValues[ 8]*w10[ 8] +PointValues[ 9]*w10[ 9]
                   +PointValues[10]*w10[10] +PointValues[11]*w10[11]
                   +PointValues[12]*w10[12] +PointValues[13]*w10[13]
                   +PointValues[14]*w10[14];
  Functionals[11]=  PointValues[ 0]*w11[ 0] +PointValues[ 1]*w11[ 1]
                   +PointValues[ 2]*w11[ 2] +PointValues[ 3]*w11[ 3]
                   +PointValues[ 4]*w11[ 4] +PointValues[ 5]*w11[ 5]
                   +PointValues[ 6]*w11[ 6] +PointValues[ 7]*w11[ 7]
                   +PointValues[ 8]*w11[ 8] +PointValues[ 9]*w11[ 9]
                   +PointValues[10]*w11[10] +PointValues[11]*w11[11]
                   +PointValues[12]*w11[12] +PointValues[13]*w11[13]
                   +PointValues[14]*w11[14];
  Functionals[12]=  PointValues[ 0]*w12[ 0] +PointValues[ 1]*w12[ 1]
                   +PointValues[ 2]*w12[ 2] +PointValues[ 3]*w12[ 3]
                   +PointValues[ 4]*w12[ 4] +PointValues[ 5]*w12[ 5]
                   +PointValues[ 6]*w12[ 6] +PointValues[ 7]*w12[ 7]
                   +PointValues[ 8]*w12[ 8] +PointValues[ 9]*w12[ 9]
                   +PointValues[10]*w12[10] +PointValues[11]*w12[11]
                   +PointValues[12]*w12[12] +PointValues[13]*w12[13]
                   +PointValues[14]*w12[14];
  Functionals[13]=  PointValues[ 0]*w13[ 0] +PointValues[ 1]*w13[ 1]
                   +PointValues[ 2]*w13[ 2] +PointValues[ 3]*w13[ 3]
                   +PointValues[ 4]*w13[ 4] +PointValues[ 5]*w13[ 5]
                   +PointValues[ 6]*w13[ 6] +PointValues[ 7]*w13[ 7]
                   +PointValues[ 8]*w13[ 8] +PointValues[ 9]*w13[ 9]
                   +PointValues[10]*w13[10] +PointValues[11]*w13[11]
                   +PointValues[12]*w13[12] +PointValues[13]*w13[13]
                   +PointValues[14]*w13[14];
  Functionals[14]=  PointValues[ 0]*w14[ 0] +PointValues[ 1]*w14[ 1]
                   +PointValues[ 2]*w14[ 2] +PointValues[ 3]*w14[ 3]
                   +PointValues[ 4]*w14[ 4] +PointValues[ 5]*w14[ 5]
                   +PointValues[ 6]*w14[ 6] +PointValues[ 7]*w14[ 7]
                   +PointValues[ 8]*w14[ 8] +PointValues[ 9]*w14[ 9]
                   +PointValues[10]*w14[10] +PointValues[11]*w14[11]
                   +PointValues[12]*w14[12] +PointValues[13]*w14[13]
                   +PointValues[14]*w14[14];
}

void NF_D_T_P4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

TNodalFunctional2D *NF_D_T_P4_2D_Obj = new TNodalFunctional2D
        (NF_D_T_P4_2D, 15, 0, 15, 0, NF_D_T_P4_2D_Xi, NF_D_T_P4_2D_Eta,
         NULL, NF_D_T_P4_2D_EvalAll, NF_D_T_P4_2D_EvalEdge);
