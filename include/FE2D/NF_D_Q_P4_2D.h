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
   
static double NF_D_Q_P4_2D_Xi[25]={
-.9061798459386639927976269,
-.5384693101056830910363144,
0,
.5384693101056830910363144,
.9061798459386639927976269,
-.9061798459386639927976269,
-.5384693101056830910363144,
0,
.5384693101056830910363144,
.9061798459386639927976269,
-.9061798459386639927976269,
-.5384693101056830910363144,
0,
.5384693101056830910363144,
.9061798459386639927976269,
-.9061798459386639927976269,
-.5384693101056830910363144,
0,
.5384693101056830910363144,
.9061798459386639927976269,
-.9061798459386639927976269,
-.5384693101056830910363144,
0,
.5384693101056830910363144,
.9061798459386639927976269
};
 
static double NF_D_Q_P4_2D_Eta[25]={
-.9061798459386639927976269,
-.9061798459386639927976269,
-.9061798459386639927976269,
-.9061798459386639927976269,
-.9061798459386639927976269,
-.5384693101056830910363144,
-.5384693101056830910363144,
-.5384693101056830910363144,
-.5384693101056830910363144,
-.5384693101056830910363144,
 0, 0, 0, 0, 0,
.5384693101056830910363144,
.5384693101056830910363144,
.5384693101056830910363144,
.5384693101056830910363144,
.5384693101056830910363144,
.9061798459386639927976269,
.9061798459386639927976269,
.9061798459386639927976269,
.9061798459386639927976269,
.9061798459386639927976269
};

static double NF_D_Q_P4_2D_W0[25]={
.1403358721560715898866278e-1,
.2835000000000000000000000e-1,
.3369626809688022577980643e-1,
.2835000000000000000000000e-1,
.1403358721560715898866278e-1,
.2835000000000000000000000e-1,
.5727135105599777928294217e-1,
.6807163313768767545476146e-1,
.5727135105599777928294217e-1,
.2835000000000000000000000e-1,
.3369626809688022577980643e-1,
.6807163313768767545476146e-1,
.8090864197530864197530862e-1,
.6807163313768767545476146e-1,
.3369626809688022577980643e-1,
.2835000000000000000000000e-1,
.5727135105599777928294217e-1,
.6807163313768767545476146e-1,
.5727135105599777928294217e-1,
.2835000000000000000000000e-1,
.1403358721560715898866278e-1,
.2835000000000000000000000e-1,
.3369626809688022577980643e-1,
.2835000000000000000000000e-1,
.1403358721560715898866278e-1
};

static double NF_D_Q_P4_2D_W1[25]={
-.01271695390100569992231234,
-.01526560494149611563087951,
0.,
.01526560494149611563087951,
.01271695390100569992231234,
-.02569019863236112419581272,
-.03083886489194350897944190,
0.,
.03083886489194350897944190,
.02569019863236112419581272,
-.03053487903273884153490688,
-.03665448533341783819951018,
0.,
.03665448533341783819951018,
.03053487903273884153490688,
-.02569019863236112419581272,
-.03083886489194350897944190,
0.,
.03083886489194350897944190,
.02569019863236112419581272,
-.01271695390100569992231234,
-.01526560494149611563087951,
0.,
.01526560494149611563087951,
.01271695390100569992231234
};

static double NF_D_Q_P4_2D_W2[25]={
-.01271695390100569992231234,
-.02569019863236112419581272,
-.03053487903273884153490688,
-.02569019863236112419581272,
-.01271695390100569992231234,
-.01526560494149611563087951,
-.03083886489194350897944190,
-.03665448533341783819951018,
-.03083886489194350897944190,
-.01526560494149611563087951,
0.,
0.,
0.,
0.,
0.,
.01526560494149611563087951,
.03083886489194350897944190,
.03665448533341783819951018,
.03083886489194350897944190,
.01526560494149611563087951,
.01271695390100569992231234,
.02569019863236112419581272,
.03053487903273884153490688,
.02569019863236112419581272,
.01271695390100569992231234
};

static double NF_D_Q_P4_2D_W3[25]={
.02053795476486015268531385,
-.003689820716420039795410648,
-.03369626809688022577980643,
-.003689820716420039795410648,
.02053795476486015268531385,
.04148982071642003979541064,
-.007454004147576202068029924,
-.06807163313768767545476146,
-.007454004147576202068029924,
.04148982071642003979541064,
.04931400783604884868032754,
-.008859686848394527692673248,
-.08090864197530864197530862,
-.008859686848394527692673248,
.04931400783604884868032754,
.04148982071642003979541064,
-.007454004147576202068029924,
-.06807163313768767545476146,
-.007454004147576202068029924,
.04148982071642003979541064,
.02053795476486015268531385,
-.003689820716420039795410648,
-.03369626809688022577980643,
-.003689820716420039795410648,
.02053795476486015268531385
};

static double NF_D_Q_P4_2D_W4[25]={
.01152384732682243722465888,
.01383338353404545781734840,
0.,
-.01383338353404545781734840,
-.01152384732682243722465888,
.01383338353404545781734840,
.01660578230280719240497075,
0.,
-.01660578230280719240497075,
-.01383338353404545781734840,
0.,
0.,
0.,
0.,
0.,
-.01383338353404545781734840,
-.01660578230280719240497075,
0.,
.01660578230280719240497075,
.01383338353404545781734840,
-.01152384732682243722465888,
-.01383338353404545781734840,
0.,
.01383338353404545781734840,
.01152384732682243722465888
};

static double NF_D_Q_P4_2D_W5[25]={
.02053795476486015268531385,
.04148982071642003979541064,
.04931400783604884868032754,
.04148982071642003979541064,
.02053795476486015268531385,
-.003689820716420039795410648,
-.007454004147576202068029924,
-.008859686848394527692673248,
-.007454004147576202068029924,
-.003689820716420039795410648,
-.03369626809688022577980643,
-.06807163313768767545476146,
-.08090864197530864197530862,
-.06807163313768767545476146,
-.03369626809688022577980643,
-.003689820716420039795410648,
-.007454004147576202068029924,
-.008859686848394527692673248,
-.007454004147576202068029924,
-.003689820716420039795410648,
.02053795476486015268531385,
.04148982071642003979541064,
.04931400783604884868032754,
.04148982071642003979541064,
.02053795476486015268531385
};

static double NF_D_Q_P4_2D_W6[25]={
-.01406252927318610548780807,
.02366556528130208110365336,
0.,
-.02366556528130208110365336,
.01406252927318610548780807,
-.02840846739823233672363567,
.04780807397404177513743252,
0.,
-.04780807397404177513743252,
.02840846739823233672363567,
-.03376576132882956896313540,
.05682376288623610930968317,
0.,
-.05682376288623610930968317,
.03376576132882956896313540,
-.02840846739823233672363567,
.04780807397404177513743252,
0.,
-.04780807397404177513743252,
.02840846739823233672363567,
-.01406252927318610548780807,
.02366556528130208110365336,
0.,
-.02366556528130208110365336,
.01406252927318610548780807
};

static double NF_D_Q_P4_2D_W7[25]={
-.01861108068471622323053472,
.003343641168346802463081287,
.03053487903273884153490688,
.003343641168346802463081287,
-.01861108068471622323053472,
-.02234099513757917699866109,
.004013752470870257898906007,
.03665448533341783819951018,
.004013752470870257898906007,
-.02234099513757917699866109,
0.,
0.,
0.,
0.,
0.,
.02234099513757917699866109,
-.004013752470870257898906007,
-.03665448533341783819951018,
-.004013752470870257898906007,
.02234099513757917699866109,
.01861108068471622323053472,
-.003343641168346802463081287,
-.03053487903273884153490688,
-.003343641168346802463081287,
.01861108068471622323053472
};

static double NF_D_Q_P4_2D_W8[25]={
-.01861108068471622323053472,
-.02234099513757917699866108,
0.,
.02234099513757917699866108,
.01861108068471622323053472,
.003343641168346802463081285,
.004013752470870257898906010,
0.,
-.004013752470870257898906010,
-.003343641168346802463081285,
.03053487903273884153490688,
.03665448533341783819951018,
0.,
-.03665448533341783819951018,
-.03053487903273884153490688,
.003343641168346802463081285,
.004013752470870257898906010,
0.,
-.004013752470870257898906010,
-.003343641168346802463081285,
-.01861108068471622323053472,
-.02234099513757917699866108,
0.,
.02234099513757917699866108,
.01861108068471622323053472
};

static double NF_D_Q_P4_2D_W9[25]={
-.01406252927318610548780807,
-.02840846739823233672363567,
-.03376576132882956896313540,
-.02840846739823233672363567,
-.01406252927318610548780807,
.02366556528130208110365336,
.04780807397404177513743252,
.05682376288623610930968317,
.04780807397404177513743252,
.02366556528130208110365336,
0.,
0.,
0.,
0.,
0.,
-.02366556528130208110365336,
-.04780807397404177513743252,
-.05682376288623610930968317,
-.04780807397404177513743252,
-.02366556528130208110365336,
.01406252927318610548780807,
.02840846739823233672363567,
.03376576132882956896313540,
.02840846739823233672363567,
.01406252927318610548780807
};

static double NF_D_Q_P4_2D_W10[25]={
.02758839997740570516854408,
-.07813280212272604383825375,
.1010888042906406773394193,
-.07813280212272604383825375,
.02758839997740570516854408,
.05573280212272604383825376,
-.1578402518292575570203960,
.2042148994130630263642844,
-.1578402518292575570203960,
.05573280212272604383825376,
.06624294328457679424142933,
-.1876059062475397572043923,
.2427259259259259259259259,
-.1876059062475397572043923,
.06624294328457679424142933,
.05573280212272604383825376,
-.1578402518292575570203960,
.2042148994130630263642844,
-.1578402518292575570203960,
.05573280212272604383825376,
.02758839997740570516854408,
-.07813280212272604383825375,
.1010888042906406773394193,
-.07813280212272604383825375,
.02758839997740570516854408
};

static double NF_D_Q_P4_2D_W11[25]={
.01274318061028373760349795,
-.02144525830066171525144593,
0.,
.02144525830066171525144593,
-.01274318061028373760349795,
.01529708784108595622151332,
-.02574318061028373760349795,
0.,
.02574318061028373760349795,
-.01529708784108595622151332,
0.,
0.,
0.,
0.,
0.,
-.01529708784108595622151332,
.02574318061028373760349795,
0.,
-.02574318061028373760349795,
.01529708784108595622151332,
-.01274318061028373760349795,
.02144525830066171525144593,
0.,
-.02144525830066171525144593,
.01274318061028373760349795
};

static double NF_D_Q_P4_2D_W12[25]={
.03005700391802442434016376,
-.005400000000000000000000006,
-.04931400783604884868032754,
-.005400000000000000000000006,
.03005700391802442434016376,
-.005399999999999999999999995,
.0009701565758027361536633800,
.008859686848394527692673248,
.0009701565758027361536633800,
-.005399999999999999999999995,
-.04931400783604884868032754,
.008859686848394527692673248,
.08090864197530864197530862,
.008859686848394527692673248,
-.04931400783604884868032754,
-.005399999999999999999999995,
.0009701565758027361536633800,
.008859686848394527692673248,
.0009701565758027361536633800,
-.005399999999999999999999995,
.03005700391802442434016376,
-.005400000000000000000000006,
-.04931400783604884868032754,
-.005400000000000000000000006,
.03005700391802442434016376
};

static double NF_D_Q_P4_2D_W13[25]={
.01274318061028373760349794,
.01529708784108595622151332,
0.,
-.01529708784108595622151332,
-.01274318061028373760349794,
-.02144525830066171525144594,
-.02574318061028373760349795,
0.,
.02574318061028373760349795,
.02144525830066171525144594,
0.,
0.,
0.,
0.,
0.,
.02144525830066171525144594,
.02574318061028373760349795,
0.,
-.02574318061028373760349795,
-.02144525830066171525144594,
-.01274318061028373760349794,
-.01529708784108595622151332,
0.,
.01529708784108595622151332,
.01274318061028373760349794
};

static double NF_D_Q_P4_2D_W14[25]={
.02758839997740570516854408,
.05573280212272604383825376,
.06624294328457679424142933,
.05573280212272604383825376,
.02758839997740570516854408,
-.07813280212272604383825375,
-.1578402518292575570203960,
-.1876059062475397572043923,
-.1578402518292575570203960,
-.07813280212272604383825375,
.1010888042906406773394193,
.2042148994130630263642844,
.2427259259259259259259259,
.2042148994130630263642844,
.1010888042906406773394193,
-.07813280212272604383825375,
-.1578402518292575570203960,
-.1876059062475397572043923,
-.1578402518292575570203960,
-.07813280212272604383825375,
.02758839997740570516854408,
.05573280212272604383825376,
.06624294328457679424142933,
.05573280212272604383825376,
.02758839997740570516854408
};

static double *NF_D_Q_P4_2D_T = NULL;

void NF_D_Q_P4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] =  NF_D_Q_P4_2D_W0[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W0[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W0[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W0[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W0[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W0[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W0[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W0[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W0[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W0[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W0[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W0[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W0[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W0[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W0[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W0[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W0[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W0[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W0[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W0[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W0[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W0[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W0[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W0[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W0[24]*PointValues[24];
  Functionals[1] =  NF_D_Q_P4_2D_W1[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W1[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W1[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W1[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W1[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W1[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W1[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W1[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W1[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W1[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W1[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W1[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W1[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W1[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W1[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W1[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W1[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W1[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W1[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W1[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W1[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W1[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W1[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W1[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W1[24]*PointValues[24];
  Functionals[2] =  NF_D_Q_P4_2D_W2[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W2[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W2[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W2[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W2[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W2[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W2[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W2[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W2[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W2[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W2[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W2[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W2[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W2[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W2[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W2[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W2[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W2[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W2[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W2[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W2[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W2[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W2[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W2[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W2[24]*PointValues[24];
  Functionals[3] =  NF_D_Q_P4_2D_W3[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W3[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W3[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W3[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W3[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W3[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W3[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W3[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W3[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W3[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W3[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W3[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W3[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W3[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W3[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W3[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W3[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W3[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W3[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W3[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W3[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W3[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W3[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W3[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W3[24]*PointValues[24];
  Functionals[4] =  NF_D_Q_P4_2D_W4[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W4[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W4[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W4[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W4[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W4[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W4[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W4[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W4[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W4[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W4[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W4[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W4[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W4[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W4[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W4[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W4[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W4[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W4[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W4[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W4[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W4[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W4[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W4[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W4[24]*PointValues[24];
  Functionals[5] =  NF_D_Q_P4_2D_W5[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W5[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W5[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W5[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W5[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W5[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W5[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W5[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W5[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W5[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W5[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W5[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W5[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W5[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W5[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W5[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W5[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W5[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W5[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W5[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W5[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W5[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W5[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W5[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W5[24]*PointValues[24];
  Functionals[6] =  NF_D_Q_P4_2D_W6[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W6[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W6[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W6[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W6[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W6[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W6[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W6[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W6[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W6[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W6[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W6[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W6[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W6[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W6[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W6[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W6[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W6[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W6[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W6[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W6[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W6[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W6[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W6[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W6[24]*PointValues[24];
  Functionals[7] =  NF_D_Q_P4_2D_W7[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W7[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W7[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W7[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W7[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W7[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W7[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W7[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W7[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W7[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W7[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W7[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W7[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W7[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W7[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W7[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W7[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W7[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W7[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W7[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W7[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W7[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W7[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W7[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W7[24]*PointValues[24];
  Functionals[8] =  NF_D_Q_P4_2D_W8[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W8[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W8[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W8[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W8[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W8[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W8[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W8[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W8[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W8[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W8[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W8[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W8[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W8[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W8[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W8[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W8[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W8[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W8[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W8[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W8[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W8[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W8[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W8[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W8[24]*PointValues[24];
  Functionals[9] =  NF_D_Q_P4_2D_W9[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W9[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W9[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W9[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W9[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W9[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W9[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W9[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W9[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W9[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W9[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W9[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W9[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W9[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W9[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W9[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W9[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W9[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W9[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W9[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W9[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W9[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W9[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W9[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W9[24]*PointValues[24];
  Functionals[10] = NF_D_Q_P4_2D_W10[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W10[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W10[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W10[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W10[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W10[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W10[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W10[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W10[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W10[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W10[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W10[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W10[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W10[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W10[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W10[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W10[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W10[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W10[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W10[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W10[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W10[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W10[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W10[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W10[24]*PointValues[24];
  Functionals[11] = NF_D_Q_P4_2D_W11[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W11[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W11[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W11[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W11[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W11[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W11[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W11[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W11[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W11[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W11[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W11[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W11[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W11[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W11[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W11[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W11[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W11[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W11[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W11[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W11[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W11[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W11[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W11[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W11[24]*PointValues[24];
  Functionals[12] = NF_D_Q_P4_2D_W12[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W12[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W12[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W12[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W12[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W12[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W12[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W12[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W12[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W12[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W12[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W12[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W12[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W12[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W12[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W12[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W12[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W12[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W12[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W12[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W12[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W12[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W12[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W12[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W12[24]*PointValues[24];
  Functionals[13] = NF_D_Q_P4_2D_W13[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W13[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W13[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W13[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W13[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W13[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W13[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W13[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W13[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W13[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W13[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W13[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W13[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W13[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W13[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W13[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W13[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W13[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W13[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W13[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W13[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W13[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W13[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W13[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W13[24]*PointValues[24];
  Functionals[14] = NF_D_Q_P4_2D_W14[0]*PointValues[0]
                   +NF_D_Q_P4_2D_W14[1]*PointValues[1]
                   +NF_D_Q_P4_2D_W14[2]*PointValues[2]
                   +NF_D_Q_P4_2D_W14[3]*PointValues[3]
                   +NF_D_Q_P4_2D_W14[4]*PointValues[4]
                   +NF_D_Q_P4_2D_W14[5]*PointValues[5]
                   +NF_D_Q_P4_2D_W14[6]*PointValues[6]
                   +NF_D_Q_P4_2D_W14[7]*PointValues[7]
                   +NF_D_Q_P4_2D_W14[8]*PointValues[8]
                   +NF_D_Q_P4_2D_W14[9]*PointValues[9]
                   +NF_D_Q_P4_2D_W14[10]*PointValues[10]
                   +NF_D_Q_P4_2D_W14[11]*PointValues[11]
                   +NF_D_Q_P4_2D_W14[12]*PointValues[12]
                   +NF_D_Q_P4_2D_W14[13]*PointValues[13]
                   +NF_D_Q_P4_2D_W14[14]*PointValues[14]
                   +NF_D_Q_P4_2D_W14[15]*PointValues[15]
                   +NF_D_Q_P4_2D_W14[16]*PointValues[16]
                   +NF_D_Q_P4_2D_W14[17]*PointValues[17]
                   +NF_D_Q_P4_2D_W14[18]*PointValues[18]
                   +NF_D_Q_P4_2D_W14[19]*PointValues[19]
                   +NF_D_Q_P4_2D_W14[20]*PointValues[20]
                   +NF_D_Q_P4_2D_W14[21]*PointValues[21]
                   +NF_D_Q_P4_2D_W14[22]*PointValues[22]
                   +NF_D_Q_P4_2D_W14[23]*PointValues[23]
                   +NF_D_Q_P4_2D_W14[24]*PointValues[24];
}

void NF_D_Q_P4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
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

TNodalFunctional2D *NF_D_Q_P4_2D_Obj = new TNodalFunctional2D
        (NF_D_Q_P4_2D, 15, 0, 25, 0, NF_D_Q_P4_2D_Xi, NF_D_Q_P4_2D_Eta,
         NF_D_Q_P4_2D_T, NF_D_Q_P4_2D_EvalAll, NULL);
