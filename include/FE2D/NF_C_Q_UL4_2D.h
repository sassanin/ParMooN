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
   
static double NF_C_Q_UL4_2D_Xi[] = {
-1., -.5000000000000000000000000, 0., .5000000000000000000000000, 1., 1., 1., 1., 1.,
.5000000000000000000000000, 0., -.5000000000000000000000000, -1., -1., -1., -1.,
-.9061798459386639927976267, -.5384693101056830910363143, 0., .5384693101056830910363143,
.9061798459386639927976267, -.9061798459386639927976267, -.5384693101056830910363143, 0.,
.5384693101056830910363143, .9061798459386639927976267, -.9061798459386639927976267,
-.5384693101056830910363143, 0., .5384693101056830910363143, .9061798459386639927976267,
-.9061798459386639927976267, -.5384693101056830910363143, 0., .5384693101056830910363143,
.9061798459386639927976267, -.9061798459386639927976267, -.5384693101056830910363143, 0.,
.5384693101056830910363143, .9061798459386639927976267
};

static double NF_C_Q_UL4_2D_Eta[] = {
-1., -1., -1., -1., -1., -.5000000000000000000000000, 0., .5000000000000000000000000, 1., 1.,
1., 1., 1., .5000000000000000000000000, 0., -.5000000000000000000000000, -.9061798459386639927976267,
-.9061798459386639927976267, -.9061798459386639927976267, -.9061798459386639927976267,
-.9061798459386639927976267, -.5384693101056830910363143, -.5384693101056830910363143,
-.5384693101056830910363143, -.5384693101056830910363143, -.5384693101056830910363143, 0., 0., 0., 0.,
0., .5384693101056830910363143, .5384693101056830910363143, .5384693101056830910363143,
.5384693101056830910363143, .5384693101056830910363143, .9061798459386639927976267,
.9061798459386639927976267, .9061798459386639927976267, .9061798459386639927976267,
.9061798459386639927976267
};

static double NF_C_Q_UL4_2D_T[] = {
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00
};

static double NF_C_Q_UL4_2D_W16[] = {
0.5613434886242863595465133e-1, .1134000000000000000000000, .1347850723875209031192266,
.1134000000000000000000001, 0.5613434886242863595465135e-1, .1134000000000000000000000,
.2290854042239911171317679, .2722865325507507018190467, .2290854042239911171317681,
.1134000000000000000000001, .1347850723875209031192266, .2722865325507507018190467,
.3236345679012345679012375, .2722865325507507018190469, .1347850723875209031192266,
.1134000000000000000000001, .2290854042239911171317681, .2722865325507507018190469,
.2290854042239911171317682, .1134000000000000000000001, 0.5613434886242863595465135e-1,
.1134000000000000000000001, .1347850723875209031192266, .1134000000000000000000001,
0.5613434886242863595465138e-1
};

static double NF_C_Q_UL4_2D_W17[] = {
-0.5086781560402279968924955e-1, -.1027607945294444967832509, -.1221395161309553661396283,
-.1027607945294444967832510, -0.5086781560402279968924957e-1, -0.6106241976598446252351804e-1,
-.1233554595677740359177672, -.1466179413336713527980412, -.1233554595677740359177673,
-0.6106241976598446252351810e-1, 0., 0., 0., 0., 0., 0.6106241976598446252351810e-1,
.1233554595677740359177673, .1466179413336713527980413, .1233554595677740359177673,
0.6106241976598446252351810e-1, 0.5086781560402279968924957e-1, .1027607945294444967832510,
.1221395161309553661396283, .1027607945294444967832510, 0.5086781560402279968924959e-1
};

static double NF_C_Q_UL4_2D_W18[] = {
0.4107590952972030537062780e-1, 0.8297964143284007959082117e-1, 0.9862801567209769736065559e-1,
0.8297964143284007959082124e-1, 0.4107590952972030537062781e-1, -0.7379641432840079590821313e-2,
-0.1490800829515240413605983e-1, -0.1771937369678905538534659e-1, -0.1490800829515240413605984e-1,
-0.7379641432840079590821319e-2, -0.6739253619376045155961330e-1, -.1361432662753753509095234,
-.1618172839506172839506188, -.1361432662753753509095234, -0.6739253619376045155961330e-1,
-0.7379641432840079590821319e-2, -0.1490800829515240413605984e-1, -0.1771937369678905538534660e-1,
-0.1490800829515240413605985e-1, -0.7379641432840079590821319e-2, 0.4107590952972030537062781e-1,
0.8297964143284007959082124e-1, 0.9862801567209769736065559e-1, 0.8297964143284007959082124e-1,
0.4107590952972030537062783e-1
};

static double NF_C_Q_UL4_2D_W19[] = {
-0.5086781560402279968924955e-1, -0.6106241976598446252351804e-1, 0.,
0.6106241976598446252351810e-1, 0.5086781560402279968924957e-1, -.1027607945294444967832509,
-.1233554595677740359177672, 0., .1233554595677740359177673, .1027607945294444967832510,
-.1221395161309553661396283, -.1466179413336713527980412, 0., .1466179413336713527980413,
.1221395161309553661396283, -.1027607945294444967832510, -.1233554595677740359177673, 0.,
.1233554595677740359177673, .1027607945294444967832510, -0.5086781560402279968924957e-1,
-0.6106241976598446252351810e-1, 0., 0.6106241976598446252351810e-1, 0.5086781560402279968924959e-1
};

static double NF_C_Q_UL4_2D_W20[] = {
0.4609538930728974889863567e-1, 0.5533353413618183126939358e-1, 0.,
-0.5533353413618183126939363e-1, -0.4609538930728974889863568e-1, 0.5533353413618183126939358e-1,
0.6642312921122876961988275e-1, 0., -0.6642312921122876961988280e-1, -0.5533353413618183126939363e-1,
-0., -0., 0., 0., 0., -0.5533353413618183126939363e-1, -0.6642312921122876961988280e-1, 0.,
0.6642312921122876961988283e-1, 0.5533353413618183126939363e-1, -0.4609538930728974889863568e-1,
-0.5533353413618183126939363e-1, 0., 0.5533353413618183126939363e-1, 0.4609538930728974889863571e-1
};

static double NF_C_Q_UL4_2D_W21[] = {
-0.3722216136943244646106953e-1, -0.4468199027515835399732213e-1, 0.,
0.4468199027515835399732217e-1, 0.3722216136943244646106954e-1, 0.6687282336693604926162598e-2,
0.8027504941740515797812026e-2, 0., -0.8027504941740515797812033e-2, -0.6687282336693604926162603e-2,
0.6106975806547768306981415e-1, 0.7330897066683567639902060e-1, 0., -0.7330897066683567639902065e-1,
-0.6106975806547768306981415e-1, 0.6687282336693604926162603e-2, 0.8027504941740515797812033e-2, 0.,
-0.8027504941740515797812037e-2, -0.6687282336693604926162603e-2, -0.3722216136943244646106954e-1,
-0.4468199027515835399732217e-1, 0., 0.4468199027515835399732217e-1, 0.3722216136943244646106956e-1
};

static double NF_C_Q_UL4_2D_W22[] = {
0.4107590952972030537062780e-1, -0.7379641432840079590821313e-2,
-0.6739253619376045155961330e-1, -0.7379641432840079590821319e-2, 0.4107590952972030537062781e-1,
0.8297964143284007959082117e-1, -0.1490800829515240413605983e-1, -.1361432662753753509095234,
-0.1490800829515240413605984e-1, 0.8297964143284007959082124e-1, 0.9862801567209769736065559e-1,
-0.1771937369678905538534659e-1, -.1618172839506172839506188, -0.1771937369678905538534660e-1,
0.9862801567209769736065559e-1, 0.8297964143284007959082124e-1, -0.1490800829515240413605984e-1,
-.1361432662753753509095234, -0.1490800829515240413605985e-1, 0.8297964143284007959082124e-1,
0.4107590952972030537062781e-1, -0.7379641432840079590821319e-2, -0.6739253619376045155961330e-1,
-0.7379641432840079590821319e-2, 0.4107590952972030537062783e-1
};

static double NF_C_Q_UL4_2D_W23[] = {
-0.3722216136943244646106951e-1, 0.6687282336693604926162588e-2,
0.6106975806547768306981415e-1, 0.6687282336693604926162594e-2, -0.3722216136943244646106952e-1,
-0.4468199027515835399732212e-1, 0.8027504941740515797812004e-2, 0.7330897066683567639902060e-1,
0.8027504941740515797812011e-2, -0.4468199027515835399732216e-1, 0., -0., -0., -0., 0.,
0.4468199027515835399732216e-1, -0.8027504941740515797812011e-2, -0.7330897066683567639902065e-1,
-0.8027504941740515797812014e-2, 0.4468199027515835399732216e-1, 0.3722216136943244646106952e-1,
-0.6687282336693604926162594e-2, -0.6106975806547768306981415e-1, -0.6687282336693604926162594e-2,
0.3722216136943244646106954e-1
};

static double NF_C_Q_UL4_2D_W24[] = {
0.3005700391802442434016381e-1, -0.5400000000000000000000013e-2,
-0.4931400783604884868032782e-1, -0.5400000000000000000000018e-2, 0.3005700391802442434016382e-1,
-0.5400000000000000000000006e-2, 0.9701565758027361536633835e-3, 0.8859686848394527692673296e-2,
0.9701565758027361536633844e-3, -0.5400000000000000000000011e-2, -0.4931400783604884868032779e-1,
0.8859686848394527692673296e-2, 0.8090864197530864197530938e-1, 0.8859686848394527692673302e-2,
-0.4931400783604884868032779e-1, -0.5400000000000000000000011e-2, 0.9701565758027361536633844e-3,
0.8859686848394527692673302e-2, 0.9701565758027361536633848e-3, -0.5400000000000000000000011e-2,
0.3005700391802442434016382e-1, -0.5400000000000000000000018e-2, -0.4931400783604884868032782e-1,
-0.5400000000000000000000018e-2, 0.3005700391802442434016383e-1
};

static double NF_C_Q_UL4_2D_W25[] = {
-0.2812505854637221097561622e-1, -0.5681693479646467344727129e-1,
-0.6753152265765913792627117e-1, -0.5681693479646467344727134e-1, -0.2812505854637221097561623e-1,
0.4733113056260416220730673e-1, 0.9561614794808355027486474e-1, .1136475257724722186193667,
0.9561614794808355027486483e-1, 0.4733113056260416220730677e-1, 0., 0., 0., 0., 0.,
-0.4733113056260416220730677e-1, -0.9561614794808355027486483e-1, -.1136475257724722186193668,
-0.9561614794808355027486487e-1, -0.4733113056260416220730677e-1, 0.2812505854637221097561623e-1,
0.5681693479646467344727134e-1, 0.6753152265765913792627117e-1, 0.5681693479646467344727134e-1,
0.2812505854637221097561625e-1
};

static double NF_C_Q_UL4_2D_W26[] = {
-0.2812505854637221097561622e-1, 0.4733113056260416220730673e-1, 0.,
-0.4733113056260416220730677e-1, 0.2812505854637221097561623e-1, -0.5681693479646467344727129e-1,
0.9561614794808355027486474e-1, 0., -0.9561614794808355027486483e-1, 0.5681693479646467344727134e-1,
-0.6753152265765913792627117e-1, .1136475257724722186193667, 0., -.1136475257724722186193668,
0.6753152265765913792627117e-1, -0.5681693479646467344727134e-1, 0.9561614794808355027486483e-1, 0.,
-0.9561614794808355027486487e-1, 0.5681693479646467344727134e-1, -0.2812505854637221097561623e-1,
0.4733113056260416220730677e-1, 0., -0.4733113056260416220730677e-1, 0.2812505854637221097561625e-1
};


void NF_C_Q_UL4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // point values on the boundary: 0 - 15
  memcpy(Functionals, PointValues, 16*SizeOfDouble);
  Functionals[16] = NF_C_Q_UL4_2D_W16[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W16[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W16[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W16[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W16[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W16[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W16[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W16[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W16[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W16[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W16[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W16[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W16[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W16[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W16[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W16[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W16[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W16[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W16[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W16[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W16[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W16[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W16[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W16[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W16[24]*PointValues[40];
  Functionals[17] = NF_C_Q_UL4_2D_W17[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W17[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W17[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W17[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W17[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W17[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W17[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W17[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W17[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W17[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W17[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W17[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W17[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W17[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W17[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W17[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W17[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W17[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W17[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W17[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W17[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W17[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W17[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W17[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W17[24]*PointValues[40];
  Functionals[18] = NF_C_Q_UL4_2D_W18[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W18[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W18[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W18[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W18[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W18[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W18[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W18[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W18[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W18[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W18[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W18[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W18[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W18[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W18[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W18[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W18[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W18[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W18[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W18[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W18[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W18[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W18[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W18[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W18[24]*PointValues[40];
  Functionals[19] = NF_C_Q_UL4_2D_W19[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W19[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W19[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W19[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W19[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W19[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W19[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W19[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W19[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W19[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W19[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W19[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W19[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W19[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W19[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W19[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W19[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W19[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W19[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W19[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W19[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W19[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W19[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W19[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W19[24]*PointValues[40];
  Functionals[20] = NF_C_Q_UL4_2D_W20[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W20[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W20[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W20[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W20[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W20[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W20[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W20[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W20[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W20[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W20[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W20[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W20[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W20[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W20[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W20[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W20[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W20[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W20[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W20[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W20[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W20[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W20[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W20[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W20[24]*PointValues[40];
  Functionals[21] = NF_C_Q_UL4_2D_W21[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W21[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W21[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W21[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W21[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W21[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W21[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W21[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W21[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W21[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W21[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W21[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W21[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W21[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W21[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W21[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W21[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W21[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W21[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W21[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W21[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W21[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W21[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W21[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W21[24]*PointValues[40];
  Functionals[22] = NF_C_Q_UL4_2D_W22[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W22[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W22[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W22[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W22[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W22[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W22[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W22[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W22[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W22[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W22[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W22[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W22[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W22[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W22[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W22[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W22[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W22[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W22[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W22[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W22[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W22[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W22[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W22[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W22[24]*PointValues[40];
  Functionals[23] = NF_C_Q_UL4_2D_W23[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W23[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W23[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W23[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W23[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W23[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W23[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W23[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W23[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W23[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W23[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W23[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W23[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W23[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W23[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W23[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W23[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W23[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W23[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W23[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W23[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W23[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W23[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W23[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W23[24]*PointValues[40];
  Functionals[24] = NF_C_Q_UL4_2D_W24[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W24[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W24[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W24[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W24[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W24[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W24[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W24[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W24[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W24[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W24[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W24[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W24[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W24[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W24[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W24[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W24[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W24[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W24[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W24[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W24[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W24[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W24[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W24[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W24[24]*PointValues[40];
  Functionals[25] = NF_C_Q_UL4_2D_W25[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W25[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W25[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W25[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W25[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W25[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W25[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W25[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W25[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W25[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W25[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W25[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W25[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W25[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W25[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W25[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W25[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W25[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W25[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W25[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W25[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W25[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W25[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W25[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W25[24]*PointValues[40];
  Functionals[26] = NF_C_Q_UL4_2D_W26[ 0]*PointValues[16]
                   +NF_C_Q_UL4_2D_W26[ 1]*PointValues[17]
                   +NF_C_Q_UL4_2D_W26[ 2]*PointValues[18]
                   +NF_C_Q_UL4_2D_W26[ 3]*PointValues[19]
                   +NF_C_Q_UL4_2D_W26[ 4]*PointValues[20]
                   +NF_C_Q_UL4_2D_W26[ 5]*PointValues[21]
                   +NF_C_Q_UL4_2D_W26[ 6]*PointValues[22]
                   +NF_C_Q_UL4_2D_W26[ 7]*PointValues[23]
                   +NF_C_Q_UL4_2D_W26[ 8]*PointValues[24]
                   +NF_C_Q_UL4_2D_W26[ 9]*PointValues[25]
                   +NF_C_Q_UL4_2D_W26[10]*PointValues[26]
                   +NF_C_Q_UL4_2D_W26[11]*PointValues[27]
                   +NF_C_Q_UL4_2D_W26[12]*PointValues[28]
                   +NF_C_Q_UL4_2D_W26[13]*PointValues[29]
                   +NF_C_Q_UL4_2D_W26[14]*PointValues[30]
                   +NF_C_Q_UL4_2D_W26[15]*PointValues[31]
                   +NF_C_Q_UL4_2D_W26[16]*PointValues[32]
                   +NF_C_Q_UL4_2D_W26[17]*PointValues[33]
                   +NF_C_Q_UL4_2D_W26[18]*PointValues[34]
                   +NF_C_Q_UL4_2D_W26[19]*PointValues[35]
                   +NF_C_Q_UL4_2D_W26[20]*PointValues[36]
                   +NF_C_Q_UL4_2D_W26[21]*PointValues[37]
                   +NF_C_Q_UL4_2D_W26[22]*PointValues[38]
                   +NF_C_Q_UL4_2D_W26[23]*PointValues[39]
                   +NF_C_Q_UL4_2D_W26[24]*PointValues[40];
};

void NF_C_Q_UL4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  memcpy(Functionals, PointValues, 5*SizeOfDouble);
};

/*
TNodalFunctional2D(NodalFunctional2D id,
         int n_allfunctionals, int n_edgefunctionals, 
         int n_pointsall, int n_pointsedge, 
         double *xi, double *eta, double *t, 
         DoubleFunctVect *evalall, 
         DoubleFunctVect *evaledge); 
*/

TNodalFunctional2D *NF_C_Q_UL4_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL4_2D, 27, 5, 41, 5, NF_C_Q_UL4_2D_Xi, NF_C_Q_UL4_2D_Eta,
         NF_C_Q_UL4_2D_T, NF_C_Q_UL4_2D_EvalAll, NF_C_Q_UL4_2D_EvalEdge);
