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

static double NF_C_T_UL4_2D_W0[27] = {
                  0.01365973100267786195871938,
                  0.01365973100267786195871938,
                  0.01365973100267786195871938,
                  0.03618454050341807933828794,
                  0.03618454050341807933828794,
                  0.03618454050341807933828794,
                 0.0009270063289606760506600850,
                 0.0009270063289606760506600850,
                 0.0009270063289606760506600850,
                  0.05932297738077407311581238,
                  0.05932297738077407311581238,
                  0.05932297738077407311581238,
                  0.07714953491481312286151106,
                  0.07714953491481312286151106,
                  0.07714953491481312286151106,
                  0.05233711196220407111713894,
                  0.05233711196220407111713894,
                  0.05233711196220407111713894,
                  0.05233711196220407111713894,
                  0.05233711196220407111713894,
                  0.05233711196220407111713894,
                  0.02070765963914068888703230,
                  0.02070765963914068888703230,
                  0.02070765963914068888703230,
                  0.02070765963914068888703230,
                  0.02070765963914068888703230,
                  0.02070765963914068888703230 };


static double NF_C_T_UL4_2D_W1[27] = {
                 -0.006166720773665447533343113,
                               0.0,
                 0.006166720773665447533343115,
                 -0.01161428337229559707660497,
                               0.0,
                  0.01161428337229559707660498,
                 0.0002798785734386036356090761,
                              0.0,
                -0.0002798785734386036356090761,
                 -0.01157006575352843784923587,
                               0.0,
                  0.01157006575352843784923588,
                 0.007595676796247792031257549,
                              0.0,
                 -0.007595676796247792031257541,
                 0.006190971763611308709464438,
                 -0.01421012974024818159828179,
                 0.008019157976636872888817360,
                  0.01421012974024818159828180,
                 -0.008019157976636872888817360,
                 -0.006190971763611308709464433,
                 0.006585036001891291212813170,
                 -0.008142945338211328594027984,
                 0.001557909336320037381214818,
                 0.008142945338211328594027987,
                 -0.001557909336320037381214816,
                 -0.006585036001891291212813166 };


static double NF_C_T_UL4_2D_W2[27] = {
                 -0.006166720773665447533343113,
                 0.006166720773665447533343115,
                               0.0,
                 -0.01161428337229559707660497,
                  0.01161428337229559707660498,
                               0.0,
                 0.0002798785734386036356090761,
                -0.0002798785734386036356090759,
                               0.0,
                 -0.01157006575352843784923587,
                  0.01157006575352843784923588,
                               0.0,
                 0.007595676796247792031257549,
                 -0.007595676796247792031257541,
                               0.0,
                 -0.008019157976636872888817360,
                 -0.006190971763611308709464433,
                  0.01421012974024818159828180,
                 0.008019157976636872888817360,
                 0.006190971763611308709464438,
                 -0.01421012974024818159828179,
                 -0.001557909336320037381214816,
                 -0.006585036001891291212813166,
                 0.008142945338211328594027987,
                 0.001557909336320037381214818,
                 0.006585036001891291212813170,
                 -0.008142945338211328594027984 };


static double NF_C_T_UL4_2D_W3[27] = {
                  0.09239439077785866069415689,
                  0.09239439077785866069415689,
                  0.2784709310156610905569492,
                  0.06110671294767265446150073,
                  0.06110671294767265446150073,
                  0.2297093574717731482432842,
                 -0.003494061049012103890782514,
                 -0.003494061049012103890782514,
                  0.01336302226341469970600685,
                  -.1137915807850695810786726,
                  -.1137915807850695810786726,
                  -.1206848185980322443408267,
                  -.4166001122518296766622062,
                  -.4166001122518296766622062,
                  -.1445649390652278462257972,
                  -.1067989118773909445492667,
                  0.2938090398456628637325732,
                  -.2737566366298440985085243,
                  -.1067989118773909445492667,
                  0.2938090398456628637325732,
                  -.2737566366298440985085243,
                  0.1976392025339522551122001,
                  0.1558497410481981330029510,
                 -0.01450456110399258628373663,
                  0.1976392025339522551122001,
                  0.1558497410481981330029510,
                 -0.01450456110399258628373663 };


static double NF_C_T_UL4_2D_W4[27] = {
                  0.09368214945994376916863499,
                  -.2784709310156610905569487,
                  -.2784709310156610905569486,
                  0.1074959315764278393202820,
                  -.2297093574717731482432833,
                  -.2297093574717731482432833,
                  0.02035114436143890748757189,
                 -0.01336302226341469970600685,
                 -0.01336302226341469970600686,
                  0.1068983429721069178165186,
                  0.1206848185980322443408265,
                  0.1206848185980322443408264,
                  0.6886352854384315070986156,
                  0.1445649390652278462257972,
                  0.1445649390652278462257975,
                  -.4607667645981160176918310,
                  -.1268513150932097097733154,
                  0.6743645883528979067903642,
                  0.6743645883528979067903639,
                  -.4607667645981160176918310,
                  -.1268513150932097097733154,
                  -.3679935046861429743988872,
                  0.05629402258974670839298525,
                 -0.02728490038176153582551210,
                 -0.02728490038176153582551199,
                  -.3679935046861429743988873,
                  0.05629402258974670839298523 };


static double NF_C_T_UL4_2D_W5[27] = {
                  0.09239439077785866069415689,
                  0.2784709310156610905569492,
                  0.09239439077785866069415689,
                  0.06110671294767265446150073,
                  0.2297093574717731482432842,
                  0.06110671294767265446150073,
                 -0.003494061049012103890782514,
                  0.01336302226341469970600685,
                 -0.003494061049012103890782514,
                  -.1137915807850695810786726,
                  -.1206848185980322443408267,
                  -.1137915807850695810786726,
                  -.4166001122518296766622062,
                  -.1445649390652278462257972,
                  -.4166001122518296766622062,
                  0.2938090398456628637325732,
                  -.2737566366298440985085243,
                  -.1067989118773909445492667,
                  -.2737566366298440985085243,
                  -.1067989118773909445492667,
                  0.2938090398456628637325732,
                  0.1558497410481981330029510,
                 -0.01450456110399258628373663,
                  0.1976392025339522551122001,
                 -0.01450456110399258628373663,
                  0.1976392025339522551122001,
                  0.1558497410481981330029510 };


static double NF_C_T_UL4_2D_W6[27] = {
                 -0.003827977430701715674726601,
                 -0.003827977430701715674726601,
                  0.01570874620567093141772023,
                 0.003776919736611775084891085,
                 0.003776919736611775084891085,
                 -0.003950847063964964077445233,
                -0.0002291354990264545112729660,
                -0.0002291354990264545112729660,
                 -0.001050097514298223711257017,
                  0.01434756519285661340167803,
                  0.01434756519285661340167803,
                 -0.01872463340534674520176998,
                 0.001709853655180828947381827,
                 0.001709853655180828947381827,
                  0.01860738754184006943097046,
                 -0.01651289350954159570792762,
                 -0.009321476403944588847586670,
                 0.005591982984107302720974849,
                 -0.01651289350954159570792762,
                 -0.009321476403944588847586670,
                 0.005591982984107302720974849,
                 0.002056165558056609179404659,
                 -0.007291490303649334283704869,
                 0.004405208138100025761779190,
                 0.002056165558056609179404659,
                 -0.007291490303649334283704869,
                 0.004405208138100025761779190 };


static double NF_C_T_UL4_2D_W7[27] = {
                 -0.004026395672133750034133483,
                  0.01551032796423889705831329,
                 -0.02356311930850639712658022,
                 -0.001801496204629293046168433,
                 -0.009529263005206032208504773,
                 0.005926270595947446116167904,
                 0.0007541842561755663669014762,
                -0.00006677775909620283308257939,
                 0.001575146271447335566885526,
                 -0.004985248490183240800793042,
                 -0.03805744708838659940424111,
                  0.02808695010802011780265495,
                 -0.01101354742610086366286695,
                 0.005883986460558376820721662,
                 -0.02791108131276010414645574,
                -0.0007160792793129837898087173,
                 0.007005963003141900918411776,
                 -0.03001237287855888907905228,
                  0.05025475980793777091359189,
                  0.02095846620869186562434822,
                  0.01323642392623698091612773,
                 -0.02256956667483210153240745,
                 -0.001525212371376798023813823,
                 -0.01557095339316957465167241,
                  0.01640107000066227399419338,
                  0.02339968328232480087492843,
                 0.002355328978869497366334809 };


static double NF_C_T_UL4_2D_W8[27] = {
                 -0.004026395672133750034133483,
                 -0.02356311930850639712658025,
                  0.01551032796423889705831329,
                 -0.001801496204629293046168434,
                 0.005926270595947446116167895,
                 -0.009529263005206032208504787,
                 0.0007541842561755663669014743,
                 0.001575146271447335566885525,
                -0.00006677775909620283308258004,
                 -0.004985248490183240800793042,
                  0.02808695010802011780265489,
                 -0.03805744708838659940424105,
                 -0.01101354742610086366286698,
                 -0.02791108131276010414645562,
                 0.005883986460558376820721597,
                  0.02095846620869186562434822,
                  0.01323642392623698091612772,
                  0.05025475980793777091359192,
                 -0.03001237287855888907905236,
                -0.0007160792793129837898086911,
                 0.007005963003141900918411787,
                  0.02339968328232480087492843,
                 0.002355328978869497366334806,
                  0.01640107000066227399419331,
                 -0.01557095339316957465167238,
                 -0.02256956667483210153240745,
                 -0.001525212371376798023813824 };


static double NF_C_T_UL4_2D_W9[27] = {
                 -0.003827977430701715674726601,
                  0.01570874620567093141772023,
                 -0.003827977430701715674726601,
                 0.003776919736611775084891085,
                 -0.003950847063964964077445233,
                 0.003776919736611775084891085,
                -0.0002291354990264545112729660,
                 -0.001050097514298223711257017,
                -0.0002291354990264545112729660,
                  0.01434756519285661340167803,
                 -0.01872463340534674520176998,
                  0.01434756519285661340167803,
                 0.001709853655180828947381827,
                  0.01860738754184006943097046,
                 0.001709853655180828947381827,
                 -0.009321476403944588847586670,
                 0.005591982984107302720974849,
                 -0.01651289350954159570792762,
                 0.005591982984107302720974849,
                 -0.01651289350954159570792762,
                 -0.009321476403944588847586670,
                 -0.007291490303649334283704869,
                 0.004405208138100025761779190,
                 0.002056165558056609179404659,
                 0.004405208138100025761779190,
                 0.002056165558056609179404659,
                 -0.007291490303649334283704869 };

static double NF_C_T_UL4_2D_Xi[39] = {
        0.0, 0.25, 0.5, 0.75,
        1.0, 0.75, 0.5, 0.25,
        0.0, 0.0, 0.0, 0.0, 
        0.0323649481112758931588480911328593,
        0.0323649481112758931588480911328593,
        0.9352701037774482136823038177342814,
        0.119350912282581309581102091581736,
        0.119350912282581309581102091581736,
        0.761298175434837380837795816836528,
        0.534611048270758309358680864963778,
        0.534611048270758309358680864963778,
       -0.069222096541516618717361729927556,
        0.203309900431282473351326588944569,
        0.203309900431282473351326588944569,
        0.593380199137435053297346822110862,
        0.398969302965855222611381867187058,
        0.398969302965855222611381867187058,
        0.202061394068289554777236265625884,
        0.593201213428212752488840882179699,
        0.0501781383104946650738269077613887,
        0.3566206482612925824373322100589123,
        0.593201213428212752488840882179699,
        0.0501781383104946650738269077613887,
        0.3566206482612925824373322100589123,
        0.807489003159792153166724890348745,
        0.0210220165361662971236385570923633,
        0.1714889803040415497096365525588917,
        0.807489003159792153166724890348745,
        0.0210220165361662971236385570923633,
        0.1714889803040415497096365525588917 };

static double NF_C_T_UL4_2D_Eta[39] = {
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.25, 0.5, 0.75,
        1.0, 0.75, 0.5, 0.25,
        0.0323649481112758931588480911328593,
        0.9352701037774482136823038177342814,
        0.0323649481112758931588480911328593,
        0.119350912282581309581102091581736,
        0.761298175434837380837795816836528,
        0.119350912282581309581102091581736,
        0.534611048270758309358680864963778,
       -0.069222096541516618717361729927556,
        0.534611048270758309358680864963778,
        0.203309900431282473351326588944569,
        0.593380199137435053297346822110862,
        0.203309900431282473351326588944569,
        0.398969302965855222611381867187058,
        0.202061394068289554777236265625884,
        0.398969302965855222611381867187058,
        0.0501781383104946650738269077613887,
        0.3566206482612925824373322100589123,
        0.593201213428212752488840882179699,
        0.3566206482612925824373322100589123,
        0.593201213428212752488840882179699,
        0.0501781383104946650738269077613887,
        0.0210220165361662971236385570923633,
        0.1714889803040415497096365525588917,
        0.807489003159792153166724890348745,
        0.1714889803040415497096365525588917,
        0.807489003159792153166724890348745,
        0.0210220165361662971236385570923633 };

static double NF_C_T_UL4_2D_T[5] = { -1, -0.5, 0.0, 0.5, 1 };

void NF_C_T_UL4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[ 0] = PointValues[ 0];
  Functionals[ 1] = PointValues[ 1];
  Functionals[ 2] = PointValues[ 2];
  Functionals[ 3] = PointValues[ 3];
  Functionals[ 4] = PointValues[ 4];
  Functionals[ 5] = PointValues[ 5];
  Functionals[ 6] = PointValues[ 6];
  Functionals[ 7] = PointValues[ 7];
  Functionals[ 8] = PointValues[ 8];
  Functionals[ 9] = PointValues[ 9];
  Functionals[10] = PointValues[10];
  Functionals[11] = PointValues[11];

  Functionals[12] =     PointValues[12]*NF_C_T_UL4_2D_W0[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W0[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W0[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W0[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W0[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W0[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W0[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W0[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W0[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W0[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W0[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W0[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W0[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W0[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W0[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W0[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W0[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W0[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W0[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W0[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W0[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W0[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W0[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W0[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W0[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W0[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W0[26];

  Functionals[13] =     PointValues[12]*NF_C_T_UL4_2D_W1[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W1[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W1[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W1[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W1[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W1[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W1[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W1[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W1[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W1[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W1[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W1[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W1[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W1[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W1[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W1[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W1[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W1[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W1[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W1[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W1[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W1[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W1[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W1[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W1[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W1[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W1[26];

  Functionals[14] =     PointValues[12]*NF_C_T_UL4_2D_W2[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W2[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W2[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W2[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W2[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W2[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W2[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W2[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W2[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W2[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W2[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W2[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W2[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W2[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W2[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W2[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W2[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W2[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W2[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W2[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W2[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W2[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W2[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W2[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W2[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W2[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W2[26];

  Functionals[15] =     PointValues[12]*NF_C_T_UL4_2D_W3[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W3[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W3[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W3[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W3[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W3[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W3[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W3[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W3[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W3[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W3[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W3[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W3[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W3[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W3[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W3[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W3[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W3[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W3[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W3[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W3[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W3[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W3[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W3[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W3[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W3[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W3[26];

  Functionals[16] =     PointValues[12]*NF_C_T_UL4_2D_W4[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W4[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W4[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W4[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W4[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W4[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W4[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W4[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W4[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W4[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W4[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W4[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W4[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W4[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W4[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W4[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W4[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W4[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W4[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W4[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W4[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W4[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W4[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W4[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W4[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W4[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W4[26];

  Functionals[17] =     PointValues[12]*NF_C_T_UL4_2D_W5[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W5[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W5[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W5[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W5[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W5[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W5[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W5[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W5[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W5[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W5[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W5[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W5[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W5[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W5[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W5[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W5[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W5[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W5[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W5[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W5[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W5[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W5[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W5[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W5[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W5[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W5[26];

  Functionals[18] =     PointValues[12]*NF_C_T_UL4_2D_W6[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W6[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W6[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W6[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W6[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W6[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W6[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W6[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W6[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W6[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W6[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W6[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W6[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W6[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W6[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W6[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W6[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W6[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W6[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W6[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W6[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W6[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W6[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W6[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W6[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W6[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W6[26];

  Functionals[19] =     PointValues[12]*NF_C_T_UL4_2D_W7[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W7[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W7[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W7[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W7[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W7[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W7[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W7[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W7[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W7[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W7[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W7[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W7[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W7[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W7[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W7[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W7[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W7[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W7[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W7[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W7[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W7[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W7[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W7[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W7[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W7[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W7[26];

  Functionals[20] =     PointValues[12]*NF_C_T_UL4_2D_W8[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W8[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W8[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W8[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W8[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W8[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W8[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W8[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W8[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W8[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W8[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W8[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W8[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W8[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W8[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W8[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W8[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W8[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W8[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W8[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W8[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W8[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W8[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W8[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W8[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W8[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W8[26];

  Functionals[21] =     PointValues[12]*NF_C_T_UL4_2D_W9[ 0]
                      + PointValues[13]*NF_C_T_UL4_2D_W9[ 1]
                      + PointValues[14]*NF_C_T_UL4_2D_W9[ 2]
                      + PointValues[15]*NF_C_T_UL4_2D_W9[ 3]
                      + PointValues[16]*NF_C_T_UL4_2D_W9[ 4]
                      + PointValues[17]*NF_C_T_UL4_2D_W9[ 5]
                      + PointValues[18]*NF_C_T_UL4_2D_W9[ 6]
                      + PointValues[19]*NF_C_T_UL4_2D_W9[ 7]
                      + PointValues[20]*NF_C_T_UL4_2D_W9[ 8]
                      + PointValues[21]*NF_C_T_UL4_2D_W9[ 9]
                      + PointValues[22]*NF_C_T_UL4_2D_W9[10]
                      + PointValues[23]*NF_C_T_UL4_2D_W9[11]
                      + PointValues[24]*NF_C_T_UL4_2D_W9[12]
                      + PointValues[25]*NF_C_T_UL4_2D_W9[13]
                      + PointValues[26]*NF_C_T_UL4_2D_W9[14]
                      + PointValues[27]*NF_C_T_UL4_2D_W9[15]
                      + PointValues[28]*NF_C_T_UL4_2D_W9[16]
                      + PointValues[29]*NF_C_T_UL4_2D_W9[17]
                      + PointValues[30]*NF_C_T_UL4_2D_W9[18]
                      + PointValues[31]*NF_C_T_UL4_2D_W9[19]
                      + PointValues[32]*NF_C_T_UL4_2D_W9[20]
                      + PointValues[33]*NF_C_T_UL4_2D_W9[21]
                      + PointValues[34]*NF_C_T_UL4_2D_W9[22]
                      + PointValues[35]*NF_C_T_UL4_2D_W9[23]
                      + PointValues[36]*NF_C_T_UL4_2D_W9[24]
                      + PointValues[37]*NF_C_T_UL4_2D_W9[25]
                      + PointValues[38]*NF_C_T_UL4_2D_W9[26];
}

void NF_C_T_UL4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[ 0] = PointValues[ 0];
  Functionals[ 1] = PointValues[ 1];
  Functionals[ 2] = PointValues[ 2];
  Functionals[ 3] = PointValues[ 3];
  Functionals[ 4] = PointValues[ 4];
}

TNodalFunctional2D *NF_C_T_UL4_2D_Obj = new TNodalFunctional2D
        (NF_C_T_UL4_2D, 22, 5, 39, 5, NF_C_T_UL4_2D_Xi, NF_C_T_UL4_2D_Eta,
         NF_C_T_UL4_2D_T, NF_C_T_UL4_2D_EvalAll, NF_C_T_UL4_2D_EvalEdge);
