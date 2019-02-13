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
   
static double NF_N_T_P4_2D_Xi[30] = 
        { 
          0.046910077030668003601186560850,
          0.230765344947158454481842789650,
          0.5,
          0.769234655052841545518157210350,
          0.953089922969331996398813439150,
          0.953089922969331996398813439150,
          0.769234655052841545518157210350,
          0.5,
          0.230765344947158454481842789650,
          0.046910077030668003601186560850,
          0, 0, 0, 0, 0,
          0.5133469206394541,
          0.3132512106717253,
          0.6517753036487957,
          0.06510199345893917,
          0.345792011168269,
          0.2810412473151104,
          0.6306214343189561,
          0.313477887523733,
          0.8701651015635631,
          0.3623168221569262E1,
          0.2056118320454355,
          0.5612735500931855E-1,
          0.3474680882747129E-1,
          0.6473290497749777E-1,
         -0.2968960232737531E1
        };

static double NF_N_T_P4_2D_Eta[30] = 
        { 
          0, 0, 0, 0, 0, 
          0.046910077030668003601186560850,
          0.230765344947158454481842789650,
          0.5,
          0.769234655052841545518157210350,
          0.953089922969331996398813439150,
          0.953089922969331996398813439150,
          0.769234655052841545518157210350,
          0.5,
          0.230765344947158454481842789650,
          0.046910077030668003601186560850,
          0.2810412473151104,
          0.6306214343189561,
          0.313477887523733,
          0.8701651015635631,
          0.3623168221569262E1,
          0.2056118320454355,
          0.5612735500931855E-1,
          0.3474680882747129E-1,
          0.6473290497749777E-1,
         -0.2968960232737531E1,
          0.5133469206394541,
          0.3132512106717253,
          0.6517753036487957,
          0.6510199345893917E-1,
          0.345792011168269
        };

static double NF_N_T_P4_2D_T[5] = { -0.906179845938663992797626878299,
                                    -0.538469310105683091036314420699,
                                     0,
                                     0.538469310105683091036314420699,
                                     0.906179845938663992797626878299 };

static double NF_N_T_P4_2D_EdgeWeight0[5] = {
        0.1184634425280945437571320205,
        0.2393143352496832340206457575,
        0.2844444444444444444444444444,
        0.2393143352496832340206457575,
        0.1184634425280945437571320205 };

static double NF_N_T_P4_2D_EdgeWeight1[5] = {
       -0.3220475522984174693134711033,
       -0.3865902750008912622604589675,
        0.0,
        0.3865902750008912622604589675,
        0.3220475522984174693134711033 };

static double NF_N_T_P4_2D_EdgeWeight2[5] = {
        0.4334238969965230841044414528,
       -0.0778683414409675285488858967,
       -0.7111111111111111111111111111,
       -0.0778683414409675285488858967,
        0.4334238969965230841044414528 };

static double NF_N_T_P4_2D_EdgeWeight3[5] = {
       -0.4154771413508325868510801940,
        0.6991986448892333762714923941,
        0.0,
       -0.6991986448892333762714923941,
        0.4154771413508325868510801940 };


static double NF_N_T_P4_2D_CellWeight0[15] = {
        0.03396003583624729958603721,
        0.001659928772460918247370114,
        0.0004242940743545037796319238,
        0.001334070902641637542676619,
        0.1698598594336262696527415e-7,
        0.2116868384521265935109703,
        0.05170419916926577994680451,
        0.1492907562798600740896240,
        0.001349327252954068718657978,
        0.2304156165773590450668980e-9,
        0.06344714185011418859601229,
        0.2095452185809266179080118,
        0.03453432245238165684485284,
        0.2410638238638489590105253,
        0.2529641609960604545839634e-7 };

static double NF_N_T_P4_2D_CellWeight1[15] = {
        0.2116868384521265935109703,
        0.05170419916926576344113832,
        0.1492907562798600695085742,
        0.001349327252954070376766655,
        0.2304156165773590450668980e-9,
        0.06344714185011418859601229,
        0.2095452185809265846796426,
        0.03453432245238165464155087,
        0.2410638238638489811731078,
        0.2529641609960604545839634e-7,
        0.03396003583624729958603721,
        0.001659928772460915289937375,
        0.0004242940743545035354114681,
        0.001334070902641639191384834,
        0.1698598594336262696527415e-7 };

static double NF_N_T_P4_2D_CellWeight2[15] = {
        0.06344714185011418859601229,
        0.2095452185809265846796426,
        0.03453432245238165464155087,
        0.2410638238638489811731078,
        0.2529641609960604545839634e-7,
        0.03396003583624729958603721,
        0.001659928772460915289937375,
        0.0004242940743545035354114681,
        0.001334070902641639191384834,
        0.1698598594336262696527415e-7,
        0.2116868384521265935109703,
        0.05170419916926576344113832,
        0.1492907562798600695085742,
        0.001349327252954070376766655,
        0.2304156165773590450668980e-9 };

static double NF_N_T_P4_2D_CellWeight3[15] = {
        0.1695746752306618684689641,
        0.01852838771810587384581863,
        0.01591768616922033054907036,
        0.002683354785567795807609552,
        -0.3956683674146445297237118e-8,
        0.2317837342616973425616102,
        0.2081765377411735185945207,
        0.1436057814508679070699259,
        0.03607070763047549758772657,
        0.4828535725367747122020117e-8,
        0.09283678604814674727640592,
        0.03730040951271670414894234,
        0.007655771255306214695467177,
        0.03586620878187657586340455,
        -0.4145766844796558689352039e-7 };

static double NF_N_T_P4_2D_CellWeight4[15] = {
        0.2317837342616973425616102,
        0.2081765377411734853661515,
        0.1436057814508679048666239,
        0.03607070763047551975030914,
        0.4828535725367747122020117e-8,
        0.09283678604814674727640592,
        0.03730040951271670119150960,
        0.007655771255306214451246721,
        0.03586620878187657751211276,
        -0.4145766844796558689352039e-7,
        0.1695746752306618684689641,
        0.01852838771810585734015244,
        0.01591768616922032596802056,
        0.002683354785567797465718229,
        -0.3956683674146445297237118e-8 };

static double NF_N_T_P4_2D_CellWeight5[15] = {
        0.09283678604814674727640592,
        0.03730040951271673441987882,
        0.007655771255306216654548697,
        0.03586620878187655534953020,
        -0.4145766844796558689352039e-7,
        0.1695746752306618684689641,
        0.01852838771810586029758518,
        0.01591768616922032621224102,
        0.002683354785567795817010014,
        -0.3956683674146445297237118e-8,
        0.2317837342616973425616102,
        0.2081765377411735018718177,
        0.1436057814508679094476737,
        0.03607070763047551809220046,
        0.4828535725367747122020117e-8 };

void NF_N_T_P4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] =( NF_N_T_P4_2D_EdgeWeight0[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight0[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight0[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight0[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight0[4]*PointValues[ 4]);
  Functionals[1] =( NF_N_T_P4_2D_EdgeWeight0[0]*PointValues[ 5]
                   +NF_N_T_P4_2D_EdgeWeight0[1]*PointValues[ 6]
                   +NF_N_T_P4_2D_EdgeWeight0[2]*PointValues[ 7]
                   +NF_N_T_P4_2D_EdgeWeight0[3]*PointValues[ 8]
                   +NF_N_T_P4_2D_EdgeWeight0[4]*PointValues[ 9]);
  Functionals[2] =( NF_N_T_P4_2D_EdgeWeight0[0]*PointValues[10]
                   +NF_N_T_P4_2D_EdgeWeight0[1]*PointValues[11]
                   +NF_N_T_P4_2D_EdgeWeight0[2]*PointValues[12]
                   +NF_N_T_P4_2D_EdgeWeight0[3]*PointValues[13]
                   +NF_N_T_P4_2D_EdgeWeight0[4]*PointValues[14]);

  Functionals[3] =( NF_N_T_P4_2D_EdgeWeight1[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight1[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight1[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight1[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight1[4]*PointValues[ 4]);
  Functionals[4] =( NF_N_T_P4_2D_EdgeWeight1[0]*PointValues[ 5]
                   +NF_N_T_P4_2D_EdgeWeight1[1]*PointValues[ 6]
                   +NF_N_T_P4_2D_EdgeWeight1[2]*PointValues[ 7]
                   +NF_N_T_P4_2D_EdgeWeight1[3]*PointValues[ 8]
                   +NF_N_T_P4_2D_EdgeWeight1[4]*PointValues[ 9]);
  Functionals[5] =( NF_N_T_P4_2D_EdgeWeight1[0]*PointValues[10]
                   +NF_N_T_P4_2D_EdgeWeight1[1]*PointValues[11]
                   +NF_N_T_P4_2D_EdgeWeight1[2]*PointValues[12]
                   +NF_N_T_P4_2D_EdgeWeight1[3]*PointValues[13]
                   +NF_N_T_P4_2D_EdgeWeight1[4]*PointValues[14]);

  Functionals[6] =( NF_N_T_P4_2D_EdgeWeight2[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight2[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight2[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight2[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight2[4]*PointValues[ 4]);
  Functionals[7] =( NF_N_T_P4_2D_EdgeWeight2[0]*PointValues[ 5]
                   +NF_N_T_P4_2D_EdgeWeight2[1]*PointValues[ 6]
                   +NF_N_T_P4_2D_EdgeWeight2[2]*PointValues[ 7]
                   +NF_N_T_P4_2D_EdgeWeight2[3]*PointValues[ 8]
                   +NF_N_T_P4_2D_EdgeWeight2[4]*PointValues[ 9]);
  Functionals[8] =( NF_N_T_P4_2D_EdgeWeight2[0]*PointValues[10]
                   +NF_N_T_P4_2D_EdgeWeight2[1]*PointValues[11]
                   +NF_N_T_P4_2D_EdgeWeight2[2]*PointValues[12]
                   +NF_N_T_P4_2D_EdgeWeight2[3]*PointValues[13]
                   +NF_N_T_P4_2D_EdgeWeight2[4]*PointValues[14]);

  Functionals[9] =( NF_N_T_P4_2D_EdgeWeight3[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight3[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight3[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight3[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight3[4]*PointValues[ 4]);
  Functionals[10]=( NF_N_T_P4_2D_EdgeWeight3[0]*PointValues[ 5]
                   +NF_N_T_P4_2D_EdgeWeight3[1]*PointValues[ 6]
                   +NF_N_T_P4_2D_EdgeWeight3[2]*PointValues[ 7]
                   +NF_N_T_P4_2D_EdgeWeight3[3]*PointValues[ 8]
                   +NF_N_T_P4_2D_EdgeWeight3[4]*PointValues[ 9]);
  Functionals[11]=( NF_N_T_P4_2D_EdgeWeight3[0]*PointValues[10]
                   +NF_N_T_P4_2D_EdgeWeight3[1]*PointValues[11]
                   +NF_N_T_P4_2D_EdgeWeight3[2]*PointValues[12]
                   +NF_N_T_P4_2D_EdgeWeight3[3]*PointValues[13]
                   +NF_N_T_P4_2D_EdgeWeight3[4]*PointValues[14]);

  Functionals[12]=( NF_N_T_P4_2D_CellWeight0[ 0]*PointValues[15]
                   +NF_N_T_P4_2D_CellWeight0[ 1]*PointValues[16]
                   +NF_N_T_P4_2D_CellWeight0[ 2]*PointValues[17]
                   +NF_N_T_P4_2D_CellWeight0[ 3]*PointValues[18]
                   +NF_N_T_P4_2D_CellWeight0[ 4]*PointValues[19]
                   +NF_N_T_P4_2D_CellWeight0[ 5]*PointValues[20]
                   +NF_N_T_P4_2D_CellWeight0[ 6]*PointValues[21]
                   +NF_N_T_P4_2D_CellWeight0[ 7]*PointValues[22]
                   +NF_N_T_P4_2D_CellWeight0[ 8]*PointValues[23]
                   +NF_N_T_P4_2D_CellWeight0[ 9]*PointValues[24]
                   +NF_N_T_P4_2D_CellWeight0[10]*PointValues[25]
                   +NF_N_T_P4_2D_CellWeight0[11]*PointValues[26]
                   +NF_N_T_P4_2D_CellWeight0[12]*PointValues[27]
                   +NF_N_T_P4_2D_CellWeight0[13]*PointValues[28]
                   +NF_N_T_P4_2D_CellWeight0[14]*PointValues[29] );
  Functionals[13]=( NF_N_T_P4_2D_CellWeight1[ 0]*PointValues[15]
                   +NF_N_T_P4_2D_CellWeight1[ 1]*PointValues[16]
                   +NF_N_T_P4_2D_CellWeight1[ 2]*PointValues[17]
                   +NF_N_T_P4_2D_CellWeight1[ 3]*PointValues[18]
                   +NF_N_T_P4_2D_CellWeight1[ 4]*PointValues[19]
                   +NF_N_T_P4_2D_CellWeight1[ 5]*PointValues[20]
                   +NF_N_T_P4_2D_CellWeight1[ 6]*PointValues[21]
                   +NF_N_T_P4_2D_CellWeight1[ 7]*PointValues[22]
                   +NF_N_T_P4_2D_CellWeight1[ 8]*PointValues[23]
                   +NF_N_T_P4_2D_CellWeight1[ 9]*PointValues[24]
                   +NF_N_T_P4_2D_CellWeight1[10]*PointValues[25]
                   +NF_N_T_P4_2D_CellWeight1[11]*PointValues[26]
                   +NF_N_T_P4_2D_CellWeight1[12]*PointValues[27]
                   +NF_N_T_P4_2D_CellWeight1[13]*PointValues[28]
                   +NF_N_T_P4_2D_CellWeight1[14]*PointValues[29] );
  Functionals[14]=( NF_N_T_P4_2D_CellWeight2[ 0]*PointValues[15]
                   +NF_N_T_P4_2D_CellWeight2[ 1]*PointValues[16]
                   +NF_N_T_P4_2D_CellWeight2[ 2]*PointValues[17]
                   +NF_N_T_P4_2D_CellWeight2[ 3]*PointValues[18]
                   +NF_N_T_P4_2D_CellWeight2[ 4]*PointValues[19]
                   +NF_N_T_P4_2D_CellWeight2[ 5]*PointValues[20]
                   +NF_N_T_P4_2D_CellWeight2[ 6]*PointValues[21]
                   +NF_N_T_P4_2D_CellWeight2[ 7]*PointValues[22]
                   +NF_N_T_P4_2D_CellWeight2[ 8]*PointValues[23]
                   +NF_N_T_P4_2D_CellWeight2[ 9]*PointValues[24]
                   +NF_N_T_P4_2D_CellWeight2[10]*PointValues[25]
                   +NF_N_T_P4_2D_CellWeight2[11]*PointValues[26]
                   +NF_N_T_P4_2D_CellWeight2[12]*PointValues[27]
                   +NF_N_T_P4_2D_CellWeight2[13]*PointValues[28]
                   +NF_N_T_P4_2D_CellWeight2[14]*PointValues[29] );
  Functionals[15]=( NF_N_T_P4_2D_CellWeight3[ 0]*PointValues[15]
                   +NF_N_T_P4_2D_CellWeight3[ 1]*PointValues[16]
                   +NF_N_T_P4_2D_CellWeight3[ 2]*PointValues[17]
                   +NF_N_T_P4_2D_CellWeight3[ 3]*PointValues[18]
                   +NF_N_T_P4_2D_CellWeight3[ 4]*PointValues[19]
                   +NF_N_T_P4_2D_CellWeight3[ 5]*PointValues[20]
                   +NF_N_T_P4_2D_CellWeight3[ 6]*PointValues[21]
                   +NF_N_T_P4_2D_CellWeight3[ 7]*PointValues[22]
                   +NF_N_T_P4_2D_CellWeight3[ 8]*PointValues[23]
                   +NF_N_T_P4_2D_CellWeight3[ 9]*PointValues[24]
                   +NF_N_T_P4_2D_CellWeight3[10]*PointValues[25]
                   +NF_N_T_P4_2D_CellWeight3[11]*PointValues[26]
                   +NF_N_T_P4_2D_CellWeight3[12]*PointValues[27]
                   +NF_N_T_P4_2D_CellWeight3[13]*PointValues[28]
                   +NF_N_T_P4_2D_CellWeight3[14]*PointValues[29] );
  Functionals[16]=( NF_N_T_P4_2D_CellWeight4[ 0]*PointValues[15]
                   +NF_N_T_P4_2D_CellWeight4[ 1]*PointValues[16]
                   +NF_N_T_P4_2D_CellWeight4[ 2]*PointValues[17]
                   +NF_N_T_P4_2D_CellWeight4[ 3]*PointValues[18]
                   +NF_N_T_P4_2D_CellWeight4[ 4]*PointValues[19]
                   +NF_N_T_P4_2D_CellWeight4[ 5]*PointValues[20]
                   +NF_N_T_P4_2D_CellWeight4[ 6]*PointValues[21]
                   +NF_N_T_P4_2D_CellWeight4[ 7]*PointValues[22]
                   +NF_N_T_P4_2D_CellWeight4[ 8]*PointValues[23]
                   +NF_N_T_P4_2D_CellWeight4[ 9]*PointValues[24]
                   +NF_N_T_P4_2D_CellWeight4[10]*PointValues[25]
                   +NF_N_T_P4_2D_CellWeight4[11]*PointValues[26]
                   +NF_N_T_P4_2D_CellWeight4[12]*PointValues[27]
                   +NF_N_T_P4_2D_CellWeight4[13]*PointValues[28]
                   +NF_N_T_P4_2D_CellWeight4[14]*PointValues[29] );
  Functionals[17]=( NF_N_T_P4_2D_CellWeight5[ 0]*PointValues[15]
                   +NF_N_T_P4_2D_CellWeight5[ 1]*PointValues[16]
                   +NF_N_T_P4_2D_CellWeight5[ 2]*PointValues[17]
                   +NF_N_T_P4_2D_CellWeight5[ 3]*PointValues[18]
                   +NF_N_T_P4_2D_CellWeight5[ 4]*PointValues[19]
                   +NF_N_T_P4_2D_CellWeight5[ 5]*PointValues[20]
                   +NF_N_T_P4_2D_CellWeight5[ 6]*PointValues[21]
                   +NF_N_T_P4_2D_CellWeight5[ 7]*PointValues[22]
                   +NF_N_T_P4_2D_CellWeight5[ 8]*PointValues[23]
                   +NF_N_T_P4_2D_CellWeight5[ 9]*PointValues[24]
                   +NF_N_T_P4_2D_CellWeight5[10]*PointValues[25]
                   +NF_N_T_P4_2D_CellWeight5[11]*PointValues[26]
                   +NF_N_T_P4_2D_CellWeight5[12]*PointValues[27]
                   +NF_N_T_P4_2D_CellWeight5[13]*PointValues[28]
                   +NF_N_T_P4_2D_CellWeight5[14]*PointValues[29] );

  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
    {
      Functionals[3] = -Functionals[3];
      Functionals[9] = -Functionals[9];
    }
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
    {
      Functionals[4] = -Functionals[4];
      Functionals[10] = -Functionals[10];
    }
    if(Cell->GetVertex(2) > Cell->GetVertex(0))
    {
      Functionals[5] = -Functionals[5];
      Functionals[11] = -Functionals[11];
    }
  }
  */

  if(Cell)
  {
    OwnNum = Coll->GetIndex(Cell);

    neigh = Cell->GetJoint(0)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[3] = -Functionals[3];
        Functionals[9] = -Functionals[9];
      }
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[ 4] = -Functionals[ 4];
        Functionals[10] = -Functionals[10];
      }
    } // endif neigh

    neigh = Cell->GetJoint(2)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[ 5] = -Functionals[ 5];
        Functionals[11] = -Functionals[11];
      }
    } // endif neigh
  } // endif Cell
}

void NF_N_T_P4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                           double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] =( NF_N_T_P4_2D_EdgeWeight0[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight0[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight0[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight0[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight0[4]*PointValues[ 4]);

  Functionals[1] =( NF_N_T_P4_2D_EdgeWeight1[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight1[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight1[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight1[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight1[4]*PointValues[ 4]);

  Functionals[2] =( NF_N_T_P4_2D_EdgeWeight2[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight2[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight2[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight2[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight2[4]*PointValues[ 4]);

  Functionals[3] =( NF_N_T_P4_2D_EdgeWeight3[0]*PointValues[ 0]
                   +NF_N_T_P4_2D_EdgeWeight3[1]*PointValues[ 1]
                   +NF_N_T_P4_2D_EdgeWeight3[2]*PointValues[ 2]
                   +NF_N_T_P4_2D_EdgeWeight3[3]*PointValues[ 3]
                   +NF_N_T_P4_2D_EdgeWeight3[4]*PointValues[ 4]);

  if(Joint != -1)
  {
    // if(Cell->GetVertex(Joint) > Cell->GetVertex((Joint+1)%3))
    // {
    //   Functionals[1] = -Functionals[1];
    //   Functionals[3] = -Functionals[3];
    // }
    neigh = Cell->GetJoint(Joint)->GetNeighbour(Cell);
    if(neigh)
    {
      OwnNum = Coll->GetIndex(Cell);
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        Functionals[1] = -Functionals[1];
        Functionals[3] = -Functionals[3];
      }
    } // endif neigh
    // */
  }
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_N_T_P4_2D_Obj = new TNodalFunctional2D
        (NF_N_T_P4_2D, 18, 4, 30, 5, NF_N_T_P4_2D_Xi, NF_N_T_P4_2D_Eta,
         NF_N_T_P4_2D_T, NF_N_T_P4_2D_EvalAll, NF_N_T_P4_2D_EvalEdge);
