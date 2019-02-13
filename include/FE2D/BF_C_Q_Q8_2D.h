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
   

// ***********************************************************************
// Q8 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_Q8_2D_Funct(double xi, double eta, double *values)
{ 

      double xi0= 0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi*xi-0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi*xi+0.1422222222222222E1*xi*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi*xi-0.3111111111111111*xi*xi*xi-0.1428571428571429E-1*xi*xi+0.1428571428571429E-1*xi;
      double xi1= -0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi*xi+0.9752380952380952E1*xi*xi*xi*xi*xi*xi*xi+0.1706666666666667E2*xi*xi*xi*xi*xi*xi-0.128E2*xi*xi*xi*xi*xi-0.4266666666666667E1*xi*xi*xi*xi+0.32E1*xi*xi*xi+0.2031746031746032*xi*xi-0.1523809523809524*xi;
      double xi2= 0.4551111111111111E2*xi*xi*xi*xi*xi*xi*xi*xi-0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi-0.7395555555555556E2*xi*xi*xi*xi*xi*xi+0.3697777777777778E2*xi*xi*xi*xi*xi+0.3004444444444444E2*xi*xi*xi*xi-0.1502222222222222E2*xi*xi*xi-0.16E1*xi*xi+0.8*xi;
      double xi3= -0.9102222222222222E2*xi*xi*xi*xi*xi*xi*xi*xi+0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi+0.1649777777777778E3*xi*xi*xi*xi*xi*xi-0.4124444444444444E2*xi*xi*xi*xi*xi-0.8675555555555556E2*xi*xi*xi*xi+0.2168888888888889E2*xi*xi*xi+0.128E2*xi*xi-0.32E1*xi;
      double xi4= 0.1137777777777778E3*xi*xi*xi*xi*xi*xi*xi*xi-0.2133333333333333E3*xi*xi*xi*xi*xi*xi+0.1213333333333333E3*xi*xi*xi*xi-0.2277777777777778E2*xi*xi+1.0;
      double xi5= -0.9102222222222222E2*xi*xi*xi*xi*xi*xi*xi*xi-0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi+0.1649777777777778E3*xi*xi*xi*xi*xi*xi+0.4124444444444444E2*xi*xi*xi*xi*xi-0.8675555555555556E2*xi*xi*xi*xi-0.2168888888888889E2*xi*xi*xi+0.128E2*xi*xi+0.32E1*xi;
      double xi6= 0.4551111111111111E2*xi*xi*xi*xi*xi*xi*xi*xi+0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi-0.7395555555555556E2*xi*xi*xi*xi*xi*xi-0.3697777777777778E2*xi*xi*xi*xi*xi+0.3004444444444444E2*xi*xi*xi*xi+0.1502222222222222E2*xi*xi*xi-0.16E1*xi*xi-0.8*xi;
      double xi7= -0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi*xi-0.9752380952380952E1*xi*xi*xi*xi*xi*xi*xi+0.1706666666666667E2*xi*xi*xi*xi*xi*xi+0.128E2*xi*xi*xi*xi*xi-0.4266666666666667E1*xi*xi*xi*xi-0.32E1*xi*xi*xi+0.2031746031746032*xi*xi+0.1523809523809524*xi;
      double xi8= 0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi*xi+0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi-0.1428571428571429E-1*xi*xi-0.1428571428571429E-1*xi;

      double eta0= 0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta*eta-0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta*eta+0.1422222222222222E1*eta*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta*eta-0.3111111111111111*eta*eta*eta-0.1428571428571429E-1*eta*eta+0.1428571428571429E-1*eta;
      double eta1= -0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta*eta+0.9752380952380952E1*eta*eta*eta*eta*eta*eta*eta+0.1706666666666667E2*eta*eta*eta*eta*eta*eta-0.128E2*eta*eta*eta*eta*eta-0.4266666666666667E1*eta*eta*eta*eta+0.32E1*eta*eta*eta+0.2031746031746032*eta*eta-0.1523809523809524*eta;
      double eta2= 0.4551111111111111E2*eta*eta*eta*eta*eta*eta*eta*eta-0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta-0.7395555555555556E2*eta*eta*eta*eta*eta*eta+0.3697777777777778E2*eta*eta*eta*eta*eta+0.3004444444444444E2*eta*eta*eta*eta-0.1502222222222222E2*eta*eta*eta-0.16E1*eta*eta+0.8*eta;
      double eta3= -0.9102222222222222E2*eta*eta*eta*eta*eta*eta*eta*eta+0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta+0.1649777777777778E3*eta*eta*eta*eta*eta*eta-0.4124444444444444E2*eta*eta*eta*eta*eta-0.8675555555555556E2*eta*eta*eta*eta+0.2168888888888889E2*eta*eta*eta+0.128E2*eta*eta-0.32E1*eta;
      double eta4= 0.1137777777777778E3*eta*eta*eta*eta*eta*eta*eta*eta-0.2133333333333333E3*eta*eta*eta*eta*eta*eta+0.1213333333333333E3*eta*eta*eta*eta-0.2277777777777778E2*eta*eta+1.0;
      double eta5= -0.9102222222222222E2*eta*eta*eta*eta*eta*eta*eta*eta-0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta+0.1649777777777778E3*eta*eta*eta*eta*eta*eta+0.4124444444444444E2*eta*eta*eta*eta*eta-0.8675555555555556E2*eta*eta*eta*eta-0.2168888888888889E2*eta*eta*eta+0.128E2*eta*eta+0.32E1*eta;
      double eta6= 0.4551111111111111E2*eta*eta*eta*eta*eta*eta*eta*eta+0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta-0.7395555555555556E2*eta*eta*eta*eta*eta*eta-0.3697777777777778E2*eta*eta*eta*eta*eta+0.3004444444444444E2*eta*eta*eta*eta+0.1502222222222222E2*eta*eta*eta-0.16E1*eta*eta-0.8*eta;
      double eta7= -0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta*eta-0.9752380952380952E1*eta*eta*eta*eta*eta*eta*eta+0.1706666666666667E2*eta*eta*eta*eta*eta*eta+0.128E2*eta*eta*eta*eta*eta-0.4266666666666667E1*eta*eta*eta*eta-0.32E1*eta*eta*eta+0.2031746031746032*eta*eta+0.1523809523809524*eta;
      double eta8= 0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta*eta+0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta-0.1428571428571429E-1*eta*eta-0.1428571428571429E-1*eta;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi7*eta0;
      values[8] = xi8*eta0;
      values[9] = xi0*eta1;
      values[10] = xi1*eta1;
      values[11] = xi2*eta1;
      values[12] = xi3*eta1;
      values[13] = xi4*eta1;
      values[14] = xi5*eta1;
      values[15] = xi6*eta1;
      values[16] = xi7*eta1;
      values[17] = xi8*eta1;
      values[18] = xi0*eta2;
      values[19] = xi1*eta2;
      values[20] = xi2*eta2;
      values[21] = xi3*eta2;
      values[22] = xi4*eta2;
      values[23] = xi5*eta2;
      values[24] = xi6*eta2;
      values[25] = xi7*eta2;
      values[26] = xi8*eta2;
      values[27] = xi0*eta3;
      values[28] = xi1*eta3;
      values[29] = xi2*eta3;
      values[30] = xi3*eta3;
      values[31] = xi4*eta3;
      values[32] = xi5*eta3;
      values[33] = xi6*eta3;
      values[34] = xi7*eta3;
      values[35] = xi8*eta3;
      values[36] = xi0*eta4;
      values[37] = xi1*eta4;
      values[38] = xi2*eta4;
      values[39] = xi3*eta4;
      values[40] = xi4*eta4;
      values[41] = xi5*eta4;
      values[42] = xi6*eta4;
      values[43] = xi7*eta4;
      values[44] = xi8*eta4;
      values[45] = xi0*eta5;
      values[46] = xi1*eta5;
      values[47] = xi2*eta5;
      values[48] = xi3*eta5;
      values[49] = xi4*eta5;
      values[50] = xi5*eta5;
      values[51] = xi6*eta5;
      values[52] = xi7*eta5;
      values[53] = xi8*eta5;
      values[54] = xi0*eta6;
      values[55] = xi1*eta6;
      values[56] = xi2*eta6;
      values[57] = xi3*eta6;
      values[58] = xi4*eta6;
      values[59] = xi5*eta6;
      values[60] = xi6*eta6;
      values[61] = xi7*eta6;
      values[62] = xi8*eta6;
      values[63] = xi0*eta7;
      values[64] = xi1*eta7;
      values[65] = xi2*eta7;
      values[66] = xi3*eta7;
      values[67] = xi4*eta7;
      values[68] = xi5*eta7;
      values[69] = xi6*eta7;
      values[70] = xi7*eta7;
      values[71] = xi8*eta7;
      values[72] = xi0*eta8;
      values[73] = xi1*eta8;
      values[74] = xi2*eta8;
      values[75] = xi3*eta8;
      values[76] = xi4*eta8;
      values[77] = xi5*eta8;
      values[78] = xi6*eta8;
      values[79] = xi7*eta8;
      values[80] = xi8*eta8;
}


// values of the derivatives in xi direction
static void C_Q_Q8_2D_DeriveXi(double xi, double eta, double *values)
{

      double xi0= 0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi-0.1137777777777778E2*xi*xi*xi*xi*xi*xi-0.8533333333333333E1*xi*xi*xi*xi*xi+0.7111111111111111E1*xi*xi*xi*xi+0.1244444444444444E1*xi*xi*xi-0.9333333333333333*xi*xi-0.2857142857142857E-1*xi+0.1428571428571429E-1;
      double xi1= -0.1040253968253968E3*xi*xi*xi*xi*xi*xi*xi+0.6826666666666667E2*xi*xi*xi*xi*xi*xi+0.1024E3*xi*xi*xi*xi*xi-0.64E2*xi*xi*xi*xi-0.1706666666666667E2*xi*xi*xi+0.96E1*xi*xi+0.4063492063492063*xi-0.1523809523809524;
      double xi2= 0.3640888888888889E3*xi*xi*xi*xi*xi*xi*xi-0.1592888888888889E3*xi*xi*xi*xi*xi*xi-0.4437333333333333E3*xi*xi*xi*xi*xi+0.1848888888888889E3*xi*xi*xi*xi+0.1201777777777778E3*xi*xi*xi-0.4506666666666667E2*xi*xi-0.32E1*xi+0.8;
      double xi3= -0.7281777777777778E3*xi*xi*xi*xi*xi*xi*xi+0.1592888888888889E3*xi*xi*xi*xi*xi*xi+0.9898666666666667E3*xi*xi*xi*xi*xi-0.2062222222222222E3*xi*xi*xi*xi-0.3470222222222222E3*xi*xi*xi+0.6506666666666667E2*xi*xi+0.256E2*xi-0.32E1;
      double xi4= 0.9102222222222222E3*xi*xi*xi*xi*xi*xi*xi-0.128E4*xi*xi*xi*xi*xi+0.4853333333333333E3*xi*xi*xi-0.4555555555555556E2*xi;
      double xi5= -0.7281777777777778E3*xi*xi*xi*xi*xi*xi*xi-0.1592888888888889E3*xi*xi*xi*xi*xi*xi+0.9898666666666667E3*xi*xi*xi*xi*xi+0.2062222222222222E3*xi*xi*xi*xi-0.3470222222222222E3*xi*xi*xi-0.6506666666666667E2*xi*xi+0.256E2*xi+0.32E1;
      double xi6= 0.3640888888888889E3*xi*xi*xi*xi*xi*xi*xi+0.1592888888888889E3*xi*xi*xi*xi*xi*xi-0.4437333333333333E3*xi*xi*xi*xi*xi-0.1848888888888889E3*xi*xi*xi*xi+0.1201777777777778E3*xi*xi*xi+0.4506666666666667E2*xi*xi-0.32E1*xi-0.8;
      double xi7= -0.1040253968253968E3*xi*xi*xi*xi*xi*xi*xi-0.6826666666666667E2*xi*xi*xi*xi*xi*xi+0.1024E3*xi*xi*xi*xi*xi+0.64E2*xi*xi*xi*xi-0.1706666666666667E2*xi*xi*xi-0.96E1*xi*xi+0.4063492063492063*xi+0.1523809523809524;
      double xi8= 0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi+0.1137777777777778E2*xi*xi*xi*xi*xi*xi-0.8533333333333333E1*xi*xi*xi*xi*xi-0.7111111111111111E1*xi*xi*xi*xi+0.1244444444444444E1*xi*xi*xi+0.9333333333333333*xi*xi-0.2857142857142857E-1*xi-0.1428571428571429E-1;

      double eta0= 0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta*eta-0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta*eta+0.1422222222222222E1*eta*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta*eta-0.3111111111111111*eta*eta*eta-0.1428571428571429E-1*eta*eta+0.1428571428571429E-1*eta;
      double eta1= -0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta*eta+0.9752380952380952E1*eta*eta*eta*eta*eta*eta*eta+0.1706666666666667E2*eta*eta*eta*eta*eta*eta-0.128E2*eta*eta*eta*eta*eta-0.4266666666666667E1*eta*eta*eta*eta+0.32E1*eta*eta*eta+0.2031746031746032*eta*eta-0.1523809523809524*eta;
      double eta2= 0.4551111111111111E2*eta*eta*eta*eta*eta*eta*eta*eta-0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta-0.7395555555555556E2*eta*eta*eta*eta*eta*eta+0.3697777777777778E2*eta*eta*eta*eta*eta+0.3004444444444444E2*eta*eta*eta*eta-0.1502222222222222E2*eta*eta*eta-0.16E1*eta*eta+0.8*eta;
      double eta3= -0.9102222222222222E2*eta*eta*eta*eta*eta*eta*eta*eta+0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta+0.1649777777777778E3*eta*eta*eta*eta*eta*eta-0.4124444444444444E2*eta*eta*eta*eta*eta-0.8675555555555556E2*eta*eta*eta*eta+0.2168888888888889E2*eta*eta*eta+0.128E2*eta*eta-0.32E1*eta;
      double eta4= 0.1137777777777778E3*eta*eta*eta*eta*eta*eta*eta*eta-0.2133333333333333E3*eta*eta*eta*eta*eta*eta+0.1213333333333333E3*eta*eta*eta*eta-0.2277777777777778E2*eta*eta+1.0;
      double eta5= -0.9102222222222222E2*eta*eta*eta*eta*eta*eta*eta*eta-0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta+0.1649777777777778E3*eta*eta*eta*eta*eta*eta+0.4124444444444444E2*eta*eta*eta*eta*eta-0.8675555555555556E2*eta*eta*eta*eta-0.2168888888888889E2*eta*eta*eta+0.128E2*eta*eta+0.32E1*eta;
      double eta6= 0.4551111111111111E2*eta*eta*eta*eta*eta*eta*eta*eta+0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta-0.7395555555555556E2*eta*eta*eta*eta*eta*eta-0.3697777777777778E2*eta*eta*eta*eta*eta+0.3004444444444444E2*eta*eta*eta*eta+0.1502222222222222E2*eta*eta*eta-0.16E1*eta*eta-0.8*eta;
      double eta7= -0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta*eta-0.9752380952380952E1*eta*eta*eta*eta*eta*eta*eta+0.1706666666666667E2*eta*eta*eta*eta*eta*eta+0.128E2*eta*eta*eta*eta*eta-0.4266666666666667E1*eta*eta*eta*eta-0.32E1*eta*eta*eta+0.2031746031746032*eta*eta+0.1523809523809524*eta;
      double eta8= 0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta*eta+0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta-0.1428571428571429E-1*eta*eta-0.1428571428571429E-1*eta;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi7*eta0;
      values[8] = xi8*eta0;
      values[9] = xi0*eta1;
      values[10] = xi1*eta1;
      values[11] = xi2*eta1;
      values[12] = xi3*eta1;
      values[13] = xi4*eta1;
      values[14] = xi5*eta1;
      values[15] = xi6*eta1;
      values[16] = xi7*eta1;
      values[17] = xi8*eta1;
      values[18] = xi0*eta2;
      values[19] = xi1*eta2;
      values[20] = xi2*eta2;
      values[21] = xi3*eta2;
      values[22] = xi4*eta2;
      values[23] = xi5*eta2;
      values[24] = xi6*eta2;
      values[25] = xi7*eta2;
      values[26] = xi8*eta2;
      values[27] = xi0*eta3;
      values[28] = xi1*eta3;
      values[29] = xi2*eta3;
      values[30] = xi3*eta3;
      values[31] = xi4*eta3;
      values[32] = xi5*eta3;
      values[33] = xi6*eta3;
      values[34] = xi7*eta3;
      values[35] = xi8*eta3;
      values[36] = xi0*eta4;
      values[37] = xi1*eta4;
      values[38] = xi2*eta4;
      values[39] = xi3*eta4;
      values[40] = xi4*eta4;
      values[41] = xi5*eta4;
      values[42] = xi6*eta4;
      values[43] = xi7*eta4;
      values[44] = xi8*eta4;
      values[45] = xi0*eta5;
      values[46] = xi1*eta5;
      values[47] = xi2*eta5;
      values[48] = xi3*eta5;
      values[49] = xi4*eta5;
      values[50] = xi5*eta5;
      values[51] = xi6*eta5;
      values[52] = xi7*eta5;
      values[53] = xi8*eta5;
      values[54] = xi0*eta6;
      values[55] = xi1*eta6;
      values[56] = xi2*eta6;
      values[57] = xi3*eta6;
      values[58] = xi4*eta6;
      values[59] = xi5*eta6;
      values[60] = xi6*eta6;
      values[61] = xi7*eta6;
      values[62] = xi8*eta6;
      values[63] = xi0*eta7;
      values[64] = xi1*eta7;
      values[65] = xi2*eta7;
      values[66] = xi3*eta7;
      values[67] = xi4*eta7;
      values[68] = xi5*eta7;
      values[69] = xi6*eta7;
      values[70] = xi7*eta7;
      values[71] = xi8*eta7;
      values[72] = xi0*eta8;
      values[73] = xi1*eta8;
      values[74] = xi2*eta8;
      values[75] = xi3*eta8;
      values[76] = xi4*eta8;
      values[77] = xi5*eta8;
      values[78] = xi6*eta8;
      values[79] = xi7*eta8;
      values[80] = xi8*eta8;
}

// values of the derivatives in eta direction
static void C_Q_Q8_2D_DeriveEta(double xi, double eta, double *values)
{

      double xi0= 0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi*xi-0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi*xi+0.1422222222222222E1*xi*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi*xi-0.3111111111111111*xi*xi*xi-0.1428571428571429E-1*xi*xi+0.1428571428571429E-1*xi;
      double xi1= -0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi*xi+0.9752380952380952E1*xi*xi*xi*xi*xi*xi*xi+0.1706666666666667E2*xi*xi*xi*xi*xi*xi-0.128E2*xi*xi*xi*xi*xi-0.4266666666666667E1*xi*xi*xi*xi+0.32E1*xi*xi*xi+0.2031746031746032*xi*xi-0.1523809523809524*xi;
      double xi2= 0.4551111111111111E2*xi*xi*xi*xi*xi*xi*xi*xi-0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi-0.7395555555555556E2*xi*xi*xi*xi*xi*xi+0.3697777777777778E2*xi*xi*xi*xi*xi+0.3004444444444444E2*xi*xi*xi*xi-0.1502222222222222E2*xi*xi*xi-0.16E1*xi*xi+0.8*xi;
      double xi3= -0.9102222222222222E2*xi*xi*xi*xi*xi*xi*xi*xi+0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi+0.1649777777777778E3*xi*xi*xi*xi*xi*xi-0.4124444444444444E2*xi*xi*xi*xi*xi-0.8675555555555556E2*xi*xi*xi*xi+0.2168888888888889E2*xi*xi*xi+0.128E2*xi*xi-0.32E1*xi;
      double xi4= 0.1137777777777778E3*xi*xi*xi*xi*xi*xi*xi*xi-0.2133333333333333E3*xi*xi*xi*xi*xi*xi+0.1213333333333333E3*xi*xi*xi*xi-0.2277777777777778E2*xi*xi+1.0;
      double xi5= -0.9102222222222222E2*xi*xi*xi*xi*xi*xi*xi*xi-0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi+0.1649777777777778E3*xi*xi*xi*xi*xi*xi+0.4124444444444444E2*xi*xi*xi*xi*xi-0.8675555555555556E2*xi*xi*xi*xi-0.2168888888888889E2*xi*xi*xi+0.128E2*xi*xi+0.32E1*xi;
      double xi6= 0.4551111111111111E2*xi*xi*xi*xi*xi*xi*xi*xi+0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi-0.7395555555555556E2*xi*xi*xi*xi*xi*xi-0.3697777777777778E2*xi*xi*xi*xi*xi+0.3004444444444444E2*xi*xi*xi*xi+0.1502222222222222E2*xi*xi*xi-0.16E1*xi*xi-0.8*xi;
      double xi7= -0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi*xi-0.9752380952380952E1*xi*xi*xi*xi*xi*xi*xi+0.1706666666666667E2*xi*xi*xi*xi*xi*xi+0.128E2*xi*xi*xi*xi*xi-0.4266666666666667E1*xi*xi*xi*xi-0.32E1*xi*xi*xi+0.2031746031746032*xi*xi+0.1523809523809524*xi;
      double xi8= 0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi*xi+0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi-0.1428571428571429E-1*xi*xi-0.1428571428571429E-1*xi;

      double eta0= 0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta-0.1137777777777778E2*eta*eta*eta*eta*eta*eta-0.8533333333333333E1*eta*eta*eta*eta*eta+0.7111111111111111E1*eta*eta*eta*eta+0.1244444444444444E1*eta*eta*eta-0.9333333333333333*eta*eta-0.2857142857142857E-1*eta+0.1428571428571429E-1;
      double eta1= -0.1040253968253968E3*eta*eta*eta*eta*eta*eta*eta+0.6826666666666667E2*eta*eta*eta*eta*eta*eta+0.1024E3*eta*eta*eta*eta*eta-0.64E2*eta*eta*eta*eta-0.1706666666666667E2*eta*eta*eta+0.96E1*eta*eta+0.4063492063492063*eta-0.1523809523809524;
      double eta2= 0.3640888888888889E3*eta*eta*eta*eta*eta*eta*eta-0.1592888888888889E3*eta*eta*eta*eta*eta*eta-0.4437333333333333E3*eta*eta*eta*eta*eta+0.1848888888888889E3*eta*eta*eta*eta+0.1201777777777778E3*eta*eta*eta-0.4506666666666667E2*eta*eta-0.32E1*eta+0.8;
      double eta3= -0.7281777777777778E3*eta*eta*eta*eta*eta*eta*eta+0.1592888888888889E3*eta*eta*eta*eta*eta*eta+0.9898666666666667E3*eta*eta*eta*eta*eta-0.2062222222222222E3*eta*eta*eta*eta-0.3470222222222222E3*eta*eta*eta+0.6506666666666667E2*eta*eta+0.256E2*eta-0.32E1;
      double eta4= 0.9102222222222222E3*eta*eta*eta*eta*eta*eta*eta-0.128E4*eta*eta*eta*eta*eta+0.4853333333333333E3*eta*eta*eta-0.4555555555555556E2*eta;
      double eta5= -0.7281777777777778E3*eta*eta*eta*eta*eta*eta*eta-0.1592888888888889E3*eta*eta*eta*eta*eta*eta+0.9898666666666667E3*eta*eta*eta*eta*eta+0.2062222222222222E3*eta*eta*eta*eta-0.3470222222222222E3*eta*eta*eta-0.6506666666666667E2*eta*eta+0.256E2*eta+0.32E1;
      double eta6= 0.3640888888888889E3*eta*eta*eta*eta*eta*eta*eta+0.1592888888888889E3*eta*eta*eta*eta*eta*eta-0.4437333333333333E3*eta*eta*eta*eta*eta-0.1848888888888889E3*eta*eta*eta*eta+0.1201777777777778E3*eta*eta*eta+0.4506666666666667E2*eta*eta-0.32E1*eta-0.8;
      double eta7= -0.1040253968253968E3*eta*eta*eta*eta*eta*eta*eta-0.6826666666666667E2*eta*eta*eta*eta*eta*eta+0.1024E3*eta*eta*eta*eta*eta+0.64E2*eta*eta*eta*eta-0.1706666666666667E2*eta*eta*eta-0.96E1*eta*eta+0.4063492063492063*eta+0.1523809523809524;
      double eta8= 0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta+0.1137777777777778E2*eta*eta*eta*eta*eta*eta-0.8533333333333333E1*eta*eta*eta*eta*eta-0.7111111111111111E1*eta*eta*eta*eta+0.1244444444444444E1*eta*eta*eta+0.9333333333333333*eta*eta-0.2857142857142857E-1*eta-0.1428571428571429E-1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi7*eta0;
      values[8] = xi8*eta0;
      values[9] = xi0*eta1;
      values[10] = xi1*eta1;
      values[11] = xi2*eta1;
      values[12] = xi3*eta1;
      values[13] = xi4*eta1;
      values[14] = xi5*eta1;
      values[15] = xi6*eta1;
      values[16] = xi7*eta1;
      values[17] = xi8*eta1;
      values[18] = xi0*eta2;
      values[19] = xi1*eta2;
      values[20] = xi2*eta2;
      values[21] = xi3*eta2;
      values[22] = xi4*eta2;
      values[23] = xi5*eta2;
      values[24] = xi6*eta2;
      values[25] = xi7*eta2;
      values[26] = xi8*eta2;
      values[27] = xi0*eta3;
      values[28] = xi1*eta3;
      values[29] = xi2*eta3;
      values[30] = xi3*eta3;
      values[31] = xi4*eta3;
      values[32] = xi5*eta3;
      values[33] = xi6*eta3;
      values[34] = xi7*eta3;
      values[35] = xi8*eta3;
      values[36] = xi0*eta4;
      values[37] = xi1*eta4;
      values[38] = xi2*eta4;
      values[39] = xi3*eta4;
      values[40] = xi4*eta4;
      values[41] = xi5*eta4;
      values[42] = xi6*eta4;
      values[43] = xi7*eta4;
      values[44] = xi8*eta4;
      values[45] = xi0*eta5;
      values[46] = xi1*eta5;
      values[47] = xi2*eta5;
      values[48] = xi3*eta5;
      values[49] = xi4*eta5;
      values[50] = xi5*eta5;
      values[51] = xi6*eta5;
      values[52] = xi7*eta5;
      values[53] = xi8*eta5;
      values[54] = xi0*eta6;
      values[55] = xi1*eta6;
      values[56] = xi2*eta6;
      values[57] = xi3*eta6;
      values[58] = xi4*eta6;
      values[59] = xi5*eta6;
      values[60] = xi6*eta6;
      values[61] = xi7*eta6;
      values[62] = xi8*eta6;
      values[63] = xi0*eta7;
      values[64] = xi1*eta7;
      values[65] = xi2*eta7;
      values[66] = xi3*eta7;
      values[67] = xi4*eta7;
      values[68] = xi5*eta7;
      values[69] = xi6*eta7;
      values[70] = xi7*eta7;
      values[71] = xi8*eta7;
      values[72] = xi0*eta8;
      values[73] = xi1*eta8;
      values[74] = xi2*eta8;
      values[75] = xi3*eta8;
      values[76] = xi4*eta8;
      values[77] = xi5*eta8;
      values[78] = xi6*eta8;
      values[79] = xi7*eta8;
      values[80] = xi8*eta8;
}

// values of the derivatives in xi-xi  direction
static void C_Q_Q8_2D_DeriveXiXi(double xi, double eta, double *values)
{

      double xi0= 0.9102222222222222E2*xi*xi*xi*xi*xi*xi-0.6826666666666667E2*xi*xi*xi*xi*xi-0.4266666666666667E2*xi*xi*xi*xi+0.2844444444444444E2*xi*xi*xi+0.3733333333333333E1*xi*xi-0.1866666666666667E1*xi-0.2857142857142857E-1;
      double xi1= -0.7281777777777778E3*xi*xi*xi*xi*xi*xi+0.4096E3*xi*xi*xi*xi*xi+0.512E3*xi*xi*xi*xi-0.256E3*xi*xi*xi-0.512E2*xi*xi+0.192E2*xi+0.4063492063492063;
      double xi2= 0.2548622222222222E4*xi*xi*xi*xi*xi*xi-0.9557333333333333E3*xi*xi*xi*xi*xi-0.2218666666666667E4*xi*xi*xi*xi+0.7395555555555556E3*xi*xi*xi+0.3605333333333333E3*xi*xi-0.9013333333333333E2*xi-0.32E1;
      double xi3= -0.5097244444444444E4*xi*xi*xi*xi*xi*xi+0.9557333333333333E3*xi*xi*xi*xi*xi+0.4949333333333333E4*xi*xi*xi*xi-0.8248888888888889E3*xi*xi*xi-0.1041066666666667E4*xi*xi+0.1301333333333333E3*xi+0.256E2;
      double xi4= 0.6371555555555556E4*xi*xi*xi*xi*xi*xi-0.64E4*xi*xi*xi*xi+0.1456E4*xi*xi-0.4555555555555556E2;
      double xi5= -0.5097244444444444E4*xi*xi*xi*xi*xi*xi-0.9557333333333333E3*xi*xi*xi*xi*xi+0.4949333333333333E4*xi*xi*xi*xi+0.8248888888888889E3*xi*xi*xi-0.1041066666666667E4*xi*xi-0.1301333333333333E3*xi+0.256E2;
      double xi6= 0.2548622222222222E4*xi*xi*xi*xi*xi*xi+0.9557333333333333E3*xi*xi*xi*xi*xi-0.2218666666666667E4*xi*xi*xi*xi-0.7395555555555556E3*xi*xi*xi+0.3605333333333333E3*xi*xi+0.9013333333333333E2*xi-0.32E1;
      double xi7= -0.7281777777777778E3*xi*xi*xi*xi*xi*xi-0.4096E3*xi*xi*xi*xi*xi+0.512E3*xi*xi*xi*xi+0.256E3*xi*xi*xi-0.512E2*xi*xi-0.192E2*xi+0.4063492063492063;
      double xi8= 0.9102222222222222E2*xi*xi*xi*xi*xi*xi+0.6826666666666667E2*xi*xi*xi*xi*xi-0.4266666666666667E2*xi*xi*xi*xi-0.2844444444444444E2*xi*xi*xi+0.3733333333333333E1*xi*xi+0.1866666666666667E1*xi-0.2857142857142857E-1;

      double eta0= 0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta*eta-0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta*eta+0.1422222222222222E1*eta*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta*eta-0.3111111111111111*eta*eta*eta-0.1428571428571429E-1*eta*eta+0.1428571428571429E-1*eta;
      double eta1= -0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta*eta+0.9752380952380952E1*eta*eta*eta*eta*eta*eta*eta+0.1706666666666667E2*eta*eta*eta*eta*eta*eta-0.128E2*eta*eta*eta*eta*eta-0.4266666666666667E1*eta*eta*eta*eta+0.32E1*eta*eta*eta+0.2031746031746032*eta*eta-0.1523809523809524*eta;
      double eta2= 0.4551111111111111E2*eta*eta*eta*eta*eta*eta*eta*eta-0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta-0.7395555555555556E2*eta*eta*eta*eta*eta*eta+0.3697777777777778E2*eta*eta*eta*eta*eta+0.3004444444444444E2*eta*eta*eta*eta-0.1502222222222222E2*eta*eta*eta-0.16E1*eta*eta+0.8*eta;
      double eta3= -0.9102222222222222E2*eta*eta*eta*eta*eta*eta*eta*eta+0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta+0.1649777777777778E3*eta*eta*eta*eta*eta*eta-0.4124444444444444E2*eta*eta*eta*eta*eta-0.8675555555555556E2*eta*eta*eta*eta+0.2168888888888889E2*eta*eta*eta+0.128E2*eta*eta-0.32E1*eta;
      double eta4= 0.1137777777777778E3*eta*eta*eta*eta*eta*eta*eta*eta-0.2133333333333333E3*eta*eta*eta*eta*eta*eta+0.1213333333333333E3*eta*eta*eta*eta-0.2277777777777778E2*eta*eta+1.0;
      double eta5= -0.9102222222222222E2*eta*eta*eta*eta*eta*eta*eta*eta-0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta+0.1649777777777778E3*eta*eta*eta*eta*eta*eta+0.4124444444444444E2*eta*eta*eta*eta*eta-0.8675555555555556E2*eta*eta*eta*eta-0.2168888888888889E2*eta*eta*eta+0.128E2*eta*eta+0.32E1*eta;
      double eta6= 0.4551111111111111E2*eta*eta*eta*eta*eta*eta*eta*eta+0.2275555555555556E2*eta*eta*eta*eta*eta*eta*eta-0.7395555555555556E2*eta*eta*eta*eta*eta*eta-0.3697777777777778E2*eta*eta*eta*eta*eta+0.3004444444444444E2*eta*eta*eta*eta+0.1502222222222222E2*eta*eta*eta-0.16E1*eta*eta-0.8*eta;
      double eta7= -0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta*eta-0.9752380952380952E1*eta*eta*eta*eta*eta*eta*eta+0.1706666666666667E2*eta*eta*eta*eta*eta*eta+0.128E2*eta*eta*eta*eta*eta-0.4266666666666667E1*eta*eta*eta*eta-0.32E1*eta*eta*eta+0.2031746031746032*eta*eta+0.1523809523809524*eta;
      double eta8= 0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta*eta+0.1625396825396825E1*eta*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta*eta-0.1422222222222222E1*eta*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta*eta+0.3111111111111111*eta*eta*eta-0.1428571428571429E-1*eta*eta-0.1428571428571429E-1*eta;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi7*eta0;
      values[8] = xi8*eta0;
      values[9] = xi0*eta1;
      values[10] = xi1*eta1;
      values[11] = xi2*eta1;
      values[12] = xi3*eta1;
      values[13] = xi4*eta1;
      values[14] = xi5*eta1;
      values[15] = xi6*eta1;
      values[16] = xi7*eta1;
      values[17] = xi8*eta1;
      values[18] = xi0*eta2;
      values[19] = xi1*eta2;
      values[20] = xi2*eta2;
      values[21] = xi3*eta2;
      values[22] = xi4*eta2;
      values[23] = xi5*eta2;
      values[24] = xi6*eta2;
      values[25] = xi7*eta2;
      values[26] = xi8*eta2;
      values[27] = xi0*eta3;
      values[28] = xi1*eta3;
      values[29] = xi2*eta3;
      values[30] = xi3*eta3;
      values[31] = xi4*eta3;
      values[32] = xi5*eta3;
      values[33] = xi6*eta3;
      values[34] = xi7*eta3;
      values[35] = xi8*eta3;
      values[36] = xi0*eta4;
      values[37] = xi1*eta4;
      values[38] = xi2*eta4;
      values[39] = xi3*eta4;
      values[40] = xi4*eta4;
      values[41] = xi5*eta4;
      values[42] = xi6*eta4;
      values[43] = xi7*eta4;
      values[44] = xi8*eta4;
      values[45] = xi0*eta5;
      values[46] = xi1*eta5;
      values[47] = xi2*eta5;
      values[48] = xi3*eta5;
      values[49] = xi4*eta5;
      values[50] = xi5*eta5;
      values[51] = xi6*eta5;
      values[52] = xi7*eta5;
      values[53] = xi8*eta5;
      values[54] = xi0*eta6;
      values[55] = xi1*eta6;
      values[56] = xi2*eta6;
      values[57] = xi3*eta6;
      values[58] = xi4*eta6;
      values[59] = xi5*eta6;
      values[60] = xi6*eta6;
      values[61] = xi7*eta6;
      values[62] = xi8*eta6;
      values[63] = xi0*eta7;
      values[64] = xi1*eta7;
      values[65] = xi2*eta7;
      values[66] = xi3*eta7;
      values[67] = xi4*eta7;
      values[68] = xi5*eta7;
      values[69] = xi6*eta7;
      values[70] = xi7*eta7;
      values[71] = xi8*eta7;
      values[72] = xi0*eta8;
      values[73] = xi1*eta8;
      values[74] = xi2*eta8;
      values[75] = xi3*eta8;
      values[76] = xi4*eta8;
      values[77] = xi5*eta8;
      values[78] = xi6*eta8;
      values[79] = xi7*eta8;
      values[80] = xi8*eta8;
}

// values of the derivatives in xi-eta direction
static void C_Q_Q8_2D_DeriveXiEta(double xi, double eta, double *values)
{

      double xi0= 0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi-0.1137777777777778E2*xi*xi*xi*xi*xi*xi-0.8533333333333333E1*xi*xi*xi*xi*xi+0.7111111111111111E1*xi*xi*xi*xi+0.1244444444444444E1*xi*xi*xi-0.9333333333333333*xi*xi-0.2857142857142857E-1*xi+0.1428571428571429E-1;
      double xi1= -0.1040253968253968E3*xi*xi*xi*xi*xi*xi*xi+0.6826666666666667E2*xi*xi*xi*xi*xi*xi+0.1024E3*xi*xi*xi*xi*xi-0.64E2*xi*xi*xi*xi-0.1706666666666667E2*xi*xi*xi+0.96E1*xi*xi+0.4063492063492063*xi-0.1523809523809524;
      double xi2= 0.3640888888888889E3*xi*xi*xi*xi*xi*xi*xi-0.1592888888888889E3*xi*xi*xi*xi*xi*xi-0.4437333333333333E3*xi*xi*xi*xi*xi+0.1848888888888889E3*xi*xi*xi*xi+0.1201777777777778E3*xi*xi*xi-0.4506666666666667E2*xi*xi-0.32E1*xi+0.8;
      double xi3= -0.7281777777777778E3*xi*xi*xi*xi*xi*xi*xi+0.1592888888888889E3*xi*xi*xi*xi*xi*xi+0.9898666666666667E3*xi*xi*xi*xi*xi-0.2062222222222222E3*xi*xi*xi*xi-0.3470222222222222E3*xi*xi*xi+0.6506666666666667E2*xi*xi+0.256E2*xi-0.32E1;
      double xi4= 0.9102222222222222E3*xi*xi*xi*xi*xi*xi*xi-0.128E4*xi*xi*xi*xi*xi+0.4853333333333333E3*xi*xi*xi-0.4555555555555556E2*xi;
      double xi5= -0.7281777777777778E3*xi*xi*xi*xi*xi*xi*xi-0.1592888888888889E3*xi*xi*xi*xi*xi*xi+0.9898666666666667E3*xi*xi*xi*xi*xi+0.2062222222222222E3*xi*xi*xi*xi-0.3470222222222222E3*xi*xi*xi-0.6506666666666667E2*xi*xi+0.256E2*xi+0.32E1;
      double xi6= 0.3640888888888889E3*xi*xi*xi*xi*xi*xi*xi+0.1592888888888889E3*xi*xi*xi*xi*xi*xi-0.4437333333333333E3*xi*xi*xi*xi*xi-0.1848888888888889E3*xi*xi*xi*xi+0.1201777777777778E3*xi*xi*xi+0.4506666666666667E2*xi*xi-0.32E1*xi-0.8;
      double xi7= -0.1040253968253968E3*xi*xi*xi*xi*xi*xi*xi-0.6826666666666667E2*xi*xi*xi*xi*xi*xi+0.1024E3*xi*xi*xi*xi*xi+0.64E2*xi*xi*xi*xi-0.1706666666666667E2*xi*xi*xi-0.96E1*xi*xi+0.4063492063492063*xi+0.1523809523809524;
      double xi8= 0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi+0.1137777777777778E2*xi*xi*xi*xi*xi*xi-0.8533333333333333E1*xi*xi*xi*xi*xi-0.7111111111111111E1*xi*xi*xi*xi+0.1244444444444444E1*xi*xi*xi+0.9333333333333333*xi*xi-0.2857142857142857E-1*xi-0.1428571428571429E-1;

      double eta0= 0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta-0.1137777777777778E2*eta*eta*eta*eta*eta*eta-0.8533333333333333E1*eta*eta*eta*eta*eta+0.7111111111111111E1*eta*eta*eta*eta+0.1244444444444444E1*eta*eta*eta-0.9333333333333333*eta*eta-0.2857142857142857E-1*eta+0.1428571428571429E-1;
      double eta1= -0.1040253968253968E3*eta*eta*eta*eta*eta*eta*eta+0.6826666666666667E2*eta*eta*eta*eta*eta*eta+0.1024E3*eta*eta*eta*eta*eta-0.64E2*eta*eta*eta*eta-0.1706666666666667E2*eta*eta*eta+0.96E1*eta*eta+0.4063492063492063*eta-0.1523809523809524;
      double eta2= 0.3640888888888889E3*eta*eta*eta*eta*eta*eta*eta-0.1592888888888889E3*eta*eta*eta*eta*eta*eta-0.4437333333333333E3*eta*eta*eta*eta*eta+0.1848888888888889E3*eta*eta*eta*eta+0.1201777777777778E3*eta*eta*eta-0.4506666666666667E2*eta*eta-0.32E1*eta+0.8;
      double eta3= -0.7281777777777778E3*eta*eta*eta*eta*eta*eta*eta+0.1592888888888889E3*eta*eta*eta*eta*eta*eta+0.9898666666666667E3*eta*eta*eta*eta*eta-0.2062222222222222E3*eta*eta*eta*eta-0.3470222222222222E3*eta*eta*eta+0.6506666666666667E2*eta*eta+0.256E2*eta-0.32E1;
      double eta4= 0.9102222222222222E3*eta*eta*eta*eta*eta*eta*eta-0.128E4*eta*eta*eta*eta*eta+0.4853333333333333E3*eta*eta*eta-0.4555555555555556E2*eta;
      double eta5= -0.7281777777777778E3*eta*eta*eta*eta*eta*eta*eta-0.1592888888888889E3*eta*eta*eta*eta*eta*eta+0.9898666666666667E3*eta*eta*eta*eta*eta+0.2062222222222222E3*eta*eta*eta*eta-0.3470222222222222E3*eta*eta*eta-0.6506666666666667E2*eta*eta+0.256E2*eta+0.32E1;
      double eta6= 0.3640888888888889E3*eta*eta*eta*eta*eta*eta*eta+0.1592888888888889E3*eta*eta*eta*eta*eta*eta-0.4437333333333333E3*eta*eta*eta*eta*eta-0.1848888888888889E3*eta*eta*eta*eta+0.1201777777777778E3*eta*eta*eta+0.4506666666666667E2*eta*eta-0.32E1*eta-0.8;
      double eta7= -0.1040253968253968E3*eta*eta*eta*eta*eta*eta*eta-0.6826666666666667E2*eta*eta*eta*eta*eta*eta+0.1024E3*eta*eta*eta*eta*eta+0.64E2*eta*eta*eta*eta-0.1706666666666667E2*eta*eta*eta-0.96E1*eta*eta+0.4063492063492063*eta+0.1523809523809524;
      double eta8= 0.130031746031746E2*eta*eta*eta*eta*eta*eta*eta+0.1137777777777778E2*eta*eta*eta*eta*eta*eta-0.8533333333333333E1*eta*eta*eta*eta*eta-0.7111111111111111E1*eta*eta*eta*eta+0.1244444444444444E1*eta*eta*eta+0.9333333333333333*eta*eta-0.2857142857142857E-1*eta-0.1428571428571429E-1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi7*eta0;
      values[8] = xi8*eta0;
      values[9] = xi0*eta1;
      values[10] = xi1*eta1;
      values[11] = xi2*eta1;
      values[12] = xi3*eta1;
      values[13] = xi4*eta1;
      values[14] = xi5*eta1;
      values[15] = xi6*eta1;
      values[16] = xi7*eta1;
      values[17] = xi8*eta1;
      values[18] = xi0*eta2;
      values[19] = xi1*eta2;
      values[20] = xi2*eta2;
      values[21] = xi3*eta2;
      values[22] = xi4*eta2;
      values[23] = xi5*eta2;
      values[24] = xi6*eta2;
      values[25] = xi7*eta2;
      values[26] = xi8*eta2;
      values[27] = xi0*eta3;
      values[28] = xi1*eta3;
      values[29] = xi2*eta3;
      values[30] = xi3*eta3;
      values[31] = xi4*eta3;
      values[32] = xi5*eta3;
      values[33] = xi6*eta3;
      values[34] = xi7*eta3;
      values[35] = xi8*eta3;
      values[36] = xi0*eta4;
      values[37] = xi1*eta4;
      values[38] = xi2*eta4;
      values[39] = xi3*eta4;
      values[40] = xi4*eta4;
      values[41] = xi5*eta4;
      values[42] = xi6*eta4;
      values[43] = xi7*eta4;
      values[44] = xi8*eta4;
      values[45] = xi0*eta5;
      values[46] = xi1*eta5;
      values[47] = xi2*eta5;
      values[48] = xi3*eta5;
      values[49] = xi4*eta5;
      values[50] = xi5*eta5;
      values[51] = xi6*eta5;
      values[52] = xi7*eta5;
      values[53] = xi8*eta5;
      values[54] = xi0*eta6;
      values[55] = xi1*eta6;
      values[56] = xi2*eta6;
      values[57] = xi3*eta6;
      values[58] = xi4*eta6;
      values[59] = xi5*eta6;
      values[60] = xi6*eta6;
      values[61] = xi7*eta6;
      values[62] = xi8*eta6;
      values[63] = xi0*eta7;
      values[64] = xi1*eta7;
      values[65] = xi2*eta7;
      values[66] = xi3*eta7;
      values[67] = xi4*eta7;
      values[68] = xi5*eta7;
      values[69] = xi6*eta7;
      values[70] = xi7*eta7;
      values[71] = xi8*eta7;
      values[72] = xi0*eta8;
      values[73] = xi1*eta8;
      values[74] = xi2*eta8;
      values[75] = xi3*eta8;
      values[76] = xi4*eta8;
      values[77] = xi5*eta8;
      values[78] = xi6*eta8;
      values[79] = xi7*eta8;
      values[80] = xi8*eta8;
}

// values of the derivatives in eta-eta direction
static void C_Q_Q8_2D_DeriveEtaEta(double xi, double eta, double *values)
{

      double xi0= 0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi*xi-0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi*xi+0.1422222222222222E1*xi*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi*xi-0.3111111111111111*xi*xi*xi-0.1428571428571429E-1*xi*xi+0.1428571428571429E-1*xi;
      double xi1= -0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi*xi+0.9752380952380952E1*xi*xi*xi*xi*xi*xi*xi+0.1706666666666667E2*xi*xi*xi*xi*xi*xi-0.128E2*xi*xi*xi*xi*xi-0.4266666666666667E1*xi*xi*xi*xi+0.32E1*xi*xi*xi+0.2031746031746032*xi*xi-0.1523809523809524*xi;
      double xi2= 0.4551111111111111E2*xi*xi*xi*xi*xi*xi*xi*xi-0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi-0.7395555555555556E2*xi*xi*xi*xi*xi*xi+0.3697777777777778E2*xi*xi*xi*xi*xi+0.3004444444444444E2*xi*xi*xi*xi-0.1502222222222222E2*xi*xi*xi-0.16E1*xi*xi+0.8*xi;
      double xi3= -0.9102222222222222E2*xi*xi*xi*xi*xi*xi*xi*xi+0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi+0.1649777777777778E3*xi*xi*xi*xi*xi*xi-0.4124444444444444E2*xi*xi*xi*xi*xi-0.8675555555555556E2*xi*xi*xi*xi+0.2168888888888889E2*xi*xi*xi+0.128E2*xi*xi-0.32E1*xi;
      double xi4= 0.1137777777777778E3*xi*xi*xi*xi*xi*xi*xi*xi-0.2133333333333333E3*xi*xi*xi*xi*xi*xi+0.1213333333333333E3*xi*xi*xi*xi-0.2277777777777778E2*xi*xi+1.0;
      double xi5= -0.9102222222222222E2*xi*xi*xi*xi*xi*xi*xi*xi-0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi+0.1649777777777778E3*xi*xi*xi*xi*xi*xi+0.4124444444444444E2*xi*xi*xi*xi*xi-0.8675555555555556E2*xi*xi*xi*xi-0.2168888888888889E2*xi*xi*xi+0.128E2*xi*xi+0.32E1*xi;
      double xi6= 0.4551111111111111E2*xi*xi*xi*xi*xi*xi*xi*xi+0.2275555555555556E2*xi*xi*xi*xi*xi*xi*xi-0.7395555555555556E2*xi*xi*xi*xi*xi*xi-0.3697777777777778E2*xi*xi*xi*xi*xi+0.3004444444444444E2*xi*xi*xi*xi+0.1502222222222222E2*xi*xi*xi-0.16E1*xi*xi-0.8*xi;
      double xi7= -0.130031746031746E2*xi*xi*xi*xi*xi*xi*xi*xi-0.9752380952380952E1*xi*xi*xi*xi*xi*xi*xi+0.1706666666666667E2*xi*xi*xi*xi*xi*xi+0.128E2*xi*xi*xi*xi*xi-0.4266666666666667E1*xi*xi*xi*xi-0.32E1*xi*xi*xi+0.2031746031746032*xi*xi+0.1523809523809524*xi;
      double xi8= 0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi*xi+0.1625396825396825E1*xi*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi*xi-0.1422222222222222E1*xi*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi*xi+0.3111111111111111*xi*xi*xi-0.1428571428571429E-1*xi*xi-0.1428571428571429E-1*xi;

      double eta0= 0.9102222222222222E2*eta*eta*eta*eta*eta*eta-0.6826666666666667E2*eta*eta*eta*eta*eta-0.4266666666666667E2*eta*eta*eta*eta+0.2844444444444444E2*eta*eta*eta+0.3733333333333333E1*eta*eta-0.1866666666666667E1*eta-0.2857142857142857E-1;
      double eta1= -0.7281777777777778E3*eta*eta*eta*eta*eta*eta+0.4096E3*eta*eta*eta*eta*eta+0.512E3*eta*eta*eta*eta-0.256E3*eta*eta*eta-0.512E2*eta*eta+0.192E2*eta+0.4063492063492063;
      double eta2= 0.2548622222222222E4*eta*eta*eta*eta*eta*eta-0.9557333333333333E3*eta*eta*eta*eta*eta-0.2218666666666667E4*eta*eta*eta*eta+0.7395555555555556E3*eta*eta*eta+0.3605333333333333E3*eta*eta-0.9013333333333333E2*eta-0.32E1;
      double eta3= -0.5097244444444444E4*eta*eta*eta*eta*eta*eta+0.9557333333333333E3*eta*eta*eta*eta*eta+0.4949333333333333E4*eta*eta*eta*eta-0.8248888888888889E3*eta*eta*eta-0.1041066666666667E4*eta*eta+0.1301333333333333E3*eta+0.256E2;
      double eta4= 0.6371555555555556E4*eta*eta*eta*eta*eta*eta-0.64E4*eta*eta*eta*eta+0.1456E4*eta*eta-0.4555555555555556E2;
      double eta5= -0.5097244444444444E4*eta*eta*eta*eta*eta*eta-0.9557333333333333E3*eta*eta*eta*eta*eta+0.4949333333333333E4*eta*eta*eta*eta+0.8248888888888889E3*eta*eta*eta-0.1041066666666667E4*eta*eta-0.1301333333333333E3*eta+0.256E2;
      double eta6= 0.2548622222222222E4*eta*eta*eta*eta*eta*eta+0.9557333333333333E3*eta*eta*eta*eta*eta-0.2218666666666667E4*eta*eta*eta*eta-0.7395555555555556E3*eta*eta*eta+0.3605333333333333E3*eta*eta+0.9013333333333333E2*eta-0.32E1;
      double eta7= -0.7281777777777778E3*eta*eta*eta*eta*eta*eta-0.4096E3*eta*eta*eta*eta*eta+0.512E3*eta*eta*eta*eta+0.256E3*eta*eta*eta-0.512E2*eta*eta-0.192E2*eta+0.4063492063492063;
      double eta8= 0.9102222222222222E2*eta*eta*eta*eta*eta*eta+0.6826666666666667E2*eta*eta*eta*eta*eta-0.4266666666666667E2*eta*eta*eta*eta-0.2844444444444444E2*eta*eta*eta+0.3733333333333333E1*eta*eta+0.1866666666666667E1*eta-0.2857142857142857E-1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi7*eta0;
      values[8] = xi8*eta0;
      values[9] = xi0*eta1;
      values[10] = xi1*eta1;
      values[11] = xi2*eta1;
      values[12] = xi3*eta1;
      values[13] = xi4*eta1;
      values[14] = xi5*eta1;
      values[15] = xi6*eta1;
      values[16] = xi7*eta1;
      values[17] = xi8*eta1;
      values[18] = xi0*eta2;
      values[19] = xi1*eta2;
      values[20] = xi2*eta2;
      values[21] = xi3*eta2;
      values[22] = xi4*eta2;
      values[23] = xi5*eta2;
      values[24] = xi6*eta2;
      values[25] = xi7*eta2;
      values[26] = xi8*eta2;
      values[27] = xi0*eta3;
      values[28] = xi1*eta3;
      values[29] = xi2*eta3;
      values[30] = xi3*eta3;
      values[31] = xi4*eta3;
      values[32] = xi5*eta3;
      values[33] = xi6*eta3;
      values[34] = xi7*eta3;
      values[35] = xi8*eta3;
      values[36] = xi0*eta4;
      values[37] = xi1*eta4;
      values[38] = xi2*eta4;
      values[39] = xi3*eta4;
      values[40] = xi4*eta4;
      values[41] = xi5*eta4;
      values[42] = xi6*eta4;
      values[43] = xi7*eta4;
      values[44] = xi8*eta4;
      values[45] = xi0*eta5;
      values[46] = xi1*eta5;
      values[47] = xi2*eta5;
      values[48] = xi3*eta5;
      values[49] = xi4*eta5;
      values[50] = xi5*eta5;
      values[51] = xi6*eta5;
      values[52] = xi7*eta5;
      values[53] = xi8*eta5;
      values[54] = xi0*eta6;
      values[55] = xi1*eta6;
      values[56] = xi2*eta6;
      values[57] = xi3*eta6;
      values[58] = xi4*eta6;
      values[59] = xi5*eta6;
      values[60] = xi6*eta6;
      values[61] = xi7*eta6;
      values[62] = xi8*eta6;
      values[63] = xi0*eta7;
      values[64] = xi1*eta7;
      values[65] = xi2*eta7;
      values[66] = xi3*eta7;
      values[67] = xi4*eta7;
      values[68] = xi5*eta7;
      values[69] = xi6*eta7;
      values[70] = xi7*eta7;
      values[71] = xi8*eta7;
      values[72] = xi0*eta8;
      values[73] = xi1*eta8;
      values[74] = xi2*eta8;
      values[75] = xi3*eta8;
      values[76] = xi4*eta8;
      values[77] = xi5*eta8;
      values[78] = xi6*eta8;
      values[79] = xi7*eta8;
      values[80] = xi8*eta8;
}


// ***********************************************************************

TBaseFunct2D *BF_C_Q_Q8_2D_Obj = new TBaseFunct2D
        (81, BF_C_Q_Q8_2D, BFUnitSquare, 
         C_Q_Q8_2D_Funct, C_Q_Q8_2D_DeriveXi,
         C_Q_Q8_2D_DeriveEta, C_Q_Q8_2D_DeriveXiXi,
         C_Q_Q8_2D_DeriveXiEta, C_Q_Q8_2D_DeriveEtaEta, 8, 8,
         0, NULL);
