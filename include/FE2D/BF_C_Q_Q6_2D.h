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
// Q6 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_Q6_2D_Funct(double xi, double eta, double *values)
{  

      double xi0= 0.10125E1*xi*xi*xi*xi*xi*xi-0.10125E1*xi*xi*xi*xi*xi-0.5625*xi*xi*xi*xi+0.5625*xi*xi*xi+0.5E-1*xi*xi-0.5E-1*xi;
      double xi1= -0.6075E1*xi*xi*xi*xi*xi*xi+0.405E1*xi*xi*xi*xi*xi+0.675E1*xi*xi*xi*xi-0.45E1*xi*xi*xi-0.675*xi*xi+0.45*xi;
      double xi2= 0.151875E2*xi*xi*xi*xi*xi*xi-0.50625E1*xi*xi*xi*xi*xi-0.219375E2*xi*xi*xi*xi+0.73125E1*xi*xi*xi+0.675E1*xi*xi-0.225E1*xi;
      double xi3= -0.2025E2*xi*xi*xi*xi*xi*xi+0.315E2*xi*xi*xi*xi-0.1225E2*xi*xi+1.0;
      double xi4= 0.151875E2*xi*xi*xi*xi*xi*xi+0.50625E1*xi*xi*xi*xi*xi-0.219375E2*xi*xi*xi*xi-0.73125E1*xi*xi*xi+0.675E1*xi*xi+0.225E1*xi;
      double xi5= -0.6075E1*xi*xi*xi*xi*xi*xi-0.405E1*xi*xi*xi*xi*xi+0.675E1*xi*xi*xi*xi+0.45E1*xi*xi*xi-0.675*xi*xi-0.45*xi;
      double xi6= 0.10125E1*xi*xi*xi*xi*xi*xi+0.10125E1*xi*xi*xi*xi*xi-0.5625*xi*xi*xi*xi-0.5625*xi*xi*xi+0.5E-1*xi*xi+0.5E-1*xi;

      double eta0= 0.10125E1*eta*eta*eta*eta*eta*eta-0.10125E1*eta*eta*eta*eta*eta-0.5625*eta*eta*eta*eta+0.5625*eta*eta*eta+0.5E-1*eta*eta-0.5E-1*eta;
      double eta1= -0.6075E1*eta*eta*eta*eta*eta*eta+0.405E1*eta*eta*eta*eta*eta+0.675E1*eta*eta*eta*eta-0.45E1*eta*eta*eta-0.675*eta*eta+0.45*eta;
      double eta2= 0.151875E2*eta*eta*eta*eta*eta*eta-0.50625E1*eta*eta*eta*eta*eta-0.219375E2*eta*eta*eta*eta+0.73125E1*eta*eta*eta+0.675E1*eta*eta-0.225E1*eta;
      double eta3= -0.2025E2*eta*eta*eta*eta*eta*eta+0.315E2*eta*eta*eta*eta-0.1225E2*eta*eta+1.0;
      double eta4= 0.151875E2*eta*eta*eta*eta*eta*eta+0.50625E1*eta*eta*eta*eta*eta-0.219375E2*eta*eta*eta*eta-0.73125E1*eta*eta*eta+0.675E1*eta*eta+0.225E1*eta;
      double eta5= -0.6075E1*eta*eta*eta*eta*eta*eta-0.405E1*eta*eta*eta*eta*eta+0.675E1*eta*eta*eta*eta+0.45E1*eta*eta*eta-0.675*eta*eta-0.45*eta;
      double eta6= 0.10125E1*eta*eta*eta*eta*eta*eta+0.10125E1*eta*eta*eta*eta*eta-0.5625*eta*eta*eta*eta-0.5625*eta*eta*eta+0.5E-1*eta*eta+0.5E-1*eta;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi0*eta1;
      values[8] = xi1*eta1;
      values[9] = xi2*eta1;
      values[10] = xi3*eta1;
      values[11] = xi4*eta1;
      values[12] = xi5*eta1;
      values[13] = xi6*eta1;
      values[14] = xi0*eta2;
      values[15] = xi1*eta2;
      values[16] = xi2*eta2;
      values[17] = xi3*eta2;
      values[18] = xi4*eta2;
      values[19] = xi5*eta2;
      values[20] = xi6*eta2;
      values[21] = xi0*eta3;
      values[22] = xi1*eta3;
      values[23] = xi2*eta3;
      values[24] = xi3*eta3;
      values[25] = xi4*eta3;
      values[26] = xi5*eta3;
      values[27] = xi6*eta3;
      values[28] = xi0*eta4;
      values[29] = xi1*eta4;
      values[30] = xi2*eta4;
      values[31] = xi3*eta4;
      values[32] = xi4*eta4;
      values[33] = xi5*eta4;
      values[34] = xi6*eta4;
      values[35] = xi0*eta5;
      values[36] = xi1*eta5;
      values[37] = xi2*eta5;
      values[38] = xi3*eta5;
      values[39] = xi4*eta5;
      values[40] = xi5*eta5;
      values[41] = xi6*eta5;
      values[42] = xi0*eta6;
      values[43] = xi1*eta6;
      values[44] = xi2*eta6;
      values[45] = xi3*eta6;
      values[46] = xi4*eta6;
      values[47] = xi5*eta6;
      values[48] = xi6*eta6;
}


// values of the derivatives in xi direction
static void C_Q_Q6_2D_DeriveXi(double xi, double eta, double *values)
{

      double xi0= 0.6075E1*xi*xi*xi*xi*xi-0.50625E1*xi*xi*xi*xi-0.225E1*xi*xi*xi+0.16875E1*xi*xi+0.1*xi-0.5E-1;
      double xi1= -0.3645E2*xi*xi*xi*xi*xi+0.2025E2*xi*xi*xi*xi+0.27E2*xi*xi*xi-0.135E2*xi*xi-0.135E1*xi+0.45;
      double xi2= 0.91125E2*xi*xi*xi*xi*xi-0.253125E2*xi*xi*xi*xi-0.8775E2*xi*xi*xi+0.219375E2*xi*xi+0.135E2*xi-0.225E1;
      double xi3= -0.1215E3*xi*xi*xi*xi*xi+0.126E3*xi*xi*xi-0.245E2*xi;
      double xi4= 0.91125E2*xi*xi*xi*xi*xi+0.253125E2*xi*xi*xi*xi-0.8775E2*xi*xi*xi-0.219375E2*xi*xi+0.135E2*xi+0.225E1;
      double xi5= -0.3645E2*xi*xi*xi*xi*xi-0.2025E2*xi*xi*xi*xi+0.27E2*xi*xi*xi+0.135E2*xi*xi-0.135E1*xi-0.45;
      double xi6= 0.6075E1*xi*xi*xi*xi*xi+0.50625E1*xi*xi*xi*xi-0.225E1*xi*xi*xi-0.16875E1*xi*xi+0.1*xi+0.5E-1;

      double eta0= 0.10125E1*eta*eta*eta*eta*eta*eta-0.10125E1*eta*eta*eta*eta*eta-0.5625*eta*eta*eta*eta+0.5625*eta*eta*eta+0.5E-1*eta*eta-0.5E-1*eta;
      double eta1= -0.6075E1*eta*eta*eta*eta*eta*eta+0.405E1*eta*eta*eta*eta*eta+0.675E1*eta*eta*eta*eta-0.45E1*eta*eta*eta-0.675*eta*eta+0.45*eta;
      double eta2= 0.151875E2*eta*eta*eta*eta*eta*eta-0.50625E1*eta*eta*eta*eta*eta-0.219375E2*eta*eta*eta*eta+0.73125E1*eta*eta*eta+0.675E1*eta*eta-0.225E1*eta;
      double eta3= -0.2025E2*eta*eta*eta*eta*eta*eta+0.315E2*eta*eta*eta*eta-0.1225E2*eta*eta+1.0;
      double eta4= 0.151875E2*eta*eta*eta*eta*eta*eta+0.50625E1*eta*eta*eta*eta*eta-0.219375E2*eta*eta*eta*eta-0.73125E1*eta*eta*eta+0.675E1*eta*eta+0.225E1*eta;
      double eta5= -0.6075E1*eta*eta*eta*eta*eta*eta-0.405E1*eta*eta*eta*eta*eta+0.675E1*eta*eta*eta*eta+0.45E1*eta*eta*eta-0.675*eta*eta-0.45*eta;
      double eta6= 0.10125E1*eta*eta*eta*eta*eta*eta+0.10125E1*eta*eta*eta*eta*eta-0.5625*eta*eta*eta*eta-0.5625*eta*eta*eta+0.5E-1*eta*eta+0.5E-1*eta;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi0*eta1;
      values[8] = xi1*eta1;
      values[9] = xi2*eta1;
      values[10] = xi3*eta1;
      values[11] = xi4*eta1;
      values[12] = xi5*eta1;
      values[13] = xi6*eta1;
      values[14] = xi0*eta2;
      values[15] = xi1*eta2;
      values[16] = xi2*eta2;
      values[17] = xi3*eta2;
      values[18] = xi4*eta2;
      values[19] = xi5*eta2;
      values[20] = xi6*eta2;
      values[21] = xi0*eta3;
      values[22] = xi1*eta3;
      values[23] = xi2*eta3;
      values[24] = xi3*eta3;
      values[25] = xi4*eta3;
      values[26] = xi5*eta3;
      values[27] = xi6*eta3;
      values[28] = xi0*eta4;
      values[29] = xi1*eta4;
      values[30] = xi2*eta4;
      values[31] = xi3*eta4;
      values[32] = xi4*eta4;
      values[33] = xi5*eta4;
      values[34] = xi6*eta4;
      values[35] = xi0*eta5;
      values[36] = xi1*eta5;
      values[37] = xi2*eta5;
      values[38] = xi3*eta5;
      values[39] = xi4*eta5;
      values[40] = xi5*eta5;
      values[41] = xi6*eta5;
      values[42] = xi0*eta6;
      values[43] = xi1*eta6;
      values[44] = xi2*eta6;
      values[45] = xi3*eta6;
      values[46] = xi4*eta6;
      values[47] = xi5*eta6;
      values[48] = xi6*eta6;
}

// values of the derivatives in eta direction
static void C_Q_Q6_2D_DeriveEta(double xi, double eta, double *values)
{

      double xi0= 0.10125E1*xi*xi*xi*xi*xi*xi-0.10125E1*xi*xi*xi*xi*xi-0.5625*xi*xi*xi*xi+0.5625*xi*xi*xi+0.5E-1*xi*xi-0.5E-1*xi;
      double xi1= -0.6075E1*xi*xi*xi*xi*xi*xi+0.405E1*xi*xi*xi*xi*xi+0.675E1*xi*xi*xi*xi-0.45E1*xi*xi*xi-0.675*xi*xi+0.45*xi;
      double xi2= 0.151875E2*xi*xi*xi*xi*xi*xi-0.50625E1*xi*xi*xi*xi*xi-0.219375E2*xi*xi*xi*xi+0.73125E1*xi*xi*xi+0.675E1*xi*xi-0.225E1*xi;
      double xi3= -0.2025E2*xi*xi*xi*xi*xi*xi+0.315E2*xi*xi*xi*xi-0.1225E2*xi*xi+1.0;
      double xi4= 0.151875E2*xi*xi*xi*xi*xi*xi+0.50625E1*xi*xi*xi*xi*xi-0.219375E2*xi*xi*xi*xi-0.73125E1*xi*xi*xi+0.675E1*xi*xi+0.225E1*xi;
      double xi5= -0.6075E1*xi*xi*xi*xi*xi*xi-0.405E1*xi*xi*xi*xi*xi+0.675E1*xi*xi*xi*xi+0.45E1*xi*xi*xi-0.675*xi*xi-0.45*xi;
      double xi6= 0.10125E1*xi*xi*xi*xi*xi*xi+0.10125E1*xi*xi*xi*xi*xi-0.5625*xi*xi*xi*xi-0.5625*xi*xi*xi+0.5E-1*xi*xi+0.5E-1*xi;

      double eta0= 0.6075E1*eta*eta*eta*eta*eta-0.50625E1*eta*eta*eta*eta-0.225E1*eta*eta*eta+0.16875E1*eta*eta+0.1*eta-0.5E-1;
      double eta1= -0.3645E2*eta*eta*eta*eta*eta+0.2025E2*eta*eta*eta*eta+0.27E2*eta*eta*eta-0.135E2*eta*eta-0.135E1*eta+0.45;
      double eta2= 0.91125E2*eta*eta*eta*eta*eta-0.253125E2*eta*eta*eta*eta-0.8775E2*eta*eta*eta+0.219375E2*eta*eta+0.135E2*eta-0.225E1;
      double eta3= -0.1215E3*eta*eta*eta*eta*eta+0.126E3*eta*eta*eta-0.245E2*eta;
      double eta4= 0.91125E2*eta*eta*eta*eta*eta+0.253125E2*eta*eta*eta*eta-0.8775E2*eta*eta*eta-0.219375E2*eta*eta+0.135E2*eta+0.225E1;
      double eta5= -0.3645E2*eta*eta*eta*eta*eta-0.2025E2*eta*eta*eta*eta+0.27E2*eta*eta*eta+0.135E2*eta*eta-0.135E1*eta-0.45;
      double eta6= 0.6075E1*eta*eta*eta*eta*eta+0.50625E1*eta*eta*eta*eta-0.225E1*eta*eta*eta-0.16875E1*eta*eta+0.1*eta+0.5E-1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi0*eta1;
      values[8] = xi1*eta1;
      values[9] = xi2*eta1;
      values[10] = xi3*eta1;
      values[11] = xi4*eta1;
      values[12] = xi5*eta1;
      values[13] = xi6*eta1;
      values[14] = xi0*eta2;
      values[15] = xi1*eta2;
      values[16] = xi2*eta2;
      values[17] = xi3*eta2;
      values[18] = xi4*eta2;
      values[19] = xi5*eta2;
      values[20] = xi6*eta2;
      values[21] = xi0*eta3;
      values[22] = xi1*eta3;
      values[23] = xi2*eta3;
      values[24] = xi3*eta3;
      values[25] = xi4*eta3;
      values[26] = xi5*eta3;
      values[27] = xi6*eta3;
      values[28] = xi0*eta4;
      values[29] = xi1*eta4;
      values[30] = xi2*eta4;
      values[31] = xi3*eta4;
      values[32] = xi4*eta4;
      values[33] = xi5*eta4;
      values[34] = xi6*eta4;
      values[35] = xi0*eta5;
      values[36] = xi1*eta5;
      values[37] = xi2*eta5;
      values[38] = xi3*eta5;
      values[39] = xi4*eta5;
      values[40] = xi5*eta5;
      values[41] = xi6*eta5;
      values[42] = xi0*eta6;
      values[43] = xi1*eta6;
      values[44] = xi2*eta6;
      values[45] = xi3*eta6;
      values[46] = xi4*eta6;
      values[47] = xi5*eta6;
      values[48] = xi6*eta6;
}

// values of the derivatives in xi-xi  direction
static void C_Q_Q6_2D_DeriveXiXi(double xi, double eta, double *values)
{

      double xi0= 0.30375E2*xi*xi*xi*xi-0.2025E2*xi*xi*xi-0.675E1*xi*xi+0.3375E1*xi+0.1;
      double xi1= -0.18225E3*xi*xi*xi*xi+0.81E2*xi*xi*xi+0.81E2*xi*xi-0.27E2*xi-0.135E1;
      double xi2= 0.455625E3*xi*xi*xi*xi-0.10125E3*xi*xi*xi-0.26325E3*xi*xi+0.43875E2*xi+0.135E2;
      double xi3= -0.6075E3*xi*xi*xi*xi+0.378E3*xi*xi-0.245E2;
      double xi4= 0.455625E3*xi*xi*xi*xi+0.10125E3*xi*xi*xi-0.26325E3*xi*xi-0.43875E2*xi+0.135E2;
      double xi5= -0.18225E3*xi*xi*xi*xi-0.81E2*xi*xi*xi+0.81E2*xi*xi+0.27E2*xi-0.135E1;
      double xi6= 0.30375E2*xi*xi*xi*xi+0.2025E2*xi*xi*xi-0.675E1*xi*xi-0.3375E1*xi+0.1;

      double eta0= 0.10125E1*eta*eta*eta*eta*eta*eta-0.10125E1*eta*eta*eta*eta*eta-0.5625*eta*eta*eta*eta+0.5625*eta*eta*eta+0.5E-1*eta*eta-0.5E-1*eta;
      double eta1= -0.6075E1*eta*eta*eta*eta*eta*eta+0.405E1*eta*eta*eta*eta*eta+0.675E1*eta*eta*eta*eta-0.45E1*eta*eta*eta-0.675*eta*eta+0.45*eta;
      double eta2= 0.151875E2*eta*eta*eta*eta*eta*eta-0.50625E1*eta*eta*eta*eta*eta-0.219375E2*eta*eta*eta*eta+0.73125E1*eta*eta*eta+0.675E1*eta*eta-0.225E1*eta;
      double eta3= -0.2025E2*eta*eta*eta*eta*eta*eta+0.315E2*eta*eta*eta*eta-0.1225E2*eta*eta+1.0;
      double eta4= 0.151875E2*eta*eta*eta*eta*eta*eta+0.50625E1*eta*eta*eta*eta*eta-0.219375E2*eta*eta*eta*eta-0.73125E1*eta*eta*eta+0.675E1*eta*eta+0.225E1*eta;
      double eta5= -0.6075E1*eta*eta*eta*eta*eta*eta-0.405E1*eta*eta*eta*eta*eta+0.675E1*eta*eta*eta*eta+0.45E1*eta*eta*eta-0.675*eta*eta-0.45*eta;
      double eta6= 0.10125E1*eta*eta*eta*eta*eta*eta+0.10125E1*eta*eta*eta*eta*eta-0.5625*eta*eta*eta*eta-0.5625*eta*eta*eta+0.5E-1*eta*eta+0.5E-1*eta;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi0*eta1;
      values[8] = xi1*eta1;
      values[9] = xi2*eta1;
      values[10] = xi3*eta1;
      values[11] = xi4*eta1;
      values[12] = xi5*eta1;
      values[13] = xi6*eta1;
      values[14] = xi0*eta2;
      values[15] = xi1*eta2;
      values[16] = xi2*eta2;
      values[17] = xi3*eta2;
      values[18] = xi4*eta2;
      values[19] = xi5*eta2;
      values[20] = xi6*eta2;
      values[21] = xi0*eta3;
      values[22] = xi1*eta3;
      values[23] = xi2*eta3;
      values[24] = xi3*eta3;
      values[25] = xi4*eta3;
      values[26] = xi5*eta3;
      values[27] = xi6*eta3;
      values[28] = xi0*eta4;
      values[29] = xi1*eta4;
      values[30] = xi2*eta4;
      values[31] = xi3*eta4;
      values[32] = xi4*eta4;
      values[33] = xi5*eta4;
      values[34] = xi6*eta4;
      values[35] = xi0*eta5;
      values[36] = xi1*eta5;
      values[37] = xi2*eta5;
      values[38] = xi3*eta5;
      values[39] = xi4*eta5;
      values[40] = xi5*eta5;
      values[41] = xi6*eta5;
      values[42] = xi0*eta6;
      values[43] = xi1*eta6;
      values[44] = xi2*eta6;
      values[45] = xi3*eta6;
      values[46] = xi4*eta6;
      values[47] = xi5*eta6;
      values[48] = xi6*eta6;
}

// values of the derivatives in xi-eta direction
static void C_Q_Q6_2D_DeriveXiEta(double xi, double eta, double *values)
{

      double xi0= 0.6075E1*xi*xi*xi*xi*xi-0.50625E1*xi*xi*xi*xi-0.225E1*xi*xi*xi+0.16875E1*xi*xi+0.1*xi-0.5E-1;
      double xi1= -0.3645E2*xi*xi*xi*xi*xi+0.2025E2*xi*xi*xi*xi+0.27E2*xi*xi*xi-0.135E2*xi*xi-0.135E1*xi+0.45;
      double xi2= 0.91125E2*xi*xi*xi*xi*xi-0.253125E2*xi*xi*xi*xi-0.8775E2*xi*xi*xi+0.219375E2*xi*xi+0.135E2*xi-0.225E1;
      double xi3= -0.1215E3*xi*xi*xi*xi*xi+0.126E3*xi*xi*xi-0.245E2*xi;
      double xi4= 0.91125E2*xi*xi*xi*xi*xi+0.253125E2*xi*xi*xi*xi-0.8775E2*xi*xi*xi-0.219375E2*xi*xi+0.135E2*xi+0.225E1;
      double xi5= -0.3645E2*xi*xi*xi*xi*xi-0.2025E2*xi*xi*xi*xi+0.27E2*xi*xi*xi+0.135E2*xi*xi-0.135E1*xi-0.45;
      double xi6= 0.6075E1*xi*xi*xi*xi*xi+0.50625E1*xi*xi*xi*xi-0.225E1*xi*xi*xi-0.16875E1*xi*xi+0.1*xi+0.5E-1;

      double eta0= 0.6075E1*eta*eta*eta*eta*eta-0.50625E1*eta*eta*eta*eta-0.225E1*eta*eta*eta+0.16875E1*eta*eta+0.1*eta-0.5E-1;
      double eta1= -0.3645E2*eta*eta*eta*eta*eta+0.2025E2*eta*eta*eta*eta+0.27E2*eta*eta*eta-0.135E2*eta*eta-0.135E1*eta+0.45;
      double eta2= 0.91125E2*eta*eta*eta*eta*eta-0.253125E2*eta*eta*eta*eta-0.8775E2*eta*eta*eta+0.219375E2*eta*eta+0.135E2*eta-0.225E1;
      double eta3= -0.1215E3*eta*eta*eta*eta*eta+0.126E3*eta*eta*eta-0.245E2*eta;
      double eta4= 0.91125E2*eta*eta*eta*eta*eta+0.253125E2*eta*eta*eta*eta-0.8775E2*eta*eta*eta-0.219375E2*eta*eta+0.135E2*eta+0.225E1;
      double eta5= -0.3645E2*eta*eta*eta*eta*eta-0.2025E2*eta*eta*eta*eta+0.27E2*eta*eta*eta+0.135E2*eta*eta-0.135E1*eta-0.45;
      double eta6= 0.6075E1*eta*eta*eta*eta*eta+0.50625E1*eta*eta*eta*eta-0.225E1*eta*eta*eta-0.16875E1*eta*eta+0.1*eta+0.5E-1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi0*eta1;
      values[8] = xi1*eta1;
      values[9] = xi2*eta1;
      values[10] = xi3*eta1;
      values[11] = xi4*eta1;
      values[12] = xi5*eta1;
      values[13] = xi6*eta1;
      values[14] = xi0*eta2;
      values[15] = xi1*eta2;
      values[16] = xi2*eta2;
      values[17] = xi3*eta2;
      values[18] = xi4*eta2;
      values[19] = xi5*eta2;
      values[20] = xi6*eta2;
      values[21] = xi0*eta3;
      values[22] = xi1*eta3;
      values[23] = xi2*eta3;
      values[24] = xi3*eta3;
      values[25] = xi4*eta3;
      values[26] = xi5*eta3;
      values[27] = xi6*eta3;
      values[28] = xi0*eta4;
      values[29] = xi1*eta4;
      values[30] = xi2*eta4;
      values[31] = xi3*eta4;
      values[32] = xi4*eta4;
      values[33] = xi5*eta4;
      values[34] = xi6*eta4;
      values[35] = xi0*eta5;
      values[36] = xi1*eta5;
      values[37] = xi2*eta5;
      values[38] = xi3*eta5;
      values[39] = xi4*eta5;
      values[40] = xi5*eta5;
      values[41] = xi6*eta5;
      values[42] = xi0*eta6;
      values[43] = xi1*eta6;
      values[44] = xi2*eta6;
      values[45] = xi3*eta6;
      values[46] = xi4*eta6;
      values[47] = xi5*eta6;
      values[48] = xi6*eta6;
}

// values of the derivatives in eta-eta direction
static void C_Q_Q6_2D_DeriveEtaEta(double xi, double eta, double *values)
{

      double xi0= 0.10125E1*xi*xi*xi*xi*xi*xi-0.10125E1*xi*xi*xi*xi*xi-0.5625*xi*xi*xi*xi+0.5625*xi*xi*xi+0.5E-1*xi*xi-0.5E-1*xi;
      double xi1= -0.6075E1*xi*xi*xi*xi*xi*xi+0.405E1*xi*xi*xi*xi*xi+0.675E1*xi*xi*xi*xi-0.45E1*xi*xi*xi-0.675*xi*xi+0.45*xi;
      double xi2= 0.151875E2*xi*xi*xi*xi*xi*xi-0.50625E1*xi*xi*xi*xi*xi-0.219375E2*xi*xi*xi*xi+0.73125E1*xi*xi*xi+0.675E1*xi*xi-0.225E1*xi;
      double xi3= -0.2025E2*xi*xi*xi*xi*xi*xi+0.315E2*xi*xi*xi*xi-0.1225E2*xi*xi+1.0;
      double xi4= 0.151875E2*xi*xi*xi*xi*xi*xi+0.50625E1*xi*xi*xi*xi*xi-0.219375E2*xi*xi*xi*xi-0.73125E1*xi*xi*xi+0.675E1*xi*xi+0.225E1*xi;
      double xi5= -0.6075E1*xi*xi*xi*xi*xi*xi-0.405E1*xi*xi*xi*xi*xi+0.675E1*xi*xi*xi*xi+0.45E1*xi*xi*xi-0.675*xi*xi-0.45*xi;
      double xi6= 0.10125E1*xi*xi*xi*xi*xi*xi+0.10125E1*xi*xi*xi*xi*xi-0.5625*xi*xi*xi*xi-0.5625*xi*xi*xi+0.5E-1*xi*xi+0.5E-1*xi;

      double eta0= 0.30375E2*eta*eta*eta*eta-0.2025E2*eta*eta*eta-0.675E1*eta*eta+0.3375E1*eta+0.1;
      double eta1= -0.18225E3*eta*eta*eta*eta+0.81E2*eta*eta*eta+0.81E2*eta*eta-0.27E2*eta-0.135E1;
      double eta2= 0.455625E3*eta*eta*eta*eta-0.10125E3*eta*eta*eta-0.26325E3*eta*eta+0.43875E2*eta+0.135E2;
      double eta3= -0.6075E3*eta*eta*eta*eta+0.378E3*eta*eta-0.245E2;
      double eta4= 0.455625E3*eta*eta*eta*eta+0.10125E3*eta*eta*eta-0.26325E3*eta*eta-0.43875E2*eta+0.135E2;
      double eta5= -0.18225E3*eta*eta*eta*eta-0.81E2*eta*eta*eta+0.81E2*eta*eta+0.27E2*eta-0.135E1;
      double eta6= 0.30375E2*eta*eta*eta*eta+0.2025E2*eta*eta*eta-0.675E1*eta*eta-0.3375E1*eta+0.1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi6*eta0;
      values[7] = xi0*eta1;
      values[8] = xi1*eta1;
      values[9] = xi2*eta1;
      values[10] = xi3*eta1;
      values[11] = xi4*eta1;
      values[12] = xi5*eta1;
      values[13] = xi6*eta1;
      values[14] = xi0*eta2;
      values[15] = xi1*eta2;
      values[16] = xi2*eta2;
      values[17] = xi3*eta2;
      values[18] = xi4*eta2;
      values[19] = xi5*eta2;
      values[20] = xi6*eta2;
      values[21] = xi0*eta3;
      values[22] = xi1*eta3;
      values[23] = xi2*eta3;
      values[24] = xi3*eta3;
      values[25] = xi4*eta3;
      values[26] = xi5*eta3;
      values[27] = xi6*eta3;
      values[28] = xi0*eta4;
      values[29] = xi1*eta4;
      values[30] = xi2*eta4;
      values[31] = xi3*eta4;
      values[32] = xi4*eta4;
      values[33] = xi5*eta4;
      values[34] = xi6*eta4;
      values[35] = xi0*eta5;
      values[36] = xi1*eta5;
      values[37] = xi2*eta5;
      values[38] = xi3*eta5;
      values[39] = xi4*eta5;
      values[40] = xi5*eta5;
      values[41] = xi6*eta5;
      values[42] = xi0*eta6;
      values[43] = xi1*eta6;
      values[44] = xi2*eta6;
      values[45] = xi3*eta6;
      values[46] = xi4*eta6;
      values[47] = xi5*eta6;
      values[48] = xi6*eta6;
}


// ***********************************************************************

TBaseFunct2D *BF_C_Q_Q6_2D_Obj = new TBaseFunct2D
        (49, BF_C_Q_Q6_2D, BFUnitSquare, 
         C_Q_Q6_2D_Funct, C_Q_Q6_2D_DeriveXi,
         C_Q_Q6_2D_DeriveEta, C_Q_Q6_2D_DeriveXiXi,
         C_Q_Q6_2D_DeriveXiEta, C_Q_Q6_2D_DeriveEtaEta, 6, 6,
         0, NULL);
