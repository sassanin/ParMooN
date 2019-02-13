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
// Q5 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_Q5_2D_Funct(double xi, double eta, double *values)
{  

     double xi0= -0.8138020833333333*xi*xi*xi*xi*xi+0.8138020833333333*xi*xi*xi*xi+0.3255208333333333*xi*xi*xi-0.3255208333333333*xi*xi-0.1171875E-1*xi+0.1171875E-1;
     double xi1= 0.4069010416666667E1*xi*xi*xi*xi*xi-0.244140625E1*xi*xi*xi*xi-0.4231770833333333E1*xi*xi*xi+0.25390625E1*xi*xi+0.1627604166666667*xi-0.9765625E-1;
     double xi2= -0.8138020833333333E1*xi*xi*xi*xi*xi+0.1627604166666667E1*xi*xi*xi*xi+0.1106770833333333E2*xi*xi*xi-0.2213541666666667E1*xi*xi-0.29296875E1*xi+0.5859375;
     double xi3= 0.8138020833333333E1*xi*xi*xi*xi*xi+0.1627604166666667E1*xi*xi*xi*xi-0.1106770833333333E2*xi*xi*xi-0.2213541666666667E1*xi*xi+0.29296875E1*xi+0.5859375;
     double xi4= -0.4069010416666667E1*xi*xi*xi*xi*xi-0.244140625E1*xi*xi*xi*xi+0.4231770833333333E1*xi*xi*xi+0.25390625E1*xi*xi-0.1627604166666667*xi-0.9765625E-1;
     double xi5= 0.8138020833333333*xi*xi*xi*xi*xi+0.8138020833333333*xi*xi*xi*xi-0.3255208333333333*xi*xi*xi-0.3255208333333333*xi*xi+0.1171875E-1*xi+0.1171875E-1;

     double eta0= -0.8138020833333333*eta*eta*eta*eta*eta+0.8138020833333333*eta*eta*eta*eta+0.3255208333333333*eta*eta*eta-0.3255208333333333*eta*eta-0.1171875E-1*eta+0.1171875E-1;
     double eta1= 0.4069010416666667E1*eta*eta*eta*eta*eta-0.244140625E1*eta*eta*eta*eta-0.4231770833333333E1*eta*eta*eta+0.25390625E1*eta*eta+0.1627604166666667*eta-0.9765625E-1;
     double eta2= -0.8138020833333333E1*eta*eta*eta*eta*eta+0.1627604166666667E1*eta*eta*eta*eta+0.1106770833333333E2*eta*eta*eta-0.2213541666666667E1*eta*eta-0.29296875E1*eta+0.5859375;
     double eta3= 0.8138020833333333E1*eta*eta*eta*eta*eta+0.1627604166666667E1*eta*eta*eta*eta-0.1106770833333333E2*eta*eta*eta-0.2213541666666667E1*eta*eta+0.29296875E1*eta+0.5859375;
     double eta4= -0.4069010416666667E1*eta*eta*eta*eta*eta-0.244140625E1*eta*eta*eta*eta+0.4231770833333333E1*eta*eta*eta+0.25390625E1*eta*eta-0.1627604166666667*eta-0.9765625E-1;
     double eta5= 0.8138020833333333*eta*eta*eta*eta*eta+0.8138020833333333*eta*eta*eta*eta-0.3255208333333333*eta*eta*eta-0.3255208333333333*eta*eta+0.1171875E-1*eta+0.1171875E-1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi0*eta1;
      values[7] = xi1*eta1;
      values[8] = xi2*eta1;
      values[9] = xi3*eta1;
      values[10] = xi4*eta1;
      values[11] = xi5*eta1;
      values[12] = xi0*eta2;
      values[13] = xi1*eta2;
      values[14] = xi2*eta2;
      values[15] = xi3*eta2;
      values[16] = xi4*eta2;
      values[17] = xi5*eta2;
      values[18] = xi0*eta3;
      values[19] = xi1*eta3;
      values[20] = xi2*eta3;
      values[21] = xi3*eta3;
      values[22] = xi4*eta3;
      values[23] = xi5*eta3;
      values[24] = xi0*eta4;
      values[25] = xi1*eta4;
      values[26] = xi2*eta4;
      values[27] = xi3*eta4;
      values[28] = xi4*eta4;
      values[29] = xi5*eta4;
      values[30] = xi0*eta5;
      values[31] = xi1*eta5;
      values[32] = xi2*eta5;
      values[33] = xi3*eta5;
      values[34] = xi4*eta5;
      values[35] = xi5*eta5;
}


// values of the derivatives in xi direction
static void C_Q_Q5_2D_DeriveXi(double xi, double eta, double *values)
{

     double xi0= -0.4069010416666667E1*xi*xi*xi*xi+0.3255208333333333E1*xi*xi*xi+0.9765625*xi*xi-0.6510416666666667*xi-0.1171875E-1;
     double xi1= 0.2034505208333333E2*xi*xi*xi*xi-0.9765625E1*xi*xi*xi-0.126953125E2*xi*xi+0.5078125E1*xi+0.1627604166666667;
     double xi2= -0.4069010416666667E2*xi*xi*xi*xi+0.6510416666666667E1*xi*xi*xi+0.33203125E2*xi*xi-0.4427083333333333E1*xi-0.29296875E1;
     double xi3= 0.4069010416666667E2*xi*xi*xi*xi+0.6510416666666667E1*xi*xi*xi-0.33203125E2*xi*xi-0.4427083333333333E1*xi+0.29296875E1;
     double xi4= -0.2034505208333333E2*xi*xi*xi*xi-0.9765625E1*xi*xi*xi+0.126953125E2*xi*xi+0.5078125E1*xi-0.1627604166666667;
     double xi5= 0.4069010416666667E1*xi*xi*xi*xi+0.3255208333333333E1*xi*xi*xi-0.9765625*xi*xi-0.6510416666666667*xi+0.1171875E-1;

     double eta0= -0.8138020833333333*eta*eta*eta*eta*eta+0.8138020833333333*eta*eta*eta*eta+0.3255208333333333*eta*eta*eta-0.3255208333333333*eta*eta-0.1171875E-1*eta+0.1171875E-1;
     double eta1= 0.4069010416666667E1*eta*eta*eta*eta*eta-0.244140625E1*eta*eta*eta*eta-0.4231770833333333E1*eta*eta*eta+0.25390625E1*eta*eta+0.1627604166666667*eta-0.9765625E-1;
     double eta2= -0.8138020833333333E1*eta*eta*eta*eta*eta+0.1627604166666667E1*eta*eta*eta*eta+0.1106770833333333E2*eta*eta*eta-0.2213541666666667E1*eta*eta-0.29296875E1*eta+0.5859375;
     double eta3= 0.8138020833333333E1*eta*eta*eta*eta*eta+0.1627604166666667E1*eta*eta*eta*eta-0.1106770833333333E2*eta*eta*eta-0.2213541666666667E1*eta*eta+0.29296875E1*eta+0.5859375;
     double eta4= -0.4069010416666667E1*eta*eta*eta*eta*eta-0.244140625E1*eta*eta*eta*eta+0.4231770833333333E1*eta*eta*eta+0.25390625E1*eta*eta-0.1627604166666667*eta-0.9765625E-1;
     double eta5= 0.8138020833333333*eta*eta*eta*eta*eta+0.8138020833333333*eta*eta*eta*eta-0.3255208333333333*eta*eta*eta-0.3255208333333333*eta*eta+0.1171875E-1*eta+0.1171875E-1;


      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi0*eta1;
      values[7] = xi1*eta1;
      values[8] = xi2*eta1;
      values[9] = xi3*eta1;
      values[10] = xi4*eta1;
      values[11] = xi5*eta1;
      values[12] = xi0*eta2;
      values[13] = xi1*eta2;
      values[14] = xi2*eta2;
      values[15] = xi3*eta2;
      values[16] = xi4*eta2;
      values[17] = xi5*eta2;
      values[18] = xi0*eta3;
      values[19] = xi1*eta3;
      values[20] = xi2*eta3;
      values[21] = xi3*eta3;
      values[22] = xi4*eta3;
      values[23] = xi5*eta3;
      values[24] = xi0*eta4;
      values[25] = xi1*eta4;
      values[26] = xi2*eta4;
      values[27] = xi3*eta4;
      values[28] = xi4*eta4;
      values[29] = xi5*eta4;
      values[30] = xi0*eta5;
      values[31] = xi1*eta5;
      values[32] = xi2*eta5;
      values[33] = xi3*eta5;
      values[34] = xi4*eta5;
      values[35] = xi5*eta5;
}

// values of the derivatives in eta direction
static void C_Q_Q5_2D_DeriveEta(double xi, double eta, double *values)
{
  
     double xi0= -0.8138020833333333*xi*xi*xi*xi*xi+0.8138020833333333*xi*xi*xi*xi+0.3255208333333333*xi*xi*xi-0.3255208333333333*xi*xi-0.1171875E-1*xi+0.1171875E-1;
     double xi1= 0.4069010416666667E1*xi*xi*xi*xi*xi-0.244140625E1*xi*xi*xi*xi-0.4231770833333333E1*xi*xi*xi+0.25390625E1*xi*xi+0.1627604166666667*xi-0.9765625E-1;
     double xi2= -0.8138020833333333E1*xi*xi*xi*xi*xi+0.1627604166666667E1*xi*xi*xi*xi+0.1106770833333333E2*xi*xi*xi-0.2213541666666667E1*xi*xi-0.29296875E1*xi+0.5859375;
     double xi3= 0.8138020833333333E1*xi*xi*xi*xi*xi+0.1627604166666667E1*xi*xi*xi*xi-0.1106770833333333E2*xi*xi*xi-0.2213541666666667E1*xi*xi+0.29296875E1*xi+0.5859375;
     double xi4= -0.4069010416666667E1*xi*xi*xi*xi*xi-0.244140625E1*xi*xi*xi*xi+0.4231770833333333E1*xi*xi*xi+0.25390625E1*xi*xi-0.1627604166666667*xi-0.9765625E-1;
     double xi5= 0.8138020833333333*xi*xi*xi*xi*xi+0.8138020833333333*xi*xi*xi*xi-0.3255208333333333*xi*xi*xi-0.3255208333333333*xi*xi+0.1171875E-1*xi+0.1171875E-1;

     double eta0= -0.4069010416666667E1*eta*eta*eta*eta+0.3255208333333333E1*eta*eta*eta+0.9765625*eta*eta-0.6510416666666667*eta-0.1171875E-1;
     double eta1= 0.2034505208333333E2*eta*eta*eta*eta-0.9765625E1*eta*eta*eta-0.126953125E2*eta*eta+0.5078125E1*eta+0.1627604166666667;
     double eta2= -0.4069010416666667E2*eta*eta*eta*eta+0.6510416666666667E1*eta*eta*eta+0.33203125E2*eta*eta-0.4427083333333333E1*eta-0.29296875E1;
     double eta3= 0.4069010416666667E2*eta*eta*eta*eta+0.6510416666666667E1*eta*eta*eta-0.33203125E2*eta*eta-0.4427083333333333E1*eta+0.29296875E1;
     double eta4= -0.2034505208333333E2*eta*eta*eta*eta-0.9765625E1*eta*eta*eta+0.126953125E2*eta*eta+0.5078125E1*eta-0.1627604166666667;
     double eta5= 0.4069010416666667E1*eta*eta*eta*eta+0.3255208333333333E1*eta*eta*eta-0.9765625*eta*eta-0.6510416666666667*eta+0.1171875E-1;
    

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi0*eta1;
      values[7] = xi1*eta1;
      values[8] = xi2*eta1;
      values[9] = xi3*eta1;
      values[10] = xi4*eta1;
      values[11] = xi5*eta1;
      values[12] = xi0*eta2;
      values[13] = xi1*eta2;
      values[14] = xi2*eta2;
      values[15] = xi3*eta2;
      values[16] = xi4*eta2;
      values[17] = xi5*eta2;
      values[18] = xi0*eta3;
      values[19] = xi1*eta3;
      values[20] = xi2*eta3;
      values[21] = xi3*eta3;
      values[22] = xi4*eta3;
      values[23] = xi5*eta3;
      values[24] = xi0*eta4;
      values[25] = xi1*eta4;
      values[26] = xi2*eta4;
      values[27] = xi3*eta4;
      values[28] = xi4*eta4;
      values[29] = xi5*eta4;
      values[30] = xi0*eta5;
      values[31] = xi1*eta5;
      values[32] = xi2*eta5;
      values[33] = xi3*eta5;
      values[34] = xi4*eta5;
      values[35] = xi5*eta5;
}
// values of the derivatives in xi-xi  direction
static void C_Q_Q5_2D_DeriveXiXi(double xi, double eta, double *values)
{
     
     double xi0= -0.1627604166666667E2*xi*xi*xi+0.9765625E1*xi*xi+0.1953125E1*xi-0.6510416666666667;
     double xi1= 0.8138020833333333E2*xi*xi*xi-0.29296875E2*xi*xi-0.25390625E2*xi+0.5078125E1;
     double xi2= -0.1627604166666667E3*xi*xi*xi+0.1953125E2*xi*xi+0.6640625E2*xi-0.4427083333333333E1;
     double xi3= 0.1627604166666667E3*xi*xi*xi+0.1953125E2*xi*xi-0.6640625E2*xi-0.4427083333333333E1;
     double xi4= -0.8138020833333333E2*xi*xi*xi-0.29296875E2*xi*xi+0.25390625E2*xi+0.5078125E1;
     double xi5= 0.1627604166666667E2*xi*xi*xi+0.9765625E1*xi*xi-0.1953125E1*xi-0.6510416666666667;

     double eta0= -0.8138020833333333*eta*eta*eta*eta*eta+0.8138020833333333*eta*eta*eta*eta+0.3255208333333333*eta*eta*eta-0.3255208333333333*eta*eta-0.1171875E-1*eta+0.1171875E-1;
     double eta1= 0.4069010416666667E1*eta*eta*eta*eta*eta-0.244140625E1*eta*eta*eta*eta-0.4231770833333333E1*eta*eta*eta+0.25390625E1*eta*eta+0.1627604166666667*eta-0.9765625E-1;
     double eta2= -0.8138020833333333E1*eta*eta*eta*eta*eta+0.1627604166666667E1*eta*eta*eta*eta+0.1106770833333333E2*eta*eta*eta-0.2213541666666667E1*eta*eta-0.29296875E1*eta+0.5859375;
     double eta3= 0.8138020833333333E1*eta*eta*eta*eta*eta+0.1627604166666667E1*eta*eta*eta*eta-0.1106770833333333E2*eta*eta*eta-0.2213541666666667E1*eta*eta+0.29296875E1*eta+0.5859375;
     double eta4= -0.4069010416666667E1*eta*eta*eta*eta*eta-0.244140625E1*eta*eta*eta*eta+0.4231770833333333E1*eta*eta*eta+0.25390625E1*eta*eta-0.1627604166666667*eta-0.9765625E-1;
     double eta5= 0.8138020833333333*eta*eta*eta*eta*eta+0.8138020833333333*eta*eta*eta*eta-0.3255208333333333*eta*eta*eta-0.3255208333333333*eta*eta+0.1171875E-1*eta+0.1171875E-1;
    

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi0*eta1;
      values[7] = xi1*eta1;
      values[8] = xi2*eta1;
      values[9] = xi3*eta1;
      values[10] = xi4*eta1;
      values[11] = xi5*eta1;
      values[12] = xi0*eta2;
      values[13] = xi1*eta2;
      values[14] = xi2*eta2;
      values[15] = xi3*eta2;
      values[16] = xi4*eta2;
      values[17] = xi5*eta2;
      values[18] = xi0*eta3;
      values[19] = xi1*eta3;
      values[20] = xi2*eta3;
      values[21] = xi3*eta3;
      values[22] = xi4*eta3;
      values[23] = xi5*eta3;
      values[24] = xi0*eta4;
      values[25] = xi1*eta4;
      values[26] = xi2*eta4;
      values[27] = xi3*eta4;
      values[28] = xi4*eta4;
      values[29] = xi5*eta4;
      values[30] = xi0*eta5;
      values[31] = xi1*eta5;
      values[32] = xi2*eta5;
      values[33] = xi3*eta5;
      values[34] = xi4*eta5;
      values[35] = xi5*eta5;
}
// values of the derivatives in xi-eta direction
static void C_Q_Q5_2D_DeriveXiEta(double xi, double eta, double *values)
{

     double xi0= -0.4069010416666667E1*xi*xi*xi*xi+0.3255208333333333E1*xi*xi*xi+0.9765625*xi*xi-0.6510416666666667*xi-0.1171875E-1;
     double xi1= 0.2034505208333333E2*xi*xi*xi*xi-0.9765625E1*xi*xi*xi-0.126953125E2*xi*xi+0.5078125E1*xi+0.1627604166666667;
     double xi2= -0.4069010416666667E2*xi*xi*xi*xi+0.6510416666666667E1*xi*xi*xi+0.33203125E2*xi*xi-0.4427083333333333E1*xi-0.29296875E1;
     double xi3= 0.4069010416666667E2*xi*xi*xi*xi+0.6510416666666667E1*xi*xi*xi-0.33203125E2*xi*xi-0.4427083333333333E1*xi+0.29296875E1;
     double xi4= -0.2034505208333333E2*xi*xi*xi*xi-0.9765625E1*xi*xi*xi+0.126953125E2*xi*xi+0.5078125E1*xi-0.1627604166666667;
     double xi5= 0.4069010416666667E1*xi*xi*xi*xi+0.3255208333333333E1*xi*xi*xi-0.9765625*xi*xi-0.6510416666666667*xi+0.1171875E-1;

     double eta0= -0.4069010416666667E1*eta*eta*eta*eta+0.3255208333333333E1*eta*eta*eta+0.9765625*eta*eta-0.6510416666666667*eta-0.1171875E-1;
     double eta1= 0.2034505208333333E2*eta*eta*eta*eta-0.9765625E1*eta*eta*eta-0.126953125E2*eta*eta+0.5078125E1*eta+0.1627604166666667;
     double eta2= -0.4069010416666667E2*eta*eta*eta*eta+0.6510416666666667E1*eta*eta*eta+0.33203125E2*eta*eta-0.4427083333333333E1*eta-0.29296875E1;
     double eta3= 0.4069010416666667E2*eta*eta*eta*eta+0.6510416666666667E1*eta*eta*eta-0.33203125E2*eta*eta-0.4427083333333333E1*eta+0.29296875E1;
     double eta4= -0.2034505208333333E2*eta*eta*eta*eta-0.9765625E1*eta*eta*eta+0.126953125E2*eta*eta+0.5078125E1*eta-0.1627604166666667;
     double eta5= 0.4069010416666667E1*eta*eta*eta*eta+0.3255208333333333E1*eta*eta*eta-0.9765625*eta*eta-0.6510416666666667*eta+0.1171875E-1;
    

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi0*eta1;
      values[7] = xi1*eta1;
      values[8] = xi2*eta1;
      values[9] = xi3*eta1;
      values[10] = xi4*eta1;
      values[11] = xi5*eta1;
      values[12] = xi0*eta2;
      values[13] = xi1*eta2;
      values[14] = xi2*eta2;
      values[15] = xi3*eta2;
      values[16] = xi4*eta2;
      values[17] = xi5*eta2;
      values[18] = xi0*eta3;
      values[19] = xi1*eta3;
      values[20] = xi2*eta3;
      values[21] = xi3*eta3;
      values[22] = xi4*eta3;
      values[23] = xi5*eta3;
      values[24] = xi0*eta4;
      values[25] = xi1*eta4;
      values[26] = xi2*eta4;
      values[27] = xi3*eta4;
      values[28] = xi4*eta4;
      values[29] = xi5*eta4;
      values[30] = xi0*eta5;
      values[31] = xi1*eta5;
      values[32] = xi2*eta5;
      values[33] = xi3*eta5;
      values[34] = xi4*eta5;
      values[35] = xi5*eta5;
}

// values of the derivatives in eta-eta direction
static void C_Q_Q5_2D_DeriveEtaEta(double xi, double eta, double *values)
{
     double xi0= -0.8138020833333333*xi*xi*xi*xi*xi+0.8138020833333333*xi*xi*xi*xi+0.3255208333333333*xi*xi*xi-0.3255208333333333*xi*xi-0.1171875E-1*xi+0.1171875E-1;
     double xi1= 0.4069010416666667E1*xi*xi*xi*xi*xi-0.244140625E1*xi*xi*xi*xi-0.4231770833333333E1*xi*xi*xi+0.25390625E1*xi*xi+0.1627604166666667*xi-0.9765625E-1;
     double xi2= -0.8138020833333333E1*xi*xi*xi*xi*xi+0.1627604166666667E1*xi*xi*xi*xi+0.1106770833333333E2*xi*xi*xi-0.2213541666666667E1*xi*xi-0.29296875E1*xi+0.5859375;
     double xi3= 0.8138020833333333E1*xi*xi*xi*xi*xi+0.1627604166666667E1*xi*xi*xi*xi-0.1106770833333333E2*xi*xi*xi-0.2213541666666667E1*xi*xi+0.29296875E1*xi+0.5859375;
     double xi4= -0.4069010416666667E1*xi*xi*xi*xi*xi-0.244140625E1*xi*xi*xi*xi+0.4231770833333333E1*xi*xi*xi+0.25390625E1*xi*xi-0.1627604166666667*xi-0.9765625E-1;
     double xi5= 0.8138020833333333*xi*xi*xi*xi*xi+0.8138020833333333*xi*xi*xi*xi-0.3255208333333333*xi*xi*xi-0.3255208333333333*xi*xi+0.1171875E-1*xi+0.1171875E-1;

     double eta0= -0.1627604166666667E2*eta*eta*eta+0.9765625E1*eta*eta+0.1953125E1*eta-0.6510416666666667;
     double eta1= 0.8138020833333333E2*eta*eta*eta-0.29296875E2*eta*eta-0.25390625E2*eta+0.5078125E1;
     double eta2= -0.1627604166666667E3*eta*eta*eta+0.1953125E2*eta*eta+0.6640625E2*eta-0.4427083333333333E1;
     double eta3= 0.1627604166666667E3*eta*eta*eta+0.1953125E2*eta*eta-0.6640625E2*eta-0.4427083333333333E1;
     double eta4= -0.8138020833333333E2*eta*eta*eta-0.29296875E2*eta*eta+0.25390625E2*eta+0.5078125E1;
     double eta5= 0.1627604166666667E2*eta*eta*eta+0.9765625E1*eta*eta-0.1953125E1*eta-0.6510416666666667;
    

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi5*eta0;
      values[6] = xi0*eta1;
      values[7] = xi1*eta1;
      values[8] = xi2*eta1;
      values[9] = xi3*eta1;
      values[10] = xi4*eta1;
      values[11] = xi5*eta1;
      values[12] = xi0*eta2;
      values[13] = xi1*eta2;
      values[14] = xi2*eta2;
      values[15] = xi3*eta2;
      values[16] = xi4*eta2;
      values[17] = xi5*eta2;
      values[18] = xi0*eta3;
      values[19] = xi1*eta3;
      values[20] = xi2*eta3;
      values[21] = xi3*eta3;
      values[22] = xi4*eta3;
      values[23] = xi5*eta3;
      values[24] = xi0*eta4;
      values[25] = xi1*eta4;
      values[26] = xi2*eta4;
      values[27] = xi3*eta4;
      values[28] = xi4*eta4;
      values[29] = xi5*eta4;
      values[30] = xi0*eta5;
      values[31] = xi1*eta5;
      values[32] = xi2*eta5;
      values[33] = xi3*eta5;
      values[34] = xi4*eta5;
      values[35] = xi5*eta5;
}


// ***********************************************************************

TBaseFunct2D *BF_C_Q_Q5_2D_Obj = new TBaseFunct2D
        (36, BF_C_Q_Q5_2D, BFUnitSquare, 
         C_Q_Q5_2D_Funct, C_Q_Q5_2D_DeriveXi,
         C_Q_Q5_2D_DeriveEta, C_Q_Q5_2D_DeriveXiXi,
         C_Q_Q5_2D_DeriveXiEta, C_Q_Q5_2D_DeriveEtaEta, 5, 5,
         0, NULL);
