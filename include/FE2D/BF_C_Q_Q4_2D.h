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
   
// ********************************************************************
// Q4 element, conforming, 2D
// ********************************************************************

// base function values
static void C_Q_Q4_2D_Funct(double xi, double eta, double *values)
{     
      
      double xi0 = 0.666666666666666666667*xi*xi*xi*xi-0.666666666666666666667*xi*xi*xi - 0.166666666666666666667*xi*xi+0.166666666666666666667*xi;
      double xi1 = -2.66666666666666666667*xi*xi*xi*xi+1.33333333333333333333*xi*xi*xi+2.66666666666666666667*xi*xi-1.33333333333333333333*xi;
      double xi2 = 4.0*xi*xi*xi*xi-5.0*xi*xi+1.0;
      double xi3 = -2.66666666666666666667*xi*xi*xi*xi-1.33333333333333333333*xi*xi*xi+2.666666666666667*xi*xi+1.33333333333333333333*xi;
      double xi4 = 0.666666666666666666667*xi*xi*xi*xi+0.666666666666666666667*xi*xi*xi-0.166666666666666666667*xi*xi-0.166666666666666666667*xi;
      double eta0 = 0.666666666666666666667*eta*eta*eta*eta-0.666666666666666666667*eta*eta*eta-0.166666666666666666667*eta*eta+0.166666666666666666667*eta;
      double eta1 = -2.66666666666666666667*eta*eta*eta*eta+1.33333333333333333333*eta*eta*eta+2.66666666666666666667*eta*eta-1.33333333333333333333*eta;
      double eta2 = 4.0*eta*eta*eta*eta-5.0*eta*eta+1.0;
      double eta3 = -2.66666666666666666667*eta*eta*eta*eta-1.33333333333333333333*eta*eta*eta+2.66666666666666666667*eta*eta+1.33333333333333333333*eta;
      double eta4 = 0.666666666666666666667*eta*eta*eta*eta+0.666666666666666666667*eta*eta*eta-0.166666666666666666667*eta*eta-0.166666666666666666667*eta;

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi0*eta1;
      values[6] = xi1*eta1;
      values[7] = xi2*eta1;
      values[8] = xi3*eta1;
      values[9] = xi4*eta1;
      values[10] = xi0*eta2;
      values[11] = xi1*eta2;
      values[12] = xi2*eta2;
      values[13] = xi3*eta2;
      values[14] = xi4*eta2;
      values[15] = xi0*eta3;
      values[16] = xi1*eta3;
      values[17] = xi2*eta3;
      values[18] = xi3*eta3;
      values[19] = xi4*eta3;
      values[20] = xi0*eta4;
      values[21] = xi1*eta4;
      values[22] = xi2*eta4;
      values[23] = xi3*eta4;
      values[24] = xi4*eta4;
}


// values of the derivatives in xi direction
static void C_Q_Q4_2D_DeriveXi(double xi, double eta, double *values)
{

      double xi0 = 2.666666666666666666667*xi*xi*xi-2*xi*xi-0.333333333333333333333*xi+0.1666666666666666666667;
      double xi1 = -10.66666666666666666667*xi*xi*xi+4*xi*xi+5.33333333333333333333*xi-1.33333333333333333333;
      double xi2 = 16.0*xi*xi*xi-10.0*xi;
      double xi3 = -10.66666666666666666667*xi*xi*xi-4*xi*xi+5.33333333333333333333*xi+1.33333333333333333333;
      double xi4  = 2.666666666666666666667*xi*xi*xi+2*xi*xi-0.333333333333333333333*xi-0.1666666666666666666667;
      double eta0 = 0.666666666666666666667*eta*eta*eta*eta-0.6666666666666667*eta*eta*eta-0.166666666666666666667*eta*eta+0.166666666666666666667*eta;
      double eta1 = -2.66666666666666666667*eta*eta*eta*eta+1.33333333333333333333*eta*eta*eta+2.66666666666666666667*eta*eta-1.33333333333333333333*eta;
      double eta2 = 4.0*eta*eta*eta*eta-5.0*eta*eta+1.0;
      double eta3 = -2.66666666666666666667*eta*eta*eta*eta-1.33333333333333333333*eta*eta*eta+2.66666666666666666667*eta*eta+1.33333333333333333333*eta;
      double eta4 = 0.666666666666666666667*eta*eta*eta*eta+0.666666666666666666667*eta*eta*eta-0.166666666666666666667*eta*eta-0.166666666666666666667*eta;

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi0*eta1;
      values[6] = xi1*eta1;
      values[7] = xi2*eta1;
      values[8] = xi3*eta1;
      values[9] = xi4*eta1;
      values[10] = xi0*eta2;
      values[11] = xi1*eta2;
      values[12] = xi2*eta2;
      values[13] = xi3*eta2;
      values[14] = xi4*eta2;
      values[15] = xi0*eta3;
      values[16] = xi1*eta3;
      values[17] = xi2*eta3;
      values[18] = xi3*eta3;
      values[19] = xi4*eta3;
      values[20] = xi0*eta4;
      values[21] = xi1*eta4;
      values[22] = xi2*eta4;
      values[23] = xi3*eta4;
      values[24] = xi4*eta4;
}

// values of the derivatives in eta direction
static void C_Q_Q4_2D_DeriveEta(double xi, double eta, double *values)
{

      double xi0 = 0.666666666666666666667*xi*xi*xi*xi-0.666666666666666666667*xi*xi*xi-0.166666666666666666667*xi*xi+0.166666666666666666667*xi;
      double xi1 = -2.66666666666666666667*xi*xi*xi*xi+1.33333333333333333333*xi*xi*xi+2.666666666666667*xi*xi-1.33333333333333333333*xi;
      double xi2 = 4.0*xi*xi*xi*xi-5.0*xi*xi+1.0;
      double xi3 = -2.66666666666666666667*xi*xi*xi*xi-1.33333333333333333333*xi*xi*xi+2.66666666666666666667*xi*xi+1.33333333333333333333*xi;
      double xi4 = 0.666666666666666666667*xi*xi*xi*xi+0.6666666666666666666667*xi*xi*xi-0.166666666666666666667*xi*xi-0.1666666666666666666667*xi;
      double eta0 = 2.666666666666666666667*eta*eta*eta-2*eta*eta-0.333333333333333333333*eta+0.1666666666666666666667;
      double eta1 = -10.6666666666666666667*eta*eta*eta+4*eta*eta+5.33333333333333333333*eta-1.33333333333333333333;
      double eta2 = 16.0*eta*eta*eta-10.0*eta;
      double eta3 = -10.66666666666666666667*eta*eta*eta-4*eta*eta+5.33333333333333333333*eta+1.33333333333333333333;
      double eta4 = 2.66666666666666666667*eta*eta*eta+2*eta*eta-0.333333333333333333333*eta-0.166666666666666666667;

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi0*eta1;
      values[6] = xi1*eta1;
      values[7] = xi2*eta1;
      values[8] = xi3*eta1;
      values[9] = xi4*eta1;
      values[10] = xi0*eta2;
      values[11] = xi1*eta2;
      values[12] = xi2*eta2;
      values[13] = xi3*eta2;
      values[14] = xi4*eta2;
      values[15] = xi0*eta3;
      values[16] = xi1*eta3;
      values[17] = xi2*eta3;
      values[18] = xi3*eta3;
      values[19] = xi4*eta3;
      values[20] = xi0*eta4;
      values[21] = xi1*eta4;
      values[22] = xi2*eta4;
      values[23] = xi3*eta4;
      values[24] = xi4*eta4;
}
// values of the derivatives in xi-xi  direction
static void C_Q_Q4_2D_DeriveXiXi(double xi, double eta, double *values)
{

      double xi0 = 8*xi*xi-4*xi-0.333333333333333333333;
      double xi1 = -32*xi*xi+8*xi+5.33333333333333333333;
      double xi2 = 48.0*xi*xi-10.0;
      double xi3 = -32*xi*xi-8*xi+5.33333333333333333333;
      double xi4 = 8*xi*xi+4*xi-0.333333333333333333333;
      double eta0 = 0.6666666666666666666667*eta*eta*eta*eta-0.666666666666666666667*eta*eta*eta-0.166666666666666666667*eta*eta+0.166666666666666666667*eta;
      double eta1 = -2.66666666666666666667*eta*eta*eta*eta+1.33333333333333333333*eta*eta*eta+2.66666666666666666667*eta*eta-1.33333333333333333333*eta;
      double eta2 = 4.0*eta*eta*eta*eta-5.0*eta*eta+1.0;
      double eta3 = -2.66666666666666666667*eta*eta*eta*eta-1.33333333333333333333*eta*eta*eta+2.66666666666666666667*eta*eta+1.33333333333333333333*eta;
      double eta4 = 0.666666666666666666667*eta*eta*eta*eta+0.666666666666666666667*eta*eta*eta-0.166666666666666666667*eta*eta-0.166666666666666666667*eta;

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi0*eta1;
      values[6] = xi1*eta1;
      values[7] = xi2*eta1;
      values[8] = xi3*eta1;
      values[9] = xi4*eta1;
      values[10] = xi0*eta2;
      values[11] = xi1*eta2;
      values[12] = xi2*eta2;
      values[13] = xi3*eta2;
      values[14] = xi4*eta2;
      values[15] = xi0*eta3;
      values[16] = xi1*eta3;
      values[17] = xi2*eta3;
      values[18] = xi3*eta3;
      values[19] = xi4*eta3;
      values[20] = xi0*eta4;
      values[21] = xi1*eta4;
      values[22] = xi2*eta4;
      values[23] = xi3*eta4;
      values[24] = xi4*eta4;
}
// values of the derivatives in xi-eta direction
static void C_Q_Q4_2D_DeriveXiEta(double xi, double eta, double *values)
{

      double xi0 = 2.666666666666666666667*xi*xi*xi-2*xi*xi-0.333333333333333333333*xi+0.166666666666666666667;
      double xi1 = -10.6666666666666666667*xi*xi*xi+4*xi*xi+5.33333333333333333333*xi-1.33333333333333333333;
      double xi2 = 16.0*xi*xi*xi-10.0*xi;
      double xi3 = -10.6666666666666666667*xi*xi*xi-4*xi*xi+5.33333333333333333333*xi+1.33333333333333333333;
      double xi4 = 2.66666666666666666667*xi*xi*xi+2*xi*xi-0.333333333333333333333*xi-0.166666666666666666667;
      double eta0 = 2.66666666666666666667*eta*eta*eta-2*eta*eta-0.333333333333333333333*eta+0.166666666666666666667;
      double eta1 = -10.6666666666666666667*eta*eta*eta+4*eta*eta+5.33333333333333333333*eta-1.33333333333333333333;
      double eta2 = 16.0*eta*eta*eta-10.0*eta;
      double eta3 = -10.6666666666666666667*eta*eta*eta-4*eta*eta+5.33333333333333333333*eta+1.33333333333333333333;
      double eta4 = 2.66666666666666666667*eta*eta*eta+2*eta*eta-0.333333333333333333333*eta-0.166666666666666666667;

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi0*eta1;
      values[6] = xi1*eta1;
      values[7] = xi2*eta1;
      values[8] = xi3*eta1;
      values[9] = xi4*eta1;
      values[10] = xi0*eta2;
      values[11] = xi1*eta2;
      values[12] = xi2*eta2;
      values[13] = xi3*eta2;
      values[14] = xi4*eta2;
      values[15] = xi0*eta3;
      values[16] = xi1*eta3;
      values[17] = xi2*eta3;
      values[18] = xi3*eta3;
      values[19] = xi4*eta3;
      values[20] = xi0*eta4;
      values[21] = xi1*eta4;
      values[22] = xi2*eta4;
      values[23] = xi3*eta4;
      values[24] = xi4*eta4;
}

// values of the derivatives in eta-eta direction
static void C_Q_Q4_2D_DeriveEtaEta(double xi, double eta, double *values)
{

      double xi0 = 0.666666666666666666667*xi*xi*xi*xi-0.666666666666666666667*xi*xi*xi-0.166666666666666666667*xi*xi+0.166666666666666666667*xi;
      double xi1 = -2.66666666666666666667*xi*xi*xi*xi+1.33333333333333333333*xi*xi*xi+2.66666666666666666667*xi*xi-1.33333333333333333333*xi;
      double xi2 = 4.0*xi*xi*xi*xi-5.0*xi*xi+1.0;
      double xi3 = -2.666666666666666666667*xi*xi*xi*xi-1.33333333333333333333*xi*xi*xi+2.66666666666666666667*xi*xi+1.33333333333333333333*xi;
      double xi4 = 0.6666666666666667*xi*xi*xi*xi+0.6666666666666666666667*xi*xi*xi-0.1666666666666666666667*xi*xi-0.166666666666666666667*xi;
      double eta0 = 8*eta*eta-4*eta-0.333333333333333333333;
      double eta1 = -32*eta*eta+8*eta+5.33333333333333333333;
      double eta2 = 48.0*eta*eta-10.0;
      double eta3 = -32*eta*eta-8*eta+5.33333333333333333333;
      double eta4 = 8*eta*eta+4*eta-0.333333333333333333333;

      values[0] = xi0*eta0;
      values[1] = xi1*eta0;
      values[2] = xi2*eta0;
      values[3] = xi3*eta0;
      values[4] = xi4*eta0;
      values[5] = xi0*eta1;
      values[6] = xi1*eta1;
      values[7] = xi2*eta1;
      values[8] = xi3*eta1;
      values[9] = xi4*eta1;
      values[10] = xi0*eta2;
      values[11] = xi1*eta2;
      values[12] = xi2*eta2;
      values[13] = xi3*eta2;
      values[14] = xi4*eta2;
      values[15] = xi0*eta3;
      values[16] = xi1*eta3;
      values[17] = xi2*eta3;
      values[18] = xi3*eta3;
      values[19] = xi4*eta3;
      values[20] = xi0*eta4;
      values[21] = xi1*eta4;
      values[22] = xi2*eta4;
      values[23] = xi3*eta4;
      values[24] = xi4*eta4;
}


// ***********************************************************************

TBaseFunct2D *BF_C_Q_Q4_2D_Obj = new TBaseFunct2D
        (25, BF_C_Q_Q4_2D, BFUnitSquare, 
         C_Q_Q4_2D_Funct, C_Q_Q4_2D_DeriveXi,
         C_Q_Q4_2D_DeriveEta, C_Q_Q4_2D_DeriveXiXi,
         C_Q_Q4_2D_DeriveXiEta, C_Q_Q4_2D_DeriveEtaEta, 4, 4,
         0, NULL);
