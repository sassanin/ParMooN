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
// P4 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_P4_2D_Funct(double xi, double eta, double *values)
{
  double t3, t5, t8, t10, t14, t22, t23, t32;

  t3 = xi*xi;
  t5 = 15.0/4.0*t3-5.0/4.0;
  t8 = eta*eta;
  t10 = 15.0/4.0*t8-5.0/4.0;
  t14 = 35.0/4.0*t3*xi-21.0/4.0*xi;
  t22 = 35.0/4.0*t8*eta-21.0/4.0*eta;
  t23 = t3*t3;
  t32 = t8*t8;

  values[0] = 1.0;
  values[1] = 3.0*xi;
  values[2] = 3.0*eta;
  values[3] = t5;
  values[4] = 9.0*xi*eta;
  values[5] = t10;
  values[6] = t14;
  values[7] = 3.0*t5*eta;
  values[8] = 3.0*xi*t10;
  values[9] = t22;
  values[10] = 315.0/64.0*t23-135.0/32.0*t3+27.0/64.0;
  values[11] = 3.0*t14*eta;
  values[12] = t5*t10;
  values[13] = 3.0*xi*t22;
  values[14] = 315.0/64.0*t32-135.0/32.0*t8+27.0/64.0;
}

// values of the derivatives in xi direction
static void D_Q_P4_2D_DeriveXi(double xi, double eta, double *values)
{
  double t3, t5, t8;

  t3 = xi*xi;
  t5 = 105.0/4.0*t3-21.0/4.0;
  t8 = eta*eta;

  values[0] = 0.0;
  values[1] = 3.0;
  values[2] = 0.0;
  values[3] = 15.0/2.0*xi;
  values[4] = 9.0*eta;
  values[5] = 0.0;
  values[6] = t5;
  values[7] = 45.0/2.0*xi*eta;
  values[8] = 45.0/4.0*t8-15.0/4.0;
  values[9] = 0.0;
  values[10] = 315.0/16.0*t3*xi-135.0/16.0*xi;
  values[11] = 3.0*t5*eta;
  values[12] = 15.0/2.0*xi*(15.0/4.0*t8-5.0/4.0);
  values[13] = 105.0/4.0*t8*eta-63.0/4.0*eta;
  values[14] = 0.0;
}

// values of the derivatives in eta direction
static void D_Q_P4_2D_DeriveEta(double xi, double eta, double *values)
{
  double t3, t8, t10;

  t3 = xi*xi;
  t8 = eta*eta;
  t10 = 105.0/4.0*t8-21.0/4.0;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 3.0;
  values[3] = 0.0;
  values[4] = 9.0*xi;
  values[5] = 15.0/2.0*eta;
  values[6] = 0.0;
  values[7] = 45.0/4.0*t3-15.0/4.0;
  values[8] = 45.0/2.0*xi*eta;
  values[9] = t10;
  values[10] = 0.0;
  values[11] = 105.0/4.0*t3*xi-63.0/4.0*xi;
  values[12] = 15.0/2.0*(15.0/4.0*t3-5.0/4.0)*eta;
  values[13] = 3.0*xi*t10;
  values[14] = 315.0/16.0*t8*eta-135.0/16.0*eta;
}

// values of the derivatives in xi-xi direction
static void D_Q_P4_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t3, t8;

  t3 = xi*xi;
  t8 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 15.0/2.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 105.0/2.0*xi;
  values[7] = 45.0/2.0*eta;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 945.0/16.0*t3-135.0/16.0;
  values[11] = 315.0/2.0*xi*eta;
  values[12] = 225.0/8.0*t8-75.0/8.0;
  values[13] = 0.0;
  values[14] = 0.0;
}

// values of the derivatives in xi-eta direction
static void D_Q_P4_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t3, t8;

  t3 = xi*xi;
  t8 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 9.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 45.0/2.0*xi;
  values[8] = 45.0/2.0*eta;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 315.0/4.0*t3-63.0/4.0;
  values[12] = 225.0/4.0*xi*eta;
  values[13] = 315.0/4.0*t8-63.0/4.0;
  values[14] = 0.0;
}

// values of the derivatives in eta-eta direction
static void D_Q_P4_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t3, t8;

  t3 = xi*xi;
  t8 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 15.0/2.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 45.0/2.0*xi;
  values[9] = 105.0/2.0*eta;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 225.0/8.0*t3-75.0/8.0;
  values[13] = 315.0/2.0*xi*eta;
  values[14] = 945.0/16.0*t8-135.0/16.0;
}

// ***********************************************************************

TBaseFunct2D *BF_D_Q_P4_2D_Obj = new TBaseFunct2D
        (15, BF_D_Q_P4_2D, BFUnitSquare, 
         D_Q_P4_2D_Funct, D_Q_P4_2D_DeriveXi,
         D_Q_P4_2D_DeriveEta, D_Q_P4_2D_DeriveXiXi,
         D_Q_P4_2D_DeriveXiEta, D_Q_P4_2D_DeriveEtaEta, 4, 4,
         0, NULL);
