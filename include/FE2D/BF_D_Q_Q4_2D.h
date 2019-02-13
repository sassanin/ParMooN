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
// Q4 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_Q4_2D_Funct(double xi, double eta, double *values)
{
  double t1, t8, t12, t14, t15, t19, t24, t25, t27, t31, t45, t55, t61;

  t1 = xi*xi;
  t8 = t1*t1;
  t12 = xi*eta;
  t14 = -1.0+3.0*t1;
  t15 = t14*eta;
  t19 = xi*(5.0*t1-3.0);
  t24 = 3.0+35.0*t8-30.0*t1;
  t25 = t24*eta;
  t27 = eta*eta;
  t31 = -1.0+3.0*t27;
  t45 = 5.0*t27-3.0;
  t55 = t27*t27;
  t61 = 3.0+35.0*t55-30.0*t27;

  values[0] = 1.0;
  values[1] = xi;
  values[2] = -1.0/2.0+3.0/2.0*t1;
  values[3] = 5.0/2.0*t1*xi-3.0/2.0*xi;
  values[4] = 3.0/8.0+35.0/8.0*t8-15.0/4.0*t1;
  values[5] = eta;
  values[6] = t12;
  values[7] = t15/2.0;
  values[8] = t19*eta/2.0;
  values[9] = t25/8.0;
  values[10] = -1.0/2.0+3.0/2.0*t27;
  values[11] = xi*t31/2.0;
  values[12] = t14*t31/4.0;
  values[13] = t19*t31/4.0;
  values[14] = t24*t31/16.0;
  values[15] = 5.0/2.0*t27*eta-3.0/2.0*eta;
  values[16] = t12*t45/2.0;
  values[17] = t15*t45/4.0;
  values[18] = t19*eta*t45/4.0;
  values[19] = t25*t45/16.0;
  values[20] = 3.0/8.0+35.0/8.0*t55-15.0/4.0*t27;
  values[21] = xi*t61/8.0;
  values[22] = t14*t61/16.0;
  values[23] = t19*t61/16.0;
  values[24] = t24*t61/64.0;
}

// values of the derivatives in xi direction
static void D_Q_Q4_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2, t9, t11, t17, t20, t24, t28, t35, t36, t40, t49, t55;

  t2 = xi*xi;
  t9 = xi*eta;
  t11 = t2*eta;
  t17 = xi*(7.0*t2-3.0);
  t20 = eta*eta;
  t24 = -1.0+3.0*t20;
  t28 = t2*t20;
  t35 = 5.0*t20-3.0;
  t36 = eta*t35;
  t40 = t20*eta;
  t49 = t20*t20;
  t55 = 3.0+35.0*t49-30.0*t20;

  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 3.0*xi;
  values[3] = 15.0/2.0*t2-3.0/2.0;
  values[4] = 35.0/2.0*t2*xi-15.0/2.0*xi;
  values[5] = 0.0;
  values[6] = eta;
  values[7] = 3.0*t9;
  values[8] = 15.0/2.0*t11-3.0/2.0*eta;
  values[9] = 5.0/2.0*t17*eta;
  values[10] = 0.0;
  values[11] = -1.0/2.0+3.0/2.0*t20;
  values[12] = 3.0/2.0*xi*t24;
  values[13] = -15.0/4.0*t2+45.0/4.0*t28+3.0/4.0-9.0/4.0*t20;
  values[14] = 5.0/4.0*t17*t24;
  values[15] = 0.0;
  values[16] = t36/2.0;
  values[17] = 3.0/2.0*t9*t35;
  values[18] = 75.0/4.0*t2*t40-45.0/4.0*t11-15.0/4.0*t40+9.0/4.0*eta;
  values[19] = 5.0/4.0*t17*t36;
  values[20] = 0.0;
  values[21] = 3.0/8.0+35.0/8.0*t49-15.0/4.0*t20;
  values[22] = 3.0/8.0*xi*t55;
  values[23] = 45.0/16.0*t2+525.0/16.0*t2*t49-225.0/8.0*t28-9.0/16.0-105.0/16.0*t49+45.0/8.0*t20;
  values[24] = 5.0/16.0*t17*t55;
}

// values of the derivatives in eta direction
static void D_Q_Q4_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t6, t8, t13, t17, t24, t26, t29, t34, t38, t57;

  t1 = xi*xi;
  t6 = xi*(5.0*t1-3.0);
  t8 = t1*t1;
  t13 = xi*eta;
  t17 = (-1.0+3.0*t1)*eta;
  t24 = (3.0+35.0*t8-30.0*t1)*eta;
  t26 = eta*eta;
  t29 = xi*t26;
  t34 = t1*t26;
  t38 = t1*xi;
  t57 = 7.0*t26-3.0;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 1.0;
  values[6] = xi;
  values[7] = -1.0/2.0+3.0/2.0*t1;
  values[8] = t6/2.0;
  values[9] = 3.0/8.0+35.0/8.0*t8-15.0/4.0*t1;
  values[10] = 3.0*eta;
  values[11] = 3.0*t13;
  values[12] = 3.0/2.0*t17;
  values[13] = 3.0/2.0*t6*eta;
  values[14] = 3.0/8.0*t24;
  values[15] = 15.0/2.0*t26-3.0/2.0;
  values[16] = 15.0/2.0*t29-3.0/2.0*xi;
  values[17] = -15.0/4.0*t26+3.0/4.0+45.0/4.0*t34-9.0/4.0*t1;
  values[18] = 75.0/4.0*t38*t26-15.0/4.0*t38-45.0/4.0*t29+9.0/4.0*xi;
  values[19] = 45.0/16.0*t26-9.0/16.0+525.0/16.0*t8*t26-105.0/16.0*t8-225.0/8.0*t34+45.0/8.0*t1;
  values[20] = 35.0/2.0*t26*eta-15.0/2.0*eta;
  values[21] = 5.0/2.0*t13*t57;
  values[22] = 5.0/4.0*t17*t57;
  values[23] = 5.0/4.0*t6*eta*t57;
  values[24] = 5.0/16.0*t24*t57;
}

// values of the derivatives in xi-xi direction
static void D_Q_Q4_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t2, t6, t8, t12, t16, t20, t22, t28, t39;

  t2 = xi*xi;
  t6 = xi*eta;
  t8 = t2*eta;
  t12 = eta*eta;
  t16 = xi*t12;
  t20 = t2*t12;
  t22 = 45.0/4.0*t12;
  t28 = t12*eta;
  t39 = t12*t12;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 3.0;
  values[3] = 15.0*xi;
  values[4] = 105.0/2.0*t2-15.0/2.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 3.0*eta;
  values[8] = 15.0*t6;
  values[9] = 105.0/2.0*t8-15.0/2.0*eta;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = -3.0/2.0+9.0/2.0*t12;
  values[13] = -15.0/2.0*xi+45.0/2.0*t16;
  values[14] = -105.0/4.0*t2+315.0/4.0*t20+15.0/4.0-t22;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 3.0/2.0*eta*(5.0*t12-3.0);
  values[18] = 75.0/2.0*xi*t28-45.0/2.0*t6;
  values[19] = 525.0/4.0*t2*t28-315.0/4.0*t8-75.0/4.0*t28+45.0/4.0*eta;
  values[20] = 0.0;
  values[21] = 0.0;
  values[22] = 9.0/8.0+105.0/8.0*t39-t22;
  values[23] = 45.0/8.0*xi+525.0/8.0*xi*t39-225.0/4.0*t16;
  values[24] = 315.0/16.0*t2+3675.0/16.0*t2*t39-1575.0/8.0*t20-45.0/16.0-525.0/16.0*t39+225.0/8.0*t12;
}

// values of the derivatives in xi-eta direction
static void D_Q_Q4_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t2, t7, t10, t12, t18, t21, t30, t37, t42;

  t2 = xi*xi;
  t7 = xi*(7.0*t2-3.0);
  t10 = xi*eta;
  t12 = t2*eta;
  t18 = eta*eta;
  t21 = xi*t18;
  t30 = t2*xi;
  t37 = t18*eta;
  t42 = 7.0*t18-3.0;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 1.0;
  values[7] = 3.0*xi;
  values[8] = 15.0/2.0*t2-3.0/2.0;
  values[9] = 5.0/2.0*t7;
  values[10] = 0.0;
  values[11] = 3.0*eta;
  values[12] = 9.0*t10;
  values[13] = 45.0/2.0*t12-9.0/2.0*eta;
  values[14] = 15.0/2.0*t7*eta;
  values[15] = 0.0;
  values[16] = 15.0/2.0*t18-3.0/2.0;
  values[17] = 45.0/2.0*t21-9.0/2.0*xi;
  values[18] = 225.0/4.0*t2*t18-45.0/4.0*t2-45.0/4.0*t18+9.0/4.0;
  values[19] = 525.0/4.0*t30*t18-105.0/4.0*t30-225.0/4.0*t21+45.0/4.0*xi;
  values[20] = 0.0;
  values[21] = 35.0/2.0*t37-15.0/2.0*eta;
  values[22] = 15.0/2.0*t10*t42;
  values[23] = 525.0/4.0*t2*t37-225.0/4.0*t12-105.0/4.0*t37+45.0/4.0*eta;
  values[24] = 25.0/4.0*t7*eta*t42;
}

// values of the derivatives in eta-eta direction
static void D_Q_Q4_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t2, t9, t11, t14, t17, t20, t30, t33, t38;

  t2 = xi*xi;
  t9 = t2*t2;
  t11 = 45.0/4.0*t2;
  t14 = xi*eta;
  t17 = t2*eta;
  t20 = t2*xi;
  t30 = eta*eta;
  t33 = xi*t30;
  t38 = t2*t30;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 3.0;
  values[11] = 3.0*xi;
  values[12] = 9.0/2.0*t2-3.0/2.0;
  values[13] = 3.0/2.0*xi*(5.0*t2-3.0);
  values[14] = 9.0/8.0+105.0/8.0*t9-t11;
  values[15] = 15.0*eta;
  values[16] = 15.0*t14;
  values[17] = -15.0/2.0*eta+45.0/2.0*t17;
  values[18] = 75.0/2.0*t20*eta-45.0/2.0*t14;
  values[19] = 45.0/8.0*eta+525.0/8.0*t9*eta-225.0/4.0*t17;
  values[20] = 105.0/2.0*t30-15.0/2.0;
  values[21] = 105.0/2.0*t33-15.0/2.0*xi;
  values[22] = -105.0/4.0*t30+15.0/4.0+315.0/4.0*t38-t11;
  values[23] = 525.0/4.0*t20*t30-75.0/4.0*t20-315.0/4.0*t33+45.0/4.0*xi;
  values[24] = 315.0/16.0*t30+3675.0/16.0*t9*t30-525.0/16.0*t9-45.0/16.0-1575.0/8.0*t38+225.0/8.0*t2;
}

// ***********************************************************************

TBaseFunct2D *BF_D_Q_Q4_2D_Obj = new TBaseFunct2D
        (25, BF_D_Q_Q4_2D, BFUnitSquare, 
         D_Q_Q4_2D_Funct, D_Q_Q4_2D_DeriveXi,
         D_Q_Q4_2D_DeriveEta, D_Q_Q4_2D_DeriveXiXi,
         D_Q_Q4_2D_DeriveXiEta, D_Q_Q4_2D_DeriveEtaEta, 4, 4,
         0, NULL);
