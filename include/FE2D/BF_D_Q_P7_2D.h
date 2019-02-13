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
// P7 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_P7_2D_Funct(double xi, double eta, double *values)
{
  double t1 = xi*xi;
  double t4 = xi*eta;
  double t5 = eta*eta;
  double t8 = t1*xi;
  double t13 = -1.0+3.0*t1;
  double t14 = t13*eta;
  double t17 = -1.0+3.0*t5;
  double t20 = t5*eta;
  double t24 = t1*t1;
  double t30 = xi*(5.0*t1-3.0);
  double t36 = 5.0*t5-3.0;
  double t39 = t5*t5;
  double t43 = t24*xi;
  double t50 = 3.0+35.0*t24-30.0*t1;
  double t51 = t50*eta;
  double t59 = 3.0+35.0*t39-30.0*t5;
  double t62 = t39*eta;
  double t67 = t24*t1;
  double t75 = xi*(63.0*t24-70.0*t1+15.0);
  double t87 = 63.0*t39-70.0*t5+15.0;
  double t90 = t39*t5;

  values[0] = 1.0;
  values[1] = xi;
  values[2] = eta;
  values[3] = -1.0/2.0+3.0/2.0*t1;
  values[4] = t4;
  values[5] = -1.0/2.0+3.0/2.0*t5;
  values[6] = 5.0/2.0*t8-3.0/2.0*xi;
  values[7] = t14/2.0;
  values[8] = xi*t17/2.0;
  values[9] = 5.0/2.0*t20-3.0/2.0*eta;
  values[10] = 3.0/8.0+35.0/8.0*t24-15.0/4.0*t1;
  values[11] = t30*eta/2.0;
  values[12] = t13*t17/4.0;
  values[13] = t4*t36/2.0;
  values[14] = 3.0/8.0+35.0/8.0*t39-15.0/4.0*t5;
  values[15] = 63.0/8.0*t43-35.0/4.0*t8+15.0/8.0*xi;
  values[16] = t51/8.0;
  values[17] = t30*t17/4.0;
  values[18] = t14*t36/4.0;
  values[19] = xi*t59/8.0;
  values[20] = 63.0/8.0*t62-35.0/4.0*t20+15.0/8.0*eta;
  values[21] = -5.0/16.0+231.0/16.0*t67-315.0/16.0*t24+105.0/16.0*t1;
  values[22] = t75*eta/8.0;
  values[23] = t50*t17/16.0;
  values[24] = t30*eta*t36/4.0;
  values[25] = t13*t59/16.0;
  values[26] = t4*t87/8.0;
  values[27] = -5.0/16.0+231.0/16.0*t90-315.0/16.0*t39+105.0/16.0*t5;
  values[28] = 429.0/16.0*t24*t8-693.0/16.0*t43+315.0/16.0*t8-35.0/16.0*xi;
  values[29] = (-5.0+231.0*t67-315.0*t24+105.0*t1)*eta/16.0;
  values[30] = t75*t17/16.0;
  values[31] = t51*t36/16.0;
  values[32] = t30*t59/16.0;
  values[33] = t14*t87/16.0;
  values[34] = xi*(-5.0+231.0*t90-315.0*t39+105.0*t5)/16.0;
  values[35] = 429.0/16.0*t39*t20-693.0/16.0*t62+315.0/16.0*t20-35.0/16.0*eta;
}

// values of the derivatives in xi direction
static void D_Q_P7_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2 = xi*xi;
  double t5 = xi*eta;
  double t7 = eta*eta;
  double t10 = t2*xi;
  double t14 = t2*eta;
  double t19 = -1.0+3.0*t7;
  double t23 = 5.0*t7-3.0;
  double t24 = eta*t23;
  double t26 = t2*t2;
  double t32 = xi*(7.0*t2-3.0);
  double t36 = t2*t7;
  double t42 = t7*t7;
  double t58 = t7*eta;
  double t72 = 63.0*t42-70.0*t7+15.0;

  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 0.0;
  values[3] = 3.0*xi;
  values[4] = eta;
  values[5] = 0.0;
  values[6] = 15.0/2.0*t2-3.0/2.0;
  values[7] = 3.0*t5;
  values[8] = -1.0/2.0+3.0/2.0*t7;
  values[9] = 0.0;
  values[10] = 35.0/2.0*t10-15.0/2.0*xi;
  values[11] = 15.0/2.0*t14-3.0/2.0*eta;
  values[12] = 3.0/2.0*xi*t19;
  values[13] = t24/2.0;
  values[14] = 0.0;
  values[15] = 315.0/8.0*t26-105.0/4.0*t2+15.0/8.0;
  values[16] = 5.0/2.0*t32*eta;
  values[17] = -15.0/4.0*t2+45.0/4.0*t36+3.0/4.0-9.0/4.0*t7;
  values[18] = 3.0/2.0*t5*t23;
  values[19] = 3.0/8.0+35.0/8.0*t42-15.0/4.0*t7;
  values[20] = 0.0;
  values[21] = 693.0/8.0*t26*xi-315.0/4.0*t10+105.0/8.0*xi;
  values[22] = 315.0/8.0*eta*t26-105.0/4.0*t14+15.0/8.0*eta;
  values[23] = 5.0/4.0*t32*t19;
  values[24] = 75.0/4.0*t2*t58-45.0/4.0*t14-15.0/4.0*t58+9.0/4.0*eta;
  values[25] = 3.0/8.0*xi*(3.0+35.0*t42-30.0*t7);
  values[26] = eta*t72/8.0;
  values[27] = 0.0;
  values[28] = 3003.0/16.0*t26*t2-3465.0/16.0*t26+945.0/16.0*t2-35.0/16.0;
  values[29] = 21.0/8.0*xi*(33.0*t26-30.0*t2+5.0)*eta;
  values[30] = -315.0/16.0*t26+945.0/16.0*t26*t7+105.0/8.0*t2-315.0/8.0*t36-15.0/16.0+45.0/16.0*t7;
  values[31] = 5.0/4.0*t32*t24;
  values[32] = 45.0/16.0*t2+525.0/16.0*t2*t42-225.0/8.0*t36-9.0/16.0-105.0/16.0*t42+45.0/8.0*t7;
  values[33] = 3.0/8.0*t5*t72;
  values[34] = -5.0/16.0+231.0/16.0*t42*t7-315.0/16.0*t42+105.0/16.0*t7;
  values[35] = 0.0;
}

// values of the derivatives in eta direction
static void D_Q_P7_2D_DeriveEta(double xi, double eta, double *values)
{
  double t2 = xi*xi;
  double t5 = xi*eta;
  double t7 = eta*eta;
  double t12 = xi*(5.0*t2-3.0);
  double t16 = (-1.0+3.0*t2)*eta;
  double t18 = xi*t7;
  double t22 = t7*eta;
  double t26 = t2*t2;
  double t33 = t2*t7;
  double t38 = 7.0*t7-3.0;
  double t41 = t7*t7;
  double t48 = xi*(63.0*t26-70.0*t2+15.0);
  double t55 = t2*xi;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 1.0;
  values[3] = 0.0;
  values[4] = xi;
  values[5] = 3.0*eta;
  values[6] = 0.0;
  values[7] = -1.0/2.0+3.0/2.0*t2;
  values[8] = 3.0*t5;
  values[9] = 15.0/2.0*t7-3.0/2.0;
  values[10] = 0.0;
  values[11] = t12/2.0;
  values[12] = 3.0/2.0*t16;
  values[13] = 15.0/2.0*t18-3.0/2.0*xi;
  values[14] = 35.0/2.0*t22-15.0/2.0*eta;
  values[15] = 0.0;
  values[16] = 3.0/8.0+35.0/8.0*t26-15.0/4.0*t2;
  values[17] = 3.0/2.0*t12*eta;
  values[18] = -15.0/4.0*t7+3.0/4.0+45.0/4.0*t33-9.0/4.0*t2;
  values[19] = 5.0/2.0*t5*t38;
  values[20] = 315.0/8.0*t41-105.0/4.0*t7+15.0/8.0;
  values[21] = 0.0;
  values[22] = t48/8.0;
  values[23] = 3.0/8.0*(3.0+35.0*t26-30.0*t2)*eta;
  values[24] = 75.0/4.0*t55*t7-15.0/4.0*t55-45.0/4.0*t18+9.0/4.0*xi;
  values[25] = 5.0/4.0*t16*t38;
  values[26] = 315.0/8.0*xi*t41-105.0/4.0*t18+15.0/8.0*xi;
  values[27] = 693.0/8.0*t41*eta-315.0/4.0*t22+105.0/8.0*eta;
  values[28] = 0.0;
  values[29] = -5.0/16.0+231.0/16.0*t26*t2-315.0/16.0*t26+105.0/16.0*t2;
  values[30] = 3.0/8.0*t48*eta;
  values[31] = 45.0/16.0*t7-9.0/16.0+525.0/16.0*t26*t7-105.0/16.0*t26-225.0/8.0*t33+45.0/8.0*t2;
  values[32] = 5.0/4.0*t12*eta*t38;
  values[33] = -315.0/16.0*t41+105.0/8.0*t7-15.0/16.0+945.0/16.0*t2*t41-315.0/8.0*t33+45.0/16.0*t2;
  values[34] = 21.0/8.0*t5*(33.0*t41-30.0*t7+5.0);
  values[35] = 3003.0/16.0*t41*t7-3465.0/16.0*t41+945.0/16.0*t7-35.0/16.0;
}

// values of the derivatives in xi-xi direction
static void D_Q_P7_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t3 = xi*xi;
  double t6 = xi*eta;
  double t8 = eta*eta;
  double t11 = t3*xi;
  double t15 = t3*eta;
  double t20 = xi*t8;
  double t27 = t3*t3;
  double t38 = 45.0/4.0*t8;
  double t40 = t8*eta;
  double t45 = t8*t8;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 3.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 15.0*xi;
  values[7] = 3.0*eta;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 105.0/2.0*t3-15.0/2.0;
  values[11] = 15.0*t6;
  values[12] = -3.0/2.0+9.0/2.0*t8;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 315.0/2.0*t11-105.0/2.0*xi;
  values[16] = 105.0/2.0*t15-15.0/2.0*eta;
  values[17] = -15.0/2.0*xi+45.0/2.0*t20;
  values[18] = 3.0/2.0*eta*(5.0*t8-3.0);
  values[19] = 0.0;
  values[20] = 0.0;
  values[21] = 3465.0/8.0*t27-945.0/4.0*t3+105.0/8.0;
  values[22] = 315.0/2.0*eta*t11-105.0/2.0*t6;
  values[23] = -105.0/4.0*t3+315.0/4.0*t3*t8+15.0/4.0-t38;
  values[24] = 75.0/2.0*xi*t40-45.0/2.0*t6;
  values[25] = 9.0/8.0+105.0/8.0*t45-t38;
  values[26] = 0.0;
  values[27] = 0.0;
  values[28] = 9009.0/8.0*t27*xi-3465.0/4.0*t11+945.0/8.0*xi;
  values[29] = 3465.0/8.0*eta*t27-945.0/4.0*t15+105.0/8.0*eta;
  values[30] = -315.0/4.0*t11+945.0/4.0*t11*t8+105.0/4.0*xi-315.0/4.0*t20;
  values[31] = 525.0/4.0*t3*t40-315.0/4.0*t15-75.0/4.0*t40+45.0/4.0*eta;
  values[32] = 45.0/8.0*xi+525.0/8.0*xi*t45-225.0/4.0*t20;
  values[33] = 3.0/8.0*eta*(63.0*t45-70.0*t8+15.0);
  values[34] = 0.0;
  values[35] = 0.0;
}

// values of the derivatives in xi-eta direction
static void D_Q_P7_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t3 = xi*xi;
  double t6 = xi*eta;
  double t8 = eta*eta;
  double t13 = xi*(7.0*t3-3.0);
  double t15 = t3*eta;
  double t19 = xi*t8;
  double t23 = t8*eta;
  double t27 = t3*t3;
  double t42 = t8*t8;
  double t56 = t3*xi;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 1.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 3.0*xi;
  values[8] = 3.0*eta;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 15.0/2.0*t3-3.0/2.0;
  values[12] = 9.0*t6;
  values[13] = 15.0/2.0*t8-3.0/2.0;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = 5.0/2.0*t13;
  values[17] = 45.0/2.0*t15-9.0/2.0*eta;
  values[18] = 45.0/2.0*t19-9.0/2.0*xi;
  values[19] = 35.0/2.0*t23-15.0/2.0*eta;
  values[20] = 0.0;
  values[21] = 0.0;
  values[22] = 315.0/8.0*t27-105.0/4.0*t3+15.0/8.0;
  values[23] = 15.0/2.0*t13*eta;
  values[24] = 225.0/4.0*t3*t8-45.0/4.0*t3-45.0/4.0*t8+9.0/4.0;
  values[25] = 15.0/2.0*t6*(7.0*t8-3.0);
  values[26] = 315.0/8.0*t42-105.0/4.0*t8+15.0/8.0;
  values[27] = 0.0;
  values[28] = 0.0;
  values[29] = 21.0/8.0*xi*(33.0*t27-30.0*t3+5.0);
  values[30] = 945.0/8.0*eta*t27-315.0/4.0*t15+45.0/8.0*eta;
  values[31] = 525.0/4.0*t56*t8-105.0/4.0*t56-225.0/4.0*t19+45.0/4.0*xi;
  values[32] = 525.0/4.0*t3*t23-225.0/4.0*t15-105.0/4.0*t23+45.0/4.0*eta;
  values[33] = 945.0/8.0*xi*t42-315.0/4.0*t19+45.0/8.0*xi;
  values[34] = 693.0/8.0*t42*eta-315.0/4.0*t23+105.0/8.0*eta;
  values[35] = 0.0;
}

// values of the derivatives in eta-eta direction
static void D_Q_P7_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t3 = xi*xi;
  double t6 = xi*eta;
  double t8 = eta*eta;
  double t16 = t3*eta;
  double t19 = xi*t8;
  double t23 = t8*eta;
  double t27 = t3*t3;
  double t29 = 45.0/4.0*t3;
  double t31 = t3*xi;
  double t44 = t8*t8;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 3.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 3.0*xi;
  values[9] = 15.0*eta;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 9.0/2.0*t3-3.0/2.0;
  values[13] = 15.0*t6;
  values[14] = 105.0/2.0*t8-15.0/2.0;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 3.0/2.0*xi*(5.0*t3-3.0);
  values[18] = -15.0/2.0*eta+45.0/2.0*t16;
  values[19] = 105.0/2.0*t19-15.0/2.0*xi;
  values[20] = 315.0/2.0*t23-105.0/2.0*eta;
  values[21] = 0.0;
  values[22] = 0.0;
  values[23] = 9.0/8.0+105.0/8.0*t27-t29;
  values[24] = 75.0/2.0*eta*t31-45.0/2.0*t6;
  values[25] = -105.0/4.0*t8+15.0/4.0+315.0/4.0*t3*t8-t29;
  values[26] = 315.0/2.0*xi*t23-105.0/2.0*t6;
  values[27] = 3465.0/8.0*t44-945.0/4.0*t8+105.0/8.0;
  values[28] = 0.0;
  values[29] = 0.0;
  values[30] = 3.0/8.0*xi*(63.0*t27-70.0*t3+15.0);
  values[31] = 45.0/8.0*eta+525.0/8.0*eta*t27-225.0/4.0*t16;
  values[32] = 525.0/4.0*t31*t8-75.0/4.0*t31-315.0/4.0*t19+45.0/4.0*xi;
  values[33] = -315.0/4.0*t23+105.0/4.0*eta+945.0/4.0*t3*t23-315.0/4.0*t16;
  values[34] = 3465.0/8.0*xi*t44-945.0/4.0*t19+105.0/8.0*xi;
  values[35] = 9009.0/8.0*t44*eta-3465.0/4.0*t23+945.0/8.0*eta;
}

// ***********************************************************************

TBaseFunct2D *BF_D_Q_P7_2D_Obj = new TBaseFunct2D
        (36, BF_D_Q_P7_2D, BFUnitSquare, 
         D_Q_P7_2D_Funct, D_Q_P7_2D_DeriveXi,
         D_Q_P7_2D_DeriveEta, D_Q_P7_2D_DeriveXiXi,
         D_Q_P7_2D_DeriveXiEta, D_Q_P7_2D_DeriveEtaEta, 7, 7,
         0, NULL);
