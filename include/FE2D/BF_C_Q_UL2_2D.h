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
// Q2 element with bubbles, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_UL2_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
  double t30, t32, t33, t34, t35, t36, t37, t38, t40, t42;

  t1 = 5.0/16.0*eta;
  t2 = eta*eta;
  t3 = t2/16.0;
  t4 = 5.0/16.0*xi;
  t6 = xi*eta/4.0;
  t7 = xi*t2;
  t8 = 9.0/16.0*t7;
  t9 = xi*xi;
  t10 = t9/16.0;
  t11 = t9*eta;
  t12 = 9.0/16.0*t11;
  t13 = t9*t2;
  t14 = 3.0/16.0*t13;
  t15 = t9*xi;
  t16 = 5.0/16.0*t15;
  t17 = t15*t2;
  t18 = 5.0/16.0*t17;
  t19 = t2*eta;
  t20 = 5.0/16.0*t19;
  t21 = t9*t19;
  t22 = 5.0/16.0*t21;
  t23 = -1.0/16.0+t1+t3+t4+t6-t8+t10-t12+t14-t16+t18-t20+t22;
  t24 = 3.0/4.0*eta;
  t25 = 3.0/4.0*t2;
  t26 = 3.0/4.0*t11;
  t27 = t9/4.0;
  t28 = 3.0/4.0*t13;
  t29 = 5.0/4.0*t19;
  t30 = 5.0/4.0*t21;
  t32 = -1.0/16.0+t1+t3-t4-t6+t8+t10-t12+t14+t16-t18-t20+t22;
  t33 = 3.0/4.0*xi;
  t34 = 3.0/4.0*t7;
  t35 = 3.0/4.0*t9;
  t36 = t2/4.0;
  t37 = 5.0/4.0*t15;
  t38 = 5.0/4.0*t17;
  t40 = -1.0/16.0-t1+t3-t4+t6+t8+t10+t12+t14+t16-t18+t20-t22;
  t42 = -1.0/16.0-t1+t3+t4-t6-t8+t10+t12+t14-t16+t18+t20-t22;

  values[0] = t23;
  values[1] = -1.0/4.0+t24+t25-t26+t27-t28-t29+t30;
  values[2] = t32;
  values[3] = -1.0/4.0-t33+t34+t35+t36-t28+t37-t38;
  values[4] = t40;
  values[5] = -1.0/4.0-t24+t25+t26+t27-t28+t29-t30;
  values[6] = t42;
  values[7] = -1.0/4.0+t33-t34+t35+t36-t28-t37+t38;
  values[8] = 9.0/16.0-9.0/16.0*t2-9.0/16.0*t9+9.0/16.0*t13;
  values[9] = 45.0/16.0*eta-45.0/16.0*t19-45.0/16.0*t11+45.0/16.0*t21;
  values[10] = 45.0/16.0*xi-45.0/16.0*t7-45.0/16.0*t15+45.0/16.0*t17;
}

// values of the derivatives in xi direction
static void C_Q_UL2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t14, t15, t17;
  double t18, t19, t20, t23, t24, t25, t26;

  t1 = eta/4.0;
  t2 = eta*eta;
  t3 = 9.0/16.0*t2;
  t4 = xi/8.0;
  t5 = xi*eta;
  t6 = 9.0/8.0*t5;
  t7 = xi*t2;
  t8 = 3.0/8.0*t7;
  t9 = xi*xi;
  t10 = 15.0/16.0*t9;
  t11 = t9*t2;
  t12 = 15.0/16.0*t11;
  t14 = xi*t2*eta;
  t15 = 5.0/8.0*t14;
  t17 = 3.0/2.0*t5;
  t18 = xi/2.0;
  t19 = 3.0/2.0*t7;
  t20 = 5.0/2.0*t14;
  t23 = 3.0/4.0*t2;
  t24 = 3.0/2.0*xi;
  t25 = 15.0/4.0*t9;
  t26 = 15.0/4.0*t11;

  values[0] = 5.0/16.0+t1-t3+t4-t6+t8-t10+t12+t15;
  values[1] = -t17+t18-t19+t20;
  values[2] = -5.0/16.0-t1+t3+t4-t6+t8+t10-t12+t15;
  values[3] = -3.0/4.0+t23+t24-t19+t25-t26;
  values[4] = -5.0/16.0+t1+t3+t4+t6+t8+t10-t12-t15;
  values[5] = t17+t18-t19-t20;
  values[6] = 5.0/16.0-t1-t3+t4+t6+t8-t10+t12-t15;
  values[7] = 3.0/4.0-t23+t24-t19-t25+t26;
  values[8] = -9.0/8.0*xi+9.0/8.0*t7;
  values[9] = -45.0/8.0*t5+45.0/8.0*t14;
  values[10] = 45.0/16.0-45.0/16.0*t2-135.0/16.0*t9+135.0/16.0*t11;
}

// values of the derivatives in eta direction
static void C_Q_UL2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t12, t13, t14, t15, t17;
  double t18, t19, t20, t21, t24, t25, t26;

  t1 = eta/8.0;
  t2 = xi/4.0;
  t3 = xi*eta;
  t4 = 9.0/8.0*t3;
  t5 = xi*xi;
  t6 = 9.0/16.0*t5;
  t7 = t5*eta;
  t8 = 3.0/8.0*t7;
  t10 = t5*xi*eta;
  t11 = 5.0/8.0*t10;
  t12 = eta*eta;
  t13 = 15.0/16.0*t12;
  t14 = t5*t12;
  t15 = 15.0/16.0*t14;
  t17 = 3.0/2.0*eta;
  t18 = 3.0/4.0*t5;
  t19 = 3.0/2.0*t7;
  t20 = 15.0/4.0*t12;
  t21 = 15.0/4.0*t14;
  t24 = 3.0/2.0*t3;
  t25 = eta/2.0;
  t26 = 5.0/2.0*t10;

  values[0] = 5.0/16.0+t1+t2-t4-t6+t8+t11-t13+t15;
  values[1] = 3.0/4.0+t17-t18-t19-t20+t21;
  values[2] = 5.0/16.0+t1-t2+t4-t6+t8-t11-t13+t15;
  values[3] = t24+t25-t19-t26;
  values[4] = -5.0/16.0+t1+t2+t4+t6+t8-t11+t13-t15;
  values[5] = -3.0/4.0+t17+t18-t19+t20-t21;
  values[6] = -5.0/16.0+t1-t2-t4+t6+t8+t11+t13-t15;
  values[7] = -t24+t25-t19+t26;
  values[8] = -9.0/8.0*eta+9.0/8.0*t7;
  values[9] = 45.0/16.0-135.0/16.0*t12-45.0/16.0*t5+135.0/16.0*t14;
  values[10] = -45.0/8.0*t3+45.0/8.0*t10;
}

// values of the derivatives in xi-xi  direction
static void C_Q_UL2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t12, t15, t16;

  t1 = 9.0/8.0*eta;
  t2 = eta*eta;
  t3 = 3.0/8.0*t2;
  t4 = 15.0/8.0*xi;
  t5 = xi*t2;
  t6 = 15.0/8.0*t5;
  t7 = t2*eta;
  t8 = 5.0/8.0*t7;
  t10 = 3.0/2.0*eta;
  t11 = 3.0/2.0*t2;
  t12 = 5.0/2.0*t7;
  t15 = 15.0/2.0*xi;
  t16 = 15.0/2.0*t5;

  values[0] = 1.0/8.0-t1+t3-t4+t6+t8;
  values[1] = -t10+1.0/2.0-t11+t12;
  values[2] = 1.0/8.0-t1+t3+t4-t6+t8;
  values[3] = 3.0/2.0-t11+t15-t16;
  values[4] = 1.0/8.0+t1+t3+t4-t6-t8;
  values[5] = t10+1.0/2.0-t11-t12;
  values[6] = 1.0/8.0+t1+t3-t4+t6-t8;
  values[7] = 3.0/2.0-t11-t15+t16;
  values[8] = -9.0/8.0+9.0/8.0*t2;
  values[9] = -45.0/8.0*eta+45.0/8.0*t7;
  values[10] = -135.0/8.0*xi+135.0/8.0*t5;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t12, t13, t14, t17, t18;

  t1 = 9.0/8.0*eta;
  t2 = 9.0/8.0*xi;
  t3 = xi*eta;
  t4 = 3.0/4.0*t3;
  t5 = xi*xi;
  t6 = t5*eta;
  t7 = 15.0/8.0*t6;
  t8 = eta*eta;
  t9 = xi*t8;
  t10 = 15.0/8.0*t9;
  t12 = 3.0/2.0*xi;
  t13 = 3.0*t3;
  t14 = 15.0/2.0*t9;
  t17 = 3.0/2.0*eta;
  t18 = 15.0/2.0*t6;

  values[0] = 1.0/4.0-t1-t2+t4+t7+t10;
  values[1] = -t12-t13+t14;
  values[2] = -1.0/4.0+t1-t2+t4-t7+t10;
  values[3] = t17-t13-t18;
  values[4] = 1.0/4.0+t1+t2+t4-t7-t10;
  values[5] = t12-t13-t14;
  values[6] = -1.0/4.0-t1+t2+t4+t7-t10;
  values[7] = -t17-t13+t18;
  values[8] = 9.0/4.0*t3;
  values[9] = -45.0/8.0*xi+135.0/8.0*t9;
  values[10] = -45.0/8.0*eta+135.0/8.0*t6;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t12, t15, t16;

  t1 = 9.0/8.0*xi;
  t2 = xi*xi;
  t3 = 3.0/8.0*t2;
  t4 = xi*t2;
  t5 = 5.0/8.0*t4;
  t6 = 15.0/8.0*eta;
  t7 = t2*eta;
  t8 = 15.0/8.0*t7;
  t10 = 3.0/2.0*t2;
  t11 = 15.0/2.0*eta;
  t12 = 15.0/2.0*t7;
  t15 = 3.0/2.0*xi;
  t16 = 5.0/2.0*t4;

  values[0] = 1.0/8.0-t1+t3+t5-t6+t8;
  values[1] = 3.0/2.0-t10-t11+t12;
  values[2] = 1.0/8.0+t1+t3-t5-t6+t8;
  values[3] = t15+1.0/2.0-t10-t16;
  values[4] = 1.0/8.0+t1+t3-t5+t6-t8;
  values[5] = 3.0/2.0-t10+t11-t12;
  values[6] = 1.0/8.0-t1+t3+t5+t6-t8;
  values[7] = -t15+1.0/2.0-t10+t16;
  values[8] = -9.0/8.0+9.0/8.0*t2;
  values[9] = -135.0/8.0*eta+135.0/8.0*t7;
  values[10] = -45.0/8.0*xi+45.0/8.0*t4;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_UL2_2D_Obj = new TBaseFunct2D
        (11, BF_C_Q_UL2_2D, BFUnitSquare, 
         C_Q_UL2_2D_Funct, C_Q_UL2_2D_DeriveXi,
         C_Q_UL2_2D_DeriveEta, C_Q_UL2_2D_DeriveXiXi,
         C_Q_UL2_2D_DeriveXiEta, C_Q_UL2_2D_DeriveEtaEta, 3, 2,
         0, NULL);
