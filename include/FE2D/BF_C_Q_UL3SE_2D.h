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
// Q3 element with bubbles, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_UL3SE_2D_Funct(double xi, double eta, double *values)
{
  double t1 = xi/32.0;
  double t2 = xi*xi;
  double t3 = t2/32.0;
  double t4 = eta/32.0;
  double t6 = xi*eta/4.0;
  double t8 = 9.0/32.0*t2*eta;
  double t9 = eta*eta;
  double t10 = t9/32.0;
  double t12 = 9.0/32.0*xi*t9;
  double t13 = t2*t9;
  double t14 = 5.0/16.0*t13;
  double t16 = 1.0-t2;
  double t18 = (1.0-eta)*t16*xi;
  double t19 = 9.0/32.0*t18;
  double t21 = 1.0-t9;
  double t23 = (1.0-xi)*t21*eta;
  double t24 = 9.0/32.0*t23;
  double t26 = 9.0/32.0*eta;
  double t27 = 9.0/32.0*t9;
  double t28 = 9.0/32.0*t13;
  double t29 = 27.0/32.0*t18;
  double t34 = (1.0+xi)*t21*eta;
  double t35 = 9.0/32.0*t34;
  double t37 = 9.0/32.0*xi;
  double t38 = 9.0/32.0*t2;
  double t39 = 27.0/32.0*t34;
  double t44 = (1.0+eta)*t16*xi;
  double t45 = 9.0/32.0*t44;
  double t47 = 27.0/32.0*t44;
  double t51 = 27.0/32.0*t23;
  values[0] = t1-t3+t4+t6-t8-t10-t12+t14+t19+t24;
  values[1] = -t26+t8+t27-t28-t29;
  values[2] = -t26+t8+t27-t28+t29;
  values[3] = -t1-t3+t4-t6-t8-t10+t12+t14-t19+t35;
  values[4] = t37+t38-t12-t28-t39;
  values[5] = t37+t38-t12-t28+t39;
  values[6] = -t1-t3-t4+t6+t8-t10+t12+t14-t45-t35;
  values[7] = t26-t8+t27-t28+t47;
  values[8] = t26-t8+t27-t28-t47;
  values[9] = t1-t3-t4-t6+t8-t10-t12+t14+t45-t24;
  values[10] = -t37+t38+t12-t28+t51;
  values[11] = -t37+t38+t12-t28-t51;
  values[12] = 1.0-t2-t9+t13;
}

// values of the derivatives in xi direction
static void C_Q_UL3SE_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1 = xi/16.0;
  double t2 = eta/4.0;
  double t4 = 9.0/16.0*xi*eta;
  double t5 = eta*eta;
  double t6 = 9.0/32.0*t5;
  double t7 = xi*t5;
  double t8 = 5.0/8.0*t7;
  double t9 = 1.0-eta;
  double t10 = xi*xi;
  double t11 = t9*t10;
  double t12 = 9.0/16.0*t11;
  double t13 = 1.0-t10;
  double t14 = t9*t13;
  double t15 = 9.0/32.0*t14;
  double t17 = (1.0-t5)*eta;
  double t18 = 9.0/32.0*t17;
  double t20 = 9.0/16.0*t7;
  double t21 = 27.0/16.0*t11;
  double t22 = 27.0/32.0*t14;
  double t26 = 9.0/16.0*xi;
  double t27 = 27.0/32.0*t17;
  double t30 = 1.0+eta;
  double t31 = t30*t10;
  double t32 = 9.0/16.0*t31;
  double t33 = t30*t13;
  double t34 = 9.0/32.0*t33;
  double t36 = 27.0/16.0*t31;
  double t37 = 27.0/32.0*t33;
  values[0] = 1.0/32.0-t1+t2-t4-t6+t8-t12+t15-t18;
  values[1] = t4-t20+t21-t22;
  values[2] = t4-t20-t21+t22;
  values[3] = -1.0/32.0-t1-t2-t4+t6+t8+t12-t15+t18;
  values[4] = 9.0/32.0+t26-t6-t20-t27;
  values[5] = 9.0/32.0+t26-t6-t20+t27;
  values[6] = -1.0/32.0-t1+t2+t4+t6+t8+t32-t34-t18;
  values[7] = -t4-t20-t36+t37;
  values[8] = -t4-t20+t36-t37;
  values[9] = 1.0/32.0-t1-t2+t4-t6+t8-t32+t34+t18;
  values[10] = -9.0/32.0+t26+t6-t20-t27;
  values[11] = -9.0/32.0+t26+t6-t20+t27;
  values[12] = -2.0*xi+2.0*t7;
}

// values of the derivatives in eta direction
static void C_Q_UL3SE_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1 = xi/4.0;
  double t2 = xi*xi;
  double t3 = 9.0/32.0*t2;
  double t4 = eta/16.0;
  double t6 = 9.0/16.0*xi*eta;
  double t7 = t2*eta;
  double t8 = 5.0/8.0*t7;
  double t10 = (1.0-t2)*xi;
  double t11 = 9.0/32.0*t10;
  double t12 = 1.0-xi;
  double t13 = eta*eta;
  double t14 = t12*t13;
  double t15 = 9.0/16.0*t14;
  double t16 = 1.0-t13;
  double t17 = t12*t16;
  double t18 = 9.0/32.0*t17;
  double t20 = 9.0/16.0*eta;
  double t21 = 9.0/16.0*t7;
  double t22 = 27.0/32.0*t10;
  double t25 = 1.0+xi;
  double t26 = t25*t13;
  double t27 = 9.0/16.0*t26;
  double t28 = t25*t16;
  double t29 = 9.0/32.0*t28;
  double t31 = 27.0/16.0*t26;
  double t32 = 27.0/32.0*t28;
  double t39 = 27.0/16.0*t14;
  double t40 = 27.0/32.0*t17;
  values[0] = 1.0/32.0+t1-t3-t4-t6+t8-t11-t15+t18;
  values[1] = -9.0/32.0+t3+t20-t21+t22;
  values[2] = -9.0/32.0+t3+t20-t21-t22;
  values[3] = 1.0/32.0-t1-t3-t4+t6+t8+t11-t27+t29;
  values[4] = -t6-t21+t31-t32;
  values[5] = -t6-t21-t31+t32;
  values[6] = -1.0/32.0+t1+t3-t4+t6+t8-t11+t27-t29;
  values[7] = 9.0/32.0-t3+t20-t21+t22;
  values[8] = 9.0/32.0-t3+t20-t21-t22;
  values[9] = -1.0/32.0-t1+t3-t4-t6+t8+t11+t15-t18;
  values[10] = t6-t21-t39+t40;
  values[11] = t6-t21+t39-t40;
  values[12] = -2.0*eta+2.0*t7;
}

// values of the derivatives in xi-xi  direction
static void C_Q_UL3SE_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1 = 9.0/16.0*eta;
  double t2 = eta*eta;
  double t3 = 5.0/8.0*t2;
  double t5 = (1.0-eta)*xi;
  double t6 = 27.0/16.0*t5;
  double t8 = 9.0/16.0*t2;
  double t9 = 81.0/16.0*t5;
  double t13 = 1.0-t2;
  double t15 = (1.0+eta)*xi;
  double t16 = 27.0/16.0*t15;
  double t18 = 81.0/16.0*t15;
  values[0] = -1.0/16.0-t1+t3-t6;
  values[1] = t1-t8+t9;
  values[2] = t1-t8-t9;
  values[3] = -1.0/16.0-t1+t3+t6;
  values[4] = 9.0/16.0*t13;
  values[5] = 9.0/16.0*t13;
  values[6] = -1.0/16.0+t1+t3+t16;
  values[7] = -t1-t8-t18;
  values[8] = -t1-t8+t18;
  values[9] = -1.0/16.0+t1+t3-t16;
  values[10] = 9.0/16.0*t13;
  values[11] = 9.0/16.0*t13;
  values[12] = -2.0*t13;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL3SE_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = 9.0/16.0*xi;
  double t2 = 9.0/16.0*eta;
  double t3 = xi*eta;
  double t4 = 5.0/4.0*t3;
  double t5 = xi*xi;
  double t6 = 27.0/32.0*t5;
  double t7 = eta*eta;
  double t8 = 27.0/32.0*t7;
  double t10 = 9.0/8.0*t3;
  double t11 = 81.0/32.0*t5;
  double t15 = 81.0/32.0*t7;
  values[0] = -5.0/16.0-t1-t2+t4+t6+t8;
  values[1] = t1-t10-t11+27.0/32.0;
  values[2] = t1-t10+t11-27.0/32.0;
  values[3] = 5.0/16.0-t1+t2+t4-t6-t8;
  values[4] = -t2-t10+t15-27.0/32.0;
  values[5] = -t2-t10-t15+27.0/32.0;
  values[6] = -5.0/16.0+t1+t2+t4+t6+t8;
  values[7] = -t1-t10-t11+27.0/32.0;
  values[8] = -t1-t10+t11-27.0/32.0;
  values[9] = 5.0/16.0+t1-t2+t4-t6-t8;
  values[10] = t2-t10+t15-27.0/32.0;
  values[11] = t2-t10-t15+27.0/32.0;
  values[12] = 4.0*t3;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL3SE_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1 = 9.0/16.0*xi;
  double t2 = xi*xi;
  double t3 = 5.0/8.0*t2;
  double t5 = (1.0-xi)*eta;
  double t6 = 27.0/16.0*t5;
  double t8 = 1.0-t2;
  double t10 = (1.0+xi)*eta;
  double t11 = 27.0/16.0*t10;
  double t13 = 9.0/16.0*t2;
  double t14 = 81.0/16.0*t10;
  double t19 = 81.0/16.0*t5;
  values[0] = -1.0/16.0-t1+t3-t6;
  values[1] = 9.0/16.0*t8;
  values[2] = 9.0/16.0*t8;
  values[3] = -1.0/16.0+t1+t3-t11;
  values[4] = -t1-t13+t14;
  values[5] = -t1-t13-t14;
  values[6] = -1.0/16.0+t1+t3+t11;
  values[7] = 9.0/16.0*t8;
  values[8] = 9.0/16.0*t8;
  values[9] = -1.0/16.0-t1+t3+t6;
  values[10] = t1-t13-t19;
  values[11] = t1-t13+t19;
  values[12] = -2.0*t8;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_UL3SE_2D_Obj = new TBaseFunct2D
        (13, BF_C_Q_UL3SE_2D, BFUnitSquare, 
         C_Q_UL3SE_2D_Funct, C_Q_UL3SE_2D_DeriveXi,
         C_Q_UL3SE_2D_DeriveEta, C_Q_UL3SE_2D_DeriveXiXi,
         C_Q_UL3SE_2D_DeriveXiEta, C_Q_UL3SE_2D_DeriveEtaEta, 3, 3,
         0, NULL);
