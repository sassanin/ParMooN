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
// P2 element with bubbles, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_UL2_2D_Funct(double xi, double eta, double *values)
{
  double t3, t4, t5, t6, t7, t8, t9, t10, t13, t14, t16, t17, t21, t22;
  double t24, t27, t34;

  t3 = xi*xi;
  t4 = 2.0*t3;
  t5 = xi*eta;
  t6 = 4.0*t5;
  t7 = eta*eta;
  t8 = 2.0*t7;
  t9 = 1.0-xi-eta;
  t10 = t5*t9;
  t13 = t3*eta*t9;
  t14 = 63.0*t13;
  t16 = xi*t7*t9;
  t17 = 63.0*t16;
  t21 = 48.0*t10;
  t22 = 84.0*t16;
  t24 = 21.0*t10;
  t27 = 84.0*t13;
  t34 = 105.0*t10;

  values[0] = 1.0-3.0*xi-3.0*eta+t4+t6+t8-42.0*t10+t14+t17;
  values[1] = 4.0*xi-4.0*t3-t6-t21+t22;
  values[2] = -xi+t4+t24-t14;
  values[3] = t6+36.0*t10-t27-t22;
  values[4] = -eta+t8+t24-t17;
  values[5] = 4.0*eta-t6-4.0*t7-t21+t27;
  values[6] = 60.0*t10;
  values[7] = -t34+210.0*t13+105.0*t16;
  values[8] = -t34+105.0*t13+210.0*t16;
}

// values of the derivatives in xi direction
static void C_T_UL2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t6, t8, t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t20, t21, t22, t23, t25, t26, t30, t31, t36, t37;

  t1 = 4.0*xi;
  t2 = 4.0*eta;
  t3 = 1.0-xi-eta;
  t4 = eta*t3;
  t6 = xi*eta;
  t8 = t6*t3;
  t9 = 126.0*t8;
  t10 = xi*xi;
  t11 = t10*eta;
  t12 = 63.0*t11;
  t13 = eta*eta;
  t14 = t13*t3;
  t15 = 63.0*t14;
  t16 = xi*t13;
  t17 = 63.0*t16;
  t20 = 48.0*t4;
  t21 = 48.0*t6;
  t22 = 84.0*t14;
  t23 = 84.0*t16;
  t25 = 21.0*t4;
  t26 = 21.0*t6;
  t30 = 168.0*t8;
  t31 = 84.0*t11;
  t36 = 105.0*t4;
  t37 = 105.0*t6;

  values[0] = -3.0+t1+t2-42.0*t4+42.0*t6+t9-t12+t15-t17;
  values[1] = 4.0-8.0*xi-t2-t20+t21+t22-t23;
  values[2] = -1.0+t1+t25-t26-t9+t12;
  values[3] = t2+36.0*t4-36.0*t6-t30+t31-t22+t23;
  values[4] = t25-t26-t15+t17;
  values[5] = -t2-t20+t21+t30-t31;
  values[6] = 60.0*t4-60.0*t6;
  values[7] = -t36+t37+420.0*t8-210.0*t11+105.0*t14-105.0*t16;
  values[8] = -t36+t37+210.0*t8-105.0*t11+210.0*t14-210.0*t16;
}

// values of the derivatives in eta direction
static void C_T_UL2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t6, t8, t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t19, t20, t21, t22, t24, t25, t29, t30, t36, t37;

  t1 = 4.0*xi;
  t2 = 4.0*eta;
  t3 = 1.0-xi-eta;
  t4 = xi*t3;
  t6 = xi*eta;
  t8 = xi*xi;
  t9 = t8*t3;
  t10 = 63.0*t9;
  t11 = t8*eta;
  t12 = 63.0*t11;
  t13 = t6*t3;
  t14 = 126.0*t13;
  t15 = eta*eta;
  t16 = xi*t15;
  t17 = 63.0*t16;
  t19 = 48.0*t4;
  t20 = 48.0*t6;
  t21 = 168.0*t13;
  t22 = 84.0*t16;
  t24 = 21.0*t4;
  t25 = 21.0*t6;
  t29 = 84.0*t9;
  t30 = 84.0*t11;
  t36 = 105.0*t4;
  t37 = 105.0*t6;

  values[0] = -3.0+t1+t2-42.0*t4+42.0*t6+t10-t12+t14-t17;
  values[1] = -t1-t19+t20+t21-t22;
  values[2] = t24-t25-t10+t12;
  values[3] = t1+36.0*t4-36.0*t6-t29+t30-t21+t22;
  values[4] = -1.0+t2+t24-t25-t14+t17;
  values[5] = 4.0-t1-8.0*eta-t19+t20+t29-t30;
  values[6] = 60.0*t4-60.0*t6;
  values[7] = -t36+t37+210.0*t9-210.0*t11+210.0*t13-105.0*t16;
  values[8] = -t36+t37+105.0*t9-105.0*t11+420.0*t13-210.0*t16;
}

// values of the derivatives in xi-xi  direction
static void C_T_UL2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t3, t4, t5, t6, t7, t8, t10, t11, t13, t16, t17, t22;

  t3 = eta*(1.0-xi-eta);
  t4 = 126.0*t3;
  t5 = xi*eta;
  t6 = 252.0*t5;
  t7 = eta*eta;
  t8 = 126.0*t7;
  t10 = 96.0*eta;
  t11 = 168.0*t7;
  t13 = 42.0*eta;
  t16 = 168.0*t3;
  t17 = 336.0*t5;
  t22 = 210.0*eta;

  values[0] = 4.0+84.0*eta+t4-t6-t8;
  values[1] = -8.0+t10-t11;
  values[2] = 4.0-t13-t4+t6;
  values[3] = -72.0*eta-t16+t17+t11;
  values[4] = -t13+t8;
  values[5] = t10+t16-t17;
  values[6] = -120.0*eta;
  values[7] = t22+420.0*t3-840.0*t5-210.0*t7;
  values[8] = t22+210.0*t3-420.0*t5-420.0*t7;
}

// values of the derivatives in xi-eta direction
static void C_T_UL2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t3, t4, t5, t6, t8, t9, t10, t11, t12, t13, t15, t16, t17, t18;
  double t19, t21, t22, t23, t27, t29, t36, t37, t39;

  t3 = 1.0-xi-eta;
  t4 = xi*t3;
  t5 = 126.0*t4;
  t6 = xi*eta;
  t8 = xi*xi;
  t9 = 63.0*t8;
  t10 = eta*t3;
  t11 = 126.0*t10;
  t12 = eta*eta;
  t13 = 63.0*t12;
  t15 = 96.0*xi;
  t16 = 96.0*eta;
  t17 = 168.0*t10;
  t18 = 84.0*t12;
  t19 = 168.0*t6;
  t21 = 42.0*xi;
  t22 = 42.0*eta;
  t23 = 126.0*t6;
  t27 = 168.0*t4;
  t29 = 84.0*t8;
  t36 = 210.0*xi;
  t37 = 210.0*eta;
  t39 = 630.0*t6;

  values[0] = -38.0+84.0*xi+84.0*eta+t5-252.0*t6-t9+t11-t13;
  values[1] = -52.0+t15+t16+t17-t18-t19;
  values[2] = 21.0-t21-t22-t5+t23+t9;
  values[3] = 40.0-72.0*xi-72.0*eta-t27+336.0*t6+t29-t17+t18;
  values[4] = 21.0-t21-t22-t11+t13+t23;
  values[5] = -52.0+t15+t16+t27-t19-t29;
  values[6] = 60.0-120.0*xi-120.0*eta;
  values[7] = -105.0+t36+t37+420.0*t4-t39-210.0*t8+210.0*t10-105.0*t12;
  values[8] = -105.0+t36+t37+210.0*t4-t39-105.0*t8+420.0*t10-210.0*t12;
}

// values of the derivatives in eta-eta direction
static void C_T_UL2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t2, t3, t5, t6, t7, t8, t10, t11, t12, t14, t17, t22;

  t2 = xi*xi;
  t3 = 126.0*t2;
  t5 = xi*(1.0-xi-eta);
  t6 = 126.0*t5;
  t7 = xi*eta;
  t8 = 252.0*t7;
  t10 = 96.0*xi;
  t11 = 168.0*t5;
  t12 = 336.0*t7;
  t14 = 42.0*xi;
  t17 = 168.0*t2;
  t22 = 210.0*xi;

  values[0] = 4.0+84.0*xi-t3+t6-t8;
  values[1] = t10+t11-t12;
  values[2] = -t14+t3;
  values[3] = -72.0*xi+t17-t11+t12;
  values[4] = 4.0-t14-t6+t8;
  values[5] = -8.0+t10-t17;
  values[6] = -120.0*xi;
  values[7] = t22-420.0*t2+210.0*t5-420.0*t7;
  values[8] = t22-210.0*t2+420.0*t5-840.0*t7;
}

// ***********************************************************************

TBaseFunct2D *BF_C_T_UL2_2D_Obj = new TBaseFunct2D
        (9, BF_C_T_UL2_2D, BFUnitTriangle, 
         C_T_UL2_2D_Funct, C_T_UL2_2D_DeriveXi,
         C_T_UL2_2D_DeriveEta, C_T_UL2_2D_DeriveXiXi,
         C_T_UL2_2D_DeriveXiEta, C_T_UL2_2D_DeriveEtaEta, 4, 2,
         0, NULL);
