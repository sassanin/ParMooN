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
static void C_Q_UL3S_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t4, t5, t7, t8, t10, t11, t12, t14, t16, t17, t19, t21, t22;
  double t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t37;
  double t38, t39, t41, t42, t45, t46, t47, t48, t49, t50, t51, t52, t53, t55;
  double t59, t60, t61, t62, t63, t66, t67, t68, t72, t73, t74;

  t1 = xi/32.0;
  t2 = eta/32.0;
  t4 = xi*eta/4.0;
  t5 = xi*xi;
  t7 = 9.0/32.0*t5*eta;
  t8 = eta*eta;
  t10 = 9.0/32.0*xi*t8;
  t11 = t5*t8;
  t12 = 9.0/32.0*t11;
  t14 = 1.0-t5;
  t16 = (1.0-eta)*t14*xi;
  t17 = 9.0/32.0*t16;
  t19 = 1.0-t8;
  t21 = (1.0-xi)*t19*eta;
  t22 = 9.0/32.0*t21;
  t23 = t14*t19;
  t24 = t23*xi;
  t25 = 27.0/64.0*t24;
  t26 = t23*eta;
  t27 = 27.0/64.0*t26;
  t28 = -1.0/32.0+t1+t2+t4-t7-t10+t12+t17+t22-t25-t27;
  t29 = 9.0/128.0*t5;
  t30 = 9.0/32.0*eta;
  t31 = 27.0/128.0*t8;
  t32 = 27.0/128.0*t11;
  t33 = 27.0/32.0*t16;
  t34 = 81.0/64.0*t24;
  t35 = 81.0/128.0*t26;
  t37 = 9.0/64.0*t5;
  t38 = 27.0/64.0*t8;
  t39 = 27.0/64.0*t11;
  t41 = 9.0/128.0*t8;
  t42 = 45.0/128.0*t11;
  t45 = (1.0+xi)*t19*eta;
  t46 = 9.0/32.0*t45;
  t47 = 27.0/128.0*t26;
  t48 = 5.0/128.0-t1-t29+t2-t4-t7-t41+t10+t42-t17+t46+t25-t47;
  t49 = 9.0/32.0*xi;
  t50 = 27.0/64.0*t5;
  t51 = 9.0/64.0*t8;
  t52 = 27.0/32.0*t45;
  t53 = 81.0/128.0*t24;
  t55 = 27.0/128.0*t5;
  t59 = (1.0+eta)*t14*xi;
  t60 = 9.0/32.0*t59;
  t61 = 27.0/128.0*t24;
  t62 = -1.0/32.0-t1-t2+t4+t7+t10+t12-t60-t46+t61+t47;
  t63 = 27.0/32.0*t59;
  t66 = 5.0/128.0+t1-t29-t2-t4+t7-t41-t10+t42+t60-t22-t61+t27;
  t67 = 27.0/32.0*t21;
  t68 = 81.0/64.0*t26;
  t72 = 81.0/128.0*t5;
  t73 = 81.0/128.0*t8;
  t74 = 81.0/128.0*t11;

  values[0] = t28;
  values[1] = 9.0/128.0-t29-t30+t7+t31-t32-t33+t34+t35;
  values[2] = -9.0/64.0+t37-t30+t7+t38-t39+t33-t34;
  values[3] = t48;
  values[4] = -9.0/64.0+t49+t50+t51-t10-t39-t52-t53+t35;
  values[5] = 9.0/128.0+t49+t55-t41-t10-t32+t52-t35;
  values[6] = t62;
  values[7] = 9.0/128.0-t29+t30-t7+t31-t32+t63-t53;
  values[8] = -9.0/64.0+t37+t30-t7+t38-t39-t63+t53-t35;
  values[9] = t66;
  values[10] = -9.0/64.0-t49+t50+t51+t10-t39+t67-t68;
  values[11] = 9.0/128.0-t49+t55-t41+t10-t32-t67+t53+t68;
  values[12] = -243.0/128.0*t24-243.0/128.0*t26;
  values[13] = 81.0/128.0-t72-t73+t74+243.0/128.0*t24;
  values[14] = 81.0/128.0-t72-t73+t74+243.0/128.0*t26;
}

// values of the derivatives in xi direction
static void C_Q_UL3S_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t21, t23, t24, t25, t26, t27, t28, t29, t30, t31;
  double t32, t34, t35, t37, t38, t39, t40, t41, t42, t43, t45, t47, t48, t49;
  double t50, t51, t52, t53, t54, t55, t56, t59, t60, t63, t64, t65, t67, t68;

  t1 = eta/4.0;
  t3 = 9.0/16.0*xi*eta;
  t4 = eta*eta;
  t5 = 9.0/32.0*t4;
  t6 = xi*t4;
  t7 = 9.0/16.0*t6;
  t8 = 1.0-eta;
  t9 = xi*xi;
  t10 = t8*t9;
  t11 = 9.0/16.0*t10;
  t12 = 1.0-t9;
  t13 = t8*t12;
  t14 = 9.0/32.0*t13;
  t15 = 1.0-t4;
  t16 = t15*eta;
  t17 = 9.0/32.0*t16;
  t18 = t9*t15;
  t19 = 27.0/32.0*t18;
  t20 = t12*t15;
  t21 = 27.0/64.0*t20;
  t23 = xi*t15*eta;
  t24 = 27.0/32.0*t23;
  t25 = 1.0/32.0+t1-t3-t5+t7-t11+t14-t17+t19-t21+t24;
  t26 = 9.0/64.0*xi;
  t27 = 27.0/64.0*t6;
  t28 = 27.0/16.0*t10;
  t29 = 27.0/32.0*t13;
  t30 = 81.0/32.0*t18;
  t31 = 81.0/64.0*t20;
  t32 = 81.0/64.0*t23;
  t34 = 9.0/32.0*xi;
  t35 = 27.0/32.0*t6;
  t37 = 45.0/64.0*t6;
  t38 = 27.0/64.0*t23;
  t39 = -1.0/32.0-t26-t1-t3+t5+t37+t11-t14+t17-t19+t21+t38;
  t40 = 27.0/32.0*xi;
  t41 = 27.0/32.0*t16;
  t42 = 81.0/64.0*t18;
  t43 = 81.0/128.0*t20;
  t45 = 27.0/64.0*xi;
  t47 = 1.0+eta;
  t48 = t47*t9;
  t49 = 9.0/16.0*t48;
  t50 = t47*t12;
  t51 = 9.0/32.0*t50;
  t52 = 27.0/64.0*t18;
  t53 = 27.0/128.0*t20;
  t54 = -1.0/32.0+t1+t3+t5+t7+t49-t51-t17-t52+t53-t38;
  t55 = 27.0/16.0*t48;
  t56 = 27.0/32.0*t50;
  t59 = 1.0/32.0-t26-t1+t3-t5+t37-t49+t51+t17+t52-t53-t24;
  t60 = 81.0/32.0*t23;
  t63 = 243.0/64.0*t18;
  t64 = 243.0/128.0*t20;
  t65 = 243.0/64.0*t23;
  t67 = 81.0/64.0*xi;
  t68 = 81.0/64.0*t6;

  values[0] = t25;
  values[1] = -t26+t3-t27+t28-t29-t30+t31-t32;
  values[2] = t34+t3-t35-t28+t29+t30-t31;
  values[3] = t39;
  values[4] = 9.0/32.0+t40-t5-t35-t41+t42-t43-t32;
  values[5] = 9.0/32.0+t45-t5-t27+t41+t32;
  values[6] = t54;
  values[7] = -t26-t3-t27-t55+t56+t42-t43;
  values[8] = t34-t3-t35+t55-t56-t42+t43+t32;
  values[9] = t59;
  values[10] = -9.0/32.0+t40+t5-t35-t41+t60;
  values[11] = -9.0/32.0+t45+t5-t27+t41-t42+t43-t60;
  values[12] = t63-t64+t65;
  values[13] = -t67+t68-t63+t64;
  values[14] = -t67+t68-t65;
}

// values of the derivatives in eta direction
static void C_Q_UL3S_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31;
  double t33, t34, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47;
  double t48, t49, t52, t53, t56, t57, t58, t59, t60, t63, t64, t65, t67, t68;

  t1 = xi/4.0;
  t2 = xi*xi;
  t3 = 9.0/32.0*t2;
  t5 = 9.0/16.0*xi*eta;
  t6 = t2*eta;
  t7 = 9.0/16.0*t6;
  t8 = 1.0-t2;
  t9 = xi*t8;
  t10 = 9.0/32.0*t9;
  t11 = 1.0-xi;
  t12 = eta*eta;
  t13 = t11*t12;
  t14 = 9.0/16.0*t13;
  t15 = 1.0-t12;
  t16 = t11*t15;
  t17 = 9.0/32.0*t16;
  t19 = t8*eta*xi;
  t20 = 27.0/32.0*t19;
  t21 = t8*t12;
  t22 = 27.0/32.0*t21;
  t23 = t8*t15;
  t24 = 27.0/64.0*t23;
  t25 = 1.0/32.0+t1-t3-t5+t7-t10-t14+t17+t20+t22-t24;
  t26 = 27.0/64.0*eta;
  t27 = 27.0/64.0*t6;
  t28 = 27.0/32.0*t9;
  t29 = 81.0/32.0*t19;
  t30 = 81.0/64.0*t21;
  t31 = 81.0/128.0*t23;
  t33 = 27.0/32.0*eta;
  t34 = 27.0/32.0*t6;
  t36 = 9.0/64.0*eta;
  t37 = 45.0/64.0*t6;
  t38 = 1.0+xi;
  t39 = t38*t12;
  t40 = 9.0/16.0*t39;
  t41 = t38*t15;
  t42 = 9.0/32.0*t41;
  t43 = 27.0/64.0*t21;
  t44 = 27.0/128.0*t23;
  t45 = 1.0/32.0-t1-t3-t36+t5+t37+t10-t40+t42-t20+t43-t44;
  t46 = 9.0/32.0*eta;
  t47 = 27.0/16.0*t39;
  t48 = 27.0/32.0*t41;
  t49 = 81.0/64.0*t19;
  t52 = 27.0/64.0*t19;
  t53 = -1.0/32.0+t1+t3+t5+t7-t10+t40-t42-t52-t43+t44;
  t56 = -1.0/32.0-t1+t3-t36-t5+t37+t10+t14-t17+t52-t22+t24;
  t57 = 27.0/16.0*t13;
  t58 = 27.0/32.0*t16;
  t59 = 81.0/32.0*t21;
  t60 = 81.0/64.0*t23;
  t63 = 243.0/64.0*t19;
  t64 = 243.0/64.0*t21;
  t65 = 243.0/128.0*t23;
  t67 = 81.0/64.0*eta;
  t68 = 81.0/64.0*t6;

  values[0] = t25;
  values[1] = -9.0/32.0+t3+t26-t27+t28-t29-t30+t31;
  values[2] = -9.0/32.0+t3+t33-t34-t28+t29;
  values[3] = t45;
  values[4] = t46-t5-t34+t47-t48+t49-t30+t31;
  values[5] = -t36-t5-t27-t47+t48+t30-t31;
  values[6] = t53;
  values[7] = 9.0/32.0-t3+t26-t27+t28+t49;
  values[8] = 9.0/32.0-t3+t33-t34-t28-t49+t30-t31;
  values[9] = t56;
  values[10] = t46+t5-t34-t57+t58+t59-t60;
  values[11] = -t36+t5-t27+t57-t58-t49-t59+t60;
  values[12] = t63+t64-t65;
  values[13] = -t67+t68-t63;
  values[14] = -t67+t68-t64+t65;
}

// values of the derivatives in xi-xi  direction
static void C_Q_UL3S_2D_DeriveXiXi(double xi, double eta, double *values)
{
 double t1, t2, t3, t5, t6, t7, t8, t9, t10, t11, t13, t14, t15, t16, t18;
 double t20, t21, t23, t27, t28, t29, t31, t35, t38, t39, t41;

  t1 = 9.0/16.0*eta;
  t2 = eta*eta;
  t3 = 9.0/16.0*t2;
  t5 = (1.0-eta)*xi;
  t6 = 27.0/16.0*t5;
  t7 = 1.0-t2;
  t8 = xi*t7;
  t9 = 81.0/32.0*t8;
  t10 = t7*eta;
  t11 = 27.0/32.0*t10;
  t13 = 27.0/64.0*t2;
  t14 = 81.0/16.0*t5;
  t15 = 243.0/32.0*t8;
  t16 = 81.0/64.0*t10;
  t18 = 27.0/32.0*t2;
  t20 = 45.0/64.0*t2;
  t21 = 27.0/64.0*t10;
  t23 = 243.0/64.0*t8;
  t27 = (1.0+eta)*xi;
  t28 = 27.0/16.0*t27;
  t29 = 81.0/64.0*t8;
  t31 = 81.0/16.0*t27;
  t35 = 81.0/32.0*t10;
  t38 = 729.0/64.0*t8;
  t39 = 243.0/64.0*t10;
  t41 = 81.0/64.0*t2;

  values[0] = -t1+t3-t6+t9+t11;
  values[1] = -9.0/64.0+t1-t13+t14-t15-t16;
  values[2] = 9.0/32.0+t1-t18-t14+t15;
  values[3] = -9.0/64.0-t1+t20+t6-t9+t21;
  values[4] = 27.0/32.0-t18+t23-t16;
  values[5] = 27.0/64.0-t13+t16;
  values[6] = t1+t3+t28-t29-t21;
  values[7] = -9.0/64.0-t1-t13-t31+t23;
  values[8] = 9.0/32.0-t1-t18+t31-t23+t16;
  values[9] = -9.0/64.0+t1+t20-t28+t29-t11;
  values[10] = 27.0/32.0-t18+t35;
  values[11] = 27.0/64.0-t13-t23-t35;
  values[12] = t38+t39;
  values[13] = -81.0/64.0+t41-t38;
  values[14] = -81.0/64.0+t41-t39;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL3S_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t12, t13, t14, t15, t17;
  double t18, t20, t21, t22, t23, t24, t25, t27, t29, t30, t31, t33, t34, t35;
  double t38, t39, t44, t45, t48, t49, t50, t51, t53;

   t1 = 9.0/16.0*xi;
   t2 = 9.0/16.0*eta;
   t3 = xi*eta;
   t4 = 9.0/8.0*t3;
   t5 = xi*xi;
   t6 = 27.0/32.0*t5;
   t7 = eta*eta;
   t8 = 27.0/32.0*t7;
   t9 = t5*eta;
   t10 = 27.0/16.0*t9;
   t12 = (1.0-t5)*eta;
   t13 = 27.0/32.0*t12;
   t14 = xi*t7;
   t15 = 27.0/16.0*t14;
   t17 = xi*(1.0-t7);
   t18 = 27.0/32.0*t17;
   t20 = 27.0/32.0*t3;
   t21 = 81.0/32.0*t5;
   t22 = 81.0/16.0*t9;
   t23 = 81.0/32.0*t12;
   t24 = 81.0/32.0*t14;
   t25 = 81.0/64.0*t17;
   t27 = 27.0/16.0*t3;
   t29 = 45.0/32.0*t3;
   t30 = 27.0/32.0*t14;
   t31 = 27.0/64.0*t17;
   t33 = 81.0/32.0*t7;
   t34 = 81.0/32.0*t9;
   t35 = 81.0/64.0*t12;
   t38 = 27.0/32.0*t9;
   t39 = 27.0/64.0*t12;
   t44 = 81.0/16.0*t14;
   t45 = 81.0/32.0*t17;
   t48 = 243.0/32.0*t9;
   t49 = 243.0/64.0*t12;
   t50 = 243.0/32.0*t14;
   t51 = 243.0/64.0*t17;
   t53 = 81.0/32.0*t3;

  values[0] = -5.0/16.0-t1-t2+t4+t6+t8-t10+t13-t15+t18;
  values[1] = t1-t20-t21+27.0/32.0+t22-t23+t24-t25;
  values[2] = t1-t27+t21-27.0/32.0-t22+t23;
  values[3] = 5.0/16.0-t1+t2+t29-t6-t8+t10-t13-t30+t31;
  values[4] = -t2-t27+t33-27.0/32.0-t34+t35+t24-t25;
  values[5] = -t2-t20-t33+27.0/32.0-t24+t25;
  values[6] = -5.0/16.0+t1+t2+t4+t6+t8+t38-t39+t30-t31;
  values[7] = -t1-t20-t21+27.0/32.0-t34+t35;
  values[8] = -t1-t27+t21-27.0/32.0+t34-t35-t24+t25;
  values[9] = 5.0/16.0+t1-t2+t29-t6-t8-t38+t39+t15-t18;
  values[10] = t2-t27+t33-27.0/32.0-t44+t45;
  values[11] = t2-t20-t33+27.0/32.0+t34-t35+t44-t45;
  values[12] = -t48+t49-t50+t51;
  values[13] = t53+t48-t49;
  values[14] = t53+t50-t51;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL3S_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t5, t6, t7, t8, t9, t10, t11, t13, t14, t15, t17, t19;
  double t21, t22, t23, t25, t26, t29, t34, t35, t38, t39, t41;

  t1 = 9.0/16.0*xi;
  t2 = xi*xi;
  t3 = 9.0/16.0*t2;
  t5 = (1.0-xi)*eta;
  t6 = 27.0/16.0*t5;
  t7 = 1.0-t2;
  t8 = xi*t7;
  t9 = 27.0/32.0*t8;
  t10 = t7*eta;
  t11 = 81.0/32.0*t10;
  t13 = 27.0/64.0*t2;
  t14 = 81.0/32.0*t8;
  t15 = 243.0/64.0*t10;
  t17 = 27.0/32.0*t2;
  t19 = 45.0/64.0*t2;
  t21 = (1.0+xi)*eta;
  t22 = 27.0/16.0*t21;
  t23 = 81.0/64.0*t10;
  t25 = 81.0/16.0*t21;
  t26 = 81.0/64.0*t8;
  t29 = 27.0/64.0*t8;
  t34 = 81.0/16.0*t5;
  t35 = 243.0/32.0*t10;
  t38 = 243.0/64.0*t8;
  t39 = 729.0/64.0*t10;
  t41 = 81.0/64.0*t2;

  values[0] = -t1+t3-t6+t9+t11;
  values[1] = 27.0/64.0-t13-t14-t15;
  values[2] = 27.0/32.0-t17+t14;
  values[3] = -9.0/64.0+t1+t19-t22-t9+t23;
  values[4] = 9.0/32.0-t1-t17+t25+t26-t15;
  values[5] = -9.0/64.0-t1-t13-t25+t15;
  values[6] = t1+t3+t22-t29-t23;
  values[7] = 27.0/64.0-t13+t26;
  values[8] = 27.0/32.0-t17-t26+t15;
  values[9] = -9.0/64.0-t1+t19+t6+t29-t11;
  values[10] = 9.0/32.0+t1-t17-t34+t35;
  values[11] = -9.0/64.0+t1-t13+t34-t26-t35;
  values[12] = t38+t39;
  values[13] = -81.0/64.0+t41-t38;
  values[14] = -81.0/64.0+t41-t39;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_UL3S_2D_Obj = new TBaseFunct2D
        (15, BF_C_Q_UL3S_2D, BFUnitSquare, 
         C_Q_UL3S_2D_Funct, C_Q_UL3S_2D_DeriveXi,
         C_Q_UL3S_2D_DeriveEta, C_Q_UL3S_2D_DeriveXiXi,
         C_Q_UL3S_2D_DeriveXiEta, C_Q_UL3S_2D_DeriveEtaEta, 3, 3,
         0, NULL);
