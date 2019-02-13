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
static void C_Q_UL3_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54;
  double t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67;
  double t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t91;

  t1 = xi/8.0;
  t2 = eta/8.0;
  t3 = xi*xi;
  t4 = eta*eta;
  t5 = t3*t4;
  t6 = 81.0/64.0*t5;
  t7 = t4*eta;
  t8 = 3.0/32.0*t7;
  t9 = t3*xi;
  t10 = 3.0/32.0*t9;
  t11 = t9*t4;
  t12 = 3.0/16.0*t11;
  t13 = xi*t4;
  t14 = 3.0/32.0*t13;
  t15 = t3*eta;
  t16 = 3.0/32.0*t15;
  t17 = t3*t7;
  t18 = 3.0/16.0*t17;
  t19 = xi*eta;
  t20 = 15.0/32.0*t19;
  t21 = t3*t3;
  t22 = 105.0/256.0*t21;
  t23 = 147.0/256.0*t3;
  t24 = 147.0/256.0*t4;
  t25 = t4*t4;
  t26 = t3*t25;
  t27 = 105.0/256.0*t26;
  t28 = 105.0/256.0*t25;
  t29 = xi*t7;
  t30 = 7.0/16.0*t29;
  t31 = t9*eta;
  t32 = 7.0/16.0*t31;
  t33 = t9*t7;
  t34 = 5.0/32.0*t33;
  t35 = t21*t4;
  t36 = 105.0/256.0*t35;
  t37 = t1+t2+t6-t8-t10-t12-t14-t16-t18-t20+t22-t23-t24-t27+t28+t30+t32-t34-t36+17.0/128.0;
  t38 = 27.0/64.0*xi;
  t39 = 27.0/64.0*eta;
  t40 = 135.0/128.0*t5;
  t41 = 45.0/64.0*t7;
  t42 = 27.0/64.0*t9;
  t43 = 81.0/64.0*t11;
  t44 = 81.0/64.0*t13;
  t45 = 27.0/64.0*t15;
  t46 = 45.0/64.0*t17;
  t47 = 81.0/64.0*t19;
  t48 = 27.0/256.0*t3;
  t49 = 135.0/128.0*t4;
  t50 = 315.0/256.0*t26;
  t51 = 315.0/256.0*t25;
  t52 = 135.0/64.0*t29;
  t53 = 81.0/64.0*t31;
  t54 = 135.0/64.0*t33;
  t55 = t38+t39+t40-t41-t42+t43-t44-t45+t46-t47-t48-t49-t50+t51+t52+t53-t54+27.0/256.0;
  t56 = -t38+t39+t40-t41+t42-t43+t44-t45+t46+t47-t48-t49-t50+t51-t52-t53+t54+27.0/256.0;
  t57 = -t1+t2+t6-t8+t10+t12+t14-t16-t18+t20+t22-t23-t24-t27+t28-t30-t32+t34-t36+17.0/128.0;
  t58 = 27.0/64.0*t7;
  t59 = 45.0/64.0*t9;
  t60 = 45.0/64.0*t11;
  t61 = 27.0/64.0*t13;
  t62 = 81.0/64.0*t15;
  t63 = 81.0/64.0*t17;
  t64 = 315.0/256.0*t21;
  t65 = 135.0/128.0*t3;
  t66 = 27.0/256.0*t4;
  t67 = 81.0/64.0*t29;
  t68 = 135.0/64.0*t31;
  t69 = 315.0/256.0*t35;
  t70 = -t38+t39+t40-t58+t59-t60+t61-t62+t63+t47+t64-t65-t66-t67-t68+t54-t69+27.0/256.0;
  t71 = -t38-t39+t40+t58+t59-t60+t61+t62-t63-t47+t64-t65-t66+t67+t68-t54-t69+27.0/256.0;
  t72 = -t1-t2+t6+t8+t10+t12+t14+t16+t18-t20+t22-t23-t24-t27+t28+t30+t32-t34-t36+17.0/128.0;
  t73 = -t38-t39+t40+t41+t42-t43+t44+t45-t46-t47-t48-t49-t50+t51+t52+t53-t54+27.0/256.0;
  t74 = t38-t39+t40+t41-t42+t43-t44+t45-t46+t47-t48-t49-t50+t51-t52-t53+t54+27.0/256.0;
  t75 = t1-t2+t6+t8-t10-t12-t14+t16+t18+t20+t22-t23-t24-t27+t28-t30-t32+t34-t36+17.0/128.0;
  t76 = t38-t39+t40+t58-t59+t60-t61+t62-t63+t47+t64-t65-t66-t67-t68+t54-t69+27.0/256.0;
  t77 = t38+t39+t40-t58-t59+t60-t61-t62+t63-t47+t64-t65-t66+t67+t68-t54-t69+27.0/256.0;
  t91 = 315.0/32.0*t5;

  values[0] = t37;
  values[1] = t55;
  values[2] = t56;
  values[3] = t57;
  values[4] = t70;
  values[5] = t71;
  values[6] = t72;
  values[7] = t73;
  values[8] = t74;
  values[9] = t75;
  values[10] = t76;
  values[11] = t77;
  values[12] = -3.0/32.0+111.0/64.0*t4+111.0/64.0*t3-27.0/8.0*t5-105.0/64.0*t21+105.0/64.0*t35-105.0/64.0*t25+105.0/64.0*t26;
  values[13] = 45.0/16.0*eta-45.0/16.0*t7-45.0/16.0*t15+45.0/16.0*t17;
  values[14] = 45.0/16.0*xi-45.0/16.0*t13-45.0/16.0*t9+45.0/16.0*t11;
  values[15] = 225.0/16.0*t19-225.0/16.0*t29-225.0/16.0*t31+225.0/16.0*t33;
  values[16] = -105.0/64.0+315.0/32.0*t4+105.0/64.0*t3-t91-525.0/64.0*t25+525.0/64.0*t26;
  values[17] = -105.0/64.0+105.0/64.0*t4+315.0/32.0*t3-t91-525.0/64.0*t21+525.0/64.0*t35;
}

// values of the derivatives in xi direction
static void C_Q_UL3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54;
  double t55, t56, t57, t58, t59, t60, t61, t79;

  t1 = eta*eta;
  t2 = xi*t1;
  t3 = 81.0/32.0*t2;
  t4 = xi*xi;
  t5 = 9.0/32.0*t4;
  t6 = t4*t1;
  t7 = 9.0/16.0*t6;
  t8 = 3.0/32.0*t1;
  t9 = xi*eta;
  t10 = 3.0/16.0*t9;
  t11 = t1*eta;
  t12 = xi*t11;
  t13 = 3.0/8.0*t12;
  t14 = 15.0/32.0*eta;
  t15 = xi*t4;
  t16 = 105.0/64.0*t15;
  t17 = 147.0/128.0*xi;
  t18 = t1*t1;
  t19 = xi*t18;
  t20 = 105.0/128.0*t19;
  t21 = 7.0/16.0*t11;
  t22 = t4*eta;
  t23 = 21.0/16.0*t22;
  t24 = t4*t11;
  t25 = 15.0/32.0*t24;
  t26 = t15*t1;
  t27 = 105.0/64.0*t26;
  t28 = 1.0/8.0+t3-t5-t7-t8-t10-t13-t14+t16-t17-t20+t21+t23-t25-t27;
  t29 = 135.0/64.0*t2;
  t30 = 81.0/64.0*t4;
  t31 = 243.0/64.0*t6;
  t32 = 81.0/64.0*t1;
  t33 = 27.0/32.0*t9;
  t34 = 45.0/32.0*t12;
  t35 = 81.0/64.0*eta;
  t36 = 27.0/128.0*xi;
  t37 = 315.0/128.0*t19;
  t38 = 135.0/64.0*t11;
  t39 = 243.0/64.0*t22;
  t40 = 405.0/64.0*t24;
  t41 = 27.0/64.0+t29-t30+t31-t32-t33+t34-t35-t36-t37+t38+t39-t40;
  t42 = -27.0/64.0+t29+t30-t31+t32-t33+t34+t35-t36-t37-t38-t39+t40;
  t43 = -1.0/8.0+t3+t5+t7+t8-t10-t13+t14+t16-t17-t20-t21-t23+t25-t27;
  t44 = 135.0/64.0*t4;
  t45 = 135.0/64.0*t6;
  t46 = 27.0/64.0*t1;
  t47 = 81.0/32.0*t9;
  t48 = 81.0/32.0*t12;
  t49 = 315.0/64.0*t15;
  t50 = 135.0/64.0*xi;
  t51 = 81.0/64.0*t11;
  t52 = 405.0/64.0*t22;
  t53 = 315.0/64.0*t26;
  t54 = -27.0/64.0+t29+t44-t45+t46-t47+t48+t35+t49-t50-t51-t52+t40-t53;
  t55 = -27.0/64.0+t29+t44-t45+t46+t47-t48-t35+t49-t50+t51+t52-t40-t53;
  t56 = -1.0/8.0+t3+t5+t7+t8+t10+t13-t14+t16-t17-t20+t21+t23-t25-t27;
  t57 = -27.0/64.0+t29+t30-t31+t32+t33-t34-t35-t36-t37+t38+t39-t40;
  t58 = 27.0/64.0+t29-t30+t31-t32+t33-t34+t35-t36-t37-t38-t39+t40;
  t59 = 1.0/8.0+t3-t5-t7-t8+t10+t13+t14+t16-t17-t20-t21-t23+t25-t27;
  t60 = 27.0/64.0+t29-t44+t45-t46+t47-t48+t35+t49-t50-t51-t52+t40-t53;
  t61 = 27.0/64.0+t29-t44+t45-t46-t47+t48-t35+t49-t50+t51+t52-t40-t53;
  t79 = 315.0/16.0*t2;

  values[0] = t28;
  values[1] = t41;
  values[2] = t42;
  values[3] = t43;
  values[4] = t54;
  values[5] = t55;
  values[6] = t56;
  values[7] = t57;
  values[8] = t58;
  values[9] = t59;
  values[10] = t60;
  values[11] = t61;
  values[12] = 111.0/32.0*xi-27.0/4.0*t2-105.0/16.0*t15+105.0/16.0*t26+105.0/32.0*t19;
  values[13] = -45.0/8.0*t9+45.0/8.0*t12;
  values[14] = 45.0/16.0-45.0/16.0*t1-135.0/16.0*t4+135.0/16.0*t6;
  values[15] = 225.0/16.0*eta-225.0/16.0*t11-675.0/16.0*t22+675.0/16.0*t24;
  values[16] = 105.0/32.0*xi-t79+525.0/32.0*t19;
  values[17] = 315.0/16.0*xi-t79-525.0/16.0*t15+525.0/16.0*t26;
}

// values of the derivatives in eta direction
static void C_Q_UL3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54;
  double t55, t56, t57, t58, t59, t60, t61, t79;

  t1 = xi*xi;
  t2 = t1*eta;
  t3 = 81.0/32.0*t2;
  t4 = eta*eta;
  t5 = 9.0/32.0*t4;
  t6 = t1*xi;
  t7 = t6*eta;
  t8 = 3.0/8.0*t7;
  t9 = xi*eta;
  t10 = 3.0/16.0*t9;
  t11 = 3.0/32.0*t1;
  t12 = t4*t1;
  t13 = 9.0/16.0*t12;
  t14 = 15.0/32.0*xi;
  t15 = 147.0/128.0*eta;
  t16 = t4*eta;
  t17 = t1*t16;
  t18 = 105.0/64.0*t17;
  t19 = 105.0/64.0*t16;
  t20 = xi*t4;
  t21 = 21.0/16.0*t20;
  t22 = 7.0/16.0*t6;
  t23 = t6*t4;
  t24 = 15.0/32.0*t23;
  t25 = t1*t1;
  t26 = t25*eta;
  t27 = 105.0/128.0*t26;
  t28 = 1.0/8.0+t3-t5-t8-t10-t11-t13-t14-t15-t18+t19+t21+t22-t24-t27;
  t29 = 135.0/64.0*t2;
  t30 = 135.0/64.0*t4;
  t31 = 81.0/32.0*t7;
  t32 = 81.0/32.0*t9;
  t33 = 27.0/64.0*t1;
  t34 = 135.0/64.0*t12;
  t35 = 81.0/64.0*xi;
  t36 = 135.0/64.0*eta;
  t37 = 315.0/64.0*t17;
  t38 = 315.0/64.0*t16;
  t39 = 405.0/64.0*t20;
  t40 = 81.0/64.0*t6;
  t41 = 405.0/64.0*t23;
  t42 = 27.0/64.0+t29-t30+t31-t32-t33+t34-t35-t36-t37+t38+t39+t40-t41;
  t43 = 27.0/64.0+t29-t30-t31+t32-t33+t34+t35-t36-t37+t38-t39-t40+t41;
  t44 = 1.0/8.0+t3-t5+t8+t10-t11-t13+t14-t15-t18+t19-t21-t22+t24-t27;
  t45 = 81.0/64.0*t4;
  t46 = 45.0/32.0*t7;
  t47 = 27.0/32.0*t9;
  t48 = 81.0/64.0*t1;
  t49 = 243.0/64.0*t12;
  t50 = 27.0/128.0*eta;
  t51 = 243.0/64.0*t20;
  t52 = 135.0/64.0*t6;
  t53 = 315.0/128.0*t26;
  t54 = 27.0/64.0+t29-t45-t46+t47-t48+t49+t35-t50-t51-t52+t41-t53;
  t55 = -27.0/64.0+t29+t45-t46+t47+t48-t49-t35-t50+t51+t52-t41-t53;
  t56 = -1.0/8.0+t3+t5+t8+t10+t11+t13-t14-t15-t18+t19+t21+t22-t24-t27;
  t57 = -27.0/64.0+t29+t30-t31+t32+t33-t34-t35-t36-t37+t38+t39+t40-t41;
  t58 = -27.0/64.0+t29+t30+t31-t32+t33-t34+t35-t36-t37+t38-t39-t40+t41;
  t59 = -1.0/8.0+t3+t5-t8-t10+t11+t13+t14-t15-t18+t19-t21-t22+t24-t27;
  t60 = -27.0/64.0+t29+t45+t46-t47+t48-t49+t35-t50-t51-t52+t41-t53;
  t61 = 27.0/64.0+t29-t45+t46-t47-t48+t49-t35-t50+t51+t52-t41-t53;
  t79 = 315.0/16.0*t2;

  values[0] = t28;
  values[1] = t42;
  values[2] = t43;
  values[3] = t44;
  values[4] = t54;
  values[5] = t55;
  values[6] = t56;
  values[7] = t57;
  values[8] = t58;
  values[9] = t59;
  values[10] = t60;
  values[11] = t61;
  values[12] = 111.0/32.0*eta-27.0/4.0*t2+105.0/32.0*t26-105.0/16.0*t16+105.0/16.0*t17;
  values[13] = 45.0/16.0-135.0/16.0*t4-45.0/16.0*t1+135.0/16.0*t12;
  values[14] = -45.0/8.0*t9+45.0/8.0*t7;
  values[15] = 225.0/16.0*xi-675.0/16.0*t20-225.0/16.0*t6+675.0/16.0*t23;
  values[16] = 315.0/16.0*eta-t79-525.0/16.0*t16+525.0/16.0*t17;
  values[17] = 105.0/32.0*eta-t79+525.0/32.0*t26;
}

// values of the derivatives in xi-xi  direction
static void C_Q_UL3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t30;
  double t31, t32, t33, t34, t35, t36, t37, t40, t43, t54;

  t1 = eta*eta;
  t2 = 81.0/32.0*t1;
  t3 = 9.0/16.0*xi;
  t4 = t1*xi;
  t5 = 9.0/8.0*t4;
  t6 = 3.0/16.0*eta;
  t7 = t1*eta;
  t8 = 3.0/8.0*t7;
  t9 = xi*xi;
  t10 = 315.0/64.0*t9;
  t11 = t1*t1;
  t12 = 105.0/128.0*t11;
  t13 = xi*eta;
  t14 = 21.0/8.0*t13;
  t15 = xi*t7;
  t16 = 15.0/16.0*t15;
  t17 = t9*t1;
  t18 = 315.0/64.0*t17;
  t19 = t2-t3-t5-t6-t8+t10-147.0/128.0-t12+t14-t16-t18;
  t20 = 135.0/64.0*t1;
  t21 = 81.0/32.0*xi;
  t22 = 243.0/32.0*t4;
  t23 = 27.0/32.0*eta;
  t24 = 45.0/32.0*t7;
  t25 = 315.0/128.0*t11;
  t26 = 243.0/32.0*t13;
  t27 = 405.0/32.0*t15;
  t30 = t2+t3+t5-t6-t8+t10-147.0/128.0-t12-t14+t16-t18;
  t31 = 135.0/32.0*xi;
  t32 = 135.0/32.0*t4;
  t33 = 81.0/32.0*eta;
  t34 = 81.0/32.0*t7;
  t35 = 945.0/64.0*t9;
  t36 = 405.0/32.0*t13;
  t37 = 945.0/64.0*t17;
  t40 = t2+t3+t5+t6+t8+t10-147.0/128.0-t12+t14-t16-t18;
  t43 = t2-t3-t5+t6+t8+t10-147.0/128.0-t12-t14+t16-t18;
  t54 = 315.0/16.0*t1;

  values[0] = t19;
  values[1] = t20-t21+t22-t23+t24-27.0/128.0-t25+t26-t27;
  values[2] = t20+t21-t22-t23+t24-27.0/128.0-t25-t26+t27;
  values[3] = t30;
  values[4] = t20+t31-t32-t33+t34+t35-135.0/64.0-t36+t27-t37;
  values[5] = t20+t31-t32+t33-t34+t35-135.0/64.0+t36-t27-t37;
  values[6] = t40;
  values[7] = t20+t21-t22+t23-t24-27.0/128.0-t25+t26-t27;
  values[8] = t20-t21+t22+t23-t24-27.0/128.0-t25-t26+t27;
  values[9] = t43;
  values[10] = t20-t31+t32+t33-t34+t35-135.0/64.0-t36+t27-t37;
  values[11] = t20-t31+t32-t33+t34+t35-135.0/64.0+t36-t27-t37;
  values[12] = 111.0/32.0-27.0/4.0*t1-315.0/16.0*t9+315.0/16.0*t17+105.0/32.0*t11;
  values[13] = -45.0/8.0*eta+45.0/8.0*t7;
  values[14] = -135.0/8.0*xi+135.0/8.0*t4;
  values[15] = -675.0/8.0*t13+675.0/8.0*t15;
  values[16] = 105.0/32.0-t54+525.0/32.0*t11;
  values[17] = 315.0/16.0-t54-1575.0/16.0*t9+1575.0/16.0*t17;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t12, t13, t14, t15, t16;
  double t17, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
  double t33, t34, t35, t36, t37, t38, t39, t40, t43, t46, t63;

  t1 = xi*eta;
  t2 = 81.0/16.0*t1;
  t3 = xi*xi;
  t4 = t3*eta;
  t5 = 9.0/8.0*t4;
  t6 = 3.0/16.0*eta;
  t7 = 3.0/16.0*xi;
  t8 = eta*eta;
  t9 = xi*t8;
  t10 = 9.0/8.0*t9;
  t12 = xi*t8*eta;
  t13 = 105.0/32.0*t12;
  t14 = 21.0/16.0*t8;
  t15 = 21.0/16.0*t3;
  t16 = t3*t8;
  t17 = 45.0/32.0*t16;
  t19 = t3*xi*eta;
  t20 = 105.0/32.0*t19;
  t21 = t2-t5-t6-t7-t10-15.0/32.0-t13+t14+t15-t17-t20;
  t22 = 135.0/32.0*t1;
  t23 = 243.0/32.0*t4;
  t24 = 81.0/32.0*eta;
  t25 = 27.0/32.0*xi;
  t26 = 135.0/32.0*t9;
  t27 = 315.0/32.0*t12;
  t28 = 405.0/64.0*t8;
  t29 = 243.0/64.0*t3;
  t30 = 1215.0/64.0*t16;
  t33 = t2+t5+t6-t7-t10+15.0/32.0-t13-t14-t15+t17-t20;
  t34 = 135.0/32.0*t4;
  t35 = 27.0/32.0*eta;
  t36 = 81.0/32.0*xi;
  t37 = 243.0/32.0*t9;
  t38 = 243.0/64.0*t8;
  t39 = 405.0/64.0*t3;
  t40 = 315.0/32.0*t19;
  t43 = t2+t5+t6+t7+t10-15.0/32.0-t13+t14+t15-t17-t20;
  t46 = t2-t5-t6+t7+t10+15.0/32.0-t13-t14-t15+t17-t20;
  t63 = 315.0/8.0*t1;

  values[0] = t21;
  values[1] = t22+t23-t24-t25+t26-81.0/64.0-t27+t28+t29-t30;
  values[2] = t22-t23+t24-t25+t26+81.0/64.0-t27-t28-t29+t30;
  values[3] = t33;
  values[4] = t22-t34+t35-t36+t37+81.0/64.0-t38-t39+t30-t40;
  values[5] = t22-t34+t35+t36-t37-81.0/64.0+t38+t39-t30-t40;
  values[6] = t43;
  values[7] = t22-t23+t24+t25-t26-81.0/64.0-t27+t28+t29-t30;
  values[8] = t22+t23-t24+t25-t26+81.0/64.0-t27-t28-t29+t30;
  values[9] = t46;
  values[10] = t22+t34-t35+t36-t37+81.0/64.0-t38-t39+t30-t40;
  values[11] = t22+t34-t35-t36+t37-81.0/64.0+t38+t39-t30-t40;
  values[12] = -27.0/2.0*t1+105.0/8.0*t19+105.0/8.0*t12;
  values[13] = -45.0/8.0*xi+135.0/8.0*t9;
  values[14] = -45.0/8.0*eta+135.0/8.0*t4;
  values[15] = 225.0/16.0-675.0/16.0*t8-675.0/16.0*t3+2025.0/16.0*t16;
  values[16] = -t63+525.0/8.0*t12;
  values[17] = -t63+525.0/8.0*t19;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t31, t32, t33, t34, t35, t36, t37, t40, t43, t54;

  t1 = xi*xi;
  t2 = 81.0/32.0*t1;
  t3 = 9.0/16.0*eta;
  t4 = t1*xi;
  t5 = 3.0/8.0*t4;
  t6 = 3.0/16.0*xi;
  t7 = t1*eta;
  t8 = 9.0/8.0*t7;
  t9 = eta*eta;
  t10 = t9*t1;
  t11 = 315.0/64.0*t10;
  t12 = 315.0/64.0*t9;
  t13 = xi*eta;
  t14 = 21.0/8.0*t13;
  t15 = t4*eta;
  t16 = 15.0/16.0*t15;
  t17 = t1*t1;
  t18 = 105.0/128.0*t17;
  t19 = t2-t3-t5-t6-t8-147.0/128.0-t11+t12+t14-t16-t18;
  t20 = 135.0/64.0*t1;
  t21 = 135.0/32.0*eta;
  t22 = 81.0/32.0*t4;
  t23 = 81.0/32.0*xi;
  t24 = 135.0/32.0*t7;
  t25 = 945.0/64.0*t10;
  t26 = 945.0/64.0*t9;
  t27 = 405.0/32.0*t13;
  t28 = 405.0/32.0*t15;
  t31 = t2-t3+t5+t6-t8-147.0/128.0-t11+t12-t14+t16-t18;
  t32 = 81.0/32.0*eta;
  t33 = 45.0/32.0*t4;
  t34 = 27.0/32.0*xi;
  t35 = 243.0/32.0*t7;
  t36 = 243.0/32.0*t13;
  t37 = 315.0/128.0*t17;
  t40 = t2+t3+t5+t6+t8-147.0/128.0-t11+t12+t14-t16-t18;
  t43 = t2+t3-t5-t6+t8-147.0/128.0-t11+t12-t14+t16-t18;
  t54 = 315.0/16.0*t1;

  values[0] = t19;
  values[1] = t20-t21+t22-t23+t24-135.0/64.0-t25+t26+t27-t28;
  values[2] = t20-t21-t22+t23+t24-135.0/64.0-t25+t26-t27+t28;
  values[3] = t31;
  values[4] = t20-t32-t33+t34+t35-27.0/128.0-t36+t28-t37;
  values[5] = t20+t32-t33+t34-t35-27.0/128.0+t36-t28-t37;
  values[6] = t40;
  values[7] = t20+t21-t22+t23-t24-135.0/64.0-t25+t26+t27-t28;
  values[8] = t20+t21+t22-t23-t24-135.0/64.0-t25+t26-t27+t28;
  values[9] = t43;
  values[10] = t20+t32+t33-t34-t35-27.0/128.0-t36+t28-t37;
  values[11] = t20-t32+t33-t34+t35-27.0/128.0+t36-t28-t37;
  values[12] = 111.0/32.0-27.0/4.0*t1+105.0/32.0*t17-315.0/16.0*t9+315.0/16.0*t10;
  values[13] = -135.0/8.0*eta+135.0/8.0*t7;
  values[14] = -45.0/8.0*xi+45.0/8.0*t4;
  values[15] = -675.0/8.0*t13+675.0/8.0*t15;
  values[16] = 315.0/16.0-t54-1575.0/16.0*t9+1575.0/16.0*t10;
  values[17] = 105.0/32.0-t54+525.0/32.0*t17;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_UL3_2D_Obj = new TBaseFunct2D
        (18, BF_C_Q_UL3_2D, BFUnitSquare, 
         C_Q_UL3_2D_Funct, C_Q_UL3_2D_DeriveXi,
         C_Q_UL3_2D_DeriveEta, C_Q_UL3_2D_DeriveXiXi,
         C_Q_UL3_2D_DeriveXiEta, C_Q_UL3_2D_DeriveEtaEta, 4, 3,
         0, NULL);
