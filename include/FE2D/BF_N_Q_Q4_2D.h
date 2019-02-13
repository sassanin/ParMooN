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
// Q4 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_Q_Q4_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15, t17;
  double t19, t20, t21, t22, t23, t24, t26, t27, t28, t29, t31, t35, t36;
  double t37, t38, t39, t40, t41, t43, t44, t45, t46, t47, t54, t61, t63;
  double t64, t65, t66, t67, t69, t70, t71, t75, t76, t77, t78, t79, t80;
  double t81, t82, t100, t101, t102;

  t1 = 3.0/4.0*eta;
  t2 = eta*eta;
  t3 = 15.0/8.0*t2;
  t4 = t2*eta;
  t5 = 5.0/4.0*t4;
  t6 = t2*t2;
  t7 = 35.0/16.0*t6;
  t8 = xi*xi;
  t9 = t8*t8;
  t10 = 35.0/8.0*t9;
  t11 = 15.0/4.0*t8;
  t12 = 3.0/8.0+t10-t11;
  t13 = 5.0/2.0*t4;
  t15 = t13-3.0/2.0*eta;
  t17 = t12*t15/2.0;
  t19 = 3.0/4.0*xi;
  t20 = 15.0/8.0*t8;
  t21 = t8*xi;
  t22 = 5.0/4.0*t21;
  t23 = 35.0/16.0*t9;
  t24 = 5.0/2.0*t21;
  t26 = t24-3.0/2.0*xi;
  t27 = 35.0/8.0*t6;
  t28 = 15.0/4.0*t2;
  t29 = 3.0/8.0+t27-t28;
  t31 = t26*t29/2.0;
  t35 = xi/4.0;
  t36 = xi*eta;
  t37 = 3.0/4.0*t36;
  t38 = xi*t2;
  t39 = 3.0/4.0*t38;
  t40 = xi*t4;
  t41 = 5.0/4.0*t40;
  t43 = -1.0/2.0+3.0/2.0*t2;
  t44 = t26*t43;
  t45 = t44/2.0;
  t46 = t26*t15;
  t47 = t46/4.0;
  t54 = (63.0/8.0*t9*xi-35.0/4.0*t21+15.0/8.0*xi)*t15/4.0;
  t61 = t26*(63.0/8.0*t6*eta-35.0/4.0*t4+15.0/8.0*eta)/4.0;
  t63 = eta/4.0;
  t64 = t8*eta;
  t65 = 3.0/4.0*t64;
  t66 = t21*eta;
  t67 = 5.0/4.0*t66;
  t69 = -1.0/2.0+3.0/2.0*t8;
  t70 = t69*t15;
  t71 = t70/2.0;
  t75 = 3.0/16.0*t8;
  t76 = 3.0/16.0*t2;
  t77 = t8*t2;
  t78 = 9.0/16.0*t77;
  t79 = t12*t43;
  t80 = t79/4.0;
  t81 = t69*t29;
  t82 = t81/4.0;
  t100 = 9.0/8.0*t77;
  t101 = t79/2.0;
  t102 = t81/2.0;

  values[0] = 3.0/16.0+t1-t3-t5+t7+t17;
  values[1] = 3.0/16.0-t19-t20+t22+t23-t31;
  values[2] = 3.0/16.0-t1-t3+t5+t7-t17;
  values[3] = 3.0/16.0+t19-t20-t22+t23+t31;
  values[4] = -t35+t37+t39-t41-t45+t47+t31+t54-t61;
  values[5] = -t63-t37+t65+t67-t71-t47+t17+t54-t61;
  values[6] = t35+t37-t39-t41+t45+t47-t31+t54-t61;
  values[7] = t63-t37-t65+t67+t71-t47-t17+t54-t61;
  values[8] = 1.0/16.0-t75-t76+t78-t71-t80+t82+t17;
  values[9] = 1.0/16.0-t75-t76+t78+t45+t80-t82-t31;
  values[10] = 1.0/16.0-t75-t76+t78+t71-t80+t82-t17;
  values[11] = 1.0/16.0-t75-t76+t78-t45+t80-t82+t31;
  values[12] = -t47+t31+t54-t61;
  values[13] = t47+t17+t54-t61;
  values[14] = -t47-t31+t54-t61;
  values[15] = t47-t17+t54-t61;
  values[16] = 1.0/4.0+t11+t28-t10-t27;
  values[17] = 3.0*xi-t24-3.0/2.0*t38+t44;
  values[18] = 3.0*eta-3.0/2.0*t64-t13+t70;
  values[19] = -1.0+45.0/8.0*t8+3.0/8.0*t2-t10-t100+t101-t102;
  values[20] = 4.0*t36-5.0/2.0*t66-5.0/2.0*t40+t46;
  values[21] = -1.0+3.0/8.0*t8+45.0/8.0*t2-t100-t27-t101+t102;
}

static void N_Q_Q4_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t5, t6, t7, t8, t10, t12, t13, t14, t15, t16, t17;
  double t18, t21, t23, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35;
  double t40, t47, t49, t50, t51, t52, t53, t54, t58, t59, t60, t61, t62;
  double t63, t64, t76, t77, t78;

  t1 = xi*xi;
  t2 = t1*xi;
  t3 = 35.0/2.0*t2;
  t5 = t3-15.0/2.0*xi;
  t6 = eta*eta;
  t7 = t6*eta;
  t8 = 5.0/2.0*t7;
  t10 = t8-3.0/2.0*eta;
  t12 = t5*t10/2.0;
  t13 = 15.0/4.0*xi;
  t14 = 15.0/4.0*t1;
  t15 = 35.0/4.0*t2;
  t16 = 15.0/2.0*t1;
  t17 = t16-3.0/2.0;
  t18 = t6*t6;
  t21 = 3.0/8.0+35.0/8.0*t18-15.0/4.0*t6;
  t23 = t17*t21/2.0;
  t26 = 3.0/4.0*eta;
  t27 = 3.0/4.0*t6;
  t28 = 5.0/4.0*t7;
  t29 = 3.0/2.0*t6;
  t30 = -1.0/2.0+t29;
  t31 = t17*t30;
  t32 = t31/2.0;
  t33 = t17*t10;
  t34 = t33/4.0;
  t35 = t1*t1;
  t40 = (315.0/8.0*t35-105.0/4.0*t1+15.0/8.0)*t10/4.0;
  t47 = t17*(63.0/8.0*t18*eta-35.0/4.0*t7+15.0/8.0*eta)/4.0;
  t49 = xi*eta;
  t50 = 3.0/2.0*t49;
  t51 = t1*eta;
  t52 = 15.0/4.0*t51;
  t53 = xi*t10;
  t54 = 3.0/2.0*t53;
  t58 = 3.0/8.0*xi;
  t59 = xi*t6;
  t60 = 9.0/8.0*t59;
  t61 = t5*t30;
  t62 = t61/4.0;
  t63 = xi*t21;
  t64 = 3.0/4.0*t63;
  t76 = 9.0/4.0*t59;
  t77 = t61/2.0;
  t78 = 3.0/2.0*t63;

  values[0] = t12;
  values[1] = -3.0/4.0-t13+t14+t15-t23;
  values[2] = -t12;
  values[3] = 3.0/4.0-t13-t14+t15+t23;
  values[4] = -1.0/4.0+t26+t27-t28-t32+t34+t23+t40-t47;
  values[5] = -t26+t50+t52-t54-t34+t12+t40-t47;
  values[6] = 1.0/4.0+t26-t27-t28+t32+t34-t23+t40-t47;
  values[7] = -t26-t50+t52+t54-t34-t12+t40-t47;
  values[8] = -t58+t60-t54-t62+t64+t12;
  values[9] = -t58+t60+t32+t62-t64-t23;
  values[10] = -t58+t60+t54-t62+t64-t12;
  values[11] = -t58+t60-t32+t62-t64+t23;
  values[12] = -t34+t23+t40-t47;
  values[13] = t34+t12+t40-t47;
  values[14] = -t34-t23+t40-t47;
  values[15] = t34-t12+t40-t47;
  values[16] = -t5;
  values[17] = 3.0-t16-t29+t31;
  values[18] = -3.0*t49+3.0*t53;
  values[19] = 45.0/4.0*xi-t3-t76+t77-t78;
  values[20] = 4.0*eta-15.0/2.0*t51-t8+t33;
  values[21] = 3.0/4.0*xi-t76-t77+t78;
}

static void N_Q_Q4_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t10, t11, t12, t14, t16, t17, t19;
  double t20, t22, t24, t26, t27, t28, t29, t30, t31, t32, t33, t34, t41;
  double t42, t47, t49, t50, t51, t52, t53, t54, t58, t59, t60, t61, t62;
  double t63, t64, t76, t77, t78;

   t1 = 15.0/4.0*eta;
   t2 = eta*eta;
   t3 = 15.0/4.0*t2;
   t4 = t2*eta;
   t5 = 35.0/4.0*t4;
   t6 = xi*xi;
   t7 = t6*t6;
   t10 = 3.0/8.0+35.0/8.0*t7-15.0/4.0*t6;
   t11 = 15.0/2.0*t2;
   t12 = t11-3.0/2.0;
   t14 = t10*t12/2.0;
   t16 = xi*t6;
   t17 = 5.0/2.0*t16;
   t19 = t17-3.0/2.0*xi;
   t20 = 35.0/2.0*t4;
   t22 = t20-15.0/2.0*eta;
   t24 = t19*t22/2.0;
   t26 = 3.0/4.0*xi;
   t27 = xi*eta;
   t28 = 3.0/2.0*t27;
   t29 = xi*t2;
   t30 = 15.0/4.0*t29;
   t31 = t19*eta;
   t32 = 3.0/2.0*t31;
   t33 = t19*t12;
   t34 = t33/4.0;
   t41 = (63.0/8.0*t7*xi-35.0/4.0*t16+15.0/8.0*xi)*t12/4.0;
   t42 = t2*t2;
   t47 = t19*(315.0/8.0*t42-105.0/4.0*t2+15.0/8.0)/4.0;
   t49 = 3.0/4.0*t6;
   t50 = 5.0/4.0*t16;
   t51 = 3.0/2.0*t6;
   t52 = -1.0/2.0+t51;
   t53 = t52*t12;
   t54 = t53/2.0;
   t58 = 3.0/8.0*eta;
   t59 = t6*eta;
   t60 = 9.0/8.0*t59;
   t61 = t10*eta;
   t62 = 3.0/4.0*t61;
   t63 = t52*t22;
   t64 = t63/4.0;
   t76 = 9.0/4.0*t59;
   t77 = 3.0/2.0*t61;
   t78 = t63/2.0;

  values[0] = 3.0/4.0-t1-t3+t5+t14;
  values[1] = -t24;
  values[2] = -3.0/4.0-t1+t3+t5-t14;
  values[3] = t24;
  values[4] = t26+t28-t30-t32+t34+t24+t41-t47;
  values[5] = -1.0/4.0-t26+t49+t50-t54-t34+t14+t41-t47;
  values[6] = t26-t28-t30+t32+t34-t24+t41-t47;
  values[7] = 1.0/4.0-t26-t49+t50+t54-t34-t14+t41-t47;
  values[8] = -t58+t60-t54-t62+t64+t14;
  values[9] = -t58+t60+t32+t62-t64-t24;
  values[10] = -t58+t60+t54-t62+t64-t14;
  values[11] = -t58+t60-t32+t62-t64+t24;
  values[12] = -t34+t24+t41-t47;
  values[13] = t34+t14+t41-t47;
  values[14] = -t34-t24+t41-t47;
  values[15] = t34-t14+t41-t47;
  values[16] = -t22;
  values[17] = -3.0*t27+3.0*t31;
  values[18] = 3.0-t51-t11+t53;
  values[19] = 3.0/4.0*eta-t76+t77-t78;
  values[20] = 4.0*xi-t17-15.0/2.0*t29+t33;
  values[21] = 45.0/4.0*eta-t76-t20-t77+t78;
}

static void N_Q_Q4_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t8, t10, t11, t12, t13, t18, t22, t23, t24;
  double t25, t26, t32, t39, t41, t42, t43, t44, t48, t49, t50, t51, t52;
  double t54, t65, t66;

  t1 = xi*xi;
  t2 = 105.0/2.0*t1;
  t3 = t2-15.0/2.0;
  t4 = eta*eta;
  t5 = t4*eta;
  t8 = 5.0/2.0*t5-3.0/2.0*eta;
  t10 = t3*t8/2.0;
  t11 = 15.0/2.0*xi;
  t12 = 105.0/4.0*t1;
  t13 = t4*t4;
  t18 = 15.0/2.0*xi*(3.0/8.0+35.0/8.0*t13-15.0/4.0*t4);
  t22 = -1.0/2.0+3.0/2.0*t4;
  t23 = xi*t22;
  t24 = 15.0/2.0*t23;
  t25 = t8*xi;
  t26 = 15.0/4.0*t25;
  t32 = (315.0/2.0*t1*xi-105.0/2.0*xi)*t8/4.0;
  t39 = 15.0/4.0*xi*(63.0/8.0*t13*eta-35.0/4.0*t5+15.0/8.0*eta);
  t41 = 15.0/4.0*eta;
  t42 = xi*eta;
  t43 = 15.0/2.0*t42;
  t44 = 15.0/4.0*t5;
  t48 = 27.0/16.0*t4;
  t49 = 9.0/4.0*eta;
  t50 = t3*t22;
  t51 = t50/4.0;
  t52 = 105.0/32.0*t13;
  t54 = 63.0/16.0*t4;
  t65 = t50/2.0;
  t66 = 105.0/16.0*t13;

  values[0] = t10;
  values[1] = -15.0/4.0+t11+t12-t18;
  values[2] = -t10;
  values[3] = -15.0/4.0-t11+t12+t18;
  values[4] = -t24+t26+t18+t32-t39;
  values[5] = t41+t43-t44-t26+t10+t32-t39;
  values[6] = t24+t26-t18+t32-t39;
  values[7] = -t41+t43+t44-t26-t10+t32-t39;
  values[8] = -3.0/32.0-t48-t44+t49-t51+t52+t10;
  values[9] = -21.0/32.0+t54+t24+t51-t52-t18;
  values[10] = -3.0/32.0-t48+t44-t49-t51+t52-t10;
  values[11] = -21.0/32.0+t54-t24+t51-t52+t18;
  values[12] = -t26+t18+t32-t39;
  values[13] = t26+t10+t32-t39;
  values[14] = -t26-t18+t32-t39;
  values[15] = t26-t10+t32-t39;
  values[16] = -t3;
  values[17] = -15.0*xi+15.0*t23;
  values[18] = -15.0/2.0*eta+15.0/2.0*t5;
  values[19] = 171.0/16.0-t2+27.0/8.0*t4+t65-t66;
  values[20] = -15.0*t42+15.0*t25;
  values[21] = 21.0/16.0-63.0/8.0*t4-t65+t66;
}

static void N_Q_Q4_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t5, t6, t7, t8, t10, t11, t12, t16, t18, t19, t20, t21, t22;
  double t23, t24, t25, t30, t31, t36, t38, t39, t40, t41, t45, t46, t47;
  double t48, t49, t50, t61, t62, t63;

  t1 = xi*xi;
  t5 = 35.0/2.0*t1*xi-15.0/2.0*xi;
  t6 = eta*eta;
  t7 = 15.0/2.0*t6;
  t8 = t7-3.0/2.0;
  t10 = t5*t8/2.0;
  t11 = 15.0/2.0*t1;
  t12 = t11-3.0/2.0;
  t16 = 35.0/2.0*t6*eta-15.0/2.0*eta;
  t18 = t12*t16/2.0;
  t19 = 3.0/2.0*eta;
  t20 = 15.0/4.0*t6;
  t21 = t12*eta;
  t22 = 3.0/2.0*t21;
  t23 = t12*t8;
  t24 = t23/4.0;
  t25 = t1*t1;
  t30 = (315.0/8.0*t25-105.0/4.0*t1+15.0/8.0)*t8/4.0;
  t31 = t6*t6;
  t36 = t12*(315.0/8.0*t31-105.0/4.0*t6+15.0/8.0)/4.0;
  t38 = 3.0/2.0*xi;
  t39 = 15.0/4.0*t1;
  t40 = t8*xi;
  t41 = 3.0/2.0*t40;
  t45 = xi*eta;
  t46 = 9.0/4.0*t45;
  t47 = t5*eta;
  t48 = 3.0/4.0*t47;
  t49 = xi*t16;
  t50 = 3.0/4.0*t49;
  t61 = 9.0/2.0*t45;
  t62 = 3.0/2.0*t47;
  t63 = 3.0/2.0*t49;

  values[0] = t10;
  values[1] = -t18;
  values[2] = -t10;
  values[3] = t18;
  values[4] = 3.0/4.0+t19-t20-t22+t24+t18+t30-t36;
  values[5] = -3.0/4.0+t38+t39-t41-t24+t10+t30-t36;
  values[6] = 3.0/4.0-t19-t20+t22+t24-t18+t30-t36;
  values[7] = -3.0/4.0-t38+t39+t41-t24-t10+t30-t36;
  values[8] = t46-t41-t48+t50+t10;
  values[9] = t46+t22+t48-t50-t18;
  values[10] = t46+t41-t48+t50-t10;
  values[11] = t46-t22+t48-t50+t18;
  values[12] = -t24+t18+t30-t36;
  values[13] = t24+t10+t30-t36;
  values[14] = -t24-t18+t30-t36;
  values[15] = t24-t10+t30-t36;
  values[16] = 0.0;
  values[17] = -3.0*eta+3.0*t21;
  values[18] = -3.0*xi+3.0*t40;
  values[19] = -t61+t62-t63;
  values[20] = 4.0-t11-t7+t23;
  values[21] = -t61-t62+t63;
}

static void N_Q_Q4_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t10, t12, t15, t16, t17, t19, t21, t22, t23;
  double t24, t25, t26, t33, t39, t42, t43, t44, t48, t49, t50, t51, t53;
  double t54, t65, t66;

  t1 = 15.0/2.0*eta;
  t2 = eta*eta;
  t3 = 105.0/4.0*t2;
  t4 = xi*xi;
  t5 = t4*t4;
  t10 = 15.0/2.0*(3.0/8.0+35.0/8.0*t5-15.0/4.0*t4)*eta;
  t12 = t4*xi;
  t15 = 5.0/2.0*t12-3.0/2.0*xi;
  t16 = 105.0/2.0*t2;
  t17 = t16-15.0/2.0;
  t19 = t15*t17/2.0;
  t21 = 15.0/4.0*xi;
  t22 = xi*eta;
  t23 = 15.0/2.0*t22;
  t24 = 15.0/4.0*t12;
  t25 = t15*eta;
  t26 = 15.0/4.0*t25;
  t33 = 15.0/4.0*(63.0/8.0*t5*xi-35.0/4.0*t12+15.0/8.0*xi)*eta;
  t39 = t15*(315.0/2.0*t2*eta-105.0/2.0*eta)/4.0;
  t42 = -1.0/2.0+3.0/2.0*t4;
  t43 = t42*eta;
  t44 = 15.0/2.0*t43;
  t48 = 63.0/16.0*t4;
  t49 = 105.0/32.0*t5;
  t50 = t42*t17;
  t51 = t50/4.0;
  t53 = 27.0/16.0*t4;
  t54 = 9.0/4.0*xi;
  t65 = 105.0/16.0*t5;
  t66 = t50/2.0;

   values[0] = -15.0/4.0-t1+t3+t10;
   values[1] = -t19;
   values[2] = -15.0/4.0+t1+t3-t10;
   values[3] = t19;
   values[4] = t21-t23-t24+t26+t19+t33-t39;
   values[5] = -t44-t26+t10+t33-t39;
   values[6] = -t21-t23+t24+t26-t19+t33-t39;
   values[7] = t44-t26-t10+t33-t39;
   values[8] = -21.0/32.0+t48-t44-t49+t51+t10;
   values[9] = -3.0/32.0-t53+t24-t54+t49-t51-t19;
   values[10] = -21.0/32.0+t48+t44-t49+t51-t10;
   values[11] = -3.0/32.0-t53-t24+t54+t49-t51+t19;
   values[12] = -t26+t19+t33-t39;
   values[13] = t26+t10+t33-t39;
   values[14] = -t26-t19+t33-t39;
   values[15] = t26-t10+t33-t39;
   values[16] = -t17;
   values[17] = -15.0/2.0*xi+15.0/2.0*t12;
   values[18] = -15.0*eta+15.0*t43;
   values[19] = 21.0/16.0-63.0/8.0*t4+t65-t66;
   values[20] = -15.0*t22+15.0*t25;
   values[21] = 171.0/16.0+27.0/8.0*t4-t16-t65+t66;
}

static int N_Q_Q4_2D_ChangeJ0[2] = { 4, 12 };
static int N_Q_Q4_2D_ChangeJ1[2] = { 5, 13 };
static int N_Q_Q4_2D_ChangeJ2[2] = { 6, 14 };
static int N_Q_Q4_2D_ChangeJ3[2] = { 7, 15 };

static int *N_Q_Q4_2D_Change[4] = { N_Q_Q4_2D_ChangeJ0, N_Q_Q4_2D_ChangeJ1,
                                    N_Q_Q4_2D_ChangeJ2, N_Q_Q4_2D_ChangeJ3 };
// ***********************************************************************

TBaseFunct2D *BF_N_Q_Q4_2D_Obj = new TBaseFunct2D
        (22, BF_N_Q_Q4_2D, BFUnitSquare, 
         N_Q_Q4_2D_Funct, N_Q_Q4_2D_DeriveXi,
         N_Q_Q4_2D_DeriveEta, N_Q_Q4_2D_DeriveXiXi,
         N_Q_Q4_2D_DeriveXiEta, N_Q_Q4_2D_DeriveEtaEta, 5, 4,
         2, N_Q_Q4_2D_Change);
