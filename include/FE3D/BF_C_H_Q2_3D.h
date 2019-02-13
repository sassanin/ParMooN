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
// Q2 element, conforming, 3D
// ***********************************************************************

static void C_H_Q2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t2, t3, t4, t5, t6, t7, t9, t10, t11, t14, t15, t17, t18, t19;
  double t20, t22, t26, t27, t28, t33, t34, t35, t40, t45, t50, t51;
  double t55, t60;

  t2 = xi*(xi-1.0);
  t3 = t2*eta;
  t4 = eta-1.0;
  t5 = t4*zeta;
  t6 = zeta-1.0;
  t7 = t5*t6;
  t9 = xi*xi;
  t10 = 1.0-t9;
  t11 = t10*eta;
  t14 = xi*(-xi-1.0);
  t15 = t14*eta;
  t17 = eta*eta;
  t18 = 1.0-t17;
  t19 = t18*zeta;
  t20 = t19*t6;
  t22 = t10*t18;
  t26 = -eta-1.0;
  t27 = t26*zeta;
  t28 = t27*t6;
  t33 = zeta*zeta;
  t34 = 1.0-t33;
  t35 = eta*t4*t34;
  t40 = t18*t34;
  t45 = eta*t26*t34;
  t50 = -zeta-1.0;
  t51 = t5*t50;
  t55 = t19*t50;
  t60 = t27*t50;

  values[0] = t3*t7/8;
  values[1] = t11*t7/4;
  values[2] = -t15*t7/8;
  values[3] = t2*t20/4;
  values[4] = t22*zeta*t6/2;
  values[5] = -t14*t20/4;
  values[6] = -t3*t28/8;
  values[7] = -t11*t28/4;
  values[8] = t15*t28/8;
  values[9] = t2*t35/4;
  values[10] = t11*t4*t34/2;
  values[11] = -t14*t35/4;
  values[12] = t2*t40/2;
  values[13] = t22*t34;
  values[14] = -t14*t40/2;
  values[15] = -t2*t45/4;
  values[16] = -t11*t26*t34/2;
  values[17] = t14*t45/4;
  values[18] = -t3*t51/8;
  values[19] = -t11*t51/4;
  values[20] = t15*t51/8;
  values[21] = -t2*t55/4;
  values[22] = -t22*zeta*t50/2;
  values[23] = t14*t55/4;
  values[24] = t3*t60/8;
  values[25] = t11*t60/4;
  values[26] = -t15*t60/8;
}

static void C_H_Q2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t8, t9, t11, t12, t15, t16, t17, t18;
  double t20, t21, t23, t26, t27, t28, t30, t34, t35, t36, t38, t43, t47;
  double t49, t53, t54, t56, t60, t62, t66, t68;

  t1 = xi-1.0;
  t2 = t1*eta;
  t3 = eta-1.0;
  t4 = t3*zeta;
  t5 = zeta-1.0;
  t6 = t4*t5;
  t8 = xi*eta;
  t9 = t8*t6;
  t11 = -xi-1.0;
  t12 = t11*eta;
  t15 = eta*eta;
  t16 = 1.0-t15;
  t17 = t1*t16;
  t18 = zeta*t5;
  t20 = xi*t16;
  t21 = t20*t18;
  t23 = t11*t16;
  t26 = -eta-1.0;
  t27 = t26*zeta;
  t28 = t27*t5;
  t30 = t8*t28;
  t34 = zeta*zeta;
  t35 = 1.0-t34;
  t36 = t3*t35;
  t38 = t8*t36;
  t43 = t20*t35;
  t47 = t26*t35;
  t49 = t8*t47;
  t53 = -zeta-1.0;
  t54 = t4*t53;
  t56 = t8*t54;
  t60 = zeta*t53;
  t62 = t20*t60;
  t66 = t27*t53;
  t68 = t8*t66;

  values[0] = t2*t6/8+t9/8;
  values[1] = -t9/2;
  values[2] = -t12*t6/8+t9/8;
  values[3] = t17*t18/4+t21/4;
  values[4] = -t21;
  values[5] = -t23*t18/4+t21/4;
  values[6] = -t2*t28/8-t30/8;
  values[7] = t30/2;
  values[8] = t12*t28/8-t30/8;
  values[9] = t2*t36/4+t38/4;
  values[10] = -t38;
  values[11] = -t12*t36/4+t38/4;
  values[12] = t17*t35/2+t43/2;
  values[13] = -2.0*t43;
  values[14] = -t23*t35/2+t43/2;
  values[15] = -t2*t47/4-t49/4;
  values[16] = t49;
  values[17] = t12*t47/4-t49/4;
  values[18] = -t2*t54/8-t56/8;
  values[19] = t56/2;
  values[20] = t12*t54/8-t56/8;
  values[21] = -t17*t60/4-t62/4;
  values[22] = t62;
  values[23] = t23*t60/4-t62/4;
  values[24] = t2*t66/8+t68/8;
  values[25] = -t68/2;
  values[26] = -t12*t66/8+t68/8;
}

static void C_H_Q2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t6, t8, t9, t10, t12, t13, t14, t15, t17, t18;
  double t21, t23, t25, t26, t27, t30, t35, t36, t37, t39, t40, t43, t46;
  double t48, t55, t56, t58, t59, t61, t63, t66, t68;

  t2 = xi*(xi-1.0);
  t3 = eta-1.0;
  t4 = t3*zeta;
  t5 = zeta-1.0;
  t6 = t4*t5;
  t8 = eta*zeta;
  t9 = t8*t5;
  t10 = t2*t9;
  t12 = xi*xi;
  t13 = 1.0-t12;
  t14 = t13*t3;
  t15 = zeta*t5;
  t17 = t13*eta;
  t18 = t17*t15;
  t21 = xi*(-xi-1.0);
  t23 = t21*t9;
  t25 = -eta-1.0;
  t26 = t25*zeta;
  t27 = t26*t5;
  t30 = t13*t25;
  t35 = zeta*zeta;
  t36 = 1.0-t35;
  t37 = t3*t36;
  t39 = eta*t36;
  t40 = t2*t39;
  t43 = t17*t36;
  t46 = t21*t39;
  t48 = t25*t36;
  t55 = -zeta-1.0;
  t56 = t4*t55;
  t58 = t8*t55;
  t59 = t2*t58;
  t61 = zeta*t55;
  t63 = t17*t61;
  t66 = t21*t58;
  t68 = t26*t55;

  values[0] = t2*t6/8+t10/8;
  values[1] = t14*t15/4+t18/4;
  values[2] = -t21*t6/8-t23/8;
  values[3] = -t10/2;
  values[4] = -t18;
  values[5] = t23/2;
  values[6] = -t2*t27/8+t10/8;
  values[7] = -t30*t15/4+t18/4;
  values[8] = t21*t27/8-t23/8;
  values[9] = t2*t37/4+t40/4;
  values[10] = t14*t36/2+t43/2;
  values[11] = -t21*t37/4-t46/4;
  values[12] = -t40;
  values[13] = -2.0*t43;
  values[14] = t46;
  values[15] = -t2*t48/4+t40/4;
  values[16] = -t30*t36/2+t43/2;
  values[17] = t21*t48/4-t46/4;
  values[18] = -t2*t56/8-t59/8;
  values[19] = -t14*t61/4-t63/4;
  values[20] = t21*t56/8+t66/8;
  values[21] = t59/2;
  values[22] = t63;
  values[23] = -t66/2;
  values[24] = t2*t68/8-t59/8;
  values[25] = t30*t61/4-t63/4;
  values[26] = -t21*t68/8+t66/8;
}

static void C_H_Q2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t6, t8, t9, t11, t12, t13, t17, t20, t22, t24;
  double t25, t26, t28, t29, t31, t33, t36, t38, t39, t40, t42, t43, t48;
  double t51, t53, t54, t62, t69;

  t2 = xi*(xi-1.0);
  t3 = eta-1.0;
  t4 = eta*t3;
  t5 = zeta-1.0;
  t6 = t4*t5;
  t8 = t4*zeta;
  t9 = t2*t8;
  t11 = xi*xi;
  t12 = 1.0-t11;
  t13 = t12*eta;
  t17 = t13*t3*zeta;
  t20 = xi*(-xi-1.0);
  t22 = t20*t8;
  t24 = eta*eta;
  t25 = 1.0-t24;
  t26 = t25*t5;
  t28 = t25*zeta;
  t29 = t2*t28;
  t31 = t12*t25;
  t33 = t31*zeta;
  t36 = t20*t28;
  t38 = -eta-1.0;
  t39 = eta*t38;
  t40 = t39*t5;
  t42 = t39*zeta;
  t43 = t2*t42;
  t48 = t13*t38*zeta;
  t51 = t20*t42;
  t53 = -zeta-1.0;
  t54 = t4*t53;
  t62 = t25*t53;
  t69 = t39*t53;

  values[0] = t2*t6/8+t9/8;
  values[1] = t13*t3*t5/4+t17/4;
  values[2] = -t20*t6/8-t22/8;
  values[3] = t2*t26/4+t29/4;
  values[4] = t31*t5/2+t33/2;
  values[5] = -t20*t26/4-t36/4;
  values[6] = -t2*t40/8-t43/8;
  values[7] = -t13*t38*t5/4-t48/4;
  values[8] = t20*t40/8+t51/8;
  values[9] = -t9/2;
  values[10] = -t17;
  values[11] = t22/2;
  values[12] = -t29;
  values[13] = -2.0*t33;
  values[14] = t36;
  values[15] = t43/2;
  values[16] = t48;
  values[17] = -t51/2;
  values[18] = -t2*t54/8+t9/8;
  values[19] = -t13*t3*t53/4+t17/4;
  values[20] = t20*t54/8-t22/8;
  values[21] = -t2*t62/4+t29/4;
  values[22] = -t31*t53/2+t33/2;
  values[23] = t20*t62/4-t36/4;
  values[24] = t2*t69/8-t43/8;
  values[25] = t13*t38*t53/4-t48/4;
  values[26] = -t20*t69/8+t51/8;
}

static void C_H_Q2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t6, t7, t8, t9, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t21, t22;

  t2 = eta*(eta-1.0);
  t3 = zeta-1.0;
  t4 = t3*zeta;
  t5 = t2*t4;
  t6 = eta*eta;
  t7 = 1.0-t6;
  t8 = t7*zeta;
  t9 = t8*t3;
  t11 = eta*(-eta-1.0);
  t12 = t11*t4;
  t13 = zeta*zeta;
  t14 = 1.0-t13;
  t15 = t2*t14;
  t16 = t7*t14;
  t17 = t11*t14;
  t18 = -zeta-1.0;
  t19 = zeta*t18;
  t20 = t2*t19;
  t21 = t8*t18;
  t22 = t11*t19;

  values[0] = t5/4;
  values[1] = -t5/2;
  values[2] = t5/4;
  values[3] = t9/2;
  values[4] = -t9;
  values[5] = t9/2;
  values[6] = -t12/4;
  values[7] = t12/2;
  values[8] = -t12/4;
  values[9] = t15/2;
  values[10] = -t15;
  values[11] = t15/2;
  values[12] = t16;
  values[13] = -2.0*t16;
  values[14] = t16;
  values[15] = -t17/2;
  values[16] = t17;
  values[17] = -t17/2;
  values[18] = -t20/4;
  values[19] = t20/2;
  values[20] = -t20/4;
  values[21] = -t21/2;
  values[22] = t21;
  values[23] = -t21/2;
  values[24] = t22/4;
  values[25] = -t22/2;
  values[26] = t22/4;
}

static void C_H_Q2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t5, t7, t8, t9, t10, t11, t12, t15, t16, t18, t19;
  double t23, t24, t26, t27, t30, t33, t34, t36, t37, t38, t42, t47, t53;
  double t55, t56, t57, t61, t66;

  t1 = xi-1.0;
  t2 = eta-1.0;
  t3 = t1*t2;
  t5 = zeta*(zeta-1.0);
  t7 = t1*eta;
  t8 = t7*t5;
  t9 = xi*t2;
  t10 = t9*t5;
  t11 = xi*eta;
  t12 = t11*t5;
  t15 = -xi-1.0;
  t16 = t15*t2;
  t18 = t15*eta;
  t19 = t18*t5;
  t23 = -eta-1.0;
  t24 = t1*t23;
  t26 = xi*t23;
  t27 = t26*t5;
  t30 = t15*t23;
  t33 = zeta*zeta;
  t34 = 1.0-t33;
  t36 = t7*t34;
  t37 = t9*t34;
  t38 = t11*t34;
  t42 = t18*t34;
  t47 = t26*t34;
  t53 = zeta*(-zeta-1.0);
  t55 = t7*t53;
  t56 = t9*t53;
  t57 = t11*t53;
  t61 = t18*t53;
  t66 = t26*t53;

  values[0] = t3*t5/8+t8/8+t10/8+t12/8;
  values[1] = -t10/2-t12/2;
  values[2] = -t16*t5/8-t19/8+t10/8+t12/8;
  values[3] = -t8/2-t12/2;
  values[4] = 2.0*t12;
  values[5] = t19/2-t12/2;
  values[6] = -t24*t5/8+t8/8-t27/8+t12/8;
  values[7] = t27/2-t12/2;
  values[8] = t30*t5/8-t19/8-t27/8+t12/8;
  values[9] = t3*t34/4+t36/4+t37/4+t38/4;
  values[10] = -t37-t38;
  values[11] = -t16*t34/4-t42/4+t37/4+t38/4;
  values[12] = -t36-t38;
  values[13] = 4.0*t38;
  values[14] = t42-t38;
  values[15] = -t24*t34/4+t36/4-t47/4+t38/4;
  values[16] = t47-t38;
  values[17] = t30*t34/4-t42/4-t47/4+t38/4;
  values[18] = -t3*t53/8-t55/8-t56/8-t57/8;
  values[19] = t56/2+t57/2;
  values[20] = t16*t53/8+t61/8-t56/8-t57/8;
  values[21] = t55/2+t57/2;
  values[22] = -2.0*t57;
  values[23] = -t61/2+t57/2;
  values[24] = t24*t53/8-t55/8+t66/8-t57/8;
  values[25] = -t66/2+t57/2;
  values[26] = -t30*t53/8+t61/8+t66/8-t57/8;
}

static void C_H_Q2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t7, t8, t9, t10, t11, t14, t15, t17, t19, t20;
  double t21, t23, t24, t25, t26, t29, t31, t33, t34, t36, t37, t38, t39;
  double t43, t51, t52, t54, t60, t65, t67;

  t1 = xi-1.0;
  t2 = t1*eta;
  t3 = eta-1.0;
  t4 = zeta-1.0;
  t5 = t3*t4;
  t7 = t3*zeta;
  t8 = t2*t7;
  t9 = xi*eta;
  t10 = t9*t5;
  t11 = t9*t7;
  t14 = -xi-1.0;
  t15 = t14*eta;
  t17 = t15*t7;
  t19 = eta*eta;
  t20 = 1.0-t19;
  t21 = t1*t20;
  t23 = t21*zeta;
  t24 = xi*t20;
  t25 = t24*t4;
  t26 = t24*zeta;
  t29 = t14*t20;
  t31 = t29*zeta;
  t33 = -eta-1.0;
  t34 = t33*t4;
  t36 = t33*zeta;
  t37 = t2*t36;
  t38 = t9*t34;
  t39 = t9*t36;
  t43 = t15*t36;
  t51 = -zeta-1.0;
  t52 = t3*t51;
  t54 = t9*t52;
  t60 = t24*t51;
  t65 = t33*t51;
  t67 = t9*t65;

  values[0] = t2*t5/8+t8/8+t10/8+t11/8;
  values[1] = -t10/2-t11/2;
  values[2] = -t15*t5/8-t17/8+t10/8+t11/8;
  values[3] = t21*t4/4+t23/4+t25/4+t26/4;
  values[4] = -t25-t26;
  values[5] = -t29*t4/4-t31/4+t25/4+t26/4;
  values[6] = -t2*t34/8-t37/8-t38/8-t39/8;
  values[7] = t38/2+t39/2;
  values[8] = t15*t34/8+t43/8-t38/8-t39/8;
  values[9] = -t8/2-t11/2;
  values[10] = 2.0*t11;
  values[11] = t17/2-t11/2;
  values[12] = -t23-t26;
  values[13] = 4.0*t26;
  values[14] = t31-t26;
  values[15] = t37/2+t39/2;
  values[16] = -2.0*t39;
  values[17] = -t43/2+t39/2;
  values[18] = -t2*t52/8+t8/8-t54/8+t11/8;
  values[19] = t54/2-t11/2;
  values[20] = t15*t52/8-t17/8-t54/8+t11/8;
  values[21] = -t21*t51/4+t23/4-t60/4+t26/4;
  values[22] = t60-t26;
  values[23] = t29*t51/4-t31/4-t60/4+t26/4;
  values[24] = t2*t65/8-t37/8+t67/8-t39/8;
  values[25] = -t67/2+t39/2;
  values[26] = -t15*t65/8+t43/8+t67/8-t39/8;
}

static void C_H_Q2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t6, t7, t8, t9, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t21, t22;

  t2 = xi*(xi-1.0);
  t3 = zeta-1.0;
  t4 = t3*zeta;
  t5 = t2*t4;
  t6 = xi*xi;
  t7 = 1.0-t6;
  t8 = t7*zeta;
  t9 = t8*t3;
  t11 = xi*(-xi-1.0);
  t12 = t11*t4;
  t13 = zeta*zeta;
  t14 = 1.0-t13;
  t15 = t2*t14;
  t16 = t7*t14;
  t17 = t11*t14;
  t18 = -zeta-1.0;
  t19 = zeta*t18;
  t20 = t2*t19;
  t21 = t8*t18;
  t22 = t11*t19;

  values[0] = t5/4;
  values[1] = t9/2;
  values[2] = -t12/4;
  values[3] = -t5/2;
  values[4] = -t9;
  values[5] = t12/2;
  values[6] = t5/4;
  values[7] = t9/2;
  values[8] = -t12/4;
  values[9] = t15/2;
  values[10] = t16;
  values[11] = -t17/2;
  values[12] = -t15;
  values[13] = -2.0*t16;
  values[14] = t17;
  values[15] = t15/2;
  values[16] = t16;
  values[17] = -t17/2;
  values[18] = -t20/4;
  values[19] = -t21/2;
  values[20] = t22/4;
  values[21] = t20/2;
  values[22] = t21;
  values[23] = -t22/2;
  values[24] = -t20/4;
  values[25] = -t21/2;
  values[26] = t22/4;
}

static void C_H_Q2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t7, t8, t9, t10, t11, t12, t14, t15, t16, t18;
  double t19, t20, t21, t24, t26, t27, t28, t33, t34, t36, t37, t39, t41;
  double t44, t52, t53, t55, t56, t59, t62, t67;

  t2 = xi*(xi-1.0);
  t3 = eta-1.0;
  t4 = zeta-1.0;
  t5 = t3*t4;
  t7 = t3*zeta;
  t8 = t2*t7;
  t9 = eta*t4;
  t10 = t2*t9;
  t11 = eta*zeta;
  t12 = t2*t11;
  t14 = xi*xi;
  t15 = 1.0-t14;
  t16 = t15*t3;
  t18 = t16*zeta;
  t19 = t15*eta;
  t20 = t19*t4;
  t21 = t19*zeta;
  t24 = xi*(-xi-1.0);
  t26 = t24*t7;
  t27 = t24*t9;
  t28 = t24*t11;
  t33 = -eta-1.0;
  t34 = t33*t4;
  t36 = t33*zeta;
  t37 = t2*t36;
  t39 = t15*t33;
  t41 = t39*zeta;
  t44 = t24*t36;
  t52 = -zeta-1.0;
  t53 = t3*t52;
  t55 = eta*t52;
  t56 = t2*t55;
  t59 = t19*t52;
  t62 = t24*t55;
  t67 = t33*t52;

  values[0] = t2*t5/8+t8/8+t10/8+t12/8;
  values[1] = t16*t4/4+t18/4+t20/4+t21/4;
  values[2] = -t24*t5/8-t26/8-t27/8-t28/8;
  values[3] = -t10/2-t12/2;
  values[4] = -t20-t21;
  values[5] = t27/2+t28/2;
  values[6] = -t2*t34/8-t37/8+t10/8+t12/8;
  values[7] = -t39*t4/4-t41/4+t20/4+t21/4;
  values[8] = t24*t34/8+t44/8-t27/8-t28/8;
  values[9] = -t8/2-t12/2;
  values[10] = -t18-t21;
  values[11] = t26/2+t28/2;
  values[12] = 2.0*t12;
  values[13] = 4.0*t21;
  values[14] = -2.0*t28;
  values[15] = t37/2-t12/2;
  values[16] = t41-t21;
  values[17] = -t44/2+t28/2;
  values[18] = -t2*t53/8+t8/8-t56/8+t12/8;
  values[19] = -t16*t52/4+t18/4-t59/4+t21/4;
  values[20] = t24*t53/8-t26/8+t62/8-t28/8;
  values[21] = t56/2-t12/2;
  values[22] = t59-t21;
  values[23] = -t62/2+t28/2;
  values[24] = t2*t67/8-t37/8-t56/8+t12/8;
  values[25] = t39*t52/4-t41/4-t59/4+t21/4;
  values[26] = -t24*t67/8+t44/8+t62/8-t28/8;
}

static void C_H_Q2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t6, t7, t8, t9, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t21, t22;

  t2 = xi*(xi-1.0);
  t3 = eta-1.0;
  t4 = eta*t3;
  t5 = t2*t4;
  t6 = xi*xi;
  t7 = 1.0-t6;
  t8 = t7*eta;
  t9 = t8*t3;
  t11 = xi*(-xi-1.0);
  t12 = t11*t4;
  t13 = eta*eta;
  t14 = 1.0-t13;
  t15 = t2*t14;
  t16 = t7*t14;
  t17 = t11*t14;
  t18 = -eta-1.0;
  t19 = eta*t18;
  t20 = t2*t19;
  t21 = t8*t18;
  t22 = t11*t19;

  values[0] = t5/4;
  values[1] = t9/2;
  values[2] = -t12/4;
  values[3] = t15/2;
  values[4] = t16;
  values[5] = -t17/2;
  values[6] = -t20/4;
  values[7] = -t21/2;
  values[8] = t22/4;
  values[9] = -t5/2;
  values[10] = -t9;
  values[11] = t12/2;
  values[12] = -t15;
  values[13] = -2.0*t16;
  values[14] = t17;
  values[15] = t20/2;
  values[16] = t21;
  values[17] = -t22/2;
  values[18] = t5/4;
  values[19] = t9/2;
  values[20] = -t12/4;
  values[21] = t15/2;
  values[22] = t16;
  values[23] = -t17/2;
  values[24] = -t20/4;
  values[25] = -t21/2;
  values[26] = t22/4;
}

TBaseFunct3D *BF_C_H_Q2_3D_Obj = 
new TBaseFunct3D(27, BF_C_H_Q2_3D, BFUnitHexahedron, 
                 C_H_Q2_3D_Funct, C_H_Q2_3D_DeriveXi,
                 C_H_Q2_3D_DeriveEta, C_H_Q2_3D_DeriveZeta,
                 C_H_Q2_3D_DeriveXiXi, C_H_Q2_3D_DeriveXiEta,
                 C_H_Q2_3D_DeriveXiZeta, C_H_Q2_3D_DeriveEtaEta,
                 C_H_Q2_3D_DeriveEtaZeta, C_H_Q2_3D_DeriveZetaZeta,
                 2, 2,
                 0, NULL);
