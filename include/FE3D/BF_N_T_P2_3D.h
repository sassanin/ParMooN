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
// P2 element, nonconforming, 3D
// ***********************************************************************

static void N_T_P2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t34, t35, t36, t37, t38, t39, t40, t41, t42;
  double t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t54, t55, t56;
  double t58, t59, t60, t61, t62, t63, t64, t65, t67, t69, t70, t71, t73;
  double t74, t75, t77, t78, t79, t82, t85, t86, t87;

  t1 = zeta*zeta;
  t3 = eta*eta;
  t4 = 135.0/2.0*t3;
  t5 = xi*xi;
  t6 = 135.0/2.0*t5;
  t7 = t1*zeta;
  t8 = 90.0*t7;
  t9 = xi*t1;
  t10 = 135.0*t9;
  t11 = eta*t1;
  t12 = 135.0*t11;
  t13 = t5*xi;
  t14 = 45.0*t13;
  t15 = t3*eta;
  t16 = 45.0*t15;
  t17 = t5*eta;
  t18 = 90.0*t17;
  t19 = t5*zeta;
  t20 = 45.0*t19;
  t21 = xi*t3;
  t22 = 90.0*t21;
  t23 = t3*zeta;
  t24 = 45.0*t23;
  t25 = 15.0*xi;
  t26 = 15.0*eta;
  t27 = xi*eta;
  t28 = 90.0*t27;
  t29 = eta*zeta;
  t30 = 75.0*t29;
  t31 = xi*zeta;
  t32 = 75.0*t31;
  t34 = 6.0+165.0*t1-t4-t6-t8-t10-t12+t14+t16+t18+t20+t22+t24+t25+t26-t28+t30+t32-78.0*zeta;
  t35 = 135.0/2.0*t1;
  t36 = 45.0*t7;
  t37 = 90.0*t9;
  t38 = 45.0*t11;
  t39 = 45.0*t17;
  t40 = 90.0*t19;
  t41 = 45.0*t21;
  t42 = 45.0*t27;
  t43 = 45.0*t29;
  t44 = 30.0*t31;
  t45 = 27.0*zeta;
  t46 = -t35+t6+t36+t37+t38-t14-t39-t40-t41+t24-t25-3.0/2.0+t42-t43-t44+t45;
  t47 = 45.0*t9;
  t48 = 90.0*t11;
  t49 = 90.0*t23;
  t50 = 30.0*t29;
  t51 = 45.0*t31;
  t52 = -t35+t4+t36+t47+t48-t16-t39+t20-t41-t49-t26-3.0/2.0+t42-t50-t51+t45;
  t54 = 90.0*t15;
  t55 = 135.0*t21;
  t56 = 135.0*t23;
  t58 = 75.0*t27;
  t59 = 90.0*t31;
  t60 = 15.0*zeta;
  t61 = 6.0-t35+165.0*t3-t6+t36+t37+t38+t14-t54+t39+t40-t55-t56+t25-78.0*eta+t58+t30-t59+t60;
  t62 = 27.0*eta;
  t63 = t35-t4-t36-t47-t48+t16+t39-t20+t41+t49+t62-3.0/2.0-t42-t50+t51-t60;
  t64 = 30.0*t27;
  t65 = -t4+t6-t47+t38-t14+t16-t18-t20+t22+t24-t25+t62-3.0/2.0-t64-t43+t51;
  t67 = 27.0*xi;
  t69 = 15.0*t27;
  t70 = 15.0*t29;
  t71 = 3.0+t35-105.0*t3+t6-t36-t37-t38-t14+t54-t39-t40+t55+t56-t67+18.0*eta-t69-t70+t59-t45;
  t73 = 90.0*t13;
  t74 = 135.0*t17;
  t75 = 135.0*t19;
  t77 = 90.0*t29;
  t78 = 15.0*t31;
  t79 = 3.0+t35+t4-105.0*t5-t36-t47-t48+t73-t16+t74+t75-t41-t49+18.0*xi-t62-t69+t77-t78-t45;
  t82 = 3.0-105.0*t1+t4+t6+t8+t10+t12-t14-t16-t18-t20-t22-t24-t67-t62+t28-t70-t78+18.0*zeta;
  t85 = 6.0-t35-t4+165.0*t5+t36+t47+t48-t73+t16-t74-t75+t41+t49-78.0*xi+t26+t58-t77+t32+t60;
  t86 = t4-t6+t47-t38+t14-t16+t18+t20-t22-t24+t67-t26-3.0/2.0-t64+t43-t51;
  t87 = t35-t6-t36-t37-t38+t14+t39+t40+t41-t24+t67-3.0/2.0-t42+t43-t44-t60;

  values[0] = t34;
  values[1] = t46;
  values[2] = t52;
  values[3] = t61;
  values[4] = t63;
  values[5] = t65;
  values[6] = t71;
  values[7] = t79;
  values[8] = t82;
  values[9] = t85;
  values[10] = t86;
  values[11] = t87;
  values[12] = -5.0-20.0*t1-20.0*t3-20.0*t5+20.0*xi+20.0*eta-20.0*t27-20.0*t29-20.0*t31+20.0*zeta;
}

static void N_T_P2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15, t16;
  double t17, t18, t19, t20, t22, t23, t24, t25, t26, t27, t29, t31, t34;
  double t35, t36, t37;

  t1 = 135.0*xi;
  t2 = zeta*zeta;
  t3 = 135.0*t2;
  t4 = xi*xi;
  t5 = 135.0*t4;
  t6 = xi*eta;
  t7 = 180.0*t6;
  t8 = xi*zeta;
  t9 = 90.0*t8;
  t10 = eta*eta;
  t11 = 90.0*t10;
  t12 = 90.0*eta;
  t13 = 75.0*zeta;
  t15 = 90.0*t2;
  t16 = 90.0*t6;
  t17 = 180.0*t8;
  t18 = 45.0*t10;
  t19 = 45.0*eta;
  t20 = 30.0*zeta;
  t22 = 45.0*t2;
  t23 = 45.0*zeta;
  t24 = t22-t16+t9-t18+t19-t23;
  t25 = 135.0*t10;
  t26 = 75.0*eta;
  t27 = 90.0*zeta;
  t29 = 30.0*eta;
  t31 = 15.0*eta;
  t34 = 270.0*t4;
  t35 = 270.0*t6;
  t36 = 270.0*t8;
  t37 = 15.0*zeta;

  values[0] = -t1-t3+t5+t7+t9+t11+15.0-t12+t13;
  values[1] = t1+t15-t5-t16-t17-t18-15.0+t19-t20;
  values[2] = t24;
  values[3] = -t1+t15+t5+t16+t17-t25+15.0+t26-t27;
  values[4] = -t24;
  values[5] = t1-t22-t5-t7-t9+t11-15.0-t29+t23;
  values[6] = t1-t15-t5-t16-t17+t25-27.0-t31+t27;
  values[7] = -210.0*xi-t22+t34+t35+t36-t18+18.0-t31-t37;
  values[8] = t1+t3-t5-t7-t9-t11-27.0+t12-t37;
  values[9] = 330.0*xi+t22-t34-t35-t36+t18-78.0+t26+t13;
  values[10] = -t1+t22+t5+t7+t9-t11+27.0-t29-t23;
  values[11] = -t1-t15+t5+t16+t17+t18+27.0-t19-t20;
  values[12] = -40.0*xi+20.0-20.0*eta-20.0*zeta;
}

static void N_T_P2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15, t16;
  double t17, t18, t19, t20, t21, t22, t23, t26, t27, t28, t29, t32, t35;
  double t36, t38, t39;

  t1 = 135.0*eta;
  t2 = zeta*zeta;
  t3 = 135.0*t2;
  t4 = eta*eta;
  t5 = 135.0*t4;
  t6 = xi*xi;
  t7 = 90.0*t6;
  t8 = xi*eta;
  t9 = 180.0*t8;
  t10 = eta*zeta;
  t11 = 90.0*t10;
  t12 = 90.0*xi;
  t13 = 75.0*zeta;
  t15 = 45.0*t2;
  t16 = 45.0*t6;
  t17 = 90.0*t8;
  t18 = 45.0*xi;
  t19 = 45.0*zeta;
  t20 = t15-t16-t17+t11+t18-t19;
  t21 = 90.0*t2;
  t22 = 180.0*t10;
  t23 = 30.0*zeta;
  t26 = 270.0*t4;
  t27 = 270.0*t8;
  t28 = 270.0*t10;
  t29 = 75.0*xi;
  t32 = 30.0*xi;
  t35 = 15.0*xi;
  t36 = 15.0*zeta;
  t38 = 135.0*t6;
  t39 = 90.0*zeta;

  values[0] = -t1-t3+t5+t7+t9+t11+15.0-t12+t13;
  values[1] = t20;
  values[2] = t1+t21-t5-t16-t17-t22-15.0+t18-t23;
  values[3] = 330.0*eta+t15-t26+t16-t27-t28-78.0+t29+t13;
  values[4] = -t1-t21+t5+t16+t17+t22+27.0-t18-t23;
  values[5] = -t1+t15+t5-t7+t9+t11+27.0-t32-t19;
  values[6] = -210.0*eta-t15+t26-t16+t27+t28+18.0-t35-t36;
  values[7] = t1-t21-t5+t38-t17-t22-27.0-t35+t39;
  values[8] = t1+t3-t5-t7-t9-t11-27.0+t12-t36;
  values[9] = -t1+t21+t5-t38+t17+t22+15.0+t29-t39;
  values[10] = t1-t15-t5+t7-t9-t11-15.0-t32+t19;
  values[11] = -t20;
  values[12] = -40.0*eta+20.0-20.0*xi-20.0*zeta;
}

static void N_T_P2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15, t16, t17;
  double t18, t19, t20, t21, t23, t24, t25, t26, t27, t29, t30, t33, t34;
  double t36, t37, t38;

  t2 = zeta*zeta;
  t3 = 270.0*t2;
  t4 = xi*zeta;
  t5 = 270.0*t4;
  t6 = eta*zeta;
  t7 = 270.0*t6;
  t8 = xi*xi;
  t9 = 45.0*t8;
  t10 = eta*eta;
  t11 = 45.0*t10;
  t12 = 75.0*eta;
  t13 = 75.0*xi;
  t15 = 135.0*zeta;
  t16 = 135.0*t2;
  t17 = 180.0*t4;
  t18 = 90.0*t6;
  t19 = 90.0*t8;
  t20 = 45.0*eta;
  t21 = 30.0*xi;
  t23 = 90.0*t4;
  t24 = 180.0*t6;
  t25 = 90.0*t10;
  t26 = 30.0*eta;
  t27 = 45.0*xi;
  t29 = 135.0*t10;
  t30 = 90.0*xi;
  t33 = -t23+t18-t9+t11-t20+t27;
  t34 = 15.0*eta;
  t36 = 135.0*t8;
  t37 = 90.0*eta;
  t38 = 15.0*xi;

  values[0] = 330.0*zeta-t3-t5-t7+t9+t11+t12+t13-78.0;
  values[1] = -t15+t16+t17+t18-t19+t11-t20-t21+27.0;
  values[2] = -t15+t16+t23+t24+t9-t25-t26-t27+27.0;
  values[3] = -t15+t16+t17+t18+t19-t29+t12-t30+15.0;
  values[4] = t15-t16-t23-t24-t9+t25-t26+t27-15.0;
  values[5] = t33;
  values[6] = t15-t16-t17-t18-t19+t29-t34+t30-27.0;
  values[7] = t15-t16-t23-t24+t36-t25+t37-t38-27.0;
  values[8] = -210.0*zeta+t3+t5+t7-t9-t11-t34-t38+18.0;
  values[9] = -t15+t16+t23+t24-t36+t25-t37+t13+15.0;
  values[10] = -t33;
  values[11] = t15-t16-t17-t18+t19-t11+t20-t21-15.0;
  values[12] = -40.0*zeta-20.0*eta-20.0*xi+20.0;
}

static void N_T_P2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t4, t7, t8, t9, t10, t11;

  t1 = 270.0*xi;
  t4 = -135.0+t1+180.0*eta+90.0*zeta;
  t7 = 135.0-t1-90.0*eta-180.0*zeta;
  t8 = -eta+zeta;
  t9 = 540.0*xi;
  t10 = 270.0*eta;
  t11 = 270.0*zeta;

  values[0] = t4;
  values[1] = t7;
  values[2] = 90.0*t8;
  values[3] = -t7;
  values[4] = -90.0*t8;
  values[5] = -t4;
  values[6] = t7;
  values[7] = -210.0+t9+t10+t11;
  values[8] = -t4;
  values[9] = 330.0-t9-t10-t11;
  values[10] = t4;
  values[11] = -t7;
  values[12] = -40.0;
}

static void N_T_P2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t11;

  t1 = 180.0*xi;
  t2 = 180.0*eta;
  t3 = t1+t2-90.0;
  t4 = 90.0*xi;
  t5 = 90.0*eta;
  t6 = -t4-t5+45.0;
  t7 = 270.0*eta;
  t11 = 270.0*xi;

  values[0] = t3;
  values[1] = t6;
  values[2] = t6;
  values[3] = t4-t7+75.0;
  values[4] = -t6;
  values[5] = -t1+t2-30.0;
  values[6] = -t4+t7-15.0;
  values[7] = t11-t5-15.0;
  values[8] = -t3;
  values[9] = -t11+t5+75.0;
  values[10] = t1-t2-30.0;
  values[11] = -t6;
  values[12] = -20.0;
}

static void N_T_P2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t5, t7, t8, t9, t10;

  t1 = 270.0*zeta;
  t2 = 90.0*xi;
  t4 = 180.0*zeta;
  t5 = 180.0*xi;
  t7 = 90.0*zeta;
  t8 = t7+t2-45.0;
  t9 = t4+t5-90.0;
  t10 = 270.0*xi;

  values[0] = -t1+t2+75.0;
  values[1] = t4-t5-30.0;
  values[2] = t8;
  values[3] = t9;
  values[4] = -t8;
  values[5] = -t8;
  values[6] = -t9;
  values[7] = -t7+t10-15.0;
  values[8] = t1-t2-15.0;
  values[9] = t7-t10+75.0;
  values[10] = t8;
  values[11] = -t4+t5-30.0;
  values[12] = -20.0;
}

static void N_T_P2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t4, t5, t8, t9, t10, t11;

  t1 = 270.0*eta;
  t4 = -135.0+t1+180.0*xi+90.0*zeta;
  t5 = -xi+zeta;
  t8 = 135.0-t1-90.0*xi-180.0*zeta;
  t9 = 540.0*eta;
  t10 = 270.0*xi;
  t11 = 270.0*zeta;

  values[0] = t4;
  values[1] = 90.0*t5;
  values[2] = t8;
  values[3] = 330.0-t9-t10-t11;
  values[4] = -t8;
  values[5] = t4;
  values[6] = -210.0+t9+t10+t11;
  values[7] = t8;
  values[8] = -t4;
  values[9] = -t8;
  values[10] = -t4;
  values[11] = -90.0*t5;
  values[12] = -40.0;
}

static void N_T_P2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t5, t6, t7, t9, t13;

  t1 = 270.0*zeta;
  t2 = 90.0*eta;
  t4 = 90.0*zeta;
  t5 = t4+t2-45.0;
  t6 = 180.0*zeta;
  t7 = 180.0*eta;
  t9 = 270.0*eta;
  t13 = -t6-t7+90.0;

  values[0] = -t1+t2+75.0;
  values[1] = t5;
  values[2] = t6-t7-30.0;
  values[3] = t4-t9+75.0;
  values[4] = -t6+t7-30.0;
  values[5] = t5;
  values[6] = -t4+t9-15.0;
  values[7] = t13;
  values[8] = t1-t2-15.0;
  values[9] = -t13;
  values[10] = -t5;
  values[11] = -t5;
  values[12] = -20.0;
}

static void N_T_P2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t5, t8, t11, t12;

  t1 = 540.0*zeta;
  t2 = 270.0*xi;
  t3 = 270.0*eta;
  t5 = 270.0*zeta;
  t8 = -135.0+t5+180.0*xi+90.0*eta;
  t11 = -135.0+t5+90.0*xi+180.0*eta;
  t12 = -xi+eta;

  values[0] = 330.0-t1-t2-t3;
  values[1] = t8;
  values[2] = t11;
  values[3] = t8;
  values[4] = -t11;
  values[5] = 90.0*t12;
  values[6] = -t8;
  values[7] = -t11;
  values[8] = -210.0+t1+t2+t3;
  values[9] = t11;
  values[10] = -90.0*t12;
  values[11] = -t8;
  values[12] = -40.0;
}

TBaseFunct3D *BF_N_T_P2_3D_Obj = 
new TBaseFunct3D(13, BF_N_T_P2_3D, BFUnitTetrahedron, 
                 N_T_P2_3D_Funct, N_T_P2_3D_DeriveXi,
                 N_T_P2_3D_DeriveEta, N_T_P2_3D_DeriveZeta,
                 N_T_P2_3D_DeriveXiXi, N_T_P2_3D_DeriveXiEta,
                 N_T_P2_3D_DeriveXiZeta, N_T_P2_3D_DeriveEtaEta,
                 N_T_P2_3D_DeriveEtaZeta, N_T_P2_3D_DeriveZetaZeta,
                 3, 2,
                 0, NULL);
