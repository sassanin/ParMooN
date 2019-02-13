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
// Q4 element with bubbles, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_M4_2D_Funct(double xi, double eta, double *values)
{
  double t1 = 4.0/9.0*xi;
  double t2 = 4.0/9.0*eta;
  double t3 = xi*xi;
  double t4 = 7.0/9.0*t3;
  double t5 = xi*eta;
  double t6 = 5.0/12.0*t5;
  double t7 = eta*eta;
  double t8 = 7.0/9.0*t7;
  double t9 = t3*xi;
  double t10 = t9/3.0;
  double t11 = t3*eta;
  double t12 = 13.0/36.0*t11;
  double t13 = xi*t7;
  double t14 = 13.0/36.0*t13;
  double t15 = t7*eta;
  double t16 = t15/3.0;
  double t17 = t3*t3;
  double t18 = 2.0/3.0*t17;
  double t19 = t9*eta;
  double t20 = t19/3.0;
  double t21 = t3*t7;
  double t22 = t21/4.0;
  double t23 = xi*t15;
  double t24 = t23/3.0;
  double t25 = t7*t7;
  double t26 = 2.0/3.0*t25;
  double t32 = (1.0-t3)*(-1.0/2.0+3.0/2.0*t3)*(1.0+eta);
  double t33 = 2.0/9.0*t32;
  double t39 = (1.0-t7)*(-1.0/2.0+3.0/2.0*t7)*(1.0+xi);
  double t40 = 2.0/9.0*t39;
  double t41 = 2.0/9.0+t1+t2-t4-t6-t8-t10-t12-t14-t16+t18+t20+t22+t24+t26+t33+t40;
  double t42 = 2.0/3.0*xi;
  double t43 = 28.0/9.0*t3;
  double t44 = 2.0/3.0*t5;
  double t45 = 2.0/3.0*t9;
  double t46 = 4.0/9.0*t11;
  double t47 = 8.0/3.0*t17;
  double t48 = 2.0/3.0*t19;
  double t49 = 8.0/9.0*t32;
  double t51 = eta/6.0;
  double t53 = t7/2.0;
  double t54 = t11/6.0;
  double t56 = t21/2.0;
  double t57 = 4.0/3.0*t32;
  double t60 = t7/9.0;
  double t61 = -t1+t2-t4+t6+t60+t10-t12+t14-t16+t18-t20+t22-t24+t33-t40;
  double t62 = 2.0/3.0*eta;
  double t63 = 4.0/9.0*t7;
  double t64 = 4.0/9.0*t13;
  double t65 = 2.0/3.0*t15;
  double t66 = 2.0/3.0*t23;
  double t67 = 8.0/9.0*t39;
  double t69 = xi/6.0;
  double t70 = t3/2.0;
  double t72 = t13/6.0;
  double t73 = 4.0/3.0*t39;
  double t76 = t3/9.0;
  double t77 = -2.0/9.0-t1-t2+t76-t6+t60+t10+t12+t14+t16+t20+t22+t24-t33-t40;
  double t78 = 4.0/9.0*t3;
  double t83 = t1-t2+t76+t6-t8-t10+t12-t14+t16-t20+t22-t24+t26-t33+t40;
  double t84 = 28.0/9.0*t7;
  double t85 = 8.0/3.0*t25;

  values[0] = t41;
  values[1] = -4.0/9.0-t42-t2+t43+t44+t45+t46-t47-t48-t49;
  values[2] = 2.0/3.0+t51-14.0/3.0*t3+t53-t54+4.0*t17-t56+t57;
  values[3] = -4.0/9.0+t42-t2+t43-t44-t45+t46-t47+t48-t49;
  values[4] = t61;
  values[5] = 4.0/9.0+t1-t62-t44-t63-t64+t65+t66+t67;
  values[6] = -2.0/3.0-t69+t70+2.0/3.0*t7+t72-t56-t73;
  values[7] = 4.0/9.0+t1+t62+t44-t63-t64-t65-t66+t67;
  values[8] = t77;
  values[9] = 4.0/9.0+t42+t2-t78+t44-t45-t46-t48+t49;
  values[10] = -2.0/3.0-t51+2.0/3.0*t3+t53+t54-t56-t57;
  values[11] = 4.0/9.0-t42+t2-t78-t44+t45-t46+t48+t49;
  values[12] = t83;
  values[13] = -4.0/9.0-t1+t62-t44+t84+t64-t65+t66-t85-t67;
  values[14] = 2.0/3.0+t69+t70-14.0/3.0*t7-t72-t56+4.0*t25+t73;
  values[15] = -4.0/9.0-t1-t62+t44+t84+t64+t65-t66-t85-t67;
  values[16] = 1.0-t3-t7+t21;
}

// values of the derivatives in xi direction
static void C_Q_M4_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1 = 14.0/9.0*xi;
  double t2 = 5.0/12.0*eta;
  double t3 = xi*xi;
  double t4 = xi*eta;
  double t5 = 13.0/18.0*t4;
  double t6 = eta*eta;
  double t7 = 13.0/36.0*t6;
  double t8 = t3*xi;
  double t9 = 8.0/3.0*t8;
  double t10 = t3*eta;
  double t11 = xi*t6;
  double t12 = t11/2.0;
  double t13 = t6*eta;
  double t14 = t13/3.0;
  double t18 = 1.0+eta;
  double t19 = xi*(-1.0/2.0+3.0/2.0*t3)*t18;
  double t20 = 4.0/9.0*t19;
  double t23 = (1.0-t3)*xi*t18;
  double t24 = 2.0/3.0*t23;
  double t28 = (1.0-t6)*(-1.0/2.0+3.0/2.0*t6);
  double t29 = 2.0/9.0*t28;
  double t30 = 4.0/9.0-t1-t2-t3-t5-t7+t9+t10+t12+t14-t20+t24+t29;
  double t31 = 56.0/9.0*xi;
  double t32 = 2.0/3.0*eta;
  double t33 = 2.0*t3;
  double t34 = 8.0/9.0*t4;
  double t35 = 32.0/3.0*t8;
  double t36 = 2.0*t10;
  double t37 = 16.0/9.0*t19;
  double t38 = 8.0/3.0*t23;
  double t41 = t4/3.0;
  double t43 = 8.0/3.0*t19;
  double t44 = 4.0*t23;
  double t47 = -4.0/9.0-t1+t2+t3-t5+t7+t9-t10+t12-t14-t20+t24-t29;
  double t48 = 4.0/9.0*t6;
  double t49 = 2.0/3.0*t13;
  double t50 = 8.0/9.0*t28;
  double t51 = 4.0/9.0-t32-t48+t49+t50;
  double t52 = t6/6.0;
  double t53 = 4.0/3.0*t28;
  double t55 = 4.0/9.0+t32-t48-t49+t50;
  double t56 = 2.0/9.0*xi;
  double t57 = -4.0/9.0+t56-t2+t3+t5+t7+t10+t12+t14+t20-t24-t29;
  double t58 = 8.0/9.0*xi;
  double t63 = 4.0/9.0+t56+t2-t3+t5-t7-t10+t12-t14+t20-t24+t29;

  values[0] = t30;
  values[1] = -2.0/3.0+t31+t32+t33+t34-t35-t36+t37-t38;
  values[2] = -28.0/3.0*xi-t41+16.0*t8-t11-t43+t44;
  values[3] = 2.0/3.0+t31-t32-t33+t34-t35+t36+t37-t38;
  values[4] = t47;
  values[5] = t51;
  values[6] = -1.0/6.0+xi+t52-t11-t53;
  values[7] = t55;
  values[8] = t57;
  values[9] = 2.0/3.0-t58+t32-t33-t34-t36-t37+t38;
  values[10] = 4.0/3.0*xi+t41-t11+t43-t44;
  values[11] = -2.0/3.0-t58-t32+t33-t34+t36-t37+t38;
  values[12] = t63;
  values[13] = -t55;
  values[14] = 1.0/6.0+xi-t52-t11+t53;
  values[15] = -t51;
  values[16] = -2.0*xi+2.0*t11;
}

// values of the derivatives in eta direction
static void C_Q_M4_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1 = 5.0/12.0*xi;
  double t2 = 14.0/9.0*eta;
  double t3 = xi*xi;
  double t4 = 13.0/36.0*t3;
  double t5 = xi*eta;
  double t6 = 13.0/18.0*t5;
  double t7 = eta*eta;
  double t8 = t3*xi;
  double t9 = t8/3.0;
  double t10 = t3*eta;
  double t11 = t10/2.0;
  double t12 = xi*t7;
  double t13 = t7*eta;
  double t14 = 8.0/3.0*t13;
  double t18 = (1.0-t3)*(-1.0/2.0+3.0/2.0*t3);
  double t19 = 2.0/9.0*t18;
  double t23 = 1.0+xi;
  double t24 = eta*(-1.0/2.0+3.0/2.0*t7)*t23;
  double t25 = 4.0/9.0*t24;
  double t28 = (1.0-t7)*eta*t23;
  double t29 = 2.0/3.0*t28;
  double t30 = 4.0/9.0-t1-t2-t4-t6-t7+t9+t11+t12+t14+t19-t25+t29;
  double t31 = 2.0/3.0*xi;
  double t32 = 4.0/9.0*t3;
  double t33 = 2.0/3.0*t8;
  double t34 = 8.0/9.0*t18;
  double t35 = -4.0/9.0+t31+t32-t33-t34;
  double t36 = t3/6.0;
  double t37 = 4.0/3.0*t18;
  double t39 = -4.0/9.0-t31+t32+t33-t34;
  double t40 = 2.0/9.0*eta;
  double t41 = 4.0/9.0+t1+t40-t4+t6-t7-t9+t11-t12+t19+t25-t29;
  double t42 = 8.0/9.0*eta;
  double t43 = 8.0/9.0*t5;
  double t44 = 2.0*t7;
  double t45 = 2.0*t12;
  double t46 = 16.0/9.0*t24;
  double t47 = 8.0/3.0*t28;
  double t50 = t5/3.0;
  double t51 = 8.0/3.0*t24;
  double t52 = 4.0*t28;
  double t55 = -4.0/9.0-t1+t40+t4+t6+t7+t9+t11+t12-t19+t25-t29;
  double t57 = -4.0/9.0+t1-t2+t4-t6+t7-t9+t11-t12+t14-t19-t25+t29;
  double t58 = 56.0/9.0*eta;
  double t59 = 32.0/3.0*t13;

  values[0] = t30;
  values[1] = t35;
  values[2] = 1.0/6.0+eta-t36-t10+t37;
  values[3] = t39;
  values[4] = t41;
  values[5] = -2.0/3.0-t31-t42-t43+t44+t45-t46+t47;
  values[6] = 4.0/3.0*eta+t50-t10+t51-t52;
  values[7] = 2.0/3.0+t31-t42-t43-t44-t45-t46+t47;
  values[8] = t55;
  values[9] = -t39;
  values[10] = -1.0/6.0+eta+t36-t10-t37;
  values[11] = -t35;
  values[12] = t57;
  values[13] = 2.0/3.0-t31+t58+t43-t44+t45-t59+t46-t47;
  values[14] = -28.0/3.0*eta-t50-t10+16.0*t13-t51+t52;
  values[15] = -2.0/3.0+t31+t58+t43+t44-t45-t59+t46-t47;
  values[16] = -2.0*eta+2.0*t10;
}

// values of the derivatives in xi-xi  direction
static void C_Q_M4_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1 = 2.0*xi;
  double t2 = 13.0/18.0*eta;
  double t3 = xi*xi;
  double t4 = 8.0*t3;
  double t5 = xi*eta;
  double t6 = 2.0*t5;
  double t7 = eta*eta;
  double t8 = t7/2.0;
  double t11 = 1.0+eta;
  double t12 = (-1.0/2.0+3.0/2.0*t3)*t11;
  double t13 = 4.0/9.0*t12;
  double t14 = t3*t11;
  double t15 = 8.0/3.0*t14;
  double t17 = (1.0-t3)*t11;
  double t18 = 2.0/3.0*t17;
  double t20 = 4.0*xi;
  double t21 = 8.0/9.0*eta;
  double t22 = 32.0*t3;
  double t23 = 4.0*t5;
  double t24 = 16.0/9.0*t12;
  double t25 = 32.0/3.0*t14;
  double t26 = 8.0/3.0*t17;
  double t28 = eta/3.0;
  double t30 = 8.0/3.0*t12;
  double t31 = 16.0*t14;
  double t32 = 4.0*t17;
  double t36 = 1.0-t7;

  values[0] = -14.0/9.0-t1-t2+t4+t6+t8-t13-t15+t18;
  values[1] = 56.0/9.0+t20+t21-t22-t23+t24+t25-t26;
  values[2] = -28.0/3.0-t28+48.0*t3-t7-t30-t31+t32;
  values[3] = 56.0/9.0-t20+t21-t22+t23+t24+t25-t26;
  values[4] = -14.0/9.0+t1-t2+t4-t6+t8-t13-t15+t18;
  values[5] = 0.0;
  values[6] = t36;
  values[7] = 0.0;
  values[8] = 2.0/9.0+t1+t2+t6+t8+t13+t15-t18;
  values[9] = -8.0/9.0-t20-t21-t23-t24-t25+t26;
  values[10] = 4.0/3.0+t28-t7+t30+t31-t32;
  values[11] = -8.0/9.0+t20-t21+t23-t24-t25+t26;
  values[12] = 2.0/9.0-t1+t2-t6+t8+t13+t15-t18;
  values[13] = 0.0;
  values[14] = t36;
  values[15] = 0.0;
  values[16] = -2.0*t36;
}

// values of the derivatives in xi-eta direction
static void C_Q_M4_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = 13.0/18.0*xi;
  double t2 = 13.0/18.0*eta;
  double t3 = xi*xi;
  double t4 = xi*eta;
  double t5 = eta*eta;
  double t8 = xi*(-1.0/2.0+3.0/2.0*t3);
  double t9 = 4.0/9.0*t8;
  double t11 = (1.0-t3)*xi;
  double t12 = 2.0/3.0*t11;
  double t15 = eta*(-1.0/2.0+3.0/2.0*t5);
  double t16 = 4.0/9.0*t15;
  double t18 = (1.0-t5)*eta;
  double t19 = 2.0/3.0*t18;
  double t21 = 8.0/9.0*xi;
  double t22 = 2.0*t3;
  double t23 = 16.0/9.0*t8;
  double t24 = 8.0/3.0*t11;
  double t25 = 2.0/3.0+t21-t22+t23-t24;
  double t26 = xi/3.0;
  double t27 = 2.0*t4;
  double t28 = 8.0/3.0*t8;
  double t29 = 4.0*t11;
  double t31 = -2.0/3.0+t21+t22+t23-t24;
  double t33 = 8.0/9.0*eta;
  double t34 = 2.0*t5;
  double t35 = 16.0/9.0*t15;
  double t36 = 8.0/3.0*t18;
  double t37 = -2.0/3.0-t33+t34-t35+t36;
  double t38 = eta/3.0;
  double t39 = 8.0/3.0*t15;
  double t40 = 4.0*t18;
  double t42 = 2.0/3.0-t33-t34-t35+t36;

  values[0] = -5.0/12.0-t1-t2+t3+t4+t5-t9+t12-t16+t19;
  values[1] = t25;
  values[2] = -t26-t27-t28+t29;
  values[3] = t31;
  values[4] = 5.0/12.0-t1+t2-t3+t4-t5-t9+t12+t16-t19;
  values[5] = t37;
  values[6] = t38-t27+t39-t40;
  values[7] = t42;
  values[8] = -5.0/12.0+t1+t2+t3+t4+t5+t9-t12+t16-t19;
  values[9] = -t31;
  values[10] = t26-t27+t28-t29;
  values[11] = -t25;
  values[12] = 5.0/12.0+t1-t2-t3+t4-t5+t9-t12-t16+t19;
  values[13] = -t42;
  values[14] = -t38-t27-t39+t40;
  values[15] = -t37;
  values[16] = 4.0*t4;
}

// values of the derivatives in eta-eta direction
static void C_Q_M4_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1 = 13.0/18.0*xi;
  double t2 = 2.0*eta;
  double t3 = xi*xi;
  double t4 = t3/2.0;
  double t5 = xi*eta;
  double t6 = 2.0*t5;
  double t7 = eta*eta;
  double t8 = 8.0*t7;
  double t11 = 1.0+xi;
  double t12 = (-1.0/2.0+3.0/2.0*t7)*t11;
  double t13 = 4.0/9.0*t12;
  double t14 = t7*t11;
  double t15 = 8.0/3.0*t14;
  double t17 = (1.0-t7)*t11;
  double t18 = 2.0/3.0*t17;
  double t20 = 1.0-t3;
  double t22 = 8.0/9.0*xi;
  double t23 = 4.0*eta;
  double t24 = 4.0*t5;
  double t25 = 16.0/9.0*t12;
  double t26 = 32.0/3.0*t14;
  double t27 = 8.0/3.0*t17;
  double t29 = xi/3.0;
  double t30 = 8.0/3.0*t12;
  double t31 = 16.0*t14;
  double t32 = 4.0*t17;
  double t37 = 32.0*t7;

  values[0] = -14.0/9.0-t1-t2+t4+t6+t8-t13-t15+t18;
  values[1] = 0.0;
  values[2] = t20;
  values[3] = 0.0;
  values[4] = 2.0/9.0+t1-t2+t4-t6+t13+t15-t18;
  values[5] = -8.0/9.0-t22+t23+t24-t25-t26+t27;
  values[6] = 4.0/3.0+t29-t3+t30+t31-t32;
  values[7] = -8.0/9.0-t22-t23-t24-t25-t26+t27;
  values[8] = 2.0/9.0+t1+t2+t4+t6+t13+t15-t18;
  values[9] = 0.0;
  values[10] = t20;
  values[11] = 0.0;
  values[12] = -14.0/9.0-t1+t2+t4-t6+t8-t13-t15+t18;
  values[13] = 56.0/9.0+t22-t23+t24-t37+t25+t26-t27;
  values[14] = -28.0/3.0-t29-t3+48.0*t7-t30-t31+t32;
  values[15] = 56.0/9.0+t22+t23-t24-t37+t25+t26-t27;
  values[16] = -2.0*t20;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_M4_2D_Obj = new TBaseFunct2D
        (17, BF_C_Q_M4_2D, BFUnitSquare, 
         C_Q_M4_2D_Funct, C_Q_M4_2D_DeriveXi,
         C_Q_M4_2D_DeriveEta, C_Q_M4_2D_DeriveXiXi,
         C_Q_M4_2D_DeriveXiEta, C_Q_M4_2D_DeriveEtaEta, 4, 4,
         0, NULL);
