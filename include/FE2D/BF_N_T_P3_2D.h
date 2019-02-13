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
// P3 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_T_P3_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t7, t8, t9, t10, t12, t13, t14, t15, t16;
  double t17, t21, t23, t26, t28, t29, t30, t31, t32, t33, t35, t36, t37;
  double t45, t46, t48, t49, t53, t55, t58, t59, t64, t70, t72, t73, t83;
  double t87, t98, t101, t105, t106, t108, t109, t112, t116, t119, t127;

  t1 = 12.0*eta;
  t2 = eta*eta;
  t3 = 30.0*t2;
  t4 = t2*eta;
  t5 = 20.0*t4;
  t7 = 12.0*xi;
  t8 = xi*xi;
  t9 = 30.0*t8;
  t10 = xi*eta;
  t12 = t8*xi;
  t13 = 20.0*t12;
  t14 = t8*eta;
  t15 = 60.0*t14;
  t16 = xi*t2;
  t17 = 60.0*t16;
  t21 = 64.0/5.0*eta;
  t23 = 6.0*t10;
  t26 = 48.0*t14;
  t28 = 64.0*t4;
  t29 = t8*t2;
  t30 = 42.0*t29;
  t31 = t12*eta;
  t32 = 28.0*t31;
  t33 = xi*t4;
  t35 = t2*t2;
  t36 = 28.0*t35;
  t37 = -4.0/5.0-2.0/5.0*xi+t21+6.0*t8-t23-48.0*t2-4.0*t12-t26+42.0*t16+t28+t30+t32-42.0*t33-t36;
  t45 = t8*t8;
  t46 = 28.0*t45;
  t48 = 16.0/5.0*xi-16.0/5.0*eta-24.0*t8+24.0*t2+48.0*t12+t15-t17-48.0*t4-70.0*t31-t46+70.0*t33+t36;
  t49 = 64.0/5.0*xi;
  t53 = 64.0*t12;
  t55 = 48.0*t16;
  t58 = 28.0*t33;
  t59 = 4.0/5.0-t49+2.0/5.0*eta+48.0*t8+t23-6.0*t2-t53-42.0*t14+t55+4.0*t4+42.0*t31+t46-t30-t58;
  t64 = 102.0/5.0*t10;
  t70 = 126.0/5.0*t29;
  t72 = 84.0/5.0*t35;
  t73 = 18.0/5.0*xi-204.0/25.0*eta+t3+336.0/5.0*t12-196.0/5.0*t4+t64+216.0/5.0*t14-414.0/5.0*t16-186.0/5.0*t8-168.0/5.0*t45-336.0/5.0*t31+t70+294.0/5.0*t33+t72+13.0/25.0;
  t83 = 84.0/5.0*t45;
  t87 = -36.0/25.0*xi-36.0/25.0*eta+66.0/5.0*t2-28.0*t12-28.0*t4-24.0/5.0*t10+48.0/5.0*t14+48.0/5.0*t16+66.0/5.0*t8+t83+42.0/5.0*t31-252.0/5.0*t29+42.0/5.0*t33+t72-1.0/25.0;
  t98 = -204.0/25.0*xi+18.0/5.0*eta-186.0/5.0*t2-196.0/5.0*t12+336.0/5.0*t4+t64-414.0/5.0*t14+216.0/5.0*t16+t9+t83+294.0/5.0*t31+t70-336.0/5.0*t33-168.0/5.0*t35+13.0/25.0;
  t101 = 144.0*t2;
  t105 = 144.0*t8;
  t106 = 56.0*t45;
  t108 = 56.0*t35;
  t109 = 224.0/5.0*xi+224.0/5.0*eta-t101+160.0*t12+160.0*t4-84.0*t10+t26+t55-t105-t106-t32+168.0*t29-t58-t108-16.0/5.0;
  t112 = 72.0*t10;
  t116 = 84.0*t29;
  t119 = -t49-32.0*eta+t101+t53-224.0*t4+t112+156.0*t14-264.0*t16+8.0/5.0-t106-196.0*t31-t116+224.0*t33+112.0*t35;
  t127 = -32.0*xi-t21-224.0*t12+t28+t112-264.0*t14+156.0*t16+t105+8.0/5.0+112.0*t45+224.0*t31-t116-196.0*t33-t108;

  values[0] = 1.0-t1+t3-t5;
  values[1] = -1.0+t7+t1-t9-60.0*t10-t3+t13+t15+t17+t5;
  values[2] = 1.0-t7+t9-t13;
  values[3] = t37;
  values[4] = t48;
  values[5] = t59;
  values[6] = t73;
  values[7] = t87;
  values[8] = t98;
  values[9] = t109;
  values[10] = t119;
  values[11] = t127;
}

// values of the derivatives in xi direction
static void N_T_P3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t3, t4, t5, t6, t7, t8, t12, t14, t16, t17, t18, t19, t20;
  double t26, t27, t31, t33, t35, t38, t44, t52, t65, t66, t69, t73;

  t1 = 60.0*xi;
  t3 = xi*xi;
  t4 = 60.0*t3;
  t5 = xi*eta;
  t6 = 120.0*t5;
  t7 = eta*eta;
  t8 = 60.0*t7;
  t12 = 6.0*eta;
  t14 = 96.0*t5;
  t16 = xi*t7;
  t17 = 84.0*t16;
  t18 = t3*eta;
  t19 = 84.0*t18;
  t20 = t7*eta;
  t26 = t3*xi;
  t27 = 112.0*t26;
  t31 = 192.0*t3;
  t33 = 48.0*t7;
  t35 = 28.0*t20;
  t38 = 102.0/5.0*eta;
  t44 = 252.0/5.0*t16;
  t52 = 336.0/5.0*t26;
  t65 = 288.0*xi;
  t66 = 224.0*t26;
  t69 = 72.0*eta;
  t73 = 168.0*t16;

  values[0] = 0.0;
  values[1] = 12.0-t1-60.0*eta+t4+t6+t8;
  values[2] = -12.0+t1-t4;
  values[3] = -2.0/5.0+12.0*xi-t12-12.0*t3-t14+42.0*t7+t17+t19-42.0*t20;
  values[4] = 16.0/5.0-48.0*xi+144.0*t3+t6-t8-210.0*t18-t27+70.0*t20;
  values[5] = -64.0/5.0+96.0*xi+t12-t31-84.0*t5+t33+126.0*t18+t27-t17-t35;
  values[6] = 18.0/5.0+1008.0/5.0*t3+t38+432.0/5.0*t5-414.0/5.0*t7-372.0/5.0*xi-672.0/5.0*t26-1008.0/5.0*t18+t44+294.0/5.0*t20;
  values[7] = -36.0/25.0-84.0*t3-24.0/5.0*eta+96.0/5.0*t5+48.0/5.0*t7+132.0/5.0*xi+t52+126.0/5.0*t18-504.0/5.0*t16+42.0/5.0*t20;
  values[8] = -204.0/25.0-588.0/5.0*t3+t38-828.0/5.0*t5+216.0/5.0*t7+t1+t52+882.0/5.0*t18+t44-336.0/5.0*t20;
  values[9] = 224.0/5.0+480.0*t3-84.0*eta+t14+t33-t65-t66-t19+336.0*t16-t35;
  values[10] = -64.0/5.0+t31+t69+312.0*t5-264.0*t7-t66-588.0*t18-t73+224.0*t20;
  values[11] = -32.0-672.0*t3+t69-528.0*t5+156.0*t7+t65+448.0*t26+672.0*t18-t73-196.0*t20;
}

// values of the derivatives in eta direction
static void N_T_P3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t6, t7, t8, t9, t11, t13, t15, t16, t17, t18, t19;
  double t20, t22, t23, t32, t35, t38, t42, t44, t63, t67, t70, t74;

  t1 = 60.0*eta;
  t2 = eta*eta;
  t3 = 60.0*t2;
  t6 = xi*xi;
  t7 = 60.0*t6;
  t8 = xi*eta;
  t9 = 120.0*t8;
  t11 = 6.0*xi;
  t13 = 48.0*t6;
  t15 = 192.0*t2;
  t16 = t6*eta;
  t17 = 84.0*t16;
  t18 = t6*xi;
  t19 = 28.0*t18;
  t20 = xi*t2;
  t22 = t2*eta;
  t23 = 112.0*t22;
  t32 = 96.0*t8;
  t35 = 84.0*t20;
  t38 = 102.0/5.0*xi;
  t42 = 252.0/5.0*t16;
  t44 = 336.0/5.0*t22;
  t63 = 288.0*eta;
  t67 = 224.0*t22;
  t70 = 72.0*xi;
  t74 = 168.0*t16;

  values[0] = -12.0+t1-t3;
  values[1] = 12.0-60.0*xi-t1+t7+t9+t3;
  values[2] = 0.0;
  values[3] = 64.0/5.0-t11-96.0*eta-t13+84.0*t8+t15+t17+t19-126.0*t20-t23;
  values[4] = -16.0/5.0+48.0*eta+t7-t9-144.0*t2-70.0*t18+210.0*t20+t23;
  values[5] = 2.0/5.0+t11-12.0*eta-42.0*t6+t32+12.0*t2+42.0*t18-t17-t35;
  values[6] = -204.0/25.0+t1-588.0/5.0*t2+t38+216.0/5.0*t6-828.0/5.0*t8-336.0/5.0*t18+t42+882.0/5.0*t20+t44;
  values[7] = -36.0/25.0+132.0/5.0*eta-84.0*t2-24.0/5.0*xi+48.0/5.0*t6+96.0/5.0*t8+42.0/5.0*t18-504.0/5.0*t16+126.0/5.0*t20+t44;
  values[8] = 18.0/5.0-372.0/5.0*eta+1008.0/5.0*t2+t38-414.0/5.0*t6+432.0/5.0*t8+294.0/5.0*t18+t42-1008.0/5.0*t20-672.0/5.0*t22;
  values[9] = 224.0/5.0-t63+480.0*t2-84.0*xi+t13+t32-t19+336.0*t16-t35-t67;
  values[10] = -32.0+t63-672.0*t2+t70+156.0*t6-528.0*t8-196.0*t18-t74+672.0*t20+448.0*t22;
  values[11] = -64.0/5.0+t15+t70-264.0*t6+312.0*t8+224.0*t18-t74-588.0*t20-t67;
}

// values of derivatives in xi-xi direction
static void N_T_P3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2, t6, t7, t8, t9, t10, t14, t15, t17, t25, t29, t38, t43;

  t1 = 120.0*xi;
  t2 = 120.0*eta;
  t6 = 96.0*eta;
  t7 = eta*eta;
  t8 = 84.0*t7;
  t9 = xi*eta;
  t10 = 168.0*t9;
  t14 = xi*xi;
  t15 = 336.0*t14;
  t17 = 384.0*xi;
  t25 = 252.0/5.0*t7;
  t29 = 1008.0/5.0*t14;
  t38 = 672.0*t14;
  t43 = 168.0*t7;

  values[0] = 0.0;
  values[1] = -60.0+t1+t2;
  values[2] = 60.0-t1;
  values[3] = 12.0-24.0*xi-t6+t8+t10;
  values[4] = -48.0+288.0*xi+t2-420.0*t9-t15;
  values[5] = 96.0-t17-84.0*eta+252.0*t9+t15-t8;
  values[6] = 2016.0/5.0*xi+432.0/5.0*eta-372.0/5.0-2016.0/5.0*t14-2016.0/5.0*t9+t25;
  values[7] = -168.0*xi+96.0/5.0*eta+132.0/5.0+t29+252.0/5.0*t9-504.0/5.0*t7;
  values[8] = -1176.0/5.0*xi-828.0/5.0*eta+60.0+t29+1764.0/5.0*t9+t25;
  values[9] = 960.0*xi+t6-288.0-t38-t10+336.0*t7;
  values[10] = t17+312.0*eta-t38-1176.0*t9-t43;
  values[11] = -1344.0*xi-528.0*eta+288.0+1344.0*t14+1344.0*t9-t43;
}

// values of derivatives in eta-eta direction
static void N_T_P3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t4, t6, t7, t8, t9, t10, t17, t19, t24, t43;

  t1 = 120.0*xi;
  t2 = 120.0*eta;
  t4 = 96.0*xi;
  t6 = xi*eta;
  t7 = 168.0*t6;
  t8 = xi*xi;
  t9 = 84.0*t8;
  t10 = eta*eta;
  t17 = 96.0*eta;
  t19 = 84.0*t10;
  t24 = 504.0/5.0*t6;
  t43 = 336.0*t6;

  values[0] = 0.0;
  values[1] = -60.0+t1+t2;
  values[2] = 0.0;
  values[3] = -6.0-t4+84.0*eta+t7+t9-126.0*t10;
  values[4] = t1-t2-210.0*t8+210.0*t10;
  values[5] = 6.0-84.0*xi+t17+126.0*t8-t7-t19;
  values[6] = 102.0/5.0+432.0/5.0*xi-828.0/5.0*eta-1008.0/5.0*t8+t24+882.0/5.0*t10;
  values[7] = -24.0/5.0+96.0/5.0*xi+96.0/5.0*eta+126.0/5.0*t8-1008.0/5.0*t6+126.0/5.0*t10;
  values[8] = 102.0/5.0-828.0/5.0*xi+432.0/5.0*eta+882.0/5.0*t8+t24-1008.0/5.0*t10;
  values[9] = -84.0+t4+t17-t9+672.0*t6-t19;
  values[10] = 72.0+312.0*xi-528.0*eta-588.0*t8-t43+672.0*t10;
  values[11] = 72.0-528.0*xi+312.0*eta+672.0*t8-t43-588.0*t10;
}

// values of derivatives in xi-eta direction
static void N_T_P3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t3, t6, t7, t8, t9, t11, t12, t17, t19, t23, t25, t39, t43;

  t1 = 120.0*eta;
  t3 = 120.0*xi;
  t6 = 384.0*eta;
  t7 = xi*xi;
  t8 = 84.0*t7;
  t9 = xi*eta;
  t11 = eta*eta;
  t12 = 336.0*t11;
  t17 = 96.0*xi;
  t19 = 168.0*t9;
  t23 = 252.0/5.0*t7;
  t25 = 1008.0/5.0*t11;
  t39 = 672.0*t11;
  t43 = 168.0*t7;

   values[0] = 60.0-t1;
   values[1] = -60.0+t3+t1;
   values[2] = 0.0;
   values[3] = -96.0+84.0*xi+t6+t8-252.0*t9-t12;
   values[4] = 48.0-t3-288.0*eta+420.0*t9+t12;
   values[5] = -12.0+t17+24.0*eta-t8-t19;
   values[6] = 60.0-1176.0/5.0*eta-828.0/5.0*xi+t23+1764.0/5.0*t9+t25;
   values[7] = 132.0/5.0-168.0*eta+96.0/5.0*xi-504.0/5.0*t7+252.0/5.0*t9+t25;
   values[8] = -372.0/5.0+2016.0/5.0*eta+432.0/5.0*xi+t23-2016.0/5.0*t9-2016.0/5.0*t11;
   values[9] = -288.0+960.0*eta+t17+336.0*t7-t19-t39;
   values[10] = 288.0-1344.0*eta-528.0*xi-t43+1344.0*t9+1344.0*t11;
   values[11] = t6+312.0*xi-t43-1176.0*t9-t39;
}

static int N_T_P3_2D_ChangeJ0[1] = { 3 };
static int N_T_P3_2D_ChangeJ1[1] = { 4 };
static int N_T_P3_2D_ChangeJ2[1] = { 5 };

static int *N_T_P3_2D_Change[3] = { N_T_P3_2D_ChangeJ0, N_T_P3_2D_ChangeJ1,
                                    N_T_P3_2D_ChangeJ2 };

// ***********************************************************************

TBaseFunct2D *BF_N_T_P3_2D_Obj = new TBaseFunct2D
        (12, BF_N_T_P3_2D, BFUnitTriangle, 
         N_T_P3_2D_Funct, N_T_P3_2D_DeriveXi,
         N_T_P3_2D_DeriveEta, N_T_P3_2D_DeriveXiXi,
         N_T_P3_2D_DeriveXiEta, N_T_P3_2D_DeriveEtaEta, 4, 3,
         1, N_T_P3_2D_Change);
