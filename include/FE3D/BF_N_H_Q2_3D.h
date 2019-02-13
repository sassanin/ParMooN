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
// Q2Rot element, nonconforming, 3D
// ***********************************************************************

static void N_H_Q2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t3, t4, t5, t6, t8, t9, t10, t11, t13, t15, t16, t17, t18;
  double t20, t22, t24, t25, t27, t29, t35, t39, t41, t45, t47, t50, t54;
  double t56, t58, t63, t65, t67;

  t1 = zeta/2.0;
  t2 = zeta*zeta;
  t3 = 3.0/4.0*t2;
  t4 = eta*eta;
  t5 = 3.0/2.0*t4;
  t6 = t5-1.0/2.0;
  t8 = t6*zeta/2.0;
  t9 = xi*xi;
  t10 = 3.0/2.0*t9;
  t11 = t10-1.0/2.0;
  t13 = t11*zeta/2.0;
  t15 = eta/2.0;
  t16 = 3.0/4.0*t4;
  t17 = 3.0/2.0*t2;
  t18 = t17-1.0/2.0;
  t20 = eta*t18/2.0;
  t22 = t11*eta/2.0;
  t24 = xi/2.0;
  t25 = 3.0/4.0*t9;
  t27 = xi*t18/2.0;
  t29 = xi*t6/2.0;
  t35 = xi*zeta/4.0;
  t39 = 5.0/2.0*t9*xi-3.0/2.0*xi;
  t41 = t39*zeta/4.0;
  t45 = 5.0/2.0*t2*zeta-3.0/2.0*zeta;
  t47 = xi*t45/4.0;
  t50 = eta*zeta/4.0;
  t54 = 5.0/2.0*t4*eta-3.0/2.0*eta;
  t56 = t54*zeta/4.0;
  t58 = eta*t45/4.0;
  t63 = xi*eta/4.0;
  t65 = xi*t54/4.0;
  t67 = t39*eta/4.0;

  values[0] = -1.0/4.0-t1+t3+t8+t13;
  values[1] = -1.0/4.0-t15+t16+t20+t22;
  values[2] = -1.0/4.0+t24+t25-t27-t29;
  values[3] = -1.0/4.0+t15+t16-t20-t22;
  values[4] = -1.0/4.0-t24+t25+t27+t29;
  values[5] = -1.0/4.0+t1+t3-t8-t13;
  values[6] = -t35+t27+t41-t47;
  values[7] = -t50+t8-t56+t58;
  values[8] = t35+t13+t41-t47;
  values[9] = t50+t8+t56-t58;
  values[10] = -t63+t22+t65-t67;
  values[11] = t50+t20-t56+t58;
  values[12] = -t50+t20+t56-t58;
  values[13] = -t63+t29-t65+t67;
  values[14] = t63+t22-t65+t67;
  values[15] = -t63-t29-t65+t67;
  values[16] = -t35+t13-t41+t47;
  values[17] = t35+t27-t41+t47;
  values[18] = 5.0/2.0-t10-t5-t17;
}

static void N_H_Q2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t4, t5, t6, t7, t8, t9, t12, t13, t15, t17, t19, t21, t23;
  double t25, t27, t29;

  t2 = 3.0/2.0*xi*zeta;
  t4 = 3.0/2.0*xi*eta;
  t5 = 3.0/2.0*xi;
  t6 = zeta*zeta;
  t7 = 3.0/4.0*t6;
  t8 = eta*eta;
  t9 = 3.0/4.0*t8;
  t12 = zeta/8.0;
  t13 = xi*xi;
  t15 = 15.0/2.0*t13-3.0/2.0;
  t17 = t15*zeta/4.0;
  t19 = 5.0/8.0*t6*zeta;
  t21 = 5.0/8.0*zeta;
  t23 = 5.0/8.0*eta;
  t25 = 5.0/8.0*t8*eta;
  t27 = t15*eta/4.0;
  t29 = eta/8.0;

  values[0] = t2;
  values[1] = t4;
  values[2] = 1.0+t5-t7-t9;
  values[3] = -t4;
  values[4] = -1.0+t5+t7+t9;
  values[5] = -t2;
  values[6] = t12+t7-1.0/4.0+t17-t19;
  values[7] = 0.0;
  values[8] = t21+t2+t17-t19;
  values[9] = 0.0;
  values[10] = -t23+t4+t25-t27;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = t29+t9-1.0/4.0-t25+t27;
  values[14] = t23+t4-t25+t27;
  values[15] = t29-t9+1.0/4.0-t25+t27;
  values[16] = -t21+t2-t17+t19;
  values[17] = -t12+t7-1.0/4.0-t17+t19;
  values[18] = -3.0*xi;
}

static void N_H_Q2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  double t2, t3, t4, t5, t6, t7, t10, t12, t13, t15, t17, t19, t22, t24;
  double t26, t28, t31;

  t2 = 3.0/2.0*eta*zeta;
  t3 = 3.0/2.0*eta;
  t4 = zeta*zeta;
  t5 = 3.0/4.0*t4;
  t6 = xi*xi;
  t7 = 3.0/4.0*t6;
  t10 = 3.0/2.0*xi*eta;
  t12 = 5.0/8.0*zeta;
  t13 = eta*eta;
  t15 = 15.0/2.0*t13-3.0/2.0;
  t17 = t15*zeta/4.0;
  t19 = 5.0/8.0*t4*zeta;
  t22 = xi/8.0;
  t24 = xi*t15/4.0;
  t26 = 5.0/8.0*t6*xi;
  t28 = zeta/8.0;
  t31 = 5.0/8.0*xi;

  values[0] = t2;
  values[1] = -1.0+t3+t5+t7;
  values[2] = -t10;
  values[3] = 1.0+t3-t5-t7;
  values[4] = t10;
  values[5] = -t2;
  values[6] = 0.0;
  values[7] = -t12+t2-t17+t19;
  values[8] = 0.0;
  values[9] = t12+t2+t17-t19;
  values[10] = t22+t7-1.0/4.0+t24-t26;
  values[11] = -t28+t5-1.0/4.0-t17+t19;
  values[12] = t28+t5-1.0/4.0+t17-t19;
  values[13] = -t31+t10-t24+t26;
  values[14] = -t22+t7-1.0/4.0-t24+t26;
  values[15] = -t31-t10-t24+t26;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = -3.0*eta;
}

static void N_H_Q2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t8, t10, t12, t14, t15, t17, t19, t21, t23;
  double t25, t27, t30;

  t1 = 3.0/2.0*zeta;
  t2 = eta*eta;
  t3 = 3.0/4.0*t2;
  t4 = xi*xi;
  t5 = 3.0/4.0*t4;
  t8 = 3.0/2.0*eta*zeta;
  t10 = 3.0/2.0*xi*zeta;
  t12 = 5.0/8.0*xi;
  t14 = 5.0/8.0*t4*xi;
  t15 = zeta*zeta;
  t17 = 15.0/2.0*t15-3.0/2.0;
  t19 = xi*t17/4.0;
  t21 = eta/8.0;
  t23 = 5.0/8.0*t2*eta;
  t25 = eta*t17/4.0;
  t27 = xi/8.0;
  t30 = 5.0/8.0*eta;

  values[0] = -1.0+t1+t3+t5;
  values[1] = t8;
  values[2] = -t10;
  values[3] = -t8;
  values[4] = t10;
  values[5] = 1.0+t1-t3-t5;
  values[6] = -t12+t10+t14-t19;
  values[7] = t21+t3-1.0/4.0-t23+t25;
  values[8] = -t27+t5-1.0/4.0+t14-t19;
  values[9] = -t21+t3-1.0/4.0+t23-t25;
  values[10] = 0.0;
  values[11] = t30+t8-t23+t25;
  values[12] = -t30+t8+t23-t25;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = t27+t5-1.0/4.0-t14+t19;
  values[17] = t12+t10-t14+t19;
  values[18] = -3.0*zeta;
}

static void N_H_Q2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t7;

  t1 = 3.0/2.0*zeta;
  t2 = 3.0/2.0*eta;
  t4 = 15.0/4.0*xi*zeta;
  t7 = 15.0/4.0*xi*eta;

  values[0] = t1;
  values[1] = t2;
  values[2] = 3.0/2.0;
  values[3] = -t2;
  values[4] = 3.0/2.0;
  values[5] = -t1;
  values[6] = t4;
  values[7] = 0.0;
  values[8] = t1+t4;
  values[9] = 0.0;
  values[10] = t2-t7;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = t7;
  values[14] = t2+t7;
  values[15] = t7;
  values[16] = t1-t4;
  values[17] = -t4;
  values[18] = -3.0;
}

static void N_H_Q2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6;

  t1 = 3.0/2.0*xi;
  t2 = 3.0/2.0*eta;
  t3 = eta*eta;
  t4 = 15.0/8.0*t3;
  t5 = xi*xi;
  t6 = 15.0/8.0*t5;

  values[0] = 0.0;
  values[1] = t1;
  values[2] = -t2;
  values[3] = -t1;
  values[4] = t2;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = -1.0/4.0+t1+t4-t6;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = -1.0/4.0+t2-t4+t6;
  values[14] = 1.0/4.0+t1-t4+t6;
  values[15] = -1.0/4.0-t2-t4+t6;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
}

static void N_H_Q2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6;

  t1 = 3.0/2.0*xi;
  t2 = 3.0/2.0*zeta;
  t3 = xi*xi;
  t4 = 15.0/8.0*t3;
  t5 = zeta*zeta;
  t6 = 15.0/8.0*t5;

  values[0] = t1;
  values[1] = 0.0;
  values[2] = -t2;
  values[3] = 0.0;
  values[4] = t2;
  values[5] = -t1;
  values[6] = -1.0/4.0+t2+t4-t6;
  values[7] = 0.0;
  values[8] = 1.0/4.0+t1+t4-t6;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = -1.0/4.0+t1-t4+t6;
  values[17] = 1.0/4.0+t2-t4+t6;
  values[18] = 0.0;
}

static void N_H_Q2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t8;

  t1 = 3.0/2.0*zeta;
  t2 = 3.0/2.0*xi;
  t4 = 15.0/4.0*eta*zeta;
  t8 = 15.0/4.0*xi*eta;

  values[0] = t1;
  values[1] = 3.0/2.0;
  values[2] = -t2;
  values[3] = 3.0/2.0;
  values[4] = t2;
  values[5] = -t1;
  values[6] = 0.0;
  values[7] = t1-t4;
  values[8] = 0.0;
  values[9] = t1+t4;
  values[10] = t8;
  values[11] = -t4;
  values[12] = t4;
  values[13] = t2-t8;
  values[14] = -t8;
  values[15] = -t2-t8;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = -3.0;
}

static void N_H_Q2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6;

  t1 = 3.0/2.0*eta;
  t2 = 3.0/2.0*zeta;
  t3 = eta*eta;
  t4 = 15.0/8.0*t3;
  t5 = zeta*zeta;
  t6 = 15.0/8.0*t5;

  values[0] = t1;
  values[1] = t2;
  values[2] = 0.0;
  values[3] = -t2;
  values[4] = 0.0;
  values[5] = -t1;
  values[6] = 0.0;
  values[7] = -1.0/4.0+t1-t4+t6;
  values[8] = 0.0;
  values[9] = 1.0/4.0+t1+t4-t6;
  values[10] = 0.0;
  values[11] = 1.0/4.0+t2-t4+t6;
  values[12] = -1.0/4.0+t2+t4-t6;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
}

static void N_H_Q2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t7;

  t1 = 3.0/2.0*eta;
  t2 = 3.0/2.0*xi;
  t4 = 15.0/4.0*xi*zeta;
  t7 = 15.0/4.0*eta*zeta;

  values[0] = 3.0/2.0;
  values[1] = t1;
  values[2] = -t2;
  values[3] = -t1;
  values[4] = t2;
  values[5] = 3.0/2.0;
  values[6] = t2-t4;
  values[7] = t7;
  values[8] = -t4;
  values[9] = -t7;
  values[10] = 0.0;
  values[11] = t1+t7;
  values[12] = t1-t7;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = t4;
  values[17] = t2+t4;
  values[18] = -3.0;
}

static int N_H_Q2_3D_ChangeU0[1] = {  6 };
static int N_H_Q2_3D_ChangeU1[1] = {  7 };
static int N_H_Q2_3D_ChangeU2[1] = {  8 };
static int N_H_Q2_3D_ChangeU3[1] = {  9 };
static int N_H_Q2_3D_ChangeU4[1] = { 10 };
static int N_H_Q2_3D_ChangeU5[1] = { 11 };

static int N_H_Q2_3D_ChangeV0[1] = { 12 };
static int N_H_Q2_3D_ChangeV1[1] = { 13 };
static int N_H_Q2_3D_ChangeV2[1] = { 14 };
static int N_H_Q2_3D_ChangeV3[1] = { 15 };
static int N_H_Q2_3D_ChangeV4[1] = { 16 };
static int N_H_Q2_3D_ChangeV5[1] = { 17 };

static int *N_H_Q2_3D_ChangeU[6] = { N_H_Q2_3D_ChangeU0, N_H_Q2_3D_ChangeU1,
                                     N_H_Q2_3D_ChangeU2, N_H_Q2_3D_ChangeU3,
                                     N_H_Q2_3D_ChangeU4, N_H_Q2_3D_ChangeU5 };

static int *N_H_Q2_3D_ChangeV[6] = { N_H_Q2_3D_ChangeV0, N_H_Q2_3D_ChangeV1,
                                     N_H_Q2_3D_ChangeV2, N_H_Q2_3D_ChangeV3,
                                     N_H_Q2_3D_ChangeV4, N_H_Q2_3D_ChangeV5 };

static int **N_H_Q2_3D_Change[2] = { N_H_Q2_3D_ChangeU, N_H_Q2_3D_ChangeV };

TBaseFunct3D *BF_N_H_Q2_3D_Obj = 
new TBaseFunct3D(19, BF_N_H_Q2_3D, BFUnitHexahedron, 
                 N_H_Q2_3D_Funct, N_H_Q2_3D_DeriveXi,
                 N_H_Q2_3D_DeriveEta, N_H_Q2_3D_DeriveZeta,
                 N_H_Q2_3D_DeriveXiXi, N_H_Q2_3D_DeriveXiEta,
                 N_H_Q2_3D_DeriveXiZeta, N_H_Q2_3D_DeriveEtaEta,
                 N_H_Q2_3D_DeriveEtaZeta, N_H_Q2_3D_DeriveZetaZeta,
                 3, 2,
                 1, N_H_Q2_3D_Change);
