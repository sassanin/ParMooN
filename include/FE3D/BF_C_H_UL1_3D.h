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
// Q1 element with bubble, conforming, 3D
// ***********************************************************************

static void C_H_UL1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t4, t5, t7, t8, t9, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33;

      t1 = zeta/8.0;
      t2 = eta/8.0;
      t4 = eta*zeta/8.0;
      t5 = xi/8.0;
      t7 = xi*zeta/8.0;
      t8 = xi*eta;
      t9 = t8/8.0;
      t11 = t8*zeta/8.0;
      t12 = zeta*zeta;
      t13 = 27.0/64.0*t12;
      t14 = eta*eta;
      t15 = 27.0/64.0*t14;
      t16 = t14*t12;
      t17 = 27.0/64.0*t16;
      t18 = xi*xi;
      t19 = 27.0/64.0*t18;
      t20 = t18*t12;
      t21 = 27.0/64.0*t20;
      t22 = t18*t14;
      t23 = 27.0/64.0*t22;
      t24 = t22*t12;
      t25 = 27.0/64.0*t24;
      t26 = -19.0/64.0-t1-t2+t4-t5+t7+t9-t11+t13+t15-t17+t19-t21-t23+t25;
      t27 = -19.0/64.0-t1-t2+t4+t5-t7-t9+t11+t13+t15-t17+t19-t21-t23+t25;
      t28 = -19.0/64.0-t1+t2-t4-t5+t7-t9+t11+t13+t15-t17+t19-t21-t23+t25;
      t29 = -19.0/64.0-t1+t2-t4+t5-t7+t9-t11+t13+t15-t17+t19-t21-t23+t25;
      t30 = -19.0/64.0+t1-t2-t4-t5-t7+t9+t11+t13+t15-t17+t19-t21-t23+t25;
      t31 = -19.0/64.0+t1-t2-t4+t5+t7-t9-t11+t13+t15-t17+t19-t21-t23+t25;
      t32 = -19.0/64.0+t1+t2+t4-t5-t7-t9-t11+t13+t15-t17+t19-t21-t23+t25;
      t33 = -19.0/64.0+t1+t2+t4+t5+t7+t9+t11+t13+t15-t17+t19-t21-t23+t25;
      values[0] = t26;
      values[1] = t27;
      values[2] = t28;
      values[3] = t29;
      values[4] = t30;
      values[5] = t31;
      values[6] = t32;
      values[7] = t33;
      values[8] = 27.0/64.0-27.0/64.0*t12-27.0/64.0*t14+27.0/64.0*t16-27.0/64.0*t18+
27.0/64.0*t20+27.0/64.0*t22-27.0/64.0*t24;

}

static void C_H_UL1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;

      t1 = zeta/8.0;
      t2 = eta/8.0;
      t4 = eta*zeta/8.0;
      t5 = 27.0/32.0*xi;
      t6 = zeta*zeta;
      t7 = xi*t6;
      t8 = 27.0/32.0*t7;
      t9 = eta*eta;
      t10 = xi*t9;
      t11 = 27.0/32.0*t10;
      t12 = t10*t6;
      t13 = 27.0/32.0*t12;
      values[0] = -1.0/8.0+t1+t2-t4+t5-t8-t11+t13;
      values[1] = 1.0/8.0-t1-t2+t4+t5-t8-t11+t13;
      values[2] = -1.0/8.0+t1-t2+t4+t5-t8-t11+t13;
      values[3] = 1.0/8.0-t1+t2-t4+t5-t8-t11+t13;
      values[4] = -1.0/8.0-t1+t2+t4+t5-t8-t11+t13;
      values[5] = 1.0/8.0+t1-t2-t4+t5-t8-t11+t13;
      values[6] = -1.0/8.0-t1-t2-t4+t5-t8-t11+t13;
      values[7] = 1.0/8.0+t1+t2+t4+t5-t8-t11+t13;
      values[8] = -27.0/32.0*xi+27.0/32.0*t7+27.0/32.0*t10-27.0/32.0*t12;
}

static void C_H_UL1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;

      t1 = zeta/8.0;
      t2 = xi/8.0;
      t4 = xi*zeta/8.0;
      t5 = 27.0/32.0*eta;
      t6 = zeta*zeta;
      t7 = eta*t6;
      t8 = 27.0/32.0*t7;
      t9 = xi*xi;
      t10 = t9*eta;
      t11 = 27.0/32.0*t10;
      t12 = t10*t6;
      t13 = 27.0/32.0*t12;
      values[0] = -1.0/8.0+t1+t2-t4+t5-t8-t11+t13;
      values[1] = -1.0/8.0+t1-t2+t4+t5-t8-t11+t13;
      values[2] = 1.0/8.0-t1-t2+t4+t5-t8-t11+t13;
      values[3] = 1.0/8.0-t1+t2-t4+t5-t8-t11+t13;
      values[4] = -1.0/8.0-t1+t2+t4+t5-t8-t11+t13;
      values[5] = -1.0/8.0-t1-t2-t4+t5-t8-t11+t13;
      values[6] = 1.0/8.0+t1-t2-t4+t5-t8-t11+t13;
      values[7] = 1.0/8.0+t1+t2+t4+t5-t8-t11+t13;
      values[8] = -27.0/32.0*eta+27.0/32.0*t7+27.0/32.0*t10-27.0/32.0*t12;
}

static void C_H_UL1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t5, t6, t7, t8, t9, t10, t11, t13, t14;

      t1 = eta/8.0;
      t2 = xi/8.0;
      t4 = xi*eta/8.0;
      t5 = 27.0/32.0*zeta;
      t6 = eta*eta;
      t7 = t6*zeta;
      t8 = 27.0/32.0*t7;
      t9 = xi*xi;
      t10 = t9*zeta;
      t11 = 27.0/32.0*t10;
      t13 = t9*t6*zeta;
      t14 = 27.0/32.0*t13;
      values[0] = -1.0/8.0+t1+t2-t4+t5-t8-t11+t14;
      values[1] = -1.0/8.0+t1-t2+t4+t5-t8-t11+t14;
      values[2] = -1.0/8.0-t1+t2+t4+t5-t8-t11+t14;
      values[3] = -1.0/8.0-t1-t2-t4+t5-t8-t11+t14;
      values[4] = 1.0/8.0-t1-t2+t4+t5-t8-t11+t14;
      values[5] = 1.0/8.0-t1+t2-t4+t5-t8-t11+t14;
      values[6] = 1.0/8.0+t1-t2-t4+t5-t8-t11+t14;
      values[7] = 1.0/8.0+t1+t2+t4+t5-t8-t11+t14;
      values[8] = -27.0/32.0*zeta+27.0/32.0*t7+27.0/32.0*t10-27.0/32.0*t13;
}

static void C_H_UL1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4;

      t1 = zeta*zeta;
      t2 = eta*eta;
      t4 = 1.0-t1-t2+t2*t1;
      values[0] = 27.0/32.0*t4;
      values[1] = 27.0/32.0*t4;
      values[2] = 27.0/32.0*t4;
      values[3] = 27.0/32.0*t4;
      values[4] = 27.0/32.0*t4;
      values[5] = 27.0/32.0*t4;
      values[6] = 27.0/32.0*t4;
      values[7] = 27.0/32.0*t4;
      values[8] = -27.0/32.0*t4;
}

static void C_H_UL1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

      t1 = zeta/8.0;
      t2 = xi*eta;
      t3 = 27.0/16.0*t2;
      t4 = zeta*zeta;
      t5 = t2*t4;
      t6 = 27.0/16.0*t5;
      t7 = 1.0/8.0-t1-t3+t6;
      t8 = -1.0/8.0+t1-t3+t6;
      t9 = 1.0/8.0+t1-t3+t6;
      t10 = -1.0/8.0-t1-t3+t6;
      values[0] = t7;
      values[1] = t8;
      values[2] = t8;
      values[3] = t7;
      values[4] = t9;
      values[5] = t10;
      values[6] = t10;
      values[7] = t9;
      values[8] = 27.0/16.0*t2-27.0/16.0*t5;
}

static void C_H_UL1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t6, t7, t8, t9, t10, t11;

      t1 = eta/8.0;
      t2 = xi*zeta;
      t3 = 27.0/16.0*t2;
      t4 = eta*eta;
      t6 = xi*t4*zeta;
      t7 = 27.0/16.0*t6;
      t8 = 1.0/8.0-t1-t3+t7;
      t9 = -1.0/8.0+t1-t3+t7;
      t10 = 1.0/8.0+t1-t3+t7;
      t11 = -1.0/8.0-t1-t3+t7;
      values[0] = t8;
      values[1] = t9;
      values[2] = t10;
      values[3] = t11;
      values[4] = t9;
      values[5] = t8;
      values[6] = t11;
      values[7] = t10;
      values[8] = 27.0/16.0*t2-27.0/16.0*t6;
}

static void C_H_UL1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4;

      t1 = zeta*zeta;
      t2 = xi*xi;
      t4 = 1.0-t1-t2+t2*t1;
      values[0] = 27.0/32.0*t4;
      values[1] = 27.0/32.0*t4;
      values[2] = 27.0/32.0*t4;
      values[3] = 27.0/32.0*t4;
      values[4] = 27.0/32.0*t4;
      values[5] = 27.0/32.0*t4;
      values[6] = 27.0/32.0*t4;
      values[7] = 27.0/32.0*t4;
      values[8] = -27.0/32.0*t4;
}

static void C_H_UL1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t6, t7, t8, t9, t10, t11;

      t1 = xi/8.0;
      t2 = eta*zeta;
      t3 = 27.0/16.0*t2;
      t4 = xi*xi;
      t6 = t4*eta*zeta;
      t7 = 27.0/16.0*t6;
      t8 = 1.0/8.0-t1-t3+t7;
      t9 = 1.0/8.0+t1-t3+t7;
      t10 = -1.0/8.0+t1-t3+t7;
      t11 = -1.0/8.0-t1-t3+t7;
      values[0] = t8;
      values[1] = t9;
      values[2] = t10;
      values[3] = t11;
      values[4] = t10;
      values[5] = t11;
      values[6] = t8;
      values[7] = t9;
      values[8] = 27.0/16.0*t2-27.0/16.0*t6;
}

static void C_H_UL1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4;

      t1 = eta*eta;
      t2 = xi*xi;
      t4 = 1.0-t1-t2+t2*t1;
      values[0] = 27.0/32.0*t4;
      values[1] = 27.0/32.0*t4;
      values[2] = 27.0/32.0*t4;
      values[3] = 27.0/32.0*t4;
      values[4] = 27.0/32.0*t4;
      values[5] = 27.0/32.0*t4;
      values[6] = 27.0/32.0*t4;
      values[7] = 27.0/32.0*t4;
      values[8] = -27.0/32.0*t4;
}

TBaseFunct3D *BF_C_H_UL1_3D_Obj = 
new TBaseFunct3D(9, BF_C_H_UL1_3D, BFUnitHexahedron, 
                 C_H_UL1_3D_Funct, C_H_UL1_3D_DeriveXi,
                 C_H_UL1_3D_DeriveEta, C_H_UL1_3D_DeriveZeta,
                 C_H_UL1_3D_DeriveXiXi, C_H_UL1_3D_DeriveXiEta,
                 C_H_UL1_3D_DeriveXiZeta, C_H_UL1_3D_DeriveEtaEta,
                 C_H_UL1_3D_DeriveEtaZeta, C_H_UL1_3D_DeriveZetaZeta,
                 2, 1,
                 0, NULL);
