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
// P3 element, conforming, 3D
// ***********************************************************************

static void C_T_P3_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
  double t13, t14, t15, t16, t17;

  t1 = zeta*zeta;
  t2 = t1*zeta;
  t3 = eta*zeta;
  t4 = eta*t1;
  t5 = eta*eta;
  t6 = t5*zeta;
  t7 = t5*eta;
  t8 = xi*zeta;
  t9 = xi*t1;
  t10 = xi*eta;
  t11 = t10*zeta;
  t12 = xi*t5;
  t13 = xi*xi;
  t14 = t13*zeta;
  t15 = t13*eta;
  t16 = t13*xi;
  t17 = 1.0-11.0/2.0*zeta+9.0*t1-9.0/2.0*t2-11.0/2.0*eta+18.0*t3
       -27.0/2.0*t4+9.0*t5-27.0/2.0*t6-9.0/2.0*t7-11.0/2.0*xi+18.0*t8
       -27.0/2.0*t9+18.0*t10-27.0*t11-27.0/2.0*t12+9.0*t13
       -27.0/2.0*t14-27.0/2.0*t15-9.0/2.0*t16;

  values[0] = t17;
  values[1] = 9.0*xi-45.0/2.0*t8+27.0/2.0*t9-45.0/2.0*t10+27.0*t11
             +27.0/2.0*t12-45.0/2.0*t13+27.0*t14+27.0*t15+27.0/2.0*t16;
  values[2] = -9.0/2.0*xi+9.0/2.0*t8+9.0/2.0*t10+18.0*t13-27.0/2.0*t14
              -27.0/2.0*t15-27.0/2.0*t16;
  values[3] = xi-9.0/2.0*t13+9.0/2.0*t16;
  values[4] = 9.0*eta-45.0/2.0*t3+27.0/2.0*t4-45.0/2.0*t5+27.0*t6
             +27.0/2.0*t7-45.0/2.0*t10+27.0*t11+27.0*t12+27.0/2.0*t15;
  values[5] = 27.0*t10-27.0*t11-27.0*t12-27.0*t15;
  values[6] = -9.0/2.0*t10+27.0/2.0*t15;
  values[7] = -9.0/2.0*eta+9.0/2.0*t3+18.0*t5-27.0/2.0*t6-27.0/2.0*t7
              +9.0/2.0*t10-27.0/2.0*t12;
  values[8] = -9.0/2.0*t10+27.0/2.0*t12;
  values[9] = eta-9.0/2.0*t5+9.0/2.0*t7;
  values[10] = 9.0*zeta-45.0/2.0*t1+27.0/2.0*t2-45.0/2.0*t3+27.0*t4
              +27.0/2.0*t6-45.0/2.0*t8+27.0*t9+27.0*t11+27.0/2.0*t14;
  values[11] = 27.0*t8-27.0*t9-27.0*t11-27.0*t14;
  values[12] = -9.0/2.0*t8+27.0/2.0*t14;
  values[13] = 27.0*t3-27.0*t4-27.0*t6-27.0*t11;
  values[14] = 27.0*t11;
  values[15] = -9.0/2.0*t3+27.0/2.0*t6;
  values[16] = -9.0/2.0*zeta+18.0*t1-27.0/2.0*t2+9.0/2.0*t3-27.0/2.0*t4
               +9.0/2.0*t8-27.0/2.0*t9;
  values[17] = -9.0/2.0*t8+27.0/2.0*t9;
  values[18] = -9.0/2.0*t3+27.0/2.0*t4;
  values[19] = zeta-9.0/2.0*t1+9.0/2.0*t2;
}

static void C_T_P3_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6;

  t1 = zeta*zeta;
  t2 = eta*zeta;
  t3 = eta*eta;
  t4 = xi*zeta;
  t5 = xi*eta;
  t6 = xi*xi;

  values[0] = -11.0/2.0+18.0*zeta-27.0/2.0*t1+18.0*eta-27.0*t2
              -27.0/2.0*t3+18.0*xi-27.0*t4-27.0*t5-27.0/2.0*t6;
  values[1] = 9.0-45.0/2.0*zeta+27.0/2.0*t1-45.0/2.0*eta+27.0*t2
             +27.0/2.0*t3-45.0*xi+54.0*t4+54.0*t5+81.0/2.0*t6;
  values[2] = -9.0/2.0+9.0/2.0*zeta+9.0/2.0*eta+36.0*xi-27.0*t4
              -27.0*t5-81.0/2.0*t6;
  values[3] = 1.0-9.0*xi+27.0/2.0*t6;
  values[4] = -45.0/2.0*eta+27.0*t2+27.0*t3+27.0*t5;
  values[5] = 27.0*eta-27.0*t2-27.0*t3-54.0*t5;
  values[6] = -9.0/2.0*eta+27.0*t5;
  values[7] = 9.0/2.0*eta-27.0/2.0*t3;
  values[8] = -9.0/2.0*eta+27.0/2.0*t3;
  values[9] = 0.0;
  values[10] = -45.0/2.0*zeta+27.0*t1+27.0*t2+27.0*t4;
  values[11] = 27.0*zeta-27.0*t1-27.0*t2-54.0*t4;
  values[12] = -9.0/2.0*zeta+27.0*t4;
  values[13] = -27.0*t2;
  values[14] = 27.0*t2;
  values[15] = 0.0;
  values[16] = 9.0/2.0*zeta-27.0/2.0*t1;
  values[17] = -9.0/2.0*zeta+27.0/2.0*t1;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void C_T_P3_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6;

  t1 = zeta*zeta;
  t2 = eta*zeta;
  t3 = eta*eta;
  t4 = xi*zeta;
  t5 = xi*eta;
  t6 = xi*xi;

  values[0] = -11.0/2.0+18.0*zeta-27.0/2.0*t1+18.0*eta-27.0*t2
              -27.0/2.0*t3+18.0*xi-27.0*t4-27.0*t5-27.0/2.0*t6;
  values[1] = -45.0/2.0*xi+27.0*t4+27.0*t5+27.0*t6;
  values[2] = 9.0/2.0*xi-27.0/2.0*t6;
  values[3] = 0.0;
  values[4] = 9.0-45.0/2.0*zeta+27.0/2.0*t1-45.0*eta+54.0*t2+81.0/2.0*t3
             -45.0/2.0*xi+27.0*t4+54.0*t5+27.0/2.0*t6;
  values[5] = 27.0*xi-27.0*t4-54.0*t5-27.0*t6;
  values[6] = -9.0/2.0*xi+27.0/2.0*t6;
  values[7] = -9.0/2.0+9.0/2.0*zeta+36.0*eta-27.0*t2-81.0/2.0*t3
              +9.0/2.0*xi-27.0*t5;
  values[8] = -9.0/2.0*xi+27.0*t5;
  values[9] = 1.0-9.0*eta+27.0/2.0*t3;
  values[10] = -45.0/2.0*zeta+27.0*t1+27.0*t2+27.0*t4;
  values[11] = -27.0*t4;
  values[12] = 0.0;
  values[13] = 27.0*zeta-27.0*t1-54.0*t2-27.0*t4;
  values[14] = 27.0*t4;
  values[15] = -9.0/2.0*zeta+27.0*t2;
  values[16] = 9.0/2.0*zeta-27.0/2.0*t1;
  values[17] = 0.0;
  values[18] = -9.0/2.0*zeta+27.0/2.0*t1;
  values[19] = 0.0;
}

static void C_T_P3_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6;

  t1 = zeta*zeta;
  t2 = eta*zeta;
  t3 = eta*eta;
  t4 = xi*zeta;
  t5 = xi*eta;
  t6 = xi*xi;

  values[0] = -11.0/2.0+18.0*zeta-27.0/2.0*t1+18.0*eta-27.0*t2
              -27.0/2.0*t3+18.0*xi-27.0*t4-27.0*t5-27.0/2.0*t6;
  values[1] = -45.0/2.0*xi+27.0*t4+27.0*t5+27.0*t6;
  values[2] = 9.0/2.0*xi-27.0/2.0*t6;
  values[3] = 0.0;
  values[4] = -45.0/2.0*eta+27.0*t2+27.0*t3+27.0*t5;
  values[5] = -27.0*t5;
  values[6] = 0.0;
  values[7] = 9.0/2.0*eta-27.0/2.0*t3;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 9.0-45.0*zeta+81.0/2.0*t1-45.0/2.0*eta+54.0*t2
              +27.0/2.0*t3-45.0/2.0*xi+54.0*t4+27.0*t5+27.0/2.0*t6;
  values[11] = 27.0*xi-54.0*t4-27.0*t5-27.0*t6;
  values[12] = -9.0/2.0*xi+27.0/2.0*t6;
  values[13] = 27.0*eta-54.0*t2-27.0*t3-27.0*t5;
  values[14] = 27.0*t5;
  values[15] = -9.0/2.0*eta+27.0/2.0*t3;
  values[16] = -9.0/2.0+36.0*zeta-81.0/2.0*t1+9.0/2.0*eta-27.0*t2
               +9.0/2.0*xi-27.0*t4;
  values[17] = -9.0/2.0*xi+27.0*t4;
  values[18] = -9.0/2.0*eta+27.0*t2;
  values[19] = 1.0-9.0*zeta+27.0/2.0*t1;
}

static void C_T_P3_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 18.0-27.0*zeta-27.0*eta-27.0*xi;
  values[1] = -45.0+54.0*zeta+54.0*eta+81.0*xi;
  values[2] = 36.0-27.0*zeta-27.0*eta-81.0*xi;
  values[3] = -9.0+27.0*xi;
  values[4] = 27.0*eta;
  values[5] = -54.0*eta;
  values[6] = 27.0*eta;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 27.0*zeta;
  values[11] = -54.0*zeta;
  values[12] = 27.0*zeta;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void C_T_P3_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 18.0-27.0*zeta-27.0*eta-27.0*xi;
  values[1] = -45.0/2.0+27.0*zeta+27.0*eta+54.0*xi;
  values[2] = 9.0/2.0-27.0*xi;
  values[3] = 0.0;
  values[4] = -45.0/2.0+27.0*zeta+54.0*eta+27.0*xi;
  values[5] = 27.0-27.0*zeta-54.0*eta-54.0*xi;
  values[6] = -9.0/2.0+27.0*xi;
  values[7] = 9.0/2.0-27.0*eta;
  values[8] = -9.0/2.0+27.0*eta;
  values[9] = 0.0;
  values[10] = 27.0*zeta;
  values[11] = -27.0*zeta;
  values[12] = 0.0;
  values[13] = -27.0*zeta;
  values[14] = 27.0*zeta;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void C_T_P3_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 18.0-27.0*zeta-27.0*eta-27.0*xi;
  values[1] = -45.0/2.0+27.0*zeta+27.0*eta+54.0*xi;
  values[2] = 9.0/2.0-27.0*xi;
  values[3] = 0.0;
  values[4] = 27.0*eta;
  values[5] = -27.0*eta;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = -45.0/2.0+54.0*zeta+27.0*eta+27.0*xi;
  values[11] = 27.0-54.0*zeta-27.0*eta-54.0*xi;
  values[12] = -9.0/2.0+27.0*xi;
  values[13] = -27.0*eta;
  values[14] = 27.0*eta;
  values[15] = 0.0;
  values[16] = 9.0/2.0-27.0*zeta;
  values[17] = -9.0/2.0+27.0*zeta;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void C_T_P3_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 18.0-27.0*zeta-27.0*eta-27.0*xi;
  values[1] = 27.0*xi;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = -45.0+54.0*zeta+81.0*eta+54.0*xi;
  values[5] = -54.0*xi;
  values[6] = 0.0;
  values[7] = 36.0-27.0*zeta-81.0*eta-27.0*xi;
  values[8] = 27.0*xi;
  values[9] = -9.0+27.0*eta;
  values[10] = 27.0*zeta;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = -54.0*zeta;
  values[14] = 0.0;
  values[15] = 27.0*zeta;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void C_T_P3_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 18.0-27.0*zeta-27.0*eta-27.0*xi;
  values[1] = 27.0*xi;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = -45.0/2.0+27.0*zeta+54.0*eta+27.0*xi;
  values[5] = -27.0*xi;
  values[6] = 0.0;
  values[7] = 9.0/2.0-27.0*eta;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = -45.0/2.0+54.0*zeta+27.0*eta+27.0*xi;
  values[11] = -27.0*xi;
  values[12] = 0.0;
  values[13] = 27.0-54.0*zeta-54.0*eta-27.0*xi;
  values[14] = 27.0*xi;
  values[15] = -9.0/2.0+27.0*eta;
  values[16] = 9.0/2.0-27.0*zeta;
  values[17] = 0.0;
  values[18] = -9.0/2.0+27.0*zeta;
  values[19] = 0.0;
}

static void C_T_P3_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 18.0-27.0*zeta-27.0*eta-27.0*xi;
  values[1] = 27.0*xi;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 27.0*eta;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = -45.0+81.0*zeta+54.0*eta+54.0*xi;
  values[11] = -54.0*xi;
  values[12] = 0.0;
  values[13] = -54.0*eta;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = 36.0-81.0*zeta-27.0*eta-27.0*xi;
  values[17] = 27.0*xi;
  values[18] = 27.0*eta;
  values[19] = -9.0+27.0*zeta;
}

TBaseFunct3D *BF_C_T_P3_3D_Obj = 
new TBaseFunct3D(20, BF_C_T_P3_3D, BFUnitTetrahedron, 
                 C_T_P3_3D_Funct, C_T_P3_3D_DeriveXi,
                 C_T_P3_3D_DeriveEta, C_T_P3_3D_DeriveZeta,
                 C_T_P3_3D_DeriveXiXi, C_T_P3_3D_DeriveXiEta,
                 C_T_P3_3D_DeriveXiZeta, C_T_P3_3D_DeriveEtaEta,
                 C_T_P3_3D_DeriveEtaZeta, C_T_P3_3D_DeriveZetaZeta,
                 3, 3,
                 0, NULL);
