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
// P2 element with 4 face and 1 cell bubbles, conforming, 3D
// 
// Author:     Sashikumaar Ganesan
//
// ***********************************************************************

static void C_T_B2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t3, t4, t5,t6,t7,t8,t9,t10,t11,t12,t13;
  double t14,t15,t18,t19,t20,t22,t23,t24,t25,t28,t31,t33,t41,t44;

      t1 = 1.0-xi-eta-zeta;
      t2 = t1*t1;
      t3 = zeta*t1;
      t4 = t1*xi;
      t5 = eta*t1;
      t6 = eta*zeta;
      t7 = t6*t1;
      t8 = 3.0*t7;
      t9 = t3*xi;
      t10 = 3.0*t9;
      t11 = t4*eta;
      t12 = 3.0*t11;
      t13 = xi*eta;
      t14 = t13*t3;
      t15 = 4.0*t14;
      t18 = 12.0*t9;
      t19 = 12.0*t11;
      t20 = 32.0*t14;
      t22 = xi*xi;
      t23 = xi*zeta;
      t24 = t13*zeta;
      t25 = 3.0*t24;
      t28 = 12.0*t7;
      t31 = 12.0*t24;
      t33 = eta*eta;
      t41 = zeta*zeta;
      t44 = 108.0*t14;

      values[0] = t2-t3-t4-t5+t8+t10+t12-t15;
      values[1] = 4.0*t4-t18-t19+t20;
      values[2] = t22-t13-t4-t23+t25+t10+t12-t15;
      values[3] = 4.0*t5-t28-t19+t20;
      values[4] = 4.0*t13-t31-t19+t20;
      values[5] = t33-t13-t6-t5+t25+t8+t12-t15;
      values[6] = 4.0*t3-t28-t18+t20;
      values[7] = 4.0*t23-t31-t18+t20;
      values[8] = 4.0*t6-t31-t28+t20;
      values[9] = t41-t6-t3-t23+t25+t8+t10-t15;
      values[10] = 27.0*t11-t44;
      values[11] = 27.0*t9-t44;
      values[12] = 27.0*t24-t44;
      values[13] = 27.0*t7-t44;
      values[14] = 256.0*t14;
}

static void C_T_B2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
      double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12,t13, t14;
      double t15, t16, t17, t18, t19, t21, t22,t23,t24,t25,t26,t29,t35,t39,t40,t45;

      t1 = 4.0*xi;
      t2 = 4.0*eta;
      t3 = 4.0*zeta;
      t4 = eta*zeta;
      t5 = 3.0*t4;
      t6 = xi*zeta;
      t7 = 3.0*t6;
      t8 = 1.0-xi-eta-zeta;
      t9 = zeta*t8;
      t10 = 3.0*t9;
      t11 = xi*eta;
      t12 = 3.0*t11;
      t13 = eta*t8;
      t14 = 3.0*t13;
      t15 = t4*t8;
      t16 = 4.0*t15;
      t17 = t11*zeta;
      t18 = 4.0*t17;
      t19 = -3.0+t1+t2+t3-t5-t7+t10-t12+t14-t16+t18;
      t21 = 12.0*t6;
      t22 = 12.0*t9;
      t23 = 12.0*t11;
      t24 = 12.0*t13;
      t25 = 32.0*t15;
      t26 = 32.0*t17;
      t29 = 12.0*t4;
      t35 = t15-t17;
      t39 = 108.0*t15;
      t40 = 108.0*t17;
      t45 = 27.0*t4;

      values[0] = t19;
      values[1] = -8.0*xi+4.0-t2-t3+t21-t22+t23-t24+t25-t26;
      values[2] = t1-1.0+t5-t7+t10-t12+t14-t16+t18;
      values[3] = -t2+t29+t23-t24+t25-t26;
      values[4] = t2-t29+t23-t24+t25-t26;
      values[5] = -t12+t14-t16+t18;
      values[6] = -t3+t29+t21-t22+t25-t26;
      values[7] = t3-t29+t21-t22+t25-t26;
      values[8] = 32.0*t35;
      values[9] = -t7+t10-t16+t18;
      values[10] = -27.0*t11+27.0*t13-t39+t40;
      values[11] = -27.0*t6+27.0*t9-t39+t40;
      values[12] = t45-t39+t40;
      values[13] = -t45-t39+t40;
      values[14] = 256.0*t35;
}

static void C_T_B2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{

 double t1 = 4.0*xi;
 double t2 = 4.0*eta;
 double t3 = 4.0*zeta;
 double t4 = 1.0-xi-eta-zeta;
 double t5 = zeta*t4;
 double t6 = 3.0*t5;
 double t7 = eta*zeta;
 double t8 = 3.0*t7;
 double t9 = xi*zeta;
 double t10 = 3.0*t9;
 double t11 = xi*eta;
 double t12 = 3.0*t11;
 double t13 = t4*xi;
 double t14 = 3.0*t13;
 double t15 = t5*xi;
 double t16 = 4.0*t15;
 double t17 = t11*zeta;
 double t18 = 4.0*t17;
 double t19 = -3.0+t1+t2+t3+t6-t8-t10-t12+t14-t16+t18;
 double t20 = 12.0*t9;
 double t21 = 12.0*t11;
 double t22 = 12.0*t13;
 double t23 = 32.0*t15;
 double t24 = 32.0*t17;
 double t28 = 12.0*t5;
 double t29 = 12.0*t7;
 double t34 = t15-t17;
 double t39 = 108.0*t15;
 double t40 = 108.0*t17;
 double t42 = 27.0*t9;

      values[0] = t19;
      values[1] = -t1+t20+t21-t22+t23-t24;
      values[2] = -t12+t14-t16+t18;
      values[3] = 4.0-t1-8.0*eta-t3-t28+t29+t21-t22+t23-t24;
      values[4] = t1-t20+t21-t22+t23-t24;
      values[5] = t2-1.0+t10+t6-t8-t12+t14-t16+t18;
      values[6] = -t3-t28+t29+t20+t23-t24;
      values[7] = 32.0*t34;
      values[8] = t3-t20-t28+t29+t23-t24;
      values[9] = t6-t8-t16+t18;
      values[10] = -27.0*t11+27.0*t13-t39+t40;
      values[11] = -t42-t39+t40;
      values[12] = t42-t39+t40;
      values[13] = 27.0*t5-27.0*t7-t39+t40;
      values[14] = 256.0*t34;
}

static void C_T_B2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{

 double t1 = 4.0*xi;
 double t2 = 4.0*eta;
 double t3 = 4.0*zeta;
 double t4 = 1.0-xi-eta-zeta;
 double t5 = eta*t4;
 double t6 = 3.0*t5;
 double t7 = eta*zeta;
 double t8 = 3.0*t7;
 double t9 = t4*xi;
 double t10 = 3.0*t9;
 double t11 = xi*zeta;
 double t12 = 3.0*t11;
 double t13 = xi*eta;
 double t14 = 3.0*t13;
 double t15 = t9*eta;
 double t16 = 4.0*t15;
 double t17 = t13*zeta;
 double t18 = 4.0*t17;
 double t19 = -3.0+t1+t2+t3+t6-t8+t10-t12-t14-t16+t18;
 double t20 = 12.0*t9;
 double t21 = 12.0*t11;
 double t22 = 12.0*t13;
 double t23 = 32.0*t15;
 double t24 = 32.0*t17;
 double t27 = 12.0*t5;
 double t28 = 12.0*t7;
 double t30 = t15-t17;
 double t37 = 27.0*t13;
 double t38 = 108.0*t15;
 double t39 = 108.0*t17;

      values[0] = t19;
      values[1] = -t1-t20+t21+t22+t23-t24;
      values[2] = t10-t12-t16+t18;
      values[3] = -t2-t27+t28+t22+t23-t24;
      values[4] = 32.0*t30;
      values[5] = t6-t8-t16+t18;
      values[6] = 4.0-t1-t2-8.0*zeta-t27+t28-t20+t21+t23-t24;
      values[7] = t1-t20+t21-t22+t23-t24;
      values[8] = t2-t22-t27+t28+t23-t24;
      values[9] = t3-1.0+t14+t6-t8+t10-t12-t16+t18;
      values[10] = -t37-t38+t39;
      values[11] = 27.0*t9-27.0*t11-t38+t39;
      values[12] = t37-t38+t39;
      values[13] = 27.0*t5-27.0*t7-t38+t39;
      values[14] = 256.0*t30;
}

static void C_T_B2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{

   double t1 = 6.0*zeta;
   double t2 = 6.0*eta;
   double t3 = eta*zeta;
   double t4 = 8.0*t3;
   double t5 = 4.0-t1-t2+t4;
   double t6 = 24.0*zeta;
   double t7 = 24.0*eta;
   double t8 = 64.0*t3;
   double t10 = t7-t8;
   double t12 = t6-t8;
   double t15 = 216.0*t3;

      values[0] = t5;
      values[1] = -8.0+t6+t7-t8;
      values[2] = t5;
      values[3] = t10;
      values[4] = t10;
      values[5] = -t2+t4;
      values[6] = t12;
      values[7] = t12;
      values[8] = -t8;
      values[9] = -t1+t4;
      values[10] = -54.0*eta+t15;
      values[11] = t15-54.0*zeta;
      values[12] = t15;
      values[13] = t15;
      values[14] = -512.0*t3;
      
    }

static void C_T_B2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
 double t2 = 6.0*xi;
 double t3 = 6.0*eta;
 double t5 = zeta*(1.0-xi-eta-zeta);
 double t6 = 4.0*t5;
 double t7 = eta*zeta;
 double t8 = 4.0*t7;
 double t9 = xi*zeta;
 double t10 = 4.0*t9;
 double t12 = 24.0*zeta;
 double t13 = 24.0*xi;
 double t14 = 24.0*eta;
 double t15 = 32.0*t5;
 double t16 = 32.0*t7;
 double t17 = 32.0*t9;
 double t18 = -16.0+t12+t13+t14+t15-t16-t17;
 double t19 = 3.0*zeta;
 double t20 = -t19-t2+3.0-t3-t6+t8+t10;
 double t23 = t5-t7-t9;
 double t27 = 27.0*zeta;
 double t28 = 108.0*t5;
 double t29 = 108.0*t7;
 double t30 = 108.0*t9;
 double t32 = -t27-t28+t29+t30;

      values[0] = 7.0-9.0*zeta-t2-t3-t6+t8+t10;
      values[1] = t18;
      values[2] = t20;
      values[3] = t18;
      values[4] = -8.0+t13+t14+t15-t16-t17;
      values[5] = t20;
      values[6] = t12+t15-t16-t17;
      values[7] = 32.0*t23;
      values[8] = 32.0*t23;
      values[9] = -t19-t6+t8+t10;
      values[10] = -54.0*xi+27.0-54.0*eta-t27-t28+t29+t30;
      values[11] = t32;
      values[12] = t27-t28+t29+t30;
      values[13] = t32;
      values[14] = 256.0*t23;
}

static void C_T_B2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{

 double t2 = 6.0*xi;
 double t3 = 6.0*zeta;
 double t5 = eta*(1.0-xi-eta-zeta);
 double t6 = 4.0*t5;
 double t7 = eta*zeta;
 double t8 = 4.0*t7;
 double t9 = xi*eta;
 double t10 = 4.0*t9;
 double t12 = 24.0*xi;
 double t13 = 24.0*eta;
 double t14 = 24.0*zeta;
 double t15 = 32.0*t5;
 double t16 = 32.0*t7;
 double t17 = 32.0*t9;
 double t18 = -16.0+t12+t13+t14+t15-t16-t17;
 double t19 = 3.0*eta;
 double t20 = -t19-t2+3.0-t3-t6+t8+t10;
 double t22 = t5-t7-t9;
 double t25 = 27.0*eta;
 double t26 = 108.0*t5;
 double t27 = 108.0*t7;
 double t28 = 108.0*t9;
 double t29 = -t25-t26+t27+t28;

      values[0] = 7.0-9.0*eta-t2-t3-t6+t8+t10;
      values[1] = t18;
      values[2] = t20;
      values[3] = t13+t15-t16-t17;
      values[4] = 32.0*t22;
      values[5] = -t19-t6+t8+t10;
      values[6] = t18;
      values[7] = -8.0+t12+t14+t15-t16-t17;
      values[8] = 32.0*t22;
      values[9] = t20;
      values[10] = t29;
      values[11] = -54.0*xi+27.0-t25-54.0*zeta-t26+t27+t28;
      values[12] = t25-t26+t27+t28;
      values[13] = t29;
      values[14] = 256.0*t22;

     }

static void C_T_B2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
 double t1 = 6.0*zeta;
 double t2 = 6.0*xi;
 double t3 = xi*zeta;
 double t4 = 8.0*t3;
 double t5 = 4.0-t1-t2+t4;
 double t6 = 24.0*xi;
 double t7 = 64.0*t3;
 double t8 = t6-t7;
 double t10 = 24.0*zeta;
 double t12 = t10-t7;
 double t15 = 216.0*t3;

      values[0] = t5;
      values[1] = t8;
      values[2] = -t2+t4;
      values[3] = -8.0+t10+t6-t7;
      values[4] = t8;
      values[5] = t5;
      values[6] = t12;
      values[7] = -t7;
      values[8] = t12;
      values[9] = -t1+t4;
      values[10] = -54.0*xi+t15;
      values[11] = t15;
      values[12] = t15;
      values[13] = -54.0*zeta+t15;
      values[14] = -512.0*t3;
      }

static void C_T_B2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{

 double t2 = 6.0*eta;
 double t3 = 6.0*zeta;
 double t5 = (1.0-xi-eta-zeta)*xi;
 double t6 = 4.0*t5;
 double t7 = xi*zeta;
 double t8 = 4.0*t7;
 double t9 = xi*eta;
 double t10 = 4.0*t9;
 double t12 = 24.0*xi;
 double t13 = 32.0*t5;
 double t14 = 32.0*t7;
 double t15 = 32.0*t9;
 double t17 = 3.0*xi;
 double t19 = 24.0*eta;
 double t20 = 24.0*zeta;
 double t21 = -16.0+t12+t19+t20+t13-t14-t15;
 double t22 = t5-t7-t9;
 double t23 = -t17+3.0-t2-t3-t6+t8+t10;
 double t25 = 27.0*xi;
 double t26 = 108.0*t5;
 double t27 = 108.0*t7;
 double t28 = 108.0*t9;
 double t29 = -t25-t26+t27+t28;

      values[0] = 7.0-9.0*xi-t2-t3-t6+t8+t10;
      values[1] = t12+t13-t14-t15;
      values[2] = -t17-t6+t8+t10;
      values[3] = t21;
      values[4] = 32.0*t22;
      values[5] = t23;
      values[6] = t21;
      values[7] = 32.0*t22;
      values[8] = -8.0+t19+t20+t13-t14-t15;
      values[9] = t23;
      values[10] = t29;
      values[11] = t29;
      values[12] = t25-t26+t27+t28;
      values[13] = 27.0-t25-54.0*eta-54.0*zeta-t26+t27+t28;
      values[14] = 256.0*t22;
    }

static void C_T_B2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{

 double t1 = 6.0*eta;
 double t2 = 6.0*xi;
 double t3 = xi*eta;
 double t4 = 8.0*t3;
 double t5 = 4.0-t1-t2+t4;
 double t6 = 24.0*xi;
 double t7 = 64.0*t3;
 double t8 = t6-t7;
 double t10 = 24.0*eta;
 double t11 = t10-t7;
 double t14 = 216.0*t3;

      values[0] = t5;
      values[1] = t8;
      values[2] = -t2+t4;
      values[3] = t11;
      values[4] = -t7;
      values[5] = -t1+t4;
      values[6] = -8.0+t10+t6-t7;
      values[7] = t8;
      values[8] = t11;
      values[9] = t5;
      values[10] = t14;
      values[11] = -54.0*xi+t14;
      values[12] = t14;
      values[13] = -54.0*eta+t14;
      values[14] = -512.0*t3;

}

TBaseFunct3D *BF_C_T_B2_3D_Obj = 
new TBaseFunct3D(15, BF_C_T_B2_3D, BFUnitTetrahedron, 
                 C_T_B2_3D_Funct, C_T_B2_3D_DeriveXi,
                 C_T_B2_3D_DeriveEta, C_T_B2_3D_DeriveZeta,
                 C_T_B2_3D_DeriveXiXi, C_T_B2_3D_DeriveXiEta,
                 C_T_B2_3D_DeriveXiZeta, C_T_B2_3D_DeriveEtaEta,
                 C_T_B2_3D_DeriveEtaZeta, C_T_B2_3D_DeriveZetaZeta,
                 3, 2,
		 0, NULL);
