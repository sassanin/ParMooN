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
static void C_Q_M3_2D_Funct(double xi, double eta, double *values)
{
  double t1 = 19.0/32.0*xi;
  double t2 = 19.0/32.0*eta;
  double t3 = xi*xi;
  double t4 = 9.0/32.0*t3;
  double t6 = xi*eta/4.0;
  double t7 = eta*eta;
  double t8 = 9.0/32.0*t7;
  double t9 = t3*xi;
  double t10 = 9.0/16.0*t9;
  double t12 = 9.0/32.0*t3*eta;
  double t14 = 9.0/32.0*xi*t7;
  double t15 = t7*eta;
  double t16 = 9.0/16.0*t15;
  double t20 = (1.0-t3)*xi*(1.0+eta);
  double t21 = 9.0/32.0*t20;
  double t25 = (1.0-t7)*eta*(1.0+xi);
  double t26 = 9.0/32.0*t25;
  double t27 = -5.0/16.0+t1+t2+t4+t6+t8-t10-t12-t14-t16-t21-t26;
  double t28 = 27.0/16.0*xi;
  double t29 = 9.0/32.0*eta;
  double t30 = 27.0/16.0*t9;
  double t31 = 27.0/32.0*t20;
  double t34 = eta/32.0;
  double t35 = -5.0/16.0-t1+t34+t4-t6+t8+t10-t12+t14+t21+t26;
  double t36 = 9.0/32.0*xi;
  double t37 = 27.0/32.0*t25;
  double t40 = xi/32.0;
  double t44 = -5.0/16.0+t40-t2+t4-t6+t8+t12-t14+t16+t21+t26;
  double t45 = 27.0/16.0*eta;
  double t46 = 27.0/16.0*t15;

  values[0] = t27;
  values[1] = 9.0/32.0-t28-t29-t4+t30+t12+t31;
  values[2] = 9.0/32.0+t28-t29-t4-t30+t12-t31;
  values[3] = t35;
  values[4] = 9.0/32.0+t36-t8-t14-t37;
  values[5] = 9.0/32.0+t36-t8-t14+t37;
  values[6] = -5.0/16.0-t40-t34+t4+t6+t8+t12+t14-t21-t26;
  values[7] = 9.0/32.0+t29-t4-t12+t31;
  values[8] = 9.0/32.0+t29-t4-t12-t31;
  values[9] = t44;
  values[10] = 9.0/32.0-t36+t45-t8+t14-t46-t37;
  values[11] = 9.0/32.0-t36-t45-t8+t14+t46+t37;
}

// values of the derivatives in xi direction
static void C_Q_M3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1 = 9.0/16.0*xi;
  double t2 = eta/4.0;
  double t3 = xi*xi;
  double t4 = 27.0/16.0*t3;
  double t6 = 9.0/16.0*xi*eta;
  double t7 = eta*eta;
  double t8 = 9.0/32.0*t7;
  double t9 = 1.0+eta;
  double t10 = t3*t9;
  double t11 = 9.0/16.0*t10;
  double t13 = (1.0-t3)*t9;
  double t14 = 9.0/32.0*t13;
  double t16 = (1.0-t7)*eta;
  double t17 = 9.0/32.0*t16;
  double t19 = 81.0/16.0*t3;
  double t20 = 27.0/16.0*t10;
  double t21 = 27.0/32.0*t13;
  double t25 = 27.0/32.0*t16;
  double t26 = 9.0/32.0-t8-t25;
  double t27 = 9.0/32.0-t8+t25;

  values[0] = 19.0/32.0+t1+t2-t4-t6-t8+t11-t14-t17;
  values[1] = -27.0/16.0-t1+t19+t6-t20+t21;
  values[2] = 27.0/16.0-t1-t19+t6+t20-t21;
  values[3] = -19.0/32.0+t1-t2+t4-t6+t8-t11+t14+t17;
  values[4] = t26;
  values[5] = t27;
  values[6] = -1.0/32.0+t1+t2+t6+t8+t11-t14-t17;
  values[7] = -t1-t6-t20+t21;
  values[8] = -t1-t6+t20-t21;
  values[9] = 1.0/32.0+t1-t2+t6-t8-t11+t14+t17;
  values[10] = -t27;
  values[11] = -t26;
}

// values of the derivatives in eta direction
static void C_Q_M3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1 = xi/4.0;
  double t2 = 9.0/16.0*eta;
  double t3 = xi*xi;
  double t4 = 9.0/32.0*t3;
  double t6 = 9.0/16.0*xi*eta;
  double t7 = eta*eta;
  double t8 = 27.0/16.0*t7;
  double t10 = (1.0-t3)*xi;
  double t11 = 9.0/32.0*t10;
  double t12 = 1.0+xi;
  double t13 = t7*t12;
  double t14 = 9.0/16.0*t13;
  double t16 = (1.0-t7)*t12;
  double t17 = 9.0/32.0*t16;
  double t19 = 27.0/32.0*t10;
  double t20 = -9.0/32.0+t4+t19;
  double t21 = -9.0/32.0+t4-t19;
  double t23 = 27.0/16.0*t13;
  double t24 = 27.0/32.0*t16;
  double t29 = 81.0/16.0*t7;

  values[0] = 19.0/32.0+t1+t2-t4-t6-t8-t11+t14-t17;
  values[1] = t20;
  values[2] = t21;
  values[3] = 1.0/32.0-t1+t2-t4+t6+t11-t14+t17;
  values[4] = -t2-t6+t23-t24;
  values[5] = -t2-t6-t23+t24;
  values[6] = -1.0/32.0+t1+t2+t4+t6-t11+t14-t17;
  values[7] = -t21;
  values[8] = -t20;
  values[9] = -19.0/32.0-t1+t2+t4-t6+t8+t11-t14+t17;
  values[10] = 27.0/16.0-t2+t6-t29+t23-t24;
  values[11] = -27.0/16.0-t2+t6+t29-t23+t24;
}

// values of the derivatives in xi-xi  direction
static void C_Q_M3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1 = 27.0/8.0*xi;
  double t2 = 9.0/16.0*eta;
  double t4 = xi*(1.0+eta);
  double t5 = 27.0/16.0*t4;
  double t7 = 81.0/8.0*xi;
  double t8 = 81.0/16.0*t4;

  values[0] = 9.0/16.0-t1-t2+t5;
  values[1] = -9.0/16.0+t7+t2-t8;
  values[2] = -9.0/16.0-t7+t2+t8;
  values[3] = 9.0/16.0+t1-t2-t5;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 9.0/16.0+t2+t5;
  values[7] = -9.0/16.0-t2-t8;
  values[8] = -9.0/16.0-t2+t8;
  values[9] = 9.0/16.0+t2-t5;
  values[10] = 0.0;
  values[11] = 0.0;
}

// values of the derivatives in xi-eta direction
static void C_Q_M3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = 9.0/16.0*xi;
  double t2 = 9.0/16.0*eta;
  double t3 = xi*xi;
  double t4 = 27.0/32.0*t3;
  double t5 = eta*eta;
  double t6 = 27.0/32.0*t5;
  double t8 = 81.0/32.0*t3;
  double t9 = t1-t8+27.0/32.0;
  double t10 = t1+t8-27.0/32.0;
  double t12 = 81.0/32.0*t5;
  double t13 = -t2+t12-27.0/32.0;
  double t14 = -t2-t12+27.0/32.0;

  values[0] = -5.0/16.0-t1-t2+t4+t6;
  values[1] = t9;
  values[2] = t10;
  values[3] = 5.0/16.0-t1+t2-t4-t6;
  values[4] = t13;
  values[5] = t14;
  values[6] = -5.0/16.0+t1+t2+t4+t6;
  values[7] = -t10;
  values[8] = -t9;
  values[9] = 5.0/16.0+t1-t2-t4-t6;
  values[10] = -t14;
  values[11] = -t13;
}

// values of the derivatives in eta-eta direction
static void C_Q_M3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1 = 9.0/16.0*xi;
  double t2 = 27.0/8.0*eta;
  double t4 = eta*(1.0+xi);
  double t5 = 27.0/16.0*t4;
  double t8 = 81.0/16.0*t4;
  double t13 = 81.0/8.0*eta;

  values[0] = 9.0/16.0-t1-t2+t5;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 9.0/16.0+t1-t5;
  values[4] = -9.0/16.0-t1+t8;
  values[5] = -9.0/16.0-t1-t8;
  values[6] = 9.0/16.0+t1+t5;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 9.0/16.0-t1+t2-t5;
  values[10] = -9.0/16.0+t1-t13+t8;
  values[11] = -9.0/16.0+t1+t13-t8;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_M3_2D_Obj = new TBaseFunct2D
        (12, BF_C_Q_M3_2D, BFUnitSquare, 
         C_Q_M3_2D_Funct, C_Q_M3_2D_DeriveXi,
         C_Q_M3_2D_DeriveEta, C_Q_M3_2D_DeriveXiXi,
         C_Q_M3_2D_DeriveXiEta, C_Q_M3_2D_DeriveEtaEta, 3, 3,
         0, NULL);
