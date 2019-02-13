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
// P2 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_T_P2_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t3, t5, t6, t7, t8, t9, t12, t14, t15, t18, t19;
  double t24, t29, t31, t32;

  t1 = 6.0*eta;
  t2 = eta*eta;
  t3 = 6.0*t2;
  t5 = 6.0*xi;
  t6 = xi*xi;
  t7 = 6.0*t6;
  t8 = xi*eta;
  t9 = 12.0*t8;
  t12 = 4.0/3.0*xi;
  t14 = 4.0*t8;
  t15 = 2.0*t2;
  t18 = 10.0/3.0*t8*(xi-eta);
  t19 = 1.0-xi-eta;
  t24 = 10.0/3.0*eta*t19*(2.0*eta-1.0+xi);
  t29 = 10.0/3.0*t19*xi*(1.0-2.0*xi-eta);
  t31 = 4.0/3.0*eta;
  t32 = 2.0*t6;

  values[0] = 1.0-t1+t3;
  values[1] = 1.0-t5-t1+t7+t9+t3;
  values[2] = 1.0-t5+t7;
  values[3] = -2.0/3.0+t12+8.0/3.0*eta-t14-t15-t18-t24-t29;
  values[4] = t12-t31-t32+t15-t18-t24-t29;
  values[5] = 2.0/3.0-8.0/3.0*xi-t31+t32+t14-t18-t24-t29;
  values[6] = -2.0+12.0*xi+12.0*eta-12.0*t6-t9-12.0*t2;
}

// values of the derivatives in xi direction
static void N_T_P2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t5, t8, t10, t14, t15, t17, t19, t21, t23, t25, t27;

  t1 = 12.0*xi;
  t2 = 12.0*eta;
  t5 = 4.0*eta;
  t8 = 10.0/3.0*eta*(xi-eta);
  t10 = 10.0/3.0*xi*eta;
  t14 = 10.0/3.0*eta*(2.0*eta-1.0+xi);
  t15 = 1.0-xi-eta;
  t17 = 10.0/3.0*eta*t15;
  t19 = 1.0-2.0*xi-eta;
  t21 = 10.0/3.0*xi*t19;
  t23 = 10.0/3.0*t15*t19;
  t25 = 20.0/3.0*t15*xi;
  t27 = 4.0*xi;

  values[0] = 0.0;
  values[1] = -6.0+t1+t2;
  values[2] = -6.0+t1;
  values[3] = 4.0/3.0-t5-t8-t10+t14-t17+t21-t23+t25;
  values[4] = 4.0/3.0-t27-t8-t10+t14-t17+t21-t23+t25;
  values[5] = -8.0/3.0+t27+t5-t8-t10+t14-t17+t21-t23+t25;
  values[6] = 12.0-24.0*xi-t2;
}

// values of the derivatives in eta direction
static void N_T_P2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t3, t5, t6, t9, t11, t12, t14, t16, t18, t20, t24, t26;

  t1 = 12.0*eta;
  t3 = 12.0*xi;
  t5 = 4.0*xi;
  t6 = 4.0*eta;
  t9 = 10.0/3.0*xi*(xi-eta);
  t11 = 10.0/3.0*xi*eta;
  t12 = 1.0-xi-eta;
  t14 = 2.0*eta-1.0+xi;
  t16 = 10.0/3.0*t12*t14;
  t18 = 10.0/3.0*eta*t14;
  t20 = 20.0/3.0*eta*t12;
  t24 = 10.0/3.0*xi*(1.0-2.0*xi-eta);
  t26 = 10.0/3.0*t12*xi;

  values[0] = -6.0+t1;
  values[1] = -6.0+t3+t1;
  values[2] = 0.0;
  values[3] = 8.0/3.0-t5-t6-t9+t11-t16+t18-t20+t24+t26;
  values[4] = -4.0/3.0+t6-t9+t11-t16+t18-t20+t24+t26;
  values[5] = -4.0/3.0+t5-t9+t11-t16+t18-t20+t24+t26;
  values[6] = 12.0-t3-24.0*eta;
}

// values of derivatives in xi-xi direction
static void N_T_P2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = 20.0*eta;
  t2 = 40.0*xi;

  values[0] = 0.0;
  values[1] = 12.0;
  values[2] = 12.0;
  values[3] = -t1+20.0-t2;
  values[4] = 16.0-t1-t2;
  values[5] = 24.0-t1-t2;
  values[6] = -24.0;
}

// values of derivatives in eta-eta direction
static void N_T_P2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = 20.0*xi;
  t2 = 20.0*eta;

  values[0] = 0.0;
  values[1] = 12.0;
  values[2] = 0.0;
  values[3] = -4.0-t1+t2;
  values[4] = -20.0*xi+20.0*eta;
  values[5] = 4.0-t1+t2;
  values[6] = -12.0;
}

// values of derivatives in xi-eta direction
static void N_T_P2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = 20.0*xi;
  t2 = 40.0*eta;

  values[0] = 12.0;
  values[1] = 12.0;
  values[2] = 0.0;
  values[3] = -24.0+t1+t2;
  values[4] = -16.0+t1+t2;
  values[5] = t1-20.0+t2;
  values[6] = -24.0;
}

static int N_T_P2_2D_ChangeJ0[1] = { 3 };
static int N_T_P2_2D_ChangeJ1[1] = { 4 };
static int N_T_P2_2D_ChangeJ2[1] = { 5 };

static int *N_T_P2_2D_Change[3] = { N_T_P2_2D_ChangeJ0, N_T_P2_2D_ChangeJ1,
                                    N_T_P2_2D_ChangeJ2 };

// ***********************************************************************

TBaseFunct2D *BF_N_T_P2_2D_Obj = new TBaseFunct2D
        (7, BF_N_T_P2_2D, BFUnitTriangle, 
         N_T_P2_2D_Funct, N_T_P2_2D_DeriveXi,
         N_T_P2_2D_DeriveEta, N_T_P2_2D_DeriveXiXi,
         N_T_P2_2D_DeriveXiEta, N_T_P2_2D_DeriveEtaEta, 3, 2,
         1, N_T_P2_2D_Change);
