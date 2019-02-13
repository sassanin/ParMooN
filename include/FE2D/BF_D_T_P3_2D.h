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
// P3 element, discontinous, 2D, triangle
// ***********************************************************************

// base function values
static void D_T_P3_2D_Funct(double xi, double eta, double *values)
{
  double t7, t9, t10, t11, t26, t27, t29, t32, t34, t40;
  
  t7 = xi*xi;
  t9 = xi*eta;
  t10 = 6.0*t9;
  t11 = eta*eta;
  t26 = 768.0*t9;
  t27 = t7*eta;
  t29 = xi*t11;
  t32 = t7*xi;
  t34 = t11*eta;
  t40 = 1152.0*t9;

  values[0] = 1.0;
  values[1] = -8.0+24.0*xi;
  values[2] = -8.0+24.0*eta;
  values[3] = 1.0-6.0*xi-2.0*eta+6.0*t7+t10+t11;
  values[4] = 1.0-4.0*xi-4.0*eta+3.0*t7+8.0*t9+3.0*t11;
  values[5] = 1.0-2.0*xi-6.0*eta+t7+t10+6.0*t11;
  values[6] = -32.0+384.0*xi+96.0*eta-960.0*t7-t26+960.0*t27+384.0*t29
                -96.0*t11+640.0*t32+32.0*t34;
  values[7] = -32.0+288.0*xi+192.0*eta-576.0*t7-t40+1152.0*t27+864.0*t29
                -288.0*t11+320.0*t32+128.0*t34;
  values[8] = -32.0+192.0*xi+288.0*eta-288.0*t7-t40+864.0*t27+1152.0*t29
                -576.0*t11+128.0*t32+320.0*t34;
  values[9] = -32.0+96.0*xi+384.0*eta-96.0*t7-t26+384.0*t27+960.0*t29
                -960.0*t11+32.0*t32+640.0*t34;
}

// values of the derivatives in xi direction
static void D_T_P3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2, t10, t11, t13, t15, t19;

  t2 = 6.0*eta;
  t10 = 768.0*eta;
  t11 = xi*eta;
  t13 = eta*eta;
  t15 = xi*xi;
  t19 = 1152.0*eta;

  values[0] = 0.0;
  values[1] = 24.0;
  values[2] = 0.0;
  values[3] = -6.0+12.0*xi+t2;
  values[4] = -4.0+6.0*xi+8.0*eta;
  values[5] = -2.0+2.0*xi+t2;
  values[6] = 384.0-1920.0*xi-t10+1920.0*t11+384.0*t13+1920.0*t15;
  values[7] = 288.0-1152.0*xi-t19+2304.0*t11+864.0*t13+960.0*t15;
  values[8] = 192.0-576.0*xi-t19+1728.0*t11+1152.0*t13+384.0*t15;
  values[9] = 96.0-192.0*xi-t10+768.0*t11+960.0*t13+96.0*t15;
}

// values of the derivatives in eta direction
static void D_T_P3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t9, t10, t12, t15, t18;

  t1 = 6.0*xi;
  t9 = 768.0*xi;
  t10 = xi*xi;
  t12 = xi*eta;
  t15 = eta*eta;
  t18 = 1152.0*xi;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 24.0;
  values[3] = -2.0+t1+2.0*eta;
  values[4] = -4.0+8.0*xi+6.0*eta;
  values[5] = -6.0+t1+12.0*eta;
  values[6] = 96.0-t9+960.0*t10+768.0*t12-192.0*eta+96.0*t15;
  values[7] = 192.0-t18+1152.0*t10+1728.0*t12-576.0*eta+384.0*t15;
  values[8] = 288.0-t18+864.0*t10+2304.0*t12-1152.0*eta+960.0*t15;
  values[9] = 384.0-t9+384.0*t10+1920.0*t12-1920.0*eta+1920.0*t15;
}

// values of the derivatives in xi-xi direction
static void D_T_P3_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 12.0;
  values[4] = 6.0;
  values[5] = 2.0;
  values[6] = -1920.0+1920.0*eta+3840.0*xi;
  values[7] = -1152.0+2304.0*eta+1920.0*xi;
  values[8] = -576.0+1728.0*eta+768.0*xi;
  values[9] = -192.0+768.0*eta+192.0*xi;
}

// values of the derivatives in xi-eta direction
static void D_T_P3_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 6.0;
  values[4] = 8.0;
  values[5] = 6.0;
  values[6] = -768.0+1920.0*xi+768.0*eta;
  values[7] = -1152.0+2304.0*xi+1728.0*eta;
  values[8] = -1152.0+1728.0*xi+2304.0*eta;
  values[9] = -768.0+768.0*xi+1920.0*eta;
}

// values of the derivatives in eta-eta direction
static void D_T_P3_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 2.0;
  values[4] = 6.0;
  values[5] = 12.0;
  values[6] = 768.0*xi-192.0+192.0*eta;
  values[7] = 1728.0*xi-576.0+768.0*eta;
  values[8] = 2304.0*xi-1152.0+1920.0*eta;
  values[9] = 1920.0*xi-1920.0+3840.0*eta;
}

// ***********************************************************************

TBaseFunct2D *BF_D_T_P3_2D_Obj = new TBaseFunct2D
        (10, BF_D_T_P3_2D, BFUnitTriangle, 
         D_T_P3_2D_Funct, D_T_P3_2D_DeriveXi,
         D_T_P3_2D_DeriveEta, D_T_P3_2D_DeriveXiXi,
         D_T_P3_2D_DeriveXiEta, D_T_P3_2D_DeriveEtaEta, 3, 3,
         0, NULL);
