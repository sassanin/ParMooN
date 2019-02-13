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
// P4 element, discontinous, 2D, triangle
// ***********************************************************************

// base function values
static void D_T_P4_2D_Funct(double xi, double eta, double *values)
{
  double t7, t9, t10, t11, t26, t27, t29, t32, t34, t40, t65, t67, t69;
  double t70, t74, t80, t82, t84, t87, t91, t99, t114, t127, t140;

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
  t65 = xi*t34;
  t67 = t32*eta;
  t69 = t7*t11;
  t70 = 4500.0*t69;
  t74 = 3000.0*t9;
  t80 = t7*t7;
  t82 = t11*t11;
  t84 = 50.0+1000.0*t65+7000.0*t67+t70-1000.0*xi-200.0*eta+4500.0*t7
                +t74-9000.0*t27-3000.0*t29+300.0*t11-7000.0*t32-200.0*t34
                +3500.0*t80+50.0*t82;
  t87 = 9000.0*t69;
  t91 = 4800.0*t9;
  t99 = 50.0+3200.0*t65+8000.0*t67+t87-800.0*xi-400.0*eta+3000.0*t7+t91
                -12000.0*t27-7200.0*t29+900.0*t11-4000.0*t32-800.0*t34
                +1750.0*t80+250.0*t82;
  t114 = 50.0+6000.0*t65+6000.0*t67+10800.0*t69-600.0*xi-600.0*eta
                +1800.0*t7+5400.0*t9-10800.0*t27-10800.0*t29+1800.0*t11
                -2000.0*t32-2000.0*t34+750.0*t80+750.0*t82;
  t127 = 50.0+8000.0*t65+3200.0*t67+t87-400.0*xi-800.0*eta+900.0*t7+t91
                -7200.0*t27-12000.0*t29+3000.0*t11-800.0*t32-4000.0*t34
                +250.0*t80+1750.0*t82;
  t140 = 50.0+7000.0*t65+1000.0*t67+t70-200.0*xi-1000.0*eta+300.0*t7+t74
                -3000.0*t27-9000.0*t29+4500.0*t11-200.0*t32-7000.0*t34
                +50.0*t80+3500.0*t82;

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
  values[10] = t84;
  values[11] = t99;
  values[12] = t114;
  values[13] = t127;
  values[14] = t140;
}

// values of the derivatives in xi direction
static void D_T_P4_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2, t10, t11, t13, t15, t19, t34, t36, t38, t39, t41;
  double t45, t50, t52;

  t2 = 6.0*eta;
  t10 = 768.0*eta;
  t11 = xi*eta;
  t13 = eta*eta;
  t15 = xi*xi;
  t19 = 1152.0*eta;
  t34 = t13*eta;
  t36 = t15*eta;
  t38 = xi*t13;
  t39 = 9000.0*t38;
  t41 = 3000.0*eta;
  t45 = t15*xi;
  t50 = 18000.0*t38;
  t52 = 4800.0*eta;

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
  values[10] = 1000.0*t34+21000.0*t36+t39-1000.0+9000.0*xi+t41-18000.0*t11-3000.0*t13-21000.0*t15+14000.0*t45;
  values[11] = 3200.0*t34+24000.0*t36+t50-800.0+6000.0*xi+t52-24000.0*t11-7200.0*t13-12000.0*t15+7000.0*t45;
  values[12] = 6000.0*t34+18000.0*t36+21600.0*t38-600.0+3600.0*xi+5400.0*eta-21600.0*t11-10800.0*t13-6000.0*t15+3000.0*t45;
  values[13] = 8000.0*t34+9600.0*t36+t50-400.0+1800.0*xi+t52-14400.0*t11-12000.0*t13-2400.0*t15+1000.0*t45;
  values[14] = 7000.0*t34+3000.0*t36+t39-200.0+600.0*xi+t41-6000.0*t11-9000.0*t13-600.0*t15+200.0*t45;
}

// values of the derivatives in eta direction
static void D_T_P4_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t9, t10, t12, t15, t18, t34, t36, t38, t39, t40;
  double t45, t50, t51;

  t1 = 6.0*xi;
  t9 = 768.0*xi;
  t10 = xi*xi;
  t12 = xi*eta;
  t15 = eta*eta;
  t18 = 1152.0*xi;
  t34 = t15*xi;
  t36 = t10*xi;
  t38 = t10*eta;
  t39 = 9000.0*t38;
  t40 = 3000.0*xi;
  t45 = t15*eta;
  t50 = 18000.0*t38;
  t51 = 4800.0*xi;

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
  values[10] = 3000.0*t34+7000.0*t36+t39-200.0+t40-9000.0*t10-6000.0*t12+600.0*eta-600.0*t15+200.0*t45;
  values[11] = 9600.0*t34+8000.0*t36+t50-400.0+t51-12000.0*t10-14400.0*t12+1800.0*eta-2400.0*t15+1000.0*t45;
  values[12] = 18000.0*t34+6000.0*t36+21600.0*t38-600.0+5400.0*xi-10800.0*t10-21600.0*t12+3600.0*eta-6000.0*t15+3000.0*t45;
  values[13] = 24000.0*t34+3200.0*t36+t50-800.0+t51-7200.0*t10-24000.0*t12+6000.0*eta-12000.0*t15+7000.0*t45;
  values[14] = 21000.0*t34+1000.0*t36+t39-1000.0+t40-3000.0*t10-18000.0*t12+9000.0*eta-21000.0*t15+14000.0*t45;
}

// values of the derivatives in xi-xi direction
static void D_T_P4_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t13, t15, t16, t19, t23;

  t13 = xi*eta;
  t15 = eta*eta;
  t16 = 9000.0*t15;
  t19 = xi*xi;
  t23 = 18000.0*t15;

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
  values[10] = 42000.0*t13+t16+9000.0-18000.0*eta-42000.0*xi+42000.0*t19;
  values[11] = 48000.0*t13+t23+6000.0-24000.0*eta-24000.0*xi+21000.0*t19;
  values[12] = 36000.0*t13+21600.0*t15+3600.0-21600.0*eta-12000.0*xi+9000.0*t19;
  values[13] = 19200.0*t13+t23+1800.0-14400.0*eta-4800.0*xi+3000.0*t19;
  values[14] = 6000.0*t13+t16+600.0-6000.0*eta-1200.0*xi+600.0*t19;
}

// values of the derivatives in xi-eta direction
static void D_T_P4_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t13, t15, t17, t18, t24;

  t13 = eta*eta;
  t15 = xi*xi;
  t17 = xi*eta;
  t18 = 18000.0*t17;
  t24 = 36000.0*t17;

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
  values[10] = 3000.0*t13+21000.0*t15+t18+3000.0-18000.0*xi-6000.0*eta;
  values[11] = 9600.0*t13+24000.0*t15+t24+4800.0-24000.0*xi-14400.0*eta;
  values[12] = 18000.0*t13+18000.0*t15+43200.0*t17+5400.0-21600.0*xi-21600.0*eta;
  values[13] = 24000.0*t13+9600.0*t15+t24+4800.0-14400.0*xi-24000.0*eta;
  values[14] = 21000.0*t13+3000.0*t15+t18+3000.0-6000.0*xi-18000.0*eta;
}

// values of the derivatives in eta-eta direction
static void D_T_P4_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t13, t15, t16, t19, t23;

  t13 = xi*eta;
  t15 = xi*xi;
  t16 = 9000.0*t15;
  t19 = eta*eta;
  t23 = 18000.0*t15;

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
  values[10] = 6000.0*t13+t16-6000.0*xi+600.0-1200.0*eta+600.0*t19;
  values[11] = 19200.0*t13+t23-14400.0*xi+1800.0-4800.0*eta+3000.0*t19;
  values[12] = 36000.0*t13+21600.0*t15+3600.0-21600.0*xi-12000.0*eta+9000.0*t19;
  values[13] = 48000.0*t13+t23+6000.0-24000.0*eta-24000.0*xi+21000.0*t19;
  values[14] = 42000.0*t13+t16-18000.0*xi+9000.0-42000.0*eta+42000.0*t19;
}

// ***********************************************************************

TBaseFunct2D *BF_D_T_P4_2D_Obj = new TBaseFunct2D
        (15, BF_D_T_P4_2D, BFUnitTriangle, 
         D_T_P4_2D_Funct, D_T_P4_2D_DeriveXi,
         D_T_P4_2D_DeriveEta, D_T_P4_2D_DeriveXiXi,
         D_T_P4_2D_DeriveXiEta, D_T_P4_2D_DeriveEtaEta, 4, 4,
         0, NULL);
