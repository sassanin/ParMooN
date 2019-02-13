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
// P3 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_P3_2D_Funct(double xi, double eta, double *values)
{
  double t1 = xi+eta-2.0/3.0;
  double t3 = xi+eta-1.0;
  double t6 = xi*(xi-1.0/3.0);
  double t16 = eta*(eta-1.0/3.0);
  values[0] = -9.0/2.0*(xi+eta-1.0/3.0)*t1*t3;
  values[1] = 27.0/2.0*xi*t1*t3;
  values[2] = -27.0/2.0*t6*t3;
  values[3] = 9.0/2.0*t6*(xi-2.0/3.0);
  values[4] = 27.0/2.0*t1*t3*eta;
  values[5] = -27.0*xi*eta*t3;
  values[6] = 27.0/2.0*t6*eta;
  values[7] = -27.0/2.0*t16*t3;
  values[8] = 27.0/2.0*t16*xi;
  values[9] = 9.0/2.0*t16*(eta-2.0/3.0);
}

// values of the derivatives in xi direction
static void C_T_P3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1 = xi+eta-2.0/3.0;
  double t2 = xi+eta-1.0;
  double t3 = t1*t2;
  double t4 = xi*t2;
  double t7 = xi-1.0/3.0;
  double t9 = xi*t7;
  double t11 = xi-2.0/3.0;
  double t15 = t2*eta;
  double t18 = xi*eta;
  double t23 = eta*(eta-1.0/3.0);
  double t24 = xi+eta-1.0/3.0;
  values[0] = -9.0/2.0*t3-9.0/2.0*t24*t2-9.0/2.0*t24*t1;
  values[1] = 27.0/2.0*t3+27.0/2.0*t4+27.0/2.0*xi*t1;
  values[2] = -27.0/2.0*t7*t2-27.0/2.0*t4-27.0/2.0*t9;
  values[3] = 9.0/2.0*t7*t11+9.0/2.0*xi*t11+9.0/2.0*t9;
  values[4] = 27.0/2.0*t15+27.0/2.0*t1*eta;
  values[5] = -27.0*t15-27.0*t18;
  values[6] = 27.0/2.0*t18+27.0/2.0*t7*eta;
  values[7] = -27.0/2.0*t23;
  values[8] = 27.0/2.0*t23;
  values[9] = 0.0;
}

// values of the derivatives in eta direction
static void C_T_P3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1 = xi+eta-1.0;
  double t2 = xi*t1;
  double t3 = xi+eta-2.0/3.0;
  double t7 = xi*(xi-1.0/3.0);
  double t8 = t1*eta;
  double t10 = t1*t3;
  double t12 = xi*eta;
  double t14 = eta-1.0/3.0;
  double t16 = eta*t14;
  double t20 = eta-2.0/3.0;
  double t24 = xi+eta-1.0/3.0;
  values[0] = -9.0/2.0*t10-9.0/2.0*t24*t1-9.0/2.0*t24*t3;
  values[1] = 27.0/2.0*t2+27.0/2.0*xi*t3;
  values[2] = -27.0/2.0*t7;
  values[3] = 0.0;
  values[4] = 27.0/2.0*t8+27.0/2.0*t3*eta+27.0/2.0*t10;
  values[5] = -27.0*t2-27.0*t12;
  values[6] = 27.0/2.0*t7;
  values[7] = -27.0/2.0*t14*t1-27.0/2.0*t8-27.0/2.0*t16;
  values[8] = 27.0/2.0*t14*xi+27.0/2.0*t12;
  values[9] = 9.0/2.0*t14*t20+9.0/2.0*eta*t20+9.0/2.0*t16;
}
// values of the derivatives in xi-xi direction
static void C_T_P3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1 = xi+eta-2.0/3.0;
  double t2 = xi+eta-1.0;
  double t7 = xi-1.0/3.0;
  double t11 = xi-2.0/3.0;
  double t24 = xi+eta-1.0/3.0;

  values[0] = -9.0*(t1+t2+t24);
  values[1] = 27.0*(t1+t2+xi);
  values[2] = -27.0*(t7+t2+xi);
  values[3] = 9.0*(t7+t11+xi);
  values[4] = 27.0*eta;
  values[5] = -54.0*eta;
  values[6] = 27.0*eta;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
}
// values of the derivatives in xi-eta direction
static void C_T_P3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = xi+eta-2.0/3.0;
  double t2 = xi+eta-1.0;
  double t24 = xi+eta-1.0/3.0;

  values[0] = -9.0*(t1+t2+t24);
  values[1] = 27.0/2.0*(t1+t2)+27.0*xi;
  values[2] = -27.0*xi+9.0/2.0;
  values[3] = 0.0;
  values[4] = 27.0/2.0*(t1+t2)+27.0*eta;
  values[5] = -54.0*(xi+eta-0.5);
  values[6] = 27.0*xi-9.0/2.0;
  values[7] = -27.0*eta+9.0/2.0;
  values[8] = 27.0*eta-9.0/2.0;
  values[9] = 0.0;
}
// values of the derivatives in eta-eta direction
static void C_T_P3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1 = xi+eta-1.0;
  double t3 = xi+eta-2.0/3.0;
  double t14 = eta-1.0/3.0;
  double t20 = eta-2.0/3.0;
  double t24 = xi+eta-1.0/3.0;

  values[0] = -9.0*(t3+t1+t24);
  values[1] = 27.0*xi;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 27.0*(t1+eta+t3);
  values[5] = -54.0*xi;
  values[6] = 0.0;
  values[7] = -27.0*(t14+t1+eta);
  values[8] = 27.0*xi;
  values[9] = 9.0*(t14+t20+eta);
}
  
// ***********************************************************************

TBaseFunct2D *BF_C_T_P3_2D_Obj = new TBaseFunct2D
        (10, BF_C_T_P3_2D, BFUnitTriangle, 
         C_T_P3_2D_Funct, C_T_P3_2D_DeriveXi,
         C_T_P3_2D_DeriveEta, C_T_P3_2D_DeriveXiXi,
         C_T_P3_2D_DeriveXiEta, C_T_P3_2D_DeriveEtaEta, 3, 3,
         0, NULL);
