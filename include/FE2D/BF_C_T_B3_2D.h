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
// conforming P3 element with cell bubbles
// ***********************************************************************

// base function values
static void C_T_B3_2D_Funct(double xi, double eta, double *values)
{
  double t1 = xi*eta;
  double t2 = xi*xi;
  double t3 = t2*eta;
  double t4 = eta*eta;
  double t5 = xi*t4;
  double t6 = t2*xi;
  double t7 = t6*eta;
  double t8 = t2*t4;
  double t9 = t4*eta;
  double t10 = xi*t9;
  double t12 = 1.0-11.0/2.0*xi-11.0/2.0*eta+9.0*t2+9.0*t1+9.0*t4
              -9.0/2.0*t6+6.0*t3+6.0*t5-9.0/2.0*t9-21.0/2.0*t7
              -21.0*t8-21.0/2.0*t10;

  values[0] = t12;
  values[1] = 9.0*xi-45.0/2.0*t2-90.0*t1+27.0/2.0*t6+189.0*t3
             +351.0/2.0*t5-189.0/2.0*t7-189.0*t8-189.0/2.0*t10;
  values[2] = -9.0/2.0*xi+18.0*t2+63.0/2.0*t1-27.0/2.0*t6-135.0*t3
              -27.0*t5+189.0/2.0*t7+189.0/2.0*t8;
  values[3] = xi-9.0/2.0*t2+3.0/2.0*t1+9.0/2.0*t6-12.0*t3-3.0/2.0*t5
             +21.0/2.0*t7+21.0/2.0*t8;
  values[4] = 9.0*eta-90.0*t1-45.0/2.0*t4+351.0/2.0*t3+189.0*t5
             +27.0/2.0*t9-189.0/2.0*t7-189.0*t8-189.0/2.0*t10;
  values[5] = 45.0/2.0*t1-108.0*t3-27.0*t5+189.0/2.0*t7+189.0/2.0*t8;
  values[6] = -9.0/2.0*eta+63.0/2.0*t1+18.0*t4-27.0*t3-135.0*t5
              -27.0/2.0*t9+189.0/2.0*t8+189.0/2.0*t10;
  values[7] = 45.0/2.0*t1-27.0*t3-108.0*t5+189.0/2.0*t8+189.0/2.0*t10;
  values[8] = eta+3.0/2.0*t1-9.0/2.0*t4-3.0/2.0*t3-12.0*t5+9.0/2.0*t9
             +21.0/2.0*t8+21.0/2.0*t10;
  values[9] = 5.0*t1-12.0*t3-12.0*t5+7.0*t7+14.0*t8+7.0*t10;
  values[10] = -7.0*t1+21.0*t3+14.0*t5-14.0*t7-21.0*t8-7.0*t10;
  values[11] = -7.0*t1+14.0*t3+21.0*t5-7.0*t7-21.0*t8-14.0*t10;
}

// values of the derivatives in xi direction
static void C_T_B3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1 = xi*eta;
  double t2 = eta*eta;
  double t3 = xi*xi;
  double t4 = t3*eta;
  double t5 = xi*t2;
  double t6 = t2*eta;

  values[0] = -11.0/2.0+18.0*xi+9.0*eta-27.0/2.0*t3+12.0*t1+6.0*t2
              -63.0/2.0*t4-42.0*t5-21.0/2.0*t6;
  values[1] = 9.0-45.0*xi-90.0*eta+81.0/2.0*t3+378.0*t1+351.0/2.0*t2
             -567.0/2.0*t4-378.0*t5-189.0/2.0*t6;
  values[2] = -9.0/2.0+36.0*xi+63.0/2.0*eta-81.0/2.0*t3-270.0*t1-27.0*t2
              +567.0/2.0*t4+189.0*t5;
  values[3] = 1.0-9.0*xi+3.0/2.0*eta+27.0/2.0*t3-24.0*t1-3.0/2.0*t2
             +63.0/2.0*t4+21.0*t5;
  values[4] = -90.0*eta+351.0*t1+189.0*t2-567.0/2.0*t4-378.0*t5-189.0/2.0*t6;
  values[5] = 45.0/2.0*eta-216.0*t1-27.0*t2+567.0/2.0*t4+189.0*t5;
  values[6] = 63.0/2.0*eta-54.0*t1-135.0*t2+189.0*t5+189.0/2.0*t6;
  values[7] = 45.0/2.0*eta-54.0*t1-108.0*t2+189.0*t5+189.0/2.0*t6;
  values[8] = 3.0/2.0*eta-3.0*t1-12.0*t2+21.0*t5+21.0/2.0*t6;
  values[9] = 5.0*eta-24.0*t1-12.0*t2+21.0*t4+28.0*t5+7.0*t6;
  values[10] = -7.0*eta+42.0*t1+14.0*t2-42.0*t4-42.0*t5-7.0*t6;
  values[11] = -7.0*eta+28.0*t1+21.0*t2-21.0*t4-42.0*t5-14.0*t6;
}

// values of the derivatives in eta direction
static void C_T_B3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1 = xi*xi;
  double t2 = xi*eta;
  double t3 = t1*xi;
  double t4 = t1*eta;
  double t5 = eta*eta;
  double t6 = xi*t5;

  values[0] = -11.0/2.0+9.0*xi+18.0*eta+6.0*t1+12.0*t2-27.0/2.0*t5
              -21.0/2.0*t3-42.0*t4-63.0/2.0*t6;
  values[1] = -90.0*xi+189.0*t1+351.0*t2-189.0/2.0*t3-378.0*t4-567.0/2.0*t6;
  values[2] = 63.0/2.0*xi-135.0*t1-54.0*t2+189.0/2.0*t3+189.0*t4;
  values[3] = 3.0/2.0*xi-12.0*t1-3.0*t2+21.0/2.0*t3+21.0*t4;
  values[4] = 9.0-90.0*xi-45.0*eta+351.0/2.0*t1+378.0*t2+81.0/2.0*t5
             -189.0/2.0*t3-378.0*t4-567.0/2.0*t6;
  values[5] = 45.0/2.0*xi-108.0*t1-54.0*t2+189.0/2.0*t3+189.0*t4;
  values[6] = -9.0/2.0+63.0/2.0*xi+36.0*eta-27.0*t1-270.0*t2-81.0/2.0*t5
              +189.0*t4+567.0/2.0*t6;
  values[7] = 45.0/2.0*xi-27.0*t1-216.0*t2+189.0*t4+567.0/2.0*t6;
  values[8] = 1.0+3.0/2.0*xi-9.0*eta-3.0/2.0*t1-24.0*t2+27.0/2.0*t5+21.0*t4
             +63.0/2.0*t6;
  values[9] = 5.0*xi-12.0*t1-24.0*t2+7.0*t3+28.0*t4+21.0*t6;
  values[10] = -7.0*xi+21.0*t1+28.0*t2-14.0*t3-42.0*t4-21.0*t6;
  values[11] = -7.0*xi+14.0*t1+42.0*t2-7.0*t3-42.0*t4-42.0*t6;
}

// values of the derivatives in xi-xi direction
static void C_T_B3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1 = xi*eta;
  double t2 = eta*eta;
  double t10 = -54.0*eta+189.0*t2;

  values[0] = 18.0-27.0*xi+12.0*eta-63.0*t1-42.0*t2;
  values[1] = -45.0+81.0*xi+378.0*eta-567.0*t1-378.0*t2;
  values[2] = 36.0-81.0*xi-270.0*eta+567.0*t1+189.0*t2;
  values[3] = -9.0+27.0*xi-24.0*eta+63.0*t1+21.0*t2;
  values[4] = 351.0*eta-567.0*t1-378.0*t2;
  values[5] = -216.0*eta+567.0*t1+189.0*t2;
  values[6] = t10;
  values[7] = t10;
  values[8] = -3.0*eta+21.0*t2;
  values[9] = -24.0*eta+42.0*t1+28.0*t2;
  values[10] = 42.0*eta-84.0*t1-42.0*t2;
  values[11] = 28.0*eta-42.0*t1-42.0*t2;
}

// values of the derivatives in xi-eta direction
static void C_T_B3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = xi*xi;
  double t2 = xi*eta;
  double t3 = eta*eta;

  values[0] = 9.0+12.0*xi+12.0*eta-63.0/2.0*t1-84.0*t2-63.0/2.0*t3;
  values[1] = -90.0+378.0*xi+351.0*eta-567.0/2.0*t1-756.0*t2-567.0/2.0*t3;
  values[2] = 63.0/2.0-270.0*xi-54.0*eta+567.0/2.0*t1+378.0*t2;
  values[3] = 3.0/2.0-24.0*xi-3.0*eta+63.0/2.0*t1+42.0*t2;
  values[4] = -90.0+351.0*xi+378.0*eta-567.0/2.0*t1-756.0*t2-567.0/2.0*t3;
  values[5] = 45.0/2.0-216.0*xi-54.0*eta+567.0/2.0*t1+378.0*t2;
  values[6] = 63.0/2.0-54.0*xi-270.0*eta+378.0*t2+567.0/2.0*t3;
  values[7] = 45.0/2.0-54.0*xi-216.0*eta+378.0*t2+567.0/2.0*t3;
  values[8] = 3.0/2.0-3.0*xi-24.0*eta+42.0*t2+63.0/2.0*t3;
  values[9] = 5.0-24.0*xi-24.0*eta+21.0*t1+56.0*t2+21.0*t3;
  values[10] = -7.0+42.0*xi+28.0*eta-42.0*t1-84.0*t2-21.0*t3;
  values[11] = -7.0+28.0*xi+42.0*eta-21.0*t1-84.0*t2-42.0*t3;
}

// values of the derivatives in eta-eta direction
static void C_T_B3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1 = xi*xi;
  double t2 = xi*eta;
  double t6 = -54.0*xi+189.0*t1;

  values[0] = 18.0+12.0*xi-27.0*eta-42.0*t1-63.0*t2;
  values[1] = 351.0*xi-378.0*t1-567.0*t2;
  values[2] = t6;
  values[3] = -3.0*xi+21.0*t1;
  values[4] = -45.0+378.0*xi+81.0*eta-378.0*t1-567.0*t2;
  values[5] = t6;
  values[6] = 36.0-270.0*xi-81.0*eta+189.0*t1+567.0*t2;
  values[7] = -216.0*xi+189.0*t1+567.0*t2;
  values[8] = -9.0-24.0*xi+27.0*eta+21.0*t1+63.0*t2;
  values[9] = -24.0*xi+28.0*t1+42.0*t2;
  values[10] = 28.0*xi-42.0*t1-42.0*t2;
  values[11] = 42.0*xi-42.0*t1-84.0*t2;
}
  
// ***********************************************************************

TBaseFunct2D *BF_C_T_B3_2D_Obj = new TBaseFunct2D
        (12, BF_C_T_B3_2D, BFUnitTriangle, 
         C_T_B3_2D_Funct, C_T_B3_2D_DeriveXi,
         C_T_B3_2D_DeriveEta, C_T_B3_2D_DeriveXiXi,
         C_T_B3_2D_DeriveXiEta, C_T_B3_2D_DeriveEtaEta, 4, 3,
         0, NULL);
