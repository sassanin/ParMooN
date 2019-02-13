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
// Q2 element with bubbles, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_UL2SE_2D_Funct(double xi, double eta, double *values)
{
  double t1, t3, t4, t5, t7, t8, t9, t13, t17;

  t1 = xi*eta;
  t3 = xi*xi;
  t4 = 1.0-t3;
  t5 = (1.0-eta)*t4;
  t7 = eta*eta;
  t8 = 1.0-t7;
  t9 = (1.0-xi)*t8;
  t13 = (1.0+xi)*t8;
  t17 = (1.0+eta)*t4;

  values[0] = 1.0/4.0-xi/4.0-eta/4.0+t1/4.0-t5/4.0-t9/4.0;
  values[1] = t5/2.0;
  values[2] = 1.0/4.0+xi/4.0-eta/4.0-t1/4.0-t5/4.0-t13/4.0;
  values[3] = t13/2.0;
  values[4] = 1.0/4.0+xi/4.0+eta/4.0+t1/4.0-t17/4.0-t13/4.0;
  values[5] = t17/2.0;
  values[6] = 1.0/4.0-xi/4.0+eta/4.0-t1/4.0-t17/4.0-t9/4.0;
  values[7] = t9/2.0;
}

// values of the derivatives in xi direction
static void C_Q_UL2SE_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t3, t4, t5, t6, t9, t11, t12;

  t1 = eta/4.0;
  t3 = (1.0-eta)*xi;
  t4 = t3/2.0;
  t5 = eta*eta;
  t6 = t5/4.0;
  t9 = 1.0-t5;
  t11 = (1.0+eta)*xi;
  t12 = t11/2.0;

  values[0] = t1+t4-t6;
  values[1] = -t3;
  values[2] = -t1+t4+t6;
  values[3] = t9/2.0;
  values[4] = t1+t12+t6;
  values[5] = -t11;
  values[6] = -t1+t12-t6;
  values[7] = -t9/2.0;
}

// values of the derivatives in eta direction
static void C_Q_UL2SE_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t5, t6, t8, t10, t11;

  t1 = xi/4.0;
  t2 = xi*xi;
  t3 = t2/4.0;
  t5 = (1.0-xi)*eta;
  t6 = t5/2.0;
  t8 = -1.0+t2;
  t10 = (1.0+xi)*eta;
  t11 = t10/2.0;

  values[0] = t1-t3+t6;
  values[1] = t8/2.0;
  values[2] = -t1-t3+t11;
  values[3] = -t10;
  values[4] = t1+t3+t11;
  values[5] = -t8/2.0;
  values[6] = -t1+t3+t6;
  values[7] = -t5;
}

// values of the derivatives in xi-xi  direction
static void C_Q_UL2SE_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = 1.0-eta;
  t2 = 1.0+eta;

  values[0] = t1/2.0;
  values[1] = -t1;
  values[2] = t1/2.0;
  values[3] = 0.0;
  values[4] = t2/2.0;
  values[5] = -t2;
  values[6] = t2/2.0;
  values[7] = 0.0;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL2SE_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = xi/2.0;
  t2 = eta/2.0;

  values[0] = 1.0/4.0-t1-t2;
  values[1] = xi;
  values[2] = -1.0/4.0-t1+t2;
  values[3] = -eta;
  values[4] = 1.0/4.0+t1+t2;
  values[5] = -xi;
  values[6] = -1.0/4.0+t1-t2;
  values[7] = eta;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL2SE_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = 1.0-xi;
  t2 = 1.0+xi;

  values[0] = t1/2.0;
  values[1] = 0.0;
  values[2] = t2/2.0;
  values[3] = -t2;
  values[4] = t2/2.0;
  values[5] = 0.0;
  values[6] = t1/2.0;
  values[7] = -t1;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_UL2SE_2D_Obj = new TBaseFunct2D
        (8, BF_C_Q_UL2SE_2D, BFUnitSquare, 
         C_Q_UL2SE_2D_Funct, C_Q_UL2SE_2D_DeriveXi,
         C_Q_UL2SE_2D_DeriveEta, C_Q_UL2SE_2D_DeriveXiXi,
         C_Q_UL2SE_2D_DeriveXiEta, C_Q_UL2SE_2D_DeriveEtaEta, 2, 2,
         0, NULL);
