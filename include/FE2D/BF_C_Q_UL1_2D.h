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
// Q1 element with bubble, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_UL1_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t4, t5, t6, t7, t8, t9, t10;

  t1 = eta/4.0;
  t2 = xi/4.0;
  t4 = xi*eta/4.0;
  t5 = eta*eta;
  t6 = 9.0/16.0*t5;
  t7 = xi*xi;
  t8 = 9.0/16.0*t7;
  t9 = t7*t5;
  t10 = 9.0/16.0*t9;

  values[0] = -5.0/16.0-t1-t2+t4+t6+t8-t10;
  values[1] = -5.0/16.0-t1+t2-t4+t6+t8-t10;
  values[2] = -5.0/16.0+t1+t2+t4+t6+t8-t10;
  values[3] = -5.0/16.0+t1-t2-t4+t6+t8-t10;
  values[4] = 9.0/16.0-9.0/16.0*t5-9.0/16.0*t7+9.0/16.0*t9;
}

// values of the derivatives in xi direction
static void C_Q_UL1_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5;

  t1 = eta/4.0;
  t2 = 9.0/8.0*xi;
  t3 = eta*eta;
  t4 = xi*t3;
  t5 = 9.0/8.0*t4;

  values[0] = -1.0/4.0+t1+t2-t5;
  values[1] = 1.0/4.0-t1+t2-t5;
  values[2] = 1.0/4.0+t1+t2-t5;
  values[3] = -1.0/4.0-t1+t2-t5;
  values[4] = -9.0/8.0*xi+9.0/8.0*t4;
}

// values of the derivatives in eta direction
static void C_Q_UL1_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5;

  t1 = xi/4.0;
  t2 = 9.0/8.0*eta;
  t3 = xi*xi;
  t4 = t3*eta;
  t5 = 9.0/8.0*t4;

  values[0] = -1.0/4.0+t1+t2-t5;
  values[1] = -1.0/4.0-t1+t2-t5;
  values[2] = 1.0/4.0+t1+t2-t5;
  values[3] = 1.0/4.0-t1+t2-t5;
  values[4] = -9.0/8.0*eta+9.0/8.0*t4;
}
// values of the derivatives in xi-xi  direction
static void C_Q_UL1_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = eta*eta;
  t2 = 1.0-t1;

  values[0] = 9.0/8.0*t2;
  values[1] = 9.0/8.0*t2;
  values[2] = 9.0/8.0*t2;
  values[3] = 9.0/8.0*t2;
  values[4] = -9.0/8.0*t2;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL1_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t2, t3, t4;

  t2 = 9.0/4.0*xi*eta;
  t3 = 1.0/4.0-t2;
  t4 = -1.0/4.0-t2;

  values[0] = t3;
  values[1] = t4;
  values[2] = t3;
  values[3] = t4;
  values[4] = t2;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL1_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2;

  t1 = xi*xi;
  t2 = 1.0-t1;

  values[0] = 9.0/8.0*t2;
  values[1] = 9.0/8.0*t2;
  values[2] = 9.0/8.0*t2;
  values[3] = 9.0/8.0*t2;
  values[4] = -9.0/8.0*t2;
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_UL1_2D_Obj = new TBaseFunct2D
        (5, BF_C_Q_UL1_2D, BFUnitSquare, 
         C_Q_UL1_2D_Funct, C_Q_UL1_2D_DeriveXi,
         C_Q_UL1_2D_DeriveEta, C_Q_UL1_2D_DeriveXiXi,
         C_Q_UL1_2D_DeriveXiEta, C_Q_UL1_2D_DeriveEtaEta, 2, 1,
         0, NULL);
