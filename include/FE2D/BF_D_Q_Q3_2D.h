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
// Q3 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_Q3_2D_Funct(double xi, double eta, double *values)
{
  double t1, t8, t10, t11, t15, t18, t22, t34;

  t1 = xi*xi;
  t8 = xi*eta;
  t10 = -1.0+3.0*t1;
  t11 = t10*eta;
  t15 = xi*(5.0*t1-3.0);
  t18 = eta*eta;
  t22 = -1.0+3.0*t18;
  t34 = 5.0*t18-3.0;

  values[0] = 1.0;
  values[1] = xi;
  values[2] = -1.0/2.0+3.0/2.0*t1;
  values[3] = 5.0/2.0*t1*xi-3.0/2.0*xi;
  values[4] = eta;
  values[5] = t8;
  values[6] = t11/2.0;
  values[7] = t15*eta/2.0;
  values[8] = -1.0/2.0+3.0/2.0*t18;
  values[9] = xi*t22/2.0;
  values[10] = t10*t22/4.0;
  values[11] = t15*t22/4.0;
  values[12] = 5.0/2.0*t18*eta-3.0/2.0*eta;
  values[13] = t8*t34/2.0;
  values[14] = t11*t34/4.0;
  values[15] = t15*eta*t34/4.0;
}

// values of the derivatives in xi direction
static void D_Q_Q3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2, t5, t7, t11, t24, t29;

  t2 = xi*xi;
  t5 = xi*eta;
  t7 = t2*eta;
  t11 = eta*eta;
  t24 = 5.0*t11-3.0;
  t29 = t11*eta;

  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 3.0*xi;
  values[3] = 15.0/2.0*t2-3.0/2.0;
  values[4] = 0.0;
  values[5] = eta;
  values[6] = 3.0*t5;
  values[7] = 15.0/2.0*t7-3.0/2.0*eta;
  values[8] = 0.0;
  values[9] = -1.0/2.0+3.0/2.0*t11;
  values[10] = 3.0/2.0*xi*(-1.0+3.0*t11);
  values[11] = -15.0/4.0*t2+45.0/4.0*t2*t11+3.0/4.0-9.0/4.0*t11;
  values[12] = 0.0;
  values[13] = eta*t24/2.0;
  values[14] = 3.0/2.0*t5*t24;
  values[15] = 75.0/4.0*t2*t29-45.0/4.0*t7-15.0/4.0*t29+9.0/4.0*eta;
}

// values of the derivatives in eta direction
static void D_Q_Q3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t6, t17, t20, t29;

  t1 = xi*xi;
  t6 = xi*(5.0*t1-3.0);
  t17 = eta*eta;
  t20 = xi*t17;
  t29 = t1*xi;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 1.0;
  values[5] = xi;
  values[6] = -1.0/2.0+3.0/2.0*t1;
  values[7] = t6/2.0;
  values[8] = 3.0*eta;
  values[9] = 3.0*xi*eta;
  values[10] = 3.0/2.0*(-1.0+3.0*t1)*eta;
  values[11] = 3.0/2.0*t6*eta;
  values[12] = 15.0/2.0*t17-3.0/2.0;
  values[13] = 15.0/2.0*t20-3.0/2.0*xi;
  values[14] = -15.0/4.0*t17+3.0/4.0+45.0/4.0*t1*t17-9.0/4.0*t1;
  values[15] = 75.0/4.0*t29*t17-15.0/4.0*t29-45.0/4.0*t20+9.0/4.0*xi;
}

// values of the derivatives in xi-xi direction
static void D_Q_Q3_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t3, t5;

  t3 = xi*eta;
  t5 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 3.0;
  values[3] = 15.0*xi;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 3.0*eta;
  values[7] = 15.0*t3;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = -3.0/2.0+9.0/2.0*t5;
  values[11] = -15.0/2.0*xi+45.0/2.0*xi*t5;
  values[12] = 0.0;
  values[13] = 0.0;
  values[14] = 3.0/2.0*eta*(5.0*t5-3.0);
  values[15] = 75.0/2.0*xi*t5*eta-45.0/2.0*t3;
}

// values of the derivatives in eta-eta direction
static void D_Q_Q3_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t2, t12;

  t2 = xi*xi;
  t12 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 1.0;
  values[6] = 3.0*xi;
  values[7] = 15.0/2.0*t2-3.0/2.0;
  values[8] = 0.0;
  values[9] = 3.0*eta;
  values[10] = 9.0*xi*eta;
  values[11] = 45.0/2.0*t2*eta-9.0/2.0*eta;
  values[12] = 0.0;
  values[13] = 15.0/2.0*t12-3.0/2.0;
  values[14] = 45.0/2.0*xi*t12-9.0/2.0*xi;
  values[15] = 225.0/4.0*t2*t12-45.0/4.0*t2-45.0/4.0*t12+9.0/4.0;
}

// values of the derivatives in xi-eta direction
static void D_Q_Q3_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t2, t10;

  t2 = xi*xi;
  t10 = xi*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 3.0;
  values[9] = 3.0*xi;
  values[10] = 9.0/2.0*t2-3.0/2.0;
  values[11] = 3.0/2.0*xi*(5.0*t2-3.0);
  values[12] = 15.0*eta;
  values[13] = 15.0*t10;
  values[14] = -15.0/2.0*eta+45.0/2.0*t2*eta;
  values[15] = 75.0/2.0*t2*xi*eta-45.0/2.0*t10;
}

// ***********************************************************************

TBaseFunct2D *BF_D_Q_Q3_2D_Obj = new TBaseFunct2D
        (16, BF_D_Q_Q3_2D, BFUnitSquare, 
         D_Q_Q3_2D_Funct, D_Q_Q3_2D_DeriveXi,
         D_Q_Q3_2D_DeriveEta, D_Q_Q3_2D_DeriveXiXi,
         D_Q_Q3_2D_DeriveXiEta, D_Q_Q3_2D_DeriveEtaEta, 3, 3,
         0, NULL);
