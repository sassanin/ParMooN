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
// Q3 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_Q3_2D_Funct(double xi, double eta, double *values)
{
  double xi0 = -0.625E-1*(3.0*xi+1.0)*(3.0*xi-1.0)*(xi-1.0);
  double xi1 =  0.5625*(xi+1.0)*(3.0*xi-1.0)*(xi-1.0);
  double xi2 = -0.5625*(xi+1.0)*(3.0*xi+1.0)*(xi-1.0);
  double xi3 = 0.625E-1*(xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0);
  double eta0 = -0.625E-1*(3.0*eta+1.0)*(3.0*eta-1.0)*(eta-1.0);
  double eta1 =  0.5625*(eta+1.0)*(3.0*eta-1.0)*(eta-1.0);
  double eta2 = -0.5625*(eta+1.0)*(3.0*eta+1.0)*(eta-1.0);
  double eta3 = 0.625E-1*(eta+1.0)*(3.0*eta+1.0)*(3.0*eta-1.0);

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi3*eta0;
  values[4]=xi0*eta1;
  values[5]=xi1*eta1;
  values[6]=xi2*eta1;
  values[7]=xi3*eta1;
  values[8]=xi0*eta2;
  values[9]=xi1*eta2;
  values[10]=xi2*eta2;
  values[11]=xi3*eta2;
  values[12]=xi0*eta3;
  values[13]=xi1*eta3;
  values[14]=xi2*eta3;
  values[15]=xi3*eta3;
}

// values of the derivatives in xi direction
static void C_Q_Q3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1 = xi*xi;
  double xi0 = -0.16875E1*t1+0.1125E1*xi+0.625E-1;
  double xi1 = 0.50625E1*t1-0.1125E1*xi-0.16875E1;
  double xi2 = -0.50625E1*t1-0.1125E1*xi+0.16875E1;
  double xi3 = 0.16875E1*t1-0.625E-1+0.1125E1*xi;
  double eta0 = -0.625E-1*(3.0*eta+1.0)*(3.0*eta-1.0)*(eta-1.0);
  double eta1 =  0.5625*(eta+1.0)*(3.0*eta-1.0)*(eta-1.0);
  double eta2 = -0.5625*(eta+1.0)*(3.0*eta+1.0)*(eta-1.0);
  double eta3 = 0.625E-1*(eta+1.0)*(3.0*eta+1.0)*(3.0*eta-1.0);

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi3*eta0;
  values[4]=xi0*eta1;
  values[5]=xi1*eta1;
  values[6]=xi2*eta1;
  values[7]=xi3*eta1;
  values[8]=xi0*eta2;
  values[9]=xi1*eta2;
  values[10]=xi2*eta2;
  values[11]=xi3*eta2;
  values[12]=xi0*eta3;
  values[13]=xi1*eta3;
  values[14]=xi2*eta3;
  values[15]=xi3*eta3;
}

// values of the derivatives in eta direction
static void C_Q_Q3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1=eta*eta;
  double xi0 = -0.625E-1*(3.0*xi+1.0)*(3.0*xi-1.0)*(xi-1.0);
  double xi1 =  0.5625*(xi+1.0)*(3.0*xi-1.0)*(xi-1.0);
  double xi2 = -0.5625*(xi+1.0)*(3.0*xi+1.0)*(xi-1.0);
  double xi3 = 0.625E-1*(xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0);
  double eta0 = -0.16875E1*t1+0.1125E1*eta+0.625E-1;
  double eta1 = 0.50625E1*t1-0.1125E1*eta-0.16875E1;
  double eta2 = -0.50625E1*t1-0.1125E1*eta+0.16875E1;
  double eta3 = 0.16875E1*t1-0.625E-1+0.1125E1*eta;

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi3*eta0;
  values[4]=xi0*eta1;
  values[5]=xi1*eta1;
  values[6]=xi2*eta1;
  values[7]=xi3*eta1;
  values[8]=xi0*eta2;
  values[9]=xi1*eta2;
  values[10]=xi2*eta2;
  values[11]=xi3*eta2;
  values[12]=xi0*eta3;
  values[13]=xi1*eta3;
  values[14]=xi2*eta3;
  values[15]=xi3*eta3;
}
// values of the derivatives in xi-xi  direction
static void C_Q_Q3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double xi0 = -0.3375E1*xi+0.1125E1;
  double xi1 = 1.0125E1*xi-0.1125E1;
  double xi2 = -1.0125E1*xi-0.1125E1;
  double xi3 = 0.3375E1*xi+0.1125E1;
  double eta0 = -0.625E-1*(3.0*eta+1.0)*(3.0*eta-1.0)*(eta-1.0);
  double eta1 =  0.5625*(eta+1.0)*(3.0*eta-1.0)*(eta-1.0);
  double eta2 = -0.5625*(eta+1.0)*(3.0*eta+1.0)*(eta-1.0);
  double eta3 = 0.625E-1*(eta+1.0)*(3.0*eta+1.0)*(3.0*eta-1.0);

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi3*eta0;
  values[4]=xi0*eta1;
  values[5]=xi1*eta1;
  values[6]=xi2*eta1;
  values[7]=xi3*eta1;
  values[8]=xi0*eta2;
  values[9]=xi1*eta2;
  values[10]=xi2*eta2;
  values[11]=xi3*eta2;
  values[12]=xi0*eta3;
  values[13]=xi1*eta3;
  values[14]=xi2*eta3;
  values[15]=xi3*eta3;
}
// values of the derivatives in xi-eta direction
static void C_Q_Q3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = xi*xi;
  double xi0 = -0.16875E1*t1+0.1125E1*xi+0.625E-1;
  double xi1 = 0.50625E1*t1-0.1125E1*xi-0.16875E1;
  double xi2 = -0.50625E1*t1-0.1125E1*xi+0.16875E1;
  double xi3 = 0.16875E1*t1-0.625E-1+0.1125E1*xi;
  double t2=eta*eta;
  double eta0 = -0.16875E1*t2+0.1125E1*eta+0.625E-1;
  double eta1 = 0.50625E1*t2-0.1125E1*eta-0.16875E1;
  double eta2 = -0.50625E1*t2-0.1125E1*eta+0.16875E1;
  double eta3 = 0.16875E1*t2-0.625E-1+0.1125E1*eta;

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi3*eta0;
  values[4]=xi0*eta1;
  values[5]=xi1*eta1;
  values[6]=xi2*eta1;
  values[7]=xi3*eta1;
  values[8]=xi0*eta2;
  values[9]=xi1*eta2;
  values[10]=xi2*eta2;
  values[11]=xi3*eta2;
  values[12]=xi0*eta3;
  values[13]=xi1*eta3;
  values[14]=xi2*eta3;
  values[15]=xi3*eta3;
}
// values of the derivatives in eta-eta direction
static void C_Q_Q3_2D_DeriveEtaEta(double xi, double eta, double *values)
{

  double xi0 = -0.625E-1*(3.0*xi+1.0)*(3.0*xi-1.0)*(xi-1.0);
  double xi1 =  0.5625*(xi+1.0)*(3.0*xi-1.0)*(xi-1.0);
  double xi2 = -0.5625*(xi+1.0)*(3.0*xi+1.0)*(xi-1.0);
  double xi3 = 0.625E-1*(xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0);
  double eta0 = -0.3375E1*eta+0.1125E1;
  double eta1 = 1.0125E1*eta-0.1125E1;
  double eta2 = -1.0125E1*eta-0.1125E1;
  double eta3 = 0.3375E1*eta+0.1125E1;

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi3*eta0;
  values[4]=xi0*eta1;
  values[5]=xi1*eta1;
  values[6]=xi2*eta1;
  values[7]=xi3*eta1;
  values[8]=xi0*eta2;
  values[9]=xi1*eta2;
  values[10]=xi2*eta2;
  values[11]=xi3*eta2;
  values[12]=xi0*eta3;
  values[13]=xi1*eta3;
  values[14]=xi2*eta3;
  values[15]=xi3*eta3;  
}

// ***********************************************************************

TBaseFunct2D *BF_C_Q_Q3_2D_Obj = new TBaseFunct2D
        (16, BF_C_Q_Q3_2D, BFUnitSquare, 
         C_Q_Q3_2D_Funct, C_Q_Q3_2D_DeriveXi,
         C_Q_Q3_2D_DeriveEta, C_Q_Q3_2D_DeriveXiXi,
         C_Q_Q3_2D_DeriveXiEta, C_Q_Q3_2D_DeriveEtaEta, 3, 3,
         0, NULL);
