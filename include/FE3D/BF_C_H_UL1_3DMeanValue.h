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
// Q1 + bubble element for LPS, conforming, 3D
// Author : Sashi
// History: 25.06.2011 
// ***********************************************************************

static void C_H_UL1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{

 double t1, t2, t3, t4, t7, t9, t12, t14, t15, t17, t18, t22,t23, t27, t31;

      t1 = 1.0-xi;
      t2 = 1.0-eta;
      t3 = t1*t2;
      t4 = 1.0-zeta;
      t7 = xi*xi;
      t9 = eta*eta;
      t12 = zeta*zeta;
      t14 = (1.0-t7)*(1.0-t9)*(1.0-t12);
      t15 = 0.4218750001384277*t14;
      t17 = 1.0+xi;
      t18 = t17*t2;
      t22 = 1.0+eta;
      t23 = t1*t22;
      t27 = t17*t22;
      t31 = 1.0+zeta;
        
      values[0] = 0.125*t3*t4-t15;
      values[1] = 0.125*t18*t4-t15;
      values[2] = 0.125*t23*t4-t15;
      values[3] = 0.125*t27*t4-t15;
      values[4] = 0.125*t3*t31-t15;
      values[5] = 0.125*t18*t31-t15;
      values[6] = 0.125*t23*t31-t15;
      values[7] = 0.125*t27*t31-t15;
      values[8] = 0.4218750002438965*t14;
}

static void C_H_UL1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t5, t8, t10, t11, t15,t16, t21, t22,  t27;
  
      t1 = 1.0-eta;
      t2 = 1.0-zeta;
      t3 = t1*t2;
      t5 = eta*eta;
      t8 = zeta*zeta;
      t10 = xi*(1.0-t5)*(1.0-t8);
      t11 = 0.8437500002768554*t10;
      t15 = 1.0+eta;
      t16 = t15*t2;
      t21 = 1.0+zeta;
      t22 = t1*t21;
      t27 = t15*t21;  
      
      values[0] = -0.125*t3+t11;
      values[1] = 0.125*t3+t11;
      values[2] = -0.125*t16+t11;
      values[3] = 0.125*t16+t11;
      values[4] = -0.125*t22+t11;
      values[5] = 0.125*t22+t11;
      values[6] = -0.125*t27+t11;
      values[7] = 0.125*t27+t11;
      values[8] = -0.8437500004877929*t10;      
  
}

static void C_H_UL1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
 double t1,t2,  t3, t5, t8,t10,  t11, t13, t14, t21, t22,t25; 
 
      t1 = 1.0-xi;
      t2 = 1.0-zeta;
      t3 = t1*t2;
      t5 = xi*xi;
      t8 = zeta*zeta;
      t10 = (1.0-t5)*eta*(1.0-t8);
      t11 = 0.8437500002768554*t10;
      t13 = 1.0+xi;
      t14 = t13*t2;
      t21 = 1.0+zeta;
      t22 = t1*t21;
      t25 = t13*t21;
      
      values[0] = -0.125*t3+t11;
      values[1] = -0.125*t14+t11;
      values[2] = 0.125*t3+t11;
      values[3] = 0.125*t14+t11;
      values[4] = -0.125*t22+t11;
      values[5] = -0.125*t25+t11;
      values[6] = 0.125*t22+t11;
      values[7] = 0.125*t25+t11;
      values[8] = -0.8437500004877929*t10; 
 
}

static void C_H_UL1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2,  t3, t5, t7, t10, t11,t13, t14, t17, t18,  t21;
 
      t1 = 1.0-xi;
      t2 = 1.0-eta;
      t3 = t1*t2;
      t5 = xi*xi;
      t7 = eta*eta;
      t10 = (1.0-t5)*(1.0-t7)*zeta;
      t11 = 0.8437500002768554*t10;
      t13 = 1.0+xi;
      t14 = t13*t2;
      t17 = 1.0+eta;
      t18 = t1*t17;
      t21 = t13*t17;
      
      values[0] = -0.125*t3+t11;
      values[1] = -0.125*t14+t11;
      values[2] = -0.125*t18+t11;
      values[3] = -0.125*t21+t11;
      values[4] = 0.125*t3+t11;
      values[5] = 0.125*t14+t11;
      values[6] = 0.125*t18+t11;
      values[7] = 0.125*t21+t11;
      values[8] = -0.8437500004877929*t10;  
  
}

static void C_H_UL1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6;
  
      t1 = eta*eta;
      t3 = zeta*zeta;
      t5 = (1.0-t1)*(1.0-t3);
      t6 = 0.8437500002768554*t5;
      
      values[0] = t6;
      values[1] = t6;
      values[2] = t6;
      values[3] = t6;
      values[4] = t6;
      values[5] = t6;
      values[6] = t6;
      values[7] = t6;
      values[8] = -0.8437500004877929*t5;
}

static void C_H_UL1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6, t7, t8, t9, t10, t11;
  
      t1 = 0.125*zeta;
      t3 = zeta*zeta;
      t5 = xi*eta*(1.0-t3);
      t6 = 0.1687500000553711E1*t5;
      t7 = 0.125-t1-t6;
      t8 = 0.125*zeta;
      t9 = -0.125+t8-t6;
      t10 = 0.125+t8-t6;
      t11 = -0.125-t1-t6;
      
      values[0] = t7;
      values[1] = t9;
      values[2] = t9;
      values[3] = t7;
      values[4] = t10;
      values[5] = t11;
      values[6] = t11;
      values[7] = t10;
      values[8] = 0.1687500000975586E1*t5;
}

static void C_H_UL1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t5, t6, t7, t8, t9,t10, t11 ;

      t1 = 0.125*eta;
      t2 = eta*eta;
      t5 = xi*(1.0-t2)*zeta;
      t6 = 0.1687500000553711E1*t5;
      t7 = 0.125-t1-t6;
      t8 = 0.125*eta;
      t9 = -0.125+t8-t6;
      t10 = 0.125+t8-t6;
      t11 = -0.125-t1-t6;
      
      
      values[0] = t7;
      values[1] = t9;
      values[2] = t10;
      values[3] = t11;
      values[4] = t9;
      values[5] = t7;
      values[6] = t11;
      values[7] = t10;
      values[8] = 0.1687500000975586E1*t5;
      
      
}

static void C_H_UL1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6;
  
      t1 = xi*xi;
      t3 = zeta*zeta;
      t5 = (1.0-t1)*(1.0-t3);
      t6 = 0.8437500002768554*t5;
      
      values[0] = t6;
      values[1] = t6;
      values[2] = t6;
      values[3] = t6;
      values[4] = t6;
      values[5] = t6;
      values[6] = t6;
      values[7] = t6;
      values[8] = -0.8437500004877929*t5;
}

static void C_H_UL1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t5, t6,t7, t8,t9, t10, t11  ;

      t1 = 0.125*xi;
      t2 = xi*xi;
      t5 = (1.0-t2)*eta*zeta;
      t6 = 0.1687500000553711E1*t5;
      t7 = 0.125-t1-t6;
      t8 = 0.125*xi;
      t9 = 0.125+t8-t6;
      t10 = -0.125+t8-t6;
      t11 = -0.125-t1-t6;
      
      values[0] = t7;
      values[1] = t9;
      values[2] = t10;
      values[3] = t11;
      values[4] = t10;
      values[5] = t11;
      values[6] = t7;
      values[7] = t9;
      values[8] = 0.1687500000975586E1*t5;
}

static void C_H_UL1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6;
  
      t1 = xi*xi;
      t3 = eta*eta;
      t5 = (1.0-t1)*(1.0-t3);
      t6 = 0.8437500002768554*t5;
      
      values[0] = t6;
      values[1] = t6;
      values[2] = t6;
      values[3] = t6;
      values[4] = t6;
      values[5] = t6;
      values[6] = t6;
      values[7] = t6;
      values[8] = -0.8437500004877929*t5;   
  
}

TBaseFunct3D *BF_C_H_UL1_3D_Obj = 
new TBaseFunct3D(9, BF_C_H_UL1_3D, BFUnitHexahedron, 
                 C_H_UL1_3D_Funct, C_H_UL1_3D_DeriveXi,
                 C_H_UL1_3D_DeriveEta, C_H_UL1_3D_DeriveZeta,
                 C_H_UL1_3D_DeriveXiXi, C_H_UL1_3D_DeriveXiEta,
                 C_H_UL1_3D_DeriveXiZeta, C_H_UL1_3D_DeriveEtaEta,
                 C_H_UL1_3D_DeriveEtaZeta, C_H_UL1_3D_DeriveZetaZeta,
                 2, 1,
                 0, NULL);
