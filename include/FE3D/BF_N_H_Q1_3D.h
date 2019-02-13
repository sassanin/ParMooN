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
// Q1Rot element, nonconforming, 3D
// ***********************************************************************

#define __POINTVALUE__

static void N_H_Q1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t3;

  t1 = xi*xi;
  t2 = eta*eta;
  t3 = zeta*zeta;

#ifdef __POINTVALUE__
  // point value oriented
  values[0] = 1.0/6.0-zeta/2-t1/6-t2/6+t3/3;
  values[1] = 1.0/6.0-eta/2-t1/6+t2/3-t3/6;
  values[2] = 1.0/6.0+xi/2+t1/3-t2/6-t3/6;
  values[3] = 1.0/6.0+eta/2-t1/6+t2/3-t3/6;
  values[4] = 1.0/6.0-xi/2+t1/3-t2/6-t3/6;
  values[5] = 1.0/6.0+zeta/2-t1/6-t2/6+t3/3;
#else
  // mean value oriented
  values[0] = 1.0/6.0-zeta/2-t1/4-t2/4+t3/2;
  values[1] = 1.0/6.0-eta/2-t1/4+t2/2-t3/4;
  values[2] = 1.0/6.0+xi/2+t1/2-t2/4-t3/4;
  values[3] = 1.0/6.0+eta/2-t1/4+t2/2-t3/4;
  values[4] = 1.0/6.0-xi/2+t1/2-t2/4-t3/4;
  values[5] = 1.0/6.0+zeta/2-t1/4-t2/4+t3/2;
#endif
}

static void N_H_Q1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
#ifdef __POINTVALUE__
  // point value oriented
  values[0] = -xi/3;
  values[1] = -xi/3;
  values[2] = 1.0/2.0+2.0/3.0*xi;
  values[3] = -xi/3;
  values[4] = -1.0/2.0+2.0/3.0*xi;
  values[5] = -xi/3;
#else
  // mean value oriented
  values[0] = -xi/2;
  values[1] = -xi/2;
  values[2] = 1.0/2.0+xi;
  values[3] = -xi/2;
  values[4] = -1.0/2.0+xi;
  values[5] = -xi/2;
#endif

}

static void N_H_Q1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
#ifdef __POINTVALUE__
  // point value oriented
  values[0] = -eta/3;
  values[1] = -1.0/2.0+2.0/3.0*eta;
  values[2] = -eta/3;
  values[3] = 1.0/2.0+2.0/3.0*eta;
  values[4] = -eta/3;
  values[5] = -eta/3;
#else
  // mean value oriented
  values[0] = -eta/2;
  values[1] = -1.0/2.0+eta;
  values[2] = -eta/2;
  values[3] = 1.0/2.0+eta;
  values[4] = -eta/2;
  values[5] = -eta/2;
#endif
}

static void N_H_Q1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
#ifdef __POINTVALUE__
  // point value oriented
  values[0] = -1.0/2.0+2.0/3.0*zeta;
  values[1] = -zeta/3;
  values[2] = -zeta/3;
  values[3] = -zeta/3;
  values[4] = -zeta/3;
  values[5] = 1.0/2.0+2.0/3.0*zeta;
#else
  // mean value oriented
  values[0] = -1.0/2.0+zeta;
  values[1] = -zeta/2;
  values[2] = -zeta/2;
  values[3] = -zeta/2;
  values[4] = -zeta/2;
  values[5] = 1.0/2.0+zeta;
#endif
}

static void N_H_Q1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
#ifdef __POINTVALUE__
  values[0] = -1.0/3.0;
  values[1] = -1.0/3.0;
  values[2] = 2.0/3.0;
  values[3] = -1.0/3.0;
  values[4] = 2.0/3.0;
  values[5] = -1.0/3.0;
#else
  values[0] = -0.5;
  values[1] = -0.5;
  values[2] = 1;
  values[3] = -0.5;
  values[4] = 1;
  values[5] = -0.5;
#endif
}

static void N_H_Q1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
}

static void N_H_Q1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
}

static void N_H_Q1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
#ifdef __POINTVALUE__
  values[0] = -1.0/3.0;
  values[1] = 2.0/3.0;
  values[2] = -1.0/3.0;
  values[3] = 2.0/3.0;
  values[4] = -1.0/3.0;
  values[5] = -1.0/3.0;
#else
  values[0] = -0.5;
  values[1] = 1;
  values[2] = -0.5;
  values[3] = 1;
  values[4] = -0.5;
  values[5] = -0.5;
#endif
}

static void N_H_Q1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
}

static void N_H_Q1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
#ifdef __POINTVALUE__
  values[0] = 2.0/3.0;
  values[1] = -1.0/3.0;
  values[2] = -1.0/3.0;
  values[3] = -1.0/3.0;
  values[4] = -1.0/3.0;
  values[5] = 2.0/3.0;
#else
  values[0] = 1;
  values[1] = -0.5;
  values[2] = -0.5;
  values[3] = -0.5;
  values[4] = -0.5;
  values[5] = 1;
#endif
}

TBaseFunct3D *BF_N_H_Q1_3D_Obj = 
new TBaseFunct3D(6, BF_N_H_Q1_3D, BFUnitHexahedron, 
                 N_H_Q1_3D_Funct, N_H_Q1_3D_DeriveXi,
                 N_H_Q1_3D_DeriveEta, N_H_Q1_3D_DeriveZeta,
                 N_H_Q1_3D_DeriveXiXi, N_H_Q1_3D_DeriveXiEta,
                 N_H_Q1_3D_DeriveXiZeta, N_H_Q1_3D_DeriveEtaEta,
                 N_H_Q1_3D_DeriveEtaZeta, N_H_Q1_3D_DeriveZetaZeta,
                 2, 1,
                 0, NULL);
