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
// P3 element, discontinuous, 3D
// ***********************************************************************

static void D_H_P3_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t3, t5, t6, t8, t9;

  t1 = xi*xi;
  t2 = 3.0*t1-1.0;
  t3 = xi*eta;
  t5 = eta*eta;
  t6 = 3.0*t5-1.0;
  t8 = zeta*zeta;
  t9 = 3.0*t8-1.0;

  values[0] = 1.0;
  values[1] = xi;
  values[2] = eta;
  values[3] = zeta;
  values[4] = t2;
  values[5] = t3;
  values[6] = xi*zeta;
  values[7] = t6;
  values[8] = eta*zeta;
  values[9] = t9;
  values[10] = 5.0*t1*xi-3.0*xi;
  values[11] = t2*eta;
  values[12] = t2*zeta;
  values[13] = xi*t6;
  values[14] = t3*zeta;
  values[15] = xi*t9;
  values[16] = 5.0*t5*eta-3.0*eta;
  values[17] = t6*zeta;
  values[18] = eta*t9;
  values[19] = 5.0*t8*zeta-3.0*zeta;
}

static void D_H_P3_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t5, t8;

  t1 = xi*xi;
  t5 = eta*eta;
  t8 = zeta*zeta;

  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 6.0*xi;
  values[5] = eta;
  values[6] = zeta;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 15.0*t1-3.0;
  values[11] = 6.0*xi*eta;
  values[12] = 6.0*xi*zeta;
  values[13] = 3.0*t5-1.0;
  values[14] = eta*zeta;
  values[15] = 3.0*t8-1.0;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void D_H_P3_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t5, t8;

  t1 = xi*xi;
  t5 = eta*eta;
  t8 = zeta*zeta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 1.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = xi;
  values[6] = 0.0;
  values[7] = 6.0*eta;
  values[8] = zeta;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 3.0*t1-1.0;
  values[12] = 0.0;
  values[13] = 6.0*xi*eta;
  values[14] = xi*zeta;
  values[15] = 0.0;
  values[16] = 15.0*t5-3.0;
  values[17] = 6.0*eta*zeta;
  values[18] = 3.0*t8-1.0;
  values[19] = 0.0;
}

static void D_H_P3_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t5, t8;

  t1 = xi*xi;
  t5 = eta*eta;
  t8 = zeta*zeta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 1.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = xi;
  values[7] = 0.0;
  values[8] = eta;
  values[9] = 6.0*zeta;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 3.0*t1-1.0;
  values[13] = 0.0;
  values[14] = xi*eta;
  values[15] = 6.0*xi*zeta;
  values[16] = 0.0;
  values[17] = 3.0*t5-1.0;
  values[18] = 6.0*eta*zeta;
  values[19] = 15.0*t8-3.0;
}

static void D_H_P3_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 6.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 30.0*xi;
  values[11] = 6.0*eta;
  values[12] = 6.0*zeta;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void D_H_P3_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 1.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 6.0*xi;
  values[12] = 0.0;
  values[13] = 6.0*eta;
  values[14] = zeta;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void D_H_P3_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 1.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 6.0*xi;
  values[13] = 0.0;
  values[14] = eta;
  values[15] = 6.0*zeta;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void D_H_P3_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 6.0;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = 6.0*xi;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = 30.0*eta;
  values[17] = 6.0*zeta;
  values[18] = 0.0;
  values[19] = 0.0;
}

static void D_H_P3_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 1.0;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = 0.0;
  values[14] = xi;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 6.0*eta;
  values[18] = 6.0*zeta;
  values[19] = 0.0;
}

static void D_H_P3_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 6.0;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 0.0;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 6.0*xi;
  values[16] = 0.0;
  values[17] = 0.0;
  values[18] = 6.0*eta;
  values[19] = 30.0*zeta;
}

TBaseFunct3D *BF_D_H_P3_3D_Obj = 
new TBaseFunct3D(20, BF_D_H_P3_3D, BFUnitHexahedron, 
                 D_H_P3_3D_Funct, D_H_P3_3D_DeriveXi,
                 D_H_P3_3D_DeriveEta, D_H_P3_3D_DeriveZeta,
                 D_H_P3_3D_DeriveXiXi, D_H_P3_3D_DeriveXiEta,
                 D_H_P3_3D_DeriveXiZeta, D_H_P3_3D_DeriveEtaEta,
                 D_H_P3_3D_DeriveEtaZeta, D_H_P3_3D_DeriveZetaZeta,
                 3, 3,
                 0, NULL);
