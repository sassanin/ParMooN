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
// P2 element, conforming, 3D
// ***********************************************************************

static void C_T_P2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t3, t4, t5, t6;

  t1 = xi*xi;
  t2 = xi*eta;
  t3 = xi*zeta;
  t4 = eta*eta;
  t5 = eta*zeta;
  t6 = zeta*zeta;

  values[0] = 1.0-3.0*xi-3.0*eta-3.0*zeta+2.0*t1+4.0*t2+4.0*t3
             +2.0*t4+4.0*t5+2.0*t6;
  values[1] = 4.0*xi-4.0*t1-4.0*t2-4.0*t3;
  values[2] = -xi+2.0*t1;
  values[3] = 4.0*eta-4.0*t2-4.0*t4-4.0*t5;
  values[4] = 4.0*t2;
  values[5] = -eta+2.0*t4;
  values[6] = 4.0*zeta-4.0*t3-4.0*t5-4.0*t6;
  values[7] = 4.0*t3;
  values[8] = 4.0*t5;
  values[9] = -zeta+2.0*t6;
}

static void C_T_P2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -3.0+4.0*xi+4.0*eta+4.0*zeta;
  values[1] = 4.0-8.0*xi-4.0*eta-4.0*zeta;
  values[2] = -1.0+4.0*xi;
  values[3] = -4.0*eta;
  values[4] = 4.0*eta;
  values[5] = 0.0;
  values[6] = -4.0*zeta;
  values[7] = 4.0*zeta;
  values[8] = 0.0;
  values[9] = 0.0;
}

static void C_T_P2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -3.0+4.0*xi+4.0*eta+4.0*zeta;
  values[1] = -4.0*xi;
  values[2] = 0.0;
  values[3] = 4.0-4.0*xi-8.0*eta-4.0*zeta;
  values[4] = 4.0*xi;
  values[5] = -1.0+4.0*eta;
  values[6] = -4.0*zeta;
  values[7] = 0.0;
  values[8] = 4.0*zeta;
  values[9] = 0.0;
}

static void C_T_P2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -3.0+4.0*xi+4.0*eta+4.0*zeta;
  values[1] = -4.0*xi;
  values[2] = 0.0;
  values[3] = -4.0*eta;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 4.0-4.0*xi-4.0*eta-8.0*zeta;
  values[7] = 4.0*xi;
  values[8] = 4.0*eta;
  values[9] = -1.0+4.0*zeta;
}

static void C_T_P2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 4.0;
  values[1] = -8.0;
  values[2] = 4.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
}

static void C_T_P2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 4.0;
  values[1] = -4.0;
  values[2] = 0.0;
  values[3] = -4.0;
  values[4] = 4.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
}

static void C_T_P2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 4.0;
  values[1] = -4.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = -4.0;
  values[7] = 4.0;
  values[8] = 0.0;
  values[9] = 0.0;
}

static void C_T_P2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 4.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = -8.0;
  values[4] = 0.0;
  values[5] = 4.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 0.0;
}

static void C_T_P2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 4.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = -4.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = -4.0;
  values[7] = 0.0;
  values[8] = 4.0;
  values[9] = 0.0;
}

static void C_T_P2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 4.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = -8.0;
  values[7] = 0.0;
  values[8] = 0.0;
  values[9] = 4.0;
}

TBaseFunct3D *BF_C_T_P2_3D_Obj = 
new TBaseFunct3D(10, BF_C_T_P2_3D, BFUnitTetrahedron, 
                 C_T_P2_3D_Funct, C_T_P2_3D_DeriveXi,
                 C_T_P2_3D_DeriveEta, C_T_P2_3D_DeriveZeta,
                 C_T_P2_3D_DeriveXiXi, C_T_P2_3D_DeriveXiEta,
                 C_T_P2_3D_DeriveXiZeta, C_T_P2_3D_DeriveEtaEta,
                 C_T_P2_3D_DeriveEtaZeta, C_T_P2_3D_DeriveZetaZeta,
                 2, 2,
                 0, NULL);
