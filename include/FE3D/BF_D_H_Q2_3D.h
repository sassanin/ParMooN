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
// Q2 element, discontinuous, 3D
// ***********************************************************************

static void D_H_Q2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0]  = 1;

  values[1]  = xi;
  values[2]  = eta;
  values[3]  = zeta;

  values[4]  = 3*xi*xi-1;
  values[5]  = xi*eta;
  values[6]  = xi*zeta;
  values[7]  = 3*eta*eta-1;
  values[8]  = eta*zeta;
  values[9]  = 3*zeta*zeta-1;

  values[10] = (3*xi*xi-1)*eta;
  values[11] = (3*xi*xi-1)*zeta;
  values[12] = (3*eta*eta-1)*xi;
  values[13] = xi*eta*zeta;
  values[14] = (3*zeta*zeta-1)*xi;
  values[15] = (3*eta*eta-1)*zeta;
  values[16] = (3*zeta*zeta-1)*eta;

  values[17] = (3*xi*xi-1)*(3*eta*eta-1);
  values[18] = (3*xi*xi-1)*eta*zeta;
  values[19] = (3*xi*xi-1)*(3*zeta*zeta-1);
  values[20] = (3*eta*eta-1)*xi*zeta;
  values[21] = (3*zeta*zeta-1)*xi*eta;
  values[22] = (3*eta*eta-1)*(3*zeta*zeta-1);

  values[23] = (3*xi*xi-1)*(3*eta*eta-1)*zeta;
  values[24] = (3*xi*xi-1)*(3*zeta*zeta-1)*eta;
  values[25] = (3*eta*eta-1)*(3*zeta*zeta-1)*xi;

  values[26] = (3*xi*xi-1)*(3*eta*eta-1)*(3*zeta*zeta-1);
}

static void D_H_Q2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 1;
  values[2]  = 0;
  values[3]  = 0;

  values[4]  = 6*xi;
  values[5]  = eta;
  values[6]  = zeta;
  values[7]  = 0;
  values[8]  = 0;
  values[9]  = 0;

  values[10] = 6*xi*eta;
  values[11] = 6*xi*zeta;
  values[12] = (3*eta*eta-1);
  values[13] = eta*zeta;
  values[14] = (3*zeta*zeta-1);
  values[15] = 0;
  values[16] = 0;

  values[17] = 6*xi*(3*eta*eta-1);
  values[18] = 6*xi*eta*zeta;
  values[19] = 6*xi*(3*zeta*zeta-1);
  values[20] = (3*eta*eta-1)*zeta;
  values[21] = (3*zeta*zeta-1)*eta;
  values[22] = 0;

  values[23] = 6*xi*(3*eta*eta-1)*zeta;
  values[24] = 6*xi*(3*zeta*zeta-1)*eta;
  values[25] = (3*eta*eta-1)*(3*zeta*zeta-1);

  values[26] = 6*xi*(3*eta*eta-1)*(3*zeta*zeta-1);
}

static void D_H_Q2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 0;
  values[2]  = 1;
  values[3]  = 0;

  values[4]  = 0;
  values[5]  = xi;
  values[6]  = 0;
  values[7]  = 6*eta;
  values[8]  = zeta;
  values[9]  = 0;

  values[10] = (3*xi*xi-1);
  values[11] = 0;
  values[12] = 6*eta*xi;
  values[13] = xi*zeta;
  values[14] = 0;
  values[15] = 6*eta*zeta;
  values[16] = (3*zeta*zeta-1);

  values[17] = (3*xi*xi-1)*6*eta;
  values[18] = (3*xi*xi-1)*zeta;
  values[19] = 0;
  values[20] = 6*eta*xi*zeta;
  values[21] = (3*zeta*zeta-1)*xi;
  values[22] = 6*eta*(3*zeta*zeta-1);

  values[23] = (3*xi*xi-1)*6*eta*zeta;
  values[24] = (3*xi*xi-1)*(3*zeta*zeta-1);
  values[25] = 6*eta*(3*zeta*zeta-1)*xi;

  values[26] = (3*xi*xi-1)*6*eta*(3*zeta*zeta-1);
}

static void D_H_Q2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 0;
  values[2]  = 0;
  values[3]  = 1;

  values[4]  = 0;
  values[5]  = 0;
  values[6]  = xi;
  values[7]  = 0;
  values[8]  = eta;
  values[9]  = 6*zeta;

  values[10] = 0;
  values[11] = (3*xi*xi-1);
  values[12] = 0;
  values[13] = xi*eta;
  values[14] = 6*zeta*xi;
  values[15] = (3*eta*eta-1);
  values[16] = 6*zeta*eta;

  values[17] = 0;
  values[18] = (3*xi*xi-1)*eta;
  values[19] = (3*xi*xi-1)*6*zeta;
  values[20] = (3*eta*eta-1)*xi;
  values[21] = 6*zeta*xi*eta;
  values[22] = (3*eta*eta-1)*6*zeta;

  values[23] = (3*xi*xi-1)*(3*eta*eta-1);
  values[24] = (3*xi*xi-1)*6*zeta*eta;
  values[25] = (3*eta*eta-1)*6*zeta*xi;

  values[26] = (3*xi*xi-1)*(3*eta*eta-1)*6*zeta;
}

static void D_H_Q2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 0;
  values[2]  = 0;
  values[3]  = 0;

  values[4]  = 6;
  values[5]  = 0;
  values[6]  = 0;
  values[7]  = 0;
  values[8]  = 0;
  values[9]  = 0;

  values[10] = 6*eta;
  values[11] = 6*zeta;
  values[12] = 0;
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;

  values[17] = 6*(3*eta*eta-1);
  values[18] = 6*eta*zeta;
  values[19] = 6*(3*zeta*zeta-1);
  values[20] = 0;
  values[21] = 0;
  values[22] = 0;

  values[23] = 6*(3*eta*eta-1)*zeta;
  values[24] = 6*(3*zeta*zeta-1)*eta;
  values[25] = 0;

  values[26] = 6*(3*eta*eta-1)*(3*zeta*zeta-1);
}

static void D_H_Q2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 0;
  values[2]  = 0;
  values[3]  = 0;

  values[4]  = 0;
  values[5]  = 1;
  values[6]  = 0;
  values[7]  = 0;
  values[8]  = 0;
  values[9]  = 0;

  values[10] = 6*xi;
  values[11] = 0;
  values[12] = 6*eta;
  values[13] = zeta;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;

  values[17] = 6*xi*6*eta;
  values[18] = 6*xi*zeta;
  values[19] = 0;
  values[20] = 6*eta*zeta;
  values[21] = (3*zeta*zeta-1);
  values[22] = 0;

  values[23] = 6*xi*6*eta*zeta;
  values[24] = 6*xi*(3*zeta*zeta-1);
  values[25] = 6*eta*(3*zeta*zeta-1);

  values[26] = 6*xi*6*eta*(3*zeta*zeta-1);
}

static void D_H_Q2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 0;
  values[2]  = 0;
  values[3]  = 0;

  values[4]  = 0;
  values[5]  = 0;
  values[6]  = 1;
  values[7]  = 0;
  values[8]  = 0;
  values[9]  = 0;

  values[10] = 0;
  values[11] = 6*xi;
  values[12] = 0;
  values[13] = eta;
  values[14] = 6*zeta;
  values[15] = 0;
  values[16] = 0;

  values[17] = 0;
  values[18] = 6*xi*eta;
  values[19] = 6*xi*6*zeta;
  values[20] = (3*eta*eta-1);
  values[21] = 6*zeta*eta;
  values[22] = 0;

  values[23] = 6*xi*(3*eta*eta-1);
  values[24] = 6*xi*6*zeta*eta;
  values[25] = (3*eta*eta-1)*6*zeta;

  values[26] = 6*xi*(3*eta*eta-1)*6*zeta;
}

static void D_H_Q2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 0;
  values[2]  = 0;
  values[3]  = 0;

  values[4]  = 0;
  values[5]  = 0;
  values[6]  = 0;
  values[7]  = 6;
  values[8]  = 0;
  values[9]  = 0;

  values[10] = 0;
  values[11] = 0;
  values[12] = 6*xi;
  values[13] = 0;
  values[14] = 0;
  values[15] = 6*zeta;
  values[16] = 0;

  values[17] = (3*xi*xi-1)*6;
  values[18] = 0;
  values[19] = 0;
  values[20] = 6*xi*zeta;
  values[21] = 0;
  values[22] = 6*(3*zeta*zeta-1);

  values[23] = (3*xi*xi-1)*6*zeta;
  values[24] = 0;
  values[25] = 6*(3*zeta*zeta-1)*xi;

  values[26] = (3*xi*xi-1)*6*(3*zeta*zeta-1);
}

static void D_H_Q2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0]  = 0;

  values[1]  = 0;
  values[2]  = 0;
  values[3]  = 0;

  values[4]  = 0;
  values[5]  = 0;
  values[6]  = 0;
  values[7]  = 0;
  values[8]  = 1;
  values[9]  = 0;

  values[10] = 0;
  values[11] = 0;
  values[12] = 0;
  values[13] = xi;
  values[14] = 0;
  values[15] = 6*eta;
  values[16] = 6*zeta;

  values[17] = 0;
  values[18] = (3*xi*xi-1);
  values[19] = 0;
  values[20] = 6*eta*xi;
  values[21] = 6*zeta*xi;
  values[22] = 6*eta*6*zeta;

  values[23] = (3*xi*xi-1)*6*eta;
  values[24] = (3*xi*xi-1)*6*zeta;
  values[25] = 6*eta*6*zeta*xi;

  values[26] = (3*xi*xi-1)*6*eta*6*zeta;
}

static void D_H_Q2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
   values[0]  = 0;

   values[1]  = 0;
   values[2]  = 0;
   values[3]  = 0;

   values[4]  = 0;
   values[5]  = 0;
   values[6]  = 0;
   values[7]  = 0;
   values[8]  = 0;
   values[9]  = 6;

   values[10] = 0;
   values[11] = 0;
   values[12] = 0;
   values[13] = 0;
   values[14] = 6*xi;
   values[15] = 0;
   values[16] = 6*eta;

   values[17] = 0;
   values[18] = 0;
   values[19] = (3*xi*xi-1)*6;
   values[20] = 0;
   values[21] = 6*xi*eta;
   values[22] = (3*eta*eta-1)*6;

   values[23] = 0;
   values[24] = (3*xi*xi-1)*6*eta;
   values[25] = (3*eta*eta-1)*6*xi;

   values[26] = (3*xi*xi-1)*(3*eta*eta-1)*6;
}

TBaseFunct3D *BF_D_H_Q2_3D_Obj =
new TBaseFunct3D(27, BF_D_H_Q2_3D, BFUnitHexahedron,
                 D_H_Q2_3D_Funct, D_H_Q2_3D_DeriveXi,
                 D_H_Q2_3D_DeriveEta, D_H_Q2_3D_DeriveZeta,
                 D_H_Q2_3D_DeriveXiXi, D_H_Q2_3D_DeriveXiEta,
                 D_H_Q2_3D_DeriveXiZeta, D_H_Q2_3D_DeriveEtaEta,
                 D_H_Q2_3D_DeriveEtaZeta, D_H_Q2_3D_DeriveZetaZeta,
                 3, 1,
                 0, NULL);
