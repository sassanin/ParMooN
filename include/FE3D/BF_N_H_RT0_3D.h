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
// Raviart-Thomas element of zero-th order on hexahedra, 3D
// ***********************************************************************

static void N_H_RT0_3D_Funct(double xi, double eta, double zeta,
                             double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0.125*(xi+1);
  values[3] = 0;
  values[4] = 0.125*(xi-1);
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0.125*(eta-1);
  values[8] = 0;
  values[9] = 0.125*(eta+1);
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0.125*(zeta-1);
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0.125*(zeta+1);
}

static void N_H_RT0_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0.125;
  values[3] = 0;
  values[4] = 0.125;
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0;
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0;
}

static void N_H_RT0_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0.125;
  values[8] = 0;
  values[9] = 0.125;
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0;
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0;
}

static void N_H_RT0_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0.125;
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0.125;
}

static void N_H_RT0_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

TBaseFunct3D *BF_N_H_RT0_3D_Obj = 
new TBaseFunct3D(6, BF_N_H_RT0_3D, BFUnitHexahedron, 
                 N_H_RT0_3D_Funct, N_H_RT0_3D_DeriveXi,
                 N_H_RT0_3D_DeriveEta, N_H_RT0_3D_DeriveZeta,
                 N_H_RT0_3D_DeriveXiXi, N_H_RT0_3D_DeriveXiEta,
                 N_H_RT0_3D_DeriveXiZeta, N_H_RT0_3D_DeriveEtaEta,
                 N_H_RT0_3D_DeriveEtaZeta, N_H_RT0_3D_DeriveZetaZeta,
                 1, 1,
                 0, NULL, 3);
