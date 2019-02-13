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
// P1 element, discontinuous, 3D Tetrahedra
// 
// Author:     Sashikumaar Ganesan
//
// ***********************************************************************

static void D_T_P1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
   values[0] =  6.0;
   values[1] = -1.0+2.0*xi+eta+zeta;
   values[2] = -1.0+xi+2.0*eta+zeta;
   values[3] = -1.0+xi+eta+2.0*zeta;
}

static void D_T_P1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 2;
  values[2] = 1;
  values[3] = 1;
}

static void D_T_P1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 2;
  values[3] = 1;
}

static void D_T_P1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 1;
  values[3] = 2;
}

static void D_T_P1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

TBaseFunct3D *BF_D_T_P1_3D_Obj =
new TBaseFunct3D(4, BF_D_T_P1_3D, BFUnitTetrahedron,
                 D_T_P1_3D_Funct, D_T_P1_3D_DeriveXi,
                 D_T_P1_3D_DeriveEta, D_T_P1_3D_DeriveZeta,
                 D_T_P1_3D_DeriveXiXi, D_T_P1_3D_DeriveXiEta,
                 D_T_P1_3D_DeriveXiZeta, D_T_P1_3D_DeriveEtaEta,
                 D_T_P1_3D_DeriveEtaZeta, D_T_P1_3D_DeriveZetaZeta,
                 1, 1,
		 0, NULL);
