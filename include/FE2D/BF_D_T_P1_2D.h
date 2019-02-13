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
// P1 element, discontinous, 2D, triangle
// ***********************************************************************

// base function values
static void D_T_P1_2D_Funct(double xi, double eta, double *values)
{
  values[0] = 1;
  values[1] = 2*xi+eta-1;
  values[2] = xi+2*eta-1;
}

// values of the derivatives in xi direction
static void D_T_P1_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 2;
  values[2] = 1;
}

// values of the derivatives in eta direction
static void D_T_P1_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 2;
}

// values of the derivatives in xi-xi direction
static void D_T_P1_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in eta-eta direction
static void D_T_P1_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in xi-eta direction
static void D_T_P1_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}
// ***********************************************************************

TBaseFunct2D *BF_D_T_P1_2D_Obj = new TBaseFunct2D
        (3, BF_D_T_P1_2D, BFUnitTriangle, 
         D_T_P1_2D_Funct, D_T_P1_2D_DeriveXi,
         D_T_P1_2D_DeriveEta, D_T_P1_2D_DeriveXiXi,
         D_T_P1_2D_DeriveXiEta, D_T_P1_2D_DeriveEtaEta, 1, 1,
         0, NULL);
