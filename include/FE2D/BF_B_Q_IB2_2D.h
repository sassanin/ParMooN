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
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************

// base function values
static void B_Q_IB2_2D_Funct(double xi, double eta, double *values)
{
  values[0]=(1.0-xi)*(1.0-eta)*(1.0+xi)*(1.0+eta);
}

// values of the derivatives in xi direction
static void B_Q_IB2_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0]=-2.*xi+2.*eta*eta*xi;
}

// values of the derivatives in eta direction
static void B_Q_IB2_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0]=-2.*eta+2.*xi*xi*eta;
}
// values of the derivatives in xi-xi  direction
static void B_Q_IB2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  values[0]=-2.+2.*eta*eta;
}
// values of the derivatives in xi-eta direction
static void B_Q_IB2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  values[0]=4.*eta*xi;
}
// values of the derivatives in eta-eta direction
static void B_Q_IB2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  values[0]=-2.+2.*xi*xi;
}
// ***********************************************************************

TBaseFunct2D *BF_B_Q_IB2_2D_Obj = new TBaseFunct2D
        (1, BF_B_Q_IB2_2D, BFUnitSquare, 
         B_Q_IB2_2D_Funct, B_Q_IB2_2D_DeriveXi,
         B_Q_IB2_2D_DeriveEta, B_Q_IB2_2D_DeriveXiXi,
         B_Q_IB2_2D_DeriveXiEta, B_Q_IB2_2D_DeriveEtaEta, 2, 2,
         0, NULL);
