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
// Q1 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_Q_Q1_2D_Funct(double xi, double eta, double *values)
{
  double h=0.375*(xi*xi-eta*eta);
  values[0]=-h       -0.5*eta+0.25;
  values[1]= h+0.5*xi        +0.25;
  values[2]=-h       +0.5*eta+0.25;
  values[3]= h-0.5*xi        +0.25;
}

// values of the derivatives in xi direction
static void N_Q_Q1_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0]=-0.75*xi;
  values[1]= 0.75*xi+0.5;
  values[2]=-0.75*xi;
  values[3]= 0.75*xi-0.5;
}

// values of the derivatives in eta direction
static void N_Q_Q1_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0]= 0.75*eta-0.5;
  values[1]=-0.75*eta;
  values[2]= 0.75*eta+0.5;
  values[3]=-0.75*eta;
}

// values of derivatives in xi-xi direction
static void N_Q_Q1_2D_DeriveXiXi(double xi, double eta, double *values)
{
  values[0]=-0.75;
  values[1]= 0.75;
  values[2]=-0.75;
  values[3]= 0.75;
}

// values of derivatives in eta-eta direction
static void N_Q_Q1_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  values[0]= 0.75;
  values[1]=-0.75;
  values[2]= 0.75;
  values[3]=-0.75;
}

// values of derivatives in xi-eta direction
static void N_Q_Q1_2D_DeriveXiEta(double xi, double eta, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}

// ***********************************************************************

TBaseFunct2D *BF_N_Q_Q1_2D_Obj = new TBaseFunct2D
        (4, BF_N_Q_Q1_2D, BFUnitSquare, 
         N_Q_Q1_2D_Funct, N_Q_Q1_2D_DeriveXi,
         N_Q_Q1_2D_DeriveEta, N_Q_Q1_2D_DeriveXiXi,
         N_Q_Q1_2D_DeriveXiEta, N_Q_Q1_2D_DeriveEtaEta, 2, 1,
         0, NULL);
