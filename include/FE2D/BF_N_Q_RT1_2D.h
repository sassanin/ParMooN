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
// Q0 Raviart-Thomas vector element, nonconforming , 2D
// History:  17.11.2011 implementation (Ulrich)
// ***********************************************************************

// base function values
// vector function, orthogonal to edges, 
// functions 1 and 2 are orthogonal to edge 1
// functions 3 and 4 are orthogonal to edge 2
// functions 5 and 6 are orthogonal to edge 3
// functions 7 and 8 are orthogonal to edge 4
// functions 9-12 correspond to inner degrees of freedom

// coefficient matrix for different choices of the degrees of freedom
// There seems to be no difference in using integral or point evaluation for 
// inner dofs. 
// The dofs on the edges determine the boundary conditions as well. Using 
// Tschebyscheff-points might improve the boundary approximation.

// using equidistant points on edges and integral evaluation for inner dofs
/*static double N_Q_RT1_2D_CM[144] = { 
  0,0,-0.0625,-0.0625,0,0,0.0625,0.0625,0.375,0,0,0,
  0.0625,0.0625,0,0,-0.0625,-0.0625,0,0,0,0.375,0,0,
  0,0,0.125,0.125,-0,0,0.125,0.125,0,0,0,0,
 -0.1875,0.1875,0,0,-0.1875,0.1875,0,0,0,0,1.125,0,
  0,0,0.1875,-0.1875,-0,0,0.1875,-0.1875,0,0,0,1.125,
  0.125,0.125,0,-0,0.125,0.125,0,0,0,0,0,0,
  0,0,-0.375,0.375,-0,0,0.375,-0.375,0,0,0,0,
 -0.375,0.375,0,0,0.375,-0.375,0,0,0,0,0,0,
  0,0,0.1875,0.1875,-0,0,-0.1875,-0.1875,-0.375,0,0,0,
 -0.1875,-0.1875,0,0,0.1875,0.1875,0,0,0,-0.375,0,0,
  0,0,-0.5625,0.5625,-0,0,-0.5625,0.5625,0,0,0,-1.125,
  0.5625,-0.5625,0,0,0.5625,-0.5625,0,0,0,0,-1.125,0
};*/
/*
// using equidistant points on edges and point evaluation for inner dofs
static double N_Q_RT1_2D_CM[144] = { 
0,0,0,0,0,0,0,0,0.5,0,0.5,0,
0,0,0,0,0,0,0,0,0,0.5,0,0.5,
0,0,0.125,0.125,0,0,0.125,0.125,0,0,0,0,
0,0,0,0,0,0,0,0,0,1.5,0,-1.5,
0,0,0,0,0,0,0,0,-1.5,0,1.5,0,
0.125,0.125,0,0,0.125,0.125,0,0,0,0,0,0,
0,0,-0.375,0.375,0,0,0.375,-0.375,0,0,0,0,
-0.375,0.375,0,0,0.375,-0.375,0,0,0,0,0,0,
0,0,0.125,0.125,0,0,-0.125,-0.125,-0.5,0,-0.5,0,
-0.125,-0.125,0,0,0.125,0.125,0,0,0,-0.5,-0,-0.5,
0,0,-0.375,0.375,0,0,-0.375,0.375,1.5,0,-1.5,0,
0.375,-0.375,0,0,0.375,-0.375,0,0,0,-1.5,0,1.5
};
*/
// using Gauss-Points on edges and integral evaluation for inner dofs
/*static double N_Q_RT1_2D_CM[144] = { 
0,0,-0.0625,-0.0625,0,0,0.0625,0.0625,0.375,0,0,0,
0.0625,0.0625,0,0,-0.0625,-0.0625,0,0,0,0.375,0,0,
0,0,0.125,0.125,0,0,0.125,0.125,0,0,0,0,
-0.10825318,0.10825318,0,0,-0.10825318,0.10825318,0,0,0,0,1.125,0,
0,0,0.10825318,-0.10825318,0,0,0.10825318,-0.10825318,0,0,0,1.125,
0.125,0.125,0,-0,0.125,0.125,0,0,0,0,0,0,
0,0,-0.21650635,0.21650635,0,0,0.21650635,-0.21650635,0,0,0,0,
-0.21650635,0.21650635,0,0,0.21650635,-0.21650635,0,0,0,0,0,0,
0,0,0.1875,0.1875,0,0,-0.1875,-0.1875,-0.375,0,0,0,
-0.1875,-0.1875,0,0,0.1875,0.1875,0,0,0,-0.375,0,0,
0,0,-0.32475953,0.32475953,0,0,-0.32475953,0.32475953,0,0,0,-1.125,
0.32475953,-0.32475953,0,0,0.32475953,-0.32475953,0,0,0,0,-1.125,0
};*/
/*
// using Gauss-Points on edges and point evaluation for inner dofs
static double N_Q_RT1_2D_CM[144] = { 
0,0,0,0,-0,0,0,0,0.5,0,0.5,-0,
0,0,0,0,0,0,0,0,0,0.5,0,0.5,
0,0,0.125,0.125,-0,0,0.125,0.125,0,0,0,-0,
0,0,0,0,0,0,0,0,0,0.8660254038,0,-0.8660254038,
0,0,0,0,0,0,0,0,-0.8660254038,0,0.8660254038,-0,
0.125,0.125,0,0,0.125,0.125,0,0,0,0,0,-0,
0,0,-0.2165063509,0.2165063509,0,0,0.2165063509,-0.2165063509,0,0,0,-0,
-0.2165063509,0.2165063509,0,0,0.2165063509,-0.2165063509,0,0,0,0,0,-0,
0,0,0.125,0.125,0,0,-0.125,-0.125,-0.5,0,-0.5,-0,
-0.125,-0.125,0,0,0.125,0.125,0,0,0,-0.5,-0,-0.5,
0,0,-0.2165063509,0.2165063509,0,0,-0.2165063509,0.2165063509,0.8660254038,0,-0.8660254038,0,
0.2165063509,-0.2165063509,0,0,0.2165063509,-0.2165063509,0,0,0,-0.8660254038,0,0.8660254038
};*/
/*
// using Tschebyscheff-points on edges and integral evaluation for inner dofs
static double N_Q_RT1_2D_CM[144] = { 
0,0,-0.0625,-0.0625,0,0,0.0625,0.0625,0.375,0,0,0,
0.0625,0.0625,0,0,-0.0625,-0.0625,0,0,0,0.375,-0,0,
0,0,0.125,0.125,0,0,0.125,0.125,0,0,0,0,
-0.08838835,0.08838835,0,0,-0.08838835,0.08838835,0,0,0,0,1.125,0,
0,0,0.08838835,-0.08838835,0,0,0.08838835,-0.08838835,0,0,0,1.125,
0.125,0.125,0,0,0.125,0.125,0,0,0,0,0,0,
0,0,-0.1767767,0.1767767,0,0,0.1767767,-0.1767767,0,0,0,0,
-0.1767767,0.1767767,0,0,0.1767767,-0.1767767,0,0,0,0,0,0,
0,-0,0.1875,0.1875,0,0,-0.1875,-0.1875,-0.375,0,0,0,
-0.1875,-0.1875,0,0,0.1875,0.1875,0,0,0,-0.375,0,0,
0,0,-0.26516504,0.26516504,0,0,-0.26516504,0.26516504,0,0,0,-1.125,
0.26516504,-0.26516504,0,0,0.26516504,-0.26516504,0,0,0,0,-1.125,0
};*/

// using Tschebyscheff-points on edges and point evaluation for inner dofs
static double N_Q_RT1_2D_CM[144] = { 
0,0,0,0,-0,0,0,0,0.5,0,0.5,0,
0,0,0,0,0,0,0,0,0,0.5,0,0.5,
0,0,0.125,0.125,-0,0,0.125,0.125,0,0,0,0,
0,0,0,0,0,0,0,0,0,0.7071067812,0,-0.7071067812,
0,0,0,0,0,0,0,0,-0.7071067812,0,0.7071067812,0,
0.125,0.125,0,0,0.125,0.125,0,0,0,0,0,0,
0,0,-0.1767766953,0.1767766953,0,0,0.1767766953,-0.1767766953,0,0,0,0,
-0.1767766953,0.1767766953,0,0,0.1767766953,-0.1767766953,0,0,0,0,0,0,
0,0,0.125,0.125,0,0,-0.125,-0.125,-0.5,0,-0.5,0,
-0.125,-0.125,0,0,0.125,0.125,0,0,0,-0.5,0,-0.5,
0,0,-0.1767766953,0.1767766953,0,0,-0.1767766953,0.1767766953,0.7071067812,0,-0.7071067812,0,
0.1767766953,-0.1767766953,0,0,0.1767766953,-0.1767766953,0,0,0,-0.7071067812,0,0.7071067812
};


static void N_Q_RT1_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={1,0,xi,0 ,eta,0  ,xi*eta,0     ,xi*xi,0      ,xi*xi*eta,0};
  double mon_y[]={0,1,0 ,xi,0  ,eta,0     ,xi*eta,0    ,eta*eta,0,xi*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_Q_RT1_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // xi-derivative of monomials x-component and y-component
  double mon_x[]={0,0,1,0,0,0,eta,0  ,2*xi,0,2*xi*eta,0};
  double mon_y[]={0,0,0,1,0,0,0  ,eta,0   ,0,0       ,eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_Q_RT1_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // eta-derivative of monomials x-component and y-component
  double mon_x[]={0,0,0,0,1,0,xi,0 ,0,0    ,xi*xi,0};
  double mon_y[]={0,0,0,0,0,1,0 ,xi,0,2*eta,0    ,2*xi*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_Q_RT1_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // xi-derivative of monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,0,0,2,0,2*eta,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0    ,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_Q_RT1_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // eta-derivative of monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,2,0,2*xi};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_Q_RT1_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // eta-derivative of monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,2*xi,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0   ,2*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// ***********************************************************************

TBaseFunct2D *BF_N_Q_RT1_2D_Obj = new TBaseFunct2D
        (12, BF_N_Q_RT1_2D, BFUnitSquare, 
         N_Q_RT1_2D_Funct, N_Q_RT1_2D_DeriveXi,
         N_Q_RT1_2D_DeriveEta, N_Q_RT1_2D_DeriveXiXi,
         N_Q_RT1_2D_DeriveXiEta, N_Q_RT1_2D_DeriveEtaEta, 3, 2,
         0, NULL, 2);
