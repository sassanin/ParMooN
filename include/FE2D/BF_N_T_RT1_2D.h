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
// P1 Raviart-Thomas vector element, nonconforming , 2D
// History:  17.10.2011 implementation (Ulrich)
// ***********************************************************************

// base function values
// vector function, orthonormal to edges, 
// functions 1 and 2 are orthogonal to edge 1
// functions 3 and 4 are orthogonal to edge 2
// functions 5 and 6 are orthogonal to edge 3
// functions 7 and 8 correspond to inner degrees of freedom

// coefficient matrix
// using equidistant points on edges and integral evaluation for inner dofs
/*static double N_T_RT1_2D_CM[64] = { 
  0, 0, 0, 0, 1,-2, 0,  0,
 -2, 1, 0, 0, 0, 0, 0,  0,
  3,-2,-2,-1,-1, 6, 16, 8,
  3,-3, 0, 0, 0, 0, 0,  0,
  0, 0, 0, 0,-3, 3, 0,  0,
  6,-1,-1,-2,-2, 3, 8,  16,
 -4, 4, 4, 0, 0,-4,-16,-8,
 -4, 0,-0, 4, 4,-4,-8, -16
};*/
/*
// using Gauss-Points on edges and integral evaluation for inner dofs
static double N_T_RT1_2D_CM[64] = { // using Gauss-Points
0,0,0,0,0.3660254,-1.3660254,0,0,
-1.3660254,0.3660254,0,0,0,0,0,0,
1.94337567,-0.94337567,-1.78867513,-1.21132487,0.47927406,4.52072594,16,8,
1.73205081,-1.73205081,0,0,0,0,0,0,
0,0,0,0,-1.73205081,1.73205081,0,0,
4.52072594,0.47927406,-1.21132487,-1.78867513,-0.94337567,1.94337567,8,16,
-2.30940108,2.30940108,3.15470054,0.84529946,-0.84529946,-3.15470054,-16,-8,
-3.15470054,-0.84529946,0.84529946,3.15470054,2.30940108,-2.30940108,-8,-16
};
*/
/*
// using Tschebyscheff-points on edges and integral evaluation for inner dofs
static double N_T_RT1_2D_CM[64] = { 
0,0,0,0,0.20710678,-1.20710678,0,0,
-1.20710678,0.20710678,-0,-0,-0,0,0,0,
1.6785113,-0.6785113,-1.73570226,-1.26429774,0.85008418,4.14991582,16,8,
1.41421356,-1.41421356,0,0,0,-0,-0,-0,
0,0,0,0,-1.41421356,1.41421356,0,0,
4.14991582,0.85008418,-1.26429774,-1.73570226,-0.6785113,1.6785113,8,16,
-1.88561808,1.88561808,2.94280904,1.05719096,-1.05719096,-2.94280904,-16,-8,
-2.94280904,-1.05719096,1.05719096,2.94280904,1.88561808,-1.88561808,-8,-16
};
*/
// using Tschebyscheff-points on edges and point evaluation for inner dofs
static double N_T_RT1_2D_CM[64] = {
0,0,0,0,0.2071067812,-1.2071067812,0,0,
-1.2071067812,0.2071067812,0,0,0,0,0,0,
1.2071067812,-0.2071067812,-1,-1,0.5857864376,3.4142135624,6,3,
1.4142135624,-1.4142135624,0,0,0,0,0,0,
0,0,0,0,-1.4142135624,1.4142135624,0,0,
3.4142135624,0.5857864376,-1,-1,-0.2071067812,1.2071067812,3,6,
-1.4142135624,1.4142135624,2.2071067812,0.7928932188,-0.7928932188,-2.2071067812,-6,-3,
-2.2071067812,-0.7928932188,0.7928932188,2.2071067812,1.4142135624,-1.4142135624,-3,-6
};
static void N_T_RT1_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={1,0,xi,0 ,eta,0  ,xi*xi ,xi*eta };
  double mon_y[]={0,1,0 ,xi,0  ,eta,xi*eta,eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_RT1_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,1,0,0,0,2*xi,eta };
  double mon_y[]={0,0,0,1,0,0,eta ,0};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_RT1_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,1,0,0 ,xi };
  double mon_y[]={0,0,0,0,0,1,xi,2*eta};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_RT1_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,2,0};
  double mon_y[]={0,0,0,0,0,0,0,0};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_RT1_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,2};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_RT1_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,0,1};
  double mon_y[]={0,0,0,0,0,0,1,0};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// ***********************************************************************

TBaseFunct2D *BF_N_T_RT1_2D_Obj = new TBaseFunct2D
        (8, BF_N_T_RT1_2D, BFUnitTriangle, 
         N_T_RT1_2D_Funct, N_T_RT1_2D_DeriveXi,
         N_T_RT1_2D_DeriveEta, N_T_RT1_2D_DeriveXiXi,
         N_T_RT1_2D_DeriveXiEta, N_T_RT1_2D_DeriveEtaEta, 2, 1,
         0, NULL, 2);
