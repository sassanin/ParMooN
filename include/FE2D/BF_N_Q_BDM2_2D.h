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
// P2 BDM vector element, nonconforming , 2D
// History:  23.09.2013 implementation (Markus Wolff)
// ***********************************************************************

// base function values
// vector function, orthonormal to edges, 
// functions 1,2, and 3 are orthogonal to edge 1
// functions 4,5, and 6 are orthogonal to edge 2
// functions 7,8, and 9 are orthogonal to edge 3
// functions 10,11, and 12 are orthogonal to edge 4
// functions 13-14 correspond to inner degrees of freedom

// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in
// NF_N_Q_BDM2_2D.h

// coefficient matrix for different choices of the degrees of freedom
// There seems to be no difference in using integral or point evaluation for 
// inner dofs. 
// The dofs on the edges determine the boundary conditions as well. Using 
// Tschebyscheff-points might improve the boundary approximation.

/*
// using equidistant points on edges and point evaluations for inner dofs
static double N_Q_BDM2_2D_CM[196] = {
0,0,0,0,0,0,0,-0,0,0,0,0,1,0,
0,0,0,0,0,-0,0,0,0,0,0,0,0,1,
0.16666667,-0.33333333,0.16666667,0,0.25,0,0.16666667,-0.33333333,0.16666667,0,0.25,0,0,0,
0.25,0,-0.25,0,0,0,0.25,0,-0.25,0,0,0,0,0,
0,0,0,0,0.25,0,0,0,0,0,-0.25,0,-1,0,
-0.5,1,-0.5,0,0,0,0.5,-1,0.5,0,0,0,0,0,
0,0,0,-0.25,0,0.25,0,-0,0,-0.25,0,0.25,0,0,
0,0.25,0,0.16666667,-0.33333333,0.16666667,0,0.25,0,0.16666667,-0.33333333,0.16666667,0,0,
0,0,0,0.5,-1,0.5,0,0,0,-0.5,1,-0.5,0,0,
0,-0.25,0,0,0,0,0,0.25,0,0,0,0,0,-1,
0,0,0,-0.25,0,0.25,0,0,0,0.25,0,-0.25,0,0,
-0.25,0,0.25,0,0,0,0.25,0,-0.25,0,0,0,0,0,
-0.16666667,0.33333333,-0.16666667,0,-0,0,-0.16666667,0.33333333,-0.16666667,0,0,0,0,0,
0,0,0,0.16666667,-0.33333333,0.16666667,0,0,0,0.16666667,-0.33333333,0.16666667,0,0
};
*/
/*
// using equidistant points on edges and integral evaluation for inner dofs
static double N_Q_BDM2_2D_CM[196] = {
-0,0,-0,-0.25,0.375,-0.25,-0,0,-0,0.25,-0.375,0.25,0.375,-0,
0.25,-0.375,0.25,0,0,0,-0.25,0.375,-0.25,0,0,0,0,0.375,
0.16666667,-0.33333333,0.16666667,0,0.25,0,0.16666667,-0.33333333,0.16666667,0,0.25,0,0,0,
0.25,-0,-0.25,0,0,0,0.25,-0,-0.25,0,0,0,0,-0,
0,-0,0,0.25,-0.125,0.25,0,0,0,-0.25,0.125,-0.25,-0.375,0,
-0.5,1,-0.5,0,-0,0,0.5,-1,0.5,0,-0,0,0,0,
0,-0,0,-0.25,0,0.25,0,0,0,-0.25,0,0.25,0,0,
0,0.25,0,0.16666667,-0.33333333,0.16666667,0,0.25,0,0.16666667,-0.33333333,0.16666667,0,0,
0,-0,0,0.5,-1,0.5,0,0,0,-0.5,1,-0.5,0,0,
-0.25,0.125,-0.25,0,0,0,0.25,-0.125,0.25,0,0,0,0,-0.375,
0,0,0,-0.25,0,0.25,-0,0,0,0.25,0,-0.25,0,0,
-0.25,-0,0.25,0,0,0,0.25,-0,-0.25,0,0,0,0,-0,
-0.16666667,0.33333333,-0.16666667,0,0,0,-0.16666667,0.33333333,-0.16666667,0,0,0,0,-0,
0,0,0,0.16666667,-0.33333333,0.16666667,0,0,0,0.16666667,-0.33333333,0.16666667,0,0
};
*/
/*
// using Gauss-Points on edges and point evaluation for inner dofs
static double N_Q_BDM2_2D_CM[196] = {
0,0,0,0,0,0,0,-0,0,0,0,0,1,0,
0,0,0,0,0,-0,0,0,0,0,0,0,0,1,
0.069444444,-0.13888889,0.069444444,0,0.25,0,0.069444444,-0.13888889,0.069444444,0,0.25,0,0,0,
0.16137431,0,-0.16137431,0,0,0,0.16137431,-0,-0.16137431,0,0,0,0,0,
0,0,0,0,0.25,0,0,0,0,0,-0.25,0,-1,0,
-0.20833333,0.41666667,-0.20833333,0,0,0,0.20833333,-0.41666667,0.20833333,0,0,0,0,0,
0,0,0,-0.16137431,0,0.16137431,0,-0,0,-0.16137431,0,0.16137431,0,0,
0,0.25,0,0.069444444,-0.13888889,0.069444444,0,0.25,0,0.069444444,-0.13888889,0.069444444,0,0,
0,0,0,0.20833333,-0.41666667,0.20833333,0,0,0,-0.20833333,0.41666667,-0.20833333,0,0,
0,-0.25,0,0,0,0,0,0.25,0,0,0,0,0,-1,
0,0,0,-0.16137431,0,0.16137431,0,0,0,0.16137431,0,-0.16137431,0,0,
-0.16137431,-0,0.16137431,0,0,0,0.16137431,-0,-0.16137431,0,0,0,0,0,
-0.069444444,0.13888889,-0.069444444,0,-0,0,-0.069444444,0.13888889,-0.069444444,0,0,0,0,0,
0,0,0,0.069444444,-0.13888889,0.069444444,0,0,0,0.069444444,-0.13888889,0.069444444,0,0
};
*/
/*
// using Gauss-Points on edges and integral evaluation for inner dofs
static double N_Q_BDM2_2D_CM[196] = {
-0,0,-0,-0.10416667,0.083333333,-0.10416667,-0,0,-0,0.10416667,-0.083333333,0.10416667,0.375,0,
0.10416667,-0.083333333,0.10416667,-0,0,-0,-0.10416667,0.083333333,-0.10416667,-0,0,-0,0,0.375,
0.069444444,-0.13888889,0.069444444,0,0.25,0,0.069444444,-0.13888889,0.069444444,0,0.25,0,0,0,
0.16137431,-0,-0.16137431,0,0,0,0.16137431,-0,-0.16137431,0,0,0,0,-0,
0,0,0,0.10416667,0.16666667,0.10416667,0,0,0,-0.10416667,-0.16666667,-0.10416667,-0.375,0,
-0.20833333,0.41666667,-0.20833333,0,0,0,0.20833333,-0.41666667,0.20833333,0,0,0,0,0,
0,-0,0,-0.16137431,0,0.16137431,0,-0,0,-0.16137431,0,0.16137431,0,0,
0,0.25,0,0.069444444,-0.13888889,0.069444444,0,0.25,0,0.069444444,-0.13888889,0.069444444,0,0,
0,-0,0,0.20833333,-0.41666667,0.20833333,0,0,0,-0.20833333,0.41666667,-0.20833333,0,0,
-0.10416667,-0.16666667,-0.10416667,0,0,0,0.10416667,0.16666667,0.10416667,0,0,0,0,-0.375,
0,-0,0,-0.16137431,0,0.16137431,0,0,0,0.16137431,0,-0.16137431,0,0,
-0.16137431,0,0.16137431,0,0,0,0.16137431,0,-0.16137431,0,0,0,0,0,
-0.069444444,0.13888889,-0.069444444,0,0,0,-0.069444444,0.13888889,-0.069444444,0,0,0,0,0,
0,0,0,0.069444444,-0.13888889,0.069444444,0,0,0,0.069444444,-0.13888889,0.069444444,0,0
};
*/
/*
// using Tschebyscheff-points on edges and point evaluations for inner dofs
static double N_Q_BDM2_2D_CM[196] = {
0,0,0,0,0,0,0,-0,0,0,0,0,1,0,
0,0,0,0,0,-0,0,0,0,0,0,0,0,1,
0.055555556,-0.11111111,0.055555556,0,0.25,0,0.055555556,-0.11111111,0.055555556,0,0.25,0,0,0,
0.14433757,-0,-0.14433757,0,0,0,0.14433757,0,-0.14433757,0,0,0,0,0,
0,0,0,0,0.25,0,0,0,0,0,-0.25,0,-1,0,
-0.16666667,0.33333333,-0.16666667,0,0,0,0.16666667,-0.33333333,0.16666667,0,0,0,0,0,
0,0,0,-0.14433757,0,0.14433757,0,-0,0,-0.14433757,0,0.14433757,0,0,
0,0.25,0,0.055555556,-0.11111111,0.055555556,0,0.25,0,0.055555556,-0.11111111,0.055555556,0,0,
0,0,0,0.16666667,-0.33333333,0.16666667,0,0,0,-0.16666667,0.33333333,-0.16666667,0,0,
0,-0.25,0,0,0,0,0,0.25,0,0,0,0,0,-1,
0,0,0,-0.14433757,0,0.14433757,0,0,0,0.14433757,0,-0.14433757,0,0,
-0.14433757,0,0.14433757,0,0,0,0.14433757,0,-0.14433757,0,0,0,0,0,
-0.055555556,0.11111111,-0.055555556,0,-0,0,-0.055555556,0.11111111,-0.055555556,0,0,0,0,0,
0,0,0,0.055555556,-0.11111111,0.055555556,0,0,0,0.055555556,-0.11111111,0.055555556,0,0
};
*/

// using Tschebyscheff-points on edges and integral evaluations for inner dofs
static double N_Q_BDM2_2D_CM[196] = {
-0,0,-0,-0.083333333,0.041666667,-0.083333333,-0,0,-0,0.083333333,-0.041666667,0.083333333,0.375,-0,
0.083333333,-0.041666667,0.083333333,-0,0,-0,-0.083333333,0.041666667,-0.083333333,-0,0,-0,0,0.375,
0.055555556,-0.11111111,0.055555556,0,0.25,0,0.055555556,-0.11111111,0.055555556,0,0.25,0,0,-0,
0.14433757,0,-0.14433757,0,0,0,0.14433757,-0,-0.14433757,0,0,0,0,-0,
0,0,0,0.083333333,0.20833333,0.083333333,0,0,0,-0.083333333,-0.20833333,-0.083333333,-0.375,0,
-0.16666667,0.33333333,-0.16666667,0,0,0,0.16666667,-0.33333333,0.16666667,0,0,0,0,-0,
0,-0,0,-0.14433757,0,0.14433757,0,0,0,-0.14433757,0,0.14433757,0,0,
0,0.25,0,0.055555556,-0.11111111,0.055555556,0,0.25,0,0.055555556,-0.11111111,0.055555556,0,0,
0,-0,0,0.16666667,-0.33333333,0.16666667,0,0,0,-0.16666667,0.33333333,-0.16666667,0,0,
-0.083333333,-0.20833333,-0.083333333,0,0,0,0.083333333,0.20833333,0.083333333,0,0,0,0,-0.375,
0,-0,0,-0.14433757,0,0.14433757,0,0,0,0.14433757,0,-0.14433757,0,0,
-0.14433757,0,0.14433757,0,0,0,0.14433757,-0,-0.14433757,0,0,0,0,0,
-0.055555556,0.11111111,-0.055555556,0,0,0,-0.055555556,0.11111111,-0.055555556,0,0,0,0,0,
0,0,0,0.055555556,-0.11111111,0.055555556,0,0,0,0.055555556,-0.11111111,0.055555556,0,0
};



static void N_Q_BDM2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={1,0,xi, 0,xi*xi,    0,eta,  0,eta*eta,      0,xi*eta,     0,    xi*xi*xi,3*xi*eta*eta};
  double mon_y[14]={0,1, 0,xi,    0,xi*xi,  0,eta,      0,eta*eta,     0,xi*eta,-3*xi*xi*eta,-eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_Q_BDM2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,1,0,2*xi,   0,0,0,0,0,eta,  0,   3*xi*xi,3*eta*eta};
  double mon_y[14]={0,0,0,1,   0,2*xi,0,0,0,0,  0,eta,-6*xi*eta,0      };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_Q_BDM2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,0,0,0,0,1,0,2*eta,    0,xi, 0,      0,6*xi*eta};
  double mon_y[14]={0,0,0,0,0,0,0,1,    0,2*eta, 0,xi,-3*xi*xi,-3*eta*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_Q_BDM2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,0,0,2,0,0,0,0,0,0,0, 6*xi,0};
  double mon_y[14]={0,0,0,0,0,2,0,0,0,0,0,0-6*eta,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_Q_BDM2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,0,0,0,0,0,0,2,0,0,0,0,6*xi};
  double mon_y[14]={0,0,0,0,0,0,0,0,0,2,0,0,0,-6*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_Q_BDM2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,0,0,0,0,0,0,0,0,1,0,    0,6*eta};
  double mon_y[14]={0,0,0,0,0,0,0,0,0,0,0,1,-6*xi,    0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}


// ***********************************************************************

TBaseFunct2D *BF_N_Q_BDM2_2D_Obj = new TBaseFunct2D
        (14, BF_N_Q_BDM2_2D, BFUnitSquare,
         N_Q_BDM2_2D_Funct, N_Q_BDM2_2D_DeriveXi,
         N_Q_BDM2_2D_DeriveEta, N_Q_BDM2_2D_DeriveXiXi,
         N_Q_BDM2_2D_DeriveXiEta, N_Q_BDM2_2D_DeriveEtaEta, 3, 3,
         0, NULL, 2);
