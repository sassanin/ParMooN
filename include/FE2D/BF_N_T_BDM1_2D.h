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
// P1 BDM vector element, nonconforming , 2D
// History:  12.09.2013 implementation (Markus Wolff)
// ***********************************************************************

// base function values
// vector function, orthonormal to edges, 
// functions 1 and 2 are orthogonal to edge 1
// functions 3 and 4 are orthogonal to edge 2
// functions 5 and 6 are orthogonal to edge 3

// coefficient matrix
// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change the evaluation points in NF_N_T_BDM1_2D.h
// using equidistant points on edges
/*
static double N_T_BDM1_2D_CM[36] = {
-0,0,0,0,1,-2,
-2,1,0,0,0,0,
-1,2,2,-1,-1,2,
3,-3,0,0,0,0,
0,0,0,0,-3,3,
2,-1,-1,2,2,-1
};*/

// using Gauss-Points on edges
/*
static double N_T_BDM1_2D_CM[36] = {
-0,0,0,0,0.3660254,-1.3660254,
-1.3660254,0.3660254,0,0,0,0,
-0.3660254,1.3660254,1.3660254,-0.3660254,-0.3660254,1.3660254,
1.7320508,-1.7320508,0,0,0,0,
0,0,0,0,-1.7320508,1.7320508,
1.3660254,-0.3660254,-0.3660254,1.3660254,1.3660254,-0.3660254
};*/

// using Tschebyscheff-points on edges
static double N_T_BDM1_2D_CM[36] = {
-0,0,0,0,0.20710678,-1.2071068,
-1.2071068,0.20710678,-0,0,0,-0,
-0.20710678,1.2071068,1.2071068,-0.20710678,-0.20710678,1.2071068,
1.4142136,-1.4142136,0,-0,-0,0,
0,0,0,0,-1.4142136,1.4142136,
1.2071068,-0.20710678,-0.20710678,1.2071068,1.2071068,-0.20710678
};

static void N_T_BDM1_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 6; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={1,0,xi,0 ,eta,0 };
  double mon_y[]={0,1,0 ,xi,0  ,eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_BDM1_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 6; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,1,0,0,0};
  double mon_y[]={0,0,0,1,0,0};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_BDM1_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 6; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,1,0};
  double mon_y[]={0,0,0,0,0,1};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_BDM1_2D_DeriveXiXi(double xi, double eta, double *values)
{
/*
  int nBF = 6; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
  */
// first component
  values[0]= 0;
  values[1]= 0;
  values[2]= 0;
  values[3]= 0;
  values[4]= 0;
  values[5]= 0;

  // second component
  values[6]= 0;
  values[7]= 0;
  values[8]= 0;
  values[9]= 0;
  values[10]= 0;
  values[11]= 0;
}

// values of derivatives in eta-eta direction
static void N_T_BDM1_2D_DeriveEtaEta(double xi, double eta, double *values)
{
/*
  int nBF = 6; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
  */
// first component
	values[0]= 0;
	values[1]= 0;
	values[2]= 0;
	values[3]= 0;
	values[4]= 0;
	values[5]= 0;

// second component
	values[6]= 0;
	values[7]= 0;
	values[8]= 0;
	values[9]= 0;
	values[10]= 0;
	values[11]= 0;
}

// values of derivatives in xi-eta direction
static void N_T_BDM1_2D_DeriveXiEta(double xi, double eta, double *values)
{
	/*
  int nBF = 6; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,0,1};
  double mon_y[]={0,0,0,0,0,0,1,0};
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
  */
// first component
	values[0]= 0;
	values[1]= 0;
	values[2]= 0;
	values[3]= 0;
	values[4]= 0;
	values[5]= 0;

	// second component
	values[6]= 0;
	values[7]= 0;
	values[8]= 0;
	values[9]= 0;
	values[10]= 0;
	values[11]= 0;
}

// ***********************************************************************
//TODO die zahlen am ende
/*
  int polynomialdegree,
  int accuracy,
  int n_bf2change,
  int **bf2change,
  int baseVectDim
 */
TBaseFunct2D *BF_N_T_BDM1_2D_Obj = new TBaseFunct2D
        (6, BF_N_T_BDM1_2D, BFUnitTriangle,
         N_T_BDM1_2D_Funct, N_T_BDM1_2D_DeriveXi,
         N_T_BDM1_2D_DeriveEta, N_T_BDM1_2D_DeriveXiXi,
         N_T_BDM1_2D_DeriveXiEta, N_T_BDM1_2D_DeriveEtaEta, 2, 1,
         0, NULL, 2);
