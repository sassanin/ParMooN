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
// P2 Raviart-Thomas vector element, nonconforming , 2D
// History:  14.11.2011 implementation (Ulrich)
// ***********************************************************************

// base function values
// vector function, orthonormal to edges, 
// functions 1,2, and 3 are orthogonal to edge 1
// functions 4,5, and 6 are orthogonal to edge 2
// functions 7,8, and 9 are orthogonal to edge 3
// functions 10,11,12,13,14, and 15 correspond to inner degrees of freedom

// coefficient matrix

// using equidistant points on edges and integral evaluation for inner dofs
static double N_T_RT2_2D_CM[225] = {
 0,0,0,0,0,0,-1,3,-3,0,0,0,0,0,0,
-3,3,-1,-0,0,0,0,0,0,0,0,0,0,0,0,
 5.5,-7,2.5,5,-2,3,3.5,-13,20.5,180,60,-270,-120,-240,-90,
 10,-16,6,-0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,6,-16,10,0,0,0,0,0,0,
 20.5,-13,3.5,3,-2,5,2.5,-7,5.5,60,180,-90,-240,-120,-270,
-16.5,28,-11.5,-21.5,8,-6.5,-7,19,-37,-450,-150,810,390,510,180,
-8,16,-8,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,-8,16,-8,0,0,0,0,0,0,
-37,19,-7,-6.5,8,-21.5,-11.5,28,-16.5,-150,-450,180,510,390,810,
-22,14,-2,-7,4,-17,-5,40,-35,-300,-300,360,420,660,540,
-35,40,-5,-17,4,-7,-2,14,-22,-300,-300,540,660,420,360,
 12,-24,12,19.5,-9,4.5,4.5,-9,19.5,270,90,-540,-270,-270,-90,
 27,-24,-3,12,6,12,-3,-24,27,360,360,-540,-720,-720,-540,
 19.5,-9,4.5,4.5,-9,19.5,12,-24,12,90,270,-90,-270,-270,-540
};

// using Gauss-Points on edges and integral evaluation for inner dofs
/*static double N_T_RT2_2D_CM[225] = { 
0,0,0,0,0,0,-0.187836,0.666667,-1.478831,0,0,0,0,0,0,
-1.478831,0.666667,-0.187836,0,0,0,0,0,0,0,0,0,0,0,0,
2.634913,-2.333333,0.698421,2.312164,2.666667,1.021169,-0.486726,1,10.486726,180,60,-270,-120,-240,-90,
4.624328,-6.666667,2.042339,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,2.042339,-6.666667,4.624328,0,0,0,0,0,0,
10.486726,1,-0.486726,1.021169,2.666667,2.312164,0.698421,-2.333333,2.634913,60,180,-90,-240,-120,-270,
-7.447076,11.666667,-4.21959,-10.674563,-8.333333,-0.992104,0.515792,-6.666667,-18.849125,-450,-150,810,390,510,180,
-3.333333,6.666667,-3.333333,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-3.333333,6.666667,-3.333333,0,0,0,0,0,0,
-18.849125,-6.666667,0.515792,-0.992104,-8.333333,-10.674563,-4.21959,11.666667,-7.447076,-150,-450,180,510,390,810,
-11.454972,-0,1.454972,-1.772514,-10,-8.227486,1.349125,16.666667,-18.015792,-300,-300,360,420,660,540,
-18.015792,16.666667,1.349125,-8.227486,-10,-1.772514,1.454972,0,-11.454972,-300,-300,540,660,420,360,
5,-10,5,9.841229,5,0.158771,0.158771,5,9.841229,270,90,-540,-270,-270,-90,
14.682458,-10,-4.682458,5,20,5,-4.682458,-10,14.682458,360,360,-540,-720,-720,-540,
9.841229,5,0.158771,0.158771,5,9.841229,5,-10,5,90,270,-90,-270,-270,-540
};*/
/*
// using Tschebyscheff-points on edges and integral evaluation for inner dofs
static double N_T_RT2_2D_CM[225] = { 
0,0,0,0,0,0,-0.0893164,0.33333333,-1.24401694,0,0,0,0,0,0,
-1.24401694,0.33333333,-0.0893164,0,0,0,0,0,0,0,0,0,0,0,0,
2.19935874,-1.66666667,0.46730793,1.9106836,3.33333333,0.75598307,-0.90747729,3,8.90747729,180,60,-270,-120,-240,-90,
3.82136721,-5.33333333,1.51196613,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,1.51196613,-5.33333333,3.82136721,0,0,0,0,0,0,
8.90747729,3,-0.90747729,0.75598306,3.33333334,1.91068361,0.46730793,-1.66666667,2.19935874,60,180,-90,-240,-120,-270,
-6.11004234,9.33333333,-3.22329099,-8.99679368,-10.66666667,-0.33653965,1.3269207,-10.33333333,-15.99358737,-450,-150,810,390,510,180,
-2.66666667,5.33333333,-2.66666667,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-2.66666667,5.33333333,-2.66666667,0,0,0,0,0,0,
-15.99358738,-10.33333334,1.3269207,-0.33653965,-10.66666667,-8.9967937,-3.223291,9.33333334,-6.11004234,-150,-450,180,510,390,810,
-9.7735027,-2,1.77350269,-1.11324865,-12,-6.88675135,1.99358737,13.33333334,-15.32692071,-300,-300,360,420,660,540,
-15.32692071,13.33333333,1.99358737,-6.88675134,-12,-1.11324866,1.77350269,-2,-9.77350269,-300,-300,540,660,420,360,
4,-8,4,8.33012702,7,-0.33012702,-0.33012702,7,8.33012702,270,90,-540,-270,-270,-90,
12.66025404,-8,-4.66025404,4,22,4,-4.66025404,-8,12.66025404,360,360,-540,-720,-720,-540,
8.33012703,7,-0.33012702,-0.33012702,7,8.33012703,4,-8,4,90,270,-90,-270,-270,-540
};
*/
/*
// using Tschebyscheff-points on edges and point evaluation for inner dofs
static double N_T_RT2_2D_CM[225] = { 
0,0,0,0,0,0,-0.0893163975,0.3333333333,-1.2440169359,0,0,0,0,0,0,
-1.2440169359,0.3333333333,-0.0893163975,0,0,0,0,0,0,0,0,0,0,0,0,
1.2440169359,-0.3333333333,0.0893163975,0.6666666667,1.1666666667,0.6666666667,-1.2854688202,1.8333333333,7.9521354869,12.5,2.5,-2.5,-2,-4,-2,
3.821367205,-5.3333333333,1.5119661283,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,1.5119661283,-5.3333333333,3.821367205,0,0,0,0,0,0,
7.9521354869,1.8333333333,-1.2854688202,0.6666666667,1.1666666667,0.6666666667,0.0893163975,-0.3333333333,1.2440169359,2.5,12.5,-2,-4,-2,-2.5,
-3.821367205,5.3333333333,-1.5119661283,-5.5534180126,-5.1666666667,0.2200846793,2.1722201661,-5.8333333333,-12.8388868328,-27.5,-5.5,17.5,14,4,2,
-2.6666666667,5.3333333333,-2.6666666667,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-2.6666666667,5.3333333333,-2.6666666667,0,0,0,0,0,0,
-12.8388868328,-5.8333333333,2.1722201661,0.2200846793,-5.1666666667,-5.5534180126,-1.5119661283,5.3333333333,-3.821367205,-5.5,-27.5,2,4,14,17.5,
-6.708118551,-2.1666666667,1.3747852177,0.2200846793,-5.1666666667,-5.5534180126,2.7495704353,10.6666666667,-13.416237102,-17.5,-15.5,2,4,26,14.5,
-13.416237102,10.6666666667,2.7495704353,-5.5534180126,-5.1666666667,0.2200846793,1.3747852177,-2.1666666667,-6.708118551,-15.5,-17.5,14.5,26,4,2,
2.6666666667,-5.3333333333,2.6666666667,6.1307682818,3.6666666667,-0.7974349485,-0.7974349485,3.6666666667,6.1307682818,15,3,-15,-12,0,0,
9.5948698969,-5.3333333333,-4.2615365636,2.6666666667,12.6666666667,2.6666666667,-4.2615365636,-5.3333333333,9.5948698969,18,18,-12,-24,-24,-12,
6.1307682818,3.6666666667,-0.7974349485,-0.7974349485,3.6666666667,6.1307682818,2.6666666667,-5.3333333333,2.6666666667,3,15,0,0,-12,-15
};
*/
static void N_T_RT2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={1,0,xi,0,eta,0,xi*xi,0,eta*eta,0,xi*eta,0,xi*xi*xi, xi*xi*eta,xi*eta*eta};
  double mon_y[15]={0,1,0,xi,0,eta,0,xi*xi,0,eta*eta,0,xi*eta,xi*xi*eta, xi*eta*eta,eta*eta*eta };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_RT2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  //Dscal(nBF*nBF,1.0,mat);
  // monomials x-component and y-component
  double mon_x[15]={0,0,1,0,0,0,2*xi,0,0,0,eta,0,3*xi*xi,2*xi*eta,eta*eta};
  double mon_y[15]={0,0,0,1,0,0,0,2*xi,0,0,0,eta,2*xi*eta,eta*eta,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_RT2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0,0,0,1,0,0,0,2*eta,0,xi,0,0,xi*xi,2*xi*eta};
  double mon_y[15]={0,0,0,0,0,1,0,0,0,2*eta,0,xi,xi*xi,2*xi*eta,3*eta*eta };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_RT2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0,0,0,0,0,2,0,0,0,0,0,6*xi,2*eta,0};
  double mon_y[15]={0,0,0,0,0,0,0,2,0,0,0,0,2*eta,0,0};
 
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_RT2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*xi};
  double mon_y[15]={0,0,0,0,0,0,0,0,0,2,0,0,0,2*xi,6*eta };
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_RT2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component and y-component
  double mon_x[15]={0,0,0,0,0,0,0,0,0,0,1,0,0,2*xi,2*eta};
  double mon_y[15]={0,0,0,0,0,0,0,0,0,0,0,1,2*xi,2*eta,0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// ***********************************************************************

TBaseFunct2D *BF_N_T_RT2_2D_Obj = new TBaseFunct2D
        (15, BF_N_T_RT2_2D, BFUnitTriangle, 
         N_T_RT2_2D_Funct, N_T_RT2_2D_DeriveXi,
         N_T_RT2_2D_DeriveEta, N_T_RT2_2D_DeriveXiXi,
         N_T_RT2_2D_DeriveXiEta, N_T_RT2_2D_DeriveEtaEta, 3, 2,
         0, NULL, 2);
