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
// functions 10,11 and 12 correspond to inner degrees of freedom

// coefficient matrix

// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in
// NF_N_T_BDM2_2D.h

/*
// using equidistant points on edges and point evaluation for inner dofs
static double N_T_BDM2_2D_CM[144] = {
0,0,-0,-0,0,-0,-1,3,-3,-0,0,0,
-3,3,-1,0,-0,0,0,-0,0,-0,-0,-0,
-0.047619048,-0.047619048,1.3809524,0.23809524,-1.7619048,1.0952381,2.8095238,-9.1904762,14.52381,9.1428571,2.2857143,-0.57142857,
10,-16,6,0,-0,0,0,-0,0,0,0,-0,
1.047619,-2.952381,1.6190476,2.7619048,-1.2380952,-0.095238095,-1.8095238,6.1904762,-11.52381,-9.1428571,-2.2857143,0.57142857,
-8,16,-8,-0,0,-0,-0,0,0,0,-0,0,
0,0,0,0,-0,0,6,-16,10,0,0,0,
5.7619048,-4.2380952,0.9047619,-0.80952381,3.1904762,-8.5238095,-3.952381,10.047619,-5.3809524,-2.2857143,-0.57142857,9.1428571,
0,0,0,0,0,0,-8,16,-8,0,0,0,
-2.7619048,1.2380952,0.095238095,1.8095238,-6.1904762,11.52381,6.952381,-13.047619,6.3809524,2.2857143,0.57142857,-9.1428571,
-0.76190476,3.2380952,-9.9047619,-4.1904762,11.809524,-6.4761905,-7.047619,20.952381,-19.619048,-9.7142857,-11.428571,2.8571429,
-8.952381,11.047619,3.6190476,0.76190476,-3.2380952,9.9047619,4.1904762,-11.809524,6.4761905,2.8571429,9.7142857,-11.428571
};
*/
/*
// using equidistant points on edges and integral evaluation for inner dofs
static double N_T_BDM2_2D_CM[144] = {
0,0,0,0,0,0,-1,3,-3,0,0,-0,
-3,3,-1,-0,0,0,0,0,0,0,0,-0,
-1.5,4,-2.5,-3,1,-1,2.5,-6,9.5,18,6,90,
10,-16,6,-0,0,0,0,0,0,0,0,-0,
2.5,-7,5.5,6,-4,2,-1.5,3,-6.5,-18,-6,-90,
-8,16,-8,0,0,0,0,0,0,0,0,-0,
0,0,0,0,0,0,6,-16,10,0,0,-0,
9.5,-6,2.5,-1,1,-3,-2.5,4,-1.5,6,18,-90,
0,0,0,0,0,0,-8,16,-8,0,0,-0,
-6.5,3,-1.5,2,-4,6,5.5,-7,2.5,-6,-18,90,
1,-2,-1,-0,4,-0,-3,14,-13,-12,-12,-180,
-13,14,-3,0,4,0,-1,-2,1,-12,-12,180
};
*/
/*
// using Gauss-Points on edges and poinit evaluation for inner dofs
static double N_T_BDM2_2D_CM[144] = {
0,0,0,0,0,0,-0.18783611,0.66666667,-1.4788306,0,0,0,
-1.4788306,0.66666667,-0.18783611,-0,0,0,0,-0,0,-0,-0,-0,
-0.18329167,0.73015873,0.73884722,0.0011361102,-0.98412698,0.55441945,-0.16965835,0.92063492,7.3918806,9.1428571,2.2857143,-0.57142857,
4.6243278,-6.6666667,2.0423389,0,0,0,0,-0,0,0,0,-0,
0.37112778,-1.3968254,0.73998333,1.4776944,0.31746032,-0.36658334,0.35749445,-1.5873016,-5.91305,-9.1428571,-2.2857143,0.57142857,
-3.3333333,6.6666667,-3.3333333,0,-0,-0,0,0,-0,0,0,0,
0,0,0,0,0,0,2.0423389,-6.6666667,4.6243278,0,0,0,
2.956525,-0.34920635,-0.17874723,0.54533056,-2.2539683,-4.4342195,-1.483375,4.6031746,-2.4055139,-2.2857143,-0.57142857,9.1428571,
0,0,0,0,0,0,-3.3333333,6.6666667,-3.3333333,0,0,-0,
-1.4776944,-0.31746032,0.36658334,-0.35749445,1.5873016,5.91305,2.9622056,-5.2698413,2.59335,2.2857143,0.57142857,-9.1428571,
0.72862223,-2.984127,-5.1730667,-1.4845111,5.5873016,-2.9599333,-1.4981444,5.3968254,-9.6129667,-9.7142857,-11.428571,2.8571429,
-5.1685222,7.9365079,2.9463,-0.72862223,2.984127,5.1730667,1.4845111,-5.5873016,2.9599333,2.8571429,9.7142857,-11.428571
};
*/
/*
// using Gauss-Points on edges and integral evaluation for inner dofs
static double N_T_BDM2_2D_CM[144] = {
0,0,0,0,0,0,-0.18783611,0.66666667,-1.4788306,0,0,-0,
-1.4788306,0.66666667,-0.18783611,0,0,0,0,-0,0,-0,-0,-0,
-0.51058472,1.6666667,-1.1560819,-1.4788306,-1.3333333,-0.18783611,0.24075971,1,4.7592403,18,6,90,
4.6243278,-6.6666667,2.0423389,0,0,0,0,-0,0,-0,-0,0,
0.69842083,-2.3333333,2.6349125,2.9576611,0.66666667,0.37567222,-0.052923606,-1.6666667,-3.2804097,-18,-6,-90,
-3.3333333,6.6666667,-3.3333333,-0,-0,-0,0,-0,-0,0,0,-0,
0,0,0,0,0,0,2.0423389,-6.6666667,4.6243278,0,0,-0,
4.7592403,1,0.24075971,-0.18783611,-1.3333333,-1.4788306,-1.1560819,1.6666667,-0.51058472,6,18,-90,
0,0,0,0,0,0,-3.3333333,6.6666667,-3.3333333,0,0,-0,
-3.2804097,-1.6666667,-0.052923606,0.37567222,0.66666667,2.9576611,2.6349125,-2.3333333,0.69842083,-6,-18,90,
0.64549722,-2,-0.64549722,-0,4,-0,-0.10584721,4.6666667,-6.5608195,-12,-12,-180,
-6.5608195,4.6666667,-0.10584721,0,4,0,-0.64549722,-2,0.64549722,-12,-12,180
};
*/

// using Tschebyscheff-points on edges and point evaluation for inner dofs
static double N_T_BDM2_2D_CM[144] = {
0,-0,0,-0,-0,0,-0.089316398,0.33333333,-1.2440169,0,0,-0,
-1.2440169,0.33333333,-0.089316398,0,-0,-0,-0,0,-0,-0,-0,0,
-0.19017083,0.84126984,0.63461527,-0.025213607,-0.87301587,0.46965805,-0.49273412,2.3650794,6.2705119,9.1428571,2.2857143,-0.57142857,
3.8213672,-5.3333333,1.5119661,-0,0,0,0,-0,0,0,0,-0,
0.27948722,-1.1746032,0.60940166,1.2692305,0.53968254,-0.38034165,0.58205051,-2.6984127,-5.026495,-9.1428571,-2.2857143,0.57142857,
-2.6666667,5.3333333,-2.6666667,0,-0,-0,-0,0,-0,-0,-0,0,
-0,-0,-0,-0,0,-0,1.5119661,-5.3333333,3.8213672,-0,-0,0,
2.5132475,0.20634921,-0.29102526,0.67136691,-3.031746,-3.782478,-1.1431625,3.8253968,-1.9679486,-2.2857143,-0.57142857,9.1428571,
0,0,0,0,-0,0,-2.6666667,5.3333333,-2.6666667,0,0,-0,
-1.2692305,-0.53968254,0.38034165,-0.58205051,2.6984127,5.026495,2.3871794,-4.1587302,2.057265,2.2857143,0.57142857,-9.1428571,
0.86153774,-3.8730159,-4.4170933,-1.1179489,4.6984127,-2.4376067,-0.81538561,3.1746032,-8.0735033,-9.7142857,-11.428571,2.8571429,
-4.5179477,7.4920635,2.7401699,-0.86153774,3.8730159,4.4170933,1.1179489,-4.6984127,2.4376067,2.8571429,9.7142857,-11.428571
};

/*
// using Tschebyscheff-points on edges and integral evaluation for inner dofs
static double N_T_BDM2_2D_CM[144] = {
0,-0,-0,-0,-0,-0,-0.089316398,0.33333333,-1.2440169,0,0,0,
-1.2440169,0.33333333,-0.089316398,-0,-0,0,0,0,-0,0,0,-0,
-0.37799153,1.3333333,-0.9553418,-1.2440169,-1.6666667,-0.089316398,-0.020725942,2,4.0207259,18,6,90,
3.8213672,-5.3333333,1.5119661,-0,0,-0,-0,-0,0,-0,-0,0,
0.46730793,-1.6666667,2.1993587,2.4880339,1.3333333,0.1786328,0.11004234,-2.3333333,-2.776709,-18,-6,-90,
-2.6666667,5.3333333,-2.6666667,-0,-0,-0,-0,0,-0,0,0,-0,
-0,-0,0,-0,0,0,1.5119661,-5.3333333,3.8213672,-0,-0,-0,
4.0207259,2,-0.020725942,-0.089316398,-1.6666667,-1.2440169,-0.9553418,1.3333333,-0.37799153,6,18,-90,
0,-0,-0,-0,-0,-0,-2.6666667,5.3333333,-2.6666667,0,0,0,
-2.776709,-2.3333333,0.11004234,0.1786328,1.3333333,2.4880339,2.1993587,-1.6666667,0.46730793,-6,-18,90,
0.57735027,-2,-0.57735027,-0,4,-0,0.22008468,3.3333333,-5.553418,-12,-12,-180,
-5.553418,3.3333333,0.22008468,0,4,0,-0.57735027,-2,0.57735027,-12,-12,180
};
*/
static void N_T_BDM2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={1, 0, xi, 0, xi*xi, 0, eta, 0, eta*eta, 0, xi*eta, 0};
  double mon_y[12]={0, 1, 0, xi, 0, xi*xi, 0, eta, 0, eta*eta, 0, xi*eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_BDM2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0, 0 , 1, 0, 2*xi, 0, 0, 0, 0, 0, eta, 0};
  double mon_y[12]={0, 0, 0 , 1, 0, 2*xi, 0, 0, 0, 0, 0, eta};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_BDM2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0, 0, 0, 0, 0, 0, 1, 0, 2*eta, 0, xi, 0};
  double mon_y[12]={0, 0, 0, 0, 0, 0, 0, 1, 0, 2*eta, 0, xi};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_BDM2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0, 0 , 0, 0, 2, 0, 0, 0, 0, 0, 0, 0};
  double mon_y[12]={0, 0, 0 , 0, 0, 2, 0, 0, 0, 0, 0, 0};
 
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_BDM2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0};
  double mon_y[12]={0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_BDM2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
  double mon_y[12]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  
  memset(values, 0.0, 2*nBF*SizeOfDouble); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// ***********************************************************************

TBaseFunct2D *BF_N_T_BDM2_2D_Obj = new TBaseFunct2D
        (12, BF_N_T_BDM2_2D, BFUnitTriangle,
         N_T_BDM2_2D_Funct, N_T_BDM2_2D_DeriveXi,
         N_T_BDM2_2D_DeriveEta, N_T_BDM2_2D_DeriveXiXi,
         N_T_BDM2_2D_DeriveXiEta, N_T_BDM2_2D_DeriveEtaEta, 2, 2,
         0, NULL, 2);
