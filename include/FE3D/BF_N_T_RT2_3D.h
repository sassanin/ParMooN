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
// Raviart-Thomas element of second order on tetrahedra, 3D
// ***********************************************************************

static double N_T_RT2_3D_CM[1296] = {
    0,-0,0,-0,0,0,0,-0,0,-0,0,-0,0,-0,0,-0,-0,0,-6,4,-1,4,-1,-1,0,-0,-0,-0,0,-0,-0,-0,0,-0,-0,-0,
    0,0,0,0,0,0,-6,4,-1,4,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -6,4,-1,4,-1,-1,0,0,0,-0,-0,0,-0,0,0,0,-0,-0,-0,0,-0,0,-0,-0,0,-0,-0,-0,0,-0,-0,-0,-0,0,0,0,
    9.5,-9,2.5,-5,2,1,9.5,-5,1,-9,2,2.5,5,-2.5,7.5,-2.5,-2.5,5,36,-16,4,-16,1,4,1080,-1512,-1368,-1368,360,-648,-504,-504,360,-648,-504,-504,
    0,0,0,0,0,0,17.5,-5,0,-25,5,7.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    17.5,-25,7.5,-5,5,-0,0,-0,0,-0,0,-0,0,-0,0,-0,-0,0,0,-0,0,-0,0,0,0,-0,-0,-0,0,-0,-0,-0,0,-0,-0,-0,
    -0,0,0,0,-0,-0,0,0,0,0,-0,-0,0,0,-0,-0,0,-0,17.5,-25,7.5,-5,5,0,-0,0,0,0,-0,0,0,0,-0,0,0,0,
    9.5,-5,1,-9,2,2.5,36,-16,4,-16,1,4,7.5,-2.5,5,-2.5,-2.5,5,9.5,-9,2.5,-5,2,1,360,-504,-648,-504,1080,-1368,-1512,-1368,360,-504,-648,-504,
    17.5,-5,0,-25,5,7.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -0,0,0,0,-0,-0,-0,0,-0,0,-0,-0,-0,0,-0,0,0,-0,17.5,-5,-0,-25,5,7.5,-0,0,0,0,-0,0,0,0,-0,0,0,0,
    0,0,0,0,0,0,17.5,-25,7.5,-5,5,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    36,-16,4,-16,1,4,9.5,-9,2.5,-5,2,1,5,-2.5,5,-2.5,-2.5,7.5,9.5,-5,1,-9,2,2.5,360,-504,-504,-648,360,-504,-504,-648,1080,-1368,-1368,-1512,
    -26,40,-14,6,-6,0,-26,6,0,40,-6,-14,-11,9,-33,7,9,-11,-61.5,22.5,-10,22.5,3.5,-10,-2592,4536,2880,2880,-864,2160,1008,1008,-864,2160,1008,1008,
    0,0,0,0,0,0,-12.5,0,0,25,-0,-12.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -12.5,25,-12.5,0,-0,-0,-0,0,-0,0,-0,-0,-0,0,-0,0,0,-0,-0,0,-0,0,-0,-0,-0,0,0,0,-0,0,0,0,-0,0,0,0,
    -26,11,0,35,-11,-9,-36,13,-3,17,-1,-2,-27,4,-12,8,8,-11,-57,60,-3,12,-12,0,-1728,2016,3744,1872,-1728,2304,3024,2016,-864,1152,2016,1008,
    -26,35,-9,11,-11,0,-57,12,0,60,-12,-3,-12,4,-27,8,8,-11,-36,17,-2,13,-1,-3,-1728,3024,2304,2016,-1728,3744,2016,1872,-864,2016,1152,1008,
    -25,25,0,25,-25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -36,17,-2,13,-1,-3,-26,35,-9,11,-11,-0,-11,8,-12,8,4,-27,-57,12,-0,60,-12,-3,-1728,2016,1872,3744,-864,1152,1008,2016,-1728,2304,2016,3024,
    0,0,0,0,0,0,-25,25,0,25,-25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -57,60,-3,12,-12,-0,-26,11,0,35,-11,-9,-11,8,-27,8,4,-12,-36,13,-3,17,-1,-2,-1728,3024,2016,2304,-864,2016,1008,1152,-1728,3744,1872,2016,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-12.5,25,-12.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -26,6,0,40,-6,-14,-61.5,22.5,-10,22.5,3.5,-10,-33,9,-11,9,7,-11,-26,40,-14,6,-6,0,-864,1008,2160,1008,-2592,2880,4536,2880,-864,1008,2160,1008,
    -12.5,0,0,25,-0,-12.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -0,-0,-0,-0,-0,0,0,-0,0,-0,0,0,0,-0,0,-0,-0,0,-25,25,0,25,-25,-0,0,-0,-0,-0,0,-0,-0,-0,0,-0,-0,-0,
    -36,13,-3,17,-1,-2,-57,60,-3,12,-12,-0,-12,8,-11,4,8,-27,-26,11,-0,35,-11,-9,-864,1008,1152,2016,-1728,1872,2016,3744,-1728,2016,2304,3024,
    -57,12,0,60,-12,-3,-36,17,-2,13,-1,-3,-27,8,-11,4,8,-12,-26,35,-9,11,-11,0,-864,1008,2016,1152,-1728,2016,3024,2304,-1728,1872,3744,2016,
    0,-0,0,-0,0,0,0,-0,0,-0,0,0,0,-0,0,-0,-0,0,-12.5,-0,0,25,-0,-12.5,0,-0,-0,-0,0,-0,-0,-0,0,-0,-0,-0,
    0,0,0,0,0,0,-12.5,25,-12.5,0,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -61.5,22.5,-10,22.5,3.5,-10,-26,40,-14,6,-6,-0,-11,7,-11,9,9,-33,-26,6,-0,40,-6,-14,-864,1008,1008,2160,-864,1008,1008,2160,-2592,2880,2880,4536,
    17.5,-35,17.5,0,0,-0,17.5,0,-0,-35,0,17.5,7,-10.5,31.5,-3.5,-10.5,7,31.5,-10.5,7,-10.5,-3.5,7,1512,-3024,-1512,-1512,504,-1512,-504,-504,504,-1512,-504,-504,
    35,-35,-0,-35,35,-0,42,-7,-0,-35,7,-7,21,14,21,-14,-14,14,42,-35,-7,-7,7,-0,2016,-3024,-4032,-2016,2016,-4032,-3024,-2016,1008,-2016,-2016,-1008,
    42,-35,-7,-7,7,-0,35,-35,-1e-08,-35,35,0,14,-14,21,-14,14,21,42,-7,0,-35,7,-7,2016,-3024,-2016,-4032,1008,-2016,-1008,-2016,2016,-4032,-2016,-3024,
    17.5,0,-0,-35,0,17.5,31.5,-10.5,7,-10.5,-3.5,7,31.5,-10.5,7,-10.5,-3.5,7,17.5,-35,17.5,0,0,-0,504,-504,-1512,-504,1512,-1512,-3024,-1512,504,-504,-1512,-504,
    42,-7,-0,-35,7,-7,42,-35,-7,-7,7,0,21,-14,14,14,-14,21,35,-35,0,-35,35,-1e-08,1008,-1008,-2016,-2016,2016,-2016,-3024,-4032,2016,-2016,-4032,-3024,
    31.5,-10.5,7,-10.5,-3.5,7,17.5,-35,17.5,-0,0,0,7,-3.5,7,-10.5,-10.5,31.5,17.5,-0,0,-35,0,17.5,504,-504,-504,-1512,504,-504,-504,-1512,1512,-1512,-1512,-3024
};

static void N_T_RT2_3D_Funct(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
                  xi*xi,0,0,xi*eta,0,0,xi*zeta,0,0,
                  eta*eta,0,0,eta*zeta,0,0,zeta*zeta,0,0,
                  xi*xi*xi,xi*xi*eta,xi*xi*zeta,
                  xi*eta*eta,xi*eta*zeta,xi*zeta*zeta};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
                  0,xi*xi,0,0,xi*eta,0,0,xi*zeta,0,
                  0,eta*eta,0,0,eta*zeta,0,0,zeta*zeta,0,
                  xi*xi*eta,xi*eta*eta,xi*eta*zeta,
                  eta*eta*eta,eta*eta*zeta,eta*zeta*zeta};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
                  0,0,xi*xi,0,0,xi*eta,0,0,xi*zeta,
                  0,0,eta*eta,0,0,eta*zeta,0,0,zeta*zeta,
                  xi*xi*zeta,xi*eta*zeta,xi*zeta*zeta,
                  eta*eta*zeta,eta*zeta*zeta,zeta*zeta*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0,
                  2*xi,0,0,eta,0,0,zeta,0,0,
                  0,0,0,0,0,0,0,0,0,
                  3*xi*xi,2*xi*eta,2*xi*zeta,
                  eta*eta,eta*zeta,zeta*zeta};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0,
                  0,2*xi,0,0,eta,0,0,zeta,0,
                  0,0,0,0,0,0,0,0,0,
                  2*xi*eta,eta*eta,eta*zeta,
                  0,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0,
                  0,0,2*xi,0,0,eta,0,0,zeta,
                  0,0,0,0,0,0,0,0,0,
                  2*xi*zeta,eta*zeta,zeta*zeta,
                  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0,
                  0,0,0,xi,0,0,0,0,0,
                  2*eta,0,0,zeta,0,0,0,0,0,
                  0,xi*xi,0,
                  xi*2*eta,xi*zeta,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0,
                  0,0,0,0,xi,0,0,0,0,
                  0,2*eta,0,0,zeta,0,0,0,0,
                  xi*xi,xi*2*eta,xi*zeta,
                  3*eta*eta,2*eta*zeta,zeta*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0,
                  0,0,0,0,0,xi,0,0,0,
                  0,0,2*eta,0,0,zeta,0,0,0,
                  0,xi*zeta,0,
                  2*eta*zeta,zeta*zeta,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0,
                  0,0,0,0,0,0,xi,0,0,
                  0,0,0,eta,0,0,2*zeta,0,0,
                  0,0,xi*xi,
                  0,xi*eta,xi*2*zeta};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0,
                  0,0,0,0,0,0,0,xi,0,
                  0,0,0,0,eta,0,0,2*zeta,0,
                  0,0,xi*eta,
                  0,eta*eta,eta*2*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1,
                  0,0,0,0,0,0,0,0,xi,
                  0,0,0,0,0,eta,0,0,2*zeta,
                  xi*xi,xi*eta,xi*2*zeta,
                  eta*eta,eta*2*zeta,3*zeta*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  2,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  3*2*xi,2*eta,2*zeta,
                  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,2,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  2*eta,0,0,
                  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,2,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  2*zeta,0,0,
                  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,1,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,2*xi,0,
                  2*eta,zeta,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,1,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  2*xi,2*eta,zeta,
                  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,1,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,zeta,0,
                  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,1,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,2*xi,
                  0,eta,2*zeta};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,1,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,eta,
                  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,1,
                  0,0,0,0,0,0,0,0,0,
                  2*xi,eta,2*zeta,
                  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  2,0,0,0,0,0,0,0,0,
                  0,0,0,
                  xi*2,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,2,0,0,0,0,0,0,0,
                  0,xi*2,0,
                  3*2*eta,2*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,2,0,0,0,0,0,0,
                  0,0,0,
                  2*zeta,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,1,0,0,0,0,0,
                  0,0,0,
                  0,xi,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,1,0,0,0,0,
                  0,0,xi,
                  0,2*eta,2*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,1,0,0,0,
                  0,xi,0,
                  2*eta,2*zeta,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,2,0,0,
                  0,0,0,
                  0,0,xi*2};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,2,0,
                  0,0,0,
                  0,0,eta*2};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,2,
                  0,0,xi*2,
                  0,eta*2,3*2*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

TBaseFunct3D *BF_N_T_RT2_3D_Obj =
new TBaseFunct3D(36, BF_N_T_RT2_3D, BFUnitTetrahedron,
                 N_T_RT2_3D_Funct, N_T_RT2_3D_DeriveXi,
                 N_T_RT2_3D_DeriveEta, N_T_RT2_3D_DeriveZeta,
                 N_T_RT2_3D_DeriveXiXi, N_T_RT2_3D_DeriveXiEta,
                 N_T_RT2_3D_DeriveXiZeta, N_T_RT2_3D_DeriveEtaEta,
                 N_T_RT2_3D_DeriveEtaZeta, N_T_RT2_3D_DeriveZetaZeta,
                 3, 1,
                 0, NULL, 3);
