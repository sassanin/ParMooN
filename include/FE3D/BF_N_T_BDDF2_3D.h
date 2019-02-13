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
// Brezzi-Douglas-Duran-Fortin element of second order on tetrahedra, 3D
// ***********************************************************************

static double N_T_BDDF2_3D_CM[900] = {
    0,-0,-0,-0,0,0,-0,-0,0,-0,0,0,-0,0,0,0,0,-0,-6,4,-1,4,-1,-1,0,0,0,0,-0,-0,
    0,0,0,0,0,0,-6,4,-1,4,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    -6,4,-1,4,-1,-1,0,0,0,-0,-0,0,-0,-0,-0,0,-0,-0,-0,0,-0,-0,0,-0,0,0,-0,-0,0,0,
    -2.25,5,-2.75,0.9999999997,-1,5e-10,-2.250000001,1,-2e-10,5,-0.9999999999,-2.75,-1.499999999,0.9999999997,-4,0.9999999996,1,-1.5,17.5,-8.000000001,3.250000001,-8,-6e-10,3.25,72,24,24,504,504,-3.56e-08,
    0,0,0,0,0,0,17.5,-5,0,-25,5,7.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    17.5,-25,7.5,-5,5,-0,0,-0,0,0,0,-0,0,0,-0,-0,0,-0,0,-0,0,0,-0,-0,0,0,0,-0,-0,-0,
    -0,0,-0,0,-0,-0,0,-0,0,-0,-0,-0,0,0,-0,-0,-0,0,17.5,-25,7.5,-5,5,-0,-0,-0,-0,-0,0,0,
    -2.249999998,0.9999999992,2e-10,4.999999999,-0.9999999996,-2.75,17.5,-7.999999999,3.25,-8,-1e-10,3.25,-4.000000001,1.000000001,-1.500000001,0.9999999999,1.000000001,-1.500000001,-2.250000001,5.000000001,-2.750000001,1.000000001,-1,-5e-10,24.00000001,72.00000001,24.00000001,-504.0000001,-6.31e-08,503.9999998,
    17.5,-5,0,-25,5,7.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    -0,0,-0,0,-0,-0,-0,0,-0,0,-0,0,-0,-0,0,-0,-0,0,17.5,-5,-0,-25,5,7.5,-0,-0,-0,0,0,0,
    0,0,0,0,0,0,17.5,-25,7.5,-5,5,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    17.5,-7.999999999,3.25,-7.999999999,-2e-10,3.25,-2.25,4.999999999,-2.75,0.9999999999,-1,-0,-1.5,0.9999999997,-1.499999999,1.000000001,0.9999999993,-3.999999999,-2.249999999,0.9999999995,-1e-10,4.999999999,-0.9999999991,-2.749999999,23.99999999,23.99999999,71.99999999,3.05e-08,-503.9999999,-503.9999998,
    3.25,-9,8.75,3e-10,-3,0.9999999995,3.250000001,-4e-10,1,-9,-3,8.75,2.499999999,-5,10,4e-10,-5,2.5,-11.5,4.000000001,-2.250000001,4,1.000000001,-2.25,-72,-24,-24,-504,-504,3.56e-08,
    0,0,0,0,0,0,-12.5,0,0,25,-0,-12.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    -12.5,25,-12.5,0,-0,-0,-0,0,-0,0,-0,-0,0,-0,0,0,-0,0,0,-0,-0,-0,0,-0,-0,-0,-0,0,0,0,
    3.249999998,-4.5,-1.375,-4.499999998,9.500000001,-1.375000002,1.500000001,-0.500000001,-0.6249999997,-2.5,1.5,-2.375,-1.5e-09,7.500000001,-1e-10,-2.499999999,-2.5,2.5,-22.25,22.5,-3.625000002,4.499999999,-3.499999998,-0.6250000001,-48,-48,-24,-1260,-252,-251.9999999,
    3.249999999,-4.499999999,-1.375,-4.500000001,9.499999999,-1.374999999,-22.25,4.5,-0.6250000001,22.5,-3.5,-3.625,1.7e-09,7.499999999,6e-10,-2.500000001,-2.5,2.500000001,1.500000001,-2.500000003,-2.374999998,-0.5,1.499999999,-0.6249999994,-48.00000001,-48.00000001,-24.00000001,1260,-252,-251.9999999,
    -25,25,0,25,-25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    1.500000001,-2.5,-2.375,-0.500000001,1.5,-0.625,3.25,-4.5,-1.375,-4.5,9.5,-1.375,2.5,-2.5,-2e-10,-2.5,7.5,-1e-10,-22.25,4.500000001,-0.6250000004,22.5,-3.5,-3.625,-48,-24,-48,-252,-1260,252,
    0,0,0,0,0,0,-25,25,0,25,-25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    -22.25,22.5,-3.625,4.5,-3.5,-0.6250000004,3.25,-4.499999999,-1.375,-4.5,9.5,-1.375,2.499999999,-2.5,-4e-10,-2.5,7.500000001,-7e-10,1.499999999,-0.4999999994,-0.6250000004,-2.499999999,1.5,-2.375000001,-48,-23.99999999,-47.99999999,-252,1260,251.9999999,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-12.5,25,-12.5,0,0,0,0,0,0,0,-0,0,
    3.249999998,8e-10,0.9999999998,-8.999999999,-3,8.75,-11.5,3.999999999,-2.25,4,1,-2.25,10,-5.000000001,2.500000001,-5,-5e-10,2.500000001,3.250000001,-9.000000001,8.750000001,-9e-10,-3,1.000000001,-24.00000001,-72.00000001,-24.00000001,504.0000001,6.31e-08,-503.9999998,
    -12.5,0,0,25,-0,-12.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    0,-0,-0,-0,-0,-0,-0,0,-0,0,0,0,0,0,0,0,0,-0,-25,25,0,25,-25,-0,0,0,0,0,-0,-0,
    1.499999994,-0.499999998,-0.6250000007,-2.499999996,1.499999999,-2.375,-22.25,22.5,-3.625,4.5,-3.5,-0.6250000001,7e-10,-2.500000001,2.500000001,7.500000001,-2.500000001,1.6e-09,3.250000002,-4.500000002,-1.375,-4.500000003,9.500000002,-1.374999999,-24.00000001,-48.00000002,-48.00000003,252.0000001,252.0000002,-1260,
    -22.24999999,4.499999998,-0.6249999992,22.5,-3.5,-3.624999999,1.499999999,-2.499999998,-2.375,-0.4999999997,1.5,-0.625,-1e-10,-2.499999999,2.499999999,7.499999998,-2.499999998,-1.6e-09,3.249999999,-4.499999999,-1.374999999,-4.499999997,9.499999997,-1.375000001,-23.99999999,-47.99999998,-47.99999997,251.9999999,251.9999998,1260,
    0,-0,-0,-0,0,0,0,-0,0,-0,0,0,-0,0,0,0,0,-0,-12.5,-0,0,25,-0,-12.5,0,0,0,0,-0,-0,
    0,0,0,0,0,0,-12.5,25,-12.5,0,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    -11.5,3.999999999,-2.25,3.999999999,1,-2.25,3.25,-8.999999999,8.75,1e-10,-3,1,2.5,3e-10,2.499999999,-5.000000001,-4.999999999,9.999999999,3.249999999,5e-10,1,-8.999999999,-3.000000001,8.749999999,-23.99999999,-23.99999999,-71.99999999,-3.05e-08,503.9999999,503.9999998
};

static void N_T_BDDF2_3D_Funct(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
                  xi*xi,0,0,xi*eta,0,0,xi*zeta,0,0,
                  eta*eta,0,0,eta*zeta,0,0,zeta*zeta,0,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
                  0,xi*xi,0,0,xi*eta,0,0,xi*zeta,0,
                  0,eta*eta,0,0,eta*zeta,0,0,zeta*zeta,0};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
                  0,0,xi*xi,0,0,xi*eta,0,0,xi*zeta,
                  0,0,eta*eta,0,0,eta*zeta,0,0,zeta*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0,
                  2*xi,0,0,eta,0,0,zeta,0,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0,
                  0,2*xi,0,0,eta,0,0,zeta,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0,
                  0,0,2*xi,0,0,eta,0,0,zeta,
                  0,0,0,0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0,
                  0,0,0,xi,0,0,0,0,0,
                  2*eta,0,0,zeta,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0,
                  0,0,0,0,xi,0,0,0,0,
                  0,2*eta,0,0,zeta,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0,
                  0,0,0,0,0,xi,0,0,0,
                  0,0,2*eta,0,0,zeta,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0,
                  0,0,0,0,0,0,xi,0,0,
                  0,0,0,eta,0,0,2*zeta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0,
                  0,0,0,0,0,0,0,xi,0,
                  0,0,0,0,eta,0,0,2*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1,
                  0,0,0,0,0,0,0,0,xi,
                  0,0,0,0,0,eta,0,0,2*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  2,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,2,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,2,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,1,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,1,0,0,0,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,1,0,0,0,
                  0,0,0,0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,1,0,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,1,0,
                  0,0,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,1,
                  0,0,0,0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  2,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,2,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,2,0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,1,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,1,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,1,0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  int nBF = 30; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,2,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,2,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,2};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

TBaseFunct3D *BF_N_T_BDDF2_3D_Obj =
new TBaseFunct3D(30, BF_N_T_BDDF2_3D, BFUnitTetrahedron,
                 N_T_BDDF2_3D_Funct, N_T_BDDF2_3D_DeriveXi,
                 N_T_BDDF2_3D_DeriveEta, N_T_BDDF2_3D_DeriveZeta,
                 N_T_BDDF2_3D_DeriveXiXi, N_T_BDDF2_3D_DeriveXiEta,
                 N_T_BDDF2_3D_DeriveXiZeta, N_T_BDDF2_3D_DeriveEtaEta,
                 N_T_BDDF2_3D_DeriveEtaZeta, N_T_BDDF2_3D_DeriveZetaZeta,
                 2, 1,
                 0, NULL, 3);
