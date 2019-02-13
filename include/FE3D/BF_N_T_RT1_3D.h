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
// Raviart-Thomas element of first order on tetrahedra, 3D
// ***********************************************************************

static double N_T_RT1_3D_CM[225] = {
    0,-0,-0,0,0,-0,-0,0,-0,-1.6666667,0.33333333,0.33333333,0,0,0,
    0,0,0,-1.6666667,0.33333333,0.33333333,0,-0,0,0,0,0,0,0,0,
    -1.6666667,0.33333333,0.33333333,0,0,0,0,0,0,0,0,0,0,0,0,
    2.1666667,-0.83333333,-0.33333333,2.1666667,-0.33333333,-0.83333333,-1.1666667,-1.6666667,-1.1666667,5,0.5,0.5,60,30,30,
    0,0,0,2,0,-2,0,0,0,0,0,0,0,0,0,
    2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,-0,0,0,0,-0,-0,-0,2,-2,-0,0,0,0,
    2.1666667,-0.33333333,-0.83333333,5,0.5,0.5,-1.6666667,-1.1666667,-1.1666667,2.1666667,-0.83333333,-0.33333333,30,60,30,
    2,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,2,0,-2,0,0,0,
    0,0,0,2,-2,0,0,0,0,0,0,0,0,0,0,
    5,0.5,0.5,2.1666667,-0.83333333,-0.33333333,-1.1666667,-1.1666667,-1.6666667,2.1666667,-0.33333333,-0.83333333,30,30,60,
    -2.5,2.5,0,-2.5,0,2.5,0.83333333,3.3333333,0.83333333,-3.3333333,-0.83333333,-0.83333333,-60,-30,-30,
    -2.5,-0,2.5,-3.3333333,-0.83333333,-0.83333333,3.3333333,0.83333333,0.83333333,-2.5,2.5,0,-30,-60,-30,
    -3.3333333,-0.83333333,-0.83333333,-2.5,2.5,0,0.83333333,0.83333333,3.3333333,-2.5,-0,2.5,-30,-30,-60
};

static void N_T_RT1_3D_Funct(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
                  xi*xi,xi*eta,xi*zeta};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
                  xi*eta,eta*eta,eta*zeta};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
                  xi*zeta,eta*zeta,zeta*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0,
                  2*xi,eta,zeta};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0,
                  eta,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0,
                  zeta,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0,
                  0,xi,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0,
                  xi,2*eta,zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0,
                  0,zeta,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0,
                  0,0,xi};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0,
                  0,0,eta};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1,
                  xi,eta,2*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  2,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,1,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  1,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0,
                  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,1};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  1,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,2,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,1};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,1,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_RT1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  int nBF = 15; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,2};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

TBaseFunct3D *BF_N_T_RT1_3D_Obj =
new TBaseFunct3D(15, BF_N_T_RT1_3D, BFUnitTetrahedron,
                 N_T_RT1_3D_Funct, N_T_RT1_3D_DeriveXi,
                 N_T_RT1_3D_DeriveEta, N_T_RT1_3D_DeriveZeta,
                 N_T_RT1_3D_DeriveXiXi, N_T_RT1_3D_DeriveXiEta,
                 N_T_RT1_3D_DeriveXiZeta, N_T_RT1_3D_DeriveEtaEta,
                 N_T_RT1_3D_DeriveEtaZeta, N_T_RT1_3D_DeriveZetaZeta,
                 2, 1,
                 0, NULL, 3);
