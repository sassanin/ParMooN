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
// Raviart-Thomas element of first order on hexahedra, 3D
// ***********************************************************************

static double N_H_RT1_3D_CM[1296] = {
//using tschebyscheff points (see NF_N_H_RT1_3D.h)
  0,0,0,0,0,0,0,0,-0.015625,-0.015625,-0.015625,-0.015625,0,0,-0,0,0.015625,0.015625,0.015625,0.015625,0,0,0,0,0.1875,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.015625,0.015625,0.015625,0.015625,0,0,0,0,-0.015625,-0.015625,-0.015625,-0.015625,0,0,0,0,0,0,0,0,0,0,0,0,0.1875,-0,-0,0,0,0,0,0,
  0.015625,0.015625,0.015625,0.015625,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,-0.015625,-0.015625,-0.015625,-0.015625,0,0,0,0,0,0,0,0,0.1875,-0,0,0,
  0,0,0,0,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,-0,0,0,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,-0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,-0,0.5625,-0,-0,0,0,0,0,
  -0.022097087,0.022097087,-0.022097087,0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.022097087,0.022097087,-0.022097087,-0.022097087,0,0,0,0,0,0,0,0,-0,0.5625,-0,-0,
  0,0,0,0,0,0,0,0,0.022097087,0.022097087,-0.022097087,-0.022097087,0,0,0,0,-0.022097087,0.022097087,-0.022097087,0.022097087,0,0,0,0,0,0.5625,0,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,-0,-0,0,0,0,0,
  -0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.022097087,-0.022097087,0.022097087,-0.022097087,0,0,0,0,0,0,0,0,0,-0,0.5625,-0,
  0,0,0,0,0,0,0,0,0.022097087,-0.022097087,0.022097087,-0.022097087,0,0,0,0,-0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,-0,0,0.5625,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.022097087,0.022097087,-0.022097087,0.022097087,0,0,0,0,0.022097087,-0.022097087,0.022097087,-0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,0.5625,-0,0,0,0,0,
  0.03125,0.03125,0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0.046875,0.046875,0.046875,0.046875,0,0,-0,0,-0.046875,-0.046875,-0.046875,-0.046875,0,0,0,0,-0.1875,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.046875,-0.046875,-0.046875,-0.046875,0,0,0,0,0.046875,0.046875,0.046875,0.046875,0,0,0,0,0,0,0,0,0,0,0,0,-0.1875,0,0,-0,0,0,0,0,
  -0.046875,-0.046875,-0.046875,-0.046875,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.046875,0.046875,0.046875,0.046875,0,0,0,0,0,0,0,0,-0.1875,0,-0,-0,
  0,0,0,0,0,0,0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,-0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0.044194174,0.044194174,-0.044194174,-0.044194174,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,0,
  0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.03125,0.03125,0.03125,-0.03125,0,0,0,0,0,0,0,0,0,-0,-0,1.6875,
  0,0,0,0,0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,-0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0,0,1.6875,0,0,0,0,
  -0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0,0,0,0,-0,0,0,0,
  0,0,0,0,0,0,0,0,-0.03125,0.03125,0.03125,-0.03125,0,0,-0,0,0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0,-0,-0,1.6875,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,0,
  -0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,0,0,0,-0,0,0,0,
  0,0,0,0,0,0,0,0,-0.066291261,-0.066291261,0.066291261,0.066291261,0,0,-0,0,0.066291261,-0.066291261,0.066291261,-0.066291261,0,0,0,0,0,-0.5625,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.5625,0,0,0,0,0,0,
  0.066291261,-0.066291261,0.066291261,-0.066291261,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,-0.066291261,-0.066291261,0.066291261,0.066291261,0,0,0,0,0,0,0,0,0,-0.5625,0,0,
  0,0,0,0,0,0,0,0,-0.066291261,0.066291261,-0.066291261,0.066291261,0,0,-0,0,0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,-0,-0,-0.5625,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.066291261,-0.066291261,0.066291261,-0.066291261,0,0,0,0,-0.066291261,0.066291261,-0.066291261,0.066291261,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.5625,0,0,0,0,0,
  0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.066291261,0.066291261,-0.066291261,0.066291261,0,0,0,0,0,0,0,0,-0,0,-0.5625,0,
  0,0,0,0,0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,-0,0,0,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,-0.0625,0.0625,0.0625,-0.0625,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,-0,-0,0,0,0,0,
  0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,-0,0,0,0,-0,0,0,0,-0,
  0,0,0,0,0,0,0,0,0.09375,-0.09375,-0.09375,0.09375,0,0,0,0,-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,0,0,0,-1.6875,0,0,0,-0,0,0,0,-0,
  0,0,0,0,-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-1.6875,0,0,0,-0,
-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.09375,-0.09375,-0.09375,0.09375,0,0,0,0,0,0,0,0,-0,0,0,-1.6875
};

static void N_H_RT1_3D_Funct(double xi, double eta, double zeta,
							 double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
  xi*xi,0,0,xi*eta,0,0,xi*zeta,0,0,eta*zeta,0,0,
  xi*xi*eta,0,0,xi*xi*zeta,0,0,xi*eta*zeta,0,0,
  xi*xi*eta*zeta,0,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
  0,eta*eta,0,0,xi*eta,0,0,xi*zeta,0,0,eta*zeta,0,
  0,xi*eta*eta,0,0,eta*eta*zeta,0,0,xi*eta*zeta,0,
  0,xi*eta*eta*zeta,0};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
  0,0,zeta*zeta,0,0,xi*eta,0,0,xi*zeta,0,0,eta*zeta,
  0,0,xi*zeta*zeta,0,0,eta*zeta*zeta,0,0,xi*eta*zeta,
  0,0,xi*eta*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0,
  2*xi,0,0,eta,0,0,zeta,0,0,0,0,0,
  2*xi*eta,0,0,2*xi*zeta,0,0,eta*zeta,0,0,
  2*xi*eta*zeta,0,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,eta,0,0,zeta,0,0,0,0,
  0,eta*eta,0,0,0,0,0,eta*zeta,0,
  0,eta*eta*zeta,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,eta,0,0,zeta,0,0,0,
  0,0,zeta*zeta,0,0,0,0,0,eta*zeta,
  0,0,eta*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,xi,0,0,0,0,0,zeta,0,0,
  xi*xi,0,0,0,0,0,xi*zeta,0,0,
  xi*xi*zeta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0,
  0,2*eta,0,0,xi,0,0,0,0,0,zeta,0,
  0,xi*2*eta,0,0,2*eta*zeta,0,0,xi*zeta,0,
  0,xi*2*eta*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,xi,0,0,0,0,0,zeta,
  0,0,0,0,0,zeta*zeta,0,0,xi*zeta,
  0,0,xi*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,xi,0,0,eta,0,0,
  0,0,0,xi*xi,0,0,xi*eta,0,0,
  xi*xi*eta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,xi,0,0,eta,0,
  0,0,0,0,eta*eta,0,0,xi*eta,0,
  0,xi*eta*eta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1,
  0,0,2*zeta,0,0,0,0,0,xi,0,0,eta,
  0,0,xi*2*zeta,0,0,eta*2*zeta,0,0,xi*eta,
  0,0,xi*eta*2*zeta};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  2,0,0,0,0,0,0,0,0,0,0,0,
  2*eta,0,0,2*zeta,0,0,0,0,0,
  2*eta*zeta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,0,
  2*xi,0,0,0,0,0,zeta,0,0,
  2*xi*zeta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,
  0,2*eta,0,0,0,0,0,zeta,0,
  0,2*eta*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,zeta,
  0,0,zeta*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,2*xi,0,0,eta,0,0,
  2*xi*eta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,eta,0,
  0,eta*eta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,2*zeta,0,0,0,0,0,eta,
  0,0,eta*2*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,2,0,0,0,0,0,0,0,0,0,0,
  0,xi*2,0,0,2*zeta,0,0,0,0,
  0,xi*2*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,xi,0,0,
  xi*xi,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,2*eta,0,0,xi,0,
  0,xi*2*eta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,1,
  0,0,0,0,0,2*zeta,0,0,xi,
  0,0,xi*2*zeta};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,2,0,0,0,0,0,0,0,0,0,
  0,0,xi*2,0,0,eta*2,0,0,0,
  0,0,xi*eta*2};

  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

TBaseFunct3D *BF_N_H_RT1_3D_Obj =
new TBaseFunct3D(36, BF_N_H_RT1_3D, BFUnitHexahedron,
                 N_H_RT1_3D_Funct, N_H_RT1_3D_DeriveXi,
                 N_H_RT1_3D_DeriveEta, N_H_RT1_3D_DeriveZeta,
                 N_H_RT1_3D_DeriveXiXi, N_H_RT1_3D_DeriveXiEta,
                 N_H_RT1_3D_DeriveXiZeta, N_H_RT1_3D_DeriveEtaEta,
                 N_H_RT1_3D_DeriveEtaZeta, N_H_RT1_3D_DeriveZetaZeta,
                 4, 1, 0, NULL, 3);
