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
   
// red refined triangle, second order, used for LPS

// ***********************************************************************
// RR2 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_RR2_2D_Funct(double xi, double eta, double *values)
{
  if(xi>=0.5)
  {
    // right triangle
    values[0] = 0;
    values[1] = 0;
    values[2] = 8 * (xi + eta - 1) * (xi + eta - 0.75);
    values[3] = -16 * (xi - 0.5) * (xi + eta - 1);
    values[4] = 8 * (xi - 0.5) * (xi - 0.75);
    values[5] = 16 * (xi - 0.5) * eta;
    values[6] = 8 * eta * (-0.25 + eta);
    values[7] = 0;
    values[8] = 0;
    values[9] = 0;
    values[10] = 0;
    values[11] = 0;
    values[12] = 0;
    values[13] = -16 * eta * (xi + eta - 1);
    values[14] = 0;
  }
  else
  {
    if(eta>=0.5)
    {
      // upper triangle
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 8 * xi * (xi - 0.25);
      values[7] = 16 * (-0.5 + eta) * xi;
      values[8] = 8 * (-0.5 + eta) * (eta - 0.75);
      values[9] = -16 * (-0.5 + eta) * (xi + eta - 1);
      values[10] = 8 * (xi + eta - 1) * (xi + eta - 0.75);
      values[11] = 0;
      values[12] = 0;
      values[13] = 0;
      values[14] = -16 * xi * (xi + eta - 1);
    }
    else
    {
      if(xi+eta>=0.5)
      {
        // middle triangle
        values[0] = 0;
        values[1] = 0;
        values[2] = 8 * (-0.5 + eta) * (-0.25 + eta);
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 8 * (xi + eta - 0.5) * (xi + eta - 0.75);
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 8 * (xi - 0.25) * (xi - 0.5);
        values[11] = 0;
        values[12] = 16 * (-0.5 + eta) * (xi - 0.5);
        values[13] = -16 * (xi + eta - 0.5) * (-0.5 + eta);
        values[14] = -16 * (xi + eta - 0.5) * (xi - 0.5);
      }
      else
      {
        // lower left triangle
        values[0] = 8 * (xi + eta - 0.5) * (xi + eta - 0.25);
        values[1] = -16 * xi * (xi + eta - 0.5);
        values[2] = 8 * xi * (xi - 0.25);
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 0;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 8 * eta * (-0.25 + eta);
        values[11] = -16 * eta * (xi + eta - 0.5);
        values[12] = 16 * eta * xi;
        values[13] = 0;
        values[14] = 0;
      }
    }
  }
}

// values of the derivatives in xi direction
static void C_T_RR2_2D_DeriveXi(double xi, double eta, double *values)
{
  if(xi>=0.5)
  {
    // right triangle
    values[0] = 0;
    values[1] = 0;
    values[2] =  (16 * xi) +  (16 * eta) - 14;
    values[3] = - (32 * xi) -  (16 * eta) + 24;
    values[4] =  (16 * xi) - 10;
    values[5] = 16 * eta;
    values[6] = 0;
    values[7] = 0;
    values[8] = 0;
    values[9] = 0;
    values[10] = 0;
    values[11] = 0;
    values[12] = 0;
    values[13] = -16 * eta;
    values[14] = 0;
  }
  else
  {
    if(eta>=0.5)
    {
      // upper triangle
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] =  (16 * xi) - 2;
      values[7] = -8 +  (16 * eta);
      values[8] = 0;
      values[9] = 8 -  (16 * eta);
      values[10] =  (16 * xi) +  (16 * eta) - 14;
      values[11] = 0;
      values[12] = 0;
      values[13] = 0;
      values[14] = - (32 * xi) -  (16 * eta) + 16;
    }
    else
    {
      if(xi+eta>=0.5)
      {
        // middle triangle
        values[0] = 0;
        values[1] = 0;
        values[2] = 0;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] =  (16 * xi) +  (16 * eta) - 10;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] =  (16 * xi) - 6;
        values[11] = 0;
        values[12] = -8 +  (16 * eta);
        values[13] = 8 -  (16 * eta);
        values[14] = - (32 * xi) -  (16 * eta) + 16;
      }
      else
      {
        // lower left triangle
        values[0] =  (16 * xi) +  (16 * eta) - 6;
        values[1] = - (32 * xi) -  (16 * eta) + 8;
        values[2] =  (16 * xi) - 2;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 0;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 0;
        values[11] = -16 * eta;
        values[12] = 16 * eta;
        values[13] = 0;
        values[14] = 0;
      }
    }
  }
}

// values of the derivatives in eta direction
static void C_T_RR2_2D_DeriveEta(double xi, double eta, double *values)
{
  if(xi>=0.5)
  {
    // right triangle
    values[0] = 0;
    values[1] = 0;
    values[2] =  (16 * xi) +  (16 * eta) - 14;
    values[3] = - (16 * xi) + 8;
    values[4] = 0;
    values[5] =  (16 * xi) - 8;
    values[6] = -2 +  (16 * eta);
    values[7] = 0;
    values[8] = 0;
    values[9] = 0;
    values[10] = 0;
    values[11] = 0;
    values[12] = 0;
    values[13] = - (16 * xi) -  (32 * eta) + 16;
    values[14] = 0;
  }
  else
  {
    if(eta>=0.5)
    {
      // upper triangle
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 0;
      values[7] = 16 * xi;
      values[8] =  (16 * eta) - 10;
      values[9] = - (16 * xi) -  (32 * eta) + 24;
      values[10] =  (16 * xi) +  (16 * eta) - 14;
      values[11] = 0;
      values[12] = 0;
      values[13] = 0;
      values[14] = -16 * xi;
    }
    else
    {
      if(xi+eta>=0.5)
      {
        // middle triangle
        values[0] = 0;
        values[1] = 0;
        values[2] = -6 +  (16 * eta);
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] =  (16 * xi) +  (16 * eta) - 10;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 0;
        values[11] = 0;
        values[12] =  (16 * xi) - 8;
        values[13] = - (16 * xi) -  (32 * eta) + 16;
        values[14] = - (16 * xi) + 0.80e1;
      }
      else
      {
        // lower left triangle
        values[0] =  (16 * xi) +  (16 * eta) - 6;
        values[1] = -16 * xi;
        values[2] = 0;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 0;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = -2 +  (16 * eta);
        values[11] = - (16 * xi) -  (32 * eta) + 8;
        values[12] = 16 * xi;
        values[13] = 0;
        values[14] = 0;
      }
    }
  }
}

// values of the derivatives in xi-xi  direction
static void C_T_RR2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  if(xi>=0.5)
  {
    // right triangle
    values[0] = 0;
    values[1] = 0;
    values[2] = 16;
    values[3] = -32;
    values[4] = 16;
    values[5] = 0;
    values[6] = 0;
    values[7] = 0;
    values[8] = 0;
    values[9] = 0;
    values[10] = 0;
    values[11] = 0;
    values[12] = 0;
    values[13] = 0;
    values[14] = 0;
  }
  else
  {
    if(eta>=0.5)
    {
      // upper triangle
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 16;
      values[7] = 0;
      values[8] = 0;
      values[9] = 0;
      values[10] = 16;
      values[11] = 0;
      values[12] = 0;
      values[13] = 0;
      values[14] = -32;
    }
    else
    {
      if(xi+eta>=0.5)
      {
        // middle triangle
        values[0] = 0;
        values[1] = 0;
        values[2] = 0;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 16;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 16;
        values[11] = 0;
        values[12] = 0;
        values[13] = 0;
        values[14] = -32;
      }
      else
      {
        // lower left triangle
        values[0] = 16;
        values[1] = -32;
        values[2] = 16;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 0;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 0;
        values[11] = 0;
        values[12] = 0;
        values[13] = 0;
        values[14] = 0;
      }
    }
  }
}

// values of the derivatives in xi-eta direction
static void C_T_RR2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  if(xi>=0.5)
  {
    // right triangle
    values[0] = 0;
    values[1] = 0;
    values[2] = 16;
    values[3] = -16;
    values[4] = 0;
    values[5] = 16;
    values[6] = 0;
    values[7] = 0;
    values[8] = 0;
    values[9] = 0;
    values[10] = 0;
    values[11] = 0;
    values[12] = 0;
    values[13] = -16;
    values[14] = 0;
  }
  else
  {
    if(eta>=0.5)
    {
      // upper triangle
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 0;
      values[7] = 16;
      values[8] = 0;
      values[9] = -16;
      values[10] = 16;
      values[11] = 0;
      values[12] = 0;
      values[13] = 0;
      values[14] = -16;
    }
    else
    {
      if(xi+eta>=0.5)
      {
        // middle triangle
        values[0] = 0;
        values[1] = 0;
        values[2] = 0;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 16;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 0;
        values[11] = 0;
        values[12] = 16;
        values[13] = -16;
        values[14] = -16;
      }
      else
      {
        // lower left triangle
        values[0] = 16;
        values[1] = -16;
        values[2] = 0;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 0;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 0;
        values[11] = -16;
        values[12] = 16;
        values[13] = 0;
        values[14] = 0;
      }
    }
  }
}

// values of the derivatives in eta-eta direction
static void C_T_RR2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  if(xi>=0.5)
  {
    // right triangle
    values[0] = 0;
    values[1] = 0;
    values[2] = 16;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] = 16;
    values[7] = 0;
    values[8] = 0;
    values[9] = 0;
    values[10] = 0;
    values[11] = 0;
    values[12] = 0;
    values[13] = -32;
    values[14] = 0;
  }
  else
  {
    if(eta>=0.5)
    {
      // upper triangle
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 0;
      values[7] = 0;
      values[8] = 16;
      values[9] = -32;
      values[10] = 16;
      values[11] = 0;
      values[12] = 0;
      values[13] = 0;
      values[14] = 0;
    }
    else
    {
      if(xi+eta>=0.5)
      {
        // middle triangle
        values[0] = 0;
        values[1] = 0;
        values[2] = 16;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 16;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 0;
        values[11] = 0;
        values[12] = 0;
        values[13] = -32;
        values[14] = 0;
      }
      else
      {
        // lower left triangle
        values[0] = 16;
        values[1] = 0;
        values[2] = 0;
        values[3] = 0;
        values[4] = 0;
        values[5] = 0;
        values[6] = 0;
        values[7] = 0;
        values[8] = 0;
        values[9] = 0;
        values[10] = 16;
        values[11] = -32;
        values[12] = 0;
        values[13] = 0;
        values[14] = 0;
      }
    }
  }
}

// ***********************************************************************

TBaseFunct2D *BF_C_T_RR2_2D_Obj = new TBaseFunct2D
        (15, BF_C_T_B2_2D, BFUnitTriangle, 
         C_T_RR2_2D_Funct, C_T_RR2_2D_DeriveXi,
         C_T_RR2_2D_DeriveEta, C_T_RR2_2D_DeriveXiXi,
         C_T_RR2_2D_DeriveXiEta, C_T_RR2_2D_DeriveEtaEta, 2, 2,
         0, NULL);
