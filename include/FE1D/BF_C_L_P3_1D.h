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
// P3 element, conforming, 1D
// ***********************************************************************

// base function values
static void C_L_P3_1D_Funct(double xi, double *values)
{
  values[0]=-0.0625*(3*xi+1)*(3*xi-1)*(xi-1);
  values[1]= 0.5625*(xi+1)*(3*xi-1)*(xi-1);
  values[2]=-0.5625*(xi+1)*(3*xi+1)*(xi-1);
  values[3]= 0.0625*(xi+1)*(3*xi+1)*(3*xi-1);
}

// values of the derivatives in xi direction
static void C_L_P3_1D_DeriveXi(double xi, double *values)
{
  values[0]=0.0625*(xi*(-27*xi+18)+1);
  values[1]=0.0625*(xi*(81*xi-18)-27);
  values[2]=0.0625*(xi*(-81*xi-18)+27);
  values[3]=0.0625*(xi*(27*xi+18)+1);
}

// values of the derivatives in xi-xi  direction
static void C_L_P3_1D_DeriveXiXi(double xi, double *values)
{
  values[0]=0.125*(-27*xi+9);
  values[1]=0.125*( 81*xi-9);
  values[2]=0.125*(-81*xi-9);
  values[3]=0.125*( 27*xi+9);
}

// ***********************************************************************
// TBaseFunct1D *BF_C_L_P3_1D_Obj = new TBaseFunct1D
//         (4, BF_C_L_P3_1D, C_L_P3_1D_Funct, C_L_P3_1D_DeriveXi, 
//          C_L_P3_1D_DeriveXiXi);
TBaseFunct1D *BF_C_L_P3_1D_Obj = new TBaseFunct1D
        (4, BF_C_L_P3_1D, C_L_P3_1D_Funct, C_L_P3_1D_DeriveXi, 
         C_L_P3_1D_DeriveXiXi, 3, 3);
