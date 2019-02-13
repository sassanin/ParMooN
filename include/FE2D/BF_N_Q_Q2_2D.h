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
// Q2 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_Q_Q2_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t8, t10, t11, t12, t15, t20, t26, t32;
  double t7, t9, t16, t17, t18, t22;
  double t6, t19;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2 // GM
      t1 = 3.0/16.0*eta;
      t2 = eta*eta;
      t3 = 3.0/4.0*t2;
      t4 = xi*xi;
      t5 = t4-t2;
      t7 = 15.0/32.0*eta*t5;
      t9 = 3.0/16.0*xi;
      t10 = 3.0/4.0*t4;
      t12 = 15.0/32.0*xi*t5;
      t16 = 5.0/16.0*xi;
      t17 = xi*eta;
      t18 = t17/4.0;
      t20 = 5.0/8.0*t17*t5;
      t22 = 5.0/16.0*eta;
  
      values[0] = -1.0/4.0-t1+t3+t7;
      values[1] = -1.0/4.0+t9+t10+t12;
      values[2] = -1.0/4.0+t1+t3-t7;
      values[3] = -1.0/4.0-t9+t10-t12;
      values[4] = t16-t18-t12+t20;
      values[5] = t22+t18+t7+t20;
      values[6] = -t16-t18+t12+t20;
      values[7] = -t22+t18-t7+t20;
      values[8] = 2.0-3.0/2.0*t4-3.0/2.0*t2;
    break;
  
    case 3: // G. Matthies, only one difference, family 3 // FS
      t1 = 3.0/4.0*eta;
      t2 = eta*eta;
      t3 = 3.0/4.0*t2;
      t5 = 5.0/4.0*t2*eta;
      t7 = 3.0/4.0*xi;
      t8 = xi*xi;
      t9 = 3.0/4.0*t8;
      t11 = 5.0/4.0*t8*xi;
      t15 = 5.0/4.0*xi;
      t16 = xi*eta;
      t17 = t16/4.0;
      t20 = 5.0/8.0*t16*(t8-t2);
      t22 = 5.0/4.0*eta;
  
      values[0] = -1.0/4.0+t1+t3-t5;
      values[1] = -1.0/4.0-t7+t9+t11;
      values[2] = -1.0/4.0-t1+t3+t5;
      values[3] = -1.0/4.0+t7+t9-t11;
      values[4] = t15-t17-t11+t20;
      values[5] = t22+t17-t5+t20;
      values[6] = -t15-t17+t11+t20;
      values[7] = -t22+t17+t5+t20;
      values[8] = 2.0-3.0/2.0*t8-3.0/2.0*t2;
    break;
  
    case 4: // Apel, Matthies, anisotrop, family 4 // Apel/Matthies
      cout << "o";
      t1 = 3.0/4.0*eta;
      t2 = eta*eta;
      t3 = 3.0/4.0*t2;
      t4 = xi*xi;
      t6 = 3.0/4.0*t4*eta;
      t8 = 3.0/4.0*xi;
      t9 = 3.0/4.0*t4;
      t10 = t4*xi;
      t11 = 5.0/4.0*t10;
      t15 = xi*eta;
      t16 = t10*eta;
      t18 = eta/4.0;
      t19 = 3.0/4.0*t15;
      t20 = 5.0/4.0*t16;
  
      values[0] = -1.0/4.0-t1+t3+t6;
      values[1] = -1.0/4.0-t8+t9+t11;
      values[2] = -1.0/4.0+t1+t3-t6;
      values[3] = -1.0/4.0+t8+t9-t11;
      values[4] = 5.0/4.0*xi-5.0/4.0*t15-5.0/4.0*t10+5.0/4.0*t16;
      values[5] = -t18-t19+t6+t20;
      values[6] = -5.0/4.0*xi-5.0/4.0*t15+5.0/4.0*t10+5.0/4.0*t16;
      values[7] = t18-t19-t6+t20;
      values[8] = 2.0-3.0/2.0*t4-3.0/2.0*t2;
    break;
  
    case 1: // Hennart, Jaffr'e, Roberts, 1988 // Hennart et al.
      t1 = eta/2.0;
      t2 = eta*eta;
      t3 = 3.0/4.0*t2;
      t4 = xi*xi;
      t5 = 3.0/2.0*t4;
      t8 = (-1.0/2.0+t5)*eta/2.0;
      t10 = xi/2.0;
      t11 = 3.0/4.0*t4;
      t12 = 3.0/2.0*t2;
      t15 = xi*(-1.0/2.0+t12)/2.0;
      t20 = xi*eta/4.0;
      t26 = (5.0/2.0*t4*xi-3.0/2.0*xi)*eta/4.0;
      t32 = xi*(5.0/2.0*t2*eta-3.0/2.0*eta)/4.0;
  
      values[0] = -1.0/4.0-t1+t3+t8;
      values[1] = -1.0/4.0+t10+t11-t15;
      values[2] = -1.0/4.0+t1+t3-t8;
      values[3] = -1.0/4.0-t10+t11+t15;
      values[4] = -t20+t15+t26-t32;
      values[5] = t20+t8+t26-t32;
      values[6] = -t20-t15+t26-t32;
      values[7] = t20-t8+t26-t32;
      values[8] = 2.0-t5-t12;
    break;
  } // end switch
}

// values of the derivatives in xi direction
static void N_Q_Q2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2, t3, t4, t5, t8, t9, t13, t15, t17;
  double t1, t6, t7, t10, t12, t14;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t2 = 15.0/16.0*xi*eta;
      t3 = 3.0/2.0*xi;
      t4 = xi*xi;
      t5 = 45.0/32.0*t4;
      t6 = eta*eta;
      t7 = 15.0/32.0*t6;
      t10 = eta/4.0;
      t13 = 5.0/8.0*eta*(t4-t6);
      t15 = 5.0/4.0*t4*eta;
    
      values[0] = t2;
      values[1] = 3.0/16.0+t3+t5-t7;
      values[2] = -t2;
      values[3] = -3.0/16.0+t3-t5+t7;
      values[4] = 5.0/16.0-t10-t5+t7+t13+t15;
      values[5] = t10+t2+t13+t15;
      values[6] = -5.0/16.0-t10+t5-t7+t13+t15;
      values[7] = t10-t2+t13+t15;
      values[8] = -3.0*xi;
    break;
    
    case 3: // G. Matthies, only one difference, family 3
      t1 = 3.0/2.0*xi;
      t2 = xi*xi;
      t3 = 15.0/4.0*t2;
      t6 = eta/4.0;
      t7 = eta*eta;
      t10 = 5.0/8.0*eta*(t2-t7);
      t12 = 5.0/4.0*t2*eta;
      t14 = t6+t10+t12;
    
      values[0] = 0.0;
      values[1] = -3.0/4.0+t1+t3;
      values[2] = 0.0;
      values[3] = 3.0/4.0+t1-t3;
      values[4] = 5.0/4.0-t6-t3+t10+t12;
      values[5] = t14;
      values[6] = -5.0/4.0-t6+t3+t10+t12;
      values[7] = t14;
      values[8] = -3.0*xi;
    break;
    
    case 4: // Apel, Matthies, anisotrop, family 4
      t2 = 3.0/2.0*xi*eta;
      t3 = 3.0/2.0*xi;
      t4 = xi*xi;
      t5 = 15.0/4.0*t4;
      t8 = 5.0/4.0*eta;
      t10 = 15.0/4.0*t4*eta;
      t12 = 3.0/4.0*eta;
    
      values[0] = t2;
      values[1] = -3.0/4.0+t3+t5;
      values[2] = -t2;
      values[3] = 3.0/4.0+t3-t5;
      values[4] = 5.0/4.0-t8-t5+t10;
      values[5] = -t12+t2+t10;
      values[6] = -5.0/4.0-t8+t5+t10;
      values[7] = -t12-t2+t10;
      values[8] = -3.0*xi;
    break;
    
    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t2 = 3.0/2.0*xi*eta;
      t3 = 3.0/2.0*xi;
      t4 = eta*eta;
      t5 = 3.0/4.0*t4;
      t8 = eta/8.0;
      t9 = xi*xi;
      t13 = (15.0/2.0*t9-3.0/2.0)*eta/4.0;
      t15 = 5.0/8.0*t4*eta;
      t17 = 5.0/8.0*eta;
    
      values[0] = t2;
      values[1] = 3.0/4.0+t3-t5;
      values[2] = -t2;
      values[3] = -3.0/4.0+t3+t5;
      values[4] = t8-1.0/4.0+t5+t13-t15;
      values[5] = t17+t2+t13-t15;
      values[6] = t8+1.0/4.0-t5+t13-t15;
      values[7] = t17-t2+t13-t15;
      values[8] = -3.0*xi;
    break;
  }
}

// values of the derivatives in eta direction
static void N_Q_Q2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t6, t8, t10, t11, t15, t17;
  double t4, t5, t7, t12, t13;
  double t9, t14, t16, t19, t20, t21, t22;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 3.0/2.0*eta;
      t2 = xi*xi;
      t3 = 15.0/32.0*t2;
      t4 = eta*eta;
      t5 = 45.0/32.0*t4;
      t8 = 15.0/16.0*xi*eta;
      t10 = xi/4.0;
      t13 = 5.0/8.0*xi*(t2-t4);
      t15 = 5.0/4.0*xi*t4;
    
      values[0] = -3.0/16.0+t1+t3-t5;
      values[1] = -t8;
      values[2] = 3.0/16.0+t1-t3+t5;
      values[3] = t8;
      values[4] = -t10+t8+t13-t15;
      values[5] = 5.0/16.0+t10+t3-t5+t13-t15;
      values[6] = -t10-t8+t13-t15;
      values[7] = -5.0/16.0+t10-t3+t5+t13-t15;
      values[8] = -3.0*eta;
    break;
    
    case 3: // G. Matthies, only one difference, family 3
      t1 = 3.0/2.0*eta;
      t2 = eta*eta;
      t3 = 15.0/4.0*t2;
      t6 = xi/4.0;
      t7 = xi*xi;
      t10 = 5.0/8.0*xi*(t7-t2);
      t12 = 5.0/4.0*xi*t2;
      t13 = -t6+t10-t12;
    
      values[0] = 3.0/4.0+t1-t3;
      values[1] = 0.0;
      values[2] = -3.0/4.0+t1+t3;
      values[3] = 0.0;
      values[4] = t13;
      values[5] = 5.0/4.0+t6-t3+t10-t12;
      values[6] = t13;
      values[7] = -5.0/4.0+t6+t3+t10-t12;
      values[8] = -3.0*eta;
    break;
    
    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = 3.0/2.0*eta;
      t2 = xi*xi;
      t3 = 3.0/4.0*t2;
      t6 = t2*xi;
      t7 = -xi+t6;
      t8 = 3.0/4.0*xi;
      t9 = 5.0/4.0*t6;
    
      values[0] = -3.0/4.0+t1+t3;
      values[1] = 0.0;
      values[2] = 3.0/4.0+t1-t3;
      values[3] = 0.0;
      values[4] = 5.0/4.0*t7;
      values[5] = -1.0/4.0-t8+t3+t9;
      values[6] = 5.0/4.0*t7;
      values[7] = 1.0/4.0-t8-t3+t9;
      values[8] = -3.0*eta;
    break;
    
    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = 3.0/2.0*eta;
      t2 = xi*xi;
      t3 = 3.0/4.0*t2;
      t6 = 3.0/2.0*xi*eta;
      t8 = 5.0/8.0*xi;
      t10 = 5.0/8.0*t2*xi;
      t11 = eta*eta;
      t15 = xi*(15.0/2.0*t11-3.0/2.0)/4.0;
      t17 = xi/8.0;
    
      values[0] = -3.0/4.0+t1+t3;
      values[1] = -t6;
      values[2] = 3.0/4.0+t1-t3;
      values[3] = t6;
      values[4] = -t8+t6+t10-t15;
      values[5] = -t17-1.0/4.0+t3+t10-t15;
      values[6] = -t8-t6+t10-t15;
      values[7] = -t17+1.0/4.0-t3+t10-t15;
      values[8] = -3.0*eta;
    break;
  }
}

// values of derivatives in xi-xi direction
static void N_Q_Q2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t3;
  double t2, t5, t6;
  double t7;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 15.0/16.0*eta;
      t2 = 45.0/16.0*xi;
      t6 = 15.0/4.0*xi*eta;
    
      values[0] = t1;
      values[1] = 3.0/2.0+t2;
      values[2] = -t1;
      values[3] = 3.0/2.0-t2;
      values[4] = -t2+t6;
      values[5] = t1+t6;
      values[6] = t2+t6;
      values[7] = -t1+t6;
      values[8] = -3.0;
    break;
    
    case 3: // G. Matthies, only one difference, family 3
      t1 = 15.0/2.0*xi;
      t5 = 15.0/4.0*xi*eta;
    
      values[0] = 0.0;
      values[1] = 3.0/2.0+t1;
      values[2] = 0.0;
      values[3] = 3.0/2.0-t1;
      values[4] = -t1+t5;
      values[5] = t5;
      values[6] = t1+t5;
      values[7] = t5;
      values[8] = -3.0;
    break;
    
    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = 3.0/2.0*eta;
      t2 = 15.0/2.0*xi;
      t5 = xi*eta;
      t7 = 15.0/2.0*t5;
    
      values[0] = t1;
      values[1] = 3.0/2.0+t2;
      values[2] = -t1;
      values[3] = 3.0/2.0-t2;
      values[4] = -15.0/2.0*xi+15.0/2.0*t5;
      values[5] = t1+t7;
      values[6] = 15.0/2.0*xi+15.0/2.0*t5;
      values[7] = -t1+t7;
      values[8] = -3.0;
    break;
    
    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = 3.0/2.0*eta;
      t3 = 15.0/4.0*xi*eta;
    
      values[0] = t1;
      values[1] = 3.0/2.0;
      values[2] = -t1;
      values[3] = 3.0/2.0;
      values[4] = t3;
      values[5] = t1+t3;
      values[6] = t3;
      values[7] = -t1+t3;
      values[8] = -3.0;
    break;
  }
}

// values of derivatives in eta-eta direction
static void N_Q_Q2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 15.0/16.0*xi;
      t2 = 15.0/16.0*eta;
      t3 = xi*xi;
      t4 = 15.0/8.0*t3;
      t5 = eta*eta;
      t6 = 15.0/8.0*t5;
    
      values[0] = t1;
      values[1] = -t2;
      values[2] = -t1;
      values[3] = t2;
      values[4] = -1.0/4.0+t2+t4-t6;
      values[5] = 1.0/4.0+t1+t4-t6;
      values[6] = -1.0/4.0-t2+t4-t6;
      values[7] = 1.0/4.0-t1+t4-t6;
      values[8] = 0.0;
    break;
    
    case 3: // G. Matthies, only one difference, family 3
      t1 = xi*xi;
      t2 = 15.0/8.0*t1;
      t3 = eta*eta;
      t4 = 15.0/8.0*t3;
      t5 = -1.0/4.0+t2-t4;
      t6 = 1.0/4.0+t2-t4;
    
      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 0.0;
      values[3] = 0.0;
      values[4] = t5;
      values[5] = t6;
      values[6] = t5;
      values[7] = t6;
      values[8] = 0.0;
    break;
    
    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = 3.0/2.0*xi;
      t2 = xi*xi;
      t3 = 15.0/4.0*t2;
      t4 = -5.0/4.0+t3;
    
      values[0] = t1;
      values[1] = 0.0;
      values[2] = -t1;
      values[3] = 0.0;
      values[4] = t4;
      values[5] = -3.0/4.0+t1+t3;
      values[6] = t4;
      values[7] = -3.0/4.0-t1+t3;
      values[8] = 0.0;
    break;
    
    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = 3.0/2.0*xi;
      t2 = 3.0/2.0*eta;
      t3 = xi*xi;
      t4 = 15.0/8.0*t3;
      t5 = eta*eta;
      t6 = 15.0/8.0*t5;
    
      values[0] = t1;
      values[1] = -t2;
      values[2] = -t1;
      values[3] = t2;
      values[4] = -1.0/4.0+t2+t4-t6;
      values[5] = 1.0/4.0+t1+t4-t6;
      values[6] = -1.0/4.0-t2+t4-t6;
      values[7] = 1.0/4.0-t1+t4-t6;
      values[8] = 0.0;
    break;
  }
}

// values of derivatives in xi-eta direction
static void N_Q_Q2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t3;
  double t5, t6;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 45.0/16.0*eta;
      t3 = 15.0/16.0*xi;
      t6 = 15.0/4.0*xi*eta;
    
      values[0] = 3.0/2.0-t1;
      values[1] = -t3;
      values[2] = 3.0/2.0+t1;
      values[3] = t3;
      values[4] = t3-t6;
      values[5] = -t1-t6;
      values[6] = -t3-t6;
      values[7] = t1-t6;
      values[8] = -3.0;
    break;
    
    case 3: // G. Matthies, only one difference, family 3
      t1 = 15.0/2.0*eta;
      t5 = 15.0/4.0*xi*eta;
    
      values[0] = 3.0/2.0-t1;
      values[1] = 0.0;
      values[2] = 3.0/2.0+t1;
      values[3] = 0.0;
      values[4] = -t5;
      values[5] = -t1-t5;
      values[6] = -t5;
      values[7] = t1-t5;
      values[8] = -3.0;
    break;
    
    case 4: // Apel, Matthies, anisotrop, family 4
      values[0] = 3.0/2.0;
      values[1] = 0.0;
      values[2] = 3.0/2.0;
      values[3] = 0.0;
      values[4] = 0.0;
      values[5] = 0.0;
      values[6] = 0.0;
      values[7] = 0.0;
      values[8] = -3.0;
    break;
    
    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = 3.0/2.0*xi;
      t3 = 15.0/4.0*xi*eta;
    
      values[0] = 3.0/2.0;
      values[1] = -t1;
      values[2] = 3.0/2.0;
      values[3] = t1;
      values[4] = t1-t3;
      values[5] = -t3;
      values[6] = -t1-t3;
      values[7] = -t3;
      values[8] = -3.0;
    break;
  }
}

static int N_Q_Q2_2D_ChangeJ0[1] = { 4 };
static int N_Q_Q2_2D_ChangeJ1[1] = { 5 };
static int N_Q_Q2_2D_ChangeJ2[1] = { 6 };
static int N_Q_Q2_2D_ChangeJ3[1] = { 7 };

static int *N_Q_Q2_2D_Change[4] = { N_Q_Q2_2D_ChangeJ0, N_Q_Q2_2D_ChangeJ1,
                                    N_Q_Q2_2D_ChangeJ2, N_Q_Q2_2D_ChangeJ3 };

// ***********************************************************************

TBaseFunct2D *BF_N_Q_Q2_2D_Obj = new TBaseFunct2D
        (9, BF_N_Q_Q2_2D, BFUnitSquare, 
         N_Q_Q2_2D_Funct, N_Q_Q2_2D_DeriveXi,
         N_Q_Q2_2D_DeriveEta, N_Q_Q2_2D_DeriveXiXi,
         N_Q_Q2_2D_DeriveXiEta, N_Q_Q2_2D_DeriveEtaEta, 3, 2,
         1, N_Q_Q2_2D_Change);
