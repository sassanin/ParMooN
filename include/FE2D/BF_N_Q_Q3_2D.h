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
// Q3 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_Q_Q3_2D_Funct(double xi, double eta, double *values)
{
  double  t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55;
  double t58, t59, t63, t65;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 3.0/4.0*eta;
      t2 = xi*xi;
      t3 = 7.0/8.0*t2;
      t4 = eta*eta;
      t5 = t4*eta;
      t6 = 5.0/4.0*t5;
      t7 = t4*t4;
      t8 = 35.0/24.0*t7;
      t10 = -1.0/2.0+3.0/2.0*t2;
      t12 = -1.0/2.0+3.0/2.0*t4;
      t13 = t10*t12;
      t14 = t13/6.0;
      t15 = t2*t2;
      t16 = 35.0/48.0*t15;
      t18 = 3.0/4.0*xi;
      t19 = 7.0/8.0*t4;
      t20 = t2*xi;
      t21 = 5.0/4.0*t20;
      t22 = 35.0/48.0*t7;
      t23 = 35.0/24.0*t15;
      t27 = xi/8.0;
      t29 = xi*eta/4.0;
      t30 = xi*t4;
      t31 = 3.0/8.0*t30;
      t32 = 5.0/2.0*t5;
      t34 = t32-3.0/2.0*eta;
      t36 = xi*t34/4.0;
      t37 = 5.0/2.0*t20;
      t39 = t37-3.0/2.0*xi;
      t41 = t39*eta/4.0;
      t45 = xi*(3.0/8.0+35.0/8.0*t7-15.0/4.0*t4);
      t46 = t45/4.0;
      t47 = t39*t12;
      t48 = t47/4.0;
      t50 = eta/8.0;
      t51 = t2*eta;
      t52 = 3.0/8.0*t51;
      t53 = t10*t34;
      t54 = t53/4.0;
      t58 = (3.0/8.0+35.0/8.0*t15-15.0/4.0*t2)*eta;
      t59 = t58/4.0;
      t63 = 7.0/4.0*t2;
      t65 = 7.0/4.0*t4;

      values[0] = -5.0/48.0+t1+t3-t4-t6+t8-t14-t16;
      values[1] = -5.0/48.0-t18-t2+t19+t21-t22-t14+t23;
      values[2] = -5.0/48.0-t1+t3-t4+t6+t8-t14-t16;
      values[3] = -5.0/48.0+t18-t2+t19-t21-t22-t14+t23;
      values[4] = -t27-t29+t31-t36+t41+t46-t48;
      values[5] = -t50+t29+t52-t36+t41-t54+t59;
      values[6] = t27-t29-t31-t36+t41-t46+t48;
      values[7] = t50+t29-t52-t36+t41+t54-t59;
      values[8] = -7.0/48.0+t50+t63-t19-t52+t22+t14-t23-t54+t59;
      values[9] = -7.0/48.0-t27-t3+t65+t31-t8+t14+t16-t46+t48;
      values[10] = -7.0/48.0-t50+t63-t19+t52+t22+t14-t23+t54-t59;
      values[11] = -7.0/48.0+t27-t3+t65-t31-t8+t14+t16+t46-t48;
      values[12] = 17.0/12.0+t2/4.0+t4/4.0-t8+2.0/3.0*t13-t23;
      values[13] = 11.0/4.0*xi-t37-3.0/4.0*t30-t45/2.0+t47/2.0;
      values[14] = 11.0/4.0*eta-3.0/4.0*t51-t32+t53/2.0-t58/2.0;
    break;

    case 3: // G. Matthies, only one difference, family 3
      t1 = 3.0/4.0*eta;
      t2 = eta*eta;
      t3 = t2*eta;
      t4 = 5.0/4.0*t3;
      t5 = t2*t2;
      t6 = 35.0/16.0*t5;
      t7 = 15.0/8.0*t2;
      t9 = 3.0/4.0*xi;
      t10 = xi*xi;
      t11 = t10*xi;
      t12 = 5.0/4.0*t11;
      t13 = t10*t10;
      t14 = 35.0/16.0*t13;
      t15 = 15.0/8.0*t10;
      t20 = xi*eta/4.0;
      t21 = 5.0/2.0*t3;
      t25 = xi*(t21-3.0/2.0*eta)/4.0;
      t26 = 5.0/2.0*t11;
      t30 = (t26-3.0/2.0*xi)*eta/4.0;
      t31 = 35.0/8.0*t5;
      t32 = 15.0/4.0*t2;
      t34 = xi*(3.0/8.0+t31-t32);
      t35 = t34/2.0;
      t37 = 35.0/8.0*t13;
      t38 = 15.0/4.0*t10;
      t40 = (3.0/8.0+t37-t38)*eta;
      t41 = t40/2.0;
      t45 = eta/4.0;
      t46 = 21.0/8.0*t10;
      t48 = 3.0/4.0*t10*eta;
      t50 = xi/4.0;
      t51 = 21.0/8.0*t2;
      t53 = 3.0/4.0*xi*t2;

      values[0] = 3.0/16.0+t1-t4+t6-t7;
      values[1] = 3.0/16.0-t9+t12+t14-t15;
      values[2] = 3.0/16.0-t1+t4+t6-t7;
      values[3] = 3.0/16.0+t9-t12+t14-t15;
      values[4] = -t20-t25+t30+t35;
      values[5] = t20-t25+t30+t41;
      values[6] = -t20-t25+t30-t35;
      values[7] = t20-t25+t30-t41;
      values[8] = -7.0/16.0+t45+t46-t48-t14+t41;
      values[9] = -7.0/16.0-t50+t51+t53-t6-t35;
      values[10] = -7.0/16.0-t45+t46+t48-t14-t41;
      values[11] = -7.0/16.0+t50+t51-t53-t6+t35;
      values[12] = 1.0/4.0-t31+t32-t37+t38;
      values[13] = 5.0/2.0*xi-t26-t34;
      values[14] = 5.0/2.0*eta-t21-t40;
    break;

    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = 3.0/4.0*eta;
      t2 = xi*xi;
      t3 = 21.0/8.0*t2;
      t4 = eta*eta;
      t5 = 3.0/4.0*t4;
      t6 = t4*eta;
      t7 = 5.0/4.0*t6;
      t8 = 3.0/2.0*t2;
      t10 = 3.0/2.0*t4;
      t11 = -1.0/2.0+t10;
      t12 = (-1.0/2.0+t8)*t11;
      t13 = t12/2.0;
      t14 = t2*t2;
      t15 = 35.0/16.0*t14;
      t17 = 3.0/4.0*xi;
      t18 = t2*xi;
      t19 = 5.0/4.0*t18;
      t20 = 15.0/8.0*t2;
      t24 = xi/4.0;
      t26 = xi*eta/2.0;
      t27 = xi*t4;
      t28 = 3.0/4.0*t27;
      t29 = 5.0/2.0*t18;
      t31 = t29-3.0/2.0*xi;
      t32 = t31*eta;
      t33 = t32/2.0;
      t34 = t31*t11;
      t35 = t34/2.0;
      t40 = (3.0/8.0+35.0/8.0*t14-15.0/4.0*t2)*eta;
      t44 = eta/4.0;
      t46 = 3.0/4.0*t2*eta;
      t47 = t40/2.0;

      values[0] = -11.0/16.0+t1+t3+t5-t7-t13-t15;
      values[1] = 3.0/16.0-t17+t19+t15-t20;
      values[2] = -11.0/16.0-t1+t3+t5+t7-t13-t15;
      values[3] = 3.0/16.0+t17-t19+t15-t20;
      values[4] = -t24-t26+t28+t33-t35;
      values[5] = t32/2.0+t40/2.0;
      values[6] = t24-t26-t28+t33+t35;
      values[7] = t32/2.0-t40/2.0;
      values[8] = -7.0/16.0+t44+t3-t46-t15+t47;
      values[9] = 7.0/16.0-t3+t13+t15+t35;
      values[10] = -7.0/16.0-t44+t3+t46-t15-t47;
      values[11] = 7.0/16.0-t3+t13+t15-t35;
      values[12] = 2.0-t8-t10+t12;
      values[13] = 3.0*xi-t29-3.0/2.0*t27+t34;
      values[14] = 5.0/2.0*eta-5.0/2.0*t6-t40;
    break;

    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = 3.0/4.0*eta;
      t2 = eta*eta;
      t3 = 3.0/4.0*t2;
      t4 = t2*eta;
      t5 = 5.0/4.0*t4;
      t6 = xi*xi;
      t7 = 3.0/2.0*t6;
      t8 = -1.0/2.0+t7;
      t9 = 3.0/2.0*t2;
      t10 = -1.0/2.0+t9;
      t11 = t8*t10;
      t12 = t11/4.0;
      t13 = t6*t6;
      t18 = (3.0/8.0+35.0/8.0*t13-15.0/4.0*t6)*t10/4.0;
      t19 = t2*t2;
      t24 = t8*(3.0/8.0+35.0/8.0*t19-15.0/4.0*t2)/4.0;
      t26 = 3.0/4.0*xi;
      t27 = 3.0/4.0*t6;
      t28 = t6*xi;
      t29 = 5.0/4.0*t28;
      t33 = xi/4.0;
      t35 = xi*eta/4.0;
      t36 = xi*t2;
      t37 = 3.0/4.0*t36;
      t38 = 5.0/2.0*t28;
      t40 = t38-3.0/2.0*xi;
      t41 = t40*t10;
      t42 = t41/2.0;
      t44 = t40*eta/4.0;
      t45 = 5.0/2.0*t4;
      t47 = t45-3.0/2.0*eta;
      t49 = xi*t47/4.0;
      t51 = eta/4.0;
      t52 = t6*eta;
      t53 = 3.0/4.0*t52;
      t54 = t8*t47;
      t55 = t54/2.0;

      values[0] = -1.0/4.0+t1+t3-t5-t12-t18+t24;
      values[1] = -1.0/4.0-t26+t27+t29-t12+t18-t24;
      values[2] = -1.0/4.0-t1+t3+t5-t12-t18+t24;
      values[3] = -1.0/4.0+t26+t27-t29-t12+t18-t24;
      values[4] = -t33-t35+t37-t42+t44-t49;
      values[5] = -t51+t35+t53-t55+t44-t49;
      values[6] = t33-t35-t37+t42+t44-t49;
      values[7] = t51+t35-t53+t55+t44-t49;
      values[8] = t12-t55-t18+t24;
      values[9] = t12+t42+t18-t24;
      values[10] = t12+t55-t18+t24;
      values[11] = t12-t42+t18-t24;
      values[12] = 2.0-t7-t9+t11;
      values[13] = 3.0*xi-t38-3.0/2.0*t36+t41;
      values[14] = 3.0*eta-3.0/2.0*t52-t45+t54;
    break;

    default:
      OutPut("unknown NC_TYPE: " << TDatabase::ParamDB->NC_TYPE << endl);
  }
}

static void N_Q_Q3_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15, t16;
  double t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
  double t30, t31, t32, t33, t34, t35, t36, t37, t40, t41, t42, t45, t47;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 7.0/4.0*xi;
      t2 = eta*eta;
      t4 = -1.0/2.0+3.0/2.0*t2;
      t5 = xi*t4;
      t6 = t5/2.0;
      t7 = xi*xi;
      t8 = t7*xi;
      t9 = 35.0/12.0*t8;
      t10 = t1-t6-t9;
      t11 = 2.0*xi;
      t12 = 15.0/4.0*t7;
      t13 = 35.0/6.0*t8;
      t16 = eta/8.0;
      t17 = 9.0/16.0*t2;
      t18 = t2*eta;
      t19 = 5.0/8.0*t18;
      t20 = 15.0/2.0*t7;
      t21 = t20-3.0/2.0;
      t23 = t21*eta/4.0;
      t24 = t2*t2;
      t25 = 35.0/32.0*t24;
      t26 = t21*t4;
      t27 = t26/4.0;
      t29 = 5.0/8.0*eta;
      t30 = xi*eta;
      t31 = 3.0/4.0*t30;
      t35 = xi*(5.0/2.0*t18-3.0/2.0*eta);
      t36 = 3.0/4.0*t35;
      t40 = (35.0/2.0*t8-15.0/2.0*xi)*eta;
      t41 = t40/4.0;
      t45 = 7.0/2.0*xi;
      t47 = 21.0/16.0*t2;

      values[0] = t10;
      values[1] = -3.0/4.0-t11+t12-t6+t13;
      values[2] = t10;
      values[3] = 3.0/4.0-t11-t12-t6+t13;
      values[4] = -1.0/32.0+t16-t17-t19+t23+t25-t27;
      values[5] = t29+t31-t19+t23-t36+t41;
      values[6] = 1.0/32.0+t16+t17-t19+t23-t25+t27;
      values[7] = t29-t31-t19+t23+t36-t41;
      values[8] = t45-t31+t6-t13-t36+t41;
      values[9] = -7.0/32.0-t1+t47+t6+t9-t25+t27;
      values[10] = t45+t31+t6-t13+t36-t41;
      values[11] = 7.0/32.0-t1-t47+t6+t9+t25-t27;
      values[12] = xi/2.0+2.0*t5-t13;
      values[13] = 41.0/16.0-t20+9.0/8.0*t2-35.0/16.0*t24+t26/2.0;
      values[14] = -3.0/2.0*t30+3.0/2.0*t35-t40/2.0;
    break;

    case 3: // G. Matthies, only one difference, family 3
      t1 = xi*xi;
      t2 = 15.0/4.0*t1;
      t3 = t1*xi;
      t4 = 35.0/4.0*t3;
      t5 = 15.0/4.0*xi;
      t8 = eta/8.0;
      t9 = eta*eta;
      t11 = 5.0/8.0*t9*eta;
      t12 = 15.0/2.0*t1;
      t15 = (t12-3.0/2.0)*eta/4.0;
      t16 = t9*t9;
      t17 = 35.0/16.0*t16;
      t18 = 15.0/8.0*t9;
      t20 = 5.0/8.0*eta;
      t23 = 35.0/2.0*t3-15.0/2.0*xi;
      t24 = t23*eta;
      t25 = t24/2.0;
      t29 = 21.0/4.0*xi;
      t31 = 3.0/2.0*xi*eta;
      t34 = -7.0/16.0+21.0/8.0*t9-t17;

      values[0] = 0.0;
      values[1] = -3.0/4.0+t2+t4-t5;
      values[2] = 0.0;
      values[3] = 3.0/4.0-t2+t4-t5;
      values[4] = t8-t11+t15+3.0/16.0+t17-t18;
      values[5] = t20-t11+t15+t25;
      values[6] = t8-t11+t15-3.0/16.0-t17+t18;
      values[7] = t20-t11+t15-t25;
      values[8] = t29-t31-t4+t25;
      values[9] = t34;
      values[10] = t29+t31-t4-t25;
      values[11] = -t34;
      values[12] = -t23;
      values[13] = 17.0/8.0-t12-35.0/8.0*t16+15.0/4.0*t9;
      values[14] = -t24;
    break;

    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = 21.0/4.0*xi;
      t2 = eta*eta;
      t3 = 3.0/2.0*t2;
      t4 = -1.0/2.0+t3;
      t5 = xi*t4;
      t6 = 3.0/2.0*t5;
      t7 = xi*xi;
      t8 = t7*xi;
      t9 = 35.0/4.0*t8;
      t10 = t1-t6-t9;
      t11 = 15.0/4.0*t7;
      t12 = 15.0/4.0*xi;
      t15 = eta/2.0;
      t16 = 3.0/4.0*t2;
      t17 = 15.0/2.0*t7;
      t18 = t17-3.0/2.0;
      t19 = t18*eta;
      t20 = t19/2.0;
      t21 = t18*t4;
      t22 = t21/2.0;
      t27 = (35.0/2.0*t8-15.0/2.0*xi)*eta;
      t32 = 3.0/2.0*xi*eta;
      t33 = t27/2.0;

      values[0] = t10;
      values[1] = -3.0/4.0+t11+t9-t12;
      values[2] = t10;
      values[3] = 3.0/4.0-t11+t9-t12;
      values[4] = -1.0/4.0-t15+t16+t20-t22;
      values[5] = t19/2.0+t27/2.0;
      values[6] = 1.0/4.0-t15-t16+t20+t22;
      values[7] = t19/2.0-t27/2.0;
      values[8] = t1-t32-t9+t33;
      values[9] = -t1+t6+t9+t22;
      values[10] = t1+t32-t9-t33;
      values[11] = -t1+t6+t9-t22;
      values[12] = -3.0*xi+3.0*t5;
      values[13] = 3.0-t17-t3+t21;
      values[14] = -t27;
    break;

    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = eta*eta;
      t2 = 3.0/2.0*t1;
      t3 = -1.0/2.0+t2;
      t4 = xi*t3;
      t5 = 3.0/4.0*t4;
      t6 = xi*xi;
      t12 = (35.0/2.0*t6*xi-15.0/2.0*xi)*t3/4.0;
      t13 = t1*t1;
      t18 = 3.0/4.0*xi*(3.0/8.0+35.0/8.0*t13-15.0/4.0*t1);
      t19 = -t5-t12+t18;
      t20 = 3.0/2.0*xi;
      t21 = 15.0/4.0*t6;
      t24 = eta/8.0;
      t25 = 3.0/4.0*t1;
      t26 = 15.0/2.0*t6;
      t27 = t26-3.0/2.0;
      t28 = t27*t3;
      t29 = t28/2.0;
      t31 = t27*eta/4.0;
      t32 = t1*eta;
      t33 = 5.0/8.0*t32;
      t35 = 5.0/8.0*eta;
      t36 = xi*eta;
      t37 = 3.0/2.0*t36;
      t41 = xi*(5.0/2.0*t32-3.0/2.0*eta);
      t42 = 3.0/2.0*t41;

      values[0] = t19;
      values[1] = -3.0/4.0+t20+t21-t5+t12-t18;
      values[2] = t19;
      values[3] = 3.0/4.0+t20-t21-t5+t12-t18;
      values[4] = -1.0/4.0+t24+t25-t29+t31-t33;
      values[5] = t35+t37-t42+t31-t33;
      values[6] = 1.0/4.0+t24-t25+t29+t31-t33;
      values[7] = t35-t37+t42+t31-t33;
      values[8] = t5-t42-t12+t18;
      values[9] = t5+t29+t12-t18;
      values[10] = t5+t42-t12+t18;
      values[11] = t5-t29+t12-t18;
      values[12] = -3.0*xi+3.0*t4;
      values[13] = 3.0-t26-t2+t28;
      values[14] = -3.0*t36+3.0*t41;
    break;

    default:
      OutPut("unknown NC_TYPE: " << TDatabase::ParamDB->NC_TYPE << endl);
  }
}

static void N_Q_Q3_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t22, t23, t24, t25, t26, t27, t28, t29;
  double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42;
  double t45, t47;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 2.0*eta;
      t2 = eta*eta;
      t3 = 15.0/4.0*t2;
      t4 = t2*eta;
      t5 = 35.0/6.0*t4;
      t6 = xi*xi;
      t8 = -1.0/2.0+3.0/2.0*t6;
      t9 = t8*eta;
      t10 = t9/2.0;
      t12 = 7.0/4.0*eta;
      t13 = 35.0/12.0*t4;
      t14 = t12-t13-t10;
      t16 = 5.0/8.0*xi;
      t17 = xi*eta;
      t18 = 3.0/4.0*t17;
      t19 = 15.0/2.0*t2;
      t20 = t19-3.0/2.0;
      t22 = xi*t20/4.0;
      t23 = t6*xi;
      t24 = 5.0/8.0*t23;
      t28 = xi*(35.0/2.0*t4-15.0/2.0*eta);
      t29 = t28/4.0;
      t33 = (5.0/2.0*t23-3.0/2.0*xi)*eta;
      t34 = 3.0/4.0*t33;
      t36 = xi/8.0;
      t37 = 9.0/16.0*t6;
      t38 = t8*t20;
      t39 = t38/4.0;
      t40 = t6*t6;
      t41 = 35.0/32.0*t40;
      t45 = 21.0/16.0*t6;
      t47 = 7.0/2.0*eta;

      values[0] = 3.0/4.0-t1-t3+t5-t10;
      values[1] = t14;
      values[2] = -3.0/4.0-t1+t3+t5-t10;
      values[3] = t14;
      values[4] = -t16+t18-t22+t24+t29-t34;
      values[5] = -1.0/32.0-t36-t37-t22+t24-t39+t41;
      values[6] = -t16-t18-t22+t24-t29+t34;
      values[7] = 1.0/32.0-t36+t37-t22+t24+t39-t41;
      values[8] = 7.0/32.0-t12-t45+t13+t10-t39+t41;
      values[9] = t47+t18-t5+t10-t29+t34;
      values[10] = -7.0/32.0-t12+t45+t13+t10+t39-t41;
      values[11] = t47-t18-t5+t10+t29-t34;
      values[12] = eta/2.0-t5+2.0*t9;
      values[13] = -3.0/2.0*t17-t28/2.0+3.0/2.0*t33;
      values[14] = 41.0/16.0+9.0/8.0*t6-t19+t38/2.0-35.0/16.0*t40;
    break;

    case 3: // G. Matthies, only one difference, family 3
      t1 = eta*eta;
      t2 = 15.0/4.0*t1;
      t3 = t1*eta;
      t4 = 35.0/4.0*t3;
      t5 = 15.0/4.0*eta;
      t8 = 5.0/8.0*xi;
      t9 = 15.0/2.0*t1;
      t12 = xi*(t9-3.0/2.0)/4.0;
      t13 = xi*xi;
      t15 = 5.0/8.0*t13*xi;
      t18 = 35.0/2.0*t3-15.0/2.0*eta;
      t19 = xi*t18;
      t20 = t19/2.0;
      t22 = xi/8.0;
      t23 = t13*t13;
      t24 = 35.0/16.0*t23;
      t25 = 15.0/8.0*t13;
      t30 = 7.0/16.0-21.0/8.0*t13+t24;
      t31 = 21.0/4.0*eta;
      t33 = 3.0/2.0*xi*eta;

      values[0] = 3.0/4.0-t2+t4-t5;
      values[1] = 0.0;
      values[2] = -3.0/4.0+t2+t4-t5;
      values[3] = 0.0;
      values[4] = -t8-t12+t15+t20;
      values[5] = -t22-t12+t15+3.0/16.0+t24-t25;
      values[6] = -t8-t12+t15-t20;
      values[7] = -t22-t12+t15-3.0/16.0-t24+t25;
      values[8] = t30;
      values[9] = t31+t33-t4-t20;
      values[10] = -t30;
      values[11] = t31-t33-t4+t20;
      values[12] = -t18;
      values[13] = -t19;
      values[14] = 17.0/8.0-t9-35.0/8.0*t23+15.0/4.0*t13;
    break;

    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = 3.0/2.0*eta;
      t2 = eta*eta;
      t3 = 15.0/4.0*t2;
      t4 = xi*xi;
      t7 = (-1.0/2.0+3.0/2.0*t4)*eta;
      t8 = 3.0/2.0*t7;
      t11 = 5.0/4.0*xi;
      t12 = xi*eta;
      t13 = 3.0/2.0*t12;
      t14 = xi*t4;
      t15 = 5.0/4.0*t14;
      t19 = (5.0/2.0*t14-3.0/2.0*xi)*eta;
      t20 = 3.0/2.0*t19;
      t22 = 3.0/4.0*xi;
      t23 = t4*t4;
      t24 = 35.0/16.0*t23;
      t25 = 15.0/8.0*t4;
      t30 = 7.0/16.0-21.0/8.0*t4+t24;

      values[0] = 3.0/4.0+t1-t3-t8;
      values[1] = 0.0;
      values[2] = -3.0/4.0+t1+t3-t8;
      values[3] = 0.0;
      values[4] = -t11+t13+t15-t20;
      values[5] = 3.0/16.0-t22+t15+t24-t25;
      values[6] = -t11-t13+t15+t20;
      values[7] = t15-t22-3.0/16.0-t24+t25;
      values[8] = t30;
      values[9] = 3.0/2.0*t7+3.0/2.0*t19;
      values[10] = -t30;
      values[11] = 3.0/2.0*t7-3.0/2.0*t19;
      values[12] = -3.0*eta+3.0*t7;
      values[13] = -3.0*t12+3.0*t19;
      values[14] = 17.0/8.0-15.0/2.0*t2-35.0/8.0*t23+15.0/4.0*t4;
    break;

    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = 3.0/2.0*eta;
      t2 = eta*eta;
      t3 = 15.0/4.0*t2;
      t4 = xi*xi;
      t5 = 3.0/2.0*t4;
      t6 = -1.0/2.0+t5;
      t7 = t6*eta;
      t8 = 3.0/4.0*t7;
      t9 = t4*t4;
      t14 = 3.0/4.0*(3.0/8.0+35.0/8.0*t9-15.0/4.0*t4)*eta;
      t20 = t6*(35.0/2.0*t2*eta-15.0/2.0*eta)/4.0;
      t22 = -t8+t14-t20;
      t24 = 5.0/8.0*xi;
      t25 = xi*eta;
      t26 = 3.0/2.0*t25;
      t27 = t4*xi;
      t31 = (5.0/2.0*t27-3.0/2.0*xi)*eta;
      t32 = 3.0/2.0*t31;
      t33 = 5.0/8.0*t27;
      t34 = 15.0/2.0*t2;
      t35 = t34-3.0/2.0;
      t37 = xi*t35/4.0;
      t39 = xi/8.0;
      t40 = 3.0/4.0*t4;
      t41 = t6*t35;
      t42 = t41/2.0;

      values[0] = 3.0/4.0+t1-t3-t8-t14+t20;
      values[1] = t22;
      values[2] = -3.0/4.0+t1+t3-t8-t14+t20;
      values[3] = t22;
      values[4] = -t24+t26-t32+t33-t37;
      values[5] = -1.0/4.0-t39+t40-t42+t33-t37;
      values[6] = -t24-t26+t32+t33-t37;
      values[7] = 1.0/4.0-t39-t40+t42+t33-t37;
      values[8] = t8-t42-t14+t20;
      values[9] = t8+t32+t14-t20;
      values[10] = t8+t42-t14+t20;
      values[11] = t8-t32+t14-t20;
      values[12] = -3.0*eta+3.0*t7;
      values[13] = -3.0*t25+3.0*t31;
      values[14] = 3.0-t5-t34+t41;
    break;

    default:
      OutPut("unknown NC_TYPE: " << TDatabase::ParamDB->NC_TYPE << endl);
  }
}

static void N_Q_Q3_2D_DeriveXiXi(double xi, double eta, double *values)
{
 double t1, t2, t3, t4, t5, t6, t7, t9, t10, t11, t12, t13, t14, t15, t16;
 double t17, t18, t19, t20, t21, t22, t24, t26, t27, t28;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = eta*eta;
      t2 = 3.0/4.0*t1;
      t3 = xi*xi;
      t4 = 35.0/4.0*t3;
      t5 = 2.0-t2-t4;
      t6 = 15.0/2.0*xi;
      t7 = 35.0/2.0*t3;
      t10 = xi*eta;
      t13 = xi*(-1.0/2.0+3.0/2.0*t1);
      t15 = 15.0/8.0*eta;
      t16 = 15.0/4.0*t10;
      t17 = t1*eta;
      t18 = 15.0/8.0*t17;
      t21 = (105.0/2.0*t3-15.0/2.0)*eta;
      t22 = t21/4.0;
      t26 = 3.0/8.0*eta;
      t28 = 15.0/4.0*t13;

      values[0] = t5;
      values[1] = -7.0/4.0+t6-t2+t7;
      values[2] = t5;
      values[3] = -7.0/4.0-t6-t2+t7;
      values[4] = 15.0/4.0*t10-15.0/4.0*t13;
      values[5] = t15+t16-t18+t22;
      values[6] = 15.0/4.0*t10+15.0/4.0*t13;
      values[7] = -t15+t16+t18-t22;
      values[8] = 13.0/4.0+t26+t2-t7-t18+t22;
      values[9] = -2.0+t2+t4+t28;
      values[10] = 13.0/4.0-t26+t2-t7+t18-t22;
      values[11] = -2.0+t2+t4-t28;
      values[12] = -1.0/2.0+3.0*t1-t7;
      values[13] = -15.0*xi+15.0/2.0*t13;
      values[14] = -15.0/4.0*eta+15.0/4.0*t17-t21/2.0;
    break;

    case 3: // G. Matthies, only one difference, family 3
      t1 = 15.0/2.0*xi;
      t2 = xi*xi;
      t3 = 105.0/4.0*t2;
      t7 = 15.0/4.0*xi*eta;
      t9 = 105.0/2.0*t2-15.0/2.0;
      t10 = t9*eta;
      t11 = t10/2.0;
      t14 = 3.0/2.0*eta;

      values[0] = 0.0;
      values[1] = t1+t3-15.0/4.0;
      values[2] = 0.0;
      values[3] = -t1+t3-15.0/4.0;
      values[4] = t7;
      values[5] = t7+t11;
      values[6] = t7;
      values[7] = t7-t11;
      values[8] = 21.0/4.0-t14-t3+t11;
      values[9] = 0.0;
      values[10] = 21.0/4.0+t14-t3-t11;
      values[11] = 0.0;
      values[12] = -t9;
      values[13] = -15.0*xi;
      values[14] = -t10;
    break;

    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = eta*eta;
      t2 = 9.0/4.0*t1;
      t3 = xi*xi;
      t4 = 105.0/4.0*t3;
      t5 = 6.0-t2-t4;
      t6 = 15.0/2.0*xi;
      t9 = xi*eta;
      t12 = xi*(-1.0/2.0+3.0/2.0*t1);
      t14 = 15.0/2.0*t9;
      t17 = (105.0/2.0*t3-15.0/2.0)*eta;
      t18 = t17/2.0;
      t22 = 3.0/2.0*eta;
      t24 = 15.0/2.0*t12;

      values[0] = t5;
      values[1] = t6+t4-15.0/4.0;
      values[2] = t5;
      values[3] = -t6+t4-15.0/4.0;
      values[4] = 15.0/2.0*t9-15.0/2.0*t12;
      values[5] = t14+t18;
      values[6] = 15.0/2.0*t9+15.0/2.0*t12;
      values[7] = t14-t18;
      values[8] = 21.0/4.0-t22-t4+t18;
      values[9] = -6.0+t2+t4+t24;
      values[10] = 21.0/4.0+t22-t4-t18;
      values[11] = -6.0+t2+t4-t24;
      values[12] = -9.0/2.0+9.0/2.0*t1;
      values[13] = -15.0*xi+15.0*t12;
      values[14] = -t17;
    break;

    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = eta*eta;
      t2 = 63.0/16.0*t1;
      t3 = xi*xi;
      t7 = -1.0/2.0+3.0/2.0*t1;
      t9 = (105.0/2.0*t3-15.0/2.0)*t7/4.0;
      t10 = t1*t1;
      t11 = 105.0/32.0*t10;
      t12 = 21.0/32.0-t2-t9+t11;
      t13 = 15.0/2.0*xi;
      t14 = 27.0/16.0*t1;
      t17 = xi*t7;
      t18 = 15.0/2.0*t17;
      t19 = xi*eta;
      t20 = 15.0/4.0*t19;
      t22 = t1*eta;
      t26 = 15.0/4.0*t22;
      t27 = 9.0/4.0*eta;

      values[0] = t12;
      values[1] = 51.0/32.0+t13+t14+t9-t11;
      values[2] = t12;
      values[3] = 51.0/32.0-t13+t14+t9-t11;
      values[4] = -t18+t20;
      values[5] = 15.0/4.0*eta-15.0/4.0*t22+15.0/4.0*t19;
      values[6] = t18+t20;
      values[7] = -15.0/4.0*eta+15.0/4.0*t22+15.0/4.0*t19;
      values[8] = -3.0/32.0-t14-t26+t27-t9+t11;
      values[9] = -21.0/32.0+t2+t18+t9-t11;
      values[10] = -3.0/32.0-t14+t26-t27-t9+t11;
      values[11] = -21.0/32.0+t2-t18+t9-t11;
      values[12] = -9.0/2.0+9.0/2.0*t1;
      values[13] = -15.0*xi+15.0*t17;
      values[14] = -15.0/2.0*eta+15.0/2.0*t22;
    break;

    default:
      OutPut("unknown NC_TYPE: " << TDatabase::ParamDB->NC_TYPE << endl);
  }
}

static void N_Q_Q3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15, t16;
  double t17, t18, t19, t20, t21, t22, t23, t24, t25, t27, t30, t31;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = xi*eta;
      t2 = 3.0/2.0*t1;
      t3 = 9.0/8.0*eta;
      t4 = eta*eta;
      t5 = 15.0/8.0*t4;
      t6 = xi*xi;
      t7 = 15.0/8.0*t6;
      t8 = t4*eta;
      t9 = 35.0/8.0*t8;
      t12 = (15.0/2.0*t6-3.0/2.0)*eta;
      t13 = 3.0/4.0*t12;
      t15 = 9.0/8.0*xi;
      t18 = xi*(15.0/2.0*t4-3.0/2.0);
      t19 = 3.0/4.0*t18;
      t20 = t6*xi;
      t21 = 35.0/8.0*t20;
      t25 = 21.0/8.0*xi;
      t27 = 21.0/8.0*eta;

      values[0] = -t2;
      values[1] = -t2;
      values[2] = -t2;
      values[3] = -t2;
      values[4] = -1.0/4.0-t3-t5+t7+t9-t13;
      values[5] = 1.0/4.0-t15-t5+t7-t19+t21;
      values[6] = -1.0/4.0+t3-t5+t7-t9+t13;
      values[7] = 1.0/4.0+t15-t5+t7+t19-t21;
      values[8] = -t25+t2-t19+t21;
      values[9] = t27+t2-t9+t13;
      values[10] = t25+t2+t19-t21;
      values[11] = -t27+t2+t9-t13;
      values[12] = 6.0*t1;
      values[13] = 9.0/4.0*eta-35.0/4.0*t8+3.0/2.0*t12;
      values[14] = 9.0/4.0*xi+3.0/2.0*t18-35.0/4.0*t20;
    break;

    case 3: // G. Matthies, only one difference, family 3
      t1 = eta*eta;
      t2 = 15.0/8.0*t1;
      t3 = xi*xi;
      t4 = 15.0/8.0*t3;
      t5 = t1*eta;
      t6 = 35.0/4.0*t5;
      t7 = 15.0/4.0*eta;
      t9 = t3*xi;
      t10 = 35.0/4.0*t9;
      t11 = 15.0/4.0*xi;
      t16 = -21.0/4.0*xi+t10;
      t18 = 21.0/4.0*eta-t6;

      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 0.0;
      values[3] = 0.0;
      values[4] = -1.0/4.0-t2+t4+t6-t7;
      values[5] = 1.0/4.0-t2+t4+t10-t11;
      values[6] = -1.0/4.0-t2+t4-t6+t7;
      values[7] = 1.0/4.0-t2+t4-t10+t11;
      values[8] = t16;
      values[9] = t18;
      values[10] = -t16;
      values[11] = -t18;
      values[12] = 0.0;
      values[13] = -35.0/2.0*t5+15.0/2.0*eta;
      values[14] = -35.0/2.0*t9+15.0/2.0*xi;
    break;

    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = xi*eta;
      t2 = 9.0/2.0*t1;
      t3 = 3.0/2.0*eta;
      t4 = xi*xi;
      t5 = 15.0/4.0*t4;
      t8 = (15.0/2.0*t4-3.0/2.0)*eta;
      t9 = 3.0/2.0*t8;
      t11 = xi*t4;
      t12 = 35.0/4.0*t11;
      t13 = 15.0/4.0*xi;
      t18 = -21.0/4.0*xi+t12;

      values[0] = -t2;
      values[1] = 0.0;
      values[2] = -t2;
      values[3] = 0.0;
      values[4] = -5.0/4.0+t3+t5-t9;
      values[5] = -3.0/4.0+t5+t12-t13;
      values[6] = -5.0/4.0-t3+t5+t9;
      values[7] = t5-3.0/4.0-t12+t13;
      values[8] = t18;
      values[9] = t2+t9;
      values[10] = -t18;
      values[11] = t2-t9;
      values[12] = 9.0*t1;
      values[13] = -3.0*eta+3.0*t8;
      values[14] = -35.0/2.0*t11+15.0/2.0*xi;
    break;

    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = xi*eta;
      t2 = 9.0/4.0*t1;
      t3 = xi*xi;
      t9 = 3.0/4.0*(35.0/2.0*t3*xi-15.0/2.0*xi)*eta;
      t10 = eta*eta;
      t16 = 3.0/4.0*xi*(35.0/2.0*t10*eta-15.0/2.0*eta);
      t17 = -t2-t9+t16;
      t18 = -t2+t9-t16;
      t19 = 3.0/2.0*eta;
      t22 = (15.0/2.0*t3-3.0/2.0)*eta;
      t23 = 3.0/2.0*t22;
      t24 = 15.0/8.0*t3;
      t25 = 15.0/8.0*t10;
      t27 = 3.0/2.0*xi;
      t30 = xi*(15.0/2.0*t10-3.0/2.0);
      t31 = 3.0/2.0*t30;

      values[0] = t17;
      values[1] = t18;
      values[2] = t17;
      values[3] = t18;
      values[4] = -1.0/4.0+t19-t23+t24-t25;
      values[5] = 1.0/4.0+t27-t31+t24-t25;
      values[6] = -1.0/4.0-t19+t23+t24-t25;
      values[7] = 1.0/4.0-t27+t31+t24-t25;
      values[8] = t2-t31-t9+t16;
      values[9] = t2+t23+t9-t16;
      values[10] = t2+t31-t9+t16;
      values[11] = t2-t23+t9-t16;
      values[12] = 9.0*t1;
      values[13] = -3.0*eta+3.0*t22;
      values[14] = -3.0*xi+3.0*t30;
    break;

    default:
      OutPut("unknown NC_TYPE: " << TDatabase::ParamDB->NC_TYPE << endl);
  }
}

static void N_Q_Q3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t14, t15, t16;
  double t17, t18, t20, t21, t22, t26, t27, t28;

  switch(TDatabase::ParamDB->NC_TYPE)
  {
    case 2: // G. Matthies, many differences, family 2
      t1 = 15.0/2.0*eta;
      t2 = eta*eta;
      t3 = 35.0/2.0*t2;
      t4 = xi*xi;
      t5 = 3.0/4.0*t4;
      t7 = 35.0/4.0*t2;
      t8 = 2.0-t7-t5;
      t10 = 15.0/8.0*xi;
      t11 = xi*eta;
      t12 = 15.0/4.0*t11;
      t15 = xi*(105.0/2.0*t2-15.0/2.0);
      t16 = t15/4.0;
      t17 = xi*t4;
      t18 = 15.0/8.0*t17;
      t22 = (-1.0/2.0+3.0/2.0*t4)*eta;
      t26 = 15.0/4.0*t22;
      t28 = 3.0/8.0*xi;

      values[0] = -7.0/4.0-t1+t3-t5;
      values[1] = t8;
      values[2] = -7.0/4.0+t1+t3-t5;
      values[3] = t8;
      values[4] = t10-t12+t16-t18;
      values[5] = -15.0/4.0*t11-15.0/4.0*t22;
      values[6] = -t10-t12-t16+t18;
      values[7] = -15.0/4.0*t11+15.0/4.0*t22;
      values[8] = -2.0+t7+t5-t26;
      values[9] = 13.0/4.0-t28-t3+t5-t16+t18;
      values[10] = -2.0+t7+t5+t26;
      values[11] = 13.0/4.0+t28-t3+t5+t16-t18;
      values[12] = -1.0/2.0-t3+3.0*t4;
      values[13] = -15.0/4.0*xi-t15/2.0+15.0/4.0*t17;
      values[14] = -15.0*eta+15.0/2.0*t22;
    break;

    case 3: // G. Matthies, only one difference, family 3
      t1 = 15.0/2.0*eta;
      t2 = eta*eta;
      t3 = 105.0/4.0*t2;
      t7 = 15.0/4.0*xi*eta;
      t9 = 105.0/2.0*t2-15.0/2.0;
      t10 = xi*t9;
      t11 = t10/2.0;
      t14 = 3.0/2.0*xi;

      values[0] = -t1+t3-15.0/4.0;
      values[1] = 0.0;
      values[2] = t1+t3-15.0/4.0;
      values[3] = 0.0;
      values[4] = -t7+t11;
      values[5] = -t7;
      values[6] = -t7-t11;
      values[7] = -t7;
      values[8] = 0.0;
      values[9] = 21.0/4.0+t14-t3-t11;
      values[10] = 0.0;
      values[11] = 21.0/4.0-t14-t3+t11;
      values[12] = -t9;
      values[13] = -t10;
      values[14] = -15.0*eta;
    break;

    case 4: // Apel, Matthies, anisotrop, family 4
      t1 = 15.0/2.0*eta;
      t2 = xi*xi;
      t3 = 9.0/4.0*t2;
      t6 = t2*xi;
      t7 = xi-t6;
      t8 = 15.0/4.0*t6;
      t9 = 9.0/4.0*xi;

      values[0] = 9.0/4.0-t1-t3;
      values[1] = 0.0;
      values[2] = 9.0/4.0+t1-t3;
      values[3] = 0.0;
      values[4] = 15.0/4.0*t7;
      values[5] = 0.0;
      values[6] = -15.0/4.0*t7;
      values[7] = 0.0;
      values[8] = 0.0;
      values[9] = -3.0/4.0+t3+t8-t9;
      values[10] = 0.0;
      values[11] = -3.0/4.0+t3-t8+t9;
      values[12] = -9.0/2.0+9.0/2.0*t2;
      values[13] = -15.0/2.0*t7;
      values[14] = -15.0*eta;
    break;

    case 1: // Hennart, Jaffr'e, Roberts, 1988
      t1 = 15.0/2.0*eta;
      t2 = xi*xi;
      t3 = 27.0/16.0*t2;
      t4 = t2*t2;
      t5 = 105.0/32.0*t4;
      t7 = -1.0/2.0+3.0/2.0*t2;
      t8 = eta*eta;
      t12 = t7*(105.0/2.0*t8-15.0/2.0)/4.0;
      t14 = 63.0/16.0*t2;
      t15 = 21.0/32.0-t14+t5-t12;
      t17 = t2*xi;
      t18 = xi*eta;
      t20 = t7*eta;
      t21 = 15.0/2.0*t20;
      t22 = 15.0/4.0*t18;
      t27 = 15.0/4.0*t17;
      t28 = 9.0/4.0*xi;

      values[0] = 51.0/32.0-t1+t3-t5+t12;
      values[1] = t15;
      values[2] = 51.0/32.0+t1+t3-t5+t12;
      values[3] = t15;
      values[4] = 15.0/4.0*xi-15.0/4.0*t17-15.0/4.0*t18;
      values[5] = -t21-t22;
      values[6] = -15.0/4.0*xi+15.0/4.0*t17-15.0/4.0*t18;
      values[7] = t21-t22;
      values[8] = -21.0/32.0+t14-t21-t5+t12;
      values[9] = -3.0/32.0-t3+t27-t28+t5-t12;
      values[10] = -21.0/32.0+t14+t21-t5+t12;
      values[11] = -3.0/32.0-t3-t27+t28+t5-t12;
      values[12] = -9.0/2.0+9.0/2.0*t2;
      values[13] = -15.0/2.0*xi+15.0/2.0*t17;
      values[14] = -15.0*eta+15.0*t20;
    break;

    default:
      OutPut("unknown NC_TYPE: " << TDatabase::ParamDB->NC_TYPE << endl);
  }
}

static int N_Q_Q3_2D_ChangeJ0[1] = { 4 };
static int N_Q_Q3_2D_ChangeJ1[1] = { 5 };
static int N_Q_Q3_2D_ChangeJ2[1] = { 6 };
static int N_Q_Q3_2D_ChangeJ3[1] = { 7 };

static int *N_Q_Q3_2D_Change[4] = { N_Q_Q3_2D_ChangeJ0, N_Q_Q3_2D_ChangeJ1,
                                    N_Q_Q3_2D_ChangeJ2, N_Q_Q3_2D_ChangeJ3 };

// ***********************************************************************

TBaseFunct2D *BF_N_Q_Q3_2D_Obj = new TBaseFunct2D
        (15, BF_N_Q_Q3_2D, BFUnitSquare, 
         N_Q_Q3_2D_Funct, N_Q_Q3_2D_DeriveXi,
         N_Q_Q3_2D_DeriveEta, N_Q_Q3_2D_DeriveXiXi,
         N_Q_Q3_2D_DeriveXiEta, N_Q_Q3_2D_DeriveEtaEta, 4, 3,
         1, N_Q_Q3_2D_Change);
