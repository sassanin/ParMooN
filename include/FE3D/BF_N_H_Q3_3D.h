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
// Q3Rot element, nonconforming, 3D
// ***********************************************************************

static void N_H_Q3_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t16, t18;
  double t19, t20, t21, t24, t26, t28, t29, t30, t31, t32, t35, t37, t39;
  double t40, t41, t43, t44, t46, t47, t48, t50, t51, t52, t53, t58, t60;
  double t62, t63, t65, t67, t68, t69, t70, t71, t72, t74, t75, t76, t78;
  double t79, t80, t81, t82, t83, t84, t86, t88, t89, t91, t93, t95, t96;
  double t97, t98, t101, t102, t103, t104, t105, t106, t107, t108, t109;
  double t110, t111, t113, t114, t115, t116, t119, t120, t121, t133, t134;
  double t135, t136, t137, t138, t139, t141, t143, t144, t145, t147, t148;

  t1 = zeta*zeta;
  t2 = 3.0/4.0*t1;
  t3 = eta*eta;
  t4 = 3.0/2.0*t3;
  t5 = -1.0/2.0+t4;
  t6 = 3.0/2.0*t1;
  t7 = -1.0/2.0+t6;
  t8 = t5*t7;
  t9 = t8/4.0;
  t10 = t1*zeta;
  t11 = 5.0/4.0*t10;
  t12 = xi*xi;
  t13 = t12*t12;
  t16 = 3.0/8.0+35.0/8.0*t13-15.0/4.0*t12;
  t18 = t16*t7/4.0;
  t19 = 3.0/2.0*t12;
  t20 = -1.0/2.0+t19;
  t21 = t1*t1;
  t24 = 3.0/8.0+35.0/8.0*t21-15.0/4.0*t1;
  t26 = t20*t24/4.0;
  t28 = t5*t24/4.0;
  t29 = t20*t7;
  t30 = t29/4.0;
  t31 = 3.0/4.0*zeta;
  t32 = t3*t3;
  t35 = 3.0/8.0+35.0/8.0*t32-15.0/4.0*t3;
  t37 = t35*t7/4.0;
  t39 = t3*eta;
  t40 = 5.0/4.0*t39;
  t41 = 3.0/4.0*eta;
  t43 = t16*t5/4.0;
  t44 = 3.0/4.0*t3;
  t46 = t20*t35/4.0;
  t47 = t20*t5;
  t48 = t47/4.0;
  t50 = 3.0/4.0*xi;
  t51 = t12*xi;
  t52 = 5.0/4.0*t51;
  t53 = 3.0/4.0*t12;
  t58 = 5.0/2.0*t51;
  t60 = t58-3.0/2.0*xi;
  t62 = t60*zeta/4.0;
  t63 = 5.0/2.0*t10;
  t65 = t63-3.0/2.0*zeta;
  t67 = xi*t65/4.0;
  t68 = t60*t7;
  t69 = t68/2.0;
  t70 = xi*t5;
  t71 = t70*zeta;
  t72 = t71/4.0;
  t74 = xi*zeta/4.0;
  t75 = xi*t7;
  t76 = t75/2.0;
  t78 = t20*eta;
  t79 = t78*zeta;
  t80 = t79/4.0;
  t81 = t5*t65;
  t82 = t81/2.0;
  t83 = t5*zeta;
  t84 = t83/2.0;
  t86 = eta*zeta/4.0;
  t88 = eta*t65/4.0;
  t89 = 5.0/2.0*t39;
  t91 = t89-3.0/2.0*eta;
  t93 = t91*zeta/4.0;
  t95 = t20*t65;
  t96 = t95/2.0;
  t97 = t20*zeta;
  t98 = t97/2.0;
  t101 = t60*eta;
  t102 = t101/4.0;
  t103 = xi*t91;
  t104 = t103/4.0;
  t105 = xi*eta;
  t106 = t105/4.0;
  t107 = t105*t7;
  t108 = t107/4.0;
  t109 = t78/2.0;
  t110 = t20*t91;
  t111 = t110/2.0;
  t113 = eta*t7;
  t114 = t113/2.0;
  t115 = t91*t7;
  t116 = t115/2.0;
  t119 = t60*t5;
  t120 = t119/2.0;
  t121 = t70/2.0;
  t133 = t105*t65;
  t134 = t133/3.0;
  t135 = t107/2.0;
  t136 = t101*zeta;
  t137 = t136/6.0;
  t138 = t103*zeta;
  t139 = t138/6.0;
  t141 = t105*zeta/6.0;
  t143 = t133/6.0;
  t144 = t71/2.0;
  t145 = t138/3.0;
  t147 = t79/2.0;
  t148 = t136/3.0;

  values[0] = t2-t9-t11-t18+t26+t28-t30-1.0/4.0+t31-t37;
  values[1] = -t40-t9+t41-t28-t43+t44+t46-t48-1.0/4.0+t37;
  values[2] = -t50+t52+t18-t26+t43-t46-t30-t48-1.0/4.0+t53;
  values[3] = t40-t9-t41-t28-t43+t44+t46-t48-1.0/4.0+t37;
  values[4] = t50-t52+t18-t26+t43-t46-t30-t48-1.0/4.0+t53;
  values[5] = t2-t9+t11-t18+t26+t28-t30-1.0/4.0-t31-t37;
  values[6] = t62-t67-t69+t72-t74+t76;
  values[7] = t80-t82+t84-t86+t88-t93;
  values[8] = t62-t67-t96-t72+t74+t98;
  values[9] = -t80-t82+t84+t86-t88+t93;
  values[10] = -t102+t104-t106+t108+t109-t111;
  values[11] = t114-t80-t116+t86+t88-t93;
  values[12] = t114+t80-t116-t86-t88+t93;
  values[13] = t102-t104-t106-t120+t108+t121;
  values[14] = t102-t104+t106-t108+t109-t111;
  values[15] = t102-t104-t106+t120+t108-t121;
  values[16] = -t62+t67-t96+t72-t74+t98;
  values[17] = -t62+t67-t69-t72+t74+t76;
  values[18] = -t96-t18+t26+t30;
  values[19] = t9-t28-t116+t37;
  values[20] = t69+t18-t26+t30;
  values[21] = t9-t28+t116+t37;
  values[22] = -t120+t43-t46+t48;
  values[23] = t9+t28+t82-t37;
  values[24] = -t134+t135+t137+t139-t141;
  values[25] = t143+t144+t137-t145-t141;
  values[26] = -t143+t147+t148-t139+t141;
  values[27] = t143-t144+t137-t145-t141;
  values[28] = t143+t147-t148+t139-t141;
  values[29] = t134+t135-t137-t139+t141;
  values[30] = t9+t28-t82-t37;
  values[31] = -t43+t46-t111+t48;
  values[32] = t120+t43-t46+t48;
  values[33] = -t43+t46+t111+t48;
  values[34] = -t69+t18-t26+t30;
  values[35] = t96-t18+t26+t30;
  values[36] = -t6+t8+5.0/2.0-t4+t29+t47-t19;
  values[37] = t119+5.0/2.0*xi-t58+t68-t75-t70;
  values[38] = -t89-t113+5.0/2.0*eta+t115-t78+t110;
  values[39] = t95-t63+t81-t97+5.0/2.0*zeta-t83;
}

static void N_H_Q3_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t5, t6, t7, t8, t10, t11, t16, t17, t18, t19, t20, t21, t22;
  double t24, t25, t30, t31, t32, t33, t34, t35, t38, t39, t41, t42, t43;
  double t44, t45, t46, t47, t48, t49, t51, t52, t53, t54, t57, t58, t59;
  double t60, t61, t63, t64, t65, t66, t67, t68, t69, t70, t73, t74, t75;
  double t77, t78, t79, t80, t89, t90, t91, t92, t93, t94, t95, t97, t99;
  double t100, t101, t103, t104;

  t1 = xi*xi;
  t5 = 35.0/2.0*t1*xi-15.0/2.0*xi;
  t6 = zeta*zeta;
  t7 = 3.0/2.0*t6;
  t8 = -1.0/2.0+t7;
  t10 = t5*t8/4.0;
  t11 = t6*t6;
  t16 = 3.0/4.0*xi*(3.0/8.0+35.0/8.0*t11-15.0/4.0*t6);
  t17 = xi*t8;
  t18 = 3.0/4.0*t17;
  t19 = -t10+t16-t18;
  t20 = eta*eta;
  t21 = 3.0/2.0*t20;
  t22 = -1.0/2.0+t21;
  t24 = t5*t22/4.0;
  t25 = t20*t20;
  t30 = 3.0/4.0*xi*(3.0/8.0+35.0/8.0*t25-15.0/4.0*t20);
  t31 = xi*t22;
  t32 = 3.0/4.0*t31;
  t33 = -t24+t30-t32;
  t34 = 15.0/4.0*t1;
  t35 = 3.0/2.0*xi;
  t38 = 15.0/2.0*t1;
  t39 = t38-3.0/2.0;
  t41 = t39*zeta/4.0;
  t42 = t6*zeta;
  t43 = 5.0/8.0*t42;
  t44 = zeta/8.0;
  t45 = t39*t8;
  t46 = t45/2.0;
  t47 = t22*zeta;
  t48 = t47/4.0;
  t49 = 3.0/4.0*t6;
  t51 = xi*eta;
  t52 = t51*zeta;
  t53 = 3.0/4.0*t52;
  t54 = 5.0/8.0*zeta;
  t57 = 5.0/2.0*t42-3.0/2.0*zeta;
  t58 = xi*t57;
  t59 = 3.0/2.0*t58;
  t60 = xi*zeta;
  t61 = 3.0/2.0*t60;
  t63 = t39*eta;
  t64 = t63/4.0;
  t65 = t20*eta;
  t66 = 5.0/8.0*t65;
  t67 = 5.0/8.0*eta;
  t68 = eta*t8;
  t69 = t68/4.0;
  t70 = 3.0/2.0*t51;
  t73 = 5.0/2.0*t65-3.0/2.0*eta;
  t74 = xi*t73;
  t75 = 3.0/2.0*t74;
  t77 = eta/8.0;
  t78 = t39*t22;
  t79 = t78/2.0;
  t80 = 3.0/4.0*t20;
  t89 = eta*t57;
  t90 = t89/3.0;
  t91 = t68/2.0;
  t92 = t63*zeta;
  t93 = t92/6.0;
  t94 = t73*zeta;
  t95 = t94/6.0;
  t97 = eta*zeta/6.0;
  t99 = t89/6.0;
  t100 = t47/2.0;
  t101 = t94/3.0;
  t103 = 3.0/2.0*t52;
  t104 = t92/3.0;

  values[0] = t19;
  values[1] = t33;
  values[2] = -3.0/4.0+t34+t10-t16+t24-t30-t18-t32+t35;
  values[3] = t33;
  values[4] = 3.0/4.0-t34+t10-t16+t24-t30-t18-t32+t35;
  values[5] = t19;
  values[6] = t41-t43+t44-t46+t48-1.0/4.0+t49;
  values[7] = t53;
  values[8] = t41-t43+t54-t59-t48+t61;
  values[9] = -t53;
  values[10] = -t64+t66-t67+t69+t70-t75;
  values[11] = -t53;
  values[12] = t53;
  values[13] = t64-t66+t77-t79+t69-1.0/4.0+t80;
  values[14] = t64-t66+t67-t69+t70-t75;
  values[15] = t64-t66+t77+t79+t69+1.0/4.0-t80;
  values[16] = -t41+t43-t54-t59+t48+t61;
  values[17] = -t41+t43-t44-t46-t48-1.0/4.0+t49;
  values[18] = -t59-t10+t16+t18;
  values[19] = 0.0;
  values[20] = t46+t10-t16+t18;
  values[21] = 0.0;
  values[22] = -t79+t24-t30+t32;
  values[23] = 0.0;
  values[24] = -t90+t91+t93+t95-t97;
  values[25] = t99+t100+t93-t101-t97;
  values[26] = -t99+t103+t104-t95+t97;
  values[27] = t99-t100+t93-t101-t97;
  values[28] = t99+t103-t104+t95-t97;
  values[29] = t90+t91-t93-t95+t97;
  values[30] = 0.0;
  values[31] = -t24+t30-t75+t32;
  values[32] = t79+t24-t30+t32;
  values[33] = -t24+t30+t75+t32;
  values[34] = -t46+t10-t16+t18;
  values[35] = t59-t10+t16+t18;
  values[36] = 3.0*t17+3.0*t31-3.0*xi;
  values[37] = t78+7.0/2.0-t38+t45-t7-t21;
  values[38] = -3.0*t51+3.0*t74;
  values[39] = 3.0*t58-3.0*t60;
}

static void N_H_Q3_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t11, t12, t16, t18, t19, t20, t21, t22;
  double t27, t28, t29, t30, t32, t33, t34, t36, t38, t39, t40, t41, t42;
  double t43, t46, t47, t48, t49, t50, t51, t52, t53, t54, t56, t59, t60;
  double t61, t62, t63, t64, t65, t66, t67, t68, t70, t71, t72, t73, t76;
  double t79, t80, t81, t82, t90, t91, t92, t93, t94, t95, t96, t98, t100;
  double t101, t102, t104, t105;

  t1 = zeta*zeta;
  t2 = 3.0/2.0*t1;
  t3 = -1.0/2.0+t2;
  t4 = t3*eta;
  t5 = 3.0/4.0*t4;
  t6 = t1*t1;
  t11 = 3.0/4.0*eta*(3.0/8.0+35.0/8.0*t6-15.0/4.0*t1);
  t12 = eta*eta;
  t16 = 35.0/2.0*t12*eta-15.0/2.0*eta;
  t18 = t16*t3/4.0;
  t19 = -t5+t11-t18;
  t20 = 15.0/4.0*t12;
  t21 = xi*xi;
  t22 = t21*t21;
  t27 = 3.0/4.0*(3.0/8.0+35.0/8.0*t22-15.0/4.0*t21)*eta;
  t28 = 3.0/2.0*eta;
  t29 = 3.0/2.0*t21;
  t30 = -1.0/2.0+t29;
  t32 = t30*t16/4.0;
  t33 = t30*eta;
  t34 = 3.0/4.0*t33;
  t36 = t27-t32-t34;
  t38 = xi*eta;
  t39 = t38*zeta;
  t40 = 3.0/4.0*t39;
  t41 = t30*zeta;
  t42 = t41/4.0;
  t43 = t1*zeta;
  t46 = 5.0/2.0*t43-3.0/2.0*zeta;
  t47 = eta*t46;
  t48 = 3.0/2.0*t47;
  t49 = eta*zeta;
  t50 = 3.0/2.0*t49;
  t51 = 5.0/8.0*zeta;
  t52 = 5.0/8.0*t43;
  t53 = 15.0/2.0*t12;
  t54 = t53-3.0/2.0;
  t56 = t54*zeta/4.0;
  t59 = t21*xi;
  t60 = 5.0/8.0*t59;
  t61 = xi/8.0;
  t62 = xi*t54;
  t63 = t62/4.0;
  t64 = xi*t3;
  t65 = t64/4.0;
  t66 = 3.0/4.0*t21;
  t67 = t30*t54;
  t68 = t67/2.0;
  t70 = 3.0/4.0*t1;
  t71 = t54*t3;
  t72 = t71/2.0;
  t73 = zeta/8.0;
  t76 = 5.0/8.0*xi;
  t79 = 5.0/2.0*t59-3.0/2.0*xi;
  t80 = t79*eta;
  t81 = 3.0/2.0*t80;
  t82 = 3.0/2.0*t38;
  t90 = xi*t46;
  t91 = t90/3.0;
  t92 = t64/2.0;
  t93 = t79*zeta;
  t94 = t93/6.0;
  t95 = t62*zeta;
  t96 = t95/6.0;
  t98 = xi*zeta/6.0;
  t100 = t90/6.0;
  t101 = 3.0/2.0*t39;
  t102 = t95/3.0;
  t104 = t41/2.0;
  t105 = t93/3.0;

  values[0] = t19;
  values[1] = -t20-t5+3.0/4.0-t11-t27+t28+t32-t34+t18;
  values[2] = t36;
  values[3] = t20-t5-3.0/4.0-t11-t27+t28+t32-t34+t18;
  values[4] = t36;
  values[5] = t19;
  values[6] = t40;
  values[7] = t42-t48+t50-t51+t52-t56;
  values[8] = -t40;
  values[9] = -t42-t48+t50+t51-t52+t56;
  values[10] = -t60+t61+t63+t65-1.0/4.0+t66-t68;
  values[11] = -1.0/4.0+t70-t42-t72-t73+t52-t56;
  values[12] = -1.0/4.0+t70+t42-t72+t73-t52+t56;
  values[13] = t60-t76-t63-t81+t65+t82;
  values[14] = t60-t61-t63-t65-1.0/4.0+t66-t68;
  values[15] = t60-t76-t63+t81+t65-t82;
  values[16] = t40;
  values[17] = -t40;
  values[18] = 0.0;
  values[19] = t5-t11-t72+t18;
  values[20] = 0.0;
  values[21] = t5-t11+t72+t18;
  values[22] = -t81+t27-t32+t34;
  values[23] = t5+t11+t48-t18;
  values[24] = -t91+t92+t94+t96-t98;
  values[25] = t100+t101+t94-t102-t98;
  values[26] = -t100+t104+t105-t96+t98;
  values[27] = t100-t101+t94-t102-t98;
  values[28] = t100+t104-t105+t96-t98;
  values[29] = t91+t92-t94-t96+t98;
  values[30] = t5+t11-t48-t18;
  values[31] = -t27+t32-t68+t34;
  values[32] = t81+t27-t32+t34;
  values[33] = -t27+t32+t68+t34;
  values[34] = 0.0;
  values[35] = 0.0;
  values[36] = 3.0*t4-3.0*eta+3.0*t33;
  values[37] = 3.0*t80-3.0*t38;
  values[38] = -t53+7.0/2.0-t2+t71-t29+t67;
  values[39] = 3.0*t47-3.0*t49;
}

static void N_H_Q3_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t15, t16, t17, t21, t23;
  double t25, t26, t27, t28, t33, t35, t36, t38, t39, t40, t41, t42, t44;
  double t47, t48, t49, t50, t51, t52, t53, t55, t56, t57, t58, t59, t60;
  double t62, t63, t64, t66, t67, t68, t69, t72, t73, t74, t75, t76, t79;
  double t80, t81, t82, t92, t93, t94, t95, t96, t97, t98, t99, t101, t102;
  double t103, t105, t106;

  t1 = 3.0/2.0*zeta;
  t2 = eta*eta;
  t3 = 3.0/2.0*t2;
  t4 = -1.0/2.0+t3;
  t5 = t4*zeta;
  t6 = 3.0/4.0*t5;
  t7 = zeta*zeta;
  t8 = 15.0/4.0*t7;
  t9 = xi*xi;
  t10 = t9*t9;
  t15 = 3.0/4.0*(3.0/8.0+35.0/8.0*t10-15.0/4.0*t9)*zeta;
  t16 = 3.0/2.0*t9;
  t17 = -1.0/2.0+t16;
  t21 = 35.0/2.0*t7*zeta-15.0/2.0*zeta;
  t23 = t17*t21/4.0;
  t25 = t4*t21/4.0;
  t26 = t17*zeta;
  t27 = 3.0/4.0*t26;
  t28 = t2*t2;
  t33 = 3.0/4.0*(3.0/8.0+35.0/8.0*t28-15.0/4.0*t2)*zeta;
  t35 = -t6-t25+t33;
  t36 = t15-t23-t27;
  t38 = t9*xi;
  t39 = 5.0/8.0*t38;
  t40 = 5.0/8.0*xi;
  t41 = 15.0/2.0*t7;
  t42 = t41-3.0/2.0;
  t44 = xi*t42/4.0;
  t47 = 5.0/2.0*t38-3.0/2.0*xi;
  t48 = t47*zeta;
  t49 = 3.0/2.0*t48;
  t50 = xi*t4;
  t51 = t50/4.0;
  t52 = xi*zeta;
  t53 = 3.0/2.0*t52;
  t55 = t17*eta;
  t56 = t55/4.0;
  t57 = t4*t42;
  t58 = t57/2.0;
  t59 = 3.0/4.0*t2;
  t60 = eta/8.0;
  t62 = eta*t42/4.0;
  t63 = t2*eta;
  t64 = 5.0/8.0*t63;
  t66 = xi/8.0;
  t67 = t17*t42;
  t68 = t67/2.0;
  t69 = 3.0/4.0*t9;
  t72 = xi*eta;
  t73 = t72*zeta;
  t74 = 3.0/4.0*t73;
  t75 = eta*zeta;
  t76 = 3.0/2.0*t75;
  t79 = 5.0/2.0*t63-3.0/2.0*eta;
  t80 = t79*zeta;
  t81 = 3.0/2.0*t80;
  t82 = 5.0/8.0*eta;
  t92 = t72*t42;
  t93 = t92/3.0;
  t94 = 3.0/2.0*t73;
  t95 = t47*eta;
  t96 = t95/6.0;
  t97 = xi*t79;
  t98 = t97/6.0;
  t99 = t72/6.0;
  t101 = t92/6.0;
  t102 = t50/2.0;
  t103 = t97/3.0;
  t105 = t55/2.0;
  t106 = t95/3.0;

  values[0] = t1-t6-t8-t15+t23+t25-t27+3.0/4.0-t33;
  values[1] = t35;
  values[2] = t36;
  values[3] = t35;
  values[4] = t36;
  values[5] = t1-t6+t8-t15+t23+t25-t27-3.0/4.0-t33;
  values[6] = t39-t40-t44-t49+t51+t53;
  values[7] = t56-t58-1.0/4.0+t59+t60+t62-t64;
  values[8] = t39-t66-t44-t68-t51-1.0/4.0+t69;
  values[9] = -t56-t58-1.0/4.0+t59-t60-t62+t64;
  values[10] = t74;
  values[11] = t76-t56-t81+t82+t62-t64;
  values[12] = t76+t56-t81-t82-t62+t64;
  values[13] = t74;
  values[14] = -t74;
  values[15] = t74;
  values[16] = -t39+t66+t44-t68+t51-1.0/4.0+t69;
  values[17] = -t39+t40+t44-t49-t51+t53;
  values[18] = -t68-t15+t23+t27;
  values[19] = t6-t25-t81+t33;
  values[20] = t49+t15-t23+t27;
  values[21] = t6-t25+t81+t33;
  values[22] = 0.0;
  values[23] = t6+t25+t58-t33;
  values[24] = -t93+t94+t96+t98-t99;
  values[25] = t101+t102+t96-t103-t99;
  values[26] = -t101+t105+t106-t98+t99;
  values[27] = t101-t102+t96-t103-t99;
  values[28] = t101+t105-t106+t98-t99;
  values[29] = t93+t94-t96-t98+t99;
  values[30] = t6+t25-t58-t33;
  values[31] = 0.0;
  values[32] = 0.0;
  values[33] = 0.0;
  values[34] = -t49+t15-t23+t27;
  values[35] = t68-t15+t23+t27;
  values[36] = -3.0*zeta+3.0*t5+3.0*t26;
  values[37] = 3.0*t48-3.0*t52;
  values[38] = -3.0*t75+3.0*t80;
  values[39] = t67-t41+t57+7.0/2.0-t16-t3;
}

static void N_H_Q3_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t4, t6, t8, t9, t10, t11, t12, t13, t15, t17, t18, t19;
  double t20, t21, t22, t23, t24, t27, t28, t29, t30, t32, t33, t34, t36;
  double t37, t39, t40, t41, t47, t48, t52, t53, t54, t55, t58, t59;

  t1 = xi*xi;
  t3 = 105.0/2.0*t1-15.0/2.0;
  t4 = zeta*zeta;
  t6 = -1.0/2.0+3.0/2.0*t4;
  t8 = t3*t6/4.0;
  t9 = t4*t4;
  t10 = 105.0/32.0*t9;
  t11 = 63.0/16.0*t4;
  t12 = -t8+21.0/32.0+t10-t11;
  t13 = eta*eta;
  t15 = -1.0/2.0+3.0/2.0*t13;
  t17 = t3*t15/4.0;
  t18 = t13*t13;
  t19 = 105.0/32.0*t18;
  t20 = 63.0/16.0*t13;
  t21 = -t17+21.0/32.0+t19-t20;
  t22 = 15.0/2.0*xi;
  t23 = 27.0/16.0*t4;
  t24 = 27.0/16.0*t13;
  t27 = xi*zeta;
  t28 = 15.0/4.0*t27;
  t29 = xi*t6;
  t30 = 15.0/2.0*t29;
  t32 = eta*zeta;
  t33 = 3.0/4.0*t32;
  t34 = t4*zeta;
  t36 = xi*eta;
  t37 = t13*eta;
  t39 = 15.0/4.0*t36;
  t40 = xi*t15;
  t41 = 15.0/2.0*t40;
  t47 = 15.0/4.0*t34;
  t48 = 9.0/4.0*zeta;
  t52 = t36*zeta;
  t53 = 5.0/2.0*t52;
  t54 = 3.0/2.0*t32;
  t55 = 5.0*t52;
  t58 = 15.0/4.0*t37;
  t59 = 9.0/4.0*eta;

  values[0] = t12;
  values[1] = t21;
  values[2] = t22+t8+27.0/16.0-t10+t23+t17-t19+t24;
  values[3] = t21;
  values[4] = -t22+t8+27.0/16.0-t10+t23+t17-t19+t24;
  values[5] = t12;
  values[6] = t28-t30;
  values[7] = t33;
  values[8] = 15.0/4.0*t27-15.0/4.0*t34+15.0/4.0*zeta;
  values[9] = -t33;
  values[10] = -15.0/4.0*t36+15.0/4.0*eta-15.0/4.0*t37;
  values[11] = -t33;
  values[12] = t33;
  values[13] = t39-t41;
  values[14] = 15.0/4.0*t36+15.0/4.0*eta-15.0/4.0*t37;
  values[15] = t39+t41;
  values[16] = -15.0/4.0*t27-15.0/4.0*t34+15.0/4.0*zeta;
  values[17] = -t28-t30;
  values[18] = -t47+t48-t8-3.0/32.0+t10-t23;
  values[19] = 0.0;
  values[20] = t30+t8-21.0/32.0-t10+t11;
  values[21] = 0.0;
  values[22] = -t41+t17-21.0/32.0-t19+t20;
  values[23] = 0.0;
  values[24] = t53;
  values[25] = t53;
  values[26] = t54+t55;
  values[27] = t53;
  values[28] = t54-t55;
  values[29] = -t53;
  values[30] = 0.0;
  values[31] = -t17-3.0/32.0+t19-t24-t58+t59;
  values[32] = t41+t17-21.0/32.0-t19+t20;
  values[33] = -t17-3.0/32.0+t19-t24+t58-t59;
  values[34] = -t30+t8-21.0/32.0-t10+t11;
  values[35] = t47-t48-t8-3.0/32.0+t10-t23;
  values[36] = -6.0+9.0/2.0*t4+9.0/2.0*t13;
  values[37] = 15.0*t40-15.0*xi+15.0*t29;
  values[38] = -15.0/2.0*eta+15.0/2.0*t37;
  values[39] = 15.0/2.0*t34-15.0/2.0*zeta;
}

static void N_H_Q3_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t7, t8, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24;
  double t25, t26, t27, t29, t30, t31, t34, t35, t36, t37, t42, t43, t44;
  double t45, t46, t47, t48, t49, t51, t52, t53, t54, t56, t57;

  t1 = xi*xi;
  t7 = 3.0/4.0*(35.0/2.0*t1*xi-15.0/2.0*xi)*eta;
  t8 = eta*eta;
  t14 = 3.0/4.0*xi*(35.0/2.0*t8*eta-15.0/2.0*eta);
  t15 = xi*eta;
  t16 = 9.0/4.0*t15;
  t17 = -t7+t14-t16;
  t18 = t7-t14-t16;
  t19 = eta*zeta;
  t20 = 3.0/4.0*t19;
  t21 = xi*zeta;
  t22 = 3.0/4.0*t21;
  t23 = 15.0/8.0*t1;
  t24 = 15.0/8.0*t8;
  t25 = zeta*zeta;
  t26 = 3.0/8.0*t25;
  t27 = 3.0/2.0*xi;
  t29 = 15.0/2.0*t8-3.0/2.0;
  t30 = xi*t29;
  t31 = 3.0/2.0*t30;
  t34 = 15.0/2.0*t1-3.0/2.0;
  t35 = t34*eta;
  t36 = 3.0/2.0*t35;
  t37 = 3.0/2.0*eta;
  t42 = t25*zeta;
  t43 = 5.0/6.0*t42;
  t44 = zeta/3.0;
  t45 = 3.0/4.0*t25;
  t46 = t34*zeta;
  t47 = t46/6.0;
  t48 = t29*zeta;
  t49 = t48/6.0;
  t51 = 5.0/12.0*t42;
  t52 = 5.0/12.0*zeta;
  t53 = 3.0/2.0*t19;
  t54 = t48/3.0;
  t56 = 3.0/2.0*t21;
  t57 = t46/3.0;

  values[0] = 0.0;
  values[1] = t17;
  values[2] = t18;
  values[3] = t17;
  values[4] = t18;
  values[5] = 0.0;
  values[6] = t20;
  values[7] = t22;
  values[8] = -t20;
  values[9] = -t22;
  values[10] = -t23-3.0/8.0+t24+t26+t27-t31;
  values[11] = -t22;
  values[12] = t22;
  values[13] = t23-3.0/8.0-t24-t36+t26+t37;
  values[14] = t23+3.0/8.0-t24-t26+t27-t31;
  values[15] = t23-3.0/8.0-t24+t36+t26-t37;
  values[16] = t20;
  values[17] = -t20;
  values[18] = 0.0;
  values[19] = 0.0;
  values[20] = 0.0;
  values[21] = 0.0;
  values[22] = -t36+t7-t14+t16;
  values[23] = 0.0;
  values[24] = -t43+t44-1.0/4.0+t45+t47+t49;
  values[25] = t51-t52+t53+t47-t54;
  values[26] = -t51+t52+t56+t57-t49;
  values[27] = t51-t52-t53+t47-t54;
  values[28] = t51-t52+t56-t57+t49;
  values[29] = t43-t44-1.0/4.0+t45-t47-t49;
  values[30] = 0.0;
  values[31] = -t7+t14-t31+t16;
  values[32] = t36+t7-t14+t16;
  values[33] = -t7+t14+t31+t16;
  values[34] = 0.0;
  values[35] = 0.0;
  values[36] = 9.0*t15;
  values[37] = 3.0*t35-3.0*eta;
  values[38] = -3.0*xi+3.0*t30;
  values[39] = 0.0;
}

static void N_H_Q3_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t7, t8, t14, t15, t16, t17, t18, t19, t20, t22, t23, t24, t25;
  double t26, t27, t29, t30, t32, t33, t34, t35, t37, t38, t43, t44, t45;
  double t46, t47, t48, t49, t50, t52, t53, t54, t55, t57, t58;

  t1 = xi*xi;
  t7 = 3.0/4.0*(35.0/2.0*t1*xi-15.0/2.0*xi)*zeta;
  t8 = zeta*zeta;
  t14 = 3.0/4.0*xi*(35.0/2.0*t8*zeta-15.0/2.0*zeta);
  t15 = xi*zeta;
  t16 = 9.0/4.0*t15;
  t17 = -t7+t14-t16;
  t18 = t7-t14-t16;
  t19 = 15.0/8.0*t1;
  t20 = 15.0/8.0*t8;
  t22 = 15.0/2.0*t1-3.0/2.0;
  t23 = t22*zeta;
  t24 = 3.0/2.0*t23;
  t25 = eta*eta;
  t26 = 3.0/8.0*t25;
  t27 = 3.0/2.0*zeta;
  t29 = xi*eta;
  t30 = 3.0/4.0*t29;
  t32 = 15.0/2.0*t8-3.0/2.0;
  t33 = xi*t32;
  t34 = 3.0/2.0*t33;
  t35 = 3.0/2.0*xi;
  t37 = eta*zeta;
  t38 = 3.0/4.0*t37;
  t43 = eta*t32;
  t44 = t43/3.0;
  t45 = 3.0/2.0*t37;
  t46 = t22*eta;
  t47 = t46/6.0;
  t48 = t25*eta;
  t49 = 5.0/12.0*t48;
  t50 = 5.0/12.0*eta;
  t52 = t43/6.0;
  t53 = 3.0/4.0*t25;
  t54 = 5.0/6.0*t48;
  t55 = eta/3.0;
  t57 = 3.0/2.0*t29;
  t58 = t46/3.0;

  values[0] = t17;
  values[1] = 0.0;
  values[2] = t18;
  values[3] = 0.0;
  values[4] = t18;
  values[5] = t17;
  values[6] = t19-3.0/8.0-t20-t24+t26+t27;
  values[7] = t30;
  values[8] = t19+3.0/8.0-t20-t34-t26+t35;
  values[9] = -t30;
  values[10] = t38;
  values[11] = -t30;
  values[12] = t30;
  values[13] = t38;
  values[14] = -t38;
  values[15] = t38;
  values[16] = -t19-3.0/8.0+t20-t34+t26+t35;
  values[17] = -t19+3.0/8.0+t20-t24-t26+t27;
  values[18] = -t34-t7+t14+t16;
  values[19] = 0.0;
  values[20] = t24+t7-t14+t16;
  values[21] = 0.0;
  values[22] = 0.0;
  values[23] = 0.0;
  values[24] = -t44+t45+t47+t49-t50;
  values[25] = t52-1.0/4.0+t53+t47-t54+t55;
  values[26] = -t52+t57+t58-t49+t50;
  values[27] = t52+1.0/4.0-t53+t47-t54+t55;
  values[28] = t52+t57-t58+t49-t50;
  values[29] = t44+t45-t47-t49+t50;
  values[30] = 0.0;
  values[31] = 0.0;
  values[32] = 0.0;
  values[33] = 0.0;
  values[34] = -t24+t7-t14+t16;
  values[35] = t34-t7+t14+t16;
  values[36] = 9.0*t15;
  values[37] = 3.0*t23-3.0*zeta;
  values[38] = 0.0;
  values[39] = 3.0*t33-3.0*xi;
}

static void N_H_Q3_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t7, t9, t11, t12, t13, t14, t15, t16, t17;
  double t18, t20, t22, t24, t25, t27, t28, t29, t30, t33, t34, t35, t36;
  double t38, t39, t40, t43, t49, t50, t52, t53, t55, t56, t57, t58;

  t1 = zeta*zeta;
  t2 = 63.0/16.0*t1;
  t3 = t1*t1;
  t4 = 105.0/32.0*t3;
  t5 = eta*eta;
  t7 = 105.0/2.0*t5-15.0/2.0;
  t9 = -1.0/2.0+3.0/2.0*t1;
  t11 = t7*t9/4.0;
  t12 = 21.0/32.0-t2+t4-t11;
  t13 = 15.0/2.0*eta;
  t14 = 27.0/16.0*t1;
  t15 = xi*xi;
  t16 = t15*t15;
  t17 = 105.0/32.0*t16;
  t18 = 27.0/16.0*t15;
  t20 = -1.0/2.0+3.0/2.0*t15;
  t22 = t20*t7/4.0;
  t24 = 63.0/16.0*t15;
  t25 = 21.0/32.0+t17-t24-t22;
  t27 = xi*zeta;
  t28 = 3.0/4.0*t27;
  t29 = t1*zeta;
  t30 = eta*zeta;
  t33 = xi*eta;
  t34 = 15.0/4.0*t33;
  t35 = t20*eta;
  t36 = 15.0/2.0*t35;
  t38 = eta*t9;
  t39 = 15.0/2.0*t38;
  t40 = 15.0/4.0*t30;
  t43 = t15*xi;
  t49 = 15.0/4.0*t43;
  t50 = 9.0/4.0*xi;
  t52 = 15.0/4.0*t29;
  t53 = 9.0/4.0*zeta;
  t55 = t33*zeta;
  t56 = 5.0/2.0*t55;
  t57 = 3.0/2.0*t27;
  t58 = 5.0*t55;

  values[0] = t12;
  values[1] = -t13+27.0/16.0+t14-t4-t17+t18+t22+t11;
  values[2] = t25;
  values[3] = t13+27.0/16.0+t14-t4-t17+t18+t22+t11;
  values[4] = t25;
  values[5] = t12;
  values[6] = t28;
  values[7] = -15.0/4.0*t29+15.0/4.0*zeta-15.0/4.0*t30;
  values[8] = -t28;
  values[9] = -15.0/4.0*t29+15.0/4.0*zeta+15.0/4.0*t30;
  values[10] = t34-t36;
  values[11] = -t39-t40;
  values[12] = -t39+t40;
  values[13] = -15.0/4.0*t33-15.0/4.0*t43+15.0/4.0*xi;
  values[14] = -t34-t36;
  values[15] = -15.0/4.0*t33+15.0/4.0*t43-15.0/4.0*xi;
  values[16] = t28;
  values[17] = -t28;
  values[18] = 0.0;
  values[19] = -21.0/32.0+t2-t4-t39+t11;
  values[20] = 0.0;
  values[21] = -21.0/32.0+t2-t4+t39+t11;
  values[22] = -t49+t50-3.0/32.0+t17-t18-t22;
  values[23] = -3.0/32.0-t14+t4+t52-t53-t11;
  values[24] = t56;
  values[25] = t57-t58;
  values[26] = -t56;
  values[27] = -t57-t58;
  values[28] = t56;
  values[29] = -t56;
  values[30] = -3.0/32.0-t14+t4-t52+t53-t11;
  values[31] = -21.0/32.0-t17+t24+t22-t36;
  values[32] = t49-t50-3.0/32.0+t17-t18-t22;
  values[33] = -21.0/32.0-t17+t24+t22+t36;
  values[34] = 0.0;
  values[35] = 0.0;
  values[36] = -6.0+9.0/2.0*t1+9.0/2.0*t15;
  values[37] = 15.0/2.0*t43-15.0/2.0*xi;
  values[38] = -15.0*eta+15.0*t38+15.0*t35;
  values[39] = 15.0/2.0*t29-15.0/2.0*zeta;
}

static void N_H_Q3_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t9, t10, t16, t17, t18, t19, t20, t21, t22, t24, t25;
  double t26, t27, t28, t29, t32, t33, t34, t36, t37, t38, t44, t45, t46;
  double t47, t48, t49, t50, t51, t53, t54, t55, t57, t58, t59;

  t1 = eta*zeta;
  t2 = 9.0/4.0*t1;
  t3 = zeta*zeta;
  t9 = 3.0/4.0*eta*(35.0/2.0*t3*zeta-15.0/2.0*zeta);
  t10 = eta*eta;
  t16 = 3.0/4.0*(35.0/2.0*t10*eta-15.0/2.0*eta)*zeta;
  t17 = -t2+t9-t16;
  t18 = -t2-t9+t16;
  t19 = xi*eta;
  t20 = 3.0/4.0*t19;
  t21 = xi*xi;
  t22 = 3.0/8.0*t21;
  t24 = 15.0/2.0*t3-3.0/2.0;
  t25 = eta*t24;
  t26 = 3.0/2.0*t25;
  t27 = 3.0/2.0*eta;
  t28 = 15.0/8.0*t3;
  t29 = 15.0/8.0*t10;
  t32 = xi*zeta;
  t33 = 3.0/4.0*t32;
  t34 = 3.0/2.0*zeta;
  t36 = 15.0/2.0*t10-3.0/2.0;
  t37 = t36*zeta;
  t38 = 3.0/2.0*t37;
  t44 = xi*t24;
  t45 = t44/3.0;
  t46 = 3.0/2.0*t32;
  t47 = t21*xi;
  t48 = 5.0/12.0*t47;
  t49 = 5.0/12.0*xi;
  t50 = xi*t36;
  t51 = t50/6.0;
  t53 = t44/6.0;
  t54 = 3.0/2.0*t19;
  t55 = t50/3.0;
  t57 = 3.0/4.0*t21;
  t58 = 5.0/6.0*t47;
  t59 = xi/3.0;

  values[0] = t17;
  values[1] = t18;
  values[2] = 0.0;
  values[3] = t18;
  values[4] = 0.0;
  values[5] = t17;
  values[6] = t20;
  values[7] = -3.0/8.0+t22-t26+t27+t28-t29;
  values[8] = -t20;
  values[9] = 3.0/8.0-t22-t26+t27-t28+t29;
  values[10] = t33;
  values[11] = t34+3.0/8.0-t22-t38+t28-t29;
  values[12] = t34-3.0/8.0+t22-t38-t28+t29;
  values[13] = t33;
  values[14] = -t33;
  values[15] = t33;
  values[16] = t20;
  values[17] = -t20;
  values[18] = 0.0;
  values[19] = t2-t9-t38+t16;
  values[20] = 0.0;
  values[21] = t2-t9+t38+t16;
  values[22] = 0.0;
  values[23] = t2+t9+t26-t16;
  values[24] = -t45+t46+t48-t49+t51;
  values[25] = t53+t54+t48-t49-t55;
  values[26] = -t53-1.0/4.0+t57+t58-t59-t51;
  values[27] = t53-t54+t48-t49-t55;
  values[28] = t53-1.0/4.0+t57-t58+t59+t51;
  values[29] = t45+t46-t48+t49-t51;
  values[30] = t2+t9-t26-t16;
  values[31] = 0.0;
  values[32] = 0.0;
  values[33] = 0.0;
  values[34] = 0.0;
  values[35] = 0.0;
  values[36] = 9.0*t1;
  values[37] = 0.0;
  values[38] = -3.0*zeta+3.0*t37;
  values[39] = 3.0*t25-3.0*eta;
}

static void N_H_Q3_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t12, t14, t16, t18, t19, t20;
  double t22, t23, t24, t25, t27, t28, t30, t31, t32, t33, t35, t36, t37;
  double t40, t41, t42, t48, t49, t51, t52, t56, t57, t58, t60;

  t1 = eta*eta;
  t2 = 27.0/16.0*t1;
  t3 = 15.0/2.0*zeta;
  t4 = xi*xi;
  t5 = t4*t4;
  t6 = 105.0/32.0*t5;
  t7 = 27.0/16.0*t4;
  t9 = -1.0/2.0+3.0/2.0*t4;
  t10 = zeta*zeta;
  t12 = 105.0/2.0*t10-15.0/2.0;
  t14 = t9*t12/4.0;
  t16 = -1.0/2.0+3.0/2.0*t1;
  t18 = t16*t12/4.0;
  t19 = t1*t1;
  t20 = 105.0/32.0*t19;
  t22 = 63.0/16.0*t1;
  t23 = 21.0/32.0-t22-t18+t20;
  t24 = 63.0/16.0*t4;
  t25 = 21.0/32.0+t6-t24-t14;
  t27 = xi*zeta;
  t28 = t4*xi;
  t30 = t16*zeta;
  t31 = 15.0/2.0*t30;
  t32 = eta*zeta;
  t33 = 15.0/4.0*t32;
  t35 = 15.0/4.0*t27;
  t36 = t9*zeta;
  t37 = 15.0/2.0*t36;
  t40 = xi*eta;
  t41 = 3.0/4.0*t40;
  t42 = t1*eta;
  t48 = 15.0/4.0*t42;
  t49 = 9.0/4.0*eta;
  t51 = 15.0/4.0*t28;
  t52 = 9.0/4.0*xi;
  t56 = t40*zeta;
  t57 = 5.0*t56;
  t58 = 3.0/2.0*t40;
  t60 = 5.0/2.0*t56;

  values[0] = 27.0/16.0+t2-t3-t6+t7+t14+t18-t20;
  values[1] = t23;
  values[2] = t25;
  values[3] = t23;
  values[4] = t25;
  values[5] = 27.0/16.0+t2+t3-t6+t7+t14+t18-t20;
  values[6] = -15.0/4.0*t27-15.0/4.0*t28+15.0/4.0*xi;
  values[7] = -t31+t33;
  values[8] = -t35-t37;
  values[9] = -t31-t33;
  values[10] = t41;
  values[11] = 15.0/4.0*eta-15.0/4.0*t42+15.0/4.0*t32;
  values[12] = 15.0/4.0*eta-15.0/4.0*t42-15.0/4.0*t32;
  values[13] = t41;
  values[14] = -t41;
  values[15] = t41;
  values[16] = t35-t37;
  values[17] = 15.0/4.0*t27-15.0/4.0*t28+15.0/4.0*xi;
  values[18] = -t37-21.0/32.0-t6+t24+t14;
  values[19] = -3.0/32.0-t2-t18-t48+t49+t20;
  values[20] = t51-t52-3.0/32.0+t6-t7-t14;
  values[21] = -3.0/32.0-t2-t18+t48-t49+t20;
  values[22] = 0.0;
  values[23] = -21.0/32.0+t22+t18+t31-t20;
  values[24] = -t57+t58;
  values[25] = t60;
  values[26] = -t60;
  values[27] = t60;
  values[28] = t60;
  values[29] = t57+t58;
  values[30] = -21.0/32.0+t22+t18-t31-t20;
  values[31] = 0.0;
  values[32] = 0.0;
  values[33] = 0.0;
  values[34] = -t51+t52-3.0/32.0+t6-t7-t14;
  values[35] = t37-21.0/32.0-t6+t24+t14;
  values[36] = -6.0+9.0/2.0*t1+9.0/2.0*t4;
  values[37] = 15.0/2.0*t28-15.0/2.0*xi;
  values[38] = -15.0/2.0*eta+15.0/2.0*t42;
  values[39] = 15.0*t36-15.0*zeta+15.0*t30;
}

static int N_H_Q3_3D_ChangeU0[2] = {  6, 24 };
static int N_H_Q3_3D_ChangeU1[2] = {  7, 25 };
static int N_H_Q3_3D_ChangeU2[2] = {  8, 26 };
static int N_H_Q3_3D_ChangeU3[2] = {  9, 27 };
static int N_H_Q3_3D_ChangeU4[2] = { 10, 28 };
static int N_H_Q3_3D_ChangeU5[2] = { 11, 29 };

static int N_H_Q3_3D_ChangeV0[2] = { 12, 24 };
static int N_H_Q3_3D_ChangeV1[2] = { 13, 25 };
static int N_H_Q3_3D_ChangeV2[2] = { 14, 26 };
static int N_H_Q3_3D_ChangeV3[2] = { 15, 27 };
static int N_H_Q3_3D_ChangeV4[2] = { 16, 28 };
static int N_H_Q3_3D_ChangeV5[2] = { 17, 29 };

static int *N_H_Q3_3D_ChangeU[6] = { N_H_Q3_3D_ChangeU0, N_H_Q3_3D_ChangeU1,
                                     N_H_Q3_3D_ChangeU2, N_H_Q3_3D_ChangeU3,
                                     N_H_Q3_3D_ChangeU4, N_H_Q3_3D_ChangeU5 };

static int *N_H_Q3_3D_ChangeV[6] = { N_H_Q3_3D_ChangeV0, N_H_Q3_3D_ChangeV1,
                                     N_H_Q3_3D_ChangeV2, N_H_Q3_3D_ChangeV3,
                                     N_H_Q3_3D_ChangeV4, N_H_Q3_3D_ChangeV5 };

static int **N_H_Q3_3D_Change[2] = { N_H_Q3_3D_ChangeU, N_H_Q3_3D_ChangeV };

TBaseFunct3D *BF_N_H_Q3_3D_Obj = 
new TBaseFunct3D(40, BF_N_H_Q3_3D, BFUnitHexahedron, 
                 N_H_Q3_3D_Funct, N_H_Q3_3D_DeriveXi,
                 N_H_Q3_3D_DeriveEta, N_H_Q3_3D_DeriveZeta,
                 N_H_Q3_3D_DeriveXiXi, N_H_Q3_3D_DeriveXiEta,
                 N_H_Q3_3D_DeriveXiZeta, N_H_Q3_3D_DeriveEtaEta,
                 N_H_Q3_3D_DeriveEtaZeta, N_H_Q3_3D_DeriveZetaZeta,
                 4, 3,
                 2, N_H_Q3_3D_Change);
