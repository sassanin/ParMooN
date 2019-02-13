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
// P6 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_P6_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t4, t6, t7, t8, t10, t12, t14, t15, t18, t20, t22, t24;
  double t26, t28, t29, t33, t34, t35, t36, t37, t40, t41, t42, t43, t44;
  double t46, t48, t50, t51, t52, t53, t54, t57, t58, t62, t63, t64, t66;
  double t74, t76, t77, t80, t81, t82, t86, t87, t95, t97, t101, t102, t103;
  double t111, t112, t113, t116, t117, t123, t125, t126, t128, t129, t132;
  double t135, t136, t143, t148, t149, t151, t159, t165, t166, t173, t174;
  double t175, t180, t184, t186, t189, t190, t194, t195, t196, t198, t199;
  double t203, t207, t210, t213, t214, t222, t231, t242, t244, t248, t249;
  double t250, t258, t263, t266, t271, t273, t278, t287, t289, t294, t307;
  double t308, t310;

  t1 = xi*xi;
  t2 = t1*t1;
  t4 = t2*xi;
  t6 = eta*eta;
  t7 = t6*t6;
  t8 = t7*eta;
  t10 = t1*xi;
  t12 = t6*eta;
  t14 = t10*t12;
  t15 = 1296.0*t14;
  t18 = xi*t12;
  t20 = t1*eta;
  t22 = xi*eta;
  t24 = t10*eta;
  t26 = t1*t6;
  t28 = 1.0+315.0*t2-1134.0/5.0*t4-1134.0/5.0*t8-441.0/2.0*t10-441.0/2.0*t12+t15-147.0/10.0*xi-147.0/10.0*eta+1260.0*t18-1323.0/2.0*t20+812.0/5.0*t22+1260.0*t24+1890.0*t26;
  t29 = xi*t7;
  t33 = t2*t1;
  t34 = 324.0/5.0*t33;
  t35 = t7*t6;
  t36 = 324.0/5.0*t35;
  t37 = t2*eta;
  t40 = t4*eta;
  t41 = 1944.0/5.0*t40;
  t42 = xi*t8;
  t43 = 1944.0/5.0*t42;
  t44 = xi*t6;
  t46 = t10*t6;
  t48 = t1*t12;
  t50 = t2*t6;
  t51 = 972.0*t50;
  t52 = t1*t7;
  t53 = 972.0*t52;
  t54 = -1134.0*t29+406.0/5.0*t6+406.0/5.0*t1+t34+t36-1134.0*t37+315.0*t7+t41+t43-1323.0/2.0*t44-2268.0*t46-2268.0*t48+t51+t53;
  t57 = 1944.0/5.0*t33;
  t58 = 1944.0*t40;
  t62 = 5184.0*t48;
  t63 = 3888.0*t50;
  t64 = 1944.0*t52;
  t66 = 1566.0/5.0*t22;
  t74 = 5022.0*t26;
  t76 = 3888.0*t14;
  t77 = -t66-1674.0*t2+1044.0*t10+2088.0*t20-1566.0/5.0*t1+1296.0*t4+1044.0*t44-5022.0*t24-t74-1674.0*t18-t76;
  t80 = 972.0*t33;
  t81 = 3888.0*t40;
  t82 = 162.0*t29;
  t86 = 5832.0*t50;
  t87 = 513.0/2.0*t22;
  t95 = 4671.0*t26;
  t97 = -45.0*xi+t80+t81-t82-9396.0*t37-9720.0*t46-3564.0*t48+t86+t53+t87+3699.0*t2-4149.0/2.0*t10-2610.0*t20+1053.0/2.0*t1-3078.0*t4-1071.0/2.0*t44+7884.0*t24+t95+486.0*t18+t76;
  t101 = 5184.0*t46;
  t102 = 648.0*t48;
  t103 = 148.0*t22;
  t111 = 1836.0*t26;
  t112 = 72.0*t18;
  t113 = 40.0*xi-1296.0*t33-t81+8424.0*t37+t101+t102-t63-t103-4356.0*t2+2232.0*t10+1692.0*t20-508.0*t1+3888.0*t4+180.0*t44-6120.0*t24-t111-t112-t15;
  t116 = 972.0*t46;
  t117 = 99.0/2.0*t22;
  t123 = 27.0*t44;
  t125 = 297.0*t26;
  t126 = -45.0/2.0*xi+t80+t58-3726.0*t37-t116+t51+t117+2889.0*t2-2763.0/2.0*t10-1197.0/2.0*t20+297.0*t1-2754.0*t4-t123+2376.0*t24+t125;
  t128 = 648.0*t37;
  t129 = 36.0/5.0*t22;
  t132 = 90.0*t20;
  t135 = 378.0*t24;
  t136 = 36.0/5.0*xi-t57-t41+t128-t129-1026.0*t2+468.0*t10+t132-486.0/5.0*t1+5184.0/5.0*t4-t135;
  t143 = 1944.0/5.0*t35;
  t148 = 1944.0*t50;
  t149 = 3888.0*t52;
  t151 = 1944.0*t42;
  t159 = -t151-t66-1566.0/5.0*t6+1044.0*t20-1674.0*t7+2088.0*t44-1674.0*t24-t74-5022.0*t18+1044.0*t12-t76;
  t165 = 7776.0*t50;
  t166 = 7776.0*t52;
  t173 = 11664.0*t14;
  t174 = t58-5832.0*t29-5832.0*t37-17496.0*t46-17496.0*t48+t165+t166+t151+540.0*t22-3078.0*t20-3078.0*t44+6426.0*t24+12852.0*t26+6426.0*t18+t173;
  t175 = 648.0*t29;
  t180 = 360.0*t22;
  t184 = 11232.0*t26;
  t186 = -t81+t175+10368.0*t37+21384.0*t46+11664.0*t48-11664.0*t50-t149-t180+3492.0*t20+1332.0*t44-9612.0*t24-t184-1620.0*t18-t173;
  t189 = 1944.0*t48;
  t190 = 180.0*t22;
  t194 = 3996.0*t26;
  t195 = 216.0*t18;
  t196 = t81-9072.0*t37-11016.0*t46-t189+t165+t190-2016.0*t20-396.0*t44+7020.0*t24+t194+t195+t76;
  t198 = 1944.0*t46;
  t199 = 54.0*t22;
  t203 = 594.0*t26;
  t207 = 972.0*t35;
  t210 = 162.0*t37;
  t213 = 5832.0*t52;
  t214 = 3888.0*t42;
  t222 = -45.0*eta+t207-9396.0*t29-3078.0*t8-t210-3564.0*t46-9720.0*t48+t51+t213+t214+t87+1053.0/2.0*t6-1071.0/2.0*t20+3699.0*t7-2610.0*t44+486.0*t24+t95+7884.0*t18-4149.0/2.0*t12+t76;
  t231 = 10368.0*t29+t128+11664.0*t46+21384.0*t48-t63-11664.0*t52-t214-t180+1332.0*t20+3492.0*t44-1620.0*t24-t184-9612.0*t18-t173;
  t242 = -972.0*t29-972.0*t37-13608.0*t46-13608.0*t48+t86+t213+135.0*t22-1107.0*t20-1107.0*t44+1944.0*t24+8748.0*t26+1944.0*t18+t173;
  t244 = 36.0*t22;
  t248 = 2484.0*t26;
  t249 = t128+6480.0*t46+t189-t63-t244+360.0*t20+252.0*t44-972.0*t24-t248-t195-t76;
  t250 = 9.0/2.0*t22;
  t258 = 648.0*t46;
  t263 = 72.0*t24;
  t266 = 40.0*eta-1296.0*t35+8424.0*t29+3888.0*t8+t258+t62-t149-t214-t103-508.0*t6+180.0*t20-4356.0*t7+1692.0*t44-t263-t111-6120.0*t18+2232.0*t12-t15;
  t271 = 216.0*t24;
  t273 = -9072.0*t29-t198-11016.0*t48+t166+t214+t190-396.0*t20-2016.0*t44+t271+t194+7020.0*t18+t76;
  t278 = t175+t198+6480.0*t48-t149-t244+252.0*t20+360.0*t44-t271-t248-972.0*t18-t76;
  t287 = 972.0*t48;
  t289 = 27.0*t20;
  t294 = -45.0/2.0*eta+t207-3726.0*t29-2754.0*t8-t287+t53+t151+t117+297.0*t6-t289+2889.0*t7-1197.0/2.0*t44+t125+2376.0*t18-2763.0/2.0*t12;
  t307 = 90.0*t44;
  t308 = 378.0*t18;
  t310 = 36.0/5.0*eta-t143+t175+5184.0/5.0*t8-t43-t129-486.0/5.0*t6-1026.0*t7+t307-t308+468.0*t12;

  values[0] = t28+t54;
  values[1] = 36.0*xi-t57-t58+1296.0*t29+5184.0*t37+7776.0*t46+t62-t63-t64-t43+t77;
  values[2] = t97;
  values[3] = t113;
  values[4] = t126;
  values[5] = t136;
  values[6] = -xi+t34+153.0*t2-135.0/2.0*t10+137.0/10.0*t1-162.0*t4;
  values[7] = 36.0*eta-t41-t143+5184.0*t29+1296.0*t8+1296.0*t37+t101+7776.0*t48-t148-t149+t159;
  values[8] = t174;
  values[9] = t186;
  values[10] = t196;
  values[11] = -t58+3888.0*t37+t198-t148-t199+648.0*t20+54.0*t44-2538.0*t24-t203;
  values[12] = t41-t128+t129-t132+t135;
  values[13] = t222;
  values[14] = t231;
  values[15] = t242;
  values[16] = t249;
  values[17] = -t210-t116+t51+t250-99.0/2.0*t20-t123+162.0*t24+t125;
  values[18] = t266;
  values[19] = t273;
  values[20] = t278;
  values[21] = -t258-t102+4.0*t22-36.0*t20-36.0*t44+t263+324.0*t26+t112+t15;
  values[22] = t294;
  values[23] = 3888.0*t29+t189-t64-t151-t199+54.0*t20+648.0*t44-t203-2538.0*t18;
  values[24] = -t82-t287+t53+t250-t289-99.0/2.0*t44+t125+162.0*t18;
  values[25] = t310;
  values[26] = -t175+t43+t129-t307+t308;
  values[27] = -eta+t36-162.0*t8+137.0/10.0*t6+153.0*t7-135.0/2.0*t12;
}

// values of the derivatives in xi direction
static void C_T_P6_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t13, t15, t16, t18, t19;
  double t20, t21, t22, t23, t28, t33, t35, t38, t39, t40, t42, t43, t47;
  double t48, t49, t51, t59, t61, t62, t64, t65, t66, t70, t71, t79, t81;
  double t84, t85, t86, t94, t95, t96, t98, t99, t105, t107, t108, t109;
  double t110, t113, t116, t117, t126, t127, t128, t133, t138, t139, t146;
  double t147, t148, t153, t157, t159, t162, t163, t167, t168, t169, t171;
  double t172, t176, t180, t183, t184, t189, t198, t209, t211, t215, t216;
  double t217, t222, t225, t227, t232, t234, t239, t246, t247, t261;

  t1 = xi*xi;
  t2 = t1*t1;
  t3 = t2*xi;
  t4 = 1944.0/5.0*t3;
  t5 = t2*eta;
  t6 = 1944.0*t5;
  t7 = eta*eta;
  t8 = t7*t7;
  t10 = t1*xi;
  t11 = t10*eta;
  t13 = t1*t7;
  t15 = t7*eta;
  t16 = xi*t15;
  t18 = t10*t7;
  t19 = 3888.0*t18;
  t20 = xi*t8;
  t21 = 1944.0*t20;
  t22 = t8*eta;
  t23 = 1944.0/5.0*t22;
  t28 = xi*eta;
  t33 = t1*eta;
  t35 = xi*t7;
  t38 = t1*t15;
  t39 = 3888.0*t38;
  t40 = 812.0/5.0*eta+1260.0*t10-1323.0/2.0*t1-1323.0*t28+812.0/5.0*xi-1134.0*t2-1323.0/2.0*t7+3780.0*t33+3780.0*t35+1260.0*t15+t39;
  t42 = 11664.0/5.0*t3;
  t43 = 9720.0*t5;
  t47 = 10368.0*t16;
  t48 = 15552.0*t18;
  t49 = 3888.0*t20;
  t51 = 1566.0/5.0*eta;
  t59 = 10044.0*t35;
  t61 = 11664.0*t38;
  t62 = -t51-6696.0*t10+3132.0*t1+4176.0*t28-3132.0/5.0*xi+6480.0*t2+1044.0*t7-15066.0*t33-t59-1674.0*t15-t61;
  t64 = 5832.0*t3;
  t65 = 19440.0*t5;
  t66 = 162.0*t8;
  t70 = 23328.0*t18;
  t71 = 513.0/2.0*eta;
  t79 = 9342.0*t35;
  t81 = -45.0+t64+t65-t66-37584.0*t11-29160.0*t13-7128.0*t16+t70+t21+t71+14796.0*t10-12447.0/2.0*t1-5220.0*t28+1053.0*xi-15390.0*t2-1071.0/2.0*t7+23652.0*t33+t79+486.0*t15+t61;
  t84 = 15552.0*t13;
  t85 = 1296.0*t16;
  t86 = 148.0*eta;
  t94 = 3672.0*t35;
  t95 = 72.0*t15;
  t96 = 40.0-7776.0*t3-t65+33696.0*t11+t84+t85-t48-t86-17424.0*t10+6696.0*t1+3384.0*t28-1016.0*xi+19440.0*t2+180.0*t7-18360.0*t33-t94-t95-t39;
  t98 = 2916.0*t13;
  t99 = 99.0/2.0*eta;
  t105 = 27.0*t7;
  t107 = 594.0*t35;
  t108 = -45.0/2.0+t64+t43-14904.0*t11-t98+t19+t99+11556.0*t10-8289.0/2.0*t1-1197.0*t28+594.0*xi-13770.0*t2-t105+7128.0*t33+t107;
  t109 = 2592.0*t11;
  t110 = 36.0/5.0*eta;
  t113 = 180.0*t28;
  t116 = 1134.0*t33;
  t117 = 36.0/5.0-t42-t6+t109-t110-4104.0*t10+1404.0*t1+t113-972.0/5.0*xi+5184.0*t2-t116;
  t126 = 7776.0*t18;
  t127 = 7776.0*t20;
  t128 = 1944.0*t22;
  t133 = -t6+5184.0*t8+5184.0*t11+t84+15552.0*t16-t126-t127-t128-t51+2088.0*t28+2088.0*t7-5022.0*t33-t59-5022.0*t15-t61;
  t138 = 31104.0*t18;
  t139 = 15552.0*t20;
  t146 = 34992.0*t38;
  t147 = t43-5832.0*t8-23328.0*t11-52488.0*t13-34992.0*t16+t138+t139+t128+540.0*eta-6156.0*t28-3078.0*t7+19278.0*t33+25704.0*t35+6426.0*t15+t146;
  t148 = 648.0*t8;
  t153 = 360.0*eta;
  t157 = 22464.0*t35;
  t159 = -t65+t148+41472.0*t11+64152.0*t13+23328.0*t16-46656.0*t18-t127-t153+6984.0*t28+1332.0*t7-28836.0*t33-t157-1620.0*t15-t146;
  t162 = 3888.0*t16;
  t163 = 180.0*eta;
  t167 = 7992.0*t35;
  t168 = 216.0*t15;
  t169 = t65-36288.0*t11-33048.0*t13-t162+t138+t163-4032.0*t28-396.0*t7+21060.0*t33+t167+t168+t61;
  t171 = 5832.0*t13;
  t172 = 54.0*eta;
  t176 = 1188.0*t35;
  t180 = 648.0*t11;
  t183 = 11664.0*t20;
  t184 = 3888.0*t22;
  t189 = -9396.0*t8-t180-10692.0*t13-19440.0*t16+t19+t183+t184+t71-1071.0*t28-2610.0*t7+1458.0*t33+t79+7884.0*t15+t61;
  t198 = 10368.0*t8+t109+34992.0*t13+42768.0*t16-t48-23328.0*t20-t184-t153+2664.0*t28+3492.0*t7-4860.0*t33-t157-9612.0*t15-t146;
  t209 = -972.0*t8-3888.0*t11-40824.0*t13-27216.0*t16+t70+t183+135.0*eta-2214.0*t28-1107.0*t7+5832.0*t33+17496.0*t35+1944.0*t15+t146;
  t211 = 36.0*eta;
  t215 = 4968.0*t35;
  t216 = t109+19440.0*t13+t162-t48-t211+720.0*t28+252.0*t7-2916.0*t33-t215-t168-t61;
  t217 = 9.0/2.0*eta;
  t222 = 1944.0*t13;
  t225 = 216.0*t33;
  t227 = 8424.0*t8+t222+t47-t127-t184-t86+360.0*t28+1692.0*t7-t225-t94-6120.0*t15-t39;
  t232 = 648.0*t33;
  t234 = -9072.0*t8-t171-22032.0*t16+t139+t184+t163-792.0*t28-2016.0*t7+t232+t167+7020.0*t15+t61;
  t239 = t148+t171+12960.0*t16-t127-t211+504.0*t28+360.0*t7-t232-t215-972.0*t15-t61;
  t246 = 1944.0*t16;
  t247 = 54.0*t28;
  t261 = t148-t23-t110+90.0*t7-378.0*t15;

  values[0] = -147.0/10.0+t4+t6-1134.0*t8-4536.0*t11-6804.0*t13-4536.0*t16+t19+t21+t23+t40;
  values[1] = 36.0-t42-t43+1296.0*t8+20736.0*t11+23328.0*t13+t47-t48-t49-t23+t62;
  values[2] = t81;
  values[3] = t96;
  values[4] = t108;
  values[5] = t117;
  values[6] = -1.0+t4+612.0*t10-405.0/2.0*t1+137.0/5.0*xi-810.0*t2;
  values[7] = t133;
  values[8] = t147;
  values[9] = t159;
  values[10] = t169;
  values[11] = -t43+15552.0*t11+t171-t126-t172+1296.0*t28+54.0*t7-7614.0*t33-t176;
  values[12] = t6-t109+t110-t113+t116;
  values[13] = t189;
  values[14] = t198;
  values[15] = t209;
  values[16] = t216;
  values[17] = -t180-t98+t19+t217-99.0*t28-t105+486.0*t33+t107;
  values[18] = t227;
  values[19] = t234;
  values[20] = t239;
  values[21] = -t222-t85+4.0*eta-72.0*t28-36.0*t7+t225+648.0*t35+t95+t39;
  values[22] = -3726.0*t8-t246+t21+t128+t99-t247-1197.0/2.0*t7+t107+2376.0*t15;
  values[23] = 3888.0*t8+t162-t49-t128-t172+108.0*t28+648.0*t7-t176-2538.0*t15;
  values[24] = -t66-t246+t21+t217-t247-99.0/2.0*t7+t107+162.0*t15;
  values[25] = t261;
  values[26] = -t261;
  values[27] = 0.0;
}

// values of the derivatives in eta direction
static void C_T_P6_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t13, t15, t16, t18, t19;
  double t20, t21, t22, t23, t28, t33, t35, t38, t39, t40, t42, t46, t47;
  double t48, t49, t53, t55, t56, t57, t58, t62, t63, t67, t69, t71, t72;
  double t73, t77, t78, t79, t81, t82, t84, t86, t88, t89, t92, t93, t98;
  double t99, t101, t109, t115, t116, t123, t124, t125, t130, t134, t136;
  double t139, t140, t144, t145, t146, t148, t149, t153, t155, t158, t161;
  double t162, t170, t179, t190, t192, t196, t197, t198, t205, t210, t213;
  double t218, t220, t225, t233, t235, t240, t252, t253, t255;

  t1 = xi*xi;
  t2 = t1*t1;
  t3 = t2*xi;
  t4 = 1944.0/5.0*t3;
  t5 = t2*eta;
  t6 = 1944.0*t5;
  t7 = eta*eta;
  t8 = t7*t7;
  t10 = t1*xi;
  t11 = t10*eta;
  t13 = t1*t7;
  t15 = t7*eta;
  t16 = xi*t15;
  t18 = t10*t7;
  t19 = 3888.0*t18;
  t20 = xi*t8;
  t21 = 1944.0*t20;
  t22 = t8*eta;
  t23 = 1944.0/5.0*t22;
  t28 = xi*eta;
  t33 = t1*eta;
  t35 = xi*t7;
  t38 = t1*t15;
  t39 = 3888.0*t38;
  t40 = 812.0/5.0*eta+1260.0*t10-1323.0/2.0*t1-1323.0*t28+812.0/5.0*xi-1134.0*t2-1323.0/2.0*t7+3780.0*t33+3780.0*t35+1260.0*t15+t39;
  t42 = 1944.0*t3;
  t46 = 15552.0*t13;
  t47 = 7776.0*t5;
  t48 = 7776.0*t38;
  t49 = 1566.0/5.0*xi;
  t53 = 10044.0*t33;
  t55 = 11664.0*t18;
  t56 = -t42+5184.0*t16+5184.0*t2+15552.0*t11+t46-t47-t48-t21-t49+2088.0*t1+2088.0*t28-5022.0*t10-t53-5022.0*t35-t55;
  t57 = 3888.0*t3;
  t58 = 648.0*t16;
  t62 = 11664.0*t5;
  t63 = 513.0/2.0*xi;
  t67 = 9342.0*t33;
  t69 = t57-t58-9396.0*t2-19440.0*t11-10692.0*t13+t62+t39+t63-2610.0*t1-1071.0*t28+7884.0*t10+t67+1458.0*t35+t55;
  t71 = 10368.0*t11;
  t72 = 1944.0*t13;
  t73 = 148.0*xi;
  t77 = 3672.0*t33;
  t78 = 216.0*t35;
  t79 = -t57+8424.0*t2+t71+t72-t47-t73+1692.0*t1+360.0*t28-6120.0*t10-t77-t78-t19;
  t81 = 1944.0*t11;
  t82 = 99.0/2.0*xi;
  t84 = 54.0*t28;
  t86 = 594.0*t33;
  t88 = 648.0*t2;
  t89 = 36.0/5.0*xi;
  t92 = -t4+t88-t89+90.0*t1-378.0*t10;
  t93 = 11664.0/5.0*t22;
  t98 = 3888.0*t5;
  t99 = 15552.0*t38;
  t101 = 9720.0*t20;
  t109 = -t101-t49-3132.0/5.0*eta+1044.0*t1-6696.0*t15+4176.0*t28-1674.0*t10-t53-15066.0*t35+3132.0*t7-t55;
  t115 = 15552.0*t5;
  t116 = 31104.0*t38;
  t123 = 34992.0*t18;
  t124 = t42-23328.0*t16-5832.0*t2-34992.0*t11-52488.0*t13+t115+t116+t101+540.0*xi-3078.0*t1-6156.0*t28+6426.0*t10+25704.0*t33+19278.0*t35+t123;
  t125 = 2592.0*t16;
  t130 = 360.0*xi;
  t134 = 22464.0*t33;
  t136 = -t57+t125+10368.0*t2+42768.0*t11+34992.0*t13-23328.0*t5-t99-t130+3492.0*t1+2664.0*t28-9612.0*t10-t134-4860.0*t35-t123;
  t139 = 5832.0*t13;
  t140 = 180.0*xi;
  t144 = 7992.0*t33;
  t145 = 648.0*t35;
  t146 = t57-9072.0*t2-22032.0*t11-t139+t115+t140-2016.0*t1-792.0*t28+7020.0*t10+t144+t145+t55;
  t148 = 3888.0*t11;
  t149 = 54.0*xi;
  t153 = 1188.0*t33;
  t155 = 5832.0*t22;
  t158 = 162.0*t2;
  t161 = 23328.0*t38;
  t162 = 19440.0*t20;
  t170 = -45.0+t155-37584.0*t16-15390.0*t8-t158-7128.0*t11-29160.0*t13+t6+t161+t162+t63+1053.0*eta-1071.0/2.0*t1+14796.0*t15-5220.0*t28+486.0*t10+t67+23652.0*t35-12447.0/2.0*t7+t55;
  t179 = 41472.0*t16+t88+23328.0*t11+64152.0*t13-t47-46656.0*t38-t162-t130+1332.0*t1+6984.0*t28-1620.0*t10-t134-28836.0*t35-t123;
  t190 = -3888.0*t16-972.0*t2-27216.0*t11-40824.0*t13+t62+t161+135.0*xi-1107.0*t1-2214.0*t28+1944.0*t10+17496.0*t33+5832.0*t35+t123;
  t192 = 36.0*xi;
  t196 = 4968.0*t33;
  t197 = t88+12960.0*t11+t139-t47-t192+360.0*t1+504.0*t28-972.0*t10-t196-t145-t55;
  t198 = 9.0/2.0*xi;
  t205 = 1296.0*t11;
  t210 = 72.0*t10;
  t213 = 40.0-7776.0*t22+33696.0*t16+19440.0*t8+t205+t46-t99-t162-t73-1016.0*eta+180.0*t1-17424.0*t15+3384.0*t28-t210-t77-18360.0*t35+6696.0*t7-t19;
  t218 = 216.0*t10;
  t220 = -36288.0*t16-t148-33048.0*t13+t116+t162+t140-396.0*t1-4032.0*t28+t218+t144+21060.0*t35+t55;
  t225 = t125+t148+19440.0*t13-t99-t192+252.0*t1+720.0*t28-t218-t196-2916.0*t35-t55;
  t233 = 2916.0*t13;
  t235 = 27.0*t1;
  t240 = -45.0/2.0+t155-14904.0*t16-13770.0*t8-t233+t39+t101+t82+594.0*eta-t235+11556.0*t15-1197.0*t28+t86+7128.0*t35-8289.0/2.0*t7;
  t252 = 180.0*t28;
  t253 = 1134.0*t35;
  t255 = 36.0/5.0-t93+t125+5184.0*t8-t21-t89-972.0/5.0*eta-4104.0*t15+t252-t253+1404.0*t7;

  values[0] = -147.0/10.0+t4+t6-1134.0*t8-4536.0*t11-6804.0*t13-4536.0*t16+t19+t21+t23+t40;
  values[1] = t56;
  values[2] = t69;
  values[3] = t79;
  values[4] = t42-3726.0*t2-t81+t6+t82-1197.0/2.0*t1-t84+2376.0*t10+t86;
  values[5] = t92;
  values[6] = 0.0;
  values[7] = 36.0-t4-t93+20736.0*t16+6480.0*t8+1296.0*t2+t71+23328.0*t13-t98-t99+t109;
  values[8] = t124;
  values[9] = t136;
  values[10] = t146;
  values[11] = -t42+3888.0*t2+t148-t98-t149+648.0*t1+108.0*t28-2538.0*t10-t153;
  values[12] = -t92;
  values[13] = t170;
  values[14] = t179;
  values[15] = t190;
  values[16] = t197;
  values[17] = -t158-t81+t6+t198-99.0/2.0*t1-t84+162.0*t10+t86;
  values[18] = t213;
  values[19] = t220;
  values[20] = t225;
  values[21] = -t205-t72+4.0*xi-36.0*t1-72.0*t28+t210+648.0*t33+t78+t19;
  values[22] = t240;
  values[23] = 15552.0*t16+t139-t48-t101-t149+54.0*t1+1296.0*t28-t153-7614.0*t35;
  values[24] = -t58-t233+t39+t198-t235-99.0*t28+t86+486.0*t35;
  values[25] = t255;
  values[26] = -t125+t21+t89-t252+t253;
  values[27] = -1.0+t23-810.0*t8+137.0/5.0*eta+612.0*t15-405.0/2.0*t7;
}

// values of the derivatives in xi-xi direction
static void C_T_P6_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t12, t14, t15, t16, t17, t22;
  double t25, t26, t27, t28, t29, t32, t33, t34, t40, t41, t42, t43, t44;
  double t48, t54, t55, t58, t59, t65, t66, t68, t74, t75, t76, t79, t81;
  double t89, t90, t97, t98, t102, t110, t114, t117, t120, t123, t126;
  double t129, t149, t154, t156, t160, t170;

  t1 = xi*xi;
  t2 = t1*t1;
  t3 = 1944.0*t2;
  t4 = t1*xi;
  t5 = t4*eta;
  t6 = 7776.0*t5;
  t7 = t1*eta;
  t9 = eta*eta;
  t10 = xi*t9;
  t12 = t9*eta;
  t14 = t1*t9;
  t15 = 11664.0*t14;
  t16 = t9*t9;
  t17 = 1944.0*t16;
  t22 = xi*eta;
  t25 = xi*t12;
  t26 = 7776.0*t25;
  t27 = t3+t6-13608.0*t7-13608.0*t10-4536.0*t12+t15+t17+3780.0*t1-1323.0*xi-1323.0*eta+812.0/5.0-4536.0*t4+7560.0*t22+3780.0*t9+t26;
  t28 = 11664.0*t2;
  t29 = 38880.0*t5;
  t32 = 10368.0*t12;
  t33 = 46656.0*t14;
  t34 = 3888.0*t16;
  t40 = 10044.0*t9;
  t41 = 23328.0*t25;
  t42 = -t28-t29+62208.0*t7+46656.0*t10+t32-t33-t34-20088.0*t1+6264.0*xi+4176.0*eta-3132.0/5.0+25920.0*t4-30132.0*t22-t40-t41;
  t43 = 29160.0*t2;
  t44 = 77760.0*t5;
  t48 = 69984.0*t14;
  t54 = 9342.0*t9;
  t55 = t43+t44-112752.0*t7-58320.0*t10-7128.0*t12+t48+t17+44388.0*t1-12447.0*xi-5220.0*eta+1053.0-61560.0*t4+47304.0*t22+t54+t41;
  t58 = 31104.0*t10;
  t59 = 1296.0*t12;
  t65 = 3672.0*t9;
  t66 = -38880.0*t2-t44+101088.0*t7+t58+t59-t33-52272.0*t1+13392.0*xi+3384.0*eta-1016.0+77760.0*t4-36720.0*t22-t65-t26;
  t68 = 5832.0*t10;
  t74 = 594.0*t9;
  t75 = t43+t29-44712.0*t7-t68+t15+34668.0*t1-8289.0*xi-1197.0*eta+594.0-55080.0*t4+14256.0*t22+t74;
  t76 = 7776.0*t7;
  t79 = 180.0*eta;
  t81 = 2268.0*t22;
  t89 = 23328.0*t14;
  t90 = 7776.0*t16;
  t97 = 93312.0*t14;
  t98 = 15552.0*t16;
  t102 = 69984.0*t25;
  t110 = 22464.0*t9;
  t114 = 3888.0*t12;
  t117 = 7992.0*t9;
  t120 = 11664.0*t10;
  t123 = 1188.0*t9;
  t126 = 1944.0*t7;
  t129 = 11664.0*t16;
  t149 = 4968.0*t9;
  t154 = 3888.0*t10;
  t156 = 432.0*t22;
  t160 = 1296.0*t22;
  t170 = -1944.0*t12+t17-54.0*eta+t74;

  values[0] = t27;
  values[1] = t42;
  values[2] = t55;
  values[3] = t66;
  values[4] = t75;
  values[5] = -t28-t6+t76-12312.0*t1+2808.0*xi+t79-972.0/5.0+20736.0*t4-t81;
  values[6] = t3+1836.0*t1-405.0*xi+137.0/5.0-3240.0*t4;
  values[7] = -t6+15552.0*t7+t58+15552.0*t12-t89-t90+2088.0*eta-10044.0*t22-t40-t41;
  values[8] = t29-69984.0*t7-104976.0*t10-34992.0*t12+t97+t98-6156.0*eta+38556.0*t22+25704.0*t9+t102;
  values[9] = -t44+124416.0*t7+128304.0*t10+23328.0*t12-139968.0*t14-t90+6984.0*eta-57672.0*t22-t110-t102;
  values[10] = t44-108864.0*t7-66096.0*t10-t114+t97-4032.0*eta+42120.0*t22+t117+t41;
  values[11] = -t29+46656.0*t7+t120-t89+1296.0*eta-15228.0*t22-t123;
  values[12] = t6-t76-t79+t81;
  values[13] = -t126-21384.0*t10-19440.0*t12+t15+t129-1071.0*eta+2916.0*t22+t54+t41;
  values[14] = t76+69984.0*t10+42768.0*t12-t33-23328.0*t16+2664.0*eta-9720.0*t22-t110-t102;
  values[15] = -11664.0*t7-81648.0*t10-27216.0*t12+t48+t129-2214.0*eta+11664.0*t22+17496.0*t9+t102;
  values[16] = t76+38880.0*t10+t114-t33+720.0*eta-5832.0*t22-t149-t41;
  values[17] = -t126-t68+t15-99.0*eta+972.0*t22+t74;
  values[18] = t154+t32-t90+360.0*eta-t156-t65-t26;
  values[19] = -t120-22032.0*t12+t98-792.0*eta+t160+t117+t41;
  values[20] = t120+12960.0*t12-t90+504.0*eta-t160-t149-t41;
  values[21] = -t154-t59-72.0*eta+t156+648.0*t9+t26;
  values[22] = t170;
  values[23] = t114-t34+108.0*eta-t123;
  values[24] = t170;
  values[25] = 0.0;
  values[26] = 0.0;
  values[27] = 0.0;
}

// values of the derivatives in xi-eta direction
static void C_T_P6_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t12, t14, t15, t16, t17, t22;
  double t25, t26, t27, t28, t32, t33, t34, t38, t40, t41, t42, t43, t47;
  double t51, t53, t55, t56, t60, t61, t62, t64, t66, t68, t70, t73, t77;
  double t78, t79, t84, t89, t90, t96, t97, t98, t106, t108, t111, t115;
  double t116, t117, t119, t123, t126, t129, t130, t135, t144, t154, t159;
  double t160, t165, t168, t170, t175, t177, t182, t188, t189, t203;

  t1 = xi*xi;
  t2 = t1*t1;
  t3 = 1944.0*t2;
  t4 = t1*xi;
  t5 = t4*eta;
  t6 = 7776.0*t5;
  t7 = t1*eta;
  t9 = eta*eta;
  t10 = xi*t9;
  t12 = t9*eta;
  t14 = t1*t9;
  t15 = 11664.0*t14;
  t16 = t9*t9;
  t17 = 1944.0*t16;
  t22 = xi*eta;
  t25 = xi*t12;
  t26 = 7776.0*t25;
  t27 = t3+t6-13608.0*t7-13608.0*t10-4536.0*t12+t15+t17+3780.0*t1-1323.0*xi-1323.0*eta+812.0/5.0-4536.0*t4+7560.0*t22+3780.0*t9+t26;
  t28 = 9720.0*t2;
  t32 = 31104.0*t10;
  t33 = 31104.0*t5;
  t34 = 15552.0*t25;
  t38 = 20088.0*t22;
  t40 = 34992.0*t14;
  t41 = -t28+5184.0*t12+20736.0*t4+46656.0*t7+t32-t33-t34-t17-1566.0/5.0+4176.0*xi+2088.0*eta-15066.0*t1-t38-5022.0*t9-t40;
  t42 = 19440.0*t2;
  t43 = 648.0*t12;
  t47 = 46656.0*t5;
  t51 = 18684.0*t22;
  t53 = t42-t43-37584.0*t4-58320.0*t7-21384.0*t10+t47+t26+513.0/2.0-5220.0*xi-1071.0*eta+23652.0*t1+t51+1458.0*t9+t40;
  t55 = 31104.0*t7;
  t56 = 3888.0*t10;
  t60 = 7344.0*t22;
  t61 = 216.0*t9;
  t62 = -t42+33696.0*t4+t55+t56-t33-148.0+3384.0*xi+360.0*eta-18360.0*t1-t60-t61-t15;
  t64 = 5832.0*t7;
  t66 = 54.0*eta;
  t68 = 1188.0*t22;
  t70 = 2592.0*t4;
  t73 = -t3+t70-36.0/5.0+180.0*xi-1134.0*t1;
  t77 = 15552.0*t5;
  t78 = 31104.0*t25;
  t79 = 9720.0*t16;
  t84 = -t3+20736.0*t12+5184.0*t4+t55+46656.0*t10-t77-t78-t79-1566.0/5.0+2088.0*xi+4176.0*eta-5022.0*t1-t38-15066.0*t9-t40;
  t89 = 62208.0*t5;
  t90 = 62208.0*t25;
  t96 = 104976.0*t14;
  t97 = t28-23328.0*t12-23328.0*t4-104976.0*t7-104976.0*t10+t89+t90+t79+540.0-6156.0*xi-6156.0*eta+19278.0*t1+51408.0*t22+19278.0*t9+t96;
  t98 = 2592.0*t12;
  t106 = 44928.0*t22;
  t108 = -t42+t98+41472.0*t4+128304.0*t7+69984.0*t10-93312.0*t5-t78-360.0+6984.0*xi+2664.0*eta-28836.0*t1-t106-4860.0*t9-t96;
  t111 = 11664.0*t10;
  t115 = 15984.0*t22;
  t116 = 648.0*t9;
  t117 = t42-36288.0*t4-66096.0*t7-t111+t89+180.0-4032.0*xi-792.0*eta+21060.0*t1+t115+t116+t40;
  t119 = 11664.0*t7;
  t123 = 2376.0*t22;
  t126 = 648.0*t4;
  t129 = 46656.0*t25;
  t130 = 19440.0*t16;
  t135 = -37584.0*t12-t126-21384.0*t7-58320.0*t10+t6+t129+t130+513.0/2.0-1071.0*xi-5220.0*eta+1458.0*t1+t51+23652.0*t9+t40;
  t144 = 41472.0*t12+t70+69984.0*t7+128304.0*t10-t33-93312.0*t25-t130-360.0+2664.0*xi+6984.0*eta-4860.0*t1-t106-28836.0*t9-t96;
  t154 = -3888.0*t12-3888.0*t4-81648.0*t7-81648.0*t10+t47+t129+135.0-2214.0*xi-2214.0*eta+5832.0*t1+34992.0*t22+5832.0*t9+t96;
  t159 = 9936.0*t22;
  t160 = t70+38880.0*t7+t111-t33-36.0+720.0*xi+504.0*eta-2916.0*t1-t159-t116-t40;
  t165 = 3888.0*t7;
  t168 = 216.0*t1;
  t170 = 33696.0*t12+t165+t32-t78-t130-148.0+360.0*xi+3384.0*eta-t168-t60-18360.0*t9-t15;
  t175 = 648.0*t1;
  t177 = -36288.0*t12-t119-66096.0*t10+t90+t130+180.0-792.0*xi-4032.0*eta+t175+t115+21060.0*t9+t40;
  t182 = t98+t119+38880.0*t10-t78-36.0+504.0*xi+720.0*eta-t175-t159-2916.0*t9-t40;
  t188 = 5832.0*t10;
  t189 = 54.0*xi;
  t203 = t98-t17-36.0/5.0+180.0*eta-1134.0*t9;

  values[0] = t27;
  values[1] = t41;
  values[2] = t53;
  values[3] = t62;
  values[4] = t28-14904.0*t4-t64+t6+99.0/2.0-1197.0*xi-t66+7128.0*t1+t68;
  values[5] = t73;
  values[6] = 0.0;
  values[7] = t84;
  values[8] = t97;
  values[9] = t108;
  values[10] = t117;
  values[11] = -t28+15552.0*t4+t119-t77-54.0+1296.0*xi+108.0*eta-7614.0*t1-t123;
  values[12] = -t73;
  values[13] = t135;
  values[14] = t144;
  values[15] = t154;
  values[16] = t160;
  values[17] = -t126-t64+t6+9.0/2.0-99.0*xi-t66+486.0*t1+t68;
  values[18] = t170;
  values[19] = t177;
  values[20] = t182;
  values[21] = -t165-t56+4.0-72.0*xi-72.0*eta+t168+1296.0*t22+t61+t15;
  values[22] = -14904.0*t12-t188+t26+t79+99.0/2.0-t189-1197.0*eta+t68+7128.0*t9;
  values[23] = 15552.0*t12+t111-t34-t79-54.0+108.0*xi+1296.0*eta-t123-7614.0*t9;
  values[24] = -t43-t188+t26+9.0/2.0-t189-99.0*eta+t68+486.0*t9;
  values[25] = t203;
  values[26] = -t203;
  values[27] = 0.0;
}

// values of the derivatives in eta-eta direction
static void C_T_P6_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t12, t14, t15, t16, t17, t22;
  double t25, t26, t27, t30, t31, t32, t34, t36, t38, t41, t43, t46, t47;
  double t49, t50, t54, t55, t56, t60, t61, t62, t67, t71, t72, t76, t78;
  double t83, t87, t89, t90, t92, t94, t96, t101, t102, t107, t124, t129;
  double t134, t149, t154, t164, t165;

  t1 = xi*xi;
  t2 = t1*t1;
  t3 = 1944.0*t2;
  t4 = t1*xi;
  t5 = t4*eta;
  t6 = 7776.0*t5;
  t7 = t1*eta;
  t9 = eta*eta;
  t10 = xi*t9;
  t12 = t9*eta;
  t14 = t1*t9;
  t15 = 11664.0*t14;
  t16 = t9*t9;
  t17 = 1944.0*t16;
  t22 = xi*eta;
  t25 = xi*t12;
  t26 = 7776.0*t25;
  t27 = t3+t6-13608.0*t7-13608.0*t10-4536.0*t12+t15+t17+3780.0*t1-1323.0*xi-1323.0*eta+812.0/5.0-4536.0*t4+7560.0*t22+3780.0*t9+t26;
  t30 = 31104.0*t7;
  t31 = 7776.0*t2;
  t32 = 23328.0*t14;
  t34 = 10044.0*t1;
  t36 = 23328.0*t5;
  t38 = 1944.0*t10;
  t41 = 11664.0*t2;
  t43 = 9342.0*t1;
  t46 = 10368.0*t4;
  t47 = 3888.0*t7;
  t49 = 3672.0*t1;
  t50 = 432.0*t22;
  t54 = 594.0*t1;
  t55 = -1944.0*t4+t3-54.0*xi+t54;
  t56 = 11664.0*t16;
  t60 = 3888.0*t2;
  t61 = 46656.0*t14;
  t62 = 38880.0*t25;
  t67 = -t56+62208.0*t10+25920.0*t12+t46+46656.0*t7-t60-t61-t62-3132.0/5.0-20088.0*t9+4176.0*xi-t34-30132.0*t22+6264.0*eta-t36;
  t71 = 15552.0*t2;
  t72 = 93312.0*t14;
  t76 = 69984.0*t5;
  t78 = 7776.0*t10;
  t83 = 22464.0*t1;
  t87 = 11664.0*t7;
  t89 = 7992.0*t1;
  t90 = 1296.0*t22;
  t92 = 3888.0*t4;
  t94 = 1188.0*t1;
  t96 = 29160.0*t16;
  t101 = 69984.0*t14;
  t102 = 77760.0*t25;
  t107 = t96-112752.0*t10-61560.0*t12-7128.0*t4-58320.0*t7+t3+t101+t102+1053.0+44388.0*t9-5220.0*xi+t43+47304.0*t22-12447.0*eta+t36;
  t124 = 4968.0*t1;
  t129 = 1296.0*t4;
  t134 = -38880.0*t16+101088.0*t10+77760.0*t12+t129+t30-t61-t102-1016.0-52272.0*t9+3384.0*xi-t49-36720.0*t22+13392.0*eta-t6;
  t149 = 5832.0*t7;
  t154 = t96-44712.0*t10-55080.0*t12-t149+t15+t62+594.0+34668.0*t9-1197.0*xi+t54+14256.0*t22-8289.0*eta;
  t164 = 180.0*xi;
  t165 = 2268.0*t22;

  values[0] = t27;
  values[1] = 15552.0*t10+15552.0*t4+t30-t31-t32-t26+2088.0*xi-t34-10044.0*t22-t36;
  values[2] = -t38-19440.0*t4-21384.0*t7+t41+t15-1071.0*xi+t43+2916.0*t22+t36;
  values[3] = t46+t47-t31+360.0*xi-t49-t50-t6;
  values[4] = t55;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = t67;
  values[8] = -69984.0*t10-34992.0*t4-104976.0*t7+t71+t72+t62-6156.0*xi+25704.0*t1+38556.0*t22+t76;
  values[9] = t78+42768.0*t4+69984.0*t7-23328.0*t2-t61+2664.0*xi-t83-9720.0*t22-t76;
  values[10] = -22032.0*t4-t87+t71-792.0*xi+t89+t90+t36;
  values[11] = t92-t60+108.0*xi-t94;
  values[12] = 0.0;
  values[13] = t107;
  values[14] = 124416.0*t10+23328.0*t4+128304.0*t7-t31-139968.0*t14-t102+6984.0*xi-t83-57672.0*t22-t76;
  values[15] = -11664.0*t10-27216.0*t4-81648.0*t7+t41+t101-2214.0*xi+17496.0*t1+11664.0*t22+t76;
  values[16] = 12960.0*t4+t87-t31+504.0*xi-t124-t90-t36;
  values[17] = t55;
  values[18] = t134;
  values[19] = -108864.0*t10-t92-66096.0*t7+t72+t102-4032.0*xi+t89+42120.0*t22+t36;
  values[20] = t78+t92+38880.0*t7-t61+720.0*xi-t124-5832.0*t22-t36;
  values[21] = -t129-t47-72.0*xi+648.0*t1+t50+t6;
  values[22] = t154;
  values[23] = 46656.0*t10+t87-t32-t62+1296.0*xi-t94-15228.0*t22;
  values[24] = -t38-t149+t15-99.0*xi+t54+972.0*t22;
  values[25] = -t56+t78+20736.0*t12-t26-972.0/5.0-12312.0*t9+t164-t165+2808.0*eta;
  values[26] = -t78+t26-t164+t165;
  values[27] = t17-3240.0*t12+137.0/5.0+1836.0*t9-405.0*eta;
}
  
// ***********************************************************************

TBaseFunct2D *BF_C_T_P6_2D_Obj = new TBaseFunct2D
        (28, BF_C_T_P6_2D, BFUnitTriangle, 
         C_T_P6_2D_Funct, C_T_P6_2D_DeriveXi,
         C_T_P6_2D_DeriveEta, C_T_P6_2D_DeriveXiXi,
         C_T_P6_2D_DeriveXiEta, C_T_P6_2D_DeriveEtaEta, 6, 6,
         0, NULL);
