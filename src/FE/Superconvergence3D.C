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
   
//***************************************************** 
//  Superconvergence subroutine for Laplace, 
//     Stokes and Navier-Strokes 3D problems
//*****************************************************
//
#include <Superconvergence3D.h>

void Transform_Q2Q3_3D(double *fine_values, double *coarse_values)
{
  double s1, s2, s3;

  coarse_values[0] = fine_values[0];
  
  coarse_values[1] = 4.0/81.0*fine_values[30]-2.0/81.0*fine_values[27]+2.0/
9.0*fine_values[2]+68.0/81.0*fine_values[1]-7.0/81.0*fine_values[0];
      coarse_values[2] = 68.0/81.0*fine_values[30]-7.0/81.0*fine_values[27]+2.0
/9.0*fine_values[2]+4.0/81.0*fine_values[1]-2.0/81.0*fine_values[0];
      coarse_values[3] = fine_values[27];
      coarse_values[4] = 2.0/9.0*fine_values[6]+68.0/81.0*fine_values[3]-7.0/
81.0*fine_values[0]+4.0/81.0*fine_values[64]-2.0/81.0*fine_values[63];
      coarse_values[5] = 272.0/6561.0*fine_values[31]+8.0/729.0*fine_values[50]
+16.0/6561.0*fine_values[49]-8.0/6561.0*fine_values[48]-4.0/729.0*fine_values
[47]-8.0/6561.0*fine_values[46]+4.0/6561.0*fine_values[45]+8.0/729.0*
fine_values[32]-28.0/6561.0*fine_values[30]-4.0/729.0*fine_values[29]-136.0/
6561.0*fine_values[28]+14.0/6561.0*fine_values[27]+4.0/81.0*fine_values[8]+
136.0/729.0*fine_values[7]-14.0/729.0*fine_values[6]+136.0/729.0*fine_values[5]
+4624.0/6561.0*fine_values[4]-476.0/6561.0*fine_values[3]-14.0/729.0*
fine_values[2]-476.0/6561.0*fine_values[1]+272.0/6561.0*fine_values[66]-136.0/
6561.0*fine_values[65]-28.0/6561.0*fine_values[64]+14.0/6561.0*fine_values[63]+
49.0/6561.0*fine_values[0];
      coarse_values[6] = 4624.0/6561.0*fine_values[31]+8.0/729.0*fine_values
[50]+272.0/6561.0*fine_values[49]-28.0/6561.0*fine_values[48]-4.0/729.0*
fine_values[47]-136.0/6561.0*fine_values[46]+14.0/6561.0*fine_values[45]+136.0/
729.0*fine_values[32]-476.0/6561.0*fine_values[30]-14.0/729.0*fine_values[29]
-476.0/6561.0*fine_values[28]+49.0/6561.0*fine_values[27]+4.0/81.0*fine_values
[8]+8.0/729.0*fine_values[7]-4.0/729.0*fine_values[6]+136.0/729.0*fine_values
[5]+272.0/6561.0*fine_values[4]-136.0/6561.0*fine_values[3]-14.0/729.0*
fine_values[2]-28.0/6561.0*fine_values[1]+16.0/6561.0*fine_values[66]-8.0/
6561.0*fine_values[65]-8.0/6561.0*fine_values[64]+4.0/6561.0*fine_values[63]+
14.0/6561.0*fine_values[0];
      coarse_values[7] = 4.0/81.0*fine_values[48]-2.0/81.0*fine_values[45]+2.0/
9.0*fine_values[29]+68.0/81.0*fine_values[28]-7.0/81.0*fine_values[27];
      coarse_values[8] = 2.0/9.0*fine_values[6]+4.0/81.0*fine_values[3]-2.0/
81.0*fine_values[0]+68.0/81.0*fine_values[64]-7.0/81.0*fine_values[63];
      coarse_values[9] = 16.0/6561.0*fine_values[31]+136.0/729.0*fine_values
[50]+272.0/6561.0*fine_values[49]-136.0/6561.0*fine_values[48]-14.0/729.0*
fine_values[47]-28.0/6561.0*fine_values[46]+14.0/6561.0*fine_values[45]+8.0/
729.0*fine_values[32]-8.0/6561.0*fine_values[30]-4.0/729.0*fine_values[29]-8.0/
6561.0*fine_values[28]+4.0/6561.0*fine_values[27]+4.0/81.0*fine_values[8]+136.0
/729.0*fine_values[7]-14.0/729.0*fine_values[6]+8.0/729.0*fine_values[5]+272.0/
6561.0*fine_values[4]-28.0/6561.0*fine_values[3]-4.0/729.0*fine_values[2]-136.0
/6561.0*fine_values[1]+4624.0/6561.0*fine_values[66]-476.0/6561.0*fine_values
[65]-476.0/6561.0*fine_values[64]+49.0/6561.0*fine_values[63]+14.0/6561.0*
fine_values[0];
      coarse_values[10] = 272.0/6561.0*fine_values[31]+136.0/729.0*fine_values
[50]+4624.0/6561.0*fine_values[49]-476.0/6561.0*fine_values[48]-14.0/729.0*
fine_values[47]-476.0/6561.0*fine_values[46]+49.0/6561.0*fine_values[45]+136.0/
729.0*fine_values[32]-136.0/6561.0*fine_values[30]-14.0/729.0*fine_values[29]
-28.0/6561.0*fine_values[28]+14.0/6561.0*fine_values[27]+4.0/81.0*fine_values
[8]+8.0/729.0*fine_values[7]-4.0/729.0*fine_values[6]+8.0/729.0*fine_values[5]+
16.0/6561.0*fine_values[4]-8.0/6561.0*fine_values[3]-4.0/729.0*fine_values[2]
-8.0/6561.0*fine_values[1]+272.0/6561.0*fine_values[66]-28.0/6561.0*fine_values
[65]-136.0/6561.0*fine_values[64]+14.0/6561.0*fine_values[63]+4.0/6561.0*
fine_values[0];
      coarse_values[11] = 68.0/81.0*fine_values[48]-7.0/81.0*fine_values[45]+
2.0/9.0*fine_values[29]+4.0/81.0*fine_values[28]-2.0/81.0*fine_values[27];
      coarse_values[12] = fine_values[63];
      coarse_values[13] = 2.0/9.0*fine_values[47]+4.0/81.0*fine_values[46]-2.0/
81.0*fine_values[45]+68.0/81.0*fine_values[65]-7.0/81.0*fine_values[63];
      coarse_values[14] = 2.0/9.0*fine_values[47]+68.0/81.0*fine_values[46]-7.0
/81.0*fine_values[45]+4.0/81.0*fine_values[65]-2.0/81.0*fine_values[63];
      coarse_values[15] = fine_values[45];
      coarse_values[16] = 2.0/9.0*fine_values[18]-2.0/81.0*fine_values[75]+68.0
/81.0*fine_values[9]-7.0/81.0*fine_values[0]+4.0/81.0*fine_values[84];
      coarse_values[17] = 136.0/729.0*fine_values[19]+136.0/729.0*fine_values
[11]-4.0/729.0*fine_values[81]+8.0/729.0*fine_values[42]-4.0/729.0*fine_values
[39]+272.0/6561.0*fine_values[36]-136.0/6561.0*fine_values[33]-28.0/6561.0*
fine_values[30]+14.0/6561.0*fine_values[27]+4.0/81.0*fine_values[20]-14.0/729.0
*fine_values[18]+4624.0/6561.0*fine_values[10]-476.0/6561.0*fine_values[9]-14.0
/729.0*fine_values[2]-476.0/6561.0*fine_values[1]+16.0/6561.0*fine_values[100]
-8.0/6561.0*fine_values[99]-8.0/6561.0*fine_values[94]+4.0/6561.0*fine_values
[93]+8.0/729.0*fine_values[90]+272.0/6561.0*fine_values[87]-28.0/6561.0*
fine_values[84]-136.0/6561.0*fine_values[78]+14.0/6561.0*fine_values[75]+49.0/
6561.0*fine_values[0];
      coarse_values[18] = 8.0/729.0*fine_values[19]+136.0/729.0*fine_values[11]
-4.0/729.0*fine_values[81]+136.0/729.0*fine_values[42]-14.0/729.0*fine_values
[39]+4624.0/6561.0*fine_values[36]-476.0/6561.0*fine_values[33]-476.0/6561.0*
fine_values[30]+49.0/6561.0*fine_values[27]+4.0/81.0*fine_values[20]-4.0/729.0*
fine_values[18]+272.0/6561.0*fine_values[10]-136.0/6561.0*fine_values[9]-14.0/
729.0*fine_values[2]-28.0/6561.0*fine_values[1]+272.0/6561.0*fine_values[100]
-28.0/6561.0*fine_values[99]-136.0/6561.0*fine_values[94]+14.0/6561.0*
fine_values[93]+8.0/729.0*fine_values[90]+16.0/6561.0*fine_values[87]-8.0/
6561.0*fine_values[84]-8.0/6561.0*fine_values[78]+4.0/6561.0*fine_values[75]+
14.0/6561.0*fine_values[0];
      coarse_values[19] = 2.0/9.0*fine_values[39]+68.0/81.0*fine_values[33]-7.0
/81.0*fine_values[27]+4.0/81.0*fine_values[99]-2.0/81.0*fine_values[93];
      coarse_values[20] = 136.0/729.0*fine_values[21]-136.0/6561.0*fine_values
[76]+4.0/81.0*fine_values[24]-14.0/729.0*fine_values[18]+136.0/729.0*
fine_values[15]+4624.0/6561.0*fine_values[12]-476.0/6561.0*fine_values[9]-14.0/
729.0*fine_values[6]-476.0/6561.0*fine_values[3]+16.0/6561.0*fine_values[123]
-8.0/6561.0*fine_values[121]-8.0/6561.0*fine_values[119]+4.0/6561.0*fine_values
[117]+8.0/729.0*fine_values[86]+272.0/6561.0*fine_values[85]-28.0/6561.0*
fine_values[84]-4.0/729.0*fine_values[77]+14.0/6561.0*fine_values[75]+8.0/729.0
*fine_values[72]-4.0/729.0*fine_values[71]+272.0/6561.0*fine_values[68]-136.0/
6561.0*fine_values[67]-28.0/6561.0*fine_values[64]+14.0/6561.0*fine_values[63]+
49.0/6561.0*fine_values[0];
      s2 = -8.0/6561.0*fine_values[41]-1904.0/531441.0*fine_values[31]-952.0/
59049.0*fine_values[21]-952.0/59049.0*fine_values[19]+272.0/6561.0*fine_values
[17]-952.0/59049.0*fine_values[11]-16.0/59049.0*fine_values[60]+544.0/59049.0*
fine_values[56]-272.0/59049.0*fine_values[53]+544.0/59049.0*fine_values[89]
-272.0/59049.0*fine_values[82]+28.0/59049.0*fine_values[81]-9248.0/531441.0*
fine_values[79]+952.0/531441.0*fine_values[76]-544.0/531441.0*fine_values[120]
-16.0/59049.0*fine_values[115]+16.0/531441.0*fine_values[107]+272.0/531441.0*
fine_values[51]-56.0/59049.0*fine_values[50]-112.0/531441.0*fine_values[49]+
56.0/531441.0*fine_values[48]+28.0/59049.0*fine_values[47]+56.0/531441.0*
fine_values[46]-28.0/531441.0*fine_values[45]+16.0/6561.0*fine_values[44]+544.0
/59049.0*fine_values[43]-56.0/59049.0*fine_values[42]-272.0/59049.0*fine_values
[40]+28.0/59049.0*fine_values[39]+544.0/59049.0*fine_values[38]+18496.0/
531441.0*fine_values[37];
      s1 = -1904.0/531441.0*fine_values[36]-272.0/59049.0*fine_values[35]
-9248.0/531441.0*fine_values[34]+952.0/531441.0*fine_values[33]-56.0/59049.0*
fine_values[32]+196.0/531441.0*fine_values[30]+28.0/59049.0*fine_values[29]+
952.0/531441.0*fine_values[28]-98.0/531441.0*fine_values[27]+8.0/729.0*
fine_values[26]+272.0/6561.0*fine_values[25]-28.0/6561.0*fine_values[24]+272.0/
6561.0*fine_values[23]+s2+9248.0/59049.0*fine_values[22]-28.0/6561.0*
fine_values[20]+98.0/59049.0*fine_values[18]+9248.0/59049.0*fine_values[16]
-952.0/59049.0*fine_values[15]+9248.0/59049.0*fine_values[14]+314432.0/531441.0
*fine_values[13]-32368.0/531441.0*fine_values[12]-32368.0/531441.0*fine_values
[10]+3332.0/531441.0*fine_values[9]-28.0/6561.0*fine_values[8]-952.0/59049.0*
fine_values[7]+98.0/59049.0*fine_values[6]-952.0/59049.0*fine_values[5]-32368.0
/531441.0*fine_values[4]+3332.0/531441.0*fine_values[3]+98.0/59049.0*
fine_values[2]+3332.0/531441.0*fine_values[1];
      s3 = 1088.0/531441.0*fine_values[124]-112.0/531441.0*fine_values[123]
-544.0/531441.0*fine_values[122]+56.0/531441.0*fine_values[121]+56.0/531441.0*
fine_values[119]+272.0/531441.0*fine_values[118]-28.0/531441.0*fine_values[117]
+32.0/59049.0*fine_values[116]+64.0/531441.0*fine_values[114]-32.0/531441.0*
fine_values[113]-32.0/531441.0*fine_values[112]+16.0/531441.0*fine_values[111]
-16.0/59049.0*fine_values[110]+8.0/59049.0*fine_values[109]-32.0/531441.0*
fine_values[108]+16.0/531441.0*fine_values[106];
      s2 = s3-8.0/531441.0*fine_values[105]+32.0/59049.0*fine_values[104]-16.0/
59049.0*fine_values[103]+1088.0/531441.0*fine_values[102]-544.0/531441.0*
fine_values[101]-112.0/531441.0*fine_values[100]+56.0/531441.0*fine_values[99]
-16.0/59049.0*fine_values[98]+8.0/59049.0*fine_values[97]-544.0/531441.0*
fine_values[96]+272.0/531441.0*fine_values[95]+56.0/531441.0*fine_values[94]
-28.0/531441.0*fine_values[93]+16.0/6561.0*fine_values[92]+544.0/59049.0*
fine_values[91]-56.0/59049.0*fine_values[90];
      s3 = s2+18496.0/531441.0*fine_values[88]-1904.0/531441.0*fine_values[87]
-56.0/59049.0*fine_values[86]-1904.0/531441.0*fine_values[85]+196.0/531441.0*
fine_values[84]-8.0/6561.0*fine_values[83]-272.0/59049.0*fine_values[80]+952.0/
531441.0*fine_values[78]+28.0/59049.0*fine_values[77]-98.0/531441.0*fine_values
[75]+544.0/59049.0*fine_values[74]-272.0/59049.0*fine_values[73]-56.0/59049.0*
fine_values[72]+28.0/59049.0*fine_values[71]+18496.0/531441.0*fine_values[70];
      coarse_values[21] = s3-9248.0/531441.0*fine_values[69]-1904.0/531441.0*
fine_values[68]+952.0/531441.0*fine_values[67]-1904.0/531441.0*fine_values[66]+
952.0/531441.0*fine_values[65]+196.0/531441.0*fine_values[64]-98.0/531441.0*
fine_values[63]+16.0/6561.0*fine_values[62]+32.0/59049.0*fine_values[61]-8.0/
6561.0*fine_values[59]-16.0/59049.0*fine_values[58]+8.0/59049.0*fine_values[57]
+1088.0/531441.0*fine_values[55]-544.0/531441.0*fine_values[54]-544.0/531441.0*
fine_values[52]+s1-343.0/531441.0*fine_values[0];
      s2 = -28.0/6561.0*fine_values[41]-32368.0/531441.0*fine_values[31]-272.0/
59049.0*fine_values[21]-56.0/59049.0*fine_values[19]+272.0/6561.0*fine_values
[17]-952.0/59049.0*fine_values[11]-56.0/59049.0*fine_values[60]+544.0/59049.0*
fine_values[56]-272.0/59049.0*fine_values[53]+32.0/59049.0*fine_values[89]
-272.0/59049.0*fine_values[82]+28.0/59049.0*fine_values[81]-544.0/531441.0*
fine_values[79]+272.0/531441.0*fine_values[76]-32.0/531441.0*fine_values[120]
-16.0/59049.0*fine_values[115]+272.0/531441.0*fine_values[107]+952.0/531441.0*
fine_values[51]-56.0/59049.0*fine_values[50]-1904.0/531441.0*fine_values[49]+
196.0/531441.0*fine_values[48]+28.0/59049.0*fine_values[47]+952.0/531441.0*
fine_values[46]-98.0/531441.0*fine_values[45]+272.0/6561.0*fine_values[44]+
9248.0/59049.0*fine_values[43]-952.0/59049.0*fine_values[42]-952.0/59049.0*
fine_values[40]+98.0/59049.0*fine_values[39]+9248.0/59049.0*fine_values[38]+
314432.0/531441.0*fine_values[37];
      s1 = -32368.0/531441.0*fine_values[36]-952.0/59049.0*fine_values[35]
-32368.0/531441.0*fine_values[34]+3332.0/531441.0*fine_values[33]-952.0/59049.0
*fine_values[32]+3332.0/531441.0*fine_values[30]+98.0/59049.0*fine_values[29]+
3332.0/531441.0*fine_values[28]-343.0/531441.0*fine_values[27]+8.0/729.0*
fine_values[26]+16.0/6561.0*fine_values[25]-8.0/6561.0*fine_values[24]+272.0/
6561.0*fine_values[23]+s2+544.0/59049.0*fine_values[22]-28.0/6561.0*fine_values
[20]+28.0/59049.0*fine_values[18]+544.0/59049.0*fine_values[16]-272.0/59049.0*
fine_values[15]+9248.0/59049.0*fine_values[14]+18496.0/531441.0*fine_values[13]
-9248.0/531441.0*fine_values[12]-1904.0/531441.0*fine_values[10]+952.0/531441.0
*fine_values[9]-28.0/6561.0*fine_values[8]-56.0/59049.0*fine_values[7]+28.0/
59049.0*fine_values[6]-952.0/59049.0*fine_values[5]-1904.0/531441.0*fine_values
[4]+952.0/531441.0*fine_values[3]+98.0/59049.0*fine_values[2]+196.0/531441.0*
fine_values[1];
      s3 = 64.0/531441.0*fine_values[124]-32.0/531441.0*fine_values[123]-32.0/
531441.0*fine_values[122]+16.0/531441.0*fine_values[121]+16.0/531441.0*
fine_values[119]+16.0/531441.0*fine_values[118]-8.0/531441.0*fine_values[117]+
32.0/59049.0*fine_values[116]+1088.0/531441.0*fine_values[114]-544.0/531441.0*
fine_values[113]-112.0/531441.0*fine_values[112]+56.0/531441.0*fine_values[111]
-16.0/59049.0*fine_values[110]+8.0/59049.0*fine_values[109]-544.0/531441.0*
fine_values[108]+56.0/531441.0*fine_values[106];
      s2 = s3-28.0/531441.0*fine_values[105]+544.0/59049.0*fine_values[104]
-56.0/59049.0*fine_values[103]+18496.0/531441.0*fine_values[102]-1904.0/
531441.0*fine_values[101]-1904.0/531441.0*fine_values[100]+196.0/531441.0*
fine_values[99]-272.0/59049.0*fine_values[98]+28.0/59049.0*fine_values[97]
-9248.0/531441.0*fine_values[96]+952.0/531441.0*fine_values[95]+952.0/531441.0*
fine_values[94]-98.0/531441.0*fine_values[93]+16.0/6561.0*fine_values[92]+544.0
/59049.0*fine_values[91]-56.0/59049.0*fine_values[90];
      coarse_values[22] = s2+1088.0/531441.0*fine_values[88]-112.0/531441.0*
fine_values[87]-16.0/59049.0*fine_values[86]-544.0/531441.0*fine_values[85]+
56.0/531441.0*fine_values[84]-8.0/6561.0*fine_values[83]-16.0/59049.0*
fine_values[80]+56.0/531441.0*fine_values[78]+8.0/59049.0*fine_values[77]-28.0/
531441.0*fine_values[75]+32.0/59049.0*fine_values[74]-16.0/59049.0*fine_values
[73]-16.0/59049.0*fine_values[72]+8.0/59049.0*fine_values[71]+1088.0/531441.0*
fine_values[70]-544.0/531441.0*fine_values[69]-544.0/531441.0*fine_values[68]+
272.0/531441.0*fine_values[67]-112.0/531441.0*fine_values[66]+56.0/531441.0*
fine_values[65]+56.0/531441.0*fine_values[64]-28.0/531441.0*fine_values[63]+
16.0/6561.0*fine_values[62]+544.0/59049.0*fine_values[61]-8.0/6561.0*
fine_values[59]-272.0/59049.0*fine_values[58]+28.0/59049.0*fine_values[57]+
18496.0/531441.0*fine_values[55]-1904.0/531441.0*fine_values[54]-9248.0/
531441.0*fine_values[52]+s1-98.0/531441.0*fine_values[0];
      coarse_values[23] = 4.0/81.0*fine_values[41]+8.0/729.0*fine_values[60]
-136.0/6561.0*fine_values[51]-28.0/6561.0*fine_values[48]+14.0/6561.0*
fine_values[45]+136.0/729.0*fine_values[40]-14.0/729.0*fine_values[39]+136.0/
729.0*fine_values[35]+4624.0/6561.0*fine_values[34]-476.0/6561.0*fine_values
[33]-14.0/729.0*fine_values[29]-476.0/6561.0*fine_values[28]+49.0/6561.0*
fine_values[27]+16.0/6561.0*fine_values[112]-8.0/6561.0*fine_values[111]-8.0/
6561.0*fine_values[106]+4.0/6561.0*fine_values[105]+8.0/729.0*fine_values[103]+
272.0/6561.0*fine_values[101]-28.0/6561.0*fine_values[99]-4.0/729.0*fine_values
[97]-136.0/6561.0*fine_values[95]+14.0/6561.0*fine_values[93]-4.0/729.0*
fine_values[57]+272.0/6561.0*fine_values[54];
      coarse_values[24] = 8.0/729.0*fine_values[21]-8.0/6561.0*fine_values[76]+
4.0/81.0*fine_values[24]-4.0/729.0*fine_values[18]+136.0/729.0*fine_values[15]+
272.0/6561.0*fine_values[12]-136.0/6561.0*fine_values[9]-14.0/729.0*fine_values
[6]-28.0/6561.0*fine_values[3]+272.0/6561.0*fine_values[123]-28.0/6561.0*
fine_values[121]-136.0/6561.0*fine_values[119]+14.0/6561.0*fine_values[117]+8.0
/729.0*fine_values[86]+16.0/6561.0*fine_values[85]-8.0/6561.0*fine_values[84]
-4.0/729.0*fine_values[77]+4.0/6561.0*fine_values[75]+136.0/729.0*fine_values
[72]-14.0/729.0*fine_values[71]+4624.0/6561.0*fine_values[68]-476.0/6561.0*
fine_values[67]-476.0/6561.0*fine_values[64]+49.0/6561.0*fine_values[63]+14.0/
6561.0*fine_values[0];
      s2 = -8.0/6561.0*fine_values[41]-112.0/531441.0*fine_values[31]-56.0/
59049.0*fine_values[21]-272.0/59049.0*fine_values[19]+272.0/6561.0*fine_values
[17]-272.0/59049.0*fine_values[11]-272.0/59049.0*fine_values[60]+9248.0/59049.0
*fine_values[56]-952.0/59049.0*fine_values[53]+544.0/59049.0*fine_values[89]
-16.0/59049.0*fine_values[82]+8.0/59049.0*fine_values[81]-544.0/531441.0*
fine_values[79]+56.0/531441.0*fine_values[76]-9248.0/531441.0*fine_values[120]
-56.0/59049.0*fine_values[115]+56.0/531441.0*fine_values[107]+952.0/531441.0*
fine_values[51]-952.0/59049.0*fine_values[50]-1904.0/531441.0*fine_values[49]+
952.0/531441.0*fine_values[48]+98.0/59049.0*fine_values[47]+196.0/531441.0*
fine_values[46]-98.0/531441.0*fine_values[45]+16.0/6561.0*fine_values[44]+32.0/
59049.0*fine_values[43]-16.0/59049.0*fine_values[42]-16.0/59049.0*fine_values
[40]+8.0/59049.0*fine_values[39]+544.0/59049.0*fine_values[38]+1088.0/531441.0*
fine_values[37];
      s1 = -544.0/531441.0*fine_values[36]-272.0/59049.0*fine_values[35]-544.0/
531441.0*fine_values[34]+272.0/531441.0*fine_values[33]-56.0/59049.0*
fine_values[32]+56.0/531441.0*fine_values[30]+28.0/59049.0*fine_values[29]+56.0
/531441.0*fine_values[28]-28.0/531441.0*fine_values[27]+8.0/729.0*fine_values
[26]+272.0/6561.0*fine_values[25]-28.0/6561.0*fine_values[24]+16.0/6561.0*
fine_values[23]+s2+544.0/59049.0*fine_values[22]-8.0/6561.0*fine_values[20]+
28.0/59049.0*fine_values[18]+9248.0/59049.0*fine_values[16]-952.0/59049.0*
fine_values[15]+544.0/59049.0*fine_values[14]+18496.0/531441.0*fine_values[13]
-1904.0/531441.0*fine_values[12]-9248.0/531441.0*fine_values[10]+952.0/531441.0
*fine_values[9]-28.0/6561.0*fine_values[8]-952.0/59049.0*fine_values[7]+98.0/
59049.0*fine_values[6]-56.0/59049.0*fine_values[5]-1904.0/531441.0*fine_values
[4]+196.0/531441.0*fine_values[3]+28.0/59049.0*fine_values[2]+952.0/531441.0*
fine_values[1];
      s3 = 18496.0/531441.0*fine_values[124]-1904.0/531441.0*fine_values[123]
-1904.0/531441.0*fine_values[122]+196.0/531441.0*fine_values[121]+952.0/
531441.0*fine_values[119]+952.0/531441.0*fine_values[118]-98.0/531441.0*
fine_values[117]+544.0/59049.0*fine_values[116]+1088.0/531441.0*fine_values
[114]-112.0/531441.0*fine_values[113]-544.0/531441.0*fine_values[112]+56.0/
531441.0*fine_values[111]-272.0/59049.0*fine_values[110]+28.0/59049.0*
fine_values[109]-544.0/531441.0*fine_values[108]+272.0/531441.0*fine_values
[106];
      s2 = s3-28.0/531441.0*fine_values[105]+32.0/59049.0*fine_values[104]-16.0
/59049.0*fine_values[103]+64.0/531441.0*fine_values[102]-32.0/531441.0*
fine_values[101]-32.0/531441.0*fine_values[100]+16.0/531441.0*fine_values[99]
-16.0/59049.0*fine_values[98]+8.0/59049.0*fine_values[97]-32.0/531441.0*
fine_values[96]+16.0/531441.0*fine_values[95]+16.0/531441.0*fine_values[94]-8.0
/531441.0*fine_values[93]+16.0/6561.0*fine_values[92]+32.0/59049.0*fine_values
[91]-16.0/59049.0*fine_values[90];
      s3 = s2+1088.0/531441.0*fine_values[88]-544.0/531441.0*fine_values[87]
-56.0/59049.0*fine_values[86]-112.0/531441.0*fine_values[85]+56.0/531441.0*
fine_values[84]-8.0/6561.0*fine_values[83]-272.0/59049.0*fine_values[80]+272.0/
531441.0*fine_values[78]+28.0/59049.0*fine_values[77]-28.0/531441.0*fine_values
[75]+9248.0/59049.0*fine_values[74]-952.0/59049.0*fine_values[73]-952.0/59049.0
*fine_values[72]+98.0/59049.0*fine_values[71]+314432.0/531441.0*fine_values[70]
;
      coarse_values[25] = s3-32368.0/531441.0*fine_values[69]-32368.0/531441.0*
fine_values[68]+3332.0/531441.0*fine_values[67]-32368.0/531441.0*fine_values
[66]+3332.0/531441.0*fine_values[65]+3332.0/531441.0*fine_values[64]-343.0/
531441.0*fine_values[63]+272.0/6561.0*fine_values[62]+544.0/59049.0*fine_values
[61]-28.0/6561.0*fine_values[59]-56.0/59049.0*fine_values[58]+28.0/59049.0*
fine_values[57]+18496.0/531441.0*fine_values[55]-9248.0/531441.0*fine_values
[54]-1904.0/531441.0*fine_values[52]+s1-98.0/531441.0*fine_values[0];
      s2 = -28.0/6561.0*fine_values[41]-1904.0/531441.0*fine_values[31]-16.0/
59049.0*fine_values[21]-16.0/59049.0*fine_values[19]+272.0/6561.0*fine_values
[17]-272.0/59049.0*fine_values[11]-952.0/59049.0*fine_values[60]+9248.0/59049.0
*fine_values[56]-952.0/59049.0*fine_values[53]+32.0/59049.0*fine_values[89]
-16.0/59049.0*fine_values[82]+8.0/59049.0*fine_values[81]-32.0/531441.0*
fine_values[79]+16.0/531441.0*fine_values[76]-544.0/531441.0*fine_values[120]
-56.0/59049.0*fine_values[115]+952.0/531441.0*fine_values[107]+3332.0/531441.0*
fine_values[51]-952.0/59049.0*fine_values[50]-32368.0/531441.0*fine_values[49]+
3332.0/531441.0*fine_values[48]+98.0/59049.0*fine_values[47]+3332.0/531441.0*
fine_values[46]-343.0/531441.0*fine_values[45]+272.0/6561.0*fine_values[44]+
544.0/59049.0*fine_values[43]-272.0/59049.0*fine_values[42]-56.0/59049.0*
fine_values[40]+28.0/59049.0*fine_values[39]+9248.0/59049.0*fine_values[38]+
18496.0/531441.0*fine_values[37];
      s1 = -9248.0/531441.0*fine_values[36]-952.0/59049.0*fine_values[35]
-1904.0/531441.0*fine_values[34]+952.0/531441.0*fine_values[33]-952.0/59049.0*
fine_values[32]+952.0/531441.0*fine_values[30]+98.0/59049.0*fine_values[29]+
196.0/531441.0*fine_values[28]-98.0/531441.0*fine_values[27]+8.0/729.0*
fine_values[26]+16.0/6561.0*fine_values[25]-8.0/6561.0*fine_values[24]+16.0/
6561.0*fine_values[23]+s2+32.0/59049.0*fine_values[22]-8.0/6561.0*fine_values
[20]+8.0/59049.0*fine_values[18]+544.0/59049.0*fine_values[16]-272.0/59049.0*
fine_values[15]+544.0/59049.0*fine_values[14]+1088.0/531441.0*fine_values[13]
-544.0/531441.0*fine_values[12]-544.0/531441.0*fine_values[10]+272.0/531441.0*
fine_values[9]-28.0/6561.0*fine_values[8]-56.0/59049.0*fine_values[7]+28.0/
59049.0*fine_values[6]-56.0/59049.0*fine_values[5]-112.0/531441.0*fine_values
[4]+56.0/531441.0*fine_values[3]+28.0/59049.0*fine_values[2]+56.0/531441.0*
fine_values[1];
      s3 = 1088.0/531441.0*fine_values[124]-544.0/531441.0*fine_values[123]
-112.0/531441.0*fine_values[122]+56.0/531441.0*fine_values[121]+272.0/531441.0*
fine_values[119]+56.0/531441.0*fine_values[118]-28.0/531441.0*fine_values[117]+
544.0/59049.0*fine_values[116]+18496.0/531441.0*fine_values[114]-1904.0/
531441.0*fine_values[113]-1904.0/531441.0*fine_values[112]+196.0/531441.0*
fine_values[111]-272.0/59049.0*fine_values[110]+28.0/59049.0*fine_values[109]
-9248.0/531441.0*fine_values[108]+952.0/531441.0*fine_values[106];
      s2 = s3-98.0/531441.0*fine_values[105]+544.0/59049.0*fine_values[104]
-56.0/59049.0*fine_values[103]+1088.0/531441.0*fine_values[102]-112.0/531441.0*
fine_values[101]-544.0/531441.0*fine_values[100]+56.0/531441.0*fine_values[99]
-272.0/59049.0*fine_values[98]+28.0/59049.0*fine_values[97]-544.0/531441.0*
fine_values[96]+56.0/531441.0*fine_values[95]+272.0/531441.0*fine_values[94]
-28.0/531441.0*fine_values[93]+16.0/6561.0*fine_values[92]+32.0/59049.0*
fine_values[91]-16.0/59049.0*fine_values[90];
      s3 = s2+64.0/531441.0*fine_values[88]-32.0/531441.0*fine_values[87]-16.0/
59049.0*fine_values[86]-32.0/531441.0*fine_values[85]+16.0/531441.0*fine_values
[84]-8.0/6561.0*fine_values[83]-16.0/59049.0*fine_values[80]+16.0/531441.0*
fine_values[78]+8.0/59049.0*fine_values[77]-8.0/531441.0*fine_values[75]+544.0/
59049.0*fine_values[74]-56.0/59049.0*fine_values[73]-272.0/59049.0*fine_values
[72]+28.0/59049.0*fine_values[71]+18496.0/531441.0*fine_values[70];
      coarse_values[26] = s3-1904.0/531441.0*fine_values[69]-9248.0/531441.0*
fine_values[68]+952.0/531441.0*fine_values[67]-1904.0/531441.0*fine_values[66]+
196.0/531441.0*fine_values[65]+952.0/531441.0*fine_values[64]-98.0/531441.0*
fine_values[63]+272.0/6561.0*fine_values[62]+9248.0/59049.0*fine_values[61]
-28.0/6561.0*fine_values[59]-952.0/59049.0*fine_values[58]+98.0/59049.0*
fine_values[57]+314432.0/531441.0*fine_values[55]-32368.0/531441.0*fine_values
[54]-32368.0/531441.0*fine_values[52]+s1-28.0/531441.0*fine_values[0];
      coarse_values[27] = 4.0/81.0*fine_values[41]+136.0/729.0*fine_values[60]
-476.0/6561.0*fine_values[51]-476.0/6561.0*fine_values[48]+49.0/6561.0*
fine_values[45]+8.0/729.0*fine_values[40]-4.0/729.0*fine_values[39]+136.0/729.0
*fine_values[35]+272.0/6561.0*fine_values[34]-136.0/6561.0*fine_values[33]-14.0
/729.0*fine_values[29]-28.0/6561.0*fine_values[28]+14.0/6561.0*fine_values[27]+
272.0/6561.0*fine_values[112]-28.0/6561.0*fine_values[111]-136.0/6561.0*
fine_values[106]+14.0/6561.0*fine_values[105]+8.0/729.0*fine_values[103]+16.0/
6561.0*fine_values[101]-8.0/6561.0*fine_values[99]-4.0/729.0*fine_values[97]
-8.0/6561.0*fine_values[95]+4.0/6561.0*fine_values[93]-14.0/729.0*fine_values
[57]+4624.0/6561.0*fine_values[54];
      coarse_values[28] = 4.0/81.0*fine_values[121]-2.0/81.0*fine_values[117]+
2.0/9.0*fine_values[71]+68.0/81.0*fine_values[67]-7.0/81.0*fine_values[63];
      coarse_values[29] = 136.0/729.0*fine_values[53]+8.0/729.0*fine_values
[115]-8.0/6561.0*fine_values[107]-136.0/6561.0*fine_values[51]-14.0/729.0*
fine_values[47]-28.0/6561.0*fine_values[46]+14.0/6561.0*fine_values[45]+272.0/
6561.0*fine_values[122]-28.0/6561.0*fine_values[121]-136.0/6561.0*fine_values
[118]+14.0/6561.0*fine_values[117]+16.0/6561.0*fine_values[113]-8.0/6561.0*
fine_values[111]-4.0/729.0*fine_values[109]+4.0/6561.0*fine_values[105]+136.0/
729.0*fine_values[73]-14.0/729.0*fine_values[71]+4624.0/6561.0*fine_values[69]
-476.0/6561.0*fine_values[67]-476.0/6561.0*fine_values[65]+49.0/6561.0*
fine_values[63]+4.0/81.0*fine_values[59]+8.0/729.0*fine_values[58]-4.0/729.0*
fine_values[57]+272.0/6561.0*fine_values[52];
      coarse_values[30] = 136.0/729.0*fine_values[53]+8.0/729.0*fine_values
[115]-136.0/6561.0*fine_values[107]-476.0/6561.0*fine_values[51]-14.0/729.0*
fine_values[47]-476.0/6561.0*fine_values[46]+49.0/6561.0*fine_values[45]+16.0/
6561.0*fine_values[122]-8.0/6561.0*fine_values[121]-8.0/6561.0*fine_values[118]
+4.0/6561.0*fine_values[117]+272.0/6561.0*fine_values[113]-28.0/6561.0*
fine_values[111]-4.0/729.0*fine_values[109]+14.0/6561.0*fine_values[105]+8.0/
729.0*fine_values[73]-4.0/729.0*fine_values[71]+272.0/6561.0*fine_values[69]
-136.0/6561.0*fine_values[67]-28.0/6561.0*fine_values[65]+14.0/6561.0*
fine_values[63]+4.0/81.0*fine_values[59]+136.0/729.0*fine_values[58]-14.0/729.0
*fine_values[57]+4624.0/6561.0*fine_values[52];
      coarse_values[31] = -7.0/81.0*fine_values[45]+4.0/81.0*fine_values[111]
-2.0/81.0*fine_values[105]+2.0/9.0*fine_values[57]+68.0/81.0*fine_values[51];
      coarse_values[32] = 2.0/9.0*fine_values[18]-7.0/81.0*fine_values[75]+4.0/
81.0*fine_values[9]-2.0/81.0*fine_values[0]+68.0/81.0*fine_values[84];
      coarse_values[33] = 136.0/729.0*fine_values[19]+8.0/729.0*fine_values[11]
-14.0/729.0*fine_values[81]+8.0/729.0*fine_values[42]-4.0/729.0*fine_values[39]
+16.0/6561.0*fine_values[36]-8.0/6561.0*fine_values[33]-8.0/6561.0*fine_values
[30]+4.0/6561.0*fine_values[27]+4.0/81.0*fine_values[20]-14.0/729.0*fine_values
[18]+272.0/6561.0*fine_values[10]-28.0/6561.0*fine_values[9]-4.0/729.0*
fine_values[2]-136.0/6561.0*fine_values[1]+272.0/6561.0*fine_values[100]-136.0/
6561.0*fine_values[99]-28.0/6561.0*fine_values[94]+14.0/6561.0*fine_values[93]+
136.0/729.0*fine_values[90]+4624.0/6561.0*fine_values[87]-476.0/6561.0*
fine_values[84]-476.0/6561.0*fine_values[78]+49.0/6561.0*fine_values[75]+14.0/
6561.0*fine_values[0];
      coarse_values[34] = 8.0/729.0*fine_values[19]+8.0/729.0*fine_values[11]
-14.0/729.0*fine_values[81]+136.0/729.0*fine_values[42]-14.0/729.0*fine_values
[39]+272.0/6561.0*fine_values[36]-28.0/6561.0*fine_values[33]-136.0/6561.0*
fine_values[30]+14.0/6561.0*fine_values[27]+4.0/81.0*fine_values[20]-4.0/729.0*
fine_values[18]+16.0/6561.0*fine_values[10]-8.0/6561.0*fine_values[9]-4.0/729.0
*fine_values[2]-8.0/6561.0*fine_values[1]+4624.0/6561.0*fine_values[100]-476.0/
6561.0*fine_values[99]-476.0/6561.0*fine_values[94]+49.0/6561.0*fine_values[93]
+136.0/729.0*fine_values[90]+272.0/6561.0*fine_values[87]-136.0/6561.0*
fine_values[84]-28.0/6561.0*fine_values[78]+14.0/6561.0*fine_values[75]+4.0/
6561.0*fine_values[0];
      coarse_values[35] = 2.0/9.0*fine_values[39]+4.0/81.0*fine_values[33]-2.0/
81.0*fine_values[27]+68.0/81.0*fine_values[99]-7.0/81.0*fine_values[93];
      coarse_values[36] = 136.0/729.0*fine_values[21]-476.0/6561.0*fine_values
[76]+4.0/81.0*fine_values[24]-14.0/729.0*fine_values[18]+8.0/729.0*fine_values
[15]+272.0/6561.0*fine_values[12]-28.0/6561.0*fine_values[9]-4.0/729.0*
fine_values[6]-136.0/6561.0*fine_values[3]+272.0/6561.0*fine_values[123]-136.0/
6561.0*fine_values[121]-28.0/6561.0*fine_values[119]+14.0/6561.0*fine_values
[117]+136.0/729.0*fine_values[86]+4624.0/6561.0*fine_values[85]-476.0/6561.0*
fine_values[84]-14.0/729.0*fine_values[77]+49.0/6561.0*fine_values[75]+8.0/
729.0*fine_values[72]-4.0/729.0*fine_values[71]+16.0/6561.0*fine_values[68]-8.0
/6561.0*fine_values[67]-8.0/6561.0*fine_values[64]+4.0/6561.0*fine_values[63]+
14.0/6561.0*fine_values[0];
      s2 = -8.0/6561.0*fine_values[41]-544.0/531441.0*fine_values[31]-952.0/
59049.0*fine_values[21]-952.0/59049.0*fine_values[19]+16.0/6561.0*fine_values
[17]-56.0/59049.0*fine_values[11]-16.0/59049.0*fine_values[60]+32.0/59049.0*
fine_values[56]-16.0/59049.0*fine_values[53]+9248.0/59049.0*fine_values[89]
-952.0/59049.0*fine_values[82]+98.0/59049.0*fine_values[81]-32368.0/531441.0*
fine_values[79]+3332.0/531441.0*fine_values[76]-1904.0/531441.0*fine_values
[120]-272.0/59049.0*fine_values[115]+56.0/531441.0*fine_values[107]+16.0/
531441.0*fine_values[51]-16.0/59049.0*fine_values[50]-32.0/531441.0*fine_values
[49]+16.0/531441.0*fine_values[48]+8.0/59049.0*fine_values[47]+16.0/531441.0*
fine_values[46]-8.0/531441.0*fine_values[45]+16.0/6561.0*fine_values[44]+544.0/
59049.0*fine_values[43]-56.0/59049.0*fine_values[42]-272.0/59049.0*fine_values
[40]+28.0/59049.0*fine_values[39]+32.0/59049.0*fine_values[38]+1088.0/531441.0*
fine_values[37];
      s1 = -112.0/531441.0*fine_values[36]-16.0/59049.0*fine_values[35]-544.0/
531441.0*fine_values[34]+56.0/531441.0*fine_values[33]-16.0/59049.0*fine_values
[32]+56.0/531441.0*fine_values[30]+8.0/59049.0*fine_values[29]+272.0/531441.0*
fine_values[28]-28.0/531441.0*fine_values[27]+8.0/729.0*fine_values[26]+272.0/
6561.0*fine_values[25]-28.0/6561.0*fine_values[24]+272.0/6561.0*fine_values[23]
+s2+9248.0/59049.0*fine_values[22]-28.0/6561.0*fine_values[20]+98.0/59049.0*
fine_values[18]+544.0/59049.0*fine_values[16]-56.0/59049.0*fine_values[15]+
544.0/59049.0*fine_values[14]+18496.0/531441.0*fine_values[13]-1904.0/531441.0*
fine_values[12]-1904.0/531441.0*fine_values[10]+196.0/531441.0*fine_values[9]
-8.0/6561.0*fine_values[8]-272.0/59049.0*fine_values[7]+28.0/59049.0*
fine_values[6]-272.0/59049.0*fine_values[5]-9248.0/531441.0*fine_values[4]+
952.0/531441.0*fine_values[3]+28.0/59049.0*fine_values[2]+952.0/531441.0*
fine_values[1];
      s3 = 18496.0/531441.0*fine_values[124]-1904.0/531441.0*fine_values[123]
-9248.0/531441.0*fine_values[122]+952.0/531441.0*fine_values[121]+196.0/
531441.0*fine_values[119]+952.0/531441.0*fine_values[118]-98.0/531441.0*
fine_values[117]+544.0/59049.0*fine_values[116]+1088.0/531441.0*fine_values
[114]-544.0/531441.0*fine_values[113]-544.0/531441.0*fine_values[112]+272.0/
531441.0*fine_values[111]-56.0/59049.0*fine_values[110]+28.0/59049.0*
fine_values[109]-112.0/531441.0*fine_values[108]+56.0/531441.0*fine_values[106]
;
      s2 = s3-28.0/531441.0*fine_values[105]+544.0/59049.0*fine_values[104]
-272.0/59049.0*fine_values[103]+18496.0/531441.0*fine_values[102]-9248.0/
531441.0*fine_values[101]-1904.0/531441.0*fine_values[100]+952.0/531441.0*
fine_values[99]-56.0/59049.0*fine_values[98]+28.0/59049.0*fine_values[97]
-1904.0/531441.0*fine_values[96]+952.0/531441.0*fine_values[95]+196.0/531441.0*
fine_values[94]-98.0/531441.0*fine_values[93]+272.0/6561.0*fine_values[92]+
9248.0/59049.0*fine_values[91]-952.0/59049.0*fine_values[90];
      s3 = s2+314432.0/531441.0*fine_values[88]-32368.0/531441.0*fine_values
[87]-952.0/59049.0*fine_values[86]-32368.0/531441.0*fine_values[85]+3332.0/
531441.0*fine_values[84]-28.0/6561.0*fine_values[83]-952.0/59049.0*fine_values
[80]+3332.0/531441.0*fine_values[78]+98.0/59049.0*fine_values[77]-343.0/
531441.0*fine_values[75]+544.0/59049.0*fine_values[74]-272.0/59049.0*
fine_values[73]-56.0/59049.0*fine_values[72]+28.0/59049.0*fine_values[71]+
1088.0/531441.0*fine_values[70];
      coarse_values[37] = s3-544.0/531441.0*fine_values[69]-112.0/531441.0*
fine_values[68]+56.0/531441.0*fine_values[67]-544.0/531441.0*fine_values[66]+
272.0/531441.0*fine_values[65]+56.0/531441.0*fine_values[64]-28.0/531441.0*
fine_values[63]+16.0/6561.0*fine_values[62]+32.0/59049.0*fine_values[61]-8.0/
6561.0*fine_values[59]-16.0/59049.0*fine_values[58]+8.0/59049.0*fine_values[57]
+64.0/531441.0*fine_values[55]-32.0/531441.0*fine_values[54]-32.0/531441.0*
fine_values[52]+s1-98.0/531441.0*fine_values[0];
      s2 = -28.0/6561.0*fine_values[41]-9248.0/531441.0*fine_values[31]-272.0/
59049.0*fine_values[21]-56.0/59049.0*fine_values[19]+16.0/6561.0*fine_values
[17]-56.0/59049.0*fine_values[11]-56.0/59049.0*fine_values[60]+32.0/59049.0*
fine_values[56]-16.0/59049.0*fine_values[53]+544.0/59049.0*fine_values[89]
-952.0/59049.0*fine_values[82]+98.0/59049.0*fine_values[81]-1904.0/531441.0*
fine_values[79]+952.0/531441.0*fine_values[76]-112.0/531441.0*fine_values[120]
-272.0/59049.0*fine_values[115]+952.0/531441.0*fine_values[107]+56.0/531441.0*
fine_values[51]-16.0/59049.0*fine_values[50]-544.0/531441.0*fine_values[49]+
56.0/531441.0*fine_values[48]+8.0/59049.0*fine_values[47]+272.0/531441.0*
fine_values[46]-28.0/531441.0*fine_values[45]+272.0/6561.0*fine_values[44]+
9248.0/59049.0*fine_values[43]-952.0/59049.0*fine_values[42]-952.0/59049.0*
fine_values[40]+98.0/59049.0*fine_values[39]+544.0/59049.0*fine_values[38]+
18496.0/531441.0*fine_values[37];
      s1 = -1904.0/531441.0*fine_values[36]-56.0/59049.0*fine_values[35]-1904.0
/531441.0*fine_values[34]+196.0/531441.0*fine_values[33]-272.0/59049.0*
fine_values[32]+952.0/531441.0*fine_values[30]+28.0/59049.0*fine_values[29]+
952.0/531441.0*fine_values[28]-98.0/531441.0*fine_values[27]+8.0/729.0*
fine_values[26]+16.0/6561.0*fine_values[25]-8.0/6561.0*fine_values[24]+272.0/
6561.0*fine_values[23]+s2+544.0/59049.0*fine_values[22]-28.0/6561.0*fine_values
[20]+28.0/59049.0*fine_values[18]+32.0/59049.0*fine_values[16]-16.0/59049.0*
fine_values[15]+544.0/59049.0*fine_values[14]+1088.0/531441.0*fine_values[13]
-544.0/531441.0*fine_values[12]-112.0/531441.0*fine_values[10]+56.0/531441.0*
fine_values[9]-8.0/6561.0*fine_values[8]-16.0/59049.0*fine_values[7]+8.0/
59049.0*fine_values[6]-272.0/59049.0*fine_values[5]-544.0/531441.0*fine_values
[4]+272.0/531441.0*fine_values[3]+28.0/59049.0*fine_values[2]+56.0/531441.0*
fine_values[1];
      s3 = 1088.0/531441.0*fine_values[124]-544.0/531441.0*fine_values[123]
-544.0/531441.0*fine_values[122]+272.0/531441.0*fine_values[121]+56.0/531441.0*
fine_values[119]+56.0/531441.0*fine_values[118]-28.0/531441.0*fine_values[117]+
544.0/59049.0*fine_values[116]+18496.0/531441.0*fine_values[114]-9248.0/
531441.0*fine_values[113]-1904.0/531441.0*fine_values[112]+952.0/531441.0*
fine_values[111]-56.0/59049.0*fine_values[110]+28.0/59049.0*fine_values[109]
-1904.0/531441.0*fine_values[108]+196.0/531441.0*fine_values[106];
      s2 = s3-98.0/531441.0*fine_values[105]+9248.0/59049.0*fine_values[104]
-952.0/59049.0*fine_values[103]+314432.0/531441.0*fine_values[102]-32368.0/
531441.0*fine_values[101]-32368.0/531441.0*fine_values[100]+3332.0/531441.0*
fine_values[99]-952.0/59049.0*fine_values[98]+98.0/59049.0*fine_values[97]
-32368.0/531441.0*fine_values[96]+3332.0/531441.0*fine_values[95]+3332.0/
531441.0*fine_values[94]-343.0/531441.0*fine_values[93]+272.0/6561.0*
fine_values[92]+9248.0/59049.0*fine_values[91]-952.0/59049.0*fine_values[90];
      coarse_values[38] = s2+18496.0/531441.0*fine_values[88]-1904.0/531441.0*
fine_values[87]-272.0/59049.0*fine_values[86]-9248.0/531441.0*fine_values[85]+
952.0/531441.0*fine_values[84]-28.0/6561.0*fine_values[83]-56.0/59049.0*
fine_values[80]+196.0/531441.0*fine_values[78]+28.0/59049.0*fine_values[77]
-98.0/531441.0*fine_values[75]+32.0/59049.0*fine_values[74]-16.0/59049.0*
fine_values[73]-16.0/59049.0*fine_values[72]+8.0/59049.0*fine_values[71]+64.0/
531441.0*fine_values[70]-32.0/531441.0*fine_values[69]-32.0/531441.0*
fine_values[68]+16.0/531441.0*fine_values[67]-32.0/531441.0*fine_values[66]+
16.0/531441.0*fine_values[65]+16.0/531441.0*fine_values[64]-8.0/531441.0*
fine_values[63]+16.0/6561.0*fine_values[62]+544.0/59049.0*fine_values[61]-8.0/
6561.0*fine_values[59]-272.0/59049.0*fine_values[58]+28.0/59049.0*fine_values
[57]+1088.0/531441.0*fine_values[55]-112.0/531441.0*fine_values[54]-544.0/
531441.0*fine_values[52]+s1-28.0/531441.0*fine_values[0];
      coarse_values[39] = 4.0/81.0*fine_values[41]+8.0/729.0*fine_values[60]
-8.0/6561.0*fine_values[51]-8.0/6561.0*fine_values[48]+4.0/6561.0*fine_values
[45]+136.0/729.0*fine_values[40]-14.0/729.0*fine_values[39]+8.0/729.0*
fine_values[35]+272.0/6561.0*fine_values[34]-28.0/6561.0*fine_values[33]-4.0/
729.0*fine_values[29]-136.0/6561.0*fine_values[28]+14.0/6561.0*fine_values[27]+
272.0/6561.0*fine_values[112]-136.0/6561.0*fine_values[111]-28.0/6561.0*
fine_values[106]+14.0/6561.0*fine_values[105]+136.0/729.0*fine_values[103]+
4624.0/6561.0*fine_values[101]-476.0/6561.0*fine_values[99]-14.0/729.0*
fine_values[97]-476.0/6561.0*fine_values[95]+49.0/6561.0*fine_values[93]-4.0/
729.0*fine_values[57]+16.0/6561.0*fine_values[54];
      coarse_values[40] = 8.0/729.0*fine_values[21]-28.0/6561.0*fine_values[76]
+4.0/81.0*fine_values[24]-4.0/729.0*fine_values[18]+8.0/729.0*fine_values[15]+
16.0/6561.0*fine_values[12]-8.0/6561.0*fine_values[9]-4.0/729.0*fine_values[6]
-8.0/6561.0*fine_values[3]+4624.0/6561.0*fine_values[123]-476.0/6561.0*
fine_values[121]-476.0/6561.0*fine_values[119]+49.0/6561.0*fine_values[117]+
136.0/729.0*fine_values[86]+272.0/6561.0*fine_values[85]-136.0/6561.0*
fine_values[84]-14.0/729.0*fine_values[77]+14.0/6561.0*fine_values[75]+136.0/
729.0*fine_values[72]-14.0/729.0*fine_values[71]+272.0/6561.0*fine_values[68]
-28.0/6561.0*fine_values[67]-136.0/6561.0*fine_values[64]+14.0/6561.0*
fine_values[63]+4.0/6561.0*fine_values[0];
      s2 = -8.0/6561.0*fine_values[41]-32.0/531441.0*fine_values[31]-56.0/
59049.0*fine_values[21]-272.0/59049.0*fine_values[19]+16.0/6561.0*fine_values
[17]-16.0/59049.0*fine_values[11]-272.0/59049.0*fine_values[60]+544.0/59049.0*
fine_values[56]-56.0/59049.0*fine_values[53]+9248.0/59049.0*fine_values[89]
-56.0/59049.0*fine_values[82]+28.0/59049.0*fine_values[81]-1904.0/531441.0*
fine_values[79]+196.0/531441.0*fine_values[76]-32368.0/531441.0*fine_values
[120]-952.0/59049.0*fine_values[115]+196.0/531441.0*fine_values[107]+56.0/
531441.0*fine_values[51]-272.0/59049.0*fine_values[50]-544.0/531441.0*
fine_values[49]+272.0/531441.0*fine_values[48]+28.0/59049.0*fine_values[47]+
56.0/531441.0*fine_values[46]-28.0/531441.0*fine_values[45]+16.0/6561.0*
fine_values[44]+32.0/59049.0*fine_values[43]-16.0/59049.0*fine_values[42]-16.0/
59049.0*fine_values[40]+8.0/59049.0*fine_values[39]+32.0/59049.0*fine_values
[38]+64.0/531441.0*fine_values[37];
      s1 = -32.0/531441.0*fine_values[36]-16.0/59049.0*fine_values[35]-32.0/
531441.0*fine_values[34]+16.0/531441.0*fine_values[33]-16.0/59049.0*fine_values
[32]+16.0/531441.0*fine_values[30]+8.0/59049.0*fine_values[29]+16.0/531441.0*
fine_values[28]-8.0/531441.0*fine_values[27]+8.0/729.0*fine_values[26]+272.0/
6561.0*fine_values[25]-28.0/6561.0*fine_values[24]+16.0/6561.0*fine_values[23]+
s2+544.0/59049.0*fine_values[22]-8.0/6561.0*fine_values[20]+28.0/59049.0*
fine_values[18]+544.0/59049.0*fine_values[16]-56.0/59049.0*fine_values[15]+32.0
/59049.0*fine_values[14]+1088.0/531441.0*fine_values[13]-112.0/531441.0*
fine_values[12]-544.0/531441.0*fine_values[10]+56.0/531441.0*fine_values[9]-8.0
/6561.0*fine_values[8]-272.0/59049.0*fine_values[7]+28.0/59049.0*fine_values[6]
-16.0/59049.0*fine_values[5]-544.0/531441.0*fine_values[4]+56.0/531441.0*
fine_values[3]+8.0/59049.0*fine_values[2]+272.0/531441.0*fine_values[1];
      s3 = 314432.0/531441.0*fine_values[124]-32368.0/531441.0*fine_values[123]
-32368.0/531441.0*fine_values[122]+3332.0/531441.0*fine_values[121]+3332.0/
531441.0*fine_values[119]+3332.0/531441.0*fine_values[118]-343.0/531441.0*
fine_values[117]+9248.0/59049.0*fine_values[116]+18496.0/531441.0*fine_values
[114]-1904.0/531441.0*fine_values[113]-9248.0/531441.0*fine_values[112]+952.0/
531441.0*fine_values[111]-952.0/59049.0*fine_values[110]+98.0/59049.0*
fine_values[109]-1904.0/531441.0*fine_values[108]+952.0/531441.0*fine_values
[106];
      s2 = s3-98.0/531441.0*fine_values[105]+544.0/59049.0*fine_values[104]
-272.0/59049.0*fine_values[103]+1088.0/531441.0*fine_values[102]-544.0/531441.0
*fine_values[101]-544.0/531441.0*fine_values[100]+272.0/531441.0*fine_values
[99]-56.0/59049.0*fine_values[98]+28.0/59049.0*fine_values[97]-112.0/531441.0*
fine_values[96]+56.0/531441.0*fine_values[95]+56.0/531441.0*fine_values[94]
-28.0/531441.0*fine_values[93]+272.0/6561.0*fine_values[92]+544.0/59049.0*
fine_values[91]-272.0/59049.0*fine_values[90];
      s3 = s2+18496.0/531441.0*fine_values[88]-9248.0/531441.0*fine_values[87]
-952.0/59049.0*fine_values[86]-1904.0/531441.0*fine_values[85]+952.0/531441.0*
fine_values[84]-28.0/6561.0*fine_values[83]-952.0/59049.0*fine_values[80]+952.0
/531441.0*fine_values[78]+98.0/59049.0*fine_values[77]-98.0/531441.0*
fine_values[75]+9248.0/59049.0*fine_values[74]-952.0/59049.0*fine_values[73]
-952.0/59049.0*fine_values[72]+98.0/59049.0*fine_values[71]+18496.0/531441.0*
fine_values[70];
      coarse_values[41] = s3-1904.0/531441.0*fine_values[69]-1904.0/531441.0*
fine_values[68]+196.0/531441.0*fine_values[67]-9248.0/531441.0*fine_values[66]+
952.0/531441.0*fine_values[65]+952.0/531441.0*fine_values[64]-98.0/531441.0*
fine_values[63]+272.0/6561.0*fine_values[62]+544.0/59049.0*fine_values[61]-28.0
/6561.0*fine_values[59]-56.0/59049.0*fine_values[58]+28.0/59049.0*fine_values
[57]+1088.0/531441.0*fine_values[55]-544.0/531441.0*fine_values[54]-112.0/
531441.0*fine_values[52]+s1-28.0/531441.0*fine_values[0];
      s2 = -28.0/6561.0*fine_values[41]-544.0/531441.0*fine_values[31]-16.0/
59049.0*fine_values[21]-16.0/59049.0*fine_values[19]+16.0/6561.0*fine_values
[17]-16.0/59049.0*fine_values[11]-952.0/59049.0*fine_values[60]+544.0/59049.0*
fine_values[56]-56.0/59049.0*fine_values[53]+544.0/59049.0*fine_values[89]-56.0
/59049.0*fine_values[82]+28.0/59049.0*fine_values[81]-112.0/531441.0*
fine_values[79]+56.0/531441.0*fine_values[76]-1904.0/531441.0*fine_values[120]
-952.0/59049.0*fine_values[115]+3332.0/531441.0*fine_values[107]+196.0/531441.0
*fine_values[51]-272.0/59049.0*fine_values[50]-9248.0/531441.0*fine_values[49]+
952.0/531441.0*fine_values[48]+28.0/59049.0*fine_values[47]+952.0/531441.0*
fine_values[46]-98.0/531441.0*fine_values[45]+272.0/6561.0*fine_values[44]+
544.0/59049.0*fine_values[43]-272.0/59049.0*fine_values[42]-56.0/59049.0*
fine_values[40]+28.0/59049.0*fine_values[39]+544.0/59049.0*fine_values[38]+
1088.0/531441.0*fine_values[37];
      s1 = -544.0/531441.0*fine_values[36]-56.0/59049.0*fine_values[35]-112.0/
531441.0*fine_values[34]+56.0/531441.0*fine_values[33]-272.0/59049.0*
fine_values[32]+272.0/531441.0*fine_values[30]+28.0/59049.0*fine_values[29]+
56.0/531441.0*fine_values[28]-28.0/531441.0*fine_values[27]+8.0/729.0*
fine_values[26]+16.0/6561.0*fine_values[25]-8.0/6561.0*fine_values[24]+16.0/
6561.0*fine_values[23]+s2+32.0/59049.0*fine_values[22]-8.0/6561.0*fine_values
[20]+8.0/59049.0*fine_values[18]+32.0/59049.0*fine_values[16]-16.0/59049.0*
fine_values[15]+32.0/59049.0*fine_values[14]+64.0/531441.0*fine_values[13]-32.0
/531441.0*fine_values[12]-32.0/531441.0*fine_values[10]+16.0/531441.0*
fine_values[9]-8.0/6561.0*fine_values[8]-16.0/59049.0*fine_values[7]+8.0/
59049.0*fine_values[6]-16.0/59049.0*fine_values[5]-32.0/531441.0*fine_values[4]
+16.0/531441.0*fine_values[3]+8.0/59049.0*fine_values[2]+16.0/531441.0*
fine_values[1];
      s3 = 18496.0/531441.0*fine_values[124]-9248.0/531441.0*fine_values[123]
-1904.0/531441.0*fine_values[122]+952.0/531441.0*fine_values[121]+952.0/
531441.0*fine_values[119]+196.0/531441.0*fine_values[118]-98.0/531441.0*
fine_values[117]+9248.0/59049.0*fine_values[116]+314432.0/531441.0*fine_values
[114]-32368.0/531441.0*fine_values[113]-32368.0/531441.0*fine_values[112]+
3332.0/531441.0*fine_values[111]-952.0/59049.0*fine_values[110]+98.0/59049.0*
fine_values[109]-32368.0/531441.0*fine_values[108]+3332.0/531441.0*fine_values
[106];
      s2 = s3-343.0/531441.0*fine_values[105]+9248.0/59049.0*fine_values[104]
-952.0/59049.0*fine_values[103]+18496.0/531441.0*fine_values[102]-1904.0/
531441.0*fine_values[101]-9248.0/531441.0*fine_values[100]+952.0/531441.0*
fine_values[99]-952.0/59049.0*fine_values[98]+98.0/59049.0*fine_values[97]
-1904.0/531441.0*fine_values[96]+196.0/531441.0*fine_values[95]+952.0/531441.0*
fine_values[94]-98.0/531441.0*fine_values[93]+272.0/6561.0*fine_values[92]+
544.0/59049.0*fine_values[91]-272.0/59049.0*fine_values[90];
      s3 = s2+1088.0/531441.0*fine_values[88]-544.0/531441.0*fine_values[87]
-272.0/59049.0*fine_values[86]-544.0/531441.0*fine_values[85]+272.0/531441.0*
fine_values[84]-28.0/6561.0*fine_values[83]-56.0/59049.0*fine_values[80]+56.0/
531441.0*fine_values[78]+28.0/59049.0*fine_values[77]-28.0/531441.0*fine_values
[75]+544.0/59049.0*fine_values[74]-56.0/59049.0*fine_values[73]-272.0/59049.0*
fine_values[72]+28.0/59049.0*fine_values[71]+1088.0/531441.0*fine_values[70];
      coarse_values[42] = s3-112.0/531441.0*fine_values[69]-544.0/531441.0*
fine_values[68]+56.0/531441.0*fine_values[67]-544.0/531441.0*fine_values[66]+
56.0/531441.0*fine_values[65]+272.0/531441.0*fine_values[64]-28.0/531441.0*
fine_values[63]+272.0/6561.0*fine_values[62]+9248.0/59049.0*fine_values[61]
-28.0/6561.0*fine_values[59]-952.0/59049.0*fine_values[58]+98.0/59049.0*
fine_values[57]+18496.0/531441.0*fine_values[55]-1904.0/531441.0*fine_values
[54]-1904.0/531441.0*fine_values[52]+s1-8.0/531441.0*fine_values[0];
      coarse_values[43] = 4.0/81.0*fine_values[41]+136.0/729.0*fine_values[60]
-28.0/6561.0*fine_values[51]-136.0/6561.0*fine_values[48]+14.0/6561.0*
fine_values[45]+8.0/729.0*fine_values[40]-4.0/729.0*fine_values[39]+8.0/729.0*
fine_values[35]+16.0/6561.0*fine_values[34]-8.0/6561.0*fine_values[33]-4.0/
729.0*fine_values[29]-8.0/6561.0*fine_values[28]+4.0/6561.0*fine_values[27]+
4624.0/6561.0*fine_values[112]-476.0/6561.0*fine_values[111]-476.0/6561.0*
fine_values[106]+49.0/6561.0*fine_values[105]+136.0/729.0*fine_values[103]+
272.0/6561.0*fine_values[101]-136.0/6561.0*fine_values[99]-14.0/729.0*
fine_values[97]-28.0/6561.0*fine_values[95]+14.0/6561.0*fine_values[93]-14.0/
729.0*fine_values[57]+272.0/6561.0*fine_values[54];
      coarse_values[44] = 68.0/81.0*fine_values[121]-7.0/81.0*fine_values[117]+
2.0/9.0*fine_values[71]+4.0/81.0*fine_values[67]-2.0/81.0*fine_values[63];
      coarse_values[45] = 8.0/729.0*fine_values[53]+136.0/729.0*fine_values
[115]-28.0/6561.0*fine_values[107]-8.0/6561.0*fine_values[51]-4.0/729.0*
fine_values[47]-8.0/6561.0*fine_values[46]+4.0/6561.0*fine_values[45]+4624.0/
6561.0*fine_values[122]-476.0/6561.0*fine_values[121]-476.0/6561.0*fine_values
[118]+49.0/6561.0*fine_values[117]+272.0/6561.0*fine_values[113]-136.0/6561.0*
fine_values[111]-14.0/729.0*fine_values[109]+14.0/6561.0*fine_values[105]+136.0
/729.0*fine_values[73]-14.0/729.0*fine_values[71]+272.0/6561.0*fine_values[69]
-28.0/6561.0*fine_values[67]-136.0/6561.0*fine_values[65]+14.0/6561.0*
fine_values[63]+4.0/81.0*fine_values[59]+8.0/729.0*fine_values[58]-4.0/729.0*
fine_values[57]+16.0/6561.0*fine_values[52];
      coarse_values[46] = 8.0/729.0*fine_values[53]+136.0/729.0*fine_values
[115]-476.0/6561.0*fine_values[107]-28.0/6561.0*fine_values[51]-4.0/729.0*
fine_values[47]-136.0/6561.0*fine_values[46]+14.0/6561.0*fine_values[45]+272.0/
6561.0*fine_values[122]-136.0/6561.0*fine_values[121]-28.0/6561.0*fine_values
[118]+14.0/6561.0*fine_values[117]+4624.0/6561.0*fine_values[113]-476.0/6561.0*
fine_values[111]-14.0/729.0*fine_values[109]+49.0/6561.0*fine_values[105]+8.0/
729.0*fine_values[73]-4.0/729.0*fine_values[71]+16.0/6561.0*fine_values[69]-8.0
/6561.0*fine_values[67]-8.0/6561.0*fine_values[65]+4.0/6561.0*fine_values[63]+
4.0/81.0*fine_values[59]+136.0/729.0*fine_values[58]-14.0/729.0*fine_values[57]
+272.0/6561.0*fine_values[52];
      coarse_values[47] = -2.0/81.0*fine_values[45]+68.0/81.0*fine_values[111]
-7.0/81.0*fine_values[105]+2.0/9.0*fine_values[57]+4.0/81.0*fine_values[51];
      coarse_values[48] = fine_values[75];
      coarse_values[49] = 2.0/9.0*fine_values[81]+68.0/81.0*fine_values[78]-7.0
/81.0*fine_values[75]+4.0/81.0*fine_values[94]-2.0/81.0*fine_values[93];
      coarse_values[50] = 2.0/9.0*fine_values[81]+4.0/81.0*fine_values[78]-2.0/
81.0*fine_values[75]+68.0/81.0*fine_values[94]-7.0/81.0*fine_values[93];
      coarse_values[51] = fine_values[93];
      coarse_values[52] = -7.0/81.0*fine_values[75]+4.0/81.0*fine_values[119]
-2.0/81.0*fine_values[117]+2.0/9.0*fine_values[77]+68.0/81.0*fine_values[76];
      coarse_values[53] = 136.0/729.0*fine_values[82]-14.0/729.0*fine_values
[81]+4624.0/6561.0*fine_values[79]-476.0/6561.0*fine_values[76]+272.0/6561.0*
fine_values[120]-8.0/6561.0*fine_values[107]-28.0/6561.0*fine_values[119]-136.0
/6561.0*fine_values[118]+14.0/6561.0*fine_values[117]+8.0/729.0*fine_values
[110]-4.0/729.0*fine_values[109]+16.0/6561.0*fine_values[108]-8.0/6561.0*
fine_values[106]+4.0/6561.0*fine_values[105]+8.0/729.0*fine_values[98]-4.0/
729.0*fine_values[97]+272.0/6561.0*fine_values[96]-136.0/6561.0*fine_values[95]
-28.0/6561.0*fine_values[94]+14.0/6561.0*fine_values[93]+4.0/81.0*fine_values
[83]+136.0/729.0*fine_values[80]-476.0/6561.0*fine_values[78]-14.0/729.0*
fine_values[77]+49.0/6561.0*fine_values[75];
      coarse_values[54] = 136.0/729.0*fine_values[82]-14.0/729.0*fine_values
[81]+272.0/6561.0*fine_values[79]-136.0/6561.0*fine_values[76]+16.0/6561.0*
fine_values[120]-136.0/6561.0*fine_values[107]-8.0/6561.0*fine_values[119]-8.0/
6561.0*fine_values[118]+4.0/6561.0*fine_values[117]+8.0/729.0*fine_values[110]
-4.0/729.0*fine_values[109]+272.0/6561.0*fine_values[108]-28.0/6561.0*
fine_values[106]+14.0/6561.0*fine_values[105]+136.0/729.0*fine_values[98]-14.0/
729.0*fine_values[97]+4624.0/6561.0*fine_values[96]-476.0/6561.0*fine_values
[95]-476.0/6561.0*fine_values[94]+49.0/6561.0*fine_values[93]+4.0/81.0*
fine_values[83]+8.0/729.0*fine_values[80]-28.0/6561.0*fine_values[78]-4.0/729.0
*fine_values[77]+14.0/6561.0*fine_values[75];
      coarse_values[55] = 4.0/81.0*fine_values[106]-2.0/81.0*fine_values[105]+
2.0/9.0*fine_values[97]+68.0/81.0*fine_values[95]-7.0/81.0*fine_values[93];
      coarse_values[56] = -2.0/81.0*fine_values[75]+68.0/81.0*fine_values[119]
-7.0/81.0*fine_values[117]+2.0/9.0*fine_values[77]+4.0/81.0*fine_values[76];
      coarse_values[57] = 8.0/729.0*fine_values[82]-4.0/729.0*fine_values[81]+
272.0/6561.0*fine_values[79]-28.0/6561.0*fine_values[76]+4624.0/6561.0*
fine_values[120]-28.0/6561.0*fine_values[107]-476.0/6561.0*fine_values[119]
-476.0/6561.0*fine_values[118]+49.0/6561.0*fine_values[117]+136.0/729.0*
fine_values[110]-14.0/729.0*fine_values[109]+272.0/6561.0*fine_values[108]
-136.0/6561.0*fine_values[106]+14.0/6561.0*fine_values[105]+8.0/729.0*
fine_values[98]-4.0/729.0*fine_values[97]+16.0/6561.0*fine_values[96]-8.0/
6561.0*fine_values[95]-8.0/6561.0*fine_values[94]+4.0/6561.0*fine_values[93]+
4.0/81.0*fine_values[83]+136.0/729.0*fine_values[80]-136.0/6561.0*fine_values
[78]-14.0/729.0*fine_values[77]+14.0/6561.0*fine_values[75];
      coarse_values[58] = 8.0/729.0*fine_values[82]-4.0/729.0*fine_values[81]+
16.0/6561.0*fine_values[79]-8.0/6561.0*fine_values[76]+272.0/6561.0*fine_values
[120]-476.0/6561.0*fine_values[107]-136.0/6561.0*fine_values[119]-28.0/6561.0*
fine_values[118]+14.0/6561.0*fine_values[117]+136.0/729.0*fine_values[110]-14.0
/729.0*fine_values[109]+4624.0/6561.0*fine_values[108]-476.0/6561.0*fine_values
[106]+49.0/6561.0*fine_values[105]+136.0/729.0*fine_values[98]-14.0/729.0*
fine_values[97]+272.0/6561.0*fine_values[96]-28.0/6561.0*fine_values[95]-136.0/
6561.0*fine_values[94]+14.0/6561.0*fine_values[93]+4.0/81.0*fine_values[83]+8.0
/729.0*fine_values[80]-8.0/6561.0*fine_values[78]-4.0/729.0*fine_values[77]+4.0
/6561.0*fine_values[75];
      coarse_values[59] = 68.0/81.0*fine_values[106]-7.0/81.0*fine_values[105]+
2.0/9.0*fine_values[97]+4.0/81.0*fine_values[95]-2.0/81.0*fine_values[93];
      coarse_values[60] = fine_values[117];
      coarse_values[61] = 68.0/81.0*fine_values[118]-7.0/81.0*fine_values[117]+
2.0/9.0*fine_values[109]+4.0/81.0*fine_values[107]-2.0/81.0*fine_values[105];
      coarse_values[62] = 4.0/81.0*fine_values[118]-2.0/81.0*fine_values[117]+
2.0/9.0*fine_values[109]+68.0/81.0*fine_values[107]-7.0/81.0*fine_values[105];
      coarse_values[63] = fine_values[105];
 
 
}


void Transform_Q2Q4_3D(double *fine_values, double *coarse_values)
{

  double s1,s2,s3;

  coarse_values[0] = fine_values[0];
  
  coarse_values[1] = 3.0/64.0*fine_values[2]
                     +31.0/32.0*fine_values[1]
                     +fine_values[0]/128.0
                     -fine_values[30]/32.0
                     +fine_values[27]/128.0;
  
  coarse_values[2] = fine_values[2];
  
  coarse_values[3] = 3.0/64.0*fine_values[2]-fine_values[1]/32.0+
fine_values[0]/128.0+31.0/32.0*fine_values[30]+fine_values[27]/128.0;
      coarse_values[4] = fine_values[27];
      coarse_values[5] = 3.0/64.0*fine_values[6]+31.0/32.0*fine_values[3]+
fine_values[0]/128.0-fine_values[64]/32.0+fine_values[63]/128.0;
      coarse_values[6] = 93.0/2048.0*fine_values[7]+9.0/4096.0*fine_values[8]+
93.0/2048.0*fine_values[5]+3.0/8192.0*fine_values[6]+961.0/1024.0*fine_values
[4]+31.0/4096.0*fine_values[3]+3.0/8192.0*fine_values[2]+31.0/4096.0*
fine_values[1]+fine_values[49]/1024.0+3.0/8192.0*fine_values[47]-fine_values
[48]/4096.0+fine_values[45]/16384.0-fine_values[46]/4096.0-31.0/1024.0*
fine_values[31]-3.0/2048.0*fine_values[32]-fine_values[30]/4096.0+31.0/4096.0*
fine_values[28]+3.0/8192.0*fine_values[29]+fine_values[27]/16384.0+fine_values
[0]/16384.0-3.0/2048.0*fine_values[50]+31.0/4096.0*fine_values[65]-31.0/1024.0*
fine_values[66]-fine_values[64]/4096.0+fine_values[63]/16384.0;
      coarse_values[7] = 3.0/64.0*fine_values[8]+31.0/32.0*fine_values[5]+
fine_values[2]/128.0+fine_values[47]/128.0-fine_values[50]/32.0;
      coarse_values[8] = -3.0/2048.0*fine_values[7]+9.0/4096.0*fine_values[8]+
93.0/2048.0*fine_values[5]+3.0/8192.0*fine_values[6]-31.0/1024.0*fine_values[4]
+31.0/4096.0*fine_values[3]+3.0/8192.0*fine_values[2]-fine_values[1]/4096.0
-31.0/1024.0*fine_values[49]+3.0/8192.0*fine_values[47]-fine_values[48]/4096.0+
fine_values[45]/16384.0+31.0/4096.0*fine_values[46]+961.0/1024.0*fine_values
[31]+93.0/2048.0*fine_values[32]+31.0/4096.0*fine_values[30]+31.0/4096.0*
fine_values[28]+3.0/8192.0*fine_values[29]+fine_values[27]/16384.0+fine_values
[0]/16384.0-3.0/2048.0*fine_values[50]-fine_values[65]/4096.0+fine_values[66]/
1024.0-fine_values[64]/4096.0+fine_values[63]/16384.0;
      coarse_values[9] = -fine_values[48]/32.0+fine_values[45]/128.0+3.0/64.0*
fine_values[29]+fine_values[27]/128.0+31.0/32.0*fine_values[28];
      coarse_values[10] = fine_values[6];
      coarse_values[11] = 3.0/64.0*fine_values[8]+fine_values[6]/128.0+31.0/
32.0*fine_values[7]-fine_values[32]/32.0+fine_values[29]/128.0;
      coarse_values[12] = fine_values[8];
      coarse_values[13] = 3.0/64.0*fine_values[8]+fine_values[6]/128.0-
fine_values[7]/32.0+31.0/32.0*fine_values[32]+fine_values[29]/128.0;
      coarse_values[14] = fine_values[29];
      coarse_values[15] = 3.0/64.0*fine_values[6]-fine_values[3]/32.0+
fine_values[0]/128.0+31.0/32.0*fine_values[64]+fine_values[63]/128.0;
      coarse_values[16] = 93.0/2048.0*fine_values[7]+9.0/4096.0*fine_values[8]
-3.0/2048.0*fine_values[5]+3.0/8192.0*fine_values[6]-31.0/1024.0*fine_values[4]
-fine_values[3]/4096.0+3.0/8192.0*fine_values[2]+31.0/4096.0*fine_values[1]
-31.0/1024.0*fine_values[49]+3.0/8192.0*fine_values[47]+31.0/4096.0*fine_values
[48]+fine_values[45]/16384.0-fine_values[46]/4096.0+fine_values[31]/1024.0-3.0/
2048.0*fine_values[32]-fine_values[30]/4096.0-fine_values[28]/4096.0+3.0/8192.0
*fine_values[29]+fine_values[27]/16384.0+fine_values[0]/16384.0+93.0/2048.0*
fine_values[50]+31.0/4096.0*fine_values[65]+961.0/1024.0*fine_values[66]+31.0/
4096.0*fine_values[64]+fine_values[63]/16384.0;
      coarse_values[17] = 3.0/64.0*fine_values[8]-fine_values[5]/32.0+
fine_values[2]/128.0+fine_values[47]/128.0+31.0/32.0*fine_values[50];
      coarse_values[18] = -3.0/2048.0*fine_values[7]+9.0/4096.0*fine_values[8]
-3.0/2048.0*fine_values[5]+3.0/8192.0*fine_values[6]+fine_values[4]/1024.0-
fine_values[3]/4096.0+3.0/8192.0*fine_values[2]-fine_values[1]/4096.0+961.0/
1024.0*fine_values[49]+3.0/8192.0*fine_values[47]+31.0/4096.0*fine_values[48]+
fine_values[45]/16384.0+31.0/4096.0*fine_values[46]-31.0/1024.0*fine_values[31]
+93.0/2048.0*fine_values[32]+31.0/4096.0*fine_values[30]-fine_values[28]/4096.0
+3.0/8192.0*fine_values[29]+fine_values[27]/16384.0+fine_values[0]/16384.0+93.0
/2048.0*fine_values[50]-fine_values[65]/4096.0-31.0/1024.0*fine_values[66]+31.0
/4096.0*fine_values[64]+fine_values[63]/16384.0;
      coarse_values[19] = 31.0/32.0*fine_values[48]+fine_values[45]/128.0+3.0/
64.0*fine_values[29]+fine_values[27]/128.0-fine_values[28]/32.0;
      coarse_values[20] = fine_values[63];
      coarse_values[21] = -fine_values[46]/32.0+3.0/64.0*fine_values[47]+
fine_values[45]/128.0+31.0/32.0*fine_values[65]+fine_values[63]/128.0;
      coarse_values[22] = fine_values[47];
      coarse_values[23] = 31.0/32.0*fine_values[46]+3.0/64.0*fine_values[47]+
fine_values[45]/128.0-fine_values[65]/32.0+fine_values[63]/128.0;
      coarse_values[24] = fine_values[45];
      coarse_values[25] = 3.0/64.0*fine_values[18]+31.0/32.0*fine_values[9]+
fine_values[0]/128.0-fine_values[84]/32.0+fine_values[75]/128.0;
      coarse_values[26] = 9.0/4096.0*fine_values[20]+93.0/2048.0*fine_values
[19]+3.0/8192.0*fine_values[18]+93.0/2048.0*fine_values[11]+31.0/4096.0*
fine_values[9]+961.0/1024.0*fine_values[10]+3.0/8192.0*fine_values[2]+31.0/
4096.0*fine_values[1]-3.0/2048.0*fine_values[42]+3.0/8192.0*fine_values[39]
-31.0/1024.0*fine_values[36]+31.0/4096.0*fine_values[33]-fine_values[30]/4096.0
+fine_values[27]/16384.0+fine_values[0]/16384.0-fine_values[99]/4096.0+
fine_values[100]/1024.0-fine_values[94]/4096.0+fine_values[93]/16384.0-3.0/
2048.0*fine_values[90]-31.0/1024.0*fine_values[87]-fine_values[84]/4096.0+3.0/
8192.0*fine_values[81]+31.0/4096.0*fine_values[78]+fine_values[75]/16384.0;
      coarse_values[27] = 3.0/64.0*fine_values[20]+31.0/32.0*fine_values[11]+
fine_values[2]/128.0-fine_values[90]/32.0+fine_values[81]/128.0;
      coarse_values[28] = 9.0/4096.0*fine_values[20]-3.0/2048.0*fine_values[19]
+3.0/8192.0*fine_values[18]+93.0/2048.0*fine_values[11]+31.0/4096.0*fine_values
[9]-31.0/1024.0*fine_values[10]+3.0/8192.0*fine_values[2]-fine_values[1]/4096.0
+93.0/2048.0*fine_values[42]+3.0/8192.0*fine_values[39]+961.0/1024.0*
fine_values[36]+31.0/4096.0*fine_values[33]+31.0/4096.0*fine_values[30]+
fine_values[27]/16384.0+fine_values[0]/16384.0-fine_values[99]/4096.0-31.0/
1024.0*fine_values[100]+31.0/4096.0*fine_values[94]+fine_values[93]/16384.0-3.0
/2048.0*fine_values[90]+fine_values[87]/1024.0-fine_values[84]/4096.0+3.0/
8192.0*fine_values[81]-fine_values[78]/4096.0+fine_values[75]/16384.0;
      coarse_values[29] = 3.0/64.0*fine_values[39]+31.0/32.0*fine_values[33]+
fine_values[27]/128.0-fine_values[99]/32.0+fine_values[93]/128.0;
      coarse_values[30] = 93.0/2048.0*fine_values[21]+3.0/8192.0*fine_values
[18]+93.0/2048.0*fine_values[15]+961.0/1024.0*fine_values[12]+31.0/4096.0*
fine_values[9]+3.0/8192.0*fine_values[6]+31.0/4096.0*fine_values[3]+9.0/4096.0*
fine_values[24]+fine_values[0]/16384.0+3.0/8192.0*fine_values[71]-3.0/2048.0*
fine_values[72]-31.0/1024.0*fine_values[68]+31.0/4096.0*fine_values[67]-
fine_values[64]/4096.0+fine_values[63]/16384.0-3.0/2048.0*fine_values[86]-31.0/
1024.0*fine_values[85]-fine_values[84]/4096.0+fine_values[75]/16384.0+31.0/
4096.0*fine_values[76]+3.0/8192.0*fine_values[77]+fine_values[123]/1024.0-
fine_values[121]/4096.0+fine_values[117]/16384.0-fine_values[119]/4096.0;
      s2 = 93.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]+93.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]+93.0/262144.0*
fine_values[15]+2883.0/65536.0*fine_values[16]+279.0/131072.0*fine_values[17]+
29791.0/32768.0*fine_values[13]+2883.0/65536.0*fine_values[14]+93.0/262144.0*
fine_values[11]+961.0/131072.0*fine_values[12]+31.0/524288.0*fine_values[9]+
961.0/131072.0*fine_values[10]+93.0/262144.0*fine_values[7]+9.0/524288.0*
fine_values[8]+93.0/262144.0*fine_values[5]+3.0/1048576.0*fine_values[6]+961.0/
131072.0*fine_values[4]+31.0/524288.0*fine_values[3]+3.0/1048576.0*fine_values
[2]+31.0/524288.0*fine_values[1]+fine_values[49]/131072.0+3.0/1048576.0*
fine_values[47]-fine_values[48]/524288.0+fine_values[45]/2097152.0-fine_values
[46]/524288.0-93.0/65536.0*fine_values[43]-9.0/131072.0*fine_values[44]+9.0/
524288.0*fine_values[41]-3.0/262144.0*fine_values[42]+93.0/262144.0*fine_values
[40];
      s1 = s2+3.0/1048576.0*fine_values[39]-93.0/65536.0*fine_values[38]-961.0/
32768.0*fine_values[37]+961.0/131072.0*fine_values[34]-31.0/131072.0*
fine_values[36]+93.0/262144.0*fine_values[35]-31.0/131072.0*fine_values[31]-3.0
/262144.0*fine_values[32]+31.0/524288.0*fine_values[33]-fine_values[30]/
524288.0+31.0/524288.0*fine_values[28]+3.0/1048576.0*fine_values[29]+27.0/
262144.0*fine_values[26]+fine_values[27]/2097152.0+9.0/524288.0*fine_values[24]
+279.0/131072.0*fine_values[25]+fine_values[0]/2097152.0+279.0/131072.0*
fine_values[23]+2883.0/65536.0*fine_values[22]+9.0/524288.0*fine_values[59]-3.0
/262144.0*fine_values[58]+3.0/1048576.0*fine_values[57]-93.0/65536.0*
fine_values[56]+93.0/262144.0*fine_values[53]+31.0/32768.0*fine_values[55]-31.0
/131072.0*fine_values[54]-31.0/131072.0*fine_values[52]-3.0/262144.0*
fine_values[50]+31.0/524288.0*fine_values[51]+93.0/262144.0*fine_values[73]
-93.0/65536.0*fine_values[74];
      s2 = s1+3.0/1048576.0*fine_values[71]-3.0/262144.0*fine_values[72]-961.0/
32768.0*fine_values[70]+961.0/131072.0*fine_values[69]-31.0/131072.0*
fine_values[68]+31.0/524288.0*fine_values[65]-31.0/131072.0*fine_values[66]+
31.0/524288.0*fine_values[67]-fine_values[64]/524288.0+3.0/65536.0*fine_values
[61]-9.0/131072.0*fine_values[62]+fine_values[63]/2097152.0-31.0/131072.0*
fine_values[101]+31.0/32768.0*fine_values[102]-fine_values[99]/524288.0+
fine_values[100]/131072.0-3.0/262144.0*fine_values[98]-31.0/131072.0*
fine_values[96]+3.0/1048576.0*fine_values[97]-fine_values[94]/524288.0+31.0/
524288.0*fine_values[95]-93.0/65536.0*fine_values[91]+fine_values[93]/2097152.0
-9.0/131072.0*fine_values[92]-93.0/65536.0*fine_values[89]-3.0/262144.0*
fine_values[90]-961.0/32768.0*fine_values[88]-31.0/131072.0*fine_values[87]-3.0
/262144.0*fine_values[86]-31.0/131072.0*fine_values[85]+93.0/262144.0*
fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]-fine_values[84]/524288.0+3.0/
1048576.0*fine_values[81]+93.0/262144.0*fine_values[80]+961.0/131072.0*
fine_values[79]+31.0/524288.0*fine_values[78]+fine_values[75]/2097152.0+31.0/
524288.0*fine_values[76]+3.0/1048576.0*fine_values[77]-3.0/262144.0*fine_values
[60]+fine_values[112]/131072.0-fine_values[111]/524288.0+3.0/1048576.0*
fine_values[109]-3.0/262144.0*fine_values[110]+fine_values[105]/2097152.0;
      coarse_values[31] = s3-fine_values[107]/524288.0+fine_values[108]/
131072.0-fine_values[106]/524288.0-3.0/262144.0*fine_values[103]+3.0/65536.0*
fine_values[104]+31.0/32768.0*fine_values[124]-31.0/131072.0*fine_values[122]+
fine_values[123]/131072.0-31.0/131072.0*fine_values[120]-fine_values[121]/
524288.0+fine_values[117]/2097152.0+3.0/65536.0*fine_values[116]-fine_values
[119]/524288.0+31.0/524288.0*fine_values[118]-3.0/262144.0*fine_values[115]+
fine_values[113]/131072.0-fine_values[114]/32768.0;
      coarse_values[32] = 3.0/8192.0*fine_values[20]+93.0/2048.0*fine_values
[17]+961.0/1024.0*fine_values[14]+31.0/4096.0*fine_values[11]+3.0/8192.0*
fine_values[8]+31.0/4096.0*fine_values[5]+fine_values[2]/16384.0+fine_values
[47]/16384.0+9.0/4096.0*fine_values[26]+93.0/2048.0*fine_values[23]+3.0/8192.0*
fine_values[59]-31.0/1024.0*fine_values[56]+31.0/4096.0*fine_values[53]-
fine_values[50]/4096.0-3.0/2048.0*fine_values[62]-31.0/1024.0*fine_values[91]
-3.0/2048.0*fine_values[92]-fine_values[90]/4096.0+31.0/4096.0*fine_values[82]+
3.0/8192.0*fine_values[83]+fine_values[81]/16384.0+fine_values[109]/16384.0-
fine_values[110]/4096.0+fine_values[116]/1024.0-fine_values[115]/4096.0;
      s2 = 93.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]-3.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]+93.0/262144.0*
fine_values[15]-93.0/65536.0*fine_values[16]+279.0/131072.0*fine_values[17]
-961.0/32768.0*fine_values[13]+2883.0/65536.0*fine_values[14]+93.0/262144.0*
fine_values[11]+961.0/131072.0*fine_values[12]+31.0/524288.0*fine_values[9]
-31.0/131072.0*fine_values[10]-3.0/262144.0*fine_values[7]+9.0/524288.0*
fine_values[8]+93.0/262144.0*fine_values[5]+3.0/1048576.0*fine_values[6]-31.0/
131072.0*fine_values[4]+31.0/524288.0*fine_values[3]+3.0/1048576.0*fine_values
[2]-fine_values[1]/524288.0-31.0/131072.0*fine_values[49]+3.0/1048576.0*
fine_values[47]-fine_values[48]/524288.0+fine_values[45]/2097152.0+31.0/
524288.0*fine_values[46]+2883.0/65536.0*fine_values[43]+279.0/131072.0*
fine_values[44]+9.0/524288.0*fine_values[41]+93.0/262144.0*fine_values[42]+93.0
/262144.0*fine_values[40];
      s1 = s2+3.0/1048576.0*fine_values[39]+2883.0/65536.0*fine_values[38]+
29791.0/32768.0*fine_values[37]+961.0/131072.0*fine_values[34]+961.0/131072.0*
fine_values[36]+93.0/262144.0*fine_values[35]+961.0/131072.0*fine_values[31]+
93.0/262144.0*fine_values[32]+31.0/524288.0*fine_values[33]+31.0/524288.0*
fine_values[30]+31.0/524288.0*fine_values[28]+3.0/1048576.0*fine_values[29]+
27.0/262144.0*fine_values[26]+fine_values[27]/2097152.0+9.0/524288.0*
fine_values[24]-9.0/131072.0*fine_values[25]+fine_values[0]/2097152.0+279.0/
131072.0*fine_values[23]-93.0/65536.0*fine_values[22]+9.0/524288.0*fine_values
[59]+93.0/262144.0*fine_values[58]+3.0/1048576.0*fine_values[57]-93.0/65536.0*
fine_values[56]+93.0/262144.0*fine_values[53]-961.0/32768.0*fine_values[55]
-31.0/131072.0*fine_values[54]+961.0/131072.0*fine_values[52]-3.0/262144.0*
fine_values[50]+31.0/524288.0*fine_values[51]-3.0/262144.0*fine_values[73]+3.0/
65536.0*fine_values[74];
      s2 = s1+3.0/1048576.0*fine_values[71]-3.0/262144.0*fine_values[72]+31.0/
32768.0*fine_values[70]-31.0/131072.0*fine_values[69]-31.0/131072.0*fine_values
[68]-fine_values[65]/524288.0+fine_values[66]/131072.0+31.0/524288.0*
fine_values[67]-fine_values[64]/524288.0-93.0/65536.0*fine_values[61]-9.0/
131072.0*fine_values[62]+fine_values[63]/2097152.0-31.0/131072.0*fine_values
[101]-961.0/32768.0*fine_values[102]-fine_values[99]/524288.0-31.0/131072.0*
fine_values[100]+93.0/262144.0*fine_values[98]+961.0/131072.0*fine_values[96]+
3.0/1048576.0*fine_values[97]+31.0/524288.0*fine_values[94]+31.0/524288.0*
fine_values[95]-93.0/65536.0*fine_values[91]+fine_values[93]/2097152.0-9.0/
131072.0*fine_values[92]+3.0/65536.0*fine_values[89]-3.0/262144.0*fine_values
[90]+31.0/32768.0*fine_values[88]+fine_values[87]/131072.0-3.0/262144.0*
fine_values[86]-31.0/131072.0*fine_values[85]+93.0/262144.0*fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]-fine_values[84]/524288.0+3.0/
1048576.0*fine_values[81]-3.0/262144.0*fine_values[80]-31.0/131072.0*
fine_values[79]-fine_values[78]/524288.0+fine_values[75]/2097152.0+31.0/
524288.0*fine_values[76]+3.0/1048576.0*fine_values[77]-3.0/262144.0*fine_values
[60]+fine_values[112]/131072.0-fine_values[111]/524288.0+3.0/1048576.0*
fine_values[109]-3.0/262144.0*fine_values[110]+fine_values[105]/2097152.0;
      coarse_values[33] = s3+31.0/524288.0*fine_values[107]-31.0/131072.0*
fine_values[108]-fine_values[106]/524288.0-3.0/262144.0*fine_values[103]-93.0/
65536.0*fine_values[104]-fine_values[124]/32768.0+fine_values[122]/131072.0+
fine_values[123]/131072.0+fine_values[120]/131072.0-fine_values[121]/524288.0+
fine_values[117]/2097152.0+3.0/65536.0*fine_values[116]-fine_values[119]/
524288.0-fine_values[118]/524288.0-3.0/262144.0*fine_values[115]-31.0/131072.0*
fine_values[113]+31.0/32768.0*fine_values[114];
      coarse_values[34] = -fine_values[48]/4096.0+fine_values[45]/16384.0+9.0/
4096.0*fine_values[41]+93.0/2048.0*fine_values[40]+3.0/8192.0*fine_values[39]+
961.0/1024.0*fine_values[34]+93.0/2048.0*fine_values[35]+31.0/4096.0*
fine_values[33]+31.0/4096.0*fine_values[28]+3.0/8192.0*fine_values[29]+
fine_values[27]/16384.0+3.0/8192.0*fine_values[57]-31.0/1024.0*fine_values[54]+
31.0/4096.0*fine_values[51]-31.0/1024.0*fine_values[101]-fine_values[99]/4096.0
+3.0/8192.0*fine_values[97]+31.0/4096.0*fine_values[95]+fine_values[93]/16384.0
-3.0/2048.0*fine_values[60]+fine_values[112]/1024.0-fine_values[111]/4096.0+
fine_values[105]/16384.0-fine_values[106]/4096.0-3.0/2048.0*fine_values[103];
      coarse_values[35] = 31.0/32.0*fine_values[15]+fine_values[6]/128.0+3.0/
64.0*fine_values[24]-fine_values[86]/32.0+fine_values[77]/128.0;
      coarse_values[36] = 31.0/4096.0*fine_values[15]+961.0/1024.0*fine_values
[16]+93.0/2048.0*fine_values[17]+31.0/4096.0*fine_values[7]+3.0/8192.0*
fine_values[8]+fine_values[6]/16384.0-3.0/2048.0*fine_values[44]+3.0/8192.0*
fine_values[41]-31.0/1024.0*fine_values[38]+31.0/4096.0*fine_values[35]-
fine_values[32]/4096.0+fine_values[29]/16384.0+9.0/4096.0*fine_values[26]+3.0/
8192.0*fine_values[24]+93.0/2048.0*fine_values[25]-fine_values[98]/4096.0+
fine_values[97]/16384.0-3.0/2048.0*fine_values[92]-31.0/1024.0*fine_values[89]-
fine_values[86]/4096.0+3.0/8192.0*fine_values[83]+31.0/4096.0*fine_values[80]+
fine_values[77]/16384.0-fine_values[103]/4096.0+fine_values[104]/1024.0;
      coarse_values[37] = 31.0/32.0*fine_values[17]+fine_values[8]/128.0+3.0/
64.0*fine_values[26]-fine_values[92]/32.0+fine_values[83]/128.0;
      coarse_values[38] = 31.0/4096.0*fine_values[15]-31.0/1024.0*fine_values
[16]+93.0/2048.0*fine_values[17]-fine_values[7]/4096.0+3.0/8192.0*fine_values
[8]+fine_values[6]/16384.0+93.0/2048.0*fine_values[44]+3.0/8192.0*fine_values
[41]+961.0/1024.0*fine_values[38]+31.0/4096.0*fine_values[35]+31.0/4096.0*
fine_values[32]+fine_values[29]/16384.0+9.0/4096.0*fine_values[26]+3.0/8192.0*
fine_values[24]-3.0/2048.0*fine_values[25]+31.0/4096.0*fine_values[98]+
fine_values[97]/16384.0-3.0/2048.0*fine_values[92]+fine_values[89]/1024.0-
fine_values[86]/4096.0+3.0/8192.0*fine_values[83]-fine_values[80]/4096.0+
fine_values[77]/16384.0-fine_values[103]/4096.0-31.0/1024.0*fine_values[104];
      coarse_values[39] = 3.0/64.0*fine_values[41]+31.0/32.0*fine_values[35]+
fine_values[29]/128.0+fine_values[97]/128.0-fine_values[103]/32.0;
      coarse_values[40] = -3.0/2048.0*fine_values[21]+3.0/8192.0*fine_values
[18]+93.0/2048.0*fine_values[15]-31.0/1024.0*fine_values[12]+31.0/4096.0*
fine_values[9]+3.0/8192.0*fine_values[6]-fine_values[3]/4096.0+9.0/4096.0*
fine_values[24]+fine_values[0]/16384.0+3.0/8192.0*fine_values[71]+93.0/2048.0*
fine_values[72]+961.0/1024.0*fine_values[68]+31.0/4096.0*fine_values[67]+31.0/
4096.0*fine_values[64]+fine_values[63]/16384.0-3.0/2048.0*fine_values[86]+
fine_values[85]/1024.0-fine_values[84]/4096.0+fine_values[75]/16384.0-
fine_values[76]/4096.0+3.0/8192.0*fine_values[77]-31.0/1024.0*fine_values[123]-
fine_values[121]/4096.0+fine_values[117]/16384.0+31.0/4096.0*fine_values[119];
      s2 = -3.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]+93.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]+93.0/262144.0*
fine_values[15]+2883.0/65536.0*fine_values[16]+279.0/131072.0*fine_values[17]
-961.0/32768.0*fine_values[13]-93.0/65536.0*fine_values[14]+93.0/262144.0*
fine_values[11]-31.0/131072.0*fine_values[12]+31.0/524288.0*fine_values[9]+
961.0/131072.0*fine_values[10]+93.0/262144.0*fine_values[7]+9.0/524288.0*
fine_values[8]-3.0/262144.0*fine_values[5]+3.0/1048576.0*fine_values[6]-31.0/
131072.0*fine_values[4]-fine_values[3]/524288.0+3.0/1048576.0*fine_values[2]+
31.0/524288.0*fine_values[1]-31.0/131072.0*fine_values[49]+3.0/1048576.0*
fine_values[47]+31.0/524288.0*fine_values[48]+fine_values[45]/2097152.0-
fine_values[46]/524288.0+3.0/65536.0*fine_values[43]-9.0/131072.0*fine_values
[44]+9.0/524288.0*fine_values[41]-3.0/262144.0*fine_values[42]-3.0/262144.0*
fine_values[40];
      s1 = s2+3.0/1048576.0*fine_values[39]-93.0/65536.0*fine_values[38]+31.0/
32768.0*fine_values[37]-31.0/131072.0*fine_values[34]-31.0/131072.0*fine_values
[36]+93.0/262144.0*fine_values[35]+fine_values[31]/131072.0-3.0/262144.0*
fine_values[32]+31.0/524288.0*fine_values[33]-fine_values[30]/524288.0-
fine_values[28]/524288.0+3.0/1048576.0*fine_values[29]+27.0/262144.0*
fine_values[26]+fine_values[27]/2097152.0+9.0/524288.0*fine_values[24]+279.0/
131072.0*fine_values[25]+fine_values[0]/2097152.0-9.0/131072.0*fine_values[23]
-93.0/65536.0*fine_values[22]+9.0/524288.0*fine_values[59]-3.0/262144.0*
fine_values[58]+3.0/1048576.0*fine_values[57]+2883.0/65536.0*fine_values[56]+
93.0/262144.0*fine_values[53]-961.0/32768.0*fine_values[55]+961.0/131072.0*
fine_values[54]-31.0/131072.0*fine_values[52]+93.0/262144.0*fine_values[50]+
31.0/524288.0*fine_values[51]+93.0/262144.0*fine_values[73]+2883.0/65536.0*
fine_values[74];
      s2 = s1+3.0/1048576.0*fine_values[71]+93.0/262144.0*fine_values[72]+
29791.0/32768.0*fine_values[70]+961.0/131072.0*fine_values[69]+961.0/131072.0*
fine_values[68]+31.0/524288.0*fine_values[65]+961.0/131072.0*fine_values[66]+
31.0/524288.0*fine_values[67]+31.0/524288.0*fine_values[64]-93.0/65536.0*
fine_values[61]+279.0/131072.0*fine_values[62]+fine_values[63]/2097152.0+
fine_values[101]/131072.0-fine_values[102]/32768.0-fine_values[99]/524288.0+
fine_values[100]/131072.0-3.0/262144.0*fine_values[98]+fine_values[96]/131072.0
+3.0/1048576.0*fine_values[97]-fine_values[94]/524288.0-fine_values[95]/
524288.0+3.0/65536.0*fine_values[91]+fine_values[93]/2097152.0-9.0/131072.0*
fine_values[92]-93.0/65536.0*fine_values[89]-3.0/262144.0*fine_values[90]+31.0/
32768.0*fine_values[88]-31.0/131072.0*fine_values[87]-3.0/262144.0*fine_values
[86]+fine_values[85]/131072.0-3.0/262144.0*fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]-fine_values[84]/524288.0+3.0/
1048576.0*fine_values[81]+93.0/262144.0*fine_values[80]-31.0/131072.0*
fine_values[79]+31.0/524288.0*fine_values[78]+fine_values[75]/2097152.0-
fine_values[76]/524288.0+3.0/1048576.0*fine_values[77]+93.0/262144.0*
fine_values[60]-31.0/131072.0*fine_values[112]-fine_values[111]/524288.0+3.0/
1048576.0*fine_values[109]+93.0/262144.0*fine_values[110]+fine_values[105]/
2097152.0;
      coarse_values[41] = s3-fine_values[107]/524288.0-31.0/131072.0*
fine_values[108]+31.0/524288.0*fine_values[106]-3.0/262144.0*fine_values[103]+
3.0/65536.0*fine_values[104]-961.0/32768.0*fine_values[124]-31.0/131072.0*
fine_values[122]-31.0/131072.0*fine_values[123]+961.0/131072.0*fine_values[120]
-fine_values[121]/524288.0+fine_values[117]/2097152.0-93.0/65536.0*fine_values
[116]+31.0/524288.0*fine_values[119]+31.0/524288.0*fine_values[118]-3.0/
262144.0*fine_values[115]+fine_values[113]/131072.0+31.0/32768.0*fine_values
[114];
      coarse_values[42] = 3.0/8192.0*fine_values[20]+93.0/2048.0*fine_values
[17]-31.0/1024.0*fine_values[14]+31.0/4096.0*fine_values[11]+3.0/8192.0*
fine_values[8]-fine_values[5]/4096.0+fine_values[2]/16384.0+fine_values[47]/
16384.0+9.0/4096.0*fine_values[26]-3.0/2048.0*fine_values[23]+3.0/8192.0*
fine_values[59]+961.0/1024.0*fine_values[56]+31.0/4096.0*fine_values[53]+31.0/
4096.0*fine_values[50]+93.0/2048.0*fine_values[62]+fine_values[91]/1024.0-3.0/
2048.0*fine_values[92]-fine_values[90]/4096.0-fine_values[82]/4096.0+3.0/8192.0
*fine_values[83]+fine_values[81]/16384.0+fine_values[109]/16384.0+31.0/4096.0*
fine_values[110]-31.0/1024.0*fine_values[116]-fine_values[115]/4096.0;
      s2 = -3.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]-3.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]+93.0/262144.0*
fine_values[15]-93.0/65536.0*fine_values[16]+279.0/131072.0*fine_values[17]+
31.0/32768.0*fine_values[13]-93.0/65536.0*fine_values[14]+93.0/262144.0*
fine_values[11]-31.0/131072.0*fine_values[12]+31.0/524288.0*fine_values[9]-31.0
/131072.0*fine_values[10]-3.0/262144.0*fine_values[7]+9.0/524288.0*fine_values
[8]-3.0/262144.0*fine_values[5]+3.0/1048576.0*fine_values[6]+fine_values[4]/
131072.0-fine_values[3]/524288.0+3.0/1048576.0*fine_values[2]-fine_values[1]/
524288.0+961.0/131072.0*fine_values[49]+3.0/1048576.0*fine_values[47]+31.0/
524288.0*fine_values[48]+fine_values[45]/2097152.0+31.0/524288.0*fine_values
[46]-93.0/65536.0*fine_values[43]+279.0/131072.0*fine_values[44]+9.0/524288.0*
fine_values[41]+93.0/262144.0*fine_values[42]-3.0/262144.0*fine_values[40];
      s1 = s2+3.0/1048576.0*fine_values[39]+2883.0/65536.0*fine_values[38]
-961.0/32768.0*fine_values[37]-31.0/131072.0*fine_values[34]+961.0/131072.0*
fine_values[36]+93.0/262144.0*fine_values[35]-31.0/131072.0*fine_values[31]+
93.0/262144.0*fine_values[32]+31.0/524288.0*fine_values[33]+31.0/524288.0*
fine_values[30]-fine_values[28]/524288.0+3.0/1048576.0*fine_values[29]+27.0/
262144.0*fine_values[26]+fine_values[27]/2097152.0+9.0/524288.0*fine_values[24]
-9.0/131072.0*fine_values[25]+fine_values[0]/2097152.0-9.0/131072.0*fine_values
[23]+3.0/65536.0*fine_values[22]+9.0/524288.0*fine_values[59]+93.0/262144.0*
fine_values[58]+3.0/1048576.0*fine_values[57]+2883.0/65536.0*fine_values[56]+
93.0/262144.0*fine_values[53]+29791.0/32768.0*fine_values[55]+961.0/131072.0*
fine_values[54]+961.0/131072.0*fine_values[52]+93.0/262144.0*fine_values[50]+
31.0/524288.0*fine_values[51]-3.0/262144.0*fine_values[73]-93.0/65536.0*
fine_values[74];
      s2 = s1+3.0/1048576.0*fine_values[71]+93.0/262144.0*fine_values[72]-961.0
/32768.0*fine_values[70]-31.0/131072.0*fine_values[69]+961.0/131072.0*
fine_values[68]-fine_values[65]/524288.0-31.0/131072.0*fine_values[66]+31.0/
524288.0*fine_values[67]+31.0/524288.0*fine_values[64]+2883.0/65536.0*
fine_values[61]+279.0/131072.0*fine_values[62]+fine_values[63]/2097152.0+
fine_values[101]/131072.0+31.0/32768.0*fine_values[102]-fine_values[99]/
524288.0-31.0/131072.0*fine_values[100]+93.0/262144.0*fine_values[98]-31.0/
131072.0*fine_values[96]+3.0/1048576.0*fine_values[97]+31.0/524288.0*
fine_values[94]-fine_values[95]/524288.0+3.0/65536.0*fine_values[91]+
fine_values[93]/2097152.0-9.0/131072.0*fine_values[92]+3.0/65536.0*fine_values
[89]-3.0/262144.0*fine_values[90]-fine_values[88]/32768.0+fine_values[87]/
131072.0-3.0/262144.0*fine_values[86]+fine_values[85]/131072.0-3.0/262144.0*
fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]-fine_values[84]/524288.0+3.0/
1048576.0*fine_values[81]-3.0/262144.0*fine_values[80]+fine_values[79]/131072.0
-fine_values[78]/524288.0+fine_values[75]/2097152.0-fine_values[76]/524288.0+
3.0/1048576.0*fine_values[77]+93.0/262144.0*fine_values[60]-31.0/131072.0*
fine_values[112]-fine_values[111]/524288.0+3.0/1048576.0*fine_values[109]+93.0/
262144.0*fine_values[110]+fine_values[105]/2097152.0;
      coarse_values[43] = s3+31.0/524288.0*fine_values[107]+961.0/131072.0*
fine_values[108]+31.0/524288.0*fine_values[106]-3.0/262144.0*fine_values[103]
-93.0/65536.0*fine_values[104]+31.0/32768.0*fine_values[124]+fine_values[122]/
131072.0-31.0/131072.0*fine_values[123]-31.0/131072.0*fine_values[120]-
fine_values[121]/524288.0+fine_values[117]/2097152.0-93.0/65536.0*fine_values
[116]+31.0/524288.0*fine_values[119]-fine_values[118]/524288.0-3.0/262144.0*
fine_values[115]-31.0/131072.0*fine_values[113]-961.0/32768.0*fine_values[114];
      coarse_values[44] = 31.0/4096.0*fine_values[48]+fine_values[45]/16384.0+
9.0/4096.0*fine_values[41]-3.0/2048.0*fine_values[40]+3.0/8192.0*fine_values
[39]-31.0/1024.0*fine_values[34]+93.0/2048.0*fine_values[35]+31.0/4096.0*
fine_values[33]-fine_values[28]/4096.0+3.0/8192.0*fine_values[29]+fine_values
[27]/16384.0+3.0/8192.0*fine_values[57]+961.0/1024.0*fine_values[54]+31.0/
4096.0*fine_values[51]+fine_values[101]/1024.0-fine_values[99]/4096.0+3.0/
8192.0*fine_values[97]-fine_values[95]/4096.0+fine_values[93]/16384.0+93.0/
2048.0*fine_values[60]-31.0/1024.0*fine_values[112]-fine_values[111]/4096.0+
fine_values[105]/16384.0+31.0/4096.0*fine_values[106]-3.0/2048.0*fine_values
[103];
      coarse_values[45] = 3.0/64.0*fine_values[71]+31.0/32.0*fine_values[67]+
fine_values[63]/128.0-fine_values[121]/32.0+fine_values[117]/128.0;
      coarse_values[46] = 3.0/8192.0*fine_values[47]+fine_values[45]/16384.0-
fine_values[46]/4096.0+9.0/4096.0*fine_values[59]-3.0/2048.0*fine_values[58]+
3.0/8192.0*fine_values[57]+93.0/2048.0*fine_values[53]-31.0/1024.0*fine_values
[52]+31.0/4096.0*fine_values[51]+93.0/2048.0*fine_values[73]+3.0/8192.0*
fine_values[71]+961.0/1024.0*fine_values[69]+31.0/4096.0*fine_values[65]+31.0/
4096.0*fine_values[67]+fine_values[63]/16384.0-fine_values[111]/4096.0+3.0/
8192.0*fine_values[109]+fine_values[105]/16384.0-fine_values[107]/4096.0-31.0/
1024.0*fine_values[122]-fine_values[121]/4096.0+fine_values[117]/16384.0+31.0/
4096.0*fine_values[118]-3.0/2048.0*fine_values[115]+fine_values[113]/1024.0;
      coarse_values[47] = fine_values[47]/128.0+31.0/32.0*fine_values[53]+3.0/
64.0*fine_values[59]+fine_values[109]/128.0-fine_values[115]/32.0;
      coarse_values[48] = 3.0/8192.0*fine_values[47]+fine_values[45]/16384.0+
31.0/4096.0*fine_values[46]+9.0/4096.0*fine_values[59]+93.0/2048.0*fine_values
[58]+3.0/8192.0*fine_values[57]+93.0/2048.0*fine_values[53]+961.0/1024.0*
fine_values[52]+31.0/4096.0*fine_values[51]-3.0/2048.0*fine_values[73]+3.0/
8192.0*fine_values[71]-31.0/1024.0*fine_values[69]-fine_values[65]/4096.0+31.0/
4096.0*fine_values[67]+fine_values[63]/16384.0-fine_values[111]/4096.0+3.0/
8192.0*fine_values[109]+fine_values[105]/16384.0+31.0/4096.0*fine_values[107]+
fine_values[122]/1024.0-fine_values[121]/4096.0+fine_values[117]/16384.0-
fine_values[118]/4096.0-3.0/2048.0*fine_values[115]-31.0/1024.0*fine_values
[113];
      coarse_values[49] = fine_values[45]/128.0+3.0/64.0*fine_values[57]+31.0/
32.0*fine_values[51]-fine_values[111]/32.0+fine_values[105]/128.0;
      coarse_values[50] = fine_values[18];
      coarse_values[51] = 3.0/64.0*fine_values[20]+31.0/32.0*fine_values[19]+
fine_values[18]/128.0-fine_values[42]/32.0+fine_values[39]/128.0;
      coarse_values[52] = fine_values[20];
      coarse_values[53] = 3.0/64.0*fine_values[20]-fine_values[19]/32.0+
fine_values[18]/128.0+31.0/32.0*fine_values[42]+fine_values[39]/128.0;
      coarse_values[54] = fine_values[39];
      coarse_values[55] = fine_values[18]/128.0+3.0/64.0*fine_values[24]+31.0/
32.0*fine_values[21]-fine_values[72]/32.0+fine_values[71]/128.0;
      coarse_values[56] = 31.0/4096.0*fine_values[21]+3.0/8192.0*fine_values
[20]+31.0/4096.0*fine_values[19]+fine_values[18]/16384.0-31.0/1024.0*
fine_values[43]-3.0/2048.0*fine_values[44]+3.0/8192.0*fine_values[41]-
fine_values[42]/4096.0+31.0/4096.0*fine_values[40]+fine_values[39]/16384.0+9.0/
4096.0*fine_values[26]+3.0/8192.0*fine_values[24]+93.0/2048.0*fine_values[25]+
93.0/2048.0*fine_values[23]+961.0/1024.0*fine_values[22]+3.0/8192.0*fine_values
[59]-fine_values[58]/4096.0+fine_values[57]/16384.0+31.0/4096.0*fine_values[73]
-31.0/1024.0*fine_values[74]+fine_values[71]/16384.0-fine_values[72]/4096.0+
fine_values[61]/1024.0-3.0/2048.0*fine_values[62]-fine_values[60]/4096.0;
      coarse_values[57] = fine_values[20]/128.0+3.0/64.0*fine_values[26]+31.0/
32.0*fine_values[23]-fine_values[62]/32.0+fine_values[59]/128.0;
      coarse_values[58] = 31.0/4096.0*fine_values[21]+3.0/8192.0*fine_values
[20]-fine_values[19]/4096.0+fine_values[18]/16384.0+961.0/1024.0*fine_values
[43]+93.0/2048.0*fine_values[44]+3.0/8192.0*fine_values[41]+31.0/4096.0*
fine_values[42]+31.0/4096.0*fine_values[40]+fine_values[39]/16384.0+9.0/4096.0*
fine_values[26]+3.0/8192.0*fine_values[24]-3.0/2048.0*fine_values[25]+93.0/
2048.0*fine_values[23]-31.0/1024.0*fine_values[22]+3.0/8192.0*fine_values[59]+
31.0/4096.0*fine_values[58]+fine_values[57]/16384.0-fine_values[73]/4096.0+
fine_values[74]/1024.0+fine_values[71]/16384.0-fine_values[72]/4096.0-31.0/
1024.0*fine_values[61]-3.0/2048.0*fine_values[62]-fine_values[60]/4096.0;
      coarse_values[59] = 31.0/32.0*fine_values[40]+3.0/64.0*fine_values[41]+
fine_values[39]/128.0+fine_values[57]/128.0-fine_values[60]/32.0;
      coarse_values[60] = fine_values[24];
      coarse_values[61] = -fine_values[44]/32.0+fine_values[41]/128.0+31.0/32.0
*fine_values[25]+3.0/64.0*fine_values[26]+fine_values[24]/128.0;
      coarse_values[62] = fine_values[26];
      coarse_values[63] = 31.0/32.0*fine_values[44]+fine_values[41]/128.0-
fine_values[25]/32.0+3.0/64.0*fine_values[26]+fine_values[24]/128.0;
      coarse_values[64] = fine_values[41];
      coarse_values[65] = fine_values[18]/128.0+3.0/64.0*fine_values[24]-
fine_values[21]/32.0+31.0/32.0*fine_values[72]+fine_values[71]/128.0;
      coarse_values[66] = -fine_values[21]/4096.0+3.0/8192.0*fine_values[20]+
31.0/4096.0*fine_values[19]+fine_values[18]/16384.0+fine_values[43]/1024.0-3.0/
2048.0*fine_values[44]+3.0/8192.0*fine_values[41]-fine_values[42]/4096.0-
fine_values[40]/4096.0+fine_values[39]/16384.0+9.0/4096.0*fine_values[26]+3.0/
8192.0*fine_values[24]+93.0/2048.0*fine_values[25]-3.0/2048.0*fine_values[23]
-31.0/1024.0*fine_values[22]+3.0/8192.0*fine_values[59]-fine_values[58]/4096.0+
fine_values[57]/16384.0+31.0/4096.0*fine_values[73]+961.0/1024.0*fine_values
[74]+fine_values[71]/16384.0+31.0/4096.0*fine_values[72]-31.0/1024.0*
fine_values[61]+93.0/2048.0*fine_values[62]+31.0/4096.0*fine_values[60];
      coarse_values[67] = fine_values[20]/128.0+3.0/64.0*fine_values[26]-
fine_values[23]/32.0+31.0/32.0*fine_values[62]+fine_values[59]/128.0;
      coarse_values[68] = -fine_values[21]/4096.0+3.0/8192.0*fine_values[20]-
fine_values[19]/4096.0+fine_values[18]/16384.0-31.0/1024.0*fine_values[43]+93.0
/2048.0*fine_values[44]+3.0/8192.0*fine_values[41]+31.0/4096.0*fine_values[42]-
fine_values[40]/4096.0+fine_values[39]/16384.0+9.0/4096.0*fine_values[26]+3.0/
8192.0*fine_values[24]-3.0/2048.0*fine_values[25]-3.0/2048.0*fine_values[23]+
fine_values[22]/1024.0+3.0/8192.0*fine_values[59]+31.0/4096.0*fine_values[58]+
fine_values[57]/16384.0-fine_values[73]/4096.0-31.0/1024.0*fine_values[74]+
fine_values[71]/16384.0+31.0/4096.0*fine_values[72]+961.0/1024.0*fine_values
[61]+93.0/2048.0*fine_values[62]+31.0/4096.0*fine_values[60];
      coarse_values[69] = -fine_values[40]/32.0+3.0/64.0*fine_values[41]+
fine_values[39]/128.0+fine_values[57]/128.0+31.0/32.0*fine_values[60];
      coarse_values[70] = fine_values[71];
      coarse_values[71] = -fine_values[58]/32.0+fine_values[57]/128.0+31.0/32.0
*fine_values[73]+fine_values[71]/128.0+3.0/64.0*fine_values[59];
      coarse_values[72] = fine_values[59];
      coarse_values[73] = 31.0/32.0*fine_values[58]+fine_values[57]/128.0-
fine_values[73]/32.0+fine_values[71]/128.0+3.0/64.0*fine_values[59];
      coarse_values[74] = fine_values[57];
      coarse_values[75] = 3.0/64.0*fine_values[18]-fine_values[9]/32.0+
fine_values[0]/128.0+31.0/32.0*fine_values[84]+fine_values[75]/128.0;
      coarse_values[76] = 9.0/4096.0*fine_values[20]+93.0/2048.0*fine_values
[19]+3.0/8192.0*fine_values[18]-3.0/2048.0*fine_values[11]-fine_values[9]/
4096.0-31.0/1024.0*fine_values[10]+3.0/8192.0*fine_values[2]+31.0/4096.0*
fine_values[1]-3.0/2048.0*fine_values[42]+3.0/8192.0*fine_values[39]+
fine_values[36]/1024.0-fine_values[33]/4096.0-fine_values[30]/4096.0+
fine_values[27]/16384.0+fine_values[0]/16384.0+31.0/4096.0*fine_values[99]-31.0
/1024.0*fine_values[100]-fine_values[94]/4096.0+fine_values[93]/16384.0+93.0/
2048.0*fine_values[90]+961.0/1024.0*fine_values[87]+31.0/4096.0*fine_values[84]
+3.0/8192.0*fine_values[81]+31.0/4096.0*fine_values[78]+fine_values[75]/16384.0
;
      coarse_values[77] = 3.0/64.0*fine_values[20]-fine_values[11]/32.0+
fine_values[2]/128.0+31.0/32.0*fine_values[90]+fine_values[81]/128.0;
      coarse_values[78] = 9.0/4096.0*fine_values[20]-3.0/2048.0*fine_values[19]
+3.0/8192.0*fine_values[18]-3.0/2048.0*fine_values[11]-fine_values[9]/4096.0+
fine_values[10]/1024.0+3.0/8192.0*fine_values[2]-fine_values[1]/4096.0+93.0/
2048.0*fine_values[42]+3.0/8192.0*fine_values[39]-31.0/1024.0*fine_values[36]-
fine_values[33]/4096.0+31.0/4096.0*fine_values[30]+fine_values[27]/16384.0+
fine_values[0]/16384.0+31.0/4096.0*fine_values[99]+961.0/1024.0*fine_values
[100]+31.0/4096.0*fine_values[94]+fine_values[93]/16384.0+93.0/2048.0*
fine_values[90]-31.0/1024.0*fine_values[87]+31.0/4096.0*fine_values[84]+3.0/
8192.0*fine_values[81]-fine_values[78]/4096.0+fine_values[75]/16384.0;
      coarse_values[79] = 3.0/64.0*fine_values[39]-fine_values[33]/32.0+
fine_values[27]/128.0+31.0/32.0*fine_values[99]+fine_values[93]/128.0;
      coarse_values[80] = 93.0/2048.0*fine_values[21]+3.0/8192.0*fine_values
[18]-3.0/2048.0*fine_values[15]-31.0/1024.0*fine_values[12]-fine_values[9]/
4096.0+3.0/8192.0*fine_values[6]+31.0/4096.0*fine_values[3]+9.0/4096.0*
fine_values[24]+fine_values[0]/16384.0+3.0/8192.0*fine_values[71]-3.0/2048.0*
fine_values[72]+fine_values[68]/1024.0-fine_values[67]/4096.0-fine_values[64]/
4096.0+fine_values[63]/16384.0+93.0/2048.0*fine_values[86]+961.0/1024.0*
fine_values[85]+31.0/4096.0*fine_values[84]+fine_values[75]/16384.0+31.0/4096.0
*fine_values[76]+3.0/8192.0*fine_values[77]-31.0/1024.0*fine_values[123]+31.0/
4096.0*fine_values[121]+fine_values[117]/16384.0-fine_values[119]/4096.0;
      s2 = 93.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]+93.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]-3.0/262144.0*fine_values
[15]-93.0/65536.0*fine_values[16]-9.0/131072.0*fine_values[17]-961.0/32768.0*
fine_values[13]-93.0/65536.0*fine_values[14]-3.0/262144.0*fine_values[11]-31.0/
131072.0*fine_values[12]-fine_values[9]/524288.0-31.0/131072.0*fine_values[10]+
93.0/262144.0*fine_values[7]+9.0/524288.0*fine_values[8]+93.0/262144.0*
fine_values[5]+3.0/1048576.0*fine_values[6]+961.0/131072.0*fine_values[4]+31.0/
524288.0*fine_values[3]+3.0/1048576.0*fine_values[2]+31.0/524288.0*fine_values
[1]+fine_values[49]/131072.0+3.0/1048576.0*fine_values[47]-fine_values[48]/
524288.0+fine_values[45]/2097152.0-fine_values[46]/524288.0-93.0/65536.0*
fine_values[43]-9.0/131072.0*fine_values[44]+9.0/524288.0*fine_values[41]-3.0/
262144.0*fine_values[42]+93.0/262144.0*fine_values[40];
      s1 = s2+3.0/1048576.0*fine_values[39]+3.0/65536.0*fine_values[38]+31.0/
32768.0*fine_values[37]-31.0/131072.0*fine_values[34]+fine_values[36]/131072.0
-3.0/262144.0*fine_values[35]-31.0/131072.0*fine_values[31]-3.0/262144.0*
fine_values[32]-fine_values[33]/524288.0-fine_values[30]/524288.0+31.0/524288.0
*fine_values[28]+3.0/1048576.0*fine_values[29]+27.0/262144.0*fine_values[26]+
fine_values[27]/2097152.0+9.0/524288.0*fine_values[24]+279.0/131072.0*
fine_values[25]+fine_values[0]/2097152.0+279.0/131072.0*fine_values[23]+2883.0/
65536.0*fine_values[22]+9.0/524288.0*fine_values[59]-3.0/262144.0*fine_values
[58]+3.0/1048576.0*fine_values[57]+3.0/65536.0*fine_values[56]-3.0/262144.0*
fine_values[53]-fine_values[55]/32768.0+fine_values[54]/131072.0+fine_values
[52]/131072.0-3.0/262144.0*fine_values[50]-fine_values[51]/524288.0+93.0/
262144.0*fine_values[73]-93.0/65536.0*fine_values[74];
      s2 = s1+3.0/1048576.0*fine_values[71]-3.0/262144.0*fine_values[72]+31.0/
32768.0*fine_values[70]-31.0/131072.0*fine_values[69]+fine_values[68]/131072.0+
31.0/524288.0*fine_values[65]-31.0/131072.0*fine_values[66]-fine_values[67]/
524288.0-fine_values[64]/524288.0+3.0/65536.0*fine_values[61]-9.0/131072.0*
fine_values[62]+fine_values[63]/2097152.0+961.0/131072.0*fine_values[101]-961.0
/32768.0*fine_values[102]+31.0/524288.0*fine_values[99]-31.0/131072.0*
fine_values[100]-3.0/262144.0*fine_values[98]-31.0/131072.0*fine_values[96]+3.0
/1048576.0*fine_values[97]-fine_values[94]/524288.0+31.0/524288.0*fine_values
[95]+2883.0/65536.0*fine_values[91]+fine_values[93]/2097152.0+279.0/131072.0*
fine_values[92]+2883.0/65536.0*fine_values[89]+93.0/262144.0*fine_values[90]+
29791.0/32768.0*fine_values[88]+961.0/131072.0*fine_values[87]+93.0/262144.0*
fine_values[86]+961.0/131072.0*fine_values[85]+93.0/262144.0*fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]+31.0/524288.0*fine_values[84]+3.0/
1048576.0*fine_values[81]+93.0/262144.0*fine_values[80]+961.0/131072.0*
fine_values[79]+31.0/524288.0*fine_values[78]+fine_values[75]/2097152.0+31.0/
524288.0*fine_values[76]+3.0/1048576.0*fine_values[77]-3.0/262144.0*fine_values
[60]-31.0/131072.0*fine_values[112]+31.0/524288.0*fine_values[111]+3.0/
1048576.0*fine_values[109]-3.0/262144.0*fine_values[110]+fine_values[105]/
2097152.0;
      coarse_values[81] = s3-fine_values[107]/524288.0+fine_values[108]/
131072.0-fine_values[106]/524288.0+93.0/262144.0*fine_values[103]-93.0/65536.0*
fine_values[104]-961.0/32768.0*fine_values[124]+961.0/131072.0*fine_values[122]
-31.0/131072.0*fine_values[123]-31.0/131072.0*fine_values[120]+31.0/524288.0*
fine_values[121]+fine_values[117]/2097152.0-93.0/65536.0*fine_values[116]-
fine_values[119]/524288.0+31.0/524288.0*fine_values[118]+93.0/262144.0*
fine_values[115]-31.0/131072.0*fine_values[113]+31.0/32768.0*fine_values[114];
      coarse_values[82] = 3.0/8192.0*fine_values[20]-3.0/2048.0*fine_values[17]
-31.0/1024.0*fine_values[14]-fine_values[11]/4096.0+3.0/8192.0*fine_values[8]+
31.0/4096.0*fine_values[5]+fine_values[2]/16384.0+fine_values[47]/16384.0+9.0/
4096.0*fine_values[26]+93.0/2048.0*fine_values[23]+3.0/8192.0*fine_values[59]+
fine_values[56]/1024.0-fine_values[53]/4096.0-fine_values[50]/4096.0-3.0/2048.0
*fine_values[62]+961.0/1024.0*fine_values[91]+93.0/2048.0*fine_values[92]+31.0/
4096.0*fine_values[90]+31.0/4096.0*fine_values[82]+3.0/8192.0*fine_values[83]+
fine_values[81]/16384.0+fine_values[109]/16384.0-fine_values[110]/4096.0-31.0/
1024.0*fine_values[116]+31.0/4096.0*fine_values[115];
      s2 = 93.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]-3.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]-3.0/262144.0*fine_values
[15]+3.0/65536.0*fine_values[16]-9.0/131072.0*fine_values[17]+31.0/32768.0*
fine_values[13]-93.0/65536.0*fine_values[14]-3.0/262144.0*fine_values[11]-31.0/
131072.0*fine_values[12]-fine_values[9]/524288.0+fine_values[10]/131072.0-3.0/
262144.0*fine_values[7]+9.0/524288.0*fine_values[8]+93.0/262144.0*fine_values
[5]+3.0/1048576.0*fine_values[6]-31.0/131072.0*fine_values[4]+31.0/524288.0*
fine_values[3]+3.0/1048576.0*fine_values[2]-fine_values[1]/524288.0-31.0/
131072.0*fine_values[49]+3.0/1048576.0*fine_values[47]-fine_values[48]/524288.0
+fine_values[45]/2097152.0+31.0/524288.0*fine_values[46]+2883.0/65536.0*
fine_values[43]+279.0/131072.0*fine_values[44]+9.0/524288.0*fine_values[41]+
93.0/262144.0*fine_values[42]+93.0/262144.0*fine_values[40];
      s1 = s2+3.0/1048576.0*fine_values[39]-93.0/65536.0*fine_values[38]-961.0/
32768.0*fine_values[37]-31.0/131072.0*fine_values[34]-31.0/131072.0*fine_values
[36]-3.0/262144.0*fine_values[35]+961.0/131072.0*fine_values[31]+93.0/262144.0*
fine_values[32]-fine_values[33]/524288.0+31.0/524288.0*fine_values[30]+31.0/
524288.0*fine_values[28]+3.0/1048576.0*fine_values[29]+27.0/262144.0*
fine_values[26]+fine_values[27]/2097152.0+9.0/524288.0*fine_values[24]-9.0/
131072.0*fine_values[25]+fine_values[0]/2097152.0+279.0/131072.0*fine_values
[23]-93.0/65536.0*fine_values[22]+9.0/524288.0*fine_values[59]+93.0/262144.0*
fine_values[58]+3.0/1048576.0*fine_values[57]+3.0/65536.0*fine_values[56]-3.0/
262144.0*fine_values[53]+31.0/32768.0*fine_values[55]+fine_values[54]/131072.0
-31.0/131072.0*fine_values[52]-3.0/262144.0*fine_values[50]-fine_values[51]/
524288.0-3.0/262144.0*fine_values[73]+3.0/65536.0*fine_values[74];
      s2 = s1+3.0/1048576.0*fine_values[71]-3.0/262144.0*fine_values[72]-
fine_values[70]/32768.0+fine_values[69]/131072.0+fine_values[68]/131072.0-
fine_values[65]/524288.0+fine_values[66]/131072.0-fine_values[67]/524288.0-
fine_values[64]/524288.0-93.0/65536.0*fine_values[61]-9.0/131072.0*fine_values
[62]+fine_values[63]/2097152.0+961.0/131072.0*fine_values[101]+29791.0/32768.0*
fine_values[102]+31.0/524288.0*fine_values[99]+961.0/131072.0*fine_values[100]+
93.0/262144.0*fine_values[98]+961.0/131072.0*fine_values[96]+3.0/1048576.0*
fine_values[97]+31.0/524288.0*fine_values[94]+31.0/524288.0*fine_values[95]+
2883.0/65536.0*fine_values[91]+fine_values[93]/2097152.0+279.0/131072.0*
fine_values[92]-93.0/65536.0*fine_values[89]+93.0/262144.0*fine_values[90]
-961.0/32768.0*fine_values[88]-31.0/131072.0*fine_values[87]+93.0/262144.0*
fine_values[86]+961.0/131072.0*fine_values[85]+93.0/262144.0*fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]+31.0/524288.0*fine_values[84]+3.0/
1048576.0*fine_values[81]-3.0/262144.0*fine_values[80]-31.0/131072.0*
fine_values[79]-fine_values[78]/524288.0+fine_values[75]/2097152.0+31.0/
524288.0*fine_values[76]+3.0/1048576.0*fine_values[77]-3.0/262144.0*fine_values
[60]-31.0/131072.0*fine_values[112]+31.0/524288.0*fine_values[111]+3.0/
1048576.0*fine_values[109]-3.0/262144.0*fine_values[110]+fine_values[105]/
2097152.0;
      coarse_values[83] = s3+31.0/524288.0*fine_values[107]-31.0/131072.0*
fine_values[108]-fine_values[106]/524288.0+93.0/262144.0*fine_values[103]+
2883.0/65536.0*fine_values[104]+31.0/32768.0*fine_values[124]-31.0/131072.0*
fine_values[122]-31.0/131072.0*fine_values[123]+fine_values[120]/131072.0+31.0/
524288.0*fine_values[121]+fine_values[117]/2097152.0-93.0/65536.0*fine_values
[116]-fine_values[119]/524288.0-fine_values[118]/524288.0+93.0/262144.0*
fine_values[115]+961.0/131072.0*fine_values[113]-961.0/32768.0*fine_values[114]
;
      coarse_values[84] = -fine_values[48]/4096.0+fine_values[45]/16384.0+9.0/
4096.0*fine_values[41]+93.0/2048.0*fine_values[40]+3.0/8192.0*fine_values[39]
-31.0/1024.0*fine_values[34]-3.0/2048.0*fine_values[35]-fine_values[33]/4096.0+
31.0/4096.0*fine_values[28]+3.0/8192.0*fine_values[29]+fine_values[27]/16384.0+
3.0/8192.0*fine_values[57]+fine_values[54]/1024.0-fine_values[51]/4096.0+961.0/
1024.0*fine_values[101]+31.0/4096.0*fine_values[99]+3.0/8192.0*fine_values[97]+
31.0/4096.0*fine_values[95]+fine_values[93]/16384.0-3.0/2048.0*fine_values[60]
-31.0/1024.0*fine_values[112]+31.0/4096.0*fine_values[111]+fine_values[105]/
16384.0-fine_values[106]/4096.0+93.0/2048.0*fine_values[103];
      coarse_values[85] = -fine_values[15]/32.0+fine_values[6]/128.0+3.0/64.0*
fine_values[24]+31.0/32.0*fine_values[86]+fine_values[77]/128.0;
      coarse_values[86] = -fine_values[15]/4096.0-31.0/1024.0*fine_values[16]
-3.0/2048.0*fine_values[17]+31.0/4096.0*fine_values[7]+3.0/8192.0*fine_values
[8]+fine_values[6]/16384.0-3.0/2048.0*fine_values[44]+3.0/8192.0*fine_values
[41]+fine_values[38]/1024.0-fine_values[35]/4096.0-fine_values[32]/4096.0+
fine_values[29]/16384.0+9.0/4096.0*fine_values[26]+3.0/8192.0*fine_values[24]+
93.0/2048.0*fine_values[25]-fine_values[98]/4096.0+fine_values[97]/16384.0+93.0
/2048.0*fine_values[92]+961.0/1024.0*fine_values[89]+31.0/4096.0*fine_values
[86]+3.0/8192.0*fine_values[83]+31.0/4096.0*fine_values[80]+fine_values[77]/
16384.0+31.0/4096.0*fine_values[103]-31.0/1024.0*fine_values[104];
      coarse_values[87] = -fine_values[17]/32.0+fine_values[8]/128.0+3.0/64.0*
fine_values[26]+31.0/32.0*fine_values[92]+fine_values[83]/128.0;
      coarse_values[88] = -fine_values[15]/4096.0+fine_values[16]/1024.0-3.0/
2048.0*fine_values[17]-fine_values[7]/4096.0+3.0/8192.0*fine_values[8]+
fine_values[6]/16384.0+93.0/2048.0*fine_values[44]+3.0/8192.0*fine_values[41]
-31.0/1024.0*fine_values[38]-fine_values[35]/4096.0+31.0/4096.0*fine_values[32]
+fine_values[29]/16384.0+9.0/4096.0*fine_values[26]+3.0/8192.0*fine_values[24]
-3.0/2048.0*fine_values[25]+31.0/4096.0*fine_values[98]+fine_values[97]/16384.0
+93.0/2048.0*fine_values[92]-31.0/1024.0*fine_values[89]+31.0/4096.0*
fine_values[86]+3.0/8192.0*fine_values[83]-fine_values[80]/4096.0+fine_values
[77]/16384.0+31.0/4096.0*fine_values[103]+961.0/1024.0*fine_values[104];
      coarse_values[89] = 3.0/64.0*fine_values[41]-fine_values[35]/32.0+
fine_values[29]/128.0+fine_values[97]/128.0+31.0/32.0*fine_values[103];
      coarse_values[90] = -3.0/2048.0*fine_values[21]+3.0/8192.0*fine_values
[18]-3.0/2048.0*fine_values[15]+fine_values[12]/1024.0-fine_values[9]/4096.0+
3.0/8192.0*fine_values[6]-fine_values[3]/4096.0+9.0/4096.0*fine_values[24]+
fine_values[0]/16384.0+3.0/8192.0*fine_values[71]+93.0/2048.0*fine_values[72]
-31.0/1024.0*fine_values[68]-fine_values[67]/4096.0+31.0/4096.0*fine_values[64]
+fine_values[63]/16384.0+93.0/2048.0*fine_values[86]-31.0/1024.0*fine_values
[85]+31.0/4096.0*fine_values[84]+fine_values[75]/16384.0-fine_values[76]/4096.0
+3.0/8192.0*fine_values[77]+961.0/1024.0*fine_values[123]+31.0/4096.0*
fine_values[121]+fine_values[117]/16384.0+31.0/4096.0*fine_values[119];
      s2 = -3.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]+93.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]-3.0/262144.0*fine_values
[15]-93.0/65536.0*fine_values[16]-9.0/131072.0*fine_values[17]+31.0/32768.0*
fine_values[13]+3.0/65536.0*fine_values[14]-3.0/262144.0*fine_values[11]+
fine_values[12]/131072.0-fine_values[9]/524288.0-31.0/131072.0*fine_values[10]+
93.0/262144.0*fine_values[7]+9.0/524288.0*fine_values[8]-3.0/262144.0*
fine_values[5]+3.0/1048576.0*fine_values[6]-31.0/131072.0*fine_values[4]-
fine_values[3]/524288.0+3.0/1048576.0*fine_values[2]+31.0/524288.0*fine_values
[1]-31.0/131072.0*fine_values[49]+3.0/1048576.0*fine_values[47]+31.0/524288.0*
fine_values[48]+fine_values[45]/2097152.0-fine_values[46]/524288.0+3.0/65536.0*
fine_values[43]-9.0/131072.0*fine_values[44]+9.0/524288.0*fine_values[41]-3.0/
262144.0*fine_values[42]-3.0/262144.0*fine_values[40];
      s1 = s2+3.0/1048576.0*fine_values[39]+3.0/65536.0*fine_values[38]-
fine_values[37]/32768.0+fine_values[34]/131072.0+fine_values[36]/131072.0-3.0/
262144.0*fine_values[35]+fine_values[31]/131072.0-3.0/262144.0*fine_values[32]-
fine_values[33]/524288.0-fine_values[30]/524288.0-fine_values[28]/524288.0+3.0/
1048576.0*fine_values[29]+27.0/262144.0*fine_values[26]+fine_values[27]/
2097152.0+9.0/524288.0*fine_values[24]+279.0/131072.0*fine_values[25]+
fine_values[0]/2097152.0-9.0/131072.0*fine_values[23]-93.0/65536.0*fine_values
[22]+9.0/524288.0*fine_values[59]-3.0/262144.0*fine_values[58]+3.0/1048576.0*
fine_values[57]-93.0/65536.0*fine_values[56]-3.0/262144.0*fine_values[53]+31.0/
32768.0*fine_values[55]-31.0/131072.0*fine_values[54]+fine_values[52]/131072.0+
93.0/262144.0*fine_values[50]-fine_values[51]/524288.0+93.0/262144.0*
fine_values[73]+2883.0/65536.0*fine_values[74];
      s2 = s1+3.0/1048576.0*fine_values[71]+93.0/262144.0*fine_values[72]-961.0
/32768.0*fine_values[70]-31.0/131072.0*fine_values[69]-31.0/131072.0*
fine_values[68]+31.0/524288.0*fine_values[65]+961.0/131072.0*fine_values[66]-
fine_values[67]/524288.0+31.0/524288.0*fine_values[64]-93.0/65536.0*fine_values
[61]+279.0/131072.0*fine_values[62]+fine_values[63]/2097152.0-31.0/131072.0*
fine_values[101]+31.0/32768.0*fine_values[102]+31.0/524288.0*fine_values[99]
-31.0/131072.0*fine_values[100]-3.0/262144.0*fine_values[98]+fine_values[96]/
131072.0+3.0/1048576.0*fine_values[97]-fine_values[94]/524288.0-fine_values[95]
/524288.0-93.0/65536.0*fine_values[91]+fine_values[93]/2097152.0+279.0/131072.0
*fine_values[92]+2883.0/65536.0*fine_values[89]+93.0/262144.0*fine_values[90]
-961.0/32768.0*fine_values[88]+961.0/131072.0*fine_values[87]+93.0/262144.0*
fine_values[86]-31.0/131072.0*fine_values[85]-3.0/262144.0*fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]+31.0/524288.0*fine_values[84]+3.0/
1048576.0*fine_values[81]+93.0/262144.0*fine_values[80]-31.0/131072.0*
fine_values[79]+31.0/524288.0*fine_values[78]+fine_values[75]/2097152.0-
fine_values[76]/524288.0+3.0/1048576.0*fine_values[77]+93.0/262144.0*
fine_values[60]+961.0/131072.0*fine_values[112]+31.0/524288.0*fine_values[111]+
3.0/1048576.0*fine_values[109]+93.0/262144.0*fine_values[110]+fine_values[105]/
2097152.0;
      coarse_values[91] = s3-fine_values[107]/524288.0-31.0/131072.0*
fine_values[108]+31.0/524288.0*fine_values[106]+93.0/262144.0*fine_values[103]
-93.0/65536.0*fine_values[104]+29791.0/32768.0*fine_values[124]+961.0/131072.0*
fine_values[122]+961.0/131072.0*fine_values[123]+961.0/131072.0*fine_values
[120]+31.0/524288.0*fine_values[121]+fine_values[117]/2097152.0+2883.0/65536.0*
fine_values[116]+31.0/524288.0*fine_values[119]+31.0/524288.0*fine_values[118]+
93.0/262144.0*fine_values[115]-31.0/131072.0*fine_values[113]-961.0/32768.0*
fine_values[114];
      coarse_values[92] = 3.0/8192.0*fine_values[20]-3.0/2048.0*fine_values[17]
+fine_values[14]/1024.0-fine_values[11]/4096.0+3.0/8192.0*fine_values[8]-
fine_values[5]/4096.0+fine_values[2]/16384.0+fine_values[47]/16384.0+9.0/4096.0
*fine_values[26]-3.0/2048.0*fine_values[23]+3.0/8192.0*fine_values[59]-31.0/
1024.0*fine_values[56]-fine_values[53]/4096.0+31.0/4096.0*fine_values[50]+93.0/
2048.0*fine_values[62]-31.0/1024.0*fine_values[91]+93.0/2048.0*fine_values[92]+
31.0/4096.0*fine_values[90]-fine_values[82]/4096.0+3.0/8192.0*fine_values[83]+
fine_values[81]/16384.0+fine_values[109]/16384.0+31.0/4096.0*fine_values[110]+
961.0/1024.0*fine_values[116]+31.0/4096.0*fine_values[115];
      s2 = -3.0/262144.0*fine_values[21]+9.0/524288.0*fine_values[20]-3.0/
262144.0*fine_values[19]+3.0/1048576.0*fine_values[18]-3.0/262144.0*fine_values
[15]+3.0/65536.0*fine_values[16]-9.0/131072.0*fine_values[17]-fine_values[13]/
32768.0+3.0/65536.0*fine_values[14]-3.0/262144.0*fine_values[11]+fine_values
[12]/131072.0-fine_values[9]/524288.0+fine_values[10]/131072.0-3.0/262144.0*
fine_values[7]+9.0/524288.0*fine_values[8]-3.0/262144.0*fine_values[5]+3.0/
1048576.0*fine_values[6]+fine_values[4]/131072.0-fine_values[3]/524288.0+3.0/
1048576.0*fine_values[2]-fine_values[1]/524288.0+961.0/131072.0*fine_values[49]
+3.0/1048576.0*fine_values[47]+31.0/524288.0*fine_values[48]+fine_values[45]/
2097152.0+31.0/524288.0*fine_values[46]-93.0/65536.0*fine_values[43]+279.0/
131072.0*fine_values[44]+9.0/524288.0*fine_values[41]+93.0/262144.0*fine_values
[42]-3.0/262144.0*fine_values[40];
      s1 = s2+3.0/1048576.0*fine_values[39]-93.0/65536.0*fine_values[38]+31.0/
32768.0*fine_values[37]+fine_values[34]/131072.0-31.0/131072.0*fine_values[36]
-3.0/262144.0*fine_values[35]-31.0/131072.0*fine_values[31]+93.0/262144.0*
fine_values[32]-fine_values[33]/524288.0+31.0/524288.0*fine_values[30]-
fine_values[28]/524288.0+3.0/1048576.0*fine_values[29]+27.0/262144.0*
fine_values[26]+fine_values[27]/2097152.0+9.0/524288.0*fine_values[24]-9.0/
131072.0*fine_values[25]+fine_values[0]/2097152.0-9.0/131072.0*fine_values[23]+
3.0/65536.0*fine_values[22]+9.0/524288.0*fine_values[59]+93.0/262144.0*
fine_values[58]+3.0/1048576.0*fine_values[57]-93.0/65536.0*fine_values[56]-3.0/
262144.0*fine_values[53]-961.0/32768.0*fine_values[55]-31.0/131072.0*
fine_values[54]-31.0/131072.0*fine_values[52]+93.0/262144.0*fine_values[50]-
fine_values[51]/524288.0-3.0/262144.0*fine_values[73]-93.0/65536.0*fine_values
[74];
      s2 = s1+3.0/1048576.0*fine_values[71]+93.0/262144.0*fine_values[72]+31.0/
32768.0*fine_values[70]+fine_values[69]/131072.0-31.0/131072.0*fine_values[68]-
fine_values[65]/524288.0-31.0/131072.0*fine_values[66]-fine_values[67]/524288.0
+31.0/524288.0*fine_values[64]+2883.0/65536.0*fine_values[61]+279.0/131072.0*
fine_values[62]+fine_values[63]/2097152.0-31.0/131072.0*fine_values[101]-961.0/
32768.0*fine_values[102]+31.0/524288.0*fine_values[99]+961.0/131072.0*
fine_values[100]+93.0/262144.0*fine_values[98]-31.0/131072.0*fine_values[96]+
3.0/1048576.0*fine_values[97]+31.0/524288.0*fine_values[94]-fine_values[95]/
524288.0-93.0/65536.0*fine_values[91]+fine_values[93]/2097152.0+279.0/131072.0*
fine_values[92]-93.0/65536.0*fine_values[89]+93.0/262144.0*fine_values[90]+31.0
/32768.0*fine_values[88]-31.0/131072.0*fine_values[87]+93.0/262144.0*
fine_values[86]-31.0/131072.0*fine_values[85]-3.0/262144.0*fine_values[82];
      s3 = s2+9.0/524288.0*fine_values[83]+31.0/524288.0*fine_values[84]+3.0/
1048576.0*fine_values[81]-3.0/262144.0*fine_values[80]+fine_values[79]/131072.0
-fine_values[78]/524288.0+fine_values[75]/2097152.0-fine_values[76]/524288.0+
3.0/1048576.0*fine_values[77]+93.0/262144.0*fine_values[60]+961.0/131072.0*
fine_values[112]+31.0/524288.0*fine_values[111]+3.0/1048576.0*fine_values[109]+
93.0/262144.0*fine_values[110]+fine_values[105]/2097152.0;
      coarse_values[93] = s3+31.0/524288.0*fine_values[107]+961.0/131072.0*
fine_values[108]+31.0/524288.0*fine_values[106]+93.0/262144.0*fine_values[103]+
2883.0/65536.0*fine_values[104]-961.0/32768.0*fine_values[124]-31.0/131072.0*
fine_values[122]+961.0/131072.0*fine_values[123]-31.0/131072.0*fine_values[120]
+31.0/524288.0*fine_values[121]+fine_values[117]/2097152.0+2883.0/65536.0*
fine_values[116]+31.0/524288.0*fine_values[119]-fine_values[118]/524288.0+93.0/
262144.0*fine_values[115]+961.0/131072.0*fine_values[113]+29791.0/32768.0*
fine_values[114];
      coarse_values[94] = 31.0/4096.0*fine_values[48]+fine_values[45]/16384.0+
9.0/4096.0*fine_values[41]-3.0/2048.0*fine_values[40]+3.0/8192.0*fine_values
[39]+fine_values[34]/1024.0-3.0/2048.0*fine_values[35]-fine_values[33]/4096.0-
fine_values[28]/4096.0+3.0/8192.0*fine_values[29]+fine_values[27]/16384.0+3.0/
8192.0*fine_values[57]-31.0/1024.0*fine_values[54]-fine_values[51]/4096.0-31.0/
1024.0*fine_values[101]+31.0/4096.0*fine_values[99]+3.0/8192.0*fine_values[97]-
fine_values[95]/4096.0+fine_values[93]/16384.0+93.0/2048.0*fine_values[60]+
961.0/1024.0*fine_values[112]+31.0/4096.0*fine_values[111]+fine_values[105]/
16384.0+31.0/4096.0*fine_values[106]+93.0/2048.0*fine_values[103];
      coarse_values[95] = 3.0/64.0*fine_values[71]-fine_values[67]/32.0+
fine_values[63]/128.0+31.0/32.0*fine_values[121]+fine_values[117]/128.0;
      coarse_values[96] = 3.0/8192.0*fine_values[47]+fine_values[45]/16384.0-
fine_values[46]/4096.0+9.0/4096.0*fine_values[59]-3.0/2048.0*fine_values[58]+
3.0/8192.0*fine_values[57]-3.0/2048.0*fine_values[53]+fine_values[52]/1024.0-
fine_values[51]/4096.0+93.0/2048.0*fine_values[73]+3.0/8192.0*fine_values[71]
-31.0/1024.0*fine_values[69]+31.0/4096.0*fine_values[65]-fine_values[67]/4096.0
+fine_values[63]/16384.0+31.0/4096.0*fine_values[111]+3.0/8192.0*fine_values
[109]+fine_values[105]/16384.0-fine_values[107]/4096.0+961.0/1024.0*fine_values
[122]+31.0/4096.0*fine_values[121]+fine_values[117]/16384.0+31.0/4096.0*
fine_values[118]+93.0/2048.0*fine_values[115]-31.0/1024.0*fine_values[113];
      coarse_values[97] = fine_values[47]/128.0-fine_values[53]/32.0+3.0/64.0*
fine_values[59]+fine_values[109]/128.0+31.0/32.0*fine_values[115];
      coarse_values[98] = 3.0/8192.0*fine_values[47]+fine_values[45]/16384.0+
31.0/4096.0*fine_values[46]+9.0/4096.0*fine_values[59]+93.0/2048.0*fine_values
[58]+3.0/8192.0*fine_values[57]-3.0/2048.0*fine_values[53]-31.0/1024.0*
fine_values[52]-fine_values[51]/4096.0-3.0/2048.0*fine_values[73]+3.0/8192.0*
fine_values[71]+fine_values[69]/1024.0-fine_values[65]/4096.0-fine_values[67]/
4096.0+fine_values[63]/16384.0+31.0/4096.0*fine_values[111]+3.0/8192.0*
fine_values[109]+fine_values[105]/16384.0+31.0/4096.0*fine_values[107]-31.0/
1024.0*fine_values[122]+31.0/4096.0*fine_values[121]+fine_values[117]/16384.0-
fine_values[118]/4096.0+93.0/2048.0*fine_values[115]+961.0/1024.0*fine_values
[113];
      coarse_values[99] = fine_values[45]/128.0+3.0/64.0*fine_values[57]-
fine_values[51]/32.0+31.0/32.0*fine_values[111]+fine_values[105]/128.0;
      coarse_values[100] = fine_values[75];
      coarse_values[101] = fine_values[93]/128.0-fine_values[94]/32.0+3.0/64.0*
fine_values[81]+31.0/32.0*fine_values[78]+fine_values[75]/128.0;
      coarse_values[102] = fine_values[81];
      coarse_values[103] = fine_values[93]/128.0+31.0/32.0*fine_values[94]+3.0/
64.0*fine_values[81]-fine_values[78]/32.0+fine_values[75]/128.0;
      coarse_values[104] = fine_values[93];
      coarse_values[105] = 3.0/64.0*fine_values[77]+fine_values[75]/128.0+31.0/
32.0*fine_values[76]-fine_values[119]/32.0+fine_values[117]/128.0;
      coarse_values[106] = -3.0/2048.0*fine_values[98]-31.0/1024.0*fine_values
[96]+3.0/8192.0*fine_values[97]-fine_values[94]/4096.0+31.0/4096.0*fine_values
[95]+fine_values[93]/16384.0+93.0/2048.0*fine_values[82]+9.0/4096.0*fine_values
[83]+3.0/8192.0*fine_values[81]+93.0/2048.0*fine_values[80]+961.0/1024.0*
fine_values[79]+31.0/4096.0*fine_values[78]+fine_values[75]/16384.0+31.0/4096.0
*fine_values[76]+3.0/8192.0*fine_values[77]+3.0/8192.0*fine_values[109]-3.0/
2048.0*fine_values[110]+fine_values[105]/16384.0-fine_values[107]/4096.0+
fine_values[108]/1024.0-fine_values[106]/4096.0-31.0/1024.0*fine_values[120]+
fine_values[117]/16384.0-fine_values[119]/4096.0+31.0/4096.0*fine_values[118];
      coarse_values[107] = fine_values[81]/128.0+31.0/32.0*fine_values[82]+3.0/
64.0*fine_values[83]-fine_values[110]/32.0+fine_values[109]/128.0;
      coarse_values[108] = 93.0/2048.0*fine_values[98]+961.0/1024.0*fine_values
[96]+3.0/8192.0*fine_values[97]+31.0/4096.0*fine_values[94]+31.0/4096.0*
fine_values[95]+fine_values[93]/16384.0+93.0/2048.0*fine_values[82]+9.0/4096.0*
fine_values[83]+3.0/8192.0*fine_values[81]-3.0/2048.0*fine_values[80]-31.0/
1024.0*fine_values[79]-fine_values[78]/4096.0+fine_values[75]/16384.0+31.0/
4096.0*fine_values[76]+3.0/8192.0*fine_values[77]+3.0/8192.0*fine_values[109]
-3.0/2048.0*fine_values[110]+fine_values[105]/16384.0+31.0/4096.0*fine_values
[107]-31.0/1024.0*fine_values[108]-fine_values[106]/4096.0+fine_values[120]/
1024.0+fine_values[117]/16384.0-fine_values[119]/4096.0-fine_values[118]/4096.0
;
      coarse_values[109] = 3.0/64.0*fine_values[97]+31.0/32.0*fine_values[95]+
fine_values[93]/128.0-fine_values[106]/32.0+fine_values[105]/128.0;
      coarse_values[110] = fine_values[77];
      coarse_values[111] = -fine_values[98]/32.0+fine_values[97]/128.0+3.0/64.0
*fine_values[83]+31.0/32.0*fine_values[80]+fine_values[77]/128.0;
      coarse_values[112] = fine_values[83];
      coarse_values[113] = 31.0/32.0*fine_values[98]+fine_values[97]/128.0+3.0/
64.0*fine_values[83]-fine_values[80]/32.0+fine_values[77]/128.0;
      coarse_values[114] = fine_values[97];
      coarse_values[115] = 3.0/64.0*fine_values[77]+fine_values[75]/128.0-
fine_values[76]/32.0+31.0/32.0*fine_values[119]+fine_values[117]/128.0;
      coarse_values[116] = -3.0/2048.0*fine_values[98]+fine_values[96]/1024.0+
3.0/8192.0*fine_values[97]-fine_values[94]/4096.0-fine_values[95]/4096.0+
fine_values[93]/16384.0-3.0/2048.0*fine_values[82]+9.0/4096.0*fine_values[83]+
3.0/8192.0*fine_values[81]+93.0/2048.0*fine_values[80]-31.0/1024.0*fine_values
[79]+31.0/4096.0*fine_values[78]+fine_values[75]/16384.0-fine_values[76]/4096.0
+3.0/8192.0*fine_values[77]+3.0/8192.0*fine_values[109]+93.0/2048.0*fine_values
[110]+fine_values[105]/16384.0-fine_values[107]/4096.0-31.0/1024.0*fine_values
[108]+31.0/4096.0*fine_values[106]+961.0/1024.0*fine_values[120]+fine_values
[117]/16384.0+31.0/4096.0*fine_values[119]+31.0/4096.0*fine_values[118];
      coarse_values[117] = fine_values[81]/128.0-fine_values[82]/32.0+3.0/64.0*
fine_values[83]+31.0/32.0*fine_values[110]+fine_values[109]/128.0;
      coarse_values[118] = 93.0/2048.0*fine_values[98]-31.0/1024.0*fine_values
[96]+3.0/8192.0*fine_values[97]+31.0/4096.0*fine_values[94]-fine_values[95]/
4096.0+fine_values[93]/16384.0-3.0/2048.0*fine_values[82]+9.0/4096.0*
fine_values[83]+3.0/8192.0*fine_values[81]-3.0/2048.0*fine_values[80]+
fine_values[79]/1024.0-fine_values[78]/4096.0+fine_values[75]/16384.0-
fine_values[76]/4096.0+3.0/8192.0*fine_values[77]+3.0/8192.0*fine_values[109]+
93.0/2048.0*fine_values[110]+fine_values[105]/16384.0+31.0/4096.0*fine_values
[107]+961.0/1024.0*fine_values[108]+31.0/4096.0*fine_values[106]-31.0/1024.0*
fine_values[120]+fine_values[117]/16384.0+31.0/4096.0*fine_values[119]-
fine_values[118]/4096.0;
      coarse_values[119] = 3.0/64.0*fine_values[97]-fine_values[95]/32.0+
fine_values[93]/128.0+31.0/32.0*fine_values[106]+fine_values[105]/128.0;
      coarse_values[120] = fine_values[117];
      coarse_values[121] = 3.0/64.0*fine_values[109]-fine_values[107]/32.0+
fine_values[105]/128.0+31.0/32.0*fine_values[118]+fine_values[117]/128.0;
      coarse_values[122] = fine_values[109];
      coarse_values[123] = 3.0/64.0*fine_values[109]+31.0/32.0*fine_values[107]
+fine_values[105]/128.0-fine_values[118]/32.0+fine_values[117]/128.0;
      coarse_values[124] = fine_values[105]; 


}

void Transform_P1P2_1_3D(double *fine_values, double *coarse_values)
{
  coarse_values[0] = fine_values[12]/8.0
                     +fine_values[28]/8.0
                     +fine_values[24]/8.0
                     +fine_values[20]/8.0
                     +fine_values[16]/8.0
                     +fine_values[8]/8.0
                     +fine_values[4]/8.0
                     +fine_values[0]/8.0;

  coarse_values[1] = -fine_values[26]/16.0
                     +fine_values[14]/16.0
                     -3.0/16.0*fine_values[12]
                     -fine_values[6]/16.0
                     +fine_values[1]/16.0
                     +fine_values[29]/16.0
                     -3.0/16.0*fine_values[28]
                     +3.0/16.0*fine_values[24]
                     -fine_values[21]/16.0
                     +3.0/16.0*fine_values[20]
                     +fine_values[18]/16.0
                     -3.0/16.0*fine_values[16]
                     -fine_values[9]/16.0
                     +3.0/16.0*fine_values[8]
                     +3.0/16.0*fine_values[4]
                     -3.0/16.0*fine_values[0];

  coarse_values[2] = 3.0/16.0*fine_values[12]
                     -fine_values[30]/16.0
                     +3.0/16.0*fine_values[28]
                     -fine_values[25]/16.0
                     +3.0/16.0*fine_values[24]
                     +fine_values[22]/16.0
                     -3.0/16.0*fine_values[20]
                     +fine_values[17]/16.0
                     -3.0/16.0*fine_values[16]
                     -fine_values[13]/16.0
                     -fine_values[10]/16.0
                     +3.0/16.0*fine_values[8]
                     +fine_values[5]/16.0
                     -3.0/16.0*fine_values[4]
                     +fine_values[2]/16.0
                     -3.0/16.0*fine_values[0];

  coarse_values[3] = -fine_values[19]/16.0
                     -3.0/16.0*fine_values[12]
                     -fine_values[31]/16.0
                     +3.0/16.0*fine_values[28]
                     -fine_values[27]/16.0
                     +3.0/16.0*fine_values[24]
                     -fine_values[23]/16.0
                     +3.0/16.0*fine_values[20]
                     +3.0/16.0*fine_values[16]
                     +fine_values[15]/16.0
                     +fine_values[11]/16.0
                     -3.0/16.0*fine_values[8]
                     +fine_values[7]/16.0
                     -3.0/16.0*fine_values[4]
                     +fine_values[3]/16.0
                     -3.0/16.0*fine_values[0];

  coarse_values[4] = -fine_values[26]/12.0
                     -fine_values[14]/12.0
                     -fine_values[6]/12.0
                     -fine_values[1]/12.0
                     -fine_values[29]/12.0
                     -fine_values[21]/12.0
                     -fine_values[18]/12.0
                     -fine_values[9]/12.0;

  coarse_values[5] = -fine_values[12]/2.0
                     -fine_values[28]/2.0
                     +fine_values[24]/2.0
                     -fine_values[20]/2.0
                     +fine_values[16]/2.0
                     +fine_values[8]/2.0
                     -fine_values[4]/2.0
                     +fine_values[0]/2.0;

  coarse_values[6] = fine_values[12]/2.0
                     -fine_values[28]/2.0
                     +fine_values[24]/2.0
                     +fine_values[20]/2.0
                     -fine_values[16]/2.0
                     -fine_values[8]/2.0
                     -fine_values[4]/2.0
                     +fine_values[0]/2.0;

  coarse_values[7] = -fine_values[30]/12.0
                     -fine_values[25]/12.0
                     -fine_values[22]/12.0
                     -fine_values[17]/12.0
                     -fine_values[13]/12.0
                     -fine_values[10]/12.0
                     -fine_values[5]/12.0
                     -fine_values[2]/12.0;

  coarse_values[8] = -fine_values[12]/2.0
                     +fine_values[28]/2.0
                     +fine_values[24]/2.0
                     -fine_values[20]/2.0
                     -fine_values[16]/2.0
                     -fine_values[8]/2.0
                     +fine_values[4]/2.0
                     +fine_values[0]/2.0;

  coarse_values[9] = -fine_values[19]/12.0
                     -fine_values[31]/12.0
                     -fine_values[27]/12.0
                     -fine_values[23]/12.0
                     -fine_values[15]/12.0
                     -fine_values[11]/12.0
                     -fine_values[7]/12.0
                     -fine_values[3]/12.0;

}


void Transform_P1P2_2_3D(double *fine_values, double *coarse_values)
{
coarse_values[0] = fine_values[0]/8.0+fine_values[20]/8.0+fine_values[24]
/8.0+fine_values[12]/8.0+fine_values[16]/8.0+fine_values[8]/8.0+fine_values[4]/
8.0+fine_values[28]/8.0;
      coarse_values[1] = -3.0/16.0*fine_values[0]+fine_values[2]/16.0+
fine_values[18]-3.0/16.0*fine_values[20]-fine_values[21]+fine_values[22]/16.0+
fine_values[23]+3.0/16.0*fine_values[24]-fine_values[10]/16.0+3.0/16.0*
fine_values[12]-fine_values[13]/16.0-3.0/16.0*fine_values[16]-15.0/16.0*
fine_values[17]+fine_values[5]/16.0+3.0/16.0*fine_values[8]-3.0/16.0*
fine_values[4]-fine_values[25]/16.0+3.0/16.0*fine_values[28]+15.0/16.0*
fine_values[30]-fine_values[31];
      coarse_values[2] = -3.0/16.0*fine_values[0]+fine_values[2]/16.0-3.0/16.0*
fine_values[20]+fine_values[22]/16.0+3.0/16.0*fine_values[24]-fine_values[10]/
16.0+3.0/16.0*fine_values[12]-fine_values[13]/16.0-3.0/16.0*fine_values[16]+
fine_values[17]/16.0+fine_values[5]/16.0+3.0/16.0*fine_values[8]-3.0/16.0*
fine_values[4]-fine_values[25]/16.0+3.0/16.0*fine_values[28]-fine_values[30]/
16.0;
      coarse_values[3] = fine_values[1]/8.0-fine_values[2]/16.0+3.0/16.0*
fine_values[4]+9.0/16.0*fine_values[20]+7.0/8.0*fine_values[21]-fine_values[22]
/16.0-fine_values[23]+3.0/16.0*fine_values[24]+fine_values[25]/16.0+3.0/16.0*
fine_values[16]-fine_values[10]/16.0+fine_values[11]/8.0-9.0/16.0*fine_values
[12]+fine_values[13]/16.0+fine_values[15]/8.0+13.0/16.0*fine_values[17]-7.0/8.0
*fine_values[18]-fine_values[6]/8.0-3.0/16.0*fine_values[8]+fine_values[5]/16.0
-fine_values[27]/8.0-3.0/16.0*fine_values[28]-13.0/16.0*fine_values[30]+7.0/8.0
*fine_values[31]-3.0/16.0*fine_values[0];
      coarse_values[4] = -fine_values[1]/6.0+fine_values[2]/12.0-7.0/6.0*
fine_values[18]+fine_values[20]/2.0+7.0/6.0*fine_values[21]-fine_values[22]/
12.0-4.0/3.0*fine_values[23]-fine_values[9]/3.0+fine_values[10]/4.0-fine_values
[12]/2.0+fine_values[13]/12.0+13.0/12.0*fine_values[17]-fine_values[5]/12.0-
fine_values[6]/6.0+fine_values[4]/2.0-fine_values[25]/12.0-fine_values[28]/2.0
-5.0/4.0*fine_values[30]+4.0/3.0*fine_values[31];
      coarse_values[5] = 3.0/8.0*fine_values[0]-fine_values[2]/8.0-3.0/8.0*
fine_values[20]+fine_values[22]/8.0+3.0/8.0*fine_values[24]-fine_values[10]/8.0
-3.0/8.0*fine_values[12]+fine_values[13]/8.0+3.0/8.0*fine_values[16]-
fine_values[17]/8.0+fine_values[5]/8.0+3.0/8.0*fine_values[8]-3.0/8.0*
fine_values[4]-fine_values[25]/8.0-3.0/8.0*fine_values[28]+fine_values[30]/8.0;
      coarse_values[6] = -7.0/16.0*fine_values[1]+7.0/32.0*fine_values[2]-45.0/
32.0*fine_values[4]+9.0/32.0*fine_values[20]-25.0/16.0*fine_values[21]-
fine_values[22]/32.0+3.0/2.0*fine_values[23]+3.0/32.0*fine_values[24]-7.0/32.0*
fine_values[25]-21.0/32.0*fine_values[16]+7.0/32.0*fine_values[10]-3.0/16.0*
fine_values[11]+15.0/32.0*fine_values[12]+fine_values[13]/32.0-3.0/16.0*
fine_values[15]-43.0/32.0*fine_values[17]+25.0/16.0*fine_values[18]+7.0/16.0*
fine_values[6]-3.0/32.0*fine_values[8]+fine_values[5]/32.0+3.0/16.0*fine_values
[27]+21.0/32.0*fine_values[28]+35.0/32.0*fine_values[30]-21.0/16.0*fine_values
[31]+21.0/32.0*fine_values[0];
      coarse_values[7] = fine_values[1]/6.0-fine_values[2]/4.0-7.0/6.0*
fine_values[18]+fine_values[20]/2.0+7.0/6.0*fine_values[21]-fine_values[22]/
12.0-4.0/3.0*fine_values[23]-fine_values[10]/12.0-fine_values[12]/2.0+
fine_values[13]/12.0+13.0/12.0*fine_values[17]-fine_values[5]/12.0-fine_values
[6]/6.0+fine_values[4]/2.0-fine_values[25]/12.0-fine_values[28]/2.0-5.0/4.0*
fine_values[30]+4.0/3.0*fine_values[31];
      coarse_values[8] = 3.0/8.0*fine_values[0]-fine_values[2]/8.0-3.0/8.0*
fine_values[20]+fine_values[22]/8.0+3.0/8.0*fine_values[24]+fine_values[10]/8.0
-3.0/8.0*fine_values[12]+fine_values[13]/8.0-3.0/8.0*fine_values[16]+
fine_values[17]/8.0-fine_values[5]/8.0-3.0/8.0*fine_values[8]+3.0/8.0*
fine_values[4]-fine_values[25]/8.0+3.0/8.0*fine_values[28]-fine_values[30]/8.0;
      coarse_values[9] = fine_values[1]/8.0-5.0/24.0*fine_values[2]-7.0/8.0*
fine_values[18]+fine_values[20]/4.0+7.0/8.0*fine_values[21]-fine_values[22]/
24.0-fine_values[23]-fine_values[15]/24.0-fine_values[10]/12.0-fine_values[11]/
24.0-fine_values[12]/4.0+fine_values[13]/24.0-fine_values[16]/8.0+5.0/6.0*
fine_values[17]-fine_values[5]/12.0-fine_values[6]/8.0+fine_values[8]/8.0+3.0/
8.0*fine_values[4]+5.0/24.0*fine_values[25]-7.0/24.0*fine_values[27]-3.0/8.0*
fine_values[28]-2.0/3.0*fine_values[30]+17.0/24.0*fine_values[31];

}





// Q1 to Q2
void Superconvergence_Q1Q2_3D(TFEFunction3D *q1_function, 
                            TFEFunction3D *q2_function) 
{
  /* build the Q2 space */  
 
  TFESpace3D *q1_space, *q2_space;
  double *q1_values, *q2_values;
  int q1_ndof, q2_ndof;
  TCollection *q1_coll, *q2_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3,*child_cell4;
  TBaseCell *child_cell5,*child_cell6,*child_cell7;

 
  double coarse_values[27];
  double fine_values[27];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3,n4,n5,n6,n7;
  int *dof;
  TAuxParam3D *aux3d;
  TFESpace3D *fesp3d[1];
  int j,k;
  
  q1_space=q1_function->GetFESpace3D();
  q2_space=q2_function->GetFESpace3D();
  
  q2_coll=q2_space->GetCollection();
  q1_coll=q1_space->GetCollection();
  
  q1_ndof=q1_space->GetN_DegreesOfFreedom();
  q1_values=q1_function->GetValues();
  q2_ndof=q2_space->GetN_DegreesOfFreedom();
  q2_values=q2_function->GetValues();
  

  coarse_GlobalNumbers=q2_space->GetGlobalNumbers();
  coarse_BeginIndex=q2_space->GetBeginIndex();
  fine_GlobalNumbers=q1_space->GetGlobalNumbers();
  fine_BeginIndex=q1_space->GetBeginIndex();
       
  /* enumerate the cells */
       
  k=q1_coll->GetN_Cells();
  for(j=0;j<k;j++) /* for macro cell */
    q1_coll->GetCell(j)->SetClipBoard(j);
       
  k=q2_coll->GetN_Cells();
  for(j=0;j<k;j++)
    q2_coll->GetCell(j)->SetClipBoard(j);
  
  for(j=0;j<k;j++)
  {
    coarse_cell=q2_coll->GetCell(j); 
    child_cell0=coarse_cell->GetChild(0);
    child_cell1=coarse_cell->GetChild(1);
    child_cell2=coarse_cell->GetChild(2);
    child_cell3=coarse_cell->GetChild(3);
    child_cell4=coarse_cell->GetChild(4);
    child_cell5=coarse_cell->GetChild(5);
    child_cell6=coarse_cell->GetChild(6);
    child_cell7=coarse_cell->GetChild(7);
    
    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
    n4=child_cell4->GetClipBoard();
    n5=child_cell5->GetClipBoard();
    n6=child_cell6->GetClipBoard();
    n7=child_cell7->GetClipBoard();
            
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n0];

    coarse_values[0]=q1_values[dof[0]];
    coarse_values[1]=q1_values[dof[1]];
    coarse_values[3]=q1_values[dof[2]];
    coarse_values[4]=q1_values[dof[3]];
    coarse_values[9]=q1_values[dof[4]];
    coarse_values[10]=q1_values[dof[5]];
    coarse_values[12]=q1_values[dof[6]];
    coarse_values[13]=q1_values[dof[7]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    coarse_values[2]=q1_values[dof[0]];
    coarse_values[5]=q1_values[dof[1]];
    coarse_values[11]=q1_values[dof[4]];
    coarse_values[14]=q1_values[dof[5]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    coarse_values[8]=q1_values[dof[0]];
    coarse_values[7]=q1_values[dof[1]];
    coarse_values[17]=q1_values[dof[4]];
    coarse_values[16]=q1_values[dof[5]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    coarse_values[6]=q1_values[dof[0]];
    coarse_values[15]=q1_values[dof[4]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n4];
    coarse_values[18]=q1_values[dof[0]];
    coarse_values[19]=q1_values[dof[2]];
    coarse_values[21]=q1_values[dof[1]];
    coarse_values[22]=q1_values[dof[3]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n5];
    coarse_values[20]=q1_values[dof[0]];
    coarse_values[23]=q1_values[dof[2]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n6];
    coarse_values[26]=q1_values[dof[0]];
    coarse_values[25]=q1_values[dof[2]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n7];
    coarse_values[24]=q1_values[dof[0]];

       
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<27;l++)
      q2_values[dof[l]]=coarse_values[l];
  }
}


// Q2 to Q3
void Superconvergence_Q2Q3_3D(TFEFunction3D *q2_function, 
                            TFEFunction3D *q3_function) 
{
  TFESpace3D *q2_space;
  TFESpace3D *q3_space;

  double *q2_values,*q3_values;
  int q2_ndof,q3_ndof;
  
  TCollection *q3_coll, *q2_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3,*child_cell4;
  TBaseCell *child_cell5,*child_cell6,*child_cell7;


  double coarse_values[64];
  double fine_values[125];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3,n4,n5,n6,n7;
  int *dof;
  TAuxParam3D *aux3d;
  TFESpace3D *fesp3d[1];
  double s1,s2,s3;
  int j,k;

  int dofc0[27]={0,1,2,3,4,5,6,7,8,
                 9,10,11,12,13,14,15,16,17,
                 18,19,20,21,22,23,24,25,26};

  int dofc1[27]={27,28,29,30,31,32,2,5,8,
                 33,34,35,36,37,38,11,14,17,
                 39,40,41,42,43,44,20,23,26};
       
  int dofc2[27]={45,46,47,48,49,50,29,32,8,
                 51,52,53,54,55,56,35,38,17,
                 57,58,59,60,61,62,41,44,26};

  int dofc3[27]={63,64,6,65,66,7,47,50,8,
                 67,68,15,69,70,16,53,56,17,
                 71,72,24,73,74,25,59,62,26};

  int dofc4[27]={75,76,77,78,79,80,81,82,83,
                 84,85,86,87,88,89,90,91,92,
                 18,21,24,19,22,25,20,23,26};

  int dofc5[27]={93,94,81,95,96,82,97,98,83,
                 99,100,90,101,102,91,103,104,92,
                 39,42,20,40,43,23,41,44,26};
    
  int dofc6[27]={105,106,97,107,108,98,109,110,83,
                 111,112,103,113,114,104,115,116,92,
                 57,60,41,58,61,44,59,62,26};

  int dofc7[27]={117,118,109,119,120,110,77,80,83,
                 121,122,115,123,124,116,86,89,92,
                 71,73,59,72,74,62,24,25,26};


  q2_space=q2_function->GetFESpace3D();
  q3_space=q3_function->GetFESpace3D();
  
  q2_coll=q2_space->GetCollection();
  q3_coll=q3_space->GetCollection();
  
  q2_ndof=q2_space->GetN_DegreesOfFreedom();
  q2_values=q2_function->GetValues();
  q3_ndof=q3_space->GetN_DegreesOfFreedom();
  q3_values=q3_function->GetValues();
  
  coarse_GlobalNumbers=q3_space->GetGlobalNumbers();
  coarse_BeginIndex=q3_space->GetBeginIndex();
  fine_GlobalNumbers=q2_space->GetGlobalNumbers();
  fine_BeginIndex=q2_space->GetBeginIndex();

  /* enumerate the cells */
 
  k=q2_coll->GetN_Cells();
  for(j=0;j<k;j++) /* for macro cell */
    q2_coll->GetCell(j)->SetClipBoard(j);
        
  k=q3_coll->GetN_Cells();
  for(j=0;j<k;j++)
    q3_coll->GetCell(j)->SetClipBoard(j);

  for(j=0;j<k;j++)
  {
    coarse_cell=q3_coll->GetCell(j); 
    child_cell0=coarse_cell->GetChild(0);
    child_cell1=coarse_cell->GetChild(1);
    child_cell2=coarse_cell->GetChild(2);
    child_cell3=coarse_cell->GetChild(3);
    child_cell4=coarse_cell->GetChild(4);
    child_cell5=coarse_cell->GetChild(5);
    child_cell6=coarse_cell->GetChild(6);
    child_cell7=coarse_cell->GetChild(7);            
    
    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
    n4=child_cell4->GetClipBoard();
    n5=child_cell5->GetClipBoard();
    n6=child_cell6->GetClipBoard();
    n7=child_cell7->GetClipBoard();        


    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    for(l=0;l<27;l++)
      fine_values[dofc0[l]]=q2_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    for(l=0;l<27;l++)
      fine_values[dofc1[l]]=q2_values[dof[l]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    for(l=0;l<27;l++)
      fine_values[dofc2[l]]=q2_values[dof[l]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    for(l=0;l<27;l++)
      fine_values[dofc3[l]]=q2_values[dof[l]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n4];
    for(l=0;l<27;l++)
      fine_values[dofc4[l]]=q2_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n5];
    for(l=0;l<27;l++)
      fine_values[dofc5[l]]=q2_values[dof[l]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n6];
    for(l=0;l<27;l++)
      fine_values[dofc6[l]]=q2_values[dof[l]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n7];
    for(l=0;l<27;l++)
      fine_values[dofc7[l]]=q2_values[dof[l]];
            
    Transform_Q2Q3_3D(fine_values, coarse_values);    
    
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<64;l++)
      q3_values[dof[l]]=coarse_values[l];                
  }
}

void Superconvergence_Q2Q4_3D(TFEFunction3D *q2_function, 
                            TFEFunction3D *q4_function) 
{
  TFESpace3D *q2_space;
  TFESpace3D *q4_space;

  double *q2_values,*q4_values;
  int q2_ndof,q4_ndof;
  
  TCollection *q4_coll, *q2_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3,*child_cell4;
  TBaseCell *child_cell5,*child_cell6,*child_cell7;


  double coarse_values[125];
  double fine_values[125];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3,n4,n5,n6,n7;
  int *dof;
  TAuxParam3D *aux3d;
  TFESpace3D *fesp3d[1];
  double s1,s2,s3;
  int j,k;

  int dofc0[27]={0,1,2,3,4,5,6,7,8,
                 9,10,11,12,13,14,15,16,17,
                 18,19,20,21,22,23,24,25,26};

  int dofc1[27]={27,28,29,30,31,32,2,5,8,
                 33,34,35,36,37,38,11,14,17,
                 39,40,41,42,43,44,20,23,26};
       
  int dofc2[27]={45,46,47,48,49,50,29,32,8,
                 51,52,53,54,55,56,35,38,17,
                 57,58,59,60,61,62,41,44,26};

  int dofc3[27]={63,64,6,65,66,7,47,50,8,
                 67,68,15,69,70,16,53,56,17,
                 71,72,24,73,74,25,59,62,26};

  int dofc4[27]={75,76,77,78,79,80,81,82,83,
                 84,85,86,87,88,89,90,91,92,
                 18,21,24,19,22,25,20,23,26};

  int dofc5[27]={93,94,81,95,96,82,97,98,83,
                 99,100,90,101,102,91,103,104,92,
                 39,42,20,40,43,23,41,44,26};
    
  int dofc6[27]={105,106,97,107,108,98,109,110,83,
                 111,112,103,113,114,104,115,116,92,
                 57,60,41,58,61,44,59,62,26};

  int dofc7[27]={117,118,109,119,120,110,77,80,83,
                 121,122,115,123,124,116,86,89,92,
                 71,73,59,72,74,62,24,25,26};

  
  q2_space=q2_function->GetFESpace3D();
  q4_space=q4_function->GetFESpace3D();
  
  q2_coll=q2_space->GetCollection();
  q4_coll=q4_space->GetCollection();
  
  q2_ndof=q2_space->GetN_DegreesOfFreedom();
  q2_values=q2_function->GetValues();
  q4_ndof=q4_space->GetN_DegreesOfFreedom();
  q4_values=q4_function->GetValues();
  
  coarse_GlobalNumbers=q4_space->GetGlobalNumbers();
  coarse_BeginIndex=q4_space->GetBeginIndex();
  fine_GlobalNumbers=q2_space->GetGlobalNumbers();
  fine_BeginIndex=q2_space->GetBeginIndex();

  /* enumerate the cells */
 
  k=q2_coll->GetN_Cells();
  for(j=0;j<k;j++) /* for macro cell */
    q2_coll->GetCell(j)->SetClipBoard(j);
        
  k=q4_coll->GetN_Cells();
  for(j=0;j<k;j++)
    q4_coll->GetCell(j)->SetClipBoard(j);

  for(j=0;j<k;j++)
  {
    coarse_cell=q4_coll->GetCell(j); 
    
    child_cell0=coarse_cell->GetChild(0);
    child_cell1=coarse_cell->GetChild(1);
    child_cell2=coarse_cell->GetChild(2);
    child_cell3=coarse_cell->GetChild(3);
    child_cell4=coarse_cell->GetChild(4);
    child_cell5=coarse_cell->GetChild(5);
    child_cell6=coarse_cell->GetChild(6);
    child_cell7=coarse_cell->GetChild(7);

    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
    n4=child_cell4->GetClipBoard();
    n5=child_cell5->GetClipBoard();
    n6=child_cell6->GetClipBoard();
    n7=child_cell7->GetClipBoard();


    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    for(l=0;l<27;l++)
      fine_values[dofc0[l]]=q2_values[dof[l]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    for(l=0;l<27;l++)
      fine_values[dofc1[l]]=q2_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    for(l=0;l<27;l++)
      fine_values[dofc2[l]]=q2_values[dof[l]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    for(l=0;l<27;l++)
      fine_values[dofc3[l]]=q2_values[dof[l]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n4];
    for(l=0;l<27;l++)
      fine_values[dofc4[l]]=q2_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n5];
    for(l=0;l<27;l++)
      fine_values[dofc5[l]]=q2_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n6];
    for(l=0;l<27;l++)
      fine_values[dofc6[l]]=q2_values[dof[l]];
            
    dof=fine_GlobalNumbers+fine_BeginIndex[n7];
    for(l=0;l<27;l++)
      fine_values[dofc7[l]]=q2_values[dof[l]];
            

    Transform_Q2Q4_3D(fine_values, coarse_values);    

    
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<125;l++)
      q4_values[dof[l]]=coarse_values[l];                

  }
}


void Superconvergence_P1P2_3D(int version, TFEFunction3D *p1_function, 
                            TFEFunction3D *p2_function) 
{
  TFESpace3D *p1_space, *p2_space;
  double *p1_values, *p2_values;
  int p1_ndof, p2_ndof;
  TCollection *p1_coll, *p2_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3,*child_cell4;
  TBaseCell *child_cell5,*child_cell6,*child_cell7;
  double coarse_values[10];
  double fine_values[32];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3,n4,n5,n6,n7;
  int *dof;
  TAuxParam3D *aux3d;
  TFESpace3D *fesp3d[1];
  int j,k;
  
  p1_space=p1_function->GetFESpace3D();
  p2_space=p2_function->GetFESpace3D();
  
  p2_coll=p2_space->GetCollection();
  p1_coll=p1_space->GetCollection();
  
  p1_ndof=p1_space->GetN_DegreesOfFreedom();
  p1_values=p1_function->GetValues();
  p2_ndof=p2_space->GetN_DegreesOfFreedom();
  p2_values=p2_function->GetValues();
  

  coarse_GlobalNumbers=p2_space->GetGlobalNumbers();
  coarse_BeginIndex=p2_space->GetBeginIndex();
  fine_GlobalNumbers=p1_space->GetGlobalNumbers();
  fine_BeginIndex=p1_space->GetBeginIndex();
       
     /* enumerate the cells */
       
  k=p1_coll->GetN_Cells();
  for(j=0;j<k;j++) /* for macro cell */
    p1_coll->GetCell(j)->SetClipBoard(j);
       
  k=p2_coll->GetN_Cells();
  for(j=0;j<k;j++)
    p2_coll->GetCell(j)->SetClipBoard(j);
  
  for(j=0;j<k;j++)
  {
    coarse_cell=p2_coll->GetCell(j); 
    child_cell0=coarse_cell->GetChild(0);
    child_cell1=coarse_cell->GetChild(1);
    child_cell2=coarse_cell->GetChild(2);
    child_cell3=coarse_cell->GetChild(3);
    child_cell4=coarse_cell->GetChild(4);
    child_cell5=coarse_cell->GetChild(5);
    child_cell6=coarse_cell->GetChild(6);
    child_cell7=coarse_cell->GetChild(7);

    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
    n4=child_cell4->GetClipBoard();
    n5=child_cell5->GetClipBoard();
    n6=child_cell6->GetClipBoard();
    n7=child_cell7->GetClipBoard();
    

    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    for(l=0;l<4;l++)
      fine_values[l]=p1_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    for(l=0;l<4;l++)
      fine_values[4+l]=p1_values[dof[l]];
 
    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    for(l=0;l<4;l++)
      fine_values[8+l]=p1_values[dof[l]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    for(l=0;l<4;l++)
      fine_values[12+l]=p1_values[dof[l]];
 
    dof=fine_GlobalNumbers+fine_BeginIndex[n4];
    for(l=0;l<4;l++)
      fine_values[16+l]=p1_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n5];
    for(l=0;l<4;l++)
      fine_values[20+l]=p1_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n6];
    for(l=0;l<4;l++)
      fine_values[24+l]=p1_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n7];
    for(l=0;l<4;l++)
      fine_values[28+l]=p1_values[dof[l]];

    if(version==0) 
      Transform_P1P2_1_3D(fine_values, coarse_values);
    else
     Transform_P1P2_2_3D(fine_values, coarse_values);

    
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<10;l++)
      p2_values[dof[l]]=coarse_values[l];
  }
}
