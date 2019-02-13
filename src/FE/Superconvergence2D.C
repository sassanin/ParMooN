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
//     Stokes and Navier-Strokes 2D problems
//*****************************************************
//
#include <Superconvergence2D.h>

void Transform_NQ1P2_2D(double *fine_values, double *coarse_values)
{


  coarse_values[0] = fine_values[4]/4.0
                           +fine_values[8]/4.0
                           +fine_values[12]/4.0
                           +fine_values[15]/4.0;

  coarse_values[1] = -fine_values[3]/12.0
                     +fine_values[5]/12.0
                     +fine_values[11]/12.0
                     -fine_values[14]/12.0;

  coarse_values[2] = -fine_values[0]/12.0
                     -fine_values[7]/12.0
                     +fine_values[9]/12.0
                     +fine_values[13]/12.0;

  coarse_values[3] = fine_values[3]/10.0
                     -fine_values[4]/10.0
                     +fine_values[5]/10.0
                     -fine_values[8]/10.0
                     +fine_values[11]/10.0
                     -fine_values[12]/10.0
                     +fine_values[14]/10.0
                     -fine_values[15]/10.0;

  coarse_values[4] = fine_values[4]/9.0
                     -fine_values[8]/9.0
                     +fine_values[12]/9.0
                     -fine_values[15]/9.0;

  coarse_values[5] = fine_values[0]/10.0
                     -fine_values[4]/10.0
                     +fine_values[7]/10.0
                     -fine_values[8]/10.0
                     +fine_values[9]/10.0
                     -fine_values[12]/10.0
                     +fine_values[13]/10.0
                     -fine_values[15]/10.0;

}


void Transform_Q2Q3_2D(double *fine_values, double *coarse_values)
{
  coarse_values[0] = fine_values[0];
  coarse_values[1] = -7.0/81.0*fine_values[0]
                     +68.0/81.0*fine_values[1]
                     +2.0/9.0*fine_values[2]
                     -2.0/81.0*fine_values[9]
                     +4.0/81.0*fine_values[12];

  coarse_values[2] = -2.0/81.0*fine_values[0]+4.0/81.0*fine_values[1]+2.0/
9.0*fine_values[2]-7.0/81.0*fine_values[9]+68.0/81.0*fine_values[12];
      coarse_values[3] = fine_values[9];
      coarse_values[4] = -2.0/81.0*fine_values[21]+4.0/81.0*fine_values[22]+2.0
/9.0*fine_values[6]-7.0/81.0*fine_values[0]+68.0/81.0*fine_values[3];
      coarse_values[5] = -8.0/6561.0*fine_values[18]+16.0/6561.0*fine_values
[19]+8.0/729.0*fine_values[20]+14.0/6561.0*fine_values[21]-28.0/6561.0*
fine_values[22]-136.0/6561.0*fine_values[23]+272.0/6561.0*fine_values[24]+136.0
/729.0*fine_values[7]-476.0/6561.0*fine_values[1]-14.0/729.0*fine_values[2]
-476.0/6561.0*fine_values[3]+4624.0/6561.0*fine_values[4]+136.0/729.0*
fine_values[5]-14.0/729.0*fine_values[6]+4.0/81.0*fine_values[8]+14.0/6561.0*
fine_values[9]-136.0/6561.0*fine_values[10]-4.0/729.0*fine_values[11]-28.0/
6561.0*fine_values[12]+272.0/6561.0*fine_values[13]+8.0/729.0*fine_values[14]+
4.0/6561.0*fine_values[15]-8.0/6561.0*fine_values[16]-4.0/729.0*fine_values[17]
+49.0/6561.0*fine_values[0];
      coarse_values[6] = -28.0/6561.0*fine_values[18]+272.0/6561.0*fine_values
[19]+8.0/729.0*fine_values[20]+4.0/6561.0*fine_values[21]-8.0/6561.0*
fine_values[22]-8.0/6561.0*fine_values[23]+16.0/6561.0*fine_values[24]+8.0/
729.0*fine_values[7]-28.0/6561.0*fine_values[1]-14.0/729.0*fine_values[2]-136.0
/6561.0*fine_values[3]+272.0/6561.0*fine_values[4]+136.0/729.0*fine_values[5]
-4.0/729.0*fine_values[6]+4.0/81.0*fine_values[8]+49.0/6561.0*fine_values[9]
-476.0/6561.0*fine_values[10]-14.0/729.0*fine_values[11]-476.0/6561.0*
fine_values[12]+4624.0/6561.0*fine_values[13]+136.0/729.0*fine_values[14]+14.0/
6561.0*fine_values[15]-136.0/6561.0*fine_values[16]-4.0/729.0*fine_values[17]+
14.0/6561.0*fine_values[0];
      coarse_values[7] = 4.0/81.0*fine_values[18]-7.0/81.0*fine_values[9]+68.0/
81.0*fine_values[10]+2.0/9.0*fine_values[11]-2.0/81.0*fine_values[15];
      coarse_values[8] = -7.0/81.0*fine_values[21]+68.0/81.0*fine_values[22]+
2.0/9.0*fine_values[6]-2.0/81.0*fine_values[0]+4.0/81.0*fine_values[3];
      coarse_values[9] = -136.0/6561.0*fine_values[18]+272.0/6561.0*fine_values
[19]+136.0/729.0*fine_values[20]+49.0/6561.0*fine_values[21]-476.0/6561.0*
fine_values[22]-476.0/6561.0*fine_values[23]+4624.0/6561.0*fine_values[24]+
136.0/729.0*fine_values[7]-136.0/6561.0*fine_values[1]-4.0/729.0*fine_values[2]
-28.0/6561.0*fine_values[3]+272.0/6561.0*fine_values[4]+8.0/729.0*fine_values
[5]-14.0/729.0*fine_values[6]+4.0/81.0*fine_values[8]+4.0/6561.0*fine_values[9]
-8.0/6561.0*fine_values[10]-4.0/729.0*fine_values[11]-8.0/6561.0*fine_values
[12]+16.0/6561.0*fine_values[13]+8.0/729.0*fine_values[14]+14.0/6561.0*
fine_values[15]-28.0/6561.0*fine_values[16]-14.0/729.0*fine_values[17]+14.0/
6561.0*fine_values[0];
      coarse_values[10] = -476.0/6561.0*fine_values[18]+4624.0/6561.0*
fine_values[19]+136.0/729.0*fine_values[20]+14.0/6561.0*fine_values[21]-136.0/
6561.0*fine_values[22]-28.0/6561.0*fine_values[23]+272.0/6561.0*fine_values[24]
+8.0/729.0*fine_values[7]-8.0/6561.0*fine_values[1]-4.0/729.0*fine_values[2]
-8.0/6561.0*fine_values[3]+16.0/6561.0*fine_values[4]+8.0/729.0*fine_values[5]
-4.0/729.0*fine_values[6]+4.0/81.0*fine_values[8]+14.0/6561.0*fine_values[9]
-28.0/6561.0*fine_values[10]-14.0/729.0*fine_values[11]-136.0/6561.0*
fine_values[12]+272.0/6561.0*fine_values[13]+136.0/729.0*fine_values[14]+49.0/
6561.0*fine_values[15]-476.0/6561.0*fine_values[16]-14.0/729.0*fine_values[17]+
4.0/6561.0*fine_values[0];
      coarse_values[11] = 68.0/81.0*fine_values[18]-2.0/81.0*fine_values[9]+4.0
/81.0*fine_values[10]+2.0/9.0*fine_values[11]-7.0/81.0*fine_values[15];
      coarse_values[12] = fine_values[21];
      coarse_values[13] = 2.0/9.0*fine_values[17]-7.0/81.0*fine_values[21]+68.0
/81.0*fine_values[23]-2.0/81.0*fine_values[15]+4.0/81.0*fine_values[16];
      coarse_values[14] = 2.0/9.0*fine_values[17]-2.0/81.0*fine_values[21]+4.0/
81.0*fine_values[23]-7.0/81.0*fine_values[15]+68.0/81.0*fine_values[16];
      coarse_values[15] = fine_values[15];
}

void Transform_Q2Q4_2D(double *fine_values, double *coarse_values)
{
 coarse_values[0] = fine_values[0];
      coarse_values[1] = 31.0/32.0*fine_values[1]+3.0/64.0*fine_values[2]+
fine_values[0]/128.0-fine_values[12]/32.0+fine_values[9]/128.0;
      coarse_values[2] = fine_values[2];
      coarse_values[3] = -fine_values[1]/32.0+3.0/64.0*fine_values[2]+
fine_values[0]/128.0+31.0/32.0*fine_values[12]+fine_values[9]/128.0;
      coarse_values[4] = fine_values[9];
      coarse_values[5] = 31.0/32.0*fine_values[3]+fine_values[0]/128.0-
fine_values[22]/32.0+fine_values[21]/128.0+3.0/64.0*fine_values[6];
      coarse_values[6] = fine_values[0]/16384.0+3.0/8192.0*fine_values[6]+9.0/
4096.0*fine_values[8]+fine_values[9]/16384.0+93.0/2048.0*fine_values[5]+3.0/
8192.0*fine_values[2]+31.0/4096.0*fine_values[3]+961.0/1024.0*fine_values[4]+
31.0/4096.0*fine_values[1]+31.0/4096.0*fine_values[23]-31.0/1024.0*fine_values
[24]-fine_values[18]/4096.0+fine_values[19]/1024.0-3.0/2048.0*fine_values[20]+
fine_values[21]/16384.0-fine_values[22]/4096.0+fine_values[15]/16384.0-
fine_values[16]/4096.0+3.0/8192.0*fine_values[17]-fine_values[12]/4096.0-31.0/
1024.0*fine_values[13]-3.0/2048.0*fine_values[14]+31.0/4096.0*fine_values[10]+
3.0/8192.0*fine_values[11]+93.0/2048.0*fine_values[7];
      coarse_values[7] = 31.0/32.0*fine_values[5]+3.0/64.0*fine_values[8]+
fine_values[2]/128.0+fine_values[17]/128.0-fine_values[20]/32.0;
      coarse_values[8] = fine_values[0]/16384.0+3.0/8192.0*fine_values[6]+9.0/
4096.0*fine_values[8]+fine_values[9]/16384.0+93.0/2048.0*fine_values[5]+3.0/
8192.0*fine_values[2]+31.0/4096.0*fine_values[3]-31.0/1024.0*fine_values[4]-
fine_values[1]/4096.0-fine_values[23]/4096.0+fine_values[24]/1024.0-fine_values
[18]/4096.0-31.0/1024.0*fine_values[19]-3.0/2048.0*fine_values[20]+fine_values
[21]/16384.0-fine_values[22]/4096.0+fine_values[15]/16384.0+31.0/4096.0*
fine_values[16]+3.0/8192.0*fine_values[17]+31.0/4096.0*fine_values[12]+961.0/
1024.0*fine_values[13]+93.0/2048.0*fine_values[14]+31.0/4096.0*fine_values[10]+
3.0/8192.0*fine_values[11]-3.0/2048.0*fine_values[7];
      coarse_values[9] = -fine_values[18]/32.0+fine_values[15]/128.0+3.0/64.0*
fine_values[11]+fine_values[9]/128.0+31.0/32.0*fine_values[10];
      coarse_values[10] = fine_values[6];
      coarse_values[11] = 31.0/32.0*fine_values[7]+3.0/64.0*fine_values[8]-
fine_values[14]/32.0+fine_values[11]/128.0+fine_values[6]/128.0;
      coarse_values[12] = fine_values[8];
      coarse_values[13] = -fine_values[7]/32.0+3.0/64.0*fine_values[8]+31.0/
32.0*fine_values[14]+fine_values[11]/128.0+fine_values[6]/128.0;
      coarse_values[14] = fine_values[11];
      coarse_values[15] = -fine_values[3]/32.0+fine_values[0]/128.0+31.0/32.0*
fine_values[22]+fine_values[21]/128.0+3.0/64.0*fine_values[6];
      coarse_values[16] = fine_values[0]/16384.0+3.0/8192.0*fine_values[6]+9.0/
4096.0*fine_values[8]+fine_values[9]/16384.0-3.0/2048.0*fine_values[5]+3.0/
8192.0*fine_values[2]-fine_values[3]/4096.0-31.0/1024.0*fine_values[4]+31.0/
4096.0*fine_values[1]+31.0/4096.0*fine_values[23]+961.0/1024.0*fine_values[24]+
31.0/4096.0*fine_values[18]-31.0/1024.0*fine_values[19]+93.0/2048.0*fine_values
[20]+fine_values[21]/16384.0+31.0/4096.0*fine_values[22]+fine_values[15]/
16384.0-fine_values[16]/4096.0+3.0/8192.0*fine_values[17]-fine_values[12]/
4096.0+fine_values[13]/1024.0-3.0/2048.0*fine_values[14]-fine_values[10]/4096.0
+3.0/8192.0*fine_values[11]+93.0/2048.0*fine_values[7];
      coarse_values[17] = -fine_values[5]/32.0+3.0/64.0*fine_values[8]+
fine_values[2]/128.0+fine_values[17]/128.0+31.0/32.0*fine_values[20];
      coarse_values[18] = fine_values[0]/16384.0+3.0/8192.0*fine_values[6]+9.0/
4096.0*fine_values[8]+fine_values[9]/16384.0-3.0/2048.0*fine_values[5]+3.0/
8192.0*fine_values[2]-fine_values[3]/4096.0+fine_values[4]/1024.0-fine_values
[1]/4096.0-fine_values[23]/4096.0-31.0/1024.0*fine_values[24]+31.0/4096.0*
fine_values[18]+961.0/1024.0*fine_values[19]+93.0/2048.0*fine_values[20]+
fine_values[21]/16384.0+31.0/4096.0*fine_values[22]+fine_values[15]/16384.0+
31.0/4096.0*fine_values[16]+3.0/8192.0*fine_values[17]+31.0/4096.0*fine_values
[12]-31.0/1024.0*fine_values[13]+93.0/2048.0*fine_values[14]-fine_values[10]/
4096.0+3.0/8192.0*fine_values[11]-3.0/2048.0*fine_values[7];
      coarse_values[19] = 31.0/32.0*fine_values[18]+fine_values[15]/128.0+3.0/
64.0*fine_values[11]+fine_values[9]/128.0-fine_values[10]/32.0;
      coarse_values[20] = fine_values[21];
      coarse_values[21] = 31.0/32.0*fine_values[23]+3.0/64.0*fine_values[17]+
fine_values[21]/128.0+fine_values[15]/128.0-fine_values[16]/32.0;
      coarse_values[22] = fine_values[17];
      coarse_values[23] = -fine_values[23]/32.0+3.0/64.0*fine_values[17]+
fine_values[21]/128.0+fine_values[15]/128.0+31.0/32.0*fine_values[16];
      coarse_values[24] = fine_values[15];
}

void Transform_P1P2_1_2D(double *fine_values, double *coarse_values) 
{
  coarse_values[0] = fine_values[0]/4.0
                     +fine_values[3]/4.0
                     +fine_values[6]/4.0
                     +fine_values[9]/4.0;

  coarse_values[1] = -fine_values[0]/8.0
                     +fine_values[1]/8.0
                     +fine_values[3]/8.0
                     -fine_values[5]/8.0
                     +fine_values[6]/8.0
                     -fine_values[7]/8.0
                     -fine_values[9]/8.0
                     +fine_values[11]/8.0;

  coarse_values[2] = -fine_values[0]/8.0
                     +fine_values[2]/8.0
                     -fine_values[3]/8.0
                     +fine_values[4]/8.0
                     +fine_values[6]/8.0
                     -fine_values[8]/8.0
                     +fine_values[9]/8.0
                     -fine_values[10]/8.0;
 
  coarse_values[3] = -2.0/5.0*fine_values[1]
                     -2.0/5.0*fine_values[5]
                     -2.0/5.0*fine_values[7]
                     -2.0/5.0*fine_values[11];

  coarse_values[4] = fine_values[0]/9.0
                     -fine_values[3]/9.0
                     +fine_values[6]/9.0
                     -fine_values[9]/9.0;

  coarse_values[5] = -2.0/5.0*fine_values[2]
                     -2.0/5.0*fine_values[4]
                     -2.0/5.0*fine_values[8]
                     -2.0/5.0*fine_values[10];

}


void Transform_P1P2_2_2D(double *fine_values, double *coarse_values) 
{
  coarse_values[0] = fine_values[0]/4.0
                     +fine_values[3]/4.0
                     +fine_values[6]/4.0
	             +fine_values[9]/4.0;
  
  coarse_values[1] = -fine_values[0]/7.0
                     +fine_values[1]/7.0
                     +fine_values[3]/7.0
	             -fine_values[5]/7.0
                     +fine_values[6]/7.0
                     -fine_values[8]/7.0
                     -fine_values[9]/7.0
	             +fine_values[10]/7.0;
  
  coarse_values[2] = -3.0/28.0*fine_values[0]
                     -fine_values[1]/56.0
	             +fine_values[2]/8.0
                     -fine_values[3]/7.0
                     +fine_values[4]/8.0
                     +fine_values[5]/56.0
	             +3.0/28.0*fine_values[6]
                     -fine_values[7]/8.0
                     +fine_values[8]/56.0
                     +fine_values[9]/7.0
                     -15.0/56.0*fine_values[10]
                     +fine_values[11]/8.0;
  
  coarse_values[3] = fine_values[0]/7.0
                     -59.0/70.0*fine_values[1]
	             +3.0/10.0*fine_values[2]
                     +2.0/35.0*fine_values[3]
                     +3.0/10.0*fine_values[4]
	             -53.0/70.0*fine_values[5]
                     -fine_values[6]/7.0
                     -3.0/10.0*fine_values[7]
	             +3.0/70.0*fine_values[8]
                     -2.0/35.0*fine_values[9]
                     -9.0/14.0*fine_values[10]
	             +3.0/10.0*fine_values[11];

  coarse_values[4] = fine_values[0]/14.0
                     -fine_values[1]/14.0
                     -fine_values[3]/14.0
                     +fine_values[5]/14.0
                     +2.0/21.0*fine_values[6]
                     -2.0/21.0*fine_values[8]
	             -2.0/21.0*fine_values[9]
                     +2.0/21.0*fine_values[10];

  coarse_values[5] = -2.0/35.0*fine_values[0]
                     +2.0/35.0*fine_values[1]
	             -2.0/5.0*fine_values[2]
                     +2.0/35.0*fine_values[3]
                     -2.0/5.0*fine_values[4]
	             -2.0/35.0*fine_values[5]
                     +2.0/35.0*fine_values[6]
                     +2.0/5.0*fine_values[7]
	             -6.0/7.0*fine_values[8]
                     -2.0/35.0*fine_values[9]
                     +2.0/35.0*fine_values[10]
	             -2.0/5.0*fine_values[11];

}


// Q1 to Q2
void Superconvergence_Q1Q2_2D(TFEFunction2D *q1_function, 
			    TFEFunction2D *q2_function) 
{
/* build the Q2 space */  
 
  TFESpace2D *q1_space, *q2_space;
  double *q1_values, *q2_values;
  int q1_ndof, q2_ndof;
  TCollection *q1_coll, *q2_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3;
  double coarse_values[9];
  double fine_values[9];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3;
  int *dof;
  TAuxParam2D *aux2d;
  TFESpace2D *fesp2d[1];
  int j,k;
  
  q1_space=q1_function->GetFESpace2D();
  q2_space=q2_function->GetFESpace2D();
  
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
    
    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    coarse_values[0]=q1_values[dof[0]];
    coarse_values[1]=q1_values[dof[1]];
    coarse_values[3]=q1_values[dof[2]];
    coarse_values[4]=q1_values[dof[3]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    coarse_values[2]=q1_values[dof[0]];
    coarse_values[5]=q1_values[dof[1]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    coarse_values[8]=q1_values[dof[0]];
    coarse_values[7]=q1_values[dof[1]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    coarse_values[6]=q1_values[dof[0]];
	  
       
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<9;l++)
      q2_values[dof[l]]=coarse_values[l];
       
  }
}


// Q2 to Q3
void Superconvergence_Q2Q3_2D(TFEFunction2D *q2_function, 
			    TFEFunction2D *q3_function) 
{
  TFESpace2D *q2_space;
  TFESpace2D *q3_space;

  double *q2_values,*q3_values;
  int q2_ndof,q3_ndof;
  
  TCollection *q3_coll, *q2_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3;
  double coarse_values[16];
  double fine_values[25];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3;
  int *dof;
  TAuxParam2D *aux2d;
  TFESpace2D *fesp2d[1];
  int j,k;

  q2_space=q2_function->GetFESpace2D();
  q3_space=q3_function->GetFESpace2D();
  
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
	    
    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
	    
    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    for(l=0;l<9;l++)
      fine_values[l]=q2_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    for(l=0;l<6;l++)
      fine_values[l+9]=q2_values[dof[l]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    for(l=0;l<6;l++)
      fine_values[l+15]=q2_values[dof[l]];
	     
    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    fine_values[21]=q2_values[dof[0]];
    fine_values[22]=q2_values[dof[1]];	
    fine_values[23]=q2_values[dof[3]];
    fine_values[24]=q2_values[dof[4]];
	    
    Transform_Q2Q3_2D(fine_values, coarse_values);    
    
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<16;l++)
      q3_values[dof[l]]=coarse_values[l];	        
  
  }
}

void Superconvergence_Q2Q4_2D(TFEFunction2D *q2_function, 
			    TFEFunction2D *q4_function) 
{
  TFESpace2D *q2_space;
  TFESpace2D *q4_space;

  double *q2_values,*q4_values;
  int q2_ndof,q4_ndof;
  
  TCollection *q4_coll, *q2_coll;
  TBaseCell *coarse_cell;
  TBaseCell *child_cell0,*child_cell1,*child_cell2,*child_cell3;
  double coarse_values[25];
  double fine_values[25];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3;
  int *dof;
  TAuxParam2D *aux2d;
  TFESpace2D *fesp2d[1];
  int j,k;
  
  q2_space=q2_function->GetFESpace2D();
  q4_space=q4_function->GetFESpace2D();
  
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
	    
    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
	    
    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    for(l=0;l<9;l++)
      fine_values[l]=q2_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    for(l=0;l<6;l++)
      fine_values[l+9]=q2_values[dof[l]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    for(l=0;l<6;l++)
      fine_values[l+15]=q2_values[dof[l]];
	     
    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    fine_values[21]=q2_values[dof[0]];
    fine_values[22]=q2_values[dof[1]];	
    fine_values[23]=q2_values[dof[3]];
    fine_values[24]=q2_values[dof[4]];
	    
    Transform_Q2Q4_2D(fine_values, coarse_values);    
    
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<25;l++)
      q4_values[dof[l]]=coarse_values[l];	        

  }
}

void Superconvergence_P1P2_2D(int version, TFEFunction2D *p1_function, 
			    TFEFunction2D *p2_function) 
{
  TFESpace2D *p1_space, *p2_space;
  double *p1_values, *p2_values;
  int p1_ndof, p2_ndof;
  TCollection *p1_coll, *p2_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3;
  double coarse_values[6];
  double fine_values[12];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3;
  int *dof;
  TAuxParam2D *aux2d;
  TFESpace2D *fesp2d[1];
  int j,k;
  
  p1_space=p1_function->GetFESpace2D();
  p2_space=p2_function->GetFESpace2D();
  
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
    
    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    for(l=0;l<3;l++)
      fine_values[l]=p1_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    for(l=0;l<3;l++)
      fine_values[3+l]=p1_values[dof[l]];
 
    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    for(l=0;l<3;l++)
      fine_values[6+l]=p1_values[dof[l]];
    
    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    for(l=0;l<3;l++)
      fine_values[9+l]=p1_values[dof[l]];
    
   if(version==0) 
     Transform_P1P2_1_2D(fine_values, coarse_values);
   else
     Transform_P1P2_2_2D(fine_values, coarse_values);

    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<6;l++)
      p2_values[dof[l]]=coarse_values[l];
  }
}


void Superconvergence_NQ1P2_2D(TFEFunction2D *q1n_function, 
			    TFEFunction2D *p2_function) 
{
  TFESpace2D *q1n_space;
  TFESpace2D *p2_space;

  double *q1n_values,*p2_values;
  int q1n_ndof,p2_ndof;
  
  TCollection *p2_coll, *q1n_coll;
  TBaseCell *coarse_cell,*child_cell0,*child_cell1;
  TBaseCell *child_cell2,*child_cell3;
  double coarse_values[6];
  double fine_values[16];
  int *coarse_GlobalNumbers,*fine_GlobalNumbers;
  int *coarse_BeginIndex,*fine_BeginIndex;
  int l,n0,n1,n2,n3;
  int *dof;
  TAuxParam2D *aux2d;
  TFESpace2D *fesp2d[1];
  int j,k;

  q1n_space=q1n_function->GetFESpace2D();
  p2_space=p2_function->GetFESpace2D();
  
  q1n_coll=q1n_space->GetCollection();
  p2_coll=p2_space->GetCollection();
  
  q1n_ndof=q1n_space->GetN_DegreesOfFreedom();
  q1n_values=q1n_function->GetValues();
  p2_ndof=p2_space->GetN_DegreesOfFreedom();
  p2_values=p2_function->GetValues();
  
  coarse_GlobalNumbers=p2_space->GetGlobalNumbers();
  coarse_BeginIndex=p2_space->GetBeginIndex();
  fine_GlobalNumbers=q1n_space->GetGlobalNumbers();
  fine_BeginIndex=q1n_space->GetBeginIndex();

  /* enumerate the cells */
 
  k=q1n_coll->GetN_Cells();
  for(j=0;j<k;j++) /* for macro cell */
    q1n_coll->GetCell(j)->SetClipBoard(j);
	
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
	    
    n0=child_cell0->GetClipBoard();
    n1=child_cell1->GetClipBoard();
    n2=child_cell2->GetClipBoard();
    n3=child_cell3->GetClipBoard();
	    
    dof=fine_GlobalNumbers+fine_BeginIndex[n0];
    for(l=0;l<5;l++)
      fine_values[l]=q1n_values[dof[l]];

    dof=fine_GlobalNumbers+fine_BeginIndex[n1];
    fine_values[5]=q1n_values[dof[0]];
    fine_values[6]=q1n_values[dof[1]];
    fine_values[7]=q1n_values[dof[3]];
    fine_values[8]=q1n_values[dof[4]];


    dof=fine_GlobalNumbers+fine_BeginIndex[n2];
    fine_values[9]=q1n_values[dof[0]];
    fine_values[10]=q1n_values[dof[1]];	     
    fine_values[11]=q1n_values[dof[3]];
    fine_values[12]=q1n_values[dof[4]];


    dof=fine_GlobalNumbers+fine_BeginIndex[n3];
    fine_values[13]=q1n_values[dof[3]];
    fine_values[14]=q1n_values[dof[0]];	
    fine_values[15]=q1n_values[dof[4]];
	    
    Transform_NQ1P2_2D(fine_values, coarse_values);    
    
    dof=coarse_GlobalNumbers+coarse_BeginIndex[j];
    for(l=0;l<6;l++)
      p2_values[dof[l]]=coarse_values[l];	        
  
  }
}
