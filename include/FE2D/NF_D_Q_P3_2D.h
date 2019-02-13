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
   
// xi
static double NF_D_Q_P3_2D_Weight1[16]={
-0.8611363115940526, -0.3399810435848563,
 0.3399810435848563,  0.8611363115940526,
-0.8611363115940526, -0.3399810435848563,
 0.3399810435848563,  0.8611363115940526,
-0.8611363115940526, -0.3399810435848563,
 0.3399810435848563,  0.8611363115940526,
-0.8611363115940526, -0.3399810435848563,
 0.3399810435848563,  0.8611363115940526
                                };
// eta
static double NF_D_Q_P3_2D_Weight2[16]={
-0.8611363115940526, -0.8611363115940526,
-0.8611363115940526, -0.8611363115940526,
-0.3399810435848563, -0.3399810435848563,
-0.3399810435848563, -0.3399810435848563,
 0.3399810435848563,  0.3399810435848563,
 0.3399810435848563,  0.3399810435848563,
 0.8611363115940526,  0.8611363115940526,
 0.8611363115940526,  0.8611363115940526
                                 };
// 3*xi**2-1
static double NF_D_Q_P3_2D_Weight3[16]={
 0.1224667241437428E1, -0.6532386700088562,
-0.6532386700088562,  0.1224667241437428E1,
 0.1224667241437428E1, -0.6532386700088562,
-0.6532386700088562,  0.1224667241437428E1,
 0.1224667241437428E1, -0.6532386700088562,
-0.6532386700088562,  0.1224667241437428E1,
 0.1224667241437428E1, -0.6532386700088562,
-0.6532386700088562,  0.1224667241437428E1
                                    };
// xi*eta
static double NF_D_Q_P3_2D_Weight4[16]={
 0.7415557471458092,  0.29277002188456,
-0.29277002188456, -0.7415557471458092,
 0.29277002188456,  0.1155871099970479,
-0.1155871099970479, -0.29277002188456,
-0.29277002188456, -0.1155871099970479,
 0.1155871099970479,  0.29277002188456,
-0.7415557471458092, -0.29277002188456,
 0.29277002188456,  0.7415557471458092
                                    };
// 3*eta**2-1
static double NF_D_Q_P3_2D_Weight5[16]={
 0.1224667241437428E1,  0.1224667241437428E1,
 0.1224667241437428E1,  0.1224667241437428E1,
-0.6532386700088562, -0.6532386700088562,
-0.6532386700088562, -0.6532386700088562,
-0.6532386700088562, -0.6532386700088562,
-0.6532386700088562, -0.6532386700088562,
 0.1224667241437428E1,  0.1224667241437428E1,
 0.1224667241437428E1,  0.1224667241437428E1
                                     };
// 5*xi**3-3*xi
static double NF_D_Q_P3_2D_Weight6[16]={
-0.6094939699104124,  0.8234559993457991,
-0.8234559993457991,  0.6094939699104124,
-0.6094939699104124,  0.8234559993457991,
-0.8234559993457991,  0.6094939699104124,
-0.6094939699104124,  0.8234559993457991,
-0.8234559993457991,  0.6094939699104124,
-0.6094939699104124,  0.8234559993457991,
-0.8234559993457991,  0.6094939699104124
                                     };
// 3*xi**2*eta-eta
static double NF_D_Q_P3_2D_Weight7[16]={
-0.1054605431221489E1,  0.5625275388820309,
 0.5625275388820309, -0.1054605431221489E1,
-0.4163636467880838,  0.2220887647395945,
 0.2220887647395945, -0.4163636467880838,
 0.4163636467880838, -0.2220887647395945,
-0.2220887647395945,  0.4163636467880838,
 0.1054605431221489E1, -0.5625275388820309,
-0.5625275388820309,  0.1054605431221489E1
                                     };
// 3*eta**2*xi-xi
static double NF_D_Q_P3_2D_Weight8[16]={
-0.1054605431221489E1, -0.4163636467880838,
 0.4163636467880838,  0.1054605431221489E1,
 0.5625275388820309,  0.2220887647395945,
-0.2220887647395945, -0.5625275388820309,
 0.5625275388820309,  0.2220887647395945,
-0.2220887647395945, -0.5625275388820309,
-0.1054605431221489E1, -0.4163636467880838,
 0.4163636467880838,  0.1054605431221489E1
                                     };
// 5*eta**3-3*eta
static double NF_D_Q_P3_2D_Weight9[16]={
-0.6094939699104124, -0.6094939699104124,
-0.6094939699104124, -0.6094939699104124,
 0.8234559993457991,  0.8234559993457991,
 0.8234559993457991,  0.8234559993457991,
-0.8234559993457991, -0.8234559993457991,
-0.8234559993457991, -0.8234559993457991,
 0.6094939699104124,  0.6094939699104124,
 0.6094939699104124,  0.6094939699104124
                                     };

static double *NF_D_Q_P3_2D_T = NULL;

void NF_D_Q_P3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  static double weights[16]={ 
                0.12100299328560200552, 0.22685185185185185185,
                0.22685185185185185185, 0.12100299328560200552,
                0.22685185185185185185, 0.42529330301069429078,
                0.42529330301069429078, 0.22685185185185185185,
                0.22685185185185185185, 0.42529330301069429078,
                0.42529330301069429078, 0.22685185185185185185,
                0.12100299328560200552, 0.22685185185185185185,
                0.22685185185185185185, 0.12100299328560200552
              };
  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2]
                   +weights[3]*PointValues[3]
                   +weights[4]*PointValues[4]
                   +weights[5]*PointValues[5]
                   +weights[6]*PointValues[6]
                   +weights[7]*PointValues[7]
                   +weights[8]*PointValues[8]
                   +weights[9]*PointValues[9]
                   +weights[10]*PointValues[10]
                   +weights[11]*PointValues[11]
                   +weights[12]*PointValues[12]
                   +weights[13]*PointValues[13]
                   +weights[14]*PointValues[14]
                   +weights[15]*PointValues[15];
  Functionals[1] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight1[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight1[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight1[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight1[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight1[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight1[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight1[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight1[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight1[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight1[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight1[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight1[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight1[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight1[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight1[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight1[15];
  Functionals[2] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight2[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight2[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight2[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight2[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight2[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight2[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight2[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight2[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight2[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight2[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight2[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight2[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight2[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight2[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight2[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight2[15];
  Functionals[3] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight3[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight3[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight3[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight3[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight3[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight3[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight3[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight3[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight3[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight3[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight3[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight3[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight3[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight3[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight3[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight3[15];
  Functionals[4] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight4[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight4[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight4[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight4[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight4[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight4[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight4[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight4[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight4[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight4[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight4[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight4[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight4[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight4[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight4[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight4[15];
  Functionals[5] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight5[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight5[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight5[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight5[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight5[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight5[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight5[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight5[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight5[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight5[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight5[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight5[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight5[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight5[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight5[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight5[15];
  Functionals[6] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight6[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight6[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight6[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight6[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight6[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight6[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight6[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight6[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight6[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight6[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight6[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight6[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight6[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight6[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight6[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight6[15];
  Functionals[7] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight7[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight7[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight7[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight7[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight7[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight7[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight7[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight7[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight7[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight7[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight7[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight7[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight7[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight7[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight7[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight7[15];
  Functionals[8] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight8[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight8[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight8[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight8[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight8[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight8[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight8[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight8[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight8[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight8[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight8[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight8[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight8[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight8[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight8[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight8[15];
  Functionals[9] =  weights[0]*PointValues[0]*NF_D_Q_P3_2D_Weight9[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P3_2D_Weight9[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P3_2D_Weight9[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P3_2D_Weight9[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P3_2D_Weight9[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P3_2D_Weight9[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P3_2D_Weight9[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P3_2D_Weight9[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P3_2D_Weight9[8]
                   +weights[9]*PointValues[9]*NF_D_Q_P3_2D_Weight9[9]
                   +weights[10]*PointValues[10]*NF_D_Q_P3_2D_Weight9[10]
                   +weights[11]*PointValues[11]*NF_D_Q_P3_2D_Weight9[11]
                   +weights[12]*PointValues[12]*NF_D_Q_P3_2D_Weight9[12]
                   +weights[13]*PointValues[13]*NF_D_Q_P3_2D_Weight9[13]
                   +weights[14]*PointValues[14]*NF_D_Q_P3_2D_Weight9[14]
                   +weights[15]*PointValues[15]*NF_D_Q_P3_2D_Weight9[15];

  Functionals[0] *= 0.25;
  Functionals[1] *= 0.25;
  Functionals[2] *= 0.25;
  Functionals[3] *= 0.25;
  Functionals[4] *= 0.25;
  Functionals[5] *= 0.25;
  Functionals[6] *= 0.25;
  Functionals[7] *= 0.25;
  Functionals[8] *= 0.25;
  Functionals[9] *= 0.25;
}

void NF_D_Q_P3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_D_Q_P3_2D_Obj = new TNodalFunctional2D
        (NF_D_Q_P3_2D, 10, 0, 16, 0, NF_D_Q_P3_2D_Weight1, NF_D_Q_P3_2D_Weight2,
         NF_D_Q_P3_2D_T, NF_D_Q_P3_2D_EvalAll, NULL);
