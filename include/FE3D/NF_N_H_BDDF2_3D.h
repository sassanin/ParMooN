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
// Brezzi-Douglas-Duran-Fortin element of second order on hexahedra, 3D
// ***********************************************************************

//tschebyscheff point
static double BDDF2_tp = 0.866025403784439;

static double NF_N_H_BDDF2_3D_Xi[] = {-BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp,
    -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp,
    1,1,1,1,1,1,
    -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp,
    -1, -1, -1, -1, -1, -1,
    -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp,
    // innere dofs
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
    0.3399810435848563, 0.8611363115940526
};
static double NF_N_H_BDDF2_3D_Eta[]  = {-BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp,
    -1, -1, -1, -1, -1, -1,
    -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp,
    1,1,1,1,1,1,
    -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp,
    -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp,
    //innere dofs
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526
};
static double NF_N_H_BDDF2_3D_Zeta[] = {-1, -1, -1, -1, -1, -1,
    -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp,
    -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp,
    -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp,
    -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp,
    1,1,1,1,1,1,
    //innere dofs
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.3399810435848563, 0.3399810435848563,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526,
    0.8611363115940526, 0.8611363115940526
};

// face 0
static double NF_N_H_BDDF2_3D_F0_Xi[]   = { -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp };
static double NF_N_H_BDDF2_3D_F0_Eta[]  = { -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp };
static double NF_N_H_BDDF2_3D_F0_Zeta[] = { -1, -1, -1, -1, -1, -1};

// face 1                               1
static double NF_N_H_BDDF2_3D_F1_Xi[]   = { -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp };
static double NF_N_H_BDDF2_3D_F1_Eta[]  = { -1, -1, -1, -1, -1, -1};
static double NF_N_H_BDDF2_3D_F1_Zeta[] = { -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp };

// face 2                               2
static double NF_N_H_BDDF2_3D_F2_Xi[]   = { 1,1,1,1,1,1 };
static double NF_N_H_BDDF2_3D_F2_Eta[]  = { -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp };
static double NF_N_H_BDDF2_3D_F2_Zeta[] = { -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp };

 //face 3                               3
static double NF_N_H_BDDF2_3D_F3_Xi[]   = { -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp };
static double NF_N_H_BDDF2_3D_F3_Eta[]  = { 1,1,1,1,1,1 };
static double NF_N_H_BDDF2_3D_F3_Zeta[] = { -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp };

// face 4                               4
static double NF_N_H_BDDF2_3D_F4_Xi[]   = { -1, -1, -1, -1, -1, -1};
static double NF_N_H_BDDF2_3D_F4_Eta[]  = { -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp };
static double NF_N_H_BDDF2_3D_F4_Zeta[] = { -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp };

// face 5                               5
static double NF_N_H_BDDF2_3D_F5_Xi[]   = { -BDDF2_tp,-BDDF2_tp,0,-BDDF2_tp,0,BDDF2_tp };
static double NF_N_H_BDDF2_3D_F5_Eta[]  = { -BDDF2_tp,0,-BDDF2_tp,BDDF2_tp,0,-BDDF2_tp };
static double NF_N_H_BDDF2_3D_F5_Zeta[] = { 1,1,1,1,1,1 };

static double *NF_N_H_BDDF2_3D_XiArray[6] = {
                        NF_N_H_BDDF2_3D_F0_Xi,
                        NF_N_H_BDDF2_3D_F1_Xi,
                        NF_N_H_BDDF2_3D_F2_Xi,
                        NF_N_H_BDDF2_3D_F3_Xi,
                        NF_N_H_BDDF2_3D_F4_Xi,
                        NF_N_H_BDDF2_3D_F5_Xi };

static double *NF_N_H_BDDF2_3D_EtaArray[6] = {
                        NF_N_H_BDDF2_3D_F0_Eta,
                        NF_N_H_BDDF2_3D_F1_Eta,
                        NF_N_H_BDDF2_3D_F2_Eta,
                        NF_N_H_BDDF2_3D_F3_Eta,
                        NF_N_H_BDDF2_3D_F4_Eta,
                        NF_N_H_BDDF2_3D_F5_Eta };

static double *NF_N_H_BDDF2_3D_ZetaArray[6] = {
                        NF_N_H_BDDF2_3D_F0_Zeta,
                        NF_N_H_BDDF2_3D_F1_Zeta,
                        NF_N_H_BDDF2_3D_F2_Zeta,
                        NF_N_H_BDDF2_3D_F3_Zeta,
                        NF_N_H_BDDF2_3D_F4_Zeta,
                        NF_N_H_BDDF2_3D_F5_Zeta };

static double NF_N_H_BDDF2_3D_T[] = {-100};//???
static double NF_N_H_BDDF2_3D_S[] = {-100};//???



void NF_N_H_BDDF2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  //face 0
  Functionals[0]  = -PointValues[200] * 4.0;
  Functionals[1]  = -PointValues[201] * 4.0;
  Functionals[2]  = -PointValues[202] * 4.0;
  Functionals[3]  = -PointValues[203] * 4.0;
  Functionals[4]  = -PointValues[204] * 4.0;
  Functionals[5]  = -PointValues[205] * 4.0;
  //face 1
  Functionals[6]  = -PointValues[106] * 4.0;
  Functionals[7]  = -PointValues[107] * 4.0;
  Functionals[8]  = -PointValues[108] * 4.0;
  Functionals[9]  = -PointValues[109] * 4.0;
  Functionals[10] = -PointValues[110] * 4.0;
  Functionals[11] = -PointValues[111] * 4.0;
  //face 2
  Functionals[12] =  PointValues[12]  * 4.0;
  Functionals[13] =  PointValues[13]  * 4.0;
  Functionals[14] =  PointValues[14]  * 4.0;
  Functionals[15] =  PointValues[15]  * 4.0;
  Functionals[16] =  PointValues[16]  * 4.0;
  Functionals[17] =  PointValues[17]  * 4.0;
  //face 3
  Functionals[18] =  PointValues[118] * 4.0;
  Functionals[19] =  PointValues[119] * 4.0;
  Functionals[20] =  PointValues[120] * 4.0;
  Functionals[21] =  PointValues[121] * 4.0;
  Functionals[22] =  PointValues[122] * 4.0;
  Functionals[23] =  PointValues[123] * 4.0;
  //face 4
  Functionals[24] = -PointValues[24]  * 4.0;
  Functionals[25] = -PointValues[25]  * 4.0;
  Functionals[26] = -PointValues[26]  * 4.0;
  Functionals[27] = -PointValues[27]  * 4.0;
  Functionals[28] = -PointValues[28]  * 4.0;
  Functionals[29] = -PointValues[29]  * 4.0;
  //face 5
  Functionals[30] =  PointValues[230] * 4.0;
  Functionals[31] =  PointValues[231] * 4.0;
  Functionals[32] =  PointValues[232] * 4.0;
  Functionals[33] =  PointValues[233] * 4.0;
  Functionals[34] =  PointValues[234] * 4.0;
  Functionals[35] =  PointValues[235] * 4.0;
  //inner dofs
  int i;
  double s;

  //x-component
  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+36] * NF_D_H_P3_3D_Weights[i];
  Functionals[36] = s;

  //y-component
  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+136] * NF_D_H_P3_3D_Weights[i];
  Functionals[37] = s;

  //z-component
  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+236] * NF_D_H_P3_3D_Weights[i];
  Functionals[38] = s;
}

void NF_N_H_BDDF2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
                           double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  #ifdef __3D__
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 4, length == {4,4,4,4}
  Cell->GetVertex(faceVertex[4*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[4*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[4*face + 2])->GetCoords(x2,y2,z2);
  #endif
  // compute measure of this face
  s = sqrt( POW((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + POW((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + POW((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  Functionals[0] = PointValues[0]*s;
  Functionals[1] = PointValues[1]*s;
  Functionals[2] = PointValues[2]*s;
  Functionals[3] = PointValues[3]*s;
  Functionals[4] = PointValues[4]*s;
  Functionals[5] = PointValues[5]*s;
}

static int NF_N_H_BDDF2_3D_N_AllFunctionals = 39;
static int NF_N_H_BDDF2_3D_N_PointsAll = 100;
static int NF_N_H_BDDF2_3D_N_FaceFunctionals[] = { 6, 6, 6, 6, 6, 6 };
static int NF_N_H_BDDF2_3D_N_PointsFace[] = { 6, 6, 6, 6, 6, 6 };

TNodalFunctional3D *NF_N_H_BDDF2_3D_Obj = new TNodalFunctional3D
        (NF_N_H_BDDF2_3D, NF_N_H_BDDF2_3D_N_AllFunctionals,
         NF_N_H_BDDF2_3D_N_FaceFunctionals, NF_N_H_BDDF2_3D_N_PointsAll,
         NF_N_H_BDDF2_3D_N_PointsFace,
         NF_N_H_BDDF2_3D_Xi, NF_N_H_BDDF2_3D_Eta, NF_N_H_BDDF2_3D_Zeta,
         NF_N_H_BDDF2_3D_XiArray, NF_N_H_BDDF2_3D_EtaArray,
         NF_N_H_BDDF2_3D_ZetaArray,
         NF_N_H_BDDF2_3D_T, NF_N_H_BDDF2_3D_S,
         NF_N_H_BDDF2_3D_EvalAll, NF_N_H_BDDF2_3D_EvalFace);
