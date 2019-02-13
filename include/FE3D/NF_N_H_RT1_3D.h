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
// Raviart-Thomas element of first order on hexahedra, 3D
// ***********************************************************************

/* for all functionals */

//tschebyscheff point
static double tscheb_point = 0.707106781186547;

static double NF_N_H_RT1_3D_Xi[] = {
    -tscheb_point, tscheb_point, -tscheb_point, tscheb_point,
    -tscheb_point, -tscheb_point, tscheb_point, tscheb_point,
    1, 1, 1, 1,
    tscheb_point, tscheb_point, -tscheb_point, -tscheb_point,
    -1, -1, -1, -1,
    -tscheb_point, -tscheb_point, tscheb_point, tscheb_point,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526};
static double NF_N_H_RT1_3D_Eta[]  = {
    -tscheb_point, -tscheb_point, tscheb_point, tscheb_point,
    -1, -1, -1, -1,
    -tscheb_point, -tscheb_point, tscheb_point, tscheb_point,
    1, 1, 1, 1,
    -tscheb_point, tscheb_point, -tscheb_point, tscheb_point,
    -tscheb_point, tscheb_point, -tscheb_point, tscheb_point,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.8611363115940526, -0.8611363115940526,
    -0.3399810435848563, -0.3399810435848563,
    -0.3399810435848563, -0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526};
static double NF_N_H_RT1_3D_Zeta[] = {
    -1, -1, -1, -1,
    -tscheb_point, tscheb_point, -tscheb_point, tscheb_point,
    -tscheb_point, tscheb_point, -tscheb_point, tscheb_point,
    -tscheb_point, tscheb_point, -tscheb_point, tscheb_point,
    -tscheb_point, -tscheb_point, tscheb_point, tscheb_point,
    1, 1, 1, 1,
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
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.3399810435848563,  0.3399810435848563,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526,
     0.8611363115940526,  0.8611363115940526};

// face 0
static double NF_N_H_RT1_3D_F0_Xi[]   = { -tscheb_point, tscheb_point, -tscheb_point, tscheb_point };
static double NF_N_H_RT1_3D_F0_Eta[]  = { -tscheb_point, -tscheb_point, tscheb_point, tscheb_point };
static double NF_N_H_RT1_3D_F0_Zeta[] = {-1, -1, -1, -1 };

// face 1                               1
static double NF_N_H_RT1_3D_F1_Xi[]   = { -tscheb_point, -tscheb_point, tscheb_point, tscheb_point };
static double NF_N_H_RT1_3D_F1_Eta[]  = {-1, -1, -1, -1 };
static double NF_N_H_RT1_3D_F1_Zeta[] = { -tscheb_point, tscheb_point, -tscheb_point, tscheb_point };

// face 2                               2
static double NF_N_H_RT1_3D_F2_Xi[]   = { 1, 1, 1, 1 };
static double NF_N_H_RT1_3D_F2_Eta[]  = { -tscheb_point, -tscheb_point, tscheb_point, tscheb_point };
static double NF_N_H_RT1_3D_F2_Zeta[] = { -tscheb_point, tscheb_point, -tscheb_point, tscheb_point };

 //face 3                               3
static double NF_N_H_RT1_3D_F3_Xi[]   = { tscheb_point, tscheb_point, -tscheb_point, -tscheb_point };
static double NF_N_H_RT1_3D_F3_Eta[]  = { 1, 1, 1, 1 };
static double NF_N_H_RT1_3D_F3_Zeta[] = { -tscheb_point, tscheb_point, -tscheb_point, tscheb_point };

// face 4                               4
static double NF_N_H_RT1_3D_F4_Xi[]   = {-1, -1, -1, -1 };
static double NF_N_H_RT1_3D_F4_Eta[]  = { -tscheb_point, tscheb_point, -tscheb_point, tscheb_point };
static double NF_N_H_RT1_3D_F4_Zeta[] = { -tscheb_point, -tscheb_point, tscheb_point, tscheb_point };

// face 5                               5
static double NF_N_H_RT1_3D_F5_Xi[]   = { -tscheb_point, -tscheb_point, tscheb_point, tscheb_point };
static double NF_N_H_RT1_3D_F5_Eta[]  = { -tscheb_point, tscheb_point, -tscheb_point, tscheb_point };
static double NF_N_H_RT1_3D_F5_Zeta[] = { 1, 1, 1, 1 };

static double *NF_N_H_RT1_3D_XiArray[6] = {
                        NF_N_H_RT1_3D_F0_Xi,
                        NF_N_H_RT1_3D_F1_Xi,
                        NF_N_H_RT1_3D_F2_Xi,
                        NF_N_H_RT1_3D_F3_Xi,
                        NF_N_H_RT1_3D_F4_Xi,
                        NF_N_H_RT1_3D_F5_Xi };

static double *NF_N_H_RT1_3D_EtaArray[6] = {
                        NF_N_H_RT1_3D_F0_Eta,
                        NF_N_H_RT1_3D_F1_Eta,
                        NF_N_H_RT1_3D_F2_Eta,
                        NF_N_H_RT1_3D_F3_Eta,
                        NF_N_H_RT1_3D_F4_Eta,
                        NF_N_H_RT1_3D_F5_Eta };

static double *NF_N_H_RT1_3D_ZetaArray[6] = {
                        NF_N_H_RT1_3D_F0_Zeta,
                        NF_N_H_RT1_3D_F1_Zeta,
                        NF_N_H_RT1_3D_F2_Zeta,
                        NF_N_H_RT1_3D_F3_Zeta,
                        NF_N_H_RT1_3D_F4_Zeta,
                        NF_N_H_RT1_3D_F5_Zeta };

static double NF_N_H_RT1_3D_T[] = {-100};  //??? initilize the correct value
static double NF_N_H_RT1_3D_S[] = {-100};  //???



void NF_N_H_RT1_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  //face 0
  Functionals[0]  = -PointValues[176] * 4.0;
  Functionals[1]  = -PointValues[177] * 4.0;
  Functionals[2]  = -PointValues[178] * 4.0;
  Functionals[3]  = -PointValues[179] * 4.0;
  //face 1
  Functionals[4]  = -PointValues[92] * 4.0;
  Functionals[5]  = -PointValues[93] * 4.0;
  Functionals[6]  = -PointValues[94] * 4.0;
  Functionals[7]  = -PointValues[95] * 4.0;
  //face 2
  Functionals[8]  =  PointValues[8]  * 4.0;
  Functionals[9]  =  PointValues[9]  * 4.0;
  Functionals[10] =  PointValues[10] * 4.0;
  Functionals[11] =  PointValues[11] * 4.0;
  //face 3
  Functionals[12] =  PointValues[100] * 4.0;
  Functionals[13] =  PointValues[101] * 4.0;
  Functionals[14] =  PointValues[102] * 4.0;
  Functionals[15] =  PointValues[103] * 4.0;
  //face 4
  Functionals[16] = -PointValues[16] * 4.0;
  Functionals[17] = -PointValues[17] * 4.0;
  Functionals[18] = -PointValues[18] * 4.0;
  Functionals[19] = -PointValues[19] * 4.0;
  //face 5
  Functionals[20] =  PointValues[196] * 4.0;
  Functionals[21] =  PointValues[197] * 4.0;
  Functionals[22] =  PointValues[198] * 4.0;
  Functionals[23] =  PointValues[199] * 4.0;

  //inner DOF
  int i;
  double s;
  //x-component
  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+24] * NF_D_H_P3_3D_Weights[i];
  Functionals[24] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+24] * NF_D_H_P3_3D_Array3[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[25] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+24] * NF_D_H_P3_3D_Array4[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[26] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+24] * NF_D_H_P3_3D_Array9[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[27] = s;
  //y-component
  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+112] * NF_D_H_P3_3D_Weights[i];
  Functionals[28] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+112] * NF_D_H_P3_3D_Array2[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[29] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+112] * NF_D_H_P3_3D_Array4[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[30] = s ;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+112] * NF_D_H_P3_3D_Array7[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[31] = s;
  //z-component
  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+200] * NF_D_H_P3_3D_Weights[i];
  Functionals[32] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+200] * NF_D_H_P3_3D_Array2[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[33] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+200] * NF_D_H_P3_3D_Array3[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[34] = s;

  s = 0;
  for(i=0;i<64;i++)
    s += PointValues[i+200] * NF_D_H_P3_3D_Array6[i] * NF_D_H_P3_3D_Weights[i];
  Functionals[35] = s;
}

void NF_N_H_RT1_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
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
}

static int NF_N_H_RT1_3D_N_AllFunctionals = 36;
static int NF_N_H_RT1_3D_N_PointsAll = 88;
static int NF_N_H_RT1_3D_N_FaceFunctionals[] = { 4, 4, 4, 4, 4, 4 };
static int NF_N_H_RT1_3D_N_PointsFace[] = { 4, 4, 4, 4, 4, 4 };

TNodalFunctional3D *NF_N_H_RT1_3D_Obj = new TNodalFunctional3D
        (NF_N_H_RT1_3D, NF_N_H_RT1_3D_N_AllFunctionals,
         NF_N_H_RT1_3D_N_FaceFunctionals, NF_N_H_RT1_3D_N_PointsAll,
         NF_N_H_RT1_3D_N_PointsFace,
         NF_N_H_RT1_3D_Xi, NF_N_H_RT1_3D_Eta, NF_N_H_RT1_3D_Zeta,
         NF_N_H_RT1_3D_XiArray, NF_N_H_RT1_3D_EtaArray,
         NF_N_H_RT1_3D_ZetaArray,
         NF_N_H_RT1_3D_T, NF_N_H_RT1_3D_S,
         NF_N_H_RT1_3D_EvalAll, NF_N_H_RT1_3D_EvalFace);
