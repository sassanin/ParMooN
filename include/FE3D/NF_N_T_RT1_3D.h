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
// Raviart-Thomas element of first order on tetrahedra, 3D
// ***********************************************************************

//point
static double RT1T_tp = 0.166666666666667;

/* for all functionals */
static double NF_N_T_RT1_3D_Xi[]  = {RT1T_tp, 0.5+RT1T_tp, RT1T_tp, //face 0
    RT1T_tp, RT1T_tp, 0.5+RT1T_tp, //face 1
    RT1T_tp, 1-2*RT1T_tp, RT1T_tp, //face 2
    0,0,0, //face 3
    //innere dofs
    0.138196601125010515179541316563,
    0.138196601125010515179541316563,
    0.138196601125010515179541316563,
    0.585410196624968454461376050311

};
static double NF_N_T_RT1_3D_Eta[] = {RT1T_tp, RT1T_tp, 0.5+RT1T_tp, //face 0
    0,0,0, //face 1
    1-2*RT1T_tp, RT1T_tp, RT1T_tp, //face 2
    RT1T_tp, 0.5+RT1T_tp, RT1T_tp, //face 3
    //innere dofs
    0.138196601125010515179541316563,
    0.138196601125010515179541316563,
    0.585410196624968454461376050311,
    0.138196601125010515179541316563
};
static double NF_N_T_RT1_3D_Zeta[]= {0,0,0, //face 0
    RT1T_tp, 0.5+RT1T_tp, RT1T_tp, //face 1
    RT1T_tp, RT1T_tp, 1-2*RT1T_tp, //face 2
    RT1T_tp, RT1T_tp, 0.5+RT1T_tp, //face 3
    //inner dofs
    0.138196601125010515179541316563,
    0.585410196624968454461376050311,
    0.138196601125010515179541316563,
    0.138196601125010515179541316563
};

static double NF_N_T_RT1_3D_Weights[]= {
    0.04166666666666666666666667, 0.04166666666666666666666667,
    0.04166666666666666666666667, 0.04166666666666666666666667
};

/* face 0                               0 */
static double NF_N_T_RT1_3D_F0_Xi[]   = { RT1T_tp, 0.5+RT1T_tp, RT1T_tp };
static double NF_N_T_RT1_3D_F0_Eta[]  = { RT1T_tp, RT1T_tp, 0.5+RT1T_tp };
static double NF_N_T_RT1_3D_F0_Zeta[] = { 0,0,0 };

/* face 1                               1 */
static double NF_N_T_RT1_3D_F1_Xi[]   = { RT1T_tp, RT1T_tp, 0.5+RT1T_tp };
static double NF_N_T_RT1_3D_F1_Eta[]  = { 0,0,0 };
static double NF_N_T_RT1_3D_F1_Zeta[] = { RT1T_tp, 0.5+RT1T_tp, RT1T_tp };

/* face 2                               2 */
static double NF_N_T_RT1_3D_F2_Xi[]   = { RT1T_tp, 1-2*RT1T_tp, RT1T_tp };
static double NF_N_T_RT1_3D_F2_Eta[]  = { 1-2*RT1T_tp, RT1T_tp, RT1T_tp };
static double NF_N_T_RT1_3D_F2_Zeta[] = { RT1T_tp, RT1T_tp, 1-2*RT1T_tp };


/* face 3                               3 */
static double NF_N_T_RT1_3D_F3_Xi[]   = { 0,0,0 };
static double NF_N_T_RT1_3D_F3_Eta[]  = { RT1T_tp, 0.5+RT1T_tp, RT1T_tp };
static double NF_N_T_RT1_3D_F3_Zeta[] = { RT1T_tp, RT1T_tp, 0.5+RT1T_tp };

static double *NF_N_T_RT1_3D_XiArray[4] = {
                        NF_N_T_RT1_3D_F0_Xi,
                        NF_N_T_RT1_3D_F1_Xi,
                        NF_N_T_RT1_3D_F2_Xi,
                        NF_N_T_RT1_3D_F3_Xi };

static double *NF_N_T_RT1_3D_EtaArray[4] = {
                        NF_N_T_RT1_3D_F0_Eta,
                        NF_N_T_RT1_3D_F1_Eta,
                        NF_N_T_RT1_3D_F2_Eta,
                        NF_N_T_RT1_3D_F3_Eta };

static double *NF_N_T_RT1_3D_ZetaArray[4] = {
                        NF_N_T_RT1_3D_F0_Zeta,
                        NF_N_T_RT1_3D_F1_Zeta,
                        NF_N_T_RT1_3D_F2_Zeta,
                        NF_N_T_RT1_3D_F3_Zeta };

static double NF_N_T_RT1_3D_T[1] = {};// ???
static double NF_N_T_RT1_3D_S[1] = {};// ???

void NF_N_T_RT1_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  // PointValues[4*i + j] means i-th component (i=0 for x, i=1 for y, i=2 for z)
  // at j-th evaluation point (see NF_N_T_RT1_3D_Xi, ...Eta, ...Zeta)
  //face 0
  Functionals[0]  = -PointValues[32];
  Functionals[1]  = -PointValues[33];
  Functionals[2]  = -PointValues[34];
  //face 1
  Functionals[3]  = -PointValues[19];
  Functionals[4]  = -PointValues[20];
  Functionals[5]  = -PointValues[21];
  //face 2
  Functionals[6]  =  PointValues[6]+PointValues[22]+PointValues[38];
  Functionals[7]  =  PointValues[7]+PointValues[23]+PointValues[39];
  Functionals[8]  =  PointValues[8]+PointValues[24]+PointValues[40];
  //face 3
  Functionals[9]  = -PointValues[9];
  Functionals[10] = -PointValues[10];
  Functionals[11] = -PointValues[11];

  //inner dofs
  int i;
  double s;

  //x-component
  s = 0;
  for(i=0;i<4;i++)
    s += PointValues[i+12] * NF_N_T_RT1_3D_Weights[i];
  Functionals[12] = s;

  //y-component
  s = 0;
  for(i=0;i<4;i++)
    s += PointValues[i+28] * NF_N_T_RT1_3D_Weights[i];
  Functionals[13] = s;

  //z-component
  s = 0;
  for(i=0;i<4;i++)
    s += PointValues[i+44] * NF_N_T_RT1_3D_Weights[i];
  Functionals[14] = s;
}

void NF_N_T_RT1_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
                           double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  #ifdef __3D__
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 3, length == {3,3,3,3}
  Cell->GetVertex(faceVertex[3*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[3*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[3*face + 2])->GetCoords(x2,y2,z2);
  #endif
  // compute measure of this face
  s = sqrt( POW((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + POW((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + POW((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  Functionals[0] = PointValues[0]*s;
  Functionals[1] = PointValues[1]*s;
  Functionals[2] = PointValues[2]*s;
}

static int NF_N_T_RT1_3D_N_AllFunctionals = 15;
static int NF_N_T_RT1_3D_N_PointsAll = 16;
static int NF_N_T_RT1_3D_N_FaceFunctionals[] = { 3, 3, 3, 3 };
static int NF_N_T_RT1_3D_N_PointsFace[] = { 3, 3, 3, 3 };

TNodalFunctional3D *NF_N_T_RT1_3D_Obj = new TNodalFunctional3D
        (NF_N_T_RT1_3D, NF_N_T_RT1_3D_N_AllFunctionals,
         NF_N_T_RT1_3D_N_FaceFunctionals, NF_N_T_RT1_3D_N_PointsAll,
         NF_N_T_RT1_3D_N_PointsFace,
         NF_N_T_RT1_3D_Xi, NF_N_T_RT1_3D_Eta, NF_N_T_RT1_3D_Zeta,
         NF_N_T_RT1_3D_XiArray, NF_N_T_RT1_3D_EtaArray,
         NF_N_T_RT1_3D_ZetaArray,
         NF_N_T_RT1_3D_T, NF_N_T_RT1_3D_S,
         NF_N_T_RT1_3D_EvalAll, NF_N_T_RT1_3D_EvalFace);
