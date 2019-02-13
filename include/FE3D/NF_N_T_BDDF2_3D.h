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
// Brezzi-Douglas-Duran-Fortin element of second order on tetrahedra, 3D
// ***********************************************************************

static double NF_N_T_BDDF2_3D_Xi[]  = {
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    0,0,0,0,0,0,
    //inner dof
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0,
    0.333333333333333333333333333333333,
    0.25,
    0.909090909090909090909090909090909e-1,
    0.909090909090909090909090909090909e-1,
    0.727272727272727272727272727272727,
    0.909090909090909090909090909090909e-1,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697
};
static double NF_N_T_BDDF2_3D_Eta[] = {
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    0,0,0,0,0,0,
    0.6, 0.4, 0.2, 0.4, 0.2, 0.2,
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    //inner dof
    0.333333333333333333333333333333333,
    0,
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0.25,
    0.909090909090909090909090909090909e-1,
    0.727272727272727272727272727272727,
    0.909090909090909090909090909090909e-1,
    0.909090909090909090909090909090909e-1,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1,
    0.433449846426335701760119535773697,
    0.665501535736642982398804642263025e-1
};
static double NF_N_T_BDDF2_3D_Zeta[]= {
    0,0,0,0,0,0,
    0.2, 0.4, 0.6, 0.2, 0.4, 0.2,
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    0.2, 0.2, 0.2, 0.4, 0.4, 0.6,
    //inner dof
    0,
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0.333333333333333333333333333333333,
    0.25,
    0.727272727272727272727272727272727,
    0.90909090909090909090909090909091e-1,
    0.909090909090909090909090909090911e-1,
    0.909090909090909090909090909090911e-1,
    0.433449846426335701760119535773699,
    0.433449846426335701760119535773699,
    0.66550153573664298239880464226304e-1,
    0.433449846426335701760119535773699,
    0.66550153573664298239880464226304e-1,
    0.665501535736642982398804642263035e-1
};
// points are x-x^2-2xy-xz
static double NF_N_T_BDDF2_3D_XiEta_Xi[]= {
    0,
    0.111111111111111,
    0,
    -0.111111111111111,
    0,
    0,
    -0.057851239669421,
    0,
    0.057851239669421,
    -0.024417230905476,
    0,
    -0.159032615520860,
    0.024417230905476,
    0,
    0.159032615520860
};
// points are -y+y^2+2xy+yz
static double NF_N_T_BDDF2_3D_XiEta_Eta[]= {
    0,
    0,
    -0.111111111111111,
    0.111111111111111,
    0,
    0,
    0,
    0.057851239669421,
    -0.057851239669421,
    0,
    0.024417230905476,
    0.159032615520860,
    -0.024417230905476,
    -0.159032615520860,
    0
};
// points are x-x^2-xy-2xz
static double NF_N_T_BDDF2_3D_XiZeta_Xi[]= {
    0.111111111111111,
    0,
    0,
    -0.111111111111111,
    0,
    -0.057851239669421,
    0,
    0,
    0.057851239669421,
    -0.024417230905476,
    -0.159032615520860,
    0,
    0,
    0.024417230905476,
    0.159032615520860
};
// points are -z+z^2+2xz+yz
static double NF_N_T_BDDF2_3D_XiZeta_Zeta[]= {
    0,
    0,
    -0.111111111111111,
    0.111111111111111,
    0,
    0,
    0,
    0.057851239669421,
    -0.057851239669421,
    0,
    0.159032615520860,
    0.024417230905476,
    -0.159032615520860,
    -0.024417230905476,
    0
};
// points are y-y^2-xy-2yz
static double NF_N_T_BDDF2_3D_EtaZeta_Eta[]= {
    0.111111111111111,
    0,
    0,
    -0.111111111111111,
    0,
    -0.057851239669421,
    0,
    0,
    0.057851239669421,
    -0.159032615520860,
    -0.024417230905476,
    -0.000000000000000,
    0,
    0.159032615520860,
    0.024417230905476
};
// points are -z+z^2+xz+2yz
static double NF_N_T_BDDF2_3D_EtaZeta_Zeta[]= {
    0,
    -0.111111111111111,
    0,
    0.111111111111111,
    0,
    0,
    0.057851239669421,
    0,
    -0.057851239669421,
    0.159032615520860,
    0,
    0.024417230905476,
    -0.159032615520860,
    0,
    -0.024417230905476
};


static double NF_N_T_BDDF2_3D_Weights[]= {
    0.602678571428571428571428571428571e-2,
    0.602678571428571428571428571428571e-2,
    0.602678571428571428571428571428571e-2,
    0.602678571428571428571428571428571e-2,
    0.302836780970891758063769725577305e-1,
    0.116452490860289694108936091443380e-1,
    0.116452490860289694108936091443380e-1,
    0.116452490860289694108936091443380e-1,
    0.116452490860289694108936091443380e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1,
    0.109491415613864593456430191124068e-1};

/* face 0                               0 */
static double NF_N_T_BDDF2_3D_F0_Xi[]   = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};
static double NF_N_T_BDDF2_3D_F0_Eta[]  = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};
static double NF_N_T_BDDF2_3D_F0_Zeta[] = { 0,0,0,0,0,0 };

/* face 1                               1 */
static double NF_N_T_BDDF2_3D_F1_Xi[]   = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};
static double NF_N_T_BDDF2_3D_F1_Eta[]  = { 0,0,0,0,0,0 };
static double NF_N_T_BDDF2_3D_F1_Zeta[] = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};

/* face 2                               2 */
static double NF_N_T_BDDF2_3D_F2_Xi[]   = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};
static double NF_N_T_BDDF2_3D_F2_Eta[]  = {0.6, 0.4, 0.2, 0.4, 0.2, 0.2};
static double NF_N_T_BDDF2_3D_F2_Zeta[] = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};

/* face 3                               3 */
static double NF_N_T_BDDF2_3D_F3_Xi[]   = { 0,0,0,0,0,0 };
static double NF_N_T_BDDF2_3D_F3_Eta[]  = {0.2, 0.4, 0.6, 0.2, 0.4, 0.2};
static double NF_N_T_BDDF2_3D_F3_Zeta[] = {0.2, 0.2, 0.2, 0.4, 0.4, 0.6};

static double *NF_N_T_BDDF2_3D_XiArray[4] = {
                        NF_N_T_BDDF2_3D_F0_Xi,
                        NF_N_T_BDDF2_3D_F1_Xi,
                        NF_N_T_BDDF2_3D_F2_Xi,
                        NF_N_T_BDDF2_3D_F3_Xi };

static double *NF_N_T_BDDF2_3D_EtaArray[4] = {
                        NF_N_T_BDDF2_3D_F0_Eta,
                        NF_N_T_BDDF2_3D_F1_Eta,
                        NF_N_T_BDDF2_3D_F2_Eta,
                        NF_N_T_BDDF2_3D_F3_Eta };

static double *NF_N_T_BDDF2_3D_ZetaArray[4] = {
                        NF_N_T_BDDF2_3D_F0_Zeta,
                        NF_N_T_BDDF2_3D_F1_Zeta,
                        NF_N_T_BDDF2_3D_F2_Zeta,
                        NF_N_T_BDDF2_3D_F3_Zeta };

static double NF_N_T_BDDF2_3D_T[] = {-100};// ???
static double NF_N_T_BDDF2_3D_S[] = {-100};// ???

void NF_N_T_BDDF2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  // PointValues[4*i + j] means i-th component (i=0 for x, i=1 for y, i=2 for z)
  // at j-th evaluation point (see NF_N_T_BDDF2_3D_Xi, ...Eta, ...Zeta)
  //face 0
  Functionals[0]  = -PointValues[78];
  Functionals[1]  = -PointValues[79];
  Functionals[2]  = -PointValues[80];
  Functionals[3]  = -PointValues[81];
  Functionals[4]  = -PointValues[82];
  Functionals[5]  = -PointValues[83];
  //face 1
  Functionals[6]  = -PointValues[45];
  Functionals[7]  = -PointValues[46];
  Functionals[8]  = -PointValues[47];
  Functionals[9]  = -PointValues[48];
  Functionals[10] = -PointValues[49];
  Functionals[11] = -PointValues[50];
  //face 2
  Functionals[12] =  PointValues[12]+PointValues[51]+PointValues[90];
  Functionals[13] =  PointValues[13]+PointValues[52]+PointValues[91];
  Functionals[14] =  PointValues[14]+PointValues[53]+PointValues[92];
  Functionals[15] =  PointValues[15]+PointValues[54]+PointValues[93];
  Functionals[16] =  PointValues[16]+PointValues[55]+PointValues[94];
  Functionals[17] =  PointValues[17]+PointValues[56]+PointValues[95];
  //face 3
  Functionals[18] = -PointValues[18];
  Functionals[19] = -PointValues[19];
  Functionals[20] = -PointValues[20];
  Functionals[21] = -PointValues[21];
  Functionals[22] = -PointValues[22];
  Functionals[23] = -PointValues[23];

  //inner dofs
  int i;
  double s;

  //int_T u \cdot \grad p
  //x-component
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+24] * NF_N_T_BDDF2_3D_Weights[i];
  Functionals[24] = s;
  //y-component
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+63] * NF_N_T_BDDF2_3D_Weights[i];
  Functionals[25] = s;
  //z-component
  s = 0;
  for(i=0;i<15;i++)
    s += PointValues[i+102] * NF_N_T_BDDF2_3D_Weights[i];
  Functionals[26] = s;

  //int_T uv and v a subspace of P_2(T)
  s = 0;
  for(i=0;i<15;i++){
      //x-component
      s += PointValues[i+24] * NF_N_T_BDDF2_3D_XiEta_Xi[i] * NF_N_T_BDDF2_3D_Weights[i];
      //y-component
      s += PointValues[i+63] * NF_N_T_BDDF2_3D_XiEta_Eta[i] * NF_N_T_BDDF2_3D_Weights[i];
  }
  Functionals[27] = s;

  s = 0;
  for(i=0;i<15;i++){
      //x-component
      s += PointValues[i+24] * NF_N_T_BDDF2_3D_XiZeta_Xi[i] * NF_N_T_BDDF2_3D_Weights[i];
      //z-component
      s += PointValues[i+102] * NF_N_T_BDDF2_3D_XiZeta_Zeta[i] * NF_N_T_BDDF2_3D_Weights[i];
  }
  Functionals[28] = s;

  s = 0;
  for(i=0;i<15;i++){
      //y-component
      s += PointValues[i+63] * NF_N_T_BDDF2_3D_EtaZeta_Eta[i] * NF_N_T_BDDF2_3D_Weights[i];
      //z-component
      s += PointValues[i+102] * NF_N_T_BDDF2_3D_EtaZeta_Zeta[i] * NF_N_T_BDDF2_3D_Weights[i];
  }
  Functionals[29] = s;

}

void NF_N_T_BDDF2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
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
  for(int i=0; i<6; i++)
    Functionals[i] = PointValues[i]*s;
}

static int NF_N_T_BDDF2_3D_N_AllFunctionals = 30;
static int NF_N_T_BDDF2_3D_N_PointsAll = 39;
static int NF_N_T_BDDF2_3D_N_FaceFunctionals[] = { 6, 6, 6, 6 };
static int NF_N_T_BDDF2_3D_N_PointsFace[] = { 6, 6, 6, 6 };

TNodalFunctional3D *NF_N_T_BDDF2_3D_Obj = new TNodalFunctional3D
        (NF_N_T_BDDF2_3D, NF_N_T_BDDF2_3D_N_AllFunctionals,
         NF_N_T_BDDF2_3D_N_FaceFunctionals, NF_N_T_BDDF2_3D_N_PointsAll,
         NF_N_T_BDDF2_3D_N_PointsFace,
         NF_N_T_BDDF2_3D_Xi, NF_N_T_BDDF2_3D_Eta, NF_N_T_BDDF2_3D_Zeta,
         NF_N_T_BDDF2_3D_XiArray, NF_N_T_BDDF2_3D_EtaArray,
         NF_N_T_BDDF2_3D_ZetaArray,
         NF_N_T_BDDF2_3D_T, NF_N_T_BDDF2_3D_S,
         NF_N_T_BDDF2_3D_EvalAll, NF_N_T_BDDF2_3D_EvalFace);
