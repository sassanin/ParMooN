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
// Raviart-Thomas element of second order on tetrahedra, 3D
// ***********************************************************************

/* for all functionals */

//tschebyscheff point
static double tscheb_point_3 = 0.866025403784439;

static double NF_N_H_RT2_3D_Xi[] = {
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  1, 1, 1,  1, 1, 1,  1, 1, 1,
  tscheb_point_3, tscheb_point_3, tscheb_point_3, 0, 0, 0, -tscheb_point_3, -tscheb_point_3, -tscheb_point_3,
  -1, -1, -1,  -1, -1, -1,  -1, -1, -1,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3
};
static double NF_N_H_RT2_3D_Eta[]  = {
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  -1, -1, -1,  -1, -1, -1,  -1, -1, -1,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  1, 1, 1,  1, 1, 1,  1, 1, 1,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3
};
static double NF_N_H_RT2_3D_Zeta[] = {
  -1, -1, -1,  -1, -1, -1,  -1, -1, -1,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3,
  -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3,
  1, 1, 1,  1, 1, 1,  1, 1, 1
};

// face 0
static double NF_N_H_RT2_3D_F0_Xi[]   = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };
static double NF_N_H_RT2_3D_F0_Eta[]  = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F0_Zeta[] = {-1, -1, -1,  -1, -1, -1,  -1, -1, -1};

// face 1                               1
static double NF_N_H_RT2_3D_F1_Xi[]   = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F1_Eta[]  = {-1, -1, -1,  -1, -1, -1,  -1, -1, -1};
static double NF_N_H_RT2_3D_F1_Zeta[] = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };

// face 2                               2
static double NF_N_H_RT2_3D_F2_Xi[]   = { 1, 1, 1,  1, 1, 1,  1, 1, 1};
static double NF_N_H_RT2_3D_F2_Eta[]  = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F2_Zeta[] = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };

 //face 3                               3
static double NF_N_H_RT2_3D_F3_Xi[]   = { tscheb_point_3, tscheb_point_3, tscheb_point_3, 0, 0, 0, -tscheb_point_3, -tscheb_point_3, -tscheb_point_3 };
static double NF_N_H_RT2_3D_F3_Eta[]  = { 1, 1, 1,  1, 1, 1,  1, 1, 1};
static double NF_N_H_RT2_3D_F3_Zeta[] = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };

// face 4                               4
static double NF_N_H_RT2_3D_F4_Xi[]   = {-1, -1, -1,  -1, -1, -1,  -1, -1, -1};
static double NF_N_H_RT2_3D_F4_Eta[]  = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };
static double NF_N_H_RT2_3D_F4_Zeta[] = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };

// face 5                               5
static double NF_N_H_RT2_3D_F5_Xi[]   = { -tscheb_point_3, -tscheb_point_3, -tscheb_point_3, 0, 0, 0, tscheb_point_3, tscheb_point_3, tscheb_point_3 };
static double NF_N_H_RT2_3D_F5_Eta[]  = { -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3, -tscheb_point_3, 0, tscheb_point_3 };
static double NF_N_H_RT2_3D_F5_Zeta[] = { 1, 1, 1,  1, 1, 1,  1, 1, 1};

static double *NF_N_H_RT2_3D_XiArray[6] = {
                        NF_N_H_RT2_3D_F0_Xi,
                        NF_N_H_RT2_3D_F1_Xi,
                        NF_N_H_RT2_3D_F2_Xi,
                        NF_N_H_RT2_3D_F3_Xi,
                        NF_N_H_RT2_3D_F4_Xi,
                        NF_N_H_RT2_3D_F5_Xi };

static double *NF_N_H_RT2_3D_EtaArray[6] = {
                        NF_N_H_RT2_3D_F0_Eta,
                        NF_N_H_RT2_3D_F1_Eta,
                        NF_N_H_RT2_3D_F2_Eta,
                        NF_N_H_RT2_3D_F3_Eta,
                        NF_N_H_RT2_3D_F4_Eta,
                        NF_N_H_RT2_3D_F5_Eta };

static double *NF_N_H_RT2_3D_ZetaArray[6] = {
                        NF_N_H_RT2_3D_F0_Zeta,
                        NF_N_H_RT2_3D_F1_Zeta,
                        NF_N_H_RT2_3D_F2_Zeta,
                        NF_N_H_RT2_3D_F3_Zeta,
                        NF_N_H_RT2_3D_F4_Zeta,
                        NF_N_H_RT2_3D_F5_Zeta };

static double NF_N_H_RT2_3D_T[] = {-100};//???
static double NF_N_H_RT2_3D_S[] = {-100};//???



void NF_N_H_RT2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
	cout << "Raviart-Thomas elements of order 2 on hexaeder: "
	     << "Nodal functionals are not fully implemented properly!" << endl;
	  for(int i=0; i<108; i++)
	    Functionals[i] = 0;
}

void NF_N_H_RT2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
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
  Functionals[6] = PointValues[6]*s;
  Functionals[7] = PointValues[7]*s;
  Functionals[8] = PointValues[8]*s;
}

static int NF_N_H_RT2_3D_N_AllFunctionals = 108;
static int NF_N_H_RT2_3D_N_PointsAll = 54;
static int NF_N_H_RT2_3D_N_FaceFunctionals[] = { 9, 9, 9, 9, 9, 9 };
static int NF_N_H_RT2_3D_N_PointsFace[] = { 9, 9, 9, 9, 9, 9 };

TNodalFunctional3D *NF_N_H_RT2_3D_Obj = new TNodalFunctional3D
        (NF_N_H_RT2_3D, NF_N_H_RT2_3D_N_AllFunctionals,
         NF_N_H_RT2_3D_N_FaceFunctionals, NF_N_H_RT2_3D_N_PointsAll,
         NF_N_H_RT2_3D_N_PointsFace,
         NF_N_H_RT2_3D_Xi, NF_N_H_RT2_3D_Eta, NF_N_H_RT2_3D_Zeta,
         NF_N_H_RT2_3D_XiArray, NF_N_H_RT2_3D_EtaArray,
         NF_N_H_RT2_3D_ZetaArray,
         NF_N_H_RT2_3D_T, NF_N_H_RT2_3D_S,
         NF_N_H_RT2_3D_EvalAll, NF_N_H_RT2_3D_EvalFace);
