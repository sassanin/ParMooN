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
   
// =======================================================================
// @(#)QuadFormulaHexa.C        1.2 05/04/99
//
// Class:      TQuadFormulaHexa
// Superclass: TQuadFormula3D
//
// Purpose:    quadrature formula for a 3D integral on unit hexahedron
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#include <QuadFormulaHexa.h>

TQuadFormulaHexa::TQuadFormulaHexa() : TQuadFormula3D()
{
}

TQuadFormulaHexa::TQuadFormulaHexa(int n_points, double* weights, 
                     double* xi, double* eta, double* zeta, int acc)
  : TQuadFormula3D(n_points, weights, xi, eta, zeta, acc)
{
}

void TQuadFormulaHexa::Vertex()
{
  if(Weights!=NULL) return;
  double w[8]={ 1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0 };
  double x[8]={ -1.0,  1.0, -1.0,  1.0,
                -1.0,  1.0, -1.0,  1.0 };
  double e[8]={ -1.0, -1.0,  1.0,  1.0,
                -1.0, -1.0,  1.0,  1.0 };
  double z[8]={ -1.0, -1.0, -1.0, -1.0,
                 1.0,  1.0,  1.0,  1.0 };

  InitObject(8,w,x,e,z, 1);
}

void TQuadFormulaHexa::Gauss2()
{
  if(Weights!=NULL) return;
  double w[8]={ 1.0, 1.0, 1.0, 1.0, 
                1.0, 1.0, 1.0, 1.0 };
  double x[8]={-0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489};
  double e[8]={-0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489};
  double z[8]={-0.5773502691896257645091489, -0.5773502691896257645091489,
               -0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489};

  InitObject(8,w,x,e,z, 3);
}

void TQuadFormulaHexa::Gauss3()
{
  if(Weights!=NULL) return;

  double oneweight[3]={ 0.5555555555555555555555558,
                        0.8888888888888888888888888, 
                        0.5555555555555555555555558 };
  double onexi[3]={ -0.7745966692414833770358530, 
                     0,
                     0.7745966692414833770358530 };

  double w[27], x[27], e[27], z[27];

  int l;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
      {
        l=9*i+3*j+k;
        w[l]=oneweight[i]*oneweight[j]*oneweight[k];
        x[l]=onexi[k];
        e[l]=onexi[j];
        z[l]=onexi[i];
      }

  InitObject(27,w,x,e,z, 5);
}

void TQuadFormulaHexa::Gauss4()
{
  if(Weights!=NULL) return;

  double oneweight[4]={ 0.652145154862546142626936051,
                        0.347854845137453857373063949,
                        0.347854845137453857373063949,
                        0.652145154862546142626936051 };

  double onexi[4]={  0.33998104358485626480266576,
                     0.86113631159405257522394649,
                    -0.86113631159405257522394649,
                    -0.33998104358485626480266576 };

  double w[64], x[64], e[64], z[64];

  int l;
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      for(int k=0;k<4;k++)
      {
        l=16*i+4*j+k;
        w[l]=oneweight[i]*oneweight[j]*oneweight[k];
        x[l]=onexi[k];
        e[l]=onexi[j];
        z[l]=onexi[i];
      }

  InitObject(64,w,x,e,z, 7);
}

void TQuadFormulaHexa::Gauss5()
{
  if(Weights!=NULL) return;
  double onexi[]={
                0,
                0.538469310105683091036314421,
                0.906179845938663992797626878,
               -0.906179845938663992797626878,
               -0.538469310105683091036314421
             }; 
  double oneweight[]={
                0.568888888888888888888888889,
                0.478628670499366468041291515,
                0.236926885056189087514264041,
                0.236926885056189087514264041,
                0.478628670499366468041291515
             };

  double w[125], x[125], e[125], z[125];

  int l;
  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
      for(int k=0;k<5;k++)
      {
        l=25*i+5*j+k;
        w[l]=oneweight[i]*oneweight[j]*oneweight[k];
        x[l]=onexi[k];
        e[l]=onexi[j];
        z[l]=onexi[i];
      }

  InitObject(125,w,x,e,z, 9);
}

void TQuadFormulaHexa::Gauss6()
{
  if(Weights!=NULL) return;
  double onexi[]={
               -0.238619186083196908630501721681, 
                0.238619186083196908630501721681,
                0.661209386466264513661399595021, 
                0.932469514203152027812301554495,
               -0.932469514203152027812301554495,
               -0.661209386466264513661399595021
             }; 
  double oneweight[]={
               0.467913934572691047389870343990, 
               0.467913934572691047389870343990,
               0.360761573048138607569833513836, 
               0.171324492379170345040296142175,
               0.171324492379170345040296142175,
               0.360761573048138607569833513836
             };

  double w[216], x[216], e[216], z[216];

  int l;
  for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
      for(int k=0;k<6;k++)
      {
        l=36*i+6*j+k;
        w[l]=oneweight[i]*oneweight[j]*oneweight[k];
        x[l]=onexi[k];
        e[l]=onexi[j];
        z[l]=onexi[i];
      }

  InitObject(216,w,x,e,z, 11);
}

void TQuadFormulaHexa::Gauss7()
{
  if(Weights!=NULL) return;
  double onexi[]={
               -0.405845151377397166906606412077, 
                0,
                0.405845151377397166906606412077,
                0.741531185599394439863864773281, 
                0.949107912342758524526189684048,
               -0.949107912342758524526189684048, 
               -0.741531185599394439863864773281
             }; 
  double oneweight[]={
                0.381830050505118944950369775484,
                0.417959183673469387755102040812,
                0.381830050505118944950369775484,
                0.279705391489276667901467771426,
                0.129484966168869693270611432679,
                0.129484966168869693270611432679,
                0.279705391489276667901467771426
             };

  double w[343], x[343], e[343], z[343];

  int l;
  for(int i=0;i<7;i++)
    for(int j=0;j<7;j++)
      for(int k=0;k<7;k++)
      {
        l=49*i+7*j+k;
        w[l]=oneweight[i]*oneweight[j]*oneweight[k];
        x[l]=onexi[k];
        e[l]=onexi[j];
        z[l]=onexi[i];
      }

  InitObject(343,w,x,e,z, 13);
}

void TQuadFormulaHexa::Gauss8()
{
  if(Weights!=NULL) return;
  double onexi[]={
                -0.525532409916328985817739049189, 
                -0.183434642495649804939476142360,
                 0.183434642495649804939476142360, 
                 0.525532409916328985817739049189,
                 0.796666477413626739591553936476, 
                 0.960289856497536231683560868569,
                -0.960289856497536231683560868569,
                -0.796666477413626739591553936476
             }; 
  double oneweight[]={
                0.313706645877887287337962201982,
                0.362683783378361982965150449268,
                0.362683783378361982965150449268,
                0.313706645877887287337962201982,
                0.222381034453374470544355994424,
                0.101228536290376259152531354312,
                0.101228536290376259152531354312,
                0.222381034453374470544355994424
             };

  double w[512], x[512], e[512], z[512];

  int l;
  for(int i=0;i<8;i++)
    for(int j=0;j<8;j++)
      for(int k=0;k<8;k++)
      {
        l=64*i+8*j+k;
        w[l]=oneweight[i]*oneweight[j]*oneweight[k];
        x[l]=onexi[k];
        e[l]=onexi[j];
        z[l]=onexi[i];
      }

  InitObject(512,w,x,e,z,15);
}

void TQuadFormulaHexa::Gauss9()
{
  if(Weights!=NULL) return;
  double onexi[]={
               -0.613371432700590397308702039341, 
               -0.324253423403808929038538014644, 
                0,
                0.324253423403808929038538014644, 
                0.613371432700590397308702039341,
                0.836031107326635794299429788070, 
                0.968160239507626089835576202904,
               -0.968160239507626089835576202904,
               -0.836031107326635794299429788070
             }; 
  double oneweight[]={
                0.260610696402935462318742869422, 
                0.312347077040002840068630406582,
                0.330239355001259763164525069284, 
                0.312347077040002840068630406582,
                0.260610696402935462318742869422, 
                0.180648160694857404058472031240,
                0.0812743883615744119718921581134,
                0.0812743883615744119718921581134,
                0.180648160694857404058472031240
             };

  double w[729], x[729], e[729], z[729];

  int l;
  for(int i=0;i<9;i++)
    for(int j=0;j<9;j++)
      for(int k=0;k<9;k++)
      {
        l=81*i+9*j+k;
        w[l]=oneweight[i]*oneweight[j]*oneweight[k];
        x[l]=onexi[k];
        e[l]=onexi[j];
        z[l]=onexi[i];
      }

  InitObject(729,w,x,e,z,17);
}

void TQuadFormulaHexa::VerticesAndOrigin()
{
  double t = 1.0/3.0;
  
  if(Weights!=NULL) return;
  double w[9]={ t, t, t, t,
                t, t, t, t, 16.0/3.0 };
  double x[9]={ -1.0,  1.0, -1.0,  1.0,
                -1.0,  1.0, -1.0,  1.0, 0.0 };
  double e[9]={ -1.0, -1.0,  1.0,  1.0,
                -1.0, -1.0,  1.0,  1.0, 0.0 };
  double z[9]={ -1.0, -1.0, -1.0, -1.0,
                 1.0,  1.0,  1.0,  1.0, 0.0 };

  InitObject(9,w,x,e,z, 1);
}
void TQuadFormulaHexa::VerticesAndOrigin15()
{
  double t = 1.0/9.0, s= 20.0/9;
  
  if(Weights!=NULL) return;
  double w[15]={ t, t, t, t,
                t, t, t, t, - 56.0/9.0,
                s, s, s, s, s, s};
  double x[15]={ -1.0,  1.0, -1.0,  1.0,
                -1.0,  1.0, -1.0,  1.0, 0.0,
                 0.63245553203367586, 0.0, 0.0,
                 -0.63245553203367586, 0.0, 0.0};
  double e[15]={ -1.0, -1.0,  1.0,  1.0,
                -1.0, -1.0,  1.0,  1.0, 0.0,
                 0.0, 0.63245553203367586, 0.0,
                 0.0, -0.63245553203367586, 0.0};
  double z[15]={ -1.0, -1.0, -1.0, -1.0,
                 1.0,  1.0,  1.0,  1.0, 0.0,
                 0.0, 0.0,  0.63245553203367586,
                 0.0, 0.0,  -0.63245553203367586};

  InitObject(15,w,x,e,z, 1);
}

void TQuadFormulaHexa::VerticesAndOrigin57()
{
  double t0 = 0.010323148691722060, t1 = 0.69052862985186805;
  double t2 = 0.071977066056033603, t3 = 0.090534233349122512;
  double t4 = 0.14514879873706936;
  double x1 = 0.49757102301189591, x2 = 0.98959984011500087;
  double x3 = 0.87103917604353433, x4 = 0.53081040829645735;
  double x5 = 0.87906774032385155;
  
  if(Weights!=NULL) return;
  double w[57]={ t0, t0, t0, t0,
                 t0, t0, t0, t0, 
                 -1.2276013348603214,
                 t1, t1, t1, t1, t1, t1,
                 t2, t2, t2, t2, t2, t2,
                 t3, t3, t3, t3, t3, t3,
                 t3, t3, t3, t3, t3, t3,
                 t4, t4, t4, t4, t4, t4,
                 t4, t4, t4, t4, t4, t4,
                 t4, t4, t4, t4, t4, t4,
                 t4, t4, t4, t4, t4, t4};

  double x[57]={ -1.0,  1.0, -1.0,  1.0,
                 -1.0,  1.0, -1.0,  1.0, 
                 0.0,
                 x1 , 0.0, 0.0, -x1, 0.0, 0.0,
                 x2 , 0.0, 0.0, -x2, 0.0, 0.0,
                 x3 , x3,  x3,  x3 , -x3, -x3,
                 -x3 , -x3, 0.0, 0.0 , 0.0, 0.0,
                 x4, x4,  x4,  x4, x4, x4,  x4,  x4,
                 -x4, -x4,  -x4,  -x4, -x4, -x4,  -x4,  -x4,
                 x5, x5, x5, x5, -x5, -x5, -x5, -x5
                 };
  double e[57]={ -1.0, -1.0,  1.0,  1.0,
                 -1.0, -1.0,  1.0,  1.0, 
                 0.0,
                 0.0, x1, 0.0, 0.0, -x1, 0.0,
                 0.0, x2, 0.0, 0.0, -x2, 0.0,
                 x3 , -x3 , 0.0, 0.0,  x3 , -x3, 
                 0.0, 0.0, x3, x3, -x3, -x3,
                 x4,  -x4, x5, -x5, x4,  -x4, x5, -x5,
                 x4,  -x4, x5, -x5, x4,  -x4, x5, -x5,
                 x4, x4, -x4, -x4, x4, x4, -x4, -x4};
  double z[57]={ -1.0, -1.0, -1.0, -1.0,
                 1.0,  1.0,  1.0,  1.0, 0.0,
                 0.0, 0.0, x1, 0.0, 0.0,  -x1,
                 0.0, 0.0, x2, 0.0, 0.0,  -x2,
                 0.0, 0.0, x3, -x3, 0.0, 0.0,
                 x3, -x3, x3, -x3, x3, -x3,                 
                 x5, x5, x4, x4,  -x5, -x5, -x4, -x4,
                 x5, x5, x4, x4,  -x5, -x5, -x4, -x4,
                 x4, -x4, x4, -x4, x4, -x4, x4, -x4
                 };

  InitObject(57,w,x,e,z, 1);
}
void TQuadFormulaHexa::Degree7_Points38()
{
  double t0 = 0.29518973826262290, t1 = 0.40405541726620058;
  double t2 =0.12485075967894408; 
  double c0,c1,c2, c3;
  
  if(Weights!=NULL) return;
  double w[38]={ t0, t0, t0, t0, t0, t0,
                 t1, t1, t1, t1, t1, t1, t1, t1,
                 t2, t2, t2, t2, t2, t2, t2, t2, t2, t2,
                 t2, t2, t2, t2, t2, t2, t2, t2, t2, t2,
                 t2, t2, t2, t2 };

  c0 = 0.90168780782129128;
  c1 = 0.40837222149947467;
  c2 = 0.85952309020105419;
  c3 = 0.41473591372798772;

  double x[38]={ c0, 0, 0, -c0, 0, 0,
                 c1, -c1, c1, c1, c1, -c1, -c1, -c1,
                 c2, c2, c3, -c2, -c2, -c3, c2, c2, c3, c2, c2, c3,
                 c2, c2, c3, -c2, -c2, -c3, -c2, -c2, -c3, -c2, -c2, -c3};
  double e[38]={ 0, c0, 0, 0, -c0, 0,
                 c1, c1, -c1, c1, -c1, c1, -c1, -c1,
                 c2, c3, c2, c2, c3, c2, -c2, -c3, -c2, c2, c3, c2,
                 -c2, -c3, -c2, c2, c3, c2, -c2, -c3, -c2,-c2, -c3, -c2};

  double z[38]={ 0, 0, c0, 0, 0, -c0,
                 c1, c1, c1, -c1, -c1, -c1, c1, -c1,
                 c3, c2, c2, c3, c2, c2, c3, c2, c2, -c3, -c2, -c2,
                 -c3, -c2, -c2, -c3, -c2, -c2, c3, c2, c2, -c3, -c2, -c2};

  InitObject(38,w,x,e,z, 1);
}
