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
// @(#)TriaIsoparametric.C        1.8 04/13/00
//
// Class:      TTriaIsoparametric
//
// Purpose:    reference transformations for triangle
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#include <TriaIsoparametric.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <IsoJointEqN.h>
#include <BoundComp.h>
#include <string.h>
#include <stdlib.h>

BaseFunct2D TTriaIsoparametric::BaseFunctFromOrder[] = { 
                BF_C_T_P0_2D, BF_C_T_P1_2D, BF_C_T_P2_2D, BF_C_T_P3_2D,
                BF_C_T_P4_2D, BF_C_T_P5_2D, BF_C_T_P6_2D, BF_C_T_P7_2D,
                BF_C_T_P8_2D, BF_C_T_P9_2D };

FEDesc2D TTriaIsoparametric::FEDescFromOrder[] = { 
                FE_C_T_P0_2D, FE_C_T_P1_2D, FE_C_T_P2_2D, FE_C_T_P3_2D,
                FE_C_T_P4_2D, FE_C_T_P5_2D, FE_C_T_P6_2D, FE_C_T_P7_2D,
                FE_C_T_P8_2D, FE_C_T_P9_2D };

/** constuctor */
TTriaIsoparametric::TTriaIsoparametric()
{
}

/** transfer from reference element to original element */
void TTriaIsoparametric::GetOrigFromRef(double xi, double eta, 
                                        double &X, double &Y)
{
/*
  int i, j;

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(fabs(xi-XI[i])<1e-8 && fabs(eta-ETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in TTriaIsoparametric::GetOrigFromRef(1)" << endl;
    return;
  }

  X = xc0 + xc1*xi + xc2*eta;
  Y = yc0 + yc1*xi + yc2*eta;

  for(i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
  }

  // cout << "(" << xi << ", " << eta << ") = ";
  // cout << "(" << X << "," << Y << ")" << endl;
*/
  int N_Points;
  double Xi[1], Eta[1], Xarray[1], Yarray[1], absdetjk[1];

  N_Points = 1;
  Xi[0] = xi; Eta[0] = eta;
  Xarray[0] = X; Yarray[0] = Y;
  GetOrigFromRef(N_Points, Xi, Eta, Xarray, Yarray, absdetjk);
  X = Xarray[0]; Y = Yarray[0];
}

/** transfer a set of points from reference to original element */
void TTriaIsoparametric::GetOrigFromRef(int N_Points, 
                                double *xi, double *eta,
                                double *X, double *Y, double *absdetjk)
{
/*
  int i, j, k;
  double Xi, Eta, a11, a12, a21, a22;

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];

    a11 = xc1;
    a21 = xc2;
    a12 = yc1;
    a22 = yc2;

    j = -1;
    for(k=0;k<N_QuadPoints;k++)
    {
      if(fabs(Xi-XI[k])<1e-8 && fabs(Eta-ETA[k])<1e-8)
      {
        j = k;
        break;
      }
    } // endfor

    if(j==-1)
    {
      cout << "error in TTriaIsoparametric::GetOrigFromRef(2)" << endl;
      return;
    }

    // cout << "j= " << j << endl;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta;

    // cout << "(" << Xi << ", " << Eta << ") = ";
    // cout << "(" << X[i] << "," << Y[i] << ") -> ";
    for(k=0;k<N_AuxPoints;k++)
    {
      // cout << "fct: " << FctValues[j][k] << endl;
      X[i] += XDistance[k] * FctValues[j][k];
      Y[i] += YDistance[k] * FctValues[j][k];

      a11 += XDistance[k] * XiDerValues[j][k];
      a21 += XDistance[k] * EtaDerValues[j][k];

      a12 += YDistance[k] * XiDerValues[j][k];
      a22 += YDistance[k] * EtaDerValues[j][k];
    }

    absdetjk[i] = fabs(a11*a22 - a12*a21);

    // cout << "(" << X[i] << "," << Y[i] << ")" << endl;
  } // endfor i
*/
  int i, j, k;
  double Xi, Eta;
  double dx1, dx2;
  double dy1, dy2;
  double detjk;
  double AuxVector[3*MaxN_BaseFunctions2D];
  TBaseFunct2D *bf;

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];

    dx1 = xc1;
    dx2 = xc2;
    
    dy1 = yc1;
    dy2 = yc2;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta;

    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);
    bf->GetDerivatives(D00, Xi, Eta, AuxVector);
    bf->GetDerivatives(D10, Xi, Eta, AuxVector+MaxN_BaseFunctions2D);
    bf->GetDerivatives(D01, Xi, Eta, AuxVector+2*MaxN_BaseFunctions2D);

    for(k=0;k<N_AuxPoints;k++)
    {
      j = IntAux[k];
      X[i] += XDistance[k] * AuxVector[j];
      Y[i] += YDistance[k] * AuxVector[j];

      dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];

      dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];
    }

    // cout << "2-x: " << X[i] << " y: " << Y[i] << " z: " << Z[i] << endl;
    // cout << "3-x: " << Xi << " " << Eta << " " << Zeta;
    // cout << " " << dx1 << " y: " << dx2 << " z: " << dx3 << endl;

    detjk = dx1*dy2- dx2*dy1;
    absdetjk[i] = fabs(detjk);

    // cout << endl;
    // cout << Xi << " " << Eta << " " << Zeta << endl;
    // cout << X[i] << " " << Y[i] << " " << Z[i] << endl;
    // cout << detjk << endl;
    // cout << dx1 << " " << dx2 << " " << dx3 << endl;
    // cout << dy1 << " " << dy2 << " " << dy3 << endl;
    // cout << dz1 << " " << dz2 << " " << dz3 << endl;
  }
}

/** transfer from reference element to original element */
void TTriaIsoparametric::GetOrigFromRef(double *ref, double *orig)
{
  int i, j;
  double X, Y, xi, eta;

  xi = ref[0];
  eta = ref[1];

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(fabs(xi-XI[i])<1e-8 && fabs(eta-ETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in TTriaIsoparametric::GetOrigFromRef(3)" << endl;
    return;
  }

  X = xc0 + xc1*xi + xc2*eta;
  Y = yc0 + yc1*xi + yc2*eta;

  for(i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
  }

  orig[0] = X;
  orig[1] = Y;

  // cout << "(" << xi << ", " << eta << ") = ";
  // cout << "(" << X << "," << Y << ")" << endl;
}

/** transfer from original element to reference element */
void TTriaIsoparametric::GetRefFromOrig(double X, double Y, 
                                        double &xi, double &eta)
{
  cout << "not implemented yet 1 " << endl;
}

/** transfer from original element to reference element */
void TTriaIsoparametric::GetRefFromOrig(double *orig, double *ref)
{
  cout << "not implemented yet 2" << endl;
}

/** calculate functions and derivatives from reference element
    to original element */
void TTriaIsoparametric::GetOrigValues(BaseFunct2D BaseFunct,
                               int N_Points, double *xi, double *eta,
                               int N_Functs, QuadFormula2D formula)
{
  cout << "not implemented yet 3" << endl;
}

/** calculate functions and derivatives from reference element
    to original element, for all given elements */
void TTriaIsoparametric::GetOrigValues(int N_Sets, BaseFunct2D *BaseFuncts,
                               int N_Points, double *xi, double *eta,
                               QuadFormula2D formula,
                               bool *Needs2ndDer)
{
  int i,j,k,N_, start, end;
  double **refvaluesD00, **origvaluesD00;
  double **refvaluesD10, **origvaluesD10;
  double **refvaluesD01, **origvaluesD01;
  double **refvaluesD20, **origvaluesD20;
  double **refvaluesD11, **origvaluesD11;
  double **refvaluesD02, **origvaluesD02;
  double *refD00, *origD00;
  double *refD10, *origD10;
  double *refD01, *origD01;
  double *refD20, *origD20;
  double *refD11, *origD11;
  double *refD02, *origD02;
  double r20, r11, r02, o20, o11, o02;
  double *aux;
  double GeoData[3][3];
  double Eye[3][3];
  BaseFunct2D BaseFunct;
  int N_Functs;
  bool SecondDer;
  double rec_detjk, a11, a12, a21, a22;

  SecondDer = FALSE;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();

    refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D00);
    if(refvaluesD00==NULL)
    {
      TFEDatabase2D::GetBaseFunct2D(BaseFunct)->MakeRefElementData(formula);
      refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                                    formula, D00);
    }
  
    origvaluesD00=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    if(origvaluesD00==NULL)
    {
      origvaluesD00 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD00[j] = aux+j*MaxN_BaseFunctions2D;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D00, origvaluesD00);
    }
  
    for(j=0;j<N_Points;j++)
    {
      refD00 = refvaluesD00[j];
      origD00 = origvaluesD00[j];

      memcpy(origD00, refD00, N_Functs*SizeOfDouble);
    } // endfor j

    refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D10);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    if(origvaluesD10==NULL)
    {
      origvaluesD10 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD10[j] = aux+j*MaxN_BaseFunctions2D;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D10, origvaluesD10);
    }
  
    refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D01);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);
    if(origvaluesD01==NULL)
    {
      origvaluesD01 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD01[j] = aux+j*MaxN_BaseFunctions2D;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D01, origvaluesD01);
    }
  
    if(Needs2ndDer[i])
    {
      SecondDer = TRUE;

      refvaluesD20=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                        formula, D20);
      origvaluesD20=TFEDatabase2D::GetOrigElementValues(BaseFunct, D20);
      if(origvaluesD20==NULL)
      {
        origvaluesD20 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD20[j] = aux+j*MaxN_BaseFunctions2D;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D20, origvaluesD20);
      }
    
      refvaluesD11=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                        formula, D11);
      origvaluesD11=TFEDatabase2D::GetOrigElementValues(BaseFunct, D11);
      if(origvaluesD11==NULL)
      {
        origvaluesD11 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD11[j] = aux+j*MaxN_BaseFunctions2D;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D11, origvaluesD11);
      }
    
      refvaluesD02=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                        formula, D02);
      origvaluesD02=TFEDatabase2D::GetOrigElementValues(BaseFunct, D02);
      if(origvaluesD02==NULL)
      {
        origvaluesD02 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD02[j] = aux+j*MaxN_BaseFunctions2D;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D02, origvaluesD02);
      }
    } // endfor Needs2ndDer[i]
  } // endfor i
  
  // D10 and D01
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();

    refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D10);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);

    refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D01);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);

    for(j=0;j<N_Points;j++)
    {
      refD10 = refvaluesD10[j];
      origD10 = origvaluesD10[j];
  
      refD01 = refvaluesD01[j];
      origD01 = origvaluesD01[j];
  
      a11 = xc1;
      a21 = xc2;
      a12 = yc1;
      a22 = yc2;
    
      for(k=0;k<N_AuxPoints;k++)
      {
        a11 += XDistance[k] * XiDerValues[j][k];
        a21 += XDistance[k] * EtaDerValues[j][k];
    
        a12 += YDistance[k] * XiDerValues[j][k];
        a22 += YDistance[k] * EtaDerValues[j][k];
      }
    
      rec_detjk = 1 / (a11*a22 - a12*a21);
    
      for(k=0;k<N_Functs;k++)
      {
        origD10[k]=( a22*refD10[k]-a12*refD01[k]) * rec_detjk;
        origD01[k]=(-a21*refD10[k]+a11*refD01[k]) * rec_detjk;
      } // endfor k
    } // endfor j
  } // endfor i

  // leave if no second derivatives are needed
  if(!SecondDer) return;
  
  // implement calculation of second derivatives later

}

/** calculate functions and derivatives from reference element
    to original element */
void TTriaIsoparametric::GetOrigValues(double xi, double eta, 
                int N_BaseFunct,
                double *uref, double *uxiref, double *uetaref,
                double *uorig, double *uxorig, double *uyorig)
{
/*
  int i, j, k;
  double a11, a12, a21, a22, rec_detjk;

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(fabs(xi-XI[i])<1e-8 && fabs(eta-ETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in TTriaIsoparametric::GetOrigFromRef(1)" << endl;
    return;
  }

  // D00
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  a11 = xc1;
  a21 = xc2;
  a12 = yc1;
  a22 = yc2;

  for(k=0;k<N_AuxPoints;k++)
  {
    a11 += XDistance[k] * XiDerValues[j][k];
    a21 += XDistance[k] * EtaDerValues[j][k];

    a12 += YDistance[k] * XiDerValues[j][k];
    a22 += YDistance[k] * EtaDerValues[j][k];
  } // endfor k

  rec_detjk = 1 / (a11*a22 - a12*a21);

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( a22*uxiref[i]-a12*uetaref[i]) * rec_detjk;
    uyorig[i]=(-a21*uxiref[i]+a11*uetaref[i]) * rec_detjk;
  } // endfor i
*/
  int i, j, k;
  double dx1, dx2;
  double dy1, dy2;
  double rec_detjk;
  double AuxVector[3*MaxN_BaseFunctions2D];
  TBaseFunct2D *bf;


  dx1 = xc1;
  dx2 = xc2;
    
  dy1 = yc1;
  dy2 = yc2;
    
  bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);
  bf->GetDerivatives(D10, xi, eta, AuxVector+MaxN_BaseFunctions2D);
  bf->GetDerivatives(D01, xi, eta, AuxVector+2*MaxN_BaseFunctions2D);

  for(k=0;k<N_AuxPoints;k++)
  {
    j = IntAux[k];

    dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
    dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];

    dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
    dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];
  }

  rec_detjk = 1/(dx1*dy2 - dx2*dy1);

  // D00
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( dy2*uxiref[i]-dy1*uetaref[i]) * rec_detjk;
    uyorig[i]=(-dx2*uxiref[i]+dx1*uetaref[i]) * rec_detjk;
  } // endfor i
}

/** calculate functions and derivatives from reference element
    to original element */
void TTriaIsoparametric::GetOrigValues(int joint, double zeta,
            int N_BaseFunct,
            double *uref, double *uxiref, double *uetaref,
            double *uorig, double *uxorig, double *uyorig)
{

  int i, j, k;
  double a11, a12, a21, a22, rec_detjk;
  double xi, eta;
  TBaseFunct2D *bf;

  bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);

  double valxi[MaxN_BaseFunctions2D];
  double valeta[MaxN_BaseFunctions2D];
  double *values[1];

  switch(joint)
  {
    case 0:
      xi = 0.5*(1+zeta); eta = 0;
    break;

    case 1:
      xi = 0.5*(1-zeta); eta = 0.5*(1+zeta);
    break;

    case 2:
      xi = 0; eta = 0.5*(1-zeta);
    break;
  }

  a11 = xc1;
  a21 = xc2;
  a12 = yc1;
  a22 = yc2;

  values[0] = valxi;
  bf->GetValues(1, &zeta, joint, D10, values);
  values[0] = valeta;
  bf->GetValues(1, &zeta, joint, D01, values);

  // add correction due to isoparamatric boundary approximation
  for(i=0;i<N_AuxPoints;i++)
  {
    k = IntAux[i];
    a11 += XDistance[i] * valxi[k];
    a21 += XDistance[i] * valeta[k];

    a12 += YDistance[i] * valxi[k];
    a22 += YDistance[i] * valeta[k];
  } // endfor i

  rec_detjk = 1. / (a11*a22 - a12*a21);

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( a22*uxiref[i]-a12*uetaref[i]) * rec_detjk;
    uyorig[i]=(-a21*uxiref[i]+a11*uetaref[i]) * rec_detjk;
    uorig[i] = uref[i];
  } // endfor i
}

void TTriaIsoparametric::SetCell(TBaseCell *cell)
{
  int i, j, k, N_;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp2D *comp;
  TInterfaceJoint *interface;
  JointType type;
  double t0, t1, t, dt;
  double xa, ya, xe, ye, xm, ym, xp, yp, dx, dy;
  int compid;
  TBaseFunct2D *bf;
  TFEDesc2D *fedesc;
  int *JointDOF;
  int N_Vertices;
  TVertex **Vertices;
#ifdef __3D__
  double z[3], zp;
#endif

  N_AuxPoints = 0;

  Cell = cell;

  fedesc = TFEDatabase2D::GetFEDesc2D(FEDescFromOrder[ApproximationOrder]);
  bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);

  TFEDatabase2D::GetQuadFormula2D(QuadFormula)->GetFormulaData(N_QuadPoints, W, XI, ETA);

  // cout << endl;
  for(i=0;i<3;i++)
  {
#ifdef __3D__
    Cell->GetVertex(i)->GetCoords(x[i], y[i], z[i]);
#else
    Cell->GetVertex(i)->GetCoords(x[i], y[i]);
#endif
    // cout << setw(20) << x[i] << setw(20) << y[i] << endl;
  }

  xc0=x[0];
  xc1=x[1]-x[0];
  xc2=x[2]-x[0];

  yc0=y[0];
  yc1=y[1]-y[0];
  yc2=y[2]-y[0];

  for(i=0;i<3;i++)
  {
    // check whether the joint i is curved
    joint = Cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryEdge || type == InterfaceJoint)
    {
//       cout << "joint: " << i << endl;
      xa = x[i];
      ya = y[i];
      xe = x[(i+1) % 3];
      ye = y[(i+1) % 3];
      if(type == BoundaryEdge)
      {
        boundedge = (TBoundEdge *)(joint);
        comp = boundedge->GetBoundComp();
        compid = comp->GetID();
        boundedge->GetParameters(t0, t1);
      }
      else
      {
        interface = (TInterfaceJoint *)(joint);
        comp = interface->GetBoundComp();
        compid = comp->GetID();

        if(Cell == interface->GetNeighbour(0))
          interface->GetParameters(t0, t1); // correct order
        else
          interface->GetParameters(t1, t0); // backward order
      }

      JointDOF = fedesc->GetJointDOF(i);

      dt = (t1-t0)/ApproximationOrder;
      dx = (xe-xa)/ApproximationOrder;
      dy = (ye-ya)/ApproximationOrder;
      for(j=1;j<ApproximationOrder;j++)
      {
        t = t0 + j * dt;
        comp->GetXYofT(t, xp, yp);
        xm = xa + j * dx;
        ym = ya + j * dy;

        // cout << "m: (" << xm << ", " << ym << ")" << endl;
        // cout << "p: (" << xp << ", " << yp << ")" << endl;
        // cout << "d: (" << xp-xm << ", " << yp-ym << ")" << endl;

        if(fabs(xp-xm) > 1e-8 || fabs(yp-ym) > 1e-8)
        {
          // curved boundary
          XDistance[N_AuxPoints] = xp - xm;
          YDistance[N_AuxPoints] = yp - ym;
          IntAux[N_AuxPoints] = JointDOF[j];
          // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
          N_AuxPoints++;
        } // endif 
      } // endfor
    } // endif
    else
    {
      if(type == IsoInterfaceJoint || type == IsoBoundEdge ||
         type == IsoJointEqN)
      {
        switch(type)
        {
          case IsoInterfaceJoint:
            N_Vertices = ((TIsoInterfaceJoint *)joint)->GetN_Vertices();
            Vertices = ((TIsoInterfaceJoint *)joint)->GetVertices();
          break;

          case IsoBoundEdge:
            N_Vertices = ((TIsoBoundEdge *)joint)->GetN_Vertices();
            Vertices = ((TIsoBoundEdge *)joint)->GetVertices();
          break;

          case IsoJointEqN:
            N_Vertices = ((TIsoJointEqN *)joint)->GetN_Vertices();
            Vertices = ((TIsoJointEqN *)joint)->GetVertices();
          break;
	  default:
	    break;
        } // endswitch

        
        // cout << "N_AuxVertices: " << N_Vertices << endl;
        ApproximationOrder = N_Vertices+1;

        fedesc = TFEDatabase2D::GetFEDesc2D(FEDescFromOrder[ApproximationOrder]);
        bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);

        JointDOF = fedesc->GetJointDOF(i);

        xa = x[i];
        ya = y[i];
        xe = x[(i+1) % 3];
        ye = y[(i+1) % 3];

//         dt = (t1-t0)/ApproximationOrder;
        dx = (xe-xa)/ApproximationOrder;
        dy = (ye-ya)/ApproximationOrder;
        for(j=1;j<ApproximationOrder;j++)
        {
//           t = t0 + j * dt;
#ifdef __3D__
          Vertices[j-1]->GetCoords(xp, yp, zp);
#else
          Vertices[j-1]->GetCoords(xp, yp);
#endif
          xm = xa + j * dx;
          ym = ya + j * dy;
  
          // cout << "(" << xm << ", " << ym << ")" << endl;
          // cout << "(" << xp << ", " << yp << ")" << endl;
  
          if(fabs(xp-xm) > 1e-8 || fabs(yp-ym) > 1e-8)
          {
            // curved boundary
            XDistance[N_AuxPoints] = xp - xm;
            YDistance[N_AuxPoints] = yp - ym;
            IntAux[N_AuxPoints] = JointDOF[j];
	    
	    //if(XDistance[N_AuxPoints]<=0)
	     // cout <<  " xa "  << xa <<    " xe " << xe <<  " xp " << xp<< " XDist " << XDistance[N_AuxPoints] << endl;
	      
            // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
            N_AuxPoints++;
          } // endif 
        } // endfor
      } // endif
    } // endelse
  } // endfor

  if(N_AuxPoints && ApproximationOrder==3)
  {
    xm = 0;
    ym = 0;
    for(i=0;i<N_AuxPoints;i++)
    {
      xm += XDistance[i];
      ym += YDistance[i];
    }
    XDistance[N_AuxPoints] = xm/4;
    YDistance[N_AuxPoints] = ym/4;
    IntAux[N_AuxPoints] = 5;
    N_AuxPoints++;
  }

  if(N_AuxPoints)
  {
    for(i=0;i<N_QuadPoints;i++)
    {
      bf->GetDerivatives(D00, XI[i], ETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        FctValues[i][j] = DoubleAux[IntAux[j]]; 

      bf->GetDerivatives(D10, XI[i], ETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        XiDerValues[i][j] = DoubleAux[IntAux[j]]; 

      bf->GetDerivatives(D01, XI[i], ETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        EtaDerValues[i][j] = DoubleAux[IntAux[j]]; 

    } // endfor i
  } // endif
} // TTriaIsoparametric::SetCell

/** return outer normal vector */
void TTriaIsoparametric::GetOuterNormal(int j, double zeta,
                                        double &n1, double &n2)
{
  double len;
  double t1, t2;

  GetTangent(j, zeta, t1, t2);
  len = sqrt(t1*t1+t2*t2);

  n1 =  t2/len;
  n2 = -t1/len;
}

/** return tangent */
void TTriaIsoparametric::GetTangent(int j, double zeta,
                                        double &t1, double &t2)
{
  TBaseFunct2D *bf;
  bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);
  double xi, eta;
  double valx[MaxN_BaseFunctions2D];
  double valy[MaxN_BaseFunctions2D];
  double *values[1];
  double xip, etap;
  double a11, a12, a21, a22;
  int i, k;

  switch(j)
  {
    case 0:
      xip = 0.5;
      etap = 0;
    break;

    case 1:
      xip = -0.5;
      etap = 0.5;
    break;

    case 2:
      xip = 0;
      etap = -0.5;
    break;

    default:
      cerr << "wrong joint number" << endl;
      t1 = -1;
      t2 = -1;
  } // endswitch

  a11 = xc1;
  a21 = xc2;
  a12 = yc1;
  a22 = yc2;

  values[0] = valx;
  bf->GetValues(1, &zeta, j, D10, values);
  values[0] = valy;
  bf->GetValues(1, &zeta, j, D01, values);

  for(i=0;i<N_AuxPoints;i++)
  {
    k = IntAux[i];
    a11 += XDistance[i] * valx[k];
    a21 += XDistance[i] * valy[k];

    a12 += YDistance[i] * valx[k];
    a22 += YDistance[i] * valy[k];
  }

  t1 = xip*a11 + etap*a21;
  t2 = xip*a12 + etap*a22;
}

/** return volume of cell according to isoparametric boundary */
double TTriaIsoparametric::GetVolume()
{
  double locvol;

  int i, j, k;
  double Xi, Eta, a11, a12, a21, a22;
  double absdetjk;

  locvol = 0;
  for(i=0;i<N_QuadPoints;i++)
  {
    Xi = XI[i];
    Eta = ETA[i];

    a11 = xc1;
    a21 = xc2;
    a12 = yc1;
    a22 = yc2;

    for(k=0;k<N_AuxPoints;k++)
    {
      a11 += XDistance[k] * XiDerValues[i][k];
      a21 += XDistance[k] * EtaDerValues[i][k];

      a12 += YDistance[k] * XiDerValues[i][k];
      a22 += YDistance[k] * EtaDerValues[i][k];
    }
    absdetjk = fabs(a11*a22 - a12*a21);

    locvol += absdetjk*W[i];
  } // endfor i
  
  return locvol;
}

void TTriaIsoparametric::GetOrigBoundFromRef( int joint, int N_Points,
                                              double *zeta, double *X, double *Y)
{
  int i, j, k;
  double xi,  eta;
  double AuxVector[3*MaxN_BaseFunctions2D];
  TBaseFunct2D *bf;

  for(i=0;i<N_Points;i++)
  {
   switch(joint)
    {
     case 0:
      xi = 0.5*(1+zeta[i]); eta = 0;
     break;

     case 1:
      xi = 0.5*(1-zeta[i]); eta = 0.5*(1+zeta[i]);
     break;

     case 2:
      xi = 0.; eta = 0.5*(1-zeta[i]);
     break;
    }



    X[i] = xc0 + xc1*xi + xc2*eta;
    Y[i] = yc0 + yc1*xi + yc2*eta;

    
//     cout << i << " "  <<joint << " xi " <<  X[i] << " eta " << Y[i] << endl;
  
    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);
    bf->GetDerivatives(D00, xi, eta, AuxVector);
    for(k=0;k<N_AuxPoints;k++)
     {
      j = IntAux[k];
      X[i] += XDistance[k] * AuxVector[j];
      Y[i] += YDistance[k] * AuxVector[j];
      //cout << i << " "  <<joint << " isoxi " <<  X[i] << " XDistance[k] " << XDistance[k] << " xi " << xi << " eta " << eta  << endl;     
     }
     

   }

}
