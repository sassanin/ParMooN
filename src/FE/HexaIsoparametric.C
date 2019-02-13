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
// @(#)HexaIsoparametric.C      1.6 02/22/00
//
// Class:      THexaIsoparametric
//
// Purpose:    trilinear reference transformations for hexahedron
//
// Author:     Daniel Quoos  
//
// History:    11.07.00 start implementation
// 
// =======================================================================

#include <HexaIsoparametric.h>
#include <BoundFace.h>
#include <IsoBoundFace.h>
#include <IsoJointEqN.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>
#include <MooNMD_Io.h>

BaseFunct3D THexaIsoparametric::BaseFunctFromOrder[] = { 
                BF_C_H_Q0_3D, BF_C_H_Q1_3D, BF_C_H_Q2_3D, BF_C_H_Q3_3D,
                BF_C_H_Q4_3D, BF_C_H_Q00_3D};

FEDesc3D THexaIsoparametric::FEDescFromOrder[] = { 
                FE_C_H_Q0_3D, FE_C_H_Q1_3D, FE_C_H_Q2_3D, FE_C_H_Q3_3D,
                FE_C_H_Q4_3D, FE_C_H_Q00_3D};

/** constuctor */
THexaIsoparametric::THexaIsoparametric()
{
}

/** transfer from reference element to original element */
void THexaIsoparametric::GetOrigFromRef(double xi, double eta, double zeta,
                                        double &X, double &Y, double &Z)
{
  int i,j;

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(fabs(xi-XI[i])<1e-8 && fabs(eta-ETA[i])<1e-8
                           && fabs(zeta-ZETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in THexaIsoparametric::GetOrigFromRef(1)" << endl;
    exit(-1);
    return;
  }

  X = xc0 + xc1*xi + xc2*eta + xc3*zeta + xc4*xi*eta
          + xc5*xi*zeta + xc6*eta*zeta + xc7*xi*eta*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta + yc4*xi*eta
          + yc5*xi*zeta + yc6*eta*zeta + yc7*xi*eta*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta + zc4*xi*eta
          + zc5*xi*zeta + zc6*eta*zeta + zc7*xi*eta*zeta;

  for(i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
    Z += ZDistance[i] * FctValues[j][i];
  }
}

/** transfer a set of point from reference to original element */
void THexaIsoparametric::GetOrigFromRef(int N_Points,
                double *xi, double *eta, double *zeta,
                double *X, double *Y, double *Z, double *absdetjk)
{
  int i, j, k;
  double Xi, Eta, Zeta;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  double AuxVector[4*MaxN_BaseFunctions3D];
  TBaseFunct3D *bf;

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    Zeta = zeta[i];

    dx1=xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
    dx2=xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
    dx3=xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
    
    dy1=yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
    dy2=yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
    dy3=yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
    
    dz1=zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    dz2=zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    dz3=zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta + xc4*Xi*Eta
               + xc5*Xi*Zeta + xc6*Eta*Zeta + xc7*Xi*Eta*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta + yc4*Xi*Eta
               + yc5*Xi*Zeta + yc6*Eta*Zeta + yc7*Xi*Eta*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta + zc4*Xi*Eta
               + zc5*Xi*Zeta + zc6*Eta*Zeta + zc7*Xi*Eta*Zeta;

    bf = TFEDatabase3D::GetBaseFunct3D(
              BaseFunctFromOrder[ApproximationOrder]);
    bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
    bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
    bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
    bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);

    for(k=0;k<N_AuxPoints;k++)
    {
      j = IntAux[k];
      X[i] += XDistance[k] * AuxVector[j];
      Y[i] += YDistance[k] * AuxVector[j];
      Z[i] += ZDistance[k] * AuxVector[j];

      dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dx3 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

      dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dy3 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

      dz1 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
      dz2 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
      dz3 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
    }

    detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
                - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
    absdetjk[i] = fabs(detjk);
  }
}

/** transfer from reference element to original element */
void THexaIsoparametric::GetOrigFromRef(double *ref, double *orig)
{
  GetOrigFromRef(ref[0], ref[1], ref[2], orig[0], orig[1], orig[2]);
}


/** transfer from original element to reference element */
void THexaIsoparametric::GetRefFromOrig(double X, double Y, double Z, double &xi, double &eta, double &zeta)
{
  double xt = X - xc0;
  double yt = Y - yc0;
  double zt = Z - zc0;
  double t0, t1, t2;
  double xi0, eta0, zeta0;
  double recdetaffine;
  double eps=1e-14;

  recdetaffine = 1/( xc1*yc2*zc3-xc1*yc3*zc2-yc1*xc2*zc3
                    +yc1*xc3*zc2+zc1*xc2*yc3-zc1*xc3*yc2 );

  t0 = (-(yc2*zc3-yc3*zc2)*xt+(xc2*zc3-xc3*zc2)*yt-(xc2*yc3-xc3*yc2)*zt)
                *recdetaffine;
  t1 = ( (yc1*zc3-yc3*zc1)*xt-(xc1*zc3-xc3*zc1)*yt+(xc1*yc3-xc3*yc1)*zt)
                *recdetaffine;
  t2 = (-(yc1*zc2-yc2*zc1)*xt+(xc1*zc2-xc2*zc1)*yt+(xc2*yc1-xc1*yc2)*zt)
                *recdetaffine;

  xi0 = t0;
  eta0 = t1;
  zeta0 = t2;

  xi = xi0+10;
  eta = eta0+10;
  zeta = zeta0+10;
  while( fabs(xi-xi0)+fabs(eta-eta0)+fabs(zeta-zeta0) > eps)
  {
    xi = xi0;
    eta = eta0;
    zeta = zeta0;

    xt = xc4*xi*eta + xc5*xi*zeta + xc6*eta*zeta + xc7*xi*eta*zeta;
    yt = yc4*xi*eta + yc5*xi*zeta + yc6*eta*zeta + yc7*xi*eta*zeta;
    zt = zc4*xi*eta + zc5*xi*zeta + zc6*eta*zeta + zc7*xi*eta*zeta;

    xi0   = -t0+(-(yc2*zc3-yc3*zc2)*xt+(xc2*zc3-xc3*zc2)*yt-(xc2*yc3-xc3*yc2)*zt)
                *recdetaffine;
    eta0  = -t1+ ( (yc1*zc3-yc3*zc1)*xt-(xc1*zc3-xc3*zc1)*yt+(xc1*yc3-xc3*yc1)*zt)
                *recdetaffine;
    zeta0  = -t2+ (-(yc1*zc2-yc2*zc1)*xt+(xc1*zc2-xc2*zc1)*yt+(xc2*yc1-xc1*yc2)*zt)
                *recdetaffine;
  }
  xi = xi0;
  eta = eta0;
  zeta = zeta0;
}

/** transfer from original element to reference element */
void THexaIsoparametric::GetRefFromOrig(double *orig, double *ref)
{
  OutPut("THexaIsoparametric::GetRefFromOrig is not implemented" << endl);
}

/** calculate functions and derivatives from reference element
    to original element */
void THexaIsoparametric::GetOrigValues(BaseFunct3D BaseFunct,
                               int N_Points, double *xi, double *eta, double *zeta,
                               int N_Functs, QuadFormula3D quadformula)
{
  cout << "THexaIsoparametric::GetOrigValues is not implemented yet" << endl;
  cout << __FILE__ << " " << __LINE__ << endl;
}

void THexaIsoparametric::GetOrigValues(int JointNr, double p1, double p2,
          int N_BaseFunct,
          double *uref, double *uxiref, double *uetaref, double *uzetaref,
          double *uorig, double *uxorig, double *uyorig, double *uzorig)
{
  int i, j, k;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  double Xi, Eta, Zeta;
  TBaseFunct3D *bf;
  bf = TFEDatabase3D::GetBaseFunct3D(BaseFunctFromOrder[ApproximationOrder]);
  double AuxVector[4*MaxN_BaseFunctions3D];

  switch(JointNr)
  {
    case 0:
      Xi = p1; Eta = p2; Zeta = -1;
    break;

    case 1:
      Xi = p2; Eta = -1; Zeta = p1;
    break;

    case 2:
      Xi = 1; Eta = p2; Zeta = p1;
    break;

    case 3:
      Xi = -p2; Eta = 1; Zeta = p1;
    break;

    case 4:
      Xi = -1; Eta = p1; Zeta = p2;
    break;

    case 5:
      Xi = p2; Eta = p1; Zeta = 1;
    break;

    default:
      Error("Wrong local joint number" << endl);
      return;
  }

  dx1=xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
  dx2=xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
  dx3=xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
    
  dy1=yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
  dy2=yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
  dy3=yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
    
  dz1=zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
  dz2=zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
  dz3=zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;

  bf = TFEDatabase3D::GetBaseFunct3D(
            BaseFunctFromOrder[ApproximationOrder]);
  bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
  bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
  bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
  bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);

  for(k=0;k<N_AuxPoints;k++)
  {
    j = IntAux[k];
    dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
    dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
    dx3 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

    dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
    dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
    dy3 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];

    dz1 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D];
    dz2 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*2];
    dz3 += ZDistance[k] * AuxVector[j+MaxN_BaseFunctions3D*3];
  } // endfor k

  detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
         -dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
  rec_detjk=1/detjk;

  for(k=0;k<N_BaseFunct;k++)
  {
    uxorig[k] = ( (dy2*dz3 - dy3*dz2)*uxiref[k]
                   +(dy3*dz1 - dy1*dz3)*uetaref[k]
                   +(dy1*dz2 - dy2*dz1)*uzetaref[k] )*rec_detjk;
    uyorig[k] = ( (dx3*dz2 - dx2*dz3)*uxiref[k]
                   +(dx1*dz3 - dx3*dz1)*uetaref[k] 
                   +(dx2*dz1 - dx1*dz2)*uzetaref[k] )*rec_detjk;
    uzorig[k] = ( (dx2*dy3 - dy2*dx3)*uxiref[k]
                   +(dx3*dy1 - dx1*dy3)*uetaref[k] 
                   +(dx1*dy2 - dx2*dy1)*uzetaref[k] )*rec_detjk;                            
    uorig[k] = uref[k];
  } // endfor k
}

/** calculate functions and derivatives from reference element
    to original element, for all given elements */
void THexaIsoparametric::GetOrigValues(int N_Sets, BaseFunct3D *BaseFuncts,
                               int N_Points, double *xi, double *eta,
                               double *zeta, QuadFormula3D quadformula,
                               bool *Needs2ndDer)
{
  int i,j,k,N_, start, end;
  double **refvaluesD000, **origvaluesD000;
  double **refvaluesD100, **origvaluesD100;
  double **refvaluesD010, **origvaluesD010;
  double **refvaluesD001, **origvaluesD001;
  double **refvaluesD200, **origvaluesD200;
  double **refvaluesD110, **origvaluesD110;
  double **refvaluesD101, **origvaluesD101;
  double **refvaluesD011, **origvaluesD011;
  double **refvaluesD020, **origvaluesD020;
  double **refvaluesD002, **origvaluesD002;
  double *refD000, *origD000;
  double *refD100, *origD100;
  double *refD010, *origD010;
  double *refD001, *origD001;
  double *refD200, *origD200;
  double *refD110, *origD110;
  double *refD101, *origD101;
  double *refD011, *origD011;
  double *refD020, *origD020;
  double *refD002, *origD002;
  double r20, r11, r02, o20, o11, o02;
  double *aux;
  double GeoData[3][3];
  double Eye[3][3];
  BaseFunct3D BaseFunct;
  int N_Functs;
  bool SecondDer;
  int ii,ij,ik;
  double tmp,Eye1[3][3];
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;

  QuadFormula = quadformula;
  
  SecondDer = FALSE;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
    int BaseVectDim = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
    if(BaseVectDim != 1)
    {
      ErrMsg("HexaIsoparametric for vector valued functions not " <<
             "implemented. (Piola transform)");
      exit(0);
    }
      
    refvaluesD000=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D000);
    if(refvaluesD000==NULL)
    {
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)
                ->MakeRefElementData(QuadFormula);
      refvaluesD000=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D000);
    }
      
    origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    if(origvaluesD000==NULL)
    {
      origvaluesD000 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD000[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D000, origvaluesD000);
    }
      
    for(j=0;j<N_Points;j++)
    {
      refD000 = refvaluesD000[j];
      origD000 = origvaluesD000[j];
          
      memcpy(origD000, refD000, N_Functs*SizeOfDouble);
    } // endfor j
      
    refvaluesD100=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D100);
    origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    if(origvaluesD100==NULL)
    {
      origvaluesD100 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD100[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D100, origvaluesD100);
    }
      
    refvaluesD010=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D010);
    origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    if(origvaluesD010==NULL)
    {
      origvaluesD010 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD010[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D010, origvaluesD010);
    }
      
    refvaluesD001=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D001);
    origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    if(origvaluesD001==NULL)
    {
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD001[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D001, origvaluesD001);
    }
  }

/*
  if(Needs2ndDer[i])
  {
    SecondDer = TRUE;
          
    refvaluesD200=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D200);
    origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
    if(origvaluesD200==NULL)
    {
      origvaluesD200 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD200[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D200, origvaluesD200);
    }
          
    refvaluesD110=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D110);
    origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
    if(origvaluesD110==NULL)
    {
      origvaluesD110 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD110[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D110, origvaluesD110);
    }
          
    refvaluesD101=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D101);
    origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
    if(origvaluesD101==NULL)
    {
      origvaluesD101 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD101[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D101, origvaluesD101);
    }
          
    refvaluesD011=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D011);
    origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
    if(origvaluesD011==NULL)
    {
      origvaluesD011 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD011[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D011, origvaluesD011);
    }
          
    refvaluesD020=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D020);
    origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
    if(origvaluesD020==NULL)
    {
      origvaluesD020 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD020[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues
                        (BaseFunct, D020, origvaluesD020);
    }
          
    refvaluesD002=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D002);
    origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
    if(origvaluesD002==NULL)
    {
      origvaluesD002 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD002[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D002, origvaluesD002);
    } // endfor Needs2ndDer[i]
  } // endfor i

*/

  // D100, D010 and D001
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
          
    refvaluesD100=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D100);
    origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
          
    refvaluesD010=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D010);
    origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
          
    refvaluesD001=TFEDatabase3D::GetRefElementValues
                        (BaseFunct, QuadFormula, D001);
    origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
          
    for(j=0;j<N_Points;j++)
    {
      refD100 = refvaluesD100[j];
      origD100 = origvaluesD100[j];
              
      refD010 = refvaluesD010[j];
      origD010 = origvaluesD010[j];
              
      refD001 = refvaluesD001[j];
      origD001 = origvaluesD001[j];

      dx1=xc1 + xc4*eta[j] + xc5*zeta[j] + xc7*eta[j]*zeta[j];
      dx2=xc2 + xc4*xi[j] + xc6*zeta[j] + xc7*xi[j]*zeta[j];
      dx3=xc3 + xc5*xi[j] + xc6*eta[j] + xc7*xi[j]*eta[j];
              
      dy1=yc1 + yc4*eta[j] + yc5*zeta[j] + yc7*eta[j]*zeta[j];
      dy2=yc2 + yc4*xi[j] + yc6*zeta[j] + yc7*xi[j]*zeta[j];
      dy3=yc3 + yc5*xi[j] + yc6*eta[j] + yc7*xi[j]*eta[j];
              
      dz1=zc1 + zc4*eta[j] + zc5*zeta[j] + zc7*eta[j]*zeta[j];
      dz2=zc2 + zc4*xi[j] + zc6*zeta[j] + zc7*xi[j]*zeta[j];
      dz3=zc3 + zc5*xi[j] + zc6*eta[j] + zc7*xi[j]*eta[j];

      for(k=0;k<N_AuxPoints;k++)
      {
        dx1 += XDistance[k] * XiDerValues[j][k];
        dx2 += XDistance[k] * EtaDerValues[j][k];
        dx3 += XDistance[k] * ZetaDerValues[j][k];
  
        dy1 += YDistance[k] * XiDerValues[j][k];
        dy2 += YDistance[k] * EtaDerValues[j][k];
        dy3 += YDistance[k] * ZetaDerValues[j][k];
  
        dz1 += ZDistance[k] * XiDerValues[j][k];
        dz2 += ZDistance[k] * EtaDerValues[j][k];
        dz3 += ZDistance[k] * ZetaDerValues[j][k];
      }
      detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
             -dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
      rec_detjk=1/detjk;
              
      for(k=0;k<N_Functs;k++)
      {
        origD100[k] = ( (dy2*dz3 - dy3*dz2)*refD100[k]
                       +(dy3*dz1 - dy1*dz3)*refD010[k]
                       +(dy1*dz2 - dy2*dz1)*refD001[k] )*rec_detjk;
        origD010[k] = ( (dx3*dz2 - dx2*dz3)*refD100[k]
                       +(dx1*dz3 - dx3*dz1)*refD010[k] 
                       +(dx2*dz1 - dx1*dz2)*refD001[k] )*rec_detjk;
        origD001[k] = ( (dx2*dy3 - dy2*dx3)*refD100[k]
                       +(dx3*dy1 - dx1*dy3)*refD010[k] 
                       +(dx1*dy2 - dx2*dy1)*refD001[k] )*rec_detjk;                            
      } // endfor k
    } // endfor j
  } // endfor i
   
  // leave if no second derivatives are needed
  if(!SecondDer) return;

  // calculate second derivatives
  OutPut("Second derivatives are not implemented yet" << endl);
}

/** calculate functions and derivatives from reference element
    to original element */
void THexaIsoparametric::GetOrigValues(double xi, double eta, double zeta,
                               int N_BaseFunct,
                               double *uref, double *uxiref, double *uetaref, double *uzetaref,
                               double *uorig, double *uxorig, double *uyorig, double *uzorig)
{
  int i,j,k;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  
  // D000
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(fabs(xi-XI[i])<1e-8 && fabs(eta-ETA[i])<1e-8
                           && fabs(zeta-ZETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in THexaIsoparametric::GetOrigFromRef(3)" << endl;
    exit(-1);
    return;
  }

  dx1=xc1 + xc4*eta + xc5*zeta + xc7*eta*zeta;
  dx2=xc2 + xc4*xi + xc6*zeta + xc7*xi*zeta;
  dx3=xc3 + xc5*xi + xc6*eta + xc7*xi*eta;
              
  dy1=yc1 + yc4*eta + yc5*zeta + yc7*eta*zeta;
  dy2=yc2 + yc4*xi + yc6*zeta + yc7*xi*zeta;
  dy3=yc3 + yc5*xi + yc6*eta + yc7*xi*eta;
  
  dz1=zc1 + zc4*eta + zc5*zeta + zc7*eta*zeta;
  dz2=zc2 + zc4*xi + zc6*zeta + zc7*xi*zeta;
  dz3=zc3 + zc5*xi + zc6*eta + zc7*xi*eta;
  
  for(k=0;k<N_AuxPoints;k++)
  {
    dx1 += XDistance[k] * XiDerValues[j][k];
    dx2 += XDistance[k] * EtaDerValues[j][k];
    dx3 += XDistance[k] * ZetaDerValues[j][k];

    dy1 += YDistance[k] * XiDerValues[j][k];
    dy2 += YDistance[k] * EtaDerValues[j][k];
    dy3 += YDistance[k] * ZetaDerValues[j][k];

    dz1 += ZDistance[k] * XiDerValues[j][k];
    dz2 += ZDistance[k] * EtaDerValues[j][k];
    dz3 += ZDistance[k] * ZetaDerValues[j][k];
  }
  detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
         -dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
  rec_detjk=1/detjk;

  // D100, D010 and D001
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((dy2*dz3 - dy3*dz2)*uxiref[i] +
                 (dy3*dz1 - dy1*dz3)*uetaref[i] +
                 (dy1*dz2 - dy2*dz1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((dx3*dz2 - dx2*dz3)*uxiref[i] +
                 (dx1*dz3 - dx3*dz1)*uetaref[i] +
                 (dx2*dz1 - dx1*dz2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((dx2*dy3 - dx3*dy2)*uxiref[i] +
                 (dx3*dy1 - dx1*dy3)*uetaref[i] + 
                 (dx1*dy2 - dx2*dy1)*uzetaref[i]) *rec_detjk;
  } // endfor i
}

void THexaIsoparametric::SetCell(TBaseCell *cell)
{
  int i, j, k, l;
  TJoint *joint;
  JointType type;
  double xp1, xp2, xp3, xp4, yp1, yp2, yp3, yp4, zp1, zp2, zp3, zp4;
  TShapeDesc *ShapeDesc;
  const int *TmpFaceVertex, *TmpLen;
  int MaxLen;
  int *JointDOF;
  double Param1[4], Param2[4];
  TBoundFace *bdface;
  TFEDesc3D *fedesc;
  TBaseFunct3D *bf;
  double dt, ds;
  double T, S, t[4],s[4], xm, ym, zm, xp, yp, zp;
  TVertex **Vertices;
  TVertex **AuxVertices;
  double x[4], y[4], z[4], factor;
  double LinComb[4];
  int CurvedJoint;

  N_AuxPoints = 0;
  Cell = cell;
  Vertices = ((TGridCell*)Cell)->GetVertices();
  ShapeDesc = Cell->GetShapeDesc();
  ShapeDesc->GetFaceVertex(TmpFaceVertex, TmpLen, MaxLen);

  TFEDatabase3D::GetQuadFormula3D(QuadFormula)
        ->GetFormulaData(N_QuadPoints, W, XI, ETA, ZETA);

  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);
  Cell->GetVertex(4)->GetCoords(x4, y4, z4);
  Cell->GetVertex(5)->GetCoords(x5, y5, z5);
  Cell->GetVertex(6)->GetCoords(x6, y6, z6);
  Cell->GetVertex(7)->GetCoords(x7, y7, z7);

  xc0 = ( x0 + x1 + x3 + x2 + x4 + x5 + x7 + x6) * 0.125;
  xc1 = (-x0 + x1 - x3 + x2 - x4 + x5 - x7 + x6) * 0.125;
  xc2 = (-x0 - x1 + x3 + x2 - x4 - x5 + x7 + x6) * 0.125;
  xc3 = (-x0 - x1 - x3 - x2 + x4 + x5 + x7 + x6) * 0.125;
  xc4 = ( x0 - x1 - x3 + x2 + x4 - x5 - x7 + x6) * 0.125;
  xc5 = ( x0 - x1 + x3 - x2 - x4 + x5 - x7 + x6) * 0.125;
  xc6 = ( x0 + x1 - x3 - x2 - x4 - x5 + x7 + x6) * 0.125;
  xc7 = (-x0 + x1 + x3 - x2 + x4 - x5 - x7 + x6) * 0.125;

  yc0 = ( y0 + y1 + y3 + y2 + y4 + y5 + y7 + y6) * 0.125;
  yc1 = (-y0 + y1 - y3 + y2 - y4 + y5 - y7 + y6) * 0.125;
  yc2 = (-y0 - y1 + y3 + y2 - y4 - y5 + y7 + y6) * 0.125;
  yc3 = (-y0 - y1 - y3 - y2 + y4 + y5 + y7 + y6) * 0.125;
  yc4 = ( y0 - y1 - y3 + y2 + y4 - y5 - y7 + y6) * 0.125;
  yc5 = ( y0 - y1 + y3 - y2 - y4 + y5 - y7 + y6) * 0.125;
  yc6 = ( y0 + y1 - y3 - y2 - y4 - y5 + y7 + y6) * 0.125;
  yc7 = (-y0 + y1 + y3 - y2 + y4 - y5 - y7 + y6) * 0.125;

  zc0 = ( z0 + z1 + z3 + z2 + z4 + z5 + z7 + z6) * 0.125;
  zc1 = (-z0 + z1 - z3 + z2 - z4 + z5 - z7 + z6) * 0.125;
  zc2 = (-z0 - z1 + z3 + z2 - z4 - z5 + z7 + z6) * 0.125;
  zc3 = (-z0 - z1 - z3 - z2 + z4 + z5 + z7 + z6) * 0.125;
  zc4 = ( z0 - z1 - z3 + z2 + z4 - z5 - z7 + z6) * 0.125;
  zc5 = ( z0 - z1 + z3 - z2 - z4 + z5 - z7 + z6) * 0.125;
  zc6 = ( z0 + z1 - z3 - z2 - z4 - z5 + z7 + z6) * 0.125;
  zc7 = (-z0 + z1 + z3 - z2 + z4 - z5 - z7 + z6) * 0.125;

  // OutPut(Cell->GetClipBoard() << endl);
  for(i=0;i<6;i++)
  {
    joint = Cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryFace || type == IsoBoundFace || type == IsoJointEqN)
    {
      fedesc = TFEDatabase3D::GetFEDesc3D(
                FEDescFromOrder[ApproximationOrder]);
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);

      JointDOF = fedesc->GetJointDOF(i);

      if(type == BoundaryFace || type == IsoBoundFace)
      {
        bdface = (TBoundFace*)joint;
        bdface->GetParameters(Param1, Param2);
      }
      
      dt = ds = 1.0/ApproximationOrder;
      Vertices[TmpFaceVertex[i*MaxLen+0]]->GetCoords(x[0],y[0],z[0]);
      Vertices[TmpFaceVertex[i*MaxLen+1]]->GetCoords(x[1],y[1],z[1]);
      Vertices[TmpFaceVertex[i*MaxLen+2]]->GetCoords(x[2],y[2],z[2]);
      Vertices[TmpFaceVertex[i*MaxLen+3]]->GetCoords(x[3],y[3],z[3]);

      if(type == IsoBoundFace)
        AuxVertices = ((TIsoBoundFace *)joint)->GetVertices();

      if(type == IsoJointEqN)
        AuxVertices = ((TIsoJointEqN *)joint)->GetVertices();

      for(j=0;j<=ApproximationOrder;j++)
      {
        for(k=0;k<=ApproximationOrder;k++)
        {
          LinComb[0] = (1-j*dt)*(1-k*ds);
          LinComb[1] = (j*dt)*(1-k*ds);
          LinComb[2] = (j*dt)*(k*ds);
          LinComb[3] = (1-j*dt)*(k*ds);

          xm = 0; ym = 0; zm = 0;
          for(l=0;l<4;l++)
          {
            factor = LinComb[l];
            // cout << "LinComb: " << LinComb[l] << endl;
            xm += factor*x[l]; 
            ym += factor*y[l]; 
            zm += factor*z[l]; 
            // cout << x[l] << " " << y[l] << " " << z[l] << endl;
          }

          if(type == BoundaryFace || ApproximationOrder == 1)
          {
            if(type == BoundaryFace || type == IsoBoundFace)
              bdface->GetBoundComp()->GetXYZandTS(4, LinComb, x, y, z,
                                                  Param1, Param2,
                                                  xp, yp, zp, T, S);
            else
            {
              xp = xm; yp = ym; zp = zm;
            }
          }
          else
          {
            if(type == IsoJointEqN)
            {
              if(joint->GetNeighbour(0) == cell)
              {
                AuxVertices[k*(ApproximationOrder+1)+j]->GetCoords(xp, yp, zp);
              }
              else
              {
                switch(joint->GetMapType())
                {
                  case 0:
                    AuxVertices[j*(ApproximationOrder+1)
                                +k]
                                ->GetCoords(xp, yp, zp);
                  break;
                  case 1:
                    AuxVertices[k*(ApproximationOrder+1)
                                +(ApproximationOrder-j)]
                                ->GetCoords(xp, yp, zp);
                  break;
                  case 2:
                    AuxVertices[(ApproximationOrder-j)*(ApproximationOrder+1)
                                +(ApproximationOrder-k)]
                                ->GetCoords(xp, yp, zp);
                  break;
                  case 3:
                    AuxVertices[(ApproximationOrder-k)*(ApproximationOrder+1)
                                +j]
                                ->GetCoords(xp, yp, zp);
                  break;
                } // endswitch
              } // endelse
            } // endelse
            else
            {
              AuxVertices[k*(ApproximationOrder+1)+j]->GetCoords(xp, yp, zp);
            }
          } // endelse

          if(fabs(xp-xm) > 1e-8 || fabs(yp-ym) > 1e-8 || fabs(zp-zm) > 1e-8)
          {
            // cout << "j: " << j << " k: " << k << endl;
            // cout << "I: " << N_AuxPoints << " T:" << T << " S:" << S << endl;
            // cout << "Straight: " << xm << " " << ym << " " << zm << endl;
            // cout << "Boundary: " << xp << " " << yp << " " << zp << endl;
            // cout << "Diff:" << xp-xm << " " << yp-ym << " " << zp-zm << endl;
            XDistance[N_AuxPoints] = xp - xm;
            YDistance[N_AuxPoints] = yp - ym;
            ZDistance[N_AuxPoints] = zp - zm;
            IntAux[N_AuxPoints] = JointDOF[k*(ApproximationOrder+1)+j];
            // cout << "jd: " << JointDOF[k*(ApproximationOrder+1)+j] << endl;
            N_AuxPoints++;

            CurvedJoint = i;
          } // endif
        } // endfor k
      } // endfor j
    } // endif
  } // endfor i

  // cout << "cell nr: " << cell->GetClipBoard() << " " << N_AuxPoints << endl;

  if(N_AuxPoints)
  {
    if(ApproximationOrder == 2)
    {
      XDistance[N_AuxPoints] = XDistance[0]*0.5;
      YDistance[N_AuxPoints] = YDistance[0]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[0]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*0.5;
      YDistance[N_AuxPoints] = YDistance[1]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[1]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*0.5;
      YDistance[N_AuxPoints] = YDistance[2]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[2]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*0.5;
      YDistance[N_AuxPoints] = YDistance[3]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[3]*0.5;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*0.5;
      YDistance[N_AuxPoints] = YDistance[4]*0.5;
      ZDistance[N_AuxPoints] = ZDistance[4]*0.5;
      N_AuxPoints++;

      N_AuxPoints -= 5;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] + 9;
          N_AuxPoints++;
        break;

        case 5:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 5] - 9;
          N_AuxPoints++;
        break;

        /*
        default:
          OutPut("This case has not been implemented yet!" << endl);
        break;
        */
      } 
    } // ApproximationOrder == 2

    if(ApproximationOrder == 3)
    {
      XDistance[N_AuxPoints] = XDistance[0]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[0]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[1]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[2]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[3]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[3]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[4]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[4]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[5]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[5]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[6]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[6]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[7]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[7]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[8]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[8]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[9]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[9]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[10]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[10]*2.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[11]*2.0/3;
      ZDistance[N_AuxPoints] = ZDistance[11]*2.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[0]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[1]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[2]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[3]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[3]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[4]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[4]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[5]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[5]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[6]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[6]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[7]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[7]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[8]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[8]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[9]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[9]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[10]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[10]*1.0/3;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[11]*1.0/3;
      ZDistance[N_AuxPoints] = ZDistance[11]*1.0/3;
      N_AuxPoints++;

      N_AuxPoints -= 24;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints - 12] + 16;
          N_AuxPoints++;
        break;
        
        default:
          OutPut("This case has not been implemented yet!" << endl);
        break;
      }
    }  // ApproximationOrder == 3

    if(ApproximationOrder == 4)
    {
      XDistance[N_AuxPoints] = XDistance[0]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[0]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[1]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[2]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[3]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[3]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[4]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[4]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[5]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[5]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[6]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[6]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[7]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[7]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[8]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[8]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[9]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[9]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[10]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[10]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[11]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[11]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[12]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[12]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[12]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[13]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[13]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[13]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[14]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[14]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[14]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[15]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[15]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[15]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[16]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[16]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[16]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[17]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[17]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[17]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[18]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[18]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[18]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[19]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[19]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[19]*3.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[20]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[20]*3.0/4;
      ZDistance[N_AuxPoints] = ZDistance[20]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[0]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[1]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[2]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[3]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[3]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[4]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[4]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[5]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[5]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[6]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[6]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[7]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[7]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[8]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[8]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[9]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[9]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[10]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[10]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[11]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[11]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[12]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[12]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[12]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[13]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[13]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[13]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[14]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[14]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[14]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[15]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[15]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[15]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[16]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[16]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[16]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[17]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[17]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[17]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[18]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[18]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[18]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[19]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[19]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[19]*2.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[20]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[20]*2.0/4;
      ZDistance[N_AuxPoints] = ZDistance[20]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[0]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[1]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[1]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[2]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[2]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[3]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[3]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[3]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[4]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[4]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[4]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[5]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[5]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[5]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[6]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[6]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[6]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[7]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[7]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[7]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[8]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[8]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[8]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[9]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[9]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[9]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[10]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[10]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[10]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[11]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[11]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[11]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[12]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[12]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[12]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[13]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[13]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[13]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[14]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[14]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[14]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[15]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[15]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[15]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[16]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[16]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[16]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[17]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[17]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[17]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[18]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[18]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[18]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[19]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[19]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[19]*1.0/4;
      N_AuxPoints++;
      XDistance[N_AuxPoints] = XDistance[20]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[20]*1.0/4;
      ZDistance[N_AuxPoints] = ZDistance[20]*1.0/4;
      N_AuxPoints++;

      N_AuxPoints -= 63;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = IntAux[N_AuxPoints-21] + 25;
          N_AuxPoints++;

        break;
        
        default:
          OutPut("This case has not been implemented yet!" << endl);
        break;
      }
    }  // ApproximationOrder == 4
  } // N_AuxPoints > 0

  if(N_AuxPoints)
  {
    for(i=0;i<N_QuadPoints;i++)
    {
      bf->GetDerivatives(D000, XI[i], ETA[i], ZETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        FctValues[i][j] = DoubleAux[IntAux[j]]; 

      bf->GetDerivatives(D100, XI[i], ETA[i], ZETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        XiDerValues[i][j] = DoubleAux[IntAux[j]]; 

      bf->GetDerivatives(D010, XI[i], ETA[i], ZETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        EtaDerValues[i][j] = DoubleAux[IntAux[j]]; 

      bf->GetDerivatives(D001, XI[i], ETA[i], ZETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        ZetaDerValues[i][j] = DoubleAux[IntAux[j]]; 
    } // endfor i
  } // endif
}

/** return outer normal unit vector */
void THexaIsoparametric::GetOuterNormal(int j, double s, double t,
                                        double &n1, double &n2, double &n3)
{
  double t11, t12, t13, t21, t22, t23;
//   double len;

  GetTangentVectors(j, s, t, t11, t12, t13, t21, t22, t23);
  
  n1 = t12*t23 - t13*t22;
  n2 = t13*t21 - t11*t23;
  n3 = t11*t22 - t12*t21;

  
  // without scaling needed for finding area, Sashi, 20. Jan 2012  
//   len = sqrt(n1*n1 + n2*n2 + n3*n3);
// 
//   n1 /= len;
//   n2 /= len;
//   n3 /= len;
}

/** return two tangent vectors */
void THexaIsoparametric::GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23)
{
  double Xi, Eta, Zeta;
  int i, k;
  double AuxVector[4*MaxN_BaseFunctions3D];
  TBaseFunct3D *bf;

  switch(j)
  {
    case 0:
      Xi = p1; Eta = p2; Zeta = -1;
      t11 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t12 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t13 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;

      t21 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t22 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t23 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);
      bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);
    
      for(k=0;k<N_AuxPoints;k++)
      {
        i = IntAux[k];
        t21 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t11 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
    
        t22 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t12 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
    
        t23 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t13 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
      } // endfor k
    break;

    case 1:
      Xi = p2; Eta = -1; Zeta = p1;
      t11 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t12 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t13 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;

      t21 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t22 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t23 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);
      bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);
    
      for(k=0;k<N_AuxPoints;k++)
      {
        i = IntAux[k];
        t11 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t21 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t12 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t22 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t13 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t23 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
      } // endfor k
    break;

    case 2:
      Xi = 1; Eta = p2; Zeta = p1;
      t11 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t12 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t13 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;

      t21 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t22 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t23 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);
      bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);
    
      for(k=0;k<N_AuxPoints;k++)
      {
        i = IntAux[k];
        t11 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
        t21 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t12 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
        t22 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t13 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
        t23 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
      } // endfor k
    break;

    case 3:
      Xi = -p2; Eta = 1; Zeta = p1;
      t11 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t12 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t13 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;

      t21 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t22 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t23 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);
      bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);
    
      for(k=0;k<N_AuxPoints;k++)
      {
        i = IntAux[k];
        t21 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t11 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t22 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t12 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t23 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t13 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
      } // endfor k
    break;

    case 4:
      Xi = -1; Eta = p1; Zeta = p2;
      t11 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t12 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t13 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
      
      t21 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t22 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t23 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);
      bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);
    
      for(k=0;k<N_AuxPoints;k++)
      {
        i = IntAux[k];
        t21 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
        t11 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t22 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
        t12 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
    
        t23 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
        t13 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*3];
      } // endfor k
    break;

    case 5:
      Xi = p2; Eta = p1; Zeta = 1;
      t11 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t12 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t13 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;

      t21 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t22 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t23 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);
      bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+3*MaxN_BaseFunctions3D);
    
      for(k=0;k<N_AuxPoints;k++)
      {
        i = IntAux[k];
        t11 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t21 += XDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
    
        t12 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t22 += YDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
    
        t13 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D];
        t23 += ZDistance[k] * AuxVector[i+MaxN_BaseFunctions3D*2];
      } // endfor k
    break;

    default:
      Error("Wrong local joint number" << endl);
      return;
  }
}
