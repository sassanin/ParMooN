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
// @(#)TetraIsoparametric.C      1.6 02/22/00
//
// Class:      TTetraIsoparametric
//
// Purpose:    trilinear reference transformations for Tetrahedron
//
// Author:     Gunar Matthies
//
// History:    24.01.01 start implementation
// 
// =======================================================================

#include <TetraIsoparametric.h>
#include <BoundFace.h>
#include <IsoBoundFace.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>
#include <MooNMD_Io.h>
#include <Database.h>
// 
//  BaseFunct3D TTetraIsoparametric::BaseFunctFromOrder[] = { 
//                 BF_C_T_P0_3D, BF_C_T_P1_3D, BF_C_T_B2_3D, BF_C_T_P3_3D };
 BaseFunct3D TTetraIsoparametric::BaseFunctFromOrder[] = { 
                BF_C_T_P0_3D, BF_C_T_P1_3D, BF_C_T_P2_3D, BF_C_T_P3_3D };

//  FEDesc3D TTetraIsoparametric::FEDescFromOrder[] = { 
//                 FE_C_T_P0_3D, FE_C_T_P1_3D, FE_C_T_B2_3D, FE_C_T_P3_3D };
 FEDesc3D TTetraIsoparametric::FEDescFromOrder[] = { 
                FE_C_T_P0_3D, FE_C_T_P1_3D, FE_C_T_P2_3D, FE_C_T_P3_3D };

/** constuctor */
TTetraIsoparametric::TTetraIsoparametric()
{
}

/** transfer a point from reference face to original element face */
void TTetraIsoparametric::GetOrigBoundFromRef(int Joint, double xi, double eta,
					      double &X, double &Y, double &Z)
{
  int l, k;
  double Xi, Eta, Zeta;
  double AuxVector[MaxN_BaseFunctions3D];
  TBaseFunct3D *bf;

  bf = TFEDatabase3D::GetBaseFunct3D(
            BaseFunctFromOrder[ApproximationOrder]);
  
  switch(Joint)
  {
    case 0:
      Xi = xi;
      Eta = eta;
      Zeta = 0.;
      break;

    case 1:
      Xi =  eta;
      Eta = 0. ;
      Zeta = xi;
      break;

    case 2:
      Xi = xi;
      Eta = 1.- xi - eta;
      Zeta = eta;
      break;

    case 3:
      Xi = 0.;
      Eta = xi;
      Zeta = eta;
      break;
      
    default:
      Error("Wrong joint number for tetrahedron! " << endl);
      return;
  }
      
  X = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta;
  Y = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta;
  Z = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta;
  
  bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);
  
  for(l=0;l<N_AuxPoints;l++)
  {
    k = IntAux[l];
    X += XDistance[l] * AuxVector[k];
    Y += YDistance[l] * AuxVector[k];
    Z += ZDistance[l] * AuxVector[k];
  }
      
}

/** transfer a set of points from reference face to original element face */
void TTetraIsoparametric::GetOrigBoundFromRef(int Joint, int N_Points, double *xi, double *eta,
                                              double *X, double *Y, double *Z)
{
  int i, k, l, temp;
  double Xi, Eta, Zeta;
  double AuxVector[MaxN_BaseFunctions3D];
  TBaseFunct3D *bf;


  bf = TFEDatabase3D::GetBaseFunct3D(
            BaseFunctFromOrder[ApproximationOrder]);

   temp = Joint*MaxN_BaseFunctions3D;
  for(i=0;i<N_Points;i++)
  {
    switch(Joint)
     {
      case 0:
                 Xi = xi[i];
                 Eta = eta[i];
                 Zeta = 0.;
      break;

      case 1:
                 Xi =  eta[i];
                 Eta = 0. ;
                 Zeta = xi[i];
    break;

    case 2:
                 Xi = xi[i];
                 Eta = 1.- xi[i] - eta[i];
                 Zeta = eta[i];
    break;

    case 3:
                 Xi = 0.;
                 Eta = xi[i];
                 Zeta = eta[i];
    break;

    default:
      Error("Wrong joint number for tetrahedron! " << endl);
      return;
  }

    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta;


    bf->GetDerivatives(D000, Xi, Eta, Zeta, AuxVector);

    for(l=0;l<N_AuxPoints;l++)
     {
       k = IntAux[l];
       X[i] += XDistance[l] * AuxVector[k];
       Y[i] += YDistance[l] * AuxVector[k];
       Z[i] += ZDistance[l] * AuxVector[k];
     }
  } // endfor i
}

/** transfer from reference element to original element */
void TTetraIsoparametric::GetOrigFromRef(double xi, double eta, double zeta,
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
    Error("error in TTetraIsoparametric::GetOrigFromRef(1)" << endl);
    exit(-1);
    return;
  }

  X = xc0 + xc1*xi + xc2*eta + xc3*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta;
	    
  for(i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
    Z += ZDistance[i] * FctValues[j][i];
  }
  
  // cout << xi << " " << eta << " " << zeta << endl;
  // cout << X << " " << Y << " " << Z << endl;
  // cout << "1-x: " << X << " y: " << Y << " z: " << Z << endl;
}

/** transfer a set of point from reference to original element */
void TTetraIsoparametric::GetOrigFromRef(int N_Points,
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

    dx1 = xc1;
    dx2 = xc2;
    dx3 = xc3;
    
    dy1 = yc1;
    dy2 = yc2;
    dy3 = yc3;
    
    dz1 = zc1;
    dz2 = zc2;
    dz3 = zc3;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta;

    if (N_AuxPoints)
    {
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
    }

    // cout << "2-x: " << X[i] << " y: " << Y[i] << " z: " << Z[i] << endl;
    // cout << "3-x: " << Xi << " " << Eta << " " << Zeta;
    // cout << " " << dx1 << " y: " << dx2 << " z: " << dx3 << endl;

    detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
                - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
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
void TTetraIsoparametric::GetOrigFromRef(double *ref, double *orig)
{
  GetOrigFromRef(ref[0], ref[1], ref[2], orig[0], orig[1], orig[2]);
}


/** transfer from original element to reference element */
void TTetraIsoparametric::GetRefFromOrig(double X, double Y, double Z, double &xi, double &eta, double &zeta)
{

  OutPut("TTetraIsoparametric::GetRefFromOrig is not implemented" << endl);
  OutPut("TTetraIsoparametric::GetRefFromOrig works only for corners of cells" << endl);
  double xt = X - xc0;
  double yt = Y - yc0;
  double zt = Z - zc0;
  double xi0, eta0, zeta0;
  double recdetaffine;
  double eps=1e-14;

  recdetaffine = 1/( xc1*yc2*zc3-xc1*yc3*zc2-yc1*xc2*zc3
                    +yc1*xc3*zc2+zc1*xc2*yc3-zc1*xc3*yc2 );

  xi = (-(yc2*zc3-yc3*zc2)*xt+(xc2*zc3-xc3*zc2)*yt-(xc2*yc3-xc3*yc2)*zt)
                *recdetaffine;
  eta = ( (yc1*zc3-yc3*zc1)*xt-(xc1*zc3-xc3*zc1)*yt+(xc1*yc3-xc3*yc1)*zt)
                *recdetaffine;
  zeta = (-(yc1*zc2-yc2*zc1)*xt+(xc1*zc2-xc2*zc1)*yt+(xc2*yc1-xc1*yc2)*zt)
                *recdetaffine;
}

/** transfer from original element to reference element */
void TTetraIsoparametric::GetRefFromOrig(double *orig, double *ref)
{
  OutPut("TTetraIsoparametric::GetRefFromOrig is not implemented" << endl);
}

/** calculate functions and derivatives from reference element
    to original element */
void TTetraIsoparametric::GetOrigValues(BaseFunct3D BaseFunct,
                               int N_Points, double *xi, double *eta, double *zeta,
                               int N_Functs, QuadFormula3D quadformula)
{
  OutPut("TTetraIsoparametric::GetOrigValues is not implemented yet" << endl);
}

/** calculate functions and derivatives from reference element
    to original element, for all given elements */
void TTetraIsoparametric::GetOrigValues(int N_Sets, BaseFunct3D *BaseFuncts,
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

      dx1 = xc1;
      dx2 = xc2;
      dx3 = xc3;
              
      dy1 = yc1;
      dy2 = yc2;
      dy3 = yc3;
              
      dz1 = zc1;
      dz2 = zc2;
      dz3 = zc3;

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
      } // endfor k
      detjk = dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2
             -dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
      rec_detjk=1/detjk;

      // cout << "derivatives: " << detjk << endl;
      // cout << "4: " << dx1 << " " << dx2 << " " << dx3 << endl;
      // cout << "4: " << dy1 << " " << dy2 << " " << dy3 << endl;
      // cout << "4: " << dz1 << " " << dz2 << " " << dz3 << endl;

      // if (N_AuxPoints)
      //  OutPut("aux " << N_AuxPoints << endl);
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

        if (fabs(origD001[k])<1e-10)
          origD001[k] = 0;
        if (fabs(origD010[k])<1e-10)
          origD010[k] = 0;
        if (fabs(origD100[k])<1e-10)
          origD100[k] = 0;
        //if (!N_AuxPoints)
        //  OutPut(k << " " << origD100[k] << " " << origD010[k] << " " <<  origD001[k] << endl);
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
void TTetraIsoparametric::GetOrigValues(double xi, double eta, double zeta,
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
    Error(__FILE__ << ":" << __LINE__ << ": error in TTetraIsoparametric::GetOrigValues(3)" << endl);
    exit(-1);
    return;
  }

  dx1 = xc1;
  dx2 = xc2;
  dx3 = xc3;
              
  dy1 = yc1;
  dy2 = yc2;
  dy3 = yc3;
  
  dz1 = zc1;
  dz2 = zc2;
  dz3 = zc3;

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

  // cout << "derivatives: " << detjk << endl;
  // cout << "5: " << dx1 << " " << dx2 << " " << dx3 << endl;
  // cout << "5: " << dy1 << " " << dy2 << " " << dy3 << endl;
  // cout << "5: " << dz1 << " " << dz2 << " " << dz3 << endl;

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


void TTetraIsoparametric::GetOrigValues(int JointNr, double p1, double p2,
          int N_BaseFunct,
          double *uref, double *uxiref, double *uetaref, double *uzetaref,
          double *uorig, double *uxorig, double *uyorig, double *uzorig)
{

  int i, j, k, temp;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  double Xi, Eta, Zeta;
  TBaseFunct3D *bf;
  double AuxVector[4*MaxN_BaseFunctions3D];

//   OutPut("GetOrigValues() on joint"<<endl);
  
  bf = TFEDatabase3D::GetBaseFunct3D(
            BaseFunctFromOrder[ApproximationOrder]);

  temp = JointNr*MaxN_BaseFunctions3D;

  switch(JointNr)
  {
    case 0:
      Xi = p1; Eta = p2; Zeta = 0;

    break;

    case 1:
      Xi = p2; Eta = 0; Zeta = p1;

    break;

    case 2:
      Xi = p1; Eta = 1.-p1- p2; Zeta = p2;

    break;

    case 3:
      Xi = 0; Eta = p1; Zeta = p2;

    break;

    default:
      Error("Wrong joint number for tetrahedron! " << endl);
      return;
  }

  dx1 = xc1;
  dx2 = xc2;
  dx3 = xc3;

  dy1 = yc1;
  dy2 = yc2;
  dy3 = yc3;

  dz1 = zc1;
  dz2 = zc2;
  dz3 = zc3;


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
                - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
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

// cout<< " uxiref[k] " << uxiref[k] << " " <<uetaref[k]<< " " << uzetaref[k]<< endl;
  } // endfor k

//   Error("Check xi, eta, zeta values for tetrahedra!" << endl);
// exit(-1);
}

void TTetraIsoparametric::SetCell(TBaseCell *cell)
{
  int i, ii, j, k, l, m, n, auxvert, N_AuxVertices=0, N_IsoFaces=0;
  TJoint *joint;
  JointType type;
  TShapeDesc *ShapeDesc;
  const int *TmpFaceVertex, *TmpLen;
  int MaxLen;
  int *JointDOF, LocIsoDOF[4];
  double Param1[4], Param2[4];
  TBoundFace *bdface;
  TFEDesc3D *fedesc;
  TBaseFunct3D *bf;
  double dt, ds;
  double T, S, t[3], s[3], xm, ym, zm, xp, yp, zp;
  TVertex **Vertices, **AuxVertices;
  double x[3], y[3], z[3], factor;
  double LinComb[3][3];
  double a,b,c, project_a, project_b, project_c;

  N_AuxPoints = 0;

  Vertices = ((TGridCell*)cell)->GetVertices();
  ShapeDesc = cell->GetShapeDesc();
  ShapeDesc->GetFaceVertex(TmpFaceVertex, TmpLen, MaxLen);

  TFEDatabase3D::GetQuadFormula3D(QuadFormula)
        ->GetFormulaData(N_QuadPoints, W, XI, ETA, ZETA);

  cell->GetVertex(0)->GetCoords(x0, y0, z0);
  cell->GetVertex(1)->GetCoords(x1, y1, z1);
  cell->GetVertex(2)->GetCoords(x2, y2, z2);
  cell->GetVertex(3)->GetCoords(x3, y3, z3);

  xc0 = x0;
  xc1 = x1-x0;
  xc2 = x2-x0;
  xc3 = x3-x0;

  yc0 = y0;
  yc1 = y1-y0;
  yc2 = y2-y0;
  yc3 = y3-y0;

  zc0 = z0;
  zc1 = z1-z0;
  zc2 = z2-z0;
  zc3 = z3-z0;

  N_IsoFaces=0;
  for(i=0;i<4;i++) // run through all faces
  {
    joint = cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryFace || type == IsoBoundFace)
     {
      fedesc = TFEDatabase3D::GetFEDesc3D(
                FEDescFromOrder[ApproximationOrder]);
      bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);

      JointDOF = fedesc->GetJointDOF(i);

      bdface = (TBoundFace*)joint;
//       bdface->GetParameters(Param1, Param2);

//       // do not apply moving of nodes on the bottom and the
//       // top of a sandwich mesh
//       if (bdface->GetBoundComp()->GetID()>=1000)
//         continue;

//       dt = ds = 1.0/ApproximationOrder;
      // coordinates of vertices of the face
      Vertices[TmpFaceVertex[i*MaxLen+0]]->GetCoords(x[0],y[0],z[0]);
      Vertices[TmpFaceVertex[i*MaxLen+1]]->GetCoords(x[1],y[1],z[1]);
      Vertices[TmpFaceVertex[i*MaxLen+2]]->GetCoords(x[2],y[2],z[2]);

       if(type == IsoBoundFace)
         {
           AuxVertices = ((TIsoBoundFace *)joint)->GetVertices();
           N_AuxVertices = ((TIsoBoundFace *)joint)->GetN_Vertices();
	       N_AuxVertices -= N_AuxVertices > 3 ? 1 : 0; // face bubble not needed
           N_IsoFaces++;
           if(N_IsoFaces>1)
           {
             OutPut("More than one isoface are in a sigle cell " << N_IsoFaces << endl);
             OutPut("Check  TTetraIsoparametric::SetCell" << endl);
             OutPut("x[0]: " << x[0] << " y[0]: " << y[0] << " z[0]: " << z[0] <<endl);
             OutPut("x[1]: " << x[1] << " y[1]: " << y[1] << " z[1]: " << z[1] <<endl);
             OutPut("x[2]: " << x[2] << " y[2]: " << y[2] << " z[2]: " << z[2] <<endl);
             exit(0);
           }
         }

     if(TDatabase::ParamDB->MOVING_BOUNDARY)
      {  
        if(ApproximationOrder==2)
         {
          if(N_AuxVertices!=3)
           {
            cout<< " ApproximationOrder " <<ApproximationOrder<< " " << N_AuxVertices << " " <<TDatabase::ParamDB->MOVING_BOUNDARY<< endl;
            OutPut("Moving IsoBoundFace, but no face vertices are generated "<<endl;);
            exit(0);
           }
              
              
           //      additional vertices should have been already added in this face in the main program
           LocIsoDOF[0] = 1;
           LocIsoDOF[1] = 3;
           LocIsoDOF[2] = 4;
           //         LocIsoDOF[3] = 6;      // local dof of face point, see C_T_B2_3D_J

           LinComb[0][0]=0.5; LinComb[1][0]=0.5; LinComb[2][0]=0.0;
           LinComb[0][1]=0.5; LinComb[1][1]=0.0; LinComb[2][1]=0.5;
           LinComb[0][2]=0.0; LinComb[1][2]=0.5; LinComb[2][2]=0.5;

           for(j=0;j<N_AuxVertices;j++)
            {
             AuxVertices[j]->GetCoords(xp, yp, zp);
             xm = 0.;  ym=0.; zm=0.;
             for(k=0;k<3;k++)    // number of vertices on the face
              {
               factor = LinComb[j][k];
               xm += factor*x[k];
               ym += factor*y[k];
               zm += factor*z[k];
              }
 
            if(fabs(xp-xm) > 1e-8 || fabs(yp-ym) > 1e-8 || fabs(zp-zm) > 1e-8)
             {

//             cout << endl;
//             cout << "xpzp " << xp*xp + yp*yp +  zp*zp<< endl;
//             cout << "xmzm " << xm*xm + ym*ym +  zm*zm<< endl;
//             OutPut("xp: " << xp << " xm: " << xm );
//             cout << "X: " << xp-xm << endl;
//             OutPut(" yp: " << yp << " ym: " << ym);
//             cout << "Y: " << yp-ym << endl;
//             OutPut(" zp: " << zp << " zm: " << zm );
//             cout << "Z: " << zp-zm << endl;
//            OutPut(" R0: " << sqrt((xp-xm)*(xp-xm) + (yp-ym)*(yp-ym)+ (zp-zm)*(zp-zm)) << endl );
              XDistance[N_AuxPoints] = xp - xm;
              YDistance[N_AuxPoints] = yp - ym;
              ZDistance[N_AuxPoints] = zp - zm;
              IntAux[N_AuxPoints] = JointDOF[LocIsoDOF[j]];
              N_AuxPoints++;
             }
            } //    for(j=0;j<=N_A
           }
          else if(ApproximationOrder>2)
           {
            cout<< " ApproximationOrder " <<ApproximationOrder<< " " << N_AuxVertices << " " <<TDatabase::ParamDB->MOVING_BOUNDARY<< endl;
            OutPut("Tetra isoparametric for approximation order other than 2 has to be implemented !!! "<<endl;);
            exit(0);
           }
          } // if(TDatabase::ParamDB->MOVING_BOUNDARY)
        else if(TDatabase::ParamDB->USE_PRM!=0)
         {
          //     /***************************************************************************/
          //     // this is for isoparametric tetrahedral finite elements in the 3d circular
        if(N_AuxPoints == 0)
         {
          for(i=0;i<4;i++) // for all vertices  
           {
            // get coordinates of vertex i
            cell->GetVertex(i)->GetCoords(xp, yp, zp);
           // point on cylinder
           if(fabs((xp-0.5)*(xp-0.5)+(yp-0.2)*(yp-0.2)-0.0025)<1e-10)
           //if(fabs(xp*xp+yp*yp+zp*zp-1)<1e-10)
            {
          // point on boundary
          // look to all other vertices
          for(j=i+1;j<4;j++)
          {
            // get coordinates of second vertex
            cell->GetVertex(j)->GetCoords(xm, ym, zm);
            // if second vertex on cylinder
            if(fabs((xm-0.5)*(xm-0.5)+(ym-0.2)*(ym-0.2)-0.0025)<1e-10)
              //if(fabs(xm*xm+ym*ym+zm*zm-1)<1e-10)
            {
              // point on boundary
              // find number of dof for P_2
              bf = TFEDatabase3D::GetBaseFunct3D(
                BaseFunctFromOrder[ApproximationOrder]);
              if (ApproximationOrder==2)
              {
                // switch number of first vertex
                switch(i)
                {
                  case 0:
                    // switch number of second vertex
                    switch(j)
                    {
                      case 1: m = 1; break;
                      case 2: m = 3; break;
                      case 3: m = 6; break;
                    } // endswitch(j)
                    break;
                    
                  case 1:
                    switch(j)
                    {
                      case 2: m = 4; break;
                      case 3: m = 7; break;
                    } // endswitch(j)
                    break;
                    
                  case 2:
                    m = 8;
                    break;
                } // endswitch(i)
                
                // centre of edge
                a = 0.5*(xp+xm);
                b = 0.5*(yp+ym);
                c = 0.5*(zp+zm);
                // projection to the boundary
                T = sqrt((a-0.5)*(a-0.5)+(b-0.2)*(b-0.2));
                project_a = 0.5 + 0.05*(a-0.5)/T;
                project_b = 0.2 + 0.05*(b-0.2)/T;    
                project_c = c;
                // T = sqrt(a*a+b*b+c*c);
                // project_a = a/T;
                // project_b = b/T;    
                //project_c = c/T;    
                if(fabs(a-project_a) > 1e-8 || fabs(b-project_b) > 1e-8  || fabs(c-project_c) > 1e-8)
                {
                  // OutPut("edge: x " << a << " project x " << project_a << " y " << b 
                  //      << " project y " << project_b << " z " << (zm+zp)/2.0 << " "  
                  //       << sqrt((project_a-0.5)*(project_a-0.5)+(project_b-0.2)*(project_b-0.2)) << endl);
                  //<< sqrt((project_a)*(project_a)+(project_b)*(project_b)+(project_c)*(project_c)) << endl);
                  // differences
                  XDistance[N_AuxPoints] = project_a - a;
                  YDistance[N_AuxPoints] = project_b - b;
                  ZDistance[N_AuxPoints] = project_c - c;
                  // dof which is moved    
                  IntAux[N_AuxPoints] = m;
                  N_AuxPoints++;
                }
              }
              // P_3
              if (ApproximationOrder==3)
              {
                switch(i)
                {
                  case 0:
                    switch(j)
                    {
                      case 1: m = 1; n = 2; break;
                      case 2: m = 4; n = 7; break;
                      case 3: m = 10; n = 16;  break;
                    } // endswitch(j)
                    break;
                    
                  case 1:
                    switch(j)
                    {
                      case 2: m = 6; n = 8; break;
                      case 3: m = 12; n = 17; break;
                    } // endswitch(j)
                    break;
                    
                  case 2:
                    m = 15; n = 18;
                    break;
                } // endswitch(i)
                
                //  first dof
                a = (2*xp+xm)/3;
                b = (2*yp+ym)/3;
                // projection to the boundary
                T = sqrt((a-0.5)*(a-0.5)+(b-0.2)*(b-0.2));
                project_a = 0.5 + 0.05*(a-0.5)/T;
                project_b = 0.2 + 0.05*(b-0.2)/T;         
                // OutPut("0: a " << a << " project a " << project_a << " b " << b 
                //       << " project b " << project_b << " "
                //       << sqrt((project_a-0.5)*(project_a-0.5)+(project_b-0.2)*(project_b-0.2)) << endl);
                // differences
                XDistance[N_AuxPoints] = project_a - a;
                YDistance[N_AuxPoints] = project_b - b;
                ZDistance[N_AuxPoints] = 0;
                IntAux[N_AuxPoints] = m;
                N_AuxPoints++;
                //  second dof
                a = (xp+2*xm)/3;
                b = (yp+2*ym)/3;
                // projection to the boundary
                project_a = 0.5 + 0.05*(a-0.5)/T;
                project_b = 0.2 + 0.05*(b-0.2)/T;         
                //OutPut("1: a " << a << " project a " << project_a << " b " << b 
                //       << " project b " << project_b << " " 
                //      << sqrt((project_a-0.5)*(project_a-0.5)+(project_b-0.2)*(project_b-0.2)) << endl);
                // differences
                XDistance[N_AuxPoints] = project_a - a;
                YDistance[N_AuxPoints] = project_b - b;
                ZDistance[N_AuxPoints] = 0;
                IntAux[N_AuxPoints] = n;
                N_AuxPoints++;
              }
              
            } // endif
          } // endfor j
        } // endif
      } // endfor i
    } // endif                      
            
            
         } // Not MovingBoundary but IsoBoundFace
        else
         {
          cout<< " ApproximationOrder " <<ApproximationOrder<< " " << N_AuxVertices << " " <<TDatabase::ParamDB->MOVING_BOUNDARY<< endl;
          OutPut("Tetra isoparametric : Neither PRM used nor additional vertices added on IsoBoundFace, Use TetraAffine !!! "<<endl;);
          exit(0);          
         }
         
     } // if(type == BoundaryFace || type == IsoBoundFace)
    } // endfor i
//  OutPut(N_AuxPoints << " : N_AuxPoints, " <<  " N_QuadPoints  :" << N_QuadPoints << endl);

//   if (N_AuxPoints)
//   {
//     OutPut("auxa " << N_AuxPoints << endl);
//     OutPut(x0 << " " << y0 << " " << z0 << endl);
//     OutPut(x1 << " " << y1 << " " << z1 << endl);
//     OutPut(x2 << " " << y2 << " " << z2 << endl);
//     OutPut(x3 << " " << y3 << " " << z3 << endl);
//     }

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
      else
       { 
      OutPut(N_AuxPoints << " : N_AuxPoints, " <<  " N_QuadPoints  >" << N_QuadPoints << endl);
      exit(0);
     }
}

/** return outer normal unit vector */
void TTetraIsoparametric::GetOuterNormal(int j, double s, double t,
                                         double &n1, double &n2, double &n3)
{
  double a1, a2, a3, b1, b2, b3;


  GetTangentVectors(j, s, t, a1, a2, a3, b1, b2, b3);

  n1 = a2*b3 - a3*b2;
  n2 = a3*b1 - a1*b3;
  n3 = a1*b2 - a2*b1;

}

/** return two tangent vectors */
void TTetraIsoparametric::GetTangentVectors(int j, double p1, double p2,
    double &t11, double &t12, double &t13,
    double &t21, double &t22, double &t23)
{

  double Xi, Eta, Zeta;
  int i, k, temp;
  double AuxVector[3*MaxN_BaseFunctions3D];
  TBaseFunct3D *bf;

   bf = TFEDatabase3D::GetBaseFunct3D(BaseFunctFromOrder[ApproximationOrder]);

  switch(j)
  {
    case 0:
      Xi = p1; Eta = p2; Zeta = 0;

      t11 = xc2;
      t12 = yc2;
      t13 = zc2;

      t21 = xc1;
      t22 = yc1;
      t23 = zc1;

      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);

      for(k=0;k<N_AuxPoints;k++)
      {
        i = IntAux[k];  // local dof in this face
        t21 += XDistance[k] * AuxVector[i];
        t11 += XDistance[k] * AuxVector[MaxN_BaseFunctions3D   +  i];

        t22 += YDistance[k] * AuxVector[i];
        t12 += YDistance[k] * AuxVector[MaxN_BaseFunctions3D   +  i];

        t23 += ZDistance[k] * AuxVector[i];
        t13 += ZDistance[k] * AuxVector[MaxN_BaseFunctions3D   +  i];
      } // endfor k

    break;

    case 1:
      Xi = p2; Eta = 0; Zeta = p1;

      t11 = xc1;
      t12 = yc1;
      t13 = zc1;

      t21 = xc3;
      t22 = yc3;
      t23 = zc3;

      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);

      for(k=0;k<N_AuxPoints;k++)
       {
        i = IntAux[k];  // local dof in this face
        t11 += XDistance[k] * AuxVector[i];
        t21 += XDistance[k] * AuxVector[MaxN_BaseFunctions3D + i];

        t12 += YDistance[k] * AuxVector[i];
        t22 += YDistance[k] * AuxVector[MaxN_BaseFunctions3D + i];

        t13 += ZDistance[k] * AuxVector[i];
        t23 += ZDistance[k] * AuxVector[MaxN_BaseFunctions3D + i];
      } // endfor k

    break;

    case 2:
      Xi = p1; Eta = 1-p1- p2; Zeta = p2;

      t11 = (xc3-xc2);
      t12 = (yc3-yc2);
      t13 = (zc3-zc2);

      t21 = (xc1-xc2);
      t22 = (yc1-yc2);
      t23 = (zc1-zc2);

      bf->GetDerivatives(D100, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+2*MaxN_BaseFunctions3D);

      for(k=0;k<N_AuxPoints;k++)
       {
        i = IntAux[k];  // local dof in this face
        t21 += XDistance[k] *(  AuxVector[i]   - AuxVector[MaxN_BaseFunctions3D + i] );
        t11 += XDistance[k] *(  AuxVector[2*MaxN_BaseFunctions3D + i]
                              - AuxVector[MaxN_BaseFunctions3D + i] );

        t22 += YDistance[k] *(  AuxVector[i]   - AuxVector[MaxN_BaseFunctions3D + i] );
        t12 += YDistance[k] *(  AuxVector[2*MaxN_BaseFunctions3D + i]
                              - AuxVector[MaxN_BaseFunctions3D + i] );

        t23 += ZDistance[k] *(  AuxVector[i]   - AuxVector[MaxN_BaseFunctions3D + i] );
        t13 += ZDistance[k] *(  AuxVector[2*MaxN_BaseFunctions3D + i]
                              - AuxVector[MaxN_BaseFunctions3D + i] );
      } // endfor k

    break;

    case 3:
      Xi = 0; Eta = p1; Zeta = p2;

      t11 = xc3;
      t12 = yc3;
      t13 = zc3;

      t21 = xc2;
      t22 = yc2;
      t23 = zc2;

      bf->GetDerivatives(D010, Xi, Eta, Zeta, AuxVector);
      bf->GetDerivatives(D001, Xi, Eta, Zeta, AuxVector+MaxN_BaseFunctions3D);

      for(k=0;k<N_AuxPoints;k++)
       {
        i = IntAux[k];  // local dof in this face
        t21 += XDistance[k] * AuxVector[i];
        t11 += XDistance[k] * AuxVector[MaxN_BaseFunctions3D + i];

        t22 += YDistance[k] * AuxVector[i];
        t12 += YDistance[k] * AuxVector[MaxN_BaseFunctions3D + i];

        t23 += ZDistance[k] * AuxVector[i];
        t13 += ZDistance[k] * AuxVector[MaxN_BaseFunctions3D + i];
      } // endfor k

    break;

    default:
      Error("Wrong joint number for tetrahedron! " << endl);
      return;
  }

}
