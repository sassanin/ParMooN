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
// @(#)HexaTrilinear.C      1.6 02/22/00
//
// Class:      THexaTrilinear
//
// Purpose:    trilinear reference transformations for hexahedron
//
// Author:     Daniel Quoos  
//
// History:    11.07.00 start implementation
// 
// =======================================================================

#include <HexaTrilinear.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>
#include <MooNMD_Io.h>

/** constuctor */
THexaTrilinear::THexaTrilinear()
{
}

/** transfer from reference element to original element */
void THexaTrilinear::GetOrigFromRef(double xi, double eta, double zeta, double &X, double &Y, double &Z)
{
  X = xc0 + xc1*xi + xc2*eta + xc3*zeta + xc4*xi*eta
          + xc5*xi*zeta + xc6*eta*zeta + xc7*xi*eta*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta + yc4*xi*eta
          + yc5*xi*zeta + yc6*eta*zeta + yc7*xi*eta*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta + zc4*xi*eta
          + zc5*xi*zeta + zc6*eta*zeta + zc7*xi*eta*zeta;
}

/** transfer a set of point from reference to original element */
void THexaTrilinear::GetOrigFromRef(
                int N_Points,
                double *xi, double *eta, double *zeta,
                double *X, double *Y, double *Z, double *absdetjk)
{
  int i;
  double Xi, Eta, Zeta;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  
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
    
    detjk=dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2 - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta + xc4*Xi*Eta
               + xc5*Xi*Zeta + xc6*Eta*Zeta + xc7*Xi*Eta*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta + yc4*Xi*Eta
               + yc5*Xi*Zeta + yc6*Eta*Zeta + yc7*Xi*Eta*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta + zc4*Xi*Eta
               + zc5*Xi*Zeta + zc6*Eta*Zeta + zc7*Xi*Eta*Zeta;
    absdetjk[i] = fabs(detjk);
  }
}

/** transfer from reference element to original element */
void THexaTrilinear::GetOrigFromRef(double *ref, double *orig)
{
  orig[0] = xc0 + xc1*ref[0] + xc2*ref[1] + xc3*ref[2] + xc4*ref[0]*ref[1]
                + xc5*ref[0]*ref[2] + xc6*ref[1]*ref[2]
                + xc7*ref[0]*ref[1]*ref[2];
  orig[1] = yc0 + yc1*ref[0] + yc2*ref[1] + yc3*ref[2] + yc4*ref[0]*ref[1]
                + yc5*ref[0]*ref[2] + yc6*ref[1]*ref[2]
                + yc7*ref[0]*ref[1]*ref[2];
  orig[2] = zc0 + zc1*ref[0] + zc2*ref[1] + zc3*ref[2] + zc4*ref[0]*ref[1]
                + zc5*ref[0]*ref[2] + zc6*ref[1]*ref[2]
                + zc7*ref[0]*ref[1]*ref[2];
}


/** transfer from original element to reference element */
void THexaTrilinear::GetRefFromOrig(double X, double Y, double Z,
                                    double &xi, double &eta, double &zeta)
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
void THexaTrilinear::GetRefFromOrig(double *orig, double *ref)
{
  GetRefFromOrig(orig[0], orig[1], orig[2], ref[0], ref[1], ref[2]);
}

/** calculate functions and derivatives from reference element
    to original element */
void THexaTrilinear::GetOrigValues(BaseFunct3D BaseFunct,
                               int N_Points, double *xi, double *eta, double *zeta,
                               int N_Functs, QuadFormula3D QuadFormula)
{
  int i,j,k;
  double **refvaluesD000, **origvaluesD000;
  double **refvaluesD100, **origvaluesD100;
  double **refvaluesD010, **origvaluesD010;
  double **refvaluesD001, **origvaluesD001;
  double **refvaluesD200, **origvaluesD200;
  double **refvaluesD110, **origvaluesD110;
  double **refvaluesD020, **origvaluesD020;
  double **refvaluesD101, **origvaluesD101;
  double **refvaluesD011, **origvaluesD011;
  double **refvaluesD002, **origvaluesD002;
  double *refD000, *origD000;
  double *refD100, *origD100;
  double *refD010, *origD010;
  double *refD001, *origD001;
  double *refD200, *origD200;
  double *refD110, *origD110;
  double *refD020, *origD020;
  double *refD101, *origD101;
  double *refD011, *origD011;
  double *refD002, *origD002;
  double *aux;
  double AllData[MaxN_BaseFunctions3D][5];
  double GeoData[5][5];
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  
  int BaseVectDim = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
  if(BaseVectDim != 1)
  {
    ErrMsg("HexaIsoparametric for vector valued functions not " <<
           "implemented. (Piola transform)");
    exit(0);
  }
 
  refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D000);
  if(refvaluesD000==NULL)
    {
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, 
                                                       QuadFormula, D000);
    }
  
  origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
  if(origvaluesD000==NULL)
    {
      origvaluesD000 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD000[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D000, origvaluesD000);
    }
  
  refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
  origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
  if(origvaluesD100==NULL)
    {
      origvaluesD100 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD100[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D100, origvaluesD100);
    }
  
  refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
  origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
  if(origvaluesD010==NULL)
    {
      origvaluesD010 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD010[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D010, origvaluesD010);
    } 
  
  refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
  origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
  if(origvaluesD001==NULL)
    { 
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD001[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D001, origvaluesD001);
    } 
  
  refvaluesD200=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D200);
  origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
  if(origvaluesD200==NULL)
    { 
      origvaluesD200 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD200[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D200, origvaluesD200);
    } 
  
  refvaluesD110=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D110);
  origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
  if(origvaluesD110==NULL)
    { 
      origvaluesD110 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD110[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D110, origvaluesD110);
    } 
  
  refvaluesD101=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D101);
  origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
  if(origvaluesD101==NULL)
    { 
      origvaluesD101 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD101[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D101, origvaluesD101);
    } 
  
  refvaluesD011=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D011);
  origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
  if(origvaluesD011==NULL)
    { 
      origvaluesD011 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD011[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D011, origvaluesD011);
    } 
  
  refvaluesD020=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D020);
  origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
  if(origvaluesD020==NULL)
    {
      origvaluesD020 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD020[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D020, origvaluesD020);
    } 
  
  refvaluesD002=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D002);
  origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
  if(origvaluesD002==NULL)
    {
      origvaluesD002 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD002[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D002, origvaluesD002);
    } 
  // D000
  for(i=0;i<N_Points;i++)
    {
      refD000 = refvaluesD000[i];
      origD000 = origvaluesD000[i];
      
      for(j=0;j<N_Functs;j++)
        {
          origD000[j] = refD000[j];
        } // endfor j
    } // endfor i
  
  // /*
  // D100, D010 and D001
  for(i=0;i<N_Points;i++)
    { 
      refD100 = refvaluesD100[i];
      origD100 = origvaluesD100[i];
      
      refD010 = refvaluesD010[i];
      origD010 = origvaluesD010[i];
      
      refD001 = refvaluesD001[i];
      origD001 = origvaluesD001[i];
      
      dx1=xc1 + xc4*eta[i] + xc5*zeta[i] + xc7*eta[i]*zeta[i];
      dx2=xc2 + xc4*xi[i] + xc6*zeta[i] + xc7*xi[i]*zeta[i];
      dx3=xc3 + xc5*xi[i] + xc6*eta[i] + xc7*xi[i]*eta[i];

      dy1=yc1 + yc4*eta[i] + yc5*zeta[i] + yc7*eta[i]*zeta[i];
      dy2=yc2 + yc4*xi[i] + yc6*zeta[i] + yc7*xi[i]*zeta[i];
      dy3=yc3 + yc5*xi[i] + yc6*eta[i] + yc7*xi[i]*eta[i];

      dz1=zc1 + zc4*eta[i] + zc5*zeta[i] + zc7*eta[i]*zeta[i];
      dz2=zc2 + zc4*xi[i] + zc6*zeta[i] + zc7*xi[i]*zeta[i];
      dz3=zc3 + zc5*xi[i] + zc6*eta[i] + zc7*xi[i]*eta[i];
      
      detjk=dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2 - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
      rec_detjk=1/detjk;

      for(j=0;j<N_Functs;j++)
        {
          origD100[j] = ((dy2*dz3 - dy3*dz2)*refD100[j] + (dy3*dz1 - dy1*dz3)*refD010[j] 
                         + (dy1*dz2 - dy2*dz1)*refD001[j])*rec_detjk;
          origD010[j] = ((dx3*dz2 - dx2*dz3)*refD100[j] + (dx1*dz3 - dx3*dz1)*refD010[j] 
                         + (dx2*dz1 - dx1*dz2)*refD001[j])*rec_detjk;
          origD001[j] = ((dx2*dy3 - dy2*dx3)*refD100[j] + (dx3*dy1 - dx1*dy3)*refD010[j] 
                         + (dx1*dy2 - dx2*dy1)*refD001[j])*rec_detjk;
          
          
          // cout << "10: " << origD10[j] << endl;
          // cout << "01: " << origD01[j] << endl;
          // cout << endl;
        } // endfor j
      // cout << "----------" << endl;
    } // endfor i
  
  /*
    for(i=0;i<N_Points;i++)
    {
    // reset matrix
    memset(GeoData, 0, 25*SizeOfDouble);
    
    refD10 = refvaluesD10[i];
    refD01 = refvaluesD01[i];
    refD20 = refvaluesD20[i];
    refD11 = refvaluesD11[i];
    refD02 = refvaluesD02[i];
    
    origD10 = origvaluesD10[i];
    origD01 = origvaluesD01[i];
    origD20 = origvaluesD20[i];
    origD11 = origvaluesD11[i];
    origD02 = origvaluesD02[i];
    
    GeoData[0][0] = xc1;
    GeoData[0][1] = yc1;

    GeoData[1][0] = xc2;
    GeoData[1][1] = yc2;

    GeoData[2][2] = xc1*xc1;
    GeoData[2][3] = 2*xc1*yc1;
    GeoData[2][4] = yc1*yc1;
    
    GeoData[3][2] = xc1*xc2;
    GeoData[3][3] = yc1*xc2+xc1*yc2; 
    GeoData[3][4] = yc1*yc2;
    
    GeoData[4][2] = xc2*xc2;
    GeoData[4][3] = 2*xc2*yc2;
    GeoData[4][4] = yc2*yc2;

    for(j=0;j<N_Functs;j++)
    {
    AllData[j][0] = refD10[j];
    AllData[j][1] = refD01[j];
    AllData[j][2] = refD20[j];
    AllData[j][3] = refD11[j];
    AllData[j][4] = refD02[j];
    } // endfor j
    
    // subroutine for solving a multiple systems of linear equations
    // void SolveMultipleSystems(double *a, double *b, int N_Eqn,
    //                        int LDA, int LDB, int N_Rhs);
    SolveMultipleSystems((double *)GeoData, (double *)AllData, 5,
    5, 5, N_Functs);
    
    for(j=0;j<N_Functs;j++)
    {
    origD10[j] = AllData[j][0];
    origD01[j] = AllData[j][1];
    origD20[j] = AllData[j][2];
    origD11[j] = AllData[j][3];
    origD02[j] = AllData[j][4];
    
    } // endfor j
    
    } // endfor i
     */
      }

/** calculate functions and derivatives from reference element
    to original element, for all given elements */
void THexaTrilinear::GetOrigValues(int N_Sets, BaseFunct3D *BaseFuncts,
                               int N_Points, double *xi, double *eta, double *zeta, 
                               QuadFormula3D QuadFormula,
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
  
  SecondDer = FALSE;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
      
    refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D000);
    if(refvaluesD000==NULL)
    {
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, 
                                                           QuadFormula, D000);
    }
      
    origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    if(origvaluesD000==NULL)
    {
      origvaluesD000 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD000[j] = aux+j*N_Functs;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D000, origvaluesD000);
    }
      
    for(j=0;j<N_Points;j++)
    {
      refD000 = refvaluesD000[j];
      origD000 = origvaluesD000[j];
          
      memcpy(origD000, refD000, N_Functs*SizeOfDouble);
    } // endfor j
      
    refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
    origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    if(origvaluesD100==NULL)
    {
      origvaluesD100 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD100[j] = aux+j*N_Functs;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D100, origvaluesD100);
    }
      
    refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
    origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    if(origvaluesD010==NULL)
    {
      origvaluesD010 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD010[j] = aux+j*N_Functs;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D010, origvaluesD010);
    }
      
    refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
    origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    if(origvaluesD001==NULL)
    {
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD001[j] = aux+j*N_Functs;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D001, origvaluesD001);
    }
  }
      /*
      if(Needs2ndDer[i])
        {
          SecondDer = TRUE;
          
          refvaluesD200=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D200);
          origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
          if(origvaluesD200==NULL)
            {
              origvaluesD200 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*N_Functs];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD200[j] = aux+j*N_Functs;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D200, origvaluesD200);
            }
          
          refvaluesD110=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D110);
          origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
          if(origvaluesD110==NULL)
            {
              origvaluesD110 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*N_Functs];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD110[j] = aux+j*N_Functs;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D110, origvaluesD110);
            }
          
          refvaluesD101=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D101);
          origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
          if(origvaluesD101==NULL)
            {
              origvaluesD101 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*N_Functs];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD101[j] = aux+j*N_Functs;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D101, origvaluesD101);
            }
          
          refvaluesD011=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D011);
          origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
          if(origvaluesD011==NULL)
            {
              origvaluesD011 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*N_Functs];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD011[j] = aux+j*N_Functs;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D011, origvaluesD011);
            }
          
          refvaluesD020=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D020);
          origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
          if(origvaluesD020==NULL)
            {
              origvaluesD020 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*N_Functs];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD020[j] = aux+j*N_Functs;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D020, origvaluesD020);
            }
          
          refvaluesD002=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D002);
          origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
          if(origvaluesD002==NULL)
            {
              origvaluesD002 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*N_Functs];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD002[j] = aux+j*N_Functs;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D002, origvaluesD002);
            } // endfor Needs2ndDer[i]
        } // endfor i
      */
      // D100, D010 and D001
      for(i=0;i<N_Sets;i++)
        {
          BaseFunct=BaseFuncts[i];
          N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
          
          refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
          origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
          
          refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
          origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
          
          refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
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

              detjk=dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2 - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
              rec_detjk=1/detjk;
              
              for(k=0;k<N_Functs;k++)
                {
                  origD100[k] = ((dy2*dz3 - dy3*dz2)*refD100[k] + (dy3*dz1 - dy1*dz3)*refD010[k]
                                 + (dy1*dz2 - dy2*dz1)*refD001[k])*rec_detjk;
                  origD010[k] = ((dx3*dz2 - dx2*dz3)*refD100[k] + (dx1*dz3 - dx3*dz1)*refD010[k] 
                                 + (dx2*dz1 - dx1*dz2)*refD001[k])*rec_detjk;
                  origD001[k] = ((dx2*dy3 - dy2*dx3)*refD100[k] + (dx3*dy1 - dx1*dy3)*refD010[k] 
                                 + (dx1*dy2 - dx2*dy1)*refD001[k])*rec_detjk;                            
                } // endfor k
            } // endfor j
        } // endfor i
   
      // leave if no second derivatives are needed
      if(!SecondDer) return;
      /*
        // find transformation matrix for second derivatives
        // do this only once since matrix is constant for a trilinear mapping
        
        // reset matrices
        memset(GeoData, 0, 9*SizeOfDouble);
        memset(Eye, 0, 9*SizeOfDouble);
        Eye[0][0] = 1;
        Eye[1][1] = 1;
        Eye[2][2] = 1;
        
        GeoData[0][0] = xc1*xc1;
        GeoData[0][1] = 2*xc1*yc1;
        GeoData[0][2] = yc1*yc1;
        
        GeoData[1][0] = xc1*xc2;
        GeoData[1][1] = yc1*xc2+xc1*yc2; 
        GeoData[1][2] = yc1*yc2;
        
        GeoData[2][0] = xc2*xc2;
        GeoData[2][1] = 2*xc2*yc2;
        GeoData[2][2] = yc2*yc2;
        
        // subroutine for solving a multiple systems of linear equations
        // void SolveMultipleSystems(double *a, double *b, int N_Eqn,
        //                        int LDA, int LDB, int N_Rhs);
        //SolveMultipleSystems((double *)GeoData, (double *)Eye, 3,
        //                 3, 3, 3);
        //SolveLinearSystem((double *)GeoData, (double *)Eye, 3, 3);
                    
        SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 3,
        3, 3, 3);
        
        for(i=0;i<N_Sets;i++)
        {
  if(Needs2ndDer[i])
  {
  BaseFunct=BaseFuncts[i];
  N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
  
  refvaluesD20=TFEDatabase3D::GetRefElementValues
  (BaseFunct, QuadFormula, D20);
  origvaluesD20=TFEDatabase3D::GetOrigElementValues
  (BaseFunct, D20);
  
  refvaluesD11=TFEDatabase3D::GetRefElementValues
  (BaseFunct, QuadFormula, D11);
  origvaluesD11=TFEDatabase3D::GetOrigElementValues
  (BaseFunct, D11);
  
  refvaluesD02=TFEDatabase3D::GetRefElementValues
  (BaseFunct, QuadFormula, D02);
  origvaluesD02=TFEDatabase3D::GetOrigElementValues
  (BaseFunct, D02);
  
  for(j=0;j<N_Points;j++)
  {
  refD20 = refvaluesD20[j];
  refD11 = refvaluesD11[j];
  refD02 = refvaluesD02[j];
  
  origD20 = origvaluesD20[j];
  origD11 = origvaluesD11[j];
  origD02 = origvaluesD02[j];
  
  for(k=0;k<N_Functs;k++)
  {
  r20 = refD20[k];
  r11 = refD11[k];
  r02 = refD02[k];
  
  o20 = Eye[0][0]*r20+Eye[0][1]*r11+Eye[0][2]*r02;
  o11 = Eye[1][0]*r20+Eye[1][1]*r11+Eye[1][2]*r02;
  o02 = Eye[2][0]*r20+Eye[2][1]*r11+Eye[2][2]*r02;
  
  origD20[k] = o20;
  origD11[k] = o11;
  origD02[k] = o02;
  } // endfor k
  } // endif
  } // endfor j
  } // endfor i
      */
}
/** calculate functions and derivatives from reference element
    to original element */
void THexaTrilinear::GetOrigValues(double xi, double eta, double zeta,
                               int N_BaseFunct,
                               double *uref, double *uxiref, double *uetaref, double *uzetaref,
                               double *uorig, double *uxorig, double *uyorig, double *uzorig, 
                               int _BaseVectDim)
{
  int i;
  double dx1, dx2, dx3;
  double dy1, dy2, dy3;
  double dz1, dz2, dz3;
  
  if(_BaseVectDim != 1)
  {
    ErrMsg("mixed finite elements with trilinear reference transformation not available");
    exit(0);
  }
  
  // D000
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  dx1=xc1 + xc4*eta + xc5*zeta + xc7*eta*zeta;
  dx2=xc2 + xc4*xi + xc6*zeta + xc7*xi*zeta;
  dx3=xc3 + xc5*xi + xc6*eta + xc7*xi*eta;
              
  dy1=yc1 + yc4*eta + yc5*zeta + yc7*eta*zeta;
  dy2=yc2 + yc4*xi + yc6*zeta + yc7*xi*zeta;
  dy3=yc3 + yc5*xi + yc6*eta + yc7*xi*eta;
  
  dz1=zc1 + zc4*eta + zc5*zeta + zc7*eta*zeta;
  dz2=zc2 + zc4*xi + zc6*zeta + zc7*xi*zeta;
  dz3=zc3 + zc5*xi + zc6*eta + zc7*xi*eta;

  detjk=dx1*dy2*dz3 + dx2*dy3*dz1 + dx3*dy1*dz2 - dx3*dy2*dz1 - dx2*dy1*dz3 - dx1*dy3*dz2;
  rec_detjk=1/detjk;

  // D100, D010 and D001
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((dy2*dz3 - dy3*dz2)*uxiref[i] + (dy3*dz1 - dy1*dz3)*uetaref[i] + (dy1*dz2 - dy2*dz1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((dx3*dz2 - dx2*dz3)*uxiref[i] + (dx1*dz3 - dx3*dz1)*uetaref[i] + (dx2*dz1 - dx1*dz2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((dx2*dy3 - dx3*dy2)*uxiref[i] + (dx3*dy1 - dx1*dy3)*uetaref[i] + (dx1*dy2 - dx2*dy1)*uzetaref[i]) *rec_detjk;
  } // endfor i
}

void THexaTrilinear::SetCell(TBaseCell *cell)
{
  int i;

  Cell = cell;

  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);
  Cell->GetVertex(4)->GetCoords(x4, y4, z4);
  Cell->GetVertex(5)->GetCoords(x5, y5, z5);
  Cell->GetVertex(6)->GetCoords(x6, y6, z6);
  Cell->GetVertex(7)->GetCoords(x7, y7, z7);

/*
  OutPut("0: " << Cell->GetVertex(0) << endl);
  OutPut("1: " << Cell->GetVertex(1) << endl);
  OutPut("2: " << Cell->GetVertex(2) << endl);
  OutPut("3: " << Cell->GetVertex(3) << endl);
  OutPut("4: " << Cell->GetVertex(4) << endl);
  OutPut("5: " << Cell->GetVertex(5) << endl);
  OutPut("6: " << Cell->GetVertex(6) << endl);
  OutPut("7: " << Cell->GetVertex(7) << endl);
*/

  xc0 = (x0 + x1 + x3 + x2 + x4 + x5 + x7 + x6) * 0.125;
  xc1 = (-x0 + x1 - x3 + x2 - x4 + x5 - x7 + x6) * 0.125;
  xc2 = (-x0 - x1 + x3 + x2 - x4 - x5 + x7 + x6) * 0.125;
  xc3 = (-x0 - x1 - x3 - x2 + x4 + x5 + x7 + x6) * 0.125;
  xc4 = (x0 - x1 - x3 + x2 + x4 - x5 - x7 + x6) * 0.125;
  xc5 = (x0 - x1 + x3 - x2 - x4 + x5 - x7 + x6) * 0.125;
  xc6 = (x0 + x1 - x3 - x2 - x4 - x5 + x7 + x6) * 0.125;
  xc7 = (-x0 + x1 + x3 - x2 + x4 - x5 - x7 + x6) * 0.125;

  // cout << "x:" << xc4 << xc5 << xc6 << xc7 << endl;

  yc0 = (y0 + y1 + y3 + y2 + y4 + y5 + y7 + y6) * 0.125;
  yc1 = (-y0 + y1 - y3 + y2 - y4 + y5 - y7 + y6) * 0.125;
  yc2 = (-y0 - y1 + y3 + y2 - y4 - y5 + y7 + y6) * 0.125;
  yc3 = (-y0 - y1 - y3 - y2 + y4 + y5 + y7 + y6) * 0.125;
  yc4 = (y0 - y1 - y3 + y2 + y4 - y5 - y7 + y6) * 0.125;
  yc5 = (y0 - y1 + y3 - y2 - y4 + y5 - y7 + y6) * 0.125;
  yc6 = (y0 + y1 - y3 - y2 - y4 - y5 + y7 + y6) * 0.125;
  yc7 = (-y0 + y1 + y3 - y2 + y4 - y5 - y7 + y6) * 0.125;

  // cout << "y:" << yc4 << yc5 << yc6 << yc7 << endl;

  zc0 = (z0 + z1 + z3 + z2 + z4 + z5 + z7 + z6) * 0.125;
  zc1 = (-z0 + z1 - z3 + z2 - z4 + z5 - z7 + z6) * 0.125;
  zc2 = (-z0 - z1 + z3 + z2 - z4 - z5 + z7 + z6) * 0.125;
  zc3 = (-z0 - z1 - z3 - z2 + z4 + z5 + z7 + z6) * 0.125;
  zc4 = (z0 - z1 - z3 + z2 + z4 - z5 - z7 + z6) * 0.125;
  zc5 = (z0 - z1 + z3 - z2 - z4 + z5 - z7 + z6) * 0.125;
  zc6 = (z0 + z1 - z3 - z2 - z4 - z5 + z7 + z6) * 0.125;
  zc7 = (-z0 + z1 + z3 - z2 + z4 - z5 - z7 + z6) * 0.125;

  // cout << "z:" << zc4 << zc5 << zc6 << zc7 << endl;
}

/** return outer normal unit vector */
void THexaTrilinear::GetOuterNormal(int j, double s, double t,
                                    double &n1, double &n2, double &n3)
{
//   double len;
  double Xi, Eta, Zeta;
  double t11, t12, t13;
  double t21, t22, t23;

  switch(j)
  {
    case 0:
      Xi = s; Eta = t; Zeta = -1;
      t11 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t12 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t13 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;

      t21 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t22 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t23 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    break;

    case 1:
      Xi = t; Eta = -1; Zeta = s;
      t11 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t12 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t13 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;

      t21 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t22 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t23 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    break;

    case 2:
      Xi = 1; Eta = t; Zeta = s;
      t11 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t12 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t13 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;

      t21 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t22 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t23 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    break;

    case 3:
      Xi = -t; Eta = 1; Zeta = s;
      t11 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t12 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t13 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;

      t21 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t22 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t23 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    break;

    case 4:
      Xi = -1; Eta = s; Zeta = t;
      t11 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t12 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t13 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
      
      t21 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t22 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t23 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    break;

    case 5:
      Xi = t; Eta = s; Zeta = 1;
      t11 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t12 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t13 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;

      t21 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t22 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t23 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    break;

    default:
      Error("Wrong local joint number" << endl);
      n1 = -1; n2 = -1; n3 = -1;
      return;
  }

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
void THexaTrilinear::GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23)
{
  double Xi, Eta, Zeta;

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
    break;

    case 1:
      Xi = p2; Eta = -1; Zeta = p1;
      t11 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t12 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t13 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;

      t21 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t22 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t23 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    break;

    case 2:
      Xi = 1; Eta = p2; Zeta = p1;
      t11 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t12 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t13 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;

      t21 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t22 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t23 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
    break;

    case 3:
      Xi = -p2; Eta = 1; Zeta = p1;
      t11 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t12 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t13 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;

      t21 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t22 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t23 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;
    break;

    case 4:
      Xi = -1; Eta = p1; Zeta = p2;
      t11 = xc3 + xc5*Xi + xc6*Eta + xc7*Xi*Eta;
      t12 = yc3 + yc5*Xi + yc6*Eta + yc7*Xi*Eta;
      t13 = zc3 + zc5*Xi + zc6*Eta + zc7*Xi*Eta;
      
      t21 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t22 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t23 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    break;

    case 5:
      Xi = p2; Eta = p1; Zeta = 1;
      t11 = xc1 + xc4*Eta + xc5*Zeta + xc7*Eta*Zeta;
      t12 = yc1 + yc4*Eta + yc5*Zeta + yc7*Eta*Zeta;
      t13 = zc1 + zc4*Eta + zc5*Zeta + zc7*Eta*Zeta;

      t21 = xc2 + xc4*Xi + xc6*Zeta + xc7*Xi*Zeta;
      t22 = yc2 + yc4*Xi + yc6*Zeta + yc7*Xi*Zeta;
      t23 = zc2 + zc4*Xi + zc6*Zeta + zc7*Xi*Zeta;
    break;

    default:
      Error("Wrong local joint number" << endl);
      return;
  }
} // end THexaTrilinear::GetTangentVectors


void THexaTrilinear::PiolaMapOrigFromRef(int N_Functs, double *refD000, 
                                         double *origD000)
{
  ErrMsg("Piola Map for HexaTrilinear reference map not yet implemented");
  exit(1);
}
