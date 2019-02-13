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
// @(#)TetraAffin.C        1.6 12/07/99
//
// Class:      TTetraAffin
//
// Purpose:    reference transformations for tetrahedron
//
// Author:     Daniel Quoos
//
// History:    07.02.00 start implementation
// 
// =======================================================================

#include <MooNMD_Io.h>
#include <FEDatabase3D.h>
#include <TetraAffin.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>

/** constuctor */
TTetraAffin::TTetraAffin()
{
}

/** transfer a point from reference face to original element face */
void TTetraAffin::GetOrigBoundFromRef(int Joint, double xi, double eta,
					      double &X, double &Y, double &Z)
{
  double Xi, Eta, Zeta;
  
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
}


/** transfer from reference element to original element */
void TTetraAffin::GetOrigFromRef(double xi, double eta, double zeta, double &X, double &Y, double &Z)
{
  X = xc0 + xc1*xi + xc2*eta + xc3*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta;
}

/** transfer a set of points from reference to original element */
void TTetraAffin::GetOrigFromRef(int N_Points, double *xi, double *eta, double *zeta,
                                double *X, double *Y, double *Z, double *absdetjk)
{
  int i;
  double Xi, Eta, Zeta;
  double absdet = fabs(detjk);

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    Zeta = zeta[i];
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta;
    absdetjk[i] = absdet;
  } // endfor i
}

/** transfer from reference element to original element */
void TTetraAffin::GetOrigFromRef(double *ref, double *orig)
{
  orig[0]=xc0 + xc1*ref[0] + xc2*ref[1] + xc3*ref[2];
  orig[1]=yc0 + yc1*ref[0] + yc2*ref[1] + yc3*ref[2];
  orig[2]=zc0 + zc1*ref[0] + zc2*ref[1] + zc3*ref[2];
}

/** transfer from original element to reference element */
void TTetraAffin::GetRefFromOrig(double X, double Y, double Z, double &xi, double &eta, double &zeta)
{
  double xt=(X - xc0)/detjk;
  double yt=(Y - yc0)/detjk;
  double zt=(Z - zc0)/detjk;

  xi  = (yc2*zc3 - yc3*zc2)*xt - (xc2*zc3 - xc3*zc2)*yt + (xc2*yc3 - xc3*yc2)*zt;
  eta = -(yc1*zc3 - yc3*zc1)*xt + (xc1*zc3 - xc3*zc1)*yt - (xc1*yc3 - xc3*yc1)*zt;
  zeta = (yc1*zc2 - yc2*zc1)*xt - (xc1*zc2 - xc2*zc1)*yt + (xc1*yc2 - xc2*yc1)*zt;


  /*  cout << " a11: " << (yc2*zc3 - yc3*zc2);
  cout << " a12: " << (xc2*zc3 - xc3*zc2);
  cout << " a13: " << (xc2*yc3 - xc3*yc2) << endl;
  cout << " a21: " << (yc1*zc3 - yc3*zc1);
  cout << " a22: " << (xc1*zc3 - xc3*zc1);
  cout << " a23: " << (xc1*yc3 - xc3*yc1) << endl;
  cout << " a31: " << (yc1*zc2 - yc2*zc1);
  cout << " a32: " << (xc1*zc2 - xc2*zc1);
  cout << " a33: " << (xc1*yc2 - xc2*yc1) << endl;
  cout << " xt,yt,zt: " << xt << " " << yt << " " << zt << endl; */
}

/** transfer from original element to reference element */
void TTetraAffin::GetRefFromOrig(double *orig, double *ref)
{
  double xt=(orig[0] - xc0)/detjk;
  double yt=(orig[1] - yc0)/detjk;
  double zt=(orig[2] - zc0)/detjk;

  ref[0]  = (yc2*zc3 - yc3*zc2)*xt - (xc2*zc3 - xc3*zc2)*yt + (xc2*yc3 - xc3*yc2)*zt;
  ref[1] = -(yc1*zc3 - yc3*zc1)*xt + (xc1*zc3 - xc3*zc1)*yt - (xc1*yc3 - xc3*yc1)*zt;
  ref[2] = (yc1*zc2 - yc2*zc1)*xt - (xc1*zc2 - xc2*zc1)*yt + (xc1*yc2 - xc2*yc1)*zt;
  //  ref[0]  = (yc2*zc3 - yc3*zc2)*xt + (yc3*zc1 - yc1*zc3)*yt + (yc1*zc2 - yc2*zc1)*zt;
  // ref[1] = (xc3*zc2 - xc2*zc3)*xt + (xc1*zc3 - xc3*zc1)*yt + (xc2*zc1 - xc1*zc2)*zt;
  // ref[2] = (xc2*yc3 - xc3*yc2)*xt + (xc3*yc1 - xc1*yc3)*yt + (xc1*yc2 - xc2*yc1)*zt;
}

/** calculate functions and derivatives from reference element
    to original element */
void TTetraAffin::GetOrigValues(BaseFunct3D BaseFunct,
                                int N_Points, double *xi, double *eta, double *zeta,
                                int N_Functs, QuadFormula3D QuadFormula)
{
  int i,j,k;
  double **refvaluesD000, **origvaluesD000;
  double **refvaluesD100, **origvaluesD100;
  double **refvaluesD010, **origvaluesD010;
  double **refvaluesD001, **origvaluesD001;
  double **refvaluesD200, **origvaluesD200;
  double **refvaluesD020, **origvaluesD020;
  double **refvaluesD002, **origvaluesD002;
  double **refvaluesD110, **origvaluesD110;
  double **refvaluesD101, **origvaluesD101;
  double **refvaluesD011, **origvaluesD011;
  double *refD000, *origD000;
  double *refD100, *origD100;
  double *refD010, *origD010;
  double *refD001, *origD001;
  double *refD200, *origD200;
  double *refD020, *origD020;
  double *refD002, *origD002;
  double *refD110, *origD110;
  double *refD101, *origD101;
  double *refD011, *origD011;
  double *aux;
  double AllData[MaxN_BaseFunctions3D][5];
  double GeoData[5][5];
  
  int BaseVectDim = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
  if(BaseVectDim != 1)
  {
    ErrMsg("TTetraAffin::GetOrigValues for mixed elements not implemented\n");
    exit(1);
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
  
  // D100, D010 and D001
  for(i=0;i<N_Points;i++)
    {
      refD100 = refvaluesD100[i];
      origD100 = origvaluesD100[i];
      
      refD010 = refvaluesD010[i];
      origD010 = origvaluesD010[i];
      
      refD001 = refvaluesD001[i];
      origD001 = origvaluesD001[i];
      
      for(j=0;j<N_Functs;j++)
        {
          origD100[j] = ((yc2*zc3 - yc3*zc2)*refD100[j] + (yc3*zc1 - yc1*zc3)*refD010[j] + (yc1*zc2 - yc2*zc1)*refD001[j]) *rec_detjk;
          origD010[j] = ((xc3*zc2 - xc2*zc3)*refD100[j] + (xc1*zc3 - xc3*zc1)*refD010[j] + (xc2*zc1 - xc1*zc2)*refD001[j]) *rec_detjk;
          origD001[j] = ((xc2*yc3 - xc3*yc2)*refD100[j] + (xc3*yc1 - xc1*yc3)*refD010[j] + (xc1*yc2 - xc2*yc1)*refD001[j]) *rec_detjk;
          // cout << "10: " << origD10[j] << endl;
          // cout << "01: " << origD01[j] << endl;
          // cout << endl;
        } // endfor j
      // cout << "----------" << endl;
    } // endfor i
}


/** calculate functions and derivatives from reference element
    to original element, for all given elements */
void TTetraAffin::GetOrigValues(int N_Sets, BaseFunct3D *BaseFuncts,
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
  BaseFunct3D BaseFunct;
  int N_Functs;
  bool SecondDer;

  SecondDer = false;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
    int BaseVectDim = 
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
    
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
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD000[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D000, origvaluesD000);
      }
    
    for(j=0;j<N_Points;j++)
      {
        refD000 = refvaluesD000[j];
        origD000 = origvaluesD000[j];
        
       if(BaseVectDim == 1)
    {
      // simply copy values
      memcpy(origD000, refD000, N_Functs*SizeOfDouble);
    }
    else
    {
      // do Piola transformation
      PiolaMapOrigFromRef(N_Functs, refD000, origD000);
    }
      } // endfor j
    
    refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
    origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    if(origvaluesD100==NULL)
      {
        origvaluesD100 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD100[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D100, origvaluesD100);
      }
  
    refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
    origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    if(origvaluesD010==NULL)
      {
        origvaluesD010 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD010[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D010, origvaluesD010);
      }
    
    refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
    origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    if(origvaluesD001==NULL)
    {
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD001[j] = aux+j*N_Functs*BaseVectDim;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D001, origvaluesD001);
    }
  

    if(Needs2ndDer[i])
    {
      SecondDer = TRUE;
      
      refvaluesD200 = TFEDatabase3D::GetRefElementValues(BaseFunct,
                                                         QuadFormula, D200);
      origvaluesD200 = TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
      if(origvaluesD200 == NULL)
      {
        origvaluesD200 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(int j = 0; j < MaxN_QuadPoints_3D; j++)
          origvaluesD200[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D200, 
                                                 origvaluesD200);
      }
      
      refvaluesD110 = TFEDatabase3D::GetRefElementValues(BaseFunct,
                                                         QuadFormula, D110);
      origvaluesD110 = TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
      if(origvaluesD110 == NULL)
      {
        origvaluesD110 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(int j = 0; j < MaxN_QuadPoints_3D; j++)
          origvaluesD110[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D110,
                                                 origvaluesD110);
      }
      
      refvaluesD101 = TFEDatabase3D::GetRefElementValues(BaseFunct,
                                                         QuadFormula, D101);
      origvaluesD101 = TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
      if(origvaluesD101==NULL)
      {
        origvaluesD101 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(int j = 0; j < MaxN_QuadPoints_3D; j++)
          origvaluesD101[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D101,
                                                 origvaluesD101);
      }
    
      refvaluesD020 = TFEDatabase3D::GetRefElementValues(BaseFunct, 
                                                         QuadFormula, D020);
      origvaluesD020 = TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
      if(origvaluesD020==NULL)
      {
        origvaluesD020 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(int j = 0; j < MaxN_QuadPoints_3D; j++)
          origvaluesD020[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D020,
                                                 origvaluesD020);
      }
      
      refvaluesD011 = TFEDatabase3D::GetRefElementValues(BaseFunct,
                                                         QuadFormula, D011);
      origvaluesD011 = TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
      if(origvaluesD011==NULL)
      {
        origvaluesD011 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(int j = 0; j < MaxN_QuadPoints_3D; j++)
          origvaluesD011[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D011,
                                                 origvaluesD011);
      }
      
      refvaluesD002 = TFEDatabase3D::GetRefElementValues(BaseFunct,
                                                         QuadFormula, D002);
      origvaluesD002 = TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
      if(origvaluesD002==NULL)
      {
        origvaluesD002 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(int j = 0; j < MaxN_QuadPoints_3D; j++)
          origvaluesD002[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D002,
                                                 origvaluesD002);
      }
    } // endfor Needs2ndDer[i]
  } // endfor i
  
  // D100, D010 and D001
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
    int BaseVectDim =
        TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
    
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
      
      // inverse of transformation (DF^{-1}) multiplied by determinant^{-1}
      // for the inverse you can use maxima and type
      // A:matrix([x1, x2, x3],[y1, y2, y3],[z1, z2, z3]);
      // B:invert(A);
      // ratsimp(A.B); // check
      double i11 = (yc2*zc3 - yc3*zc2) * rec_detjk;
      double i12 = (yc3*zc1 - yc1*zc3) * rec_detjk;
      double i13 = (yc1*zc2 - yc2*zc1) * rec_detjk;
      double i21 = (xc3*zc2 - xc2*zc3) * rec_detjk;
      double i22 = (xc1*zc3 - xc3*zc1) * rec_detjk;
      double i23 = (xc2*zc1 - xc1*zc2) * rec_detjk;
      double i31 = (xc2*yc3 - xc3*yc2) * rec_detjk;
      double i32 = (xc3*yc1 - xc1*yc3) * rec_detjk;
      double i33 = (xc1*yc2 - xc2*yc1) * rec_detjk;
      
      if(BaseVectDim == 1)
      {
        for(int k = 0; k < N_Functs; k++)
        {
          origD100[k] = i11 * refD100[k] + i12 * refD010[k] + i13 * refD001[k];
          origD010[k] = i21 * refD100[k] + i22 * refD010[k] + i23 * refD001[k];
          origD001[k] = i31 * refD100[k] + i32 * refD010[k] + i33 * refD001[k];
        } // endfor k
      }
      else
      {
        // Piola transformation
        // this is DF divided by determinant
        double a11 = xc1 * rec_detjk;
        double a12 = xc2 * rec_detjk;
        double a13 = xc3 * rec_detjk;
        double a21 = yc1 * rec_detjk;
        double a22 = yc2 * rec_detjk;
        double a23 = yc3 * rec_detjk;
        double a31 = zc1 * rec_detjk;
        double a32 = zc2 * rec_detjk;
        double a33 = zc3 * rec_detjk;
        for(int k = 0; k < N_Functs; k++)
        {
          double p11 = a11 * refD100[k             ] 
                    + a12 * refD100[k + N_Functs  ] 
                    + a13 * refD100[k + 2*N_Functs];
          double p21 = a11 * refD010[k             ] 
                    + a12 * refD010[k + N_Functs  ] 
                    + a13 * refD010[k + 2*N_Functs];
          double p31 = a11 * refD001[k             ] 
                    + a12 * refD001[k + N_Functs  ] 
                    + a13 * refD001[k + 2*N_Functs];
          origD100[k] = i11*p11 + i12*p21 + i13*p31;
          origD010[k] = i21*p11 + i22*p21 + i23*p31;
          origD001[k] = i31*p11 + i32*p21 + i33*p31;
          
          double p12 = a21 * refD100[k             ] 
                    + a22 * refD100[k + N_Functs  ] 
                    + a23 * refD100[k + 2*N_Functs];
          double p22 = a21 * refD010[k             ] 
                    + a22 * refD010[k + N_Functs  ] 
                    + a23 * refD010[k + 2*N_Functs];
          double p32 = a21 * refD001[k             ] 
                    + a22 * refD001[k + N_Functs  ] 
                    + a23 * refD001[k + 2*N_Functs];
          origD100[k + N_Functs] = i11*p12 + i12*p22 + i13*p32;
          origD010[k + N_Functs] = i21*p12 + i22*p22 + i23*p32;
          origD001[k + N_Functs] = i31*p12 + i32*p22 + i33*p32;
          
          double p13 = a31 * refD100[k             ] 
                    + a32 * refD100[k + N_Functs  ] 
                    + a33 * refD100[k + 2*N_Functs];
          double p23 = a31 * refD010[k             ] 
                    + a32 * refD010[k + N_Functs  ] 
                    + a33 * refD010[k + 2*N_Functs];
          double p33 = a31 * refD001[k             ] 
                    + a32 * refD001[k + N_Functs  ] 
                    + a33 * refD001[k + 2*N_Functs];
          origD100[k + 2*N_Functs] = i11*p13 + i12*p23 + i13*p33;
          origD010[k + 2*N_Functs] = i21*p13 + i22*p23 + i23*p33;
          origD001[k + 2*N_Functs] = i31*p13 + i32*p23 + i33*p33;
        } // end for k
      }
    } // endfor j
  } // endfor i
  
  // leave if no second derivatives are needed
  if(!SecondDer) return;
  // find transformation matrix for second derivatives
  // do this only once since matrix is constant for an affine mapping
  
  // reset matrices
  double GeoData[6][6];
  double Eye[6][6];
  // reset matrices
  memset(GeoData, 0, 36*SizeOfDouble);
  memset(Eye, 0, 36*SizeOfDouble);
  Eye[0][0] = 1;
  Eye[1][1] = 1;
  Eye[2][2] = 1;
  Eye[3][3] = 1;
  Eye[4][4] = 1;
  Eye[5][5] = 1;
  
  GeoData[0][0] = xc1*xc1;
  GeoData[0][1] = 2*xc1*yc1;
  GeoData[0][2] = 2*xc1*zc1;
  GeoData[0][3] = yc1*yc1;
  GeoData[0][4] = 2*yc1*zc1;
  GeoData[0][5] = zc1*zc1;
  
  GeoData[1][0] = xc1*xc2;
  GeoData[1][1] = yc1*xc2 + xc1*yc2; 
  GeoData[1][2] = xc1*zc2 + xc2*zc1;
  GeoData[1][3] = yc1*yc2;
  GeoData[1][4] = yc1*zc2 + yc2*zc1;
  GeoData[1][5] = zc1*zc2;
  
  GeoData[2][0] = xc1*xc3;
  GeoData[2][1] = xc1*yc3 + xc3*yc1;
  GeoData[2][2] = xc1*zc3 + xc3*zc1;
  GeoData[2][3] = yc1*yc3;
  GeoData[2][4] = yc1*zc3 + yc3*zc1;
  GeoData[2][5] = zc1*zc3;
  
  GeoData[3][0] = xc2*xc2;
  GeoData[3][1] = 2*xc2*yc2;
  GeoData[3][2] = 2*xc2*zc2;
  GeoData[3][3] = yc2*yc2;
  GeoData[3][4] = 2*yc2*zc2;
  GeoData[3][5] = zc2*zc2;
  
  GeoData[4][0] = xc2*xc3;
  GeoData[4][1] = xc2*yc3 + xc3*yc2;
  GeoData[4][2] = xc2*zc3 + xc3*zc2;
  GeoData[4][3] = yc2*yc3;
  GeoData[4][4] = yc2*zc3 + yc3*zc2;
  GeoData[4][5] = zc2*zc3;
  
  GeoData[5][0] = xc3*xc3;
  GeoData[5][1] = 2*xc3*yc3;
  GeoData[5][2] = 2*xc3*zc3;
  GeoData[5][3] = yc3*yc3;
  GeoData[5][4] = 2*yc3*zc3;
  GeoData[5][5] = zc3*zc3;
  
  // invert GeoData and store result in Eye
  SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 6, 6, 6, 6);

  for(int i = 0; i < N_Sets; i++)
  {
    if(Needs2ndDer[i])
    {
      BaseFunct=BaseFuncts[i];
      N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
      int BaseVectDim =
        TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
      if(BaseVectDim != 1)
      {
        ErrMsg("second derivatives of vector valued basis function not " <<
               "supported");
        exit(0);
      }
      
      refvaluesD200 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D200);
      origvaluesD200 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D200);
      
      refvaluesD110 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D110);
      origvaluesD110 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D110);
      
      refvaluesD101 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D101);
      origvaluesD101 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D101);
      
      refvaluesD020 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D020);
      origvaluesD020 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D020);
      
      refvaluesD011 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D011);
      origvaluesD011 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D011);
      
      refvaluesD002 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D002);
      origvaluesD002 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D002);
      
      for(int j = 0; j < N_Points; j++)
      {
        double *refD200 = refvaluesD200[j];
        double *refD110 = refvaluesD110[j];
        double *refD101 = refvaluesD101[j];
        double *refD020 = refvaluesD020[j];
        double *refD011 = refvaluesD011[j];
        double *refD002 = refvaluesD002[j];
        
        double *origD200 = origvaluesD200[j];
        double *origD110 = origvaluesD110[j];
        double *origD101 = origvaluesD101[j];
        double *origD020 = origvaluesD020[j];
        double *origD011 = origvaluesD011[j];
        double *origD002 = origvaluesD002[j];
        
        for(int k = 0; k < N_Functs; k++)
        {
          double r200 = refD200[k];
          double r110 = refD110[k];
          double r101 = refD101[k];
          double r020 = refD020[k];
          double r011 = refD011[k];
          double r002 = refD002[k];
          
          double o200, o110, o101, o020, o011, o002;
          o200 =  Eye[0][0]*r200 + Eye[0][1]*r110 + Eye[0][2]*r101 
                + Eye[0][3]*r020 + Eye[0][4]*r011 + Eye[0][5]*r002;
          o110 =  Eye[1][0]*r200 + Eye[1][1]*r110 + Eye[1][2]*r101 
                + Eye[1][3]*r020 + Eye[1][4]*r011 + Eye[1][5]*r002;
          o101 =  Eye[2][0]*r200 + Eye[2][1]*r110 + Eye[2][2]*r101 
                + Eye[2][3]*r020 + Eye[2][4]*r011 + Eye[2][5]*r002;
          o020 =  Eye[3][0]*r200 + Eye[3][1]*r110 + Eye[3][2]*r101 
                + Eye[3][3]*r020 + Eye[3][4]*r011 + Eye[3][5]*r002;
          o011 =  Eye[4][0]*r200 + Eye[4][1]*r110 + Eye[4][2]*r101 
                + Eye[4][3]*r020 + Eye[4][4]*r011 + Eye[4][5]*r002;
          o002 =  Eye[5][0]*r200 + Eye[5][1]*r110 + Eye[5][2]*r101 
                + Eye[5][3]*r020 + Eye[5][4]*r011 + Eye[5][5]*r002;
          
          origD200[k] = o200;
          origD110[k] = o110;
          origD101[k] = o101;
          origD020[k] = o020;
          origD011[k] = o011;
          origD002[k] = o002;
        } // endfor k
      } // endif
    } // endfor j
  } // endfor i
}

/** calculate functions and derivatives from reference element
    to original element */
void TTetraAffin::GetOrigValues(double xi, double eta, double zeta,
                int N_BaseFunct,
                double *uref, double *uxiref, double *uetaref, double *uzetaref,
                double *uorig, double *uxorig, double *uyorig, double *uzorig,
                int _BaseVectDim)
{
  int i;

  if(_BaseVectDim == 1)
  {
    // D000
    for(i=0;i<N_BaseFunct;i++)
      uorig[i] = uref[i];
    
    // D100, D010 and D001
    for(i=0;i<N_BaseFunct;i++)
    {
      uxorig[i] = ((yc2*zc3 - yc3*zc2)*uxiref[i] + (yc3*zc1 - yc1*zc3)*uetaref[i] + (yc1*zc2 - yc2*zc1)*uzetaref[i]) *rec_detjk;
      uyorig[i] = ((xc3*zc2 - xc2*zc3)*uxiref[i] + (xc1*zc3 - xc3*zc1)*uetaref[i] + (xc2*zc1 - xc1*zc2)*uzetaref[i]) *rec_detjk;
      uzorig[i] = ((xc2*yc3 - xc3*yc2)*uxiref[i] + (xc3*yc1 - xc1*yc3)*uetaref[i] + (xc1*yc2 - xc2*yc1)*uzetaref[i]) *rec_detjk;
    } // endfor i
  }
  else
  {
    // D000
    PiolaMapOrigFromRef(N_BaseFunct, uref, uorig);
    
    // D100, D010, D001
    // not yet implemented
  }
}

/** calculate functions and derivatives from reference element
        to original element on joint, parameters on joint are p1, p2 */
void TTetraAffin::GetOrigValuesJoint(int JointNr, double p1, double p2, int N_BaseFunct,
          double *uref, double *uxiref, double *uetaref, double *uzetaref,
          double *uorig, double *uxorig, double *uyorig, double *uzorig)
{
  int i;

//   cout << "right!" << endl;
  // D000
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  // D100, D010 and D001
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((yc2*zc3 - yc3*zc2)*uxiref[i] + (yc3*zc1 - yc1*zc3)*uetaref[i] + (yc1*zc2 - yc2*zc1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((xc3*zc2 - xc2*zc3)*uxiref[i] + (xc1*zc3 - xc3*zc1)*uetaref[i] + (xc2*zc1 - xc1*zc2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((xc2*yc3 - xc3*yc2)*uxiref[i] + (xc3*yc1 - xc1*yc3)*uetaref[i] + (xc1*yc2 - xc2*yc1)*uzetaref[i]) *rec_detjk;
 } // endfor i
}

// for compatibility
void TTetraAffin::GetOrigValues(int JointNr, double p1, double p2, int N_BaseFunct,
          double *uref, double *uxiref, double *uetaref, double *uzetaref,
          double *uorig, double *uxorig, double *uyorig, double *uzorig)
{
  int i;

//   cout << "right!" << endl;
  // D000
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  // D100, D010 and D001
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((yc2*zc3 - yc3*zc2)*uxiref[i] + (yc3*zc1 - yc1*zc3)*uetaref[i] + (yc1*zc2 - yc2*zc1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((xc3*zc2 - xc2*zc3)*uxiref[i] + (xc1*zc3 - xc3*zc1)*uetaref[i] + (xc2*zc1 - xc1*zc2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((xc2*yc3 - xc3*yc2)*uxiref[i] + (xc3*yc1 - xc1*yc3)*uetaref[i] + (xc1*yc2 - xc2*yc1)*uzetaref[i]) *rec_detjk;
 } // endfor i
}

void TTetraAffin::SetCell(TBaseCell *cell)
{
  int i;

  Cell = cell;

  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);

  xc0=x0;
  xc1=x1-x0;
  xc2=x2-x0;
  xc3=x3-x0;

  yc0=y0;
  yc1=y1-y0;
  yc2=y2-y0;
  yc3=y3-y0;

  zc0=z0;
  zc1=z1-z0;
  zc2=z2-z0;
  zc3=z3-z0;

  detjk= xc1*yc2*zc3 + xc2*yc3*zc1 + xc3*yc1*zc2
        -xc3*yc2*zc1 - xc2*yc1*zc3 - xc1*yc3*zc2;
  
  rec_detjk = 1/detjk;
}

/** return outer normal unit vector */
void TTetraAffin::GetOuterNormal(int j, double s, double t,
                                 double &n1, double &n2, double &n3)
{
//   double len;

  switch(j)
  {
    case 0:
      n1 = yc2*zc1 - zc2*yc1;
      n2 = zc2*xc1 - xc2*zc1;
      n3 = xc2*yc1 - yc2*xc1;
    break;

    case 1:
      n1 = yc1*zc3 - zc1*yc3;
      n2 = zc1*xc3 - xc1*zc3;
      n3 = xc1*yc3 - yc1*xc3;
    break;

    case 2:
      n1 = (yc3-yc2)*(zc1-zc2) - (zc3-zc2)*(yc1-yc2);
      n2 = (zc3-zc2)*(xc1-xc2) - (xc3-xc2)*(zc1-zc2);
      n3 = (xc3-xc2)*(yc1-yc2) - (yc3-yc2)*(xc1-xc2);
    break;

    case 3:
      n1 = yc3*zc2 - zc3*yc2;
      n2 = zc3*xc2 - xc3*zc2;
      n3 = xc3*yc2 - yc3*xc2;
    break;

    default:
      Error("Wrong local joint number" << endl);
      n1 = -1; n2 = -1; n3 = -1;
      return;
  }
  // without scaling needed for finding area, Sashi, 20. Jan 2012
//   len = sqrt(n1*n1 + n2*n2 + n3*n3);

//   n1 /= len;
//   n2 /= len;
//   n3 /= len;
}

/** return two tangent vectors */
void TTetraAffin::GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23)
{
  switch(j)
  {
    case 0:
      t11 = xc2; t12 = yc2; t13 = zc2;
      t21 = xc1; t22 = yc1; t23 = zc1;
    break;

    case 1:
      t11 = xc1; t12 = yc1; t13 = zc1;
      t21 = xc3; t22 = yc3; t23 = zc3;
    break;

    case 2:
      t11 = (xc3-xc2); t12 = (yc3-yc2); t13 = (zc3-zc2);
      t21 = (xc1-xc2); t22 = (yc1-yc2); t23 = (zc1-zc2);
    break;

    case 3:
      t11 = xc3; t12 = yc3; t13 = zc3;
      t21 = xc2; t22 = yc2; t23 = zc2;
    break;

    default:
      Error("Wrong joint number for tetrahedron! " << endl);
      return;
  }
} // end TTetraAffin::GetTangentVectors


/** Piola transformation for vectorial basis functions */
void TTetraAffin::PiolaMapOrigFromRef(int N_Functs, double *refD000, 
                                      double *origD000 )
{
  double a11 = xc1 * rec_detjk;
  double a12 = xc2 * rec_detjk;
  double a13 = xc3 * rec_detjk;
  double a21 = yc1 * rec_detjk;
  double a22 = yc2 * rec_detjk;
  double a23 = yc3 * rec_detjk;
  double a31 = zc1 * rec_detjk;
  double a32 = zc2 * rec_detjk;
  double a33 = zc3 * rec_detjk;
  for(int k = 0; k < N_Functs; k++)
  {
    // three components:
    origD000[k             ] = a11 * refD000[k             ] 
                             + a12 * refD000[k +   N_Functs] 
                             + a13 * refD000[k + 2*N_Functs]; 
    origD000[k + N_Functs  ] = a21 * refD000[k             ] 
                             + a22 * refD000[k +   N_Functs] 
                             + a23 * refD000[k + 2*N_Functs]; 
    origD000[k + 2*N_Functs] = a31 * refD000[k             ] 
                             + a32 * refD000[k +   N_Functs] 
                             + a33 * refD000[k + 2*N_Functs]; 
  }
}
