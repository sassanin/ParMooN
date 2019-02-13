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
// @(#)QuadAffin.C        1.8 04/13/00
//
// Class:      TQuadAffin
//
// Purpose:    affin reference transformations for parallelogram
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#include <QuadAffin.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>

/** constuctor */
TQuadAffin::TQuadAffin()
{
}

/** transfer from reference element to original element */
void TQuadAffin::GetOrigFromRef(double xi, double eta, double &X, double &Y)
{
  X = xc0 + xc1*xi + xc2*eta;
  Y = yc0 + yc1*xi + yc2*eta;
}

/** transfer a set of point from reference to original element */
void TQuadAffin::GetOrigFromRef(int N_Points, double *xi, double *eta,
                                double *X, double *Y, double *absdetjk)
{
  int i;
  double Xi, Eta;
  double absdet=fabs(detjk);
  
  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    X[i] = xc0 + xc1*Xi + xc2*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta;
    absdetjk[i] = absdet;
  }
}

/** transfer from reference element to original element */
void TQuadAffin::GetOrigFromRef(double *ref, double *orig)
{
  orig[0]=xc0 + xc1*ref[0] + xc2*ref[1];
  orig[1]=yc0 + yc1*ref[0] + yc2*ref[1];
}

/** transfer from original element to reference element */
void TQuadAffin::GetRefFromOrig(double X, double Y, double &xi, double &eta)
{
  double xt=(X - xc0)/detjk;
  double yt=(Y - yc0)/detjk;

  xi  =  yc2*xt - xc2*yt;
  eta = -yc1*xt + xc1*yt;
}

/** transfer from original element to reference element */
void TQuadAffin::GetRefFromOrig(double *orig, double *ref)
{
  double xt=(orig[0] - xc0)/detjk;
  double yt=(orig[1] - yc0)/detjk;

  ref[0] =  yc2*xt - xc2*yt;
  ref[1] = -yc1*xt + xc1*yt;
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadAffin::GetOrigValues(BaseFunct2D BaseFunct,
                               int N_Points, double *xi, double *eta,
                               int N_Functs, QuadFormula2D QuadFormula)
{
  int i,j,k, BaseVectDim;
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
  double *aux;
//   double AllData[MaxN_BaseFunctions2D][5];
  double GeoData[5][5];

  BaseVectDim = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();

  refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D00);
  if(refvaluesD00==NULL)
  {
    TFEDatabase2D::GetBaseFunct2D(BaseFunct)->MakeRefElementData(QuadFormula);
    refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                                  QuadFormula, D00);
  }

  origvaluesD00=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
  if(origvaluesD00==NULL)
  {
    origvaluesD00 = new double* [MaxN_QuadPoints_2D];
    aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D*BaseVectDim];
    for(i=0;i<MaxN_QuadPoints_2D;i++)
      origvaluesD00[i] = aux+i*MaxN_BaseFunctions2D*BaseVectDim;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D00, origvaluesD00);
  }

  refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D10);
  origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
  if(origvaluesD10==NULL)
  {
    origvaluesD10 = new double* [MaxN_QuadPoints_2D];
    aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D*BaseVectDim];
    for(i=0;i<MaxN_QuadPoints_2D;i++)
      origvaluesD10[i] = aux+i*MaxN_BaseFunctions2D*BaseVectDim;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D10, origvaluesD10);
  }

  refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D01);
  origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);
  if(origvaluesD01==NULL)
  {
    origvaluesD01 = new double* [MaxN_QuadPoints_2D];
    aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D*BaseVectDim];
    for(i=0;i<MaxN_QuadPoints_2D;i++)
      origvaluesD01[i] = aux+i*MaxN_BaseFunctions2D*BaseVectDim;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D01, origvaluesD01);
  }

  refvaluesD20=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D20);
  origvaluesD20=TFEDatabase2D::GetOrigElementValues(BaseFunct, D20);
  if(origvaluesD20==NULL)
  {
    origvaluesD20 = new double* [MaxN_QuadPoints_2D];
    aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D*BaseVectDim];
    for(i=0;i<MaxN_QuadPoints_2D;i++)
      origvaluesD20[i] = aux+i*MaxN_BaseFunctions2D*BaseVectDim;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D20, origvaluesD20);
  }

  refvaluesD11=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D11);
  origvaluesD11=TFEDatabase2D::GetOrigElementValues(BaseFunct, D11);
  if(origvaluesD11==NULL)
  {
    origvaluesD11 = new double* [MaxN_QuadPoints_2D];
    aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D*BaseVectDim];
    for(i=0;i<MaxN_QuadPoints_2D;i++)
      origvaluesD11[i] = aux+i*MaxN_BaseFunctions2D*BaseVectDim;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D11, origvaluesD11);
  }

  refvaluesD02=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D02);
  origvaluesD02=TFEDatabase2D::GetOrigElementValues(BaseFunct, D02);
  if(origvaluesD02==NULL)
  {
    origvaluesD02 = new double* [MaxN_QuadPoints_2D];
    aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D*BaseVectDim];
    for(i=0;i<MaxN_QuadPoints_2D;i++)
      origvaluesD02[i] = aux+i*MaxN_BaseFunctions2D*BaseVectDim;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D02, origvaluesD02);
  }

  // D00
  for(i=0;i<N_Points;i++)
  {
    refD00 = refvaluesD00[i];
    origD00 = origvaluesD00[i];

    if(BaseVectDim == 1)
    {
      for(j=0;j<N_Functs*BaseVectDim;j++)
      {
        origD00[j] = refD00[j];
      } // endfor j
    }
    else
    {
      PiolaMapOrigFromRef(N_Functs, refD00, origD00);
    }
  } // endfor i

// /*
  // D10 and D01
  for(i=0;i<N_Points;i++)
  {
    refD10 = refvaluesD10[i];
    origD10 = origvaluesD10[i];

    refD01 = refvaluesD01[i];
    origD01 = origvaluesD01[i];

    if(BaseVectDim == 1)
    {
      for(j=0;j<N_Functs*BaseVectDim;j++)
      {
        origD10[j]=(yc2*refD10[j]-yc1*refD01[j]) *rec_detjk;
        origD01[j]=(-xc2*refD10[j]+xc1*refD01[j]) *rec_detjk;
  
        // cout << "10: " << origD10[j] << endl;
        // cout << "01: " << origD01[j] << endl;
        // cout << endl;
      } // endfor j
    }
    else
    {
      PiolaMapOrigFromRef(N_Functs, refD10, refD01, origD10, origD01);
    }
    // cout << "----------" << endl;
  } // endfor i
// */
  
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
void TQuadAffin::GetOrigValues(int N_Sets, BaseFunct2D *BaseFuncts,
                               int N_Points, double *xi, double *eta,
                               QuadFormula2D QuadFormula,
                               bool *Needs2ndDer)
{
  int i,j,k,N_, start, end, BaseVectDim;
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
  int ii,ij,ik;
  double tmp,Eye1[3][3];

  SecondDer = false;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    BaseVectDim = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();

    N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();

    refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D00);
    if(refvaluesD00==NULL)
    {
      TFEDatabase2D::GetBaseFunct2D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                                    QuadFormula, D00);
    }
  
    origvaluesD00=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    if(origvaluesD00==NULL)
    {
      origvaluesD00 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*N_Functs*BaseVectDim];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD00[j] = aux+j*N_Functs*BaseVectDim;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D00, origvaluesD00);
    }
  
    for(j=0;j<N_Points;j++)
    {
      refD00 = refvaluesD00[j];
      origD00 = origvaluesD00[j];
       
      if(BaseVectDim == 1)
        memcpy(origD00, refD00, N_Functs*BaseVectDim*SizeOfDouble);
      else
        this->PiolaMapOrigFromRef(N_Functs, refD00, origD00);
    } // endfor j

    refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D10);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    if(origvaluesD10==NULL)
    {
      origvaluesD10 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*N_Functs*BaseVectDim];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD10[j] = aux+j*N_Functs*BaseVectDim;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D10, origvaluesD10);
    }
  
    refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D01);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);
    if(origvaluesD01==NULL)
    {
      origvaluesD01 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*N_Functs*BaseVectDim];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD01[j] = aux+j*N_Functs*BaseVectDim;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D01, origvaluesD01);
    }
  
    if(Needs2ndDer[i])
    {
      SecondDer = true;

      refvaluesD20=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D20);
      origvaluesD20=TFEDatabase2D::GetOrigElementValues(BaseFunct, D20);
      if(origvaluesD20==NULL)
      {
        origvaluesD20 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD20[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D20, origvaluesD20);
      }
    
      refvaluesD11=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D11);
      origvaluesD11=TFEDatabase2D::GetOrigElementValues(BaseFunct, D11);
      if(origvaluesD11==NULL)
      {
        origvaluesD11 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD11[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D11, origvaluesD11);
      }
    
      refvaluesD02=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D02);
      origvaluesD02=TFEDatabase2D::GetOrigElementValues(BaseFunct, D02);
      if(origvaluesD02==NULL)
      {
        origvaluesD02 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD02[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D02, origvaluesD02);
      }
    } // endfor Needs2ndDer[i]
  } // endfor i
  
  // D10 and D01
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();
    BaseVectDim = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();
    
    refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D10);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);

    refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D01);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);

    for(j=0;j<N_Points;j++)
    {
      refD10 = refvaluesD10[j];
      origD10 = origvaluesD10[j];
  
      refD01 = refvaluesD01[j];
      origD01 = origvaluesD01[j];
  
      if(BaseVectDim == 1)
      {
        for(k=0;k<N_Functs;k++)
        {
          origD10[k]=(yc2*refD10[k]-yc1*refD01[k]) *rec_detjk;
          origD01[k]=(-xc2*refD10[k]+xc1*refD01[k]) *rec_detjk;
        } // endfor k
      }
      else
      {
        this->PiolaMapOrigFromRef(N_Functs, refD10, refD01, origD10, origD01);
      }
    } // endfor j
  } // endfor i

  // leave if no second derivatives are needed
  if(!SecondDer) return;
  
  // find transformation matrix for second derivatives
  // do this only once since matrix is constant for an affine mapping

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
      N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();
      BaseVectDim = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();
    
      refvaluesD20=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D20);
      origvaluesD20=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D20);
  
      refvaluesD11=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D11);
      origvaluesD11=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D11);
  
      refvaluesD02=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D02);
      origvaluesD02=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D02);

      for(j=0;j<N_Points;j++)
      {
        refD20 = refvaluesD20[j];
        refD11 = refvaluesD11[j];
        refD02 = refvaluesD02[j];
    
        origD20 = origvaluesD20[j];
        origD11 = origvaluesD11[j];
        origD02 = origvaluesD02[j];
    
        for(k=0;k<N_Functs*BaseVectDim;k++)
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
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadAffin::GetOrigValues(double xi, double eta, int N_BaseFunct,
            double *uref, double *uxiref, double *uetaref,
            double *uorig, double *uxorig, double *uyorig, int _BaseVectDim)
{
  if(_BaseVectDim == 1)
  {
    for(int i = 0; i < N_BaseFunct; i++) 
    {
      // D00
      uorig[i] = uref[i];
      // D10 and D01
      if (uxiref!=NULL && uetaref!=NULL && uxorig!=NULL && uyorig!=NULL)
      {
        uxorig[i]=(yc2*uxiref[i]-yc1*uetaref[i]) * rec_detjk;
        uyorig[i]=(-xc2*uxiref[i]+xc1*uetaref[i]) * rec_detjk;
      }
    }
  }
  else if(_BaseVectDim == 2)
  {
    // D00
    this->PiolaMapOrigFromRef(N_BaseFunct, uref, uorig);
    // D10 and D01
    if (uxiref!=NULL && uetaref!=NULL && uxorig!=NULL && uyorig!=NULL)
    {
      this->PiolaMapOrigFromRef(N_BaseFunct, uxiref, uetaref, uxorig, uyorig);
    }
  }
}

void TQuadAffin::GetOrigValues(int joint, double zeta,
             int N_BaseFunct,
             double *uref, double *uxiref, double *uetaref,
             double *uorig, double *uxorig, double *uyorig,
             int _BaseVectDim)
{
  double xi, eta;
  switch(joint)
  {
    case 0:
      xi = zeta; eta = -1;
    break;

    case 1:
      xi = 1; eta = zeta;
    break;

    case 2:
      xi  = -zeta; eta = 1;
    break;

    case 3:
      xi = -1 ; eta = -zeta;
    break;
   }
   this->GetOrigValues(xi, eta, N_BaseFunct, uref, uxiref, uetaref, 
                       uorig, uxorig, uyorig, _BaseVectDim);
}

void TQuadAffin::SetCell(TBaseCell *cell)
{
  int i;
#ifdef __3D__
  double z0, z1, z2, z3;
#endif

  Cell = cell;

#ifdef __3D__
  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);
#else
  Cell->GetVertex(0)->GetCoords(x0, y0);
  Cell->GetVertex(1)->GetCoords(x1, y1);
  Cell->GetVertex(2)->GetCoords(x2, y2);
  Cell->GetVertex(3)->GetCoords(x3, y3);
#endif

  xc0 = (x1 + x3) * 0.5;
  xc1 = (x1 - x0) * 0.5;
  xc2 = (x3 - x0) * 0.5;

  yc0 = (y1 + y3) * 0.5;
  yc1 = (y1 - y0) * 0.5;
  yc2 = (y3 - y0) * 0.5;

  detjk=xc1*yc2-xc2*yc1;
  rec_detjk=1/detjk;
}

/** return outer normal vector */
void TQuadAffin::GetOuterNormal(int j, double zeta,
                                double &n1, double &n2)
{
  double len;

  switch(j)
  {
    case 0:
      n1 = y1-y0;
      n2 = x0-x1;
    break;

    case 1:
      n1 = y2-y1;
      n2 = x1-x2;
    break;

    case 2:
      n1 = y3-y2;
      n2 = x2-x3;
    break;

    case 3:
      n1 = y0-y3;
      n2 = x3-x0;
    break;

    default:
      cerr << "wrong joint number" << endl;
      n1 = -1;
      n2 = -1;
  } // endswitch

  len = sqrt(n1*n1+n2*n2);

  n1 /= len;
  n2 /= len;
}

/** return derivative of boundary parametrization */
void TQuadAffin::GetTangent(int j, double zeta,
                                double &t1, double &t2)
{
  // factor 0.5 since edge parameter runs from -1 to +1
  switch(j)
  {
    case 0:
      t1 = 0.5*(x1-x0);
      t2 = 0.5*(y1-y0);
    break;

    case 1:
      t1 = 0.5*(x2-x1);
      t2 = 0.5*(y2-y1);
    break;

    case 2:
      t1 = 0.5*(x3-x2);
      t2 = 0.5*(y3-y2);
    break;

    case 3:
      t1 = 0.5*(x0-x3);
      t2 = 0.5*(y0-y3);
    break;

    default:
      cerr << "wrong joint number" << endl;
      t1 = -1;
      t2 = -1;
  } // endswitch
}

/** return volume of cell */
double TQuadAffin::GetVolume()
{
  return 4*detjk;
}



/** transfer a  set of boundary points from reference to original element */
void TQuadAffin::GetOrigBoundFromRef(int joint, int N_Points, double *zeta, double *X, double *Y)
{
  int i;
  double Xi, Eta;


 for(i=0;i<N_Points;i++)
 {
  switch(joint)
  {
    case 0:
      Xi  = zeta[i]; Eta = -1;
    break;

    case 1:
      Xi = 1; Eta = zeta[i];
    break;

    case 2:
      Xi  = -zeta[i]; Eta = 1;
    break;

    case 3:
      Xi = -1 ; Eta = -zeta[i];
    break;
   }
   
   X[i] = xc0 + xc1*Xi + xc2*Eta;
   Y[i] = yc0 + yc1*Xi + yc2*Eta;    
  }

}


/** Piola transformation for vectorial basis functions */
void TQuadAffin::PiolaMapOrigFromRef(int N_Functs, double *refD00, double *origD00 )
{
  double a11 = xc1*rec_detjk;
  double a12 = xc2*rec_detjk;
  double a21 = yc1*rec_detjk;
  double a22 = yc2*rec_detjk;
  for(int k = 0; k < N_Functs; k++)
  {
    // Piola transformation
    // phi = 1/|J| DF phi_hat
    // def.gradient/detjk
    origD00[k] = a11*refD00[k] + a12*refD00[N_Functs+k];
    origD00[N_Functs+k] = a21*refD00[k] + a22*refD00[N_Functs+k];

  }
}
   
void TQuadAffin::PiolaMapOrigFromRef(int N_Functs, double *refD10, 
                                     double *refD01, double *origD10, 
                                     double *origD01)
{
  double a11 = xc1*rec_detjk;
  double a12 = xc2*rec_detjk;
  double a21 = yc1*rec_detjk;
  double a22 = yc2*rec_detjk;
  for(int k = 0; k < N_Functs; k++)
  {
    // Piola transformation
    // phi = 1/|J| DF phi_hat
    // def.gradient/detjk
    
    // x-component (k=0,N_Functs-1)
    double refPiolaD10_k = a11*refD10[k] + a12*refD10[N_Functs+k];
    double refPiolaD01_k = a11*refD01[k] + a12*refD01[N_Functs+k];
    origD10[k] = ( yc2*refPiolaD10_k - yc1*refPiolaD01_k) * rec_detjk;
    origD01[k] = (-xc2*refPiolaD10_k + xc1*refPiolaD01_k) * rec_detjk;
    
    // y-component (k=N_Functs,BaseVectDim*N_Functs-1)
    refPiolaD10_k = a21*refD10[k] + a22*refD10[N_Functs+k];
    refPiolaD01_k = a21*refD01[k] + a22*refD01[N_Functs+k];
    origD10[k+N_Functs] = ( yc2*refPiolaD10_k - yc1*refPiolaD01_k) * rec_detjk;
    origD01[k+N_Functs] = (-xc2*refPiolaD10_k + xc1*refPiolaD01_k) * rec_detjk;
  } // endfor k
}
