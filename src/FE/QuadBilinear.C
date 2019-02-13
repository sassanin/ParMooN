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
// @(#)QuadBilinear.C        1.11 04/13/00
//
// Class:      TQuadBilinear
//
// Purpose:    Bilinear reference transformations for parallelogram
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
//        :    10.06.2010 modified for vector basis functions (Sashikumaar Ganesan)
// =======================================================================

#include <QuadBilinear.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>

/** constuctor */
TQuadBilinear::TQuadBilinear()
{
}

/** transfer from reference element to original element */
void TQuadBilinear::GetOrigFromRef(double xi, double eta, double &X, double &Y)
{
  X = xc0 + xc1*xi + xc2*eta + xc3*xi*eta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*xi*eta;
}

/** transfer a set of point from reference to original element */
void TQuadBilinear::GetOrigFromRef(int N_Points, double *xi, double *eta,
                                double *X, double *Y, double *absdetjk)
{
  int i;
  double Xi, Eta;
  
  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    X[i] = xc0 + xc1*Xi + xc2*Eta +xc3*Eta*Xi;
    Y[i] = yc0 + yc1*Xi + yc2*Eta +yc3*Eta*Xi;
    absdetjk[i] = ABS( (xc1+xc3*Eta)*(yc2+yc3*Xi)
                      -(xc2+xc3*Xi)*(yc1+yc3*Eta));
  }
}

/** transfer from reference element to original element */
void TQuadBilinear::GetOrigFromRef(double *ref, double *orig)
{
  orig[0] = xc0 + xc1*ref[0] + xc2*ref[1] + xc3*ref[0]*ref[1];
  orig[1] = yc0 + yc1*ref[0] + yc2*ref[1] + yc3*ref[0]*ref[1];
}

/** transfer from original element to reference element */
void TQuadBilinear::GetRefFromOrig(double x, double y, double &xi, double &eta)
{
  double z0, z1, xi0, eta0;
  double recdetaffine;
  double eps=1e-14;

  recdetaffine = 1/(xc1*yc2-xc2*yc1);

  z0 = ( yc2 * (x-xc0) - xc2 * (y-yc0)) * recdetaffine;
  z1 = (-yc1 * (x-xc0) + xc1 * (y-yc0)) * recdetaffine;
  
  xi0 = z0;
  eta0 = z1;

  xi = xi0+10;
  eta = eta0+10;
  while( fabs(xi-xi0)+fabs(eta-eta0) > eps )
  {
    xi  = xi0;
    eta = eta0;
    xi0  = z0 - ( yc2 * xc3*xi*eta - xc2 * yc3*xi*eta) * recdetaffine;
    eta0 = z1 - (-yc1 * xc3*xi*eta + xc1 * yc3*xi*eta) * recdetaffine;
  }
  xi = xi0;
  eta = eta0;

  return;
}

/** transfer from original element to reference element */
void TQuadBilinear::GetRefFromOrig(double *orig, double *ref)
{
  GetRefFromOrig(orig[0], orig[1], ref[0], ref[1]);
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadBilinear::GetOrigValues(BaseFunct2D BaseFunct,
                               int N_Points, double *xi, double *eta,
                               int N_Functs, QuadFormula2D QuadFormula)
{
  int i,j,k, BaseVectDim, M;
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
  double rec_detjk, Xi, Eta;
  double AllData[MaxN_BaseFunctions2D][5];
  double GeoData[5][5];

  BaseVectDim = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();
  M = BaseVectDim*MaxN_QuadPoints_2D;
  
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
    origvaluesD00 = new double* [M];
    aux = new double [M*MaxN_BaseFunctions2D];
    for(i=0;i<M;i++)
      origvaluesD00[i] = aux+i*MaxN_BaseFunctions2D;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D00, origvaluesD00);
  }

  refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D10);
  origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
  if(origvaluesD10==NULL)
  {
    origvaluesD10 = new double* [M];
    aux = new double [M*MaxN_BaseFunctions2D];
    for(i=0;i<M;i++)
      origvaluesD10[i] = aux+i*MaxN_BaseFunctions2D;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D10, origvaluesD10);
  }

  refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D01);
  origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);
  if(origvaluesD01==NULL)
  {
    origvaluesD01 = new double* [M];
    aux = new double [M*MaxN_BaseFunctions2D];
    for(i=0;i<M;i++)
      origvaluesD01[i] = aux+i*MaxN_BaseFunctions2D;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D01, origvaluesD01);
  }

  refvaluesD20=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D20);
  origvaluesD20=TFEDatabase2D::GetOrigElementValues(BaseFunct, D20);
  if(origvaluesD20==NULL)
  {
    origvaluesD20 = new double* [M];
    aux = new double [M*MaxN_BaseFunctions2D];
    for(i=0;i<M;i++)
      origvaluesD20[i] = aux+i*MaxN_BaseFunctions2D;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D20, origvaluesD20);
  }

  refvaluesD11=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D11);
  origvaluesD11=TFEDatabase2D::GetOrigElementValues(BaseFunct, D11);
  if(origvaluesD11==NULL)
  {
    origvaluesD11 = new double* [M];
    aux = new double [M*MaxN_BaseFunctions2D];
    for(i=0;i<M;i++)
      origvaluesD11[i] = aux+i*MaxN_BaseFunctions2D;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D11, origvaluesD11);
  }

  refvaluesD02=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D02);
  origvaluesD02=TFEDatabase2D::GetOrigElementValues(BaseFunct, D02);
  if(origvaluesD02==NULL)
  {
    origvaluesD02 = new double* [M];
    aux = new double [M*MaxN_BaseFunctions2D];
    for(i=0;i<M;i++)
      origvaluesD02[i] = aux+i*MaxN_BaseFunctions2D;
    TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D02, origvaluesD02);
  }

  if(BaseVectDim == 1)
  {
    // D00
    for(i=0;i<N_Points;i++)
    {
      refD00 = refvaluesD00[i];
      origD00 = origvaluesD00[i];
      
      for(j=0;j<N_Functs;j++)
      {
        origD00[j] = refD00[j];
      } // endfor j
    } // endfor i
    
    // D10 and D01
    for(i=0;i<N_Points;i++)
    {
      refD10 = refvaluesD10[i];
      origD10 = origvaluesD10[i];
      
      refD01 = refvaluesD01[i];
      origD01 = origvaluesD01[i];
      
      for(j=0;j<N_Functs;j++)
      {
        Xi = xi[i];
        Eta = eta[i];
        rec_detjk = 1/( (xc1+xc3*Eta)*(yc2+yc3*Xi)
                      -(xc2+xc3*Xi)*(yc1+yc3*Eta));
        origD10[j]=( (yc2+yc3*Xi) * refD10[j]
                    -(yc1+yc3*Eta) * refD01[j]) * rec_detjk;
        origD01[j]=(-(xc2+xc3*Xi) * refD10[j]
                    +(xc1+xc3*Eta) * refD01[j]) * rec_detjk;
      } // endfor j
    } // endfor i
  }
  else 
  {
    // D00
    for(int j = 0; j < N_Points; j++)
    {
      refD00 = refvaluesD00[j];
      origD00 = origvaluesD00[j];
      
      refD10 = refvaluesD10[i];
      origD10 = origvaluesD10[i];

      refD01 = refvaluesD01[i];
      origD01 = origvaluesD01[i];
      
      this->PiolaMapOrigFromRefNotAffine(N_Functs, refD00, origD00, xi[j], 
                                           eta[j]);
      this->PiolaMapOrigFromRefNotAffine(N_Functs, refD00, refD10, refD01, 
                                         origD10, origD01, xi[j], eta[j]);
    }
  }
// */

/*
  for(i=0;i<N_Points;i++)
  {
    // reset matrix
    memset(GeoData, 0, 25*SizeOfDouble);

    Xi = xi[i];
    Eta = eta[i];

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

    GeoData[0][0] = xc1+xc3*Eta;
    GeoData[0][1] = yc1+yc3*Eta;

    GeoData[1][0] = xc2+xc3*Xi;
    GeoData[1][1] = yc2+yc3*Xi;

    GeoData[2][0] = 0;
    GeoData[2][1] = 0;
    GeoData[2][2] = GeoData[0][0]*GeoData[0][0];
    GeoData[2][3] = 2*GeoData[0][0]*GeoData[0][1];
    GeoData[2][4] = GeoData[0][1]*GeoData[0][1];
    
    GeoData[3][0] = xc3;
    GeoData[3][1] = yc3;
    GeoData[3][2] = GeoData[0][0]*GeoData[1][0];
    GeoData[3][3] = GeoData[0][1]*GeoData[1][0]+GeoData[0][0]*GeoData[1][1];
    GeoData[3][4] = GeoData[0][1]*GeoData[1][1];
    
    GeoData[4][0] = 0;
    GeoData[4][1] = 0;
    GeoData[4][2] = GeoData[1][0]*GeoData[1][0];
    GeoData[4][3] = 2*GeoData[1][0]*GeoData[1][1];
    GeoData[4][4] = GeoData[1][1]*GeoData[1][1];

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
void TQuadBilinear::GetOrigValues(int N_Sets, BaseFunct2D *BaseFuncts,
                               int N_Points, double *xi, double *eta,
                               QuadFormula2D QuadFormula,
                               bool *Needs2ndDer)
{
  int i,j,k,N_, start, end, BaseVectDim, M;
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
  double r10, r01, r20, r11, r02;
  double o10, o01, o20, o11, o02;
  double *aux, Xi, Eta, rec_detjk;
  double GeoData[5][5], GeoData1[5][5];
  double Eye[5][5];
  BaseFunct2D BaseFunct;
  int N_Functs;
  bool SecondDer;

  SecondDer = false;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    BaseVectDim = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();
    M=BaseVectDim*MaxN_QuadPoints_2D;
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
      origvaluesD00 = new double* [M];
      aux = new double [M*N_Functs];
      for(j=0;j<M;j++)
        origvaluesD00[j] = aux+j*N_Functs;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D00, origvaluesD00);
    }
  
    for(j=0;j<N_Points;j++)
    {
      refD00 = refvaluesD00[j];
      origD00 = origvaluesD00[j];

      if(BaseVectDim == 1)
        memcpy(origD00, refD00, N_Functs*SizeOfDouble);
      else
        this->PiolaMapOrigFromRefNotAffine(N_Functs, refD00, origD00, xi[j], 
                                           eta[j]);
    } // endfor j

    refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D10);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    if(origvaluesD10==NULL)
    {
      origvaluesD10 = new double* [M];
      aux = new double [M*N_Functs];
      for(j=0;j<M;j++)
        origvaluesD10[j] = aux+j*N_Functs;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D10, origvaluesD10);
    }
  
    refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D01);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);
    if(origvaluesD01==NULL)
    {
      origvaluesD01 = new double* [M];
      aux = new double [M*N_Functs];
      for(j=0;j<M;j++)
        origvaluesD01[j] = aux+j*N_Functs;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D01, origvaluesD01);
    }
  
    if(Needs2ndDer[i])
    {
      SecondDer = true;

      refvaluesD20=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D20);
      origvaluesD20=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D20);
      if(origvaluesD20==NULL)
      {
        origvaluesD20 = new double* [M];
        aux = new double [M*N_Functs];
        for(j=0;j<M;j++)
          origvaluesD20[j] = aux+j*N_Functs;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D20, origvaluesD20);
      }
    
      refvaluesD11=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D11);
      origvaluesD11=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D11);
      if(origvaluesD11==NULL)
      {
        origvaluesD11 = new double* [M];
        aux = new double [M*N_Functs];
        for(j=0;j<M;j++)
          origvaluesD11[j] = aux+j*N_Functs;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D11, origvaluesD11);
      }
    
      refvaluesD02=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D02);
      origvaluesD02=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D02);
      if(origvaluesD02==NULL)
      {
        origvaluesD02 = new double* [M];
        aux = new double [M*N_Functs];
        for(j=0;j<M;j++)
          origvaluesD02[j] = aux+j*N_Functs;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D02, origvaluesD02);
      }
    } // endfor Needs2ndDer[i]
  } // endfor i
  
  if(!SecondDer)
  {
    // no element needs second derivatives

    // D10 and D01
    for(i=0;i<N_Sets;i++)
    {
      BaseFunct=BaseFuncts[i];
      BaseVectDim = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();
      N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();
  
      refvaluesD10=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D10);
      origvaluesD10=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D10);
  
      refvaluesD01=TFEDatabase2D::GetRefElementValues
                        (BaseFunct, QuadFormula, D01);
      origvaluesD01=TFEDatabase2D::GetOrigElementValues
                        (BaseFunct, D01);
  
      for(j=0;j<N_Points;j++)
      {
        refD10 = refvaluesD10[j];
        origD10 = origvaluesD10[j];
    
        refD01 = refvaluesD01[j];
        origD01 = origvaluesD01[j];
    
        Xi = xi[j];
        Eta = eta[j];
        if(BaseVectDim == 1)
        {
          rec_detjk = 1/( (xc1+xc3*Eta)*(yc2+yc3*Xi)
                        -(xc2+xc3*Xi)*(yc1+yc3*Eta));
          for(k=0;k<N_Functs;k++)
          {
            origD10[k]=( (yc2+yc3*Xi) * refD10[k]
                        -(yc1+yc3*Eta) * refD01[k]) * rec_detjk;
            origD01[k]=(-(xc2+xc3*Xi) * refD10[k]
                        +(xc1+xc3*Eta) * refD01[k]) * rec_detjk;
          } // endfor k
        }
        else // BaseVectDim == 2
        {
          refD00 = refvaluesD00[j];
          this->PiolaMapOrigFromRefNotAffine(N_Functs, refD00, refD10, refD01, 
                                             origD10, origD01, Xi, Eta);
        }
      } // endfor j
        
    } // endfor i
  }
  else
  {
    // at least one element needs second derivatives

    // find transformation matrix for second derivatives
    // do this only once since matrix is constant for an affine mapping
  
    for(j=0;j<N_Points;j++)
    {
      Xi = xi[j];
      Eta = eta[j];
    
      // reset matrices
      memset(GeoData, 0, 25*SizeOfDouble);
      memset(Eye, 0, 25*SizeOfDouble);
      Eye[0][0] = 1;
      Eye[1][1] = 1;
      Eye[2][2] = 1;
      Eye[3][3] = 1;
      Eye[4][4] = 1;
      
      GeoData[0][0] = xc1+xc3*Eta;
      GeoData[0][1] = yc1+yc3*Eta;
  
      GeoData[1][0] = xc2+xc3*Xi;
      GeoData[1][1] = yc2+yc3*Xi;
  
      GeoData[2][0] = 0;
      GeoData[2][1] = 0;
      GeoData[2][2] = GeoData[0][0]*GeoData[0][0];
      GeoData[2][3] = 2*GeoData[0][0]*GeoData[0][1];
      GeoData[2][4] = GeoData[0][1]*GeoData[0][1];
      
      GeoData[3][0] = xc3;
      GeoData[3][1] = yc3;
      GeoData[3][2] = GeoData[0][0]*GeoData[1][0];
      GeoData[3][3] = GeoData[0][1]*GeoData[1][0]+GeoData[0][0]*GeoData[1][1];
      GeoData[3][4] = GeoData[0][1]*GeoData[1][1];
    
      GeoData[4][0] = 0;
      GeoData[4][1] = 0;
      GeoData[4][2] = GeoData[1][0]*GeoData[1][0];
      GeoData[4][3] = 2*GeoData[1][0]*GeoData[1][1];
      GeoData[4][4] = GeoData[1][1]*GeoData[1][1];

      // subroutine for solving a multiple systems of linear equations
      // void SolveMultipleSystems(double *a, double *b, int N_Eqn,
      //                        int LDA, int LDB, int N_Rhs);
      SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 5,
                           5, 5, 5);
      //SolveMultipleSystems((double *)GeoData, (double *)Eye, 5,
      //                    5, 5, 5);
      for(i=0;i<N_Sets;i++)
      {
        BaseFunct=BaseFuncts[i];
        BaseVectDim = 
          TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim();
        if(BaseVectDim > 1)
        {
          ErrMsg("Piola transform for bilinear reference transformation for " <<
                 "second derivatives not yet implemented");
          exit(1);
        }

        N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();
    
        refvaluesD10=TFEDatabase2D::GetRefElementValues
                                (BaseFunct, QuadFormula, D10);
        origvaluesD10=TFEDatabase2D::GetOrigElementValues
                                (BaseFunct, D10);
    
        refvaluesD01=TFEDatabase2D::GetRefElementValues
                                (BaseFunct, QuadFormula, D01);
        origvaluesD01=TFEDatabase2D::GetOrigElementValues
                                (BaseFunct, D01);
    
        refD10 = refvaluesD10[j];
        refD01 = refvaluesD01[j];

        origD10 = origvaluesD10[j];
        origD01 = origvaluesD01[j];

        if(Needs2ndDer[i])
        {
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
  
          refD20 = refvaluesD20[j];
          refD11 = refvaluesD11[j];
          refD02 = refvaluesD02[j];

          origD20 = origvaluesD20[j];
          origD11 = origvaluesD11[j];
          origD02 = origvaluesD02[j];
      
          for(k=0;k<N_Functs;k++)
          {
            r10 = refD10[k];
            r01 = refD01[k];
            r20 = refD20[k];
            r11 = refD11[k];
            r02 = refD02[k];

            o10 = Eye[0][0]*r10+Eye[0][1]*r01+Eye[0][2]*r20
                 +Eye[0][3]*r11+Eye[0][4]*r02;
            o01 = Eye[1][0]*r10+Eye[1][1]*r01+Eye[1][2]*r20
                 +Eye[1][3]*r11+Eye[1][4]*r02;
            o20 = Eye[2][0]*r10+Eye[2][1]*r01+Eye[2][2]*r20
                 +Eye[2][3]*r11+Eye[2][4]*r02;
            o11 = Eye[3][0]*r10+Eye[3][1]*r01+Eye[3][2]*r20
                 +Eye[3][3]*r11+Eye[3][4]*r02;
            o02 = Eye[4][0]*r10+Eye[4][1]*r01+Eye[4][2]*r20
                 +Eye[4][3]*r11+Eye[4][4]*r02;
      
            origD10[k] = o10;
            origD01[k] = o01;
            origD20[k] = o20;
            origD11[k] = o11;
            origD02[k] = o02;
          } // endfor k
        }
        else
        {
          for(k=0;k<N_Functs;k++)
          {
            r10 = refD10[k];
            r01 = refD01[k];

            o10 = Eye[0][0]*r10+Eye[0][1]*r01;
            o01 = Eye[1][0]*r10+Eye[1][1]*r01;

            origD10[k] = o10;
            origD01[k] = o01;
          } // endfor k
        } // end SecondDer
      } // endfor i
    } // endfor j
  } // end SecondDer
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadBilinear::GetOrigValues(double xi, double eta, 
                int N_BaseFunct,
                double *uref, double *uxiref, double *uetaref,
                double *uorig, double *uxorig, double *uyorig, int _BaseVectDim)
{
  int i;
  double rec_detjk;

  if(_BaseVectDim==1) // standard case
  {
    // D00
    for(i=0;i<N_BaseFunct;i++)
    {
      uorig[i] = uref[i];
    } // endfor i
    
    // D10 and D01
    if (uxiref!=NULL && uetaref!=NULL && uxorig!=NULL && uyorig!=NULL)
    {
      rec_detjk = 1/( (xc1+xc3*eta)*(yc2+yc3*xi) - (xc2+xc3*xi)*(yc1+yc3*eta) );
      for(i=0;i<N_BaseFunct;i++)
      {
        uxorig[i]= ( (yc2+yc3*xi) * uxiref[i] - (yc1+yc3*eta) * uetaref[i] )
                  *rec_detjk;
        uyorig[i]= (-(xc2+xc3*xi) * uxiref[i] + (xc1+xc3*eta) * uetaref[i] )
                  *rec_detjk;
      } // endfor i
    }
  }
  else
  {
    // Piola transformation
    // D00
    this->PiolaMapOrigFromRefNotAffine(N_BaseFunct, uref, uorig, xi, eta);
    // D10, D01
    if (uxiref!=NULL && uetaref!=NULL && uxorig!=NULL && uyorig!=NULL)
    {
      this->PiolaMapOrigFromRefNotAffine(N_BaseFunct, uref, uxiref, uetaref,
                                         uxorig, uyorig, xi, eta);
    }
  }
}

void TQuadBilinear::GetOrigValues(int joint, double zeta, int N_BaseFunct,
                                  double *uref, double *uxiref, double *uetaref,
                                  double *uorig, double *uxorig, double *uyorig,
                                  int _BaseVectDim)
{
  double xi, eta;
  switch(joint)
  {
    case 0:
      xi  = zeta; eta = -1;
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


void TQuadBilinear::SetCell(TBaseCell *cell)
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

  xc0=( x0 + x1 + x2 + x3) * 0.25;
  xc1=(-x0 + x1 + x2 - x3) * 0.25;
  xc2=(-x0 - x1 + x2 + x3) * 0.25;
  xc3=( x0 - x1 + x2 - x3) * 0.25;

  yc0=( y0 + y1 + y2 + y3) * 0.25;
  yc1=(-y0 + y1 + y2 - y3) * 0.25;
  yc2=(-y0 - y1 + y2 + y3) * 0.25;
  yc3=( y0 - y1 + y2 - y3) * 0.25;
}

/** return outer normal vector */
void TQuadBilinear::GetOuterNormal(int j, double zeta,
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

/** return tangent */
void TQuadBilinear::GetTangent(int j, double zeta,
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
double TQuadBilinear::GetVolume()
{
  double a1, a2, b1, b2, c1, c2, d1, d2;
  double locvol;
#ifdef __3D__
  double e;
#endif

#ifdef __3D__
  Cell->GetVertex(0)->GetCoords(a1, a2, e);
  Cell->GetVertex(1)->GetCoords(b1, b2, e);
  Cell->GetVertex(2)->GetCoords(c1, c2, e);
  Cell->GetVertex(3)->GetCoords(d1, d2, e);
#else
  Cell->GetVertex(0)->GetCoords(a1, a2);
  Cell->GetVertex(1)->GetCoords(b1, b2);
  Cell->GetVertex(2)->GetCoords(c1, c2);
  Cell->GetVertex(3)->GetCoords(d1, d2);
#endif

  locvol = 0.5*((b1-a1)*(c2-a2)-(b2-a2)*(c1-a1)
               +(c1-a1)*(d2-a2)-(c2-a2)*(d1-a1));

  return locvol;
}

void TQuadBilinear::PiolaMapOrigFromRefNotAffine(int N_Functs, double *refD00, 
                                                 double *origD00, double xi, 
                                                 double eta)
{
  double rec_detjk = 1/( (xc1+xc3*eta)*(yc2+yc3*xi)-(xc2+xc3*xi)*(yc1+yc3*eta));
  double a11 = (xc1+xc3*eta)*rec_detjk;
  double a12 = (xc2+xc3*xi )*rec_detjk;
  double a21 = (yc1+yc3*eta)*rec_detjk;
  double a22 = (yc2+yc3*xi )*rec_detjk;
  
  for(int k = 0; k < N_Functs; k++)
  {
    // Piola transformation
    // phi = 1/|J| DF phi_hat
    // def.gradient/detjk
    origD00[k] = a11*refD00[k]+a12*refD00[N_Functs+k];
    origD00[N_Functs+k] = a21*refD00[k]+a22*refD00[N_Functs+k];
  }
}

void TQuadBilinear::PiolaMapOrigFromRefNotAffine(int N_Functs, double *refD00, 
                                                 double *refD10, 
                                                 double *refD01, 
                                                 double *origD10, 
                                                 double *origD01, double xi,
                                                 double eta)
{
  // Piola transformation
  // phi = 1/|J| DF phi_hat
  // 
  double rec_detjk = 1/((xc1+xc3*eta)*(yc2+yc3*xi ) 
                       -(xc2+xc3*xi )*(yc1+yc3*eta));
  // this is \hat{D}F
  double a11 = (xc1+xc3*eta);
  double a12 = (xc2+xc3* xi);
  double a21 = (yc1+yc3*eta);
  double a22 = (yc2+yc3* xi);
  // this is DF^{-1}
  double z11 =  a22*rec_detjk;
  double z12 = -a12*rec_detjk;
  double z21 = -a21*rec_detjk;
  double z22 =  a11*rec_detjk;
  
  for(int k = 0; k < N_Functs; k++)
  {
    // the derivative has three parts (since the Piola transform is a 
    // product of three factors, namely 1/|J|, DF, and phi_hat) 
    // according to the product rule
  
    // refPiola_(xy)_D.. is the derivative with respect to the 
    // reference variables xi and eta. In the end we apply the chain 
    // rule to get the derivative with respect to x and y.
    double refPiola_x_D10,refPiola_x_D01,refPiola_y_D10,refPiola_y_D01;
    
    // first part (differentiate 1/|J|)
    // \hat{D}F \hat{v}
    double DFv1 = a11*refD00[k] + a12*refD00[k+N_Functs];
    double DFv2 = a21*refD00[k] + a22*refD00[k+N_Functs];
    // D(det(DF)) (row vector)
    double b1 = xc1*yc3 - yc1*xc3;
    double b2 = yc2*xc3 - xc2*yc3;
    b1 *= -SIGN(rec_detjk)*rec_detjk*rec_detjk;
    b2 *= -SIGN(rec_detjk)*rec_detjk*rec_detjk;
    // Dfv (-sign(J)/(J^2))DJ
    refPiola_x_D10 = DFv1*b1;
    refPiola_x_D01 = DFv1*b2;
    refPiola_y_D10 = DFv2*b1;
    refPiola_y_D01 = DFv2*b2;
    
    // second part (differentiate DF)
    // D^2F v /det(DF)
    refPiola_x_D10 += refD00[k+N_Functs]*xc3*rec_detjk;
    refPiola_x_D01 += refD00[k]         *xc3*rec_detjk;
    refPiola_y_D10 += refD00[k+N_Functs]*yc3*rec_detjk;
    refPiola_y_D01 += refD00[k]         *yc3*rec_detjk;
    
    // third part (differentiate phi_hat), similar to the affine linear 
    // case
    // \hat{D}F \hat{D}\hat{v} /|J|
    refPiola_x_D10 += (a11*refD10[k]+a12*refD10[N_Functs+k])*rec_detjk;
    refPiola_x_D01 += (a11*refD01[k]+a12*refD01[N_Functs+k])*rec_detjk;
    refPiola_y_D10 += (a21*refD10[k]+a22*refD10[N_Functs+k])*rec_detjk;
    refPiola_y_D01 += (a21*refD01[k]+a22*refD01[N_Functs+k])*rec_detjk;
    
    // apply DF^{-1} from right (chain rule)
    origD10[k] = z11*refPiola_x_D10 + z21*refPiola_x_D01;
    origD01[k] = z12*refPiola_x_D10 + z22*refPiola_x_D01;
    origD10[k+N_Functs] = z11*refPiola_y_D10 + z21*refPiola_y_D01;
    origD01[k+N_Functs] = z12*refPiola_y_D10 + z22*refPiola_y_D01;
    
  } // endfor k
}
