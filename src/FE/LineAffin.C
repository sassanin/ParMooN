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
// @(#)LineAffin.C
//
// Class:      LineAffin
//
// Purpose:    reference transformations for Line
//
// Author:     Sashikumaar Ganesan
//
// History:    17.05.2007 start implementation
// 
// =======================================================================

#include <LineAffin.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>

/** constuctor */
TLineAffin::TLineAffin()
{
}

/** transfer from reference element to original element */
void TLineAffin::GetOrigFromRef(double xi, double &X
#ifdef __2D__
, double &Y
#endif
                                )
{
    X = xc0 + xc1*xi;
#ifdef __2D__
    Y = 0.;
#endif
}

/** transfer a set of points from reference to original element */
void TLineAffin::GetOrigFromRef(int N_Points, double *xi, double *X, 
                                double *Y,  double *absdetjk)
{
  int i;
  double Xi;
  double absdet = fabs(detjk);

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    X[i] = xc0 + xc1*Xi;
    absdetjk[i] = absdet;
    
    #ifdef __2D__
    Y[i] = 0;
    #endif
  } // endfor i
}

/** transfer a set of points from reference to original element */
void TLineAffin::GetOrigFromRef(int N_Points, double *xi, double *X, double *absdetjk)
{
  int i;
  double Xi;
  double absdet = fabs(detjk);

  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    X[i] = xc0 +  xc1*Xi;
    absdetjk[i] = absdet;
  } // endfor i

}


void TLineAffin::GetOrigValues(int N_Sets, BaseFunct1D *BaseFuncts, int N_Points, double *zeta,
                               QuadFormula1D QuadFormula, bool *Needs2ndDer)
{
  int i, j, k, l, N_Functs;

  double **refvaluesD0, **origvaluesD0;
  double *refD0, *origD0;
  double **refvaluesD1, **origvaluesD1;
  double *refD1, *origD1;
  double **refvaluesD2, **origvaluesD2;
  double *refD2, *origD2;
  double *aux;

  bool SecondDer=FALSE;

  BaseFunct1D BaseFunct;


  for(i=0;i<N_Sets;i++)
   {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct1D(BaseFunct)->GetDimension();

    refvaluesD0 = TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D0);
    if(refvaluesD0==NULL)
     {
      TFEDatabase2D::GetBaseFunct1D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD0=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D0);
     }

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct, D0);
    if(origvaluesD0==NULL)
     {
      origvaluesD0 = new double* [MaxN_QuadPoints_1D];
      aux = new double [MaxN_QuadPoints_1D*MaxN_BaseFunctions1D];

      for(j=0;j<MaxN_QuadPoints_1D;j++)
        origvaluesD0[j] = aux+j*MaxN_BaseFunctions1D;

      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D0, origvaluesD0);
     }

    for(j=0;j<N_Points;j++)
     {
      refD0 = refvaluesD0[j];
      origD0 = origvaluesD0[j];

      memcpy(origD0, refD0, N_Functs*SizeOfDouble);
     } // endfor j

    refvaluesD1=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D1);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct, D1);

    if(origvaluesD1==NULL)
     {
      origvaluesD1 = new double* [MaxN_QuadPoints_1D];
      aux = new double [MaxN_QuadPoints_1D*MaxN_BaseFunctions1D];

      for(j=0;j<MaxN_QuadPoints_1D;j++)
        origvaluesD1[j] = aux+j*MaxN_BaseFunctions1D;

      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D1, origvaluesD1);
     }

    if(Needs2ndDer[i])
     {
      SecondDer = TRUE;
      
      refvaluesD2=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D2);
      origvaluesD2=TFEDatabase2D::GetOrigElementValues(BaseFunct, D2);    
    
     if(origvaluesD2==NULL)
     {
      origvaluesD2 = new double* [MaxN_QuadPoints_1D];
      aux = new double [MaxN_QuadPoints_1D*MaxN_BaseFunctions1D];

      for(j=0;j<MaxN_QuadPoints_1D;j++)
        origvaluesD2[j] = aux+j*MaxN_BaseFunctions1D;

      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D2, origvaluesD2);
     }     
    } //     if(Needs2ndDer[i])
   } //  for(i=0;i<N_Sets;i++)

  // D1
  for(i=0;i<N_Sets;i++)
   {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct1D(BaseFunct)->GetDimension();

    refvaluesD1=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D1);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct, D1);

    for(j=0;j<N_Points;j++)
     {
      refD1 = refvaluesD1[j];
      origD1 = origvaluesD1[j];
      for(k=0;k<N_Functs;k++)
       {
        origD1[k]=refD1[k]*rec_detjk;
       } // endfor k
     }// endfor j
   } //  for(i=0;i<N_Sets;i++)
  
  // leave if no second derivatives are needed
  if(!SecondDer) return;
  
   //D2
    for(i=0;i<N_Sets;i++)
     {
      BaseFunct=BaseFuncts[i];
      N_Functs = TFEDatabase2D::GetBaseFunct1D(BaseFunct)->GetDimension(); 
      
      refvaluesD2=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D2);
      origvaluesD2=TFEDatabase2D::GetOrigElementValues(BaseFunct, D2);
      
      for(j=0;j<N_Points;j++)
       {
        refD2 = refvaluesD2[j];
        origD2 = origvaluesD2[j];
        for(k=0;k<N_Functs;k++)
          origD2[k]=refD2[k]*rec_detjk*rec_detjk;
       }// endfor j         
     } // for(i=0;i<N   
     
}

void TLineAffin::SetCell(TBaseCell *cell)
{
  int i;
// #ifdef __2D__
  double z0, z1;
// #endif
  x0 = 0.;   x1 = 0.; 
  y0 = 0.;   y1 = 0.; 
  Cell = cell;

#ifdef __3D__
  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
#else
  Cell->GetVertex(0)->GetCoords(x0, y0);
  Cell->GetVertex(1)->GetCoords(x1, y1);
//   y0 = 0.;   y1 = 0.; 
#endif

// cout<< " y0 " << y0<< " y1 " << y1<<endl;
// cout<< " x0 " << x0<< " x1 " << x1<<endl;

  xc0 = (x0 + x1) * 0.5;
  xc1 = (x1 - x0) * 0.5;

  yc0 = 0;
  yc1 = 0;

  detjk=xc1;

  rec_detjk = 1./detjk;
}

 
