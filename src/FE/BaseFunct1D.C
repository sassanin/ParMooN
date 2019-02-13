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
// @(#)BaseFunct1D.C        1.1 10/30/98
//
// Class:      TBaseFunct1D
//
// Purpose:    represents the set of base functions for a finite element
//             in one dimensions
//
// Author:     Gunar Matthies
//
// History:    08.10.97 start implementation
// 
// =======================================================================

#include <Constants.h>
#include <BaseFunct1D.h>
#include <FEDatabase2D.h>
#include <stdlib.h>

/** constructor, fill in all information */
TBaseFunct1D::TBaseFunct1D(int dimension, BaseFunct1D basefunct,
                           DoubleFunct1D* functions, 
                           DoubleFunct1D* derivativesxi,
                           DoubleFunct1D* derivativesxixi)
{
  Dimension=dimension;

  BaseFunct = basefunct;

  Functions[D0]=functions;
  Functions[D1]=derivativesxi;
  Functions[D2]=derivativesxixi;

  changable = false;
}

/** constructor, fill in all information */
TBaseFunct1D::TBaseFunct1D(int dimension, BaseFunct1D basefunct,
                           DoubleFunct1D* functions, 
                           DoubleFunct1D* derivativesxi,
                           DoubleFunct1D* derivativesxixi,
                           int polynomialdegree,
                           int accuracy)
{
  Dimension=dimension;

  BaseFunct = basefunct;

  Functions[D0]=functions;
  Functions[D1]=derivativesxi;
  Functions[D2]=derivativesxixi;

  changable = false;

  PolynomialDegree = polynomialdegree;
  Accuracy = accuracy;
}


/** constructor without filling data structure */
TBaseFunct1D::TBaseFunct1D(int dimension)
{
  Dimension = dimension;

  changable = false;
}

/** return the values for derivative MultiIndex at all
    quadrature points */
void TBaseFunct1D::GetDerivatives(MultiIndex1D MultiIndex, 
                        TQuadFormula1D *formula, double **values)
{
  int i, N_;
  double *w, *xi;

  formula->GetFormulaData(N_, w, xi);

  for(i=0;i<N_;i++)
    GetDerivatives(MultiIndex, xi[i], values[i]);

}

/** set function for derivative MultiIndex */
void TBaseFunct1D::SetFunction(MultiIndex1D MultiIndex, 
                               DoubleFunct1D* function)
{
  if(changable)
    Functions[MultiIndex] = function;
}

/** make data on reference element */
void TBaseFunct1D::MakeRefElementData(QuadFormula1D QuadFormula)
{
  int j;
  double **Values, *AllValues;
  TQuadFormula1D *qf;

  qf = TFEDatabase2D::GetQuadFormula1D(QuadFormula);

  // D0
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D0);
  if( Values==NULL)
  {
    Values = new double* [MaxN_QuadPoints_1D];
    AllValues = new double [MaxN_QuadPoints_1D*MaxN_BaseFunctions1D];
    for(j=0;j<MaxN_QuadPoints_1D;j++)
      Values[j] = AllValues+j*MaxN_BaseFunctions1D;
    GetDerivatives(D0, qf, Values);

    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D0, Values);
  }

  // D1
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D1);
  if( Values==NULL)
  {
    Values = new double* [MaxN_QuadPoints_1D];
    AllValues = new double [MaxN_QuadPoints_1D*MaxN_BaseFunctions1D];
    for(j=0;j<MaxN_QuadPoints_1D;j++)
      Values[j] = AllValues+j*MaxN_BaseFunctions1D;
    GetDerivatives(D1, qf, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D1, Values);
  }

  // D2
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D2);
  if( Values==NULL)
  {
    Values = new double* [MaxN_QuadPoints_1D];
    AllValues = new double [MaxN_QuadPoints_1D*MaxN_BaseFunctions1D];
    for(j=0;j<MaxN_QuadPoints_1D;j++)
      Values[j] = AllValues+j*MaxN_BaseFunctions1D;
    GetDerivatives(D2, qf, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D2, Values);
  }

}
