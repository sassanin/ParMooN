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
// @(#)QuadFormula2D.C        1.3 05/04/99
//
// Class:      TQuadFormula2D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 2D integral
// Author:     Gunar Matthies
//
// History:    29.08.1997 start implementation
// 
// =======================================================================

#include <QuadFormula2D.h>
#include <FEDatabase2D.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>

TQuadFormula2D::TQuadFormula2D() : TQuadFormula()
{
  Xi=NULL;
  Eta=NULL;
}

TQuadFormula2D::TQuadFormula2D(int n_points, double* weights, 
                               double* xi, double* eta, int acc)
{
  InitObject(n_points, weights, xi, eta, acc);
}

void TQuadFormula2D::InitObject(int n, double* w, double* xi, 
                                 double* eta, int acc)
{
  int i;

  N_QuadPoints=n;

  Weights=new double[n];
  Xi=new double[n];
  Eta=new double[n];

  for(i=0;i<n;i++)
  {
    Weights[i]=w[i];
    Xi[i]=xi[i];
    Eta[i]=eta[i];
  }

  Accuracy = acc;
}

double *TQuadFormula2D::GetCoords(int i)
{
  double *ret=NULL;
  if(i==0) 
    ret=Xi;
  else
    if(i==1)
      ret=Eta;
    
  return ret;
}

void TQuadFormula2D::
  GetFormulaData(int &n_points, double* &weights, double* &xi, double* &eta)
{
  n_points=N_QuadPoints;
  weights=Weights;
  xi=Xi;
  eta=Eta;
}

// ########################################################################
// modify later
// ########################################################################
void TQuadFormula2D::FindQuadFormula2D(FE2D *UsedElements,
        QuadFormula1D &qf1, QuadFormula2D &qf2)
{
  switch(TFEDatabase2D::GetFE2D(UsedElements[0])->GetRefTransID())
  {
    case TriaAffin:
    case TriaIsoparametric:
      qf1 = Gauss3Line;
      qf2 = Gauss3Tria;
    break;

    case QuadAffin:
    case QuadBilinear:
    case QuadIsoparametric:  
      qf1 = Gauss3Line;
      qf2 = Gauss5Quad;
    break;
   default: cerr << "unknown RefTransID()" << endl;
          exit(-1);
     break; 
  } // endswitch
    
}

void TQuadFormula2D::FindQF_2D(FE2D CurrentElement,
        QuadFormula1D &qf1, QuadFormula2D &qf2)
{
  switch(TFEDatabase2D::GetFE2D(CurrentElement)->GetRefTransID())
  {
    case TriaAffin:
    case TriaIsoparametric:      
      qf1 = Gauss3Line;
      qf2 = Gauss3Tria;
    break;

    case QuadAffin:
    case QuadBilinear:
    case QuadIsoparametric:        
      qf1 = Gauss3Line;
      qf2 = Gauss5Quad;
    break;
   default: cerr << "unknown RefTransID()" << endl;
          exit(-1);
     break;     
  } // endswitch
}

void TQuadFormula2D::FindLocalQuadFormula2D
        (int N_LocalUsedElements, FE2D *LocalUsedElements,
         QuadFormula1D &qf1, QuadFormula2D &qf2)
{
  int i,j, MaxPolynomialDegree, PolynomialDegree;
  BF2DRefElements RefElement;
  int *PolynomialDegreeFromFE2D;
  RefTrans2D RefTrans, *RefTransArray, CurrentRefTrans;
  TRefTrans2D *rt;
  BaseFunct2D BaseFuncts[N_FEs2D];

  // find adequate quadrature formula for all elements
  // and find needed reference transformation
  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();
  RefTrans = RefTransArray[LocalUsedElements[0]];
  PolynomialDegreeFromFE2D = 
        TFEDatabase2D::GetPolynomialDegreeFromFE2D();
  MaxPolynomialDegree = 0;
  RefElement = 
        TFEDatabase2D::GetRefElementFromFE2D(LocalUsedElements[0]);
  for(i=0;i<N_LocalUsedElements;i++)
  {
    BaseFuncts[i] = 
        TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(LocalUsedElements[i]);
    PolynomialDegree = PolynomialDegreeFromFE2D[LocalUsedElements[i]];
    if(PolynomialDegree > MaxPolynomialDegree) 
      MaxPolynomialDegree = PolynomialDegree;

    CurrentRefTrans = RefTransArray[LocalUsedElements[i]];
    if(CurrentRefTrans > RefTrans)
      RefTrans = CurrentRefTrans;
  }

  switch(RefElement)
  {
    case BFUnitSquare:
      qf2 = TFEDatabase2D::GetQFQuadFromDegree
                                        (2*MaxPolynomialDegree);
    break;

    case BFUnitTriangle:
      qf2 = TFEDatabase2D::GetQFTriaFromDegree
                                        (2*MaxPolynomialDegree);
    break;
  } // endswitch

  qf1 = TFEDatabase2D::GetQFLineFromDegree
                                        (2*MaxPolynomialDegree);
}

std::ostream& operator << (std::ostream &s, TQuadFormula2D *qf)
{
  int i,N_;

  s << endl;
  s << "instance of TQuadFormula2D" << endl;
  N_=qf->N_QuadPoints;
  s << "N_QuadPoints: " << N_ << endl;
  s << setw(3) << "No";
  s << setw(11) << "weights";
  s << setw(11) << "xi";
  s << setw(11) << "eta";
  s << endl;
  s.setf(std::ios::fixed);
  for(i=0;i<N_;i++)
  {
    s << setw(3) << i;
    s << setprecision(6) << setw(11) << qf->Weights[i];
    s << setprecision(6) << setw(11) << qf->Xi[i];
    s << setprecision(6) << setw(11) << qf->Eta[i];
    s << endl;
  }

  return s << endl;
}
