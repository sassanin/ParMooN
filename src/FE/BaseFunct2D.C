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
// @(#)BaseFunct2D.C        1.8 02/08/00
//
// Class:      TBaseFunct2D
//
// Purpose:    represents the set of base functions for a finite element
//             in two dimensions
//
// Author:     Gunar Matthies
//
// History:    26.11.97 start implementation
// 
//       :     10.06.2010 methods for vector basis functions (Sashikumaar Ganesan)
// =======================================================================

#include <Constants.h>
#include <BaseFunct2D.h>
#include <FEDatabase2D.h>
#include <Database.h>
#include <JointEqN.h>
#include <stdlib.h>
/** constructor, fill in all information */
TBaseFunct2D::TBaseFunct2D(int dimension, BaseFunct2D basefunct,
                           BF2DRefElements refelement, 
                           DoubleFunct2D* functions, 
                           DoubleFunct2D* derivativesxi,
                           DoubleFunct2D* derivativeseta,
                           DoubleFunct2D* derivativesxixi,
                           DoubleFunct2D* derivativesxieta,
                           DoubleFunct2D* derivativesetaeta,
                           int polynomialdegree,
                           int accuracy,
                           int n_bf2change,
                           int **bf2change)
{
  Dimension=dimension;

  RefElement = refelement;

  BaseFunct = basefunct;

  Functions[D00]=functions;
  Functions[D10]=derivativesxi;
  Functions[D01]=derivativeseta;
  Functions[D20]=derivativesxixi;
  Functions[D11]=derivativesxieta;
  Functions[D02]=derivativesetaeta;

  changable = false;

  PolynomialDegree = polynomialdegree;
  Accuracy = accuracy;

  N_BF2Change = n_bf2change;
  BF2Change = bf2change;
  
  //default is a scalar basis function
  BaseVectDim = 1;
  SpaceDeptBasis = false;
  
}

/** constructor, fill in all information with scalar basis function dimension*/
TBaseFunct2D::TBaseFunct2D(int dimension, BaseFunct2D basefunct,
                           BF2DRefElements refelement, 
                           DoubleFunct2D* functions, 
                           DoubleFunct2D* derivativesxi,
                           DoubleFunct2D* derivativeseta,
                           DoubleFunct2D* derivativesxixi,
                           DoubleFunct2D* derivativesxieta,
                           DoubleFunct2D* derivativesetaeta,
                           int polynomialdegree,
                           int accuracy,
                           int n_bf2change,
                           int **bf2change, int baseVectDim)
{
  Dimension=dimension;

  RefElement = refelement;

  BaseFunct = basefunct;

  Functions[D00]=functions;
  Functions[D10]=derivativesxi;
  Functions[D01]=derivativeseta;
  Functions[D20]=derivativesxixi;
  Functions[D11]=derivativesxieta;
  Functions[D02]=derivativesetaeta;

  changable = false;

  PolynomialDegree = polynomialdegree;
  Accuracy = accuracy;

  N_BF2Change = n_bf2change;
  BF2Change = bf2change;
  
  BaseVectDim = baseVectDim;
  SpaceDeptBasis = false;
}

/** constructor, fill in all information with space dept. basis functions*/
TBaseFunct2D::TBaseFunct2D(int dimension, BaseFunct2D basefunct,
                 BF2DRefElements refelement,
                 DoubleFunct2D* functions, 
                 DoubleFunct2D* derivativesxi,
                 DoubleFunct2D* derivativeseta,
                 DoubleFunct2D* derivativesxixi,
                 DoubleFunct2D* derivativesxieta,
                 DoubleFunct2D* derivativesetaeta,
                 int polynomialdegree,
                 int accuracy,
                 int n_bf2change,
                 int **bf2change,
                 bool spaceDeptBasis
                )
{
  Dimension=dimension;

  RefElement = refelement;

  BaseFunct = basefunct;

  Functions[D00]=functions;
  Functions[D10]=derivativesxi;
  Functions[D01]=derivativeseta;
  Functions[D20]=derivativesxixi;
  Functions[D11]=derivativesxieta;
  Functions[D02]=derivativesetaeta;

  changable = false;

  PolynomialDegree = polynomialdegree;
  Accuracy = accuracy;

  N_BF2Change = n_bf2change;
  BF2Change = bf2change;
  
  //default is a scalar basis function
  BaseVectDim = 1;
  SpaceDeptBasis = spaceDeptBasis;  
}


/** constructor without filling data structure */
TBaseFunct2D::TBaseFunct2D(int dimension)
{
  Dimension = dimension;

  changable = true;
  SpaceDeptBasis = false;
}

/** return the values for derivative MultiIndex at all
    quadrature points */
void TBaseFunct2D::GetDerivatives(MultiIndex2D MultiIndex, 
                        TQuadFormula2D *formula, double **values)
{
  int i, N_;
  double *w, *xi, *eta;

  formula->GetFormulaData(N_, w, xi, eta);

  for(i=0;i<N_;i++)
    GetDerivatives(MultiIndex, xi[i], eta[i], values[i]);
}

/** return values on joint */
void TBaseFunct2D::GetValues(int N_Points, double *zeta, 
                             int joint, double **Values)
{
  int i;
  double xi[N_Points], eta[N_Points];

  switch(RefElement)
  {
    case BFUnitTriangle:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 0.5*(zeta[i]+1);
            eta[i] = 0;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            eta[i] = 0.5*(zeta[i]+1);
            xi[i] = 1-eta[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            eta[i] = 0.5*(-zeta[i]+1);
            xi[i] = 0;
          }
          break;
      } // endswitch
      break;
    case BFUnitSquare:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = zeta[i];
            eta[i] = -1;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 1;
            eta[i] = zeta[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -zeta[i];
            eta[i] = 1;
          }
          break;
        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -1;
            eta[i] = -zeta[i];
          }
          break;
      } // endswitch
      break;
  } // endswitch

  for(i=0;i<N_Points;i++)
  {
    GetDerivatives(D00, xi[i], eta[i], Values[i]);
  }
}

/** return derivatives on joint */
void TBaseFunct2D::GetDerivatives(MultiIndex2D MultiIndex, int N_Points,
                                  double *zeta, int joint, double **Values)
{
  int i;
  double xi[N_Points], eta[N_Points];

  switch(RefElement)
  {
    case BFUnitTriangle:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 0.5*(zeta[i]+1);
            eta[i] = 0;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            eta[i] = 0.5*(zeta[i]+1);
            xi[i] = 1-eta[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            eta[i] = 0.5*(-zeta[i]+1);
            xi[i] = 0;
          }
          break;
      } // endswitch
      break;
    case BFUnitSquare:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = zeta[i];
            eta[i] = -1;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 1;
            eta[i] = zeta[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -zeta[i];
            eta[i] = 1;
          }
          break;
        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -1;
            eta[i] = -zeta[i];
          }
          break;
      } // endswitch
      break;
  } // endswitch

  for(i=0;i<N_Points;i++)
  {
    GetDerivatives(MultiIndex, xi[i], eta[i], Values[i]);
  }
}

/** return values of derivative index on joint */
void TBaseFunct2D::GetValues(int N_Points, double *zeta, 
                             int joint, MultiIndex2D index, 
                             double **Values)
{
  int i;
  double xi[N_Points], eta[N_Points];

  switch(RefElement)
  {
    case BFUnitTriangle:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 0.5*(zeta[i]+1);
            eta[i] = 0;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            eta[i] = 0.5*(zeta[i]+1);
            xi[i] = 1-eta[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            eta[i] = 0.5*(-zeta[i]+1);
            xi[i] = 0;
          }
          break;
      } // endswitch
      break;
    case BFUnitSquare:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = zeta[i];
            eta[i] = -1;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 1;
            eta[i] = zeta[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -zeta[i];
            eta[i] = 1;
          }
          break;
        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -1;
            eta[i] = -zeta[i];
          }
          break;
      } // endswitch
      break;
  } // endswitch

  for(i=0;i<N_Points;i++)
  {
    GetDerivatives(index, xi[i], eta[i], Values[i]);
  }
}

/** set function for derivative MultiIndex */
void TBaseFunct2D::SetFunction(MultiIndex2D MultiIndex, 
                               DoubleFunct2D* function)
{
  if(changable)
    Functions[MultiIndex] = function;
}

/** make data on reference element */
/** added methods for vector basis function */
void TBaseFunct2D::MakeRefElementData(QuadFormula1D LineQuadFormula)
{
  int i, j, N_Joints;
  double **Values, *AllValues;
  TQuadFormula1D *qf1;
  double *w, *zeta;
  int N_Points;

  qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  qf1->GetFormulaData(N_Points, w, zeta);

  switch(RefElement)
  {
    case BFUnitTriangle:
      N_Joints = 3;
      break;
    case BFUnitSquare:
      N_Joints = 4;
      break;
    default:
      N_Joints = 0;
  }

  // joint values
  for(i=0;i<N_Joints;i++)
  {
    Values=TFEDatabase2D::GetJointValues2D(BaseFunct, LineQuadFormula, i);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension*BaseVectDim];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension*BaseVectDim;
      GetValues(N_Points, zeta, i, Values);
      TFEDatabase2D::RegisterJointValues2D(BaseFunct, LineQuadFormula, 
                       i, Values);
    } // endif

    Values=TFEDatabase2D::GetJointDerivatives2D(BaseFunct, LineQuadFormula, i,
                                                D10);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension*BaseVectDim];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension*BaseVectDim;
      GetDerivatives(D10, N_Points, zeta, i, Values);
      TFEDatabase2D::RegisterJointDerivatives2D(BaseFunct, LineQuadFormula, 
                       i, D10, Values);
    } // endif

    Values=TFEDatabase2D::GetJointDerivatives2D(BaseFunct, LineQuadFormula, i,
                                                D01);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension*BaseVectDim];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension*BaseVectDim;
      GetDerivatives(D01, N_Points, zeta, i, Values);
      TFEDatabase2D::RegisterJointDerivatives2D(BaseFunct, LineQuadFormula, 
                       i, D01, Values);
    } // endif
    
    // second derivatives
    // D20
    Values=TFEDatabase2D::GetJointDerivatives2D(BaseFunct, LineQuadFormula, i,
                                                D20);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension*BaseVectDim];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension*BaseVectDim;
      GetDerivatives(D20, N_Points, zeta, i, Values);
      TFEDatabase2D::RegisterJointDerivatives2D(BaseFunct, LineQuadFormula, 
            i, D20, Values);
    }
    
    // D11
    Values=TFEDatabase2D::GetJointDerivatives2D(BaseFunct, LineQuadFormula, i,
                                                D11);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension*BaseVectDim];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension*BaseVectDim;
      GetDerivatives(D11, N_Points, zeta, i, Values);
      TFEDatabase2D::RegisterJointDerivatives2D(BaseFunct, LineQuadFormula, 
            i, D11, Values);
    }
  
    // D02
    Values=TFEDatabase2D::GetJointDerivatives2D(BaseFunct, LineQuadFormula, i,
                                                D02);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension*BaseVectDim];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension*BaseVectDim;
      GetDerivatives(D02, N_Points, zeta, i, Values);
      TFEDatabase2D::RegisterJointDerivatives2D(BaseFunct, LineQuadFormula, 
            i, D02, Values);
    }

  } // endfor
}

/** make data on reference element */
void TBaseFunct2D::MakeRefElementData(QuadFormula2D QuadFormula)
{
  int j;
  double **Values, *AllValues;
  TQuadFormula2D *qf2;

  qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
  int N_Points =TFEDatabase2D::GetQuadFormula2D(QuadFormula)->GetN_QuadPoints();
  
  // D00
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D00);
  if( Values==NULL)
  {
    Values = new double* [N_Points];
    AllValues = new double [N_Points*Dimension*BaseVectDim];
    for(j=0;j<N_Points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D00, qf2, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D00, Values);
  }

  // D10
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D10);
  if( Values==NULL)
  {
    Values = new double* [N_Points];
    AllValues = new double [N_Points*Dimension*BaseVectDim];
    for(j=0;j<N_Points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D10, qf2, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D10, Values);
  }

  // D01
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D01);
  if( Values==NULL)
  {
    Values = new double* [N_Points];
    AllValues = new double [N_Points*Dimension*BaseVectDim];
    for(j=0;j<N_Points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D01, qf2, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D01, Values);
  }

  // D20
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D20);
  if( Values==NULL)
  {
    Values = new double* [N_Points];
    AllValues = new double [N_Points*Dimension*BaseVectDim];
    for(j=0;j<N_Points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D20, qf2, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D20, Values);
  }

  // D11
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D11);
  if( Values==NULL)
  {
    Values = new double* [N_Points];
    AllValues = new double [N_Points*Dimension*BaseVectDim];
    for(j=0;j<N_Points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D11, qf2, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D11, Values);
  }

  // D02
  Values=TFEDatabase2D::GetRefElementValues(BaseFunct, QuadFormula, D02);
  if( Values==NULL)
  {
    Values = new double* [N_Points];
    AllValues = new double [N_Points*Dimension*BaseVectDim];
    for(j=0;j<N_Points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D02, qf2, Values);
    TFEDatabase2D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D02, Values);
  }

}

TGridCell *TBaseFunct2D::GenerateRefElement()
{
  TGridCell *Cell;
  TVertex *v[4];
  TJointEqN *b[4];
  int i;

  switch(RefElement)
  {
    case BFUnitSquare:
      Cell = new TGridCell(TDatabase::RefDescDB[Rectangle], 0);
#ifdef __2D__
      v[0] = new TVertex(-1.0, -1.0);
      v[1] = new TVertex( 1.0, -1.0);
      v[2] = new TVertex( 1.0,  1.0);
      v[3] = new TVertex(-1.0,  1.0);
#else
      v[0] = new TVertex(-1.0, -1.0, 0.0);
      v[1] = new TVertex( 1.0, -1.0, 0.0);
      v[2] = new TVertex( 1.0,  1.0, 0.0);
      v[3] = new TVertex(-1.0,  1.0, 0.0);
#endif
    
      Cell->SetClipBoard(Refinement);

      for(i=0;i<4;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetVertex(i, v[i]);
        Cell->SetJoint(i, b[i]);
      } // endfor i
      break;
    case BFUnitTriangle:
      Cell = new TGridCell(TDatabase::RefDescDB[Triangle], 0);
#ifdef __2D__
      v[0] = new TVertex(0.0, 0.0);
      v[1] = new TVertex(1.0, 0.0);
      v[2] = new TVertex(0.0, 1.0);
#else
      v[0] = new TVertex(0.0, 0.0, 0.0);
      v[1] = new TVertex(1.0, 0.0, 0.0);
      v[2] = new TVertex(0.0, 1.0, 0.0);
#endif
    
      Cell->SetClipBoard(Refinement);

      for(i=0;i<3;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetVertex(i, v[i]);
        Cell->SetJoint(i, b[i]);
      } // endfor i
      break;
    default:
      cout << "unknown reference element!" << endl;
  }

  return Cell;
}

/** change basis functions on cell if needed */
void TBaseFunct2D::ChangeBF(TCollection *Coll, TBaseCell *Cell, double *Values)
{
  int *JointArray, i, j;
  TJoint *joint;
  TBaseCell *neigh;
  int OwnNum, NeighNum;
  int N_Joints;

  if(BF2Change == NULL) return;

  // /*
  OwnNum = Coll->GetIndex(Cell);
  N_Joints = Cell->GetN_Joints();
  
  for(j=0;j<N_Joints;j++)
  {
    joint = Cell->GetJoint(j);
    neigh = joint->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        JointArray = BF2Change[j];
        for(i=0;i<N_BF2Change;i++)
          Values[JointArray[i]] = -Values[JointArray[i]];
      } // endif NeighNum < OwnNum
    } // endif neigh
  } // endfor j
  // */

  /*
  switch(RefElement)
  {
    case BFUnitSquare:
      if(Cell->GetVertex(0) > Cell->GetVertex(1))
      {
        JointArray = BF2Change[0];
        for(i=0;i<N_BF2Change;i++)
          Values[JointArray[i]] = -Values[JointArray[i]];
      }
      if(Cell->GetVertex(1) > Cell->GetVertex(2))
      {
        JointArray = BF2Change[1];
        for(i=0;i<N_BF2Change;i++)
        {
          Values[JointArray[i]] = -Values[JointArray[i]];
        }
      }
      if(Cell->GetVertex(2) > Cell->GetVertex(3))
      {
        JointArray = BF2Change[2];
        for(i=0;i<N_BF2Change;i++)
          Values[JointArray[i]] = -Values[JointArray[i]];
      }
      if(Cell->GetVertex(3) > Cell->GetVertex(0))
      {
        JointArray = BF2Change[3];
        for(i=0;i<N_BF2Change;i++)
          Values[JointArray[i]] = -Values[JointArray[i]];
      }
    break;

    case BFUnitTriangle:
      if(Cell->GetVertex(0) > Cell->GetVertex(1))
      {
        JointArray = BF2Change[0];
        for(i=0;i<N_BF2Change;i++)
          Values[JointArray[i]] = -Values[JointArray[i]];
      }
      if(Cell->GetVertex(1) > Cell->GetVertex(2))
      {
        JointArray = BF2Change[1];
        for(i=0;i<N_BF2Change;i++)
          Values[JointArray[i]] = -Values[JointArray[i]];
      }
      if(Cell->GetVertex(2) > Cell->GetVertex(0))
      {
        JointArray = BF2Change[2];
        for(i=0;i<N_BF2Change;i++)
          Values[JointArray[i]] = -Values[JointArray[i]];
      }
    break;
  } // end switch
  */
}

/** change basis functions on cell in all points if needed */
void TBaseFunct2D::ChangeBF(TCollection *Coll, TBaseCell *Cell, int N_Points, double **Values)
{
  int *JointArray, i, j, k;
  double *Array;
  TJoint *joint;
  TBaseCell *neigh;
  int OwnNum, NeighNum;
  int N_Joints;

  if(BF2Change == NULL) return;

  // /*
  OwnNum = Coll->GetIndex(Cell);
  N_Joints = Cell->GetN_Joints();

  for(j=0;j<N_Joints;j++)
  {
    joint = Cell->GetJoint(j);
    neigh = joint->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
      {
        JointArray = BF2Change[j];
        for(k=0;k<N_Points;k++)
        {
          Array = Values[k];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor k
      } // endif NeighNum < OwnNum
    } // endif neigh
  } // endfor j
  // */

  /*
  switch(RefElement)
  {
    case BFUnitSquare:
      if(Cell->GetVertex(0) > Cell->GetVertex(1))
      {
        JointArray = BF2Change[0];
        for(j=0;j<N_Points;j++)
        {
          Array = Values[j];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor j
      }
      if(Cell->GetVertex(1) > Cell->GetVertex(2))
      {
        JointArray = BF2Change[1];
        for(j=0;j<N_Points;j++)
        {
          Array = Values[j];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor j
      }
      if(Cell->GetVertex(2) > Cell->GetVertex(3))
      {
        JointArray = BF2Change[2];
        for(j=0;j<N_Points;j++)
        {
          Array = Values[j];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor j
      }
      if(Cell->GetVertex(3) > Cell->GetVertex(0))
      {
        JointArray = BF2Change[3];
        for(j=0;j<N_Points;j++)
        {
          Array = Values[j];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor j
      }
    break;

    case BFUnitTriangle:
      if(Cell->GetVertex(0) > Cell->GetVertex(1))
      {
        JointArray = BF2Change[0];
        for(j=0;j<N_Points;j++)
        {
          Array = Values[j];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor j
      }
      if(Cell->GetVertex(1) > Cell->GetVertex(2))
      {
        JointArray = BF2Change[1];
        for(j=0;j<N_Points;j++)
        {
          Array = Values[j];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor j
      }
      if(Cell->GetVertex(2) > Cell->GetVertex(0))
      {
        JointArray = BF2Change[2];
        for(j=0;j<N_Points;j++)
        {
          Array = Values[j];
          for(i=0;i<N_BF2Change;i++)
            Array[JointArray[i]] = -Array[JointArray[i]];
        } // endfor j
      }
    break;
  } // end switch
  */
}
