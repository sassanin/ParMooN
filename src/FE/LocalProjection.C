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
// LocalProjection.C
//
// Purpose:   routines for local projection stabilization
//
// Author:    Gunar Matthies  2007/03/06
//
// =======================================================================

#include <Database.h>
#include <MooNMD_Io.h>

#ifdef __2D__
  #include <FEDatabase2D.h>
  #include <FEFunction2D.h>
  #include <NodalFunctional2D.h>
#else  
  #include <FEDatabase3D.h>
  #include <FEFunction3D.h>
  #include <NodalFunctional3D.h>
#endif

#include <ConvDiff.h>
#include <LinAlg.h>

#include <string.h>
#include <stdlib.h>

#ifdef __2D__
FE2D GetElement2D(TBaseCell *cell, int CoarseOrder)
{
  FE2D ele = (FE2D)0;
  Shapes shapetype;

  shapetype = cell->GetType();
  switch(shapetype)
  {
    // regularly refined quadrilateral
    case Quadrangle:
    case Parallelogram:
    case Rectangle:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_Q0_2D_Q_M;
        break;

        case 1:
          ele = D_P1_2D_Q_M;
        break;

        case 2:
          ele = D_P2_2D_Q_M;
        break;

        case 3:
          ele = D_P3_2D_Q_M;
        break;

        case 4:
          ele = D_P4_2D_Q_M;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_Q00_2D_Q_M;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 4" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break; // end regularly refined quadrilateral

    case Triangle:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_P0_2D_T_A;
        break;

        case 1:
          ele = D_P1_2D_T_A;
        break;

        case 2:
          ele = D_P2_2D_T_A;
        break;

        case 3:
          ele = D_P3_2D_T_A;
        break;

        case 4:
          ele = D_P4_2D_T_A;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_P00_2D_T_A;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 4" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break;
    default:
      OutPut("Invalid shape" << endl);
      exit(-1);
  } // end switch reftype
  return ele;
}

/** Navier--Stokes type 2 (NSTYPE==2) with C*/
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;
  double *CEntries;
  int *CRowPtr, *CKCol;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;

  N_Active = A->GetActiveBound();

  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
      t -= value * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s -= value1 * u1[index] + value2 * u2[index];
    } // endfor j
    r3[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    }
    r3[i] -= s;
  } // endfor i
}

void Defect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], B[2], B[3], B[4], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF, TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;
  double *CEntries;
  int *CRowPtr, *CKCol;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];

  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
      t += value * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    v1[i] += s;
    v2[i] += t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], B[3], B[4], x, y);
  return;
}

// stabilisation of full gradient (velocity or pressure)
void UltraLocalProjection(void* A, bool ForPressure)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    OutPut("Local projection stabilization is implemented only for full gradient!" << endl);
    exit(-1);
  }

  if(ForPressure)
  {
    lpcoeff = -(TDatabase::ParamDB->LP_PRESSURE_COEFF);
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  if(ForPressure)
  {
    fespace = (TFESpace2D*)(((TMatrix2D *)A)->GetStructure()->GetTestSpace());
    ActiveBound = -1;
    RowPtr = ((TMatrix2D *)A)->GetRowPtr();
    KCol = ((TMatrix2D *)A)->GetKCol();
    Entries = ((TMatrix2D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix2D *)A)->GetFESpace();
    ActiveBound = fespace->GetActiveBound();
    RowPtr = ((TSquareMatrix2D *)A)->GetRowPtr();
    KCol = ((TSquareMatrix2D *)A)->GetKCol();
    Entries = ((TSquareMatrix2D *)A)->GetEntries();
//     cout << "not for pressure" << endl;
  }

  
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ]
                                           +H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // 


void AddDeformationTensorTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixA11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s22;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;

  fespace = A11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();
  EntriesA11 = A11->GetEntries();
  EntriesA12 = A12->GetEntries();
  EntriesA21 = A21->GetEntries();
  EntriesA22 = A22->GetEntries();
  // cout << "" << endl;


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixA11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA22, 0, N_DOF*N_DOF*SizeOfDouble);
    // since A21=transpose(A12) A21 is not explicitly needed

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixA11[k*N_DOF+l] += w*((ChildValueX[k]*ChildValueX[l]) + (0.5*ChildValueY[k]*ChildValueY[l]));
          LocMatrixA12[k*N_DOF+l] += w*0.5*ChildValueX[k]*ChildValueY[l];
          LocMatrixA22[k*N_DOF+l] += w*((ChildValueY[k]*ChildValueY[l]) + (0.5*ChildValueX[k]*ChildValueX[l]));
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *((H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m])
		   + (0.5*H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]));
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *0.5*H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *((H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF])
		  + (0.5*H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m]));
          }
        LocMatrixA11[l*N_DOF+m] += s11;
        LocMatrixA12[l*N_DOF+m] += s12;
        LocMatrixA22[l*N_DOF+m] += s22;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += ((P[i2*2*N_DOF+l] * H[i2*2*N_DOF+m]) + (0.5*P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]));
          s12 += 0.5*P[i2*2*N_DOF+l] * H[i2*2*N_DOF+m+N_DOF];
          s22 += ((P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]) + (0.5*P[i2*2*N_DOF+l] * H[i2*2*N_DOF+m]));
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += ((P[i1*2*N_DOF+m] * H[i1*2*N_DOF+l]) + (0.5*P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]));
          s12 += 0.5*P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l];
          s22 += ((P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]) + (0.5*P[i1*2*N_DOF+m] * H[i1*2*N_DOF+l]));
        }
        LocMatrixA11[l*N_DOF+m] -= s11;
        LocMatrixA12[l*N_DOF+m] -= s12;
        LocMatrixA22[l*N_DOF+m] -= s22;
//         LocMatrixA21[m*N_DOF+l] = LocMatrixA12[l*N_DOF+m];
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesA11[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA11[l*N_DOF+m];
              EntriesA12[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[l*N_DOF+m];
              // since A21=transpose(A12)
              EntriesA21[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[m*N_DOF+l];
              EntriesA22[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA22[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddDeformationTensorTerm



void AddDeformationTensorTerm_2PhaseOrImpDropFlow(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points, Phase_No;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixA11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s22;
  int i1, i2;
  double hK, lpcoeff=0.0;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;

  fespace = A11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();
  EntriesA11 = A11->GetEntries();
  EntriesA12 = A12->GetEntries();
  EntriesA21 = A21->GetEntries();
  EntriesA22 = A22->GetEntries();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    
    if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
    {
    Phase_No = cell->GetRegionID();
    if(Phase_No == 0)
      lpcoeff = TDatabase::ParamDB->DELTA1 * (1.0 - TDatabase::ParamDB->P13);
    else if(Phase_No == 1)
      lpcoeff = TDatabase::ParamDB->DELTA1 * (1.0 - TDatabase::ParamDB->P14);
    else
    {
     cout<<"Invalid Phase No. !!!!!\n";
     exit(4711);
    }
    }
    else if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
    {
      lpcoeff = TDatabase::ParamDB->DELTA1 * (1.0 - TDatabase::ParamDB->P13);
    }
    else
   {
    OutPut("Implemented only for two-phase and impinging droplet flows \n");
    OutPut("Change TWO_PHASE_FLOW or FREE_SURFACE_FLOW  to 1 in dat file !!!!! " << endl);
    exit(4711);
   }
    
    if(lpcoeff == 0.0)
    {
     continue; 
    }
    else
    {
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixA11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA22, 0, N_DOF*N_DOF*SizeOfDouble);
    // since A21=transpose(A12) A21 is not explicitly needed

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixA11[k*N_DOF+l] += w*((ChildValueX[k]*ChildValueX[l]) + (0.5*ChildValueY[k]*ChildValueY[l]));
          LocMatrixA12[k*N_DOF+l] += w*0.5*ChildValueX[k]*ChildValueY[l];
          LocMatrixA22[k*N_DOF+l] += w*((ChildValueY[k]*ChildValueY[l]) + (0.5*ChildValueX[k]*ChildValueX[l]));
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *((H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m])
		   + (0.5*H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]));
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *0.5*H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *((H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF])
		  + (0.5*H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m]));
          }
        LocMatrixA11[l*N_DOF+m] += s11;
        LocMatrixA12[l*N_DOF+m] += s12;
        LocMatrixA22[l*N_DOF+m] += s22;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += ((P[i2*2*N_DOF+l] * H[i2*2*N_DOF+m]) + (0.5*P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]));
          s12 += 0.5*P[i2*2*N_DOF+l] * H[i2*2*N_DOF+m+N_DOF];
          s22 += ((P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]) + (0.5*P[i2*2*N_DOF+l] * H[i2*2*N_DOF+m]));
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += ((P[i1*2*N_DOF+m] * H[i1*2*N_DOF+l]) + (0.5*P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]));
          s12 += 0.5*P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l];
          s22 += ((P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]) + (0.5*P[i1*2*N_DOF+m] * H[i1*2*N_DOF+l]));
        }
        LocMatrixA11[l*N_DOF+m] -= s11;
        LocMatrixA12[l*N_DOF+m] -= s12;
        LocMatrixA22[l*N_DOF+m] -= s22;
//         LocMatrixA21[m*N_DOF+l] = LocMatrixA12[l*N_DOF+m];
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesA11[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA11[l*N_DOF+m];
              EntriesA12[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[l*N_DOF+m];
              // since A21=transpose(A12)
              EntriesA21[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[m*N_DOF+l];
              EntriesA22[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA22[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
    
    } // if(lpcoeff == 0.0)
    
    
  } // endfor i
} // AddDeformationTensorTerm_2PhaseFlow



void AddDeformationTensorTerm_2PhaseOrImpDropFlow_3DAxial(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points, Phase_No;
  double *xi, *eta, *weights, r;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave1[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D], Gsave2[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[3*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[3*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixA11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA21[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s21, s22;
  int i1, i2;
  double hK, lpcoeff=0.0;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;

  fespace = A11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();
  EntriesA11 = A11->GetEntries();
  EntriesA12 = A12->GetEntries();
  EntriesA21 = A21->GetEntries();
  EntriesA22 = A22->GetEntries();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    
    if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
    {
    Phase_No = cell->GetRegionID();
    if(Phase_No == 0)
      lpcoeff = TDatabase::ParamDB->DELTA1 * (1.0 - TDatabase::ParamDB->P13);
    else if(Phase_No == 1)
      lpcoeff = TDatabase::ParamDB->DELTA1 * (1.0 - TDatabase::ParamDB->P14);
    else
    {
     cout<<"Invalid Phase No. !!!!!\n";
     exit(4711);
    }
    }
    else if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
    {
      lpcoeff = TDatabase::ParamDB->DELTA1 * (1.0 - TDatabase::ParamDB->P13);
    }
    else
   {
    OutPut("Implemented only for two-phase and impinging droplet flows \n");
    OutPut("Change TWO_PHASE_FLOW or FREE_SURFACE_FLOW  to 1 in dat file !!!!! " << endl);
    exit(4711);
   }
    
    if(lpcoeff == 0.0)
    {
     continue; 
    }
    else
    {
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    memset(Gsave1, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    memset(Gsave2, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      
      r  = fabs(X[j]);
    if(r<1e-12)
    {
     OutPut("check Local Projection Stabilization  r value zero !!!!! "<< X[j] <<endl);
     OutPut("Quad formula: Change all integral points as internal points"<<endl);
     exit(0);
    }
    
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l]*r;
	  Gsave1[k*N_CoarseDOF+l] += val*CoarseValue[l]*r;
	  Gsave2[k*N_CoarseDOF+l] += val*CoarseValue[l]/r;
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;


    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*3*N_DOF*SizeOfDouble);
    memset(P, 0, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    memset(LocMatrixA11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA21, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA22, 0, N_DOF*N_DOF*SizeOfDouble);
    

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues [j];
      w = AbsDetjk[j]*weights[j];
      
      r  = fabs(X[j]);
    if(r<1e-12)
    {
     OutPut("check Local Projection Stabilization  r value zero !!!!! "<< X[j] <<endl);
     OutPut("Quad formula: Change all integral points as internal points"<<endl);
     exit(0);
    }
    
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*3*N_DOF+l      ] += val*ChildValueX[l]*r;
          H[k*3*N_DOF+l+N_DOF] += val*ChildValueY[l]*r;
	  H[k*3*N_DOF+l+2*N_DOF] += val*ChildValue[l]*r;
	  
	  P[k*3*N_DOF+l      ] += val*r*ChildValueX[l];
          P[k*3*N_DOF+l+N_DOF] += val*r*ChildValueY[l];
	  P[k*3*N_DOF+l+2*N_DOF] += val*ChildValue[l]/r; 
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixA11[k*N_DOF+l] += w*r*((ChildValueX[k]*ChildValueX[l]) + (0.5*ChildValueY[k]*ChildValueY[l]) + (ChildValue[k]*ChildValue[l])/(r*r));
          LocMatrixA12[k*N_DOF+l] += w*r*0.5*ChildValueY[k]*ChildValueX[l];
	  LocMatrixA21[k*N_DOF+l] += w*r*0.5*ChildValueX[k]*ChildValueY[l];
          LocMatrixA22[k*N_DOF+l] += w*r*((ChildValueY[k]*ChildValueY[l]) + (0.5*ChildValueX[k]*ChildValueX[l]));
        }
      }
    } // end for j
//     memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

 
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3*N_DOF, 3*N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
	s21 = 0;
        s22 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave1[i1*N_CoarseDOF+i2] *(H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m]);
	    s11 += Gsave1[i1*N_CoarseDOF+i2] *(0.5*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF]);  
	    s11 += Gsave2[i1*N_CoarseDOF+i2] *(H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF]);
	    
            s12 += Gsave1[i1*N_CoarseDOF+i2]*0.5*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m];
	    
	    s21 += Gsave1[i1*N_CoarseDOF+i2]*0.5*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m+N_DOF];
		  
            s22 += Gsave1[i1*N_CoarseDOF+i2]*((H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF]) + (0.5*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m]));
          }
        LocMatrixA11[l*N_DOF+m] += s11;
        LocMatrixA12[l*N_DOF+m] += s12;
	LocMatrixA21[l*N_DOF+m] += s21;
        LocMatrixA22[l*N_DOF+m] += s22;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
	s21 = 0;
        s22 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += (P[i2*3*N_DOF+l] * H[i2*3*N_DOF+m]) + (0.5*P[i2*3*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF]) + (P[i2*3*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF]);
          s12 += 0.5*P[i2*3*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m];
	  s21 += 0.5*P[i2*3*N_DOF+l] * H[i2*3*N_DOF+m+N_DOF];
          s22 += ((P[i2*3*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF]) + (0.5*P[i2*3*N_DOF+l] * H[i2*3*N_DOF+m]));
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += ((P[i1*3*N_DOF+m] * H[i1*3*N_DOF+l]) + (0.5*P[i1*3*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF]) + (P[i1*3*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF]));
          s12 += 0.5*P[i1*3*N_DOF+m] * H[i1*3*N_DOF+l+N_DOF];
	  s21 += 0.5*P[i1*3*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l];
          s22 += ((P[i1*3*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF]) + (0.5*P[i1*3*N_DOF+m] * H[i1*3*N_DOF+l]));
        }
        LocMatrixA11[l*N_DOF+m] -= s11;
        LocMatrixA12[l*N_DOF+m] -= s12;
	LocMatrixA21[l*N_DOF+m] -= s21;
        LocMatrixA22[l*N_DOF+m] -= s22;
      } // end for m
    } // end for l

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesA11[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA11[l*N_DOF+m];
              EntriesA12[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[l*N_DOF+m];
              EntriesA21[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA21[l*N_DOF+m];
              EntriesA22[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA22[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
    
    } // if(lpcoeff == 0.0)
    
    
  } // endfor i
} // AddDeformationTensorTerm_2PhaseFlow_3DAxial



double UltraLocalError(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double error, locerror;
  double exactval[4];

  error = 0.0;

  fespace = uh->GetFESpace2D();
  Values = uh->GetValues();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    // only two right-hand sides (x and y derivative)
    memset(H, 0, N_CoarseDOF*2*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        val = Values[DOF[k]];
        valx += ChildValueX[k]*val;
        valy += ChildValueY[k]*val;
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k*2  ] += val*valx;
        H[k*2+1] += val*valy;
      } // end for k

      // grad-grad term
      locerror += w*(valx*valx + valy*valy);

    } // end for j
    memcpy(P, H, N_CoarseDOF*2*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2, 2);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2  ]*H[i2*2  ]
                                       +H[i1*2+1]*H[i2*2+1]);
    locerror += s;

    // grad-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2*2  ] * H[i2*2  ];
      s += P[i2*2+1] * H[i2*2+1];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1*2  ] * H[i1*2  ];
      s += P[i1*2+1] * H[i1*2+1];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

void AddStreamlineTerm(TSquareMatrix2D* A, TFEFunction2D *uh1,
                       TFEFunction2D *uh2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2;
  double BValue[MaxN_BaseFunctions2D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm

void AddStreamlineTermPWConst(TSquareMatrix2D* A, TFEFunction2D *uh1,
                              TFEFunction2D *uh2,
                              double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2;
  double BValue[MaxN_BaseFunctions2D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    // calculate pw constant approximation of velocity field (uh1, uh2)
    val  = 0.0;
    valx = 0.0;
    valy = 0.0;
    for(j=0;j<N_Points;j++)
    {
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        val  += w*ChildValue[k]*1;
        valx += w*ChildValue[k]*Values1[l];
        valy += w*ChildValue[k]*Values2[l];
      }
    }
    valx /= val;
    valy /= val;

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      /*
      valx = 0.0;
      valy = 0.0;
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
      }
      */
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTermPWConst

void AddDivergenceTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixA11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s22;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;

  fespace = A11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();
  EntriesA11 = A11->GetEntries();
  EntriesA12 = A12->GetEntries();
  EntriesA21 = A21->GetEntries();
  EntriesA22 = A22->GetEntries();
  // cout << "" << endl;


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixA11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA22, 0, N_DOF*N_DOF*SizeOfDouble);
    // since A21=transpose(A12) A21 is not explicitly needed

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixA11[k*N_DOF+l] += w*ChildValueX[k]*ChildValueX[l];
          LocMatrixA12[k*N_DOF+l] += w*ChildValueX[k]*ChildValueY[l];
          LocMatrixA22[k*N_DOF+l] += w*ChildValueY[k]*ChildValueY[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ];
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF];
          }
        LocMatrixA11[l*N_DOF+m] += s11;
        LocMatrixA12[l*N_DOF+m] += s12;
        LocMatrixA22[l*N_DOF+m] += s22;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s12 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m+N_DOF];
          s22 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s12 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l      ];
          s22 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrixA11[l*N_DOF+m] -= s11;
        LocMatrixA12[l*N_DOF+m] -= s12;
        LocMatrixA22[l*N_DOF+m] -= s22;
//         LocMatrixA21[m*N_DOF+l] = LocMatrixA12[l*N_DOF+m];
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesA11[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA11[l*N_DOF+m];
              EntriesA12[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[l*N_DOF+m];
              // since A21=transpose(A12)
              EntriesA21[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[m*N_DOF+l];
              EntriesA22[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA22[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddDivergenceTerm

/** Navier--Stokes type 4 (NSTYPE==4) with C*/
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    TMatrix *C,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2,value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
  double *CEntries;
  int *CRowPtr, *CKCol;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  N_Active = A11->GetActiveBound();
  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A21Entries[j];
      value3 = A22Entries[j];
      s += value * u1[index] + value1 * u2[index];
      t += value2* u1[index] + value3 * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i
  // Dirichlet and hanging nodes
  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    v1[i] += s;
    v2[i] += t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], B[4], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   TMatrix *C,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
  int N_Active;
  double *CEntries;
  int *CRowPtr, *CKCol;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;

  N_Active = A11->GetActiveBound();

  j = ARowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A21Entries[j];
      value3 = A22Entries[j];
      s -= value * u1[index] + value1 * u2[index];
      t -= value2* u1[index] + value3 * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s -= value1 * u1[index] + value2 * u2[index];
    } // endfor j
    r3[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    }
    r3[i] -= s;
  } // endfor i
}

void Defect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], B[4], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

double UltraLocalErrorDivergence(TFEFunction2D *uh1, TFEFunction2D *uh2,
                       DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values1, *Values2;
  double error, locerror;
  double exactval1[4], exactval2[4];
  double div;

  error = 0.0;

  fespace = uh1->GetFESpace2D();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    // only one right-hand side (div (u-uh))
    memset(H, 0, N_CoarseDOF*1*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate x derivative of uh1 and y derivative of uh2 in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValueX[k]*Values1[l];
        valy += ChildValueY[k]*Values2[l];
      }

      // get x derivative of u1 and y derivative of u2
      ExactU1(X[j], Y[j], exactval1);
      ExactU2(X[j], Y[j], exactval2);

      valx -= exactval1[1];
      valy -= exactval2[2];

      div = valx+valy;

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*div;
      } // end for k

      // id-id term
      locerror += w*div*div;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

double UltraLocalErrorStreamline(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double *BValues1, *BValues2;
  double error, locerror;
  double exactval[4];
  double valb1, valb2, StreamlineDerivative;

  error = 0.0;

  fespace = uh->GetFESpace2D();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values = uh->GetValues();

    BValues1 = b1->GetValues();
    BValues2 = b2->GetValues();
    
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    // only one right-hand side (b.grad(ui))
    memset(H, 0, N_CoarseDOF*1*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];

      // calculate derivative of uh and components of b in this quadrature point
      valx = 0.0;
      valy = 0.0;
      valb1 = 0.0;
      valb2 = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valb1 += ChildValue[k]*BValues1[l];
        valb2 += ChildValue[k]*BValues2[l];
        valx += ChildValueX[k]*Values[l];
        valy += ChildValueY[k]*Values[l];
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      StreamlineDerivative = valb1*valx + valb2*valy;

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*StreamlineDerivative;
      } // end for k

      // id-id term
      locerror += w*StreamlineDerivative*StreamlineDerivative;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

double UltraLocalErrorStreamlinePWConst(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double *BValues1, *BValues2;
  double error, locerror;
  double exactval[4];
  double valb1, valb2, StreamlineDerivative;

  error = 0.0;

  fespace = uh->GetFESpace2D();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values = uh->GetValues();

    BValues1 = b1->GetValues();
    BValues2 = b2->GetValues();
    
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    // only one right-hand side (b.grad(ui))
    memset(H, 0, N_CoarseDOF*1*SizeOfDouble);

    // calculate pw constant approximation of velocity field (uh1, uh2)
    val  = 0.0;
    valx = 0.0;
    valy = 0.0;
    for(j=0;j<N_Points;j++)
    {
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        val  += w*ChildValue[k]*1;
        valx += w*ChildValue[k]*BValues1[l];
        valy += w*ChildValue[k]*BValues2[l];
      }
    }
    valb1 = valx/val;
    valb2 = valy/val;

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];

      // calculate derivative of uh and components of b in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValueX[k]*Values[l];
        valy += ChildValueY[k]*Values[l];
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      StreamlineDerivative = valb1*valx + valb2*valy;

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*StreamlineDerivative;
      } // end for k

      // id-id term
      locerror += w*StreamlineDerivative*StreamlineDerivative;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalErrorPWConst



// for conformation stress tensor equation
void AddStreamlineTerm(TSquareMatrix2D *G11,TSquareMatrix2D *G22,
                       TSquareMatrix2D *G33, TFEFunction2D *u1, TFEFunction2D *u2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *GlobalNumbers_vel, *BeginIndex, *BeginIndex_vel, *DOF, *DOF_vel;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF, N_DOF_vel;
  TCollection *Coll;
  TFESpace2D *fespace, *fespace_vel;
  FE2D CurrEleID,CurrEleID_vel, UsedElements[3];
  int N_UsedElements;
  TFE2D *CurrentElement, *CurrentElement_vel, *CoarseElement;
  TBaseFunct2D *BF, *BF_vel, *CoarseBF;
  BaseFunct2D BF_ID, BF_ID_vel, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[3] = { false, false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrixG[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG22, *EntriesG33;
  double *Values1, *Values2;
    double BValue[MaxN_BaseFunctions2D];
  double delta, norm_u;
  fespace = G11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = G11->GetRowPtr();
  KCol = G11->GetKCol();
  EntriesG11 = G11->GetEntries();
  EntriesG22 = G22->GetEntries();
  EntriesG33 = G33->GetEntries();

  fespace_vel = u1->GetFESpace2D();
  
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  BeginIndex_vel = fespace_vel->GetBeginIndex();
  GlobalNumbers_vel = fespace_vel->GetGlobalNumbers();
  
  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CurrEleID_vel = fespace_vel->GetFE2D(i, cell);
    CurrentElement_vel = TFEDatabase2D::GetFE2D(CurrEleID_vel);

    BF_vel = CurrentElement_vel->GetBaseFunct2D();
    BF_ID_vel = BF_vel->GetID();
    N_DOF_vel = BF_vel->GetDimension();
    
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1), proj (index 0) space, velocity (index 2) space
    N_UsedElements = 3;
    UsedElements[1] = CurrEleID;
    UsedElements[2] = CurrEleID_vel;
    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];
    DOF_vel = GlobalNumbers_vel + BeginIndex_vel[i];
      Values1 = u1->GetValues();
      Values2 = u2->GetValues();
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID_vel, D00);
    
    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrixG, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;
      
        // compute components of uh in j
      for(k=0;k<N_DOF_vel;k++)
      {
        l = DOF_vel[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);
    
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;

        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s += Gsave[i1*N_CoarseDOF+i2]*H[i1*N_DOF+l]*H[i2*N_DOF+m];
          }
        LocMatrixG[l*N_DOF+m] += s;

      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;

        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l] * H[i2*N_DOF+m];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m] * H[i1*N_DOF+l];
        }
        LocMatrixG[l*N_DOF+m] -= s;

      } // end for m
    } // end for l
    

   // norm_u = sqrt((valx*valx) + (valy*valy));
   norm_u = 1.0;
    delta =  lpcoeff*pow(hK,lpexponent)/norm_u;

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] += delta*LocMatrixG[l*N_DOF+m];
              EntriesG22[p] += delta*2.0*LocMatrixG[l*N_DOF+m];
              EntriesG33[p] += delta*LocMatrixG[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
  
} 



// for conformation stress tensor equation
void AddStreamlineTerm_2PhaseOrImpDropFlow_3DAxial(TSquareMatrix2D *G11,TSquareMatrix2D *G22,
                       TSquareMatrix2D *G33, TFEFunction2D *u1, TFEFunction2D *u2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *GlobalNumbers_vel, *BeginIndex, *BeginIndex_vel, *DOF, *DOF_vel;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF, N_DOF_vel;
  TCollection *Coll;
  TFESpace2D *fespace, *fespace_vel;
  FE2D CurrEleID,CurrEleID_vel, UsedElements[3];
  int N_UsedElements, Phase_No;
  TFE2D *CurrentElement, *CurrentElement_vel, *CoarseElement;
  TBaseFunct2D *BF, *BF_vel, *CoarseBF;
  BaseFunct2D BF_ID, BF_ID_vel, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[3] = { false, false, false };
  int N_Points;
  double *xi, *eta, *weights, r;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesV, *ChildValueV;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrixG[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG22, *EntriesG33;
  double *Values1, *Values2;
    double BValue[MaxN_BaseFunctions2D];
  double delta, norm_u;
  fespace = G11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = G11->GetRowPtr();
  KCol = G11->GetKCol();
  EntriesG11 = G11->GetEntries();
  EntriesG22 = G22->GetEntries();
  EntriesG33 = G33->GetEntries();

  fespace_vel = u1->GetFESpace2D();
  
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  BeginIndex_vel = fespace_vel->GetBeginIndex();
  GlobalNumbers_vel = fespace_vel->GetGlobalNumbers();
  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    
      if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
    {
        if(Phase_No == 0)
    { 
      if(TDatabase::ParamDB->PHASE1_TYPE == 1)
      {
	lpcoeff = 0.0;
      }
    } 
      else if(Phase_No == 1)
     {
      if(TDatabase::ParamDB->PHASE2_TYPE == 1)
      {
	lpcoeff = 0.0;
      }
     }
      else
    {
     cout<<"Invalid Phase No. !!!!!\n";
     exit(4711);
    }
    }
    else if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
    {
      if(TDatabase::ParamDB->PHASE1_TYPE == 1)
      lpcoeff = 0.0;
    }
    else
   {
    OutPut("Implemented only for two-phase and impinging droplet flows \n");
    OutPut("Change TWO_PHASE_FLOW or FREE_SURFACE_FLOW  to 1 in dat file !!!!! " << endl);
    exit(4711);
   }
    
    if(lpcoeff == 0.0)
    {
     continue; 
    }
    else
    {
    
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CurrEleID_vel = fespace_vel->GetFE2D(i, cell);
    CurrentElement_vel = TFEDatabase2D::GetFE2D(CurrEleID_vel);

    BF_vel = CurrentElement_vel->GetBaseFunct2D();
    BF_ID_vel = BF_vel->GetID();
    N_DOF_vel = BF_vel->GetDimension();
    
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1), proj (index 0) space, velocity (index 2) space
    N_UsedElements = 3;
    UsedElements[1] = CurrEleID;
    UsedElements[2] = CurrEleID_vel;
    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      r  = fabs(X[j]);
     if(r<1e-12)
     {
     OutPut("check Local Projection Stabilization  r value zero !!!!! "<< X[j] <<endl);
     OutPut("Quad formula: Change all integral points as internal points"<<endl);
     exit(0);
     }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l]*r;
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];
    DOF_vel = GlobalNumbers_vel + BeginIndex_vel[i];
      Values1 = u1->GetValues();
      Values2 = u2->GetValues();
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValuesV  = TFEDatabase2D::GetOrigElementValues(BF_ID_vel, D00);
    
    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrixG, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueV  = ChildValuesV[j];
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;
      
        // compute components of uh in j
      for(k=0;k<N_DOF_vel;k++)
      {
        l = DOF_vel[k];
        valx += ChildValueV[k]*Values1[l];
        valy += ChildValueV[k]*Values2[l];
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      
     r  = fabs(X[j]);
      if(r<1e-12)
     {
     OutPut("check Local Projection Stabilization  r value zero !!!!! "<< X[j] <<endl);
     OutPut("Quad formula: Change all integral points as internal points"<<endl);
     exit(0);
     }
      
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l] += val*BValue[l]*r;
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG[k*N_DOF+l] += w*BValue[k]*BValue[l]*r;
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);
    
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;

        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s += Gsave[i1*N_CoarseDOF+i2]*H[i1*N_DOF+l]*H[i2*N_DOF+m];
          }
        LocMatrixG[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l] * H[i2*N_DOF+m];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m] * H[i1*N_DOF+l];
        }
        LocMatrixG[l*N_DOF+m] -= s;

      } // end for m
    } // end for l
    

    delta =  lpcoeff*pow(hK,lpexponent);

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] += delta*LocMatrixG[l*N_DOF+m];
              EntriesG22[p] += delta*2.0*LocMatrixG[l*N_DOF+m];
              EntriesG33[p] += delta*LocMatrixG[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
    
    }
    
  } // endfor i
  
} 


// for conformation stress tensor equation
void AddStretchingTerm(TSquareMatrix2D **SQMATRICES, TFEFunction2D *u1, TFEFunction2D *u2,
                       double lpcoeff, double lpexponent, int OrderDiff)

{
 int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *GlobalNumbers_vel, *BeginIndex, *BeginIndex_vel, *DOF, *DOF_vel;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF, N_DOF_vel;
  TCollection *Coll;
    TFESpace2D *fespace, *fespace_vel;
  FE2D CurrEleID,CurrEleID_vel, UsedElements[3];
  int N_UsedElements;
  TFE2D *CurrentElement, *CurrentElement_vel, *CoarseElement;
  TBaseFunct2D *BF, *BF_vel, *CoarseBF;
  BaseFunct2D BF_ID, BF_ID_vel, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[3] = { false, false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[6*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[6*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **ChildValues_Tau, *ChildValue_Tau;
  double **PCValues;
  double *PCValue;
  double w, val, valu1, valu1x, valu1y, valu2, valu2x, valu2y;
  double LocMatrixG11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG21[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG23[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG32[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG33[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s21, s23, s32, s33;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG23, *EntriesG32, *EntriesG33;
  double *Values1, *Values2;
  double delta, norm_u;
  fespace = SQMATRICES[0]->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = SQMATRICES[0]->GetRowPtr();
  KCol = SQMATRICES[0]->GetKCol();
  
  EntriesG11 = SQMATRICES[0]->GetEntries();
  EntriesG12 = SQMATRICES[1]->GetEntries();
  EntriesG21 = SQMATRICES[2]->GetEntries();
  EntriesG23 = SQMATRICES[4]->GetEntries();
  EntriesG32 = SQMATRICES[5]->GetEntries();
  EntriesG33 = SQMATRICES[6]->GetEntries();
  
 
  fespace_vel = u1->GetFESpace2D();
  
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  BeginIndex_vel = fespace_vel->GetBeginIndex();
  GlobalNumbers_vel = fespace_vel->GetGlobalNumbers();
  
  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
   for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CurrEleID_vel = fespace_vel->GetFE2D(i, cell);
    CurrentElement_vel = TFEDatabase2D::GetFE2D(CurrEleID_vel);

    BF_vel = CurrentElement_vel->GetBaseFunct2D();
    BF_ID_vel = BF_vel->GetID();
    N_DOF_vel = BF_vel->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;
   
    UsedElements[0] = GetElement2D(cell, CoarseOrder);

// approx (index 1), proj (index 0) space, velocity (index 2) space
    N_UsedElements = 3;
    UsedElements[1] = CurrEleID;
    UsedElements[2] = CurrEleID_vel;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);
  //  cout<<N_CoarseDOF<<"\n";
    
    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];
     DOF_vel = GlobalNumbers_vel + BeginIndex_vel[i];
      Values1 = u1->GetValues();
      Values2 = u2->GetValues();
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID_vel, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID_vel, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID_vel, D00);
    ChildValues_Tau  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);
    
    memset(H, 0, N_CoarseDOF*6*N_DOF*SizeOfDouble);

    memset(LocMatrixG11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG21, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG23, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG32, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG33, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      ChildValue_Tau  = ChildValues_Tau[j];
      w = AbsDetjk[j]*weights[j];
      valu1 = 0.0;
      valu1x = 0.0;
      valu1y = 0.0;
      valu2 = 0.0;
      valu2x = 0.0;
      valu2y = 0.0;
      
        // compute components of uh in j
      for(k=0;k<N_DOF_vel;k++)
      {
        l = DOF_vel[k];
        valu1 += ChildValue[k]*Values1[l];
	valu1x += ChildValueX[k]*Values1[l];
	valu1y += ChildValueY[k]*Values1[l];
        valu2 += ChildValue[k]*Values2[l];
	valu2x += ChildValueX[k]*Values2[l];
	valu2y += ChildValueY[k]*Values2[l];
      }

      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*6*N_DOF+l         ] += val*2*ChildValue_Tau[l]*valu1x;
	  H[k*6*N_DOF+l+   N_DOF] += val*2*ChildValue_Tau[l]*valu1y;
	  H[k*6*N_DOF+l+ 2*N_DOF] += val*ChildValue_Tau[l]*valu2x;
	  H[k*6*N_DOF+l+ 3*N_DOF] += val*ChildValue_Tau[l]*valu1y;
	  H[k*6*N_DOF+l+ 4*N_DOF] += val*2*ChildValue_Tau[l]*valu2x;
	  H[k*6*N_DOF+l+ 5*N_DOF] += val*2*ChildValue_Tau[l]*valu2y;
	  
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG11[k*N_DOF+l] += w*4*ChildValue_Tau[k]*valu1x*ChildValue_Tau[l]*valu1x;
          LocMatrixG12[k*N_DOF+l] += w*4*ChildValue_Tau[k]*valu1y*ChildValue_Tau[l]*valu1y;
          LocMatrixG21[k*N_DOF+l] += w*ChildValue_Tau[k]*valu2x*ChildValue_Tau[l]*valu2x;
	  LocMatrixG23[k*N_DOF+l] += w*ChildValue_Tau[k]*valu1y*ChildValue_Tau[l]*valu1y;
	  LocMatrixG32[k*N_DOF+l] += w*4*ChildValue_Tau[k]*valu2x*ChildValue_Tau[l]*valu2x;
	  LocMatrixG33[k*N_DOF+l] += w*4*ChildValue_Tau[k]*valu2y*ChildValue_Tau[l]*valu2y;
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*6*N_DOF*SizeOfDouble);
    
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 6*N_DOF, 6*N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s23 = 0;
	s32 = 0;
	s33 = 0;
	
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]*H[i1*6*N_DOF+l]*H[i2*6*N_DOF+m];
            s12 += Gsave[i1*N_CoarseDOF+i2]*H[i1*6*N_DOF+l+N_DOF]*H[i2*6*N_DOF+m+N_DOF];
            s21 += Gsave[i1*N_CoarseDOF+i2]*H[i1*6*N_DOF+l+2*N_DOF]*H[i2*6*N_DOF+m+2*N_DOF];
	    s23 += Gsave[i1*N_CoarseDOF+i2]*H[i1*6*N_DOF+l+3*N_DOF]*H[i2*6*N_DOF+m+3*N_DOF];
            s32 += Gsave[i1*N_CoarseDOF+i2]*H[i1*6*N_DOF+l+4*N_DOF]*H[i2*6*N_DOF+m+4*N_DOF];
            s33 += Gsave[i1*N_CoarseDOF+i2]*H[i1*6*N_DOF+l+5*N_DOF]*H[i2*6*N_DOF+m+5*N_DOF];
          }
        LocMatrixG11[l*N_DOF+m] += s11;
        LocMatrixG12[l*N_DOF+m] += s12;
        LocMatrixG21[l*N_DOF+m] += s21;
	LocMatrixG23[l*N_DOF+m] += s23;
        LocMatrixG32[l*N_DOF+m] += s32;
        LocMatrixG33[l*N_DOF+m] += s33;
	
	
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s23 = 0;
	s32 = 0;
	s33 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*6*N_DOF+l] * H[i2*6*N_DOF+m];
          s12 += P[i2*6*N_DOF+l+N_DOF] * H[i2*6*N_DOF+m+N_DOF];
          s21 += P[i2*6*N_DOF+l+2*N_DOF] * H[i2*6*N_DOF+m+2*N_DOF];
	  s23 += P[i2*6*N_DOF+l+3*N_DOF] * H[i2*6*N_DOF+m+3*N_DOF];
          s32 += P[i2*6*N_DOF+l+4*N_DOF] * H[i2*6*N_DOF+m+4*N_DOF];
          s33 += P[i2*6*N_DOF+l+5*N_DOF] * H[i2*6*N_DOF+m+5*N_DOF];
	  
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*6*N_DOF+m] * H[i1*6*N_DOF+l];
          s12 += P[i1*6*N_DOF+m+N_DOF] * H[i1*6*N_DOF+l+N_DOF];
          s21 += P[i1*6*N_DOF+m+2*N_DOF] * H[i1*6*N_DOF+l+2*N_DOF];
          s23 += P[i1*6*N_DOF+m+3*N_DOF] * H[i1*6*N_DOF+l+3*N_DOF];
          s32 += P[i1*6*N_DOF+m+4*N_DOF] * H[i1*6*N_DOF+l+4*N_DOF];
          s33 += P[i1*6*N_DOF+m+5*N_DOF] * H[i1*6*N_DOF+l+5*N_DOF];
        }
        LocMatrixG11[l*N_DOF+m] -= s11;
        LocMatrixG12[l*N_DOF+m] -= s12;
        LocMatrixG21[l*N_DOF+m] -= s21;
	LocMatrixG23[l*N_DOF+m] -= s23;
        LocMatrixG32[l*N_DOF+m] -= s32;
        LocMatrixG33[l*N_DOF+m] -= s33;
      } // end for m
    } // end for l

  //  norm_u = sqrt((valu1*valu1) + (valu2*valu2));
  norm_u = 1.0;
    delta =  lpcoeff*pow(hK,lpexponent)/norm_u;

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] -= delta*LocMatrixG11[l*N_DOF+m];
              EntriesG12[p] -= delta*LocMatrixG12[l*N_DOF+m];
              EntriesG21[p] -= delta*LocMatrixG21[l*N_DOF+m];
	      EntriesG23[p] -= delta*LocMatrixG23[l*N_DOF+m];
              EntriesG32[p] -=delta*LocMatrixG32[l*N_DOF+m];
              EntriesG33[p] -= delta*LocMatrixG33[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
  
  
} 
  
// for conformation stress tensor equation (div of stress tensor)
void AddDivergenceTerm(TSquareMatrix2D **SQMATRICES,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixG11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG21[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG33[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s21, s22, s33;
  int i1, i2;
  double hK, delta;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG22, *EntriesG23, *EntriesG32, *EntriesG33;

  fespace = SQMATRICES[0]->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = SQMATRICES[0]->GetRowPtr();
  KCol = SQMATRICES[0]->GetKCol();
  EntriesG11 = SQMATRICES[0]->GetEntries();
  EntriesG12 = SQMATRICES[1]->GetEntries();
  EntriesG21 = SQMATRICES[2]->GetEntries();
  EntriesG22 = SQMATRICES[3]->GetEntries();
  EntriesG23 = SQMATRICES[4]->GetEntries();
  EntriesG32 = SQMATRICES[5]->GetEntries();
  EntriesG33 = SQMATRICES[6]->GetEntries();


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixG11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG21, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG22, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG33, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG11[k*N_DOF+l] += w*ChildValueX[k]*ChildValueX[l];
          LocMatrixG12[k*N_DOF+l] += w*ChildValueX[k]*ChildValueY[l];
          LocMatrixG21[k*N_DOF+l] += w*ChildValueY[k]*ChildValueX[l];
	  LocMatrixG22[k*N_DOF+l] += w*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
          LocMatrixG33[k*N_DOF+l] += w*ChildValueY[k]*ChildValueY[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ];
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s21 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
                     H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m]        );  
            s33 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF];
          }
        LocMatrixG11[l*N_DOF+m] += s11;
        LocMatrixG12[l*N_DOF+m] += s12;
        LocMatrixG21[l*N_DOF+m] += s21;
	LocMatrixG22[l*N_DOF+m] += s22;
        LocMatrixG33[l*N_DOF+m] += s33;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s12 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m+N_DOF];
	  s21 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m      ];
	  s22 += (P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]);
          s33 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s12 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l      ];
	  s21 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l+N_DOF];
	  s22 += (P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]);
          s33 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrixG11[l*N_DOF+m] -= s11;
        LocMatrixG12[l*N_DOF+m] -= s12;
        LocMatrixG21[l*N_DOF+m] -= s21;
	LocMatrixG22[l*N_DOF+m] -= s22;
        LocMatrixG33[l*N_DOF+m] -= s33;

      } // end for m
    } // end for l

    delta =  lpcoeff*pow(hK,lpexponent);
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] += delta*LocMatrixG11[l*N_DOF+m];
              EntriesG12[p] += delta*LocMatrixG12[l*N_DOF+m];
              EntriesG21[p] += delta*LocMatrixG21[l*N_DOF+m];
	      EntriesG22[p] += delta*LocMatrixG22[l*N_DOF+m];
	      EntriesG23[p] += delta*LocMatrixG12[l*N_DOF+m];
              EntriesG32[p] += delta*LocMatrixG21[l*N_DOF+m];
              EntriesG33[p] += delta*LocMatrixG33[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddDivergenceTerm

// for conformation stress tensor equation (grad of stress tensor)
void AddGradTauTerm(TSquareMatrix2D *G11,TSquareMatrix2D *G22,
                       TSquareMatrix2D *G33, double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixG11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG33[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s22, s33;
  int i1, i2;
  double hK, delta;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG22, *EntriesG23, *EntriesG32, *EntriesG33;

  fespace = G11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = G11->GetRowPtr();
  KCol = G11->GetKCol();
  EntriesG11 = G11->GetEntries();
  EntriesG22 = G22->GetEntries();
  EntriesG33 = G33->GetEntries();


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixG11, 0, N_DOF*N_DOF*SizeOfDouble);
//     memset(LocMatrixG22, 0, N_DOF*N_DOF*SizeOfDouble);
//     memset(LocMatrixG33, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG11[k*N_DOF+l] += w*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
// 	  LocMatrixG22[k*N_DOF+l] += w*2.0*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
//           LocMatrixG33[k*N_DOF+l] += w*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
// 	s22 = 0;
// 	s33 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
                     H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m] );  
//             s22 += Gsave[i1*N_CoarseDOF+i2]
//                   *2.0*(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
//                      H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m] );  
//             s33 += Gsave[i1*N_CoarseDOF+i2]
//                   *(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
//                      H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m] );  
          }
        LocMatrixG11[l*N_DOF+m] += s11;
// 	LocMatrixG22[l*N_DOF+m] += s22;
//         LocMatrixG33[l*N_DOF+m] += s33;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
// 	s22 = 0;
// 	s33 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += (P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]);
// 	  s22 += 2.0*((P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]));
//           s33 += (P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]);
        
	}
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += (P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]);
// 	  s22 += 2.0*((P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]));
//           s33 += (P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]);
        }
        LocMatrixG11[l*N_DOF+m] -= s11;
// 	LocMatrixG22[l*N_DOF+m] -= s22;
//         LocMatrixG33[l*N_DOF+m] -= s33;

      } // end for m
    } // end for l

    delta =  lpcoeff*pow(hK,lpexponent);
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] += delta*LocMatrixG11[l*N_DOF+m];
	      EntriesG22[p] += delta*2.0*LocMatrixG11[l*N_DOF+m];
              EntriesG33[p] += delta*LocMatrixG11[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddGradTauTerm


// for conformation stress tensor equation (grad and div of stress tensor)
void AddTauTerm(TSquareMatrix2D *G11, TSquareMatrix2D *G12, TSquareMatrix2D *G21,
		TSquareMatrix2D *G22, TSquareMatrix2D *G23, TSquareMatrix2D *G32,
                       TSquareMatrix2D *G33, double lpcoeff_grad, double lpcoeff_div, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixG11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG21[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG33[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG11_Grad[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s21, s22, s33, s11_Grad;
  int i1, i2;
  double hK, delta1, delta2;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG22, *EntriesG23, *EntriesG32, *EntriesG33;

  fespace = G11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = G11->GetRowPtr();
  KCol = G11->GetKCol();
  EntriesG11 = G11->GetEntries();
  EntriesG12 = G12->GetEntries();
  EntriesG21 = G21->GetEntries();
  EntriesG22 = G22->GetEntries();
  EntriesG23 = G23->GetEntries();
  EntriesG32 = G32->GetEntries();
  EntriesG33 = G33->GetEntries();


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixG11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG21, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG22, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG33, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG11_Grad, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG11[k*N_DOF+l] += w*ChildValueX[k]*ChildValueX[l];
          LocMatrixG12[k*N_DOF+l] += w*ChildValueX[k]*ChildValueY[l];
          LocMatrixG21[k*N_DOF+l] += w*ChildValueY[k]*ChildValueX[l];
	  LocMatrixG22[k*N_DOF+l] += w*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
          LocMatrixG33[k*N_DOF+l] += w*ChildValueY[k]*ChildValueY[l];
	  LocMatrixG11_Grad[k*N_DOF+l] += w*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
	}
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
	s11_Grad = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ];
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s21 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
                     H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m]        );  
            s33 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF];
	    s11_Grad += Gsave[i1*N_CoarseDOF+i2]
                  *(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
                     H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m] );  
          }
        LocMatrixG11[l*N_DOF+m] += s11;
        LocMatrixG12[l*N_DOF+m] += s12;
        LocMatrixG21[l*N_DOF+m] += s21;
	LocMatrixG22[l*N_DOF+m] += s22;
        LocMatrixG33[l*N_DOF+m] += s33;
	LocMatrixG11_Grad[l*N_DOF+m] += s11_Grad;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
	s11_Grad = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s12 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m+N_DOF];
	  s21 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m      ];
	  s22 += (P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]);
          s33 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
	  s11_Grad += (P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]);
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s12 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l      ];
	  s21 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l+N_DOF];
	  s22 += (P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]);
          s33 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
	  s11_Grad += (P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]);
        }
        LocMatrixG11[l*N_DOF+m] -= s11;
        LocMatrixG12[l*N_DOF+m] -= s12;
        LocMatrixG21[l*N_DOF+m] -= s21;
	LocMatrixG22[l*N_DOF+m] -= s22;
        LocMatrixG33[l*N_DOF+m] -= s33;
        LocMatrixG11_Grad[l*N_DOF+m] -= s11_Grad;
      } // end for m
    } // end for l

    delta1 =  lpcoeff_grad*pow(hK,lpexponent);
    delta2 =  lpcoeff_div*pow(hK,lpexponent);
     
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] += delta2*LocMatrixG11[l*N_DOF+m] + delta1*LocMatrixG11_Grad[l*N_DOF+m];
              EntriesG12[p] += delta2*LocMatrixG12[l*N_DOF+m];
              EntriesG21[p] += delta2*LocMatrixG21[l*N_DOF+m];
	      EntriesG22[p] += delta2*LocMatrixG22[l*N_DOF+m]+ delta1*2.0*LocMatrixG11_Grad[l*N_DOF+m];
	      EntriesG23[p] += delta2*LocMatrixG12[l*N_DOF+m];
              EntriesG32[p] += delta2*LocMatrixG21[l*N_DOF+m];
              EntriesG33[p] += delta2*LocMatrixG33[l*N_DOF+m]+ delta1*LocMatrixG11_Grad[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddTauTerm


// for conformation stress tensor equation (grad and div of stress tensor)
void AddTauTerm_2PhaseOrImpDropFlow(TSquareMatrix2D *G11, TSquareMatrix2D *G12, TSquareMatrix2D *G21,
		TSquareMatrix2D *G22, TSquareMatrix2D *G23, TSquareMatrix2D *G32, TSquareMatrix2D *G33, 
		       double lpcoeff_grad, double lpcoeff_div, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF, Phase_No;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixG11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG21[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG33[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG11_Grad[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s21, s22, s33, s11_Grad;
  int i1, i2;
  double hK, delta1, delta2;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG22, *EntriesG23, *EntriesG32, *EntriesG33;

  fespace = G11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = G11->GetRowPtr();
  KCol = G11->GetKCol();
  EntriesG11 = G11->GetEntries();
  EntriesG12 = G12->GetEntries();
  EntriesG21 = G21->GetEntries();
  EntriesG22 = G22->GetEntries();
  EntriesG23 = G23->GetEntries();
  EntriesG32 = G32->GetEntries();
  EntriesG33 = G33->GetEntries();


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    
   if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
   {
    Phase_No = cell->GetRegionID();
    if(Phase_No == 0)
    { 
      if(TDatabase::ParamDB->PHASE1_TYPE == 1)
      {
	lpcoeff_grad = 0.0;
	lpcoeff_div = 0.0; 
      }
    } 
      else if(Phase_No == 1)
     {
      if(TDatabase::ParamDB->PHASE2_TYPE == 1)
      {
	lpcoeff_grad = 0.0;
	lpcoeff_div = 0.0; 
      }
     }
      else
    {
     cout<<"Invalid Phase No. !!!!!\n";
     exit(4711);
    }
   }
  else if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
 {
    if(TDatabase::ParamDB->PHASE1_TYPE == 1)
      {
	lpcoeff_grad = 0.0;
	lpcoeff_div = 0.0; 
      }
 }
  else
{
    OutPut("Implemented only for two-phase and free surface flows \n");
    OutPut("Change TWO_PHASE_FLOW or FREE_SURFACE_FLOW  to 1 in dat file !!!!! " << endl);
    exit(4711);
  
}
    
    
    
    if(lpcoeff_grad == 0.0 && lpcoeff_div == 0.0)
    {
     continue; 
    }
    else
    {
    
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixG11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG21, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG22, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG33, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG11_Grad, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG11[k*N_DOF+l] += w*ChildValueX[k]*ChildValueX[l];
          LocMatrixG12[k*N_DOF+l] += w*ChildValueX[k]*ChildValueY[l];
          LocMatrixG21[k*N_DOF+l] += w*ChildValueY[k]*ChildValueX[l];
	  LocMatrixG22[k*N_DOF+l] += w*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
          LocMatrixG33[k*N_DOF+l] += w*ChildValueY[k]*ChildValueY[l];
	  LocMatrixG11_Grad[k*N_DOF+l] += w*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
	}
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
	s11_Grad = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ];
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s21 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
                     H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m]        );  
            s33 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF];
	    s11_Grad += Gsave[i1*N_CoarseDOF+i2]
                  *(H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF] + 
                     H[i1*2*N_DOF+l]*H[i2*2*N_DOF+m] );  
          }
        LocMatrixG11[l*N_DOF+m] += s11;
        LocMatrixG12[l*N_DOF+m] += s12;
        LocMatrixG21[l*N_DOF+m] += s21;
	LocMatrixG22[l*N_DOF+m] += s22;
        LocMatrixG33[l*N_DOF+m] += s33;
	LocMatrixG11_Grad[l*N_DOF+m] += s11_Grad;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
	s11_Grad = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s12 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m+N_DOF];
	  s21 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m      ];
	  s22 += (P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]);
          s33 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
	  s11_Grad += (P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ]) + (P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF]);
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s12 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l      ];
	  s21 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l+N_DOF];
	  s22 += (P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]);
          s33 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
	  s11_Grad += (P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ])+(P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF]);
        }
        LocMatrixG11[l*N_DOF+m] -= s11;
        LocMatrixG12[l*N_DOF+m] -= s12;
        LocMatrixG21[l*N_DOF+m] -= s21;
	LocMatrixG22[l*N_DOF+m] -= s22;
        LocMatrixG33[l*N_DOF+m] -= s33;
        LocMatrixG11_Grad[l*N_DOF+m] -= s11_Grad;
      } // end for m
    } // end for l

    delta1 =  lpcoeff_grad*pow(hK,lpexponent);
    delta2 =  lpcoeff_div*pow(hK,lpexponent);
     
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] += delta2*LocMatrixG11[l*N_DOF+m] + delta1*LocMatrixG11_Grad[l*N_DOF+m];
              EntriesG12[p] += delta2*LocMatrixG12[l*N_DOF+m];
              EntriesG21[p] += delta2*LocMatrixG21[l*N_DOF+m];
	      EntriesG22[p] += delta2*LocMatrixG22[l*N_DOF+m]+ delta1*2.0*LocMatrixG11_Grad[l*N_DOF+m];
	      EntriesG23[p] += delta2*LocMatrixG12[l*N_DOF+m];
              EntriesG32[p] += delta2*LocMatrixG21[l*N_DOF+m];
              EntriesG33[p] += delta2*LocMatrixG33[l*N_DOF+m]+ delta1*LocMatrixG11_Grad[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
    
  } // if(lpcoeff_grad == 0.0 && lpcoeff_div == 0.0)
    
  } // endfor i
} // AddTauTerm_2PhaseFlow



// for conformation stress tensor equation in 3D-Axis-symmetric form (grad and div of stress tensor)
void AddTauTerm_2PhaseOrImpDropFlow_3DAxial(TSquareMatrix2D *G11, TSquareMatrix2D *G12, TSquareMatrix2D *G21,
		TSquareMatrix2D *G22, TSquareMatrix2D *G23, TSquareMatrix2D *G32, TSquareMatrix2D *G33, 
		       double lpcoeff_grad, double lpcoeff_div, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF, Phase_No;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights, r;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave1[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave2[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave3[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[3*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[6*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixG11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG21[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG33[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG11_Grad[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG22_Grad[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixG33_Grad[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  
  double s11, s12, s21, s22, s33, s11_Grad, s22_Grad, s33_Grad;
  int i1, i2;
  double hK, delta1, delta2;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG22, *EntriesG23, *EntriesG32, *EntriesG33;

  fespace = G11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = G11->GetRowPtr();
  KCol = G11->GetKCol();
  EntriesG11 = G11->GetEntries();
  EntriesG12 = G12->GetEntries();
  EntriesG21 = G21->GetEntries();
  EntriesG22 = G22->GetEntries();
  EntriesG23 = G23->GetEntries();
  EntriesG32 = G32->GetEntries();
  EntriesG33 = G33->GetEntries();


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    
   if(TDatabase::ParamDB->TWO_PHASE_FLOW==1)
   {
    Phase_No = cell->GetRegionID();
    if(Phase_No == 0)
    { 
      if(TDatabase::ParamDB->PHASE1_TYPE == 1)
      {
	lpcoeff_grad = 0.0;
	lpcoeff_div = 0.0; 
      }
    } 
      else if(Phase_No == 1)
     {
      if(TDatabase::ParamDB->PHASE2_TYPE == 1)
      {
	lpcoeff_grad = 0.0;
	lpcoeff_div = 0.0; 
      }
     }
      else
    {
     cout<<"Invalid Phase No. !!!!!\n";
     exit(4711);
    }
   }
  else if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
 {
    if(TDatabase::ParamDB->PHASE1_TYPE == 1)
      {
	lpcoeff_grad = 0.0;
	lpcoeff_div = 0.0; 
      }
 }
  else
{
    OutPut("Implemented only for two-phase and free surface flows \n");
    OutPut("Change TWO_PHASE_FLOW or FREE_SURFACE_FLOW  to 1 in dat file !!!!! " << endl);
    exit(4711);
  
}
     
    if(lpcoeff_grad == 0.0 && lpcoeff_div == 0.0)
    {
     continue; 
    }
    else
    {
    
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    memset(Gsave1, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    memset(Gsave2, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    memset(Gsave3, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      r  = fabs(X[j]);
    if(r<1e-12)
    {
     OutPut("check Local Projection Stabilization  r value zero !!!!! "<< X[j] <<endl);
     OutPut("Quad formula: Change all integral points as internal points"<<endl);
     exit(0);
    }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l]*r;
	  Gsave1[k*N_CoarseDOF+l] += val*CoarseValue[l]*r;
	  Gsave2[k*N_CoarseDOF+l] += val*CoarseValue[l]/r;
	  Gsave3[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;


    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*3*N_DOF*SizeOfDouble);
    memset(P, 0, N_CoarseDOF*6*N_DOF*SizeOfDouble);

    memset(LocMatrixG11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG21, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG22, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG33, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG11_Grad, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG22_Grad, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixG33_Grad, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues [j];
      w = AbsDetjk[j]*weights[j];
      r  = fabs(X[j]);
    if(r<1e-12)
    {
     OutPut("check Local Projection Stabilization  r value zero !!!!! "<< X[j] <<endl);
     OutPut("Quad formula: Change all integral points as internal points"<<endl);
     exit(0);
    }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*3*N_DOF+l      ] += val*ChildValueX[l]*r;
          H[k*3*N_DOF+l+N_DOF] += val*ChildValueY[l]*r;
	  H[k*3*N_DOF+l+2*N_DOF] += val*ChildValue[l]*r;
	  
	  P[k*6*N_DOF+l      ] += val*r*ChildValueX[l];
          P[k*6*N_DOF+l+N_DOF] += val*r*ChildValueY[l];
	  P[k*6*N_DOF+l+2*N_DOF] += val*ChildValue[l]/r; 
	  P[k*6*N_DOF+l+3*N_DOF] += val*ChildValueX[l];
	  P[k*6*N_DOF+l+4*N_DOF] += val*ChildValueY[l];
	  P[k*6*N_DOF+l+5*N_DOF] += val*ChildValue[l];
	  
        } // end for l
      } // end for k
      

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixG11[k*N_DOF+l] += w*r*(ChildValueX[k] + ChildValue[k]/r)*(ChildValueX[l] + ChildValue[l]/r);
          LocMatrixG12[k*N_DOF+l] += w*r*(ChildValueX[k] + ChildValue[k]/r)*ChildValueY[l];
          LocMatrixG21[k*N_DOF+l] += w*r*ChildValueY[k]*(ChildValueX[l] + ChildValue[l]/r);
	  LocMatrixG22[k*N_DOF+l] += w*r*( (ChildValueX[k] + ChildValue[k]/r)*(ChildValueX[l] + ChildValue[l]/r) + ChildValueY[k]*ChildValueY[l]);
          LocMatrixG33[k*N_DOF+l] += w*r*ChildValueY[k]*ChildValueY[l];
	  LocMatrixG11_Grad[k*N_DOF+l] += w*r*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l] + 2.0*ChildValue[k]*ChildValue[l]/(r*r));
	  LocMatrixG22_Grad[k*N_DOF+l] += w*r*2.0*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l] + ChildValue[k]*ChildValue[l]/(r*r));
	  LocMatrixG33_Grad[k*N_DOF+l] += w*r*(ChildValueX[k]*ChildValueX[l] + ChildValueY[k]*ChildValueY[l]);
	}
      }
    } // end for j
//     memcpy(P, H, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3*N_DOF, 3*N_DOF);

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
	s11_Grad = 0;
	s22_Grad = 0;
	s33_Grad = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m];
	    s11 += Gsave3[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m+2*N_DOF];
	    s11 += Gsave3[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m];
	    s11 += Gsave2[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
	    
            s12 += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m+N_DOF];
	    s12 += Gsave3[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+N_DOF];
	    
            s21 += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m]; 
	    s21 += Gsave3[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
	    
            s22 += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF];
	    s22 += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m];
	    s22 += Gsave3[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m];
	    s22 += Gsave3[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m+2*N_DOF];
	    s22 += Gsave2[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];

            s33 += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF];
	    
	    s11_Grad += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF];
	    s11_Grad += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m]; 
            s11_Grad += Gsave2[i1*N_CoarseDOF+i2]*2.0*H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
	    
            s22_Grad += Gsave1[i1*N_CoarseDOF+i2]*2.0*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF];
	    s22_Grad += Gsave1[i1*N_CoarseDOF+i2]*2.0*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m];
            s22_Grad += Gsave2[i1*N_CoarseDOF+i2]*2.0*H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF]; 
	    
	    s33_Grad += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF];
            s33_Grad += Gsave1[i1*N_CoarseDOF+i2]*H[i1*3*N_DOF+l]*H[i2*3*N_DOF+m]; 
		  
          }
        LocMatrixG11[l*N_DOF+m] += s11;
        LocMatrixG12[l*N_DOF+m] += s12;
        LocMatrixG21[l*N_DOF+m] += s21;
	LocMatrixG22[l*N_DOF+m] += s22;
        LocMatrixG33[l*N_DOF+m] += s33;
	LocMatrixG11_Grad[l*N_DOF+m] += s11_Grad;
	LocMatrixG22_Grad[l*N_DOF+m] += s22_Grad;
	LocMatrixG33_Grad[l*N_DOF+m] += s33_Grad;
	
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s21 = 0;
	s22 = 0;
	s33 = 0;
	s11_Grad = 0;
        s22_Grad = 0;
	s33_Grad = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*6*N_DOF+l]*H[i2*3*N_DOF+m];
	  s11 += P[i2*6*N_DOF+l+3*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
          s11 += P[i2*6*N_DOF+l+5*N_DOF]*H[i2*3*N_DOF+m];
	  s11 += P[i2*6*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
          
	  s12 += P[i2*6*N_DOF+l] * H[i2*3*N_DOF+m+N_DOF];
	  s12 += P[i2*6*N_DOF+l+5*N_DOF]* H[i2*3*N_DOF+m+N_DOF];
	  
	  s21 += P[i2*6*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m]; 
          s21 += P[i2*6*N_DOF+l+4*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
	  
	  s22 += P[i2*6*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF];
	  s22 += P[i2*6*N_DOF+l]*H[i2*3*N_DOF+m];
	  s22 += P[i2*6*N_DOF+l+5*N_DOF]*H[i2*3*N_DOF+m];
	  s22 += P[i2*6*N_DOF+l+3*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
	  s22 += P[i2*6*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF];
	  
	  s33 += P[i2*6*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF];
	  
	  s11_Grad += (P[i2*6*N_DOF+l] * H[i2*3*N_DOF+m]) + (P[i2*6*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF]) + 2.0*(P[i2*6*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF]);
          s22_Grad += 2.0*((P[i2*6*N_DOF+l] * H[i2*3*N_DOF+m]) + (P[i2*6*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF]) + (P[i2*6*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF]));
	  s33_Grad += (P[i2*6*N_DOF+l] * H[i2*3*N_DOF+m]) + (P[i2*6*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF]);
	}
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*6*N_DOF+m]*H[i1*3*N_DOF+l];
	  s11 += P[i1*6*N_DOF+m+5*N_DOF]*H[i1*3*N_DOF+l];
	  s11 += P[i1*6*N_DOF+m+3*N_DOF]*H[i1*3*N_DOF+l+2*N_DOF];
	  s11 += P[i1*6*N_DOF+m+2*N_DOF]*H[i1*3*N_DOF+l+2*N_DOF];
	  
          s12 += P[i1*6*N_DOF+m+N_DOF]*H[i1*3*N_DOF+l] ;
	  s12 += P[i1*6*N_DOF+m+4*N_DOF]*H[i1*3*N_DOF+l+2*N_DOF];
	  
	  s21 += (P[i1*6*N_DOF+m] + P[i1*6*N_DOF+m+5*N_DOF] )* H[i1*3*N_DOF+l+N_DOF];
	  
	  s22 += P[i1*6*N_DOF+m+N_DOF]*H[i1*3*N_DOF+l+N_DOF];
	  s22 += P[i1*6*N_DOF+m]*H[i1*3*N_DOF+l];
	  s22 += P[i1*6*N_DOF+m+3*N_DOF]*H[i1*3*N_DOF+l+2*N_DOF];
	  s22 += P[i1*6*N_DOF+m+5*N_DOF]*H[i1*3*N_DOF+l];
	  s22 += P[i1*6*N_DOF+m+2*N_DOF]*H[i1*3*N_DOF+l+2*N_DOF];
	  
	  s33 += P[i1*6*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF];
	  
	  s11_Grad += (P[i1*6*N_DOF+m] * H[i1*3*N_DOF+l])+(P[i1*6*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF]) + 2.0*(P[i1*6*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF]);
	  s22_Grad += 2.0*((P[i1*6*N_DOF+m] * H[i1*3*N_DOF+l])+(P[i1*6*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF]) + (P[i1*6*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF]));
	  s33_Grad += (P[i1*6*N_DOF+m] * H[i1*3*N_DOF+l])+(P[i1*6*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF]);
        }
        LocMatrixG11[l*N_DOF+m] -= s11;
        LocMatrixG12[l*N_DOF+m] -= s12;
        LocMatrixG21[l*N_DOF+m] -= s21;
	LocMatrixG22[l*N_DOF+m] -= s22;
        LocMatrixG33[l*N_DOF+m] -= s33;
        LocMatrixG11_Grad[l*N_DOF+m] -= s11_Grad;
	LocMatrixG22_Grad[l*N_DOF+m] -= s22_Grad;
	LocMatrixG33_Grad[l*N_DOF+m] -= s33_Grad;
      } // end for m
    } // end for l

    delta1 =  lpcoeff_grad*pow(hK,lpexponent);
    delta2 =  lpcoeff_div*pow(hK,lpexponent);
     
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesG11[p] += delta2*LocMatrixG11[l*N_DOF+m] + delta1*LocMatrixG11_Grad[l*N_DOF+m];
              EntriesG12[p] += delta2*LocMatrixG12[l*N_DOF+m];
              EntriesG21[p] += delta2*LocMatrixG21[l*N_DOF+m];
	      EntriesG22[p] += delta2*LocMatrixG22[l*N_DOF+m]+ delta1*LocMatrixG22_Grad[l*N_DOF+m];
	      EntriesG23[p] += delta2*LocMatrixG12[l*N_DOF+m];
              EntriesG32[p] += delta2*LocMatrixG21[l*N_DOF+m];
              EntriesG33[p] += delta2*LocMatrixG33[l*N_DOF+m]+ delta1*LocMatrixG33_Grad[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
    
  } // if(lpcoeff_grad == 0.0 && lpcoeff_div == 0.0)
    
  } // endfor i
} // AddTauTerm_2PhaseOrImpDropFlow_3DAxial



//**************************************************************
//  UltraLocalProjectionStreamlinePLaplacian
//  ultra LPS for scalar equation
//  projection of streamline derivative
//  p-Laplacian term can be added
//**************************************************************

void UltraLocalProjectionStreamlinePLaplacian(TSquareMatrix2D* A, 
                TFEFunction2D *uh,
                CoeffFct2D *Coeff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double rhs_l2_proj[MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s, loc_proj_L2;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end, ij, N_Edges;
  int *RowPtr, *KCol;
  double *Entries;
  double *Values, Values1, Values2;
  double *coeffs, *params, sx, sy, crosswind_uh[MaxN_QuadPoints_2D], loc_proj;
  double BValue[MaxN_BaseFunctions2D];
  int OrderDiff;
  double lpcoeff, lpexponent, stab_coeff, norm_b, lpcoeff_crosswind, lpexponent_crosswind;
  double vall[3];

  if (TDatabase::ParamDB->SC_VERBOSE>1)
    OutPut("LPS streamline" << endl);

  coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * SizeOfDouble);

  lpcoeff = TDatabase::ParamDB->LP_STREAMLINE_COEFF;
  lpexponent = TDatabase::ParamDB->LP_STREAMLINE_EXPONENT;
  OrderDiff = TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE;
  lpcoeff_crosswind = TDatabase::ParamDB->LP_CROSSWIND_COEFF;
  lpexponent_crosswind = TDatabase::ParamDB->LP_CROSSWIND_EXPONENT;
  
  // get fespace and matrices
  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();

  // get values of fe function if available
  if (uh != NULL)
  {
      Values = uh->GetValues();     
  }
  else
  {
      if (TDatabase::ParamDB->LP_CROSSWIND)
      {
    OutPut("for crosswind LPS the current finite element function has to be given to the routine"
     << endl);
    exit(4711);
      }
  }

  // get collection
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // get arrays for dofs  
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    // number of edges
    N_Edges = cell->GetN_Edges();
    // get diameter of the cell
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
  case 0: // diameter
      hK = cell->GetDiameter();
      break;
  case 1: // with reference map
      hK = cell->GetLengthWithReferenceMap();
      break;
  case 2: // shortest edge
      hK = cell->GetShortestEdge();
      break;
  case 3: // measure
      hK = cell->GetMeasure();
      hK = pow(hK,1.0/3.0);
      break;
  case 4: // mesh cell in convection direction
            TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
            sx = sy = 0;
      for (ij=0;ij<N_Edges;ij++)
      {
    TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
    sx += TDatabase::ParamDB->INTERNAL_VERTEX_X[ij];
    TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
    sy += TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij];
      }
      if (N_Edges==3)
    TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      // center of mesh cell
      sx /= N_Edges;
      sy /= N_Edges;
      hK = cell->GetDiameter(); 
      // get coefficients in center of mesh cell 
      Coeff(1, &sx ,&sy, &params, &coeffs);
      hK = Mesh_size_in_convection_direction(hK, coeffs[1], coeffs[2]); 
      break;
  default: // diameter
      hK = cell->GetDiameter();
      break;
    }
    
    // get finite element in the cell
    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);
    // get basis functions of the finite element
    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    // compute order of the local coarse space
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // get type of the mesh cell
    shapetype = cell->GetType();

    // determine finite element of the local coarse space
    switch(shapetype)
    {
      // regularly refined quadrilateral
      // discontinuous finite elements  
      case Quadrangle:
      case Parallelogram:
      case Rectangle:
        switch(CoarseOrder)
        {
          case 0:
            UsedElements[0] = C_Q0_2D_Q_M;
          break;

          case 1:
            UsedElements[0] = D_P1_2D_Q_M;
          break;

          case 2:
            UsedElements[0] = D_P2_2D_Q_M;
          break;

          case 3:
            UsedElements[0] = D_P3_2D_Q_M;
          break;

          case 4:
            UsedElements[0] = D_P4_2D_Q_M;
          break;

          default:
            OutPut("Projection space is defined up to order 4" << endl);
            exit(-1);
        } // end switch CoarseOrder
      break; // end regularly refined quadrilateral

      case Triangle:
        switch(CoarseOrder)
        {
          case 0:
            UsedElements[0] = C_P0_2D_T_A;
          break;

          case 1:
            UsedElements[0] = D_P1_2D_T_A;
          break;

          case 2:
            UsedElements[0] = D_P2_2D_T_A;
          break;

          case 3:
            UsedElements[0] = D_P3_2D_T_A;
          break;

          case 4:
            UsedElements[0] = D_P4_2D_T_A;
          break;

          default:
            OutPut("Projection space is defined up to order 4" << endl);
            exit(-1);
        }
      break;
      default: 
    OutPut("No coarse finite elements for mesh cell type "
     << shapetype << " implemented !!!" << endl);
    exit(4711);
    } // end switch reftype

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // number of dof in the local coarse space
    N_CoarseDOF = CoarseBF->GetDimension();
    // get function values for the basis functions of the coarse space
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);
    // initialize array G, stores mass matrix of coarse space
    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
      // values of the coarse basis functions in the quad points
      CoarseValue = CoarseValues[j];
      // factor for numerical quadrature
      w = AbsDetjk[j]*weights[j];
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
    // first factor of integrand
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
      // update integral 
           G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;
    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    // degrees of freedom of the fine space 
    DOF = GlobalNumbers + BeginIndex[i];
    // this should be the same as CoarseValues ???    
    //PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);
    PCValues = CoarseValues;
    // get derivatives of fe functions of approximation space
    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);
    // initialize array H, holds mixed products of coarse basis 
    // functions and streamline derivatives of fine basis functions
    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);
    // initialize array LocMatrix, holds products of derivatives
    // of fine basis functions
    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
      // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = AbsDetjk[j]*weights[j];
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      // compute streamline derivative
      for(k=0;k<N_DOF;k++)
      {
    BValue[k] = coeffs[1]*ChildValueX[k] + coeffs[2]*ChildValueY[k];
    //BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
      }
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
  // compute products of coarse basis function and
  // streamline derivatives of fine basis functions
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for kcoeffs[2]

      // fine-fine couplings
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
      LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    // save H for later use
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);
     // solve G * X = H, solution X stored on H
    // the right hand side and the solution are stored column wise
    // right hand side: a finecoeffs[2] fct. tested with all coarse fcts. 
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    // update LocMatrix
    // proj-proj coupling (coarse-coarse coupling)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
      // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling (fine-coarse couplings)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
      // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l
    // compute stabilzation coefficient
    switch(TDatabase::ParamDB->LP_COEFF_TYPE)
    {
  case 0:
      stab_coeff = lpcoeff*pow(hK,lpexponent);
      break;
  case 1:
      // h^2/epsilon
      stab_coeff = hK*hK/coeffs[0];
      norm_b = sqrt(coeffs[1]*coeffs[1]+coeffs[2]*coeffs[2]);
      if (norm_b > 1e-10)
      {
    // h/|b|
    val = hK/norm_b;
    if (val < stab_coeff)
        stab_coeff = val;
      }
      stab_coeff *= lpcoeff;
      break;
  default:
      stab_coeff = lpcoeff*pow(hK,lpexponent);
      break;
    }
 
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
    // parameter is lpcoeff*pow(hK,lpexponent)
               Entries[p] += stab_coeff*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l

    //***************************************************
    // LPS with SOLD term
    //***************************************************
    if (TDatabase::ParamDB->LP_CROSSWIND)
    {
  // recover mass matrix of coarse space
  memcpy(G, Gsave, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
  memset(rhs_l2_proj, 0, N_CoarseDOF*SizeOfDouble);
  // compute L2 projection of crosswind derivative
  // this gives just an array (one value for each quad point)
  // mass matrix of coarse space already computed: G
  // compute right hand side
  // loop over the quadrature points
  for(j=0;j<N_Points;j++)
  {
      // get all values in the quad point j 
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // values of the coarse basis functions in the quad points
      CoarseValue = CoarseValues[j];
      // factor for numerical quadrature
      w = AbsDetjk[j]*weights[j];
      // compute crosswind derivative in this point
      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
    l = DOF[k];
    val = Values[l];
    valx += ChildValueX[k]*val;
    valy += ChildValueY[k]*val;
      }
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      
      // calculate crosswind derivative in this quadrature point
      crosswind_uh[j] = -coeffs[2] * valx + coeffs[1] * valy;
      //OutPut(crosswind_uh[j] << " ");
      w *= crosswind_uh[j];
      // loop over the test functions of the coarse space
      // right hand side of linear system for local projection
      for(k=0;k<N_CoarseDOF;k++)
      {
    rhs_l2_proj[k] += w*CoarseValue[k];
      } // end for k
      //OutPut(rhs_l2_proj[0] << " :: ");
  } // end for j  
  // solution is on rhs_l2_proj
  // these are coefficients wrt to coarse space basis
  //OutPut("G " << G[0] << " " <<  rhs_l2_proj[0] << endl);
  SolveLinearSystemLapack(G, rhs_l2_proj, N_CoarseDOF, N_CoarseDOF);
  //OutPut("G " << G[0] << " " <<  rhs_l2_proj[0] << endl);
  // recover mass matrix of coarse space
  memcpy(G, Gsave, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
 
  // initialize array H, holds mixed products of coarse basis 
  // functions and crosswind derivatives of fine basis functions
  memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);
  // initialize array LocMatrix, holds products of derivatives
  // of fine basis functions
  memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

  // loop over the quadrature points
  // this computes the projections of the crosswind derivatives 
  // for each basis function of the fine space
  for(j=0;j<N_Points;j++)
  {
      // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = AbsDetjk[j]*weights[j];
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      // compute crosswind derivative
      for(k=0;k<N_DOF;k++)
      {
    BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
        }
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
    val = w*PCValue[k];
    // compute products of coarse basis function and
    // crosswind derivatives of fine basis functions
    for(l=0;l<N_DOF;l++)
    {
        H[k*N_DOF+l] += val*BValue[l];
    } // end for l
      } // end for k
  } // end for j

  // solve G * X = H, solution X stored on H
  // the right hand side and the solution are stored column wise
  // right hand side: a fine fct. tested with all coarse fcts. 
  SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

  // compute L^2 norm of fluctuations 
  if (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE >= 200)
  {
      loc_proj_L2 = 0.0;
      for(j=0;j<N_Points;j++)
      {
    // get all values in the quad point j 
    PCValue = PCValues[j];
    ChildValueX = ChildValuesX[j];
    ChildValueY = ChildValuesY[j];
    // factor in the integral
    w = AbsDetjk[j]*weights[j];
    // values of convection in this quadrature point
    Coeff(1, &X[j] , &Y[j], &params, &coeffs);
    // compute crosswind derivative
    for(k=0;k<N_DOF;k++)
    {
        BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
    }
    // compute projection of crosswind derivative for basis functions
    for(k=0;k<N_DOF;k++)
    {
        P[k] = 0;
        for (l=0; l < N_CoarseDOF; l++)
      P[k] += H[l*N_DOF+k] * PCValue[l];
    }
    // calculate projection of crosswind derivative discrete uh in this quadrature point
    val = 0.0;
    for(k=0;k<N_CoarseDOF;k++)
    {
        val += rhs_l2_proj[k] * PCValue[k]; 
    }
    // local projection of crosswind derivative of current solution in this quadrature point
    // this is the coefficient
    loc_proj = val - crosswind_uh[j];
    // update norm
    loc_proj_L2 += w * loc_proj * loc_proj;
      } // end for j
      // L2 norm and scale such that independent of h
      loc_proj_L2 = sqrt(loc_proj_L2/cell->GetMeasure());
  }

  // loop over the quadrature points
  // this computes the LPS term
  // for each basis function of the fine space
  for(j=0;j<N_Points;j++)
  {
      // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = AbsDetjk[j]*weights[j];
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      // compute crosswind derivative
      for(k=0;k<N_DOF;k++)
      {
    BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
      }
      // compute projection of crosswind derivative for basis functions
      for(k=0;k<N_DOF;k++)
      {
    P[k] = 0;
    for (l=0; l < N_CoarseDOF; l++)
        P[k] += H[l*N_DOF+k] * PCValue[l];
      }
      // calculate projection of crosswind derivative discrete uh in this quadrature point
      if (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE < 100)
      {
    val = 0.0;
    for(k=0;k<N_CoarseDOF;k++)
    {
        val += rhs_l2_proj[k] * PCValue[k]; 
    }
    // local projection of crosswind derivative of current solution in this quadrature point
    // this is the coefficient
    loc_proj = fabs(val - crosswind_uh[j]);
      }
      // linear case
      if ((TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE >= 100) && (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE <200))
    loc_proj = 1.0;
      // local average
      if (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE >= 200)
      {
    loc_proj = loc_proj_L2;
      }

      // update factor in integral
      w *= loc_proj;
      // update local matrix
      for(k=0;k<N_DOF;k++)
      {
    for(l=0;l<N_DOF;l++)
    {
        // fine-fine
        val = BValue[k]*BValue[l];
        // coarse-fine
        val -= BValue[k]*P[l];
        val -= BValue[l]*P[k];
        // coarse-coarse
        val += P[l]*P[k];       
        LocMatrix[k*N_DOF+l] += w*val;
    }
      }
  } // end for j


  // compute stabilzation coefficient
  switch(TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE)
  {
      case 0:
      case 100:
      case 200:
    stab_coeff = lpcoeff_crosswind*pow(hK,lpexponent_crosswind);
    break;
      case 1:
      case 101:
      case 201:
    // norm squared of convection
    norm_b = coeffs[1]*coeffs[1]+coeffs[2]*coeffs[2];
    if (norm_b > 1e-20)
    {
        // h/|b|
        stab_coeff = lpcoeff_crosswind*hK/norm_b;
    }
    else
        stab_coeff = 0.0;
    break;
      default:
    stab_coeff = lpcoeff_crosswind*pow(hK,lpexponent_crosswind);
    break;
  }
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
    // parameter is lpcoeff*pow(hK,lpexponent)
               Entries[p] += stab_coeff*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
    } // end LP_CROSSWIND

  } // endfor i
  delete params;
  delete coeffs;
 if (TDatabase::ParamDB->SC_VERBOSE>1)
    OutPut("LPS streamline done" << endl);
} // UltraLocalProjectionStreamlinePLaplacian


//**************************************************************
// computes the local projection to Q0 on a coarse grid
// uh      -  current fe function
// uh_proj -  fe function for projection, pw constant
//**************************************************************

void LocalProjectionCoarseGridQ0(TFEFunction2D *uh,
         TFEFunction2D *uh_proj,
               CoeffFct2D *Coeff,
               int convection_flag)
{
  int i, j, iq, index;
  int N_Cells, N_Edges;
  int *GlobalNumbers_fine, *BeginIndex_fine, *GlobalNumbers_coarse, *BeginIndex_coarse, *DOF;
  double sx, sy, b1, b2, area, detJK, weight_det, value;
  double *Values_fine, *coeffs, *params, *Values_coarse,  x[4], y[4], val[4];
  double gauss2_x[4]=
  {
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  double gauss2_y[4]=
  {
    -0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  TCollection *Coll;
  TFESpace2D *fespace_fine, *fespace_coarse;
  TBaseCell *cell, *parent_cell, *child_cell;

  OutPut("compute local projection to Q0 on coarse grid"<<endl);
  coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * SizeOfDouble);

  fespace_fine = uh->GetFESpace2D();
  fespace_coarse = uh_proj->GetFESpace2D();
  // get collection, same for fine and coarse fe space
  Coll = fespace_fine->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
   // get arrays for dofs  
  BeginIndex_fine = fespace_fine->GetBeginIndex();
  GlobalNumbers_fine = fespace_fine->GetGlobalNumbers();
 
  // get values of fe function 
  Values_fine = uh->GetValues();     
  Values_coarse = uh_proj->GetValues();     
  // pw constant space
  memset(Values_coarse, 0.0, N_Cells * SizeOfDouble);
  
  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);   
    // number of edges
    N_Edges = cell->GetN_Edges();
    // get barycenter
    sx = sy = 0;
    for (j=0;j<N_Edges;j++)
    {
      x[j] = cell->GetVertex(j)->GetX(); 
      sx += x[j];
      y[j] = cell->GetVertex(j)->GetY();
      sy += y[j];
    }
    sx /= N_Edges;
    sy /= N_Edges;
    // area of parallelogramm with vector product
    area = fabs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
    // functional determinant
    detJK = area/4.0;
    // get coefficients in center of mesh cell 
    Coeff(1, &sx ,&sy, &params, &coeffs);
    // take the desired direction
    switch(convection_flag)
    {
      case 0:
      // convection
        b1 = coeffs[1];
        b2 = coeffs[2];
        break;
      case 1:
      // normalized crosswind convection
        b2 = sqrt(coeffs[1]*coeffs[1] + coeffs[2] * coeffs[2]);
  b1 = -coeffs[2]/b2;
  b2 = coeffs[1]/b2;
  break;
    }
    // compute integral of b times grad_uh
    // loop over the quadrature points
     for (iq = 0;iq < 4; iq++)
    {
       sx = (x[1]-x[0])* gauss2_x[iq] + (x[3]-x[0])*gauss2_y[iq] + x[1] + x[3];
       sx /= 2.0;
       sy = (y[1]-y[0])* gauss2_x[iq] + (y[3]-y[0])*gauss2_y[iq] + y[1] + y[3];
       sy /= 2.0;
       uh->FindGradientLocal(cell,i,sx,sy,val);
       value = b1 * val[1] + b2 * val[2];
       weight_det = detJK;
       Values_coarse[i] += value * weight_det;
    }   
  }
   
  // compute the mean values on the macros
    // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    // cell not yet child of a macro cell
    if (cell->GetClipBoard()>=0)
    {
      parent_cell = cell->GetParent();
      if (parent_cell == NULL)
      {
  OutPut("no coarse grid to project !!!" << endl);
        delete params;
        delete coeffs;
  return;
      }
      area = 0.0;
      value = 0.0;
      for (j=0;j<4;j++)
      {
  child_cell = parent_cell->GetChild(j);
  // get index of pw constant space
  index = child_cell->GetClipBoard();
  // integral of projected quantity
  value += Values_coarse[index];
  // area of macro cell
  area += child_cell->GetMeasure();
      }
      value /= area;
      for (j=0;j<4;j++)
      {
  child_cell = parent_cell->GetChild(j);
  // get index of pw constant space
  index = child_cell->GetClipBoard();
        // projection
  Values_coarse[index] = value;
  child_cell->SetClipBoard(-1);
      }
    }
  }

  delete[] params;
  delete[] coeffs;
}


void LocalProjectionCrossWindCoarseGridQ0(TDomain *Domain, int mg_level,
            TFEFunction2D *uh,
            TFEFunction2D *uh_proj,
            CoeffFct2D *Coeff,
            double *rhs,
                        int convection_flag) 
{
  int i, j, k, iq, N_Cells, N_Edges, N_Unknowns, index, dof_index;
  int *global_numbers, *begin_index, *dof;
  double  sx, sy, b1, b2, area, detJK, weight_det, value, value1, value_proj, hK;
  double *Values_fine, *coeffs, *params, *Values_coarse,  x[4], y[4], val[4];
  double *sol, *sol_copy, *proj_test, area_coarse, norm_b, stab_coeff_cw;
  double lpcoeff_crosswind = TDatabase::ParamDB->LP_CROSSWIND_COEFF;
  double gauss2_x[4]=
  {
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  double gauss2_y[4]=
  {
    -0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  TCollection *coll_coarse, *coll;
  TBaseCell *cell, *child_cell;
  TFESpace2D *fespace;
 
   OutPut("update rhs of crosswind local projection to Q0 on coarse grid"<<endl);
  // get coarse grid
  coll_coarse=Domain->GetCollection(It_EQ, mg_level+TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR-1);
  if (coll_coarse == NULL)
  {
    OutPut("No coarse grid !!!" << endl);
    return;
  }
   coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * SizeOfDouble);

  // get fine grid
  coll = Domain->GetCollection(It_EQ, mg_level+TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR);
  N_Cells = coll->GetN_Cells();
  // initialise ClipBoard for fine grid
  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);
  
  // fespace on the fine grid
  fespace = uh->GetFESpace2D();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // copy solution
  N_Unknowns = uh->GetLength();
  sol = uh->GetValues();
  sol_copy = new double[N_Unknowns];
  memcpy(sol_copy,sol, N_Unknowns*SizeOfDouble);

  proj_test = new double[N_Unknowns];
  // the already computed projections of the current solution 
  Values_coarse = uh_proj->GetValues();  
  // loop over coarse grid
  N_Cells = coll_coarse->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = coll_coarse->GetCell(i);  
    area_coarse = cell->GetMeasure();
    hK = cell->GetDiameter();
    memset(proj_test,0,N_Unknowns*SizeOfDouble);
    for (j=0;j<4;j++)
    {
        // get child cell
  child_cell = cell->GetChild(j);
  // get degrees of freedom in the cell on the fine grid
  index = child_cell->GetClipBoard();
  dof = global_numbers+begin_index[index];
        // number of edges
        N_Edges = child_cell->GetN_Edges();
  // for each degree of freedom
       for (k=0;k<N_Edges;k++)
       {
         x[k] = child_cell->GetVertex(k)->GetX(); 
         sx += x[k];
         y[k] = child_cell->GetVertex(k)->GetY();
         sy += y[k];
       }
       sx /= N_Edges;
       sy /= N_Edges;
       // area of parallelogramm with vector product
       area = fabs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
       // functional determinant
       detJK = area/4.0;
       // get coefficients in center of mesh cell 
       Coeff(1, &sx ,&sy, &params, &coeffs);
       // take the desired direction
       switch(convection_flag)
       {
         case 0:
         // convection
           b1 = coeffs[1];
           b2 = coeffs[2];
           break;
        case 1:
        // normalized crosswind convection
           b2 = sqrt(coeffs[1]*coeffs[1] + coeffs[2] * coeffs[2]);
           b1 = -coeffs[2]/b2;
     b2 = coeffs[1]/b2;
     break;
        }
        // loop over the local dofs
  for (k=0;k<N_Edges;k++)
  {
    memset(sol,0,N_Unknowns*SizeOfDouble);
    dof_index = dof[k];
    sol[dof_index] = 1.0;
          // compute integral of b times grad_vh
          // loop over the quadrature points
          for (iq = 0;iq < 4; iq++)
          {
             sx = (x[1]-x[0])* gauss2_x[iq] + (x[3]-x[0])*gauss2_y[iq] + x[1] + x[3];
             sx /= 2.0;
             sy = (y[1]-y[0])* gauss2_x[iq] + (y[3]-y[0])*gauss2_y[iq] + y[1] + y[3];
             sy /= 2.0;
             uh->FindGradientLocal(child_cell,index,sx,sy,val);
             value = b1 * val[1] + b2 * val[2];
             weight_det = detJK;
             proj_test[dof_index] += value * weight_det/area_coarse;
           }   
           OutPut(dof_index << " " << proj_test[dof_index] << ":");
  } 
    } // end j
    // projection of test function for macro cell i computed
    // now compute contribution to rhs
    for (j=0;j<4;j++)
    {
        // get child cell
  child_cell = cell->GetChild(j);
  // get degrees of freedom in the cell on the fine grid
  index = child_cell->GetClipBoard();
  dof = global_numbers+begin_index[index];
        // number of edges
        N_Edges = child_cell->GetN_Edges();
  // for each degree of freedom
       for (k=0;k<N_Edges;k++)
       {
         x[k] = child_cell->GetVertex(k)->GetX(); 
         sx += x[k];
         y[k] = child_cell->GetVertex(k)->GetY();
         sy += y[k];
       }
       sx /= N_Edges;
       sy /= N_Edges;
       // area of parallelogramm with vector product
       area = fabs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
       // functional determinant
       detJK = area/4.0;
       // get coefficients in center of mesh cell 
       Coeff(1, &sx ,&sy, &params, &coeffs);
       // take the desired direction
       switch(convection_flag)
       {
         case 0:
         // convection
           b1 = coeffs[1];
           b2 = coeffs[2];
           break;
        case 1:
        // normalized crosswind convection
           b2 = sqrt(coeffs[1]*coeffs[1] + coeffs[2] * coeffs[2]);
           b1 = -coeffs[2]/b2;
     b2 = coeffs[1]/b2;
     break;
        }
        // value of projection
  value_proj = Values_coarse[index];
        // loop over the local dofs
  for (k=0;k<N_Edges;k++)
  {
    dof_index = dof[k];
          // compute integral of b times grad_vh
          // loop over the quadrature points
    value1 = 0;
          for (iq = 0;iq < 4; iq++)
          {
             sx = (x[1]-x[0])* gauss2_x[iq] + (x[3]-x[0])*gauss2_y[iq] + x[1] + x[3];
             sx /= 2.0;
             sy = (y[1]-y[0])* gauss2_x[iq] + (y[3]-y[0])*gauss2_y[iq] + y[1] + y[3];
             sy /= 2.0;
             uh->FindGradientLocal(child_cell,index,sx,sy,val);
             value = b1 * val[1] + b2 * val[2] - value_proj;
       value = fabs(value) * value * proj_test[dof_index];
             weight_det = detJK;
             value1 += value * weight_det;
           }
           norm_b = sqrt(b1*b1+b2*b2);
           if (norm_b > 1e-20)
           {
           // h/|b|
             stab_coeff_cw = 2*lpcoeff_crosswind*hK/norm_b;
           }
           else
             stab_coeff_cw = 0.0;
     
           rhs[dof_index] += stab_coeff_cw * value1;
     //OutPut( -hK*value1 << " ");
  } 
  //OutPut(endl);
    } // end j
   
  }
  
  // recover solution
  memcpy(sol,sol_copy, N_Unknowns*SizeOfDouble);

  delete[] proj_test;
  delete[] sol_copy; 
  delete[] params;
  delete[] coeffs;

}

// Adaptive post processing basd on Friedhelm Schiweck talk at MAFELAP 09 - Sashi
void AdaptivePostProcess(TFEFunction2D *FeFunction, double *PostSol, bool DirichletBC)
{
  int i, j, k, l, CoarseOrder, N_Points, N_U;
  int N_Cells, N_CoarseDOF, N_DOF;
  int N_UsedElements;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int OrderDiff;

  double **CoarseValues, *CoarseValue, *LPS_sol;
  double val, w, *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double rhs[MaxN_BaseFunctions2D];
  double Values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double PCValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double *Value, *CurrentValue, sol;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double *W, maxbubble=-1E8, minbubble=1E8;


  bool SecondDer[1] = { false };

  TFESpace2D *fespace;
  TCollection *coll;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  FE2D CurrEleID, UsedElements[1];
  TNodalFunctional2D *nf;

  fespace = FeFunction->GetFESpace2D();
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();
  LPS_sol = FeFunction->GetValues();
  N_U = FeFunction->GetLength();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  W = new double[N_U];
  memset(W, 0, N_U*SizeOfDouble);


  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);
    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;
    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    N_UsedElements = 1;
    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();
    N_CoarseDOF = CoarseBF->GetDimension();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements, 
                           coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     BF->GetDerivatives(D00, xi[j], eta[j], Values[j]);


    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    memset(rhs, 0, N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
     {
      CoarseValue = CoarseValues[j];
      Value = Values[j];
      w = AbsDetjk[j]*weights[j];

//    find the lps solution at this quadrature point
      sol=0.;
      for(k=0;k<N_DOF;k++)
       sol+= LPS_sol[DOF[k]]*Value[k];

      for(k=0;k<N_CoarseDOF;k++)
       {
        val = w*CoarseValue[k];
        rhs[k] += sol*val;

        for(l=0;l<N_CoarseDOF;l++)
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];

       } // end for k
     } // for(j=0;j<N_Points;j++)

//     for(j=0;j<N_CoarseDOF;j++)
//       for(k=0;k<N_CoarseDOF;k++)
//         cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
//     cout  << endl;
//     for(j=0;j<N_CoarseDOF;j++)
//         cout << j << "  " << rhs[j] << endl;

    SolveLinearSystemLapack(G, rhs, N_CoarseDOF, N_CoarseDOF);

// interpolate the discont solution to the original LPS space

    nf = CurrentElement->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);

    for(j=0;j<N_Points;j++)
     CoarseBF->GetDerivatives(D00, xi[j], eta[j], PCValues[j]);

    memset(PointValues, 0, N_Points*SizeOfDouble);

    for(j=0;j<N_Points;j++)
     for(k=0;k<N_CoarseDOF;k++)
      PointValues[j] +=  rhs[k]*PCValues[j][k];

    nf->GetAllFunctionals(coll, cell, PointValues,
                          FunctionalValues);

    for(j=0;j<N_DOF;j++)
     {
      PostSol[DOF[j]] += FunctionalValues[j];
      W[DOF[j]] += 1.;
      if(j==N_DOF-1)
        {
//          cout<< j<< " PostSol " << PostSol[DOF[j]] << endl;
          if (maxbubble< PostSol[DOF[j]]) maxbubble = PostSol[DOF[j]];
          if (minbubble > PostSol[DOF[j]]) minbubble = PostSol[DOF[j]];
        }
     }
   } // for(i=0;i<N_Cells;i++)

  for(i=0;i<N_U;i++)
   PostSol[i] /=W[i];

//   cout<<" maxbubble " << maxbubble  <<" minbubble " << minbubble  << endl;
//     exit(0);
}


//Fefunction can be a different FEspace - Sashi
void AddALEStreamlineLPS(TSquareMatrix2D* A, int N_FeFunct, TFEFunction2D **FeFunct,
                         double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int *U_GlobalNumbers, *U_BeginIndex, *U_DOF;  
  int *W_GlobalNumbers, *W_BeginIndex, *W_DOF;  
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF, N_UDOF, N_WDOF;
  TCollection *Coll;
  TFESpace2D *fespace, *U_fespace, *W_fespace;
  FE2D CurrEleID, U_CurrEleID, W_CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement, *U_CurrentElement, *W_CurrentElement;
  TBaseFunct2D *BF, *CoarseBF, *U_BF, *W_BF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  BaseFunct2D  U_BF_ID, W_BF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  TFEFunction2D *uh1, *uh2, *wh1, *wh2;  
  
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
//   double **ChildValues, *ChildValue;
  double **U_Values, *U_Value; 
  double **W_Values, *W_Value;
  
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2, *Values3, *Values4;
  double BValue[MaxN_BaseFunctions2D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  if(N_FeFunct==2)
   {
    uh1 = FeFunct[0];
    uh2 = FeFunct[1];      
   }
  else if(N_FeFunct==4)
   {
    uh1 = FeFunct[0];
    uh2 = FeFunct[1];   
    wh1 = FeFunct[2];
    wh2 = FeFunct[3]; 
   }
  else
   {
    cout << "N_FeFunct must be 2 or 4 " <<endl;
    exit(0);    
   }
   
  U_fespace = uh1->GetFESpace2D();
  U_BeginIndex = U_fespace->GetBeginIndex();
  U_GlobalNumbers = U_fespace->GetGlobalNumbers(); 
  Values1 = uh1->GetValues();
  Values2 = uh2->GetValues();
    
  if(N_FeFunct==4) 
   {
    W_fespace = wh1->GetFESpace2D();  
    W_BeginIndex = W_fespace->GetBeginIndex();
    W_GlobalNumbers = W_fespace->GetGlobalNumbers(); 
    Values3 = wh1->GetValues();
    Values4 = wh2->GetValues();       
   }
  
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);
    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    U_CurrEleID = U_fespace->GetFE2D(i, cell);
    U_CurrentElement = TFEDatabase2D::GetFE2D(U_CurrEleID);
    U_BF = U_CurrentElement->GetBaseFunct2D();
    U_BF_ID = U_BF->GetID();
    N_UDOF = U_BF->GetDimension();    
    
    W_CurrEleID = W_fespace->GetFE2D(i, cell);
    W_CurrentElement = TFEDatabase2D::GetFE2D(W_CurrEleID);
    W_BF = W_CurrentElement->GetBaseFunct2D();
    W_BF_ID = W_BF->GetID();
    N_WDOF = W_BF->GetDimension();   
    
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];    
    U_DOF = U_GlobalNumbers + U_BeginIndex[i];  
    W_DOF = W_GlobalNumbers + W_BeginIndex[i];  

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    U_Values  = TFEDatabase2D::GetOrigElementValues(U_BF_ID, D00);
    W_Values  = TFEDatabase2D::GetOrigElementValues(W_BF_ID, D00);    

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);
    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      U_Value  = U_Values[j];
      W_Value  = W_Values[j];      
      
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;     
      // compute components of uh in j
      for(k=0;k<N_UDOF;k++)
      {
        l = U_DOF[k];
        valx += U_Value[k]*Values1[l];
        valy += U_Value[k]*Values2[l];
      }
      // sub mesh velo (uh-wh)
      for(k=0;k<N_WDOF;k++)
      {
        l = W_DOF[k];
        valx -= W_Value[k]*Values3[l];
        valy -= W_Value[k]*Values4[l];
      }

      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm
#else


FE3D GetElement3D(TBaseCell *cell, int CoarseOrder)
{
  FE3D ele = (FE3D)0;
  Shapes shapetype;

  shapetype = cell->GetType();
  switch(shapetype)
  {
    // regularly refined Hexahedron
    case Hexahedron:
    case Brick:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_Q0_3D_H_M;
        break;

        case 1:
          ele = D_P1_3D_H_M;
        break;

        case 2:
          ele = D_P2_3D_H_M;
        break;

        case 3:
          ele = D_P3_3D_H_M;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_Q00_3D_H_M;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 3" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break; // end regularly refined quadrilateral

    case Tetrahedron:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_P0_3D_T_A;
        break;

        case 1:
          ele = D_P1_3D_T_A;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_P00_3D_T_A;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 1" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break;
    default:
      OutPut("Invalid shape" << endl);
      exit(-1);
  } // end switch reftype
  return ele;
}

// ADDED ON 17.06.2011 BY SASHI
void AddStreamlineTerm(TSquareMatrix3D* A, TFEFunction3D *uh1,
                       TFEFunction3D *uh2, TFEFunction3D *uh3,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  int N_UsedElements, N_Points;  
  int i1, i2;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;

  TCollection *Coll;
  TFESpace3D *fespace;
  FE3D CurrEleID, UsedElements[2];
  TFE3D *CurrentElement, *CoarseElement;
  TBaseFunct3D *BF, *CoarseBF;
  BaseFunct3D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };

  double *xi, *eta, *zeta, *weights;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[2*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[2*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;  
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy, valz;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s, hK;
  double *Entries;
  double *Values1, *Values2, *Values3;
  double BValue[MaxN_BaseFunctions3D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE3D(i, cell);
    CurrentElement = TFEDatabase3D::GetFE3D(CurrEleID);

    BF = CurrentElement->GetBaseFunct3D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement3D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase3D::GetFE3D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct3D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase3D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();
    Values3 = uh3->GetValues();
    
    PCValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    ChildValuesX = TFEDatabase3D::GetOrigElementValues(BF_ID, D100);
    ChildValuesY = TFEDatabase3D::GetOrigElementValues(BF_ID, D010);
    ChildValuesZ = TFEDatabase3D::GetOrigElementValues(BF_ID, D001);   
    ChildValues  = TFEDatabase3D::GetOrigElementValues(BF_ID, D000);

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];     
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;
      valz = 0.0;
      
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
        valz += ChildValue[k]*Values3[l];        
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k] + valz*ChildValueZ[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm


// stabilisation of full gradient (scalar)
void UltraLocalProjection(TSquareMatrix3D* A, 
                          double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  int N_UsedElements, N_Points;  
  int i1, i2;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;

  TCollection *Coll;
  TFESpace3D *fespace;
  FE3D CurrEleID, UsedElements[2];
  TFE3D *CurrentElement, *CoarseElement;
  TBaseFunct3D *BF, *CoarseBF;
  BaseFunct3D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };

  double *xi, *eta, *zeta, *weights;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue; 
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;  
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy, valz;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s, hK;
  double *Entries;
  double BValue[MaxN_BaseFunctions3D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE3D(i, cell);
    CurrentElement = TFEDatabase3D::GetFE3D(CurrEleID);

    BF = CurrentElement->GetBaseFunct3D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement3D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase3D::GetFE3D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct3D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase3D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];
    
    PCValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    ChildValuesX = TFEDatabase3D::GetOrigElementValues(BF_ID, D100);
    ChildValuesY = TFEDatabase3D::GetOrigElementValues(BF_ID, D010);
    ChildValuesZ = TFEDatabase3D::GetOrigElementValues(BF_ID, D001);   
 
    memset(H, 0, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];     
      w = AbsDetjk[j]*weights[j];

      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*3*N_DOF+l      ] += val*ChildValueX[l];
          H[k*3*N_DOF+l+N_DOF] += val*ChildValueY[l];
          H[k*3*N_DOF+l+2*N_DOF] += val*ChildValueZ[l];  
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*(  ChildValueX[k]*ChildValueX[l]
                                     + ChildValueY[k]*ChildValueY[l]
                                     + ChildValueZ[k]*ChildValueZ[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3*N_DOF, 3*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */
    
    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*3*N_DOF+l      ]*H[i2*3*N_DOF+m      ]
                                           +H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF]
                                           +H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*3*N_DOF+l      ] * H[i2*3*N_DOF+m      ];
          s += P[i2*3*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF];
          s += P[i2*3*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*3*N_DOF+m      ] * H[i1*3*N_DOF+l      ];
          s += P[i1*3*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF];
          s += P[i1*3*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm


//**************************************************************
//  UltraLocalProjection
//  checked: Volker John 08/02/19
//**************************************************************
#ifdef __2D__
void UltraLocalProjection(void* A, bool ForPressure, CoeffFct2D *Coeff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s, sx, sy;
  int i1, i2, ij, N_Edges;
  double hK, *coeffs, *params;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  OutPut("LPS full gradient" << endl);

  coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * SizeOfDouble);

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    OutPut("Local projection stabilization is implemented only for full gradient!" << endl);
    exit(-1);
  }

  if(ForPressure)
  {
    lpcoeff = TDatabase::ParamDB->LP_PRESSURE_COEFF;
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  // get fespace and matrices
  if(ForPressure)
  {
    fespace = (TFESpace2D*)(((TMatrix2D *)A)->GetStructure()->GetTestSpace());
    ActiveBound = -1;
    RowPtr = ((TMatrix2D *)A)->GetRowPtr();
    KCol = ((TMatrix2D *)A)->GetKCol();
    Entries = ((TMatrix2D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix2D *)A)->GetFESpace();
    ActiveBound = fespace->GetActiveBound();
    RowPtr = ((TSquareMatrix2D *)A)->GetRowPtr();
    KCol = ((TSquareMatrix2D *)A)->GetKCol();
    Entries = ((TSquareMatrix2D *)A)->GetEntries();
    // cout << "not for pressure" << endl;
  }
  // get collection
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();
  // get arrays for dofs
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);
  // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
      // get cell
    cell = coll->GetCell(i);
    // get diameter of the cell
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
  case 0: // diameter
      hK = cell->GetDiameter();
      break;
  case 1: // with reference map
      hK = cell->GetLengthWithReferenceMap();
      break;
  case 2: // shortest edge
      hK = cell->GetShortestEdge();
      break;
  case 3: // measure
      hK = cell->GetMeasure();
      hK = pow(hK,1.0/3.0);
      break;
  case 4: // mesh cell in convection direction
      N_Edges = cell->GetN_Edges();
      sx = sy = 0;
      for (ij=0;ij<N_Edges;ij++)
      {
    TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
    sx += TDatabase::ParamDB->INTERNAL_VERTEX_X[ij];
    TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
    sy += TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij];
      }
      if (N_Edges==3)
    TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      // center of mesh cell
      sx /= N_Edges;
      sy /= N_Edges;
      hK = cell->GetDiameter(); 
      // get coefficients in center of mesh cell 
      Coeff(1, &sx ,&sy, &params, &coeffs);
      hK = Mesh_size_in_convection_direction(hK, coeffs[1], coeffs[2]); 
      break;
  default: // diameter
      hK = cell->GetDiameter();
      break;
    }
    // get finite element in the cell
    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);
    // get basis functions of the finite element
    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();
    // compute order of the local coarse space
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // get type of the mesh cell
    shapetype = cell->GetType();

    // determine finite element of the local coarse space
    switch(shapetype)
    {
      // regularly refined quadrilateral
      // discontinuous finite elements  
      case Quadrangle:
      case Parallelogram:
      case Rectangle:
        switch(CoarseOrder)
        {
          case 0:
            UsedElements[0] = C_Q0_2D_Q_M;
          break;

          case 1:
            UsedElements[0] = D_P1_2D_Q_M;
          break;

          case 2:
            UsedElements[0] = D_P2_2D_Q_M;
          break;

          case 3:
            UsedElements[0] = D_P3_2D_Q_M;
          break;

          case 4:
            UsedElements[0] = D_P4_2D_Q_M;
          break;

          default:
            OutPut("Projection space is defined upto order 4" << endl);
            exit(-1);
        } // end switch CoarseOrder
      break; // end regularly refined quadrilateral

      case Triangle:
        switch(CoarseOrder)
        {
          case 0:
            UsedElements[0] = C_P0_2D_T_A;
          break;

          case 1:
            UsedElements[0] = D_P1_2D_T_A;
          break;

          case 2:
            UsedElements[0] = D_P2_2D_T_A;
          break;

          case 3:
            UsedElements[0] = D_P3_2D_T_A;
          break;

          case 4:
            UsedElements[0] = D_P4_2D_T_A;
          break;

          default:
            OutPut("Projection space is defined upto order 4" << endl);
            exit(-1);
        }
      break;
      default: 
    OutPut("No coarse finite elements for mesh cell type "
     << shapetype << " implemented !!!" << endl);
    exit(4711);
    } // end switch reftype

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);
    // number of dof in the local coarse space
    N_CoarseDOF = CoarseBF->GetDimension();
    // get function values for the basis functions of the coarse space
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);
    // initialize array G, stores mass matrix of coarse space
    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
      // values of the coarse basis functions in the quad points
      CoarseValue = CoarseValues[j];
      // factor for numerical quadrature
      w = AbsDetjk[j]*weights[j];
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
    // first factor of integrand
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
      // update integral 
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    // this should be the same as CoarseValues ???
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);
    // get derivatives of fe functions of approximation space
    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    // initialize array H, holds mixed products of coarse basis 
    // functions and derivatives of fine basis functions
    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);
    // initialize array LocMatrix, holds products of derivatives
    // of fine basis functions
    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
  // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = AbsDetjk[j]*weights[j];
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
  // compute products of coarse basis function and
  // derivatives of fine basis functions
  // values for the same quad point are stored after
  // each other (x deriv of all fine fcts., y deriv
  // of all fine fcts.) tested with the k-th coarse fct. 
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix (fine-fine coupling)
      // LocMatrix will store the local update of the global
      // system matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    // copy array H to array P
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // solve G * X = H, solution X stored on H
    // the right hand side and the solution are stored column wise
    // right hand side: a fine fct. tested with all coarse fcts. 
    // X - coefficients of local L^2 projection of grad of fine function
    //     ((first fct.)_x, first coarse fct.), ((second fct.)_x, first coarse fct.), 
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);
    //SolveMultipleSystems(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*// checked 08/02/19
    double val;
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
    {
  for(l=0;l<2*N_DOF;l++)
  {
      val = 0;
      
      for(k=0;k<N_CoarseDOF;k++)
      {
    val += G[j*N_CoarseDOF+k] * H[l+k*2*N_DOF];
    OutPut(l+2*k*N_DOF << " ");
      }
      OutPut(j << " " << l << " "  << P[j*2*N_DOF+l]<<  " " << val << endl);
  }
    }
    exit(1);
    */

    // update LocMatrix
    // proj-proj coupling (coarse-coarse coupling)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
  // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ]
                                           +H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling (fine-coarse couplings)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
  // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // get numbers of global dof belonging to this mesh cell
    DOF = GlobalNumbers + BeginIndex[i];

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
    // parameter is lpcoeff*pow(hK,lpexponent)
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
        break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
  delete params;
  delete coeffs;
} // UltraLocalProjection



bool TestCell(TBaseCell *cell)
{
  int i, N_;
  double x, y;
  bool val = true;
  
  N_ = cell->GetN_Vertices();
  
  for(i=0;i<N_;i++)
  {
    cell->GetVertex(i)->GetCoords(x,y);
    if(x<TDatabase::ParamDB->P6 && y<TDatabase::ParamDB->P6)
      val = val && true;
    else 
      val = val && false;    
  }
  return val;
}

// stabilisation of function velocity
void UltraLocalProjectionFunction(void* A, bool ForPressure)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    OutPut("Local projection stabilization is implemented only for full gradient!" << endl);
    exit(-1);
  }

  if(ForPressure)
  {
    lpcoeff = -(TDatabase::ParamDB->LP_PRESSURE_COEFF);
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  if(ForPressure)
  {
    fespace = (TFESpace2D*)(((TMatrix2D *)A)->GetStructure()->GetTestSpace());
    ActiveBound = -1;
    RowPtr = ((TMatrix2D *)A)->GetRowPtr();
    KCol = ((TMatrix2D *)A)->GetKCol();
    Entries = ((TMatrix2D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix2D *)A)->GetFESpace();
    ActiveBound = fespace->GetActiveBound();
    RowPtr = ((TSquareMatrix2D *)A)->GetRowPtr();
    KCol = ((TSquareMatrix2D *)A)->GetKCol();
    Entries = ((TSquareMatrix2D *)A)->GetEntries();
    // cout << "not for pressure" << endl;
  }

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // cout << "N_Points: " << N_Points << endl;

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    // cout << "N_CoarseDOF: " << N_CoarseDOF << endl;

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // (u,v)--matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // (\pi u, \pi v) proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ]
                                           +H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // function-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += 1000.*hK*hK*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // UltraLocalProjection

double UltraLocalErrorSmooth(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double error, locerror;
  double exactval[4];

  error = 0.0;

  fespace = uh->GetFESpace2D();
  Values = uh->GetValues();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();
    
    if((TestCell(cell)) == false) continue;
    
    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    // only two right-hand sides (x and y derivative)
    memset(H, 0, N_CoarseDOF*2*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        val = Values[DOF[k]];
        valx += ChildValueX[k]*val;
        valy += ChildValueY[k]*val;
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k*2  ] += val*valx;
        H[k*2+1] += val*valy;
      } // end for k

      // grad-grad term
      locerror += w*(valx*valx + valy*valy);

    } // end for j
    memcpy(P, H, N_CoarseDOF*2*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2, 2);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2  ]*H[i2*2  ]
                                       +H[i1*2+1]*H[i2*2+1]);
    locerror += s;

    // grad-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2*2  ] * H[i2*2  ];
      s += P[i2*2+1] * H[i2*2+1];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1*2  ] * H[i1*2  ];
      s += P[i1*2+1] * H[i1*2+1];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalProjection

#endif // __2D__

// stabilisation of full gradient (velocity or pressure)
void UltraLocalProjection3D(void* A, bool ForPressure)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace3D *fespace;
  FE3D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE3D *CurrentElement, *CoarseElement;
  TBaseFunct3D *BF, *CoarseBF;
  BaseFunct3D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *zeta, *weights;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s;
  int i1, i2, i3;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    OutPut("Local projection stabilization is implemented only for full gradient!" << endl);
    exit(-1);
  }

  if(ForPressure)
  {
    lpcoeff = -(TDatabase::ParamDB->LP_PRESSURE_COEFF);
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  if(ForPressure)
  {
    fespace = (TFESpace3D*)(((TMatrix3D *)A)->GetStructure()->GetTestSpace());
    ActiveBound = -1;
    RowPtr = ((TMatrix3D *)A)->GetRowPtr();
    KCol = ((TMatrix3D *)A)->GetKCol();
    Entries = ((TMatrix3D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix3D *)A)->GetFESpace();
    ActiveBound = fespace->GetActiveBound();
    RowPtr = ((TSquareMatrix3D *)A)->GetRowPtr();
    KCol = ((TSquareMatrix3D *)A)->GetKCol();
    Entries = ((TSquareMatrix3D *)A)->GetEntries();
    // cout << "not for pressure" << endl;
  }

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE3D(i, cell);
    CurrentElement = TFEDatabase3D::GetFE3D(CurrEleID);

    BF = CurrentElement->GetBaseFunct3D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement3D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase3D::GetFE3D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct3D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase3D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    PCValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    ChildValuesX = TFEDatabase3D::GetOrigElementValues(BF_ID, D100);
    ChildValuesY = TFEDatabase3D::GetOrigElementValues(BF_ID, D010);
    ChildValuesZ = TFEDatabase3D::GetOrigElementValues(BF_ID, D001);

    memset(H, 0, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*3*N_DOF+l        ] += val*ChildValueX[l];
          H[k*3*N_DOF+l+N_DOF  ] += val*ChildValueY[l];
          H[k*3*N_DOF+l+2*N_DOF] += val*ChildValueZ[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]
                                     +ChildValueZ[k]*ChildValueZ[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3*N_DOF, 3*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*3*N_DOF+l        ]*H[i2*3*N_DOF+m        ]
                                           +H[i1*3*N_DOF+l+N_DOF  ]*H[i2*3*N_DOF+m+N_DOF  ]
                                           +H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
//        for(i3=0;i3<N_CoarseDOF;i3++)
//        {
//          s += P[i3*3*N_DOF+l        ] * H[i3*3*N_DOF+m        ];
//          s += P[i3*3*N_DOF+l+N_DOF  ] * H[i3*3*N_DOF+m+N_DOF  ];
//          s += P[i3*3*N_DOF+l+2*N_DOF] * H[i3*3*N_DOF+m+2*N_DOF];
//        }
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*3*N_DOF+l        ] * H[i2*3*N_DOF+m        ];
          s += P[i2*3*N_DOF+l+N_DOF  ] * H[i2*3*N_DOF+m+N_DOF  ];
          s += P[i2*3*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*3*N_DOF+m        ] * H[i1*3*N_DOF+l        ];
          s += P[i1*3*N_DOF+m+N_DOF  ] * H[i1*3*N_DOF+l+N_DOF  ];
          s += P[i1*3*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
  cout << "\n end of lps \n" << endl;
} // UltraLocalProjection

double UltraLocalError3D(TFEFunction3D *uh, DoubleFunct3D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff)
{
  // cout << "\n start of lps error \n" << endl;
  
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace3D *fespace;
  FE3D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE3D *CurrentElement, *CoarseElement;
  TBaseFunct3D *BF, *CoarseBF;
  BaseFunct3D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double *xi, *eta, *zeta, *weights;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy, valz;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s;
  int i1, i2, i3;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double error, locerror;
  double exactval[5];

  error = 0.0;

  fespace = uh->GetFESpace3D();
  Values = uh->GetValues();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE3D(i, cell);
    CurrentElement = TFEDatabase3D::GetFE3D(CurrEleID);

    BF = CurrentElement->GetBaseFunct3D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement3D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase3D::GetFE3D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct3D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase3D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    ChildValuesX = TFEDatabase3D::GetOrigElementValues(BF_ID, D100);
    ChildValuesY = TFEDatabase3D::GetOrigElementValues(BF_ID, D010);
    ChildValuesZ = TFEDatabase3D::GetOrigElementValues(BF_ID, D001);

    // only two right-hand sides (x and y derivative)
    memset(H, 0, N_CoarseDOF*3*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];

      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      valz = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        val = Values[DOF[k]];
        valx += ChildValueX[k]*val;
        valy += ChildValueY[k]*val;
        valz += ChildValueZ[k]*val;
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], Z[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];
      valz -= exactval[3];

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k*3  ] += val*valx;
        H[k*3+1] += val*valy;
        H[k*3+2] += val*valz;
      } // end for k

      // grad-grad term
      locerror += w*(valx*valx + valy*valy + valz*valz);

    } // end for j
    memcpy(P, H, N_CoarseDOF*3*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3, 3);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*3  ]*H[i2*3  ]
                                       +H[i1*3+1]*H[i2*3+1]
                                       +H[i1*3+2]*H[i2*3+2]);
    locerror += s;

    // grad-proj coupling
    s = 0;
//    for(i3=0;i3<N_CoarseDOF;i3++)
//    {
//      s += P[i3*3  ] * H[i3*3  ];
//      s += P[i3*3+1] * H[i3*3+1];
//      s += P[i3*3+2] * H[i3*3+2];
//    }
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2*3  ] * H[i2*3  ];
      s += P[i2*3+1] * H[i2*3+1];
      s += P[i2*3+2] * H[i2*3+2];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1*3  ] * H[i1*3  ];
      s += P[i1*3+1] * H[i1*3+1];
      s += P[i1*3+2] * H[i1*3+2];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  // error = 0;
  cout << " end of local error " << endl;
  return error;
} // UltraLocalProjection


#endif
