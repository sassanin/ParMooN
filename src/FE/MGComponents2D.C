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
// MGComponents2D.C
//
// Purpose:     components for multigrid in 2d
//
// Author:      Gunar Matthies          27.01.1999
//              Volker John             27.10.1999  
//
//              parallel methods  (Sashikumaar Ganesan) 13.10.2009
// =======================================================================

#ifdef _MPI
#  include "mpi.h"
#endif

#include <LinAlg.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <MooNMD_Io.h>

#include <stdlib.h>
#include <string.h>

#ifdef __2D__

#include <SquareMatrixNSE2D.h>

#define AT(i,j) (a[(j)*LDA+(i)])
#define A(i,j) (a[(i)*LDA+(j)])

/* solve diagonal Vanka system */
void SolveDiagonalVanka2D(double *a, double *b, int N_U, int N_P, int LDA)
// Arguments:
//    a         double array which contains the system matrix
//              A(i,j) = A[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_U       number of velocity unknowns
//    N_P       number of presure unknowns
//    LDA       leading dimension of matrix a
{
  int i,j,k,l,m, row;
  int N_Eqn;
  double pp, dp, ffp, tmp;
  double Ai[MaxN_BaseFunctions2D];

  int ii, jj;

  N_Eqn = 2*N_U+N_P;
  row = 2*N_U;

  /* for(ii=0;ii<N_Eqn;ii++)
  {
    cout << ii;
    for(jj=0;jj<N_Eqn;jj++)
      cout << setw(8) << A(ii,jj);
    cout << endl;
  }
  cout << endl;*/


  if(N_P==1)
  {
    for(i=0;i<N_U;i++)
      Ai[i] = 1/AT(i,i);
  
    dp = 0;
    ffp = b[row];
    for(i=0,j=N_U; i<N_U; i++,j++)
    {
      tmp = Ai[i];
      dp  -= tmp * (AT(row, i) * AT(i, row) + AT(row, j) * AT(j, row));
      ffp -= tmp * (AT(row, i) * b[i] + AT(row, j) * b[j]);
    }
  
    pp = ffp / dp;
    b[row] = pp;
  
    for(i=0, j=N_U; i<N_U; i++, j++)
    {
      tmp = Ai[i];
      b[i] = tmp * (b[i] - AT(i, row) * pp);
      b[j] = tmp * (b[j] - AT(j, row) * pp);
    }
  }
  else
    SolveLinearSystemLapack(a, b, N_Eqn, N_Eqn);

}
// determine L2 and H1 error
void L1Int(int N_Points, double *X, double *Y, double *AbsDetjk, 
                double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *Loc)
{
  int i;
  double *deriv, *exactval, w, t;

  Loc[0] = 0.0;
  Loc[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    w = Weights[i]*AbsDetjk[i];

    // int(f)
    t = deriv[0];
    Loc[0] += w*t;

    // int(1)
    Loc[1] += w;
  } // endfor i
}

/** prolongate */
void Prolongate(TFESpace2D *CoarseSpace, 
        TFESpace2D *FineSpace, double *CoarseFunction, 
        double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int *DOF, Index;
  double *entry;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  //cout << "N_FineCells: " << N_FineCells << endl;
  //cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(FineFunction, 0, SizeOfDouble*N_FineDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      BaseFunctions->ChangeBF(CoarseColl, parent, Val);

      for(j=0;j<N_Children;j++)
      {
        // cout << "child: " << j << endl;
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do prolongation
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          entry = QQ+k*MaxN_BaseFunctions2D;
          for(l=0;l<N_Coarse;l++)
          {
            // s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
            s += entry[l] * Val[l];
            // cout << k << " " << l << " " << entry[l] << endl;
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val2);

        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
          FineFunction[Index] += Val2[k];
          aux[Index] += 1;
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      // do prolongation
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      BaseFunctions->ChangeBF(CoarseColl, cell, Val);

      for(k=0;k<N_Fine;k++)
      {
        s = 0;
        for(l=0;l<N_Coarse;l++)
        {
          s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase2D::GetBaseFunct2D(FineBF)
                      ->ChangeBF(FineColl, cell, Val2);
      for(k=0;k<N_Fine;k++)
      {
        Index = FineDOF[k];
        FineFunction[Index] += Val2[k];
        aux[Index] += 1;
      }
    } // endelse
  } // endfor i

  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }
}

void Prolongate(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int *DOF, Index;
  double *entry;
  int CoarseOffset, FineOffset, IFunct;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(FineFunction, 0, SizeOfDouble*N_Functions*N_FineDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }


  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do prolongation
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          CoarseOffset = IFunct*N_CoarseDOFs;
          FineOffset = IFunct*N_FineDOFs;

          for(l=0;l<N_Coarse;l++)
            Val[l] = CoarseFunction[CoarseOffset+CoarseDOF[l]];

          BaseFunctions->ChangeBF(CoarseColl, parent, Val);

          for(k=0;k<N_Fine;k++)
          {
            s = 0;
            entry = QQ+k*MaxN_BaseFunctions2D;
            for(l=0;l<N_Coarse;l++)
            {
              // s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
              s += entry[l] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val2);

          for(k=0;k<N_Fine;k++)
          {
            Index = FineDOF[k];
            FineFunction[FineOffset+Index] += Val2[k];
            aux[Index] += 1;
          } // endfor k
        } // endfor IFunct
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      // do prolongation
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        CoarseOffset = IFunct*N_CoarseDOFs;
        FineOffset = IFunct*N_FineDOFs;

        for(l=0;l<N_Coarse;l++)
          Val[l] = CoarseFunction[CoarseOffset + CoarseDOF[l]];

        BaseFunctions->ChangeBF(CoarseColl, cell, Val);

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          for(l=0;l<N_Coarse;l++)
          {
            s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val2);

        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
          FineFunction[FineOffset + Index] += Val2[k];
          aux[Index] += 1;
        } // endfor k
      } // endfor IFunct
    } // endelse
  } // endfor i

  Dscal(N_FineDOFs, 1.0/N_Functions, aux);
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    FineOffset = IFunct*N_FineDOFs;
    for(i=0;i<N_FineDOFs;i++)
    {
      FineFunction[FineOffset + i] /= aux[i];
    }
  }
}

/** defect restriction from level+1 to level */
void DefectRestriction(TFESpace2D *CoarseSpace,
        TFESpace2D *FineSpace, double *CoarseFunction,
        double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int *DOF, Index;
  double *entry;
  
  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  //cout << "N_CoarseCells: " << N_CoarseCells << endl;
  
  //for parallel, use the aux array from the input
#ifndef _MPI
  memset(aux, 0, SizeOfDouble*N_FineDOFs);
#endif    
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
#ifndef _MPI
    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE2D(i, cell);
    FineElement = TFEDatabase2D::GetFE2D(FineId);
    FineBF = FineElement->GetBaseFunct2D_ID();
    N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
      aux[DOF[j]] += 1;
#endif    
  }

  // modify fine function values, will be repaired at end
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions2D+k] * Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          Index = CoarseDOF[k];
          CoarseFunction[Index] += Val2[k];
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[l*MaxN_BaseFunctions2D+k]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
        CoarseFunction[CoarseDOF[k]] += Val2[k];
    } // endelse
  } // endfor i

  // repair fine function values since they are modified at beginning
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] *= aux[i];
  }
}

/** defect restriction from level+1 to level */
void DefectRestriction(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int *DOF, Index;
  double *entry;
  int FineOffset, CoarseOffset, IFunct;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs*N_Functions);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);

    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE2D(i, cell);
    FineElement = TFEDatabase2D::GetFE2D(FineId);
    FineBF = FineElement->GetBaseFunct2D_ID();
    N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
      aux[DOF[j]] += 1;
  }

  // modify fine function values, will be repaired at end
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    FineOffset = IFunct*N_FineDOFs;
    for(i=0;i<N_FineDOFs;i++)
    {
      FineFunction[FineOffset + i] /= aux[i];
    }
  } // endfor IFunct

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          CoarseOffset = IFunct*N_CoarseDOFs;

          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase2D::GetBaseFunct2D(FineBF)
                            ->ChangeBF(FineColl, cell, Val);

          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[l*MaxN_BaseFunctions2D+k] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                          ->ChangeBF(CoarseColl, parent, Val2);

          for(k=0;k<N_Coarse;k++)
          {
            Index = CoarseDOF[k];
            CoarseFunction[CoarseOffset + Index] += Val2[k];
          } // endfor k
        } // endfor IFunct
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions2D+k]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, cell, Val2);

        for(k=0;k<N_Coarse;k++)
          CoarseFunction[CoarseOffset + CoarseDOF[k]] += Val2[k];
      } // endfor IFunct
    } // endelse
  } // endfor i

  // repair fine function values since they are modified at beginning
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    FineOffset = IFunct*N_FineDOFs;
    for(i=0;i<N_FineDOFs;i++)
    {
      FineFunction[FineOffset + i] *= aux[i];
    }
  } // endfor IFunct
}

/** function restriction from level+1 to level */
void RestrictFunction(TFESpace2D *CoarseSpace, 
    TFESpace2D *FineSpace,
    double *CoarseFunction, double *FineFunction,
    double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int *DOF, Index;
  double *entry;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_CoarseDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    // cout << "i= " << i << "    ";
    // cout << "k= " << k << endl;
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, MaxN_BaseFunctions2D*SizeOfDouble);

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetRestrictionMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions2D+l] * Val[l];
          } // endfor l
          Val2[k] += s;
        } // endfor k
      } // endfor j

      TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                      ->ChangeBF(CoarseColl, parent, Val2);

      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        aux[l] += 1;
        CoarseFunction[l] += Val2[k];
      } // endfor k
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetRestrictionMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase2D::GetBaseFunct2D(FineBF)
                      ->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        CoarseFunction[l] += Val2[k];
        aux[l] += 1;
      } // endfor k
    } // endelse
  } // endfor i

  for(i=0;i<N_CoarseDOFs;i++)
    CoarseFunction[i] /= aux[i];

/*
  for(i=0;i<N_CoarseDOFs;i++)
    cout << "CoarseFunction[" << i << "]: " << CoarseFunction[i] << endl;

  for(i=0;i<N_FineDOFs;i++)
    cout << "FineFunction[" << i << "]: " << FineFunction[i] << endl;
*/

} // RestrictFunction

/** function restriction from level+1 to level */
void RestrictFunction(TFESpace2D *CoarseSpace, TFESpace2D *FineSpace,
    int N_Functions,
    double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double *Val2;
  int *DOF, Index;
  double *entry;
  int FineOffset, CoarseOffset, Offset, IFunct;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_CoarseDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs*N_Functions);
  Val2 = new double[MaxN_BaseFunctions2D*N_Functions];

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    // cout << "i= " << i << "    ";
    // cout << "k= " << k << endl;
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, N_Functions*MaxN_BaseFunctions2D*SizeOfDouble);

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetRestrictionMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

          Offset = IFunct*N_Coarse;
          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[k*MaxN_BaseFunctions2D+l] * Val[l];
            } // endfor l
            Val2[Offset + k] += s;
          } // endfor k
        } // endfor IFunct
      } // endfor j

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        CoarseOffset = IFunct*N_CoarseDOFs;
        Offset = IFunct*N_Coarse;

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          l=CoarseDOF[k];
          aux[l] += 1;
          CoarseFunction[CoarseOffset + l] += Val2[Offset + k];
        } // endfor k
      } // endfor IFunct
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetRestrictionMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, cell, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          l=CoarseDOF[k];
          CoarseFunction[CoarseOffset + l] += Val2[k];
          aux[l] += 1;
        } // endfor k
      } // endfor IFunct

    } // endelse
  } // endfor i

  Dscal(N_CoarseDOFs, 1.0/N_Functions, aux);
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    CoarseOffset = IFunct*N_CoarseDOFs;
    for(i=0;i<N_CoarseDOFs;i++)
      CoarseFunction[CoarseOffset + i] /= aux[i];
  } // endfor IFunct

  delete Val2;

} // RestrictFunction

/** project vector v into L20 */
void IntoL20Vector2D(double *v, int Length, int order)
{
  double s;
  int i;

  switch(order)
  {
    case -11:
      s=0;
      for(i=0;i<Length;i+=3)
        s += v[i];

      s /= Length;
      s *= 3;

      for(i=0;i<Length;i+=3)
        v[i] -= s;
      break;

    case -12:
      s=0;
      for(i=0;i<Length;i+=6)
        s += v[i];

      s /= Length;
      s *= 6;

      for(i=0;i<Length;i+=6)
        v[i] -= s;
     
      break;

    case -13:
      s=0;
      for(i=0;i<Length;i+=10)
        s += v[i];

      s /= Length;
      s *= 10;

      for(i=0;i<Length;i+=10)
        v[i] -= s;
      
        break;

    case -14:
      s=0;
      for(i=0;i<Length;i+=15)
        s += v[i];

      s /= Length;
      s *= 15;

      for(i=0;i<Length;i+=15)
        v[i] -= s;
      
        break;

    case 201:
      break;

    default :
      s=0;
      for(i=0;i<Length;i++)
        s += v[i];

      s /= Length;
//      OutPut("VECTOR " << Length << " " << s << endl);
      for(i=0;i<Length;i++)
        v[i] -= s;
      break;

  }
} // IntoL20

/** project function v into L20 */
void IntoL20FEFunction_OLD(double *v, int Length, TFESpace2D *FESpace,
                       int velocity_space, int pressure_space)
{
  double s;
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs2D], *N_BaseFunct;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE2D *ele;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double *Derivatives[MaxN_QuadPoints_2D], der[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, number;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1;

  if (pressure_space==-4711)
  {
    switch (velocity_space)
    {
        case 1:
          number = 0;
          break;
        case 2:
        case 3:
        case 4:
        case 5:
          number = velocity_space;
          break;
        case 12:
        case 13:
        case 14:
          number = -velocity_space+1;
          break;
        case 22:
        case 23:
          number = -velocity_space+11;
          break;

        case -1:
          number = -TDatabase::ParamDB->VELOCITY_SPACE -1;
        break;

        case -2:
        case -3:
        case -4:
        case -5:
          number = TDatabase::ParamDB->VELOCITY_SPACE - 9;
        break;
    }
  }
  else
    number = pressure_space;

//  OutPut("NUMBER " << number << " " << Length << " ");

  switch(number)
  {
    // pw constant
    case 0:
      s=0;
      for(i=0;i<Length;i++)
        s += v[i];
//      OutPut(s << endl);
      s /= Length;
      for(i=0;i<Length;i++)
        v[i] -= s;
      break;
      
    // pw linear discontinuous   
    case -11:
      s=0;
      for(i=0;i<Length;i+=3)
        s += v[i];

      s /= Length;
      s *= 3;

      for(i=0;i<Length;i+=3)
        v[i] -= s;
      break;

    // pw quadratics discontinuous   
    case -12:
      s=0;
      for(i=0;i<Length;i+=6)
        s += v[i];

      s /= Length;
      s *= 6;

      for(i=0;i<Length;i+=6)
        v[i] -= s;
      break;

    // pw cubics discontinuous   
    case -13:
      s=0;
      for(i=0;i<Length;i+=10)
        s += v[i];

      s /= Length;
      s *= 10;

      for(i=0;i<Length;i+=10)
        v[i] -= s;
      
      break;
      
    // discontinuous, pw P4
    case -14:
      s=0;
      for(i=0;i<Length;i+=15)
        s += v[i];

      s /= Length;
      s *= 15;

      for(i=0;i<Length;i+=15)
        v[i] -= s;
      
      break;
      
    // pw polynomials, continuous
    // with and without bubbles   
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 12:
    case 13:
    case 14:
    case 22:
    case 23:
      BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
      N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

      for(i=0;i<MaxN_QuadPoints_2D;i++)
        Derivatives[i] = der+i;
    
      SecondDer[0] = FALSE;
    
      GlobalNumbers = FESpace->GetGlobalNumbers();
      BeginIndex = FESpace->GetBeginIndex();
    
      error0 = 0.0;
      error1 = 0.0;
    
    // ########################################################################
    // loop over all cells
    // ########################################################################
      Coll = FESpace->GetCollection(); // all spaces use same Coll
      N_Cells = Coll->GetN_Cells();
      for(i=0;i<N_Cells;i++)
      {
        cell = Coll->GetCell(i);
    
        hK = cell->GetDiameter();
    
        // ####################################################################
        // find local used elements on this cell
        // ####################################################################
        N_LocalUsedElements = 1;
        CurrentElement = FESpace->GetFE2D(i, cell);
        LocalUsedElements[0] = CurrentElement;
    
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                             Coll, cell, SecondDer,
                             N_Points, xi, eta, weights, X, Y, AbsDetjk);
    
        // calculate all needed derivatives of this FE function
        BaseFunct = BaseFuncts[CurrentElement];
        N_ = N_BaseFunct[CurrentElement];
    
        DOF = GlobalNumbers + BeginIndex[i];
        for(l=0;l<N_;l++)
          FEFunctValues[l] = v[DOF[l]];
    
        OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
        for(j=0;j<N_Points;j++)
        {
          Orig = OrigFEValues[j];
          value = 0;
          for(l=0;l<N_;l++)
          {
            value += FEFunctValues[l] * Orig[l];
          } // endfor l
          Derivatives[j][0] = value;
        } // endfor j
    
        L1Int(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives, 
                  ExactVal, AuxArray, LocError);
    
        error0 += LocError[0];
        error1 += LocError[1];
    
      } // endfor i
    
      s = error0/error1;

      for(i=0;i<Length;i++)
        v[i] -= s; 
    
    break;

    default:
      cout << "The L^2_0 projection does not ";
      cout << "work for this type of elements" << endl;
      exit(-1);
  }
} // IntoL20Function

/** project function v into L20 */
void IntoL20FEFunction(double *v, int Length, TFESpace2D *FESpace,
                       int velocity_space, int pressure_space
#ifdef _MPI
                        , MPI_Comm comm
#endif
                      )
{
  double s;
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs2D], *N_BaseFunct;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE2D *ele;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double *Derivatives[MaxN_QuadPoints_2D], der[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, number;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1, temp;
  double *interpol;
  TNodalFunctional2D *nf;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];

  #ifdef _MPI
  int ID, rank;
  MPI_Comm_rank(comm, &rank);
  #endif

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  for(i=0;i<MaxN_QuadPoints_2D;i++)
    Derivatives[i] = der+i;

  SecondDer[0] = FALSE;

  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();

  interpol = new double[Length];

  error0 = 0.0;
  error1 = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = FESpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    #ifdef _MPI
    ID = cell->GetSubDomainNo();
    if(ID!=rank) 
      continue; // this cell is a halo cell
    #endif
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 1;
    CurrentElement = FESpace->GetFE2D(i, cell);
    LocalUsedElements[0] = CurrentElement;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // calculate all needed derivatives of this FE function
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = v[DOF[l]];

    OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
      {
        value += FEFunctValues[l] * Orig[l];
      } // endfor l
      Derivatives[j][0] = value;
    } // endfor j

    L1Int(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    error0 += LocError[0];
    error1 += LocError[1];

  } // endfor i

  #ifdef _MPI
  temp = error0;
  MPI_Allreduce(&temp, &error0, 1, MPI_DOUBLE, MPI_SUM, comm);
  temp = error1;
  MPI_Allreduce(&temp, &error1, 1, MPI_DOUBLE, MPI_SUM, comm);
  #endif

  s = error0/error1;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    CurrentElement = FESpace->GetFE2D(i, cell);
    N_ = N_BaseFunct[CurrentElement];
    nf = TFEDatabase2D::GetNodalFunctional2DFromFE2D(CurrentElement);
    nf->GetPointsForAll(N_Points, xi, eta);
    for(j=0;j<N_Points;j++)
      PointValues[j] = s;
    nf->GetAllFunctionals(Coll, cell, PointValues, FunctionalValues);
    DOF = GlobalNumbers+BeginIndex[i];
    for(j=0;j<N_;j++)
      interpol[DOF[j]] = FunctionalValues[j];
  } // endfor i

  for(i=0;i<Length;i++)
    v[i] -= interpol[i];

  delete [] interpol;
} // IntoL20Function

/** Navier--Stokes type 1 (NSTYPE==1) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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
      //if (i>=N_Active)
//	  OutPut(i << " " << index << " matv " << value << endl);
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

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
      }
    } // endfor j
    r3[i] = s;
  } // endfor i
}
void Defect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}
 
/** Navier--Stokes type 2 (NSTYPE==2) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T,
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

  return;
}
void MatVect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], B[3], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T,
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
}

void Defect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], B[2], B[3], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();



  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** Navier--Stokes type 3 (NSTYPE==3) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *A12Entries, *A21Entries, *A22Entries;
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

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], B[0], B[1], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12,
                   TSquareMatrix *A21, TSquareMatrix *A22, 
                   TMatrix *B1, TMatrix *B2,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *A12Entries, *A21Entries, *A22Entries;
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

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
      }
    } // endfor j
    r3[i] = s;
  } // endfor i
}
void Defect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], B[0], B[1], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** Navier--Stokes type 4 (NSTYPE==4) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21, 
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
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

  return;
}
void MatVect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
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
}
void Defect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{ 
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}


/** Navier--Stokes type 14 (NSTYPE==14) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
/** with equal order interpolation */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21, 
                    TSquareMatrix *A22, TSquareMatrix *C, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2,value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries, *CEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

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
      index = CKCol[j];
      value1 = CEntries[j];
      s += value1 * p[index];
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_EquOrd_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  // A[4] is the pressure-pressure matrix  
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], B[0], B[1], B[2], B[3], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TSquareMatrix *C, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries, *CEntries;
  int N_Active;

  // Aij and C have the same structure
  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

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

  // velocity-velocity block, active dofs
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

  // velocity-velocity block, non-active dofs, like Dirichlet
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

  // divergence matrix with velocity
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

  // gradient matrix with pressure
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

  // pressure-pressure block
  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      index = CKCol[j];
      value1 = CEntries[j];
      value = p[index];
      s += value1 * value;
    }
    r3[i] -= s;
   } // endfor i
}

void Defect_EquOrd_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;
  // A[4] is the pressure-pressure matrix
  CoupledDefect(A[0], A[1], A[2], A[3], A[4], B[0], B[1], B[2], B[3], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

void MatVect_Raviart_Thomas_NSE4(TSquareMatrix **A, TMatrix **B, double *x, 
    double *y)
{
  CoupledMatVect(A[0],B[0],x,y);
  return;
}

/**  Darcy Raviart-Thomas
 ( A B' )
 ( B 0  )
 block A: flux x flux
*/
void CoupledMatVect(TSquareMatrix *A, TMatrix *B, double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value;
  double *u, *p;
  double *v, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *BEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();

  BEntries = B->GetEntries();
  
  N_UDOF = A->GetN_Rows();
  N_PDOF = B->GetN_Rows();

  u = x;
  p = u+N_UDOF;

  v = y;
  q = v+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u[index];
    }
    v[i] = s;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s += value * u[index];

      if(index<N_Active)
      {
        t = p[i];
        v[index] += value * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B, 
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value;
  double *u, *p;
  double *v, *q;
  double *r1, *r2;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *BEntries;
  int N_Active;
  
  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();

  BEntries = B->GetEntries();
  
  N_UDOF = A->GetN_Rows();
  N_PDOF = B->GetN_Rows();

  u = x;
  p  = u+N_UDOF;

  v = b;
  q  = v+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  
  // ( r1 ) = ( v ) _ ( A B' ) ( u )
  // ( r2 )   ( q )   ( B 0  ) ( p )
  
  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u[index];
    }
    r1[i] = s;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s -= value * u[index];

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value * t;
      }
    } // endfor j
    r2[i] = s;
  } // endfor i
}

void Defect_Raviart_Thomas_NSE4(TSquareMatrix **A, TMatrix **B, 
                                double *x, double *b, double *r)
{
  CoupledDefect(A[0], B[0], x, b, r);
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    int N_UDOF,N_PDOF;
    N_UDOF = A[0]->GetN_Rows();
    N_PDOF = B[0]->GetN_Rows();
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  }
}



/** Navier--Stokes type 5 (NSTYPE==5) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, double *x, double *y)
{
  int i,j,k,l,m;
  TSquareMatrixNSE2D *NSE_matrix;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *RowPtr, *KCol;
  double *Entries;
  int N_UDOF, N_PDOF, N_Active;
  TFESpace2D *USpace;
  TCollection *coll;
  double *u1, *u2, *p, *v1, *v2, *q, *r1, *r2, *r3;
  double s, t, value1, value2;
  double value11, value12, value21, value22;
  int index;
  int N_Cells, N_Joints;
  TBaseCell *cell;
  int *GlobalNumbers, *BeginIndex, *DOF, *LocJb;
  int **JointDOFs, *EdgeDOF;
  double *LocAlpha;
  FE2D UElement;
  TFE2D *ele;
  TFEDesc2D *FEDesc_Obj;
  int N_U;
  
  NSE_matrix = (TSquareMatrixNSE2D*)A;

  RowPtr = NSE_matrix->GetRowPtr();
  KCol = NSE_matrix->GetKCol();
  Entries = NSE_matrix->GetEntries();

  BeginJb = NSE_matrix->GetBeginJb();
  jb = NSE_matrix->GetJb();
  N_DOFperJoint = NSE_matrix->GetN_DOFperJoint();
  Alpha = NSE_matrix->GetAlpha();

  USpace = NSE_matrix->GetFESpace();
  GlobalNumbers = USpace->GetGlobalNumbers();
  BeginIndex = USpace->GetBeginIndex();

  coll = USpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  N_UDOF = USpace->GetN_DegreesOfFreedom();
  N_Active = USpace->GetActiveBound();
  N_PDOF = N_Cells;

  u1 = x;
  u2 = u1 + N_UDOF;
  p  = u2 + N_UDOF;

  r1 = y;
  r2 = r1 + N_UDOF;
  r3 = r2 + N_UDOF;

  memset(r3, 0, N_PDOF*SizeOfDouble);

  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = RowPtr[i+1];
    if(k == j)
    {
      // bubble row => no defect
      r1[i] = 0;
      r2[i] = 0;
    }
    else
    {
      // skeleton dof
      s = 0;
      t = 0;
      for(;j<k;j++)
      {
        index = KCol[j];
        value11 = Entries[4*j + 0];
        value12 = Entries[4*j + 1];
        value21 = Entries[4*j + 2];
        value22 = Entries[4*j + 3];

        value1 = u1[index];
        value2 = u2[index];

        s += value11*value1 + value12*value2;
        t += value21*value1 + value22*value2;
      } // endfor j
      r1[i] = s;
      r2[i] = t;
    } // endelse
  } // endfor i

  // Dirichlet rows
  for(i=N_Active;i<N_UDOF;i++)
  {
    r1[i] = u1[i];
    r2[i] = u2[i];
  }

  // use data from B-blocks
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(UElement);
    FEDesc_Obj = ele->GetFEDesc2D();
    JointDOFs = FEDesc_Obj->GetJointDOF();
    N_U = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    //  B entries
    LocJb = jb+BeginJb[i];

    // NOTE: Alpha stores the negative B block entry
    for(j=0;j<N_Joints;j++)
    {
      // get local data
      LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
      EdgeDOF = JointDOFs[j]; 
        
      l = LocJb[j];
      if(l<N_U)
      {
        // jb in first component
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        m = DOF[l];
        r3[i] -= LocAlpha[k]*u1[m];
        if(m<N_Active)
          r1[m] -= LocAlpha[k]*p[i];
      }
      else
      {
        // jb in second component
        l -= N_U;
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        m = DOF[l];
        r3[i] -= LocAlpha[k+N_DOFperJoint]*u2[m];
        if(m<N_Active)
          r2[m] -= LocAlpha[k+N_DOFperJoint]*p[i];
      }
    } // endfor j
  } // endfor i
}

void MatVect_NSE5(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, double *x, double *b, double *r)
{
  int i,j,k,l,m;
  TSquareMatrixNSE2D *NSE_matrix;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *RowPtr, *KCol;
  double *Entries;
  int N_UDOF, N_PDOF, N_Active;
  TFESpace2D *USpace;
  TCollection *coll;
  double *u1, *u2, *p, *v1, *v2, *q, *r1, *r2, *r3;
  double s, t, value1, value2;
  double value11, value12, value21, value22;
  int index;
  int N_Cells, N_Joints;
  TBaseCell *cell;
  int *GlobalNumbers, *BeginIndex, *DOF, *LocJb;
  int **JointDOFs, *EdgeDOF;
  double *LocAlpha;
  FE2D UElement;
  TFE2D *ele;
  TFEDesc2D *FEDesc_Obj;
  int N_U;
  
  NSE_matrix = (TSquareMatrixNSE2D*)A;

  RowPtr = NSE_matrix->GetRowPtr();
  KCol = NSE_matrix->GetKCol();
  Entries = NSE_matrix->GetEntries();

  BeginJb = NSE_matrix->GetBeginJb();
  jb = NSE_matrix->GetJb();
  N_DOFperJoint = NSE_matrix->GetN_DOFperJoint();
  Alpha = NSE_matrix->GetAlpha();

  USpace = NSE_matrix->GetFESpace();
  GlobalNumbers = USpace->GetGlobalNumbers();
  BeginIndex = USpace->GetBeginIndex();

  coll = USpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  N_UDOF = USpace->GetN_DegreesOfFreedom();
  N_Active = USpace->GetActiveBound();
  N_PDOF = N_Cells;

  u1 = x;
  u2 = u1 + N_UDOF;
  p  = u2 + N_UDOF;

  v1 = b;
  v2 = v1 + N_UDOF;
  q  = v2 + N_UDOF;

  r1 = r;
  r2 = r1 + N_UDOF;
  r3 = r2 + N_UDOF;

  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = RowPtr[i+1];
    if(k == j)
    {
      // bubble row => no defect
      r1[i] = 0;
      r2[i] = 0;
    }
    else
    {
      // skeleton dof
      s = v1[i];
      t = v2[i];
      for(;j<k;j++)
      {
        index = KCol[j];
        value11 = Entries[4*j + 0];
        value12 = Entries[4*j + 1];
        value21 = Entries[4*j + 2];
        value22 = Entries[4*j + 3];

        value1 = u1[index];
        value2 = u2[index];

        s -= value11*value1 + value12*value2;
        t -= value21*value1 + value22*value2;
      } // endfor j
      r1[i] = s;
      r2[i] = t;
    } // endelse
  } // endfor i

  // Dirichlet rows
  for(i=N_Active;i<N_UDOF;i++)
  {
    r1[i] = v1[i] - u1[i];
    r2[i] = v2[i] - u2[i];
  }

  memcpy(r3, q, N_PDOF*SizeOfDouble);

  // use data from B-blocks
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(UElement);
    FEDesc_Obj = ele->GetFEDesc2D();
    JointDOFs = FEDesc_Obj->GetJointDOF();
    N_U = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    //  B entries
    LocJb = jb+BeginJb[i];

    // NOTE: Alpha stores the negative B block entry
    for(j=0;j<N_Joints;j++)
    {
      // get local data
      LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
      EdgeDOF = JointDOFs[j]; 
        
      l = LocJb[j];
      if(l<N_U)
      {
        // jb in first component
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        m = DOF[l];
        r3[i] += LocAlpha[k]*u1[m];
        if(m<N_Active)
          r1[m] += LocAlpha[k]*p[i];

        // cout << " " << DOF[l] << " " << i;
        // cout << " " << LocAlpha[k] << endl;
      }
      else
      {
        // jb in second component
        l -= N_U;
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        m = DOF[l];
        r3[i] += LocAlpha[k+N_DOFperJoint]*u2[m];
        if(m<N_Active)
          r2[m] += LocAlpha[k+N_DOFperJoint]*p[i];

        // cout << " " << DOF[l]+N_UDOF << "(" << DOF[l] << ")";
        // cout << " " << i;
        // cout << " " << LocAlpha[k+N_DOFperJoint] << endl;
      }
    } // endfor j
  } // endfor i
} // end CoupledDefect

void Defect_NSE5(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = ((TSquareMatrixNSE2D *)A[0])->GetFESpace()->GetN_Cells();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,0);
  return;
}

/** Convection-diffusion problem with VMM */
void MatVectCD_VMM(TSquareMatrix *A, TMatrix *B, TMatrix *C, TSquareMatrix *D, 
                   double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *p;
  double *v1, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *CRowPtr, *DRowPtr, *CKCol, *DKCol;
  double *AEntries, *BEntries, *CEntries, *DEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();
  BEntries = B->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p  = u1+N_UDOF;

  v1 = y;
  q  = v1+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
    }
    v1[i] = s;
  } // endfor i

  // block B 
  j = BRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s += value * p[index];
    }
    v1[i] += s;
  } // endfor i

  // block C
  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      index = CKCol[j];
      value1 = CEntries[j];
      s += value1 * u1[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p[index];
    } // endfor j
    q[i] += s;
  } // endfor i
  return;
}
void MatVect_CD_VMM(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  MatVectCD_VMM(A[0], B[0], B[1], A[1], x, y);
  return;
}

/** r := b - A * x */
void DefectCD_VMM(TSquareMatrix *A, TMatrix *B, TMatrix *C, TSquareMatrix *D, 
                   double *x, double *y, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *p;
  double *v1, *q, *r1, *r2;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *CRowPtr, *DRowPtr, *CKCol, *DKCol;
  double *AEntries, *BEntries, *CEntries, *DEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();
  BEntries = B->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p  = u1+N_UDOF;

  v1 = y;
  q  = v1+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
    }
    r1[i] = s;
  } // endfor i

  // block B 
  j = BRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s += value * p[index];
    }
    r1[i] -= s;
  } // endfor i

  // block C
  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      index = CKCol[j];
      value1 = CEntries[j];
      s -= value1 * u1[index];
    } // endfor j
    r2[i] = s;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p[index];
    } // endfor j
    r2[i] -= s;
  } // endfor i
  return;
}

void Defect_CD_VMM(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  DefectCD_VMM(A[0], B[0], B[1], A[1], x, b, r);
  return;
}
 
/** Convection-diffusion problem with VMM [KL02] */
void MatVectCD_VMM_KL02(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, 
                        TMatrix *C1, TMatrix *C2, TSquareMatrix *D, 
                        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *p1, *p2;
  double *v1, *q1, *q2;
  int *ARowPtr, *B1RowPtr, *AKCol, *B1KCol;
  int *C1RowPtr, *DRowPtr, *C1KCol, *DKCol;
  double *AEntries, *B1Entries, *C1Entries, *DEntries, *B2Entries, *C2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  B1RowPtr = B1->GetRowPtr();
  B1KCol = B1->GetKCol();
  B1Entries = B1->GetEntries();

  B2Entries = B2->GetEntries();

  C1RowPtr = C1->GetRowPtr();
  C1KCol = C1->GetKCol();
  C1Entries = C1->GetEntries();

  C2Entries = C2->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p1  = u1+N_UDOF;
  p2  = p1+N_PDOF;

  v1 = y;
  q1  = v1+N_UDOF;
  q2  = q1+N_PDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
    }
    v1[i] = s;
  } // endfor i

  // block B1, B2 
  j = B1RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    t = 0;
    k = B1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = B1KCol[j];
      value = B1Entries[j];
      s += value * p1[index];
      value = B2Entries[j];
      t += value * p2[index];      
    }
    v1[i] += s+t;
  } // endfor i

  // block C1, C2
  j = C1RowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    t = 0;
    k = C1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = C1KCol[j];
      value1 = C1Entries[j];
      s += value1 * u1[index];
      value2 = C2Entries[j];
      t += value2 * u1[index];
   } // endfor j
    q1[i] = s;
    q2[i] = t;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    t = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p1[index];
      t += value1 * p2[index];
    } // endfor j
    q1[i] += s;
    q2[i] += t;
  } // endfor i
  return;
}
void MatVect_CD_VMM_KL02(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  MatVectCD_VMM_KL02(A[0], B[0], B[1], B[2], B[3], A[1], x, y);
  return;
}
void DefectCD_VMM_KL02(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, 
                        TMatrix *C1, TMatrix *C2, TSquareMatrix *D, 
                        double *x, double *y, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *p1, *p2;
  double *v1, *q1, *q2, *r1, *r2, *r3;
  int *ARowPtr, *B1RowPtr, *AKCol, *B1KCol;
  int *C1RowPtr, *DRowPtr, *C1KCol, *DKCol;
  double *AEntries, *B1Entries, *C1Entries, *DEntries, *B2Entries, *C2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  B1RowPtr = B1->GetRowPtr();
  B1KCol = B1->GetKCol();
  B1Entries = B1->GetEntries();

  B2Entries = B2->GetEntries();

  C1RowPtr = C1->GetRowPtr();
  C1KCol = C1->GetKCol();
  C1Entries = C1->GetEntries();

  C2Entries = C2->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p1  = u1+N_UDOF;
  p2  = p1+N_PDOF;

  v1 = y;
  q1  = v1+N_UDOF;
  q2  = q1+N_PDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_PDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
    }
    r1[i] = s;
  } // endfor i

  // block B1, B2 
  j = B1RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    k = B1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = B1KCol[j];
      value = B1Entries[j];
      value1 = B2Entries[j];
      s += value * p1[index] + value1 * p2[index];
    }
    r1[i] -= s;
  } // endfor i

  // block C1, C2
  j = C1RowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q1[i];
    t = q2[i];
    k = C1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = C1KCol[j];
      value1 = C1Entries[j];
      s -= value1 * u1[index];
      value2 = C2Entries[j];
      t -= value2 * u1[index];
   } // endfor j
    r2[i] = s;
    r3[i] = t;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    t = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p1[index];
      t += value1 * p2[index];
    } // endfor j
    r2[i] -= s;
    r3[i] -= t;
  } // endfor i
  return;
}

void Defect_CD_VMM_KL02(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  DefectCD_VMM_KL02(A[0], B[0], B[1], B[2], B[3], A[1], x, b, r);
  return;
}

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVectLV96(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *y, double delta)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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
 
  for(i=0;i<N_UDOF;i++) // A*u
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
  for(i=0;i<N_PDOF;i++) // B^T*p and B*u
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t; 
        v2[index] += value2 * t;
      }
    } // endfor j
    q[i] = s;               // B*u
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++) // B^T * B * u = B^T * q
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];

      if(index<N_Active)
      {
        t = q[i];
        v1[index] += value1 * t/delta;
        v2[index] += value2 * t/delta;
      }
    } // endfor j
  } // endfor i
  return;
}

/** r := b - A * x */
void CoupledDefectLV96(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *b, double *r, double delta)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q, *help;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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
  for(i=0;i<N_UDOF;i++) // -A*u
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

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
      }
    } // endfor j
    r3[i] = s;
  } // endfor i
  
  help = new double [N_PDOF];
  memset(help,0,N_PDOF*SizeOfDouble);
  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++) // B*u
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
    help[i] = s;               // help = B*u
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++) 
  {
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];

      if(index<N_Active)
      {
        t = help[i];
        r1[index] -= value1 * t/delta;
        r2[index] -= value2 * t/delta;
      }
    } // endfor j
  } // endfor i
  delete help;
}

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVectMortar(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *matrix_mortar, 
        double *x, double *y)
{
  int N_UDOF, N_PDOF,N_MORTDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p, *m1, *m2;
  double *v1, *v2, *q, *n1, *n2;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol, *MortarRowPtr, *MortarKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries, *MortarEntries;
  int N_Active;

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

  MortarRowPtr = matrix_mortar->GetRowPtr();
  MortarKCol =  matrix_mortar->GetKCol();
  MortarEntries = matrix_mortar->GetEntries();
  
  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();
  N_MORTDOF = matrix_mortar->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;
  m1 = p+N_PDOF;
  m2 = m1+N_MORTDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;
  n1 = q+N_PDOF;
  n2 = n1+N_MORTDOF;

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

  for (i=0;i<N_MORTDOF;i++)
    {
      l=MortarRowPtr[i];
      k=MortarRowPtr[i+1];
      s=t=0;
      for(j=l;j<k;j++)
        {
          index = MortarKCol[j];
          s+=MortarEntries[j]*u1[index];
          t+=MortarEntries[j]*u2[index];
        }
      n1[i] = s;
      n2[i] = t;
    } 

  for (i=0;i<N_MORTDOF;i++)
    {
      l=MortarRowPtr[i];
      k=MortarRowPtr[i+1];
      for(j=l;j<k;j++)
        {
          index = MortarKCol[j];
          v1[index]+=MortarEntries[j]*m1[i];
          v2[index]+=MortarEntries[j]*m2[i];
        }
    } 
 return;
}

#endif
