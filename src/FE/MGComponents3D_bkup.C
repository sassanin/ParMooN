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
// MGComponents3D.C
//
// Purpose:     components for multigrid in 3d
//
// Author:      Gunar Matthies          27.01.1999
//              Volker John             27.10.1999  
// 
//              parallel methods  (Sashikumaar Ganesan) 19.09.2010
// =======================================================================

#ifdef _MPI
# include "mpi.h"
#endif

#include <LinAlg.h>
#include <Database.h>
#include <FEDatabase3D.h>

#include <stdlib.h>
#include <string.h>

#ifdef __3D__

#include <SquareMatrixNSE3D.h>

#define AT(i,j) (a[(j)*LDA+(i)])
#define A(i,j) (a[(i)*LDA+(j)])

void SolveDiagonalVanka3D(double *a, double *b, int N_U, int N_P, int LDA)
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
  double Ai[8*MaxN_BaseFunctions3D];

  int ii, jj;

  double S[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double r[MaxN_BaseFunctions3D];
  double q[MaxN_BaseFunctions3D];
  int NU3;

  N_Eqn = 3*N_U+N_P;
  row = 3*N_U;

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
    for(i=0,j=N_U,k=2*N_U; i<N_U; i++,j++,k++)
    {
      tmp = Ai[i];
      dp  -= tmp * (AT(row, i) * AT(i, row) + AT(row, j) * AT(j, row)
        + AT(row, k) * AT(k, row));
      ffp -= tmp * (AT(row, i) * b[i] + AT(row, j) * b[j]+ AT(row, k) * b[k]);
    }
 
    pp = ffp / dp;
    b[row] = pp;
  
    for(i=0, j=N_U,k=2*N_U; i<N_U; i++, j++,k++)
    {
      tmp = Ai[i];
      b[i] = tmp * (b[i] - AT(i, row) * pp);
      b[j] = tmp * (b[j] - AT(j, row) * pp);
      b[k] = tmp * (b[k] - AT(k, row) * pp);
    }
  }
  else
  {
    // OutPut("SolveDiagonalVanka3D" << endl);
    //if(fabs(TDatabase::ParamDB->REACTOR_P25 - 8)<1e-10)
    {
      SolveLinearSystemLapack(a, b, N_Eqn, N_Eqn);
    }
    /* else
    {
      NU3 = 3*N_U;
      for(i=0;i<N_P;i++)
      {
        for(j=0;j<N_P;j++)
        {
          pp = 0;
          for(k=0;k<NU3;k++)
          {
            pp += A(i+NU3,k)*A(k,j+NU3)/A(k,k);
          } // endfor k
          S[i*N_P+j] = pp;
        } // endfor j
      } // endfor i

      for(i=0;i<N_P;i++)
      {
        pp = -b[NU3+i];
        for(k=0;k<NU3;k++)
        {
          pp += A(i+NU3,k)*b[k]/A(k,k);
        } // endfor k
        r[i] = pp;
      } // endfor i

      for(i=0;i<N_P;i++)
      {
        for(j=0;j<N_P;j++)
        {
          cout << setw(4) << i << setw(4) << j << setw(25) << S[i*N_P+j] << endl;
        }
      }

      SolveLinearSystem(S, r, N_P, N_P);

      for(i=0;i<N_P;i++)
        b[NU3+i] = r[i];

      for(i=0;i<NU3;i++)
      {
        pp = b[i];
        for(j=0;j<N_P;j++)
        {
          pp -= A(i,j+NU3)*r[j];
        } // endfor j
        b[i] = pp/A(i,i);
      } // endfor i

      for(i=0;i<N_Eqn;i++)
        b[i] *= TDatabase::ParamDB->REACTOR_P24;
	}*/
  }
}

// determine L2 and H1 error
void L1Int3D(int N_Points, double *X, double *Y, double *Z, 
             double *AbsDetjk, 
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

double tP=0.0,tR=0.0;
/** prolongate */
void Prolongate(TFESpace3D *CoarseSpace, 
        TFESpace3D *FineSpace, double *CoarseFunction, 
        double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
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
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int *DOF, Index;
  double *entry;
  double t1,t2;
#ifdef _MPI
  t1 = MPI_Wtime();
#else
  t1 = GetTime();
#endif
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

#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,j,k,l,cell,DOF,FineId,FineElement,FineBF,N_Fine, \
                                             parent, N_Children, CoarseNumber, CoarseId, CoarseElement, \
                                             CoarseBF, BaseFunctions, Ref, FineNumber,\
                                             QQ, FineDOF, CoarseDOF, Val, s, Val2, \
                                             Index, N_Coarse, entry)
{
#pragma omp for schedule(guided) 
#endif
  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
    
    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE3D(i, cell);
    FineElement = TFEDatabase3D::GetFE3D(FineId);
    FineBF = FineElement->GetBaseFunct3D_ID();
    N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
#ifdef _HYBRID      
      #pragma omp atomic 
#endif  
      aux[DOF[j]] += 1;
  }
  
#ifdef _HYBRID
#pragma omp for schedule(guided) 
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }


#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }
  
#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())   continue;
#endif
    
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
// #ifdef _HYBRID
// 	#pragma omp critical
// #endif
        {
         k = cell->GetClipBoard();
         cell->SetClipBoard(-2);
        }
        
#ifdef _HYBRID
	  if(k==-2) continue;
#endif
	FineNumber = -(k+10);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
                (CoarseId, Ref, FineId, j);
	}

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          entry = QQ+k*MaxN_BaseFunctions3D;
          for(l=0;l<N_Coarse;l++)
          {
            // s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
            s += entry[l] * Val[l];
            // cout << k << " " << l << " " << entry[l] << endl;
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val2);
			
        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
	 {
#ifdef _HYBRID
       #pragma omp atomic
#endif  
          FineFunction[Index] += Val2[k];
	 }
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
              (CoarseId, Ref, FineId, 0);
	}

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
          s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase3D::GetBaseFunct3D(FineBF)
                      ->ChangeBF(FineColl, cell, Val2);
      
      for(k=0;k<N_Fine;k++)
      {
        Index = FineDOF[k];
	{
#ifdef _HYBRID
       #pragma omp atomic
#endif  
         FineFunction[Index] += Val2[k];
	}
      }
    } // endelse
  } // endfor i

#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif 
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }
#ifdef _HYBRID
} 
#endif

#ifdef _MPI
  t2 = MPI_Wtime();
#else
  t2 = GetTime();
#endif
  tP += (t2-t1);
}

void Prolongate(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
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
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do prolongation
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetProlongationMatrix3D 
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
            entry = QQ+k*MaxN_BaseFunctions3D;
            for(l=0;l<N_Coarse;l++)
            {
              // s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
              s += entry[l] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase3D::GetBaseFunct3D(FineBF)
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
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      // do prolongation
      QQ = TFEDatabase3D::GetProlongationMatrix3D 
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
            s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(FineBF)
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
void DefectRestriction(TFESpace3D *CoarseSpace,
        TFESpace3D *FineSpace, double *CoarseFunction,
        double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
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
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int *DOF, Index;
  double *entry;
  double t1,t2;
#ifdef _MPI
  t1 = MPI_Wtime();
#else
  t1 = GetTime();
#endif
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

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs);
  
  
  
  
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,j,k,l,cell,DOF,FineId,FineElement,FineBF,N_Fine, \
                                             parent, N_Children, CoarseNumber, CoarseId, CoarseElement, \
                                             CoarseBF, BaseFunctions, Ref, FineNumber,\
                                             QQ, FineDOF, CoarseDOF, Val, s, Val2, \
                                             Index, N_Coarse)
{
#pragma omp for schedule(guided) 
#endif
//   set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);

    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE3D(i, cell);
    FineElement = TFEDatabase3D::GetFE3D(FineId);
    FineBF = FineElement->GetBaseFunct3D_ID();
    N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
#ifdef _HYBRID      
      #pragma omp atomic 
#endif      
      aux[DOF[j]] += 1;
  }

#ifdef _HYBRID
#pragma omp for schedule(guided) nowait 
#endif
  // modify fine function values, will be repaired at end
  for(i=0;i<N_FineDOFs;i++)
  {    
    FineFunction[i] /= aux[i];
  }
   
#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())	continue;
#endif
    
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();
  
      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
// #ifdef _HYBRID
// 	#pragma omp critical
// #endif
	{
         k = cell->GetClipBoard();
         cell->SetClipBoard(-2);
	}
	
#ifdef _HYBRID
	  if(k==-2) continue;
#endif
	  
	FineNumber = -(k+10);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
                (CoarseId, Ref, FineId, j);
	}

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions3D+k] * Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                        ->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          Index = CoarseDOF[k];
#ifdef _HYBRID
       #pragma omp atomic
#endif
           CoarseFunction[Index] += Val2[k];
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;
#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
              (CoarseId, Ref, FineId, 0);
	}

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[l*MaxN_BaseFunctions3D+k]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
      {
#ifdef _HYBRID
       #pragma omp atomic
#endif
	CoarseFunction[CoarseDOF[k]] += Val2[k];
      }
    } // endelse
  } // endfor i

  // repair fine function values since they are modified at beginning
#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif  
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] *= aux[i];
  }
#ifdef _HYBRID
} 
#endif
  
#ifdef _MPI
  t2 = MPI_Wtime();
#else
  t2 = GetTime();
#endif
  tR+=(t2-t1);
  
}

/** defect restriction from level+1 to level */
void DefectRestriction(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
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
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
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
    FineId = FineSpace->GetFE3D(i, cell);
    FineElement = TFEDatabase3D::GetFE3D(FineId);
    FineBF = FineElement->GetBaseFunct3D_ID();
    N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetProlongationMatrix3D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          CoarseOffset = IFunct*N_CoarseDOFs;

          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase3D::GetBaseFunct3D(FineBF)
                            ->ChangeBF(FineColl, cell, Val);

          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[l*MaxN_BaseFunctions3D+k] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
      QQ = TFEDatabase3D::GetProlongationMatrix3D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions3D+k]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
void RestrictFunction(TFESpace3D *CoarseSpace, 
    TFESpace3D *FineSpace,
    double *CoarseFunction, double *FineFunction,
    double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
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
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
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
  
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,cell,k)
{
#pragma omp for schedule(static) nowait 
#endif
  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }
#ifdef _HYBRID
#pragma omp for schedule(static) nowait 
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }
#ifdef _HYBRID
#pragma omp for schedule(static) nowait 
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }
#ifdef _HYBRID
}
#endif

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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, MaxN_BaseFunctions3D*SizeOfDouble);
#ifdef _HYBRID
// #pragma omp parallel default(shared) private(j,k,s,l,cell,FineNumber,FineId,FineElement,FineBF,N_Fine,QQ,FineDOF,Val2,Index)
// #pragma omp for schedule(guided) nowait 
#endif
      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetRestrictionMatrix3D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions3D+l] * Val[l];
          } // endfor l
          Val2[k] += s;
        } // endfor k
      } // endfor j

      TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                      ->ChangeBF(CoarseColl, parent, Val2);
		      
#ifdef _HYBRID
// #pragma omp parallel default(shared) private(k,l)
// #pragma omp for schedule(guided) nowait 
#endif
      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        aux[l] += 1;
// 	 #pragma omp critical
        CoarseFunction[l] += Val2[k];
      } // endfor k
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
      QQ = TFEDatabase3D::GetRestrictionMatrix3D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

#ifdef _HYBRID
// #pragma omp parallel default(shared) private(k,s,l,Val2,Index)
{
// #pragma omp for schedule(guided) nowait 
#endif
      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase3D::GetBaseFunct3D(FineBF)
                      ->ChangeBF(FineColl, cell, Val);

#ifdef _HYBRID
// #pragma omp for schedule(guided) nowait 
#endif
      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

#ifdef _HYBRID
// #pragma omp for schedule(guided) nowait 
#endif
      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
#ifdef _HYBRID
	#pragma omp critical
#endif
	{
         CoarseFunction[l] += Val2[k];
         aux[l] += 1;
	}
      } // endfor k
#ifdef _HYBRID
}
#endif
    } // endelse
  } // endfor i

#ifdef _HYBRID
#pragma omp parallel default(shared) private(i)
#pragma omp for schedule(static) nowait 
#endif
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
void RestrictFunction(TFESpace3D *CoarseSpace, TFESpace3D *FineSpace,
    int N_Functions,
    double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
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
  double Val[MaxN_BaseFunctions3D];
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
  Val2 = new double[MaxN_BaseFunctions3D*N_Functions];

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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, N_Functions*MaxN_BaseFunctions3D*SizeOfDouble);

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetRestrictionMatrix3D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase3D::GetBaseFunct3D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

          Offset = IFunct*N_Coarse;
          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[k*MaxN_BaseFunctions3D+l] * Val[l];
            } // endfor l
            Val2[Offset + k] += s;
          } // endfor k
        } // endfor IFunct
      } // endfor j

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        CoarseOffset = IFunct*N_CoarseDOFs;
        Offset = IFunct*N_Coarse;

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
      QQ = TFEDatabase3D::GetRestrictionMatrix3D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
void IntoL20Vector3D(double *v, int Length, int order)
{
  double s;
  int i;

#ifdef _MPI
  int totallength;
  double temp;
#endif

  switch(order)
  {
    case -11:
      s=0;
      for(i=0;i<Length;i+=4)
        s += v[i];

#ifdef _MPI
      MPI_Allreduce(&Length, &totallength, 1, MPI_INT, MPI_SUM, TDatabase::ParamDB->Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);

      s = temp/(double)totallength;
#else
      s /= Length;
#endif

      s *= 4;

      for(i=0;i<Length;i+=4)
        v[i] -= s;
      break;

    case -12:
      s=0;
      for(i=0;i<Length;i+=10)
        s += v[i];

#ifdef _MPI
      MPI_Allreduce(&Length, &totallength, 1, MPI_INT, MPI_SUM, TDatabase::ParamDB->Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
      s = temp/(double)totallength;
#else
      s /= Length;
#endif

      s *= 10;

      for(i=0;i<Length;i+=10)
        v[i] -= s;

      break;

    default :
      s=0;
      for(i=0;i<Length;i++)
        s += v[i];

#ifdef _MPI
      MPI_Allreduce(&Length, &totallength, 1, MPI_INT, MPI_SUM, TDatabase::ParamDB->Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
      s = temp/(double)totallength;
#else
      s /= Length;
#endif

      for(i=0;i<Length;i++)
        v[i] -= s;
      break;

  }
} // IntoL20


/** project vector v into L20 */
void IntoL20FEFunction3D(double *v, int Length, TFESpace3D *FESpace)
{
  double s;
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs3D], *N_BaseFunct;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D], der[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last, number;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1;
  double *interpol;
  TNodalFunctional3D *nf;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D], temp;

// #ifdef _MPI
//   int rank, ID;
//   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
// #endif

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  for(i=0;i<MaxN_QuadPoints_3D;i++)
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
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 1;
    CurrentElement = FESpace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    // calculate all needed derivatives of this FE function
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = v[DOF[l]];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
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

    L1Int3D(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    error0 += LocError[0];
    error1 += LocError[1];

  } // endfor i

#ifdef _MPI
  temp = error0;
  MPI_Allreduce(&temp, &error0, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  temp = error1;
  MPI_Allreduce(&temp, &error1, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
#endif

  s = error0/error1;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    CurrentElement = FESpace->GetFE3D(i, cell);
    N_ = N_BaseFunct[CurrentElement];
    nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(CurrentElement);
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    for(j=0;j<N_Points;j++)
      PointValues[j] = s;
    nf->GetAllFunctionals(Coll, cell, PointValues, FunctionalValues);
    DOF = GlobalNumbers+BeginIndex[i];
    for(j=0;j<N_;j++)
      interpol[DOF[j]] = FunctionalValues[j];
  } // endfor i

  for(i=0;i<Length;i++)
    v[i] -= interpol[i];

  delete interpol;
} // IntoL20Function

/** project function v into L20 */
void IntoL20FEFunction3D(double *v, int Length, TFESpace3D *FESpace,
                       int velocity_space, int pressure_space)
{
  double s;
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs3D], *N_BaseFunct;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D], der[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last, number;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1;
  double *interpol;
  TNodalFunctional3D *nf;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];

#ifdef _MPI
//   MPI_Comm Comm;
//   int rank, ID;
  double temp;
// 
//   Comm = TDatabase::ParamDB->Comm;
//   MPI_Comm_rank(Comm, &rank);
#endif

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  for(i=0;i<MaxN_QuadPoints_3D;i++)
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

// #ifdef _MPI
//     ID = cell->GetSubDomainNo();
//     if(ID!=rank) 
//       continue; // this cell not belongs to this processor
// #endif

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 1;
    CurrentElement = FESpace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    // calculate all needed derivatives of this FE function
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = v[DOF[l]];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
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

    L1Int3D(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    error0 += LocError[0];
    error1 += LocError[1];

  } // endfor i

#ifdef _MPI
  temp = error0;
  MPI_Allreduce(&temp, &error0, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  temp = error1;
  MPI_Allreduce(&temp, &error1, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
#endif


  s = error0/error1;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    CurrentElement = FESpace->GetFE3D(i, cell);
    N_ = N_BaseFunct[CurrentElement];
    nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(CurrentElement);
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
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
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
      t += value * u2[index];
      u += value * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
        v3[index] += value3 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], x, y);
  return;
}
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  N_Active = A->GetActiveBound();

  //OutPut("r0 " << Ddot(N_UDOF+N_PDOF,r,r) << endl);
  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      //cout << "a " << value << endl;
      s -= value * u1[index];
      t -= value * u2[index];
      u -= value * u3[index];
      //OutPut(k<< " " << value << " " << endl);
    }
    //OutPut(s << " " << t << " " << u << endl);
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i
  //OutPut("r1 " << Ddot(N_UDOF+N_PDOF,r,r) << endl);

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
      value3 = B3Entries[j];
      //cout << "b " << i << " " << index << " " << value1 << " " << value2 << " " << value3 << endl;
      s -= value1 * u1[index] + value2 * u2[index] + 
        value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
        r3[index] -= value3 * t;
      }
    } // endfor j
    r4[i] = s;
  } // endfor i
  // OutPut("r2 " << Ddot(N_UDOF+N_PDOF,r,r) << endl);
}
void Defect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

// #ifdef _MPI
//  N_UDOF = 0;
//  N_PDOF = 0;
// 
//  if(TDatabase::ParamDB->ActiveProcess)
// #endif
  {
   CoupledDefect(A[0], B[0], B[1], B[2], x, b, r);
   N_UDOF = A[0]->GetN_Rows();
   N_PDOF = B[0]->GetN_Rows();
  }

  cout<<"TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE :: "<<TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE<<endl;
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

 
/** Navier--Stokes type 2 (NSTYPE==2) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T, 
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
      t += value * u2[index];
      u += value * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    q[i] = s;
  } // endfor i
 
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    v1[i] += s;
    v2[i] += t;
    v3[i] += u;
  } // endfor i

  return;
}
void MatVect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], B[3], B[4], B[5], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  N_Active = A->GetActiveBound();

  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
      t -= value * u2[index];
      u -= value * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    r4[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
    r3[i] -= u;
  } // endfor i
}
void Defect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], B[2], B[3], B[4], B[5], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  
//   cout<<"TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE :: "<<TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE<<endl;
  
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);

  return;
}

/** Navier--Stokes type 3 (NSTYPE==3) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  N_Active = A11->GetActiveBound();
  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      s += value * u1[index] + value1 * u2[index] + value2 * u3[index];
      t += value3* u1[index] + value4 * u2[index] + value5 * u3[index];
      u += value6* u1[index] + value7 * u2[index] + value8 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
  } // endfor i
  // Dirichlet and hanging nodes
  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
      u += value2 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        +  value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
        v3[index] += value3 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
                 B[0], B[1], B[2], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3;
  double value4, value5, value6, value7, value8;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  N_Active = A11->GetActiveBound();

  j = ARowPtr[0];
 
  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      s -= value * u1[index] + value1 * u2[index] + value2 * u3[index];
      t -= value3* u1[index] + value4 * u2[index] + value5 * u3[index];
      u -= value6* u1[index] + value7 * u2[index] + value8 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
      u -= value2 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
        r3[index] -= value3 * t;
      }
    } // endfor j
    r4[i] = s;
  } // endfor i
}
void Defect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

// #ifdef _MPI
//  N_UDOF = 0;
//  N_PDOF = 0;
// 
//  if(TDatabase::ParamDB->ActiveProcess)
// #endif
  {
   CoupledDefect(A[0], A[1], A[2], A[3], A[4], A[5], A[6],A[7], A[8],
                B[0], B[1], B[2], x, b, r);
   N_UDOF = A[0]->GetN_Rows();
   N_PDOF = B[0]->GetN_Rows();
  }

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** Navier--Stokes type 4 (NSTYPE==4) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s += value * value9 + value1 * value10 + value2 * value11;
      t += value3* value9 + value4 * value10 + value5 * value11;
      u += value6* value9 + value7 * value10 + value8 * value11;
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
  } // endfor i
  // Dirichlet and hanging nodes

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
      u += value2 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    q[i] = s;
  } // endfor i
 
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    v1[i] += s;
    v2[i] += t;
    v3[i] += u;
  } // endfor i

  return;
}
void MatVect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
                 B[0], B[1], B[2], B[3], B[4], B[5], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();
 
  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  j = ARowPtr[0];

  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s -= value * value9 + value1 * value10+ value2 * value11;
      t -= value3* value9 + value4 * value10+ value5 * value11;
      u -= value6* value9 + value7 * value10 + value8 * value11;
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
      u -= value2 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    r4[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
    r3[i] -= u;
  } // endfor i
}
void Defect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

// #ifdef _MPI
//  N_UDOF = 0;
//  N_PDOF = 0;
// 
//  if(TDatabase::ParamDB->ActiveProcess)
// #endif
  {
   CoupledDefect(A[0], A[1], A[2], A[3], A[4], A[5], A[6],A[7], A[8],
                B[0], B[1], B[2], B[3], B[4], B[5], x, b, r);
   N_UDOF = A[0]->GetN_Rows();
   N_PDOF = B[0]->GetN_Rows();
  }


  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** Navier--Stokes type 5 (NSTYPE==5) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, double *x, double *y)
{
  int i,j,k,l,m, jj;
  TSquareMatrixNSE3D *NSE_matrix;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *RowPtr, *KCol;
  double *Entries;
  int N_UDOF, N_PDOF, N_Active;
  TFESpace3D *USpace;
  TCollection *coll;
  double *u1, *u2, *u3, *p, *v1, *v2, *v3, *q, *r1, *r2, *r3, *r4;
  double s, t, r, value1, value2, value3;
  double value11, value12, value13, value21, value22, value23;
  double value31, value32, value33;
  int index;
  int N_Cells, N_Joints;
  TBaseCell *cell;
  int *GlobalNumbers, *BeginIndex, *DOF, *LocJb;
  int **JointDOFs, *EdgeDOF;
  double *LocAlpha;
  FE3D UElement;
  TFE3D *ele;
  TFEDesc3D *FEDesc_Obj;
  int N_U;
  
  NSE_matrix = (TSquareMatrixNSE3D*)A;

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
  u3 = u2 + N_UDOF;
  p  = u3 + N_UDOF;

  r1 = y;
  r2 = r1 + N_UDOF;
  r3 = r2 + N_UDOF;
  r4 = r3 + N_UDOF;

  memset(r4, 0, N_PDOF*SizeOfDouble);

  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = RowPtr[i+1];
    if(k == j)
    {
      // bubble row => no defect
      r1[i] = 0;
      r2[i] = 0;
      r3[i] = 0;
    }
    else
    {
      // skeleton dof
      s = 0;
      t = 0;
      r = 0;
      for(;j<k;j++)
      {
	  jj = 9*j;  
        index = KCol[j];
        value11 = Entries[jj + 0];
        value12 = Entries[jj + 1];
        value13 = Entries[jj + 2];
        value21 = Entries[jj + 3];
        value22 = Entries[jj + 4];
        value23 = Entries[jj + 5];
        value31 = Entries[jj + 6];
        value32 = Entries[jj + 7];
        value33 = Entries[jj + 8];

        value1 = u1[index];
        value2 = u2[index];
        value3 = u3[index];

        s += value11*value1 + value12*value2 + value13*value3;
        t += value21*value1 + value22*value2 + value23*value3;
        r += value31*value1 + value32*value2 + value33*value3;
      } // endfor j
      r1[i] = s;
      r2[i] = t;
      r3[i] = r;
    } // endelse
  } // endfor i

  // Dirichlet rows
  for(i=N_Active;i<N_UDOF;i++)
  {
    r1[i] = u1[i];
    r2[i] = u2[i];
    r3[i] = u3[i];
  }

  // use data from B-blocks
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(UElement);
    FEDesc_Obj = ele->GetFEDesc3D();
    JointDOFs = FEDesc_Obj->GetJointDOF();
    N_U = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    //  B entries
    LocJb = jb+BeginJb[i];

    // NOTE: Alpha stores the negative B block entry
    for(j=0;j<N_Joints;j++)
    {
      // get local data
      LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
      EdgeDOF = JointDOFs[j]; 
        
      l = LocJb[j];
      if(l<N_U)
      {
        // jb in first component
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        m = DOF[l];
        r4[i] -= LocAlpha[k]*u1[m];
        if(m<N_Active)
          r1[m] -= LocAlpha[k]*p[i];
      }
      else
      {
        if(l<2*N_U)
        {
          // jb in second component
          l -= N_U;
          for(k=0;k<N_DOFperJoint;k++)
            if(EdgeDOF[k] == l)
              break;

          m = DOF[l];
          r4[i] -= LocAlpha[k+N_DOFperJoint]*u2[m];
          if(m<N_Active)
            r2[m] -= LocAlpha[k+N_DOFperJoint]*p[i];
        }
        else
        {
          // jb in third component
          l -= 2*N_U;
          for(k=0;k<N_DOFperJoint;k++)
            if(EdgeDOF[k] == l)
              break;

          m = DOF[l];
          r4[i] -= LocAlpha[k+2*N_DOFperJoint]*u3[m];
          if(m<N_Active)
            r3[m] -= LocAlpha[k+2*N_DOFperJoint]*p[i];
        }
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
  TSquareMatrixNSE3D *NSE_matrix;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *RowPtr, *KCol;
  double *Entries;
  int N_UDOF, N_PDOF, N_Active;
  TFESpace3D *USpace;
  TCollection *coll;
  double *u1, *u2, *u3, *p, *v1, *v2, *v3, *q, *r1, *r2, *r3, *r4;
  double s, t, rr, value1, value2, value3;
  double value11, value12, value13, value21, value22, value23;
  double value31, value32, value33;
  int index;
  int N_Cells, N_Joints;
  TBaseCell *cell;
  int *GlobalNumbers, *BeginIndex, *DOF, *LocJb;
  int **JointDOFs, *EdgeDOF;
  double *LocAlpha;
  FE3D UElement;
  TFE3D *ele;
  TFEDesc3D *FEDesc_Obj;
  int N_U;
  
  NSE_matrix = (TSquareMatrixNSE3D*)A;

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
  u3 = u2 + N_UDOF;
  p  = u3 + N_UDOF;

  v1 = b;
  v2 = v1 + N_UDOF;
  v3 = v2 + N_UDOF;
  q  = v3 + N_UDOF;

  r1 = r;
  r2 = r1 + N_UDOF;
  r3 = r2 + N_UDOF;
  r4 = r3 + N_UDOF;

  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = RowPtr[i+1];
    if(k == j)
    {
      // bubble row => no defect
      r1[i] = 0;
      r2[i] = 0;
      r3[i] = 0;
    }
    else
    {
      // skeleton dof
      s = v1[i];
      t = v2[i];
      rr = v3[i];
      for(;j<k;j++)
      {
        index = KCol[j];
        value11 = Entries[9*j + 0];
        value12 = Entries[9*j + 1];
        value13 = Entries[9*j + 2];
        value21 = Entries[9*j + 3];
        value22 = Entries[9*j + 4];
        value23 = Entries[9*j + 5];
        value31 = Entries[9*j + 6];
        value32 = Entries[9*j + 7];
        value33 = Entries[9*j + 8];

        value1 = u1[index];
        value2 = u2[index];
        value3 = u3[index];

        s -= value11*value1 + value12*value2 + value13*value3;
        t -= value21*value1 + value22*value2 + value23*value3;
        rr -= value31*value1 + value32*value2 + value33*value3;
      } // endfor j
      r1[i] = s;
      r2[i] = t;
      r3[i] = rr;
    } // endelse
  } // endfor i

  // Dirichlet rows
  for(i=N_Active;i<N_UDOF;i++)
  {
    r1[i] = v1[i] - u1[i];
    r2[i] = v2[i] - u2[i];
    r3[i] = v3[i] - u3[i];
  }

  memcpy(r4, q, N_PDOF*SizeOfDouble);

  // use data from B-blocks
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(UElement);
    FEDesc_Obj = ele->GetFEDesc3D();
    JointDOFs = FEDesc_Obj->GetJointDOF();
    N_U = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    //  B entries
    LocJb = jb+BeginJb[i];

    // NOTE: Alpha stores the negative B block entry
    for(j=0;j<N_Joints;j++)
    {
      // get local data
      LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
      EdgeDOF = JointDOFs[j]; 
        
      l = LocJb[j];
      if(l<N_U)
      {
        // jb in first component
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        m = DOF[l];
        r4[i] += LocAlpha[k]*u1[m];
        if(m<N_Active)
          r1[m] += LocAlpha[k]*p[i];

        // cout << " " << DOF[l] << " " << i;
        // cout << " " << LocAlpha[k] << endl;
      }
      else
      {
        if(l<2*N_U)
        {
          // jb in second component
          l -= N_U;
          for(k=0;k<N_DOFperJoint;k++)
            if(EdgeDOF[k] == l)
              break;

          m = DOF[l];
          r4[i] += LocAlpha[k+N_DOFperJoint]*u2[m];
          if(m<N_Active)
            r2[m] += LocAlpha[k+N_DOFperJoint]*p[i];

          // cout << " " << DOF[l]+N_UDOF << "(" << DOF[l] << ")";
          // cout << " " << i;
          // cout << " " << LocAlpha[k+N_DOFperJoint] << endl;
        }
        else
        {
          // jb in third component
          l -= 2*N_U;
          for(k=0;k<N_DOFperJoint;k++)
            if(EdgeDOF[k] == l)
              break;

          m = DOF[l];
          r4[i] += LocAlpha[k+2*N_DOFperJoint]*u3[m];
          if(m<N_Active)
            r3[m] += LocAlpha[k+2*N_DOFperJoint]*p[i];

          // cout << " " << DOF[l]+2*N_UDOF << "(" << DOF[l] << ")";
          // cout << " " << i;
          // cout << " " << LocAlpha[k+2*N_DOFperJoint] << endl;
        }
      }
    } // endfor j
  } // endfor i
} // end CoupledDefect

void Defect_NSE5(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

// #ifdef _MPI
//  N_UDOF = 0;
//  N_PDOF = 0;
// 
//  if(TDatabase::ParamDB->ActiveProcess)
// #endif
  {
   CoupledDefect(A[0], x, b, r);
   N_UDOF = A[0]->GetN_Rows();
   N_PDOF = ((TSquareMatrixNSE3D *)A[0])->GetFESpace()->GetN_Cells();
  }

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,0);
  return;
}

/** Navier--Stokes type 4 (NSTYPE==4) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
		    TSquareMatrix *C,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *CEntries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s += value * value9 + value1 * value10 + value2 * value11;
      t += value3* value9 + value4 * value10 + value5 * value11;
      u += value6* value9 + value7 * value10 + value8 * value11;
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
  } // endfor i
  // Dirichlet and hanging nodes

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
      u += value2 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    q[i] = s;
  } // endfor i
 
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    v1[i] += s;
    v2[i] += t;
    v3[i] += u;
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
void MatVect_NSE14(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
                 A[9], B[0], B[1], B[2], B[3], B[4], B[5], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
		   TSquareMatrix *C,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *CEntries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();
 
  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  j = ARowPtr[0];

  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s -= value * value9 + value1 * value10+ value2 * value11;
      t -= value3* value9 + value4 * value10+ value5 * value11;
      u -= value6* value9 + value7 * value10 + value8 * value11;
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
      u -= value2 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    r4[i] = s;
    
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
    r3[i] -= u;
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
	  s += value1 * p[index];
      } // endfor j
      q[i] -= s;
  } // endfor i
}
void Defect_NSE14(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], A[4], A[5], A[6],A[7], A[8],
                A[9], B[0], B[1], B[2], B[3], B[4], B[5], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}




#endif 
