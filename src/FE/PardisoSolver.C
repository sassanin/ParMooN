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
   
#include <PardisoSolver.h>
#include <string.h>
#include <stdlib.h>

#include <stdio.h>
#include <sstream>

static const char *stringtable_iparam[] = 
{"Use default values",
 "Fill-In reduction reorderings",
 "Number of processors",
 "Preconditioned CGS",
 "User permutation",
 "Write solution on X",
 "Number of performed iterative refinement steps",
 "Iterative refinement steps",
 "Unused",
 "Piviot perturbation",
 "Scaling vector",
 "Solving with transposed Matrix",
 "Improved accuracy using (non-)symmetric weighted matchings",
 "Number of perturbed pivots",
 "Peak memory symbolic factorization",
 "Permanent memory symbolic factorization",
 "Memory numerical factorization and solution",
 "Number nonzeros in factors",
 "MFlops of factorization",
 "CG/CGS diagnostics",
 "Pivoting for symmetric indefinit matrices",
 "Inertia: Number of positive eigenvalues",
 "Inertia: Number of negative eigenvalues",
 "Parallel numerical Factorization",
 "Parallel Forward/Backward Solve",
 "Splitting Forward/Backward Solve",
 "Unused",
 "Parallel Reordering for METIS",
 "Switch between 32-bit and 64-bit factorization",
 "Control the size of supernodes",
 "Partial solve for sparse right-hand sides and sparse solution",
 "Use the multi-recursive iterative linear solver",
 "Determinant of real symmetric indefinite matrix",
 "Identical solution independent on the number of processors",
 "Unused",};
			

extern "C"
{
   // PARDISO 4.0.0

  int pardisoinit_(void**, int*, int*, int*, double*, int*);

  int pardiso_(void **, int *, int *, int *, int *, int *,
               double *, int *, int *, int *, int *, int *,
               int *, double *, double *, int *, double*);
} 


// Constructor

TPardisoSolver::TPardisoSolver() : TDirectSolver()
{
  int solver = 0; // solver method (0 = sparse direct solver)
 
  int err=0;
  char *env;

  KCol = NULL;
  RowPtr = NULL;
  Entries = NULL;
  Entries_size = 0;
  KCol_size = 0;
  RowPtr_size = 0;
  
  N_Entries = 0;
  N_Eq = 0;
  mtype = 11; // type of matrix (11 = real and nonsymmetric)

  iparam = new int [64];
  dparam = new double [64];
  pt = new void* [64];

  memset(iparam, 0, sizeof(int)*64);
  memset(dparam, 0, sizeof(double)*64);  

  env = getenv("OMP_NUM_THREADS");
  if ( env ) 
  {
    sscanf(env, "%d", &num_prc);
    OutPut("OMP_NUM_THREADS = "<<num_prc<<endl);
  }
  else 
  {
     setenv("OMP_NUM_THREADS", "1", 1);
     OutPut("Set environment OMP_NUM_THREADS to 1"<<endl);
     exit(0);
  }
    
  iparam[2] = num_prc; //there is no default value for number of processors  

  pardisoinit_(pt, &mtype, &solver, iparam, dparam, &err);

  ErrorMsg(err); 
//   cout << "init complete!" << endl;

  iparam[1] = 2; // fill-in reduction reordering (2 - metis)
//   iparam[3] = 82;
  iparam[7] = 100; // iterative refinement steps
  iparam[9] = 6;  //  Pivot perturbation
  iparam[10] = 1; // scaling
//   iparam[11] = 1; // solve transpose system
  iparam[12] = 2; // weighted matching
//   iparam[26] = 1; // check array for ordering
  iparam[23] = 0; // parallel factorization
  iparam[31] = 0;

//   PrintIparam();

//   std::ostringstream os;
// 
//   os << "bench_" << num_prc << ends;
// 
//   dat = new std::ofstream(os.str().c_str());
//   time_a = time_f = time_s = 0.0;
//   runs_a = runs_f = runs_s = 0;
//   exit(0);

  pp = false;
}

// Destructor

TPardisoSolver::~TPardisoSolver()
{
  int maxfct = 1;
  int mnum = 1;
  int phase = -1; // release all internal memory
  int msglvl = SetMsgLvl();
  int nrhs = 0;
  int err;

  pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eq, Entries, RowPtr,
	   KCol, NULL, &nrhs, iparam, &msglvl, NULL, NULL, &err, dparam);
    
  ErrorMsg(err);

  FreeMemory();

  delete [] iparam;
  delete [] dparam;
  delete [] pt;

//   dat->close();

//   delete dat;
}

/** Methods */

#ifdef __2D__
void TPardisoSolver::SetMatrix(TSquareMatrix2D *Matrix)
{
  double *entries = Matrix->GetEntries();
  int *kcol = Matrix->GetKCol();
  int *rowptr = Matrix->GetRowPtr();

  N_Eq = Matrix->GetN_Rows();
  N_Entries = Matrix->GetN_Entries();

  FreeMemory();

  Entries = new double [N_Entries];
  KCol = new int [N_Entries];
  RowPtr = new int [N_Eq+1];

  memcpy(Entries, entries, sizeof(double)*N_Entries);
  memcpy(KCol, kcol, sizeof(int)*N_Entries);
  memcpy(RowPtr, rowptr, sizeof(int)*(N_Eq+1));

  if( Matrix->GetColOrder() != 1)
    Sort();

  ShiftIndicies();  
}

// nstype 1
void TPardisoSolver::SetMatrix(TSquareMatrix2D *A, TMatrix *B1, TMatrix *B2)
{
  int i, j, k, l, begin, end, pos, len;
  int N_U, N_P, N_, N_Active;

  int *KColA = A->GetKCol();
  int *RowPtrA = A->GetRowPtr();
  int *KColB = B1->GetKCol();
  int *RowPtrB = B1->GetRowPtr();

  double *EntriesA = A->GetEntries();
  double *EntriesB1 = B1->GetEntries();
  double *EntriesB2 = B2->GetEntries();

  double t1, t2;

  FreeMemory();

  N_U = A->GetN_Rows();
  N_P = B1->GetN_Rows();
  N_Active = A->GetActiveBound();

  N_Eq = 2*N_U + N_P;
  N_Entries = 2*RowPtrA[N_U] + 4*RowPtrB[N_P];

  AllocMemory();

  RowPtr[0] = 0;
  pos = 0;
  t1 = GetTime();

  // fill combined matrix
  for(i=0;i<N_U;i++)
  {
    // first velocity component
    begin = RowPtrA[i];
    end = RowPtrA[i+1];

//     pos  = begin;
//     if(i<N_Active) pos += RowPtr

    for(j=begin;j<end;j++)
    {
      // A11
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
      pos++;
    }
    // B1T
    if(i<N_Active)
    {
      // this is quit inefficient, think about more efficient solutions
      // later
      // loop over column indices of matrix B1
      for (k=0;k< N_P; k++)
      {
        begin = RowPtrB[k];
        end = RowPtrB[k+1];

        // if column index equal to i
        for(l=begin;l<end;l++)
        {
          if (KColB[l] == i)
          {
            Entries[pos] = EntriesB1[l];
            KCol[pos] = k+2*N_U;
            pos++;
          }
        }
      }
    }
    RowPtr[i+1] = pos;
  }

  // second velocity component
  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A22
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }
    // B2T
    if(i<N_Active)
    {
      // this is quit inefficient, think about more efficient solutions
      // later
      // loop over column indices of matrix B1
      for (k=0;k< N_P; k++)
      {
        begin = RowPtrB[k];
        end = RowPtrB[k+1];

        // if column index equal to i
        for(l=begin;l<end;l++)
        {
          if (KColB[l] == i)
          {
            Entries[pos] = EntriesB2[l];
            KCol[pos] = k+2*N_U;
            pos++;
          }
        }
      }
    }
    RowPtr[N_U+i+1] = pos;
  }

  // pressure
  for(i=0;i<N_P;i++)
  {
    // B1
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB1[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    // B2
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB2[j];
      KCol[pos] = KColB[j]+N_U;
      pos++;
    }
    RowPtr[2*N_U+i+1] = pos;
  }

  t2 = GetTime();
  cout << "Time for building system: " << (t2-t1) << " sec" << endl;

  Sort();
  
  ShiftIndicies(); 

//   exit(0);
}

void TPardisoSolver::SetMatrixPar(TSquareMatrix2D *A, TMatrix *B1, TMatrix *B2)
{
  int i, j, k, l, begin, end, pos1, pos2;
  int N_U, N_P, N_, N_Active;

  int *ColPtrB, *KRowB, *MapBT;
  int *KColA = A->GetKCol();
  int *RowPtrA = A->GetRowPtr();
  int *KColB = B1->GetKCol();
  int *RowPtrB = B1->GetRowPtr();

  double *EntriesA = A->GetEntries();
  double *EntriesB1 = B1->GetEntries();
  double *EntriesB2 = B2->GetEntries();

// cout << "enter SetMatrixPar()" << endl;

  double t1, t2;

  FreeMemory();

  N_U = A->GetN_Rows();
  N_P = B1->GetN_Rows();
  N_Active = A->GetActiveBound();
// cout << "N_Active = " << N_Active << endl;

  N_Eq = 2*N_U + N_P;
  N_Entries = 2*RowPtrA[N_U] + 4*RowPtrB[N_P];

  AllocMemory();

  t1 = GetTime();

  GetTransposedArrays(B1, KRowB, ColPtrB, MapBT);
  FillRowPtr(N_U, N_P, N_Active, RowPtrA, RowPtrB, ColPtrB);

  // fill combined matrix
//   #pragma omp parallel for private(i,j,k,l,pos1,pos2,begin,end) \
// 			   firstprivate(RowPtrA, RowPtrB) \
// 			   firstprivate(EntriesA, EntriesB1, KColA, KColB) \
// 			   firstprivate(N_U, N_P, N_Active, EntriesB2)
  for(i=0;i<N_U;++i)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    pos1 = RowPtr[i];
    pos2 = RowPtr[i+N_U];

    for(j=begin;j<end;++j)
    {
      // A11
      Entries[pos1] = EntriesA[j];
      KCol[pos1] = KColA[j];
      ++pos1;

      //A22
      Entries[pos2] = EntriesA[j];
      KCol[pos2] = KColA[j]+N_U;
      ++pos2;
    }
    // B1T
    if(i<N_Active)
    {
      begin = ColPtrB[i];
      end = ColPtrB[i+1];

      for(l=begin;l<end;++l)
      {
	Entries[pos1] = EntriesB1[MapBT[l]];
	KCol[pos1] = KRowB[MapBT[l]] + 2*N_U;
	++pos1;

	Entries[pos2] = EntriesB2[MapBT[l]];
	KCol[pos2] = KRowB[MapBT[l]] + 2*N_U;
	++pos2;
      }
    }
//     if( RowPtr[i+1] != pos1 )
//     {
//       cout << i << ": " << RowPtr[i+1] << "/" << pos1 << endl;exit(0);
//     }
//     cout << pos1 << endl;
  }

//   #pragma omp parallel for private(i,j,pos1,begin,end)
  for(i=0;i<N_P;i++)
  {
    // B1
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    pos1 = RowPtr[i+2*N_U];

    for(j=begin;j<end;j++)
    {
      Entries[pos1] = EntriesB1[j];
      KCol[pos1] = KColB[j];
      pos1++;
    }
    // B2
    for(j=begin;j<end;j++)
    {
      Entries[pos1] = EntriesB2[j];
      KCol[pos1] = KColB[j]+N_U;
      pos1++;
    }
//     RowPtr[2*N_U+i+1] = pos;
  }

  t2 = GetTime();
  cout << "Time for building system: " << (t2-t1) << " sec" << endl;

// cout << N_Entries << "/" << RowPtr[N_Eq] << endl;

  Sort();
  
  ShiftIndicies(); 

  delete [] ColPtrB;
  delete [] KRowB;
  delete [] MapBT;

// cout << "leave SetMatrixPar()" << endl;
//   exit(0);
}

// nstype 2
void TPardisoSolver::SetMatrix(TSquareMatrix2D *sqmatrixA,
				 TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
				 TMatrix2D *matrixB1,  TMatrix2D *matrixB2)
{
  int i, j, pos, begin, end, N_U, N_P, N_Active;
  int *KColA, *KColB, *KColBT, *RowPtrA, *RowPtrB, *RowPtrBT;
  double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;

  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_Eq = 2*N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();

  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  EntriesA = sqmatrixA->GetEntries();
  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();

  N_Entries = 2*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U];

  FreeMemory();
  AllocMemory();

  RowPtr[0] = 0;
  pos = 0;
  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
      pos++;
    }

    if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1T[j];
        KCol[pos] = KColBT[j]+2*N_U;
        pos++;
      }
    }
    RowPtr[i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2T[j];
        KCol[pos] = KColBT[j]+2*N_U;
        pos++;
      }
    }
    RowPtr[N_U+i+1] = pos;
  }

  for(i=0;i<N_P;i++)
  {
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB1[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB2[j];
      KCol[pos] = KColB[j]+N_U;
      pos++;
    }
    RowPtr[2*N_U+i+1] = pos;
  }

  Sort();

  ShiftIndicies();
}

// nstype 4
void TPardisoSolver::SetMatrix(TSquareMatrix2D *sqmatrixA11,
      TSquareMatrix2D *sqmatrixA12, TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22, TMatrix2D *matrixB1T,
      TMatrix2D *matrixB2T, TMatrix2D *matrixB1,  TMatrix2D *matrixB2)
{
  int i, j, pos, begin, end, N_U, N_P, N_Active;
  int *KColA, *KColB, *KColBT, *RowPtrA, *RowPtrB, *RowPtrBT;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  
cout << "NSTYPE 4" << endl;


  N_U = sqmatrixA11->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_Eq = 2*N_U + N_P;
  N_Active = sqmatrixA11->GetActiveBound();

  KColA = sqmatrixA11->GetKCol();
  RowPtrA = sqmatrixA11->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  EntriesA11 = sqmatrixA11->GetEntries();
  EntriesA12 = sqmatrixA12->GetEntries();
  EntriesA21 = sqmatrixA21->GetEntries();
  EntriesA22 = sqmatrixA22->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();

  N_Entries = 4*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U];

  FreeMemory();
  AllocMemory();

  RowPtr[0] = 0;
  pos = 0;
  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesA11[j];
      KCol[pos] = KColA[j];
      pos++;

      Entries[pos] = (i<N_Active)?EntriesA12[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1T[j];
        KCol[pos] = KColBT[j]+2*N_U;
        pos++;
      }
    }
    RowPtr[i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active)?EntriesA21[j]:0;
      KCol[pos] = KColA[j];
      pos++;

      Entries[pos] = EntriesA22[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2T[j];
        KCol[pos] = KColBT[j]+2*N_U;
        pos++;
      }
    }
    RowPtr[N_U+i+1] = pos;
  }

  for(i=0;i<N_P;i++)
  {
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB1[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB2[j];
      KCol[pos] = KColB[j]+N_U;
      pos++;
    }
    RowPtr[2*N_U+i+1] = pos;
  }

  Sort();
  ShiftIndicies();
}
#endif

#ifdef __3D__
/** NSTYPE 2 */
void TPardisoSolver::SetMatrix(TSquareMatrix3D *sqmatrixA,
		      TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		      TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3)
{
  int N_U, N_P, N_Active, pos;
  int *KColA, *KColB, *KColBT;
  int *RowPtrA, *RowPtrB, *RowPtrBT;
  int i, j, begin, end;
  
  double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB3;
  double *EntriesB1T, *EntriesB2T, *EntriesB3T;
  
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_Eq = 3*N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();
  
  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  EntriesA = sqmatrixA->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB3 = matrixB3->GetEntries();

  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();
  EntriesB3T = matrixB3T->GetEntries();
  
  N_Entries = 3*RowPtrA[N_U] + 3*RowPtrB[N_P] + 3*RowPtrBT[N_U];
  
//   FreeMemory();
  AllocMemory();
  
  RowPtr[0] = 0;
  pos = 0;
  
  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
      pos++;
    }

    if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[N_U+i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB3T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[2*N_U+i+1] = pos;
  }

  for(i=0;i<N_P;i++)
  {
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB1[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB2[j];
      KCol[pos] = KColB[j]+N_U;
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB3[j];
      KCol[pos] = KColB[j]+2*N_U;
      pos++;
    }
    RowPtr[3*N_U+i+1] = pos;
  }

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[3*N_U];
    end = RowPtr[3*N_U+1];
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 3*N_U;
//     rhs[3*N_U] = 0;

    pp = true;
    rhs_index = 3*N_U;
  }

  Sort();

  ShiftIndicies();
}

void TPardisoSolver::SetMatrix(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
		      TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
		      TSquareMatrix3D *sqmatrixA22, TSquareMatrix3D *sqmatrixA23,
		      TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
		      TSquareMatrix3D *sqmatrixA33,
		      TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		      TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3)
{
  double *EntriesA11, *EntriesA12, *EntriesA13, *EntriesA21, *EntriesA22, *EntriesA23;
  double *EntriesA31, *EntriesA32, *EntriesA33, *EntriesB1, *EntriesB2, *EntriesB3;
  double *EntriesB1T, *EntriesB2T, *EntriesB3T;
  int *KColA, *RowPtrA, *KColB, *RowPtrB, *KColBT, *RowPtrBT;
  int N_U, N_P, N_, N_Active;
  int i, j, begin, end, pos;
  
  N_U = sqmatrixA11->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 3*N_U + N_P;
  N_Eq = N_;
  N_Active = sqmatrixA11->GetActiveBound();
  
  KColA = sqmatrixA11->GetKCol();
  RowPtrA = sqmatrixA11->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();
  
  EntriesA11 = sqmatrixA11->GetEntries();
  EntriesA12 = sqmatrixA12->GetEntries();
  EntriesA13 = sqmatrixA13->GetEntries();
  EntriesA21 = sqmatrixA21->GetEntries();
  EntriesA22 = sqmatrixA22->GetEntries();
  EntriesA23 = sqmatrixA23->GetEntries();
  EntriesA31 = sqmatrixA31->GetEntries();
  EntriesA32 = sqmatrixA32->GetEntries();
  EntriesA33 = sqmatrixA33->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB3 = matrixB3->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();
  EntriesB3T = matrixB3T->GetEntries();

  N_Entries = 9*RowPtrA[N_U] + 3*RowPtrB[N_P] + 3*RowPtrBT[N_U];
  
//   FreeMemory();
  AllocMemory();
  
  RowPtr[0] = 0;

  pos = 0;
  
  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A11
      Entries[pos] = EntriesA11[j];
      KCol[pos] = KColA[j];
      pos++;
      // A12
      Entries[pos] = (i<N_Active)?EntriesA12[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
      // A13
      Entries[pos] = (i<N_Active)?EntriesA13[j]:0;
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
      // B1T
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A21
      Entries[pos] = (i<N_Active)?EntriesA21[j]:0;
      KCol[pos] = KColA[j];
      pos++;
      // A22
      Entries[pos] = EntriesA22[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
      // A23
      Entries[pos] = (i<N_Active)?EntriesA23[j]:0;
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
      // B2T
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[N_U+i+1] = pos;
  }

  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A31
      Entries[pos] = (i<N_Active)?EntriesA31[j]:0;
      KCol[pos] = KColA[j];
      pos++;
      // A32
      Entries[pos] = (i<N_Active)?EntriesA32[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
      // A33
      Entries[pos] = EntriesA33[j];
      KCol[pos] = KColA[j]+2*N_U;
      pos++;
    }

    if(i<N_Active)
    {
      // B3T
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB3T[j];
        KCol[pos] = KColBT[j]+3*N_U;
        pos++;
      }
    }
    RowPtr[2*N_U+i+1] = pos;
  }

  for(i=0;i<N_P;i++)
  {
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      // B1
      Entries[pos] = EntriesB1[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      // B2
      Entries[pos] = EntriesB2[j];
      KCol[pos] = KColB[j]+N_U;
      pos++;
    }
    for(j=begin;j<end;j++)
    {
      // B3
      double _tmp = EntriesB3[j];
      Entries[pos] = _tmp;
      KCol[pos] = KColB[j]+2*N_U;
      pos++;
    }
    RowPtr[3*N_U+i+1] = pos;
  }

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[3*N_U];
    end = RowPtr[3*N_U+1];
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 3*N_U;
//     rhs[3*N_U] = 0;

    pp = true;
    rhs_index = 3*N_U;
  }
  
  Sort();
  ShiftIndicies();
}

#endif

void TPardisoSolver::Solve(double *sol, double *rhs)
{
  int maxfct = 1;
  int mnum = 1;
  int phase = 33; // solve
  int msglvl = SetMsgLvl();
  int nrhs = 1;
  int err;
  double t1, t2;

// cout << "enter TPardisoSolver::Solve()" << endl;
  t1 = GetTime();

  if ( pp ) 
    rhs[rhs_index] = 0;
  
  pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eq, Entries, RowPtr,
	   KCol, NULL, &nrhs, iparam, &msglvl, rhs, sol, &err, dparam);
    
  t2 = GetTime();

  OutPut("Time for solve: " << (t2-t1) << " sec" << endl);

  ErrorMsg(err);

  OutPut("Number of iterative refinement steps done: "<<iparam[6]<<endl);
//   PrintIparam();
// cout << "leave ParDirectSolver::Solve()" << endl;

  runs_s++;
  time_s += (t2-t1);

//   iparam[3] = 101;
//   return (t2-t1);
}

void TPardisoSolver::FactorizeSolve(double *sol, double *rhs)
{
  int maxfct = 1;
  int mnum = 1;
  int phase = 23; // factorize and solve
  int msglvl = SetMsgLvl();
  int nrhs = 1;
  int err;
  double t1, t2;

// cout << "enter TPardisoSolver::Solve()" << endl;
  t1 = GetTime();

  if ( pp ) 
  {
    cout << "pp" << endl;
    rhs[rhs_index] = 0;
  }
  
  pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eq, Entries, RowPtr,
	   KCol, NULL, &nrhs, iparam, &msglvl, rhs, sol, &err, dparam);
    
  t2 = GetTime();

  OutPut("Time for factoize and solve: " << (t2-t1) << " sec" << endl);

  ErrorMsg(err);

  OutPut("Number of iterative refinement steps done: "<<iparam[6]<<endl);
//   PrintIparam();
// cout << "leave ParDirectSolver::Solve()" << endl;

  runs_s++;
  time_s += (t2-t1);

//   iparam[3] = 82;
//   return (t2-t1);
}

void TPardisoSolver::Factorize()
{
  int maxfct = 1;
  int mnum = 1;
  int phase = 22; // factorize
  int msglvl = SetMsgLvl();
  int nrhs = 0;
  int err;
  double t1, t2;

// cout << "enter TPardisoSolver::Factorization()" << endl;
  t1 = GetTime();

  pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eq, Entries, RowPtr,
	   KCol, NULL, &nrhs, iparam, &msglvl, NULL, NULL, &err, dparam);
    
  t2 = GetTime();

  OutPut("Time for factorization: " << (t2-t1) << " sec" << endl);

  ErrorMsg(err);

//   OutPut("Number of iterative refinement steps done: "<<iparam[6]<<endl);
//   PrintIparam();
// cout << "leave TPardisoSolver::Factorization()" << endl;

  runs_f++;
  time_f += (t2-t1);

  
  iparam[3] = 101;
//   return (t2-t1);
}


void TPardisoSolver::Analyse()
{
  int maxfct = 1;
  int mnum = 1;
  int phase = 11; // analyse
  int msglvl = SetMsgLvl();
  int nrhs = 0;
  int err;
  double t1, t2;

  iparam[3] = 0;
  
// cout << "enter TPardisoSolver::Analysis()" << endl;
  t1 = GetTime();

  pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eq, Entries, RowPtr,
	   KCol, NULL, &nrhs, iparam, &msglvl, NULL, NULL, &err, dparam);
    
  t2 = GetTime();

  OutPut("Time for analyse: " << (t2-t1) << " sec" << endl);

  ErrorMsg(err);

  runs_a++;
  time_a += (t2-t1);

// cout << "leave TPardisoSolver::Analysis()" << endl;

  cout << "Pardiso: Peek memory: " << iparam[14] << " KByte (";
  cout << iparam[14]/1024 << " MByte)" << endl;

//   return (t2-t1);
}

void TPardisoSolver::ShiftIndicies()
{
  for(int i=0;i<N_Entries;i++)
    ++KCol[i];

  for(int i=0;i<N_Eq+1;i++)
    ++RowPtr[i];
}

void TPardisoSolver::Sort()
{
  int i, j, k, l, begin, end;
  double value;

    for(i=0;i<N_Eq;i++)
    {
      begin=RowPtr[i];
      end=RowPtr[i+1];

      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Entries[j];
            KCol[j] = KCol[k]; Entries[j] = Entries[k];
            KCol[k] = l;       Entries[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
}

void TPardisoSolver::PrintIparam()
{
  for(int i=0;i<34;i++) {
    OutPut("IPARAM["<<i+1<<"]: "<<iparam[i]);
    OutPut(" - "<<stringtable_iparam[i]<<endl);
  }

  OutPut("IPARAM[35] - IPARAM[64] - Unused"<<endl);
}

void TPardisoSolver::ErrorMsg(int err)
{
  switch(err)
  {
    case 0: break;
    case -4:
      OutPut("Zero pivot, numerical fact. or iterative reﬁnement problem."<<endl);
      break;
    case -10:
      OutPut("No license file \"pardiso.lic\" found"<<endl); break;
    case -11:
      OutPut("License is expired"<<endl);break;
    case -12:
      OutPut("Wrong username or hostname"<<endl);break;

    default:
      OutPut("Unknown error = "<< err << endl);
  }

  if(err != 0) exit(-1);
}

void TPardisoSolver::AllocMemory()
{
  if ( !Entries || N_Entries > Entries_size)
  {
    delete [] Entries;
    Entries = new double [N_Entries];
  }
  if ( !KCol || N_Entries > Entries_size )
  {
    delete [] KCol;
    KCol = new int [N_Entries];
  }
  if ( !RowPtr || (N_Eq+1) > RowPtr_size )
  {
    delete [] RowPtr;
    RowPtr = new int [N_Eq+1];
  }
}

void TPardisoSolver::FreeMemory()
{
  if ( Entries ) delete [] Entries;
  if ( KCol ) delete [] KCol;
  if ( RowPtr ) delete [] RowPtr;

  Entries = NULL;
  KCol = NULL;
  RowPtr = NULL;
}

#ifdef __2D__
// nstype 1
void TPardisoSolver::FillRowPtr(int N_U, int N_P, int N_Active,
				  int *RowPtrA, int *RowPtrB, int *ColPtrB)
{
  int N_EntriesA, N_EntriesB, N_;

// cout << "Enter FillRowPtr()" << endl;

  N_EntriesA = RowPtrA[N_U];
  N_EntriesB = RowPtrB[N_P];

  N_ = N_EntriesA + ColPtrB[N_Active];

  RowPtr[0] = 0;

  //  A 0 B1T
  //  0 A B2T
  for(int i=0;i<N_U;++i)
  {
    RowPtr[i+1] = RowPtrA[i+1];
    RowPtr[i+1+N_U] = RowPtrA[i+1] + N_;

    if(i<N_Active)
    {
      RowPtr[i+1] += ColPtrB[i+1];
      RowPtr[i+1+N_U] += ColPtrB[i+1];
    }
    else
    {
      RowPtr[i+1] += ColPtrB[N_Active];
      RowPtr[i+1+N_U] += ColPtrB[N_Active];
    }
  }

  // B1 B2 0
  for(int i=0;i<N_P;i++)
  {
    RowPtr[i+1+2*N_U] = 2 * ( RowPtrB[i+1] + N_);
  }
//  cout << "Leave FillRowPtr()" << endl;
// exit(0);
}

void TPardisoSolver::GetTransposedArrays(TMatrix *B, int *&KRowB,
					   int *&ColPtrB, int *&MapB)
{
  int N_C, N_R, N_EntriesB, begin, end, len, index;
  int *KColB , *RowPtrB;

// cout << "Enter GetTransposedArrays()" << endl;

  KColB = B->GetKCol();
  RowPtrB = B->GetRowPtr();

  N_C = B->GetN_Columns();
  N_R = B->GetN_Rows();
  N_EntriesB = B->GetN_Entries();

  KRowB = new int [N_EntriesB];
  MapB = new int [N_EntriesB];
  ColPtrB = new int [N_C+2];
  memset(ColPtrB, 0, sizeof(int)*(N_C+1));

  // count 
  for(int i=0;i<N_EntriesB;++i)
  {
    index = KColB[i];
    ++ColPtrB[index+2];
  }

  for(int i=1;i<N_C+1;++i)
  {
    ColPtrB[i+1] += ColPtrB[i];
  }

  // row
  for(int i=0;i<N_R;++i)
  {
    begin = RowPtrB[i];
    end = RowPtrB[i+1];

    for(int j=begin;j<end;++j)
    {
      KRowB[j] = i;
    }
  }

  // map
  for(int i=0;i<N_EntriesB;++i)
  {
    index = ColPtrB[KColB[i]+1];
    ++ColPtrB[KColB[i]+1];

    MapB[index] = i;
  }
  
// cout << "Leave GetTransposedArrays()" << endl;
// exit(0);
}
#endif

void TPardisoSolver::BenchReset()
{
  if(runs_s > 0)
  {
    time_a /= runs_a;
	time_f /= runs_f;
	time_s /= runs_s;
    *dat << num_prc << "\t" <<  N_Eq;
	*dat << "\t" << time_a << "\t" << time_f << "\t" << time_s << endl;
    runs_a = runs_f = runs_s = 0;
	time_a = time_f = time_s = 0.0;
  }
}
